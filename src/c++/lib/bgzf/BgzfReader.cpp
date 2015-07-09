/**
 ** Isaac Genome Alignment Software
 ** Copyright (c) 2010-2014 Illumina, Inc.
 ** All rights reserved.
 **
 ** This software is provided under the terms and conditions of the
 ** GNU GENERAL PUBLIC LICENSE Version 3
 **
 ** You should have received a copy of the GNU GENERAL PUBLIC LICENSE Version 3
 ** along with this program. If not, see
 ** <https://github.com/illumina/licenses/>.
 **
 ** \file BgzfReader.cpp
 **
 ** Component to read bgzf blocks.
 **
 ** \author Roman Petrovski
 **/

#include <boost/format.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>


#include "bgzf/BgzfReader.hh"

namespace isaac
{
namespace bgzf
{

inline void validateHeader(const bgzf::Header &header)
{
    ISAAC_ASSERT_MSG(header.ID1 == 31U, " got " << unsigned(header.ID1));
    ISAAC_ASSERT_MSG(header.ID2 == 139U, " got " << unsigned(header.ID2));
    ISAAC_ASSERT_MSG(header.CM == 8U, " got " << header.CM);
    ISAAC_ASSERT_MSG(header.xfield.XLEN[0] == 6U, " got " << unsigned(header.xfield.XLEN[0]));
    ISAAC_ASSERT_MSG(header.xfield.XLEN[1] == 0U, " got " << unsigned(header.xfield.XLEN[1]));
    ISAAC_ASSERT_MSG(header.xfield.SI1 == 66U, " got " << unsigned(header.xfield.SI1));
    ISAAC_ASSERT_MSG(header.xfield.SI2 == 67U, " got " << unsigned(header.xfield.SI2));
}

inline bool isValidHeader(const bgzf::Header &header)
{
    return
        header.ID1 == 31U && header.ID2 == 139U &&
        header.CM == 8U &&
        header.xfield.XLEN[0] == 6U && header.xfield.XLEN[1] == 0U &&
        header.xfield.SI1 == 66U && header.xfield.SI2 == 67U;
}

unsigned BgzfReader::readNextBlock(std::istream &is)
{
    compressedBlockBuffer_.clear();

    unsigned ret = 0;
    for (unsigned i = 0; i < blocksAtOnce_; ++i)
    {
        compressedBlockBuffer_.resize(compressedBlockBuffer_.size() + sizeof(bgzf::Header));
        is.read(&compressedBlockBuffer_.front() + compressedBlockBuffer_.size() - sizeof(bgzf::Header), sizeof(bgzf::Header));
        if (is.eof())
        {
            break;
        }
        if (!is)
        {
            BOOST_THROW_EXCEPTION(common::IoException(errno, (boost::format("Failed to read bgzf file header: %s") %
                strerror(errno)).str()));
        }
        const bgzf::Header &header = *reinterpret_cast<bgzf::Header *>(
            &compressedBlockBuffer_.front() + compressedBlockBuffer_.size() - sizeof(bgzf::Header));
        validateHeader(header);

        compressedBlockBuffer_.resize(compressedBlockBuffer_.size() + header.getCDATASize() + sizeof(bgzf::Footer));
        if (!is.read(&compressedBlockBuffer_.front() + compressedBlockBuffer_.size() - header.getCDATASize() - sizeof(bgzf::Footer), header.getCDATASize() + sizeof(bgzf::Footer)))
        {
            BOOST_THROW_EXCEPTION(common::IoException(errno, (boost::format("Failed to read %d bytes of bgzf CDATA : %s") %
                header.getCDATASize() % strerror(errno)).str()));
        }
    //    ISAAC_THREAD_CERR << "header.getCDATASize() " << header.getCDATASize() << " bytes" << std::endl;

        if (header.getCDATASize())
        {
            const bgzf::Footer &footer = *reinterpret_cast<bgzf::Footer*>(
                &compressedBlockBuffer_.front() + compressedBlockBuffer_.size() - sizeof(bgzf::Footer));
    //        ISAAC_ASSERT_MSG(footer.getISIZE(), "Unexpected ISIZE 0 for CDATA size " << header.getCDATASize());
            ret += footer.getISIZE();
        }
    }

    return ret;
}

/**
 * \brief verifies if the stream is bgzf-compressed by reading the header at the current position. Restores the position
 *        before returning.
 */
bool BgzfReader::isBgzfCompressed(std::istream &is)
{
    const std::streampos currentPos = is.tellg();
    bgzf::Header header;
    is.read(reinterpret_cast<char *>(&header), sizeof(header));
    if (is.eof())
    {
        return false;
    }
    if (!is)
    {
        BOOST_THROW_EXCEPTION(common::IoException(errno, (boost::format("Failed to read bgzf file header: %s") %
            strerror(errno)).str()));
    }
    is.seekg(currentPos);
    return isValidHeader(header);
}


void BgzfReader::uncompressCurrentBlock(char* p, std::size_t size)
{
    strm_.next_in = reinterpret_cast<Bytef *>(&compressedBlockBuffer_.front());
    strm_.avail_in = compressedBlockBuffer_.size();
    while (size)
    {
        reset();
        strm_.next_out = reinterpret_cast<Bytef *>(p);
        strm_.avail_out = size;

        int err = inflate(&strm_, Z_SYNC_FLUSH);
        if (Z_OK != err && Z_STREAM_END != err)
        {
            BOOST_THROW_EXCEPTION(BgzfInflateException(err, strm_));
        }
        const std::size_t decompressedBytes = size - strm_.avail_out;

//        if (decompressedBytes != size)
//        {
//            BOOST_THROW_EXCEPTION(common::IoException(EINVAL, (boost::format("Unexpected number of BGZF bytes uncompressed. "
//                "Expected %d, uncompressed %d") % size % decompressedBytes).str()));
//        }
        size -= decompressedBytes;
        p+= decompressedBytes;
    }
}

void ParallelBgzfReader::open(const boost::filesystem::path &filePath)
{
    fileBuffer_.reopen(filePath.c_str(), io::FileBufWithReopen::SequentialOnce);
    if (!fileBuffer_.is_open())
    {
        BOOST_THROW_EXCEPTION(common::IoException(errno, (boost::format("Failed to open bgzf file: %s") % filePath).str()));
    }
    pendingBlockSize_ = 0;
    is_.rdbuf(&fileBuffer_);
    ISAAC_THREAD_CERR << "Opened bgzf stream on " << filePath << std::endl;
}

void ParallelBgzfReader::waitForLoadSlot(boost::unique_lock<boost::mutex> &lock)
{
    while (!loadSlotAvailable_)
    {
        stateChangedCondition_.wait(lock);
    }
    loadSlotAvailable_ = false;
}

void ParallelBgzfReader::releaseLoadSlot()
{
    ISAAC_ASSERT_MSG(!loadSlotAvailable_, "Invalid attempt to release a load slot that is not locked");
    loadSlotAvailable_ = true;
    stateChangedCondition_.notify_all();
}

void ParallelBgzfReader::waitForComputeSlot(boost::unique_lock<boost::mutex> &lock)
{
    while (!computeSlotsAvailable_)
    {
        stateChangedCondition_.wait(lock);
    }
    --computeSlotsAvailable_;
}

void ParallelBgzfReader::releaseComputeSlot()
{
    ++computeSlotsAvailable_;
    stateChangedCondition_.notify_all();
}

//void ParallelBgzfReader::readMoreDataParallel(const unsigned threadNumber, std::vector<char> &buffer)
void ParallelBgzfReader::readMoreDataParallel(const unsigned threadNumber, std::istream &is, char *buffer, const std::size_t bufferMax, std::size_t &decompressed)
{
    BgzfReader &ourThreadReader = readers_.at(threadNumber);
    unsigned long &ourThreadOffset = threadOffsets_.at(threadNumber);
    unsigned long &ourThreadBlockSize = threadBlockSizes_.at(threadNumber);
    boost::unique_lock<boost::mutex> lock(stateMutex_);
    while (true)
    {

        // first get rid of any available data that is pending delivery
        if (-1UL != ourThreadOffset)
        {
            ISAAC_ASSERT_MSG(decompressed >= ourThreadOffset + ourThreadBlockSize,
                             "Result buffer is unexpectedly insufficient to fit the uncompressed block."
                             " ourThreadOffset: " << ourThreadOffset <<
                             " ourThreadBlockSize" << ourThreadBlockSize <<
                             " buffer.size()" << decompressed);

            waitForComputeSlot(lock);
            ISAAC_BLOCK_WITH_CLENAUP(boost::bind(&ParallelBgzfReader::releaseComputeSlot, this))
            {
                common::unlock_guard<boost::unique_lock<boost::mutex> > unlock(lock);
                ourThreadReader.uncompressCurrentBlock(buffer + ourThreadOffset, ourThreadBlockSize);
            }
            ourThreadOffset = -1UL;
        }

        if (is.eof())
        {
            ISAAC_THREAD_CERR << "Thread " << threadNumber << " terminating due to eof" << " buffer.size()" << decompressed << std::endl;
            break;
        }

        // Load next bgzf block
        waitForLoadSlot(lock);
        ourThreadBlockSize = 0;
        ISAAC_BLOCK_WITH_CLENAUP(boost::bind(&ParallelBgzfReader::releaseLoadSlot, this))
        {
            if (nextUncompressedOffset_ < bufferMax)
            {
                {
                    common::unlock_guard<boost::unique_lock<boost::mutex> > unlock(lock);
                    while (!is.eof() && !ourThreadBlockSize)
                    {
                        ourThreadBlockSize = ourThreadReader.readNextBlock(is);
                    }
                }
                if (!ourThreadBlockSize)
                {
                    ISAAC_ASSERT_MSG(is.eof(), "Unexpectedly stopped reading bgzf before the end of the input stream");
                    ISAAC_THREAD_CERR << "Thread " << threadNumber << " reached eof while reading empty blocks" << std::endl;
                    break;
                }
                else if (is.eof())
                {
                    ISAAC_THREAD_CERR << "Thread " << threadNumber << " reached eof" << std::endl;
                }
//                ISAAC_THREAD_CERR << "Thread " << threadNumber << " read " << ourThreadBlockSize << " bytes" << std::endl;
                ourThreadOffset = nextUncompressedOffset_;
                nextUncompressedOffset_ += ourThreadBlockSize;

                if (nextUncompressedOffset_ <= bufferMax)
                {
                    decompressed = nextUncompressedOffset_;
                    pendingBlockSize_ = 0;
                }
                else
                {
                    ourThreadOffset = 0;
                    pendingBlockSize_ = ourThreadBlockSize;
//                    ISAAC_THREAD_CERR << "Thread " << threadNumber << " has no room to place " << ourThreadBlockSize << " bytes in buffer of size " << bufferMax << std::endl;
                    break;
                }
            }
            else
            {
//                ISAAC_THREAD_CERR << "Thread " << threadNumber <<
//                    " nextUncompressedOffset_ " << nextUncompressedOffset_ <<
//                    " buffer.capacity()" << buffer.capacity() << std::endl;
                break;
            }
        }
    }

}

bool ParallelBgzfReader::readMoreData(std::vector<char> &buffer)
{
    const std::size_t oldSize = buffer.size();
    buffer.resize(buffer.capacity());
    const std::size_t ret = readMoreData(&buffer.front() + oldSize, buffer.size() - oldSize);
    buffer.resize(oldSize + ret);

    return ret;
}

std::size_t ParallelBgzfReader::readMoreData(char *buffer, const std::size_t capacity)
{
    return readMoreData(is_, buffer, capacity);
}

std::size_t ParallelBgzfReader::readMoreData(std::istream &is, char *buffer, const std::size_t capacity)
{
    std::vector<unsigned long>::iterator pendingBlockOffsetIt =
        std::find(threadOffsets_.begin(), threadOffsets_.end(), 0UL);
    if (threadOffsets_.end() != pendingBlockOffsetIt)
    {
        // preserve existing buffer.size() as there is some data the client wants to merge with what we're about
        // to decompress
        *pendingBlockOffsetIt = 0;
    }
    nextUncompressedOffset_ = pendingBlockSize_;
    pendingBlockSize_ = 0;
    std::size_t newSize = nextUncompressedOffset_;

    threads_.execute(boost::bind(&ParallelBgzfReader::readMoreDataParallel, this, _1,
                                 boost::ref(is), buffer, capacity, boost::ref(newSize)), coresMax_);
    if (threadOffsets_.end() != pendingBlockOffsetIt && !newSize)
    {
        BOOST_THROW_EXCEPTION(common::IoException(errno, (boost::format("Insufficient buffer capacity to uncompress even one bgzf block. pending %d, capacity: %d") %
            pendingBlockSize_ % capacity).str()));
    }
    return newSize;
}

void ParallelBgzfReader::initializeReaderThread(const int threadNumber)
{
    readers_.at(threadNumber).reserveBuffers();
}

} // namespace bgzf
} // namespace isaac
