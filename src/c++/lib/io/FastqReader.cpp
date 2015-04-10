/**
 ** Isaac Genome Alignment Software
 ** Copyright (c) 2010-2014 Illumina, Inc.
 ** All rights reserved.
 **
 ** This software is provided under the terms and conditions of the
 ** Illumina Public License 1
 **
 ** You should have received a copy of the Illumina Public License 1
 ** along with this program. If not, see
 ** <https://github.com/sequencing/licenses/>.
 **
 ** \file FastqReader.cpp
 **
 ** Component to read FASTQ files.
 **
 ** \author Roman Petrovski
 **/
//#include <xmmintrin.h>
#include <boost/bind.hpp>

#include "common/Debug.hh"
#include "common/FileSystem.hh"
#include "io/FastqReader.hh"

namespace isaac
{
namespace io
{

const oligo::Translator FastqReader::translator_ = oligo::getTranslator(true, FastqReader::INCORRECT_FASTQ_BASE);

FastqReader::FastqReader(const bool allowVariableLength, const unsigned threadsMax, const std::size_t maxPathLength) :
    // The uncompressed buffer can be fairly small for flat and gzipped fastq, however we need a decent amount of
    // space for parallel decompression to be effective with bgzf-compressed fastq.
    uncompressedBufferSize_(std::size_t(bgzf::BgzfReader::UNCOMPRESSED_BGZF_BLOCK_SIZE) * threadsMax * BGZF_BLOCKS_PER_THREAD),
    allowVariableLength_(allowVariableLength),
    fileBuffer_(std::ios_base::in),
    gzReader_(),
    // give the threads a chance to rebalance if some of the blocks take longer to uncompress
    bgzfReader_(threadsMax, BGZF_BLOCKS_PER_THREAD / 8),
    fastqPath_(std::string(maxPathLength, ' ')),
    compressed_(false),
    bgzfCompressed_(false),
    reachedEof_(false),
    filePos_(0),
    zeroLengthRead_(false)
{
    ISAAC_THREAD_CERR << "FastqReader uncompressedBufferSize_=" << uncompressedBufferSize_ << std::endl;
    buffer_.reserve(uncompressedBufferSize_);
    resetBuffer();
}

void FastqReader::resetBuffer()
{
    buffer_.resize(uncompressedBufferSize_);
    headerBegin_ = buffer_.end();
    headerEnd_ = buffer_.end();
    baseCallsBegin_ = buffer_.end();
    baseCallsEnd_ = buffer_.end();
    qScoresBegin_ = buffer_.end();
    endIt_ = buffer_.end();
}

void FastqReader::open(const boost::filesystem::path &fastqPath)
{
    if (fastqPath.string() != fastqPath_)
    {
        resetBuffer();
        // ensure actual copying, prevent path buffer sharing
        fastqPath_ = fastqPath.c_str();
        compressed_ = common::isDotGzPath(fastqPath_);
        fileBuffer_.reopen(fastqPath_.c_str(), FileBufWithReopen::SequentialOnce);
        buffer_.resize(uncompressedBufferSize_);
        gzReader_.reset();
        filePos_ = 0;

        if (fileBuffer_.is_open())
        {
            std::istream is(&fileBuffer_);
            bgzfCompressed_ = compressed_ ? bgzf::BgzfReader::isBgzfCompressed(is) : false;
            reachedEof_ = false;
            next();
        }
        else
        {
            BOOST_THROW_EXCEPTION(common::IoException(errno, (boost::format("Failed to open fastq file: %s") %
                fastqPath_).str()));
        }
        ISAAC_THREAD_CERR << "Opened " << (bgzfCompressed_ ? "bgzf" : compressed_ ? "gz" : "") <<
            " fastq stream on " << fastqPath_ << std::endl;
    }
    else
    {
        ISAAC_THREAD_CERR << "Keeping open fastq stream on " << fastqPath_ << std::endl;
    }
}

std::size_t FastqReader::getOffset(BufferType::const_iterator it) const
{
    const BufferType::const_iterator begin = buffer_.begin();
    return filePos_ - buffer_.size() + std::distance(begin, it);
}

template <typename IteratorT>
IteratorT findNotNewLine(IteratorT itBegin, IteratorT itEnd)
{
    return std::find_if(
        itBegin, itEnd,
         boost::bind(std::not_equal_to<char>(), '\r', _1) &&
         boost::bind(std::not_equal_to<char>(), '\n', _1));
}

template <typename IteratorT>
IteratorT findNewLine(IteratorT itBegin, IteratorT itEnd)
{
    static const char *rn = "\n\r";
    return std::find_first_of(itBegin, itEnd, rn, rn+2);
}

void FastqReader::findHeader()
{
//    ISAAC_THREAD_CERR << "FastqReader::findHeader " << std::string(endIt_, endIt_ + std::min(300L, std::distance<BufferType::const_iterator>(endIt_, buffer_.end()))) << std::endl;
    headerBegin_ = findNotNewLine(endIt_, BufferType::const_iterator(buffer_.end()));
    if (buffer_.end() == headerBegin_)
    {
        // We've reached the end of the buffer before we reached the beginning of the next record
        if (!fetchMore())
        {
            return;
        }
        // TODO: this allows more than one newline which is against fastq format
        headerBegin_ = findNotNewLine(headerBegin_, BufferType::const_iterator(buffer_.end()));
        if (buffer_.end() == headerBegin_)
        {
            if (reachedEof_)
            {
                return;
            }
            BOOST_THROW_EXCEPTION(FastqFormatException((boost::format(
                "Too many newline characters in fastq to fit in the buffer: %s, offset %u") %
                getPath() % getOffset(headerBegin_)).str()));
        }
    }
    headerEnd_ = findNewLine(headerBegin_, BufferType::const_iterator(buffer_.end()));
    if (buffer_.end() == headerEnd_)
    {
        // We've reached the end of the buffer before we reached the end of the header
        if (!fetchMore())
        {
            BOOST_THROW_EXCEPTION(FastqFormatException((boost::format("Fastq file end while reading the header line: %s, offset %u") %
                getPath() % getOffset(headerEnd_)).str()));
        }
        headerEnd_ = findNewLine(headerEnd_, BufferType::const_iterator(buffer_.end()));
        if (buffer_.end() == headerEnd_)
        {
            BOOST_THROW_EXCEPTION(FastqFormatException((boost::format("Fastq header too long to fit in the buffer: %s, offset %u") %
                getPath() % getOffset(baseCallsBegin_)).str()));
        }
    }
}

void FastqReader::findSequence()
{
    baseCallsBegin_ = findNotNewLine(headerEnd_, BufferType::const_iterator(buffer_.end()));
    if (buffer_.end() == baseCallsBegin_)
    {
        // We've reached the end of the buffer before we reached the beginning of the sequence
        if (!fetchMore())
        {
            BOOST_THROW_EXCEPTION(FastqFormatException((boost::format("Fastq file end while looking for sequence start: %s, offset %u") %
                getPath() % getOffset(baseCallsBegin_)).str()));
        }
        baseCallsBegin_ = findNotNewLine(baseCallsBegin_, BufferType::const_iterator(buffer_.end()));
        if (buffer_.end() == baseCallsBegin_)
        {
            BOOST_THROW_EXCEPTION(FastqFormatException((boost::format(
                "Too many newline characters in fastq to fit in the buffer while looking for sequence start: %s, offset %u") %
                getPath() % getOffset(baseCallsBegin_)).str()));
        }
    }

    // special case for zero-length reads
    if ('+' == *baseCallsBegin_)
    {
        zeroLengthRead_ = true;
        baseCallsEnd_ = baseCallsBegin_;
    }
    else
    {
        zeroLengthRead_ = false;
        baseCallsEnd_ = findNewLine(baseCallsBegin_, BufferType::const_iterator(buffer_.end()));
    }
    if (buffer_.end() == baseCallsEnd_)
    {
        // We've reached the end of the buffer before we reached the end of the sequence
        if (!fetchMore())
        {
            BOOST_THROW_EXCEPTION(FastqFormatException((boost::format("Fastq file end while reading the sequence line: %s, offset %u") %
                getPath() % getOffset(baseCallsEnd_)).str()));
        }
        baseCallsEnd_ = findNewLine(baseCallsEnd_, BufferType::const_iterator(buffer_.end()));
        if (buffer_.end() == baseCallsEnd_)
        {
//            ISAAC_THREAD_CERR << " findSequence " << std::string(buffer_.begin(), buffer_.end()) << " buffer_.size()=" << buffer_.size() << std::endl;
            BOOST_THROW_EXCEPTION(FastqFormatException((boost::format("Fastq sequence too long to fit in the buffer: %s, offset %u") %
                getPath() % getOffset(baseCallsEnd_)).str()));
        }
    }
}

void FastqReader::findQScores()
{
    qScoresBegin_ = findNotNewLine(baseCallsEnd_, BufferType::const_iterator(buffer_.end()));
    if (buffer_.end() == qScoresBegin_)
    {
        // We've reached the end of the buffer before we reached the beginning of the sequence
        if (!fetchMore())
        {
            BOOST_THROW_EXCEPTION(FastqFormatException((boost::format("Fastq file end while looking for + sign: %s, offset %u") %
                getPath() % getOffset(qScoresBegin_)).str()));
        }
        qScoresBegin_ = findNotNewLine(qScoresBegin_, BufferType::const_iterator(buffer_.end()));
        if (buffer_.end() == qScoresBegin_)
        {
            BOOST_THROW_EXCEPTION(FastqFormatException((boost::format(
                "Too many newline characters in fastq to fit in the buffer while looking for + sign: %s, offset %u") %
                getPath() % getOffset(qScoresBegin_)).str()));
        }
    }

    if ('+' != *qScoresBegin_)
    {
        BOOST_THROW_EXCEPTION(FastqFormatException((boost::format(
            "+ sign not found where expected: %s, offset %u") %
            getPath() % getOffset(qScoresBegin_)).str()));
    }
    // in some fastq files (like the ones produced by sra tools) + is followed by the header string. Just skip to the newline...
    qScoresBegin_ = findNewLine(qScoresBegin_, BufferType::const_iterator(buffer_.end()));

    qScoresBegin_ = findNotNewLine(qScoresBegin_, BufferType::const_iterator(buffer_.end()));
    if (buffer_.end() == qScoresBegin_)
    {
        // We've reached the end of the buffer before we reached the beginning of the qscores
        if (!fetchMore())
        {
            BOOST_THROW_EXCEPTION(FastqFormatException((boost::format("Fastq file end while looking for qscores: %s, offset %u") %
                getPath() % getOffset(qScoresBegin_)).str()));
        }
        qScoresBegin_ = findNotNewLine(qScoresBegin_, BufferType::const_iterator(buffer_.end()));
        if (buffer_.end() == qScoresBegin_)
        {
            BOOST_THROW_EXCEPTION(FastqFormatException((boost::format(
                "Too many newline characters in fastq to fit in the buffer while looking for qscores: %s, offset %u") %
                getPath() % getOffset(qScoresBegin_)).str()));
        }
    }
}

void FastqReader::findQScoresEnd()
{
    if (zeroLengthRead_)
    {
        endIt_ = qScoresBegin_;
    }
    else
    {
        endIt_ = findNewLine(qScoresBegin_, BufferType::const_iterator(buffer_.end()));
        if (buffer_.end() == endIt_)
        {
            // We've reached the end of the buffer before we reached the newline...
            if (!fetchMore())
            {
                return;
            }
            endIt_ = findNewLine(endIt_, BufferType::const_iterator(buffer_.end()));
        }
    }
}

std::size_t FastqReader::readCompressedFastq(std::istream &is, char *buffer, std::size_t amount)
{
    const std::size_t decompressedBytes = gzReader_.read(is, 0, buffer, amount);
    reachedEof_ = is.eof();
    return decompressedBytes;

//    const std::size_t amountOri = amount;
//    while (amount)
//    {
//        const std::streamsize decompressedBytes =
//            decompressor_.read(is, buffer, amount);
//
//        if (-1 == decompressedBytes)
//        {
////            ISAAC_THREAD_CERR << "FastqReader::readCompressedFastq std::size_t(-1) == decompressedBytes still need: " << amount << std::endl;
///*
//            if (!is.eof())
//            {
//                BOOST_THROW_EXCEPTION(common::IoException(errno, (boost::format(
//                    "readCompressedFastq failed: %s, uncompressed offset %u") %
//                    getPath() % getOffset(headerBegin_)).str()));
//            }
//*/
//            reachedEof_ = true;
//            return amountOri - amount;
//        }
//        amount -= decompressedBytes;
//        buffer += decompressedBytes;
////        ISAAC_THREAD_CERR << "FastqReader::readCompressedFastq read " << decompressedBytes <<
////            " decompressedBytes, still need: " << amount << std::endl;
//    }
//
//    return amountOri;
}

std::size_t FastqReader::readBgzfFastq(std::istream &is, char *buffer, std::size_t amount)
{
    const std::size_t ret = bgzfReader_.readMoreData(is, buffer, amount);
    reachedEof_ = is.eof();
    return ret;
}

std::size_t FastqReader::readFlatFastq(std::istream &is, char *buffer, std::size_t amount)
{
    is.read(buffer, amount);
    if (!is.good() && !is.eof())
    {
        BOOST_THROW_EXCEPTION(common::IoException(errno, (boost::format(
            "readFlatFastq failed: %s, offset %u") %
            getPath() % getOffset(headerBegin_)).str()));
    }
    reachedEof_ = is.eof();
    return is.gcount();

}

bool FastqReader::fetchMore()
{
    if (reachedEof_)
    {
        return false;
    }

    // move the remaining data to the start of the buffer
    std::copy(headerBegin_, BufferType::const_iterator(buffer_.end()),  buffer_.begin());
    const std::size_t moved = std::distance(headerBegin_, BufferType::const_iterator(buffer_.end()));
//    ISAAC_THREAD_CERR << "fetchMore moved=" << moved << " in buffer of size " << buffer_.size() << std::endl;
    const std::size_t distance = std::distance(BufferType::const_iterator(buffer_.begin()), headerBegin_);
    if (!distance)
    {
        BOOST_THROW_EXCEPTION(common::IoException(EINVAL, (boost::format(
            "Fastq fetchMore impossible. No more space in the buffer: %s, offset %u") %
            getPath() % getOffset(headerBegin_)).str()));
    }

    headerBegin_ -= distance;
    headerEnd_ -= distance;
    baseCallsBegin_ -= distance;
    baseCallsEnd_ -= distance;
    qScoresBegin_ -= distance;
    endIt_ -= distance;

    BufferType::iterator firstUnreadByte = buffer_.end() - distance;


    try
    {
        std::istream is(&fileBuffer_);
        buffer_.resize(uncompressedBufferSize_);
        const std::size_t availableSpace = buffer_.size() - moved;
        const std::size_t readBytes = bgzfCompressed_ ?
            readBgzfFastq(is, &*firstUnreadByte, availableSpace) :
            compressed_ ?
            readCompressedFastq(is, &*firstUnreadByte, availableSpace) :
            readFlatFastq(is, &*firstUnreadByte, availableSpace);

        filePos_ += readBytes;
        buffer_.resize(moved + readBytes);
    }
    catch (boost::exception &e)
    {
        e << errmsg_info(" While reading from " + getPath());
        throw;
    }
    return true;
}

void FastqReader::next()
{
    findHeader();
    if (buffer_.end() == headerBegin_)
    {
        return;
    }
    findSequence();
    findQScores();
    findQScoresEnd();
}


} // namespace io
} // namespace isaac
