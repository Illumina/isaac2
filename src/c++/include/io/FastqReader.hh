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
 ** \file FastqReader.hh
 **
 ** Component to read FASTQ files.
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_IO_FASTQ_READER_HH
#define iSAAC_IO_FASTQ_READER_HH

#include <vector>
#include <iostream>
#include <fstream>
#include <boost/filesystem.hpp>
#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <boost/ptr_container/ptr_vector.hpp>

#include "../common/StaticVector.hh"
#include "bgzf/BgzfReader.hh"
#include "common/Debug.hh"
#include "common/Exceptions.hh"
#include "flowcell/ReadMetadata.hh"
#include "io/InflateGzipDecompressor.hh"
#include "io/FileBufCache.hh"
#include "oligo/Nucleotides.hh"

namespace isaac
{
namespace io
{

/**
 ** \brief Exception thrown when fastq violation is encountered.
 **
 **/
class FastqFormatException: public common::IsaacException
{
public:
    FastqFormatException(const std::string &message):
        common::IsaacException(EINVAL, message){}
};

class FastqReader
{
    const std::size_t uncompressedBufferSize_;
    const bool allowVariableLength_;
    char q0Base_;

    FileBufWithReopen fileBuffer_;
    typedef std::vector<char> BufferType;
    io::InflateGzipDecompressor<BufferType> gzReader_;
    bgzf::ParallelBgzfReader bgzfReader_;

    //boost::filesystem::path forces intermediate string construction during reassignment...
    std::string fastqPath_;
    bool compressed_;
    bool bgzfCompressed_;
    bool reachedEof_;
    std::size_t filePos_;

    BufferType buffer_;
    BufferType::const_iterator headerBegin_;
    BufferType::const_iterator headerEnd_;
    BufferType::const_iterator baseCallsBegin_;
    BufferType::const_iterator baseCallsEnd_;
    BufferType::const_iterator qScoresBegin_;
    BufferType::const_iterator endIt_;
    bool zeroLengthRead_;

    static const oligo::Translator translator_;

public:
    static const unsigned INCORRECT_FASTQ_BASE = 5;
    // experimentally determined to provide reasonable balance between thread stop/start overhead and RAM requirements
    // give each thread a chance to unpack a sensible number of bgzf blocks. Otherwise thread synchronization
    // slows the whole thing down
    static const unsigned BGZF_BLOCKS_PER_THREAD = 1024;


    FastqReader(const bool allowVariableLength, const unsigned threadsMax, const std::size_t maxPathLength);

    void open(const boost::filesystem::path &fastqPath, const char q0Base);

    void next();

    template <typename InsertIt>
    InsertIt extractBcl(const flowcell::ReadMetadata &readMetadata, InsertIt it) const;

    template <typename InsertIt>
    InsertIt extractReadName(const unsigned nameLengthMax, InsertIt it) const;

    const std::string &getPath() const
    {
        return fastqPath_;
    }

    std::size_t getRecordOffset() const
    {
        return getOffset(headerBegin_);
    }

    bool hasData() const
    {
        return (!reachedEof_ || buffer_.end() != headerBegin_);
    }

    typedef std::pair<BufferType::const_iterator, BufferType::const_iterator > IteratorPair;
    IteratorPair getHeader() const
    {
        return std::make_pair(headerBegin_, headerEnd_);
    }

    unsigned getReadLength() const
    {
        return std::distance(baseCallsBegin_, baseCallsEnd_);
    }

private:
    typedef boost::error_info<struct tag_errmsg, std::string> errmsg_info;

    void resetBuffer();
    std::size_t getOffset(BufferType::const_iterator it) const;
    void findHeader();
    void findSequence();
    void findQScores();
    void findQScoresEnd();
    bool fetchMore();

    std::size_t readCompressedFastq(std::istream &is, char *buffer, std::size_t amount);
    std::size_t readBgzfFastq(std::istream &is, char *buffer, std::size_t amount);
    std::size_t readFlatFastq(std::istream &is, char *buffer, std::size_t amount);
};

template <typename InsertIt>
InsertIt FastqReader::extractBcl(const flowcell::ReadMetadata &readMetadata, InsertIt it) const
{
    const InsertIt start = it;
    BufferType::const_iterator baseCallsIt = baseCallsBegin_;
    BufferType::const_iterator qScoresIt = qScoresBegin_;
    std::vector<unsigned>::const_iterator cycleIterator = readMetadata.getCycles().begin();
    unsigned currentCycle = readMetadata.getFirstReadCycle();
    for(;endIt_ != qScoresIt && readMetadata.getCycles().end() != cycleIterator; ++baseCallsIt, ++qScoresIt, ++currentCycle)
    {
//        ISAAC_THREAD_CERR << "cycle " << *cycleIterator << std::endl;
        if (*cycleIterator != currentCycle)
        {
            continue;
        }
        const unsigned char baseValue = translator_[*baseCallsIt];
        if (oligo::INVALID_OLIGO == baseValue)
        {
            *it = 0;
        }
        else if(INCORRECT_FASTQ_BASE == baseValue)
        {
            BOOST_THROW_EXCEPTION(FastqFormatException((boost::format("Invalid oligo %c found in %s at offset %u") %
                *baseCallsIt % getPath() % getOffset(baseCallsIt)).str()));
        }
        else
        {
            const unsigned char baseQuality = (*qScoresIt - q0Base_);
            if ((1 << 6) <= baseQuality)
            {
                BOOST_THROW_EXCEPTION(FastqFormatException((boost::format("Invalid quality %d found in %s at offset %u. Base quality scores [0-63] supported only. baseq0=%d") %
                    unsigned(baseQuality) % getPath() % getOffset(baseCallsIt) % q0Base_).str()));
            }
            *it = baseValue | (baseQuality << 2);
        }
//        ISAAC_THREAD_CERR << "stored cycle " << *cycleIterator << std::endl;
        ++it;
        ++cycleIterator;
    }

    const unsigned extracted = std::distance<std::vector<unsigned>::const_iterator>(readMetadata.getCycles().begin(), cycleIterator);
    if (!allowVariableLength_)
    {
        if (readMetadata.getCycles().size() != extracted)
        {
    //            ISAAC_ASSERT_MSG(false, "Fastq read shorter than expected." << ret);
            BOOST_THROW_EXCEPTION(common::IoException(EINVAL, (boost::format("Read length (%d) "
                " is different from expected %d in %s:%u. Record %s") %
                extracted % readMetadata.getCycles().size() %
                getPath() % getRecordOffset() %
                std::string(headerBegin_, endIt_)).str()));
        }
        else
        {
//            ISAAC_THREAD_CERR << "length ok " << std::endl;30=.cilprw
        }
    }
    else if (extracted < readMetadata.getCycles().size())
    {
        it = std::fill_n(it, readMetadata.getCycles().size() - extracted, 0);
    }

    ISAAC_ASSERT_MSG(readMetadata.getCycles().size() == std::size_t(std::distance(start, it)),
        "unexpected number of cycles read: " << std::distance(start, it) << " expected: " << readMetadata);

    return it;
}

/**
 * \brief retrieve name of the fastq read to the first whitespace. Pad up to nameLengthMax with 0
 *
 * \return it + nameLengthMax
 */
template <typename InsertIt>
InsertIt FastqReader::extractReadName(const unsigned nameLengthMax, InsertIt it) const
{
    static const char whitespace[] = {' ', '\t'};
    const std::size_t nameAvailable = std::distance(
        headerBegin_+1, std::find_first_of(headerBegin_+1, headerEnd_, whitespace, whitespace + sizeof(whitespace)));
    if (nameLengthMax > nameAvailable)
    {
        it = std::copy(headerBegin_+1, headerBegin_ + 1 + nameAvailable, it);
        it = std::fill_n(it, nameLengthMax - nameAvailable, 0);
    }
    else
    {
        it = std::copy(headerBegin_ + 1 + nameAvailable - nameLengthMax, headerBegin_ + 1 + nameAvailable, it);
    }
    return it;
}

} // namespace io
} // namespace isaac

#endif // #ifndef iSAAC_IO_FASTQ_READER_HH
