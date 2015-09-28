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
 ** \file Cigar.hh
 **
 ** \brief Tools for creation, handling, management of BAM CIGAR
 ** 
 ** \author Come Raczy
 **/

#ifndef iSAAC_ALIGNMENT_CIGAR_HH
#define iSAAC_ALIGNMENT_CIGAR_HH

#include <string>
#include <vector>
#include <utility>
#include <algorithm>
#include <numeric>
#include <stdint.h>
#include <boost/bind.hpp>
#include <boost/foreach.hpp>

#include "common/Debug.hh"
#include "common/FastIo.hh"
#include "flowcell/Layout.hh"
#include "flowcell/ReadMetadata.hh"
#include "reference/ReferencePosition.hh"

namespace isaac
{
namespace alignment
{

class Cigar: std::vector<uint32_t>
{
    typedef std::vector<uint32_t> BaseT;
public:
    Cigar (const std::size_t n, const uint32_t &init) : std::vector<uint32_t>(n, init)
    {
    }

    Cigar(const std::vector<uint32_t> &that) : std::vector<uint32_t>(that)
    {
    }

    Cigar(unsigned reservedSize = 0)
    {
        reserve(reservedSize);
    }
    enum OpCode {
        ALIGN = 0, // 'M'
        INSERT = 1, // 'I'
        DELETE = 2, // 'D'
        SKIP = 3, // 'N' Not being used in iSAAC.
                 // Use DELETE whenever there is a jump in the reference between two bases of a read.
        SOFT_CLIP = 4, // 'S'
        HARD_CLIP = 5, // 'H'
        PAD = 6, // 'P'
        MATCH = 7, // '='
        MISMATCH = 8, // 'X'

        CONTIG = 9, // Absolute contig change. No equivalent in Bam CIGAR.
                    // Simply changes the contig while leaving the position unchanged.
                    // Normally followed by DELETE or BACK for position adjustment.
        BACK = 10,  // same as DELETE, but backwards. No equivalent in Bam CIGAR.
        FLIP = 11,  // Changes the orientation from F to R and vice versa.
                    // The argument indicates the number of sequence bases covered by subsequent CIGAR
                    // No equivalent in Bam CIGAR.
        UNKNOWN = 12 // '?'
    };

    using BaseT::at;
    using BaseT::back;
    using BaseT::begin;
    using BaseT::capacity;
    using BaseT::clear;
    using BaseT::iterator;
    using BaseT::const_iterator;
    using BaseT::empty;
    using BaseT::end;
    using BaseT::erase;
    using BaseT::front;
    using BaseT::pop_back;
    using BaseT::reserve;
    using BaseT::resize;
    using BaseT::size;
    using BaseT::value_type;
    using BaseT::operator[];

    void reserve(size_type n)
    {
        BaseT::reserve(n);
    }

    // we're in read length of hundreds. assume read length of thousands plus the operation code char
    static const unsigned OPERATION_CHARS_MAX = 5;
    template <typename IteratorT>
    void addOperations(IteratorT b, IteratorT e)
    {
        // reallocation of Cigar buffer can cause lots of bad memory access. Just prohibit it
        ISAAC_ASSERT_MSG(capacity() >= size() + std::distance(b, e), "CIGAR buffer is out of capacity:" << capacity() << " needed:" << std::distance(b, e) << " more");
        insert(end(), b, e);
    }

    void addOperation(long length, OpCode opCode)
    {
        unsigned long ulength = std::abs(length);
        if (DELETE == opCode && 0 > length)
        {
            opCode = BACK;
        }
        else
        {
            ISAAC_ASSERT_MSG(0 <= length, "Negative length is only allowed for DELETE. opCode:" << opCodeToChar(opCode) << " length:" << length);
        }

        do
        {
            // reallocation of Cigar buffer can cause lots of bad memory access. Just prohibit it
            ISAAC_ASSERT_MSG(capacity() >= size() + 1, "CIGAR buffer is out of capacity:" << capacity());
            push_back(encode(ulength, opCode));
        } while (ulength);
    }

    void updateOperation(const unsigned offset, unsigned long len, const Cigar::OpCode op)
    {
        at(offset) = Cigar::encode(len, op);
        ISAAC_ASSERT_MSG(!len, "Expected the length to fit in one component");
    }

    std::string toString() const;
    std::string toString(unsigned offset, unsigned length) const;
    static std::string toString(const Cigar &cigarBuffer, unsigned offset, unsigned length);
    static char opCodeToChar(const Cigar::OpCode opCode)
    {
        static const char opCodes[] = {'M','I','D','N','S','H','P','=','X','C','B','F','?'};
        ISAAC_ASSERT_MSG(sizeof(opCodes) > std::size_t(opCode), "Unexpected CIGAR op code: " << opCode);
        return opCodes[opCode];
    }
    /**
     * \brief Serializes cigar to a container. Does not push terminating zero
     *
     * \return const reference to result
     */
    template <typename IteratorT, typename ContainerT>
    static const ContainerT &toString(IteratorT begin, IteratorT end, ContainerT &result)
    {
        for (IteratorT v = begin; end != v; ++v)
        {
            const Component d = Cigar::decode(*v);
            common::appendUnsignedInteger(result, d.first);
            result.push_back(opCodeToChar(d.second));
        }
        return result;
    }

    template <typename IteratorT>
    static std::string toString(IteratorT begin, IteratorT end)
    {
        std::string result;
        return toString(begin, end, result);
    }

    template <typename IteratorT>
    static std::ostream& toStream(IteratorT begin, IteratorT end, std::ostream &os)
    {
        for (IteratorT v = begin; end != v; ++v)
        {
            const Component d = Cigar::decode(*v);
            os << d.first << opCodeToChar(d.second);
        }
        return os;
    }

    template <typename IteratorT>
    static unsigned getReadLength(IteratorT begin, IteratorT end)
    {
        unsigned ret = 0;
        BOOST_FOREACH(const value_type v, std::make_pair(begin, end))
        {
            const Component d = Cigar::decode(v);
            switch (d.second)
            {
            case ALIGN:
            case INSERT:
            case SOFT_CLIP:
                ret += d.first;
                break;
            default:
                break;
            }
        }
        return ret;
    }

    template <typename IteratorT>
    static unsigned getMappedLength(IteratorT begin, IteratorT end)
    {
        unsigned ret = 0;
        BOOST_FOREACH(const value_type v, std::make_pair(begin, end))
        {
            const Component d = Cigar::decode(v);
            switch (d.second)
            {
            case ALIGN:
                ret += d.first;
                break;
            default:
                break;
            }
        }
        return ret;
    }

    typedef std::pair<uint32_t, OpCode> Component;

    struct Encoder
    {
        static const uint32_t LENGTH_MAX = ~(~uint32_t(0) << 28);
        Encoder(unsigned long length, OpCode opCode):
            opCode_(opCode), length_(LENGTH_MAX < length ? LENGTH_MAX : length)
        {
            ISAAC_ASSERT_MSG(ALIGN <= opCode && UNKNOWN >= opCode, "Invalid CIGAR code " << opCode);
//            ISAAC_ASSERT_MSG(length == length_, "Supplied length value does not fit the length_ data field: " << length << " length_:" << length_);
        }

        Encoder(uint32_t value) : value_(value)
        {
        }

        union
        {
            struct
            {
                OpCode opCode_: 4;
                uint32_t length_ : 28;
            };
            uint32_t value_;
        };

        uint32_t getValue() const {return value_;}
        uint32_t getLength() const {return length_;}
        Component unpack() const {return Component(Component::first_type(length_), Component::second_type(opCode_));}
    };

    static Component decode(const uint32_t value)
    {
        return Encoder(value).unpack();
//        const unsigned length = value>>4;
//        const unsigned code = std::min(value & 0x0F, static_cast<unsigned>(UNKNOWN));
//        return std::pair<unsigned, OpCode>(length, static_cast<OpCode>(code));
    }


    static unsigned getMaxLength(const unsigned readLength)
    {
        return getMaxOpeations(readLength) * sizeof(value_type);
    }

    static unsigned getMinLength()
    {
        return 1 * sizeof(value_type);
    }

    static unsigned getMaxOpeations(const unsigned readLength)
    {
        const unsigned minBasesToFindAnIndel = 10; //assuming at least 10 bases are needed to identify an indel
        const unsigned maxCigarIndels = (readLength / minBasesToFindAnIndel);
        const unsigned cigarOpsPerIndel = 2; //assuming one indel requries one match
        const unsigned oneMatchOp = 1;
        const unsigned maxHardClipOps = 2; //one at each end
        const unsigned maxSoftClipOps = 2; //one at each end
        const unsigned maxCigarOperations =
            maxSoftClipOps + maxHardClipOps + oneMatchOp + maxCigarIndels * cigarOpsPerIndel;

        return maxCigarOperations;
    }

    static unsigned getMaxOpeationsForRead(const flowcell::ReadMetadata &readMetadata)
    {
        return getMaxLength(readMetadata.getLength());
    }

    static unsigned getMaxOperationsForReads(const flowcell::ReadMetadataList &readMetadataList)
    {
        return std::accumulate(readMetadataList.begin(), readMetadataList.end(), 0,
                               boost::bind(std::plus<unsigned>(), _1,
                                           boost::bind(&Cigar::getMaxOpeationsForRead, _2)));
    }

    static unsigned getMaxOperationsForReads(const flowcell::FlowcellLayoutList &flowcellLayoutList)
    {
        unsigned ret = 0;
        BOOST_FOREACH(const flowcell::Layout &flowcell, flowcellLayoutList)
        {
            ret = std::max(ret, getMaxOperationsForReads(flowcell.getReadMetadataList()));
        }
        return ret;
    }
private:
    static uint32_t encode(unsigned long &length, OpCode opCode)
    {
        Encoder c(length, opCode);
        length -= c.getLength();
        return c.getValue();
//        assert(ALIGN <= opCode);
//        assert(UNKNOWN >= opCode);
//        return (length<<4) | opCode;
    }
};

inline std::ostream &operator <<(std::ostream &os, const Cigar::Component &cigarComponent)
{
    return os << "CigarComponent(" << cigarComponent.first << "," << cigarComponent.second << ")";
}


template <typename IteratorT>
struct CigarPosition
{
    CigarPosition(const IteratorT cigarIt, const IteratorT cigarEnd,
        const reference::ReferencePosition referencePos, const bool reverse,
        const unsigned readLength):
        referencePos_(referencePos), cigarIt_(cigarIt), cigarEnd_(cigarEnd),
        sequenceOffset_(0), reverse_(reverse), readLength_(readLength){}

    reference::ReferencePosition referencePos_;
    IteratorT cigarIt_;
    IteratorT cigarEnd_;
    unsigned sequenceOffset_;
    bool reverse_;
    unsigned readLength_;

    bool operator ==(const CigarPosition &other) const
    {
        ISAAC_ASSERT_MSG(cigarEnd_ == other.cigarEnd_, "Attempt to compare CigarPosition objects of difference CIGARs");
        if (cigarIt_ == other.cigarIt_)
        {
            ISAAC_ASSERT_MSG(
                referencePos_ == other.referencePos_ ||
                isaac::reference::ReferencePosition(isaac::reference::ReferencePosition::NoMatch) == referencePos_ ||
                isaac::reference::ReferencePosition(isaac::reference::ReferencePosition::NoMatch) == other.referencePos_,
                "Invalid comparison of CigarPositions. When cigar iterators match, positions must match or be not set");
            return true;
        }
        return false;
    }

    bool operator != (const CigarPosition &other) const
    {
        return !(other == *this);
    }

    Cigar::Component component() const {return Cigar::decode(*cigarIt_); }

    CigarPosition & operator ++()
    {
        ISAAC_ASSERT_MSG(cigarEnd_ != cigarIt_, "Attempt to advance past the end of the CIGAR. sequenceOffset_:" << sequenceOffset_);
        const Cigar::Component component = Cigar::decode(*cigarIt_++);
        if (Cigar::ALIGN == component.second)
        {
            referencePos_ += component.first;
            sequenceOffset_ += component.first;
        }
        else if (Cigar::INSERT == component.second)
        {
            sequenceOffset_ += component.first;
        }
        else if (Cigar::DELETE == component.second)
        {
            referencePos_ += component.first;
        }
        else if (Cigar::BACK == component.second)
        {
            referencePos_ -= component.first;
            skipToNextAlign();
        }
        else if (Cigar::CONTIG == component.second)
        {
            referencePos_ = reference::ReferencePosition(component.first, referencePos_.getPosition());
            skipToNextAlign();
        }
        else if (Cigar::FLIP == component.second)
        {
            sequenceOffset_ = readLength_ - sequenceOffset_ - component.first;
            reverse_ = !reverse_;
            skipToNextAlign();
        }
        else if (Cigar::SOFT_CLIP == component.second)
        {
            sequenceOffset_ += component.first;
        }
        else
        {
            using boost::format;
            const format message = format("Unexpected Cigar OpCode: %d") % component.second;
            BOOST_THROW_EXCEPTION(common::PostConditionException(message.str()));
        }
        return *this;
    }

    bool end() const {return cigarEnd_ == cigarIt_;}
private:
    void skipToNextAlign()
    {
        while (!end())
        {
            const Cigar::Component nextComponent = Cigar::decode(*cigarIt_);
            if (Cigar::DELETE == nextComponent.second)
            {
                referencePos_ += nextComponent.first;
            }
            else if (Cigar::BACK == nextComponent.second)
            {
                referencePos_ -= nextComponent.first;
            }
            else if (Cigar::CONTIG == nextComponent.second)
            {
                referencePos_ = reference::ReferencePosition(nextComponent.first, referencePos_.getPosition());
            }
            else
            {
                ISAAC_ASSERT_MSG(
                    Cigar::ALIGN == nextComponent.second || Cigar::SOFT_CLIP == nextComponent.second,
                    "Unexpected cigar operation " << nextComponent.second << " while looking for Cigar::ALIGN or Cigar::SOFT_CLIP in " << *this);
                break;
            }
            ++cigarIt_;
        }
    }

    friend std::ostream & operator <<(std::ostream &os, const CigarPosition& pos)
    {
        return os << "CigarPosition(" <<
            pos.referencePos_ << "," <<
            alignment::Cigar::toString(pos.cigarIt_, pos.cigarEnd_) << "," <<
            pos.reverse_ <<  "," <<
            pos.sequenceOffset_ <<  "," <<
            pos.readLength_ << ")";
    }
};

template <typename IteratorT>
unsigned computeObservedLength(const IteratorT cigarBegin, const IteratorT cigarEnd)
{
//    ISAAC_THREAD_CERR << Cigar::toString(cigarBegin, cigarEnd) << std::endl;

    reference::ReferencePosition lastRefPos(0,0);
    CigarPosition<IteratorT> cp(cigarBegin, cigarEnd, lastRefPos, false, 0);
    while(!cp.end())
    {
//        ISAAC_THREAD_CERR << lastRefPos << cp.referencePos_ << std::endl;
        ISAAC_ASSERT_MSG(lastRefPos.getContigId() == cp.referencePos_.getContigId(),
                         "Only regular CIGARs are supported. lastRefPos:" << lastRefPos <<
                         " cp.referencePos_:" << cp.referencePos_ << Cigar::toString(cigarBegin, cigarEnd));
        ISAAC_ASSERT_MSG(!cp.reverse_, "Only regular CIGARs are supported. lastRefPos:" << lastRefPos <<
                         " cp.referencePos_:" << cp.referencePos_ << Cigar::toString(cigarBegin, cigarEnd));
        ISAAC_ASSERT_MSG(lastRefPos <= cp.referencePos_, "Only regular CIGARs are supported");
        lastRefPos = cp.referencePos_;
        ++cp;
    }

    return cp.referencePos_ - reference::ReferencePosition(0,0);
}

} // namespace alignment
} // namespace isaac

#endif // #ifndef iSAAC_ALIGNMENT_CIGAR_HH
