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
 ** \file ReferencePosition.hh
 **
 ** Representation of a position in a reference genome.
 **
 ** \author Come Raczy
 **/

#ifndef iSAAC_REFERENCE_REFERENCE_POSITION_HH
#define iSAAC_REFERENCE_REFERENCE_POSITION_HH

#include <iostream>

#include "common/Debug.hh"

namespace isaac
{
namespace reference
{

/**
 ** \brief Representation of a position in a reference genome.
 **
 ** The reference position is identified by its contig id ant the actual
 ** position on the contig. It also encapsulates some metadata about the
 ** neighborhood of that position -- currently a flag indicating if there are
 ** any neighbors with 1 or 2 mismatches in the suffix.
 **
 ** There are two special ReferencePosition values at the moment:
 **     TooManyMatch - position to be used as an indicator of number of matches exceeding the threshold
 **     NoMatch      - position to indicate that no matching sequence found in the reference
 **
 ** The bits are packed so that the natural sort order of the value will be as follows:
 **     1. TooManyMatch positions are least-significant values so that it is easy to skip matches of seeds
 **        that have TooManyMatch in the list
 **     2. Regular positions ordered by their contig, position and neighbor flag
 **     3. NoMatch closes the ordered list, no more positions, except for more NoMatch are expected after
 **        a NoMatch is found
 **/
class ReferencePosition
{
public:
    typedef unsigned long value_type;
    enum SpecialPosition
    {
        TooManyMatch, // Excessive number of repeat or neighbor matches
        NoMatch       // No matches
    };

    explicit ReferencePosition(SpecialPosition special) :
        value_(((((TooManyMatch == special ? 0 : MAX_CONTIG_ID) << POSITION_BITS) | 0) << NEIGHBORS_BITS) | 0) {}

    /**
     * \brief Creates TooManyMatches position which normally should not be used
     */
    ReferencePosition() : value_(0UL) {}
    /// Encodes a contig Id and a position into a ReferencePosition
    ReferencePosition(
        const unsigned long contigId,
        const unsigned long position,
        const bool neighbors = false,
        const bool reverse = false)
        : value_((((((contigId + 1) << POSITION_BITS) | position) << NEIGHBORS_BITS) | neighbors) << REVERSE_BITS | reverse)
    {
        ISAAC_ASSERT_MSG(0 == (position >> POSITION_BITS), "Position exceeds maximum allowed:" << position);
        ISAAC_ASSERT_MSG(0 == ((contigId + 1) >> CONTIG_ID_BITS), "Contig exceeds maximum allowed: " << contigId);
        //ISAAC_ASSERT_MSG(0 == (neighbors >> NEIGHBORS_BITS), "Neighbors exceeds maximum allowed");
    }
    /// Direct conversion from an integer value
    explicit ReferencePosition(value_type value) : value_(value) {}

    /**
     * \brief translate position into a system where contigs are ordered differently
     */
    const ReferencePosition translateContig(
            const std::vector<unsigned> &contigTranslationTable) const
    {
        const unsigned contigValue = (value_ >> (POSITION_BITS + REVERSE_BITS + NEIGHBORS_BITS));
        if (contigValue)
        {
            const unsigned long translatedContig = contigTranslationTable.at(contigValue - 1) + 1;

            return ReferencePosition((translatedContig << (POSITION_BITS + REVERSE_BITS + NEIGHBORS_BITS)) |
                                     (value_ & POSITION_REVERSE_NEIGBORS_MASK));
        }
        // special postions don't get translated
        return *this;
    }
    unsigned long getLocation() const
    {
        ISAAC_ASSERT_MSG(!isNoMatch(), "Location cannot be requested from no match positions");
        ISAAC_ASSERT_MSG(!isTooManyMatch(), "Location cannot be requested from too many match positions");
        return (value_ >> (REVERSE_BITS + NEIGHBORS_BITS)) - (1UL << POSITION_BITS);
    }
    unsigned long getContigId() const
    {
        ISAAC_ASSERT_MSG(!isNoMatch(), "Contig id cannot be requested from no match positions");
        ISAAC_ASSERT_MSG(!isTooManyMatch(), "Contig id cannot be requested from too many match positions");
        return (value_ >> (POSITION_BITS + REVERSE_BITS + NEIGHBORS_BITS)) - 1;
    }
    unsigned long getPosition() const {return (value_ & POSITION_MASK) >> (REVERSE_BITS + NEIGHBORS_BITS);}
    bool hasNeighbors() const {return value_ & NEIGHBORS_MASK;}
    bool reverse() const {return value_ & REVERSE_MASK;}
    bool isNoMatch() const {return ReferencePosition(ReferencePosition::NoMatch) == *this;}
    bool isTooManyMatch() const
        {return (ReferencePosition(ReferencePosition::TooManyMatch).value_ >> (REVERSE_BITS + NEIGHBORS_BITS)) == (value_ >> (REVERSE_BITS + NEIGHBORS_BITS));}
    ReferencePosition(const ReferencePosition &p) {value_ = p.value_;}
    ReferencePosition &operator=(const ReferencePosition &p)
    {
        if (this != &p)
        {
            value_ = p.value_;
        }
        return *this;
    }
    ReferencePosition & setNeighbors(const bool neighbors)
    {
        value_ = neighbors ? (value_ | NEIGHBORS_MASK) : (value_ & ~NEIGHBORS_MASK);
        return *this;
    }
    ReferencePosition & setReverse(const bool reverse)
    {
        value_ = reverse ? (value_ | REVERSE_MASK) : (value_ & ~REVERSE_MASK);
        return *this;
    }
    value_type getValue() const {return value_;}
    bool operator<(const ReferencePosition &p) const {return value_ < p.value_;}
    bool operator>(const ReferencePosition &p) const {return p < *this;}
    bool operator>=(const ReferencePosition &p) const {return !(*this < p);}
    bool operator<=(const ReferencePosition &p) const {return !(p < *this);}
    bool operator==(const ReferencePosition &p) const {return value_ == p.value_;}
    bool operator!=(const ReferencePosition &p) const {return !(*this == p);}

    friend std::ostream &operator<<(std::ostream &os, const ReferencePosition r);
    ReferencePosition & operator += (const long offset)
    {
        ISAAC_ASSERT_MSG(0 == ((getPosition() + offset) >> POSITION_BITS),
                        "New position is negative or exceeds maximum allowed. " << *this <<
                        ", offset " << offset);
        value_ += offset << (REVERSE_BITS + NEIGHBORS_BITS);
        return *this;
    }

    ReferencePosition & operator -= (const long offset)
    {
        return *this += -offset;
    }

    long operator - (const ReferencePosition right) const
    {
        ISAAC_ASSERT_MSG(getContigId() == right.getContigId(), "Contigs must match");
        return getPosition() - right.getPosition();
    }

    ReferencePosition operator + (const long offset) const
    {
        ReferencePosition ret(*this);
        return ret += offset;
    }

    ReferencePosition operator - (const long offset) const
    {
        return (*this) + -offset;
    }

//    // it is not unusual for a metagenomics reference to have 16M contigs.
    static const unsigned long CONTIG_ID_BITS = 25;
    static const unsigned long POSITION_BITS = 37;
    static const unsigned long REVERSE_BITS = 1;
    static const unsigned long NEIGHBORS_BITS = 1;
    static const unsigned long MAX_CONTIG_ID = ((~static_cast<value_type>(0)) >> (POSITION_BITS + REVERSE_BITS + NEIGHBORS_BITS));
    static const unsigned long POSITION_MASK = (((~static_cast<value_type>(0)) >> (CONTIG_ID_BITS + REVERSE_BITS + NEIGHBORS_BITS)) << (REVERSE_BITS + NEIGHBORS_BITS));
    static const unsigned long POSITION_REVERSE_NEIGBORS_MASK = ((~static_cast<value_type>(0)) >> CONTIG_ID_BITS);
    static const unsigned long REVERSE_MASK = ((~static_cast<value_type>(0)) >> (CONTIG_ID_BITS + POSITION_BITS + NEIGHBORS_BITS));
    static const unsigned long NEIGHBORS_MASK = ((~static_cast<value_type>(0)) >> (CONTIG_ID_BITS + POSITION_BITS + REVERSE_BITS)) << REVERSE_BITS;
private:
    value_type value_;

    /*
     * \brief enable serialization
     */
    template <class Archive> friend void serialize(Archive &ar, ReferencePosition &bm, const unsigned int version);
};

inline std::ostream &operator<<(std::ostream &os, const ReferencePosition r)
{
    if (r.isNoMatch())
    {
        return os << "ReferencePosition(nomatch)";
    }
    else if (r.isTooManyMatch())
    {
        return os << "ReferencePosition(toomanymatch)";
    }
    else
    {
        return os << "ReferencePosition(" << r.getContigId() << ":" << r.getPosition() << ":" << r.hasNeighbors() << (r.reverse() ? 'r' : 'f') << ")";
    }
}

} // namespace reference
} // namespace isaac

#endif // #ifndef iSAAC_REFERENCE_REFERENCE_POSITION_HH
