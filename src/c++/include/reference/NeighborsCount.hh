/**
 ** Isaac Genome Alignment Software
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
 ** \file NeighborsCount.hh
 **
 ** \brief Basic declarations for counting neighbors
 **
 ** \author Roman Petrovski
 **/

#ifndef ISAAC_REFERENCE_NEIGHBORS_COUNT_HH
#define ISAAC_REFERENCE_NEIGHBORS_COUNT_HH

namespace isaac
{
namespace reference
{

struct NeighborsCount
{
    NeighborsCount(unsigned short count) : second(count){}
    unsigned short second;
}__attribute__ ((packed)); // this significantly reduces memory requirement for finder especially considering the fact
// that parallel sort needs twice the memory for processing.

typedef std::vector<NeighborsCount> NeighborCounts;

inline NeighborsCount
operator -(NeighborsCount const& l, NeighborsCount const& r)
{
    return NeighborsCount(l.second - r.second);
}

static const NeighborsCount NEIGHBORS_TOO_MANY = NeighborsCount(0)-1;

inline bool
operator <(NeighborsCount const& l, NeighborsCount const& r)
{
    return l.second < r.second;
}

inline bool
operator <=(NeighborsCount const& l, NeighborsCount const& r)
{
    return l.second <= r.second;
}

inline bool
operator !=(NeighborsCount const& l, NeighborsCount const& r)
{
    return l.second != r.second;
}

inline bool
operator ==(NeighborsCount const& l, NeighborsCount const& r)
{
    return l.second == r.second;
}

inline bool
operator >=(NeighborsCount const& l, NeighborsCount const& r)
{
    return l.second >= r.second;
}


inline NeighborsCount
operator +(NeighborsCount const& l, NeighborsCount const& r)
{
    if (NEIGHBORS_TOO_MANY - l >= r)
    {
        return NeighborsCount(l.second + r.second);
    }
    else
    {
        return NEIGHBORS_TOO_MANY;
    }
}

inline NeighborsCount &
operator +=(NeighborsCount &l, NeighborsCount const& r)
{
    if (NEIGHBORS_TOO_MANY - l >= r)
    {
        l.second += r.second;
    }
    else
    {
        l = NEIGHBORS_TOO_MANY;
    }
    return l;
}

inline NeighborsCount &
operator ++(NeighborsCount &l)
{
    if (NEIGHBORS_TOO_MANY != l)
    {
        ++l.second;
    }
    return l;
}

inline std::ostream &operator <<(std::ostream &os, const NeighborsCount &nc)
{
    return os << nc.second;
}

} // namespace reference
} //namespace isaac

#endif // #ifndef ISAAC_REFERENCE_NEIGHBORS_COUNT_HH
