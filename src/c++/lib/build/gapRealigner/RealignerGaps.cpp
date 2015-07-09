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
 ** \file RealignerGaps.cpp
 **
 ** Attempts to reduce read mismatches by introducing gaps found on other reads.
 **
 ** \author Roman Petrovski
 **/

#include "build/gapRealigner/RealignerGaps.hh"

namespace isaac
{
namespace build
{
namespace gapRealigner
{

void RealignerGaps::reserve(const std::size_t gaps)
{
    deletionEndGroups_.reserve(gaps);
    gapGroups_.reserve(gaps);
}

void RealignerGaps::unreserve()
{
    gapRealigner::Gaps().swap(gapGroups_);
    gapRealigner::Gaps().swap(deletionEndGroups_);
}


inline bool orderByGapStartAndTypeLength(const gapRealigner::Gap &left, const gapRealigner::Gap &right)
{
    // ordering by signed length is required to ensure that intermixing of insertions and deletions does not prevent
    // duplicate gaps from being removed in finalizeGaps.
    return left.getBeginPos() < right.getBeginPos() ||
        (left.getBeginPos() == right.getBeginPos() && left.length_ < right.length_);
}

inline bool orderByDeletionGapEnd(const gapRealigner::Gap &left, const gapRealigner::Gap &right)
{
    return left.getDeletionEndPos() < right.getDeletionEndPos();
}

/*
inline std::ostream &operator <<(std::ostream &os, const std::pair<Gaps::const_iterator, Gaps::const_iterator> &gaps)
{
    if (gaps.second == gaps.first)
    {
        return os << "(no gaps)";
    }

    BOOST_FOREACH(const gapRealigner::Gap &gap, gaps)
    {
        os << gap << ",";
    }
    return os;
}
*/

void RealignerGaps::finalizeGaps()
{
    std::sort(gapGroups_.begin(), gapGroups_.end(), orderByGapStartAndTypeLength);
    gapGroups_.erase(std::unique(gapGroups_.begin(), gapGroups_.end()), gapGroups_.end());

    std::remove_copy_if(gapGroups_.begin(), gapGroups_.end(), std::back_inserter(deletionEndGroups_),
                        !boost::bind(&gapRealigner::Gap::isDeletion, _1));
    std::sort(deletionEndGroups_.begin(), deletionEndGroups_.end(), orderByDeletionGapEnd);
}

gapRealigner::GapsRange RealignerGaps::allGaps() const
{
    return gapRealigner::GapsRange(gapGroups_.begin(), gapGroups_.end());
}

/**
 * \brief Find gaps given the position range.
 */
gapRealigner::GapsRange RealignerGaps::findGaps(
    const unsigned long clusterId,
    const reference::ReferencePosition binStartPos,
    const reference::ReferencePosition rangeBegin,
    const reference::ReferencePosition rangeEnd,
    gapRealigner::Gaps &foundGaps) const
{
    foundGaps.clear();

    //ISAAC_THREAD_CERR << "findGaps effective range: [" << earliestPossibleGapBegin << ";" << rangeEnd << ")" << std::endl;
    gapRealigner::GapsRange gapStarts;
//    ISAAC_THREAD_CERR_DEV_TRACE("findGaps all gaps: " << gapRealigner::GapsRange(sampleGaps.begin(), sampleGaps.end()));
    // the first one that begins on or after the rangeBegin
    gapStarts.first = std::lower_bound(gapGroups_.begin(), gapGroups_.end(),
                                 gapRealigner::Gap(rangeBegin, -1000000), orderByGapStartAndTypeLength);
    // the first one that ends on or after the rangeEnd
    gapStarts.second = std::lower_bound(gapStarts.first, gapGroups_.end(), gapRealigner::Gap(rangeEnd, 0), orderByGapStartAndTypeLength);

    gapRealigner::GapsRange gapEnds;
    gapEnds.first = std::lower_bound(deletionEndGroups_.begin(), deletionEndGroups_.end(),
                                     // gap length 1 is simply to prevent orderByDeletionGapEnd for failing an assertion
                                     gapRealigner::Gap(rangeBegin, 1), orderByDeletionGapEnd);
    gapEnds.second = std::lower_bound(gapEnds.first, deletionEndGroups_.end(),
                                      // gap length 1 is simply to prevent orderByDeletionGapEnd for failing an assertion
                                      gapRealigner::Gap(rangeEnd, 1), orderByDeletionGapEnd);

    if (foundGaps.capacity() < gapStarts.size() + gapEnds.size())
    {
        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(clusterId, "GapRealigner::findGaps: Too many gaps (" << gapStarts.size() << "+" << gapEnds.size() << ")");
    }
    else
    {
        foundGaps.insert(foundGaps.end(), gapStarts.first, gapStarts.second);
        foundGaps.insert(foundGaps.end(), gapEnds.first, gapEnds.second);

        if (!gapEnds.empty() && !gapStarts.empty())
        {
            // if both ranges contain some gaps, we need to consolidate...
            std::sort(foundGaps.begin(), foundGaps.end(), orderByGapStartAndTypeLength);
            foundGaps.erase(std::unique(foundGaps.begin(), foundGaps.end()), foundGaps.end());
        }
    }

    return gapRealigner::GapsRange(foundGaps.begin(), foundGaps.end());
}

} // namespace gapRealigner
} // namespace build
} // namespace isaac
