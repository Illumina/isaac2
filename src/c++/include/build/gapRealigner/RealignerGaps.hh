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
 ** \file RealignerGaps.hh
 **
 ** Container for gaps required by gap realigner.
 ** 
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_BUILD_GAP_REALIGNER_REALIGNER_GAPS_HH
#define iSAAC_BUILD_GAP_REALIGNER_REALIGNER_GAPS_HH

#include "build/gapRealigner/Gap.hh"

namespace isaac
{
namespace build
{
namespace gapRealigner
{

class RealignerGaps
{
    // All the gaps in the sample
    gapRealigner::Gaps gapGroups_;
    // Deletion gaps sorted by their end position
    gapRealigner::Gaps deletionEndGroups_;

public:
    typedef gapRealigner::Gap GapType;

    void reserve(const std::size_t gaps);
    void unreserve();

    template<typename IteratorT>
    void addGaps(
        const reference::ReferencePosition fStrandPosition,
        const IteratorT cigarBegin, const IteratorT cigarEnd)
    {
        using alignment::Cigar;
        reference::ReferencePosition pos = fStrandPosition;
        IteratorT cigarIterator = cigarBegin;
        for(;cigarEnd != cigarIterator; ++cigarIterator)
        {
            const Cigar::Component decoded = Cigar::decode(*cigarIterator);
            if (decoded.second == Cigar::ALIGN)
            {
                pos += decoded.first;
            }
            else if (decoded.second == Cigar::INSERT)
            {
                addGap(gapRealigner::Gap(pos, -decoded.first));
            }
            else if (decoded.second == Cigar::DELETE)
            {
                addGap(gapRealigner::Gap(pos, decoded.first));
                pos += decoded.first;
            }
            else if (decoded.second == Cigar::SOFT_CLIP)
            {
                if (fStrandPosition == pos)
                {
                    ISAAC_ASSERT_MSG(cigarBegin == cigarIterator || cigarEnd == cigarIterator + 1,
                                     "First soft clip can be only the first or the last component of the cigar");
                    // not advancing pos as it does not include the soft-clipped area
                }
                else
                {
                    ISAAC_ASSERT_MSG(cigarEnd == cigarIterator + 1,
                                     "At most two soft-clips are expected with second one being the last component of the cigar");
                }
            }
            else
            {
                ISAAC_ASSERT_MSG(false, "Unexpected Cigar OpCode: " << decoded.second << " only canonical cigar codes are allowed");
            }
        }
    }

    void clear()
    {
        gapGroups_.clear();
    }

    void addGap(
        const gapRealigner::Gap &gap)
    {
        gapGroups_.push_back(gap);
//        ISAAC_THREAD_CERR_DEV_TRACE("Added " << gaps.back());
    }

    void finalizeGaps();

    std::size_t getGapsCount() const
    {
        return gapGroups_.size();
    }

    gapRealigner::GapsRange findGaps(
        const unsigned long clusterId,
        const reference::ReferencePosition binStartPos,
        const reference::ReferencePosition rangeBegin,
        const reference::ReferencePosition rangeEnd,
        gapRealigner::Gaps &foundGaps_) const;

    gapRealigner::GapsRange allGaps() const;

};

} // namespace gapRealigner
} // namespace build
} // namespace isaac

#endif // #ifndef iSAAC_BUILD_GAP_REALIGNER_REALIGNER_GAPS_HH
