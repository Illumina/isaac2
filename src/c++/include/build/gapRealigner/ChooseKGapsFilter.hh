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
 ** \file OverlappingGapsFilter.hh
 **
 ** Gap realigner implementation details.
 ** 
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_BUILD_GAP_REALIGNER_CHOOSE_K_GAPS_FILTER_HH
#define iSAAC_BUILD_GAP_REALIGNER_CHOOSE_K_GAPS_FILTER_HH

#include "build/gapRealigner/Gap.hh"
#include "common/BitHacks.hh"

namespace isaac
{
namespace build
{
namespace gapRealigner
{

template <typename BitTrackingType>
class ChooseKGapsFilter
{
    static const BitTrackingType ALL_ONES = BitTrackingType(0) - 1;
    static const unsigned BITS_IN_BYTE = 8;

    const unsigned gaps_;
    const unsigned maxK_;
    unsigned currentK_;
public:
    ChooseKGapsFilter(
        const gapRealigner::GapsRange &gapsRange,
        const unsigned maxK)
    : gaps_(gapsRange.size()), maxK_(std::min(maxK, gaps_)), currentK_(0)
    {
        ISAAC_ASSERT_MSG(gaps_ <= sizeof(BitTrackingType) * BITS_IN_BYTE, "Too many gaps: " << gaps_);
    }

    BitTrackingType next(BitTrackingType combination)
    {
        const BitTrackingType lastWithKBitsSet = (~(ALL_ONES >> currentK_)) >> (sizeof(BitTrackingType) * BITS_IN_BYTE - gaps_);

        if (!combination || lastWithKBitsSet == combination)
        {
            ++currentK_;
            if (maxK_ < currentK_)
            {
                return 0;
            }
            return ~(ALL_ONES << currentK_);
        }

//        ISAAC_THREAD_CERR << "before snoob combination: " << std::hex << combination << std::endl;
//        ISAAC_THREAD_CERR << "before snoob lastWithKBitsSet: " << std::hex << lastWithKBitsSet << std::endl;
//        ISAAC_THREAD_CERR << "before snoob currentK_: " << std::hex << currentK_ << std::endl;
        return common::snoob(combination);
    }
};

} // namespace gapRealigner
} // namespace build
} // namespace isaac

#endif // #ifndef iSAAC_BUILD_GAP_REALIGNER_CHOOSE_K_GAPS_FILTER_HH
