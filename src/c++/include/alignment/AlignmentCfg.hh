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
 ** \file AlignmentConfig.hh
 **
 ** \brief Structure for various details of how alignment gets penalized for mismatches, gaps, etc
 ** 
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_ALIGNMENT_ALIGNMENT_CFG_HH
#define iSAAC_ALIGNMENT_ALIGNMENT_CFG_HH

namespace isaac
{
namespace alignment
{

struct AlignmentCfg
{
    AlignmentCfg(
        const int gapMatchScore,
        const int gapMismatchScore,
        const int gapOpenScore,
        const int gapExtendScore,
        const int minGapExtendScore,
        const unsigned splitGapLength)
    : gapMatchScore_(gapMatchScore)
    , gapMismatchScore_(gapMismatchScore)
    , gapOpenScore_(gapOpenScore)
    , gapExtendScore_(gapExtendScore)
    , minGapExtendScore_(minGapExtendScore)
    , normalizedMismatchScore_(gapMatchScore - gapMismatchScore)
    , normalizedGapOpenScore_(gapMatchScore - gapOpenScore)
    , normalizedGapExtendScore_(gapMatchScore - gapExtendScore)
    , normalizedMaxGapExtendScore_(-minGapExtendScore)
    , splitGapLength_(splitGapLength)
    {
    }

    const int gapMatchScore_;
    const int gapMismatchScore_;
    const int gapOpenScore_;
    const int gapExtendScore_;
    const int minGapExtendScore_;

    const unsigned normalizedMismatchScore_;
    const unsigned normalizedGapOpenScore_;
    const unsigned normalizedGapExtendScore_;
    const unsigned normalizedMaxGapExtendScore_;
    const unsigned splitGapLength_;
};

} // namespace alignment
} // namespace isaac

#endif // #ifndef iSAAC_ALIGNMENT_ALIGNMENT_CFG_HH
