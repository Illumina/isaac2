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
 ** \file GapRealigner.hh
 **
 ** Attempts to reduce read mismatches by introducing gaps found on other reads.
 ** 
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_BUILD_GAP_REALIGNER_HH
#define iSAAC_BUILD_GAP_REALIGNER_HH

#include <boost/math/special_functions/binomial.hpp>

#include "alignment/Cigar.hh"
#include "alignment/TemplateLengthStatistics.hh"
#include "build/gapRealigner/RealignerGaps.hh"
#include "build/BarcodeBamMapping.hh"
#include "build/PackedFragmentBuffer.hh"
#include "flowcell/BarcodeMetadata.hh"
#include "reference/Contig.hh"
#include "reference/ReferencePosition.hh"

namespace isaac
{
namespace build
{

/**
 * \brief Attempts to insert gaps found on other fragments while preserving the ones that
 *        are already there.
 */
class GapRealigner
{
public:
    typedef unsigned long GapChoiceBitmask;
private:
    // number of bits that can represent the on/off state for each gap.
    // Currently unsigned is used to hold the choice
    static const unsigned MAX_GAPS_AT_A_TIME = 64;

    const bool realignGapsVigorously_;
    const bool realignDodgyFragments_;
    const unsigned gapsPerFragmentMax_;
    const unsigned combinationsLimit_;
    // Recommended value to be lower than gapOpenCost_ in a way that
    // no less than two mismatches would warrant adding a gap
    const unsigned mismatchCost_;// = 3;
    const unsigned gapOpenCost_;// = 4;
    // Recommended 0 as it does not matter how long the introduced gap is for realignment
    const unsigned gapExtendCost_;// = 0;
    static const int mismatchPercentReductionMin_ = 20;

    const bool clipSemialigned_;

    const flowcell::BarcodeMetadataList &barcodeMetadataList_;
    const std::vector<alignment::TemplateLengthStatistics> &barcodeTemplateLengthStatistics_;
    const std::vector<std::vector<reference::Contig> > &contigList_;

    gapRealigner::Gaps currentAttemptGaps_;

    gapRealigner::RealignerGaps fragmentGaps_;

public:
    typedef gapRealigner::Gap GapType;
    GapRealigner(
        const bool realignGapsVigorously,
        const bool realignDodgyFragments,
        const unsigned gapsPerFragmentMax,
        const unsigned mismatchCost,
        const unsigned gapOpenCost,
        const unsigned gapExtendCost,
        const bool clipSemialigned,
        const flowcell::BarcodeMetadataList &barcodeMetadataList,
        const std::vector<alignment::TemplateLengthStatistics> &barcodeTemplateLengthStatistics,
        const std::vector<std::vector<reference::Contig> > &contigList):
            realignGapsVigorously_(realignGapsVigorously),
            realignDodgyFragments_(realignDodgyFragments),
            gapsPerFragmentMax_(gapsPerFragmentMax),
            combinationsLimit_(boost::math::binomial_coefficient<double>(MAX_GAPS_AT_A_TIME, gapsPerFragmentMax_)),
            mismatchCost_(mismatchCost),
            gapOpenCost_(gapOpenCost),
            gapExtendCost_(gapExtendCost),
            clipSemialigned_(clipSemialigned),
            barcodeMetadataList_(barcodeMetadataList),
            barcodeTemplateLengthStatistics_(barcodeTemplateLengthStatistics),
            contigList_(contigList)
    {
        reserve();
    }

    void reserve()
    {
        currentAttemptGaps_.reserve(MAX_GAPS_AT_A_TIME * 10);
        // number of existing gaps to be expected in one fragment. No need to be particularly precise.
        fragmentGaps_.reserve(currentAttemptGaps_.capacity());
    }

    bool realign(
        const gapRealigner::RealignerGaps &realignerGaps,
        const reference::ReferencePosition binStartPos,
        const reference::ReferencePosition binEndPos,
        PackedFragmentBuffer::Index &index,
        io::FragmentAccessor &fragment,
        PackedFragmentBuffer &dataBuffer,
        alignment::Cigar &realignedCigars);

private:

    struct RealignmentBounds
    {
        /*
         * \brief Position of the first non soft-clipped base of the read
         */
        reference::ReferencePosition beginPos_;
        /*
         * \breif   Position of the first insertion base or the first base before the first deletion.
         *          If there are no indels, equals to endPos.
         */
        reference::ReferencePosition firstGapStartPos_;
        /*
         * \brief   Position of the first base following the last insertion or the first base
         *          that is not part of the last deletion. If there are no indels, equals to beginPos_
         */
        reference::ReferencePosition lastGapEndPos_;
        /*
         * \brief   Position of the base that follows the last non soft-clipped base of the read
         */
        reference::ReferencePosition endPos_;
    };
    friend std::ostream & operator << (std::ostream &os, const GapRealigner::RealignmentBounds &fragmentGaps);

    const gapRealigner::GapsRange findMoreGaps(
        gapRealigner::GapsRange range,
        const gapRealigner::Gaps &gaps,
        const reference::ReferencePosition binStartPos,
        const reference::ReferencePosition binEndPos);

    const gapRealigner::GapsRange findGaps(
        const unsigned sampleId,
        const reference::ReferencePosition binStartPos,
        const reference::ReferencePosition binEndPos,
        const reference::ReferencePosition rangeBegin,
        const reference::ReferencePosition rangeEnd);

    bool applyChoice(
        const GapChoiceBitmask &choice,
        const gapRealigner::GapsRange &gaps,
        const reference::ReferencePosition binEndPos,
        const reference::ReferencePosition contigEndPos,
        PackedFragmentBuffer::Index &index,
        const io::FragmentAccessor &fragment,
        alignment::Cigar &realignedCigars);


    struct GapChoice
    {
        GapChoice() : choice_(0), editDistance_(0), mismatches_(0), mismatchesPercent_(0), cost_(0), mappedLength_(0){}
        GapChoiceBitmask choice_;
        unsigned editDistance_;
        unsigned mismatches_;
        unsigned mismatchesPercent_;
        unsigned cost_;
        unsigned mappedLength_;
        reference::ReferencePosition startPos_;

        friend std::ostream & operator <<(std::ostream &os, const GapChoice &gapChoice)
        {
            return os << "GapChoice(" << gapChoice.choice_ << "," <<
                gapChoice.editDistance_ << "ed," <<
                gapChoice.mismatches_ << "mm," <<
                gapChoice.cost_ << "c," <<
                gapChoice.mappedLength_ << "ml," <<
                gapChoice.startPos_ << ")";
        }
    };


    GapChoice verifyGapsChoice(
        const GapChoiceBitmask &choice,
        const gapRealigner::GapsRange &gaps,
        const reference::ReferencePosition newBeginPos,
        const io::FragmentAccessor &fragment,
        const std::vector<reference::Contig> &reference);

    bool isBetterChoice(
        const GapChoice &choice,
        const unsigned maxMismatchesPercent,
        const GapChoice &bestChoice) const;

    const RealignmentBounds extractRealignmentBounds(const PackedFragmentBuffer::Index &index) const;

    bool findStartPos(
        const GapChoiceBitmask &choice,
        const gapRealigner::GapsRange &gaps,
        const reference::ReferencePosition binStartPos,
        const reference::ReferencePosition binEndPos,
        const unsigned pivotGapIndex,
        const reference::ReferencePosition pivotPos,
        long alignmentPos,
        reference::ReferencePosition &ret);

    bool compactCigar(
        const std::vector<reference::Contig> &reference,
        const reference::ReferencePosition binEndPos,
        PackedFragmentBuffer::Index &index,
        io::FragmentAccessor &fragment,
        alignment::Cigar &realignedCigars);

    GapChoice getAlignmentCost(
        const io::FragmentAccessor &fragment,
        const PackedFragmentBuffer::Index &index) const;

    void compactRealignedCigarBuffer(
        std::size_t bufferSizeBeforeRealignment,
        PackedFragmentBuffer::Index &index,
        alignment::Cigar &realignedCigars);

    void updatePairDetails(
        const PackedFragmentBuffer::Index &index,
        io::FragmentAccessor &fragment,
        PackedFragmentBuffer &dataBuffer);

    GapChoice findBetterGapsChoice(
        const gapRealigner::GapsRange& gaps,
        const reference::ReferencePosition& binStartPos,
        const reference::ReferencePosition& binEndPos,
        const std::vector<reference::Contig>& reference,
        const io::FragmentAccessor& fragment,
        const PackedFragmentBuffer::Index& index,
        unsigned &leftToEvaluate);

    long undoExistingGaps(const PackedFragmentBuffer::Index& index,
                          const reference::ReferencePosition& pivotPos);

    void verifyGapsChoice(
        const GapChoiceBitmask &choice,
        const gapRealigner::GapsRange& gaps,
        const reference::ReferencePosition& binStartPos,
        const reference::ReferencePosition& binEndPos,
        const io::FragmentAccessor& fragment,
        const std::vector<reference::Contig>& reference,
        const int originalMismatchesPercent,
        const long undoneAlignmentPos,
        GapChoice& bestChoice);
};

} // namespace build
} // namespace isaac

#endif // #ifndef iSAAC_BUILD_GAP_REALIGNER_HH
