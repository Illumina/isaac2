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
 ** \file TemplateBuilder.hh
 **
 ** \brief Construction of BamTemplate instances
 ** 
 ** \author Come Raczy
 **/

#ifndef iSAAC_ALIGNMENT_TEMPLATE_BUILDER_HH
#define iSAAC_ALIGNMENT_TEMPLATE_BUILDER_HH

#include <boost/noncopyable.hpp>

#include "alignment/BamTemplate.hh"
#include "alignment/BandedSmithWaterman.hh"
#include "alignment/Cigar.hh"
#include "alignment/Cluster.hh"
#include "alignment/Match.hh"
#include "alignment/RestOfGenomeCorrection.hh"
#include "alignment/SeedMetadata.hh"
#include "alignment/ShadowAligner.hh"
#include "alignment/TemplateLengthStatistics.hh"
#include "flowcell/ReadMetadata.hh"
#include "reference/Contig.hh"
#include "templateBuilder/BestPairInfo.hh"

namespace isaac
{
namespace alignment
{

/**
 ** \brief Utility component creating Template instances from Seed Matches.
 **
 ** The intended use is to create an instance of a TemplateBuilder for each
 ** thread and delegate to that instance the identification of the most likely
 ** template for each cluster. This is done by invoking the build method on the
 ** complete set of matches identified so far for the cluster. In the build
 ** method, the TemplateBuilder will do the alignment (first gapped, then
 ** ungapped), calculate the alignment quality of the individual fragments,
 ** select the most likely combination of fragments, resolve repeats and try
 ** aligning orphans.
 **
 ** TODO: add Template statistics
 **/
class TemplateBuilder: boost::noncopyable
{
public:
	typedef short DodgyAlignmentScore;
	static const DodgyAlignmentScore DODGY_ALIGNMENT_SCORE_UNKNOWN=255;
	static const DodgyAlignmentScore DODGY_ALIGNMENT_SCORE_UNALIGNED=-1;

    /**
     ** \brief Construct a template builder for a reference genome and a given
     ** set of reads.
     **
     ** TODO: Add support for clipping (variable length genomes)
     **/
    TemplateBuilder(
        const flowcell::FlowcellLayoutList &flowcellLayoutList,
        const unsigned repeatThreshold,
        const unsigned maxSeedsPerRead,
        const bool scatterRepeats,
        const bool rescueShadows,
        const bool trimPEAdaptors,
        const bool anchorMate,
        const unsigned gappedMismatchesMax,
        const unsigned smitWatermanGapsMax,
        const bool smartSmithWaterman,
        const bool noSmithWaterman,
        const bool splitAlignments,
        const int gapMatchScore,
        const int gapMismatchScore,
        const int gapOpenScore,
        const int gapExtendScore,
        const int minGapExtendScore,
        const unsigned splitGapLength,
        const DodgyAlignmentScore dodgyAlignmentScore);

    bool buildFragments(
        const reference::ContigList &contigList,
        const isaac::reference::ContigAnnotations &kUniqenessAnnotation,
        const flowcell::ReadMetadataList &readMetadataList,
        const SeedMetadataList &seedMetadataList,
        const matchSelector::SequencingAdapterList &sequencingAdapters,
        const TemplateLengthStatistics &templateLengthStatistics,
        const std::vector<Match>::const_iterator matchBegin,
        const std::vector<Match>::const_iterator matchEnd,
        const Cluster &cluster,
        const bool withGaps)
    {
        return fragmentBuilder_.build(contigList, kUniqenessAnnotation, readMetadataList, seedMetadataList, sequencingAdapters,
                                      templateLengthStatistics, matchBegin, matchEnd, cluster, withGaps);
    }

    const std::vector<std::vector<FragmentMetadata> > &getFragments() const {return fragmentBuilder_.getFragments();}

    /**
     ** \brief Build the most likely template for a single cluster, givena set of fragments
     **
     ** \return false means template ended up not having a single read aligned anywhere.
     **
     ** This method will initialize the internal template of the builder.
     **
     ** Precondition: the input list of fragments is partitioned by readIndex
     ** and sorted genomic position. This means that the order is (tileId,
     ** clusterId, seedIndex, reverse, contig, position) where 'tileId' and
     ** 'clusterId' are constant.
     **/
    bool buildTemplate(
        const std::vector<reference::Contig> &contigList,
        const isaac::reference::ContigAnnotations &kUniqeness,
        const RestOfGenomeCorrection &restOfGenomeCorrection,
        const flowcell::ReadMetadataList &readMetadataList,
        const matchSelector::SequencingAdapterList &sequencingAdapters,
        const Cluster &cluster,
        const TemplateLengthStatistics &templateLengthStatistics,
        const unsigned mapqThreshold);

    /**
     * \brief Same as above but unit testing friendly.
     */
    bool buildTemplate(
        const std::vector<reference::Contig> &contigList,
        const isaac::reference::ContigAnnotations &kUniqeness,
        const RestOfGenomeCorrection &restOfGenomeCorrection,
        const flowcell::ReadMetadataList &readMetadataList,
        const matchSelector::SequencingAdapterList &sequencingAdapters,
        const std::vector<std::vector<FragmentMetadata> > &fragments,
        const Cluster &cluster,
        const TemplateLengthStatistics &templateLengthStatistics);

    /**
     ** \brief Getter for the BamTemplate
     **/
    const BamTemplate &getBamTemplate() const {return bamTemplate_;}
    BamTemplate &getBamTemplate() {return bamTemplate_;}
private:
    // when considering orphans for shadow alignment, don't look at those that are further than
    // orphanLogProbabilitySlack_ away from the best orphan
    static const double ORPHAN_LOG_PROBABILITY_SLACK_;
    static const unsigned READS_IN_A_PAIR = 2;
    static const unsigned SKIP_ORPHAN_HAMMING_DISTANCE = 3;

    // Maximum alignment score given to fragments and templates that are not well anchored but don't have any mismatches
    static const unsigned DODGY_BUT_CLEAN_ALIGNMENT_SCORE = 10U;

    const bool scatterRepeats_;
    const bool rescueShadows_;
    const bool trimPEAdaptors_;
    const bool anchorMate_;
    const DodgyAlignmentScore dodgyAlignmentScore_;
    const AlignmentCfg alignmentCfg_;

    /// Helper component to align fragments individually
    FragmentBuilder fragmentBuilder_;
    /// Cached storage for iterative template building
    BamTemplate bamTemplate_;
    /// Helper component to rescue shadow reads (or poorly aligned fragments)
    ShadowAligner shadowAligner_;

    /// Buffer for the cigar strings of rescued clusters
    Cigar cigarBuffer_;
    /// Buffer for the list of shadows rescued by the shadow aligner
    std::vector<FragmentMetadata> shadowList_;

    /// temporary buffer for modified alignments
    std::vector<FragmentMetadata> trimmedAlignments_;

    /// Holds the information about pairs obtained via combining alignments from MatchFinder
    templateBuilder::BestPairInfo bestCombinationPairInfo_;
    /// Holds the information about the pairs rescued via rescueShadow or buildDisjoinedTemplate
    templateBuilder::BestPairInfo bestRescuedPair_;

    /// Helper method to select the best fragment for single-ended runs
    bool pickBestFragment(
        const RestOfGenomeCorrection &restOfGenomeCorrection,
        const TemplateLengthStatistics &templateLengthStatistics,
        const std::vector<FragmentMetadata> &fragmentList,
        BamTemplate &result);
    bool rescueShadow(
        const std::vector<reference::Contig> &contigList,
        const isaac::reference::ContigAnnotations &kUniqeness,
        const RestOfGenomeCorrection &restOfGenomeCorrection,
        const flowcell::ReadMetadataList &readMetadataList,
        const matchSelector::SequencingAdapterList &sequencingAdapters,
        const std::vector<std::vector<FragmentMetadata> > &fragments,
        const FragmentMetadataList::const_iterator bestOrphanIterator,
        const unsigned orphanIndex,
        const unsigned shadowIndex,
        const TemplateLengthStatistics &templateLengthStatistics,
        templateBuilder::BestPairInfo &ret);

    /// Helper method to find the best pair of fragments for paired-end runs
    bool pickBestPair(
        const std::vector<reference::Contig> &contigList,
        const isaac::reference::ContigAnnotations &kUniqenessAnnotation,
        const RestOfGenomeCorrection &restOfGenomeCorrection,
        const flowcell::ReadMetadataList &readMetadataList,
        const matchSelector::SequencingAdapterList &sequencingAdapters,
        const std::vector<FragmentMetadataList> &fragments,
        const TemplateLengthStatistics &templateLengthStatistics);

    bool locateBestAnchoredPair(
        const std::vector<reference::Contig> &contigList,
        const isaac::reference::ContigAnnotations &kUniqenessAnnotation,
        const RestOfGenomeCorrection &restOfGenomeCorrection,
        const flowcell::ReadMetadataList &readMetadataList,
        const std::vector<std::vector<FragmentMetadata> > &fragments,
        const TemplateLengthStatistics &templateLengthStatistics,
        templateBuilder::BestPairInfo &ret);

    bool buildPairedEndTemplate(
        const std::vector<reference::Contig> &contigList,
        const isaac::reference::ContigAnnotations &kUniqenessAnnotation,
        const flowcell::ReadMetadataList &readMetadataList,
        const RestOfGenomeCorrection &restOfGenomeCorrection,
        const TemplateLengthStatistics &templateLengthStatistics,
        const std::vector<std::vector<FragmentMetadata> > &fragments,
        templateBuilder::BestPairInfo &bestPairInfo);

    bool rescueDisjointedTemplate(
        const std::vector<reference::Contig> &contigList,
        const isaac::reference::ContigAnnotations &kUniqenessAnnotation,
        const RestOfGenomeCorrection &restOfGenomeCorrection,
        const flowcell::ReadMetadataList &readMetadataList,
        const matchSelector::SequencingAdapterList &sequencingAdapters,
        const std::vector<std::vector<FragmentMetadata> > &fragments,
        const TemplateLengthStatistics &templateLengthStatistics,
        const templateBuilder::BestPairInfo &knownBestPair,
        templateBuilder::BestPairInfo &ret);

    void pickRandomRepeatAlignment(
        const unsigned clusterId,
        const templateBuilder::BestPairInfo &bestPair,
        BamTemplate &bamTemplate) const;

    void scoreRescuedTemplate(
        const RestOfGenomeCorrection &restOfGenomeCorrection,
        const TemplateLengthStatistics &templateLengthStatistics,
        templateBuilder::BestPairInfo &bestPair,
        BamTemplate &bamTemplate) const;

    bool scoreDisjointedTemplate(
        const std::vector<std::vector<FragmentMetadata> > &fragments,
        const RestOfGenomeCorrection &restOfGenomeCorrection,
        const TemplateLengthStatistics &templateLengthStatistics,
        const FragmentMetadata *bestDisjointedFragments[READS_IN_A_PAIR],
        BamTemplate &bamTemplate) const;

    bool buildSingletonShadowTemplate(
        const RestOfGenomeCorrection &restOfGenomeCorrection,
        const TemplateLengthStatistics &templateLengthStatistics,
        const std::vector<std::vector<FragmentMetadata> > &fragments,
        const FragmentIterator bestOrphanIterator,
        const unsigned orphanIndex,
        const unsigned shadowIndex,
        BamTemplate &bamTemplate);

    FragmentMetadata cloneWithCigar(const FragmentMetadata &right);
    std::pair<FragmentMetadata &, FragmentMetadata &> checkTrimPEAdaptor(
        const std::vector<reference::Contig> &contigList,
        const isaac::reference::ContigAnnotations &kUniqenessAnnotation,
        const flowcell::ReadMetadataList &readMetadataList,
        const FragmentMetadata &r1Fragment,
        const FragmentMetadata &r2Fragment);

    std::pair<FragmentMetadata &, FragmentMetadata &> trimPEAdaptor(
        const std::vector<reference::Contig> &contigList,
        const isaac::reference::ContigAnnotations &kUniqenessAnnotation,
        const flowcell::ReadMetadataList &readMetadataList,
        const FragmentMetadata & forwardFragment,
        const FragmentMetadata & reverseFragment);

    FragmentMetadata trimForwardPEAdaptor(
        const std::vector<reference::Contig> &contigList,
        const isaac::reference::ContigAnnotations &kUniqenessAnnotation,
        const flowcell::ReadMetadataList &readMetadataList,
        const FragmentMetadata & forwardFragment,
        const reference::ReferencePosition &adaptorPosition);

    FragmentMetadata trimReversePEAdaptor(
        const std::vector<reference::Contig> &contigList,
        const isaac::reference::ContigAnnotations &kUniqenessAnnotation,
        const flowcell::ReadMetadataList &readMetadataList,
        const FragmentMetadata &reverseFragment,
        const reference::ReferencePosition &adaptorPosition);

    /// Helper method to select the fragment with the highest logProbability
    std::vector<FragmentMetadata>::const_iterator getBestFragment(
        const std::vector<FragmentMetadata> &fragmentList) const;

    /// Helper function to calculate the alignment score of a fragment
    bool updateMappingScore(
        FragmentMetadata &fragment,
        const RestOfGenomeCorrection &restOfGenomeCorrection,
        const TemplateLengthStatistics &templateLengthStatistics,
        const FragmentMetadata &listFragment,
        const std::vector<FragmentMetadata> &fragmentList,
        const bool forceWellAnchored) const;

    bool flagDodgyTemplate(BamTemplate &bamTemplate) const;
    void decideOnAsGoodPair(
        const std::vector<std::vector<FragmentMetadata> > &fragments,
        const bool isNewKUnique,
        const templateBuilder::PairInfo& pairInfo,
        const std::pair<FragmentMetadata&, FragmentMetadata&>& updatedAlignments,
        templateBuilder::BestPairInfo& ret);
};

} // namespace alignment
} // namespace isaac

#endif // #ifndef iSAAC_ALIGNMENT_TEMPLATE_BUILDER_HH
