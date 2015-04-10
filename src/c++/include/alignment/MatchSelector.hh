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
 ** \file MatchSelector.hh
 **
 ** \brief Selection the best matches among all possible candidates.
 **
 ** \author Come Raczy
 **/

#ifndef iSAAC_ALIGNMENT_MATCH_SELECTOR_HH
#define iSAAC_ALIGNMENT_MATCH_SELECTOR_HH

#include <string>
#include <vector>
#include <boost/filesystem.hpp>
#include <boost/noncopyable.hpp>
#include <boost/ptr_container/ptr_vector.hpp>

#include "flowcell/ReadMetadata.hh"
#include "flowcell/TileMetadata.hh"
#include "alignment/BclClusters.hh"
#include "alignment/Cluster.hh"
#include "alignment/FragmentBuilder.hh"
#include "alignment/TemplateBuilder.hh"
#include "alignment/Match.hh"
#include "alignment/MatchTally.hh"
#include "alignment/SeedMetadata.hh"
#include "alignment/TemplateLengthStatistics.hh"
#include "alignment/matchSelector/BufferingFragmentStorage.hh"
#include "alignment/matchSelector/MatchSelectorStats.hh"
#include "alignment/matchSelector/ParallelMatchLoader.hh"
#include "alignment/matchSelector/SemialignedEndsClipper.hh"
#include "alignment/matchSelector/OverlappingEndsClipper.hh"
#include "common/Threads.hpp"
#include "flowcell/BarcodeMetadata.hh"
#include "reference/AnnotationLoader.hh"
#include "reference/Contig.hh"

namespace isaac
{
namespace alignment
{

namespace bfs = boost::filesystem;

class MatchSelector: boost::noncopyable
{
public:
    typedef flowcell::TileMetadataList TileMetadataList;
    typedef flowcell::ReadMetadataList ReadMetadataList;
    /// Construction of an instance for a given reference
    MatchSelector(
        matchSelector::FragmentStorage &fragmentStorage,
        unsigned int maxThreadCount,
        const TileMetadataList &tileMetadataList,
        const flowcell::BarcodeMetadataList &barcodeMetadataList,
        const flowcell::FlowcellLayoutList &flowcellLayoutList,
        const reference::ContigLists &contigLists,
        const isaac::reference::ContigAnnotationsList &kUniquenessAnnotations,
        const unsigned repeatThreshold,
        const int mateDriftRange,
        const TemplateLengthStatistics &defaultTemplateLengthStatistics,
        const unsigned mapqThreshold,
        const bool perTileTls,
        const bool pfOnly,
        const unsigned baseQualityCutoff,
        const bool keepUnaligned,
        const bool clipSemialigned,
        const bool clipOverlapping,
        const bool scatterRepeats,
        const bool rescueShadows,
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
        const TemplateBuilder::DodgyAlignmentScore dodgyAlignmentScore);

    /**
     * \brief frees the major memory reservations to make it safe to use dynamic memory allocations again
     */
    void unreserve()
    {
        templateLengthDistribution_.unreserve();
        threadTemplateBuilders_.clear();
        std::vector<Cluster>().swap(threadCluster_);
        fragmentStorage_.unreserve();
        std::vector<matchSelector::MatchSelectorStats>().swap(threadStats_);
    }

    void dumpStats(const boost::filesystem::path &statsXmlPath);

    void parallelSelect(
        const MatchTally &matchTally,
        std::vector<alignment::TemplateLengthStatistics> &barcodeTemplateLengthStatistics,
        const flowcell::TileMetadata &tileMetadata,
        std::vector<Match> &matchList,
        const BclClusters &bclData);

private:
    // The threading code in selectTileMatches can't deal with exception cleanup. Let it just crash for now.
    common::UnsafeThreadVector computeThreads_;

    const TileMetadataList tileMetadataList_;
    /**
     * \brief threadBclFilePaths_ gets resized for every tile total readlength. If the tile read lengths
     *        changes from lower to bigger, more threadBclFilePaths_ strings get allocated which breaks the whole
     *        concept of allocating things once. For now this list contains tiles in the processing order so
     *        that the total read length goes only down. TODO: cleanup this mess for example by creating
     *        MatchSelector only for the group of tiles that have the same geometry.
     */
    const flowcell::BarcodeMetadataList &barcodeMetadataList_;
    const flowcell::FlowcellLayoutList &flowcellLayoutList_;
    const reference::ContigLists &contigLists_;
    const isaac::reference::ContigAnnotationsList &kUniquenessAnnotations_;
    const unsigned repeatThreshold_;

    const TemplateLengthStatistics userTemplateLengthStatistics_;
    const unsigned mapqThreshold_;
    const bool perTileTls_;
    const bool pfOnly_;
    const unsigned baseQualityCutoff_;
    const bool keepUnaligned_;
    const bool clipSemialigned_;
    const bool clipOverlapping_;
    const std::vector<matchSelector::SequencingAdapterList> barcodeSequencingAdapters_;

    std::vector<matchSelector::MatchSelectorStats> allStats_;
    std::vector<matchSelector::MatchSelectorStats> threadStats_;

    matchSelector::FragmentStorage &fragmentStorage_;

    std::vector<Cluster> threadCluster_;
    boost::ptr_vector<TemplateBuilder> threadTemplateBuilders_;
    std::vector<matchSelector::SemialignedEndsClipper> threadSemialignedEndsClippers_;
    std::vector<matchSelector::OverlappingEndsClipper> threadOverlappingEndsClippers_;
    TemplateLengthDistribution templateLengthDistribution_;

    mutable boost::mutex mutex_;

    void processMatchList(
        const RestOfGenomeCorrection &restOfGenomeCorrection,
        const flowcell::BarcodeMetadata &barcodeMetadata,
        std::pair<std::vector<Match>::const_iterator, std::vector<Match>::const_iterator> &matchRange,
        const flowcell::TileMetadata & tileMetadata,
        const BclClusters &bclData,
        const TemplateLengthStatistics & templateLengthStatistics,
        const unsigned threadNumber);


    /**
     * \brief Construct the contig list from the SortedReference XML
     */
    std::vector<reference::Contig> getContigList(
        const reference::SortedReferenceMetadata &sortedReferenceMetadata) const;

    TemplateLengthStatistics determineTemplateLength(
        const flowcell::TileMetadata &tileMetadata,
        const flowcell::BarcodeMetadata &barcodeMetadata,
        const std::vector<Match>::const_iterator barcodeMatchListBegin,
        const std::vector<Match>::const_iterator barcodeMatchListEnd,
        const BclClusters &bclData,
        const unsigned threadNumber);

    /**
     ** \brief Helper method to generate the 'rest of the genome' correction for
     ** uniquely aligned reads and fragments.
     **
     ** There is one value for each individual reads in the readMetadataList (at
     ** the corresponding location) and one additional value for cases whare all
     ** the reads match uniquely.
     **/
    std::vector<double> getRestOfGenomeCorrectionList(
        const std::vector<flowcell::ReadMetadata> &readMetadataList) const;

    /**
     * \return Iterator to the first match of the next cluster
     */
    static std::vector<Match>::const_iterator findNextCluster(
        std::vector<Match>::const_iterator currentClusterIt,
        std::vector<Match>::const_iterator endIt)
    {
        if (currentClusterIt == endIt)
        {
            return endIt;
        }
        const unsigned long currentClusterId = currentClusterIt->getCluster();
        const unsigned long currentTileBarcode = currentClusterIt->getTileBarcode();
        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(currentClusterId, "    match: " << *currentClusterIt);
        while ((++currentClusterIt != endIt) && (currentClusterId == currentClusterIt->getCluster()))
        {
            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(currentClusterId, "    match: " << *currentClusterIt);
            ISAAC_ASSERT_MSG(currentTileBarcode == currentClusterIt->getTileBarcode(), "Matches of the same cluster expected to have the same barcode and tile.");
        }
        return currentClusterIt;
    }

    void processAlignedCluster(
        const reference::ContigList& barcodeContigList,
        const isaac::reference::ContigAnnotations& barcodeKUniqeness,
        const flowcell::ReadMetadataList& tileReads,
        const SeedMetadataList& tileSeeds,
        const matchSelector::SequencingAdapterList& sequencingAdapters,
        const TemplateLengthStatistics& templateLengthStatistics,
        const std::vector<Match>::const_iterator matchBegin,
        const std::vector<Match>::const_iterator matchEnd,
        const RestOfGenomeCorrection& restOfGenomeCorrection,
        const unsigned threadNumber, TemplateBuilder& templateBuilder,
        Cluster& cluster, BamTemplate& bamTemplate,
        matchSelector::MatchSelectorStats& stats);
    void processUnalignedCluster(
        const flowcell::ReadMetadataList& tileReads,
        const TemplateLengthStatistics& templateLengthStatistics,
        std::vector<Match>::const_iterator matchBegin,
        const BclClusters& bclData, unsigned clusterId,
        BamTemplate& ourThreadBamTemplate, Cluster& ourThreadCluster,
        matchSelector::MatchSelectorStats& ourThreadStats);
};

} // namespace alignemnt
} // namespace isaac

#endif // #ifndef iSAAC_ALIGNMENT_MATCH_SELECTOR_HH
