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
 ** \file FragmentBuilder.cpp
 **
 ** \brief See FragmentBuilder.hh
 ** 
 ** \author Come Raczy
 **/

#include <numeric>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/foreach.hpp>
#include <boost/format.hpp>

#include "common/Exceptions.hh"
#include "alignment/Mismatch.hh"
#include "alignment/FragmentBuilder.hh"
#include "alignment/Quality.hh"
#include "common/Debug.hh"

namespace isaac
{
namespace alignment
{

FragmentBuilder::FragmentBuilder(
    const flowcell::FlowcellLayoutList &flowcellLayoutList,
    const unsigned repeatThreshold,
    const unsigned maxSeedsPerRead,
    const unsigned gappedMismatchesMax,
    const unsigned smitWatermanGapsMax,
    const bool smartSmithWaterman,
    const bool noSmithWaterman,
    const bool splitAlignments,
    const AlignmentCfg &alignmentCfg)
    : repeatThreshold_(repeatThreshold)
    , gappedMismatchesMax_(gappedMismatchesMax)
    , smitWatermanGapsMax_(smitWatermanGapsMax)
    // one seed generates up to repeat threshold matches
    // Assuming split read aligner in worst case will make pair each read no more than once
    //, alignmentsMax_((repeatThreshold_ * maxSeedsPerRead * READS_MAX) + (repeatThreshold_ * maxSeedsPerRead * READS_MAX + 1) / 2)
    , alignmentsMax_((repeatThreshold_ * maxSeedsPerRead) * (repeatThreshold_ * maxSeedsPerRead))
    , noSmithWaterman_(noSmithWaterman)
    , splitAlignments_(splitAlignments)
    , alignmentCfg_(alignmentCfg)
    , seedMatchCounts_(maxSeedsPerRead * READS_MAX)
    , fragments_(READS_MAX) // max number of read ends we ever have to deal with
    , cigarBuffer_(Cigar::getMaxOperationsForReads(flowcellLayoutList) * alignmentsMax_)
    , ungappedAligner_(alignmentCfg_)
    , gappedAligner_(flowcellLayoutList, smartSmithWaterman, alignmentCfg_)
    , splitReadAligner_(alignmentCfg_, splitAlignments_)
{
    std::for_each(fragments_.begin(), fragments_.end(), boost::bind(&std::vector<FragmentMetadata>::reserve, _1, alignmentsMax_));
}

void FragmentBuilder::clear()
{
    BOOST_FOREACH(std::vector<FragmentMetadata> &fragmentList, fragments_)
    {
        fragmentList.clear();
    }
    cigarBuffer_.clear();
    std::fill(seedMatchCounts_.begin(), seedMatchCounts_.end(), 0);
}

/**
 * \return true if at least one fragment was built.
 */
bool FragmentBuilder::build(
    const std::vector<reference::Contig> &contigList,
    const isaac::reference::ContigAnnotations &kUniqenessAnnotation,
    const flowcell::ReadMetadataList &readMetadataList,
    const SeedMetadataList &seedMetadataList,
    const matchSelector::SequencingAdapterList &sequencingAdapters,
    const TemplateLengthStatistics &templateLengthStatistics,
    std::vector<Match>::const_iterator matchBegin,
    const std::vector<Match>::const_iterator matchEnd,
    const Cluster &cluster,
    const bool withGaps)
{
    clear();
    if (matchBegin < matchEnd)
    {
        bool repeatSeeds = false;
        ISAAC_ASSERT_MSG(!matchBegin->isNoMatch(), "Fake match lists must be dealt with outside");
        ISAAC_ASSERT_MSG(cluster.getNonEmptyReadsCount() == readMetadataList.size(), "cluster geometry must match");
        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(cluster.getId(), "FragmentBuilder::build: cluster " << matchBegin->seedId.getCluster() << " (" << cluster.getId() << ")");

        for(;matchEnd != matchBegin && !matchBegin->isNoMatch(); ++matchBegin)
        {
            if (matchBegin->isTooManyMatch())
            {
                repeatSeeds = true;
                seedMatchCounts_[matchBegin->getSeedId().getSeed()] = repeatThreshold_;
                const unsigned matchReadIndex = seedMetadataList[matchBegin->getSeedId().getSeed()].getReadIndex();
                ISAAC_ASSERT_MSG(fragments_[matchReadIndex].empty(), "Too-many matches are expected to sort to the top");
                ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(cluster.getId(), "FragmentBuilder::build: toomanymatch read: " << matchReadIndex  << ", seed " << matchBegin->getSeedId().getSeed());
            }
            else
            {
                if (repeatThreshold_ == seedMatchCounts_[matchBegin->getSeedId().getSeed()]++)
                {
                    repeatSeeds = true;
                }
                addMatch(readMetadataList, seedMetadataList, *matchBegin, cluster);
            }
        }

        if (repeatSeeds)
        {
            // if some of the seeds produced more than repeatThreshold_ matches, we need to remove the seed anchors
            std::for_each(fragments_.begin(), fragments_.end(),
                          boost::bind(&FragmentBuilder::resetHighRepeatSeedsAnchors, this, _1));
        }

        return alignFragments(contigList, kUniqenessAnnotation, readMetadataList, seedMetadataList, sequencingAdapters, templateLengthStatistics, withGaps);
    }
    return false;
}

bool FragmentBuilder::alignFragments(
    const std::vector<reference::Contig> &contigList,
    const isaac::reference::ContigAnnotations &kUniqenessAnnotation,
    const flowcell::ReadMetadataList &readMetadataList,
    const SeedMetadataList &seedMetadataList,
    const matchSelector::SequencingAdapterList &sequencingAdapters,
    const TemplateLengthStatistics &templateLengthStatistics,
    const bool withGaps)
{
    bool alignmentsExist = false;
    unsigned fragmentIndex = 0;
    BOOST_FOREACH(std::vector<FragmentMetadata> &fragmentList, fragments_)
    {
        if (!fragmentList.empty())
        {
            alignmentsExist = true;
            consolidateDuplicateAlignments(fragmentList, false);

            ISAAC_DEV_TRACE_BLOCK(const Cluster & cluster = fragmentList.at(0).getCluster();)
            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(cluster.getId(), "Ungapped alignment for cluster " << cluster.getTile() <<
                ":" << cluster.getId() << " [" << fragmentIndex << "] " << sequencingAdapters.size() << " adapters");

            matchSelector::FragmentSequencingAdapterClipper adapterClipper(sequencingAdapters);

            ungappedAligner_.alignCandidates(
                contigList, kUniqenessAnnotation, readMetadataList, fragmentList, adapterClipper, cigarBuffer_);

            // make sure all positions that are in the list are unique.
            consolidateDuplicateAlignments(fragmentList, true);

            if (withGaps)
            {
                if (splitAlignments_)
                {
                    splitReadAligner_.alignSimpleSv(cigarBuffer_, contigList, kUniqenessAnnotation, readMetadataList, templateLengthStatistics, fragmentList);
                    consolidateDuplicateAlignments(fragmentList, true);
                }

                ISAAC_THREAD_CERR_DEV_TRACE("Gapped alignment for cluster " << cluster.getTile() <<
                                            ":" << cluster.getId() << " [" << fragmentIndex << "] " << sequencingAdapters.size() << " adapters");

                if (withGaps && !noSmithWaterman_)
                {
                    // If there are still bad alignments, try to do expensive smith-waterman on them.
                    gappedAligner_.realignBadUngappedAlignments(
                        gappedMismatchesMax_, smitWatermanGapsMax_, contigList, kUniqenessAnnotation, readMetadataList, fragmentList, adapterClipper, cigarBuffer_);
                }

                // gapped alignment and adapter trimming may have adjusted the alignment position
                consolidateDuplicateAlignments(fragmentList, true);
            }
        }
        ++fragmentIndex;
    }
    return alignmentsExist;
}

void FragmentBuilder::addMatch(
    const flowcell::ReadMetadataList &readMetadataList,
    const SeedMetadataList &seedMetadataList, const Match &match,
    const Cluster &cluster)
{
    const unsigned seedIndex = match.seedId.getSeed();
    const SeedMetadata &seedMetadata = seedMetadataList.at(seedIndex);
    const unsigned readIndex = seedMetadata.getReadIndex();
    const flowcell::ReadMetadata &readMetadata = readMetadataList.at(readIndex);
    const reference::ReferencePosition &seedLocation = match.location;
    const unsigned contigId = seedLocation.getContigId();
    const bool reverse = match.seedId.isReverse() != match.location.reverse();
    const long readPosition = getReadPosition(readMetadata, seedMetadata, seedLocation.getPosition(), reverse);

    fragments_[readIndex].push_back(FragmentMetadata(&cluster, readIndex, 0, 0, reverse, contigId, readPosition));
    FragmentMetadata &fragment = fragments_[readIndex].back();
    fragment.firstSeedIndex = seedIndex;
    if (!match.location.hasNeighbors())
    {
        // set seed anchors only for alignments where we've visited all candidates.
        fragments_[readIndex].back().setSeedAnchor(readMetadata, seedMetadata);
    }
    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(cluster.getId(), "    adding r" << readIndex <<
                                           " match " << match << " " << fragments_[readIndex].back());
}

/**
 * \brief Remove anchors introduced by the seeds which later on have been found to have over repeatThreshold_ matches
 */
void FragmentBuilder::resetHighRepeatSeedAnchors(FragmentMetadata &fragment) const
{
    if (seedMatchCounts_[fragment.getFirstSeedIndex()] >= repeatThreshold_)
    {
        fragment.highRepeat = true;
        if (fragment.isWellAnchored())
        {
            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment.getCluster().getId(), "resetting seed anchor for: " << fragment);
            fragment.resetAnchors();
        }
    }
}


/**
 * \brief Remove anchors introduced by the seeds which later on have been found to have over repeatThreshold_ matches
 */
void FragmentBuilder::resetHighRepeatSeedsAnchors(FragmentMetadataList &fragmentList) const
{
    std::for_each(fragmentList.begin(), fragmentList.end(), boost::bind(&FragmentBuilder::resetHighRepeatSeedAnchors, this, _1));
}

/**
 * \brief Removes entries designating the same alignment location and strand
 *
 * This is important not only to avoid repetitive processing but also to avoid good alignments
 * scored poorly because they are present multiple times in the list.
 *
 * \param removeUnaligned if true, entries with !isAlignment() are removed
 *
 * \postcondition the fragmentList is ordered
 * \postcondition The fragmentList contains unique alignments only
 */
void FragmentBuilder::consolidateDuplicateAlignments(FragmentMetadataList &fragmentList, const bool removeUnaligned)
{
    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragmentList[0].getCluster().getId(), "consolidateDuplicateFragments initial size: " << fragmentList.size());
    // although initially matches arrive ordered by location, gapped alignment might have moved the
    // start position of some.
    std::sort(fragmentList.begin(), fragmentList.end());
    std::vector<FragmentMetadata>::iterator lastFragment = fragmentList.begin();
    while(fragmentList.end() != lastFragment && removeUnaligned && !lastFragment->isAligned())
    {
        ++lastFragment;
    }

    lastFragment = fragmentList.erase(fragmentList.begin(), lastFragment);

    if (2 > fragmentList.size())
    {
        return;
    }

    ISAAC_ASSERT_MSG(fragmentList.end() != lastFragment, "Unexpected end in list of size " << fragmentList.size());
    ISAAC_ASSERT_MSG(fragmentList.end() != lastFragment + 1, "Unexpected end - 1 in list of size " << fragmentList.size());

    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragmentList[0].getCluster().getId(), "    keeping   " << *lastFragment);

    for (std::vector<FragmentMetadata>::iterator currentFragment = lastFragment + 1; fragmentList.end() != currentFragment; ++currentFragment)
    {
        if (removeUnaligned && !currentFragment->isAligned())
        {
            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragmentList[0].getCluster().getId(), "    ignoring " << *currentFragment);
        }
        else if (*lastFragment == *currentFragment)
        {
            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragmentList[0].getCluster().getId(), "    skipping " << *currentFragment);
            lastFragment->mergeAnchors(*currentFragment);
        }
        else
        {
            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragmentList[0].getCluster().getId(), "    keeping   " << *currentFragment);
            ++lastFragment;
            if (lastFragment != currentFragment)
            {
                *lastFragment = *currentFragment;
            }
        }
    }
    fragmentList.resize(1 + lastFragment - fragmentList.begin());
    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragmentList[0].getCluster().getId(), "consolidateDuplicateFragments consolidated size: " << fragmentList.size());
}

inline long FragmentBuilder::getReadPosition(
    const flowcell::ReadMetadata &readMetadata,
    const SeedMetadata &seedMetadata, const long seedPosition, const bool reverse) const
{
    const int seedOffset = seedMetadata.getOffset();
    if (reverse)
    {
        const unsigned readLength = readMetadata.getLength();
        const unsigned seedLength = seedMetadata.getLength();
        // 'seedPosition + seedLength + seedOffset' is the first position past the end of the read
        return seedPosition + seedLength + seedOffset - readLength;
    }
    else
    {
        return seedPosition - seedOffset;
    }
}

} // namespace alignment
} // namespace isaac
