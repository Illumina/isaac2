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
 ** \file TemplateBuilder.cpp
 **
 ** \brief See TemplateBuilder.hh
 **
 ** \author Come Raczy
 **/

#include <numeric>
#include <limits>

#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/foreach.hpp>
#include <boost/format.hpp>

#include "alignment/Mismatch.hh"
#include "alignment/FragmentBuilder.hh"
#include "alignment/Quality.hh"
#include "alignment/TemplateBuilder.hh"
#include "oligo/KmerGenerator.hpp"

namespace isaac
{
namespace alignment
{

const unsigned TemplateBuilder::DODGY_BUT_CLEAN_ALIGNMENT_SCORE;
const double TemplateBuilder::ORPHAN_LOG_PROBABILITY_SLACK_ = 100.0;

/**
 * \brief checks if the alignment simply does not make sense regardless
 *        of how unique it is.
 *
 * A bad alignment is either:
 * - less than 32 consecutive matches and either
 * - an alignment with too many mismatches (more than 1/8 the of the unclipped length) or
 * - an alignment with a low log probability (arbitrarily set to the equivalent of a read with
 *   all bases being Q40 and mismatches on 1/4 of them)
 */
static const double LOG_MISMATCH_Q40 = Quality::getLogMismatchSlow(40);
inline bool isVeryBadAlignment(const FragmentMetadata &fragment)
{
    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment.getCluster().getId(), "fragment.mismatchCount > fragment.getMappedLength() / 8 " << fragment.mismatchCount << ">" << fragment.getMappedLength() / 8);
    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment.getCluster().getId(), "fragment.logProbability < LOG_MISMATCH_Q40 / 4 * fragment.getMappedLength() " << fragment.logProbability << "<" << LOG_MISMATCH_Q40 / 4 * fragment.getMappedLength());
    return fragment.matchesInARow < 32 &&
        (fragment.mismatchCount > fragment.getMappedLength() / 8 || fragment.logProbability < LOG_MISMATCH_Q40 / 4 * fragment.getMappedLength());
}

TemplateBuilder::TemplateBuilder(
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
    const DodgyAlignmentScore dodgyAlignmentScore)
    : scatterRepeats_(scatterRepeats)
    , rescueShadows_(rescueShadows)
    , trimPEAdaptors_(trimPEAdaptors)
    , anchorMate_(anchorMate)
    , dodgyAlignmentScore_(dodgyAlignmentScore)
    , alignmentCfg_(gapMatchScore, gapMismatchScore, gapOpenScore, gapExtendScore, minGapExtendScore, splitGapLength)
    , fragmentBuilder_(flowcellLayoutList, repeatThreshold, maxSeedsPerRead, gappedMismatchesMax, smitWatermanGapsMax,
                       smartSmithWaterman, noSmithWaterman, splitAlignments,
                       alignmentCfg_)
    , bamTemplate_()
    , shadowAligner_(flowcellLayoutList,
                     gappedMismatchesMax, smitWatermanGapsMax, smartSmithWaterman, noSmithWaterman, alignmentCfg_)
    , cigarBuffer_(10000)
    , shadowList_(templateBuilder::TRACKED_REPEATS_MAX_ONE_READ)
    , bestCombinationPairInfo_(repeatThreshold, maxSeedsPerRead)
    , bestRescuedPair_(repeatThreshold, maxSeedsPerRead)

{
    trimmedAlignments_.reserve(READS_IN_A_PAIR);
}
bool TemplateBuilder::buildTemplate(
    const std::vector<reference::Contig> &contigList,
    const isaac::reference::ContigAnnotations &kUniqeness,
    const RestOfGenomeCorrection &restOfGenomeCorrection,
    const flowcell::ReadMetadataList &readMetadataList,
    const matchSelector::SequencingAdapterList &sequencingAdapters,
    const Cluster &cluster,
    const TemplateLengthStatistics &templateLengthStatistics,
    const unsigned mapqThreshold)
{
    const std::vector<std::vector<FragmentMetadata> > &fragments = fragmentBuilder_.getFragments();
    bool ret = buildTemplate(
        contigList, kUniqeness, restOfGenomeCorrection, readMetadataList, sequencingAdapters, fragments, cluster, templateLengthStatistics);

    if (ret && bamTemplate_.hasAlignmentScore())
    {
        if (!bamTemplate_.isProperPair())
        {
            // for improper pairs, mark fragment individually unaligned if they are below threshold
            ret = bamTemplate_.filterLowQualityFragments(mapqThreshold);
        }
        else if (mapqThreshold > bamTemplate_.getAlignmentScore())
        {
            // mark the whole thing unaligned if the pair is below threshold
            bamTemplate_.filterLowQualityFragments(-1U);
            ret = false;
        }
        else
        {
            if (anchorMate_ && bamTemplate_.getFragmentCount() == READS_IN_A_PAIR && bamTemplate_.isProperPair())
            {
                FragmentMetadata &read1 = bamTemplate_.getFragmentMetadata(0);
                FragmentMetadata &read2 = bamTemplate_.getFragmentMetadata(1);
                // when anchorMate_ is set, rescue mates of proper pairs if one of the mates is k-unique
                if (!read1.dodgy || !read2.dodgy || (read1.isWellAnchored() && read2.isWellAnchored()))
                {
                    read1.dodgy = false;
                    read2.dodgy = false;
                }
            }
        }
    }
    return ret;
}
bool TemplateBuilder::buildTemplate(
    const std::vector<reference::Contig> &contigList,
    const isaac::reference::ContigAnnotations &kUniqenessAnnotation,
    const RestOfGenomeCorrection &restOfGenomeCorrection,
    const flowcell::ReadMetadataList &readMetadataList,
    const matchSelector::SequencingAdapterList &sequencingAdapters,
    const std::vector<std::vector<FragmentMetadata> > &fragments,
    const Cluster &cluster,
    const TemplateLengthStatistics &templateLengthStatistics)
{
    cigarBuffer_.clear();
    // Initialize the bamTemplate_ with unaligned fragments for the cluster
    bamTemplate_.initialize(readMetadataList, cluster);
    if (2 == readMetadataList.size() && 2 == fragments.size())
    {
        if (!fragments[0].empty() && !fragments[1].empty())
        {
            return pickBestPair(
                contigList, kUniqenessAnnotation, restOfGenomeCorrection, readMetadataList, sequencingAdapters, fragments, templateLengthStatistics);
        }
        else if (!fragments[0].empty() || !fragments[1].empty())
        {
            const unsigned orphanIndex = fragments[0].empty() ? 1 : 0;
            const unsigned shadowIndex = (orphanIndex + 1) % READS_IN_A_PAIR;
            const FragmentMetadataList::const_iterator bestOrphanIterator = getBestFragment(fragments[orphanIndex]);
            if (rescueShadows_ &&
                rescueShadow(
                    contigList, kUniqenessAnnotation, restOfGenomeCorrection, readMetadataList, sequencingAdapters,
                    fragments, bestOrphanIterator, orphanIndex, shadowIndex, templateLengthStatistics, bestRescuedPair_))
            {
                pickRandomRepeatAlignment(fragments[orphanIndex][0].getCluster().getId(), bestRescuedPair_, bamTemplate_);
                scoreRescuedTemplate(restOfGenomeCorrection, templateLengthStatistics, bestRescuedPair_, bamTemplate_);

                ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragments[orphanIndex][0].getCluster().getId(),"rescueShadow: rescued  template: " << bamTemplate_);

                return true;
            }
            else
            {
                return buildSingletonShadowTemplate(
                    restOfGenomeCorrection, templateLengthStatistics, fragments,
                    bestOrphanIterator, orphanIndex, shadowIndex, bamTemplate_);
            }
        }
        else
        {
            // TODO: implement a recovery mechanism for unaligned clusters
            return false;
        }
    }
    else if(1 == readMetadataList.size() || 1 == fragments.size())
    {
        ISAAC_ASSERT_MSG(fragments[1].empty(), "With single-ended data expecting the fragment to be placed at index 0");
        if (!fragments[0].empty())
        {
            return pickBestFragment(restOfGenomeCorrection, templateLengthStatistics, fragments[0], bamTemplate_);
        }
        else
        {
            // TODO: implement a recovery mechanism for unaligned clusters
            return false;
        }
    }
    else
    {
        ISAAC_ASSERT_MSG(false, (boost::format("TemplateBuilder supports at most 2 reads: %d reads found") % fragments.size()).str().c_str());
        return false;
    }
}


std::vector<FragmentMetadata>::const_iterator TemplateBuilder::getBestFragment(const std::vector<FragmentMetadata> &fragmentList) const
{
    common::FiniteCapacityVector<std::vector<FragmentMetadata>::const_iterator, templateBuilder::TRACKED_REPEATS_MAX_ONE_READ> bestFragments;

    unsigned bestFragmentScore = -1U;
    double bestFragmentLogProbability = -std::numeric_limits<double>::max();
    for(FragmentIterator fragmentIterator = fragmentList.begin();
        fragmentList.end() != fragmentIterator; ++fragmentIterator)
    {
        if (bestFragmentScore > fragmentIterator->smithWatermanScore ||
            (bestFragmentScore == fragmentIterator->smithWatermanScore &&
                ISAAC_LP_LESS(bestFragmentLogProbability, fragmentIterator->logProbability)))
        {
            bestFragmentScore = fragmentIterator->smithWatermanScore;
            bestFragmentLogProbability = fragmentIterator->logProbability;
            bestFragments.clear();
            bestFragments.push_back(fragmentIterator);
        }
        else if (bestFragmentScore == fragmentIterator->smithWatermanScore &&
            ISAAC_LP_EQUALS(bestFragmentLogProbability, fragmentIterator->logProbability))
        {
            bestFragments.push_back(fragmentIterator);
        }

    }

    const unsigned clusterId = fragmentList[0].getCluster().getId();
    const unsigned repeatIndex = scatterRepeats_ ? (clusterId % bestFragments.size()) : 0;

    ISAAC_THREAD_CERR_DEV_TRACE("TemplateBuilder::getBestFragment returning repeat " <<
                                repeatIndex << " out of " << bestFragments.size() << " cluster id: " << clusterId << " " <<
                                *bestFragments[repeatIndex]);

    return bestFragments[repeatIndex];

}

unsigned computeAlignmentScore(
    const double restOfGenomeCorrection,
    const double alignmentProbability,
    const double otherAlignmentsProbability)
{
    return floor(-10.0 * log10((otherAlignmentsProbability + restOfGenomeCorrection) / (otherAlignmentsProbability + alignmentProbability + restOfGenomeCorrection)));
}
/**
 * \brief Set mapping score given the probability of the best choice and probabilities of alternative choices
 *
 * \return true if the end is considered well anchored.
 */
bool TemplateBuilder::updateMappingScore(
    FragmentMetadata &fragment,
    const RestOfGenomeCorrection &restOfGenomeCorrection,
    const TemplateLengthStatistics &templateLengthStatistics,
    const FragmentMetadata &bestAlignment,
    const std::vector<FragmentMetadata> &fragmentList,
    const bool forceWellAnchored) const
{
    bool bestIgnored = false;
    ISAAC_ASSERT_MSG(fragmentList.end() == std::adjacent_find(fragmentList.begin(), fragmentList.end()), "Expecting unique alignments in fragment list");
    if (forceWellAnchored || fragment.isWellAnchored())
    {
        double neighborProbability = 0;
        for (std::vector<FragmentMetadata>::const_iterator i = fragmentList.begin(); fragmentList.end() != i; ++i)
        {
            if (!bestIgnored && i->logProbability == bestAlignment.logProbability)
            {
                bestIgnored = true;
            }
            else
            {
                neighborProbability += exp(i->logProbability);
            }
        }

        fragment.alignmentScore = computeAlignmentScore(
            restOfGenomeCorrection.getReadRogCorrection(bestAlignment.getReadIndex()),
            exp(bestAlignment.logProbability),
            neighborProbability);
        ISAAC_THREAD_CERR_DEV_TRACE("TemplateBuilder::updateMappingScore: " << " sm=" << fragment.getAlignmentScore());
        return true;
    }
    else
    {
        ISAAC_THREAD_CERR_DEV_TRACE("TemplateBuilder::updateMappingScore: none of the candidate alignments anchor well");
        return false;
    }
}

FragmentMetadata TemplateBuilder::cloneWithCigar(const FragmentMetadata &right)
{
    FragmentMetadata ret = right;
    ret.cigarBuffer = &cigarBuffer_;
    ret.cigarOffset = cigarBuffer_.size();
    cigarBuffer_.insert(cigarBuffer_.end(), right.cigarBegin(), right.cigarEnd());
    return ret;
}

FragmentMetadata TemplateBuilder::trimForwardPEAdaptor(
    const std::vector<reference::Contig> &contigList,
    const isaac::reference::ContigAnnotations &kUniqenessAnnotation,
    const flowcell::ReadMetadataList &readMetadataList,
    const FragmentMetadata & forwardFragment,
    const reference::ReferencePosition &adaptorPosition)
{
    // trim forward fragment
    const unsigned cigarOffset = cigarBuffer_.size();
    cigarBuffer_.insert(cigarBuffer_.end(), forwardFragment.cigarBegin(), forwardFragment.cigarEnd());
    reference::ReferencePosition revPos = forwardFragment.getRStrandReferencePosition();
    unsigned long clip = 0;
    unsigned long align = 0;
    while (cigarOffset != cigarBuffer_.size())
    {
        const std::pair<int, Cigar::OpCode> cigar = Cigar::decode(cigarBuffer_.back());
        cigarBuffer_.pop_back();
        if(Cigar::SOFT_CLIP == cigar.second)
        {
            clip += cigar.first;
        }
        else if (Cigar::INSERT == cigar.second)
        {
            clip += cigar.first;
        }
        else if (Cigar::DELETE == cigar.second)
        {
            const long rightOverhang = revPos - adaptorPosition;
            if (rightOverhang < cigar.first)
            {
                break;
            }
            revPos -= cigar.first;
        }
        else if (Cigar::ALIGN == cigar.second)
        {
            const long rightOverhang = revPos - adaptorPosition;
            if (rightOverhang < cigar.first)
            {
                clip += rightOverhang;
                align = cigar.first - rightOverhang;
                break;
            }
            clip += cigar.first;
            revPos -= cigar.first;
        }
        else
        {
            // unexpected cigar operation. Just cancel
            cigarBuffer_.resize(cigarOffset);
            return forwardFragment;
        }
    }
    if (align)
    {
        cigarBuffer_.push_back(Cigar::encode(align, Cigar::ALIGN));
    }
    if (clip)
    {
        cigarBuffer_.push_back(Cigar::encode(clip, Cigar::SOFT_CLIP));
    }
    FragmentMetadata ret = forwardFragment;
    ret.resetAlignment();
    ret.rightClipped() = std::max<unsigned short>(clip, ret.rightClipped());
    ret.updateAlignment(alignmentCfg_, readMetadataList, contigList, kUniqenessAnnotation, forwardFragment.contigId, forwardFragment.position, cigarBuffer_, cigarOffset);
    return ret;
}

FragmentMetadata TemplateBuilder::trimReversePEAdaptor(
    const std::vector<reference::Contig> &contigList,
    const isaac::reference::ContigAnnotations &kUniqenessAnnotation,
    const flowcell::ReadMetadataList &readMetadataList,
    const FragmentMetadata & reverseFragment,
    const reference::ReferencePosition &adaptorPosition)
{
    const unsigned cigarOffset = cigarBuffer_.size();
    unsigned long clip = 0;
    unsigned long align = 0;
    reference::ReferencePosition fwPos = reverseFragment.getFStrandReferencePosition();

    Cigar::const_iterator current = reverseFragment.cigarBegin();
    for(;
        reverseFragment.cigarEnd() != current;
        ++current)
    {
        const std::pair<int, Cigar::OpCode> cigar = Cigar::decode(*current);
        if(Cigar::SOFT_CLIP == cigar.second)
        {
            clip += cigar.first;
        }
        else if (Cigar::INSERT == cigar.second)
        {
            clip += cigar.first;
        }
        else if (Cigar::DELETE == cigar.second)
        {
            const long leftOverhang = adaptorPosition - fwPos;
            if (leftOverhang < cigar.first)
            {
                break;
            }
            fwPos += cigar.first;
        }
        else if (Cigar::ALIGN == cigar.second)
        {
            const long leftOverhang = adaptorPosition - fwPos;
            if (leftOverhang < cigar.first)
            {
                clip += leftOverhang;
                align = cigar.first - leftOverhang;
                fwPos += leftOverhang;
                break;
            }
            clip += cigar.first;
            fwPos += cigar.first;
        }
        else
        {
            // unexpected cigar operation. Just cancel
            return reverseFragment;
        }
    }

    if (clip)
    {
        cigarBuffer_.push_back(Cigar::encode(clip, Cigar::SOFT_CLIP));
    }
    if (align)
    {
        cigarBuffer_.push_back(Cigar::encode(align, Cigar::ALIGN));
    }
    cigarBuffer_.insert(cigarBuffer_.end(), current + 1, reverseFragment.cigarEnd());

    FragmentMetadata ret = reverseFragment;
    ret.resetAlignment();
    ret.updateAlignment(alignmentCfg_, readMetadataList, contigList, kUniqenessAnnotation, reverseFragment.contigId, fwPos.getPosition(), cigarBuffer_, cigarOffset);
    ret.leftClipped() = std::max<unsigned short>(clip, ret.leftClipped());
    return ret;
}

std::pair<FragmentMetadata &, FragmentMetadata &> TemplateBuilder::trimPEAdaptor(
    const std::vector<reference::Contig> &contigList,
    const isaac::reference::ContigAnnotations &kUniqenessAnnotation,
    const flowcell::ReadMetadataList &readMetadataList,
    const FragmentMetadata & forwardFragment,
    const FragmentMetadata & reverseFragment)
{
    trimmedAlignments_.clear();
    if (!forwardFragment.splitAlignment &&
        reverseFragment.getRStrandReferencePosition() >= forwardFragment.getFStrandReferencePosition() &&
        reverseFragment.getRStrandReferencePosition() < forwardFragment.getRStrandReferencePosition())
    {
        trimmedAlignments_.push_back(trimForwardPEAdaptor(
            contigList, kUniqenessAnnotation, readMetadataList, forwardFragment, reverseFragment.getRStrandReferencePosition()));
        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(forwardFragment.getCluster().getId(), " trim fw:" << forwardFragment << " to " << trimmedAlignments_.back());
    }
    else
    {
        trimmedAlignments_.push_back(forwardFragment);
    }

    if (!reverseFragment.splitAlignment &&
        forwardFragment.getFStrandReferencePosition() > reverseFragment.getFStrandReferencePosition() &&
        forwardFragment.getFStrandReferencePosition() <= reverseFragment.getRStrandReferencePosition())
    {
        trimmedAlignments_.push_back(trimReversePEAdaptor(
            contigList, kUniqenessAnnotation, readMetadataList, reverseFragment, forwardFragment.getFStrandReferencePosition()));
        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(reverseFragment.getCluster().getId(), " trim rv:" << reverseFragment << " to " << trimmedAlignments_.back());
    }
    else
    {
        trimmedAlignments_.push_back(reverseFragment);
    }

    return std::pair<FragmentMetadata &, FragmentMetadata &>(trimmedAlignments_.front(), trimmedAlignments_.back());
}

std::pair<FragmentMetadata &, FragmentMetadata &> TemplateBuilder::checkTrimPEAdaptor(
    const std::vector<reference::Contig> &contigList,
    const isaac::reference::ContigAnnotations &kUniqenessAnnotation,
    const flowcell::ReadMetadataList &readMetadataList,
    const FragmentMetadata &r1Fragment,
    const FragmentMetadata &r2Fragment)
{
    if (trimPEAdaptors_)
    {
        if (r1Fragment.getContigId() == r1Fragment.getContigId())
        {
            if (r1Fragment.isReverse() != r2Fragment.isReverse())
            {
                if (!r1Fragment.isReverse())
                {
                    return trimPEAdaptor(contigList, kUniqenessAnnotation, readMetadataList, r1Fragment, r2Fragment);
                }
                else
                {
                    std::pair<FragmentMetadata &, FragmentMetadata &> ret =
                        trimPEAdaptor(contigList, kUniqenessAnnotation, readMetadataList, r2Fragment, r1Fragment);
                    return std::pair<FragmentMetadata &, FragmentMetadata &>(ret.second, ret.first);
                }
            }
        }
    }
    trimmedAlignments_.clear();
    trimmedAlignments_.push_back(r1Fragment);
    trimmedAlignments_.push_back(r2Fragment);
    return std::pair<FragmentMetadata &, FragmentMetadata &>(trimmedAlignments_.front(), trimmedAlignments_.back());
}

/**
 * \brief   Heuristics that handles all combinations of discovered and a new pair when they are equivalent in terms of
 *          mismatches. The main goal is to accumulate the truly best equivalent ones whilst skipping candidates that
 *          are worse than those we've already accumulated If the accumulated are worse than the new ones, reset the
 *          accumulated
 *
 * old      old match   new         new match
 * k-uniq   model       k-unique    model       action
 *
 * 0        0           0           0           skip
 * 0        0           0           1           reset
 * 0        0           1           0           reset
 * 0        0           1           1           reset
 * 0        1           0           0           skip
 * 0        1           0           1           append
 * 0        1           1           0           reset
 * 0        1           1           1           reset
 * 1        0           0           0           skip
 * 1        0           0           1           skip
 * 1        0           1           0           append
 * 1        0           1           1           reset
 * 1        1           0           0           skip
 * 1        1           0           1           skip
 * 1        1           1           0           skip
 * 1        1           1           1           append
 */
void TemplateBuilder::decideOnAsGoodPair(
    const std::vector<std::vector<FragmentMetadata> > &fragments,
    const bool isNewKUnique,
    const templateBuilder::PairInfo& pairInfo,
    const std::pair<FragmentMetadata&, FragmentMetadata&>& updatedAlignments,
    templateBuilder::BestPairInfo& ret)
{
    if (!isNewKUnique)
    {
        if (!pairInfo.matchModel_)
        {
            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragments[0].front().getCluster().getId(), "asgood skip1 " << pairInfo);
        }
        else if (ret.isKUnique())
        {
            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragments[0].front().getCluster().getId(), "asgood skip2 " << pairInfo);
        }
        else if (ret.info_.matchModel_)
        {
            ret.appendBest(updatedAlignments.first, updatedAlignments.second);
            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragments[0].front().getCluster().getId(), " append1 " << ret);
        }
        else
        {
            ret.resetBest(pairInfo, updatedAlignments.first, updatedAlignments.second);
            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragments[0].front().getCluster().getId(), " reset as good1 " << ret);
        }
    }
    else
    {
        if (!ret.isKUnique())
        {
            ret.resetBest(pairInfo, updatedAlignments.first, updatedAlignments.second);
            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragments[0].front().getCluster().getId(), " reset as good2 " << ret);
        }
        else if (!ret.info_.matchModel_)
        {
            if (pairInfo.matchModel_)
            {
                ret.resetBest(pairInfo, updatedAlignments.first, updatedAlignments.second);
                ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragments[0].front().getCluster().getId(), " reset as good3 " << ret);
            }
            else
            {
                ret.appendBest(updatedAlignments.first, updatedAlignments.second);
                ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragments[0].front().getCluster().getId(), " append2 " << ret);
            }
        }
        else
        {
            if (pairInfo.matchModel_)
            {
                ret.appendBest(updatedAlignments.first, updatedAlignments.second);
                ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragments[0].front().getCluster().getId(), " append3 " << ret);
            }
            else
            {
                ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragments[0].front().getCluster().getId(), "asgood skip3 " << pairInfo);
            }
        }
    }
}

bool TemplateBuilder::locateBestAnchoredPair(
    const std::vector<reference::Contig> &contigList,
    const isaac::reference::ContigAnnotations &kUniqenessAnnotation,
    const RestOfGenomeCorrection &restOfGenomeCorrection,
    const flowcell::ReadMetadataList &readMetadataList,
    const std::vector<std::vector<FragmentMetadata> > &fragments,
    const TemplateLengthStatistics &templateLengthStatistics,
    templateBuilder::BestPairInfo &ret)
{
    ret.clear();

    std::size_t pairsCounted = 0;
    for(FragmentIterator r1Fragment = fragments[0].begin(); fragments[0].end() != r1Fragment; ++r1Fragment)
    {
        for(FragmentIterator r2Fragment = fragments[1].begin(); fragments[1].end() != r2Fragment; ++r2Fragment)
        {
            const std::pair<FragmentMetadata &, FragmentMetadata &> updatedAlignments =
                checkTrimPEAdaptor(contigList, kUniqenessAnnotation, readMetadataList, *r1Fragment, *r2Fragment);
            const templateBuilder::PairInfo pairInfo(
                updatedAlignments.first, updatedAlignments.second,
                templateLengthStatistics.matchModel(updatedAlignments.first, updatedAlignments.second));
            const bool isNewKUnique = anchorMate_ ?
                updatedAlignments.first.isKUnique() || updatedAlignments.second.isKUnique() :
                updatedAlignments.first.isKUnique() && updatedAlignments.second.isKUnique();
            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(
                fragments[0].front().getCluster().getId(), "locateBestAnchoredPair pair: " << updatedAlignments.first << "-" << updatedAlignments.second << " " << pairInfo);

            if (ret.isWorseThan(pairInfo))
            {
                // avoid picking non-k-unique pair if k-unique is available
                if (isNewKUnique || !ret.isKUnique())
                {
                    ret.resetBest(pairInfo, updatedAlignments.first, updatedAlignments.second);
                    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragments[0].front().getCluster().getId(), " reset better " << ret);
                }
                else
                {
                    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragments[0].front().getCluster().getId(), " broken pair " << pairInfo);
                }
            }
            else if (ret.isAsGood(pairInfo))
            {
                decideOnAsGoodPair(fragments, isNewKUnique, pairInfo, updatedAlignments, ret);
            }
            else
            {
                if (isNewKUnique || pairInfo.matchModel_)
                {
                    ret.appendProbabilities(updatedAlignments.first, updatedAlignments.second);
                    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragments[0].front().getCluster().getId(), " appendProbabilities " << pairInfo);
                }
                else
                {
                    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragments[0].front().getCluster().getId(), " !KU and !matchModel " << pairInfo);
                }
                ++pairsCounted;
            }
        }
    }

    if (!ret.empty())
    {
        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragments[0].front().getCluster().getId(), "locateBestAnchoredPair " << ret);
        return true;
    }
    else
    {
        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragments[0].front().getCluster().getId(), "locateBestAnchoredPair: nothing good...");
        return false;
    }

    return !ret.empty();
}

/**
 * \return false if template needs to go into unaligned bin
 */
bool TemplateBuilder::flagDodgyTemplate(BamTemplate &bamTemplate) const
{
    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(bamTemplate.getFragmentMetadata(0).getCluster().getId(), "flagDodgyTemplate with " << dodgyAlignmentScore_ << " : " << bamTemplate);

    if (DODGY_ALIGNMENT_SCORE_UNALIGNED == dodgyAlignmentScore_)
    {
        // both must sort into unaligned bin. setUnaligned will not do it.
        bamTemplate.getFragmentMetadata(0).setNoMatch();
        if (2 == bamTemplate.getFragmentCount())
        {
            bamTemplate.getFragmentMetadata(1).setNoMatch();
        }
        bamTemplate.resetAlignmentScore();
        return false;
    }
    else
    {
        bamTemplate.getFragmentMetadata(0).dodgy = true;
        if (2 == bamTemplate.getFragmentCount())
        {
            bamTemplate.getFragmentMetadata(1).dodgy = true;
        }
    }
    return true;
}

bool TemplateBuilder::rescueShadow(
    const std::vector<reference::Contig> &contigList,
    const isaac::reference::ContigAnnotations &kUniqenessAnnotation,
    const RestOfGenomeCorrection &restOfGenomeCorrection,
    const flowcell::ReadMetadataList &readMetadataList,
    const matchSelector::SequencingAdapterList &sequencingAdapters,
    const std::vector<FragmentMetadataList > &fragments,
    const FragmentMetadataList::const_iterator bestOrphanIterator,
    const unsigned orphanIndex,
    const unsigned shadowIndex,
    const TemplateLengthStatistics &templateLengthStatistics,
    templateBuilder::BestPairInfo &ret)
{
    ret.clear();
    for(FragmentMetadataList::const_iterator orphanIterator = fragments[orphanIndex].begin();
        fragments[orphanIndex].end() != orphanIterator; ++orphanIterator)
    {
        const FragmentMetadata &orphan = *orphanIterator;
        if (orphan.isHighRepeat())
        {
            continue;
        }

        shadowList_.clear();
        ISAAC_THREAD_CERR_DEV_TRACE("TemplateBuilder::rescueShadow Orphan: " << orphan);

        std::size_t skipOneShadow = 0;

        if (ISAAC_LP_LESS(orphan.logProbability + ORPHAN_LOG_PROBABILITY_SLACK_, bestOrphanIterator->logProbability))
        {
            ISAAC_THREAD_CERR_DEV_TRACE("TemplateBuilder::rescueShadow orphan too bad to try rescuing shadows");
        }
        else if (shadowAligner_.rescueShadow(
            contigList, kUniqenessAnnotation, orphan, shadowList_, readMetadataList, sequencingAdapters,
            templateLengthStatistics))
        {
            // the best shadow for this orphan is the first in the list
            const FragmentMetadata &bestRescued = shadowList_.front();

            if(isVeryBadAlignment(bestRescued))
            {
                ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(orphan.getCluster().getId(), "    Rescued shadow too bad: " << bestRescued);
            }
            else
            {
                ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(orphan.getCluster().getId(), "    Rescued shadow not too bad: " << bestRescued);
                const std::pair<FragmentMetadata &, FragmentMetadata &> updatedAlignments =
                    checkTrimPEAdaptor(contigList, kUniqenessAnnotation, readMetadataList, *orphanIterator, *shadowList_.begin());
                const templateBuilder::PairInfo pairInfo(updatedAlignments.first, updatedAlignments.second, true);
                if (ret.isWorseThan(pairInfo))
                {
                    // Shadow may be the original from shadowAligner_ which means we have to clone the CIGAR
                    ret.resetBest(pairInfo, updatedAlignments.first, cloneWithCigar(updatedAlignments.second));
                    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(orphan.getCluster().getId(), " reset " << ret);
                }
                else if (ret.isAsGood(pairInfo))
                {
                    // Shadow may be the original from shadowAligner_ which means we have to clone the CIGAR
                    ret.appendBest(updatedAlignments.first, cloneWithCigar(updatedAlignments.second));
                    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(orphan.getCluster().getId(), " append " << ret);
                }
                else
                {
                    ret.appendProbabilities(updatedAlignments.first, updatedAlignments.second);
                }
                // In some cases, not registering clipped and unclipped version of pair might improve alignment score
                skipOneShadow = 1;
            }
        }
        else if (!shadowList_.empty())
        {
            // The shadow hits a repetitive region next to one of the orphans
            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(orphan.getCluster().getId(), (boost::format("Shadow rescue hits a repeat. Orphan: %s") % orphan).str());
        }
        else
        {
            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(orphan.getCluster().getId(), "No shadows rescued");
        }

        for (FragmentIterator shadow = shadowList_.begin() + skipOneShadow; shadowList_.end() != shadow; ++shadow)
        {
            ret.appendProbabilities(orphan, *shadow);
        }
    }

    if (!ret.empty())
    {
        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragments[orphanIndex][0].getCluster().getId(),"rescueShadow: rescued  template: " << ret);
        return true;
    }

    return false;
}

bool TemplateBuilder::buildSingletonShadowTemplate(
    const RestOfGenomeCorrection &restOfGenomeCorrection,
    const TemplateLengthStatistics &templateLengthStatistics,
    const std::vector<std::vector<FragmentMetadata> > &fragments,
    const FragmentIterator bestOrphanIterator,
    const unsigned orphanIndex,
    const unsigned shadowIndex,
    BamTemplate &bamTemplate)
{
    FragmentMetadata &orphan = bamTemplate.getFragmentMetadata(orphanIndex);
    FragmentMetadata &shadow = bamTemplate.getFragmentMetadata(shadowIndex);

    orphan = *bestOrphanIterator;
    if (isVeryBadAlignment(orphan))
    {
        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(orphan.getCluster().getId(), "    Singleton too bad: " << orphan);
        // both must sort into unaligned bin. setUnaligned will not do it.
        orphan.setNoMatch();
        shadow.setNoMatch();
        return false;
    }
    else
    {
        // mark shadow as 'shadow' (singleton's position, etc)
        shadow.contigId = orphan.contigId;
        shadow.position = orphan.position;
        shadow.readIndex = shadowIndex;
        shadow.alignmentScore = 0;
        shadow.cigarLength = 0;
        updateMappingScore(orphan, restOfGenomeCorrection, templateLengthStatistics, *bestOrphanIterator, fragments[orphanIndex], true);
    }
    return true;
}

bool TemplateBuilder::rescueDisjointedTemplate(
    const std::vector<reference::Contig> &contigList,
    const isaac::reference::ContigAnnotations &kUniqenessAnnotation,
    const RestOfGenomeCorrection &restOfGenomeCorrection,
    const flowcell::ReadMetadataList &readMetadataList,
    const matchSelector::SequencingAdapterList &sequencingAdapters,
    const std::vector<std::vector<FragmentMetadata> > &fragments,
    const TemplateLengthStatistics &templateLengthStatistics,
    const templateBuilder::BestPairInfo &knownBestPair,
    templateBuilder::BestPairInfo &ret)
{
    const FragmentMetadata *bestDisjointedFragments[READS_IN_A_PAIR] =
    {
        knownBestPair.empty() ? &*getBestFragment(fragments[0]) : &knownBestPair.repeats_.front().getFragmentMetadata(0),
        knownBestPair.empty() ? &*getBestFragment(fragments[1]) : &knownBestPair.repeats_.front().getFragmentMetadata(1)
    };

    ret.clear();
    unsigned bestOrphanIndex = 0;

    for (unsigned orphanIndex = 0; READS_IN_A_PAIR > orphanIndex; ++orphanIndex)
    {
        for(FragmentIterator orphanIterator = fragments[orphanIndex].begin();
            fragments[orphanIndex].end() != orphanIterator;
            ++orphanIterator)
        {
            const FragmentMetadata &orphan = *orphanIterator;
            if (orphan.isHighRepeat())
            {
                continue;
            }

            const bool skipThisOrphan = (!knownBestPair.empty() ? //when ed cutoff is set, no point in trying orphans with greater edit distance
                    orphan.getMismatchCount() > (knownBestPair.bestPairMismatchCount() + SKIP_ORPHAN_HAMMING_DISTANCE) :
                    ISAAC_LP_LESS(orphan.logProbability + ORPHAN_LOG_PROBABILITY_SLACK_, bestDisjointedFragments[orphanIndex]->logProbability));

            std::size_t skipOneShadow = 0;
            shadowList_.clear();
            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(orphan.getCluster().getId(), "TemplateBuilder::rescueDisjointedTemplate Orphan: " << orphan);
            if (!skipThisOrphan &&
                shadowAligner_.rescueShadow(
                    contigList, kUniqenessAnnotation, orphan, shadowList_,
                    readMetadataList, sequencingAdapters, templateLengthStatistics))
            {
                const FragmentMetadata &bestRescued = shadowList_.front();

                const unsigned rescuedMismatches = orphan.getMismatchCount() + bestRescued.getMismatchCount();

                if (isVeryBadAlignment(bestRescued))
                {
                    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(orphan.getCluster().getId(), "    Rescued shadow too bad: " << bestRescued);
                }
                else if (knownBestPair.empty() || (knownBestPair.bestPairMismatchCount() + SKIP_ORPHAN_HAMMING_DISTANCE) >= rescuedMismatches)
                {
                    std::pair<FragmentMetadata &, FragmentMetadata&> updatedAlignments =
                        checkTrimPEAdaptor(contigList, kUniqenessAnnotation, readMetadataList, *orphanIterator, *shadowList_.begin());

                    if (1 == knownBestPair.repeats_.size())
                    {
                        // carry over the non k-unique anchor if we happen to rescue the shadow which we've discovered through seed matches
                        const FragmentMetadata &knownShadowRead = knownBestPair.repeats_.front().getFragmentMetadata(updatedAlignments.second.getReadIndex());
                        if (knownShadowRead == updatedAlignments.second)
                        {
                            updatedAlignments.second.highRepeat = knownShadowRead.highRepeat;
                            updatedAlignments.second.mergeAnchors(knownShadowRead);
                        }
                    }

                    const templateBuilder::PairInfo pairInfo(updatedAlignments.first, updatedAlignments.second, true);
                    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(orphan.getCluster().getId(), "    comparing : " << ret << " vs " << pairInfo);
                    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(orphan.getCluster().getId(), "    comparing : " << ret.info_.logProbability_ - pairInfo.logProbability_);
                    if (ret.isWorseThan(pairInfo))
                    {
                        bestOrphanIndex = orphanIndex;
                        // Shadow may be the original from shadowAligner_ which means we have to clone the CIGAR
                        ret.resetBest(pairInfo, updatedAlignments.first, cloneWithCigar(updatedAlignments.second));
                        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(orphan.getCluster().getId(), " reset " << ret/* << *orphanIterator << *shadowList_.begin()*/);
                    }
                    else if (bestOrphanIndex == orphanIndex && // Avoid tracking the same pair twice. This would cause assertion failure in scoreRescuedTemplate
                        ret.isAsGood(pairInfo))
                    {
                        // Shadow may be the original from shadowAligner_ which means we have to clone the CIGAR
                        ret.appendBest(updatedAlignments.first, cloneWithCigar(updatedAlignments.second));
                        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(orphan.getCluster().getId(), " appendBest " << updatedAlignments.first << "-" << updatedAlignments.second);
                    }
                    else
                    {
                        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(orphan.getCluster().getId(), " append " << updatedAlignments.first << "-" << updatedAlignments.second);
                        ret.appendProbabilities(updatedAlignments.first, updatedAlignments.second);
                    }
                    // In some cases, not registering clipped and unclipped version of pair might improve alignment score
                    skipOneShadow = 1;
                }
                else
                {
                    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(orphan.getCluster().getId(), "rescueDisjointedTemplate: ignoring template: " << orphan << "-" << bestRescued << " for mismatches " << rescuedMismatches << ">" << (knownBestPair.bestPairMismatchCount() + SKIP_ORPHAN_HAMMING_DISTANCE));
                }
            }
            else if (!shadowList_.empty())
            {
                // The shadow hits a repetitive region next to one of the orphans
                ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(orphan.getCluster().getId(), "Mate rescue hits a repeat. Orphan: " << orphan);
            }

            for (FragmentIterator shadow = shadowList_.begin() + skipOneShadow; shadowList_.end() != shadow; ++shadow)
            {
                ret.appendProbabilities(orphan, *shadow);
            }
        }
    }

    if (!ret.empty())
    {
        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragments[0][0].getCluster().getId(), "rescueDisjointedTemplate: rescued  template: " << ret);
        return true;
    }
//    else if (!knownBestPair.empty())
//    {
//        // don't build disjointed template if there is a reasonable pair with knownBestPair.bestPairEditDistance already
//        // However, since we have not been able to rediscover this pair, it means it is either outside of
//        // consensus template boundaries or keeps hitting highly repetitive  locations with the shadow  all the time.
//        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragments[0][0].getCluster().getId(), "Disjointed template: nothing better but the rediscovery did not happen");
//        return flagDodgyTemplate(bamTemplate_);
//    }
    else
    {
        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragments[0][0].getCluster().getId(), "Disjointed template: Couldn't rescue anything");
        return false;
    }
}

bool TemplateBuilder::scoreDisjointedTemplate(
    const std::vector<std::vector<FragmentMetadata> > &fragments,
    const RestOfGenomeCorrection &restOfGenomeCorrection,
    const TemplateLengthStatistics &templateLengthStatistics,
    const FragmentMetadata *bestDisjointedFragments[READS_IN_A_PAIR],
    BamTemplate &bamTemplate) const
{
    FragmentMetadata &read1 = bamTemplate.getFragmentMetadata(0);
    FragmentMetadata &read2 = bamTemplate.getFragmentMetadata(1);
    read1 = *bestDisjointedFragments[0];
    read2 = *bestDisjointedFragments[1];
    bamTemplate.resetAlignmentScore();
    bamTemplate.setProperPair(false);

    updateMappingScore(read1, restOfGenomeCorrection, templateLengthStatistics, *bestDisjointedFragments[0], fragments[0], true);
    updateMappingScore(read2, restOfGenomeCorrection, templateLengthStatistics, *bestDisjointedFragments[1], fragments[1], true);
    return true;
}

void TemplateBuilder::pickRandomRepeatAlignment(
    const unsigned clusterId,
    const templateBuilder::BestPairInfo &bestPair,
    BamTemplate &bamTemplate) const
{
    const unsigned repeatIndex = scatterRepeats_ ? clusterId % bestPair.repeats_.size() : 0;

    bamTemplate = bestPair.repeats_[repeatIndex];

    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(clusterId, "TemplateBuilder::pickRandomRepeatAlignment: Picked repeat " << repeatIndex <<
                                " out of " << bestPair.repeats_.size() << " " << bamTemplate);
}

void TemplateBuilder::scoreRescuedTemplate(
    const RestOfGenomeCorrection &restOfGenomeCorrection,
    const TemplateLengthStatistics &templateLengthStatistics,
    templateBuilder::BestPairInfo &bestPair,
    BamTemplate &bamTemplate) const
{
    FragmentMetadata &r1Alignment = bamTemplate.getFragmentMetadata(0);
    FragmentMetadata &r2Alignment = bamTemplate.getFragmentMetadata(1);

    const double otherR1Probability = bestPair.sumUniqueReadProbabilities(0, r1Alignment.logProbability);
    const double otherR2Probability = bestPair.sumUniqueReadProbabilities(1, r2Alignment.logProbability);

    const double otherTemplateProbability = bestPair.sumUniquePairProbabilities(r1Alignment.logProbability + r2Alignment.logProbability);

    r1Alignment.alignmentScore = computeAlignmentScore(
        restOfGenomeCorrection.getReadRogCorrection(0),
        exp(r1Alignment.logProbability),
        otherR1Probability);

    r2Alignment.alignmentScore = computeAlignmentScore(
        restOfGenomeCorrection.getReadRogCorrection(1),
        exp(r2Alignment.logProbability),
        otherR2Probability);

    bamTemplate.setAlignmentScore(computeAlignmentScore(
        restOfGenomeCorrection.getRogCorrection(),
        bestPair.info_.probability(),
        otherTemplateProbability));

    bamTemplate.setProperPair(
        !r1Alignment.splitAlignment && !r2Alignment.splitAlignment &&
        alignment::TemplateLengthStatistics::Nominal == templateLengthStatistics.checkModel(r1Alignment, r2Alignment));

    ISAAC_ASSERT_MSG(
        bamTemplate.getAlignmentScore() < 4 || 1 == bestPair.removeRepeatDuplicates() ,
        "alignment score too high for a repeat of " << bestPair.repeats_.size() <<
        ":" << bamTemplate <<
        common::makeFastIoString(
            r1Alignment.getRead().getForwardSequence().begin(), r1Alignment.getRead().getForwardSequence().end()) <<
        "-" <<
        common::makeFastIoString(
            r2Alignment.getRead().getForwardSequence().begin(), r2Alignment.getRead().getForwardSequence().end()));
}

bool TemplateBuilder::pickBestFragment(
    const RestOfGenomeCorrection &restOfGenomeCorrection,
    const TemplateLengthStatistics &templateLengthStatistics,
    const std::vector<FragmentMetadata> &fragmentList,
    BamTemplate &result)
{
    if (!fragmentList.empty())
    {
        typedef std::vector<FragmentMetadata>::const_iterator FragmentIterator;
        const FragmentIterator bestFragment = getBestFragment(fragmentList);
        result.getFragmentMetadata(0) = *bestFragment;
        if (!updateMappingScore(result.getFragmentMetadata(0), restOfGenomeCorrection, templateLengthStatistics, *bestFragment, fragmentList, false))
        {
            return flagDodgyTemplate(result);
        }
        else
        {
            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(bamTemplate_.getFragmentMetadata(0).getCluster().getId(), "Single-ended template: " << result.getFragmentMetadata(0));
            return true;
        }
    }

    ISAAC_THREAD_CERR_DEV_TRACE("Single-ended template not aligned: ");
    return false;
}

bool TemplateBuilder::pickBestPair(
    const std::vector<reference::Contig> &contigList,
    const isaac::reference::ContigAnnotations &kUniqenessAnnotation,
    const RestOfGenomeCorrection &restOfGenomeCorrection,
    const flowcell::ReadMetadataList &readMetadataList,
    const matchSelector::SequencingAdapterList &sequencingAdapters,
    const std::vector<FragmentMetadataList> &fragments,
    const TemplateLengthStatistics &templateLengthStatistics)
{
    ISAAC_ASSERT_MSG(READS_IN_A_PAIR == fragments.size(), "TemplateBuilder::pickBestPair must be called for paired templates only");
    ISAAC_ASSERT_MSG(!fragments[0].empty() && !fragments[1].empty(), "TemplateBuilder::pickBestPair must be called for paired templates only");

    const bool anchoredPairExists = locateBestAnchoredPair(
        contigList, kUniqenessAnnotation, restOfGenomeCorrection, readMetadataList, fragments,
        templateLengthStatistics, bestCombinationPairInfo_);

    if (anchoredPairExists)
    {
        pickRandomRepeatAlignment(fragments[0][0].getCluster().getId(), bestCombinationPairInfo_, bamTemplate_);
        if (bamTemplate_.isHighRepeat())
        {
            return flagDodgyTemplate(bamTemplate_);
        }
        scoreRescuedTemplate(restOfGenomeCorrection, templateLengthStatistics, bestCombinationPairInfo_, bamTemplate_);

        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(bamTemplate_.getFragmentMetadata(0).getCluster().getId(), "Pair-end  template well anchored: " << bamTemplate_);
    }

    const bool anchoredPairMatchesModel = anchoredPairExists && 
        templateLengthStatistics.matchModel(bamTemplate_.getFragmentMetadata(0), bamTemplate_.getFragmentMetadata(1));

    if (anchoredPairMatchesModel &&
        (bestCombinationPairInfo_.isKUnique() ||
        (!bamTemplate_.getEditDistance() && bamTemplate_.isRepeat())))
    {
        //proper pair is k-unique or both ends have multiple perfect alignments. No point to spend any more time.
        return true;
    }

    if (!rescueShadows_)
    {
        if (!anchoredPairExists)
        {
            return false;
        }
        else if (!anchoredPairMatchesModel)
        {
            // if rescue shadows is not allowed, flag all abnormal pairs as dodgy even if they are k-unique.
            // Otherwise there are too many random picks
            return flagDodgyTemplate(bamTemplate_);
        }

        // anchored pair matches model and shadow rescuing is not allowed, just hope the pair is true...
        return true;
    }

    // nothing resolved from match combinations or the resolved pair does not have enough unique matches to
    // anchor it in a trustworthy way. Give rescuing a chance to find something better or lower-scored
    if (rescueDisjointedTemplate(contigList, kUniqenessAnnotation, restOfGenomeCorrection, readMetadataList, sequencingAdapters, fragments,
                                  templateLengthStatistics, bestCombinationPairInfo_, bestRescuedPair_))
    {
        BamTemplate rescuedTemplate;
        pickRandomRepeatAlignment(fragments[0][0].getCluster().getId(), bestRescuedPair_, rescuedTemplate);
        scoreRescuedTemplate(restOfGenomeCorrection, templateLengthStatistics, bestRescuedPair_, rescuedTemplate);
        bamTemplate_ = rescuedTemplate;
        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragments[0][0].getCluster().getId(), "rescueDisjointedTemplate: rescued  template: " << bamTemplate_);
        return true;
    }

    // rescuing did not give anything useful
    if (anchoredPairExists && bestCombinationPairInfo_.isKUnique())
    {
        // both seed-discovered ends are anchored but the pair is anomalous. Hope this is an actual SV
        return true;
    }

    // Couldn't rescue anything within the dominant template and and best combination pair is not k-unique.
    // Make sure this is flagged as dodgy
    const FragmentMetadata *bestDisjointedFragments[READS_IN_A_PAIR] = {&*getBestFragment(fragments[0]), &*getBestFragment(fragments[1])};
    return scoreDisjointedTemplate(
        fragments,
        restOfGenomeCorrection,
        templateLengthStatistics,
        bestDisjointedFragments,
        bamTemplate_);
    return flagDodgyTemplate(bamTemplate_);

}

} // namespace alignment
} // namespace isaac
