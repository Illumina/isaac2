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
 ** \file GappedAligner.cpp
 **
 ** \brief See GappedAligner.hh
 ** 
 ** \author Come Raczy
 **/
#include "alignment/fragmentBuilder/GappedAligner.hh"

namespace isaac
{
namespace alignment
{
namespace fragmentBuilder
{

const unsigned short GappedAligner::UNINITIALIZED_OFFSET_MAGIC;
const unsigned short GappedAligner::REPEAT_OFFSET_MAGIC;

GappedAligner::GappedAligner(
    const flowcell::FlowcellLayoutList &flowcellLayoutList,
    const bool smartSmithWaterman,
    const AlignmentCfg &alignmentCfg)
    : AlignerBase(alignmentCfg)
    , smartSmithWaterman_(smartSmithWaterman)
    , bandedSmithWaterman_(alignmentCfg.gapMatchScore_, alignmentCfg.gapMismatchScore_, -alignmentCfg.gapOpenScore_, -alignmentCfg.gapExtendScore_,
                           flowcell::getMaxTotalReadLength(flowcellLayoutList))
    , hashedQueryTile_(2, -1U)
    , hashedQueryCluster_(2, -1U)
    , hashedQueryReadIndex_(2, -1U)
    , queryKmerOffsets_(oligo::MaxKmer<HASH_KMER_LENGTH, unsigned short>::value + 1, UNINITIALIZED_OFFSET_MAGIC)
{
}

/// calculate the left and right flanks of the database WRT the query
static std::pair<unsigned, unsigned> getFlanks(
    const long strandPosition,
    const unsigned readLength,
    const unsigned long referenceSize,
    const unsigned widestGapSize)
{
    assert(widestGapSize >= 2);
    if (strandPosition >= widestGapSize/2)
    {
        // take into account the rounding
        if (strandPosition + readLength + (widestGapSize - widestGapSize/2) < static_cast<long>(referenceSize))
        {
            const unsigned left = widestGapSize/2;
            const unsigned right = widestGapSize - left - 1;
            return std::pair<unsigned, unsigned>(left, right);
        }
        else
        {
            assert((long)referenceSize >= readLength + strandPosition);
            const unsigned right = referenceSize - readLength - strandPosition;
            const unsigned left = widestGapSize - right - 1;
            return std::pair<unsigned, unsigned>(left, right);
        }
    }
    else
    {
        assert(strandPosition >= 0);
        const unsigned left = strandPosition;
        const unsigned right = widestGapSize - left - 1;
        return std::pair<unsigned, unsigned>(left, right);
    }
}

/**
 * \brief as smith-waterman takes about 50% of the time with only 2% of results accepted, it is worth
 *        having a cheap way of checking whether the expensive gapped alignment is going to provide a meaningful outcome
 */
bool GappedAligner::makesSenseToGapAlign(
    const unsigned tile, const unsigned cluster, const unsigned read, const bool reverse,
    const std::vector<char>::const_iterator queryBegin,
    const std::vector<char>::const_iterator queryEnd,
    const std::vector<char>::const_iterator databaseBegin,
    const std::vector<char>::const_iterator databaseEnd)
{
    if (hashedQueryTile_[reverse] != tile || hashedQueryCluster_[reverse] != cluster ||
        hashedQueryReadIndex_[reverse] != read)
    {
        oligo::KmerGenerator<HASH_KMER_LENGTH, unsigned> queryKmerGenerator(queryBegin, queryEnd);
        std::fill(queryKmerOffsets_.begin(), queryKmerOffsets_.end(), UNINITIALIZED_OFFSET_MAGIC);
        unsigned uniqueCount = 0;
        unsigned repeatCount = 0;
        unsigned kmer;
        std::vector<char>::const_iterator queryIt;
        while (queryKmerGenerator.next(kmer, queryIt))
        {
            if (UNINITIALIZED_OFFSET_MAGIC == queryKmerOffsets_[kmer])
            {
                queryKmerOffsets_[kmer] = (queryIt - queryBegin);
                ++uniqueCount;
            }
            else
            {
                queryKmerOffsets_[kmer] = REPEAT_OFFSET_MAGIC;
                ++repeatCount;
                --uniqueCount;
            }
        }
        hashedQueryTile_[reverse] = tile;
        hashedQueryCluster_[reverse] = cluster;
        hashedQueryReadIndex_[reverse] = read;
    }

    // counts of corresponding offsets of the first base of the data from the first base of the reference
    // offset of 0 is in the middle of the trackedOffsets.
    common::FiniteCapacityVector<unsigned char, QUERY_LENGTH_MAX> trackedOffsets(QUERY_LENGTH_MAX, 0);
    const int queryLength = std::distance(queryBegin, queryEnd);

    const unsigned maxOffset = queryLength * 2 + std::distance(databaseBegin, databaseEnd);
    ISAAC_ASSERT_MSG(maxOffset < trackedOffsets.size(),
                     "number of tracked offsets is too small: " << trackedOffsets.size() << " required: " << maxOffset);

    oligo::KmerGenerator<HASH_KMER_LENGTH, unsigned> databaseKmerGenerator(databaseBegin, databaseEnd);
    std::vector<char>::const_iterator databaseIt;
    int lastConfirmedOffset = std::numeric_limits<int>::max();
    unsigned kmer;
    while (databaseKmerGenerator.next(kmer, databaseIt))
    {
        const int queryOffset = queryKmerOffsets_[kmer];
        if (REPEAT_OFFSET_MAGIC == queryOffset)
        {
            // Ambiguous mapping between reference and data,
            // Ignore, assume this is not informative...
        }
        else if (UNINITIALIZED_OFFSET_MAGIC != queryOffset)
        {
            const int databaseOffset = std::distance(databaseBegin, databaseIt);
            const int firstBaseOffset = databaseOffset - queryOffset + queryLength;
            // look for two sufficiently confirmed anchoring points that will disagree about the read offset.
            // If not found, sequence is unlikely to produce gaps with smith-waterman.
            if (++trackedOffsets[firstBaseOffset] == SUFFICIENT_NUMBER_OF_HITS)
            {
                if (std::numeric_limits<int>::max() == lastConfirmedOffset)
                {
                    lastConfirmedOffset = firstBaseOffset;
                }
                else if (lastConfirmedOffset != firstBaseOffset)
                {
                    return true;
                }
            }
        }
    }

    return false;
}

unsigned GappedAligner::alignGapped(
    FragmentMetadata &fragmentMetadata,
    Cigar &cigarBuffer,
    const flowcell::ReadMetadataList &readMetadataList,
    const matchSelector::FragmentSequencingAdapterClipper &adapterClipper,
    const reference::ContigList &contigList,
    const isaac::reference::ContigAnnotations &kUniqenessAnnotation)
{
    const unsigned cigarOffset = cigarBuffer.size();
    fragmentMetadata.resetAlignment();
    fragmentMetadata.resetClipping();

    const Read &read = fragmentMetadata.getRead();
    const bool reverse = fragmentMetadata.reverse;
    const std::vector<char> &sequence = read.getStrandSequence(reverse);
    const reference::Contig &contig = contigList[fragmentMetadata.contigId];
    const std::vector<char> &reference = contig.forward_;

    std::vector<char>::const_iterator sequenceBegin = sequence.begin();
    std::vector<char>::const_iterator sequenceEnd = sequence.end();

    adapterClipper.clip(contig, fragmentMetadata, sequenceBegin, sequenceEnd);
    clipReadMasking(read, fragmentMetadata, sequenceBegin, sequenceEnd);

    clipReference(reference.size(), fragmentMetadata, sequenceBegin, sequenceEnd);

    const unsigned firstMappedBaseOffset = std::distance(sequence.begin(), sequenceBegin);
    if (firstMappedBaseOffset)
    {
        cigarBuffer.addOperation(firstMappedBaseOffset, Cigar::SOFT_CLIP);
    }

    const unsigned sequenceLength = std::distance(sequenceBegin, sequenceEnd);

    // position of the fragment on the strand
    long strandPosition = fragmentMetadata.position;
    ISAAC_ASSERT_MSG(0 <= strandPosition, "alignUngapped should have clipped reads beginning before the reference");

    // no gapped alignment if the reference is too short
    if (static_cast<long>(reference.size()) < sequenceLength + strandPosition + BandedSmithWaterman::WIDEST_GAP_SIZE)
    {
        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragmentMetadata.getCluster().getId(), "alignGapped: reference too short!");
        return 0;
    }
    // find appropriate beginning and end for the database
    const std::pair<unsigned, unsigned> flanks = getFlanks(strandPosition, sequenceLength, reference.size(), BandedSmithWaterman::WIDEST_GAP_SIZE);
    assert(flanks.first + flanks.second == BandedSmithWaterman::WIDEST_GAP_SIZE - 1);
    assert(flanks.first <= strandPosition);
    assert(strandPosition + sequenceLength + flanks.second <= (long)reference.size());
    const std::vector<char>::const_iterator databaseBegin = reference.begin() + strandPosition - flanks.first;
    const std::vector<char>::const_iterator databaseEnd = databaseBegin + flanks.first + sequenceLength + flanks.second;

    if (smartSmithWaterman_ && !makesSenseToGapAlign(
        fragmentMetadata.getCluster().getTile(), fragmentMetadata.getCluster().getId(),
        fragmentMetadata.getReadIndex(), fragmentMetadata.isReverse(),
        sequenceBegin, sequenceEnd, databaseBegin, databaseEnd))
    {
        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragmentMetadata.getCluster().getId(), "Gap-aligning does not make sense" << common::makeFastIoString(sequenceBegin, sequenceEnd) <<
            " against " << common::makeFastIoString(databaseBegin, databaseEnd));
        return 0;
    }


    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragmentMetadata.getCluster().getId(), "Gap-aligning " << common::makeFastIoString(sequenceBegin, sequenceEnd) <<
        " against " << common::makeFastIoString(databaseBegin, databaseEnd) << " strandPosition:"<<strandPosition);
    strandPosition += bandedSmithWaterman_.align(sequenceBegin, sequenceEnd, databaseBegin, databaseEnd, cigarBuffer);

    if (firstMappedBaseOffset)
    {
        const Cigar::Component firstComponent = Cigar::decode(cigarBuffer.at(cigarOffset + 1));
        // avoid two softclips in a row
        if (Cigar::SOFT_CLIP == firstComponent.second)
        {
            cigarBuffer.erase(cigarBuffer.begin() + cigarOffset);
            cigarBuffer.updateOperation(cigarOffset, firstMappedBaseOffset + firstComponent.first, Cigar::SOFT_CLIP);
        }
    }

    const unsigned clipEndBases = std::distance(sequenceEnd, sequence.end());
    if (clipEndBases)
    {
        const Cigar::Component lastComponent = Cigar::decode(cigarBuffer.back());
        if (Cigar::SOFT_CLIP == lastComponent.second)
        {
            cigarBuffer.updateOperation(cigarBuffer.size() - 1, clipEndBases + lastComponent.first, Cigar::SOFT_CLIP);
        }
        else
        {
            cigarBuffer.addOperation(clipEndBases, Cigar::SOFT_CLIP);
        }
    }

    // adjust the start position of the fragment
    strandPosition -= flanks.first;

    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragmentMetadata.getCluster().getId(), "gapped CIGAR: " <<
                                           alignment::Cigar::toString(cigarBuffer.begin() + cigarOffset, cigarBuffer.end()) << " strandPosition:"<<strandPosition);

    const unsigned matchCount = updateFragmentCigar(readMetadataList, contigList, kUniqenessAnnotation, fragmentMetadata,
                                                    fragmentMetadata.contigId, strandPosition, cigarBuffer, cigarOffset);

    return matchCount;
}

void GappedAligner::realignBadUngappedAlignments(
    const unsigned gappedMismatchesMax,
    const unsigned smitWatermanGapsMax,
    const reference::ContigList &contigList,
    const isaac::reference::ContigAnnotations &kUniqenessAnnotation,
    const flowcell::ReadMetadataList &readMetadataList,
    FragmentMetadataList &fragmentList,
    matchSelector::FragmentSequencingAdapterClipper &adapterClipper,
    Cigar &cigarBuffer)
{
    BOOST_FOREACH(FragmentMetadata &fragmentMetadata, fragmentList)
    {
        const unsigned contigId = fragmentMetadata.contigId;

        adapterClipper.checkInitStrand(fragmentMetadata, contigList[contigId]);
        ISAAC_THREAD_CERR_DEV_TRACE("    Original    : " << fragmentMetadata);
        if (BandedSmithWaterman::mismatchesCutoff < fragmentMetadata.mismatchCount)
        {
            FragmentMetadata tmp = fragmentMetadata;
            const unsigned matchCount = alignGapped(tmp, cigarBuffer, readMetadataList, adapterClipper, contigList, kUniqenessAnnotation);
            ISAAC_THREAD_CERR_DEV_TRACE("    Gap-aligned: " << tmp);
//            if (matchCount && matchCount + BandedSmithWaterman::WIDEST_GAP_SIZE > fragmentMetadata.getObservedLength() &&
//                (tmp.mismatchCount <= gappedMismatchesMax) &&
//                (fragmentMetadata.mismatchCount > tmp.mismatchCount) &&
//                ISAAC_LP_LESS(fragmentMetadata.logProbability, tmp.logProbability))
            if (matchCount && tmp.gapCount <= smitWatermanGapsMax && (
                tmp.smithWatermanScore < fragmentMetadata.smithWatermanScore ||
                (tmp.smithWatermanScore == fragmentMetadata.smithWatermanScore &&
                    ISAAC_LP_LESS(fragmentMetadata.logProbability, tmp.logProbability))))
            {
                ISAAC_THREAD_CERR_DEV_TRACE("    Using gap-aligned: " << tmp);
                fragmentMetadata = tmp;
            }
        }
    }
}

} // namespace fragmentBuilder
} // namespace alignment
} // namespace isaac
