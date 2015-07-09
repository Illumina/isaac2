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
 ** \file UngappedAligner.cpp
 **
 ** \brief See UngappedAligner.hh
 ** 
 ** \author Roman Petrovski
 **/
#include "alignment/fragmentBuilder/UngappedAligner.hh"
#include "alignment/Mismatch.hh"

namespace isaac
{
namespace alignment
{
namespace fragmentBuilder
{

UngappedAligner::UngappedAligner(
    const AlignmentCfg &alignmentCfg)
    : AlignerBase(alignmentCfg)
{
}

unsigned UngappedAligner::alignUngapped(
    FragmentMetadata &fragmentMetadata,
    Cigar &cigarBuffer,
    const flowcell::ReadMetadataList &readMetadataList,
    const matchSelector::FragmentSequencingAdapterClipper &adapterClipper,
    const reference::ContigList &contigList,
    const isaac::reference::ContigAnnotations &contigAnnotations) const
{
    const unsigned cigarOffset = cigarBuffer.size();

// Don't reset alignment to preserve the seed-based anchors.
//    fragmentMetadata.resetAlignment();
    ISAAC_ASSERT_MSG(!fragmentMetadata.isAligned(), "alignUngapped is expected to be performend on a clean fragment");
    fragmentMetadata.resetClipping();

    const reference::Contig &contig = contigList[fragmentMetadata.contigId];

    const Read &read = fragmentMetadata.getRead();
    const bool reverse = fragmentMetadata.reverse;
    const std::vector<char> &sequence = read.getStrandSequence(reverse);
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

    const unsigned mappedBases = std::distance(sequenceBegin, sequenceEnd);
    if (mappedBases)
    {
        const Cigar::OpCode opCode = Cigar::ALIGN;
        cigarBuffer.addOperation(mappedBases, opCode);
    }

    const unsigned clipEndBases = std::distance(sequenceEnd, sequence.end());
    if (clipEndBases)
    {
        cigarBuffer.addOperation(clipEndBases, Cigar::SOFT_CLIP);
    }

    const unsigned ret = updateFragmentCigar(
        readMetadataList, contigList, contigAnnotations, fragmentMetadata,
        fragmentMetadata.contigId, fragmentMetadata.position, cigarBuffer, cigarOffset);

    if (!ret)
    {
        fragmentMetadata.setUnaligned();
    }

    return ret;
}

void UngappedAligner::alignCandidates(
    const reference::ContigList &contigList,
    const isaac::reference::ContigAnnotations &kUniqenessAnnotation,
    const flowcell::ReadMetadataList &readMetadataList,
    FragmentMetadataList &fragmentList,
    matchSelector::FragmentSequencingAdapterClipper &adapterClipper,
    Cigar &cigarBuffer)
{
    BOOST_FOREACH(FragmentMetadata &fragmentMetadata, fragmentList)
    {
        ISAAC_THREAD_CERR_DEV_TRACE("    Original   : " << fragmentMetadata);
        const unsigned contigId = fragmentMetadata.contigId;
        ISAAC_ASSERT_MSG(contigList.size() > contigId, "Unexpected contig id");
        const reference::Contig &contig = contigList[contigId];

        adapterClipper.checkInitStrand(fragmentMetadata, contig);
        alignUngapped(fragmentMetadata, cigarBuffer, readMetadataList, adapterClipper, contigList, kUniqenessAnnotation);
        ISAAC_THREAD_CERR_DEV_TRACE("    Aligned    : " << fragmentMetadata);
    }
}

} // namespace fragmentBuilder
} // namespace alignment
} // namespace isaac
