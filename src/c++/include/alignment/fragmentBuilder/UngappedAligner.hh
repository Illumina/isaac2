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
 ** \file UngappedAligner.hh
 **
 ** \brief Aligns read to the reference without gaps
 ** 
 ** \author Come Raczy
 **/

#ifndef iSAAC_ALIGNMENT_FRAGMENT_BUILDER_UNGAPPED_ALIGNER_HH
#define iSAAC_ALIGNMENT_FRAGMENT_BUILDER_UNGAPPED_ALIGNER_HH

#include "alignment/fragmentBuilder/AlignerBase.hh"

namespace isaac
{
namespace alignment
{

namespace fragmentBuilder
{

class Cluster;

/**
 ** \brief Utility component creating and scoring all Fragment instances from a
 ** list Seed Matches for a single Cluster (each Read independently).
 **/
class UngappedAligner: AlignerBase
{
public:
    UngappedAligner(
        const AlignmentCfg &alignmentCfg);

    void alignCandidates(
        const reference::ContigList &contigList,
        const isaac::reference::ContigAnnotations &kUniqenessAnnotation,
        const flowcell::ReadMetadataList &readMetadataList,
        FragmentMetadataList &fragmentList,
        matchSelector::FragmentSequencingAdapterClipper &adapterClipper,
        Cigar &cigarBuffer);

    /**
     ** \brief Calculate the ungapped alignment of a fragment
     **
     ** \return number of mapped matches. If 0, read is marked as unaligned
     **/
    unsigned alignUngapped(
        FragmentMetadata &fragmentMetadata,
        Cigar &cigarBuffer,
        const flowcell::ReadMetadataList &readMetadataList,
        const matchSelector::FragmentSequencingAdapterClipper  &adapterClipper,
        const reference::ContigList &contigList,
        const isaac::reference::ContigAnnotations &contigAnnotations) const;

private:
};

} // namespace fragmentBuilder
} // namespace alignment
} // namespace isaac

#endif // #ifndef iSAAC_ALIGNMENT_FRAGMENT_BUILDER_UNGAPPED_ALIGNER_HH
