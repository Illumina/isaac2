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
 ** \file SplitReadAligner.hh
 **
 ** \brief Uses seed match discrepancies to detect insertions, deletions, translocations, inversions in the reads that span single event
 ** 
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_ALIGNMENT_FRAGMENT_BUILDER_SIMPLE_INDEL_ALIGNER_HH
#define iSAAC_ALIGNMENT_FRAGMENT_BUILDER_SIMPLE_INDEL_ALIGNER_HH


#include "alignment/fragmentBuilder/AlignerBase.hh"
#include "alignment/TemplateLengthStatistics.hh"

namespace isaac
{
namespace alignment
{

namespace fragmentBuilder
{

class SplitReadAligner: public AlignerBase
{
    const bool splitAlignments_;
public:
    SplitReadAligner(
        const AlignmentCfg &alignmentCfg,
        const unsigned splitAlignments);

    void alignSimpleSv(
        Cigar &cigarBuffer,
        const std::vector<reference::Contig> &contigList,
        const isaac::reference::ContigAnnotations &kUniqenessAnnotation,
        const flowcell::ReadMetadataList &readMetadataList,
        const TemplateLengthStatistics &templateLengthStatistics,
        FragmentMetadataList &fragmentList) const;

private:

    bool alignIndel(
        Cigar &cigarBuffer,
        const std::vector<reference::Contig> &contigList,
        const isaac::reference::ContigAnnotations &kUniqenessAnnotation,
        const flowcell::ReadMetadataList &readMetadataList,
        FragmentMetadata &head,
        const FragmentMetadata &tail) const;

    bool alignTranslocation(
        Cigar &cigarBuffer,
        const std::vector<reference::Contig> &contigList,
        const isaac::reference::ContigAnnotations &kUniqenessAnnotation,
        const flowcell::ReadMetadataList &readMetadataList,
        FragmentMetadata &head,
        const FragmentMetadata &tail) const;

    bool alignSimpleDeletion(
        Cigar &cigarBuffer,
        FragmentMetadata &headFragment,
        const unsigned headSeedOffset,
        const FragmentMetadata &tailFragment,
        const unsigned tailSeedOffset,
        const std::vector<reference::Contig> &contigList,
        const isaac::reference::ContigAnnotations &kUniqenessAnnotation,
        const flowcell::ReadMetadataList &readMetadataList) const;

    bool alignLeftAnchoredInversion(
        Cigar &cigarBuffer,
        FragmentMetadata &headFragment,
        const unsigned headSeedOffset,
        const FragmentMetadata &tailFragment,
        const unsigned tailSeedOffset,
        const std::vector<reference::Contig> &contigList,
        const isaac::reference::ContigAnnotations &kUniqenessAnnotation,
        const flowcell::ReadMetadataList &readMetadataList) const;

    bool alignRightAnchoredInversion(
        Cigar &cigarBuffer,
        FragmentMetadata &headFragment,
        const unsigned headSeedOffset,
        const FragmentMetadata &tailFragment,
        const unsigned tailSeedOffset,
        const std::vector<reference::Contig> &contigList,
        const isaac::reference::ContigAnnotations &kUniqenessAnnotation,
        const flowcell::ReadMetadataList &readMetadataList) const;

    bool alignSimpleInsertion(
        Cigar &cigarBuffer,
        const FragmentMetadata &headAlignment,
        const unsigned headSeedOffset,
        const unsigned headSeedLength,
        FragmentMetadata &tailAlignment,
        const unsigned tailSeedOffset,
        const unsigned tailSeedLength,
        const std::vector<reference::Contig> &contigList,
        const isaac::reference::ContigAnnotations &kUniqenessAnnotation,
        const flowcell::ReadMetadataList &readMetadataList) const;

    bool mergeDeletionAlignments(
        Cigar &cigarBuffer,
        FragmentMetadata &headAlignment,
        const FragmentMetadata &tailAlignment,
        const unsigned bestOffset,
        const reference::ContigList &contigList,
        const isaac::reference::ContigAnnotations &kUniqenessAnnotation,
        const unsigned bestMismatches,
        const int deletionLength,
        const flowcell::ReadMetadataList &readMetadataList) const;

    bool mergeInversionAlignments(
        Cigar &cigarBuffer,
        FragmentMetadata &headAlignment,
        const FragmentMetadata &tailAlignment,
        const unsigned bestOffset,
        const reference::ContigList &contigList,
        const isaac::reference::ContigAnnotations &kUniqenessAnnotation,
        const unsigned bestMismatches,
        const int deletionLength,
        const flowcell::ReadMetadataList &readMetadataList) const;

    bool mergeRightAnchoredInversionAlignments(
        Cigar &cigarBuffer,
        FragmentMetadata &headAlignment,
        const FragmentMetadata &tailAlignment,
        const unsigned bestOffset,
        const reference::ContigList &contigList,
        const isaac::reference::ContigAnnotations &kUniqenessAnnotation,
        const unsigned bestMismatches,
        const int deletionLength,
        const flowcell::ReadMetadataList &readMetadataList) const;

    bool mergeInsertionAlignments(
        Cigar &cigarBuffer,
        const FragmentMetadata &headAlignment,
        FragmentMetadata &tailAlignment,
        const unsigned bestOffset,
        const reference::ContigList &contigList,
        const isaac::reference::ContigAnnotations &kUniqenessAnnotation,
        const unsigned bestMismatches,
        const unsigned insertionLength,
        const flowcell::ReadMetadataList &readMetadataList) const;

};

} // namespace fragmentBuilder
} // namespace alignment
} // namespace isaac

#endif // #ifndef iSAAC_ALIGNMENT_FRAGMENT_BUILDER_SIMPLE_INDEL_ALIGNER_HH
