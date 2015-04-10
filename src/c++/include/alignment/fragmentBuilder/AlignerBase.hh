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
 ** \file AlignerBase.hh
 **
 ** \brief Helper functionality for various read aligners
 ** 
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_ALIGNMENT_FRAGMENT_BUILDER_ALIGNER_BASE_HH
#define iSAAC_ALIGNMENT_FRAGMENT_BUILDER_ALIGNER_BASE_HH

#include "alignment/Cigar.hh"
#include "alignment/FragmentMetadata.hh"
#include "alignment/matchSelector/FragmentSequencingAdapterClipper.hh"
#include "reference/Contig.hh"
#include "reference/KUniqueness.hh"

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
class AlignerBase: public boost::noncopyable
{
public:
    AlignerBase(const AlignmentCfg &alignmentCfg);

protected:
    const AlignmentCfg &alignmentCfg_;

    unsigned updateFragmentCigar(
        const flowcell::ReadMetadataList &readMetadataList,
        const reference::ContigList &contigList,
        const isaac::reference::ContigAnnotations &contigAnnotations,
        FragmentMetadata &fragmentMetadata,
        const unsigned contigId,
        const long strandPosition,
        const Cigar &cigarBuffer,
        const unsigned cigarOffset) const;

    static void clipReference(
        const long referenceSize,
        FragmentMetadata &fragment,
        std::vector<char>::const_iterator &sequenceBegin,
        std::vector<char>::const_iterator &sequenceEnd);

    /**
     * \brief Sets the sequence iterators according to the masking information stored in the read.
     *        Adjusts fragment.position to point at the first non-clipped base.
     *
     */
    static void clipReadMasking(
        const alignment::Read &read,
        FragmentMetadata &fragment,
        std::vector<char>::const_iterator &sequenceBegin,
        std::vector<char>::const_iterator &sequenceEnd);
};

} // namespace fragmentBuilder
} // namespace alignment
} // namespace isaac

#endif // #ifndef iSAAC_ALIGNMENT_FRAGMENT_BUILDER_ALIGNER_BASE_HH
