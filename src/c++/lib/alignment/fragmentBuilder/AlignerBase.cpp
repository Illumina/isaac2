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
 ** \file AlignerBase.cpp
 **
 ** \brief See AlignerBase.hh
 ** 
 ** \author Roman Petrovski
 **/
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>

#include "alignment/fragmentBuilder/AlignerBase.hh"
#include "alignment/Mismatch.hh"

//#pragma GCC push_options
//#pragma GCC optimize ("0")


namespace isaac
{
namespace alignment
{
namespace fragmentBuilder
{

AlignerBase::AlignerBase(const AlignmentCfg &alignmentCfg)
    : alignmentCfg_(alignmentCfg)
{
}

/**
 * \brief Adjusts the sequence iterators to stay within the reference. Adjusts sequenceBeginReferencePosition
 *        to point at the first not clipped base.
 *
 */
void AlignerBase::clipReference(
    const long referenceSize,
    FragmentMetadata &fragment,
    std::vector<char>::const_iterator &sequenceBegin,
    std::vector<char>::const_iterator &sequenceEnd)
{
    const long referenceLeft = referenceSize - fragment.position;
    if (referenceLeft >= 0)
    {
        if (referenceLeft < std::distance(sequenceBegin, sequenceEnd))
        {
            sequenceEnd = sequenceBegin + referenceLeft;
        }

        if (0 > fragment.position)
        {
            sequenceBegin -= fragment.position;
            fragment.position = 0L;
        }

        // in some cases other clipping can end the sequence before the reference even begins
        // or begin after it ends...
        sequenceEnd = std::max(sequenceEnd, sequenceBegin);
    }
    else
    {
        // the picard sam ValidateSamFile does not like it when alignment position points to the next base after the end of the contig.
        fragment.position += referenceLeft - 1;
        sequenceBegin += referenceLeft - 1;
        --sequenceBegin;
        sequenceEnd = sequenceBegin;
    }
}

/**
 * \brief Sets the sequence iterators according to the masking information stored in the read.
 *        Adjusts fragment.position to point at the first non-clipped base.
 *
 */
void AlignerBase::clipReadMasking(
    const alignment::Read &read,
    FragmentMetadata &fragment,
    std::vector<char>::const_iterator &sequenceBegin,
    std::vector<char>::const_iterator &sequenceEnd)
{
    std::vector<char>::const_iterator maskedBegin;
    std::vector<char>::const_iterator maskedEnd;
    if (fragment.reverse)
    {
        maskedBegin = read.getReverseSequence().begin() + read.getEndCyclesMasked();
        maskedEnd = read.getReverseSequence().end() - read.getBeginCyclesMasked();
    }
    else
    {
        maskedBegin = read.getForwardSequence().begin() + read.getBeginCyclesMasked();
        maskedEnd = read.getForwardSequence().end() - read.getEndCyclesMasked();
    }

    if (maskedBegin > sequenceBegin)
    {
        fragment.incrementClipLeft(std::distance(sequenceBegin, maskedBegin));
        sequenceBegin = maskedBegin;
    }

    if (maskedEnd < sequenceEnd)
    {
        fragment.incrementClipRight(std::distance(maskedEnd, sequenceEnd));
        sequenceEnd = maskedEnd;
    }
}

unsigned AlignerBase::updateFragmentCigar(
    const flowcell::ReadMetadataList &readMetadataList,
    const reference::ContigList &contigList,
    const isaac::reference::ContigAnnotations &contigAnnotations,
    FragmentMetadata &fragmentMetadata,
    unsigned contigId,
    const long strandPosition,
    const Cigar &cigarBuffer,
    const unsigned cigarOffset) const
{
    return fragmentMetadata.updateAlignment(
        alignmentCfg_,
        readMetadataList,
        contigList, contigAnnotations,
        contigId, strandPosition,
        cigarBuffer, cigarOffset);
}

} // namespace fragmentBuilder
} // namespace alignment
} // namespace isaac

//#pragma GCC pop_options

