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
 ** \file FragmentAccessorBamAdapter.hh
 **
 ** Generates bam SA tag out of iSAAC internal CIGAR read and reference
 ** 
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_BUILD_SA_TAG_MAKER_HH
#define iSAAC_BUILD_SA_TAG_MAKER_HH

#include "alignment/Cigar.hh"
#include "alignment/Mismatch.hh"
#include "common/FastIo.hh"
#include "io/Fragment.hh"
#include "reference/Contig.hh"
#include "reference/ReferencePosition.hh"

namespace isaac
{
namespace build
{

template<typename ContainerT>
std::size_t beginSaPart(const isaac::reference::ReferencePosition& pos,
                 const bool reverse,
                 const isaac::reference::ContigList& contigList,
                 ContainerT& result)
{
    //    beginSAPart(pos, reverse, contigList);
    const std::string& contigName = contigList.at(pos.getContigId()).name_;
    const std::size_t ret = result.size();
    std::copy(contigName.begin(), contigName.end(), std::back_inserter(result));
    result.push_back(',');
    common::appendUnsignedNumber(result, pos.getPosition() + 1);
    result.push_back(',');
    result.push_back(reverse ? '-' : '+');
    result.push_back(',');
//    ISAAC_THREAD_CERR << "beginSaPart:" << std::string(result.begin(), result.end()) << std::endl;
    return ret;
}

template<typename IteratorT, typename ContainerT>
void endSaPart(IteratorT compoundCigarBegin, IteratorT compoundCigarEnd,
                         const unsigned softClipBack, const unsigned mapQ,
                         const unsigned currentEditDistance, ContainerT& result)
{
    alignment::Cigar::toString(compoundCigarBegin, compoundCigarEnd, result);
    if (softClipBack)
    {
        common::appendUnsignedNumber(result, softClipBack);
        result.push_back('S');
    }
    result.push_back(',');
    common::appendUnsignedNumber(result, mapQ);
    result.push_back(',');
    common::appendUnsignedNumber(result, currentEditDistance);
    result.push_back(';');
}

/**
 * \param sequenceReverse direction of sequence
 * \param alignedReverse direction of this alignment
 */
inline unsigned countMismatches(
    const reference::ContigList& contigList,
    const bool sequenceReverse,
    const unsigned char* alignedBegin,
    const unsigned char* alignedRBegin,
    const bool alignedReverse,
    const unsigned alignedBases,
    const reference::ReferencePosition& pos)
{
    if (sequenceReverse == alignedReverse)
    {
        return alignment::countEditDistanceMismatches(
            contigList,
            alignedBegin,
            pos,
            alignedBases);
    }
    else
    {
        return alignment::countEditDistanceMismatches(
            contigList,
            boost::make_transform_iterator(boost::make_reverse_iterator(alignedRBegin), &oligo::getReverseBcl),
            pos,
            alignedBases);
    }
}

template<typename IteratorT>
unsigned computeEditDistance(
    const isaac::reference::ContigList& contigList,
    const alignment::CigarPosition<IteratorT>& last,
    const alignment::CigarPosition<IteratorT>& current,
    const bool sequenceReverse,
    const unsigned char* sequenceBegin,
    const unsigned char* sequenceEnd)
{
    const unsigned sequenceMoved = current.sequenceOffset_ - last.sequenceOffset_;
    const unsigned referenceMoved = current.referencePos_ - last.referencePos_;
    if (sequenceMoved == referenceMoved)
    {
        return countMismatches(
            contigList,
            sequenceReverse,
            sequenceBegin + last.sequenceOffset_, sequenceEnd - last.sequenceOffset_,
            current.reverse_,
            sequenceMoved,
            last.referencePos_);
    }
    else
    {
        ISAAC_ASSERT_MSG(
            !sequenceMoved || !referenceMoved,
            "Unexpected unequal transition in CIGAR between reference and sequence");
        if (alignment::Cigar::SOFT_CLIP
            != alignment::Cigar::decode(*last.cigarIt_).second)
        {
            return sequenceMoved + referenceMoved;
        }
    }
    return 0;
}

template<typename ExcludeCigarIteratorT, typename CurrentCigarIteratorT>
bool notTheExcludeAlignment(
    const reference::ReferencePosition& excludePos, const bool excludeReverse,
    const ExcludeCigarIteratorT excludeCigarBegin,
    const ExcludeCigarIteratorT excludeCigarEnd,
    alignment::Cigar::value_type excludeCigarAutoSoftClipComponent,
    const reference::ReferencePosition& currentPos, const bool currentReverse,
    const CurrentCigarIteratorT currentCigarBegin,
    const CurrentCigarIteratorT currentCigarEnd,
    const unsigned int autoSoftClip)
{
    const std::size_t excludeCigarLength = std::distance(excludeCigarBegin, excludeCigarEnd);
    if (excludePos != currentPos || excludeReverse != currentReverse)
    {
        return true;
    }
    const unsigned currentCigarLength = std::distance(currentCigarBegin, currentCigarEnd);
    if (!autoSoftClip)
    {
        return currentCigarLength != excludeCigarLength ||
            currentCigarEnd != std::mismatch(currentCigarBegin, currentCigarEnd, excludeCigarBegin).first;
    }
    else if (currentCigarLength + 1 != excludeCigarLength)
    {
        return true;
    }

    const alignment::Cigar::Component decoded = alignment::Cigar::decode(excludeCigarAutoSoftClipComponent);
    return alignment::Cigar::SOFT_CLIP != decoded.second || autoSoftClip != decoded.first;
}

/**
 * \brief Serializes cigar to a container in bam SA tag format
 *
 * \return true if the alignment had any splits in it
 */
template <typename CigarIteratorT, typename ContainerT>
static bool makeSaTagString(
    const io::FragmentAccessor &fragment,
    const reference::ReferencePosition &excludePos,
    const bool excludeReverse,
    const CigarIteratorT excludeCigarBegin, const CigarIteratorT excludeCigarEnd,
    const unsigned mapQ,
    const isaac::reference::ContigList &contigList,
    ContainerT &resultBuffer,
    const unsigned splitGapLength,
    unsigned &excludeEditDistance)
{
//    ISAAC_THREAD_CERR << "makeSaTagString excludeCigar:" << alignment::Cigar::toString(excludeCigarBegin, excludeCigarEnd) << std::endl;

    excludeEditDistance = -1U;
    std::size_t currentBufferOffset = beginSaPart(fragment.getFStrandReferencePosition(), fragment.isReverse(), contigList, resultBuffer);

    // First component in the compound cigar to denote the current alignment
    alignment::CigarPosition<const unsigned *> currentBegin(
        fragment.cigarBegin(), fragment.cigarEnd(),
        fragment.getFStrandReferencePosition(), fragment.isReverse(), fragment.readLength_);
    alignment::CigarPosition<const unsigned *> current = currentBegin;
    // Last processed component in thecompound cigar
    alignment::CigarPosition<const unsigned *> last = current;
    bool splitsDetected = false;
    unsigned currentEditDistance = 0;
    unsigned lastFrontAutoSoftClip = 0;
    for (; !current.end(); ++current)
    {
        if (current.referencePos_.getContigId() != last.referencePos_.getContigId() ||
            current.referencePos_ < last.referencePos_ ||
            current.reverse_ != last.reverse_ ||
            ((current.sequenceOffset_ - last.sequenceOffset_) != (current.referencePos_ - last.referencePos_) &&
                (current.referencePos_ - last.referencePos_) > splitGapLength))
        {
            ISAAC_ASSERT_MSG(!splitsDetected, "Only one split per CIGAR is supported at the moment: " <<
                             alignment::Cigar::toString(fragment.cigarBegin(), fragment.cigarEnd()) << " " <<
                             oligo::bclToString(fragment.basesBegin(), fragment.readLength_));

            if (notTheExcludeAlignment(
                excludePos, excludeReverse, excludeCigarBegin, excludeCigarEnd, *(excludeCigarEnd - 1),
                currentBegin.referencePos_, currentBegin.reverse_,
                currentBegin.cigarIt_, last.cigarIt_, fragment.readLength_ - last.sequenceOffset_))
            {
                endSaPart(currentBegin.cigarIt_, last.cigarIt_,
                          fragment.readLength_ - last.sequenceOffset_, mapQ, currentEditDistance, resultBuffer);
                splitsDetected = true;
            }
            else
            {
                // exclude current SA component from the tag.
                ISAAC_ASSERT_MSG(-1U == excludeEditDistance, "Excluding more than once:" << fragment);
                excludeEditDistance = currentEditDistance;
                resultBuffer.resize(currentBufferOffset);
            }
            currentEditDistance = 0;

            currentBufferOffset = beginSaPart(current.referencePos_, current.reverse_, contigList, resultBuffer);
            if (current.reverse_ == last.reverse_)
            {
                common::appendUnsignedNumber(resultBuffer, last.sequenceOffset_);
                resultBuffer.push_back('S');
                lastFrontAutoSoftClip = last.sequenceOffset_;
            }
            else
            {
                lastFrontAutoSoftClip = 0;
            }
            currentBegin = current;
        }
        else
        {
            currentEditDistance += computeEditDistance(
                contigList, last, current, fragment.isReverse(), fragment.basesBegin(), fragment.basesBegin() + fragment.readLength_);

        }
        last = current;
    }

    currentEditDistance += computeEditDistance(
        contigList, last, current, fragment.isReverse(), fragment.basesBegin(), fragment.basesBegin() + fragment.readLength_);
    if (notTheExcludeAlignment(
        excludePos, excludeReverse, excludeCigarBegin, excludeCigarEnd, *excludeCigarBegin,
        currentBegin.referencePos_, currentBegin.reverse_, currentBegin.cigarIt_, current.cigarIt_, lastFrontAutoSoftClip))
    {

        endSaPart(currentBegin.cigarIt_, current.cigarIt_,
                  // compound cigar must cover all the sequence bases to the end. Automatic soft clipping
                  // is only performed at the points where compound cigar gets broken up.
                  0,//current.reverse_ != prevEnd.reverse_ ?  prevEnd.sequenceOffset_ : 0,
                  mapQ, currentEditDistance, resultBuffer);
        splitsDetected = true;
    }
    else
    {
        // exclude current SA component from the tag.
        ISAAC_ASSERT_MSG(-1U == excludeEditDistance, "Excluding more than once:" << fragment);
        excludeEditDistance = currentEditDistance;
        resultBuffer.resize(currentBufferOffset);
    }



    return splitsDetected;
}

} // namespace build
} // namespace isaac

#endif // #ifndef iSAAC_BUILD_SA_TAG_MAKER_HH
