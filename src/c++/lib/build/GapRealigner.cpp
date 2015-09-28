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
 ** \file GapRealigner.cpp
 **
 ** Attempts to reduce read mismatches by introducing gaps found on other reads.
 **
 ** \author Roman Petrovski
 **/
#include <boost/foreach.hpp>
#include <boost/format.hpp>

#include "alignment/BandedSmithWaterman.hh"
#include "build/GapRealigner.hh"
#include "build/gapRealigner/ChooseKGapsFilter.hh"

#include "SemialignedEndsClipper.hh"

namespace isaac
{
namespace build
{
inline std::ostream & operator << (std::ostream &os, const GapRealigner::RealignmentBounds &fragmentGaps)
{
    return os << "FragmentGaps(" <<
        fragmentGaps.beginPos_ << "," <<
        fragmentGaps.firstGapStartPos_ << "," <<
        fragmentGaps.lastGapEndPos_ << "," <<
        fragmentGaps.endPos_ << ")";
}

const GapRealigner::RealignmentBounds GapRealigner::extractRealignmentBounds(
    const PackedFragmentBuffer::Index &index) const
{
    using alignment::Cigar;
    RealignmentBounds ret = {index.pos_, index.pos_, index.pos_, index.pos_};

    const unsigned* cigarIterator = index.cigarBegin_;
    for(;index.cigarEnd_ != cigarIterator; ++cigarIterator)
    {
        const Cigar::Component decoded = Cigar::decode(*cigarIterator);
        if (decoded.second == Cigar::ALIGN)
        {
            ret.firstGapStartPos_ += decoded.first;
            ret.endPos_ += decoded.first;
        }
        else if (decoded.second == Cigar::INSERT)
        {
            ret.lastGapEndPos_ = ret.endPos_;
            ++cigarIterator;
            break;
        }
        else if (decoded.second == Cigar::DELETE)
        {
            ret.endPos_ += decoded.first;
            ret.lastGapEndPos_ = ret.endPos_;
            ++cigarIterator;
            break;
        }
        else if (decoded.second == Cigar::SOFT_CLIP)
        {
            if (index.cigarBegin_ == cigarIterator)
            {
                ISAAC_ASSERT_MSG(ret.firstGapStartPos_ == ret.beginPos_, "first soft clip must happen before anything");
                ISAAC_ASSERT_MSG(ret.lastGapEndPos_ == ret.beginPos_, "first soft clip must happen before anything");
                ISAAC_ASSERT_MSG(ret.endPos_ == ret.beginPos_, "first soft clip must happen before anything");
                ret.lastGapEndPos_ = ret.firstGapStartPos_ = (ret.beginPos_ -= decoded.first);
            }
            else
            {
                ISAAC_ASSERT_MSG(index.cigarEnd_ == cigarIterator + 1,
                                 "At most two soft-clips are expected with second one being the last component of the cigar");
                ret.endPos_ += decoded.first;
            }
        }
        else
        {
            const boost::format message = boost::format("Unexpected Cigar OpCode: %d") % decoded.second;
            ISAAC_ASSERT_MSG(false, message.str().c_str());
        }
    }

    for(;index.cigarEnd_ != cigarIterator; ++cigarIterator)
    {
        const Cigar::Component decoded = Cigar::decode(*cigarIterator);
        if (decoded.second == Cigar::ALIGN)
        {
            ret.endPos_ += decoded.first;
        }
        else if (decoded.second == Cigar::INSERT)
        {
            ret.lastGapEndPos_ = ret.endPos_;
        }
        else if (decoded.second == Cigar::DELETE)
        {
            ret.endPos_ += decoded.first;
            ret.lastGapEndPos_ = ret.endPos_;
        }
        else if (decoded.second == Cigar::SOFT_CLIP)
        {
            ISAAC_ASSERT_MSG(index.cigarEnd_ == cigarIterator + 1,
                             "Last soft-clip has to be the last component of the cigar");
            ret.endPos_ += decoded.first;
        }
        else
        {
            const boost::format message = boost::format("Unexpected Cigar OpCode: %d") % decoded.second;
            ISAAC_ASSERT_MSG(false, message.str().c_str());
        }
    }

    ISAAC_ASSERT_MSG(ret.endPos_ >= ret.lastGapEndPos_, "Fragment gap cannot end outside the fragment");
    return ret;
}

/**
 * \brief Adjusts template length and mate fields
 *
 * \param beginPosShift change in the fragment alignment position. negative for insertions, positive for deletions
 * \param endPosShift   change in the fragment end alignment position. negative for insertions, positive for deletions
 */
void GapRealigner::updatePairDetails(
    const std::vector<alignment::TemplateLengthStatistics> &barcodeTemplateLengthStatistics,
    const PackedFragmentBuffer::Index &index,
    io::FragmentAccessor &fragment,
    PackedFragmentBuffer &dataBuffer)
{
    fragment.flags_.realigned_ = true;
    if (!index.hasMate() || fragment.flags_.mateUnmapped_)
    {
        const int oldBamTlen = fragment.bamTlen_;
        fragment.bamTlen_ = fragment.bamTlen_ < 0 ?
            (-fragment.observedLength_ + 1) : (fragment.observedLength_ - 1);
        if (index.hasMate())
        {
            io::FragmentAccessor &mate  = dataBuffer.getMate(index);
            ISAAC_ASSERT_MSG(mate.bamTlen_ == -oldBamTlen,
                             "Incoherent BAM template length between fragment and unmapped mate " << index <<
                             " " << fragment << " mate " << mate);
            mate.bamTlen_ = -fragment.bamTlen_;
            fragment.mateFStrandPosition_ = fragment.fStrandPosition_;
            mate.mateFStrandPosition_ = fragment.fStrandPosition_;
            mate.fStrandPosition_ = fragment.fStrandPosition_;
        }
        return;
    }

    io::FragmentAccessor &mate  = dataBuffer.getMate(index);

    ISAAC_ASSERT_MSG(mate.bamTlen_ == -fragment.bamTlen_,
                     "Incoherent BAM template length between fragment and mate " << index <<
                     " " << fragment <<
                     " mate " << mate);

    const reference::ReferencePosition fragmentBeginPos = fragment.fStrandPosition_;
    const reference::ReferencePosition fragmentEndPos = fragmentBeginPos + fragment.observedLength_;
    const reference::ReferencePosition mateBeginPos = fragment.mateFStrandPosition_;
    const reference::ReferencePosition mateEndPos = mateBeginPos + mate.observedLength_;

    fragment.bamTlen_ = io::FragmentAccessor::getTlen(fragmentBeginPos, fragmentEndPos, mateBeginPos, mateEndPos, !fragment.flags_.secondRead_);
    mate.bamTlen_ = -fragment.bamTlen_;
    mate.mateFStrandPosition_ = fragment.fStrandPosition_;
    mate.flags_.properPair_ = fragment.flags_.properPair_ =
        alignment::TemplateLengthStatistics::Nominal ==
            barcodeTemplateLengthStatistics.at(fragment.barcode_).checkModel(fragment, mate);
}

/**
 * \breif Collapses gaps on the ends of CIGAR into soft-clips
 *        This keeps CASAVA happy and probably in general is a right thing to do.
 *        Also adjusts fragment position and observed length.
 *
 * \return false, when the compacting the cigar would produce invalid fragment. Currently two cases are considered invalid:
 *          a) move the read past the binEndPos
 *          b) fragment gets completely soft-clipped
 *         If false is returned index and fragment are guaranteed to be unchanged.
 */
bool GapRealigner::compactCigar(
    const std::vector<reference::Contig> &reference,
    const reference::ReferencePosition binEndPos,
    PackedFragmentBuffer::Index &index,
    io::FragmentAccessor &fragment,
    alignment::Cigar &realignedCigars)
{
    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment.clusterId_, " Compacting " << fragment << index);

    ISAAC_ASSERT_MSG(alignment::Cigar::getReadLength(index.cigarBegin_, index.cigarEnd_) == fragment.readLength_,
                     "Broken CIGAR before compacting: " << alignment::Cigar::toString(index.cigarBegin_, index.cigarEnd_) <<
                     " " << index << " " << fragment);

    // at this point the fragment.editDistance_ is not adjusted for soft clipped bases
    using alignment::Cigar;
    const uint32_t *cigarIterator = index.cigarBegin_;
    unsigned softClipStart = 0;
    bool needCompacting = false;
    reference::ReferencePosition newPos = index.pos_;
    for (; index.cigarEnd_ != cigarIterator; ++cigarIterator)
    {
        Cigar::Component decoded = Cigar::decode(*cigarIterator);
        if (Cigar::ALIGN == decoded.second)
        {
            break;
        }
        else if (Cigar::SOFT_CLIP == decoded.second)
        {
            if (index.cigarBegin_ != cigarIterator)
            {
                ISAAC_ASSERT_MSG(index.cigarEnd_ == cigarIterator + 1,
                                 "At most two soft-clips are expected with second one being the last component of the cigar");
            }
            softClipStart += decoded.first;
        }
        else if (Cigar::INSERT == decoded.second)
        {
            needCompacting = true;
            softClipStart += decoded.first;
        }
        else if (Cigar::DELETE == decoded.second)
        {
            needCompacting = true;
            if (binEndPos <= newPos + decoded.first)
            {
                ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment.clusterId_, " Cannot compact CIGAR as read will move to next bins " << fragment << index);
                return false;
            }
            newPos += decoded.first;
        }
        else
        {
            ISAAC_ASSERT_MSG(false, "Unexpected CIGAR operation");
        }
    }

    if (index.cigarEnd_ == cigarIterator)
    {
        // fragment gets completely soft-clipped
        return false;
    }

    const uint32_t *cigarBackwardsIterator = index.cigarEnd_ - 1;
    unsigned softClipEnd = 0;
    for (;cigarIterator != cigarBackwardsIterator; --cigarBackwardsIterator)
    {
        Cigar::Component decoded = Cigar::decode(*cigarBackwardsIterator);
        if (Cigar::ALIGN == decoded.second)
        {
            break;
        }
        else if (Cigar::SOFT_CLIP == decoded.second)
        {
            ISAAC_ASSERT_MSG(index.cigarEnd_ == cigarBackwardsIterator + 1,
                             "At most two soft-clips are expected with second one being the last component of the cigar" << index);
            softClipEnd += decoded.first;
        }
        else if (Cigar::INSERT == decoded.second)
        {
            needCompacting = true;
            softClipEnd += decoded.first;
        }
        else if (Cigar::DELETE == decoded.second)
        {
            needCompacting = true;
            // nothing to be done for trailing deletes;
        }
        else
        {
            ISAAC_ASSERT_MSG(false, "Unexpected CIGAR operation");
        }
    }

    if (needCompacting)
    {
        std::size_t before = realignedCigars.size();

        if (softClipStart)
        {
            realignedCigars.addOperation(softClipStart, Cigar::SOFT_CLIP);
        }
        ISAAC_ASSERT_MSG(std::distance(cigarIterator, cigarBackwardsIterator + 1) < 200 &&
                         std::distance(cigarIterator, cigarBackwardsIterator + 1) > 0,
                         "Suspiciously long (" << std::distance(cigarIterator, cigarBackwardsIterator + 1) <<
                         ") CIGAR in the middle of compacting: " <<
                         " " << alignment::Cigar::toString(index.cigarBegin_, index.cigarEnd_) <<
                         " " << index << " " << fragment);

        ISAAC_ASSERT_MSG(alignment::Cigar::getReadLength(cigarIterator, cigarBackwardsIterator + 1) + softClipStart + softClipEnd == fragment.readLength_,
                         "Broken CIGAR in the middle of compacting: " << alignment::Cigar::toString(cigarIterator, cigarBackwardsIterator + 1) <<
                         " " << index << " " << fragment);
        realignedCigars.addOperations(cigarIterator, cigarBackwardsIterator + 1);
        if (softClipEnd)
        {
            realignedCigars.addOperation(softClipEnd, Cigar::SOFT_CLIP);
        }

        index.cigarBegin_ = &realignedCigars.at(before);
        index.cigarEnd_ = &realignedCigars.back() + 1;
        index.pos_ = newPos;
    }

    // recompute editDistance and observed length
    unsigned short newEditDistance = 0;
    const unsigned char *basesIterator = fragment.basesBegin() + softClipStart;
    std::vector<char>::const_iterator referenceIterator =
        reference.at(index.pos_.getContigId()).forward_.begin() + index.pos_.getPosition();

    reference::ReferencePosition newEndPos = index.pos_;
    for (;cigarIterator != cigarBackwardsIterator + 1; ++cigarIterator)
    {
        Cigar::Component decoded = Cigar::decode(*cigarIterator);
        if (Cigar::ALIGN == decoded.second)
        {
            newEditDistance += alignment::countEditDistanceMismatches(reference, basesIterator, newEndPos, decoded.first);
            newEndPos += decoded.first;
            basesIterator += decoded.first;
            referenceIterator += decoded.first;
        }
        else if (Cigar::SOFT_CLIP == decoded.second)
        {
            ISAAC_ASSERT_MSG(false, "At most two soft-clips are expected. Both at the ends of the CIGAR " << index << " " << fragment);
        }
        else if (Cigar::INSERT == decoded.second)
        {
            newEditDistance += decoded.first;
            basesIterator += decoded.first;
        }
        else if (Cigar::DELETE == decoded.second)
        {
            newEditDistance += decoded.first;
            newEndPos += decoded.first;
            referenceIterator += decoded.first;
        }
        else
        {
            ISAAC_ASSERT_MSG(false, "Unexpected CIGAR operation");
        }
    }

    ISAAC_ASSERT_MSG(alignment::Cigar::getReadLength(index.cigarBegin_, index.cigarEnd_) == fragment.readLength_,
                      "Broken CIGAR after compacting: " << alignment::Cigar::toString(index.cigarBegin_, index.cigarEnd_) <<
                      " " << index << " " << fragment);

    fragment.editDistance_ = newEditDistance;
    fragment.fStrandPosition_ = index.pos_;
    fragment.observedLength_ = newEndPos - fragment.fStrandPosition_;

    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment.clusterId_, " Compacted " << fragment << index);

    return true;
}

inline int calculateMismatchesPercent(
    unsigned mismatches,
    unsigned mappedLength)
{
//    ISAAC_THREAD_CERR << "calculateMismatchesPercent mismatches=" << mismatches << " mappedLength=" << mappedLength << std::endl;
    return mappedLength ? mismatches * 100 / mappedLength : 100;
}


/**
 * \brief bits in choice determine whether the corresponding gaps are on or off
 *
 * \return cost of the new choice or -1U if choice is inapplicable.
 */
GapRealigner::GapChoice GapRealigner::verifyGapsChoice(
    const GapChoiceBitmask &choice,
    const gapRealigner::GapsRange &gaps,
    const reference::ReferencePosition newBeginPos,
    const io::FragmentAccessor &fragment,
    const std::vector<reference::Contig> &reference)
{
    GapChoice ret;
    ret.choice_ = choice;
    ret.startPos_ = newBeginPos;

    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment.clusterId_, "newBeginPos " << newBeginPos);

    // keeping as int to allow debug checks for running into negative
    int basesLeft = fragment.readLength_;
    int leftClippedLeft = fragment.leftClipped();

//    newBeginPos += fragment.leftClipped();
//    basesLeft -= fragment.leftClipped();

    reference::ReferencePosition lastGapEndPos = newBeginPos;
    reference::ReferencePosition lastGapBeginPos; // initially set to an invalid position which would not match any gap pos
    unsigned currentGapIndex = 0;
    BOOST_FOREACH(const gapRealigner::Gap& gap, std::make_pair(gaps.first, gaps.second))
    {
        if (choice & (GapChoiceBitmask(1) << currentGapIndex))
        {
//            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment.clusterId_, "testing " << gap);
            if (gap.getEndPos(true) <= lastGapEndPos)// || gap.getBeginPos() > lastGapEndPos + basesLeft)
            {
                // the choice requires a gap that cannot be applied.
                // just bail out. there will be another choice just like
                // this one but without the useless gap
                ret.cost_ = -1U;
                ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment.clusterId_, ret << " does not fit in range [" << lastGapEndPos << ";" << lastGapEndPos + basesLeft << ")");
                return ret;
            }
//            ISAAC_THREAD_CERR << " lastGapEndPos=" << lastGapEndPos << std::endl;

            if (gap.getBeginPos() < lastGapEndPos)
            {
                ret.cost_ = -1U;
                ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment.clusterId_, ret << " contains overlapping deletions ");
                return ret;
                // Allowing overlapping deletions is tricky because it is hard to track back the
                // newBeginPos from the pivot see findStartPos.
            }

            if (gap.getBeginPos() == lastGapBeginPos)
            {
                // The only case where it makes sense to allow two or more gaps starting at the same
                // position is when we want to combine multiple insertions into a larger
                // one. Unfortunately, with enough gaps, it consumes the read into one single insertion...
                // Other cases:
                // deletion/deletion - is disallowed above
                // insertion/deletion - does not make sense (and cause trouble SAAC-253)
                ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment.clusterId_, ret << " contains overlapping gaps ");
                ret.cost_ = -1U;
                return ret;
            }

//            ISAAC_THREAD_CERR << " lastGapEndPos=" << lastGapEndPos << " lastGapBeginPos=" << lastGapBeginPos <<
//                " gap.isInsertion()=" << gap.isInsertion() << std::endl;

            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment.clusterId_, "basesLeft - fragment.rightClipped(): " << basesLeft - fragment.rightClipped());
            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment.clusterId_, "gap.getBeginPos() - lastGapEndPos: " << gap.getBeginPos() - lastGapEndPos);
            const int mappedBases = std::min<int>(basesLeft - fragment.rightClipped(), gap.getBeginPos() - lastGapEndPos);
//            ISAAC_THREAD_CERR << " mappedBases=" << mappedBases << " basesLeft=" << basesLeft << std::endl;

            const unsigned length = mappedBases - std::min(mappedBases, leftClippedLeft);
            const unsigned mm = alignment::countEditDistanceMismatches(reference,
                                                fragment.basesBegin() + (fragment.readLength_ - basesLeft) + leftClippedLeft,
                                                lastGapEndPos + leftClippedLeft, length);

            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment.clusterId_, "length: " << length);
//            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment.clusterId_, "countMismatches: " << mm);
//            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment.clusterId_, "leftClippedLeft: " << leftClippedLeft);

            ret.mappedLength_ += length;
            ret.editDistance_ += mm;
            ret.mismatches_ += mm;
            ret.cost_ += mm * mismatchCost_;
            basesLeft -= mappedBases;
            leftClippedLeft -= std::min(leftClippedLeft, mappedBases);
//            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment.clusterId_, "leftClippedLeft: " << leftClippedLeft);
            if (!basesLeft)
            {
                // gap begins after the read ends.
                ret.cost_ = -1U;
                return ret;
            }
            unsigned clippedGapLength = 0;
            if (gap.isInsertion())
            {
                clippedGapLength = std::min<int>(basesLeft - fragment.rightClipped(), gap.getLength());
                // insertions reduce read length
                basesLeft -= clippedGapLength;
                leftClippedLeft -= std::min<int>(leftClippedLeft, gap.getLength());
            }
            else
            {
                clippedGapLength = leftClippedLeft ? 0 : gap.getLength();
            }
            ret.editDistance_ += clippedGapLength;
            ret.cost_ += clippedGapLength ? (gapOpenCost_ + (clippedGapLength - 1) * gapExtendCost_) : 0;
            lastGapEndPos = gap.getEndPos(false);
            lastGapBeginPos = gap.getBeginPos();

            if (basesLeft == leftClippedLeft + fragment.rightClipped())
            {
                break;
            }
            ISAAC_ASSERT_MSG(basesLeft > leftClippedLeft + fragment.rightClipped(), "Was not supposed to run into the clipping");

        }
        ++currentGapIndex;
    }

    if(basesLeft > leftClippedLeft + fragment.rightClipped())
    {
        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment.clusterId_, "leftClippedLeft: " << leftClippedLeft);

        const unsigned length = basesLeft - std::min<unsigned>(basesLeft, leftClippedLeft) - fragment.rightClipped();

        const reference::ReferencePosition firstUnclippedPos = lastGapEndPos + leftClippedLeft;
        if (firstUnclippedPos.getPosition() > reference.at(firstUnclippedPos.getContigId()).forward_.size())
        {
            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment.clusterId_, ret << " gap pushes part of the read outside the reference " << firstUnclippedPos << " " << basesLeft);
            ret.cost_ = -1U;
            return ret;
        }
        else
        {
            const unsigned mm = alignment::countEditDistanceMismatches(
                reference, fragment.basesBegin() + (fragment.readLength_ - basesLeft) + leftClippedLeft, firstUnclippedPos, length);
            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment.clusterId_, "final countMismatches: " << mm);
            ret.mappedLength_ += length;
            ret.editDistance_ += mm;
            ret.mismatches_ += mm;
            ret.cost_ += mm * mismatchCost_;
        }
    }
    else
    {
        ISAAC_ASSERT_MSG(leftClippedLeft + fragment.rightClipped() == basesLeft, "Spent more than readLength. basesLeft: " << basesLeft <<
            " choice: " << int(choice) <<
            " " << fragment <<
            ", newBeginPos " << newBeginPos <<
            " lastGapEndPos " << lastGapEndPos <<
            " leftClippedLeft " << leftClippedLeft <<
            " gaps: " << gaps);
    }

    ret.mismatchesPercent_ = calculateMismatchesPercent(ret.mismatches_, ret.mappedLength_);
    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment.clusterId_, "verifyGapsChoice:" << ret);
    return ret;

}

/**
 * \return false if choice cannot be applied. Currently this can happen due to left-side clipping having to move the
 *         read into the next bin
 *
 * \postcondition When false is returned, all state variables and writable inputs retain their original values.
 */
bool GapRealigner::applyChoice(
    const GapChoiceBitmask &choice,
    const gapRealigner::GapsRange &gaps,
    const reference::ReferencePosition binEndPos,
    const reference::ReferencePosition contigEndPos,
    PackedFragmentBuffer::Index &index,
    const io::FragmentAccessor &fragment,
    alignment::Cigar &realignedCigars)
{
//    ISAAC_THREAD_CERR << "GapRealigner::applyChoice index=" << index << std::endl;
    reference::ReferencePosition newBeginPos = index.pos_;
    using alignment::Cigar;

    const std::size_t before = realignedCigars.size();

    int basesLeft = fragment.readLength_;

    int leftClippedLeft = fragment.leftClipped();
    // since insertions overlapped by left clipping don't move the alignment position, we need to count the overlaps for the final alignment position adjustment
    int leftClippedInsertionBases = 0;

    if (fragment.leftClipped())
    {
        Cigar::Component decoded = Cigar::decode(*fragment.cigarBegin());
        ISAAC_ASSERT_MSG(Cigar::SOFT_CLIP == decoded.second, "Original CIGAR is expected to have soft clip at the start " << fragment);
        ISAAC_ASSERT_MSG(fragment.leftClipped() <= decoded.first, "Original CIGAR soft clip at the start is shorter than left-clipped bases" << fragment);
        realignedCigars.addOperation(fragment.leftClipped(), Cigar::SOFT_CLIP);
//        newBeginPos += fragment.leftClipped();
//        basesLeft -= fragment.leftClipped();
    }

    if (fragment.rightClipped())
    {
        Cigar::Component decoded = Cigar::decode(*(fragment.cigarEnd() - 1));
        ISAAC_ASSERT_MSG(Cigar::SOFT_CLIP == decoded.second, "Original CIGAR is expected to have soft clip at the end " << fragment);
        ISAAC_ASSERT_MSG(fragment.rightClipped() <= decoded.first, "Original CIGAR soft clip at the end is shorter than right-clipped bases" << fragment);
//        basesLeft -= fragment.rightClipped();
        // actual cigar is appended at the end of the function
    }

    reference::ReferencePosition lastGapEndPos = newBeginPos;
    unsigned currentGapIndex = 0;
    Cigar::OpCode lastOperation = Cigar::UNKNOWN;
    BOOST_FOREACH(const gapRealigner::Gap& gap, std::make_pair(gaps.first, gaps.second))
    {
        if (choice & (GapChoiceBitmask(1) << currentGapIndex))
        {
            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment.clusterId_, " Applying choice " << choice <<
                " gap index " << currentGapIndex << " " << gap << " to " << fragment);

//            ISAAC_THREAD_CERR << "lastGapEndPos=" << lastGapEndPos << std::endl;

            const reference::ReferencePosition gapClippedBeginPos = std::max(gap.getBeginPos(), newBeginPos);
//            ISAAC_THREAD_CERR << "gapClippedBeginPos=" << gapClippedBeginPos << std::endl;

            if (gapClippedBeginPos >= lastGapEndPos)
            {
                const int mappedBases = std::min<int>(basesLeft - fragment.rightClipped(), gapClippedBeginPos - lastGapEndPos);
                const unsigned softClippedMappedLength = mappedBases - std::min(mappedBases, leftClippedLeft);
//                ISAAC_THREAD_CERR << "basesLeft=" << basesLeft << std::endl;
//                ISAAC_THREAD_CERR << "lastGapEndPos=" << lastGapEndPos << std::endl;
//                ISAAC_THREAD_CERR << "mappedBases=" << mappedBases << std::endl;
//                ISAAC_THREAD_CERR << "softClippedMappedLength=" << softClippedMappedLength << std::endl;

                if (softClippedMappedLength)
                {
                    // avoid 0M in CIGAR
                    ISAAC_ASSERT_MSG(realignedCigars.capacity() > realignedCigars.size(), "Realigned CIGAR buffer is out of capacity");
                    realignedCigars.addOperation(softClippedMappedLength, Cigar::ALIGN);
                }
                basesLeft -= mappedBases;
                leftClippedLeft -= std::min(mappedBases, leftClippedLeft);

                if (gap.isInsertion())
                {
                    const int clippedGapLength =
                        std::min<int>(basesLeft - fragment.rightClipped(), gap.getEndPos(true) - gapClippedBeginPos);
                    const int softClippedGapLength = clippedGapLength - std::min(clippedGapLength, leftClippedLeft);
//                    ISAAC_THREAD_CERR << "clippedGapLength=" << clippedGapLength << std::endl;
                    if (softClippedGapLength)
                    {
                        if (alignment::Cigar::INSERT == lastOperation && !mappedBases)
                        {
                            const alignment::Cigar::Component old = alignment::Cigar::decode(realignedCigars.back());
                            realignedCigars.pop_back();
                            //GATK does not like 2I1I type cigars
                            realignedCigars.addOperation(old.first + softClippedGapLength, alignment::Cigar::INSERT);
                        }
                        else
                        {
                            realignedCigars.addOperation(softClippedGapLength, gap.getOpCode());
                            ISAAC_ASSERT_MSG(realignedCigars.capacity() > realignedCigars.size(), "Realigned CIGAR buffer is out of capacity");
                            lastOperation = alignment::Cigar::INSERT;
                        }
                    }

                    basesLeft -= clippedGapLength;
                    lastGapEndPos = gapClippedBeginPos;
                    leftClippedLeft -= std::min(clippedGapLength, leftClippedLeft);
                    leftClippedInsertionBases += clippedGapLength - softClippedGapLength;

//                    ISAAC_THREAD_CERR << "2nd lastGapEndPos=" << lastGapEndPos << std::endl;
                }
                else
                {
                    // if left-side alignment-indepndent soft-clipping is in place, the deletions that we put in will simply move the alignment position forward in compactCigar
                    const int clippedGapLength = gap.getEndPos(true) - gapClippedBeginPos;
//                    ISAAC_THREAD_CERR << "clippedGapLength=" << clippedGapLength << std::endl;
                    if (!leftClippedLeft)
                    {
                        if (alignment::Cigar::DELETE == lastOperation && !mappedBases)
                        {
                            const alignment::Cigar::Component old = alignment::Cigar::decode(realignedCigars.back());
                            realignedCigars.pop_back();
                            //assuming GATK does not like 2D1D type cigars either...
                            realignedCigars.addOperation(old.first + clippedGapLength, alignment::Cigar::DELETE);
                        }
                        else
                        {
                            ISAAC_ASSERT_MSG(realignedCigars.capacity() > realignedCigars.size(), "Realigned CIGAR buffer is out of capacity");
                            realignedCigars.addOperation(clippedGapLength, gap.getOpCode());
                            lastOperation = alignment::Cigar::DELETE;
                        }
                    }
                    else
                    {
                        newBeginPos += clippedGapLength;
                    }

                    lastGapEndPos = gap.getEndPos(false);
//                    ISAAC_THREAD_CERR << "2nd lastGapEndPos=" << lastGapEndPos << std::endl;
                }
            }
            else
            {
                ISAAC_ASSERT_MSG(false, "Overlapping gaps are not allowed");
//                unsigned appliedGapLength = gap.getEndPos(true) - lastGapEndPos;
//                if (gap.isInsertion())
//                {
//                    appliedGapLength = std::min<unsigned>(appliedGapLength, basesLeft);
//                    basesLeft -= appliedGapLength;
//                }
//                realignedCigars_.push_back(Cigar::encode(appliedGapLength, gap.getOpCode()));
//                lastOperation = gap.getOpCode();
            }

            if (basesLeft == leftClippedLeft + fragment.rightClipped())
            {
                break;
            }
            ISAAC_ASSERT_MSG(basesLeft > leftClippedLeft + fragment.rightClipped(), "Was not supposed to run into the clipping");
        }
        ++currentGapIndex;
    }
    if(basesLeft > leftClippedLeft + fragment.rightClipped())
    {
        ISAAC_ASSERT_MSG(realignedCigars.capacity() > realignedCigars.size(), "Realigned CIGAR buffer is out of capacity");
        const int basesToTheEndOfContig = contigEndPos - lastGapEndPos - leftClippedLeft;
        const int mappedBases = std::min(basesToTheEndOfContig, basesLeft - leftClippedLeft - fragment.rightClipped());
        if (mappedBases)
        {
            // avoid 0M in CIGAR
            realignedCigars.addOperation(mappedBases, Cigar::ALIGN);
        }
        basesLeft -= leftClippedLeft + mappedBases;
        leftClippedLeft = 0;
    }
    else
    {
        ISAAC_ASSERT_MSG(fragment.rightClipped() + leftClippedLeft == basesLeft, "Spent more than readLength (" << basesLeft <<
            " left) bases on the CIGAR: " << alignment::Cigar::toString(&realignedCigars.at(before), &realignedCigars.back() + 1));
    }


    if (basesLeft)
    {
        realignedCigars.addOperation(basesLeft, Cigar::SOFT_CLIP);
    }

    newBeginPos += fragment.leftClipped() - leftClippedInsertionBases;
    if (newBeginPos >= binEndPos)
    {
        // all things considered, we can't apply this choice as it would place read into the next bin.
        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment.clusterId_, "applying would move start " << newBeginPos << " past binEndPos " << binEndPos << index);
        realignedCigars.resize(before);
        return false;
    }
    index.pos_ = newBeginPos;
    index.cigarBegin_ = &realignedCigars.at(before);
    index.cigarEnd_ = &realignedCigars.back() + 1;
    return true;
}

/**
 * \brief Find the start position such that the base at pivotPos does not move when all existing fragment gaps are removed
 *
 * \return alignment start position with all existing gaps removed so that base at pivotPos stays at pivotPos
 */
long GapRealigner::undoExistingGaps(
    const PackedFragmentBuffer::Index& index,
    const reference::ReferencePosition& pivotPos)
{
    using alignment::Cigar;
    reference::ReferencePosition lastGapEndPos = index.getUnclippedPosition();
    // initialize with the distance between the (possibly clipped) read start and pivotPos,
    // then offset by the gaps and soft clipping from the original CIGAR.
    long position = index.pos_.getPosition();
    //    ISAAC_THREAD_CERR << " restoring startPos offset=" << offset << " from pivot=" << pivotPos << "and original=" << index.pos_ << std::endl;
    for (const uint32_t *it = index.cigarBegin_; index.cigarEnd_ != it; ++it)
    {
        const uint32_t cigarOp = *it;
        if (lastGapEndPos > pivotPos)
        {
            break;
        }
        const Cigar::Component decoded = Cigar::decode(cigarOp);
        if (Cigar::ALIGN == decoded.second)
        {
            lastGapEndPos += decoded.first;
        }
        else if (Cigar::INSERT == decoded.second)
        {
            position -= decoded.first;
        }
        else if (Cigar::DELETE == decoded.second)
        {
            lastGapEndPos += decoded.first;
//            if (lastGapEndPos > pivotPos)
//            {
//                // existing deletion overlaps pivot position like this:
//                //"A-----C",
//                //"AGATCAG",
//                //"   **");
//                offset = -1;
//                break;
//            }
            position += decoded.first;
        }
        else if (Cigar::SOFT_CLIP == decoded.second)
        {
            if (index.cigarBegin_ == it)
            {
                // soft clip at the start eats bases just like an insertion, but also moves the position same way as deletion does
                position -= decoded.first;
            }
            lastGapEndPos += decoded.first;
            // soft clip at the end is treated same way as mapped bases
        }
        else
        {
            ISAAC_ASSERT_MSG(false, "Unexpected CIGAR operation " << decoded.second << " in " << cigarOp);
        }
    }
    return position;
}

/**
 * \brief Find the start position such that the base that would be the read base
 *        at pivotPos if read originally had no gaps, would still be at pivotPos
 *        given all gaps in this choice are applied
 */
 

bool GapRealigner::findStartPos(
    const GapChoiceBitmask &choice,
    const gapRealigner::GapsRange &gaps,
    const reference::ReferencePosition binStartPos,
    const reference::ReferencePosition binEndPos,
    const unsigned pivotGapIndex,
    const reference::ReferencePosition pivotPos,
    long alignmentPos,
    reference::ReferencePosition &ret)
{
    // Shift found start position by the gaps included in the mask while keeping the base at pivotPos steady
    unsigned gapIndex = pivotGapIndex - 1;
    reference::ReferencePosition overlapPos = pivotPos;
    unsigned basesLeft = pivotPos.getPosition() - alignmentPos;
    BOOST_REVERSE_FOREACH(const gapRealigner::Gap& gap, std::make_pair(gaps.first, gaps.first + pivotGapIndex))
    {
        if (choice & (GapChoiceBitmask(1) << gapIndex))
        {
            if (gap.getEndPos(false) > overlapPos)
            {
//                ISAAC_THREAD_CERR << " overlapping gaps are not allowed " << overlapPos << "-" << gap << std::endl;
                // overlapping gaps are not allowed
                return false;
            }
            if (gap.isInsertion())
            {
                const unsigned insertionBases = std::min(basesLeft, gap.getLength());
                alignmentPos += insertionBases;
                basesLeft -= insertionBases;
                if (!basesLeft)
                {
                    // insertions have eaten all the bases
                    break;
                }
            }
            else
            {
                //ISAAC_THREAD_CERR << " thinking of shift" << pivotPos << "-" << ret << std::endl;
                alignmentPos -= gap.getLength();
                //ISAAC_THREAD_CERR << "shift" << pivotPos << "-" << ret << std::endl;

                overlapPos = gap.getBeginPos();
            }
        }
        --gapIndex;
    }

    if (long(binStartPos.getPosition()) > alignmentPos)
    {
//        ISAAC_THREAD_CERR << " gap places read before bin start" << std::endl;
        // this combination of gaps will have the read start alignmentPos moved before the binStartPos.
        // Don't realign this way.
        return false;
    }

    if (long(binEndPos.getPosition()) < alignmentPos)
    {
//         ISAAC_THREAD_CERR << " gap places read after bin end" << std::endl;
         // this combination of gaps will have the read start position moved at or after the binEndPos.
         // Don't realign this way.
         return false;
    }

    ret = reference::ReferencePosition(pivotPos.getContigId(), alignmentPos);
//    ISAAC_THREAD_CERR << " new startPos offset=" << offset << " startPos=" << ret << " from pivot=" << pivotPos << "and original=" << index.pos_ << std::endl;

    return true;
}

static unsigned short getTotalGapsLength(
    const uint32_t *cigarIterator,
    const uint32_t * const cigarEnd,
    unsigned &gapsCount,
    unsigned &mappedLength)
{
    gapsCount = 0;
    using alignment::Cigar;
    unsigned short ret = 0;
    for (;cigarIterator != cigarEnd; ++cigarIterator)
    {
        Cigar::Component decoded = Cigar::decode(*cigarIterator);
        if (Cigar::ALIGN == decoded.second)
        {
            mappedLength += decoded.first;
        }
        else if (Cigar::SOFT_CLIP == decoded.second)
        {
        }
        else if (Cigar::INSERT == decoded.second)
        {
            ret += decoded.first;
            ++gapsCount;
        }
        else if (Cigar::DELETE == decoded.second)
        {
            ret += decoded.first;
            ++gapsCount;
        }
        else
        {
            ISAAC_ASSERT_MSG(false, "Unexpected CIGAR operation");
        }
    }

    return ret;
}

GapRealigner::GapChoice GapRealigner::getAlignmentCost(
    const io::FragmentAccessor &fragment,
    const PackedFragmentBuffer::Index &index) const
{
    unsigned gapsCount = 0;
    GapChoice ret;
    const unsigned totalGapsLength = getTotalGapsLength(index.cigarBegin_, index.cigarEnd_, gapsCount, ret.mappedLength_);
    ret.mismatches_ = fragment.editDistance_ - totalGapsLength;
//    ISAAC_THREAD_CERR << "getAlignmentCost " << fragment << std::endl;
    ret.mismatchesPercent_ = calculateMismatchesPercent(ret.mismatches_, ret.mappedLength_);
    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment.clusterId_, "Initial mm:" << ret.mismatches_ << " gaps:" << gapsCount << " gapslength:" << totalGapsLength);

    ret.cost_ = ret.mismatches_ * mismatchCost_ + gapsCount * gapOpenCost_ + gapExtendCost_ * (totalGapsLength - gapsCount);
    ret.editDistance_ = fragment.editDistance_;
    ret.startPos_ = index.pos_;

    return ret;
}

bool GapRealigner::isBetterChoice(
    const GapChoice &choice,
    const unsigned maxMismatchesPercent,
    const GapChoice &bestChoice) const
{
//    ISAAC_THREAD_CERR << "isBetterChoice " << choice << std::endl;
    const bool ret =
        choice.mappedLength_ &&
        choice.mismatchesPercent_ <= maxMismatchesPercent &&
        (choice.cost_ < bestChoice.cost_ || (choice.cost_ == bestChoice.cost_ && choice.editDistance_ < bestChoice.editDistance_));

    return ret;
}

class TraceGapsChoice
{
    const GapRealigner::GapChoiceBitmask choice_;
    const gapRealigner::GapsRange& gaps_;

public:
    TraceGapsChoice(const GapRealigner::GapChoiceBitmask &choice,
                    const gapRealigner::GapsRange& gaps) :
                        choice_(choice), gaps_(gaps)
    {}

    friend std::ostream &operator <<(std::ostream &os, const TraceGapsChoice &trace)
    {
        unsigned gapIndex = 0;
        BOOST_FOREACH(const gapRealigner::Gap &gap, std::make_pair(trace.gaps_.first, trace.gaps_.second) )
        {
            if (trace.choice_ & (GapRealigner::GapChoiceBitmask(1) << gapIndex))
            {
                os << gap << " ";
            }
            ++gapIndex;
        }
        return os;
    }
};

void GapRealigner::verifyGapsChoice(
    const GapChoiceBitmask &choice,
    const gapRealigner::GapsRange& gaps,
    const reference::ReferencePosition& binStartPos,
    const reference::ReferencePosition& binEndPos,
    const io::FragmentAccessor& fragment,
    const std::vector<reference::Contig>& reference,
    const int originalMismatchesPercent,
    const long undoneAlignmentPos,
    GapChoice& bestChoice)
{
    unsigned pivotGapIndex = 0;
    BOOST_FOREACH(const gapRealigner::Gap& pivotGap, std::make_pair(gaps.first, gaps.second))
    //                for (unsigned pivotGapIndex = 0; gaps.size() != pivotGapIndex; ++pivotGapIndex)
    {
        if (choice & (GapChoiceBitmask(1) << pivotGapIndex))
        {
            //TODO: check binBorder overrun
            reference::ReferencePosition newStarPos;

            // verify case when anchoring occurs before pivot gap
            if (pivotGap.getBeginPos() >= binStartPos)
            {
                ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment.clusterId_, "Testing pivot before " << pivotGap);
                if (findStartPos(choice, gaps, binStartPos, binEndPos,
                                 pivotGapIndex, pivotGap.getBeginPos(), undoneAlignmentPos, newStarPos))
                {
                    const GapChoice thisChoice = verifyGapsChoice(choice, gaps, newStarPos, fragment, reference);
                    if (isBetterChoice(thisChoice, originalMismatchesPercent, bestChoice))
                    {
                        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment.clusterId_, thisChoice << "better than " << bestChoice);
                        bestChoice = thisChoice;
                    }
                    else
                    {
                        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment.clusterId_, thisChoice << "no better than " << bestChoice);
                    }
                }
            }

            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment.clusterId_, "Testing pivot after " << pivotGap);

            // verify case when anchoring occurs after pivot gap
            if (findStartPos(choice, gaps, binStartPos, binEndPos,
                             pivotGapIndex + 1, pivotGap.getEndPos(false), undoneAlignmentPos, newStarPos))
            {
                const GapChoice thisChoice = verifyGapsChoice(choice, gaps, newStarPos, fragment, reference);
                if (isBetterChoice(thisChoice, originalMismatchesPercent, bestChoice))
                {
                    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment.clusterId_, thisChoice << "better than " << bestChoice);
                    bestChoice = thisChoice;
                }
                else
                {
                    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment.clusterId_, thisChoice << "no better than " << bestChoice);
                }
            }
        }
        ++pivotGapIndex;
    }
}

GapRealigner::GapChoice GapRealigner::findBetterGapsChoice(
    const gapRealigner::GapsRange& gaps,
    const reference::ReferencePosition& binStartPos,
    const reference::ReferencePosition& binEndPos,
    const std::vector<reference::Contig>& reference,
    const io::FragmentAccessor& fragment,
    const PackedFragmentBuffer::Index& index,
    unsigned &leftToEvaluate)
{
//    if ("C5DKDANXX:1:2115:702297:0" == std::string(fragment.nameBegin(), fragment.nameEnd()))
//    {
//        ISAAC_THREAD_CERR << "C5DKDANXX:1:2115:702297:0 id=" << fragment << std::endl;
//        exit(1);
//    }
//            const gapRealigner::OverlappingGapsFilter overlappingGapsFilter(gaps);
    gapRealigner::ChooseKGapsFilter<GapChoiceBitmask> gapsFilter(gaps, gapsPerFragmentMax_);

    GapChoice bestChoice = getAlignmentCost(fragment, index);
    const int originalMismatchesPercent = bestChoice.mismatchesPercent_;
    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment.clusterId_, "Initial bestChoice " << bestChoice);

    fragmentGaps_.clear();
    fragmentGaps_.addGaps(fragment.fStrandPosition_, fragment.cigarBegin(), fragment.cigarEnd());
    const gapRealigner::GapsRange fragmentGapsRange = fragmentGaps_.allGaps();

    // 0 means none of the gaps apply. It also means the original alignment should be kept.
    for (GapChoiceBitmask choice = 0; (choice = gapsFilter.next(choice));)
    {
        if (!--leftToEvaluate)
        {
            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(
                fragment.clusterId_,
                "GapRealigner::realign: Too many gaps (" << combinationsLimit_ << " checked so far). " << fragment);
            // We've spent too much time already. Just go with what we've got.
            break;
        }

        long undoneAlignmentPos = undoExistingGaps(index, index.pos_);
        // existing deletion overlaps pivot position like this:
        //A-----C
        //AGATCAG
        //   ^pp
        ISAAC_ASSERT_MSG(undoneAlignmentPos <= long(index.pos_.getPosition()), "undoPivotPos pos " << index.pos_ << " overlapped by an existing deletion " << index);
        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment.clusterId_, "Testing choice " << int(choice) << ":" << TraceGapsChoice(choice, gaps) << " undoneAlignmentPos:" << undoneAlignmentPos);


        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID( fragment.clusterId_, "undo pivot: " << index.pos_);
        verifyGapsChoice(choice, gaps, binStartPos, binEndPos, fragment, reference, originalMismatchesPercent, undoneAlignmentPos, bestChoice);

        long lastUndoneAlignmentPos = undoneAlignmentPos;
        BOOST_FOREACH(const gapRealigner::Gap& undoPivotGap, std::make_pair(fragmentGapsRange.first, fragmentGapsRange.second))
        {
            undoneAlignmentPos = undoExistingGaps(index, undoPivotGap.getEndPos(false));
            if (lastUndoneAlignmentPos != undoneAlignmentPos)
            {
                // existing deletion overlaps pivot position like this:
                //A-----C
                //AGATCAG
                //   ^pp
                ISAAC_ASSERT_MSG(undoneAlignmentPos <= long(undoPivotGap.getEndPos(false).getPosition()), "undoPivotPos pos " << index.pos_ << " overlapped by an existing deletion " << index);
                ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment.clusterId_, "Testing choice " << int(choice) << ":" << TraceGapsChoice(choice, gaps) << " undoneAlignmentPos:" << undoneAlignmentPos);

                ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID( fragment.clusterId_, "undo pivot: " << undoPivotGap.getEndPos(false));
                verifyGapsChoice(choice, gaps, binStartPos, binEndPos, fragment, reference, originalMismatchesPercent, undoneAlignmentPos, bestChoice);
                lastUndoneAlignmentPos = undoneAlignmentPos;
            }
        }
    }
    return bestChoice;
}

/**
 * \brief Perform full realignment discarding all the existing gaps
 */
bool GapRealigner::realign(
    const gapRealigner::RealignerGaps &realignerGaps,
    const reference::ReferencePosition binStartPos,
    reference::ReferencePosition binEndPos,
    PackedFragmentBuffer::Index &index,
    io::FragmentAccessor &fragment,
    PackedFragmentBuffer &dataBuffer,
    alignment::Cigar &realignedCigars)
{
//    if ("HCCCLCCXX:7:1203:5402837:0" == std::string(fragment.nameBegin(), fragment.nameEnd()))
//    {
//        ISAAC_THREAD_CERR << fragment << std::endl;
//    }

    std::size_t bufferSizeBeforeRealignment = realignedCigars.size();
    bool makesSenseToTryAgain = false;
    bool ret = false;
//    bool firstAttempt = true;
//    lastAttemptGaps_.clear();

    if (fragment.flags_.unmapped_)
    {
        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment.clusterId_, "GapRealigner::realign: Will not realign unmapped " << fragment);
        return false;
    }

    // when realignGapsVigorously_ is set, realignment can continue forever. Make sure it does not go over the arbitrary hard limit we set
    unsigned leftToEvaluate = combinationsLimit_;

    do
    {
        makesSenseToTryAgain = false;
        using alignment::Cigar;
        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment.clusterId_, "GapRealigner::realign " << index << fragment);
        const std::vector<reference::Contig> &reference = contigList_.at(barcodeMetadataList_.at(fragment.barcode_).getReferenceIndex());
        binEndPos = reference::ReferencePosition(
            binEndPos.getContigId(), std::min<unsigned long>(binEndPos.getPosition(), reference.at(binEndPos.getContigId()).getLength()));
        if (!fragment.editDistance_)
        {
            // don't bother to look at perfectly matching ones
            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment.clusterId_, "GapRealigner::realign: Will not realign " << fragment << " because it aligns perfectly already");
            break;
        }

        if (fragment.flags_.splitAlignment_)
        {
            // if a template contains a split read, even the read that is not split cannot be realigned as then
            // the mate pos updates might not get applied to all the split parts.
            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment.clusterId_, "GapRealigner::realign: Will not realign " << fragment << " because it is a split alignment");
            break;
        }

        if (fragment.flags_.mateSplit_ )
        {
            // if a template contains a split read, even the read that is not split cannot be realigned as then
            // the mate pos updates might not get applied to all the split parts.
            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment.clusterId_, "GapRealigner::realign: Will not realign " << fragment << " because it has a split mate alignment");
            break;
        }

        if (fragment.flags_.paired_ && fragment.flags_.mateUnmapped_)
        {
            // CASAVA IndelFinder bam reader skips reads containing soft-clip operations. This becomes
            // lethal with singleton-shadow pairs as then it pairs the shadow with something else (following read I guess)
            // In cases when it pairs such a orphaned shadow with an end of a chimeric read it then produces a monster pair
            // that even the almighty ClusterMerger cannot swallow.
            // Avoid realigning singletons for now. TODO: don't realign singletons only if doing so produces soft-clipping
            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment.clusterId_, "GapRealigner::realign: Will not realign " << fragment << " because it is a singleton");
            break;
        }

        if (fragment.flags_.paired_ && (binStartPos > fragment.mateFStrandPosition_ || binEndPos <= fragment.mateFStrandPosition_))
        {
            // TODO: we can't update template lengths if mate is in a different bin. If realignment requires
            // template update and mate is in a different bin, the realignment should not be done. At the
            // moment just not realigning bin-spanning pairs at all.
            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment.clusterId_, "GapRealigner::realign: Will not realign " << fragment << " because mate is in a different bin");
            break;
        }

        if (!realignDodgyFragments_ && fragment.flags_.dodgy_)
        {
            // Normal alignments don't hit gaps by large. Zero-scored templates tend to pile up around single locations.
            // Although with gcc -O3 this passes unnoticed, it causes enormous time waste trying to realign heap of
            // misplaced reads against the gaps they happen to have.
            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment.clusterId_, "GapRealigner::realign: Will not realign " << fragment << " because it is dodgy");
            break;
        }

        if (index.pos_.getPosition() < index.getBeginClippedLength())
        {
            // reference-clipped reads are difficult to compute because ReferencePosition is not allowed to go negative.
            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment.clusterId_, "GapRealigner::realign: Will not realign " << fragment << " because begins before the reference does");
            break;
        }

        // Mate realignment might have updated our fragment.fStrandPosition_. Make sure index.pos_ is up to date
        index.pos_ = fragment.fStrandPosition_;
        RealignmentBounds bounds = extractRealignmentBounds(index);
        const gapRealigner::GapsRange gaps = realignerGaps.findGaps(fragment.clusterId_, binStartPos, bounds.beginPos_, bounds.endPos_, currentAttemptGaps_);

//            if (!firstAttempt && lastAttemptGaps_.size() == currentAttemptGaps_.size() &&
//                std::equal(lastAttemptGaps_.begin(), lastAttemptGaps_.end(), currentAttemptGaps_.begin()))
//            {
//                // no new gaps found. stop trying.
//                ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment.clusterId_, "no new gaps found. stop trying ");
//                break;
//            }
//            firstAttempt = false;
//            lastAttemptGaps_ = currentAttemptGaps_;

        if (!realignGapsVigorously_ && MAX_GAPS_AT_A_TIME < gaps.size())
        {
            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment.clusterId_, "GapRealigner::realign: Too many gaps (" << gaps.size() << "). " << fragment);
            break;
        }

        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment.clusterId_, "Found Gaps " << gaps << " for bounds " << bounds);

        GapChoice bestChoice = findBetterGapsChoice(gaps, binStartPos, binEndPos, reference, fragment, index, leftToEvaluate);

        // 0 means none of the gaps apply. It also means the original alignment should be kept.
        if (bestChoice.choice_)
        {
            if (binEndPos > bestChoice.startPos_)
            {
                ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment.clusterId_, "Applying choice for bin end pos:" << binEndPos <<
                                                       " " << bestChoice << " to gaps " << gaps << index << fragment);
                PackedFragmentBuffer::Index tmp = index;
                tmp.pos_ = bestChoice.startPos_;

                const reference::ReferencePosition contigEndPos(binEndPos.getContigId(), reference.at(binEndPos.getContigId()).forward_.size());
                if (applyChoice(bestChoice.choice_, gaps, binEndPos, contigEndPos, tmp, fragment, realignedCigars))
                {
        //            ISAAC_THREAD_CERR << " before compactCigar=" << index << fragment << std::endl;
                    if (compactCigar(reference, binEndPos, tmp, fragment, realignedCigars))
                    {
                        if (clipSemialigned_)
                        {
                            // Note! this has to be called after compactCigar as otherwise the fragment.observedLength_ is incorrect
                            SemialignedEndsClipper clipper(realignedCigars);
                            clipper.clip(reference, binEndPos, tmp, fragment);
                        }

                        index = tmp;
        //                ISAAC_THREAD_CERR << " before updatePairDetails=" << index << fragment << std::endl;
//                        updatePairDetails(index, fragment, dataBuffer);
        //                ISAAC_THREAD_CERR << "Applying choice done " << int(bestChoice) << " to gaps " << gaps << index << fragment << std::endl;
                        makesSenseToTryAgain = realignGapsVigorously_;
                        ret = true;
                    }
                    else
                    {
                        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment.clusterId_, "Ignoring choice" << bestChoice << " to gaps " << gaps << index << fragment <<
                                "due to cigar compacting having to move read into the next bin.");
                    }
                }
                else
                {
                    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment.clusterId_, "Ignoring choice" << bestChoice << " to gaps " << gaps << index << fragment <<
                            "due to applyChoice having to move read into the next bin.");
                }
            }
            else
            {
                ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment.clusterId_, "Ignoring choice" << bestChoice << " to gaps " << gaps << index << fragment <<
                        "due to realigned read having to move to the next bin.");
            }
        }
        else
        {
            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment.clusterId_, "Ignoring choice" << bestChoice << " to gaps " << gaps << index << fragment <<
                    "due to no choice being available.");
        }
    } while(makesSenseToTryAgain);

    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment.clusterId_, "before compactRealignedCigarBuffer:" << index);
    if (ret)
    {
        compactRealignedCigarBuffer(bufferSizeBeforeRealignment, index, realignedCigars);
        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment.clusterId_, "after compactRealignedCigarBuffer:" << index);
        return true;
    }
    return false;
}

void GapRealigner::compactRealignedCigarBuffer(
    const std::size_t bufferSizeBeforeRealignment,
    PackedFragmentBuffer::Index &index,
    alignment::Cigar &realignedCigars)
{
    const std::size_t cigarLength = std::distance(index.cigarBegin_, index.cigarEnd_);
    const std::size_t expectedBufferSize = bufferSizeBeforeRealignment + cigarLength;
    ISAAC_ASSERT_MSG(expectedBufferSize <= realignedCigars.size(), "Unexpected buffer increase needed. realignedCigars.size():" << realignedCigars.size() << " expectedBufferSize:" << expectedBufferSize)
    if (expectedBufferSize  != realignedCigars.size())
    {
        std::copy(index.cigarBegin_, index.cigarEnd_, realignedCigars.begin() + bufferSizeBeforeRealignment);
        index.cigarBegin_ = &*(realignedCigars.begin() + bufferSizeBeforeRealignment);
        index.cigarEnd_ = index.cigarBegin_ + cigarLength;
        realignedCigars.resize(expectedBufferSize);
    }
}


} // namespace build
} // namespace isaac
