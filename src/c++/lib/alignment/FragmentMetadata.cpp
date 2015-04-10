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
 ** \file FragmentBuilder.cpp
 **
 ** \brief See FragmentMetadata.hh
 ** 
 ** \author Roman Petrovski
 **/

#include "alignment/FragmentMetadata.hh"
#include "common/Debug.hh"

namespace isaac
{
namespace alignment
{

static long skipBackDels(const Cigar &cigarBuffer, const unsigned cigarOffset, const unsigned cigarLength, unsigned &i)
{
    long ret = 0;
    for (;i < cigarLength; ++i)
    {
        const std::pair<int, Cigar::OpCode> cigar = Cigar::decode(cigarBuffer[cigarOffset + i]);
        switch (cigar.second)
        {
            case Cigar::DELETE:
            {
                ret += cigar.first;
                break;
            }

            case Cigar::BACK:
            {
                ret -= cigar.first;
                break;
            }

            case Cigar::SOFT_CLIP:
            {
                return ret;
            }

            case Cigar::ALIGN:
            {
                return ret;
            }

            default:
            {
                ISAAC_ASSERT_MSG(false, "Reached unexpected CIGAR element before reaching ALIGN element." <<
                                 Cigar::toString(cigarBuffer, cigarOffset, cigarLength) << " got " << cigar.second << " at offset " << i);
            }
        }
    }
    ISAAC_ASSERT_MSG(false, "Reached the end of CIGAR before reaching ALIGN element." << Cigar::toString(cigarBuffer, cigarOffset, cigarLength));
    return 0;
}

unsigned FragmentMetadata::scanForAnchors(
    unsigned sequenceOffset,
    const unsigned length,
    std::vector<char>::const_iterator currentReference,
    std::vector<char>::const_iterator currentSequence,
    isaac::reference::Annotation::const_iterator currentAnnotation,
    Anchor& firstAnchor, Anchor& lastAnchor) const
{
    unsigned ret = 0;
    unsigned matchesInARow = 0;
    // scan backwards so that it's easier to check for k-uniqueness
    currentReference += length;
    sequenceOffset += length;
    currentSequence += length;
    currentAnnotation += length;
    for (unsigned j = 0; length > j; ++j)
    {
        --currentReference;
        --sequenceOffset;
        --currentSequence;
        --currentAnnotation;
        if (isMatch(*currentSequence, *currentReference))
        {
            ++matchesInARow;
            if (reference::K_UNIQUE_TOO_FAR != *currentAnnotation
                && matchesInARow >= *currentAnnotation)
            {
                firstAnchor.first = sequenceOffset;
                firstAnchor.second = sequenceOffset + *currentAnnotation;
                if (lastAnchor.empty())
                {
                    lastAnchor = firstAnchor;
                }
            }
        }
        else
        {
            ret = std::max(ret, matchesInARow);
            matchesInARow = 0;
        }
    }
    return ret;
}

double FragmentMetadata::calculateLogProbability(
    unsigned length,
    std::vector<char>::const_iterator currentReference,
    std::vector<char>::const_iterator currentSequence,
    std::vector<char>::const_iterator currentQuality) const
{
    double ret = 0.0;
    while (length--)
    {
        if (isMatch(*currentSequence, *currentReference))
        {
            ret += Quality::getLogMatch(*currentQuality);
        }
        else
        {
            ret += Quality::getLogMismatch(*currentQuality);
        }
        ++currentReference;
        ++currentSequence;
        ++currentQuality;
    }
    return ret;
}

//double FragmentMetadata::calculateInsertionLogProbability(
//    unsigned length,
//    std::vector<char>::const_iterator currentQuality) const
//{
//    double ret  = std::accumulate(currentQuality, currentQuality + length, 0.0,
//                           bind(std::plus<double>(), _1, boost::bind(&Quality::getLogMatch, _2)));
//    const char qualityMax = *std::max(currentQuality, currentQuality + length);
//    // assume one highest-quality base mismatches
//    ret -= Quality::getLogMatch(qualityMax);
//    ret += Quality::getLogMismatch(qualityMax);
//    return ret;
//}

unsigned FragmentMetadata::scanMismatches(
    unsigned sequenceOffset,
    unsigned length, bool reverse,
    const unsigned lastCycle,
    const unsigned firstCycle,
    std::vector<char>::const_iterator currentReference,
    std::vector<char>::const_iterator currentSequence)
{
    unsigned ret = 0;
    while (length--)
    {
        if (!isMatch(*currentSequence, *currentReference))
        {
            addMismatchCycle(reverse ? lastCycle - sequenceOffset : firstCycle + sequenceOffset);
            ++ret;
        }
        ++currentReference;
        ++sequenceOffset;
        ++currentSequence;
    }
    return ret;
}

unsigned FragmentMetadata::updateAlignment(
    const AlignmentCfg &cfg,
    const flowcell::ReadMetadataList &readMetadataList,
    const reference::ContigList &contigList,
    const isaac::reference::ContigAnnotations &contigAnnotations,
    unsigned contigId,
    const long strandPosition,
    const Cigar &cigarBuffer,
    const unsigned cigarOffset)
{
    const Read &read = this->getRead();
    bool reverse = this->reverse;
    std::vector<char>::const_iterator sequenceBegin = read.getStrandSequence(reverse).begin();
    std::vector<char>::const_iterator qualityBegin = read.getStrandQuality(reverse).begin();

    ISAAC_ASSERT_MSG(!contigAnnotations.at(contigId).empty(), "Empty annotation unexpected for " << contigList.at(contigId));
    ISAAC_ASSERT_MSG(!contigList.at(contigId).forward_.empty(), "Reference contig was not loaded for " << *this);

    ISAAC_ASSERT_MSG(0 <= strandPosition, "position must be positive for CIGAR update " << *this << " strandPosition:" << strandPosition <<
                     " CIGAR: " << Cigar::toString(cigarBuffer.begin() + cigarOffset, cigarBuffer.end()));
    std::vector<char>::const_iterator referenceBegin = contigList.at(contigId).forward_.begin();
    isaac::reference::Annotation::const_iterator annotationBegin = contigAnnotations.at(contigId).begin();

//    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(this->getCluster().getId(), " updateFragmentCigar : " <<
//                                           Cigar::toString(cigarBuffer.begin() + cigarOffset, cigarBuffer.end()) << " for " << *this);


    const unsigned firstCycle = readMetadataList[this->readIndex].getFirstCycle();
    const unsigned lastCycle = readMetadataList[this->readIndex].getLastCycle();

    this->cigarBuffer = &cigarBuffer;
    this->cigarOffset = cigarOffset;
    this->cigarLength = cigarBuffer.size() - this->cigarOffset;
    // adjust cigarOffset and cigarLength
    ISAAC_ASSERT_MSG(cigarBuffer.size() > this->cigarOffset, "Expecting the new cigar is not empty");

    long currentPosition = strandPosition;
    unsigned currentBase = 0;
    unsigned matchCount = 0;
    Anchor firstAnchor(true);
    Anchor lastAnchor(true);
    for (unsigned i = 0; this->cigarLength > i; ++i)
    {
        const std::pair<int, Cigar::OpCode> cigar = Cigar::decode(cigarBuffer[this->cigarOffset + i]);
        const unsigned arg = cigar.first;
        const Cigar::OpCode opCode = cigar.second;
        if (opCode == Cigar::ALIGN)
        {
            // scan backwards so that it's easier to check for k-uniqueness
            this->matchesInARow = std::max(
                this->matchesInARow,
                scanForAnchors(
                    currentBase,arg, referenceBegin + currentPosition, sequenceBegin + currentBase,
                    annotationBegin + currentPosition, firstAnchor, lastAnchor));

            this->logProbability += calculateLogProbability(
                arg, referenceBegin + currentPosition, sequenceBegin + currentBase, qualityBegin + currentBase);

            this->smithWatermanScore += cfg.normalizedMismatchScore_ * scanMismatches(
                currentBase, arg, reverse, lastCycle, firstCycle,
                referenceBegin + currentPosition, sequenceBegin + currentBase);

            matchCount += std::inner_product(
                sequenceBegin + currentBase, sequenceBegin + currentBase + arg, referenceBegin + currentPosition,
                0, std::plus<unsigned>(), &isMatch);

            // the edit distance includes all mismatches and ambiguous bases (Ns)
            this->editDistance += std::inner_product(
                sequenceBegin + currentBase, sequenceBegin + currentBase + arg, referenceBegin + currentPosition,
                0, std::plus<unsigned>(), std::not_equal_to<char>());

            currentPosition += arg;
            currentBase += arg;
        }
        else if (opCode == Cigar::INSERT)
        {
//            this->logProbability += calculateInsertionLogProbability(arg, qualityBegin + currentBase);
            this->logProbability += calculateLogProbability(
                arg, referenceBegin + currentPosition, sequenceBegin + currentBase, qualityBegin + currentBase);

            currentBase += arg;
            this->editDistance += arg;
            ++this->gapCount;
            this->smithWatermanScore += cfg.normalizedGapOpenScore_ + std::min(cfg.normalizedMaxGapExtendScore_, (arg - 1) * cfg.normalizedGapExtendScore_);
        }
        else if (opCode == Cigar::DELETE)
        {
            currentPosition += arg;
            this->editDistance += arg;
            ++this->gapCount;
            this->smithWatermanScore += cfg.normalizedGapOpenScore_ + std::min(cfg.normalizedMaxGapExtendScore_, (arg - 1) * cfg.normalizedGapExtendScore_);
            this->splitAlignment |= (arg > cfg.splitGapLength_);
        }
        else if (opCode == Cigar::BACK)
        {
            currentPosition -= arg;
            this->editDistance += arg;
            ++this->gapCount;
            this->smithWatermanScore += cfg.normalizedGapOpenScore_ + std::min(cfg.normalizedMaxGapExtendScore_, (arg - 1) * cfg.normalizedGapExtendScore_);
            this->splitAlignment = true;
        }
        else if (opCode == Cigar::FLIP)
        {
            reverse = !reverse;
            sequenceBegin = read.getStrandSequence(reverse).begin();
            qualityBegin = read.getStrandQuality(reverse).begin();
            currentBase = read.getLength() - currentBase - arg;
            this->splitAlignment = true;
            // Notice, this will count flips followed by CONTIG or position adjustment as multiple gaps which is probably not
            // ideal, but so far the gapCount is not being used for anything that requires precise value.
            ++this->gapCount;
        }
        else if (opCode == Cigar::CONTIG)
        {
            ++i;
            currentPosition += skipBackDels(cigarBuffer, this->cigarOffset, this->cigarLength, i);
            ISAAC_ASSERT_MSG(0 <= currentPosition, "Unexpected negative position adjustment: " << currentPosition <<
                             " " << Cigar::toString(cigarBuffer.begin() + cigarOffset, cigarBuffer.end()))
            ISAAC_ASSERT_MSG(contigList.at(arg).forward_.size() > std::size_t(currentPosition), "Position adjustment outside the contig bounds: " << currentPosition <<
                             " " << Cigar::toString(cigarBuffer.begin() + cigarOffset, cigarBuffer.end()))
            referenceBegin = contigList.at(arg).forward_.begin();
            annotationBegin = contigAnnotations.at(arg).begin();
            contigId = arg;
            ++this->gapCount;
            --i;
            this->splitAlignment = true;
        }
        else if (opCode == Cigar::SOFT_CLIP)
        {
            // With inversions, soft clipping can occur in the middle of CIGAR
//            ISAAC_ASSERT_MSG(0 == i || i + 1 == this->cigarLength, "Soft clippings are expected to be "
//                "found only at the ends of cigar string");
            this->logProbability =
                std::accumulate(qualityBegin + currentBase, qualityBegin + currentBase + arg,
                                this->logProbability,
                                boost::bind(std::plus<double>(), _1, boost::bind(Quality::getLogMatch, _2)));

            // NOTE! Not advancing the reference for soft clips
            currentBase += arg;
        }
        else
        {
            using boost::format;
            const format message = format("Unexpected Cigar OpCode: %d") % opCode;
            BOOST_THROW_EXCEPTION(common::PostConditionException(message.str()));
        }
    }
    this->observedLength = currentPosition - strandPosition;
    this->position = strandPosition;

    // The new anchors are best possible k-unique. Update if they are available.
    if (!firstAnchor.empty())
    {
        firstAnchor_ = firstAnchor;
        lastAnchor_ = lastAnchor;
        dodgy = false;
    }
    else if (!gapCount)
    {
        // Keep anchors and make sure they are valid only for ungapped alignments.

        // Make sure empty anchors are not placed inside the soft-clipped ends (including quality trimming).
        // This will mess up split read alignment
        const unsigned endClipOffset = getReadLength() - getEndClippedLength();
        if (firstAnchor_.empty() || firstAnchor_.first < getBeginClippedLength() || firstAnchor_.second > endClipOffset)
        {
            firstAnchor_ = Anchor(getBeginClippedLength(), getBeginClippedLength(), false);
        }
        if (lastAnchor_.empty() || lastAnchor_.first < getBeginClippedLength() || lastAnchor_.second > endClipOffset)
        {
            lastAnchor_ = Anchor(endClipOffset, endClipOffset, false);
        }
        if (lastAnchor_.empty() && !firstAnchor_.empty())
        {
            lastAnchor_ = firstAnchor_;
        }
        if (firstAnchor_.empty() && !lastAnchor_.empty())
        {
            firstAnchor_ = lastAnchor_;
        }
        // if Seed or empty anchors must be outside of clipped areas
        ISAAC_ASSERT_MSG(firstAnchor_.first >= getBeginClippedLength(), "left anchor is clipped " << *this);
        ISAAC_ASSERT_MSG(lastAnchor_.second <= endClipOffset, "right anchor is clipped " << *this);
    }

//  With FLIP command, the currentBase ends in the middle of the read.
//    ISAAC_ASSERT_MSG(currentBase == read.getLength(),
//                     "Unexpected discrepancy between cigar and sequence" << fragmentMetadata);

    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(this->getCluster().getId(), " updateFragmentCigar : " << *this);

    return matchCount;
}

} // namespace alignment
} // namespace isaac
