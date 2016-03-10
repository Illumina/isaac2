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
 ** \file FragmentMetadata.hh
 **
 ** \brief Component to encapsulate all metadata associated to a fragment.
 ** 
 ** \author Come Raczy
 **/

#ifndef iSAAC_ALIGNMENT_FRAGMENT_METADATA_HH
#define iSAAC_ALIGNMENT_FRAGMENT_METADATA_HH

#include <vector>
#include <algorithm>
#include <iostream>

#include "../common/StaticVector.hh"
#include "alignment/AlignmentCfg.hh"
#include "alignment/Mismatch.hh"
#include "alignment/Cigar.hh"
#include "alignment/Cluster.hh"
#include "alignment/SeedId.hh"
#include "reference/KUniqueness.hh"
#include "reference/ReferencePosition.hh"

namespace isaac
{
namespace alignment
{


struct Anchor : public std::pair<unsigned short, unsigned short>
{
    Anchor(unsigned short f, unsigned short s, bool kUnique) : std::pair<unsigned short, unsigned short>(f,s), kUnique_(kUnique){}
    Anchor(bool kUnique = false) : std::pair<unsigned short, unsigned short>(0,0), kUnique_(kUnique){}
    unsigned short length() const {return second - first;}
    bool empty() const {return second == first;}
    bool kUnique_;

    friend std::ostream &operator <<(std::ostream &os, const Anchor &a)
    {
        return os << "Anchor(" << a.first << "," << a.second << "," << a.kUnique_ << ")";
    }
};

/**
 ** \brief Alignment information for a fragment (as defined by the SAM Format
 ** Specification (v1.4-r962)
 **
 ** This component is the building block for the FragmentBuilder. It is
 ** designed for efficiency, and does not involve any memory allocation (which
 ** is the reason why it does not and should not store the CIGAR).
 **/
struct FragmentMetadata
{
    static const unsigned MAX_CYCLES = 1024;
    /**
     ** \brief Cluster associated to the fragment
     **
     ** Note, it is the responsibility of the calling code to ensure the
     ** validity of the cluster.
     **/
    const Cluster *cluster;

    /// Id of the contig where the fragment is located
    unsigned contigId;
    /**
     ** \brief 0-based leftmost position of the fragment on the forward strand
     ** of the contig.
     **
     ** Even though the position can become negative while building the fragment
     ** (before calculating the cigar), the final position will be guaranteed to
     ** be positive (or 0). The final position is the position of the first
     ** aligned base ('M', '=' or 'X'). If the read extends outside the contig,
     ** this will be reflected by appropriate insertions or clipping operations
     ** in the cigar.
     **/
    long position;

    /// number of bases clipped from the lowest read cycle irrespective of alignment
    unsigned short lowClipped;
    /// number of bases clipped from the highest read cycle irrespective of alignment
    unsigned short highClipped;

    /**
     ** \brief observed length of the fragment on the contig
     **
     ** If there are no indels and no clipping, this would be the length of the
     ** read. With indels this is the read length minus the insertions plus the
     ** deletions (resp. to and from the reference). Clipped areas of the read
     ** are also subtracted from the length of the read.
     **/
    unsigned observedLength;
    /// 0-based index of the read in the list of ReadMetadata
    unsigned readIndex;
    /// Orientation of the read. False is forward, true is reverse
    bool reverse:1;
    bool splitAlignment:1;
    bool dodgy:1;
    bool highRepeat:1;
    /// Cigar offset in the associated cigar buffer (see FragmentBuilder)
    unsigned cigarOffset;
    /// Number of operations in the cigar
    unsigned cigarLength;
    /// Buffer containing the cigar data
    const Cigar *cigarBuffer;

    /// Number of mismatches in the alignment (can't be more than read length)
    unsigned mismatchCount;

    /// Longest stretch of matches
    unsigned matchesInARow;

    /// Number of breakpoints (indels and other splits) in the fragment
    unsigned gapCount;

    /// Edit distance from the alignment (including indels and ambiguous bases)
    unsigned editDistance;

    // set fo cycles containing mismatches (outside indels).
    typedef std::bitset<MAX_CYCLES> MismatchCycles;
    std::bitset<MAX_CYCLES> mismatchCycles;

    class ConstMismatchCycleIterator :
        public boost::iterator_facade<
        ConstMismatchCycleIterator
        , unsigned
        , boost::forward_traversal_tag
        , unsigned const &
        >
    {
        unsigned cycle_;
        boost::reference_wrapper<const MismatchCycles> mismatchCycles_;
     public:
        explicit ConstMismatchCycleIterator(const MismatchCycles &mismatchCycles, const std::size_t cycle = 1)
          : cycle_(cycle), mismatchCycles_(mismatchCycles)
        {
            ISAAC_ASSERT_MSG(cycle_, "Cycle must be a 1-based integer.");
            while (cycle_ != (mismatchCycles_.get().size() + 1) && !mismatchCycles_.get().test(cycle_ - 1)) ++cycle_;
        }

     private:
        friend class boost::iterator_core_access;

        void increment() { ++cycle_; while (cycle_ != (mismatchCycles_.get().size() + 1) && !mismatchCycles_.get().test(cycle_ - 1)) ++cycle_; }

        bool equal(ConstMismatchCycleIterator const& other) const
        {
            ISAAC_ASSERT_MSG(mismatchCycles_.get_pointer() == other.mismatchCycles_.get_pointer(), "Illegal compare for iterators of different containers.");
            return this->cycle_ == other.cycle_;
        }

        const unsigned &dereference() const {
            return cycle_; }
    };

    /**
     ** \brief natural logarithm of the probability that the fragment is correct
     **
     ** This is intrinsic to the fragment and depends only of the quality of all
     ** matching base and and the quality of all mismatching bases (indel are
     ** counted as matching bases). See AlignmentQuality for the detail of the
     ** probabilities used in each case.
     **/
    double logProbability;

    /// The id of the seed that produced the alignment candidate. Valid only prior to consolidation.
    int firstSeedIndex;

    /// lowest subsequence in the direction of the reference that matches a k-unique region in the reference
    Anchor firstAnchor_;
    /// highest subsequence in the direction of the reference that matches a k-unique region in the reference
    Anchor lastAnchor_;

    /**
     ** \brief Alignment score in the global context of the reference
     **
     ** This depends on the all the pLogCorrect values for all the possible
     ** alignments for this read across the whole reference. It also takes into
     ** account the rest-of-genome correction. Value of -1U indicates unknown alignment score.
     **/
    unsigned alignmentScore;

    // Weighted sum of mismatch and gap penalties similar to what's used for Smith-Waterman alignment
    unsigned smithWatermanScore;

    /**
     ** \brief Comparison of FragmentMetadata by reference position
     **/
    bool operator<(const FragmentMetadata &f) const
    {
        ISAAC_ASSERT_MSG(cluster == f.cluster && readIndex == f.readIndex,
                         "Comparison makes sense only for metadata representing the same fragment " <<
                         *this << " vs " << f);
        return
            (this->contigId < f.contigId ||
                (this->contigId == f.contigId && (this->position < f.position ||
                    (this->position == f.position && (this->reverse < f.reverse ||
                        (reverse == f.reverse && (observedLength < f.observedLength ||
                            (observedLength == f.observedLength && logProbability < f.logProbability))))))));
    }

    bool operator == (const FragmentMetadata &that) const
    {
        ISAAC_ASSERT_MSG(cluster == that.cluster && readIndex == that.readIndex,
                         "Comparison makes sense only for metadata representing the same fragment " <<
                         *this << " vs " << that);
        return position == that.position && contigId == that.contigId && reverse == that.reverse &&
            observedLength ==that.observedLength &&
            ISAAC_LP_EQUALS(logProbability, that.logProbability);
    }

    bool operator != (const FragmentMetadata &that) const
    {
        return !(*this == that);
    }

    friend std::ostream &operator<<(std::ostream &os, const FragmentMetadata &f);

    /**
     * \param readLength is needed to preallocate buffer to avoid memory operations during the processing.
     *        auxiliary code, such as unit tests, does not need to supply it.
     */
    FragmentMetadata():
        cluster(0), contigId(reference::ReferencePosition::MAX_CONTIG_ID), position(0), lowClipped(0), highClipped(0),
        observedLength(0), readIndex(0), reverse(false), splitAlignment(false), dodgy(true), highRepeat(false), cigarOffset(0),
        cigarLength(0), cigarBuffer(0), mismatchCount(0), matchesInARow(0), gapCount(0), editDistance(0), logProbability(0.0),
        firstSeedIndex(-1),
        firstAnchor_(0, 0, false),
        lastAnchor_(0, 0, false),
        alignmentScore(-1U),
        smithWatermanScore(0)
    {
    }

    FragmentMetadata(const Cluster *cluster, unsigned readIndex,
                     const unsigned short lowClipped = 0, const unsigned short highClipped = 0,
                     const bool reverse = false,
                     const unsigned contigId = reference::ReferencePosition::MAX_CONTIG_ID, const long position = 0):
        cluster(cluster), contigId(contigId), position(position), lowClipped(lowClipped), highClipped(highClipped),
        observedLength(0), readIndex(readIndex), reverse(reverse), splitAlignment(false), dodgy(true), highRepeat(false), cigarOffset(0),
        cigarLength(0), cigarBuffer(0), mismatchCount(0), matchesInARow(0), gapCount(0), editDistance(0), logProbability(0.0),
        firstSeedIndex(-1),
        firstAnchor_(0, 0, false),
        lastAnchor_(getReadLength(), getReadLength(), false),
        alignmentScore(-1U),
        smithWatermanScore(0)
    {
    }

    bool isReverse() const {return reverse;}
    unsigned getReadLength() const
    {
        assert(0 != cluster);
        assert(cluster->size() > readIndex);
        return (*cluster)[readIndex].getLength();
    }
    unsigned getReadIndex() const {return readIndex;}
    unsigned getObservedLength() const {return isAligned() ? observedLength : 0;}
    unsigned getAlignmentScore() const {return alignmentScore;}
    void setAlignmentScore(unsigned as) {alignmentScore = as;}
    unsigned getCigarLength() const {return cigarLength;}
    /// Position of the first base of the fragment
    reference::ReferencePosition getFStrandReferencePosition() const
    {
        return !isNoMatch() ?
            reference::ReferencePosition(contigId, position) :
            reference::ReferencePosition(reference::ReferencePosition::NoMatch);
    }
    /// Position of the last base of the fragment
    reference::ReferencePosition getRStrandReferencePosition() const
    {
        // ensure that the position is positive!
        return !isNoMatch() ?
            // observedLength can be zero if the CIGAR is soft-clipped to death or in case of
            // split alignment with alignment position immediately following the alignment position of the last base
            // think local translocations (most of them are fake though)
            reference::ReferencePosition(contigId, position + std::max(observedLength, 1U) - 1) :
            reference::ReferencePosition(reference::ReferencePosition::NoMatch);
    }
    /// Position of the fragment
    reference::ReferencePosition getStrandReferencePosition() const {
        return isReverse() ? getRStrandReferencePosition() : getFStrandReferencePosition();
    }

    /// Same as f-strand position
    reference::ReferencePosition getBeginReferencePosition() const
    {
        return getFStrandReferencePosition();
    }

    /// Different from r-strand position in that it always points to the base following the last unclipped base of the fragment
    reference::ReferencePosition getEndReferencePosition() const
    {
        return !isNoMatch() ?
            reference::ReferencePosition(contigId, position + observedLength) :
            reference::ReferencePosition(reference::ReferencePosition::NoMatch);
    }

/// First cycle of fragment bcl data
    std::vector<char>::const_iterator getBclData() const {
        return cluster->getBclData(getReadIndex());
    }
    /// Cluster of the fragment
    const Cluster &getCluster() const {
        return *cluster;
    }

    const Read &getRead() const {
        return getCluster()[getReadIndex()];
    }

    Cigar::const_iterator cigarBegin() const
    {
        ISAAC_ASSERT_MSG(isAligned(), "Requesting CIGAR of unaligned fragment is not allowed");
        return cigarBuffer->begin() + cigarOffset;
    }

    Cigar::const_iterator cigarEnd() const
    {
        ISAAC_ASSERT_MSG(isAligned(), "Requesting CIGAR of unaligned fragment is not allowed");
        return cigarBuffer->begin() + cigarOffset + cigarLength;
    }

    unsigned getBeginClippedLength() const
    {
        if (cigarBuffer && cigarLength)
        {
            Cigar::Component operation = Cigar::decode(cigarBuffer->at(cigarOffset));
            if (Cigar::SOFT_CLIP == operation.second)
            {
                return operation.first;
            }
        }
        return 0;
    }

    long getEndClippedLength() const
    {
        if (cigarBuffer && cigarLength)
        {
            Cigar::Component operation = Cigar::decode(cigarBuffer->at(cigarOffset + cigarLength - 1));
            if (Cigar::SOFT_CLIP == operation.second)
            {
                return operation.first;
            }
        }
        return 0;
    }

    /// Unlike the observed length, excludes gaps (deletions and insertion bases)
    unsigned getMappedLength() const
    {
        ISAAC_ASSERT_MSG(cigarBuffer && cigarLength, "Read must have a valid CIGAR");
        return Cigar::getMappedLength(cigarBuffer->begin() + cigarOffset,
                                      cigarBuffer->begin() + cigarOffset + cigarLength);
    }

    /**
     * \brief Returns unadjusted position if it is adjusted due to a soft clipping
     */
    long getUnclippedPosition() const
    {
        return position - getBeginClippedLength();
    }

    unsigned getMismatchCount() const {
        return mismatchCount;
    }

    unsigned getGapCount() const {
        return gapCount;
    }

    unsigned getEditDistance() const {
        return editDistance;
    }


    ConstMismatchCycleIterator getMismatchCyclesBegin() const {return ConstMismatchCycleIterator(mismatchCycles, 1);}
    ConstMismatchCycleIterator getMismatchCyclesEnd() const {return ConstMismatchCycleIterator(mismatchCycles, MAX_CYCLES + 1);}

    void addMismatchCycle(const unsigned cycle)
    {
        ISAAC_ASSERT_MSG(cycle > 0, "Cycle numbers expected to be 1-based." << *this);
        ISAAC_ASSERT_MSG(MAX_CYCLES >= cycle, "Cycle number is too high. Check MAX_CYCLES." << *this);
        mismatchCycles.set(cycle - 1);
        ++mismatchCount;
    }


    std::string getCigarString() const
    {
        if (cigarBuffer && cigarLength)
        {
            return Cigar::toString(*cigarBuffer, cigarOffset, cigarLength);
        }
        else
        {
            return "";
        }
    }

    std::ostream &serializeCigar(std::ostream &os) const
    {
        if (cigarBuffer && cigarLength)
        {
            return Cigar::toStream(cigarBuffer->begin() + cigarOffset, cigarBuffer->begin() + cigarOffset + cigarLength, os);
        }
        else
        {
            return os;
        }
    }

    /**
     ** \brief The cigarLength can be used to identify if a fragment has been
     ** successfully aligned
     **/
    bool isAligned() const {return 0 != cigarLength;}
    /**
     ** \brief Marks read as unaligned.
     **/
    void setUnaligned() {cigarBuffer = 0; cigarLength = 0; alignmentScore = -1U; reverse = false, splitAlignment = false;}

    /**
     ** \brief Marks read as something that has no match position. This is different from setUnaligned
     **        as unaligned shadows still have a position of their orphan
     **/
    void setNoMatch()
        {setUnaligned(); contigId = reference::ReferencePosition::MAX_CONTIG_ID; position = 0; dodgy = true;}
    bool isNoMatch() const {return reference::ReferencePosition::MAX_CONTIG_ID == contigId;}

    /**
     * \brief Notion of uniquely aligned in CASAVA means that a fragment was
     *        seen to have only one candidate alignment position. As this is
     *        highly dependent on the choice of seeds, alignment score based
     *        approximation should do.
     */
    bool isUniquelyAligned() const {return isAligned() && hasAlignmentScore() && 3 < getAlignmentScore();}

    bool isRepeat() const {return isAligned() && hasAlignmentScore() && 3 >= getAlignmentScore();}
    bool isHighRepeat() const {return highRepeat;}

    unsigned getQuality() const
    {
        const std::vector<char> &quality = getRead().getForwardQuality();

        return std::accumulate(quality.begin(), quality.end(), 0U);
    }

    const std::vector<char> &getStrandSequence() const {return getRead().getStrandSequence(reverse);}
    const std::vector<char> &getStrandQuality() const {return getRead().getStrandQuality(reverse);}

    int getFirstSeedIndex() const {return firstSeedIndex;}

    bool hasAlignmentScore() const {return -1U != alignmentScore;}

    void incrementClipLeft(const unsigned short bases) {position += bases; if (reverse) {highClipped += bases;} else {lowClipped += bases;}}
    void incrementClipRight(const unsigned short bases) {if (reverse) {lowClipped += bases;} else {highClipped += bases;}}

    /// \return number of bases clipped on the left side of the fragment in respect to the reference
    unsigned short leftClipped() const {return reverse ? highClipped : lowClipped;}
    /// \return number of bases clipped on the right side of the fragment in respect to the reference
    unsigned short rightClipped() const {return reverse ? lowClipped : highClipped;}

    /// \return number of bases clipped on the left side of the fragment in respect to the reference
    unsigned short &leftClipped() {return reverse ? highClipped : lowClipped;}
    /// \return number of bases clipped on the right side of the fragment in respect to the reference
    unsigned short &rightClipped() {return reverse ? lowClipped : highClipped;}

    void resetAlignment()
    {
        *this = FragmentMetadata(cluster, readIndex, lowClipped, highClipped, reverse, contigId, getUnclippedPosition());
    }

    void resetClipping()
    {
        ISAAC_ASSERT_MSG(!isAligned(), "Alignment must be reset before clipping");
        lowClipped = 0;
        highClipped = 0;
    }

    bool isWellAnchored() const
    {
        return  startAnchored();
    }

    bool isKUnique() const
    {
        // normally empty anchors are not k-unique.
        // however in some unit tests it is important to be able to have an empty k-unique anchor
        return firstAnchor_.kUnique_;
    }

    unsigned getContigId() const {return contigId;}

    long getPosition() const {return position;}

    bool startAnchored() const {return !firstAnchor_.empty();}
    bool endAnchored() const {return !lastAnchor_.empty();}

    /*
     * \brief When alignment candidate is created based on seed matches, we can assume seeds anchor the alignment if they don't
     *        have the neighbor flag set. No neighbor flag means that there are no neighbors in genome and we have visited all repeats
     */
    void setSeedAnchor(
        const flowcell::ReadMetadata &readMetadata,
        const SeedMetadata &seedMetadata)
    {
        ISAAC_ASSERT_MSG(!isWellAnchored(), "Unexpected setSeedAnchors on an anchored fragment");
        firstAnchor_.kUnique_ = false;
        if (reverse)
        {
            // 'seedPosition + seedLength + seedOffset' is the first position past the end of the read
            firstAnchor_.first = readMetadata.getLength() - seedMetadata.getOffset() - seedMetadata.getLength();
        }
        else
        {
            firstAnchor_.first = seedMetadata.getOffset();
        }
        firstAnchor_.second = firstAnchor_.first + seedMetadata.getLength();
        lastAnchor_ = firstAnchor_;
    }

    void setAnchors(const Anchor &firstAnchor, const Anchor &lastAnchor)
    {
        firstAnchor_ = firstAnchor;
        lastAnchor_ = lastAnchor;
    }
    /*
     * \brief make sure best anchors are kept in the fragment. Best anchors are those that are most further apart
     */
    void mergeAnchors(const FragmentMetadata &fragment)
    {
        highRepeat &= fragment.highRepeat;
        const Anchor &firstAnchor = fragment.firstAnchor_;
        const Anchor &lastAnchor = fragment.lastAnchor_;
        // endAnchored if startAnchored
        if (!firstAnchor.empty())
        {
            if (startAnchored())
            {
                if (firstAnchor_.second > firstAnchor.second)
                {
                    firstAnchor_ = firstAnchor;
                }
            }
            else
            {
                firstAnchor_ = firstAnchor;
            }

            if (lastAnchor_.empty())
            {
                lastAnchor_ = firstAnchor_;
            }
            else if (lastAnchor_.first < lastAnchor.first)
            {
                lastAnchor_ = lastAnchor;
            }
        }
        else
        {
            // otherwise new anchors are empty, so, nothing to merge.
            ISAAC_ASSERT_MSG(lastAnchor.empty(), "Invalid combination of anchors");
        }
    }

    void resetAnchors()
    {
        firstAnchor_ = Anchor(false);
        lastAnchor_ = Anchor(false);
    }


    unsigned updateAlignment(
        const AlignmentCfg &cfg,
        const flowcell::ReadMetadataList &readMetadataList,
        const reference::ContigList &contigList,
        const isaac::reference::ContigAnnotations &contigAnnotations,
        unsigned contigId,
        const long strandPosition,
        const Cigar &cigarBuffer,
        const unsigned cigarOffset);

private:
    unsigned scanForAnchors(
        unsigned sequenceOffset,
        const unsigned length,
        std::vector<char>::const_iterator currentReference,
        std::vector<char>::const_iterator currentSequence,
        isaac::reference::Annotation::const_iterator currentAnnotation,
        Anchor& firstAnchor, Anchor& lastAnchor) const;

    double calculateLogProbability(
        unsigned length,
        std::vector<char>::const_iterator currentReference,
        std::vector<char>::const_iterator currentSequence,
        std::vector<char>::const_iterator currentQuality) const;

    double calculateInsertionLogProbability(
        unsigned length,
        std::vector<char>::const_iterator currentQuality) const;

    unsigned scanMismatches(
        unsigned sequenceOffset,
        unsigned length, bool reverse, const unsigned lastCycle,
        const unsigned firstCycle,
        std::vector<char>::const_iterator currentReference,
        std::vector<char>::const_iterator currentSequence);
};

#ifndef _GLIBCXX_DEBUG
BOOST_STATIC_ASSERT(sizeof(FragmentMetadata) <= 232);
#endif

typedef std::vector<FragmentMetadata> FragmentMetadataList;
typedef FragmentMetadataList::const_iterator FragmentIterator;


inline std::ostream &operator<<(std::ostream &os, const FragmentMetadata &f)
{
    os << "FragmentMetadata("
              << (f.cluster ? f.getCluster().getId() : 0UL) << "id "
              << f.contigId << ":"
              << f.position << ", "
              << f.observedLength << "bp r"
              << f.readIndex << " "
              << (f.reverse ? 'R' : 'F') << " "
              << f.mismatchCount << "mm "
              << f.matchesInARow << "mir "
              << f.gapCount << "g "
              << f.editDistance << "ed ";
    return f.serializeCigar(os) << " "
              << f.dodgy << "dod "
              << f.splitAlignment << "sa "
              << f.logProbability << "lp "
              << f.alignmentScore << "sm "
              << f.smithWatermanScore << "sws "
              << f.isWellAnchored() << "wa "
              << f.firstAnchor_<< "fa "
              << f.lastAnchor_<< "la "
              << f.leftClipped() << "lc "
              << f.rightClipped() << "rc ";
    if (f.cluster)
    {
        return os << common::makeFastIoString(f.getCluster().nameBegin(), f.getCluster().nameEnd())
            << ")";
    }
    return os << ")";
}


} // namespace alignment
} // namespace isaac

#endif // #ifndef iSAAC_ALIGNMENT_FRAGMENT_METADATA_HH
