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
 ** Implements a translation interface required for serializing FragmentAccessor into bam
 ** 
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_BUILD_FRAGMENT_ACCESSOR_BAM_ADAPTER_HH
#define iSAAC_BUILD_FRAGMENT_ACCESSOR_BAM_ADAPTER_HH

#include <bitset>
#include <iterator>

#include "bam/Bam.hh"
#include "build/BuildContigMap.hh"
#include "build/PackedFragmentBuffer.hh"
#include "build/SaTagMaker.hh"
#include "io/Fragment.hh"
#include "reference/Contig.hh"

namespace isaac
{
namespace build
{
struct IncludeTags
{
    IncludeTags(
        bool includeAS,
        bool includeBC,
        bool includeNM,
        bool includeOC,
        bool includeRG,
        bool includeSM,
        bool includeZX,
        bool includeZY) :
            includeAS_(includeAS),
            includeBC_(includeBC),
            includeNM_(includeNM),
            includeOC_(includeOC),
            includeRG_(includeRG),
            includeSM_(includeSM),
            includeZX_(includeZX),
            includeZY_(includeZY){}
    bool includeAS_;
    bool includeBC_;
    bool includeNM_;
    bool includeOC_;
    bool includeRG_;
    bool includeSM_;
    bool includeZX_;
    bool includeZY_;
} ;


class FragmentAccessorBamAdapter
{
    const unsigned maxReadLength_;
    const flowcell::TileMetadataList &tileMetadataList_;
    const flowcell::BarcodeMetadataList &barcodeMetadataList_;
    const BuildContigMap &contigMap_;
    const isaac::reference::ContigLists &contigLists_;
    reference::ReferencePosition pos_;
    // if read got realigned pFragment->cigarBegin() points at the original cigar
    const io::FragmentAccessor* pFragment_;
    common::StaticVector<char, 100> readGroupNameBuffer_;
    common::StaticVector<char, 1024> barcodeNameBuffer_;
    common::StaticVector<char, 1024> readNameBuffer_;
    common::StaticVector<char, 10240> originalCigarBuffer_;
    const unsigned *cigarBegin_;
    const unsigned *cigarEnd_;
    bool reverse_;
    typedef common::StaticVector<char, 10240> SeqBufferType;
    SeqBufferType seqBuffer_;
    typedef common::StaticVector<char, 10240> QualBufferType;
    QualBufferType qualBuffer_;
    unsigned char forcedDodgyAlignmentScore_;
    flowcell::FlowcellLayoutList const &flowCellLayoutList_;
    const IncludeTags includeTags_;
    const bool pessimisticMapQ_;
    const unsigned splitGapLength_;

    common::StaticVector<char, 10240> saTagBuffer_;
    // If this is part of split alignment, conventional NM represents edit distance for the entire alignment. saNM_ is
    // the edit distance of this part
    unsigned saNM_;

public:
    FragmentAccessorBamAdapter(
        unsigned maxReadLength,
        const flowcell::TileMetadataList &tileMetadataList,
        const flowcell::BarcodeMetadataList &barcodeMetadataList,
        const BuildContigMap &contigMap,
        const isaac::reference::ContigLists &contigLists,
        const unsigned char forcedDodgyAlignmentScore,
        const flowcell::FlowcellLayoutList &flowCellLayoutList,
        const IncludeTags includeTags,
        const bool pessimisticMapQ,
        const unsigned splitGapLength
        ) :
        maxReadLength_(maxReadLength), tileMetadataList_(tileMetadataList),
        barcodeMetadataList_(barcodeMetadataList),
        contigMap_(contigMap), contigLists_(contigLists),
        pos_(0ul), pFragment_(0), cigarBegin_(0), cigarEnd_(0), reverse_(false),
        forcedDodgyAlignmentScore_(forcedDodgyAlignmentScore),
        flowCellLayoutList_(flowCellLayoutList),
        includeTags_(includeTags),
        pessimisticMapQ_(pessimisticMapQ),
        splitGapLength_(splitGapLength),
        saNM_(-1U)
    {
    }

    /// For serialization of aligned fragments and shadows
    FragmentAccessorBamAdapter &operator()(
        const PackedFragmentBuffer::Index &index,
        const io::FragmentAccessor& fragment)
    {
        pos_ = index.pos_;
        pFragment_ = &fragment;
        cigarBegin_ = index.cigarBegin_;
        cigarEnd_ = index.cigarEnd_;
        reverse_ = index.reverse_;

        saTagBuffer_.clear();
        saNM_ = -1U;
        ISAAC_ASSERT_MSG(
            !pFragment_->isAligned() || !pFragment_->flags_.splitAlignment_ ||
            makeSaTagString(
                *pFragment_,
                pos_,
                reverse_,
                cigarBegin_ ,cigarEnd_,
                mapq(),
                contigLists_.at(barcodeMetadataList_.at(pFragment_->barcode_).getReferenceIndex()),
                saTagBuffer_, splitGapLength_, saNM_), "Split alignment did not produce SA tag:" << *pFragment_);

        return *this;
    }

    /// For serialization of unaligned fragments
    FragmentAccessorBamAdapter &operator()(
        const io::FragmentAccessor& fragment)
    {
        pos_ = reference::ReferencePosition(reference::ReferencePosition::NoMatch);
        pFragment_ = &fragment;
        cigarBegin_ = 0;
        cigarEnd_ = 0;
        reverse_ = fragment.isReverse();
        saTagBuffer_.clear();
        saNM_ = -1U;
        return *this;
    }

    const char *readName() {
        if (pFragment_->hasName())
        {
            return pFragment_->nameBegin();
        }
        const flowcell::TileMetadata &tileMetadata = tileMetadataList_[pFragment_->tile_];

        readNameBuffer_.clear();
        const std::string &flowcellId = tileMetadata.getFlowcellId();
        std::copy(flowcellId.begin(), flowcellId.end(), std::back_inserter(readNameBuffer_));
        readNameBuffer_.push_back(':');
        const std::string &laneString = tileMetadata.getLaneString();
        std::copy(laneString.begin(), laneString.end(), std::back_inserter(readNameBuffer_));
        readNameBuffer_.push_back(':');
        const std::string &tileString = tileMetadata.getTileString();
        std::copy(tileString.begin(), tileString.end(), std::back_inserter(readNameBuffer_));
        readNameBuffer_.push_back(':');
        common::appendUnsignedNumber(readNameBuffer_, pFragment_->clusterId_);
        readNameBuffer_.push_back(':');
        readNameBuffer_.push_back('0');
        readNameBuffer_.push_back('\0');
        return &readNameBuffer_.front();
    }

    typedef std::pair<const unsigned *, const unsigned *> CigarBeginEnd;

    /**
     * \return true if the cigar changed from the original. Currently there are two situations where this can happen:
     *      a) when read gets realigned
     *      b) when split read gets split. In this case the original cigar has non-compliant with bam spec components
     */
    bool cigarChanged() const {
        // cigarBegin_ == 0 for unaligned reads. This prevents from bogus empty OC:Z being included with unaligned data
        return cigarBegin_ && cigarBegin_ != pFragment_->cigarBegin();
    }

    bam::zTag getFragmentOC() {
        static const char OC[2] = {'O', 'C'};

        if (!includeTags_.includeOC_ || !cigarChanged())
        {
            return bam::zTag();
        }
        originalCigarBuffer_.clear();
        alignment::Cigar::toString(pFragment_->cigarBegin(), pFragment_->cigarEnd(), originalCigarBuffer_);
        // push terminating zero
        originalCigarBuffer_.push_back(0);
        return bam::zTag(OC, &originalCigarBuffer_.front(), &originalCigarBuffer_.back() + 1);
    }

    bam::zTag getFragmentSA() {
        static const char SA[2] = {'S', 'A'};

        if (saTagBuffer_.empty())
        {
            return bam::zTag();
        }

        // push terminating zero
        saTagBuffer_.push_back(0);
        return bam::zTag(SA, &saTagBuffer_.front(), &saTagBuffer_.back() + 1);
    }

    CigarBeginEnd cigar() const {
        return std::make_pair(cigarBegin_, cigarEnd_);
    }

    int seqLen() const {
        return pFragment_->readLength_;
    }

    static unsigned char bamBaseFromBclByte(unsigned char bclByte){
        return !oligo::isBclN(bclByte) ? 1 << (bclByte & 0x03) : 15;
    }

    static unsigned char bamBasesFromBclShort(unsigned short bclShort){
        unsigned char *pBclShort(reinterpret_cast<unsigned char*>(&bclShort));
        return (bamBaseFromBclByte(pBclShort[0]) << 4) | bamBaseFromBclByte(pBclShort[1]);
    }

    static unsigned char reverseBamBasesFromBclShort(unsigned short bclShort){
        unsigned char *pBclShort(reinterpret_cast<unsigned char*>(&bclShort));
        return bamBaseFromBclByte(oligo::getReverseBcl(pBclShort[0])) | (bamBaseFromBclByte(oligo::getReverseBcl(pBclShort[1])) << 4);
    }

    static unsigned char bamQualFromBclByte(unsigned char bclByte){
        return bclByte >> 2;
    }

    typedef std::pair<SeqBufferType::const_iterator, SeqBufferType::const_iterator> SeqBeginEnd;
    SeqBeginEnd seq() {
        seqBuffer_.resize((seqLen()+1)/2, 15);
        //todo: verify if endianness is not an issue here

        if (pFragment_->isReverse() == reverse_)
        {
            const unsigned short *twoBasesBegin(reinterpret_cast<const unsigned short*>(pFragment_->basesBegin()));
            std::transform(twoBasesBegin, twoBasesBegin + seqLen() / 2, seqBuffer_.begin(), bamBasesFromBclShort);
            if (seqLen() % 2)
            {
                seqBuffer_.back() = bamBaseFromBclByte(*(pFragment_->basesEnd() - 1)) << 4;
            }
        }
        else
        {
            const std::reverse_iterator<const unsigned short *> twoBasesBegin(reinterpret_cast<const unsigned short*>(pFragment_->basesEnd()));
            std::transform(twoBasesBegin, twoBasesBegin + seqLen() / 2, seqBuffer_.begin(), reverseBamBasesFromBclShort);
            if (seqLen() % 2)
            {
                seqBuffer_.back() = reverseBamBasesFromBclShort(*(pFragment_->basesBegin())) << 4;
            }
        }

        return SeqBeginEnd(seqBuffer_.begin(), seqBuffer_.end());
    }

    typedef std::pair<QualBufferType::const_iterator, QualBufferType::const_iterator> QualBeginEnd;
    QualBeginEnd qual() {
        qualBuffer_.clear();
        if (pFragment_->isReverse() == reverse_)
        {
            std::transform(pFragment_->basesBegin(), pFragment_->basesEnd(), std::back_inserter(qualBuffer_), bamQualFromBclByte);
        }
        else
        {
            std::transform(
                std::reverse_iterator<const unsigned char *>(pFragment_->basesEnd()),
                std::reverse_iterator<const unsigned char *>(pFragment_->basesBegin()),
                std::back_inserter(qualBuffer_), bamQualFromBclByte);
        }
        return SeqBeginEnd(qualBuffer_.begin(), qualBuffer_.end());
    }

    int refId() const {
        return pos_.isNoMatch() ?
            -1 :
            contigMap_.getMappedContigIndex(barcodeMetadataList_.at(pFragment_->barcode_).getReferenceIndex(),
                                            pos_.getContigId());
    }

    int pos() const {
        return pos_.isNoMatch() ? -1 : pos_.getPosition();
    }

    unsigned char mapq() const {
        if (pFragment_->flags_.properPair_)
        {
            return  pFragment_->flags_.dodgy_ ?
                forcedDodgyAlignmentScore_:
                std::min<unsigned>(60U,
                                   pessimisticMapQ_ ?
                                       std::min(pFragment_->alignmentScore_, pFragment_->templateAlignmentScore_) :
                                       std::max(pFragment_->alignmentScore_, pFragment_->templateAlignmentScore_));
        }
        return pFragment_->flags_.dodgy_ ?
            (forcedDodgyAlignmentScore_) :
            std::min<unsigned>(60U, pFragment_->alignmentScore_);
    }

    bam::iTag getFragmentSM() const {
        static const char SM[2] = {'S', 'M'};
        return !includeTags_.includeSM_ ?
            bam::iTag() : bam::iTag(SM, pFragment_->alignmentScore_);
    }

    bam::iTag getFragmentAS() const {
        static const char AS[2] = {'A', 'S'};
        return (!includeTags_.includeAS_ || !pFragment_->flags_.properPair_) ?
            bam::iTag() :  bam::iTag(AS, pFragment_->templateAlignmentScore_);
    }

    /**
     * \brief Read group IDs are
     */
    bam::zTag getFragmentRG() {
        static const char RG[2] = {'R', 'G'};

        if (!includeTags_.includeRG_)
        {
            return bam::zTag();
        }

        readGroupNameBuffer_.clear();
        const flowcell::BarcodeMetadata &barcode = barcodeMetadataList_[pFragment_->barcode_];
        // barcode index is unique within the data analysis
        common::appendUnsignedInteger(readGroupNameBuffer_, barcode.getIndex());
        readGroupNameBuffer_.push_back('\0');

        return bam::zTag(RG, &readGroupNameBuffer_.front());
    }

    bam::iTag getFragmentNM() const {
        static const char NM[2] = {'N', 'M'};
        return includeTags_.includeNM_ ?
            bam::iTag(NM, pFragment_->flags_.splitAlignment_ ? saNM_ : pFragment_->editDistance_) : bam::iTag();
    }

    bam::zTag getFragmentBC() {
        static const char BC[2] = {'B', 'C'};

        if (!includeTags_.includeBC_)
        {
            return bam::zTag();
        }

        barcodeNameBuffer_.clear();
        const std::string &barcodeName =
            barcodeMetadataList_[pFragment_->barcode_].getName();

        if(!flowCellLayoutList_[tileMetadataList_[pFragment_->tile_].getFlowcellIndex()].getBarcodeCycles().size())
        {
            // Use barcode from the sample sheet
            std::copy(barcodeName.begin(), barcodeName.end(), std::back_inserter(barcodeNameBuffer_));
        }
        else
        {
            // Use barcode from fragment
            oligo::unpackKmer(
                pFragment_->barcodeSequence_,
                flowCellLayoutList_[tileMetadataList_[pFragment_->tile_].getFlowcellIndex()].getBarcodeCycles().size(),
                back_inserter(barcodeNameBuffer_));
        }

        // Null terminate barcode
        barcodeNameBuffer_.push_back('\0');
        return bam::zTag(BC, &barcodeNameBuffer_.front());
    }

    bam::iTag getFragmentZX() const {
        static const char ZX[2] = {'Z', 'X'};
        return (includeTags_.includeZX_ && pFragment_->isClusterXySet()) ? bam::iTag(ZX, pFragment_->clusterX_) :  bam::iTag();
    }

    bam::iTag getFragmentZY() const {
        static const char ZY[2] = {'Z', 'Y'};
        return (includeTags_.includeZY_ && pFragment_->isClusterXySet()) ? bam::iTag(ZY, pFragment_->clusterY_) :  bam::iTag();
    }

    unsigned flag() const {
        std::bitset<12> bs;
        bs.set(0, pFragment_->flags_.paired_);
        bs.set(1, pFragment_->flags_.properPair_);
        bs.set(2, pFragment_->flags_.unmapped_);
        bs.set(3, pFragment_->flags_.paired_ && pFragment_->flags_.mateUnmapped_);
        // for parts of a split read, the reverse status might not match the pFragment_->flags_.reverse_
        bs.set(4, reverse_);
        bs.set(5, pFragment_->flags_.mateReverse_);
        bs.set(6, pFragment_->flags_.paired_ && !pFragment_->flags_.secondRead_);
        bs.set(7, pFragment_->flags_.paired_ && pFragment_->flags_.secondRead_);

        bs.set(9, pFragment_->flags_.failFilter_);
        bs.set(10, pFragment_->flags_.duplicate_);
        bs.set(11, pFragment_->flags_.splitAlignment_ && isSupplementaryAlignment());
        return bs.to_ulong();
    }

    bool unmapped() const
    {
        return pFragment_->flags_.unmapped_;
    }

    int nextRefId() const {
        return pFragment_->flags_.paired_ ?
            (pFragment_->flags_.unmapped_ && pFragment_->flags_.mateUnmapped_ ?
                -1 : contigMap_.getMappedContigIndex(barcodeMetadataList_.at(pFragment_->barcode_).getReferenceIndex(),
                                                     pFragment_->mateFStrandPosition_.getContigId())) : -1;
    }

    int nextPos() const {
        return pFragment_->flags_.paired_ ?
            (pFragment_->flags_.unmapped_ && pFragment_->flags_.mateUnmapped_ ?
                            -1 : pFragment_->mateFStrandPosition_.getPosition()) : -1;
    }

    int tlen() const { return pFragment_->bamTlen_; }

    int observedLength() const { return alignment::computeObservedLength(cigarBegin_, cigarEnd_); }
private:

    /**
     * \brief Alignment designated by the beginning of the compound CIGAR is primary. All others are supplementary
     */
    bool isSupplementaryAlignment() const
    {
        if (reverse_ != pFragment_->isReverse() || pos_ != pFragment_->fStrandPosition_ ||
            pFragment_->cigarLength_ <= std::distance(cigarBegin_ ,cigarEnd_))
        {
            return true;
        }

        // compare up to one less than current CIGAR. Assumption is that current CIGAR contains a soft-clip for the
        // rest of the sequence at the end.
        if (cigarEnd_ - 1 == std::mismatch(cigarBegin_, cigarEnd_ - 1, pFragment_->cigarBegin()).first)
        {
            // our cigar matches beginning of the original unsplit cigar, we're the primary alignment
            return false;
        }

        return true;
    }
};

} // namespace build
} // namespace isaac

#endif // #ifndef iSAAC_BUILD_FRAGMENT_ACCESSOR_BAM_ADAPTER_HH
