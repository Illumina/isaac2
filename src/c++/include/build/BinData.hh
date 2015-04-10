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
 ** \file BinData.hh
 **
 ** Performs sorting and duplicate marking on a single alignment bin.
 ** 
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_BUILD_BIN_DATA_HH
#define iSAAC_BUILD_BIN_DATA_HH

#include <numeric>

#include <boost/filesystem.hpp>
#include <boost/foreach.hpp>
#include <boost/ptr_container/ptr_vector.hpp>

#include "alignment/BinMetadata.hh"
#include "build/BarcodeBamMapping.hh"
#include "build/FragmentAccessorBamAdapter.hh"
#include "build/FragmentIndex.hh"
#include "build/PackedFragmentBuffer.hh"
#include "build/gapRealigner/RealignerGaps.hh"
#include "io/FileBufCache.hh"

namespace isaac
{
namespace build
{

enum GapRealignerMode
{
    /// don't realign
    REALIGN_NONE,
    /// Realign against gaps found within the sample
    REALIGN_SAMPLE,
    /// Realign against gaps found in all samples of the same project
    REALIGN_PROJECT,
    /// Realign against all gaps present in the data
    REALIGN_ALL
};

struct BinData : public std::vector<PackedFragmentBuffer::Index>
{
    typedef std::vector<PackedFragmentBuffer::Index> BaseType;

public:
    BinData(
        const unsigned realignedGapsPerFragment,
        const BarcodeBamMapping &barcodeBamMapping,
        const flowcell::BarcodeMetadataList &barcodeMetadataList,
        const build::GapRealignerMode realignGaps,
        const alignment::BinMetadata &bin,
        const unsigned binStatsIndex,
        const flowcell::TileMetadataList &tileMetadataList,
        const BuildContigMap &contigMap,
        const isaac::reference::ContigLists &contigLists,
        const unsigned maxReadLength,
        const unsigned char forcedDodgyAlignmentScore,
        const flowcell::FlowcellLayoutList &flowCellLayoutList,
        const IncludeTags includeTags,
        const bool pessimisticMapQ,
        const unsigned splitGapLength) :
            bin_(bin),
            binStatsIndex_(binStatsIndex),
            barcodeBamMapping_(barcodeBamMapping),
            realignGaps_(realignGaps),
            realignerGaps_(getGapGroupsCount()),
            dataDistribution_(bin_.getDataDistribution()),
            inputFileBuf_(1, std::ios_base::binary|std::ios_base::in),
            bamAdapter_(
                maxReadLength, tileMetadataList, barcodeMetadataList,
                contigMap, contigLists, forcedDodgyAlignmentScore, flowCellLayoutList, includeTags, pessimisticMapQ, splitGapLength)
    {
        data_.resize(bin_);

        std::vector<PackedFragmentBuffer::Index>::reserve(bin_.getTotalElements() * 2);
        seIdx_.reserve(bin_.getSeIdxElements());
        rIdx_.reserve(bin_.getRIdxElements());
        fIdx_.reserve(bin_.getFIdxElements());
        if (REALIGN_NONE != realignGaps_)
        {
            reserveGaps(bin_, barcodeMetadataList);
        }
        // assume each split read will not have any noticeable amount of extra elements above of what's listed below
        // translocation gaps with inversion require up to 7 CIGAR components: SOFT_CLIP,ALIGN,FLIP,CONTIG,DELETE,ALIGN,SOFTCLIP
        // so, when this gets broken up, each leftover will have no more than 5 components
        static const unsigned SINGLE_SPLIT_LEFTOVER_COMPONENTS = 5;
        additionalCigars_.reserve(
            // Assuming each split read will result in two separate fragments
            SINGLE_SPLIT_LEFTOVER_COMPONENTS * bin.getTotalSplitCount() * 2 +
            // assume each existing cigar gets realignedGapsPerFragment_ gaps introduced...
            (bin.getTotalCigarLength() + bin.getTotalElements() * realignedGapsPerFragment));

        // summarize chunk sizes to get offsets
        dataDistribution_.tallyOffsets();
        inputFileBuf_.reservePathBuffers(bin_.getPathString().size());
    }

    void finalize();

    static unsigned long getMemoryRequirements(const alignment::BinMetadata& bin)
    {
        return PackedFragmentBuffer::getMemoryRequirements(bin) +
            bin.getSeIdxElements() * sizeof(SeFragmentIndex) +
            bin.getRIdxElements() * sizeof(RStrandOrShadowFragmentIndex) +
            bin.getFIdxElements() * sizeof(FStrandFragmentIndex) +
            bin.getTotalElements() * sizeof(PackedFragmentBuffer::Index);
    }

    void unreserveIndexes()
    {
        std::vector<SeFragmentIndex>().swap(seIdx_);
        std::vector<RStrandOrShadowFragmentIndex>().swap(rIdx_);
        std::vector<FStrandFragmentIndex>().swap(fIdx_);
    }

    unsigned getBinIndex() const
    {
        return bin_.getIndex();
    }

    bool isUnalignedBin() const {return bin_.isUnalignedBin();}
    unsigned long getUniqueRecordsCount() const {return isUnalignedBin() ? bin_.getTotalElements() : size();}

    BaseType::iterator indexBegin() {return begin();}
    BaseType::iterator indexEnd() {return end();}

    const gapRealigner::RealignerGaps &getRealignerGaps(const unsigned barcode) const;
//private:
    const alignment::BinMetadata &bin_;
    const unsigned binStatsIndex_;
    const BarcodeBamMapping &barcodeBamMapping_;
    std::vector<SeFragmentIndex> seIdx_;
    std::vector<RStrandOrShadowFragmentIndex> rIdx_;
    std::vector<FStrandFragmentIndex> fIdx_;
    PackedFragmentBuffer data_;
    const GapRealignerMode realignGaps_;
    alignment::Cigar additionalCigars_;
    std::vector<gapRealigner::RealignerGaps> realignerGaps_;
    alignment::BinDataDistribution dataDistribution_;
    io::FileBufCache<io::FileBufWithReopen> inputFileBuf_;
    FragmentAccessorBamAdapter bamAdapter_;

private:
    void reserveGaps(
        const alignment::BinMetadata& bin,
        const flowcell::BarcodeMetadataList &barcodeMetadataList);

    unsigned getGapGroupIndex(const unsigned barcode) const;
    unsigned getGapGroupsCount() const;

};


} // namespace build
} // namespace isaac

#endif // #ifndef iSAAC_BUILD_BIN_SORTER_HH
