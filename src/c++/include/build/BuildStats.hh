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
 ** \file BuildStats.hh
 **
 ** \brief Build statistics helper.
 **
 ** \author Roman Petrovski
 **/

#ifndef ISAAC_BUILD_BUILD_STATS_H
#define ISAAC_BUILD_BUILD_STATS_H

#include "alignment/BinMetadata.hh"

namespace isaac
{
namespace build
{

struct BinBarcodeStats
{
    BinBarcodeStats() : totalFragments_(0), uniqueFragments_(0){}
    unsigned long totalFragments_;
    unsigned long uniqueFragments_;

    BinBarcodeStats &operator +=(const BinBarcodeStats &right)
    {
        totalFragments_ += right.totalFragments_;
        uniqueFragments_ += right.uniqueFragments_;
        return *this;
    }
};

inline const BinBarcodeStats operator +(BinBarcodeStats left, const BinBarcodeStats &right)
{
    left += right;
    return left;
}

inline unsigned highestBinIndex(const alignment::BinMetadataCRefList &binMetadataList)
{
    unsigned ret = 0;
    BOOST_FOREACH(const alignment::BinMetadata &bin, binMetadataList)
    {
        ret = std::max(ret, bin.getIndex());
    }
    return ret;
}

class BuildStats
{
public:

    BuildStats(
        const alignment::BinMetadataCRefList &binMetadataList,
        const flowcell::BarcodeMetadataList &barcodeMetadataList) :
            barcodeMetadataList_(barcodeMetadataList),
            binBarcodeStats_(barcodeMetadataList_.size() * binMetadataList.size())
    {
    }

    void incrementTotalFragments(
        const unsigned binIndex,
        const unsigned barcodeIndex)
    {
        ++binBarcodeStats_.at(binBarcodeIndex(binIndex, barcodeIndex)).totalFragments_;
    }

    void incrementUniqueFragments(
        const unsigned binIndex,
        const unsigned barcodeIndex)
    {
        ++binBarcodeStats_.at(binBarcodeIndex(binIndex, barcodeIndex)).uniqueFragments_;
    }

    unsigned long getTotalFragments(
        const unsigned binIndex,
        const unsigned barcodeIndex) const
    {
        return binBarcodeStats_.at(binBarcodeIndex(binIndex, barcodeIndex)).totalFragments_;
    }

    unsigned long getUniqueFragments(
        const unsigned binIndex,
        const unsigned barcodeIndex) const
    {
        return binBarcodeStats_.at(binBarcodeIndex(binIndex, barcodeIndex)).uniqueFragments_;
    }

    BuildStats &operator +=(const BuildStats &right)
    {
        std::transform(binBarcodeStats_.begin(), binBarcodeStats_.end(),
                       right.binBarcodeStats_.begin(), binBarcodeStats_.begin(), std::plus<BinBarcodeStats>());
        return *this;
    }

    BuildStats & operator =(const BuildStats &that) {
        binBarcodeStats_ = that.binBarcodeStats_;
        return *this;
    }

private:
    const flowcell::BarcodeMetadataList &barcodeMetadataList_;
    std::vector<BinBarcodeStats>  binBarcodeStats_;

    unsigned binBarcodeIndex(const unsigned binIndex, const unsigned barcodeIndex) const
    {
        return binIndex * barcodeMetadataList_.size() + barcodeIndex;
    }
};

inline const BuildStats operator +(BuildStats left, const BuildStats &right)
{
    left += right;
    return left;
}


} //namespace build
} //namespace isaac

#endif //ISAAC_BUILD_BUILD_STATS_H
