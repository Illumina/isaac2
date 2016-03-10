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
 ** \file FragmentBinner.cpp
 **
 ** \author Roman Petrovski
 **/

#include <cerrno>
#include <fstream>
#include <boost/foreach.hpp>
#include <boost/format.hpp>

#include "common/Debug.hh"
#include "common/Exceptions.hh"
#include "alignment/BinMetadata.hh"
#include "alignment/matchSelector/BinningFragmentStorage.hh"

namespace isaac
{
namespace alignment
{
namespace matchSelector
{

const unsigned FragmentBinner::FRAGMENT_BINS_MAX;
const unsigned FragmentBinner::UNMAPPED_BIN;

FragmentBinner::FragmentBinner(
    const bool keepUnaligned,
    const unsigned maxSavers,
    const BinIndexMap &binIndexMap,
    const alignment::BinMetadataList &binPathList,
    const unsigned long maxTileClusters,
    const unsigned long expectedBinSize):
        keepUnaligned_(keepUnaligned),
        maxTileReads_(maxTileClusters * READS_MAX),
        expectedBinSize_(expectedBinSize),
        binIndexMap_(binIndexMap),
        binZeroClustersBinned_(0),
        files_(maxSavers, io::FileBufWithReopen(std::ios_base::out | std::ios_base::app | std::ios_base::binary)),
        binFiles_(binPathList.size())
{
    ISAAC_ASSERT_MSG(binPathList.back().getIndex() == binPathList.size() - 1, "Expected contiguous bin indexes " << binPathList.back() << " for " << binPathList.size());
}

void FragmentBinner::storeFragment(
    const io::FragmentAccessor &fragment,
    const bool splitRead,
    BinMetadata &binMetadata)
{
    ISAAC_ASSERT_MSG(fragment.flags_.initialized_, "Attempt to store an uninitialized " << fragment);
    if (!fragment.isAligned() && !fragment.isMateAligned())
    {
        // Bin 0 gets split during bam generation. It is important bin 0 chunks reflect distribution in the order
        // in which the records are stored in bin 0. Notice that this is not the case
        // with aligned bins where distribution reflects alignment position
        binMetadata.incrementDataSize(binZeroClustersBinned_, fragment.getTotalLength());
        binMetadata.incrementNmElements(binZeroClustersBinned_, 1, fragment.barcode_);
        ++binZeroClustersBinned_;
    }
    else
    {
        binMetadata.incrementDataSize(fragment.fStrandPosition_, fragment.getTotalLength());
        if (!fragment.flags_.paired_)
        {
            binMetadata.incrementSeIdxElements(fragment.fStrandPosition_, 1, fragment.barcode_);
        }
        else if (fragment.flags_.reverse_ || fragment.flags_.unmapped_)
        {
            binMetadata.incrementRIdxElements(fragment.fStrandPosition_, 1, fragment.barcode_);
        }
        else
        {
            binMetadata.incrementFIdxElements(fragment.fStrandPosition_, 1, fragment.barcode_);
        }

        if (splitRead)
        {
            binMetadata.incrementSplitCount(fragment.fStrandPosition_, fragment.gapCount_, fragment.barcode_);
        }
        else
        {
            binMetadata.incrementGapCount(fragment.fStrandPosition_, fragment.gapCount_, fragment.barcode_);
        }
        binMetadata.incrementCigarLength(fragment.fStrandPosition_, fragment.cigarLength_, fragment.barcode_);
    }

    std::ostream osData(&files_.at(binFiles_.at(binMetadata.getIndex())));
    if (!osData.write(reinterpret_cast<const char*>(&fragment), fragment.getTotalLength())) {
        BOOST_THROW_EXCEPTION(common::IoException(errno, "Failed to write into " + binMetadata.getPathString()));
    }
}

void FragmentBinner::storePaired(
    const io::FragmentAccessor &fragment0,
    const io::FragmentAccessor &fragment1)
{
    FragmentBins bins;
    getFragmentStorageBins(fragment0, bins);
    getFragmentStorageBins(fragment1, bins);

    BOOST_FOREACH(const unsigned binIndex, bins)
    {
        if (!binIndex && (fragment0.isAligned() || fragment1.isAligned()))
        {
            // only when both reads are unaligned, the pair goes into bin 0
            continue;
        }
        const unsigned fileIndex = binFiles_.at(binIndex);
        if (UNMAPPED_BIN != fileIndex)
        {
            alignment::BinMetadata &bin = *(binsBegin_ + fileIndex);
            boost::unique_lock<boost::mutex> lock(binMutex_[fileIndex % binMutex_.size()]);
            // make sure orphan cometh first in shadow/orphan pair
            if (fragment0.isAligned())
            {
                storeFragment(fragment0, fragment0.flags_.splitAlignment_, bin);
                ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment0.clusterId_, "BinningFragmentStorage::storePaired: " << fragment0);
            }
            storeFragment(fragment1, 0 != binIndex && fragment1.flags_.splitAlignment_, bin);
            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment1.clusterId_, "BinningFragmentStorage::storePaired: " << fragment1);
            if (!fragment0.isAligned())
            {
                storeFragment(fragment0, false, bin);
                ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment0.clusterId_, "BinningFragmentStorage::storePaired: " << fragment0);
            }
        }
    }
}

void FragmentBinner::storeSingle(
    const io::FragmentAccessor &fragment)
{
    FragmentBins bins;
    getFragmentStorageBins(fragment, bins);

    BOOST_FOREACH(const unsigned binIndex, bins)
    {
        const unsigned fileIndex = binFiles_.at(binIndex);
        if (UNMAPPED_BIN != fileIndex)
        {
            boost::unique_lock<boost::mutex> lock(binMutex_[fileIndex % binMutex_.size()]);
            storeFragment(fragment, 0 != binIndex && fragment.flags_.splitAlignment_, *(binsBegin_ + fileIndex));
        }
    }
}

void FragmentBinner::getFragmentStorageBins(const io::FragmentAccessor &fragment, FragmentBins &bins)
{
    if (!fragment.isAligned())
    {
        bins.push_back(0);
    }
    else
    {
        for (CigarPosition<const unsigned *> it(
            fragment.cigarBegin(), fragment.cigarEnd(), fragment.getFStrandReferencePosition(), fragment.isReverse(), fragment.readLength_);
            !it.end(); ++it)
        {
            if (Cigar::ALIGN == it.component().second)
            {
                const unsigned startBinIndex = binIndexMap_.getBinIndex(it.referencePos_);
                if (bins.empty() || bins.back() != startBinIndex)
                {
                    bins.push_back(startBinIndex);
                }
                // push the last position as well. This is important for duplicate detection of r-stranded alignments that end in the same bin but begin in different ones
                const unsigned endBinIndex = binIndexMap_.getBinIndex(it.referencePos_ + it.component().first - 1);
                if (bins.back() != endBinIndex)
                {
                    bins.push_back(endBinIndex);
                }
            }
        }
    }
    std::sort(bins.begin(), bins.end());
    bins.erase(std::unique(bins.begin(), bins.end()), bins.end());
}

void FragmentBinner::open(
    const alignment::BinMetadataList::iterator binsBegin,
    const alignment::BinMetadataList::iterator binsEnd)
{
    ISAAC_ASSERT_MSG(std::size_t(std::distance(binsBegin, binsEnd)) <= files_.size(),
                     "Requested more than maximum number of files to be written at the same time:" << std::distance(binsBegin, binsEnd) << " max: " << files_.size());
    ISAAC_THREAD_CERR << "Reopening output files for " << std::distance(binsBegin, binsEnd) << " bins" << std::endl;

    std::size_t file = 0;
    std::fill(binFiles_.begin(), binFiles_.end(), UNMAPPED_BIN);
    BOOST_FOREACH(const BinMetadata &binMetadata, std::make_pair(binsBegin, binsEnd))
    {
        if (binMetadata.isEmpty())
        {
            // make sure file is empty first time we decide to put data in it.
            // boost::filesystem::remove for some stupid reason needs to allocate strings for this...
            unlink(binMetadata.getPath().c_str());
            files_[file].reopen(binMetadata.getPath().c_str(), expectedBinSize_, io::FileBufWithReopen::SequentialOnce);
        }
        else
        {
            files_[file].reopen(binMetadata.getPath().c_str(), io::FileBufWithReopen::SequentialOnce);
        }

        if (!files_[file].is_open())
        {
            BOOST_THROW_EXCEPTION(common::IoException(errno, "Failed to open bin file " + binMetadata.getPathString()));
        }
        binFiles_.at(binMetadata.getIndex()) = file;
        ++file;
    }

    binsBegin_ = binsBegin;
    binsEnd_ = binsEnd;

    ISAAC_THREAD_CERR << "Reopening output files done for " << std::distance(binsBegin, binsEnd) << " bins" << std::endl;
}

void FragmentBinner::close(
    const alignment::BinMetadataList::iterator binsBegin,
    const alignment::BinMetadataList::iterator binsEnd)
{
    ISAAC_THREAD_CERR << "truncating output files for " << std::distance(binsBegin, binsEnd) << " bins" << std::endl;

    // make sure everything is written out for those that are open
    std::for_each(files_.begin(), files_.end(), boost::bind(&io::FileBufWithReopen::close, _1));

    std::fill(binFiles_.begin(), binFiles_.end(), UNMAPPED_BIN);

    // if some files were open with fallocate, trim them back to their sizes to avoid wasted storage
    BOOST_FOREACH(const BinMetadata &binMetadata, std::make_pair(binsBegin, binsEnd))
    {
        truncate(binMetadata.getPath().c_str(), binMetadata.getDataSize());
    }

    ISAAC_THREAD_CERR << "truncating output files done for " << std::distance(binsBegin, binsEnd) << " bins" << std::endl;
}

} //namespace matchSelector
} // namespace alignment
} // namespace isaac
