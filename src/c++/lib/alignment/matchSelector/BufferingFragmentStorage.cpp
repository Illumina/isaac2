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
 ** \file BufferingFragmentStorage.cpp
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
#include "alignment/matchSelector/BufferingFragmentStorage.hh"
#include "io/Fragment.hh"

namespace isaac
{
namespace alignment
{
namespace matchSelector
{

BufferingFragmentStorage::BufferingFragmentStorage(
    const bool keepUnaligned,
    const unsigned maxSavers,
    const BinIndexMap &binIndexMap,
    alignment::BinMetadataList &binMetadataList,
    const flowcell::FlowcellLayoutList &flowcellLayoutList,
    const unsigned long maxTileClusters,
    const unsigned long expectedBinSize)
    : FragmentBinner(keepUnaligned, maxSavers, binIndexMap, binMetadataList, maxTileClusters, expectedBinSize)
    , keepUnaligned_(keepUnaligned)
    , maxTileReads_(maxTileClusters * READS_MAX)
    , binIndexMap_(binIndexMap)
    , maxSavers_(maxSavers)
    , binMetadataList_(binMetadataList)
    , flushBuffer_(maxTileClusters, flowcellLayoutList)
    , storeBuffer_(maxTileClusters, flowcellLayoutList)
{
}

void BufferingFragmentStorage::store(const BamTemplate &bamTemplate, const unsigned barcodeIdx)
{
    const alignment::FragmentMetadata &fragment = bamTemplate.getFragmentMetadata(0);
    if (2 == bamTemplate.getFragmentCount())
    {
        FragmentPacker::packPairedFragment(bamTemplate, 0, barcodeIdx, binIndexMap_,
                                           storeBuffer_.getRecordInsertIterator(fragment.getCluster().getId(), 0));
        FragmentPacker::packPairedFragment(bamTemplate, 1, barcodeIdx, binIndexMap_,
                                           storeBuffer_.getRecordInsertIterator(fragment.getCluster().getId(), 1));
    }
    else
    {
        FragmentPacker::packSingleFragment(bamTemplate, barcodeIdx,
                                           storeBuffer_.getRecordInsertIterator(fragment.getCluster().getId(), 0));
    }
}

void BufferingFragmentStorage::prepareFlush()
{
    ISAAC_ASSERT_MSG(flushBuffer_.empty(), "Buffer must be flushed at this point");
    storeBuffer_.swap(flushBuffer_);
}

void BufferingFragmentStorage::flush()
{
    ISAAC_THREAD_CERR << "Flushing buffer" << std::endl;

    for (alignment::BinMetadataList::iterator binIterator = binMetadataList_.begin();
        binMetadataList_.end() != binIterator; )
    {
        const unsigned binsToFlush = std::min<unsigned>(maxSavers_, std::distance(binIterator, binMetadataList_.end()));
        open(binIterator, binIterator + binsToFlush);
        ISAAC_THREAD_CERR << "Flushing bins " << std::distance(binMetadataList_.begin(), binIterator) << "-" << std::distance(binMetadataList_.begin(), binIterator) + binsToFlush << std::endl;

        for (FragmentBuffer::const_iterator it = flushBuffer_.dataBegin(); flushBuffer_.dataEnd() != it;
            it = flushBuffer_.nexCluster(it))
        {
            const io::FragmentAccessor &fragment0 = *reinterpret_cast<const io::FragmentAccessor*>(&*it);
            if(fragment0.flags_.paired_)
            {
                const FragmentBuffer::const_iterator r2It = flushBuffer_.nextRead(it);
                ISAAC_ASSERT_MSG(flushBuffer_.dataEnd() != r2It, "Unexpected end of buffer reached");
                const io::FragmentAccessor &fragment1 = *reinterpret_cast<const io::FragmentAccessor*>(&*r2It);
                ISAAC_ASSERT_MSG(fragment0.flags_.initialized_ == fragment1.flags_.initialized_,
                                 "Both reads have to be either initialized or not: " << fragment0 << " " << fragment1);
                if (fragment0.flags_.initialized_)
                {
                    storePaired(fragment0, fragment1);
                }
            }
            else
            {
                if (fragment0.flags_.initialized_)
                {
                    storeSingle(fragment0);
                }
            }
        }

        binIterator += binsToFlush;
    }

    flushBuffer_.clear();

    ISAAC_THREAD_CERR << "Flushing buffer done " << std::endl;
}


} //namespace matchSelector
} // namespace alignment
} // namespace isaac
