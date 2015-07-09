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
 ** \file BufferingFragmentStorage.hh
 **
 ** \brief Fragment buffer flushing and output file management.
 ** 
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_ALIGNMENT_MATCH_SELECTOR_BUFFERING_FRAGMENT_STORAGE_HH
#define iSAAC_ALIGNMENT_MATCH_SELECTOR_BUFFERING_FRAGMENT_STORAGE_HH

#include <boost/noncopyable.hpp>
#include <boost/filesystem.hpp>
#include <boost/foreach.hpp>
#include <boost/thread.hpp>

#include "alignment/MatchDistribution.hh"
#include "BinIndexMap.hh"
#include "common/Threads.hpp"
#include "FragmentBinner.hh"
#include "FragmentBuffer.hh"
#include "FragmentStorage.hh"
#include "io/FileBufCache.hh"


namespace isaac
{
namespace alignment
{
namespace matchSelector
{

namespace bfs = boost::filesystem;

class BufferingFragmentStorage: FragmentBinner, public FragmentStorage
{
public:
    BufferingFragmentStorage(
        const bool keepUnaligned,
        const unsigned maxSavers,
        const BinIndexMap &binIndexMap,
        alignment::BinMetadataList &binMetadataList,
        const flowcell::FlowcellLayoutList &flowcellLayoutList,
        const unsigned long maxTileClusters,
        const unsigned long expectedBinSize);

    using FragmentBinner::close;

    virtual void store(const BamTemplate &bamTemplate, const unsigned barcodeIdx);
    virtual void prepareFlush();
    virtual void flush();
    virtual void resize(const unsigned long clusters)
    {
        storeBuffer_.clear();
        storeBuffer_.resize(clusters);
    }
    virtual void unreserve()
    {
        flushBuffer_.unreserve();
        storeBuffer_.swap(flushBuffer_);
        flushBuffer_.unreserve();
    }

private:
    static const unsigned READS_MAX = 2;

    const bool keepUnaligned_;
    const unsigned long maxTileReads_;

    const BinIndexMap &binIndexMap_;
    const unsigned maxSavers_;

    /// association of a bin index to a path
    alignment::BinMetadataList &binMetadataList_;
    matchSelector::FragmentBuffer flushBuffer_;
    matchSelector::FragmentBuffer storeBuffer_;
};

} // namespace matchSelector
} // namespace alignment
} // namespace isaac

#endif // #ifndef iSAAC_ALIGNMENT_MATCH_SELECTOR_BUFFERING_FRAGMENT_STORAGE_HH
