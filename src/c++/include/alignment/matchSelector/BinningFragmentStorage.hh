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
 ** \file BinningFragmentStorage.hh
 **
 ** \brief Stores fragments in bin files without buffering.
 ** 
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_ALIGNMENT_MATCH_SELECTOR_BINNING_FRAGMENT_STORAGE_HH
#define iSAAC_ALIGNMENT_MATCH_SELECTOR_BINNING_FRAGMENT_STORAGE_HH

#include <boost/noncopyable.hpp>
#include <boost/filesystem.hpp>
#include <boost/foreach.hpp>
#include <boost/thread.hpp>

#include "alignment/MatchDistribution.hh"
#include "BinIndexMap.hh"
#include "common/Threads.hpp"
#include "FragmentStorage.hh"
#include "FragmentBinner.hh"
#include "io/FileBufCache.hh"
#include "io/Fragment.hh"


namespace isaac
{
namespace alignment
{
namespace matchSelector
{

namespace bfs = boost::filesystem;

class BinningFragmentStorage: FragmentPacker, FragmentBinner, public FragmentStorage
{
public:
    BinningFragmentStorage(
        const bool keepUnaligned,
        const unsigned maxSavers,
        const BinIndexMap &binIndexMap,
        const alignment::BinMetadataList &binMetadataList,
        const unsigned long maxTileClusters,
        const unsigned long expectedBinSize);

    using FragmentBinner::open;
    using FragmentBinner::close;

    virtual void store(
        const BamTemplate &bamTemplate,
        const unsigned barcodeIdx);

    virtual void prepareFlush()
    {
    }
    virtual void flush()
    {
    }
    virtual void resize(const unsigned long clusters)
    {
    }
    virtual void unreserve()
    {
    }

private:
    /// Maximum number of bytes a packed fragment is expected to take. Change and recompile when needed
    static const unsigned FRAGMENT_BYTES_MAX = 10*1024;
    static const unsigned READS_MAX = 2;
    const BinIndexMap &binIndexMap_;
};

} // namespace matchSelector
} // namespace alignment
} // namespace isaac

#endif // #ifndef iSAAC_ALIGNMENT_MATCH_SELECTOR_BINNING_FRAGMENT_STORAGE_HH
