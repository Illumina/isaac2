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
 ** \file FragmentBinner.hh
 **
 ** \brief Stores fragments in bin files.
 ** 
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_ALIGNMENT_MATCH_SELECTOR_FRAGMENT_BINNER_HH
#define iSAAC_ALIGNMENT_MATCH_SELECTOR_FRAGMENT_BINNER_HH

#include <boost/noncopyable.hpp>
#include <boost/filesystem.hpp>
#include <boost/foreach.hpp>

#include "alignment/BinMetadata.hh"
#include "BinIndexMap.hh"
#include "io/FileBufCache.hh"
#include "io/Fragment.hh"


namespace isaac
{
namespace alignment
{
namespace matchSelector
{

namespace bfs = boost::filesystem;

class FragmentPacker: boost::noncopyable
{
public:
    template <typename InsertIT>
    static InsertIT packPairedFragment(
        const alignment::BamTemplate &bamTemplate,
        const unsigned fragmentIndex,
        const unsigned barcodeIdx,
        const BinIndexMap &binIndexMap,
        InsertIT insertIt)
    {
        ISAAC_ASSERT_MSG(READS_MAX == bamTemplate.getFragmentCount(), "Expected paired data");

        const alignment::FragmentMetadata &fragment = bamTemplate.getFragmentMetadata(fragmentIndex);
        const alignment::FragmentMetadata &mate = bamTemplate.getMateFragmentMetadata(fragment);

        const unsigned mateStorageBin = mate.isNoMatch() && fragment.isNoMatch() ?
            0 :
            binIndexMap.getBinIndex(mate.getFStrandReferencePosition());

        const io::FragmentHeader header(bamTemplate, fragment, mate, barcodeIdx, mateStorageBin);
        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(header.clusterId_, "FragmentPacker::packPairedFragment " << header << " " << fragment);

        insertIt = std::copy(header.bytesBegin(), header.bytesEnd(), insertIt);
        insertIt = storeBclAndCigar(fragment, insertIt);
        insertIt = std::copy(bamTemplate.nameBegin(), bamTemplate.nameEnd(), insertIt);
        *insertIt++ = 0;

        return insertIt;
    }

    template <typename InsertIT>
    static InsertIT packSingleFragment(
        const alignment::BamTemplate &bamTemplate,
        const unsigned barcodeIdx,
        InsertIT insertIt)
    {
        ISAAC_ASSERT_MSG(1 == bamTemplate.getFragmentCount(), "Expected single-ended data");

        const alignment::FragmentMetadata &fragment = bamTemplate.getFragmentMetadata(0);
        const io::FragmentHeader header(bamTemplate, fragment, barcodeIdx);
        insertIt = std::copy(header.bytesBegin(), header.bytesEnd(), insertIt);
        insertIt = storeBclAndCigar(fragment, insertIt);
        insertIt = std::copy(bamTemplate.nameBegin(), bamTemplate.nameEnd(), insertIt);
        *insertIt++ = 0;

        return insertIt;
    }

private:
    static const unsigned READS_MAX = 2;

    template <typename InsertIT>
    static InsertIT storeBclAndCigar(
        const alignment::FragmentMetadata & fragment,
        InsertIT insertIt)
    {
        // copy the bcl data (reverse-complement the sequence if the fragment is reverse-aligned)
        std::vector<char>::const_iterator bclData = fragment.getBclData();
        if (fragment.isReverse())
        {
            insertIt = std::transform(std::reverse_iterator<std::vector<char>::const_iterator>(bclData + fragment.getReadLength()),
                                          std::reverse_iterator<std::vector<char>::const_iterator>(bclData),
                                          insertIt, oligo::getReverseBcl);
        }
        else
        {
            insertIt = std::copy(bclData, bclData + fragment.getReadLength(), insertIt);

        }

        if (fragment.isAligned())
        {
            const alignment::Cigar::const_iterator cigarBegin = fragment.cigarBuffer->begin() + fragment.cigarOffset;
            const alignment::Cigar::const_iterator cigarEnd = cigarBegin + fragment.cigarLength;
            BOOST_FOREACH(const unsigned cig, std::make_pair(cigarBegin, cigarEnd))
            {
                *insertIt++ =cig;
                *insertIt++ =(cig >> 8);
                *insertIt++ =(cig >> 16);
                *insertIt++ =(cig >> 24);
            }
        }
        return insertIt;
    }
};

class FragmentBinner: boost::noncopyable
{
public:
    FragmentBinner(
        const bool keepUnaligned,
        const unsigned maxSavers,
        const BinIndexMap &binIndexMap,
        const alignment::BinMetadataList &binMetadataList,
        const unsigned long maxTileClusters,
        const unsigned long expectedBinSize);

    void open(
        const alignment::BinMetadataList::iterator binsBegin,
        const alignment::BinMetadataList::iterator binsEnd);

    /**
     * \brief reclaims any unused storage.
     */
    void close(
        const alignment::BinMetadataList::iterator binsBegin,
        const alignment::BinMetadataList::iterator binsEnd);

    void storeSingle(
        const io::FragmentAccessor &fragment);

    void storePaired(
        const io::FragmentAccessor &fragment0,
        const io::FragmentAccessor &fragment1);

private:
    /// Maximum number of bins a fragment is expected to cover. In theory this can be up to total number of bins.
    static const unsigned FRAGMENT_BINS_MAX = 10*1024;
    static const unsigned UNMAPPED_BIN = -1U;

    static const unsigned READS_MAX = 2;
    const bool keepUnaligned_;
    const unsigned long maxTileReads_;
    const unsigned long expectedBinSize_;

    const BinIndexMap &binIndexMap_;

    unsigned long binZeroClustersBinned_;
    // number of mutexes is arbitrary. It allows reducing the synchronization collisions between threads requesting write
    boost::array<boost::mutex, 8> binMutex_;
    std::vector<io::FileBufWithReopen> files_;
    std::vector<unsigned> binFiles_;

    alignment::BinMetadataList::iterator binsBegin_;
    alignment::BinMetadataList::iterator binsEnd_;

    void storeFragment(
        const io::FragmentAccessor &fragment,
        const bool splitRead,
        BinMetadata &binMetadata);

    typedef common::FiniteCapacityVector<unsigned, FRAGMENT_BINS_MAX * 2> FragmentBins;
    void getFragmentStorageBins(const io::FragmentAccessor &fragment, FragmentBins &bins);
};

} // namespace matchSelector
} // namespace alignment
} // namespace isaac

#endif // #ifndef iSAAC_ALIGNMENT_MATCH_SELECTOR_FRAGMENT_BINNER_HH
