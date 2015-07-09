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
 ** \file NeighborsFinder.hh
 **
 ** \brief Top level component to find neighbors.
 **
 ** \author Come Raczy
 **/

#ifndef ISAAC_REFERENCE_NEIGHBORS_FINDER_HH
#define ISAAC_REFERENCE_NEIGHBORS_FINDER_HH

#include <vector>
#include <boost/filesystem.hpp>
#include <boost/noncopyable.hpp>
#include <boost/ptr_container/ptr_list.hpp>

#include "oligo/Kmer.hh"
#include "oligo/Permutate.hh"
#include "reference/ContigLoader.hh"
#include "reference/SortedReferenceMetadata.hh"
#include "reference/neighborsFinder/AnnotatedKmer.hh"


namespace isaac
{
namespace reference
{

template <typename KmerT, typename AnnotatorT>
class NeighborsFinder
{
public:
    NeighborsFinder(
        AnnotatorT &annotator,
        const unsigned maskWidth,
        const unsigned neighborhoodWidth,
        common::ThreadVector &threads);

    template <typename ReferenceKmerT>
    const typename AnnotatorT::Annotation &annotate();
private:
    AnnotatorT &annotator_;
    const unsigned maskWidth_;
    const unsigned neighborhoodWidth_;
    const KmerT permutedX_;

    typedef typename AnnotatorT::KmerList KmerList;
    typedef typename AnnotatorT::AnnotatedKmer AnnotatedKmer;
    common::ThreadVector &threads_;
    // permutation list for the whole genome
    const std::vector<oligo::Permutate> permutateList_;
    // local permutation list to be used when the same suffix block length is too big to run n^2 neighbor finding on it
    const std::vector<oligo::Permutate> allSuffixPermutateList_;

    mutable boost::mutex mutex_;
    boost::condition_variable stateChangedCondition_;

    void annotateNeighbors(
        const std::vector<std::vector<KmerT> > &mismatchBlockMasks,
        const std::vector<oligo::Permutate> &suffixPermutateList,
        KmerList &kmerList);

    void annotateBlock(
        const std::vector<std::vector<KmerT> > &mismatchBlockMasks,
        const std::vector<oligo::Permutate> &suffixPermutateList,
        typename KmerList::iterator blockBegin,
        typename KmerList::iterator blockEnd) const;


    struct AnnotateTask
    {
        typename KmerList::iterator begin_;
        typename KmerList::iterator end_;
        AnnotateTask(typename KmerList::iterator begin, typename KmerList::iterator end) : begin_(begin), end_(end){}
    };
    typedef std::vector<AnnotateTask> Tasks;


    void processTasks(
        Tasks &tasks,
        const std::vector<std::vector<KmerT> > &mismatchBlockMasks,
        const std::vector<oligo::Permutate> &suffixPermutateList,
        const unsigned suffixLength,
        const std::size_t threadNum);

    void annotateBlocks(
        const std::vector<std::vector<KmerT> > &mismatchBlockMasks,
        const std::vector<oligo::Permutate> &suffixPermutateList,
        const unsigned suffixLength,
        typename KmerList::iterator ourBegin,
        typename KmerList::iterator ourEnd) const;

    static typename KmerList::iterator findBlockEnd(
        typename KmerList::iterator randomPosition,
        const unsigned prefixBits,
        const typename KmerList::iterator kmerListEnd);

    std::size_t markNeighbors(
        const typename KmerList::iterator blockBegin,
        const typename KmerList::iterator blockEnd) const;

    std::size_t markNeighbors(
        const std::vector<std::vector<KmerT> > &mismatchBlockMasks,
        const typename KmerList::iterator blockBegin,
        const typename KmerList::iterator blockEnd) const;


    static void permuateAnnotatedKmerSuffix(
        const oligo::Permutate &permutate,
        AnnotatedKmer &kmer)
    {
        kmer.kmer_ = permutate(typename oligo::KmerTraits<KmerT>::SuffixType(kmer.kmer_));
    }

    std::size_t annotatedSubblock(
        const typename KmerList::iterator& blockBegin,
        const typename KmerList::iterator& blockEnd) const;

    void annotate(
        const oligo::Permutate &permutate,
        KmerList &kmerList);

    void reorderKmers(
        const oligo::Permutate &permutate, KmerList &kmerList, const std::size_t threadNum) const;
    void permutateKmers(
        const oligo::Permutate &permutate, KmerList &kmerList, const std::size_t threadNum) const;
};

template <typename KmerT, typename AnnotatorT>
NeighborsFinder<KmerT, AnnotatorT>::NeighborsFinder(
    AnnotatorT &annotator,
    const unsigned maskWidth,
    const unsigned neighborhoodWidth,
    common::ThreadVector &threads)
    : annotator_(annotator)
    , maskWidth_(maskWidth)
    , neighborhoodWidth_(neighborhoodWidth)
    , permutedX_(~(~KmerT(0) << (neighborhoodWidth_ * oligo::BITS_PER_BASE)))
    , threads_(threads)
    // neighborhoodWidth can occur with 2,4,...,neighborhoodWidth*2 blocks.
    // Make sure all relevant partitioning is present in such way that no two
    // permutations detect the same pair of mismatching neighbors.
    , permutateList_(oligo::getPermutateList<KmerT>(neighborhoodWidth, maskWidth_, true))
    , allSuffixPermutateList_(oligo::getPermutateList(
        oligo::KmerTraits<typename oligo::KmerTraits<KmerT>::SuffixType>::KMER_BASES, 1,
        neighborhoodWidth, 0, false))
{
}

template <typename KmerT, typename PermutaeT>
KmerT permutateSuffix(const KmerT &ori, const PermutaeT &permutate)
{
    const typename oligo::KmerTraits<KmerT>::SuffixType suffix(ori);
    const typename oligo::KmerTraits<KmerT>::SuffixType permutedSuffix = permutate(suffix);
    static const unsigned SUFFIX_BITS = oligo::KmerTraits<KmerT>::SUFFIX_BITS;
    const KmerT permuted = (ori & (~KmerT(0) << SUFFIX_BITS)) | permutedSuffix;

    return permuted;
}

/**
 * \brief check if the every block that must have a mismatch has it and none of those that must match do
 */
template <typename KmerT>
bool mismatchesInCorrectBlocks(const KmerT x, const KmerT matchMask, const std::vector<KmerT> &blockMasks)
{
    BOOST_FOREACH(const KmerT blockMask, blockMasks)
    {
        // if block must have a mismatch, make sure it does. Exception is a "full match" entry which
        // can't have any mismatches anywhere
        if (!!blockMask && !(x & blockMask))
        {
            return false;
        }
    }

    return !(matchMask & x);
}

template <typename KmerT, typename AnnotatorT>
void NeighborsFinder<KmerT, AnnotatorT>::annotate(
    const oligo::Permutate &permutate,
    KmerList &kmerList)
{
//    ISAAC_ASSERT_MSG(permutate.getMismatchMask(), "Permutation is expected to have mismatch mask");
    ISAAC_THREAD_CERR << " Sorting all k-mers (" << kmerList.size() << " k-mers)" << std::endl;
    common::parallelSort(kmerList, &AnnotatedKmer::kmerLess);
//            std::sort(kmerList.begin(), kmerList.end(), &AnnotatedKmer::kmerLess);
    ISAAC_THREAD_CERR << " Sorting all k-mers done " << std::endl;

    // skip the expensive matching for permutation that does not allow neighbors
    if (!permutate.repeatsOnly())
    {
        // one mask per block so that each block can be tested independently
        std::vector<std::vector<KmerT> > mismatchBlockMasks;

        BOOST_FOREACH(unsigned mismatchMask, permutate.mismatchMasks())
        {
            ISAAC_THREAD_CERR << "mismatchMask:" << mismatchMask << std::endl;
            KmerT blockMask = ~(~KmerT(0) << permutate.getBlockLength() * oligo::BITS_PER_BASE);
            mismatchBlockMasks.resize(mismatchBlockMasks.size() + 1);
            while (mismatchMask)
            {
                if (mismatchMask & 0x01)
                {
                    mismatchBlockMasks.back().push_back(blockMask);
                    ISAAC_THREAD_CERR << "blockMask:" << blockMask << std::endl;
                }
                mismatchMask >>=1;
                blockMask <<= permutate.getBlockLength() * oligo::BITS_PER_BASE;
            }
        }

        std::vector<oligo::Permutate> suffixPermutateList;
        oligo::Permutate::Order from = allSuffixPermutateList_.front().getFromOrder();

        BOOST_FOREACH(const oligo::Permutate &suffixPermutate, allSuffixPermutateList_)
        {
            // check if this permutation is compatible with outer permutation
            const KmerT x = permutateSuffix(
                permutedX_,
                boost::bind(&oligo::Permutate::reorder<typename oligo::KmerTraits<KmerT>::SuffixType>,
                            suffixPermutate, _1));

            bool skip = true;
            BOOST_FOREACH(const std::vector<KmerT> &blockMasks, mismatchBlockMasks)
            {
                KmerT allMismatchesMask(0);
                BOOST_FOREACH(const KmerT blockMask, blockMasks)
                {
                    allMismatchesMask |= blockMask;
                }
                const KmerT matchMask = ~allMismatchesMask;

                if (mismatchesInCorrectBlocks(x, matchMask, blockMasks))
                {
                    // local permutation is compatible with the outer permutation. use it to count neighbors
                    skip = false;
                    break;
                }
            }
            if (!skip)
            {
                suffixPermutateList.push_back(
                    oligo::Permutate(suffixPermutate.getBlockLength(),
                                     from, suffixPermutate.getToOrder(), ~(~unsigned(0) << neighborhoodWidth_)));
                from = suffixPermutate.getToOrder();
            }
        }
//
//        longPermutateList_(oligo::getPermutateList(
//            oligo::KmerTraits<typename oligo::KmerTraits<KmerT>::SuffixType>::KMER_BASES, 1,
//            neighborhoodWidth, true, false))

        annotateNeighbors(mismatchBlockMasks, suffixPermutateList, kmerList);
    }
}

template <typename KmerT, typename AnnotatorT>
void NeighborsFinder<KmerT, AnnotatorT>::reorderKmers(
    const oligo::Permutate &permutate, KmerList &kmerList, const std::size_t threadNum) const
{
    const std::size_t blockLength = ((kmerList.size() + threads_.size() - 1) / threads_.size());
    if (kmerList.size() > blockLength * threadNum)
    {
        typename KmerList ::iterator begin = kmerList.begin() + blockLength * threadNum;
        typename KmerList ::iterator end = std::size_t(std::distance(begin, kmerList.end())) < blockLength ?
            kmerList.end() : begin + blockLength;
        BOOST_FOREACH(AnnotatedKmer &ak, std::make_pair(begin, end))
        {
            ak.kmer_ = permutate.reorder(ak.kmer_);
        }
    }
}

template <typename KmerT, typename AnnotatorT>
void NeighborsFinder<KmerT, AnnotatorT>::permutateKmers(
    const oligo::Permutate &permutate, KmerList &kmerList, const std::size_t threadNum) const
{
    const std::size_t blockLength = ((kmerList.size() + threads_.size() - 1) / threads_.size());
    if (kmerList.size() > blockLength * threadNum)
    {
        typename KmerList ::iterator begin = kmerList.begin() + blockLength * threadNum;
        typename KmerList ::iterator end = std::size_t(std::distance(begin, kmerList.end())) < blockLength ?
            kmerList.end() : begin + blockLength;
        BOOST_FOREACH(AnnotatedKmer &ak, std::make_pair(begin, end))
        {
            ak.kmer_ = permutate(ak.kmer_);
        }
    }
}

template <typename KmerT, typename AnnotatorT>
template <typename ReferenceKmerT>
const typename AnnotatorT::Annotation &NeighborsFinder<KmerT, AnnotatorT>::annotate()
{
    // iterate over all possible permutations
    KmerList kmerList;
    unsigned permutation = 0;
    BOOST_FOREACH(const oligo::Permutate &permutate, permutateList_)
    {
        ++permutation;
        ISAAC_THREAD_CERR << "Processing permutation " << permutate << " (" << permutation << "/" << permutateList_.size() << ")" << std::endl;

        const bool cantChain = (1 == permutation) || !permutate.isSamePrefix(maskWidth_, permutateList_.at(permutation - 2));
        if (kmerList.empty() || cantChain)
        {
            kmerList.clear();
            annotator_.getKmers(permutate, kmerList, threads_);
        }
        else
        {
            // when noMask_ is set, we have all the data, just apply next permuatation instead of regenerating kmers from linear reference
            ISAAC_THREAD_CERR << "Permutating all k-mers (" << kmerList.size() << " k-mers)" << std::endl;
            threads_.execute(boost::bind(&NeighborsFinder::permutateKmers, this, boost::ref(permutate), boost::ref(kmerList), _1));
            ISAAC_THREAD_CERR << "Permutating all k-mers done " << std::endl;
        }

        annotate(permutate, kmerList);

        // dump counts on last permutation or if the next permutation will use different set of kmers
        if (permutateList_.size() == permutation || !permutate.isSamePrefix(maskWidth_, permutateList_.at(permutation)))
        {
            if (!permutate.repeatsOnly())
            {
                ISAAC_THREAD_CERR << "Reordering all k-mers (" << kmerList.size() << " k-mers)" << std::endl;
                threads_.execute(boost::bind(&NeighborsFinder::reorderKmers, this, boost::ref(permutate), boost::ref(kmerList), _1));
                ISAAC_THREAD_CERR << "Reordering all k-mers done " << std::endl;

                ISAAC_THREAD_CERR << "Sorting all k-mers (" << kmerList.size() << " k-mers)" << std::endl;
                common::parallelSort(kmerList, &AnnotatedKmer::kmerLess);
                ISAAC_THREAD_CERR << "Sorting all k-mers done " << std::endl;
            }

            annotator_.template updateAnnotation<ReferenceKmerT>(kmerList, permutate.repeatsOnly());
        }
    }

    return annotator_.getAnnotation();
}

template <typename KmerT>
bool isReverseComplement(const KmerT &kmer)
{
    return kmer.isReverseComplement();
}


template <typename KmerT, typename AnnotatorT>
void NeighborsFinder<KmerT, AnnotatorT>::processTasks(
    Tasks &tasks,
    const std::vector<std::vector<KmerT> > &mismatchBlockMasks,
    const std::vector<oligo::Permutate> &suffixPermutateList,
    const unsigned suffixLength,
    const std::size_t threadNum)
{
    boost::unique_lock<boost::mutex> lock(mutex_);
    std::size_t tasksProcessed = 0;
    while (!tasks.empty())
    {
        AnnotateTask task = tasks.back();
        tasks.pop_back();
        {
            common::unlock_guard<boost::unique_lock<boost::mutex> > unlock(lock);
            annotateBlocks(mismatchBlockMasks, suffixPermutateList, suffixLength, task.begin_, task.end_);
        }
        ++tasksProcessed;
//        if (tasks.size() && !(tasks.size() % int(pow(10.0, std::max(0, int(log10(tasks.size())))))))
        if (!(tasks.size() % 100))
        {
            ISAAC_THREAD_CERR << "Tasks left " << tasks.size() << std::endl;
        }
    }
//    ISAAC_THREAD_CERR << "processed " << tasksProcessed << " tasks on thread " << threadNum << std::endl;
}

template <typename KmerT, typename AnnotatorT>
void NeighborsFinder<KmerT, AnnotatorT>::annotateNeighbors(
    const std::vector<std::vector<KmerT> > &mismatchBlockMasks,
    const std::vector<oligo::Permutate> &suffixPermutateList,
    KmerList &kmerList)
{
    ISAAC_THREAD_CERR << " Finding neighbors" << std::endl;
    Tasks tasks;

    const long BIG_ENOUGH_BLOCK = 1000000;

    typename KmerList::iterator ourBegin = kmerList.begin();
    typename KmerList::iterator current = ourBegin;
    while (kmerList.end() != current)
    {
        typename KmerList::iterator ourEnd =
            findBlockEnd(current, oligo::KmerTraits<KmerT>::SUFFIX_BITS, kmerList.end());
        if (std::distance(ourBegin, ourEnd) > BIG_ENOUGH_BLOCK || kmerList.end() == ourEnd)
        {
            tasks.push_back(AnnotateTask(ourBegin, ourEnd));
            ourBegin = ourEnd;
        }
        current = ourEnd;
    }

    ISAAC_THREAD_CERR << "Executing " << tasks.size() << " tasks" << std::endl;
    threads_.execute(boost::bind(&NeighborsFinder::processTasks, this, boost::ref(tasks),
                                 boost::ref(mismatchBlockMasks), boost::ref(suffixPermutateList),
                                 oligo::KmerTraits<KmerT>::SUFFIX_BITS, _1));
    ISAAC_THREAD_CERR << "Executing " << tasks.size() << " tasks done" << std::endl;
    ISAAC_THREAD_CERR << " Finding neighbors done" << std::endl;
}


template <typename KmerT, typename AnnotatorT>
typename NeighborsFinder<KmerT, AnnotatorT>::KmerList::iterator NeighborsFinder<KmerT, AnnotatorT>::findBlockEnd(
    typename KmerList::iterator randomPosition,
    const unsigned suffixBits,
    const typename KmerList::iterator kmerListEnd)
{
    if (randomPosition != kmerListEnd)
    {
        const KmerT currentPrefix = (randomPosition->kmer_) >> suffixBits;
        while (randomPosition != kmerListEnd && currentPrefix == (randomPosition->kmer_) >> suffixBits)
        {
            ++randomPosition;
        }
    }
    return randomPosition;
}

template<typename KmerT, typename AnnotatorT>
std::size_t NeighborsFinder<KmerT, AnnotatorT>::annotatedSubblock(
    const typename KmerList::iterator& blockBegin,
    const typename KmerList::iterator& blockEnd) const
{
    std::size_t neighborPairs = 0;
    //            ISAAC_THREAD_CERR << "long block marking neighbors" << std::endl;
    typename KmerList::iterator subblockBegin = blockBegin;
    while (blockEnd != subblockBegin)
    {
        typename KmerList::iterator subblockEnd = findBlockEnd(
            subblockBegin, neighborhoodWidth_ * oligo::BITS_PER_BASE, blockEnd);
        for (typename KmerList::iterator i = subblockBegin; subblockEnd != i;
            ++i)
        {
            neighborPairs += markNeighbors(i, subblockEnd);
        }
        subblockBegin = subblockEnd;
    }
    return neighborPairs;
}

template <typename KmerT, typename AnnotatorT>
void NeighborsFinder<KmerT, AnnotatorT>::annotateBlock(
    const std::vector<std::vector<KmerT> > &mismatchBlockMasks,
    const std::vector<oligo::Permutate> &suffixPermutateList,
    typename KmerList::iterator blockBegin,
    typename KmerList::iterator blockEnd) const
{
    std::size_t neighborPairs = 0;
    if (1000 < std::distance(blockBegin, blockEnd))
    {
        BOOST_FOREACH(const oligo::Permutate &permutate, suffixPermutateList)
        {
            BOOST_FOREACH(AnnotatedKmer &ak, std::make_pair(blockBegin, blockEnd))
            {
                ak.kmer_ = permutateSuffix(ak.kmer_, permutate);
            }

            std::sort(blockBegin, blockEnd, &AnnotatedKmer::kmerLess);

            neighborPairs += annotatedSubblock(blockBegin, blockEnd);
        }
        // restore original order so that kmers can be used to match up with their genomic positions
        BOOST_FOREACH(AnnotatedKmer &ak, std::make_pair(blockBegin, blockEnd))
        {
            ak.kmer_ = permutateSuffix(
                ak.kmer_,
                boost::bind(&oligo::Permutate::reorder<typename oligo::KmerTraits<KmerT>::SuffixType>,
                            boost::ref(suffixPermutateList.back()), _1));
        }
    }
    else
    {
        for (typename KmerList::iterator i = blockBegin; blockEnd != i; ++i)
        {
            neighborPairs += markNeighbors(mismatchBlockMasks, i, blockEnd);
        }
    }
//    if (10000 < std::distance(blockBegin, blockEnd))
//    {
//        ISAAC_THREAD_CERR << "annotateBlock done " << std::distance(blockBegin, blockEnd) << " pairs:" << neighborPairs << std::endl;
//    }
}

template <typename KmerT, typename AnnotatorT>
void NeighborsFinder<KmerT, AnnotatorT>::annotateBlocks(
    const std::vector<std::vector<KmerT> > &mismatchBlockMasks,
    const std::vector<oligo::Permutate> &suffixPermutateList,
    const unsigned suffixLength,
    const typename KmerList::iterator ourBegin,
    const typename KmerList::iterator ourEnd) const
{
    typename KmerList::iterator blockBegin = ourBegin;
    while (ourEnd != blockBegin)
    {
        typename KmerList::iterator blockEnd = findBlockEnd(blockBegin, suffixLength, ourEnd);

        annotateBlock(mismatchBlockMasks, suffixPermutateList, blockBegin, blockEnd);
        blockBegin = blockEnd;
    }
}

static boost::array<char, 256> generateBasesDifferentArray()
{
    boost::array<char, 256> ret;
    unsigned pos = 0;
    BOOST_FOREACH(char &c, ret)
    {
        unsigned x = pos;
        while (x)
        {
            if (oligo::BITS_PER_BASE_MASK & x)
            {
                ++c;
            }
            x >>= oligo::BITS_PER_BASE;
        }
        ++pos;
    }
    return ret;
}

static const boost::array<char, 256> BASES_DIFFERENT = generateBasesDifferentArray();

/**
 * \brief Marks all kmers within neighborhoodWidth mismatches of *blockBegin.
 *        Also marks *blockBegin if it has any neighbors.
 */
template <typename KmerT, typename AnnotatorT>
std::size_t NeighborsFinder<KmerT, AnnotatorT>::markNeighbors(
    const typename KmerList::iterator blockBegin,
    const typename KmerList::iterator blockEnd) const
{
    std::size_t ret = 0;
    ISAAC_ASSERT_MSG(blockBegin <= blockEnd, "Improper range");

    if (blockEnd == blockBegin)
    {
        return ret;
    }

    typename KmerList::iterator current = blockBegin + 1;
    AnnotatedKmer &kmer = *blockBegin;
    while (blockEnd != current)
    {
        AnnotatedKmer &currentKmer = *current;
        // at least one must be f-stranded and at least one must have a position associated
        if ((!isReverseComplement(kmer) || !isReverseComplement(currentKmer)))
        {
            KmerT x = kmer.kmer_ ^ currentKmer.kmer_;
            // x has a bitmask of mismatching bases. Since it is more than 1 bit per base, we still have to
            // count now many bases mismatch. This should not take longer than one check as here the prefix
            // contains all but neighborhoodWidth last bases.
            unsigned mismatchCount = 0;
            while (!!x)
            {
                mismatchCount += BASES_DIFFERENT[x & 0xff];
                x >>= 8;
            }

            if (kmer.isReverseComplement())
            {
                ret += annotator_.update(mismatchCount, kmer, currentKmer);
            }
            else
            {
                ret += annotator_.update(mismatchCount, currentKmer, kmer);
            }
        }
        ++current;
    }
    return ret;
}


/**
 * \brief Marks all kmers within neighborhoodWidth mismatches of *blockBegin.
 *        Also marks *blockBegin if it has any neighbors.
 *
 *  This one is more complex as the mismatches can be in any bases of the suffix. mismatchBlockMasks define
 *  proper combinations of mismatches.
 */
template <typename KmerT, typename AnnotatorT>
std::size_t NeighborsFinder<KmerT, AnnotatorT>::markNeighbors(
    const std::vector<std::vector<KmerT> > &mismatchBlockMasks,
    const typename KmerList::iterator blockBegin,
    const typename KmerList::iterator blockEnd) const
{
    std::size_t ret = 0;
    ISAAC_ASSERT_MSG(blockBegin <= blockEnd, "Improper range");

    if (blockEnd == blockBegin)
    {
        return ret;
    }

    // each blockMasks represent the combination of blocks in which each block is required to have a mismatch
    // For example, when looking for distance-2 neighbors, both mismatches can be in the last block or spread over
    // the two suffix blocks.
    BOOST_FOREACH(const std::vector<KmerT> &blockMasks, mismatchBlockMasks)
    {
        KmerT allMismatchesMask(0);
        BOOST_FOREACH(const KmerT blockMask, blockMasks)
        {
            allMismatchesMask |= blockMask;
        }
        const KmerT matchMask = ~allMismatchesMask;


        typename KmerList::iterator current = blockBegin + 1;
        AnnotatedKmer &kmer = *blockBegin;
        while (blockEnd != current)
        {
            AnnotatedKmer &currentKmer = *current;
            // at least one must be f-stranded and at least one must have a position associated
            if ((!isReverseComplement(kmer) || !isReverseComplement(currentKmer)))
            {
                KmerT x = kmer.kmer_ ^ currentKmer.kmer_;
                // to avoid double counting, make sure there are no mismatches where we don't need them
                if (mismatchesInCorrectBlocks(x, matchMask, blockMasks))
                {
                    unsigned mismatchCount = 0;
                    while (!!x)
                    {
                        mismatchCount += BASES_DIFFERENT[x & 0xff];
                        x >>= 8;
                    }

                    if (kmer.isReverseComplement())
                    {
                        ret += annotator_.update(mismatchCount, kmer, currentKmer);
                    }
                    else
                    {
                        ret += annotator_.update(mismatchCount, currentKmer, kmer);
                    }
                }
            }
            ++current;
        }

    }
    return ret;
}

} // namespace reference
} //namespace isaac

#endif // #ifndef ISAAC_REFERENCE_NEIGHBORS_FINDER_HH
