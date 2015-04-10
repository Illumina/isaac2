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
 ** \file PermutatedKmerGenerator.hh
 **
 ** \brief Generates kmers for a given permutation that match the mask
 **
 ** \author Roman Petrovski
 **/

#ifndef ISAAC_REFERENCE_PERMUTATED_KMER_GENERATOR_HH
#define ISAAC_REFERENCE_PERMUTATED_KMER_GENERATOR_HH

#include "common/Threads.hpp"
#include "oligo/KmerGenerator.hpp"
#include "oligo/Permutate.hh"
#include "reference/Contig.hh"
#include "reference/ReferencePosition.hh"
#include "reference/SortedReferenceMetadata.hh"

namespace isaac
{
namespace reference
{
namespace permutatedKmerGenerator
{

enum Mode
{
    BothStrands, // Produce forward and reverse kmers
    ForwardOnly, // Produce forward kmers only
    Min // Produce the minimum of forward and reverse kmers. Produces both forward and reverse kmers for palindromes.
};

} // namespace permutatedKmerGenerator

template <typename KmerT, permutatedKmerGenerator::Mode mode>
class PermutatedKmerGenerator: boost::noncopyable
{
public:
    /**
     * \param contigList index-ordered list of contigs
     */
    PermutatedKmerGenerator(
        const unsigned int maskWidth,
        const unsigned mask,
        const reference::ContigList &contigList,
        const SortedReferenceMetadata::Contigs &metadataContigs) :
            maskWidth_(maskWidth),
            msbMask_(oligo::safeShl(KmerT(mask), oligo::KmerTraits<KmerT>::KMER_BITS - maskWidth)),
            msbHighlight_(~(~KmerT(0) >> maskWidth_ )),
            contigList_(contigList),
            metadataContigs_(metadataContigs),
            contigOffsets_(computeContigOffsets(metadataContigs_))
    {
        ISAAC_THREAD_CERR << "maskWidth_:" << maskWidth_ << std::endl;
        ISAAC_THREAD_CERR << "oligo::KmerTraits<KmerT>::KMER_BITS - maskWidth:" << (oligo::KmerTraits<KmerT>::KMER_BITS - maskWidth) << std::endl;
        ISAAC_THREAD_CERR << "mask:" << mask << std::endl;
        ISAAC_THREAD_CERR << "msbMask_:" << oligo::bases(msbMask_) << std::endl;
    }

    template<typename ConstrucT, typename KmerLisT, typename AnnotationT>
    void generate(
        const oligo::Permutate &permutate,
        AnnotationT initV,
        KmerLisT &kmerList,
        const ConstrucT &construct,
        common::ThreadVector &threads) const;

private:
    const unsigned int maskWidth_;
    const KmerT msbMask_;
    const KmerT msbHighlight_;
    /// index-ordered list of contigs
    const reference::ContigList &contigList_;
    const SortedReferenceMetadata::Contigs &metadataContigs_;
    const std::vector<unsigned long> contigOffsets_;

    bool matchesMask(
        const KmerT &k) const;

    template<typename KmerLisT, typename ConstrucT>
    void generateThread(
        const oligo::Permutate &permutate,
        KmerLisT &kmerList,
        const ConstrucT &construct,
        const unsigned threadNumber,
        const std::size_t threads,
        std::vector<std::size_t> &threadKmerOffsets) const;

    template<typename KmerLisT>
    void countThread(
        const oligo::Permutate &permutate,
        KmerLisT &kmerList,
        const unsigned threadNumber,
        const std::size_t threads,
        std::vector<std::size_t> &threadKmerCounts) const;

    template <typename KmerLisT, typename CallbackT>
    void kmerThread(
        const oligo::Permutate &permutate,
        KmerLisT &kmerList,
        const unsigned threadNumber,
        const std::size_t threads,
        const CallbackT &callback) const;

    template <typename CallbackT>
    void generate(
        const KmerT &kmer,
        const oligo::Permutate &permutate,
        const CallbackT &callback) const;

    template <typename KmerLisT, typename ConstrucT>
    static void storeThreadKmer(
        const unsigned thread,
        const KmerT &kmer,
        const int contigId,
        const unsigned long kmerPosition,
        const bool reverse,
        KmerLisT &kmerList,
        const ConstrucT &construct,
        std::vector<std::size_t> &threadKmerOffsets);

};

template <typename KmerT, permutatedKmerGenerator::Mode mode>
bool PermutatedKmerGenerator<KmerT, mode>::matchesMask(
    const KmerT &k) const
{
    return (k & msbHighlight_) == msbMask_;
}

template <typename KmerT, permutatedKmerGenerator::Mode mode>
template <typename CallbackT>
void PermutatedKmerGenerator<KmerT, mode>::generate(
    const KmerT &kmer,
    const oligo::Permutate &permutate,
    const CallbackT &callback) const
{
    using namespace permutatedKmerGenerator;
    const KmerT k = permutate(kmer);
    if (ForwardOnly == mode)
    {
        if (matchesMask(k))
        {
            callback(k, false);
        }
        return;
    }

    const KmerT reverseK = permutate(oligo::reverseComplement(kmer));
    if (BothStrands == mode)
    {
        if (matchesMask(k))
        {
            callback(k, false);
        }
        if (matchesMask(reverseK))
        {
            callback(reverseK, true);
        }
        return;
    }

    if (Min == mode)
    {
        if(k < reverseK)
        {
            if (matchesMask(k))
            {
                callback(k, false);
            }
        }
        else if(reverseK < k)
        {
            if (matchesMask(reverseK))
            {
                callback(reverseK, true);
            }
        }
        else // palindrome, register twice
        {
            if (matchesMask(k))
            {
                callback(k, false);
                callback(k, true);
            }
        }
    }
}


template <typename KmerT, permutatedKmerGenerator::Mode mode>
template <typename KmerLisT, typename CallbackT>
void PermutatedKmerGenerator<KmerT, mode>::kmerThread(
    const oligo::Permutate &permutate,
    KmerLisT &kmerList,
    const unsigned threadNumber,
    const std::size_t threads,
    const CallbackT &callback) const
{
//    ISAAC_THREAD_CERR << "genomeLength:" << genomeLength << std::endl;
//    ISAAC_THREAD_CERR << "threadSectionLength:" << threadSectionLength << std::endl;
//    ISAAC_THREAD_CERR << "beginGenomicOffset:" << beginGenomicOffset << std::endl;
//    ISAAC_THREAD_CERR << "endGenomicOffset:" << endGenomicOffset << std::endl;

    std::size_t generated = 0;
    BOOST_FOREACH(const reference::Contig &contig, contigList_)
    {
        ISAAC_ASSERT_MSG(metadataContigs_.at(contig.index_).totalBases_ == contig.getLength(),
                         "Discrepancy between metadata and data sizes:" << contig << " " << metadataContigs_.at(contig.index_))

        const std::size_t threadSectionLength = (contig.getLength() + threads - 1) / threads;

        const std::size_t beginGenomicOffset = threadSectionLength * threadNumber;
        if (contig.getLength() > beginGenomicOffset)
        {
            ISAAC_THREAD_CERR << "Generating kmers for " << contig << std::endl;

            const std::size_t endGenomicOffset =
                (contig.getLength() - beginGenomicOffset < threadSectionLength + oligo::KmerTraits<KmerT>::KMER_BASES - 1) ?
                    contig.getLength() : (beginGenomicOffset + threadSectionLength + oligo::KmerTraits<KmerT>::KMER_BASES - 1);

            oligo::KmerGenerator<oligo::KmerTraits<KmerT>::KMER_BASES, KmerT, std::vector<char>::const_iterator> kmerGenerator(
                contig.forward_.begin() + beginGenomicOffset,
                contig.forward_.begin() + endGenomicOffset);

    //            ISAAC_THREAD_CERR << " contig:" << contig << std::endl;
    //            ISAAC_THREAD_CERR << "  genomicOffset:" << "  generatorOffset:" << generatorOffset << "  generatorEnd:" << generatorEnd << std::endl;

            KmerT kmer(0);
            for (std::vector<char>::const_iterator it = kmerGenerator.next(kmer);
                contig.forward_.begin() + endGenomicOffset != it;
                it = kmerGenerator.next(kmer))
            {
    //                ISAAC_ASSERT_MSG(genomicOffset < endGenomicOffset,
    //                                 "Passed the end of expected range genomicOffset=" << genomicOffset << " endGenomicOffset=" << endGenomicOffset);
    //                genomicOffset += std::distance(prevIt, it);
    //                prevIt = it;
                const unsigned long kmerPosition = std::distance(contig.forward_.begin(), it);
                generate(kmer, permutate, boost::bind(callback, threadNumber, _1, contig.index_, kmerPosition, _2));
                ++generated;
            }
        }
    }

//    ISAAC_THREAD_CERR << " thread" << threadNumber << " generated " << generated << " kmers" << std::endl;
}

inline void incrementThreadKmerCounts(const unsigned thread, std::vector<std::size_t> &threadKmerCounts)
{
    ++threadKmerCounts.at(thread);
}

template <typename KmerT, permutatedKmerGenerator::Mode mode>
template <typename KmerLisT, typename ConstrucT>
void PermutatedKmerGenerator<KmerT, mode>::storeThreadKmer(
    const unsigned thread,
    const KmerT &kmer,
    const int contigId,
    const unsigned long kmerPosition,
    const bool reverse,
    KmerLisT &kmerList,
    const ConstrucT &construct,
    std::vector<std::size_t> &threadKmerOffsets)
{
    const std::size_t threadOffset = threadKmerOffsets.at(thread)++;
    kmerList.at(threadOffset) = construct(kmer, contigId, kmerPosition, reverse);
}

template <typename KmerT, permutatedKmerGenerator::Mode mode>
template <typename KmerLisT>
void PermutatedKmerGenerator<KmerT, mode>::countThread(
    const oligo::Permutate &permutate,
    KmerLisT &kmerList,
    const unsigned threadNumber,
    const std::size_t threads,
    std::vector<std::size_t> &threadKmerCounts) const
{
    kmerThread(
        permutate, kmerList, threadNumber, threads,
        boost::bind(&incrementThreadKmerCounts, _1, boost::ref(threadKmerCounts)));
}

template <typename KmerT, permutatedKmerGenerator::Mode mode>
template<typename KmerLisT, typename ConstrucT>
void PermutatedKmerGenerator<KmerT, mode>::generateThread(
    const oligo::Permutate &permutate,
    KmerLisT &kmerList,
    const ConstrucT &construct,
    const unsigned threadNumber,
    const std::size_t threads,
    std::vector<std::size_t> &threadKmerOffsets) const
{
    kmerThread(
        permutate, kmerList, threadNumber, threads,
        boost::bind(&PermutatedKmerGenerator::storeThreadKmer<KmerLisT, ConstrucT>, _1, _2, _3, _4, _5,
                    boost::ref(kmerList),
                    boost::ref(construct),
                    boost::ref(threadKmerOffsets)));
}


template <typename KmerT, permutatedKmerGenerator::Mode mode>
template<typename ConstrucT, typename KmerLisT, typename AnnotationT>
void PermutatedKmerGenerator<KmerT, mode>::generate(
    const oligo::Permutate &permutate,
    AnnotationT initV,
    KmerLisT &kmerList,
    const ConstrucT &construct,
    common::ThreadVector &threads) const
{
    std::vector<std::size_t> threadKmerCounts(threads.size(), 0);

    threads.execute(boost::bind(
        &PermutatedKmerGenerator::countThread<KmerLisT>,
        this,
        boost::ref(permutate),
        boost::ref(kmerList),
        _1, threads.size(),
        boost::ref(threadKmerCounts)));

    const std::size_t estimatedKmers = std::accumulate(threadKmerCounts.begin(), threadKmerCounts.end(), 0UL);

    ISAAC_THREAD_CERR << " found " << estimatedKmers << " kmers on " << threads.size() << " threads" << std::endl;
    // if we don't deallocate, resize will attempt to allocate memory separately and then copy the old contents (which we don't need)
    // in some cases this will obviously cause bad_alloc despite the total amount of ram needed being available.
    {KmerLisT().swap(kmerList);}
    kmerList.resize(estimatedKmers, typename KmerLisT::value_type(initV));
    ISAAC_THREAD_CERR << " reserving memory done for " << kmerList.size() << " kmers" << std::endl;

    unsigned long offset = 0;
    BOOST_FOREACH(std::size_t &count, threadKmerCounts)
    {
        using std::swap; swap(offset, count);
        offset += count;
    }
    threads.execute(boost::bind(
        &PermutatedKmerGenerator::generateThread<KmerLisT, ConstrucT>,
        this,
        boost::ref(permutate),
        boost::ref(kmerList),
        boost::ref(construct),
        _1, threads.size(),
        boost::ref(threadKmerCounts)));

    ISAAC_THREAD_CERR << " generated " << estimatedKmers << " kmers" << std::endl;
}

} // namespace reference
} // namespace isaac

#endif // #ifndef ISAAC_REFERENCE_PERMUTATED_KMER_GENERATOR_HH
