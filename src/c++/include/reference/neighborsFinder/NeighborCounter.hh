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
 ** \file NeighborCounter.hh
 **
 ** \brief Counts the number of same-k neighbors for each genomic position.
 **
 ** \author Roman Petrovski
 **/

#ifndef ISAAC_REFERENCE_NEIGHBORS_FINDER_NEIGHBOR_COUNTER_HH
#define ISAAC_REFERENCE_NEIGHBORS_FINDER_NEIGHBOR_COUNTER_HH

#include <boost/format.hpp>

#include "common/ParallelSort.hpp"
#include "oligo/KmerGenerator.hpp"
#include "reference/Contig.hh"
#include "reference/NeighborsCount.hh"
#include "reference/ReferenceKmer.hh"
#include "reference/ReferencePosition.hh"
#include "reference/SortedReferenceMetadata.hh"

#include "reference/neighborsFinder/AnnotatedKmer.hh"
#include "reference/PermutatedKmerGenerator.hh"

namespace isaac
{
namespace reference
{
namespace neighborsFinder
{

template <typename KmerT>
class NeighborCounter: boost::noncopyable
{
public:
    NeighborCounter(
        const unsigned int maskWidth,
        const unsigned mask,
        const unsigned neighborhoodWidth,
        const reference::ContigList &contigList,
        const reference::SortedReferenceMetadata &sortedReferenceMetadata
        ) :
            neighborhoodWidth_(neighborhoodWidth),
            contigList_(contigList),
            sortedReferenceMetadata_(sortedReferenceMetadata),
            contigOffsets_(computeContigOffsets(sortedReferenceMetadata.getContigs())),
            kmerGenerator_(maskWidth, mask, contigList, sortedReferenceMetadata.getContigs()),
            // as we're not going to try all positions against all other positions, set 0 for pairs that will never match by prefix,
            annotation_(reference::genomeLength(sortedReferenceMetadata.getContigs()), 0)
    {
    }

    typedef neighborsFinder::AnnotatedKmer<KmerT, NeighborsCount> AnnotatedKmer;
    typedef std::vector<AnnotatedKmer > KmerList;
    void getKmers(const oligo::Permutate &permutate, KmerList &kmerList, common::ThreadVector &threads);

    bool update(
        const unsigned kmerMismatchCount,
        AnnotatedKmer &one,
        AnnotatedKmer &another);

    typedef std::vector<NeighborsCount> Annotation;
    template <typename ReferenceKmerT>
    void updateAnnotation(const KmerList &kmerList, const bool repeatsOnly);
    const Annotation &getAnnotation() const {return annotation_;}

private:
    const unsigned neighborhoodWidth_;
    const reference::ContigList &contigList_;
    const reference::SortedReferenceMetadata &sortedReferenceMetadata_;
    const std::vector<unsigned long> contigOffsets_;
    const PermutatedKmerGenerator<KmerT, permutatedKmerGenerator::BothStrands> kmerGenerator_;

    Annotation annotation_;

    void markRepeats(
        typename KmerList::iterator begin,
        typename KmerList::iterator end) const;

    void collapseRepeats(KmerList& kmerList);
};


template <typename KmerT>
template <typename ReferenceKmerT>
void NeighborCounter<KmerT>::updateAnnotation(
    const KmerList &kmerList, const bool repeatsOnly)
{
    if (kmerList.empty())
    {
        return;
    }
    const reference::SortedReferenceMetadata::MaskFiles &maskFileList =
        sortedReferenceMetadata_.getMaskFileList(oligo::KmerTraits<ReferenceKmerT>::KMER_BASES);
    typename KmerList::const_iterator it = kmerList.begin();
    BOOST_FOREACH(const SortedReferenceMetadata::MaskFile &maskFile, maskFileList)
    {
        ISAAC_THREAD_CERR << "Scanning " << maskFile.path << std::endl;
        if (!boost::filesystem::exists(maskFile.path))
        {
            const boost::format message = boost::format("Mask file %s does not exist: %s") % maskFile.path % strerror(errno);
            BOOST_THROW_EXCEPTION(common::IoException(ENOENT, message.str()));
        }

        std::ifstream maskInput(maskFile.path.string().c_str());
        if (!maskInput)
        {
            const boost::format message = boost::format("Failed to open mask file %s for reading: %s") % maskFile.path % strerror(errno);
            BOOST_THROW_EXCEPTION(common::IoException(errno, message.str()));
        }

        while(maskInput)
        {
            ReferenceKmer<ReferenceKmerT> referenceKmer;
            if (maskInput.read(reinterpret_cast<char *>(&referenceKmer), sizeof(referenceKmer)))
            {
//                const ReferenceKmerT maskedReferenceKmer = referenceKmer.getKmer() >> (oligo::KmerTraits<ReferenceKmerT>::KMER_BITS - oligo::KmerTraits<KmerT>::KMER_BITS);
                while (kmerList.end() != it && it->kmer_ < referenceKmer.getKmer())
                {
                    ++it;
                }

                if (kmerList.end() == it)
                {
                    break;
                }

                const AnnotatedKmer &currentNeighbor = *it;
                if (currentNeighbor.kmer_ == referenceKmer.getKmer())
                {
                    // The only reverse complements that exist at this point are those that don't have f-strand
                    // duplicates. Since reference kmers are all f-strand, reverse complements should not match to
                    // any reference kmer.
                    ISAAC_ASSERT_MSG(!currentNeighbor.isReverseComplement(), "Found match to unique reverse complement " <<
                                     currentNeighbor << " in reference kmers:" << referenceKmer);
                    const std::size_t offset = contigOffsets_.at(referenceKmer.getReferencePosition().getContigId()) + referenceKmer.getReferencePosition().getPosition();
                    if (repeatsOnly)
                    {
                        // don't count the kmer at the position itself as a repeat
                        annotation_.at(offset) += (currentNeighbor.repeats_ - 1);
                    }
                    else
                    {
                        annotation_.at(offset) += currentNeighbor.annotation_;
                    }
                }
            }
        }
        if (!maskInput.eof() && kmerList.end() != it)
        {
            const boost::format message = boost::format("Failed to update %s with neighbors information: %s") % maskFile.path % strerror(errno);
            BOOST_THROW_EXCEPTION(common::IoException(errno, message.str()));
        }
        ISAAC_THREAD_CERR << "Scanning done for " << maskFile.path << std::endl;
    }
}

template <typename KmerT>
bool NeighborCounter<KmerT>::update(
    const unsigned kmerMismatchCount,
    AnnotatedKmer &one,
    AnnotatedKmer &another)
{
    ISAAC_ASSERT_MSG(!another.isReverseComplement(), "second of the two must be a forward-facing kmer: " << one << " " << another);

    // Exact equivalence is required as we produce each k annotation separately from others
    if (neighborhoodWidth_ == kmerMismatchCount)
    {
        if (!one.isReverseComplement())
        {
            one.annotation_ += another.repeats_;
//            ISAAC_THREAD_CERR << "TADA! one:" << kmerMismatchCount << " " <<
//                AnnotatedKmer(currentPermutate->reorder(one.kmer_), one.annotation_) << " " <<
//                AnnotatedKmer(currentPermutate->reorder(another.kmer_), another.annotation_) << std::endl;
        }

        another.annotation_ += one.repeats_;
//        ISAAC_THREAD_CERR << "TADA! another:" << kmerMismatchCount << " " <<
//            AnnotatedKmer(currentPermutate->reorder(one.kmer_), one.annotation_) << " " <<
//            AnnotatedKmer(currentPermutate->reorder(another.kmer_), another.annotation_) << std::endl;
        return true;
    }

    return false;
}

template <typename KmerT>
void NeighborCounter<KmerT>::markRepeats(
    typename KmerList::iterator begin,
    typename KmerList::iterator end) const
{
    ISAAC_ASSERT_MSG(end != begin, "Empty range unexpected");
    for (typename KmerList::iterator it = begin + 1; end != it; ++it)
    {
        begin->incrementRepeats();

        // Make sure that out of two repeats, the forward-stranded one stays.
        // Otherwise the neighbor counts will not get accumulated.
        if (begin->isReverseComplement() && !it->isReverseComplement())
        {
            it->repeats_ = begin->repeats_;
            begin->setDiscard();
            begin = it;
        }
        else
        {
            it->setDiscard();
        }

    }
}

/**
 * \brief replaces stretches of kmers of high repeats with a single toomanymatch entry, updates corresponding positions
 *        int the annotation.
 */
template <typename KmerT>
void NeighborCounter<KmerT>::collapseRepeats(
    KmerList& kmerList)
{

// no need to remove uninitialized as kmer generator does not produce those anymore.
//    // remove uninitialized k-mers first
//    kmerList.erase(std::remove_if(kmerList.begin(), kmerList.end(), boost::bind(&AnnotatedKmer::isDiscard, _1)), kmerList.end());
//    ISAAC_THREAD_CERR << " remove uninitialized done for " << kmerList.size() << " forward and reverse kmers" << std::endl;

    common::parallelSort(kmerList, &AnnotatedKmer::kmerLess);
//    std::sort(kmerList.begin(), kmerList.end(), &AnnotatedKmer::kmerLess);
    ISAAC_THREAD_CERR << " reduce repeats sorting done for " << kmerList.size() << " forward and reverse kmers" << std::endl;

    for (typename KmerList::iterator it = kmerList.begin(); kmerList.end() != it;)
    {
        typename KmerList::iterator equalEnd = it + 1;
        while (kmerList.end() != equalEnd && equalEnd->kmer_ == it->kmer_)
        {
            ++equalEnd;
        }

        if (std::distance<typename KmerList::iterator>(it, equalEnd) > 1)
        {
            markRepeats(it, equalEnd);
        }
        it = equalEnd;
    }

    kmerList.erase(std::remove_if(kmerList.begin(), kmerList.end(), boost::bind(&AnnotatedKmer::isDiscard, _1)), kmerList.end());
}

template <typename KmerT>
neighborsFinder::AnnotatedKmer<KmerT, NeighborsCount> construct(const KmerT kmer, const bool reverse)
{
    neighborsFinder::AnnotatedKmer<KmerT, NeighborsCount> ret = neighborsFinder::AnnotatedKmer<KmerT, NeighborsCount>(kmer, reverse, 0);
    ret.repeats_ = 1;
    return ret;
}

template <typename KmerT>
void NeighborCounter<KmerT>::getKmers(
    const oligo::Permutate &permutate, KmerList &kmerList, common::ThreadVector &threads)
{
        kmerGenerator_.generate(permutate, 0, kmerList, boost::bind(&construct<KmerT>, _1, _4), threads);
        collapseRepeats(kmerList);
}


} // namespace neighborsFinder
} // namespace reference
} // namespace isaac

#endif // #ifndef ISAAC_REFERENCE_NEIGHBORS_FINDER_NEIGHBOR_COUNTER_HH
