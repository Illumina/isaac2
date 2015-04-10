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
 ** \file ReferenceSorter.cpp
 **
 ** Top level component to produce a sorted reference.
 **
 ** \author Come Raczy
 **/

#include <boost/assert.hpp>
#include <boost/foreach.hpp>
#include <boost/io/ios_state.hpp>

#include "common/Exceptions.hh"
#include "common/ParallelSort.hpp"
#include "common/SystemCompatibility.hh"
#include "io/BitsetLoader.hh"
#include "io/BitsetSaver.hh"
#include "io/FastaReader.hh"
#include "oligo/Nucleotides.hh"
#include "oligo/Mask.hh"
#include "reference/ContigLoader.hh"
#include "reference/ReferencePosition.hh"
#include "reference/ReferenceSorter.hh"
#include "reference/SortedReferenceXml.hh"

namespace isaac
{
namespace reference
{

template <typename KmerT>
ReferenceSorter<KmerT>::ReferenceSorter (
    const unsigned int maskWidth,
    const unsigned mask,
    const boost::filesystem::path &contigsXmlPath,
    const boost::filesystem::path &genomeNeighborsFile,
    const boost::filesystem::path &outputFile,
    const unsigned repeatThreshold
    )
    : repeatThreshold_(repeatThreshold)
    , maskWidth_(maskWidth)
    , mask_(mask)
    , msbMask_(~((~(KmerT(0)) >> maskWidth)))
    // avoid shifting by number of bits in the type. It does not do anything on Intel and some of the custom wigth types trigger assertion failure.
    , maskBits_(mask_ ? (KmerT(mask_) << (oligo::KmerTraits<KmerT>::KMER_BASES * oligo::BITS_PER_BASE - maskWidth_)) : KmerT(mask_))
    , contigsXmlPath_(contigsXmlPath)
    , genomeNeighborsFile_(genomeNeighborsFile)
    , outputFile_(outputFile)
    , sortedReferenceMetadata_(reference::loadSortedReferenceXml(contigsXmlPath_))
    , threads_(boost::thread::hardware_concurrency())
    , contigList_(reference::loadContigs(sortedReferenceMetadata_.getContigs(), threads_))
{
    boost::io::ios_flags_saver svr(std::cerr);
    ISAAC_THREAD_CERR <<
            "Constructing ReferenceSorter: for " << oligo::KmerTraits<KmerT>::KMER_BASES << "-mers " <<
            " mask width: " << maskWidth_ <<
            " msbMask_: " << msbMask_ <<
            " maskBits_: " << maskBits_ <<
            " genomeFile_: " << contigsXmlPath_ <<
            " outputFile_: " << outputFile_ <<
            std::endl;

    BOOST_ASSERT(mask_ < isaac::oligo::getMaskCount(maskWidth_) && "Mask value cannot exceed the allowed bit width");
}


template <typename KmerT>
void ReferenceSorter<KmerT>::run()
{
    std::vector<unsigned long> contigOffsets(computeContigOffsets(sortedReferenceMetadata_.getContigs()));

    unsigned long genomeLength = reference::genomeLength(sortedReferenceMetadata_.getContigs());
    {
        loadReference();
        ISAAC_THREAD_CERR << "Loaded genome from " << contigsXmlPath_ << " found " << genomeLength << " bases" << std::endl;
    }

    sortReference();

    markRepeats();

    // remove marked repeats and reverse complements
    reference_.erase(std::remove_if(reference_.begin(), reference_.end(), boost::bind(&ReferenceKmer<KmerT>::hasNeighbors, _1)), reference_.end());

    std::vector<bool> neighbors;
    if (!genomeNeighborsFile_.empty())
    {
        io::BitsetLoader loader(genomeNeighborsFile_);
        const unsigned long neighborsCount = loader.load(genomeLength, neighbors);
        ISAAC_THREAD_CERR << "Scanning " << genomeNeighborsFile_ << " found " << neighborsCount << " neighbors among " << genomeLength << " bases" << std::endl;

        BOOST_FOREACH(ReferenceKmer<KmerT> &referenceKmer, reference_)
        {
            const bool kmerHasNeighbors = referenceKmer.isTooManyMatch() ? false :
                neighbors.at(contigOffsets.at(referenceKmer.getReferencePosition().getContigId()) +
                             referenceKmer.getReferencePosition().getPosition());
            if (kmerHasNeighbors)
            {
                referenceKmer.setNeighbors(true);
            }
        }
    }
    saveReference();
}

template <typename KmerT>
ReferenceKmer<KmerT> construct(
    const KmerT kmer, const int contigId, const unsigned long kmerPosition, const bool reverse)
{
    return ReferenceKmer<KmerT>(kmer, ReferencePosition(contigId, kmerPosition, false, reverse));
}

template <typename KmerT>
void ReferenceSorter<KmerT>::loadReference()
{
    ISAAC_THREAD_CERR << "Loading " << oligo::KmerTraits<KmerT>::KMER_BASES << "-mers" << std::endl;

    const clock_t start = clock();
    if (!repeatThreshold_)
    {
        const PermutatedKmerGenerator<KmerT, permutatedKmerGenerator::ForwardOnly>
            neighborPositionGenerator(maskWidth_, mask_, contigList_, sortedReferenceMetadata_.getContigs());
        neighborPositionGenerator.generate(oligo::Permutate(oligo::KmerTraits<KmerT>::KMER_BASES), KmerT(0), reference_,
                                           boost::bind(&construct<KmerT>, _1, _2, _3, _4), threads_);
    }
    else
    {
        const PermutatedKmerGenerator<KmerT, permutatedKmerGenerator::Min>
            referenceKmerGenerator(maskWidth_, mask_, contigList_, sortedReferenceMetadata_.getContigs());
        referenceKmerGenerator.generate(oligo::Permutate(oligo::KmerTraits<KmerT>::KMER_BASES), KmerT(0), reference_,
                                        boost::bind(&construct<KmerT>, _1, _2, _3, _4), threads_);
    }


    ISAAC_THREAD_CERR << "Loading " << oligo::KmerTraits<KmerT>::KMER_BASES << "-mers" << " done in " << (clock() - start) / 1000 << "ms" << std::endl;
}

template <typename KmerT>
void ReferenceSorter<KmerT>::sortReference()
{
    ISAAC_THREAD_CERR << "Sorting " << reference_.size() << " " << oligo::KmerTraits<KmerT>::KMER_BASES << "-mers" << std::endl;
    const clock_t start = clock();
    common::parallelSort(reference_, &compareKmerThenPosition<KmerT>);
    ISAAC_THREAD_CERR << "Sorting " << reference_.size() << " " << oligo::KmerTraits<KmerT>::KMER_BASES << "-mers" << " done in " << (clock() - start) / 1000 << "ms" << std::endl;
}

template <typename KmerT>
void ReferenceSorter<KmerT>::markRepeats()
{
    if (!repeatThreshold_)
    {
        return;
    }

    typename std::vector<ReferenceKmer<KmerT> >::iterator current(reference_.begin());

    while(reference_.end() != current)
    {
        std::pair<typename std::vector<ReferenceKmer<KmerT> >::iterator,
                  typename std::vector<ReferenceKmer<KmerT> >::iterator> sameKmerRange =
            std::equal_range(current, reference_.end(), *current, &compareKmer<KmerT>);
        const std::size_t kmerMatches = std::distance(sameKmerRange.first, sameKmerRange.second);

        if (repeatThreshold_ < kmerMatches)
        {
            ISAAC_THREAD_CERR << "Skipping kmer " << oligo::bases(current->getKmer()) << " as it generates " << kmerMatches << "matches\n";

            static const ReferencePosition tooManyMatchPosition(ReferencePosition::TooManyMatch);
            *sameKmerRange.first = ReferenceKmer<KmerT>(sameKmerRange.first->getKmer(), tooManyMatchPosition);
            std::for_each(sameKmerRange.first + 1, sameKmerRange.second, boost::bind(&ReferenceKmer<KmerT>::setNeighbors, _1, true));
        }
        current = sameKmerRange.second;
    }
}

template <typename KmerT>
void ReferenceSorter<KmerT>::saveReference()
{
    ISAAC_THREAD_CERR << "Saving " << reference_.size() << " " << oligo::KmerTraits<KmerT>::KMER_BASES << "-mers" << std::endl;
    const clock_t start = clock();

    std::ofstream os(outputFile_.c_str());
    if (!os)
    {
        BOOST_THROW_EXCEPTION(common::IoException(errno,"Failed to create file " + outputFile_.string()));
    }

    if (!reference_.empty())
    {
        if (!os.write(reinterpret_cast<const char*>(&reference_.front()), sizeof(reference_.front()) * reference_.size()))
        {
            BOOST_THROW_EXCEPTION(common::IoException(errno,"Failed to write reference kmer into " + outputFile_.string()));
        }
    }

    os.flush();
    os.close();
    ISAAC_THREAD_CERR << "Saving " << reference_.size() << " " << oligo::KmerTraits<KmerT>::KMER_BASES << "-mers done in " << (clock() - start) / 1000 << "ms" << std::endl;

    SortedReferenceMetadata sortedReference;
    sortedReference.addMaskFile(oligo::KmerTraits<KmerT>::KMER_BASES, maskWidth_, mask_, outputFile_, reference_.size());
    saveSortedReferenceXml(std::cout, sortedReference);
}

template class ReferenceSorter<oligo::VeryShortKmerType>;
template class ReferenceSorter<oligo::BasicKmerType<12> >;
//template class ReferenceSorter<oligo::BasicKmerType<14> >;
template class ReferenceSorter<oligo::ShortKmerType>;
//template class ReferenceSorter<oligo::BasicKmerType<18> >;
template class ReferenceSorter<oligo::BasicKmerType<20> >;
//template class ReferenceSorter<oligo::BasicKmerType<22> >;
template class ReferenceSorter<oligo::BasicKmerType<24> >;
//template class ReferenceSorter<oligo::BasicKmerType<26> >;
template class ReferenceSorter<oligo::BasicKmerType<28> >;
//template class ReferenceSorter<oligo::BasicKmerType<30> >;
template class ReferenceSorter<oligo::KmerType>;
//template class ReferenceSorter<oligo::BasicKmerType<34> >;
template class ReferenceSorter<oligo::BasicKmerType<36> >;
//template class ReferenceSorter<oligo::BasicKmerType<38> >;
template class ReferenceSorter<oligo::BasicKmerType<40> >;
//template class ReferenceSorter<oligo::BasicKmerType<42> >;
template class ReferenceSorter<oligo::BasicKmerType<44> >;
//template class ReferenceSorter<oligo::BasicKmerType<46> >;
template class ReferenceSorter<oligo::BasicKmerType<48> >;
//template class ReferenceSorter<oligo::BasicKmerType<50> >;
template class ReferenceSorter<oligo::BasicKmerType<52> >;
//template class ReferenceSorter<oligo::BasicKmerType<54> >;
template class ReferenceSorter<oligo::BasicKmerType<56> >;
//template class ReferenceSorter<oligo::BasicKmerType<58> >;
template class ReferenceSorter<oligo::BasicKmerType<60> >;
//template class ReferenceSorter<oligo::BasicKmerType<62> >;
template class ReferenceSorter<oligo::LongKmerType>;
//template class ReferenceSorter<oligo::BasicKmerType<66> >;
template class ReferenceSorter<oligo::BasicKmerType<68> >;
//template class ReferenceSorter<oligo::BasicKmerType<70> >;
template class ReferenceSorter<oligo::BasicKmerType<72> >;
//template class ReferenceSorter<oligo::BasicKmerType<74> >;
template class ReferenceSorter<oligo::BasicKmerType<76> >;
//template class ReferenceSorter<oligo::BasicKmerType<78> >;
template class ReferenceSorter<oligo::BasicKmerType<80> >;

} // namespace reference
} // namespace isaac
