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
 ** \file AnnotatedKmer.hh
 **
 ** \brief Combination of kmer and reference position with special meaning assigned to a number of
 **        reference position magic states
 **
 ** \author Roman Petrovski
 **/

#ifndef ISAAC_REFERENCE_NEIGHBORS_FINDER_ANNOTATED_KMER_HH
#define ISAAC_REFERENCE_NEIGHBORS_FINDER_ANNOTATED_KMER_HH

#include <boost/utility/value_init.hpp>

#include "common/ParallelSort.hpp"
#include "oligo/KmerGenerator.hpp"
#include "reference/Contig.hh"
#include "reference/KUniqueness.hh"
#include "reference/ReferencePosition.hh"
#include "reference/SortedReferenceMetadata.hh"

namespace isaac
{
namespace reference
{
namespace neighborsFinder
{


template <typename KmerT, typename AnnotationT>
struct AnnotatedKmer
{
    typedef KmerT KmerType;
    KmerT kmer_;
    struct
    {
        bool reverse_ : 1;
        /// number of repeats including this one. When 0, indicates a repeat that needs to be discarded during repeat reduction
        unsigned char repeats_ : 7;
    };
    static const unsigned char REPEATS_MAX = 0x7F;

    AnnotationT annotation_;

    // constructs discarded AnnotatedKmer
    AnnotatedKmer(
        const AnnotationT initV) :
            kmer_(KmerT(0)), reverse_(false), repeats_(0), annotation_(initV)
    {
    }

    // constructs non-discarded AnnotatedKmer
    AnnotatedKmer(
        const KmerT kmer,
        const AnnotationT initV) :
            kmer_(kmer), reverse_(false), repeats_(1), annotation_(initV)
    {
    }

    // constructs non-discarded AnnotatedKmer
    AnnotatedKmer(
        const KmerT kmer,
        bool reverseComplement,
        const AnnotationT initV) :
            kmer_(kmer), reverse_(reverseComplement), repeats_(1), annotation_(initV)
    {
    }
    bool isReverseComplement() const
        {return reverse_;}
    void setReverseComplement(const bool rc)
        {reverse_ = rc;}
    void setReverseComplement()
        {reverse_ = true;}
    bool isDiscard() const
        {return !repeats_;}
    unsigned repeats() const {return repeats_;}
    void setDiscard()
        {repeats_ = 0;}
    void incrementRepeats()
    {
        repeats_ += (REPEATS_MAX != repeats_);
    }

    friend std::ostream &operator <<(std::ostream &os, const AnnotatedKmer &kmer)
    {
        os << "AnnotatedKmer(" << oligo::Bases<2, KmerT>(kmer.kmer_, oligo::KmerTraits<KmerT>::KMER_BASES) << " " << int(kmer.repeats_) << "r " << kmer.annotation_ << "a ";
        if (kmer.isReverseComplement())
        {
            return os << "rc)";
        }
        else
        {
            return os << ")";
        }
    }

    static inline bool kmerLess(
        const AnnotatedKmer &lhs,
        const AnnotatedKmer &rhs)
    {
        return lhs.kmer_ < rhs.kmer_;
    }
}__attribute__ ((packed)); // this significantly reduces memory requirement for finder especially considering the fact
                           // that parallel sort needs twice the memory for processing.
} // namespace neighborsFinder
} // namespace reference
} // namespace isaac

#endif // #ifndef ISAAC_REFERENCE_NEIGHBORS_FINDER_ANNOTATED_KMER_HH
