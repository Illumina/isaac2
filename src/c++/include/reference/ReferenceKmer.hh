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
 ** \file ReferenceKmer.hh
 **
 ** Representation of a k-mer at a given position in a reference genome.
 **
 ** \author Come Raczy
 **/

#ifndef iSAAC_REFERENCE_REFERENCE_KMER_HH
#define iSAAC_REFERENCE_REFERENCE_KMER_HH

#include <utility>

#include <boost/mpl/back.hpp>
#include <boost/mpl/begin_end.hpp>
#include <boost/mpl/copy_if.hpp>
#include <boost/mpl/deref.hpp>
#include <boost/mpl/front.hpp>
#include <boost/mpl/modulus.hpp>

#include "oligo/Kmer.hh"
#include "oligo/Permutate.hh"
#include "reference/ReferencePosition.hh"

namespace isaac
{
namespace reference
{

typedef boost::mpl::copy_if<
    oligo::SUPPORTED_KMERS
    , boost::mpl::equal_to<boost::mpl::modulus<boost::mpl::_1, boost::mpl::int_<4> >, boost::mpl::int_<0> >
    , boost::mpl::back_inserter< boost::mpl::vector<> >
    >::type SUPPORTED_KMERS;

static const unsigned FIRST_SUPPORTED_KMER = boost::mpl::front<SUPPORTED_KMERS>::type::value;
static const unsigned LAST_SUPPORTED_KMER = boost::mpl::back<SUPPORTED_KMERS>::type::value;

template<class It,class End>
bool isSupportedKmerLength(int seedLength,boost::mpl::true_ endofvec)
{
    return false;
}

template<class It,class End>
bool isSupportedKmerLength(int seedLength, boost::mpl::false_)
{
    if(seedLength == boost::mpl::deref<It>::type::value)
    {
        return true;
    }
    else
    {
        typedef typename boost::mpl::next<It>::type Next;
        return isSupportedKmerLength<Next,End>(seedLength, typename boost::is_same<Next,End>::type());
    }
}

inline bool isSupportedKmerLength(int length)
{
    typedef boost::mpl::begin<SUPPORTED_KMERS>::type begin;
    typedef boost::mpl::end<SUPPORTED_KMERS>::type end;

    return isSupportedKmerLength<begin,end>(length, boost::is_same<begin,end>::type());
}

template <typename KmerT>
class ReferenceKmer
{
public:
    KmerT first;
    ReferencePosition::value_type second;
    /// Encodes a contig Id and a position into a ReferencePosition
    ReferenceKmer(
        const KmerT &kmer = KmerT(0),
        const ReferencePosition &referencePosition = ReferencePosition(0))
        : first(kmer), second(referencePosition.getValue())
    {
    }
    KmerT getKmer() const {return first;}
    const ReferencePosition getTranslatedPosition(
        const std::vector<unsigned> &contigTranslationTable) const {return ReferencePosition(second).translateContig(contigTranslationTable);}
    const ReferencePosition getReferencePosition() const {return ReferencePosition(second);}
    void setKmer(const KmerT kmer) {first = kmer;}
    void setNeighbors(bool set){second = ReferencePosition(second).setNeighbors(set).getValue();}
    void setTooManyMatch(){second = ReferencePosition(ReferencePosition::TooManyMatch).getValue();}
    bool isTooManyMatch() const {return second == ReferencePosition(ReferencePosition::TooManyMatch).getValue();}
    bool hasNeighbors() const {return ReferencePosition(second).hasNeighbors();}
    bool hasNoNeighbors() const {return !ReferencePosition(second).hasNeighbors();}
    bool reverse() const {return ReferencePosition(second).reverse();}
} __attribute__ ((packed));

BOOST_STATIC_ASSERT(sizeof(ReferenceKmer<oligo::KmerType>) == (sizeof(oligo::KmerType) + sizeof(ReferencePosition)));
BOOST_STATIC_ASSERT(sizeof(ReferenceKmer<oligo::LongKmerType>) == (sizeof(oligo::LongKmerType) + sizeof(ReferencePosition)));

template <typename KmerT>
inline std::ostream &operator<<(std::ostream &os, const ReferenceKmer<KmerT> &rk)
{
    return os << "ReferenceKmer(" << rk.getKmer() << ":" << isaac::oligo::bases(rk.getKmer()) << "," << rk.getReferencePosition() << ")";
}


template <typename KmerT>
inline bool compareKmer(const ReferenceKmer<KmerT> &lhs, const ReferenceKmer<KmerT> &rhs)
{
    return lhs.getKmer() < rhs.getKmer();
}

template <typename KmerT>
inline bool comparePosition(const ReferenceKmer<KmerT> &lhs, const ReferenceKmer<KmerT> &rhs)
{
    return lhs.getReferencePosition() < rhs.getReferencePosition();
}

template <typename KmerT>
inline bool compareKmerThenPosition(const ReferenceKmer<KmerT> &lhs, const ReferenceKmer<KmerT> &rhs)
{
    return lhs.getKmer() < rhs.getKmer() ? true : lhs.getKmer() != rhs.getKmer() ? false : comparePosition(lhs, rhs);
}

} // namespace reference
} // namespace isaac

#endif // #ifndef iSAAC_REFERENCE_REFERENCE_KMER_HH
