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
 ** \file Permutate.cpp
 **
 ** \brief See Permutate.hh.
 **
 ** \author Come Raczy
 **/

#include <iostream>
#include <algorithm>
#include <boost/format.hpp>
#include <boost/foreach.hpp>

#include "oligo/Permutate.hh"

namespace isaac
{
namespace oligo
{

Permutate::Permutate(const unsigned blockLength):
    blockLength_(blockLength)
    , count_(1)
    , order_(1, 0)
    , absoluteForwardOrder_(1, 0)
    , absoluteReverseOrder_(1, 0)
    , from_(1, 0)
    , to_(1, 0)
{
    mismatchMasks_.push_back(0);
}

/**
 * \brief  Create a permutation from a list of positions.
 *
 * \param from          The un-permuted order is
 * \param to            target order
 *
 */
Permutate::Permutate(
    const unsigned blockLength,
    const std::vector<unsigned char> &from, const std::vector<unsigned char> &to,
    const unsigned mismatchMask)
    : blockLength_(blockLength)
    , count_(from.size())
    , order_(encode(from, to))
    , absoluteForwardOrder_(encode(to))
    , absoluteReverseOrder_(decode(to))
    , from_(from)
    , to_(to)
{
    mismatchMasks_.push_back(mismatchMask);
}

std::string Permutate::fromToString() const
{
    std::string ret = "from " + orderToString(from_) + " to " + orderToString(to_);
    return ret;
}

std::string Permutate::abcdToString() const
{
    unsigned i = 0;
    using boost::lambda::var;
    using namespace boost::assign_detail;
    const std::vector<unsigned char> origin = generic_list<unsigned char>().repeat_fun(to_.size(), var(i)++);

    std::string ret = "from " + orderToString(origin) + " to " + orderToString(to_);
    return ret;
}

std::string Permutate::orderToString(const std::vector<unsigned char> &order)
{
    static std::string S = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    std::string ret;
    BOOST_FOREACH(unsigned v, order)
    {
        ret.push_back(v < S.size() ? S[v] : '?');
    }
    return ret;
}


bool Permutate::isSamePrefix(const unsigned prefixLength, const Permutate &that) const
{
    std::pair<std::vector<unsigned char>::const_iterator, std::vector<unsigned char>::const_iterator> mismatch =
        std::mismatch(to_.begin(), to_.end(), that.to_.begin());

    if (std::distance(to_.begin(), mismatch.first) * blockLength_ * oligo::BITS_PER_BASE >= prefixLength)
    {
        return true;
    }
    return false;
}


std::vector<unsigned char> Permutate::encode(const std::vector<unsigned char> &from, const std::vector<unsigned char> &to)
{
    //std::cerr << "\nENCODE" << std::endl;
    assert(from.size() == to.size());
    std::vector<unsigned char> ret(count_ - from.size(), 0);
    std::vector<bool> checkFrom(from.size(), false);
    std::vector<bool> checkTo(to.size(), false);
    for (unsigned origin = 0; from.size() > origin; ++origin)
    {
        assert(!checkFrom[origin]);
        assert(!checkTo[origin]);
        checkFrom[origin] = true;
        checkTo[origin] = true;
        const std::vector<unsigned char>::const_iterator found = std::find(to.begin(), to.end(), from[origin]);
        ISAAC_ASSERT_MSG(to.end() != found, unsigned(from[origin]) <<
                         " not found in the 'to' permutation to: " << orderToString(to) <<
                         " from: " << orderToString(from));
        const unsigned char target = found - to.begin();
        assert(target < from.size());
        ret.push_back(target);
        //std::cerr << (boost::format("%d:%d:%d:%08x") % origin % target % from[origin] % ret).str() << std::endl;
    }

    assert(0 == std::count(checkFrom.begin(), checkFrom.end(), false));
    assert(0 == std::count(checkTo.begin(), checkTo.end(), false));
    return ret;
}

/**
 * \brief Produces order that decodes a permutated k-mer into original unpermutated state
 */
std::vector<unsigned char> Permutate::decode(const std::vector<unsigned char> &from)
{
    unsigned i = 0;
    using boost::lambda::var;
    using namespace boost::assign_detail;
    const std::vector<unsigned char> origin = generic_list<unsigned char>().repeat_fun(from.size(), var(i)++);
    return encode(from, origin);
}

/**
 * \brief Produces order that encodes an unpermutated k-mer directly into target state
 */
std::vector<unsigned char> Permutate::encode(const std::vector<unsigned char> &to)
{
    unsigned i = 0;
    using boost::lambda::var;
    using namespace boost::assign_detail;
    const std::vector<unsigned char> origin = generic_list<unsigned char>().repeat_fun(to.size(), var(i)++);
    return encode(origin, to);
}

template <typename IteratorT>
std::string vectoToString(IteratorT b, IteratorT e)
{
    std::string ret;
    BOOST_FOREACH(const unsigned ui, std::make_pair(b, e))
    {
        const std::string str = boost::lexical_cast<std::string>(ui);
        ret.insert(ret.end(), str.begin(), str.end());
        ret.push_back(' ');
    }
    return ret;
}

/**
 * \brief Produce permutations of suffix such that second half represents all choices of n elements from the stuffix.
 *
 * \param prefix Initially must be empty. Used as part of recursion
 * \sufix sequence for which the permutations must be built;
 *
 * This implementation seems to be limited in to n == suffix.size() / 2
 */
void buildPermutationList(
    const std::vector<unsigned char> &prefix,
    const std::vector<unsigned char> &suffix,
    const unsigned n,
    std::vector<std::vector<unsigned char> > &permutationList)
{
    ISAAC_THREAD_CERR << "buildPermutationList prefix: " <<
        vectoToString(prefix.begin(), prefix.end()) << " suffix: " <<
        vectoToString(suffix.begin(), suffix.end()) << std::endl;
    assert(2 * n == prefix.size() + suffix.size());
    if (prefix.size() == n)
    {
        permutationList.push_back(prefix);
        std::vector<unsigned char> &back = permutationList.back();
        back.insert(back.end(), suffix.begin(), suffix.end());

        ISAAC_THREAD_CERR << "buildPermutationList generated: " <<
            vectoToString(back.begin(), back.end()) << std::endl;
    }
    else
    {
        for (unsigned i = 0; suffix.size() > i; ++i)
        {
            if (prefix.empty() || suffix[i] > prefix.back())
            {
                std::vector<unsigned char> newPrefix = prefix;
                newPrefix.push_back(suffix[i]);
                std::vector<unsigned char> newSuffix = suffix;
                newSuffix.erase(newSuffix.begin() + i);
                buildPermutationList(newPrefix, newSuffix, n, permutationList);
            }
        }
    }
}

template <typename T>
bool orderByPrefix(const std::vector<T> &left, const std::vector<T> &right)
{
    const std::pair<
        typename std::vector<T>::const_iterator,
        typename std::vector<T>::const_iterator> mismatch = std::mismatch(left.begin(), left.end(), right.begin());
    if (left.end() == mismatch.first)
    {
        return false;
    }

    return *mismatch.first < *mismatch.second;
}

/**
 * \brief Produce permutations of suffix such that second half represents all choices of n elements from the stuffix.
 *
 * \original sequence for which the permutations must be built;
 *
 * This implementation should accept 32 > original.size() >= n
 */
template <typename T>
void buildPermutationList(
    const std::vector<T> &original,
    const unsigned n,
    std::vector<std::vector<T> > &permutationList)
{
    ISAAC_ASSERT_MSG(original.size() >= n, "Input sequence must be at least " << n << " elements long");
    typedef unsigned long BitTrackingType;
    static const BitTrackingType ALL_ONES = BitTrackingType(0) - 1;
    ISAAC_ASSERT_MSG(original.size() <= common::countBitsSet(ALL_ONES), "Input sequence must be at most " << common::countBitsSet(ALL_ONES) << " elements long");
    BitTrackingType firstWithNBitsSet = ~(ALL_ONES << n);
    static const unsigned BITS_IN_BYTE = 8;
    BitTrackingType lastWithNBitsSet = (~(ALL_ONES >> n)) >> (sizeof(BitTrackingType) * BITS_IN_BYTE - original.size());
    ISAAC_ASSERT_MSG(common::countBitsSet(firstWithNBitsSet) == n, "Expected to produce value with " << n << " bits set. Got: " << firstWithNBitsSet);

    for (BitTrackingType current = firstWithNBitsSet; current <= lastWithNBitsSet; current = common::snoob(current))
    {
        std::vector<T> permutation = original;
        ISAAC_ASSERT_MSG(!permutation.empty(), "Unexpected empty permutation");
        typename std::vector<T>::iterator lastElement = permutation.end() - 1;
        typename std::vector<T>::iterator currentElement = permutation.end() - 1;
        for (BitTrackingType u = current; 0 != u; u >>= 1)
        {
            if (1 & u)
            {
                std::swap(*lastElement, *currentElement);
                if (1 != u)
                {
                    --lastElement;
                }
            }
            if (1 != u)
            {
                --currentElement;
            }
        }
        permutationList.push_back(permutation);

//        ISAAC_THREAD_CERR << "n=" << n << " " << current << " " << Permutate::orderToString(permutation)<< std::endl;

        if (!n)
        {
            // when n == 0, there are no bits to move. Just terminated after having pushed the original permutation once.
            break;
        }
    }

    std::sort(permutationList.begin(), permutationList.end(), &orderByPrefix<T>);

//    if (n)
//    {
//        std::vector<std::vector<T> > nMinusOnePermutations;
//        buildPermutationList(original, n - 1, nMinusOnePermutations);
//        BOOST_FOREACH(const std::vector<T> &perm, nMinusOnePermutations)
//        {
//            if(permutationList.end() == std::find(permutationList.begin(), permutationList.end(), perm))
//            {
//                ISAAC_THREAD_CERR << "Lower n permutation not found n: " << n << " " << Permutate::orderToString(perm) << std::endl;
//            }
//        }
//    }

}

std::vector<Permutate> getPermutateList(
    const unsigned blocksCount, const unsigned blockLength,
    const unsigned errorCount, const unsigned chainPrefixLength, const bool allErrorBlockCounts)
{
    using boost::assign_detail::generic_list;
    using boost::lambda::var;
    unsigned char i = 0;
    const std::vector<unsigned char> suffix = generic_list<unsigned char>().repeat_fun(blocksCount, var(i)++);

    std::vector<oligo::Permutate> ret;
    // include k-0 permutation only when explicitly requested. 
    // It will not do anything useful as repeats are removed from the data before any
    // pairs get analyzed.
    std::vector<std::vector<unsigned char> > permutationList;
    buildPermutationList(suffix, errorCount, permutationList);
    const std::vector<std::vector<unsigned char> >::const_iterator origin = permutationList.begin();
    ret.push_back(Permutate(blockLength, *origin, *origin, ~(~unsigned(0) << errorCount)));

    for (std::vector<std::vector<unsigned char> >::const_iterator to = permutationList.begin() + 1; permutationList.end() != to; ++to)
    {
        Permutate chained(blockLength, *(to - 1), *to, ~(~unsigned(0) << errorCount));
        // attempt to chain permutations that have the same prefix.
        if (chained.isSamePrefix(chainPrefixLength, ret.back()))
        {
            ret.push_back(chained);
        }
        else // otherwise make the origin-based permutaiton
        {
            ret.push_back(Permutate(blockLength, *origin, *to, ~(~unsigned(0) << errorCount)));
        }
    }

    if (errorCount && allErrorBlockCounts)
    {
        unsigned ec = errorCount;
        while(--ec)
        {
            std::vector<std::vector<unsigned char> > subPermutations;
            buildPermutationList(suffix, ec, subPermutations);

            BOOST_FOREACH(std::vector<unsigned char> &subPermutation, subPermutations)
            {
                std::size_t permutationIdx = 0;
                BOOST_FOREACH(std::vector<unsigned char> &permutation, permutationList)
                {
                    unsigned mismatchMask = 0;
                    for (std::vector<unsigned char>::const_iterator it = subPermutation.end() - ec;
                        subPermutation.end() != it; ++it)
                    {
                        std::vector<unsigned char>::const_iterator pos = std::find(permutation.end() - errorCount, permutation.end(), *it);
                        if (permutation.end() == pos)
                        {
                            // permutation does not contain one of our blocks, ignore it
                            mismatchMask = 0;
                            break;
                        }
                        else
                        {
                            BOOST_FOREACH(unsigned char uc, std::make_pair(permutation.end() - errorCount, permutation.end()))
                            {
                                ISAAC_THREAD_CERR << int(uc) << ",";
                            }
                            ISAAC_THREAD_CERR << "found: " << int(*it) << " at pos " << std::distance<std::vector<unsigned char>::const_iterator>(permutation.begin(), pos) << std::endl;
                        }
                        mismatchMask |= 1 << (std::distance<std::vector<unsigned char>::const_iterator>(pos, permutation.end()) - 1);
                    }
                    if (mismatchMask)
                    {
                        ISAAC_THREAD_CERR << "add mask: " << std::hex << mismatchMask << " to " << ret.at(permutationIdx) << std::endl;
                        ret.at(permutationIdx).addMismatchMask(mismatchMask);
                        // one permutation for each subPermutation
                        break;
                    }
                    permutationIdx++;
                }
            }
        }
    }

    return ret;
}

//template<class It,class End> static void instantiateTemplates(boost::mpl::true_ endofvec){}
//template<class It,class End> static void instantiateTemplates(boost::mpl::false_);
//static void instantiateTemplates();
//
//template<class It,class End>
//void Permutate::instantiateTemplates(boost::mpl::false_)
//{
//    typedef oligo::BasicKmerType<boost::mpl::deref<It>::type::value>  KmerType;
//    const std::vector<unsigned char> order;
//
//    Permutate p(0);
//    p.transform<KmerType>(KmerType(0), order);
//
//    typedef typename boost::mpl::next<It>::type Next;
//    instantiateTemplates<Next,End>(typename boost::is_same<Next,End>::type());
//}
//
//void Permutate::instantiateTemplates()
//{
//    typedef boost::mpl::begin<oligo::SUPPORTED_KMERS>::type begin;
//    typedef boost::mpl::end<oligo::SUPPORTED_KMERS>::type end;
//
//    instantiateTemplates<begin,end>(boost::is_same<begin,end>::type());
//}

} //namespace oligo
} // namespace isaac
