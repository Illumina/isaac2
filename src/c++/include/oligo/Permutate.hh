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
 ** \file Permutate.hh
 **
 ** \brief Utility class to permutate blocks in a kmer.
 **
 ** \author Come Raczy
 **/

#ifndef iSAAC_OLIGO_PERMUTATE_HH
#define iSAAC_OLIGO_PERMUTATE_HH

#include <vector>
#include <string>

#include <boost/lambda/lambda.hpp>

#include "common/BitHacks.hh"
#include "oligo/Kmer.hh"

namespace isaac
{
namespace oligo
{

/**
 ** \brief Utility class to reorder the blocks in a kmer, from an origin permutation
 ** to a target permutation.
 **
 ** The kmer is split into 'count' (<=16) blocks, each of them at positions [0,
 ** count). The number of bases in each block is constant.
 **/
class Permutate
{
public:
    typedef std::vector<unsigned char> Order;

    Permutate(const unsigned blockLength);

    Permutate(
        unsigned blockLength,
        const std::vector<unsigned char> &from, const std::vector<unsigned char> &to,
        const unsigned mismatchMask);

    /// Apply the permutation to the kmer
    template <typename KmerT>
    KmerT operator()(const KmerT kmer) const
    {
        return transform(kmer, order_);
    }

    /// Reorder the blocks in their natural sequence (0, 1, 2...)
    template <typename KmerT>
    KmerT reorder(const KmerT kmer) const
    {
        return transform(kmer, absoluteReverseOrder_);
    }
    /// Shows a readable representation of the permutation from natural sequence
    std::string abcdToString() const;
    /// Shows a readable representation of the permutation
    std::string fromToString() const;

    const Order &getFromOrder() const {return from_;}
    const Order &getToOrder() const {return to_;}

    const std::vector<unsigned> &mismatchMasks() const
    {
        return mismatchMasks_;
    }

    void addMismatchMask(const unsigned mask)
    {
        mismatchMasks_.push_back(mask);
    }

    bool repeatsOnly() const
    {
        return 1 == mismatchMasks_.size() && 0 == mismatchMasks_.front();
    }

    static std::string orderToString(const std::vector<unsigned char> &order);

    unsigned getBlockLength() const {return blockLength_;}

    bool isSamePrefix(const unsigned prefixLength, const Permutate &that) const;
private:
   /// the length of the blocks in bases
    unsigned blockLength_;
    /// the count of blocks
    unsigned count_;
    /// the encoded order of the blocks relatively to 'from'
    Order order_;
    /// the encoded absolute order of the blocks relatively to the 'ABCD' order
    Order absoluteForwardOrder_;
    /// the encoded absolute order of the blocks relatively to the 'to' order
    Order  absoluteReverseOrder_;
    /// expected order of input data
    std::vector<unsigned char> from_;
    /// order in which the input data comes out after permutation is applied
    std::vector<unsigned char> to_;

    // bitmasks with 1 bit per block. Each entry indicating which blocks must have mismatches at the same time
    std::vector<unsigned> mismatchMasks_;

    /**
     ** \brief encode the re-ordering of the blocks from a permutation to an other.
     **
     ** Example
     **
     ** Assuming 4 blocks ABCD, numbered respectively 0, 1, 2 and 3. If the
     ** original permutation is ABDC and the targeted permutation ACBD, the
     ** 'from' vector would be [0, 1, 3, 2] and the 'to' vector [0, 2, 1,
     ** 3]. The encoded value would be 0x0231 (the block at position 0 stays at
     ** position 0, position 1 goes to position 2, position 2 goes to position 3
     ** and position 3 goes to position 1).
     **
     ** \param from the list of block number in the original permutation
     **
     ** \param to the list of block numbers in the targetted permutation
     **
     ** \return the encoded transformation from/to
     **/
    std::vector<unsigned char> encode(const std::vector<unsigned char> &from, const std::vector<unsigned char> &to);
    /// encode the re-ordering of the blocks from the permutated order into natural order (0, 1, 2, ...)
    std::vector<unsigned char> decode(const std::vector<unsigned char> &from);
    /// encode the re-ordering of the blocks from the natural order (0, 1, 2, ...)
    std::vector<unsigned char> encode(const std::vector<unsigned char> &to);
    /// apply the permutation encoded in the given order to the kmer
    template<typename KmerT>
    KmerT transform(const KmerT kmer, const std::vector<unsigned char> &order) const;

    static const unsigned ENCODING_BITS = 8;
//    static const unsigned long ENCODING_MASK = 0x0F;

    friend std::ostream &operator <<(std::ostream &os, const Permutate& permutation)
    {
        return os << "Permutate(" <<
            permutation.blockLength_ << "bl " <<
            permutation.count_ << "c " <<
            permutation.fromToString() << ")";
    }
};

std::vector<Permutate> getPermutateList(
    const unsigned blocksCount, const unsigned blockLength,
    const unsigned errorCount, const unsigned chainPrefixLength, const bool allErrorBlockCounts);


template<typename KmerT>
KmerT Permutate::transform(const KmerT kmer, const std::vector<unsigned char> &order) const
{
    ISAAC_ASSERT_MSG(blockLength_ * count_ <= oligo::KmerTraits<KmerT>::KMER_BASES, "Permutation is incompatible with Kmer type");
    ISAAC_ASSERT_MSG(count_ <= (1 << ENCODING_BITS), "Block count should fit in " << ENCODING_BITS << " bits");

    const unsigned blockBits = BITS_PER_BASE * blockLength_;
    const KmerT blockMask = ~oligo::safeShl(~KmerT(0), blockBits);
//    std::cerr << (boost::format("\n%d:%016x:%016x:%d:%016x") % blockBits % blockMask % order_ % count_ %kmer).str() << std::endl;
    KmerT ret(0);
    for (unsigned origin = 0; count_ > origin; ++origin)
    {
        const unsigned orderIndex = (count_ - origin - 1U);
        const unsigned target = order.at(order.size() - orderIndex - 1);
        const unsigned kmerOriginShift = (count_ - origin - 1U) * blockBits;
        const unsigned kmerTargetShift = (count_ - target - 1U) * blockBits;
        ret |= (((kmer >> kmerOriginShift) & blockMask) << kmerTargetShift);
//        std::cerr << (boost::format("    %016x:%d:%d:%d:%d") % ret % origin % target % kmerOriginShift % kmerTargetShift).str() << std::endl;
    }
    return ret;
}

/**
 ** \brief Generate the list of permutations to detect a given number of errors
 **
 ** \param errorCount               the number of errors to detect
 ** \param chained                  if true, the list of permutations is such that each permutation expects its input
 **                                 sequence to be a result of the previous permutation from the list.
 **                                 Otherwise, all permutations expect the input to be the original sequence
 ** \param allErrorBlockCounts      Include permutations that allow detecting errors that are spread over the less than
 **                                 errorCount blocks.
 **
 ** \return the list of Permutate components in the order where they should be
 ** applied (starting from the natural order).
 **/
template <typename KmerT>
std::vector<Permutate> getPermutateList(const unsigned errorCount, const unsigned chainPrefixLength, const bool allErrorBlockCounts)
{
    // all blocks must not overlap and have the same length. Since k of k-mers is power of 2,
    // this means that for non power of 2 errorCount we will have more permutations to go through
    // than theoretically needed
    const unsigned blocksCount = errorCount ? 2 * common::upperPowerOfTwo(errorCount) : 1;
    const unsigned blockLength = oligo::KmerTraits<KmerT>::KMER_BASES / blocksCount;
    ISAAC_ASSERT_MSG(blocksCount * blockLength == oligo::KmerTraits<KmerT>::KMER_BASES,
                     "Kmer length must be divisible by block length."
                     " oligo::KmerTraits<KmerT>::KMER_BASES=" << oligo::KmerTraits<KmerT>::KMER_BASES <<
                     " blocksCount=" << blocksCount <<
                     " blockLength=" << blockLength);
    return getPermutateList(blocksCount, blockLength, errorCount, chainPrefixLength, allErrorBlockCounts);
}

} // namespace oligo
} // namespace isaac

#endif // #ifndef iSAAC_OLIGO_PERMUTATE_HH
