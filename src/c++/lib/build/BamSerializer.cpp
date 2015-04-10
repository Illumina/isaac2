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
 ** \file BamSerializer.cpp
 **
 ** Reorders aligned data and stores results in bam file.
 ** 
 ** \author Roman Petrovski
 **/

#include "build/BamSerializer.hh"

namespace isaac
{
namespace build
{

/**
 * \brief split into multiple bam segments if it cannot be represented as single one.
 */
void BamSerializer::splitIfNeeded(
    PackedFragmentBuffer &data,
    PackedFragmentBuffer::Index &index,
    std::vector<PackedFragmentBuffer::Index> &splitIndexEntries,
    alignment::Cigar &splitCigars)
{
    // update all index record positions before sorting as they may have been messed up by gap realignment
    io::FragmentAccessor &fragment = data.getFragment(index);
    index.pos_ = fragment.fStrandPosition_;
    alignment::CigarPosition<PackedFragmentBuffer::Index::CigarIterator> last(index.cigarBegin_, index.cigarEnd_, index.pos_, fragment.isReverse(), fragment.readLength_);
    for (alignment::CigarPosition<PackedFragmentBuffer::Index::CigarIterator> current = last;
        !current.end(); ++current)
    {
        if (current.referencePos_.getContigId() != last.referencePos_.getContigId() ||
            current.referencePos_ < last.referencePos_ ||
            current.reverse_ != last.reverse_ ||
            ((current.sequenceOffset_ - last.sequenceOffset_) != (current.referencePos_ - last.referencePos_) &&
                (current.referencePos_ - last.referencePos_) > splitGapLength_))
        {
            ISAAC_ASSERT_MSG(splitIndexEntries.size() < splitIndexEntries.capacity(), "New entries must not cause buffer reallocation");
            ISAAC_ASSERT_MSG(splitCigars.size() + std::distance(index.cigarBegin_, index.cigarEnd_) * 2 <= splitCigars.capacity(),
                "New entries must not cause cigar buffer reallocation");

            // create new entry
            PackedFragmentBuffer::Index secondPart(index);

            const std::size_t before = splitCigars.size();
            if (current.reverse_ == last.reverse_)
            {
                // soft clip sequence that belongs to the previous part of the split
                splitCigars.addOperation(last.sequenceOffset_, alignment::Cigar::SOFT_CLIP);
            }

            secondPart.pos_ = current.referencePos_;
            std::copy(current.cigarIt_, index.cigarEnd_, std::back_inserter(splitCigars));
            secondPart.cigarBegin_ = &splitCigars.front() + before;
            secondPart.cigarEnd_ = &splitCigars.back() + 1;
            secondPart.reverse_ = current.reverse_;
            splitIndexEntries.push_back(secondPart);

            //patch the old one
            const PackedFragmentBuffer::Index::CigarIterator oldBegin = index.cigarBegin_;
            index.cigarBegin_ = &splitCigars.back() + 1;
            std::copy(oldBegin, last.cigarIt_, std::back_inserter(splitCigars));

            const unsigned sequenceLeftover = data.getFragment(index).readLength_ - last.sequenceOffset_;
            if (sequenceLeftover)
            {
                // soft clip sequence that belongs to the next part of the split
                splitCigars.addOperation(sequenceLeftover, alignment::Cigar::SOFT_CLIP);
            }
            index.cigarEnd_ = &splitCigars.back() + 1;

            //change iterator to travel over the CIGAR of the new entry
            current = alignment::CigarPosition<PackedFragmentBuffer::Index::CigarIterator>(secondPart.cigarBegin_, secondPart.cigarEnd_, secondPart.pos_, current.reverse_, fragment.readLength_);

            fragment.flags_.properPair_ = false;
            io::FragmentAccessor &mate = data.getMate(index);
            mate.flags_.properPair_ = false;

            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment.clusterId_, "Split into " << index << " and " << secondPart);
        }
        last = current;
    }
}

void BamSerializer::prepareForBam(
    PackedFragmentBuffer &data,
    std::vector<PackedFragmentBuffer::Index> &dataIndex,
    alignment::Cigar &splitCigars)
{
    ISAAC_THREAD_CERR << "Sorting offsets for bam" << std::endl;
    BOOST_FOREACH(PackedFragmentBuffer::Index &index, dataIndex)
    {
        splitIfNeeded(data, index, dataIndex, splitCigars);
    }

    const clock_t startSortOffsets = clock();
    std::sort(dataIndex.begin(), dataIndex.end(), boost::bind(&PackedFragmentBuffer::orderForBam, boost::ref(data), _1, _2));
    ISAAC_THREAD_CERR << "Sorting offsets for bam" << " done in " << (clock() - startSortOffsets) / 1000 << "ms" << std::endl;
}

} // namespace build
} // namespace isaac
