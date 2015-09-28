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
 ** \file BinSorter.cpp
 **
 ** Reorders aligned data and stores results in bam file.
 ** 
 ** \author Roman Petrovski
 **/

#include <boost/foreach.hpp>
#include <boost/function_output_iterator.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/ptr_container/ptr_vector.hpp>

#include "bam/Bam.hh"

#include "build/BinSorter.hh"
#include "common/Memory.hh"

namespace isaac
{
namespace build
{

unsigned long BinSorter::serialize(
    BinData &binData,
    boost::ptr_vector<boost::iostreams::filtering_ostream> &bgzfStreams,
    boost::ptr_vector<bam::BamIndexPart> &bamIndexParts)
{
    if (!binData.getUniqueRecordsCount())
    {
        return 0;
    }
    ISAAC_THREAD_CERR << "Sorting offsets for bam " << binData.bin_ << std::endl;

    bamSerializer_.prepareForBam(binData.data_, binData, binData.additionalCigars_);

    ISAAC_THREAD_CERR << "Sorting offsets for bam done " << binData.bin_ << std::endl;

    ISAAC_THREAD_CERR << "Serializing records: " << binData.getUniqueRecordsCount() <<  " of them for bin " << binData.bin_ << std::endl;

    std::time_t serTimeStart = common::time();

    if (binData.isUnalignedBin())
    {
        unsigned long offset = 0;
        while(binData.data_.size() != offset)
        {
            const io::FragmentAccessor &fragment = binData.data_.getFragment(offset);
            bamSerializer_.storeUnaligned(fragment, bgzfStreams, bamIndexParts, binData.bamAdapter_(fragment));
            offset += fragment.getTotalLength();
        }
    }
    else
    {
        BOOST_FOREACH(const PackedFragmentBuffer::Index& idx, binData)
        {
            // realigning reads that don't belong to the bin is not very useful
            // also, it can move the read position and cause more than one copy of the
            // read to be stored in the bam file.
            if (binData.bin_.hasPosition(idx.pos_))
            {
                downgradeAlignmentScores(idx, binData.data_);
                const io::FragmentAccessor &fragment = binData.data_.getFragment(idx);
                bamSerializer_.storeAligned(fragment, bgzfStreams, bamIndexParts, binData.bamAdapter_(idx, fragment));
            }
            //else the fragment got split into a bit that does not belong to the current bin. it will get stored by another bin BinSorter.
        }
    }

    BOOST_FOREACH(boost::iostreams::filtering_ostream &bgzfStream, bgzfStreams)
    {
        ISAAC_ASSERT_MSG(bgzfStream.strict_sync(), "Expecting the compressor to flush all the data");
    }

    std::time_t serTimeEnd = common::time();

    ISAAC_THREAD_CERR << "Serializing records done: " << binData.getUniqueRecordsCount() <<  " of them for bin " << binData.bin_ << " in " << ::std::difftime(serTimeEnd, serTimeStart) << "seconds." << std::endl;
    return binData.size();
}

void BinSorter::resolveDuplicates(
    BinData &binData,
    BuildStats &buildStats)
{
    ISAAC_THREAD_CERR << "Resolving duplicates for bin " << binData.bin_ << std::endl;

    NotAFilter().filterInput(binData.data_, binData.seIdx_.begin(), binData.seIdx_.end(), buildStats, binData.binStatsIndex_, std::back_inserter(binData));
    if (keepDuplicates_ && !markDuplicates_)
    {
        NotAFilter().filterInput(binData.data_, binData.rIdx_.begin(), binData.rIdx_.end(), buildStats, binData.binStatsIndex_, std::back_inserter(binData));
        NotAFilter().filterInput(binData.data_, binData.fIdx_.begin(), binData.fIdx_.end(), buildStats, binData.binStatsIndex_, std::back_inserter(binData));
    }
    else
    {
        if (singleLibrarySamples_)
        {
            DuplicatePairEndFilter(keepDuplicates_).filterInput(
                RSDuplicateFilter<true>(binData.barcodeBamMapping_.getSampleIndexMap()),
                binData.data_, binData.rIdx_.begin(), binData.rIdx_.end(),
                buildStats, binData.binStatsIndex_, std::back_inserter(binData));
            DuplicatePairEndFilter(keepDuplicates_).filterInput(
                FDuplicateFilter<true>(binData.barcodeBamMapping_.getSampleIndexMap()),
                binData.data_, binData.fIdx_.begin(), binData.fIdx_.end(),
                buildStats, binData.binStatsIndex_, std::back_inserter(binData));
        }
        else
        {
            DuplicatePairEndFilter(keepDuplicates_).filterInput(
                RSDuplicateFilter<false>(binData.barcodeBamMapping_.getSampleIndexMap()),
                binData.data_, binData.rIdx_.begin(), binData.rIdx_.end(),
                buildStats, binData.binStatsIndex_, std::back_inserter(binData));
            DuplicatePairEndFilter(keepDuplicates_).filterInput(
                FDuplicateFilter<false>(binData.barcodeBamMapping_.getSampleIndexMap()),
                binData.data_, binData.fIdx_.begin(), binData.fIdx_.end(),
                buildStats, binData.binStatsIndex_, std::back_inserter(binData));
        }
    }

    // we will not be needing these anymore. Free up some memory so that other bins get a chance to start earlier
    binData.unreserveIndexes();

    ISAAC_THREAD_CERR << "Resolving duplicates done for bin " << binData.bin_ << std::endl;
}

template <typename BclIteratorT>
bool BinSorter::isAnchored(
    const isaac::reference::ContigList &reference,
    const isaac::reference::ContigAnnotation &contigAnnotation,
    const isaac::reference::ReferencePosition &fpos,
    BclIteratorT currentBase,
    const uint32_t *cigarBegin,
    const uint32_t *cigarEnd) const
{
    std::vector<char>::const_iterator currentReference = reference.at(fpos.getContigId()).forward_.begin() + fpos.getPosition();
    isaac::reference::Annotation::const_iterator currentAnnotation = contigAnnotation.begin() + fpos.getPosition();

    for (const uint32_t *cigarIterator = cigarBegin; cigarIterator != cigarEnd; ++cigarIterator)
    {
        const alignment::Cigar::Component cigar = alignment::Cigar::decode(*cigarIterator);
        const unsigned length = cigar.first;
        const alignment::Cigar::OpCode opCode = cigar.second;
        if (opCode == alignment::Cigar::ALIGN)
        {
            currentReference += length;
            currentAnnotation += length;
            currentBase += length;

            unsigned matchesInARow = 0;
            for (unsigned i = 0; i < length; ++i)
            {
                --currentReference;
                --currentAnnotation;
                --currentBase;

                if (!alignment::isMatch(oligo::getSequenceBaseFromBcl(*currentBase), *currentReference))
                {
                    matchesInARow = 0;
                }
                else
                {
                    ++matchesInARow;
                    if (isaac::reference::K_UNIQUE_TOO_FAR != *currentAnnotation &&
                        matchesInARow >= *currentAnnotation)
                    {
                        return true;
                    }
                }
            }

            currentReference += length;
            currentAnnotation += length;
            currentBase += length;
        }
        else if (opCode == alignment::Cigar::INSERT)
        {
            currentBase += length;
        }
        else if (opCode == alignment::Cigar::DELETE)
        {
            currentReference += length;
            currentAnnotation += length;
        }
        else if (opCode == alignment::Cigar::SOFT_CLIP)
        {
            // NOTE! Not advancing the reference for soft clips
            currentBase += length;
        }
        else
        {
            BOOST_THROW_EXCEPTION(common::PostConditionException(
                (boost::format("Only canonical CIGAR codes are allowed. Unexpected Cigar OpCode: %d") % opCode).str()));
        }

    }
    return false;
}

void BinSorter::downgradeAlignmentScores(const PackedFragmentBuffer::Index& idx, PackedFragmentBuffer &data)
{
    io::FragmentAccessor &fragment = data.getFragment(idx);
    // only rescan the CIGAR for realigned reads. Otherwise, trust the MatchSelector to do the decent job. This
    // preserves seed-anchored read alignment scores.
    if (fragment.isAligned() && fragment.flags_.realigned_)
    {
        const isaac::reference::ContigAnnotations &contigAnnotations =
            annotations_.at(barcodeMetadataList_.at(fragment.barcode_).getReferenceIndex());
        if (!contigAnnotations.empty())
        {
            const isaac::reference::ContigList &reference = contigLists_.at(barcodeMetadataList_.at(fragment.barcode_).getReferenceIndex());
            const isaac::reference::ContigAnnotation &contigAnnotation = contigAnnotations.at(idx.pos_.getContigId());
            bool anchored = idx.reverse_ == fragment.isReverse() ?
                isAnchored(reference, contigAnnotation, idx.pos_, fragment.basesBegin(), idx.cigarBegin_, idx.cigarEnd_) :
                isAnchored(reference, contigAnnotation, idx.pos_,
                           boost::make_transform_iterator(boost::make_reverse_iterator(fragment.basesEnd()), &oligo::getReverseBcl), idx.cigarBegin_, idx.cigarEnd_);

            if (!anchored && anchorMate_)
            {
                // if we're dealing with a sensible pair and mate is in the same bin
                if (fragment.flags_.properPair_)
                {
                    ISAAC_ASSERT_MSG(idx.hasMate(), "All paired reads expected to have mate copy stored in the same bin");
                    const io::FragmentAccessor &mate = data.getMate(idx);
                    if (!mate.flags_.realigned_)
                    {
                        anchored = isAnchored(reference, contigAnnotation,
                                              mate.getFStrandReferencePosition(), mate.basesBegin(), mate.cigarBegin(), mate.cigarEnd());
                    }
                }
            }

            if (!anchored)
            {
                ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment.clusterId_, "Downgrading " << fragment << " " << idx);
                fragment.flags_.dodgy_ = true;
            }
            else
            {
                fragment.flags_.dodgy_ = false;
            }
        }
    }
}

} // namespace build
} // namespace isaac
