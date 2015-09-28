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
 ** \file ParallelGapRealigner.cpp
 **
 ** Parallelizes cap realignment.
 ** 
 ** \author Roman Petrovski
 **/

#include <boost/foreach.hpp>

#include "build/ParallelGapRealigner.hh"

namespace isaac
{
namespace build
{

void ParallelGapRealigner::realign(
    isaac::build::GapRealigner& realigner,
    io::FragmentAccessor& fragment, PackedFragmentBuffer::Index &index,
    BinData& binData, isaac::alignment::Cigar& cigars)
{
    if (realigner.realign(binData.getRealignerGaps(fragment.barcode_), binData.bin_.getBinStart(), binData.bin_.getBinEnd(), index, fragment, binData.data_, cigars))
    {
        boost::unique_lock<boost::mutex> lock(cigarBufferMutex_);
        {
            const std::size_t before = binData.additionalCigars_.size();
            binData.additionalCigars_.addOperations(index.cigarBegin_, index.cigarEnd_);
            index.cigarBegin_ = &binData.additionalCigars_.at(before);
            index.cigarEnd_ = &binData.additionalCigars_.back() + 1;
            // realignment affects both reads. We must make sure realignment updates on one read don't
            // collide with post-realignmnet pair updates from another read.
            realigner.updatePairDetails(barcodeTemplateLengthStatistics_, index, fragment, binData.data_);
        }
    }
}

void ParallelGapRealigner::threadRealignGaps(boost::unique_lock<boost::mutex> &lock, BinData &binData, BinData::iterator &nextUnprocessed, unsigned long threadNumber)
{
//    ISAAC_THREAD_CERR << "threadRealignGaps this " << this  << std::endl;

    isaac::build::GapRealigner &realigner = threadGapRealigners_.at(threadNumber);
    isaac::alignment::Cigar &cigars = threadCigars_.at(threadNumber);
    static const std::size_t READS_AT_A_TIME = 1024;

//    int blockCount = 0;
    while (binData.indexEnd() != nextUnprocessed)
    {
        BinData::iterator ourBegin = nextUnprocessed;
        const std::size_t readsToProcess = std::min<std::size_t>(READS_AT_A_TIME, std::distance(ourBegin, binData.indexEnd()));
        nextUnprocessed += readsToProcess;
        {
            common::unlock_guard<boost::unique_lock<boost::mutex> > unlock(lock);
            for (const BinData::iterator ourEnd = ourBegin + readsToProcess; ourEnd != ourBegin; ++ourBegin)
            {
                PackedFragmentBuffer::Index &index = *ourBegin;
                io::FragmentAccessor &fragment = binData.data_.getFragment(index);
                if (binData.bin_.hasPosition(fragment.fStrandPosition_))
                {
                    cigars.clear();
                    realign(realigner, fragment, index, binData, cigars);
                }
            }
        }
//        ++blockCount;
    }

//    ISAAC_THREAD_CERR << "Thread " << threadNumber << " realigned " << blockCount << " blocks for " << binData.bin_ << std::endl;
}

} // namespace build
} // namespace isaac
