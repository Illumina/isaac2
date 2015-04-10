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

inline std::size_t getTotalGapsCount(const std::vector<gapRealigner::RealignerGaps> &allGaps)
{
    return std::accumulate(allGaps.begin(), allGaps.end(), 0,
                           boost::bind(std::plus<std::size_t>(), _1,
                                       boost::bind(&gapRealigner::RealignerGaps::getGapsCount, _2)));
}

void ParallelGapRealigner::threadRealignGaps(BinData &binData, unsigned long threadNumber)
{
    BOOST_FOREACH(PackedFragmentBuffer::Index &index, std::make_pair(binData.indexBegin(), binData.indexEnd()))
    {
        io::FragmentAccessor &fragment = binData.data_.getFragment(index);
        // realignment affects both reads. We must make sure each read of the same cluster is processed on the same thread
        // otherwise realignment updates on one read can collide with post-realignmnet pair updates from another read.
        if (fragment.clusterId_ % threads_.size() == threadNumber)
        {
            if (binData.bin_.hasPosition(fragment.fStrandPosition_))
            {
                if (threadGapRealigners_.at(threadNumber).realign(
                    binData.getRealignerGaps(fragment.barcode_),
                    binData.bin_.getBinStart(),
                    binData.bin_.getBinEnd(), index, fragment,
                    binData.data_,
                    threadCigars_.at(threadNumber)))
                {
                    boost::unique_lock<boost::mutex> lock(cigarBufferMutex_);
                    const std::size_t before = binData.additionalCigars_.size();
                    binData.additionalCigars_.insert(binData.additionalCigars_.end(), index.cigarBegin_, index.cigarEnd_);
                    index.cigarBegin_ = &binData.additionalCigars_.at(before);
                    index.cigarEnd_ = &binData.additionalCigars_.back() + 1;
                    threadCigars_.at(threadNumber).clear();
                }
            }
        }
    }
}

void ParallelGapRealigner::realignGaps(BinData &binData)
{
    ISAAC_THREAD_CERR << "Realigning against " << getTotalGapsCount(binData.realignerGaps_) <<
        " unique gaps. " << binData.bin_ << std::endl;

    threads_.execute(boost::bind(&ParallelGapRealigner::threadRealignGaps, this, boost::ref(binData), _1));

    ISAAC_THREAD_CERR << "Realigning gaps done" << std::endl;
}


} // namespace build
} // namespace isaac
