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
 ** \file BinSorter.hh
 **
 ** Performs sorting and duplicate marking on a single alignment bin.
 ** 
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_BUILD_PARALLEL_GAP_REALIGNER_HH
#define iSAAC_BUILD_PARALLEL_GAP_REALIGNER_HH

#include "alignment/TemplateLengthStatistics.hh"
#include "build/BinData.hh"
#include "build/GapRealigner.hh"
#include "common/Threads.hpp"

namespace isaac
{
namespace build
{

class ParallelGapRealigner
{
public:
    ParallelGapRealigner(
        const unsigned threads,
        const bool realignGapsVigorously,
        const bool realignDodgyFragments,
        const unsigned realignedGapsPerFragment,
        const bool clipSemialigned,
        const flowcell::BarcodeMetadataList &barcodeMetadataList,
        const std::vector<alignment::TemplateLengthStatistics> &barcodeTemplateLengthStatistics,
        const isaac::reference::ContigLists &contigLists) :
            threads_(threads),
            threadCigars_(threads_.size()),
            threadGapRealigners_(
                threads_.size(),
                GapRealigner(realignGapsVigorously, realignDodgyFragments, realignedGapsPerFragment, 3, 4, 0,
                             clipSemialigned, barcodeMetadataList, barcodeTemplateLengthStatistics, contigLists))
    {
        std::for_each(threadCigars_.begin(), threadCigars_.end(), boost::bind(&alignment::Cigar::reserve, _1, THREAD_CIGAR_MAX));
        std::for_each(threadGapRealigners_.begin(), threadGapRealigners_.end(), boost::bind(&GapRealigner::reserve, _1));
    }

    void realignGaps(BinData &binData);
private:
    static const std::size_t THREAD_CIGAR_MAX =1024;
    common::ThreadVector threads_;
    std::vector<alignment::Cigar> threadCigars_;
    std::vector<GapRealigner> threadGapRealigners_;
    boost::mutex cigarBufferMutex_;

    void threadRealignGaps(BinData &binData, unsigned long threadNumber);

};


} // namespace build
} // namespace isaac

#endif // #ifndef iSAAC_BUILD_PARALLEL_GAP_REALIGNER_HH
