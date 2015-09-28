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

inline std::size_t getTotalGapsCount(const std::vector<gapRealigner::RealignerGaps> &allGaps)
{
    return std::accumulate(allGaps.begin(), allGaps.end(), 0,
                           boost::bind(std::plus<std::size_t>(), _1,
                                       boost::bind(&gapRealigner::RealignerGaps::getGapsCount, _2)));
}


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
            barcodeTemplateLengthStatistics_(barcodeTemplateLengthStatistics),
            threadCigars_(threads),
            threadGapRealigners_(
                threads,
                GapRealigner(realignGapsVigorously, realignDodgyFragments, realignedGapsPerFragment, 3, 4, 0,
                             clipSemialigned, barcodeMetadataList, contigLists))
    {
        std::for_each(threadCigars_.begin(), threadCigars_.end(), boost::bind(&alignment::Cigar::reserve, _1, THREAD_CIGAR_MAX));
        std::for_each(threadGapRealigners_.begin(), threadGapRealigners_.end(), boost::bind(&GapRealigner::reserve, _1));
    }

    void threadRealignGaps(boost::unique_lock<boost::mutex> &lock, BinData &binData, BinData::iterator &nextUnprocessed, unsigned long threadNumber);
private:
    static const std::size_t THREAD_CIGAR_MAX =1024;
    const std::vector<alignment::TemplateLengthStatistics> &barcodeTemplateLengthStatistics_;
    std::vector<alignment::Cigar> threadCigars_;
    std::vector<GapRealigner> threadGapRealigners_;
    boost::mutex cigarBufferMutex_;

    void realign(
        isaac::build::GapRealigner& realigner,
        io::FragmentAccessor& fragment, PackedFragmentBuffer::Index &index,
        BinData& binData, isaac::alignment::Cigar& cigars);
};


} // namespace build
} // namespace isaac

#endif // #ifndef iSAAC_BUILD_PARALLEL_GAP_REALIGNER_HH
