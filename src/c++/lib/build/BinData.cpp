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
 ** \file BinData.cpp
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

#include "build/BinData.hh"
#include "common/Memory.hh"

namespace isaac
{
namespace build
{

unsigned BinData::getGapGroupsCount() const
{
    switch (realignGaps_)
    {
    case REALIGN_SAMPLE:
        return barcodeBamMapping_.getTotalSamples();
    case REALIGN_PROJECT:
        return barcodeBamMapping_.getMaxProjectIndex() + 1;
    case REALIGN_ALL:
        return 1;
    case REALIGN_NONE:
        return 0;
    default:
        ISAAC_ASSERT_MSG(false, "Unknown gap realignment mode: " << realignGaps_);
        break;
    }
    return -1;
}

unsigned BinData::getGapGroupIndex(const unsigned barcode) const
{
    switch (realignGaps_)
    {
    case REALIGN_SAMPLE:
        return barcodeBamMapping_.getSampleIndex(barcode);
    case REALIGN_PROJECT:
        return barcodeBamMapping_.getProjectIndex(barcode);
    case REALIGN_ALL:
        return 0;
    default:
        ISAAC_ASSERT_MSG(false, "Unknown gap realignment mode: " << realignGaps_);
        break;
    }
    return -1;
}

const gapRealigner::RealignerGaps &BinData::getRealignerGaps(const unsigned barcode) const
{
    return realignerGaps_.at(getGapGroupIndex(barcode));
}


void BinData::reserveGaps(
    const alignment::BinMetadata& bin,
    const flowcell::BarcodeMetadataList &barcodeMetadataList)
{
    std::vector<std::size_t> gapsByGroup(getGapGroupsCount(), 0);
    BOOST_FOREACH(const flowcell::BarcodeMetadata &barcode, barcodeMetadataList)
    {
        gapsByGroup.at(getGapGroupIndex(barcode.getIndex())) +=
            bin.getBarcodeGapCount(barcode.getIndex());
    }
    unsigned gapGroupId = 0;
    BOOST_FOREACH(const std::size_t gaps, gapsByGroup)
    {
        realignerGaps_.at(gapGroupId++).reserve(gaps);
    }
}

void BinData::finalize()
{
    if (REALIGN_NONE != realignGaps_)
    {
        ISAAC_THREAD_CERR << "Collecting gaps."  << std::endl;
        for(std::vector<char>::const_iterator p = data_.begin(); p != data_.end();)
        {
            const io::FragmentAccessor &fragment = *reinterpret_cast<const io::FragmentAccessor *>(&*p);

            if (!fragment.flags_.splitAlignment_ && fragment.gapCount_)
            {
                const unsigned gapGroupIndex = getGapGroupIndex(fragment.barcode_);
                realignerGaps_.at(gapGroupIndex).addGaps(fragment.fStrandPosition_, fragment.cigarBegin(), fragment.cigarEnd());
            }
            p += fragment.getTotalLength();
        }
        ISAAC_THREAD_CERR << "Finalizing gaps."  << std::endl;

        std::for_each(realignerGaps_.begin(), realignerGaps_.end(), boost::bind(&gapRealigner::RealignerGaps::finalizeGaps, _1));
    }
}


} // namespace build
} // namespace isaac
