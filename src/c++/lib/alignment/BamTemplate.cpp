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
 ** \file BamTemplate.cpp
 **
 ** \brief See BamTemplate.hh
 ** 
 ** \author Come Raczy
 **/

#include <boost/foreach.hpp>

#include "alignment/BamTemplate.hh"

namespace isaac
{
namespace alignment
{

BamTemplate::BamTemplate()
    : alignmentScore_(-1U)
    , properPair_(false)
{
}

BamTemplate::BamTemplate(
    const FragmentMetadata &read1,
    const FragmentMetadata &read2)
    : alignmentScore_(-1U)
    , properPair_(false)
{
    fragmentMetadataList_.push_back(read1);
    fragmentMetadataList_.push_back(read2);
}


void BamTemplate::initialize(const flowcell::ReadMetadataList &tileReads, const Cluster &cluster)
{
    fragmentMetadataList_.clear();
    alignmentScore_ = -1U;
    properPair_ = false;
    BOOST_FOREACH(const flowcell::ReadMetadata &read, tileReads)
    {
        fragmentMetadataList_.push_back(FragmentMetadata(&cluster, read.getIndex()));
    }
}

bool BamTemplate::filterLowQualityFragments(const unsigned mapqThreshold)
{
    bool ret = false;
    unsigned alignmentScore = 0;
    for (unsigned i = 0; getFragmentCount() > i; ++i)
    {
        FragmentMetadata &fragment = getFragmentMetadata(i);
        if (mapqThreshold > fragment.getAlignmentScore())
        {
            fragment.cigarLength = 0;
            fragment.cigarOffset = 0;
            fragment.alignmentScore = 0;
            const FragmentMetadata &mate = getFragmentMetadata((i+1) % getFragmentCount());
            fragment.position = mate.position;
            fragment.contigId = mate.contigId;
        }
        else if (fragment.isAligned())
        {
            ret = true;
        }
        // update the alignment score of the templates that didn't resolve into a proper pair
        alignmentScore += fragment.alignmentScore;
    }
    setAlignmentScore(alignmentScore);
    return ret;
}

} // namespace alignment
} // namespace isaac
