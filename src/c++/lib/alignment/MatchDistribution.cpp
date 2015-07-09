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
 ** \file MatchDistribution.cpp
 **
 ** Component to keep track of the distribution of the matches across the
 ** reference genome.
 **
 ** \author Come Raczy
 **/

#include <functional>
#include <boost/foreach.hpp>

#include "alignment/MatchDistribution.hh"
#include "common/Debug.hh"

namespace isaac
{
namespace alignment
{

MatchDistribution::MatchDistribution(const isaac::reference::SortedReferenceMetadataList &sortedReferenceMetadataList) :
    binSize_(std::max(1UL, isaac::reference::getLongestGenomeLength(sortedReferenceMetadataList) / BINS_PER_GENOME))
{
    ISAAC_THREAD_CERR << "MatchDistribution binSize_=" << binSize_ << std::endl;
    BOOST_FOREACH(const isaac::reference::SortedReferenceMetadata &sortedReferenceMetadata, sortedReferenceMetadataList)
    {
        const std::vector<isaac::reference::SortedReferenceMetadata::Contig> xmlContigs = sortedReferenceMetadata.getContigs();
        resize(std::max(size(), xmlContigs.size()));
        BOOST_FOREACH(const isaac::reference::SortedReferenceMetadata::Contig &xmlContig, xmlContigs)
        {
            const unsigned long binCount = (xmlContig.totalBases_ + getBinSize() - 1) / getBinSize();
            at(xmlContig.karyotypeIndex_).resize(std::max(binCount, at(xmlContig.karyotypeIndex_).size()), 0);
        }
    }
}

/**
 * \brief Sum up contents of two MatchDistribution objects
 */
void MatchDistribution::consolidate(const MatchDistribution &matchDistribution)
{
    ISAAC_ASSERT_MSG(size() == matchDistribution.size(), "MatchDistribution geometries must match for consolidation");
    for (size_t contig = 0; size() > contig; ++contig)
    {
        std::vector<unsigned> &ourContig = (*this)[contig];
        ISAAC_ASSERT_MSG(ourContig.size() == matchDistribution[contig].size(), "MatchDistribution geometries must match for consolidation");
        const std::vector<unsigned> &theirConting = matchDistribution[contig];
        std::transform(theirConting.begin(), theirConting.end(), ourContig.begin(),
                       ourContig.begin(), std::plus<unsigned>());
    }
}

void MatchDistribution::consolidate(const std::vector<MatchDistribution> &matchDistributionList)
{
    BOOST_FOREACH(const MatchDistribution &matchDistribution, matchDistributionList)
    {
        consolidate(matchDistribution);
    }
}


} // namespace alignemnt
} // namespace isaac
