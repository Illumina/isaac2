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
 ** \file NeighborsFinderWorkflow.hh
 **
 ** \brief Glue component to find neighbors.
 **
 ** \author Come Raczy
 **/

#ifndef ISAAC_WORKFLOW_NEIGHBORS_FINDER_WORKFLOW_HH
#define ISAAC_WORKFLOW_NEIGHBORS_FINDER_WORKFLOW_HH

#include <vector>
#include <boost/filesystem.hpp>
#include <boost/noncopyable.hpp>

#include "reference/NeighborsFinder.hh"

namespace isaac
{
namespace workflow
{

class NeighborsFinderWorkflow: boost::noncopyable
{
public:

    NeighborsFinderWorkflow(
        const unsigned seedLength,
        const boost::filesystem::path &referenceGenome,
        const boost::filesystem::path &outputFile,
        const unsigned jobs);

    void run(
        const unsigned int maskWidth,
        const unsigned mask,
        const unsigned neighborhoodWidth);

    template <typename KmerT, typename ReferenceKmerT>
    void findNeighborsT(
        const unsigned int maskWidth,
        const unsigned mask,
        const unsigned neighborhoodWidth);
private:
    const unsigned seedLength_;
    const boost::filesystem::path outputFile_;
    const unsigned jobs_;
    const reference::SortedReferenceMetadata sortedReferenceMetadata_;
    common::ThreadVector threads_;
    const reference::ContigList contigList_;


    template <unsigned seedLength, unsigned referenceSeedLength>
    void resolveSeedLength(
        const unsigned int maskWidth,
        const unsigned mask,
        const unsigned neighborhoodWidth,
        const std::vector<unsigned> &permutations,
        const unsigned repeatThreshold,
        const bool estimateKul);
};

} // namespace workflow
} //namespace isaac

#endif // #ifndef ISAAC_WORKFLOW_NEIGHBORS_FINDER_WORKFLOW_HH
