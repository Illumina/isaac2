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
 ** \file MergeReferenceWorkflow.cpp
 **
 ** \brief see MergeReferenceWorkflow.hh
 **
 ** \author Roman Petrovski
 **/

#include <boost/foreach.hpp>
#include <boost/format.hpp>

#include "common/Debug.hh"
#include "common/Exceptions.hh"
#include "reference/AnnotationsMerger.hh"
#include "workflow/MergeAnnotationsWorkflow.hh"

namespace isaac
{
namespace workflow
{

MergeAnnotationsWorkflow::MergeAnnotationsWorkflow(
    const boost::filesystem::path referenceGenome,
    const std::vector<bfs::path> &filesToMerge,
    const bfs::path &outputFilePath,
    const std::vector<std::string> &aggregateFunction,
    const std::string &mergedType)
    : referenceGenome_(referenceGenome),
      filesToMerge_(filesToMerge),
      outputFilePath_(outputFilePath),
      aggregateFunction_(aggregateFunction),
      mergedType_(mergedType)
{
}


reference::DistanceToBeNeighborless minDistanceToBeNeighborless(
    const reference::DistanceToBeNeighborless merged,
    const reference::DistanceToBeNeighborless toMerge)
{
    return std::min(merged, toMerge);
}

reference::DistanceToBeNeighborless maxDistanceToBeNeighborless(
    const reference::DistanceToBeNeighborless merged,
    const reference::DistanceToBeNeighborless toMerge)
{
    return std::max(merged, toMerge);
}

reference::NeighborsCount sumNeighborsCount(
    const reference::NeighborsCount merged,
    const reference::NeighborsCount toMerge)
{
    return merged + toMerge;
}

reference::NeighborsCount minNeighborsCount(
    const reference::NeighborsCount merged,
    const reference::NeighborsCount toMerge)
{
    return std::min(merged, toMerge);
}

reference::NeighborsCount zeroAndLe10(
    const reference::NeighborsCount merged,
    const reference::NeighborsCount toMerge)
{
    reference::NeighborsCount ret = (reference::NeighborsCount(0) == merged && reference::NeighborsCount(10) >= toMerge) ?
        reference::NeighborsCount(0) : reference::NEIGHBORS_TOO_MANY;

    return ret;
}

reference::NeighborsCount zeroAndLe100(
    const reference::NeighborsCount merged,
    const reference::NeighborsCount toMerge)
{
    reference::NeighborsCount ret = (reference::NeighborsCount(0) == merged && reference::NeighborsCount(100) >= toMerge) ?
        reference::NeighborsCount(0) : reference::NEIGHBORS_TOO_MANY;

    return ret;
}

reference::DistanceToBeNeighborless mask32(
        const reference::DistanceToBeNeighborless merged,
        const reference::NeighborsCount toMerge)
{
    if (reference::NeighborsCount(0) == toMerge)
    {
        return std::min(merged, reference::DistanceToBeNeighborless(32));
    }
    return merged;
}

reference::DistanceToBeNeighborless mask(
        const reference::DistanceToBeNeighborless merged,
        const reference::NeighborsCount toMerge,
        reference::DistanceToBeNeighborless mask)
{
    if (reference::NeighborsCount(0) == toMerge)
    {
        return std::min(merged, mask);
    }
    return merged;
}

reference::DistanceToBeNeighborless mask64(
        const reference::DistanceToBeNeighborless merged,
        const reference::NeighborsCount toMerge)
{
    if (reference::NeighborsCount(0) == toMerge)
    {
        return std::min(merged, reference::DistanceToBeNeighborless(64));
    }
    return merged;
}

void MergeAnnotationsWorkflow::run()
{
    reference::AnnotationsMerger merger(
        referenceGenome_,
        outputFilePath_,
        false
        );

    if ("kul" == mergedType_)
    {
        std::vector<reference::DistanceToBeNeighborless> annotation;
        merger.start(filesToMerge_.front(), annotation);
        std::vector<std::string>::const_iterator agIt = aggregateFunction_.begin();
        for (std::vector<bfs::path>::const_iterator it = filesToMerge_.begin() + 1;
            filesToMerge_.end() != it; ++it, ++agIt)
        {
            ISAAC_THREAD_CERR << *agIt << std::endl;
            if ("min" == *agIt)
            {
                merger.merge<reference::DistanceToBeNeighborless>(
                    *it, &minDistanceToBeNeighborless, annotation);
            }
            else if ("max" == *agIt)
            {
                merger.merge<reference::DistanceToBeNeighborless>(
                    *it, &maxDistanceToBeNeighborless, annotation);
            }
            else if ("mask" == agIt->substr(0, 4))
            {
                const int value = atoi(agIt->substr(5, 2).c_str());
                ISAAC_ASSERT_MSG(value, "Unexpected mask value in " << *agIt);
                merger.merge<reference::NeighborsCount>(
                    *it, boost::bind(&mask, _1, _2, value), annotation);
            }
            else
            {
                const boost::format message = boost::format("Unexpected aggregate function %s for meged type %s") % *agIt % mergedType_;
                BOOST_THROW_EXCEPTION(common::InvalidOptionException(message.str()));
            }
        }
        merger.finish(annotation);
    }
    else if ("neighbor-counts" == mergedType_)
    {
        std::vector<reference::NeighborsCount> annotation;
        merger.start(filesToMerge_.front(), annotation);

        std::vector<std::string>::const_iterator agIt = aggregateFunction_.begin();
        for (std::vector<bfs::path>::const_iterator it = filesToMerge_.begin() + 1;
            filesToMerge_.end() != it; ++it, ++agIt)
        {
            if ("sum" == *agIt)
            {
                merger.merge<reference::NeighborsCount>(
                    *it, &sumNeighborsCount, annotation);
            }
            else if ("zero-and-leq10" == *agIt)
            {
                merger.merge<reference::NeighborsCount>(
                    *it, boost::bind(&zeroAndLe10, _1, _2), annotation);
            }
            else if ("zero-and-leq100" == *agIt)
            {
                merger.merge<reference::NeighborsCount>(
                    *it, boost::bind(&zeroAndLe100, _1, _2), annotation);
            }
            else
            {
                const boost::format message = boost::format("Unexpected aggregate function %s") % *agIt;
                BOOST_THROW_EXCEPTION(common::InvalidOptionException(message.str()));
            }
        }
        merger.finish(annotation);
    }
    else
    {
        const boost::format message = boost::format("Unexpected merged type %s") % mergedType_;
        BOOST_THROW_EXCEPTION(common::InvalidOptionException(message.str()));
    }


}

} // namespace workflow
} // namespace isaac
