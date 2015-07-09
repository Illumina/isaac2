/**
 ** Isaac Genome Alignment Software
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
 ** \file AnnotationsMerger.cpp
 **
 ** \brief see AnnotationsMerger.hh
 **
 ** \author Roman Petrovski
 **/

#include <boost/foreach.hpp>
#include <boost/format.hpp>

#include "common/Debug.hh"
#include "common/Exceptions.hh"
#include "reference/AnnotationsMerger.hh"

namespace isaac
{
namespace reference
{
AnnotationsMerger::AnnotationsMerger(
    const bfs::path &referenceGenome,
    const bfs::path &outputFilePath,
    const bool outputReadLengths)
    : outputFilePath_(outputFilePath),
      outputReadLengths_(outputReadLengths),
      sortedReferenceMetadata_(reference::loadSortedReferenceXml(referenceGenome)),
      contigOffsets_(reference::computeContigOffsets(sortedReferenceMetadata_.getContigs()))
{
}

/**
 * \brief some positions with TOO_MANY or TOO_FAR annotations can be rescued by their adjacent positions that have
 *        sensible length stored
 */
template <>
void AnnotationsMerger::postProcess(std::vector<reference::DistanceToBeNeighborless> &annotation) const
{
    unsigned contigId = 0;
    BOOST_FOREACH(const reference::SortedReferenceMetadata::Contig &contig, sortedReferenceMetadata_.getContigs())
    {
        std::vector<reference::DistanceToBeNeighborless>::iterator contigBegin = annotation.begin() + contigOffsets_.at(contigId);
        std::vector<reference::DistanceToBeNeighborless>::iterator contigEnd = annotation.begin() + contigOffsets_.at(contigId) + contig.totalBases_;

        typedef std::reverse_iterator<std::vector<reference::DistanceToBeNeighborless>::iterator> RIter;
        RIter currentRIt(contigEnd);
        const RIter contigREnd(contigBegin);
        unsigned lastNeighborlessLength = reference::K_UNIQUE_TOO_FAR;
        while (currentRIt != contigREnd)
        {
            if (reference::K_UNIQUE_TOO_FAR != *currentRIt)
            {
                if (lastNeighborlessLength && *currentRIt > lastNeighborlessLength)
                {
                    *currentRIt = lastNeighborlessLength;
                }
                else
                {
                    lastNeighborlessLength = *currentRIt;
                }
            }
            else
            {
                if (lastNeighborlessLength < reference::K_UNIQUE_TOO_FAR)
                {
                    *currentRIt = lastNeighborlessLength;
                }
            }

            ++currentRIt;
            if (reference::K_UNIQUE_TOO_FAR != lastNeighborlessLength)
            {
                ++lastNeighborlessLength;
            }
        }

        ++contigId;
    }

}

template <>
void AnnotationsMerger::initAnnotation(const std::string &initFunction, std::vector<reference::DistanceToBeNeighborless> &annotation)
{
    if ("init-kul-65535" == initFunction)
    {
        ISAAC_THREAD_CERR << "initializing: " << initFunction << std::endl;
        annotation.clear();
        annotation.resize(reference::genomeLength(sortedReferenceMetadata_.getContigs()), reference::DistanceToBeNeighborless(65535));
    }
    else
    {
        const boost::format message = boost::format("Unexpected annotation init function %s for DistanceToBeNeighborless annotation") % initFunction;
        BOOST_THROW_EXCEPTION(common::InvalidOptionException(message.str()));
    }
}

template <>
void AnnotationsMerger::initAnnotation(const std::string &initFunction, std::vector<reference::NeighborsCount> &annotation)
{
    const boost::format message = boost::format("Unexpected annotation init function %s for NeighborsCount annotation") % initFunction;
    BOOST_THROW_EXCEPTION(common::InvalidOptionException(message.str()));
}

template <>
void AnnotationsMerger::postProcess(std::vector<reference::NeighborsCount> &annotation) const
{

}

} // namespace reference
} // namespace isaac
