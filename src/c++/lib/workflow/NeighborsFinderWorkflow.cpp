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
 ** \file NeighborsFinderWorkflow.cpp
 **
 ** \brief See NeighborsFinderWorkflow.hh.
 **
 ** \author Roman Petrovski
 **/

#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include "reference/AnnotationLoader.hh"
#include "reference/SortedReferenceXml.hh"
#include "reference/neighborsFinder/NeighborCounter.hh"
#include "workflow/NeighborsFinderWorkflow.hh"

namespace isaac
{
namespace workflow
{

namespace bfs = boost::filesystem;

NeighborsFinderWorkflow::NeighborsFinderWorkflow(
    const unsigned seedLength,
    const bfs::path &referenceGenome,
    const bfs::path &outputFile,
    const unsigned jobs)
    : seedLength_(seedLength)
    , outputFile_(outputFile)
    , jobs_(jobs)
    , sortedReferenceMetadata_(reference::loadSortedReferenceXml(referenceGenome))
    , threads_(jobs)
    , contigList_(reference::loadContigs(sortedReferenceMetadata_.getContigs(), threads_))
{
}

template <typename KmerT, typename ReferenceKmerT>
void NeighborsFinderWorkflow::findNeighborsT(
    const unsigned maskWidth,
    const unsigned mask,
    const unsigned neighborhoodWidth)
{
    typedef isaac::reference::NeighborsFinder<KmerT, reference::neighborsFinder::NeighborCounter<KmerT> >  NeighborsFinderT;
    reference::neighborsFinder::NeighborCounter<KmerT> neighborsCounter(
        maskWidth, mask, neighborhoodWidth, contigList_, sortedReferenceMetadata_);
    NeighborsFinderT neighborsFinder(
        neighborsCounter,
        maskWidth,
        neighborhoodWidth,
        threads_);
    neighborsFinder.template annotate<ReferenceKmerT>();
    const std::vector<reference::NeighborsCount> &annotation = neighborsCounter.getAnnotation();
    boost::iostreams::filtering_ostream os;
    os.push(boost::iostreams::gzip_compressor(boost::iostreams::gzip_params(1)));
    os.push(boost::iostreams::file_sink(outputFile_.string()));
    if (!os.write(reinterpret_cast<const char*>(&annotation.front()),
                  annotation.size() * sizeof(reference::NeighborsCount)))
    {
        const boost::format message = boost::format("Failed to write neighbor counts into %s: %s") % outputFile_.string() % strerror(errno);
        BOOST_THROW_EXCEPTION(common::IoException(errno, message.str()));
    }
}

template <typename It,typename End, bool endTrueFalse>
struct NeighborsFinderWorkflowSeedLengthResolver
{
    NeighborsFinderWorkflowSeedLengthResolver(
        const int seedLength,
        const unsigned int maskWidth,
        const unsigned mask,
        const unsigned neighborhoodWidth,
        const reference::SortedReferenceMetadata &sortedReferenceMetadata,
        NeighborsFinderWorkflow &workflow)
    {
        if (boost::mpl::deref<It>::type::value == seedLength)
        {
            ISAAC_ASSERT_MSG(sortedReferenceMetadata.supportsSeedLength(seedLength), "Sorted reference does not support seedLength " << seedLength);
            workflow.findNeighborsT<isaac::oligo::BasicKmerType<boost::mpl::deref<It>::type::value>, isaac::oligo::BasicKmerType<boost::mpl::deref<It>::type::value> >(
                maskWidth, mask, neighborhoodWidth);
        }
        else
        {
            typedef typename boost::mpl::next<It>::type Next;
            NeighborsFinderWorkflowSeedLengthResolver<Next, End, boost::is_same<Next,End>::type::value>(
                seedLength, maskWidth, mask, neighborhoodWidth, sortedReferenceMetadata, workflow);
        }
    }
};

template <typename It,typename End>
struct NeighborsFinderWorkflowSeedLengthResolver<It, End, true>
{
    NeighborsFinderWorkflowSeedLengthResolver(
        const int seedLength,
        const unsigned int maskWidth,
        const unsigned mask,
        const unsigned neighborhoodWidth,
        const reference::SortedReferenceMetadata &sortedReferenceMetadata,
        NeighborsFinderWorkflow &workflow)
    {
        ISAAC_ASSERT_MSG(false, "Unexpected seedLength requested: " << seedLength);
    }
};

void NeighborsFinderWorkflow::run(
    const unsigned int maskWidth,
    const unsigned mask,
    const unsigned neighborhoodWidth)
{

    // Skip first supported kmer. It does not have suffix type defined
    typedef boost::mpl::next<boost::mpl::begin<reference::SUPPORTED_KMERS>::type>::type Begin;
    typedef boost::mpl::end<reference::SUPPORTED_KMERS>::type End;

    NeighborsFinderWorkflowSeedLengthResolver<Begin, End, boost::is_same<Begin,End>::type::value>(
        seedLength_, maskWidth, mask, neighborhoodWidth, sortedReferenceMetadata_, *this);
}

} // namespace workflow
} //namespace isaac
