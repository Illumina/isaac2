/**
 ** Isaac Genome Alignment Software
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
 ** \file AnnotationsMerger.hh
 **
 ** \brief merge multiple annotation data files into one
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_REFERENCE_ANNOTATIONS_MERGER_HH
#define iSAAC_REFERENCE_ANNOTATIONS_MERGER_HH

#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>

#include "common/FileSystem.hh"
#include "reference/KUniqueness.hh"
#include "reference/NeighborsCount.hh"
#include "reference/SortedReferenceXml.hh"

namespace isaac
{
namespace reference
{

namespace bfs = boost::filesystem;

class AnnotationsMerger: boost::noncopyable
{
private:
    const bfs::path &outputFilePath_;
    const bool outputReadLengths_;
    const reference::SortedReferenceMetadata sortedReferenceMetadata_;
    const std::vector<unsigned long> contigOffsets_;

public:

    AnnotationsMerger(
        const bfs::path &referenceGenome,
        const bfs::path &outputFilePath,
        const bool outputReadLengths);

    template <typename MergedT>
    void start(const boost::filesystem::path &firstFile, std::vector<MergedT> &annotation);

    template <typename ToMergeT, typename MergedT, typename AggregaT>
    void merge(const boost::filesystem::path &path, AggregaT agg, std::vector<MergedT> &annotation);

    template <typename MergedT>
    void finish(std::vector<MergedT> &annotation);

    template <typename MergedT, typename AggregaT>
    void run(
        const std::vector<boost::filesystem::path> mergeFiles,
        AggregaT agg);

private:
    template <typename AnnotationT>
    void initAnnotation(const std::string &initFunction, AnnotationT &annotation);

    template <typename AnnotationT>
    void postProcess(AnnotationT &annotation) const;
//    void computeMinAnchoringReadLength(std::vector<reference::DistanceToBeNeighborless> &annotation) const;
    template <typename AnnotationT>
    void loadFile(
        const boost::filesystem::path &path,
        std::vector<AnnotationT> &annotation);
};


template <typename AnnotationT>
void AnnotationsMerger::loadFile(
    const boost::filesystem::path &path,
    std::vector<AnnotationT> &annotation)
{
    const std::size_t annotationSize = reference::genomeLength(sortedReferenceMetadata_.getContigs());
    annotation.resize(annotationSize, 0);
    boost::iostreams::filtering_istream is;
    if (common::isDotGzPath(path))
    {
        is.push(boost::iostreams::gzip_decompressor());
    }
    is.push(boost::iostreams::file_source(path.string()));
    if (!is.read(reinterpret_cast<char*>(&annotation.front()), annotation.size() * sizeof(AnnotationT)))
    {
        const boost::format message = boost::format("Failed to read neighbor counts from %s: %s") % path.string() % strerror(errno);
        BOOST_THROW_EXCEPTION(common::IoException(errno, message.str()));
    }
}

template <typename MergedT>
void AnnotationsMerger::start(const boost::filesystem::path &firstFile, std::vector<MergedT> &annotation)
{
    if (!boost::filesystem::exists(firstFile))
    {
        initAnnotation(firstFile.string(), annotation);
    }
    else
    {
        ISAAC_THREAD_CERR << "loading: " << firstFile << std::endl;
        loadFile(firstFile, annotation);
    }
}

template <typename ToMergeT, typename MergedT, typename AggregaT>
void AnnotationsMerger::merge(const boost::filesystem::path &path, AggregaT agg, std::vector<MergedT> &annotation)
{
    ISAAC_THREAD_CERR << "merging: " << path << std::endl;
    const std::size_t annotationSize = reference::genomeLength(sortedReferenceMetadata_.getContigs());
    if (annotation.size() != annotationSize)
    {
        const boost::format message = boost::format("All annotations must be the same byte length. Expected %d, got %d") % annotation.size() % annotationSize;
        BOOST_THROW_EXCEPTION(common::IoException(errno, message.str()));
    }
    else
    {
        std::size_t readSoFar = 0;
        ToMergeT toMerge = 0;
        boost::iostreams::filtering_istream is;
        if (common::isDotGzPath(path))
        {
            is.push(boost::iostreams::gzip_decompressor());
        }
        is.push(boost::iostreams::file_source(path.string()));
        while (is.read(reinterpret_cast<char*>(&toMerge), sizeof(toMerge)))
        {
            ISAAC_ASSERT_MSG(annotationSize > readSoFar, "Annotation size is greater than expected " << annotationSize);
            const MergedT mergedValue = annotation.at(readSoFar);
            annotation.at(readSoFar) = agg(mergedValue, toMerge);

            ++readSoFar;
        }
        if (annotation.size() != readSoFar)
        {
            const boost::format message = boost::format("Failed to read %d values from %s. Read so far %d. Error %s") %
                annotation.size() % path.string() % readSoFar % strerror(errno);
            BOOST_THROW_EXCEPTION(common::IoException(errno, message.str()));
        }
    }
}
template <typename MergedT>
void AnnotationsMerger::finish(std::vector<MergedT> &annotation)
{
    postProcess(annotation);
//
//    if (outputReadLengths_)
//    {
//        computeMinAnchoringReadLength(annotation);
//    }

    boost::iostreams::filtering_ostream os;
    os.push(boost::iostreams::gzip_compressor(boost::iostreams::gzip_params(1)));
    os.push(boost::iostreams::file_sink(outputFilePath_.string()));
    if (!os.write(reinterpret_cast<const char*>(&annotation.front()),
                  annotation.size() * sizeof(MergedT)))
    {
        const boost::format message = boost::format("Failed to write merged annotation into %s: %s") % outputFilePath_.string() % strerror(errno);
        BOOST_THROW_EXCEPTION(common::IoException(errno, message.str()));
    }
}

template <typename MergedT, typename AggregaT>
void AnnotationsMerger::run(
    const std::vector<boost::filesystem::path> mergeFiles,
    AggregaT agg)
{
    std::vector<MergedT> annotation;
    start(mergeFiles.front(), annotation);
    BOOST_FOREACH(const boost::filesystem::path &mergeFile, std::make_pair(mergeFiles.begin() + 1, mergeFiles.end()))
    {
        merge<MergedT>(mergeFile, agg, annotation);
    }
    finish(annotation);
}
} // namespace reference
} // namespace isaac

#endif // #ifndef iSAAC_REFERENCE_ANNOTATIONS_MERGER_HH
