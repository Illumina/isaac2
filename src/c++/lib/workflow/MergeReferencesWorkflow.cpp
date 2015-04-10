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
#include "workflow/MergeReferencesWorkflow.hh"

namespace isaac
{
namespace workflow
{

reference::DistanceToBeNeighborless minDistanceToBeNeighborless(
    const reference::DistanceToBeNeighborless merged,
    const reference::DistanceToBeNeighborless toMerge)
{
    return std::min(merged, toMerge);
}



MergeReferencesWorkflow::MergeReferencesWorkflow(
    const std::vector<bfs::path> &filesToMerge,
    const bfs::path &outputFilePath,
    const bool mergeAnnotations,
    const bool makeAbsolutePaths)
    : filesToMerge_(filesToMerge),
      outputFilePath_(outputFilePath),
      mergeAnnotations_(mergeAnnotations),
      makeAbsolutePaths_(makeAbsolutePaths)
{
}

void MergeReferencesWorkflow::run()
{
    reference::SortedReferenceMetadata result;
    reference::SortedReferenceMetadata::AnnotationFiles annotationFiles;
    BOOST_FOREACH(const boost::filesystem::path &path, filesToMerge_)
    {
        reference::SortedReferenceMetadata referenceToMerge = reference::loadSortedReferenceXml(path, makeAbsolutePaths_);
        if (mergeAnnotations_)
        {
            annotationFiles.push_back(referenceToMerge.getKUniquenessAnnotation());
            referenceToMerge.clearAnnotations();
        }
        result.merge(referenceToMerge);
    }
    const reference::SortedReferenceMetadata::Contigs karyotypeOrderedContigs = result.getKaryotypeOrderedContigs();
    const reference::SortedReferenceMetadata::Contigs::const_iterator collision = std::adjacent_find(
        karyotypeOrderedContigs.begin(), karyotypeOrderedContigs.end(),
        boost::bind(&reference::SortedReferenceMetadata::Contig::karyotypeIndex_, _1) ==
            boost::bind(&reference::SortedReferenceMetadata::Contig::karyotypeIndex_, _2));

    if (karyotypeOrderedContigs.end() != collision)
    {
        const boost::format message = boost::format("\n   *** Karyotype index collision detected in %s and %s ***\n") % *collision % *(collision + 1);
        BOOST_THROW_EXCEPTION(common::PostConditionException(message.str()));
    }

    if (mergeAnnotations_)
    {
        reference::SortedReferenceMetadata::AnnotationFiles::const_iterator oddAnnotation = std::find_if(
                        annotationFiles.begin(), annotationFiles.end(),
                        boost::bind(&reference::SortedReferenceMetadata::AnnotationFile::k_, _1) != annotationFiles.front().k_);
        if (annotationFiles.end() != oddAnnotation)
        {
            const boost::format message = boost::format("\n   *** Annotations expected for the same k (%d) found %s ***\n") % annotationFiles.front().k_ % *oddAnnotation;
            BOOST_THROW_EXCEPTION(common::PreConditionException(message.str()));
        }

        std::vector<boost::filesystem::path> paths;
        std::transform(annotationFiles.begin(), annotationFiles.end(), std::back_inserter(paths),
                       boost::bind(&reference::SortedReferenceMetadata::AnnotationFile::path_, _1));

        const boost::filesystem::path outputAnnotationPath = (outputFilePath_.parent_path() / annotationFiles.front().path_.filename());

        reference::AnnotationsMerger merger(filesToMerge_.front(), outputAnnotationPath, false);
        merger.run<reference::DistanceToBeNeighborless>(paths, &minDistanceToBeNeighborless);
        result.setKUniquenessAnnotation(outputAnnotationPath, annotationFiles.front().k_);
    }
    reference::saveSortedReferenceXml(outputFilePath_, result);
}

} // namespace workflow
} // namespace isaac
