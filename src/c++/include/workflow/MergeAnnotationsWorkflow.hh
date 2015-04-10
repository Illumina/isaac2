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
 ** \file MergeAnnotationsWorkflow.hh
 **
 ** \brief merge multiple reference annotations into one
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_WORKFLOW_MERGE_ANNOTATIONS_WORKFLOW_HH
#define iSAAC_WORKFLOW_MERGE_ANNOTATIONS_WORKFLOW_HH

#include "reference/SortedReferenceXml.hh"

namespace isaac
{
namespace workflow
{

namespace bfs = boost::filesystem;

class MergeAnnotationsWorkflow: boost::noncopyable
{
private:
    const boost::filesystem::path referenceGenome_;
    const std::vector<bfs::path> &filesToMerge_;
    const bfs::path &outputFilePath_;
    const std::vector<std::string> aggregateFunction_;
    const std::string &mergedType_;

public:
    MergeAnnotationsWorkflow(
        const boost::filesystem::path referenceGenome,
        const std::vector<bfs::path> &filesToMerge,
        const bfs::path &outputFilePath,
        const std::vector<std::string> &aggregateFunction,
        const std::string &mergedType);

    void run();
};
} // namespace workflow
} // namespace isaac

#endif // #ifndef iSAAC_WORKFLOW_MERGE_ANNOTATIONS_WORKFLOW_HH
