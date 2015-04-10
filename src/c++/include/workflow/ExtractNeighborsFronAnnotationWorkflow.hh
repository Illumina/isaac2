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
 ** \file ExtractNeighborsFronAnnotationWorkflow.hh
 **
 ** \brief Top level component to control the neighbor extraction process.
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_WORKFLOW_EXTRACT_NEIGHBORS_FROM_ANNOTATION_WORKFLOW_HH
#define iSAAC_WORKFLOW_EXTRACT_NEIGHBORS_FROM_ANNOTATION_WORKFLOW_HH

#include "common/Threads.hpp"
#include "reference/Contig.hh"
#include "reference/SortedReferenceXml.hh"

namespace isaac
{
namespace workflow
{

namespace bfs = boost::filesystem;

class ExtractNeighborsFromAnnotationWorkflow: boost::noncopyable
{
private:
    const bfs::path annotationFilePath_;
    const unsigned bitsPerValue_;
    const bfs::path neighborsFilePath_;
    const reference::SortedReferenceMetadata contigsXml_;

public:
    ExtractNeighborsFromAnnotationWorkflow(
        const bfs::path &contigsXmlPath,
        const bfs::path &highAnnotationFilePath,
        const unsigned bitsPerValue,
        const bfs::path &neighborsFilePath
        );

    void run();

private:
    template <typename ReadType>
    void extractNeighbors(
        std::istream &bitsetFile,
        std::vector<bool> &neighbors );
    void dumpResults(const std::vector<bool> &neighbors);
};
} // namespace workflow
} // namespace isaac

#endif // #ifndef iSAAC_WORKFLOW_EXTRACT_NEIGHBORS_FROM_ANNOTATION_WORKFLOW_HH
