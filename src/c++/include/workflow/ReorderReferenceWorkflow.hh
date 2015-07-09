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
 ** \file ReorderReferenceWorkflow.hh
 **
 ** \brief Top level component to control the reference reordering process.
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_WORKFLOW_REORDER_REFERENCE_WORKFLOW_HH
#define iSAAC_WORKFLOW_REORDER_REFERENCE_WORKFLOW_HH

#include "common/Threads.hpp"
#include "reference/Contig.hh"
#include "reference/SortedReferenceXml.hh"

namespace isaac
{
namespace workflow
{

namespace bfs = boost::filesystem;

class ReorderReferenceWorkflow: boost::noncopyable
{
public:
    ReorderReferenceWorkflow(
        const bfs::path &sortedReferenceMetadata,
        const bfs::path &newXmlPath,
        const bfs::path &newDataFileDirectory,
        const std::vector<std::string> &newOrder,
        const unsigned basesPerLine
        );

    void run();

private:
    const bfs::path sortedReferenceMetadata_;
    const bfs::path newXmlPath_;
    const bfs::path newDataFileDirectory_;
    const std::vector<std::string> &newOrder_;
    const unsigned basesPerLine_;
    bool sameOrder_;
    common::ThreadVector threads_;

    reference::SortedReferenceMetadata xml_;
    // translation array from new karyotype indexes to the original ones
    std::vector<unsigned> originalKaryotypeIndexes_;
    std::vector<unsigned long> originalContigOffsets_;
    bool orderByKaryotypeIndex(const reference::Contig& left, const reference::Contig& right);
    void storeContig(
        std::ostream &os,
        const reference::Contig &contig,
        const boost::filesystem::path &fastaPath);
    void writeBase(std::ostream &os, const char base, const bool writeNewline);
    void reorderContigs();
    void reorderAnnotation();
};
} // namespace workflow
} // namespace isaac

#endif // #ifndef iSAAC_WORKFLOW_REORDER_REFERENCE_WORKFLOW_HH
