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
 ** \file BpbToWigWorkflow.hh
 **
 ** \brief prints wig file given the bitset and contigs map
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_WORKFLOW_BPB_TO_WIG_WORKFLOW_HH
#define iSAAC_WORKFLOW_BPB_TO_WIG_WORKFLOW_HH

#include "common/Threads.hpp"
#include "reference/ContigLoader.hh"
#include "reference/SortedReferenceXml.hh"

namespace isaac
{
namespace workflow
{

namespace bfs = boost::filesystem;

class BpbToWigWorkflow: boost::noncopyable
{
    const unsigned wigDefaultValue_;
    const bfs::path sortedReferenceMetadata_;
    const bfs::path inputFilePath_;
    const std::string outputFormatString_;
    const unsigned bitsPerValue_;
    const bool knownSitesOnly_;
    const bool bedPrintAllPositions_;

    const reference::SortedReferenceMetadata xml_;
    common::ThreadVector threads_;
    const reference::ContigList contigs_;

public:
    BpbToWigWorkflow(
        const unsigned wigDefaultValue,
        const bfs::path &sortedReferenceMetadata,
        const bfs::path &inputFilePath,
        const std::string &outputFormatString,
        const unsigned bitsPerValue,
        const bool knownSitesOnly,
        const bool bedPrintAllPositions
        );

    void run();
private:

    template <typename ReadType>
    void printWig(
        std::istream &bitsetFile,
        const unsigned bitsPerValue,
        const reference::SortedReferenceMetadata::Contigs &contigs );

    template <typename ReadType>
    void printBed(
        std::istream &bitsetFile,
        const unsigned bitsPerValue,
        const reference::SortedReferenceMetadata::Contigs &contigs );

};
} // namespace workflow
} // namespace isaac

#endif // #ifndef iSAAC_WORKFLOW_BPB_TO_WIG_WORKFLOW_HH
