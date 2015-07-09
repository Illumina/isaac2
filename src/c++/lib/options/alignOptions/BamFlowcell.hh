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
 ** \file FastqFlowcell.hh
 **
 ** Generate flowcell object out of BaseCalls/laneX.bam files
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_OPTIONS_ALIGN_OPTIONS_BAM_FLOWCELL_HH
#define iSAAC_OPTIONS_ALIGN_OPTIONS_BAM_FLOWCELL_HH

#include "flowcell/Layout.hh"
#include "reference/ReferenceMetadata.hh"

namespace isaac
{
namespace options
{
namespace alignOptions
{

struct BamFlowcellInfo
{
    BamFlowcellInfo() : readNameLength_(0){}

    std::string flowcellId_;
    std::pair<unsigned, unsigned> readLengths_;
    std::vector<unsigned> lanes_;
    unsigned readNameLength_;

    const std::vector<unsigned> &getLanes() const
    {
        return lanes_;
    }
};


class BamFlowcell : boost::noncopyable
{
    static const unsigned MAX_LANE_NUMBER = 8;

public:
    static flowcell::Layout createFilteredFlowcell(
        const std::string &tilesFilter,
        const boost::filesystem::path &baseCallsDirectory,
        const unsigned laneNumberMax,
        const unsigned readNameLength,
        std::string useBasesMask,
        const bool allowVariableReadLength,
        const std::string &seedDescriptor,
        const unsigned seedLength,
        const reference::ReferenceMetadataList &referenceMetadataList,
        unsigned &firstPassSeeds);

private:
    struct BamPath
    {
        unsigned lane_;
        boost::filesystem::path path_;
    };
    static BamPath findBamPath(
        const boost::filesystem::path &baseCallsDirectory);
    static BamFlowcellInfo parseBamFlowcellInfo(
        const BamPath &laneFilePaths,
        const bool allowVariableReadLength,
        const bool allowMixedFlowcells,
        const unsigned readNameLength);

};

inline std::ostream &operator<< (std::ostream &os, const BamFlowcellInfo &fcInfo)
{
    os << "BamFlowcellInfo("<<
        fcInfo.flowcellId_ << "," <<
        fcInfo.readLengths_.first << ":" <<
        fcInfo.readLengths_.second << ",[";

    BOOST_FOREACH(const unsigned lane, fcInfo.getLanes())
    {
        os << lane << " ";
    }
    return os << "])";
}

} // namespace alignOptions
} // namespace options
} // namespace isaac

#endif // #ifndef iSAAC_OPTIONS_ALIGN_OPTIONS_BAM_FLOWCELL_HH
