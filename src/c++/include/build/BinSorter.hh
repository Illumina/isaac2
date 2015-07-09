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
 ** \file BinSorter.hh
 **
 ** Performs sorting and duplicate marking on a single alignment bin.
 ** 
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_BUILD_BIN_SORTER_HH
#define iSAAC_BUILD_BIN_SORTER_HH

#include <numeric>

#include <boost/filesystem.hpp>
#include <boost/foreach.hpp>
#include <boost/ptr_container/ptr_vector.hpp>

#include "alignment/BinMetadata.hh"
#include "build/BamSerializer.hh"
#include "build/BinData.hh"
#include "build/BinLoader.hh"
#include "build/DuplicateFragmentIndexFiltering.hh"
#include "build/DuplicatePairEndFilter.hh"
#include "build/FragmentIndex.hh"
#include "build/GapRealigner.hh"
#include "build/NotAFilter.hh"
#include "build/ParallelGapRealigner.hh"
#include "flowcell/TileMetadata.hh"
#include "io/FileBufCache.hh"
#include "io/WigLoader.hh"
#include "reference/AnnotationLoader.hh"

namespace isaac
{
namespace build
{

class BinSorter
{
public:
    BinSorter(
        const bool singleLibrarySamples,
        const bool keepDuplicates,
        const bool markDuplicates,
        const bool anchorMate,
        const BarcodeBamMapping &barcodeBamMapping,
        const flowcell::BarcodeMetadataList &barcodeMetadataList,
        const isaac::reference::ContigLists &contigLists,
        const unsigned splitGapLength,
        const isaac::reference::ContigAnnotationsList &annotations) :
            singleLibrarySamples_(singleLibrarySamples),
            keepDuplicates_(keepDuplicates),
            markDuplicates_(markDuplicates),
            anchorMate_(anchorMate),
            barcodeMetadataList_(barcodeMetadataList),
            contigLists_(contigLists),
            bamSerializer_(barcodeBamMapping.getSampleIndexMap(), splitGapLength),
            annotations_(annotations)
    {
    }

    void resolveDuplicates(
        BinData &binData,
        BuildStats &buildStats);

    std::size_t serialize(
        BinData &binData,
        boost::ptr_vector<boost::iostreams::filtering_ostream> &bgzfStreams,
        boost::ptr_vector<bam::BamIndexPart> &bamIndexParts);

private:
    const bool singleLibrarySamples_;
    const bool keepDuplicates_;
    const bool markDuplicates_;
    const bool anchorMate_;
    const flowcell::BarcodeMetadataList &barcodeMetadataList_;
    const std::vector<std::vector<isaac::reference::Contig> > &contigLists_;
    BamSerializer bamSerializer_;
    const isaac::reference::ContigAnnotationsList &annotations_;

    typedef boost::iterator_range<const unsigned char *> AnchorRange;

    template <typename BclIteratorT>
    bool isAnchored(
        const isaac::reference::ContigList &reference,
        const isaac::reference::Annotation &annotation,
        const isaac::reference::ReferencePosition &fpos,
        BclIteratorT currentBase,
        const uint32_t *cigarBegin,
        const uint32_t *cigarEnd) const;
    void downgradeAlignmentScores(
        const PackedFragmentBuffer::Index& idx,
        PackedFragmentBuffer &data);

};


} // namespace build
} // namespace isaac

#endif // #ifndef iSAAC_BUILD_BIN_SORTER_HH
