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
 ** \file BarcodeResolver.hh
 **
 ** Performs translation from barcode sequences to barcode indexes. Allows for sequence mismatches.
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_DEMULTIPLEXING_BARCODE_RESOLVER_HH
#define iSAAC_DEMULTIPLEXING_BARCODE_RESOLVER_HH

#include "demultiplexing/Barcode.hh"
#include "demultiplexing/DemultiplexingStats.hh"
#include "flowcell/BarcodeMetadata.hh"
#include "flowcell/TileMetadata.hh"

namespace isaac
{
namespace demultiplexing
{

class BarcodeResolver: boost::noncopyable
{
public:
    BarcodeResolver(
        const flowcell::TileMetadataList &allTilesMetadata,
        const flowcell::BarcodeMetadataList &allBarcodeMetadata,
        const flowcell::BarcodeMetadataList &barcodeGroup);

    /**
     * \brief updates the tile information in 'result' with the corresponding barcodeMetadataList_ indexes
     */
    void resolve(
        std::vector<Barcode> &barcodes,
        demultiplexing::DemultiplexingStats &demultiplexingStats);

    static unsigned getMismatchKmersCount(const unsigned kmerLength, const unsigned maxMismatches);

    static void generateBarcodeMismatches(
        const flowcell::BarcodeMetadata &barcodeMetadata,
        std::vector<Barcode> &result);

    static std::vector<Barcode> generateMismatches(
        const flowcell::BarcodeMetadataList &allBarcodeMetadata,
        const flowcell::BarcodeMetadataList &barcodeGroup);

    static std::pair<Kmer, unsigned> get1MismatchKmer(const Kmer original, const unsigned kmerLength,
                                                      const unsigned componentOffset,
                                                      const unsigned iteration);
    static std::pair<Kmer, unsigned> get2MismatchKmer(const Kmer original, const unsigned kmerLength,
                                                      const unsigned componentOffset,
                                                      const unsigned iteration);

private:
    const flowcell::TileMetadataList &allTilesMetadata_;
    const flowcell::BarcodeMetadataList &allBarcodeMetadata_;
    const std::vector<Barcode> mismatchBarcodes_;
    const unsigned unknownBarcodeIndex_;
    std::vector<unsigned long> barcodeHits_;
};

} // namespace demultiplexing
} // namespace isaac

#endif // #ifndef iSAAC_DEMULTIPLEXING_BARCODE_RESOLVER_HH
