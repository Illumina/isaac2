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
 ** \file ClusterSeedGenerator.hh
 **
 ** \brief Component to generate the seeds from a block of sequentially-stored bcl clusters
 ** 
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_ALIGNMENT_CLUSTER_SEED_GENERATOR_HH
#define iSAAC_ALIGNMENT_CLUSTER_SEED_GENERATOR_HH

#include <vector>
#include <boost/filesystem.hpp>
#include <boost/noncopyable.hpp>
#include <boost/thread/mutex.hpp>

#include "alignment/Seed.hh"
#include "alignment/SeedGeneratorBase.hh"
#include "alignment/SeedMetadata.hh"
#include "common/Threads.hpp"
#include "alignment/matchFinder/TileClusterInfo.hh"
#include "flowcell/BarcodeMetadata.hh"
#include "flowcell/TileMetadata.hh"
#include "flowcell/ReadMetadata.hh"
#include "oligo/Kmer.hh"
#include "reference/SortedReferenceMetadata.hh"

namespace isaac
{
namespace alignment
{

/**
 ** \brief Encapsulates the variables that are shared by all the threads while
 ** loading the seeds.
 **/
template <typename KmerT>
class ClusterSeedGenerator : private SeedGeneratorBase<KmerT>
{
    typedef SeedGeneratorBase<KmerT> BaseT;
public:
    /**
     ** \brief constructs an instance with all the required shorthands.
     **
     ** Note: all parameters are kept as references and it is the
     ** responsibility of the caller to ensure appropriate life time for the
     ** referenced variables.
     **/
    ClusterSeedGenerator(
        common::ThreadVector &threads,
        const unsigned computeThreadsMax,
        const unsigned seedBaseQualityMin,
        const flowcell::BarcodeMetadataList &barcodeMetadataList,
        const flowcell::Layout &flowcellLayout,
        const std::vector<SeedMetadata> &seedMetadataList,
        const reference::SortedReferenceMetadataList &sortedReferenceMetadataList,
        const flowcell::TileMetadataList &willRequestTiles,
        const BclClusters &clusters,
        const flowcell::TileMetadataList &loadedTiles);

    void generateSeeds(const flowcell::TileMetadataList &tiles,
                   const matchFinder::TileClusterInfo &tileClusterBarcode,
                   std::vector<Seed<KmerT> > &seeds,
                   common::ScopedMallocBlock  &mallocBlock);

    /**
     * \brief returned iterators of the vector point past the last tile for each of the references
     */
    const std::vector<typename std::vector<Seed<KmerT> >::iterator> &getReferenceSeedBounds() const
    {
        return BaseT::nextTileSeedBegins_;
    }

private:
    // The mutex guards acquisition of the next tile and the destination of the seeds
    boost::mutex mutex_;
    const BclClusters &clusters_;
    const flowcell::TileMetadataList &loadedTiles_;

    const unsigned computeThreadsMax_;
    const unsigned seedBaseQualityMin_;
    /**
     * \brief Geometry: [thread][reference]
     */
    std::vector<std::vector<typename std::vector<Seed<KmerT> >::iterator> > threadDestinations_;

    common::ThreadVector &threads_;

    void generateThread(
        const matchFinder::TileClusterInfo &tileClusterBarcode,
        std::vector<flowcell::TileMetadata>::const_iterator &nextTile,
        const std::vector<flowcell::TileMetadata>::const_iterator tilesBegin,
        const std::vector<flowcell::TileMetadata>::const_iterator tilesEnd,
        const unsigned threadNumber);
};

} // namespace alignment
} // namespace isaac

#endif // #ifndef iSAAC_ALIGNMENT_CLUSTER_SEED_GENERATOR_HH
