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
 ** \file ReferenceAnnotator.hh
 **
 ** \brief Top level component to annotate reference.
 **
 ** \author Roman Petrovski
 **/

#ifndef ISAAC_REFERENCE_REFERENCE_ANNOTATOR_HH
#define ISAAC_REFERENCE_REFERENCE_ANNOTATOR_HH

#include <vector>
#include <boost/filesystem.hpp>
#include <boost/noncopyable.hpp>

#include "oligo/Kmer.hh"
#include "reference/SortedReferenceMetadata.hh"

namespace isaac
{
namespace reference
{

template <typename KmerT>
class ReferenceAnnotator: boost::noncopyable
{
public:
    ReferenceAnnotator(
        const boost::filesystem::path &inputFile,
        const boost::filesystem::path &outputDirectory,
        const boost::filesystem::path &outputFile,
        const boost::filesystem::path &tempFile);
    void run() const;
private:
    const boost::filesystem::path inputFile_;
    const boost::filesystem::path outputDirectory_;
    const boost::filesystem::path outputFile_;
    const boost::filesystem::path tempFile_;
    const SortedReferenceMetadata sortedReferenceMetadata_;
    const std::vector<unsigned long> contigOffsets_;
    typedef unsigned NeighborsCount;

    void updateSortedReference(
        std::vector<SortedReferenceMetadata::MaskFile> &maskFileList,
        const std::vector<NeighborsCount> &neighborsAtBaseMismatch) const;
};

} // namespace reference
} //namespace isaac

#endif // #ifndef ISAAC_REFERENCE_REFERENCE_ANNOTATOR_HH
