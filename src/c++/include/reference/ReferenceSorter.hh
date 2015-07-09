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
 ** \file ReferenceSorter.hh
 **
 ** Top level component to produce a sorted reference.
 **
 ** \author Come Raczy
 **/

#ifndef iSAAC_REFERENCE_REFERENCE_SORTER_HH
#define iSAAC_REFERENCE_REFERENCE_SORTER_HH

#include <string>
#include <boost/noncopyable.hpp>
#include <boost/filesystem.hpp>

#include "oligo/Kmer.hh"
#include "PermutatedKmerGenerator.hh"
#include "reference/ReferenceKmer.hh"
#include "reference/ReferencePosition.hh"

namespace isaac
{
namespace reference
{

template <typename KmerT>
class ReferenceSorter: boost::noncopyable
{
public:
    ReferenceSorter(
        const unsigned int maskWidth,
        const unsigned mask,
        const boost::filesystem::path &contigsXmlPath,
        const boost::filesystem::path &genomeNeighborsFile,
        const boost::filesystem::path &outputFile,
        const unsigned repeatThreshold);
    void run();
private:
    const unsigned repeatThreshold_;
    const unsigned int maskWidth_;
    const unsigned mask_;

    // the mask highlight bits in original kmer (ABCD)
    const KmerT msbMask_;
    // the mask value in the original kmer (ABCD) shifted to the topmost position
    const KmerT maskBits_;

    const boost::filesystem::path contigsXmlPath_;
    const boost::filesystem::path genomeNeighborsFile_;

    const boost::filesystem::path outputFile_;
    const reference::SortedReferenceMetadata sortedReferenceMetadata_;
    common::ThreadVector threads_;
    const reference::ContigList contigList_;

    std::vector<ReferenceKmer<KmerT> > reference_;

    void loadReference();
    void sortReference();
    void markRepeats();
    void saveReference();
};

} // namespace reference
} // namespace isaac

#endif // #ifndef iSAAC_REFERENCE_REFERNECE_SORTER_HH
