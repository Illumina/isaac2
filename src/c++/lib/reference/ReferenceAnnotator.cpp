/**
 ** Isaac Genome Alignment Software
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
 ** \file ReferenceAnnotator.cpp
 **
 ** \brief See ReferenceAnnotator.hh.
 **
 ** \author Roman Petrovski
 **/

#include <fstream>
#include <cerrno>
#include <cstring>
#include <ctime>
#include <boost/thread.hpp>
#include <boost/bind.hpp>
#include <boost/format.hpp>
#include <boost/foreach.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include "common/Debug.hh"
#include "reference/ReferenceAnnotator.hh"
#include "reference/SortedReferenceXml.hh"
#include "reference/ReferenceKmer.hh"
#include "common/Exceptions.hh"
#include "oligo/Kmer.hh"
#include "oligo/Permutate.hh"

namespace isaac
{
namespace reference
{

namespace bfs = boost::filesystem;

template <typename KmerT>
ReferenceAnnotator<KmerT>::ReferenceAnnotator(
    const bfs::path &inputFile,
    const bfs::path &outputDirectory,
    const bfs::path &outputFile,
    const bfs::path &tempFile)
    : inputFile_(inputFile)
    , outputDirectory_(outputDirectory)
    , outputFile_(outputFile)
    , tempFile_(tempFile)
    , sortedReferenceMetadata_(loadSortedReferenceXml(inputFile_))
    , contigOffsets_(computeContigOffsets(sortedReferenceMetadata_.getContigs()))
{
}

template <typename KmerT>
void ReferenceAnnotator<KmerT>::run() const
{
    std::vector<NeighborsCount> neighborsAtBaseMismatch(reference::genomeLength(sortedReferenceMetadata_.getContigs()));

    using common::IoException;
    using boost::format;
    std::ifstream is(tempFile_.string().c_str());
    if (!is)
    {
        const format message = format("Failed to open neighbor counts file %s for reading: %s") % tempFile_ % strerror(errno);
        BOOST_THROW_EXCEPTION(IoException(errno, message.str()));
    }

    if (!is.read(reinterpret_cast<char*>(&neighborsAtBaseMismatch.front()),
                        neighborsAtBaseMismatch.size() * sizeof(NeighborsCount)))
    {
        const boost::format message = boost::format("Failed to read neighbor counts from %s: %s") % outputFile_.string() % strerror(errno);
        BOOST_THROW_EXCEPTION(common::IoException(errno, message.str()));
    }

    SortedReferenceMetadata sortedReferenceMetadata = sortedReferenceMetadata_;
    updateSortedReference(sortedReferenceMetadata.getMaskFileList(oligo::KmerTraits<KmerT>::KMER_BASES), neighborsAtBaseMismatch);
    saveSortedReferenceXml(outputFile_, sortedReferenceMetadata);
}


template <typename KmerT>
void ReferenceAnnotator<KmerT>::updateSortedReference(
    std::vector<SortedReferenceMetadata::MaskFile> &maskFileList,
    const std::vector<NeighborsCount> &neighborsAtBaseMismatch) const
{
    BOOST_FOREACH(SortedReferenceMetadata::MaskFile &maskFile, maskFileList)
    {
        const bfs::path oldMaskFile = maskFile.path; //bfs::path(maskFile.path).replace_extension(".orig");
        ISAAC_THREAD_CERR << "Annotating " << oldMaskFile << std::endl;
        if (!exists(oldMaskFile))
        {
            const boost::format message = boost::format("Mask file %s does not exist: %s") % oldMaskFile;
            BOOST_THROW_EXCEPTION(common::IoException(ENOENT, message.str()));
        }

        maskFile.path = outputDirectory_ / maskFile.path.filename();
        std::ifstream maskInput(oldMaskFile.string().c_str());
        if (!maskInput)
        {
            const boost::format message = boost::format("Failed to open mask file %s for reading: %s") % oldMaskFile % strerror(errno);
            BOOST_THROW_EXCEPTION(common::IoException(errno, message.str()));
        }
        std::ofstream maskOutput(maskFile.path.string().c_str());
        if (!maskOutput)
        {
            const boost::format message = boost::format("Failed to open mask file %s for writing: %s") % maskFile.path % strerror(errno);
            BOOST_THROW_EXCEPTION(common::IoException(errno, message.str()));
        }
        KmerT lastTooManyMatchKmer = ~KmerT(0);
        while(maskInput && maskOutput)
        {
            ReferenceKmer<KmerT> referenceKmer;
            if (maskInput.read(reinterpret_cast<char *>(&referenceKmer), sizeof(referenceKmer)))
            {
                if (!referenceKmer.isTooManyMatch())
                {
                    const std::size_t referenceOffset =
                        contigOffsets_.at(referenceKmer.getReferencePosition().getContigId()) +
                        referenceKmer.getReferencePosition().getPosition();

//                    const NeighborsCount neighbors = *std::max_element(
//                        neighborsAtBaseMismatch.begin() + referenceOffset,
//                        neighborsAtBaseMismatch.begin() + referenceOffset + oligo::KmerTraits<KmerT>::KMER_BASES);
                    const NeighborsCount neighbors = neighborsAtBaseMismatch.at(referenceOffset);
                    // don't store those that overlap high mutant positions
                    if (neighbors < 10)
                    {
                        referenceKmer.setNeighbors(neighbors);
                    }
                    else
                    {
                        referenceKmer.setTooManyMatch();
                    }
                }
                else
                {
                    referenceKmer.setTooManyMatch();
                }

                // only one entry per too-many-match kmer is allowed
                if (!referenceKmer.isTooManyMatch() || referenceKmer.getKmer() != lastTooManyMatchKmer)
                {
                    if (!maskOutput.write(reinterpret_cast<const char *>(&referenceKmer), sizeof(referenceKmer)))
                    {
                        const boost::format message = boost::format("Failed to write toomanymatch k-mer into %s: %s") % maskFile.path % strerror(errno);
                        BOOST_THROW_EXCEPTION(common::IoException(errno, message.str()));
                    }
                    lastTooManyMatchKmer = referenceKmer.getKmer();
                }
            }
        }
        if (!maskInput.eof())
        {
            const boost::format message = boost::format("Failed to update %s with neighbors information: %s") % maskFile.path % strerror(errno);
            BOOST_THROW_EXCEPTION(common::IoException(errno, message.str()));
        }
        ISAAC_THREAD_CERR << "Adding neighbors information done for " << maskFile.path << std::endl;
    }
}

template class ReferenceAnnotator<isaac::oligo::LongKmerType>;

} // namespace reference
} //namespace isaac
