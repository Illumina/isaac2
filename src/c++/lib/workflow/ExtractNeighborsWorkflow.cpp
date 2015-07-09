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
 ** \file ExtractNeighborsWorkflow.cpp
 **
 ** \brief see ExtractNeighborsWorkflow.hh
 **
 ** \author Roman Petrovski
 **/

#include <boost/foreach.hpp>
#include <boost/format.hpp>

#include "common/Debug.hh"
#include "common/Exceptions.hh"
#include "io/BitsetSaver.hh"
#include "reference/ContigLoader.hh"
#include "reference/ReferenceKmer.hh"
#include "workflow/ExtractNeighborsWorkflow.hh"

namespace isaac
{
namespace workflow
{

ExtractNeighborsWorkflow::ExtractNeighborsWorkflow(
    const bfs::path &sortedReferenceMetadata,
    const bfs::path &neighborsFilePath,
    const bfs::path &highRepeatsFilePath
    )
    : sortedReferenceMetadata_(sortedReferenceMetadata),
      neighborsFilePath_(neighborsFilePath),
      highRepeatsFilePath_(highRepeatsFilePath),
      threads_(boost::thread::hardware_concurrency()),
      xml_(reference::loadSortedReferenceXml(sortedReferenceMetadata_))
{
}

template <typename KmerT>
void ExtractNeighborsWorkflow::run()
{
    const reference::SortedReferenceMetadata::MaskFiles &maskFiles =
        xml_.getMaskFileList(oligo::KmerTraits<KmerT>::KMER_BASES);
    if (maskFiles.empty())
    {
        BOOST_THROW_EXCEPTION(isaac::common::PreConditionException("No mask files in " + sortedReferenceMetadata_.string()));
    }

    const reference::SortedReferenceMetadata::Contigs contigs = xml_.getContigs();
    const std::vector<unsigned long> contigOffsets = reference::computeContigOffsets(contigs);

    std::vector<unsigned> karyotypes;
    karyotypes.reserve(xml_.getContigs().size());
    std::transform(xml_.getContigs().begin(), xml_.getContigs().end(), std::back_inserter(karyotypes),
                   boost::bind(&reference::SortedReferenceMetadata::Contig::karyotypeIndex_, _1));


    std::vector<bool> neighbors(reference::genomeLength(contigs), false);
    std::vector<bool> highRepeats(highRepeatsFilePath_.empty() ? 0 : reference::genomeLength(contigs), true);

    // there could be multiple mask widths in the xml. Just fail if there are.
    unsigned maskWidth = -1U;

    BOOST_FOREACH(const reference::SortedReferenceMetadata::MaskFile &maskFile, maskFiles)
    {
        if (-1U == maskWidth)
        {
            maskWidth = maskFile.maskWidth;
        }
        ISAAC_ASSERT_MSG(maskWidth == maskFile.maskWidth, "Mixed mask widths are not supported");
        scanMaskFile<KmerT>(maskFile, contigOffsets, karyotypes, neighbors, highRepeats);
    }

    dumpResults(neighbors, highRepeats);
}

template <typename KmerT>
void ExtractNeighborsWorkflow::scanMaskFile(
    const reference::SortedReferenceMetadata::MaskFile &maskFile,
    const std::vector<unsigned long> &contigOffsets,
    const std::vector<unsigned> &karyotypes,
    std::vector<bool> &neighbors,
    std::vector<bool> &highRepeats)
{
    if (!exists(maskFile.path))
    {
        const boost::format message = boost::format("Mask file %s does not exist: %s") % maskFile.path;
        BOOST_THROW_EXCEPTION(common::IoException(ENOENT, message.str()));
    }

    std::ifstream maskInput(maskFile.path.c_str());
    if (!maskInput)
    {
        const boost::format message = boost::format("Failed to open mask file %s for reading: %s") % maskFile.path % strerror(errno);
        BOOST_THROW_EXCEPTION(common::IoException(errno, message.str()));
    }

    std::size_t scannedKmers = 0, maskNeighbors = 0, maskNonHighRepeats = 0;
    while(maskInput)
    {
        reference::ReferenceKmer<KmerT> referenceKmer;
        if (maskInput.read(reinterpret_cast<char *>(&referenceKmer), sizeof(referenceKmer)))
        {
            ++scannedKmers;
            const reference::ReferencePosition pos = referenceKmer.getReferencePosition();
            if (!pos.isTooManyMatch())
            {
                const reference::ReferencePosition translatedPos = pos.translateContig(karyotypes);
                if (translatedPos.hasNeighbors())
                {
                    neighbors.at(contigOffsets.at(translatedPos.getContigId()) + translatedPos.getPosition()) = true;
                    ++maskNeighbors;
                }

                if (!highRepeats.empty())
                {
                    highRepeats.at(contigOffsets.at(translatedPos.getContigId()) + translatedPos.getPosition()) = false;
                }
                ++maskNonHighRepeats;
            }
        }
    }
    if (!maskInput.eof())
    {
        const boost::format message = boost::format("Failed to scan %s to the end") % maskFile.path % strerror(errno);
        BOOST_THROW_EXCEPTION(common::IoException(errno, message.str()));
    }
    else
    {
        ISAAC_THREAD_CERR << "Scanning " << maskFile.path << " found " << scannedKmers << " kmers of which " <<
            maskNeighbors << " neighbors and " << (scannedKmers - maskNonHighRepeats) << " high repeats" << std::endl;
    }
}

void ExtractNeighborsWorkflow::dumpResults(const std::vector<bool> &neighbors, const std::vector<bool> &highRepeats)
{
    io::BitsetSaver neighborsSaver(neighborsFilePath_);
    neighborsSaver.save(neighbors);
    ISAAC_THREAD_CERR << "Stored " << neighbors.size() << " neighbor locations in " << neighborsFilePath_ << std::endl;
    if (!highRepeatsFilePath_.empty())
    {
        io::BitsetSaver highRepeatsSaver(highRepeatsFilePath_);
        highRepeatsSaver.save(highRepeats);
        ISAAC_THREAD_CERR << "Stored " << highRepeats.size() << " high repeats locations in " << highRepeatsFilePath_ << std::endl;
    }
}

template void  ExtractNeighborsWorkflow::run<oligo::ShortKmerType>();
template void  ExtractNeighborsWorkflow::run<oligo::BasicKmerType<28> >();
template void  ExtractNeighborsWorkflow::run<oligo::BasicKmerType<30> >();
template void  ExtractNeighborsWorkflow::run<oligo::KmerType>();
template void  ExtractNeighborsWorkflow::run<oligo::BasicKmerType<34> >();
template void  ExtractNeighborsWorkflow::run<oligo::BasicKmerType<36> >();
template void  ExtractNeighborsWorkflow::run<oligo::LongKmerType>();

} // namespace workflow
} // namespace isaac
