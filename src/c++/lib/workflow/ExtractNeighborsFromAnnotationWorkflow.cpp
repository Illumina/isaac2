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
 ** \file ExtractNeighborsFronAnnotationWorkflow.cpp
 **
 ** \brief see ExtractNeighborsFronAnnotationWorkflow.hh
 **
 ** \author Roman Petrovski
 **/

#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>

#include "common/Debug.hh"
#include "common/Exceptions.hh"
#include "common/FileSystem.hh"
#include "io/BitsetSaver.hh"
#include "reference/ContigLoader.hh"
#include "reference/ReferenceKmer.hh"
#include "workflow/ExtractNeighborsFronAnnotationWorkflow.hh"

namespace isaac
{
namespace workflow
{

ExtractNeighborsFromAnnotationWorkflow::ExtractNeighborsFromAnnotationWorkflow(
    const bfs::path &contigsXmlPath,
    const bfs::path &annotationFilePath,
    const unsigned bitsPerValue,
    const bfs::path &neighborsFilePath
    )
    : annotationFilePath_(annotationFilePath),
      bitsPerValue_(bitsPerValue),
      neighborsFilePath_(neighborsFilePath),
      contigsXml_(reference::loadSortedReferenceXml(contigsXmlPath))
{
}

void ExtractNeighborsFromAnnotationWorkflow::run()
{
    boost::iostreams::filtering_istream is;
    if (common::isDotGzPath(annotationFilePath_))
    {
        is.push(boost::iostreams::gzip_decompressor());
    }
    is.push(boost::iostreams::file_source(annotationFilePath_.string()));

    if (!is)
    {
        const boost::format message = boost::format("Failed to open bitset file %s for reading: %s") %
            annotationFilePath_ % strerror(errno);
        BOOST_THROW_EXCEPTION(common::IoException(errno, message.str()));
    }

    std::vector<bool> neighbors(reference::genomeLength(contigsXml_.getContigs()), false);

    if (16 == bitsPerValue_)
    {
        extractNeighbors<unsigned short>(is, neighbors);
    }
    else
    {
        ISAAC_ASSERT_MSG(8 == bitsPerValue_, "Invalid bits per value, only 16 and 8 are allowed: " << bitsPerValue_);
        extractNeighbors<unsigned char>(is, neighbors);
    }

    dumpResults(neighbors);
}

template <typename ReadType>
void ExtractNeighborsFromAnnotationWorkflow::extractNeighbors(
    std::istream &bitsetFile,
    std::vector<bool> &neighbors )
{
    ReadType bits = 0;

    std::size_t position = 0;
    std::size_t neighborPositions = 0;
    while(bitsetFile.read(reinterpret_cast<char*>(&bits), sizeof(ReadType)))
    {
        const bool hasNeighbors = !!bits;
        neighbors.at(position++) = hasNeighbors;
        neighborPositions += hasNeighbors;
    }
    ISAAC_ASSERT_MSG(neighbors.size() == position, "Read only " << position << " values out of " << neighbors.size() << " expected");
    ISAAC_THREAD_CERR << "Found " << neighborPositions << " positions that have neighbors " << std::endl;
}

void ExtractNeighborsFromAnnotationWorkflow::dumpResults(const std::vector<bool> &neighbors)
{
    io::BitsetSaver neighborsSaver(neighborsFilePath_);
    neighborsSaver.save(neighbors);
    ISAAC_THREAD_CERR << "Stored " << neighbors.size() << " neighbor locations in " << neighborsFilePath_ << std::endl;
}

} // namespace workflow
} // namespace isaac
