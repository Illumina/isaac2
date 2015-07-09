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
 ** \file AnnotationLoader.cpp
 **
 ** Helper utility for loading multiple contigs of a fasta file.
 **
 ** \author Roman Petrovski
 **/
#include <errno.h>

#include <boost/format.hpp>

#include "common/Exceptions.hh"
#include "reference/AnnotationLoader.hh"

namespace isaac
{
namespace reference
{

/**
 * \brief split single reference annotation blob into per-contig annotations
 */
ContigAnnotations loadContigAnnotations(const reference::SortedReferenceMetadata &sortedReferenceMetadata)
{
    const boost::filesystem::path &path = sortedReferenceMetadata.getKUniquenessAnnotation().path_;

    boost::iostreams::filtering_istream is;
    if (common::isDotGzPath(path))
    {
        is.push(boost::iostreams::gzip_decompressor());
    }
    is.push(boost::iostreams::file_source(path.string()));

    if (!is)
    {
        BOOST_THROW_EXCEPTION(common::IoException(errno, "Failed to open annotation file " + path.string()));
    }

    const SortedReferenceMetadata::Contigs contigs = sortedReferenceMetadata.getKaryotypeOrderedContigs();
    //const std::vector<unsigned long> contigOffsets = reference::computeContigOffsets(sortedReferenceMetadataList_.at(0).getKaryotypeOrderedContigs());

    ContigAnnotations ret(contigs.size());//(genomeLength(sortedReferenceMetadata.getContigs()));
    BOOST_FOREACH(const SortedReferenceMetadata::Contig &contig, contigs)
    {
        ContigAnnotation &contigAnnotation = ret.at(contig.karyotypeIndex_);
        contigAnnotation.resize(contig.totalBases_);
        std::size_t currentOffset = 0;
        if (!is.read(reinterpret_cast<char *>(&contigAnnotation.front()), contigAnnotation.size() * sizeof(AnnotationValue)) ||
            contigAnnotation.size() * sizeof(AnnotationValue) != std::size_t(is.gcount()))
        {
            const boost::format message = boost::format("Failed to read %d bytes at offset %d from annotation file %s") %
                (contigAnnotation.size() * sizeof(AnnotationValue)) % currentOffset % path.string();
            BOOST_THROW_EXCEPTION(common::IoException(errno, message.str()));
        }
        currentOffset += contigAnnotation.size();
    }

    return ret;
}

/**
 * \brief load contig annotations for each reference in the list
 */
ContigAnnotationsList loadAnnotations(const reference::SortedReferenceMetadataList &sortedReferenceMetadataList)
{
    ISAAC_TRACE_STAT("loadAnnotations ");

    ContigAnnotationsList ret(sortedReferenceMetadataList.size());
    std::size_t referenceIndex = 0;
    BOOST_FOREACH(const reference::SortedReferenceMetadata &ref, sortedReferenceMetadataList)
    {
        if (ref.hasKUniquenessAnnotation())
        {
            ret.at(referenceIndex) = reference::loadContigAnnotations(ref);
        }
        else
        {
            ISAAC_THREAD_CERR << "WARNING: No annotation is available for reference genome. Alignment scores could be misleading. " <<
                ref.getContigs().front().filePath_ << std::endl;
        }
        ++referenceIndex;
    }

    ISAAC_TRACE_STAT("loadAnnotations done ");

    return ret;
}

} // namespace reference
} // namespace isaac
