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
 ** \file ContigsLoader.hh
 **
 ** Helper utility for loading multiple contigs of a fasta file.
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_REFERENCE_ANNOTATION_LOADER_HH
#define iSAAC_REFERENCE_ANNOTATION_LOADER_HH

#include <boost/filesystem.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>

#include "common/Exceptions.hh"
#include "common/FileSystem.hh"
#include "reference/KUniqueness.hh"
#include "reference/SortedReferenceMetadata.hh"

namespace isaac
{
namespace reference
{

/**
 * \brief load annotation file as a single blob
 */
template <typename AnnotationType>
AnnotationType loadAnnotationBlob(const boost::filesystem::path &path, const std::size_t genomeLength)
{
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

    const std::size_t annotationBytes = genomeLength * sizeof(typename AnnotationType::value_type);
    AnnotationType ret(genomeLength);
    if (!is.read(reinterpret_cast<char *>(&ret.front()), annotationBytes) ||
        annotationBytes != std::size_t(is.gcount()))
    {
        BOOST_THROW_EXCEPTION(common::IoException(
            errno, "Failed to read " +
            boost::lexical_cast<std::string>(annotationBytes) + " bytes from annotation file " + path.string()));
    }

    return ret;
}

ContigAnnotationsList loadAnnotations(const reference::SortedReferenceMetadataList &sortedReferenceMetadataList);


} // namespace reference
} // namespace isaac

#endif // #ifndef iSAAC_REFERENCE_ANNOTATION_LOADER_HH
