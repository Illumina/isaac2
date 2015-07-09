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
 ** \file SortedReferenceXml.hh
 **
 ** sorted-reference.xml io helper.
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_REFERENCE_SORTED_REFERENCE_XML_HH
#define iSAAC_REFERENCE_SORTED_REFERENCE_XML_HH

#include <boost/filesystem.hpp>

#include "reference/SortedReferenceMetadata.hh"

namespace isaac
{
namespace reference
{

SortedReferenceMetadata loadSortedReferenceXml(
    std::istream &is);

SortedReferenceMetadata loadSortedReferenceXml(
    const boost::filesystem::path &xmlPath,
    const bool makeAbsolutePaths = true);

void saveSortedReferenceXml(
    std::ostream &os,
    const SortedReferenceMetadata &sortedReferenceMetadata);

void saveSortedReferenceXml(
    const boost::filesystem::path &xmlPath,
    const SortedReferenceMetadata &sortedReferenceMetadata);

} // namespace reference
} // namespace isaac

#endif // #ifndef iSAAC_REFERENCE_SORTED_REFERENCE_XML_HH
