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
 ** \file SeedDescriptorOption.hh
 **
 ** seeds option parsing
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_OPTIONS_SEED_DESCRIPTOR_OPTION_HH
#define iSAAC_OPTIONS_SEED_DESCRIPTOR_OPTION_HH

#include "flowcell/ReadMetadata.hh"

namespace isaac
{
namespace options
{
namespace alignOptions
{

alignment::SeedMetadataList parseSeedDescriptor(
    const std::vector<flowcell::ReadMetadata> &readMetadataList,
    const std::string &seedDescriptor,
    const unsigned seedLength,
    unsigned &firstPassSeeds);

} // namespace alignOptions
} // namespace options
} // namespace isaac

#endif // #ifndef iSAAC_OPTIONS_SEED_DESCRIPTOR_OPTION_HH
