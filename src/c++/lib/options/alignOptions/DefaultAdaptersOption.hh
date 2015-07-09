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
 ** \file DefaultAdaptersOption.hh
 **
 ** default-adapters option parsing
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_OPTIONS_DEFAULT_ADAPTERS_OPTION_HH
#define iSAAC_OPTIONS_DEFAULT_ADAPTERS_OPTION_HH

#include "flowcell/ReadMetadata.hh"

namespace isaac
{
namespace options
{

const flowcell::SequencingAdapterMetadataList parseDefaultAdapters (const std::string &defaultAdapters);

} // namespace options
} // namespace isaac

#endif // #ifndef iSAAC_OPTIONS_DEFAULT_ADAPTERS_OPTION_HH
