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
 ** \file InstallationPaths.hh
 **
 ** \brief Path resolution for installed components
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_PACKAGE_INSTALLATION_PATHS_HH
#define iSAAC_PACKAGE_INSTALLATION_PATHS_HH

#include <boost/filesystem.hpp>

namespace isaac
{
namespace package
{

void initialize(const boost::filesystem::path &modulePath, const char *homeOffset);
boost::filesystem::path expandPath(const char *path);

} // namespace reference
} // namespace isaac

#endif // #ifndef iSAAC_PACKAGE_INSTALLATION_PATHS_HH
