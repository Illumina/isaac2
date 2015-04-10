/**
 ** Isaac Genome Alignment Software
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
 ** \file WigLoader.cpp
 **
 ** Loader for genomic position data stored in Wig format http://genome.ucsc.edu/goldenPath/help/wiggle.html
 **
 ** \author Roman Petrovski
 **/

#include <boost/format.hpp>

#include "common/Debug.hh"
#include "common/Exceptions.hh"
#include "io/WigLoader.hh"

namespace isaac
{
namespace io
{

WigLoader::WigLoader(
    const boost::filesystem::path &filePath) :
    filePath_(filePath),
    is_(filePath_.c_str())
{
    if (!is_)
    {
        const boost::format message = boost::format("Failed to open wig file %s for reading: %s") % filePath_ % strerror(errno);
        BOOST_THROW_EXCEPTION(common::IoException(errno, message.str()));
    }
}


} // namespace reference
} // namespace isaac
