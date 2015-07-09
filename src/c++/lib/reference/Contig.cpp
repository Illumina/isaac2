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
 ** \file Contig.
 **
 ** \brief See Contig.hh
 **
 ** \author Come Raczy
 **/

#include <numeric>
#include <boost/bind.hpp>

#include "reference/Contig.hh"

namespace isaac
{
namespace reference
{

std::size_t genomeLength(const std::vector<Contig> &contigList)
{
    return std::accumulate(
        contigList.begin(), contigList.end(),
        size_t(0), boost::bind<size_t>(std::plus<size_t>(), _1, boost::bind(&Contig::getLength, _2)));
}

} // namespace reference
} // namespace isaac
