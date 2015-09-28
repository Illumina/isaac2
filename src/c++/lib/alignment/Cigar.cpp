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
 ** \file Cigar.cpp
 **
 ** \brief See Cigar.hh
 ** 
 ** \author Come Raczy
 **/

#include "alignment/Cigar.hh"
#include "common/Debug.hh"

namespace isaac
{
namespace alignment
{

std::string Cigar::toString() const
{
    return toString(begin(), end());
}

std::string Cigar::toString(const unsigned offset, const unsigned length) const
{
    ISAAC_ASSERT_MSG(this->size() >= offset + length, "Requested end is outside of cigarBuffer");
    return toString(begin() + offset, begin() + offset + length);
}

std::string Cigar::toString(const Cigar &cigarBuffer, unsigned offset, unsigned length)
{
    ISAAC_ASSERT_MSG(cigarBuffer.size() >= offset + length, "Requested end is outside of cigarBuffer");
    return toString(cigarBuffer.begin() + offset, cigarBuffer.begin() + offset + length);
}

} // namespace alignment
} // namespace isaac
