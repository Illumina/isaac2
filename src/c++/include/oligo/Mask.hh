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
 ** \file Mask.hh
 **
 ** Masking operations on oligomers.
 **
 ** \author Come Raczy
 **/

#ifndef iSAAC_OLIGO_MASK_HH
#define iSAAC_OLIGO_MASK_HH

#include <string>
#include <boost/format.hpp>
#include <boost/filesystem.hpp>
#include <boost/numeric/conversion/cast.hpp>

#include "oligo/Kmer.hh"

namespace isaac
{
namespace oligo
{

inline unsigned int getMaskCount(const unsigned int maskWidth) {return (1 << maskWidth);}

} // namespace oligo
} // namespace isaac

#endif // #ifndef iSAAC_OLIGO_MASK_HH
