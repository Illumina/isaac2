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
 ** \file FastIo.hh
 **
 ** \brief Fast IO routines for integers and fixed width floating points.
 **
 ** \author Come Raczy
 **/

#ifndef iSAAC_COMMON_ENDIANNESS_HH
#define iSAAC_COMMON_ENDIANNESS_HH


namespace isaac
{
namespace common
{

template <typename T, typename IterT> IterT extractLittleEndian(const IterT p, T &result)
{
#ifdef ISAAC_IS_BIG_ENDIAN
    // TODO: untested
    BOOST_STATIC_ASSERT(sizeof(T) / 8 == 4);
    // this should work on any endian machine but if we're on little-endian already the alternative is probably faster
    result = T(*p) + T(*(p + 1)) * 256 + T(*(p + 2)) * 256 * 256 + T(*(p + 3)) * 256 * 256 * 256;
#else
    result = *reinterpret_cast<const T*>(&*p);
#endif
    return p + sizeof(T);
}

template <typename T, typename IterT> T extractLittleEndian(const IterT p)
{
#ifdef ISAAC_IS_BIG_ENDIAN
    // TODO: untested
    BOOST_STATIC_ASSERT(sizeof(T) / 8 == 4);
    // this should work on any endian machine but if we're on little-endian already the alternative is probably faster
    return T(*p) + T(*(p + 1)) * 256 + T(*(p + 2)) * 256 * 256 + T(*(p + 3)) * 256 * 256 * 256;
#endif
    return *reinterpret_cast<const T*>(&*p);
}

} // namespace common
} // namespace isaac

#endif // #ifndef iSAAC_COMMON_ENDIANNESS_HH
