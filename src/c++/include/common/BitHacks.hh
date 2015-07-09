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
 ** \file BitHacks.hh
 **
 ** \brief http://graphics.stanford.edu/~seander/bithacks.html
 **
 ** \author Roman Petrovski
 **/

#include <stdint.h>

#include "common/Debug.hh"

#ifndef iSAAC_COMMON_BIT_HACKS_HH
#define iSAAC_COMMON_BIT_HACKS_HH

namespace isaac
{
namespace common {


/**
 * \return next power of two value that is >= v
 */
inline unsigned long upperPowerOfTwo(unsigned long v)
{
    v--;
    v |= v >> 1;
    v |= v >> 2;
    v |= v >> 4;
    v |= v >> 8;
    v |= v >> 16;
    v |= v >> 32;
    v++;
    return v;
}

inline unsigned countBitsSet(const unsigned &ui)
{
    unsigned v = ui - ((ui >> 1) & 0x55555555);                    // reuse input as temporary
    v = (v & 0x33333333) + ((v >> 2) & 0x33333333);     // temp
    const unsigned c = (((v + (v >> 4)) & 0xF0F0F0F) * 0x1010101) >> 24; // count
    return c;
}

inline unsigned countBitsSet(const unsigned long &ul)
{
    static const int BITS_IN_BYTE = 8;
    const unsigned left = ul >> (sizeof(unsigned) * BITS_IN_BYTE);
    const unsigned right = ul;
    return countBitsSet(left) + countBitsSet(right);
}

// find the number of trailing zeros in 32-bit v
inline int lsbSet(const unsigned v)
{
    static const int MultiplyDeBruijnBitPosition[32] =
    {
      0, 1, 28, 2, 29, 14, 24, 3, 30, 22, 20, 15, 25, 17, 4, 8,
      31, 27, 13, 23, 21, 19, 16, 7, 26, 12, 18, 6, 11, 5, 10, 9
    };
    const int r = MultiplyDeBruijnBitPosition[((uint32_t)((v & -v) * 0x077CB531U)) >> 27];
    return r;
}

/**
 * \brief this function returns next higher number with same number of set bits as x.
 * http://www.geeksforgeeks.org/next-higher-number-with-same-number-of-set-bits/
 */
template <typename T>
inline T snoob(const T &x)
{
    T rightOne;
    T nextHigherOneBit;
    T rightOnesPattern;

    T next = 0;

    if(x)
    {
        // right most set bit
        rightOne = x & -x;

        // reset the pattern and set next higher bit
        // left part of x will be here
        nextHigherOneBit = x + rightOne;

        // nextHigherOneBit is now part [D] of the above explanation.

        // isolate the pattern
        rightOnesPattern = x ^ nextHigherOneBit;

        // right adjust pattern
        rightOnesPattern = (rightOnesPattern)/rightOne;

        // correction factor
        rightOnesPattern >>= 2;

        // rightOnesPattern is now part [A] of the above explanation.

        // integrate new pattern (Add [D] and [A])
        next = nextHigherOneBit | rightOnesPattern;
    }

    ISAAC_ASSERT_MSG(common::countBitsSet(x) == common::countBitsSet(next), "Incorrect number of bits set in " << next << " original: " << x);

    return next;
}

} //namespace common
} //namespace isaac

#endif //iSAAC_COMMON_BIT_HACKS_HH
