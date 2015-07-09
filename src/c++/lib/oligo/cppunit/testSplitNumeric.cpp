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
 **/

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <boost/assign.hpp>
#include <boost/utility/binary.hpp>

using namespace std;

#include "RegistryName.hh"
#include "testSplitNumeric.hh"

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION( TestSplitNumeric, registryName("SplitNumeric"));

void TestSplitNumeric::setUp()
{
}

void TestSplitNumeric::tearDown()
{
}


void TestSplitNumeric::testAll()
{
    testUint256();
    testUint136();
    testUint160();
}

void TestSplitNumeric::testUint256()
{
    const isaac::common::Uint256 k00000002 =
                isaac::common::Uint256(
                    (__uint128_t(0x0000000000000000UL) << 64) | __uint128_t(0x0000000000000000UL),
                    (__uint128_t(0x0000000000000000UL) << 64) | __uint128_t(0x0000000022222222UL));

    const isaac::common::Uint256 k00000123 =
                isaac::common::Uint256(
                    (__uint128_t(0x0000000000000000UL) << 64) | __uint128_t(0x0000000000000000UL),
                    (__uint128_t(0x0000000011111111UL) << 64) | __uint128_t(0x2222222233333333UL));

    const isaac::common::Uint256 k00001234 =
                isaac::common::Uint256(
                    (__uint128_t(0x0000000000000000UL) << 64) | __uint128_t(0x0000000000000000UL),
                    (__uint128_t(0x1111111122222222UL) << 64) | __uint128_t(0x3333333344444444UL));

    const isaac::common::Uint256 k12345678 =
                isaac::common::Uint256(
                    (__uint128_t(0x1111111122222222UL) << 64) | __uint128_t(0x3333333344444444UL),
                    (__uint128_t(0x5555555566666666UL) << 64) | __uint128_t(0x7777777788888888UL));
    const isaac::common::Uint256 k23456780 =
                isaac::common::Uint256(
                    (__uint128_t(0x2222222233333333UL) << 64) | __uint128_t(0x4444444455555555UL),
                    (__uint128_t(0x6666666677777777UL) << 64) | __uint128_t(0x8888888800000000UL));

    const isaac::common::Uint256 k67800000 =
                isaac::common::Uint256(
                    (__uint128_t(0x6666666677777777UL) << 64) | __uint128_t(0x8888888800000000UL),
                    (__uint128_t(0x0000000000000000UL) << 64) | __uint128_t(0x0000000000000000UL));

    const isaac::common::Uint256 k80000000 =
                isaac::common::Uint256(
                    (__uint128_t(0x8888888800000000UL) << 64) | __uint128_t(0x0000000000000000UL),
                    (__uint128_t(0x0000000000000000UL) << 64) | __uint128_t(0x0000000000000000UL));

    const isaac::common::Uint256 k02345678 =
                isaac::common::Uint256(
                    (__uint128_t(0x0000000022222222UL) << 64) | __uint128_t(0x3333333344444444UL),
                    (__uint128_t(0x5555555566666666UL) << 64) | __uint128_t(0x7777777788888888UL));

    const isaac::common::Uint256 k0fffffff =
                isaac::common::Uint256(
                    (__uint128_t(0x00000000ffffffffUL) << 64) | __uint128_t(0xffffffffffffffffUL),
                    (__uint128_t(0xffffffffffffffffUL) << 64) | __uint128_t(0xffffffffffffffffUL));
    const isaac::common::Uint256 kf0000000 =
                isaac::common::Uint256(
                    (__uint128_t(0xffffffff00000000UL) << 64) | __uint128_t(0x0000000000000000UL),
                    (__uint128_t(0x0000000000000000UL) << 64) | __uint128_t(0x0000000000000000UL));

    CPPUNIT_ASSERT_EQUAL(k02345678, k23456780 >> 32);
    CPPUNIT_ASSERT_EQUAL(k23456780, k12345678 << 32);
    CPPUNIT_ASSERT_EQUAL(k02345678, k12345678 & k0fffffff);
    CPPUNIT_ASSERT_EQUAL(k00001234, k12345678 >> 128);
    CPPUNIT_ASSERT_EQUAL(k00000123, k12345678 >> 160);
    CPPUNIT_ASSERT_EQUAL(k67800000, k12345678 << 160);
    CPPUNIT_ASSERT_EQUAL(k0fffffff, ~kf0000000);
    CPPUNIT_ASSERT_EQUAL(k00000002, k23456780 >> 224);
    CPPUNIT_ASSERT_EQUAL(k80000000, k12345678 << 224);
    CPPUNIT_ASSERT_EQUAL(k12345678, k12345678 << 0);
    CPPUNIT_ASSERT_EQUAL(k12345678, k12345678 >> 0);

}


void TestSplitNumeric::testUint136()
{
    const isaac::common::Uint136 k12345 =
        isaac::common::Uint136(
            0x11, (__uint128_t(0x1111112222222233UL) << 64) | __uint128_t(0x3333334444444455UL));

    const isaac::common::Uint136 k23450 =
        isaac::common::Uint136(
            0x22, (__uint128_t(0x2222223333333344UL) << 64) | __uint128_t(0x4444445500000000UL));

    const isaac::common::Uint136 k02345 =
        isaac::common::Uint136(
            0x00, (__uint128_t(0x0000002222222233UL) << 64) | __uint128_t(0x3333334444444455UL));

    const isaac::common::Uint136 k0ffff =
        isaac::common::Uint136(
            0x00, (__uint128_t(0x000000ffffffffffUL) << 64) | __uint128_t(0xffffffffffffffffUL));

    const isaac::common::Uint136 k00001 =
        isaac::common::Uint136(
            0x00, (__uint128_t(0x0000000000000000UL) << 64) | __uint128_t(0x0000000000000011UL));

    const isaac::common::Uint136 kf0000 =
        isaac::common::Uint136(
            0xff, (__uint128_t(0xffffff0000000000UL) << 64) | __uint128_t(0x0000000000000000UL));

    const isaac::common::Uint136 k50000 =
        isaac::common::Uint136(
            0x55, (__uint128_t(0x0000000000000000UL) << 64) | __uint128_t(0x0000000000000000UL));

    CPPUNIT_ASSERT_EQUAL(k23450, k12345 << 32);
    CPPUNIT_ASSERT_EQUAL(k02345, k23450 >> 32);
    CPPUNIT_ASSERT_EQUAL(k02345, k12345 & k0ffff);
    CPPUNIT_ASSERT_EQUAL(k00001, k12345 >> 128);
    CPPUNIT_ASSERT_EQUAL(k50000, k12345 << 128);
    CPPUNIT_ASSERT_EQUAL(k0ffff, ~kf0000);
    CPPUNIT_ASSERT_EQUAL(k12345, k12345 << 0);
    CPPUNIT_ASSERT_EQUAL(k12345, k12345 >> 0);

}

void TestSplitNumeric::testUint160()
{
    {
        isaac::common::Uint160 before =
            isaac::common::Uint160(
                0x11, (__uint128_t(0x3f82374404afebffUL) << 64) | __uint128_t(0x82bc81ebadc883d4UL));

        const isaac::common::Uint160 after =
            isaac::common::Uint160(
                0x44, (__uint128_t(0xfe08dd1012bfaffeUL) << 64) | __uint128_t(0x0af207aeb7220f50UL));


        CPPUNIT_ASSERT_EQUAL(after, before << 2);

        before <<= 2;
        CPPUNIT_ASSERT_EQUAL(after, before);
    }

    {
        isaac::common::Uint160 before =
            isaac::common::Uint160(
                0x11, (__uint128_t(0x3f82374404afebffUL) << 64) | __uint128_t(0x82bc81ebadc883d5UL));

        const isaac::common::Uint160 after =
            isaac::common::Uint160(
                0x44, (__uint128_t(0xfe08dd1012bfaffeUL) << 64) | __uint128_t(0x0af207aeb7220f54UL));


        CPPUNIT_ASSERT_EQUAL(after, before << 2);

        before <<= 2;
        CPPUNIT_ASSERT_EQUAL(after, before);
    }

}
