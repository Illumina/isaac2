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
#include <boost/foreach.hpp>
#include <boost/mpl/bool.hpp>
#include <boost/math/tools/big_constant.hpp>
using namespace std;

#include "RegistryName.hh"
#include "testPermutate.hh"

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION( TestPermutate, registryName("Permutate"));

void TestPermutate::setUp()
{
}

void TestPermutate::tearDown()
{
}

void TestPermutate::testFourBlocks()
{
    using namespace isaac;
    using boost::assign::list_of;
    const oligo::KmerType kmer(0xFEDCBA9876543210UL);
    const oligo::KmerType abcd(0xFEDCBA9876543210UL);
    const oligo::KmerType adbc(0xFEDC3210BA987654UL);
    const oligo::KmerType dbca(0x3210BA987654FEDCUL);
    const unsigned blockLength = 8;
    const std::vector<unsigned char> ABCD = list_of(0)(1)(2)(3);
    const std::vector<unsigned char> ADBC = list_of(0)(3)(1)(2);
    const std::vector<unsigned char> DBCA = list_of(3)(1)(2)(0);
    const oligo::Permutate ABCD_ABCD(blockLength, ABCD, ABCD, 0);
    const oligo::Permutate ABCD_ADBC(blockLength, ABCD, ADBC, 0);
    const oligo::Permutate ABCD_DBCA(blockLength, ABCD, DBCA, 0);
    const oligo::Permutate ADBC_DBCA(blockLength, ADBC, DBCA, 0);
    CPPUNIT_ASSERT_EQUAL(ABCD_ABCD(kmer), kmer);
    CPPUNIT_ASSERT_EQUAL(ABCD_ADBC(kmer), oligo::KmerType(0xFEDC3210BA987654UL));
    CPPUNIT_ASSERT_EQUAL(ABCD_DBCA(kmer), oligo::KmerType(0x3210BA987654FEDCUL));
    CPPUNIT_ASSERT_EQUAL(ADBC_DBCA(kmer), oligo::KmerType(0xBA9876543210FEDCUL));
    CPPUNIT_ASSERT_EQUAL(ABCD_ABCD(abcd), abcd);
    CPPUNIT_ASSERT_EQUAL(ABCD_ADBC(abcd), adbc);
    CPPUNIT_ASSERT_EQUAL(ABCD_DBCA(abcd), dbca);
    CPPUNIT_ASSERT_EQUAL(ADBC_DBCA(adbc), dbca);
    CPPUNIT_ASSERT_EQUAL(ABCD_ABCD.reorder(abcd), abcd);
    CPPUNIT_ASSERT_EQUAL(ABCD_ADBC.reorder(adbc), abcd);
    CPPUNIT_ASSERT_EQUAL(ABCD_DBCA.reorder(dbca), abcd);
    CPPUNIT_ASSERT_EQUAL(ADBC_DBCA.reorder(dbca), abcd);
}

void TestPermutate::testEightBlocks()
{
    using namespace isaac;
    using boost::assign::list_of;
    // A  B  C  D  E  F  G  H
    //FE DC BA 98 76 54 32 10
    const oligo::KmerType kmer = oligo::KmerType(0xFEDCBA9876543210UL);
    const oligo::KmerType abcdefgh = oligo::KmerType(0xFEDCBA9876543210UL);
    const oligo::KmerType adbcefgh = oligo::KmerType(0xFE98DCBA76543210UL);
    const oligo::KmerType dghbcafe = oligo::KmerType(0x983210DCBAFE5476UL);
    const unsigned blockLength = 4;
    const std::vector<unsigned char> ABCDEFGH = list_of(0)(1)(2)(3)(4)(5)(6)(7);
    const std::vector<unsigned char> ADBCEFGH = list_of(0)(3)(1)(2)(4)(5)(6)(7);
    const std::vector<unsigned char> DGHBCAFE = list_of(3)(6)(7)(1)(2)(0)(5)(4);
    const oligo::Permutate ABCDEFGH_ABCDEFGH(blockLength, ABCDEFGH, ABCDEFGH, 0);
    const oligo::Permutate ABCDEFGH_ADBCEFGH(blockLength, ABCDEFGH, ADBCEFGH, 0);
    const oligo::Permutate ABCDEFGH_DGHBCAFE(blockLength, ABCDEFGH, DGHBCAFE, 0);
    const oligo::Permutate ADBCEFGH_DGHBCAFE(blockLength, ADBCEFGH, DGHBCAFE, 0);
    CPPUNIT_ASSERT_EQUAL(ABCDEFGH_ABCDEFGH(kmer), kmer);
    CPPUNIT_ASSERT_EQUAL(ABCDEFGH_ADBCEFGH(kmer), oligo::KmerType(0xFE98DCBA76543210UL));
    CPPUNIT_ASSERT_EQUAL(ABCDEFGH_DGHBCAFE(kmer), oligo::KmerType(0x983210DCBAFE5476UL));
    CPPUNIT_ASSERT_EQUAL(ADBCEFGH_DGHBCAFE(kmer), oligo::KmerType(0xDC3210BA98FE5476UL));
    CPPUNIT_ASSERT_EQUAL(ABCDEFGH_ABCDEFGH(abcdefgh), abcdefgh);
    CPPUNIT_ASSERT_EQUAL(ABCDEFGH_ADBCEFGH(abcdefgh), adbcefgh);
    CPPUNIT_ASSERT_EQUAL(ABCDEFGH_DGHBCAFE(abcdefgh), dghbcafe);
    CPPUNIT_ASSERT_EQUAL(ADBCEFGH_DGHBCAFE(adbcefgh), dghbcafe);
    CPPUNIT_ASSERT_EQUAL(ABCDEFGH_ABCDEFGH.reorder(abcdefgh), abcdefgh);
    CPPUNIT_ASSERT_EQUAL(ABCDEFGH_ADBCEFGH.reorder(adbcefgh), abcdefgh);
    CPPUNIT_ASSERT_EQUAL(ABCDEFGH_DGHBCAFE.reorder(dghbcafe), abcdefgh);
    CPPUNIT_ASSERT_EQUAL(ADBCEFGH_DGHBCAFE.reorder(dghbcafe), abcdefgh);
}

template <typename KmerT>
void testPermutate(KmerT original, const KmerT expected, const std::vector<isaac::oligo::Permutate> permutateList)
{
    CPPUNIT_ASSERT_EQUAL(original, permutateList.front()(original));
    CPPUNIT_ASSERT_EQUAL(original, permutateList.front().reorder(original));
    KmerT permuted = original;
    BOOST_FOREACH(const isaac::oligo::Permutate &permutate, permutateList)
    {
        permuted = permutate(permuted);
        CPPUNIT_ASSERT_EQUAL(original, permutate.reorder(permuted));
    }
    CPPUNIT_ASSERT_EQUAL(expected, permuted);
    CPPUNIT_ASSERT_EQUAL(original, permutateList.back().reorder(permuted));
}

// This is needed for cases when isaac::oligo::Kmer is defined as __uint128_t or else CPPUNIT_ASSERT_EQUAL fails to compile
inline std::ostream & operator <<(std::ostream &os, const isaac::oligo::LongKmerType &kmer)
{
    return os << isaac::oligo::bases(kmer);
}


void TestPermutate::testNoErrors()
{
    const isaac::oligo::ShortKmerType ORIGINAL16(0x76543210U);
    const isaac::oligo::KmerType ORIGINAL(0xFEDCBA9876543210UL);
    const isaac::oligo::LongKmerType ORIGINAL64(isaac::oligo::LongKmerType(0x1111222233334444UL) << 64 | isaac::oligo::LongKmerType(0x5555666677778888UL));

    using namespace isaac;
    {
        const std::vector<oligo::Permutate> permutateList = oligo::getPermutateList<oligo::ShortKmerType>(0, true, false);
        CPPUNIT_ASSERT_EQUAL(1UL, permutateList.size());
        testPermutate(ORIGINAL16, ORIGINAL16, permutateList);
    }

    {
        const std::vector<oligo::Permutate> permutateList = oligo::getPermutateList<oligo::KmerType>(0, true, false);
        CPPUNIT_ASSERT_EQUAL(1UL, permutateList.size());
        testPermutate(ORIGINAL, ORIGINAL, permutateList);
    }

    {
        const std::vector<oligo::Permutate> permutateList = oligo::getPermutateList<oligo::LongKmerType>(0, true, false);
        CPPUNIT_ASSERT_EQUAL(1UL, permutateList.size());
        testPermutate(ORIGINAL64, ORIGINAL64, permutateList);
    }
}


void TestPermutate::testTwoErrors()
{
    const isaac::oligo::ShortKmerType ORIGINAL16(0x76543210U);
    const isaac::oligo::ShortKmerType EXPECTED16(0x10547632U);
    const isaac::oligo::KmerType ORIGINAL(0xFEDCBA9876543210UL);
    const isaac::oligo::KmerType EXPECTED(0x3210ba98fedc7654UL);
    const isaac::oligo::BasicKmerType<48> ORIGINAL48(isaac::oligo::BasicKmerType<48>(0x33334444UL) << 64 | isaac::oligo::BasicKmerType<48>(0x5555666677778888UL));
    const isaac::oligo::BasicKmerType<48> EXPECTED48(isaac::oligo::BasicKmerType<48>(0x77888844UL) << 64 | isaac::oligo::BasicKmerType<48>(0x5555333344666677UL));
    const isaac::oligo::LongKmerType ORIGINAL64(isaac::oligo::LongKmerType(0x1111222233334444UL) << 64 | isaac::oligo::LongKmerType(0x5555666677778888UL));
    const isaac::oligo::LongKmerType EXPECTED64(isaac::oligo::LongKmerType(0x7777888833334444UL) << 64 | isaac::oligo::LongKmerType(0x1111222255556666UL));

    using namespace isaac;
    {
        const std::vector<oligo::Permutate> permutateList = oligo::getPermutateList<oligo::ShortKmerType>(2, 0, false);
        CPPUNIT_ASSERT_EQUAL(6UL, permutateList.size());
        testPermutate(ORIGINAL16, EXPECTED16, permutateList);
    }

    {
        const std::vector<oligo::Permutate> permutateList = oligo::getPermutateList<oligo::KmerType>(2, 0, false);
        CPPUNIT_ASSERT_EQUAL(6UL, permutateList.size());
        testPermutate(ORIGINAL, EXPECTED, permutateList);
    }

    {
        const std::vector<oligo::Permutate> permutateList = oligo::getPermutateList<oligo::BasicKmerType<48> >(2, 0, false);
        CPPUNIT_ASSERT_EQUAL(6UL, permutateList.size());
        testPermutate(ORIGINAL48, EXPECTED48, permutateList);
    }

    {
        const std::vector<oligo::Permutate> permutateList = oligo::getPermutateList<oligo::LongKmerType>(2, 0, false);
        CPPUNIT_ASSERT_EQUAL(6UL, permutateList.size());
        testPermutate(ORIGINAL64, EXPECTED64, permutateList);
    }
}

void TestPermutate::testThreeErrors()
{
    const isaac::oligo::ShortKmerType ABCDEFGH16(0x76543210U);
    const isaac::oligo::ShortKmerType FGHDEABC16(0x01543762U);
    const isaac::oligo::KmerType ABCDEFGH(0xFEDCBA9876543210UL);
    const isaac::oligo::KmerType FGHDEABC(0x1032ba9876fedc54UL);
    const isaac::oligo::BasicKmerType<48> ABCDEFGH48(isaac::oligo::BasicKmerType<48>(0x00011122UL) << 64 | isaac::oligo::BasicKmerType<48>(0x2333444555666777UL));
    const isaac::oligo::BasicKmerType<48> FGHDEABC48(isaac::oligo::BasicKmerType<48>(0x77766622UL) << 64 | isaac::oligo::BasicKmerType<48>(0x2333444000111555UL));
    const isaac::oligo::LongKmerType ABCDEFGH64(isaac::oligo::LongKmerType(0x1111222233334444UL) << 64 | isaac::oligo::LongKmerType(0x5555666677778888UL));
    const isaac::oligo::LongKmerType FGHDEABC64(isaac::oligo::LongKmerType(0x8888777733334444UL) << 64 | isaac::oligo::LongKmerType(0x5555111122226666UL));

    using namespace isaac;
    {
        const std::vector<oligo::Permutate> permutateList = oligo::getPermutateList<oligo::ShortKmerType>(3, 0, false);
        CPPUNIT_ASSERT_EQUAL(56UL, permutateList.size());
        testPermutate(ABCDEFGH16, FGHDEABC16, permutateList);
    }

    {
        const std::vector<oligo::Permutate> permutateList = oligo::getPermutateList<oligo::KmerType>(3, 0, false);
        CPPUNIT_ASSERT_EQUAL(56UL, permutateList.size());
        testPermutate(ABCDEFGH, FGHDEABC, permutateList);
    }

    {
        const std::vector<oligo::Permutate> permutateList = oligo::getPermutateList<oligo::BasicKmerType<48> >(3, 0, false);
        CPPUNIT_ASSERT_EQUAL(56UL, permutateList.size());
        testPermutate(ABCDEFGH48, FGHDEABC48, permutateList);
    }

    {
        const std::vector<oligo::Permutate> permutateList = oligo::getPermutateList<oligo::LongKmerType>(3, 0, false);
        CPPUNIT_ASSERT_EQUAL(56UL, permutateList.size());
        testPermutate(ABCDEFGH64, FGHDEABC64, permutateList);
    }
}


void TestPermutate::testFourErrors()
{
    using namespace isaac;
    {
        const isaac::oligo::ShortKmerType ORIGINAL16(0x76543210U);
        const isaac::oligo::ShortKmerType EXPECTED16(0x02147653U);
        const std::vector<oligo::Permutate> permutateList = oligo::getPermutateList<oligo::ShortKmerType>(4, 0, false);
        CPPUNIT_ASSERT_EQUAL(70UL, permutateList.size());
        testPermutate(ORIGINAL16, EXPECTED16, permutateList);
    }

    {
        const isaac::oligo::KmerType ORIGINAL(0xFEDCBA9876543210UL);
        const isaac::oligo::KmerType EXPECTED(0x10543298fedcba76UL);
        const std::vector<oligo::Permutate> permutateList = oligo::getPermutateList<oligo::KmerType>(4, 0, false);
        CPPUNIT_ASSERT_EQUAL(70UL, permutateList.size());
        testPermutate(ORIGINAL, EXPECTED, permutateList);
    }

    {
        const isaac::oligo::BasicKmerType<48> ORIGINAL48(isaac::oligo::BasicKmerType<48>(0x33334444UL) << 64 | isaac::oligo::BasicKmerType<48>(0x5555666677778888UL));
        const isaac::oligo::BasicKmerType<48> EXPECTED48(isaac::oligo::BasicKmerType<48>(0x88867777UL) << 64 | isaac::oligo::BasicKmerType<48>(0x8555333344445666UL));
        const std::vector<oligo::Permutate> permutateList = oligo::getPermutateList<oligo::BasicKmerType<48> >(4, 0, false);
        CPPUNIT_ASSERT_EQUAL(70UL, permutateList.size());
        testPermutate(ORIGINAL48, EXPECTED48, permutateList);
    }

    {
        const isaac::oligo::LongKmerType ORIGINAL64(isaac::oligo::LongKmerType(0x1111222233334444UL) << 64 | isaac::oligo::LongKmerType(0x5555666677778888UL));
        const isaac::oligo::LongKmerType EXPECTED64(isaac::oligo::LongKmerType(0x8888666677774444UL) << 64 | isaac::oligo::LongKmerType(0x1111222233335555UL));
        const std::vector<oligo::Permutate> permutateList = oligo::getPermutateList<oligo::LongKmerType>(4, 0, false);
        CPPUNIT_ASSERT_EQUAL(70UL, permutateList.size());
        testPermutate(ORIGINAL64, EXPECTED64, permutateList);
    }

    {
        const isaac::oligo::BasicKmerType<128> ORIGINAL128(
            isaac::common::Uint256(
                __uint128_t(0x1111111122222222UL) << 64 | __uint128_t(0x3333333344444444UL),
                __uint128_t(0x5555555566666666UL) << 64 | __uint128_t(0x7777777788888888UL)));
        const isaac::oligo::BasicKmerType<128> EXPECTED128(
            isaac::common::Uint256(
                __uint128_t(0x8888888866666666UL) << 64 | __uint128_t(0x7777777744444444UL),
                __uint128_t(0x1111111122222222UL) << 64 | __uint128_t(0x3333333355555555UL)));

        const std::vector<oligo::Permutate> permutateList = oligo::getPermutateList<isaac::oligo::BasicKmerType<128> >(4, 0, false);
        CPPUNIT_ASSERT_EQUAL(70UL, permutateList.size());
        testPermutate(ORIGINAL128, EXPECTED128, permutateList);
    }

    {
        const isaac::oligo::BasicKmerType<72> ORIGINAL72(
            isaac::common::Uint144(
                0x1111, __uint128_t(0x1111222222223333UL) << 64 | __uint128_t(0x3333444444445555UL)));
        const isaac::oligo::BasicKmerType<72> EXPECTED72(
            isaac::common::Uint144(
                0x1555, __uint128_t(0x7344444446223311UL) << 64 | __uint128_t(0x111111222220ccccUL)));

        const std::vector<oligo::Permutate> permutateList = oligo::getPermutateList<isaac::oligo::BasicKmerType<72> >(4, 0, false);
        CPPUNIT_ASSERT_EQUAL(70UL, permutateList.size());
        testPermutate(ORIGINAL72, EXPECTED72, permutateList);
    }

}
