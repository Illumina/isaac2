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
#include <algorithm>
#include <cstdlib>
#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <boost/assign.hpp>
#include <boost/assign/std/vector.hpp> 

using namespace std;

#include "RegistryName.hh"
#include "testBandedSmithWaterman.hh"
#include "BuilderInit.hh"
#include "alignment/Cigar.hh"

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION( TestBandedSmithWaterman, registryName("BandedSmithWaterman"));

std::string getGenome(unsigned size = 1000)
{
    static const std::string bases = "ACGT";
    std::string genome;
    genome.reserve(size);
    unsigned int seed = 4;
    while (size--)
    {
        genome.push_back(bases[rand_r(&seed) % bases.size()]);
    }
    return genome;
}

TestBandedSmithWaterman::TestBandedSmithWaterman()
    : bsw(2, -1, 15, 3, 300)
//    : bsw(2, -1, 15, 3, 300)
    , genome(getGenome())
{
}

void TestBandedSmithWaterman::setUp()
{
}

void TestBandedSmithWaterman::tearDown()
{
}


void TestBandedSmithWaterman::testCustom()
{
    {
        const std::vector<char> query =    vectorFromString("        TTTTTTTTTTTTTTTCCCC            AATAAAAGTAAATGCTGTAATATTTAAGGTAAGTAAGCCAACAGCTGAGGAAAAGGGAAATCCATTAT", true);
        const std::vector<char> database = vectorFromString("     GGGGGCAATAACCTAAGATAAAAAAAACATTGTAAACAAAAGTAAATGCTGTAATATTTAAGGTAAGTAAGCCAACAGCTGAGGAAAAGGGAAATCCATTAT  ", true);
        isaac::alignment::Cigar cigar; cigar.reserve(1024);
        isaac::alignment::BandedSmithWaterman bsw(0, -4, 6, 1, 300);
        // earlier versions would introduce a soft-clip at the beginning. This has been corrected with SAAC-770
        // CPPUNIT_ASSERT_EQUAL(29U, bsw.align(query, database.begin(), database.end(), cigar));
        // CPPUNIT_ASSERT_EQUAL(std::string("14S73M"), isaac::alignment::Cigar::toString(cigar.begin(), cigar.end()));
        CPPUNIT_ASSERT_EQUAL(15U, bsw.align(query, database.begin(), database.end(), cigar));
        CPPUNIT_ASSERT_EQUAL(std::string("87M"), isaac::alignment::Cigar::toString(cigar.begin(), cigar.end()));
    }

    {
        const std::vector<char> query =    vectorFromString("        GGCAATAACCTAAGAT               AATAAAAGTAAATGCTGTAATATTTAAGGTAAGTAAGCCAACAGCTGAGGAAAAGGGAAATCCATTAT", true);
        const std::vector<char> database = vectorFromString("        GGCAATAACCTAAGATAAAAAAAACATTGTAAACAAAAGTAAATGCTGTAATATTTAAGGTAAGTAAGCCAACAGCTGAGGAAAAGGGAAATCCATTAT  ", true);
        isaac::alignment::Cigar cigar; cigar.reserve(1024);
        isaac::alignment::BandedSmithWaterman bsw(0, -4, 6, 1, 300);
        bsw.align(query, database.begin(), database.end(), cigar);
        CPPUNIT_ASSERT_EQUAL(std::string("16M15D68M"), isaac::alignment::Cigar::toString(cigar.begin(), cigar.end()));
    }

    {
        const std::vector<char> query = vectorFromString(           "GGCAATAACCTAAGATAAAAAAAAACATTGTAAATAAAAGTAAATGCTGTAATATTTAAGGTAAGTAAGCCAACAGCTGAGGAAAAGGGAAATCCATTAT");
        const std::vector<char> database = vectorFromString("CTTTAATCAGGCAATAACCTAAGATAAAAAAAACATTGTAAACAAAAGTAAATGCTGTAATATTTAAGGTAAGTAAGCCAACAGCTGAGGAAAAGGGAAATCCATTATCCTGGGT");
        isaac::alignment::Cigar cigar; cigar.reserve(1024);
        bsw.align(query, database.begin(), database.end(), cigar);
        CPPUNIT_ASSERT_EQUAL(std::string("16M1I83M"), isaac::alignment::Cigar::toString(cigar.begin(), cigar.end()));
    }

    {
        const std::vector<char> query =    vectorFromString(       "AGGCAATAACCTAAGATAAAAAAAACATTGTAAACAAAAGTAAATGCTGTAATATTTAAGGTAAGTAAGCCAACAGCTGAGGAAAAGGGAAATCCATTA");
        const std::vector<char> database = vectorFromString("NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNAATATTTAAGGTAAGTAAGCCAACAGCTGAGGAAAAGGGAAATCCATTATCCTGGGT");
        isaac::alignment::Cigar cigar; cigar.reserve(1024);
        isaac::alignment::BandedSmithWaterman bsw(0, -4, 6, 1, 300);
        // earlier versions would introduce a soft-clip at the beginning. This has been corrected with SAAC-770
        //CPPUNIT_ASSERT_EQUAL(std::string("50S49M"), isaac::alignment::Cigar::toString(cigar.begin(), cigar.end()));
        CPPUNIT_ASSERT_EQUAL(7U, bsw.align(query, database.begin(), database.end(), cigar));
        CPPUNIT_ASSERT_EQUAL(std::string("99M"), isaac::alignment::Cigar::toString(cigar.begin(), cigar.end()));
    }

    {
        const std::vector<char> query =    vectorFromString(       "NNNNNNNNNNNNNNNNNNNNAAAAACATTGTAAACAAAAGTAAATGCTGTAATNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN");
        const std::vector<char> database = vectorFromString("TTTAATCAGGCAATAACCTAAGATAAAAAAAACATTGTAAACAAAAGTAAATGCTGTAATATTTAAGGTAAGTAAGCCAACAGCTGAGGAAAAGGGAAATCCATTATCCTGGGT");
        isaac::alignment::Cigar cigar; cigar.reserve(1024);
        isaac::alignment::BandedSmithWaterman bsw(1, -4, 6, 1, 300);
        // earlier versions would introduce a soft-clip at the beginning. This has been corrected with SAAC-770
        // CPPUNIT_ASSERT_EQUAL(27U, bsw.align(query, database.begin(), database.end(), cigar));
        // CPPUNIT_ASSERT_EQUAL(std::string("20S33M7I39M"), isaac::alignment::Cigar::toString(cigar.begin(), cigar.end()));
        CPPUNIT_ASSERT_EQUAL(7U, bsw.align(query, database.begin(), database.end(), cigar));
        CPPUNIT_ASSERT_EQUAL(std::string("53M7I39M"), isaac::alignment::Cigar::toString(cigar.begin(), cigar.end()));
    }

    {
        const std::vector<char> query =    vectorFromString(       "TTCGTTGGAAACGGGATTTCTTCATATTCTGCTAGACAGAAGAATTCTCAGTAACTTCCTTGTGTTGTGTGTATTCAACTCACAGAGTTGAACGATCCTT");
        const std::vector<char> database = vectorFromString("TTGAGGCCTTCGTTGGAAACGGGATTTCTTCATATTCTGCTAGACAGAAGAATTCNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN");

        isaac::alignment::Cigar cigar; cigar.reserve(1024);
        isaac::alignment::BandedSmithWaterman bsw(1, -4, 6, 1, 300);
        CPPUNIT_ASSERT_EQUAL(8U, bsw.align(query, database.begin(), database.end(), cigar));
        // earlier versions would introduce a soft-clip at the beginning. This has been corrected with SAAC-770
        // CPPUNIT_ASSERT_EQUAL(std::string("46M8I1M45S"), isaac::alignment::Cigar::toString(cigar.begin(), cigar.end()));
        CPPUNIT_ASSERT_EQUAL(std::string("46M8I46M"), isaac::alignment::Cigar::toString(cigar.begin(), cigar.end()));
    }

}

void TestBandedSmithWaterman::testUngapped()
{
    const std::vector<char> database = subv(genome, 100, 115);
    std::vector<char> query = subv(database, 0, 100);
    CPPUNIT_ASSERT_EQUAL(115UL, database.size());
    isaac::alignment::Cigar cigar; cigar.reserve(1024);
    bsw.align(query, database.begin(), database.end(), cigar);
    CPPUNIT_ASSERT_EQUAL(1UL, cigar.size());
    CPPUNIT_ASSERT_EQUAL(1600U, cigar[0]);
    CPPUNIT_ASSERT_EQUAL(242U, cigar[1]);
    cigar.clear();
    for (unsigned i = 1; i < 15; ++i)
    {
        query = subv(database, i, 100);
        bsw.align(query, database.begin(), database.end(), cigar);
        CPPUNIT_ASSERT_EQUAL(1UL, cigar.size());
        CPPUNIT_ASSERT_EQUAL(1600U, cigar[0]);
        cigar.clear();
    }
    query = subv(database, 15, 100);
    bsw.align(query, database.begin(), database.end(), cigar);
    CPPUNIT_ASSERT_EQUAL(1UL, cigar.size());
    CPPUNIT_ASSERT_EQUAL(1600U, cigar[0]);
    cigar.clear();
}

void TestBandedSmithWaterman::testSingleDeletion()
{
    // Note: Deletions in homopolymer won't show at the right place
    // Note: the gap penalties make it impractical to have inserts greater than 9
    const unsigned left = 40;
    const unsigned right = 40;
    const std::string deletion = "AGAGCAGCGAGCGACAGCAGCAGCAAA"; // no Ts to ensure constant location
    for (unsigned deletionLength = 1; 13 >= deletionLength; ++deletionLength)
    {
        const unsigned dl = 7 - (deletionLength / 2);
        const std::vector<char> dlS = subv(genome, 100, dl);
        const std::vector<char> leftS = subv(genome, 100 + dl, left - 1) + "T";
        const std::vector<char> rightS = subv(genome, 100 + dl + left, right);
        const std::vector<char> deletionS1 = subv(deletion, 0, deletionLength);
        const unsigned drD = 15 - dl - deletionS1.size();
        const std::vector<char> drDS = subv(genome, 100 + dl + left + right, drD);
        const std::vector<char> databaseD = dlS + leftS + deletionS1 + rightS + drDS;
        const std::vector<char> queryD = leftS + rightS;
        isaac::alignment::Cigar cigar; cigar.reserve(1024);
        bsw.align(queryD, databaseD.begin(), databaseD.end(), cigar);
        CPPUNIT_ASSERT_EQUAL(3UL, cigar.size());
        CPPUNIT_ASSERT_EQUAL(( left << 4) | 0U, cigar[0]);
        CPPUNIT_ASSERT_EQUAL((unsigned)( deletionS1.size() << 4) | 2U, cigar[1]);
        CPPUNIT_ASSERT_EQUAL((right << 4) | 0U, cigar[2]);
        cigar.clear();
    }
}

void TestBandedSmithWaterman::testSingleInsertion()
{
    // Note: Insertions require appropriate database length at the beginning
    const std::vector<char> database = subv(genome, 100, 220);
    unsigned queryLength = database.size() - 15;
    // Note: the gap penalties make it impractical to have inserts greater than 9
    for (unsigned insertLength = 1; 9 >= insertLength; ++insertLength)
    {
        unsigned left = 100;
        unsigned right = queryLength - left - insertLength;
        unsigned dl = 9;
        std::vector<char> query = subv(database, dl, left) + std::string(insertLength, 'T') + subv(database, left + dl, right);
        //CPPUNIT_ASSERT_EQUAL(115UL, database.size());
        isaac::alignment::Cigar cigar; cigar.reserve(1024);
        bsw.align(query, database.begin(), database.end(), cigar);
//        CPPUNIT_ASSERT_EQUAL_MESSAGE(isaac::alignment::Cigar::toString(cigar.begin(), cigar.end()) +
//                             " db:" + string(database.begin(), database.end()) +
//                             " qw:" + string(query.begin(), query.end()),
//                             std::string("99M1I105M"), isaac::alignment::Cigar::toString(cigar.begin(), cigar.end()));
        CPPUNIT_ASSERT_EQUAL(3UL, cigar.size());
        CPPUNIT_ASSERT_EQUAL_MESSAGE(
            isaac::alignment::Cigar::toString(cigar.begin(), cigar.end()) +
            " db:" + string(database.begin(), database.end()) +
            " qw:" + string(query.begin(), query.end()), ( left << 4) | 0U, cigar[0]);
        CPPUNIT_ASSERT_EQUAL(( insertLength << 4) | 1U, cigar[1]);
        CPPUNIT_ASSERT_EQUAL((right << 4) | 0U, cigar[2]);
        cigar.clear();
    }
}

void TestBandedSmithWaterman::testMultipleIndels()
{
    // Note: the gap penalties make it impractical to have inserts greater than 9
    const unsigned left = 20;
    const unsigned center = 20;
    const unsigned right = 20;
    const unsigned dl = 6;
    const std::string dlS = genome.substr(100, dl);
    const std::string leftS = genome.substr(100 + dl, left - 1) + "T";
    const std::string insertS1 = "A";
    const std::string insertS2 = "CG";
    const std::string centerS = genome.substr(100 + dl + left, center - 1) + "T";
    const std::string deletionS1 = "AAG";
    const std::string deletionS2 = "ACAG";
    const std::string rightS = genome.substr(100 + dl + left + center, right);
    // 1 Insertion and 1 Deletion
    const unsigned drID = 15 - dl + insertS1.length() - deletionS2.length();
    const std::string drIDS = genome.substr(100 + dl + left + center + right, drID);
    const std::vector<char> databaseID = vectorFromString(dlS + leftS + centerS + deletionS2 + rightS + drIDS);
    const std::vector<char> queryID = vectorFromString(leftS + insertS1 + centerS + rightS);
    //CPPUNIT_ASSERT_EQUAL(115UL, database.size());
    isaac::alignment::Cigar cigar; cigar.reserve(1024);
    bsw.align(queryID, databaseID.begin(), databaseID.end(), cigar);
    CPPUNIT_ASSERT_EQUAL(5UL, cigar.size());
    CPPUNIT_ASSERT_EQUAL(( left << 4) | 0U, cigar[0]);
    CPPUNIT_ASSERT_EQUAL((unsigned)( insertS1.length() << 4) | 1U, cigar[1]);
    CPPUNIT_ASSERT_EQUAL((center << 4) | 0U, cigar[2]);
    CPPUNIT_ASSERT_EQUAL((unsigned)(deletionS2.length() << 4) | 2U, cigar[3]);
    CPPUNIT_ASSERT_EQUAL((right << 4) | 0U, cigar[4]);
    cigar.clear();
    // 2 Insertions
    const unsigned drI2 = 15 - dl + insertS1.length() +insertS2.length();
    const std::string drI2S = genome.substr(100 + dl + left + center + right, drI2);
    const std::vector<char> databaseI2 = vectorFromString(dlS + leftS + centerS + rightS + drI2S);
    const std::vector<char> queryI2 = vectorFromString(leftS + insertS1 + centerS + insertS2 + rightS);
    bsw.align(queryI2, databaseI2.begin(), databaseI2.end(), cigar);
    CPPUNIT_ASSERT_EQUAL(5UL, cigar.size());
    CPPUNIT_ASSERT_EQUAL(( left << 4) | 0U, cigar[0]);
    CPPUNIT_ASSERT_EQUAL((unsigned)(insertS1.length() << 4) | 1U, cigar[1]);
    CPPUNIT_ASSERT_EQUAL((center << 4) | 0U, cigar[2]);
    CPPUNIT_ASSERT_EQUAL((unsigned)(insertS2.length() << 4) | 1U, cigar[3]);
    CPPUNIT_ASSERT_EQUAL((right << 4) | 0U, cigar[4]);
    cigar.clear();
    // 2 deletions
    const unsigned drD2 = 15 - dl - deletionS1.length() - deletionS2.length();
    const std::string drD2S = genome.substr(100 + dl + left + center + right, drD2);
    const std::vector<char> databaseD2 = vectorFromString(dlS + leftS + deletionS1 + centerS + deletionS2 + rightS + drD2S);
    const std::vector<char> queryD2 = vectorFromString(leftS + centerS + rightS);
    bsw.align(queryD2, databaseD2.begin(), databaseD2.end(), cigar);
    CPPUNIT_ASSERT_EQUAL(5UL, cigar.size());
    CPPUNIT_ASSERT_EQUAL(( left << 4) | 0U, cigar[0]);
    CPPUNIT_ASSERT_EQUAL((unsigned)(deletionS1.length() << 4) | 2U, cigar[1]);
    CPPUNIT_ASSERT_EQUAL((center << 4) | 0U, cigar[2]);
    CPPUNIT_ASSERT_EQUAL((unsigned)(deletionS2.length() << 4) | 2U, cigar[3]);
    CPPUNIT_ASSERT_EQUAL((right << 4) | 0U, cigar[4]);
    cigar.clear();
}

void TestBandedSmithWaterman::testOverflow()
{
    // -32767
    // bsw(2, -1, 15, 3, 300)
    // exactly the limit minus 1 -- should pass
    CPPUNIT_ASSERT_NO_THROW(isaac::alignment::BandedSmithWaterman(2, -1, 6, 3, 5460));
    // exactly at the limit -- should throw to differentiate just initialized to a proper score
    CPPUNIT_ASSERT_THROW(isaac::alignment::BandedSmithWaterman(2, -1, 7, 3, 4681), isaac::common::InvalidParameterException);
    // a couple more to verify that even large values will work
    CPPUNIT_ASSERT_THROW(isaac::alignment::BandedSmithWaterman(2, -1, 17, 3, 3681), isaac::common::InvalidParameterException);
    CPPUNIT_ASSERT_THROW(isaac::alignment::BandedSmithWaterman(2, -1, 11, 3, 13681), isaac::common::InvalidParameterException);
}
