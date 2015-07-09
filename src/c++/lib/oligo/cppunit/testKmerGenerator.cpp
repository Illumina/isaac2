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
#include "testKmerGenerator.hh"

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION( TestKmerGenerator, registryName("KmerGenerator"));

void TestKmerGenerator::setUp()
{
}

void TestKmerGenerator::tearDown()
{
}


void TestKmerGenerator::testUnsigned()
{
    typedef isaac::oligo::KmerGenerator<7, unsigned> KmerGenerator;
    unsigned kmer;
    std::vector<char>::const_iterator position;
    {
        const std::string s("ANAACGTA");
        const std::vector<char> v(s.begin(), s.end());
        KmerGenerator kmerGenerator(v.begin(), v.end());
        CPPUNIT_ASSERT(!kmerGenerator.next(kmer, position));
    }
    {
        const std::string s = std::string("ANAACGTAA");
        const std::vector<char> v(s.begin(), s.end());
        KmerGenerator kmerGenerator(v.begin(), v.end());
        CPPUNIT_ASSERT(kmerGenerator.next(kmer, position));
        CPPUNIT_ASSERT_EQUAL(kmer, 0x01B0U);
        CPPUNIT_ASSERT_EQUAL(position - v.begin(), 2L);
        CPPUNIT_ASSERT(!kmerGenerator.next(kmer, position));
    }
    {
        const std::string s = std::string("ANAACGTAANAAAAAACGT");
        const std::vector<char> v(s.begin(), s.end());
        KmerGenerator kmerGenerator(v.begin(), v.end());
        CPPUNIT_ASSERT(kmerGenerator.next(kmer, position));
        CPPUNIT_ASSERT_EQUAL(kmer, 0x01B0U);
        CPPUNIT_ASSERT_EQUAL(position - v.begin(), 2L);
        CPPUNIT_ASSERT(kmerGenerator.next(kmer, position));
        CPPUNIT_ASSERT_EQUAL(kmer, 0x01U);
        CPPUNIT_ASSERT_EQUAL(position - v.begin(), 10L);
        CPPUNIT_ASSERT(kmerGenerator.next(kmer, position));
        CPPUNIT_ASSERT_EQUAL(kmer, 0x06U);
        CPPUNIT_ASSERT_EQUAL(position - v.begin(), 11L);
        CPPUNIT_ASSERT(kmerGenerator.next(kmer, position));
        CPPUNIT_ASSERT_EQUAL(kmer, 0x01BU);
        CPPUNIT_ASSERT_EQUAL(position - v.begin(), 12L);
        CPPUNIT_ASSERT(!kmerGenerator.next(kmer, position));
    }
}

void TestKmerGenerator::testUnsignedIteratorReturn()
{
    typedef isaac::oligo::KmerGenerator<7, unsigned> KmerGenerator;
    unsigned kmer;
    std::vector<char>::const_iterator position;
    {
        const std::string s("ANAACGTA");
        const std::vector<char> v(s.begin(), s.end());
        KmerGenerator kmerGenerator(v.begin(), v.end());
        CPPUNIT_ASSERT(v.end() == kmerGenerator.next(kmer));
    }
    {
        const std::string s = std::string("ANAACGTAA");
        const std::vector<char> v(s.begin(), s.end());
        KmerGenerator kmerGenerator(v.begin(), v.end());
        CPPUNIT_ASSERT(v.end() != (position = kmerGenerator.next(kmer)));
        CPPUNIT_ASSERT_EQUAL(kmer, 0x01B0U);
        CPPUNIT_ASSERT_EQUAL(position - v.begin(), 2L);
        CPPUNIT_ASSERT(v.end() == kmerGenerator.next(kmer));
    }
    {
        const std::string s = std::string("ANAACGTAANAAAAAACGT");
        const std::vector<char> v(s.begin(), s.end());
        KmerGenerator kmerGenerator(v.begin(), v.end());
        CPPUNIT_ASSERT(v.end() != (position = kmerGenerator.next(kmer)));
        CPPUNIT_ASSERT_EQUAL(kmer, 0x01B0U);
        CPPUNIT_ASSERT_EQUAL(position - v.begin(), 2L);
        CPPUNIT_ASSERT(v.end() != (position = kmerGenerator.next(kmer)));
        CPPUNIT_ASSERT_EQUAL(kmer, 0x01U);
        CPPUNIT_ASSERT_EQUAL(position - v.begin(), 10L);
        CPPUNIT_ASSERT(v.end() != (position = kmerGenerator.next(kmer)));
        CPPUNIT_ASSERT_EQUAL(kmer, 0x06U);
        CPPUNIT_ASSERT_EQUAL(position - v.begin(), 11L);
        CPPUNIT_ASSERT(v.end() != (position = kmerGenerator.next(kmer)));
        CPPUNIT_ASSERT_EQUAL(kmer, 0x01BU);
        CPPUNIT_ASSERT_EQUAL(position - v.begin(), 12L);
        CPPUNIT_ASSERT(v.end() == (position = kmerGenerator.next(kmer)));
    }
}


void TestKmerGenerator::testConstMethods()
{
    unsigned kmer;
    std::vector<char>::const_iterator position;
    {
        const std::string s("AGAACGTA");
        isaac::oligo::KmerGenerator<9, unsigned, std::string::const_iterator> kmerGenerator(s.begin(), s.end());
        CPPUNIT_ASSERT_EQUAL(262143UL, isaac::oligo::getMaxKmer<unsigned long>(9));
        CPPUNIT_ASSERT(!isaac::oligo::generateKmer(9, kmer, s.begin(), s.end()));
    }
    {
        const std::string s = std::string("AACGTAA");
        isaac::oligo::KmerGenerator<7, unsigned, std::string::const_iterator> kmerGenerator(s.begin(), s.end());
        CPPUNIT_ASSERT_EQUAL(16383UL, isaac::oligo::getMaxKmer<unsigned long>(7));
        CPPUNIT_ASSERT(isaac::oligo::generateKmer(7, kmer, s.begin(), s.end()));
        CPPUNIT_ASSERT_EQUAL(kmer, unsigned(BOOST_BINARY(00 00 01 10 11 00 00)));
    }
}
