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

#ifndef iSAAC_ALIGNMENT_TEST_KMER_GENERATOR_HH
#define iSAAC_ALIGNMENT_TEST_KMER_GENERATOR_HH

#include <cppunit/extensions/HelperMacros.h>

#include <string>

#include "oligo/KmerGenerator.hpp"

class TestKmerGenerator : public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE( TestKmerGenerator );
    CPPUNIT_TEST( testUnsigned );
    CPPUNIT_TEST( testUnsignedIteratorReturn );
    CPPUNIT_TEST( testConstMethods );
    CPPUNIT_TEST_SUITE_END();
private:
public:
    void setUp();
    void tearDown();
    void testUnsigned();
    void testUnsignedIteratorReturn();
    void testConstMethods();
};

#endif // #ifndef iSAAC_ALIGNMENT_TEST_KMER_GENERATOR_HH
