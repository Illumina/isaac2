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

#ifndef iSAAC_ALIGNMENT_TEST_SPLIT_NUMERIC_HH
#define iSAAC_ALIGNMENT_TEST_SPLIT_NUMERIC_HH

#include <cppunit/extensions/HelperMacros.h>

#include <string>

#include "oligo/Kmer.hh"

class TestSplitNumeric : public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE( TestSplitNumeric );
    CPPUNIT_TEST( testAll );
    CPPUNIT_TEST_SUITE_END();
private:
public:
    void setUp();
    void tearDown();
    void testAll();
    void testUint256();
    void testUint136();
    void testUint160();
};

#endif // #ifndef iSAAC_ALIGNMENT_TEST_SPLIT_NUMERIC_HH
