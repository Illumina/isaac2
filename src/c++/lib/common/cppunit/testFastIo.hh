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
 ** \file testFastIo.hh
 **
 ** Unit tests for FastIo.hh
 **
 ** \author Come Raczy
 **/

#ifndef iSAAC_COMMON_CPPUNIT_TEST_FAST_IO
#define iSAAC_COMMON_CPPUNIT_TEST_FAST_IO

#include <cppunit/extensions/HelperMacros.h>

#include <string>

#include "common/FastIo.hh"

class TestFastIo : public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE( TestFastIo );
    CPPUNIT_TEST( testSprintFloatZeros );
    CPPUNIT_TEST( testSprintFloatZerosPadding );
    CPPUNIT_TEST( testSprintFloatSmallPositive );
    CPPUNIT_TEST( testSprintFloatSmallNegative );
    CPPUNIT_TEST( testSprintFloatMediumPositive );
    CPPUNIT_TEST( testSprintFloatMediumNegative );
    CPPUNIT_TEST( testSprintFloatLargePositive );
    CPPUNIT_TEST( testSprintFloatLargeNegative );
    CPPUNIT_TEST( testPutUnsignedInteger );
    CPPUNIT_TEST( testAppendUnsignedInteger );
    CPPUNIT_TEST( testPutInteger );
    CPPUNIT_TEST( testGetUnsignedInteger );
    CPPUNIT_TEST( testGetInteger );
    CPPUNIT_TEST( testBoolIo );
    CPPUNIT_TEST_SUITE_END();
private:
    static const unsigned int BUFFER_SIZE = 1024;
    char buffer_[BUFFER_SIZE+1];
public:
    void setUp();
    void tearDown();
    void testSprintFloatZeros();
    void testSprintFloatZerosPadding();
    void testSprintFloatSmallPositive();
    void testSprintFloatSmallNegative();
    void testSprintFloatMediumPositive();
    void testSprintFloatMediumNegative();
    void testSprintFloatLargePositive();
    void testSprintFloatLargeNegative();
    void testPutUnsignedInteger();
    void testAppendUnsignedInteger();
    void testPutInteger();
    void testGetUnsignedInteger();
    void testGetInteger();
    void testBoolIo();
};

#endif // #ifndef iSAAC_COMMON_CPPUNIT_TEST_FAST_IO

