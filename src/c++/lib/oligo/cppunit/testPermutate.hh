/**
 ** Isaac Genome Alignment Software
 ** Copyright (c) 2010-2014 Illumina, Inc.
 ** All rights reserved.
 **
 ** This software is provided under the terms and conditions of the
 ** Illumina Public License 1
 **
 ** You should have received a copy of the Illumina Public License 1
 ** along with this program. If not, see
 ** <https://github.com/sequencing/licenses/>.
 **/

#ifndef iSAAC_ALIGNMENT_TEST_PERMUTATE_HH
#define iSAAC_ALIGNMENT_TEST_PERMUTATE_HH

#include <cppunit/extensions/HelperMacros.h>

#include <string>

#include "oligo/Permutate.hh"

class TestPermutate : public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE( TestPermutate );
    CPPUNIT_TEST( testEverything );
    CPPUNIT_TEST_SUITE_END();
private:
public:
    void setUp();
    void tearDown();
    void testFourBlocks();
    void testEightBlocks();
    void testNoErrors();
    void testTwoErrors();
    void testThreeErrors();
    void testFourErrors();

    void testEverything()
    {
        testFourBlocks();
        testEightBlocks();
        testNoErrors();
        testTwoErrors();
        testThreeErrors();
        testFourErrors();
    }

};

#endif // #ifndef iSAAC_ALIGNMENT_TEST_PERMUTATE_HH
