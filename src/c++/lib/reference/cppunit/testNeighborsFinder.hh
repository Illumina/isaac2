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

#ifndef iSAAC_REFERENCE_TEST_NEIGHBORS_FINDER_HH
#define iSAAC_REFERENCE_TEST_NEIGHBORS_FINDER_HH

#include <cppunit/extensions/HelperMacros.h>

#include <string>

class TestNeighborsFinder : public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE( TestNeighborsFinder );
    CPPUNIT_TEST( testFindNeighbors );
    CPPUNIT_TEST_SUITE_END();
private:
public:
    void setUp();
    void tearDown();
    void testFindNeighbors();
};

#endif // #ifndef iSAAC_REFERENCE_TEST_NEIGHBORS_FINDER_HH
