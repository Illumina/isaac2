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

#ifndef iSAAC_ALIGNMENT_TEST_SA_TAG_MAKER_HH
#define iSAAC_ALIGNMENT_TEST_SA_TAG_MAKER_HH

#include <cppunit/extensions/HelperMacros.h>

class TestSaTagMaker : public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE( TestSaTagMaker );
    CPPUNIT_TEST( testFull );
    CPPUNIT_TEST_SUITE_END();
private:

public:
    TestSaTagMaker();
    void setUp();
    void tearDown();

    void testFull();
};

#endif // #ifndef iSAAC_ALIGNMENT_TEST_SA_TAG_MAKER_HH

