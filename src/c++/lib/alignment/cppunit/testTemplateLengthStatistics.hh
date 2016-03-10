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

#ifndef iSAAC_ALIGNMENT_TEST_TEMPLATE_LENGTH_STATISTICS_HH
#define iSAAC_ALIGNMENT_TEST_TEMPLATE_LENGTH_STATISTICS_HH

#include <cppunit/extensions/HelperMacros.h>

#include <string>
#include <vector>

#include "alignment/TemplateLengthStatistics.hh"
#include "flowcell/ReadMetadata.hh"

class TestTemplateLengthStatistics : public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE( TestTemplateLengthStatistics );
//    CPPUNIT_TEST( testAlignmentModels );
//    CPPUNIT_TEST( testAlignmentClassNames );
//    CPPUNIT_TEST( testStatistics );
//    CPPUNIT_TEST( testMateDriftRange );
//    CPPUNIT_TEST( testNoMateDriftRange );
//    CPPUNIT_TEST( testMateOrientation );
//    CPPUNIT_TEST( testMateMinPosition );
//    CPPUNIT_TEST( testMateMaxPosition );
    CPPUNIT_TEST( testAll );
    CPPUNIT_TEST_SUITE_END();
private:
    void addTemplates(isaac::alignment::TemplateLengthDistribution & tls);
public:
    TestTemplateLengthStatistics();
    void setUp();
    void tearDown();
    void testAll()
    {
        testAlignmentModels();
        testAlignmentClassNames();
        testMateOrientation();
        testMateMinPosition();
        testMateMaxPosition();
        testStatistics();
        testMateDriftRange();
        testNoMateDriftRange();
    }
    void testAlignmentModels();
    void testAlignmentClassNames();
    void testMateOrientation();
    void testMateMinPosition();
    void testMateMaxPosition();
    void testStatistics();
    void testMateDriftRange();
    void testNoMateDriftRange();
};

#endif // #ifndef iSAAC_ALIGNMENT_TEST_TEMPLATE_LENGTH_STATISTICS_HH

