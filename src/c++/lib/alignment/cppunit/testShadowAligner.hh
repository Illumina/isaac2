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

#ifndef iSAAC_ALIGNMENT_TEST_SHADOW_ALIGNER_HH
#define iSAAC_ALIGNMENT_TEST_SHADOW_ALIGNER_HH

#include <cppunit/extensions/HelperMacros.h>

#include <vector>

#include "alignment/ShadowAligner.hh"
#include "reference/Contig.hh"

class TestShadowAligner : public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE( TestShadowAligner );
    CPPUNIT_TEST( testRescueShadowShortest );
    CPPUNIT_TEST( testRescueShadowLongest );
    CPPUNIT_TEST_SUITE_END();
private:
    const std::vector<isaac::flowcell::ReadMetadata> readMetadataList;
    const isaac::flowcell::FlowcellLayoutList flowcells;
    const std::vector<isaac::reference::Contig> contigList;
    isaac::reference::ContigAnnotations contigAnnotations;

public:
    TestShadowAligner();
    void setUp();
    void tearDown();
    void testRescueShadowShortest();
    void testRescueShadowLongest();
};

#endif // #ifndef iSAAC_ALIGNMENT_TEST_SHADOW_ALIGNER_HH
