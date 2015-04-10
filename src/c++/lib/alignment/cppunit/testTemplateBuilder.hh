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

#ifndef iSAAC_ALIGNMENT_TEST_TEMPLATE_BUILDER_HH
#define iSAAC_ALIGNMENT_TEST_TEMPLATE_BUILDER_HH

#include <cppunit/extensions/HelperMacros.h>

#include <string>
#include <vector>

#include "alignment/TemplateBuilder.hh"
#include "flowcell/ReadMetadata.hh"

class DummyTemplateLengthStatistics: public isaac::alignment::TemplateLengthStatistics
{
public:
    DummyTemplateLengthStatistics();
};

class TestTemplateBuilder : public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE( TestTemplateBuilder );
/**
 * WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING
 *  Some of these tests use rand() to initialize their data. If you change the order of execution, they most
 *  likely start failing. If you are really motivated to do that, please also change rand to something else!
 */
    CPPUNIT_TEST( testConstructor );
    CPPUNIT_TEST( testEmptyMatchList );
    CPPUNIT_TEST( testOrphan );
    CPPUNIT_TEST( testUnique );
    CPPUNIT_TEST( testMultiple );
//    CPPUNIT_TEST( testAll );
    CPPUNIT_TEST_SUITE_END();
private:
    const std::vector<isaac::flowcell::ReadMetadata> readMetadataList;
    const isaac::flowcell::FlowcellLayoutList flowcells;
    const std::vector<isaac::reference::Contig> contigList;
    isaac::reference::ContigAnnotations contigAnnotations;
    const isaac::alignment::RestOfGenomeCorrection restOfGenomeCorrection;
    const isaac::alignment::Cigar cigarBuffer;
    const DummyTemplateLengthStatistics tls;
    const std::vector<char> bcl0;
    const std::vector<char> bcl2;
    const unsigned tile0;
    const unsigned tile2;
    const unsigned clusterId0;
    const unsigned clusterId2;
    isaac::alignment::Cluster cluster0;
    isaac::alignment::Cluster cluster2;
    const isaac::alignment::FragmentMetadata f0_0;
    const isaac::alignment::FragmentMetadata f0_1;
private:
    void checkUnalignedTemplate(
        const isaac::alignment::BamTemplate &bamTemplate,
        const isaac::alignment::Cluster &cluster) const;
    void checkUnalignedFragment(
        const isaac::alignment::BamTemplate &bamTemplate,
        const isaac::alignment::Cluster &cluster,
        unsigned i,
        unsigned readIndex) const;
public:
    TestTemplateBuilder();
    void setUp();
    void tearDown();
    void testConstructor();
    void testEmptyMatchList();
    void testOrphan();
    void testUnique();
    void testMultiple();
    void testAll();
};

#endif // #ifndef iSAAC_ALIGNMENT_TEST_TEMPLATE_BUILDER_HH
