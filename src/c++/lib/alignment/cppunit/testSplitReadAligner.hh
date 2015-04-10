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
 **
 ** \file testSplitReadAligner.hh
 **
 ** Tests for split read aligner.
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_ALIGNMENT_TEST_SPLIT_READ_ALIGNER_HH
#define iSAAC_ALIGNMENT_TEST_SPLIT_READ_ALIGNER_HH

#include <cppunit/extensions/HelperMacros.h>

#include <string>

#include "alignment/Cigar.hh"
#include "alignment/FragmentMetadata.hh"
#include "alignment/SeedMetadata.hh"
#include "flowcell/Layout.hh"
#include "flowcell/ReadMetadata.hh"
#include "alignment/matchSelector/SequencingAdapter.hh"

class TestSplitReadAligner : public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE( TestSplitReadAligner );
    CPPUNIT_TEST( testEverything );
    CPPUNIT_TEST_SUITE_END();
private:
    isaac::alignment::Cigar cigarBuffer_;

public:
    TestSplitReadAligner();
    void setUp();
    void tearDown();
    void testEverything();

private:
    void align(
        const std::string &read,
        const std::string &reference,
        isaac::alignment::FragmentMetadataList &fragmentMetadataList);

    void align(
        const std::string &readAlignment1,
        const std::string &readAlignment2,
        const std::string &reference,
        isaac::alignment::FragmentMetadataList &fragmentMetadataList);

    void align(
        const std::string &readAlignment1,
        const std::string &readAlignment2,
        const std::string &reference,
        const std::string &reference2,
        isaac::alignment::FragmentMetadataList &fragmentMetadataList);

    void align(
        const isaac::alignment::Cluster &cluster,
        const isaac::flowcell::ReadMetadataList &readMetadatList,
        const std::vector<char> reference1WithoutSpaces,
        const std::vector<char> reference2WithoutSpaces,
        isaac::alignment::FragmentMetadataList &fragmentMetadataList);
};

#endif // #ifndef iSAAC_ALIGNMENT_TEST_SPLIT_READ_ALIGNER_HH

