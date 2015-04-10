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
 ** \file testSplitReadAligner.cpp
 **
 ** Tests for split read aligner.
 **
 ** \author Roman Petrovski
 **/

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <boost/assign.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/lambda/construct.hpp>

using namespace std;

#include "RegistryName.hh"
#include "testSplitReadAligner.hh"

#include "alignment/fragmentBuilder/SplitReadAligner.hh"
#include "alignment/Cluster.hh"
#include "oligo/Nucleotides.hh"

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION( TestSplitReadAligner, registryName("SplitReadAligner"));

static isaac::flowcell::ReadMetadataList getReadMetadataList(const unsigned readLength)
{
    std::vector<isaac::flowcell::ReadMetadata> ret =
        boost::assign::list_of
            (isaac::flowcell::ReadMetadata(1, readLength, 0, 0))
            ;
    return ret;
}

TestSplitReadAligner::TestSplitReadAligner()
{

}

void TestSplitReadAligner::setUp()
{
}

void TestSplitReadAligner::tearDown()
{
}

inline std::vector<char> vectorFromString(const std::string &str)
{
    return std::vector<char>(str.begin(), str.end());
}

namespace testSimpleIndelAligner
{

static const std::string irrelevantQualities("CFCBBBFFFFFFFFFFJJJJJFJJJJJJJJJJJJJJJJJJJJBFFJJFJFFJJFJJJJFJJJJJJJJFFFFFFFFFFFBFFFFFFFFFFFFFFFFFFFFFFFFEEBFHEHDGBDBEDDEGEHHFHEGBHHDDDB<F>FGGBFGGFGCGGGDGGDDFHHHFEGGBGDGGBGGBEGEGGBGEHDHHHGGGGGDGGGG?GGGGCFCEEBFHEHDGBDBEDDEGEHHFHEGBHHDDDBCFCEEBFHEHDGBDBEDDEGEHHFHEGBHHDDDBFFFFFFFEEBFHEHDGBDBEDDEGEHHFHEGBHHDDDB<F>FGGBFGGFGCGGGDGGDDFHHHFEGGBGDGGBGGBEGEGGBGEHDHHHGGGGGDGGGG?GGGGCFCEEBFHEHDGBDBEDDEGEHHFHEGBHHDDDBCFCEEBFHEHDGBDBEDDEGEHHFHEGBH");

struct ReadInit : public std::pair<std::string, std::string>
{
    bool reverse_;
    typedef std::pair<std::string, std::string> BaseType;
    ReadInit(const std::string &read, const bool reverse) :
        BaseType(read, irrelevantQualities), reverse_(reverse)
    {
        ISAAC_ASSERT_MSG(irrelevantQualities.length() >= read.length(), "Size matters");
    }
};

} // namespace TestSimpleIndelAligner

namespace isaac
{
namespace alignment
{

inline void phredToBcl(char &qual)
{
    qual -= 33;
}

template<class InpuT> InpuT& operator >>(InpuT &input, isaac::alignment::Read &read);

template<> testSimpleIndelAligner::ReadInit& operator >><testSimpleIndelAligner::ReadInit >(
    testSimpleIndelAligner::ReadInit &input,
    isaac::alignment::Read &read)
{
    ISAAC_ASSERT_MSG(input.first.length() <= input.second.length(), "sequence and quality must be of equal lengths");

    read.forwardSequence_ = vectorFromString(input.first);
    read.forwardQuality_ = vectorFromString(input.second);
    read.forwardQuality_.resize(read.forwardSequence_.size());

    std::for_each(read.forwardQuality_.begin(), read.forwardQuality_.end(), &phredToBcl);

    read.reverseSequence_ = read.forwardSequence_;
    std::transform(read.reverseSequence_.begin(), read.reverseSequence_.end(), read.reverseSequence_.begin(), &isaac::oligo::reverseBase);
    std::reverse(read.reverseSequence_.begin(), read.reverseSequence_.end());

    read.reverseQuality_ = read.forwardQuality_;
    std::reverse(read.reverseQuality_.begin(), read.reverseQuality_.end());

    if (input.reverse_)
    {
        using std::swap; swap(read.reverseSequence_, read.forwardSequence_);
        using std::swap; swap(read.reverseQuality_, read.forwardQuality_);
    }


    return input;
}

}
}

static const isaac::reference::Contig makeContig(const std::vector<char> forward)
{
    isaac::reference::Contig ret(0, "vasja");
    ret.forward_ = forward;
    return ret;
}

static const int MATCH_SCORE = 0;
static const int MISMATCH_SCORE = -1;
static const int GAP_OPEN_SCORE = -2;
static const int GAP_EXTEND_SCORE = -1;
static const int MIN_GAP_EXTEND_SCORE = -5;


void TestSplitReadAligner::align(
    const std::string &readAlignment1,
    const std::string &readAlignment2,
    const std::string &reference,
    isaac::alignment::FragmentMetadataList &fragmentMetadataList)
{
    align(readAlignment1, readAlignment2, reference, std::string(), fragmentMetadataList);
}

static isaac::alignment::AlignmentCfg alignmentCfg(MATCH_SCORE, MISMATCH_SCORE, GAP_OPEN_SCORE, GAP_EXTEND_SCORE, MIN_GAP_EXTEND_SCORE, 20000);

class TestAligner : public isaac::alignment::fragmentBuilder::SplitReadAligner
{
public:
    TestAligner() : isaac::alignment::fragmentBuilder::SplitReadAligner(alignmentCfg, true){}


    unsigned updateFragmentCigar(
        const isaac::flowcell::ReadMetadataList &readMetadataList,
        const isaac::reference::ContigList &contigList,
        isaac::alignment::FragmentMetadata &fragmentMetadata,
        const unsigned contigId,
        const long strandPosition,
        isaac::alignment::Cigar &cigarBuffer,
        const unsigned cigarOffset) const
    {
        isaac::reference::ContigAnnotations kUniqueAnnotations;
        std::transform(
            contigList.begin(), contigList.end(),
            std::back_inserter(kUniqueAnnotations),
            boost::lambda::bind<isaac::reference::ContigAnnotation>(
                boost::lambda::constructor<isaac::reference::ContigAnnotation>(),
                boost::lambda::bind(&std::vector<char>::size, boost::lambda::bind(&isaac::reference::Contig::forward_, boost::lambda::_1)),
                isaac::reference::AnnotationValue(32)
            ));

        const isaac::alignment::Anchor firstAnchor = fragmentMetadata.firstAnchor_;
        const isaac::alignment::Anchor lastAnchor = fragmentMetadata.lastAnchor_;
        unsigned ret = isaac::alignment::fragmentBuilder::SplitReadAligner::updateFragmentCigar(
            readMetadataList, contigList, kUniqueAnnotations, fragmentMetadata, contigId, strandPosition, cigarBuffer, cigarOffset);

        if (!firstAnchor.empty() || firstAnchor.kUnique_)
        {
            fragmentMetadata.firstAnchor_ = firstAnchor;
        }
        if (!lastAnchor.empty() || lastAnchor.kUnique_)
        {
            fragmentMetadata.lastAnchor_ = lastAnchor;
        }

        return ret;
    }
};


void TestSplitReadAligner::align(
    const std::string &read,
    const std::string &reference,
    isaac::alignment::FragmentMetadataList &fragmentMetadataList)
{
    fragmentMetadataList.resize(2);

    const long long referenceOffset = reference.find_first_not_of(' ');
    const std::string referenceWithoutSpaces = reference.substr(referenceOffset);
    const long long alignment1Pos = read.find_first_not_of(' ');
    const std::string readWithoutSpaces = read.substr(alignment1Pos);

    //cluster is referenced from fragmentMetadataList by pointer. Don't destroy it.
    static isaac::alignment::Cluster cluster(1000);
    testSimpleIndelAligner::ReadInit init(readWithoutSpaces, fragmentMetadataList[0].reverse);
    init >> cluster.at(0);

    static isaac::flowcell::ReadMetadataList readMetadatList;
    readMetadatList = getReadMetadataList(cluster.at(0).getLength());

    fragmentMetadataList[0].readIndex = 0;
    fragmentMetadataList[0].contigId = 0;
    fragmentMetadataList[0].position = alignment1Pos - referenceOffset;

    fragmentMetadataList[1].readIndex = 0;
    fragmentMetadataList[1].contigId = 0;
    fragmentMetadataList[1].position = reference.length() - readWithoutSpaces.length() - referenceOffset;

    const std::vector<char> referenceV(referenceWithoutSpaces.begin(), referenceWithoutSpaces.end());

    align(cluster, readMetadatList, referenceV, std::vector<char>(), fragmentMetadataList);
}

void TestSplitReadAligner::align(
    const std::string &readAlignment1,
    const std::string &readAlignment2,
    const std::string &reference,
    const std::string &reference2,
    isaac::alignment::FragmentMetadataList &fragmentMetadataList)
{
    fragmentMetadataList.resize(2);
    const long long reference1Offset = reference.find_first_not_of(' ');
    const std::string reference1WithoutSpaces = reference.substr(reference1Offset);
    const long long reference2Offset = reference2.find_first_not_of(' ');
    const std::string reference2WithoutSpaces = reference2.empty() ? std::string() : reference2.substr(reference2Offset);
    const long long alignment1Pos = readAlignment1.find_first_not_of(' ');
    const long long alignment2Pos = readAlignment2.find_first_not_of(' ');
    const std::string readWithoutSpaces = readAlignment1.substr(alignment1Pos);

    //cluster is referenced from fragmentMetadataList by pointer. Don't destroy it.
    static isaac::alignment::Cluster cluster(1000);
    testSimpleIndelAligner::ReadInit init(readWithoutSpaces, fragmentMetadataList[0].reverse);
    init >> cluster.at(0);

    ISAAC_ASSERT_MSG(
        cluster.at(0).getStrandSequence(fragmentMetadataList[1].reverse) == vectorFromString(readAlignment2.substr(alignment2Pos)),
        "R/F mismatch readAlignment1:" <<
        std::string(cluster.at(0).getStrandSequence(fragmentMetadataList[1].reverse).begin(),
                    cluster.at(0).getStrandSequence(fragmentMetadataList[1].reverse).end()) <<
        " readAlignment2:" << readAlignment2.substr(alignment2Pos));

    static isaac::flowcell::ReadMetadataList readMetadatList;
    readMetadatList = getReadMetadataList(cluster.at(0).getLength());

    fragmentMetadataList[0].readIndex = 0;
    fragmentMetadataList[0].contigId = 0;
    fragmentMetadataList[0].position = alignment1Pos - reference1Offset;

    fragmentMetadataList[1].readIndex = 0;
    fragmentMetadataList[1].contigId = reference2.empty() ? 0 : 1;
    fragmentMetadataList[1].position = alignment2Pos - (reference2.empty() ? reference1Offset : reference2Offset);

    align(cluster, readMetadatList,
          std::vector<char>(reference1WithoutSpaces.begin(), reference1WithoutSpaces.end()),
          std::vector<char>(reference2WithoutSpaces.begin(), reference2WithoutSpaces.end()),
          fragmentMetadataList);
}

void TestSplitReadAligner::align(
    const isaac::alignment::Cluster &cluster,
    const isaac::flowcell::ReadMetadataList &readMetadatList,
    const std::vector<char> reference1WithoutSpaces,
    const std::vector<char> reference2WithoutSpaces,
    isaac::alignment::FragmentMetadataList &fragmentMetadataList)
{
    fragmentMetadataList.reserve(fragmentMetadataList.size() + 1);
    const TestAligner aligner;
    isaac::reference::ContigList contigList;
    contigList.push_back(makeContig(reference1WithoutSpaces));
    if (!reference2WithoutSpaces.empty())
    {
        contigList.push_back(makeContig(reference2WithoutSpaces));
    }

    isaac::reference::ContigList::const_iterator currentContig = contigList.begin();

    BOOST_FOREACH(isaac::alignment::FragmentMetadata &fragmentMetadata, fragmentMetadataList)
    {
        fragmentMetadata.cluster = &cluster;
        fragmentMetadata.cigarBuffer = &cigarBuffer_;
        fragmentMetadata.cigarOffset = cigarBuffer_.size();
        fragmentMetadata.observedLength = cluster.at(0).getLength();

        if (0 > fragmentMetadata.position)
        {
            cigarBuffer_.addOperation(std::max<long>(-fragmentMetadata.position, fragmentMetadata.leftClipped()), isaac::alignment::Cigar::SOFT_CLIP);
            ++fragmentMetadata.cigarLength;
            fragmentMetadata.observedLength -= std::max<long>(-fragmentMetadata.position, fragmentMetadata.leftClipped());
            fragmentMetadata.position += std::max<long>(-fragmentMetadata.position, fragmentMetadata.leftClipped());
        }
        else if (fragmentMetadata.leftClipped())
        {
            cigarBuffer_.addOperation(fragmentMetadata.leftClipped(), isaac::alignment::Cigar::SOFT_CLIP);
            ++fragmentMetadata.cigarLength;
            fragmentMetadata.observedLength -= fragmentMetadata.leftClipped();
            fragmentMetadata.position += fragmentMetadata.leftClipped();
        }

        long rightClip = 0;
        if (fragmentMetadata.rightClipped() || (fragmentMetadata.position + fragmentMetadata.observedLength > long(currentContig->forward_.size())))
        {
            rightClip = std::max<long>(fragmentMetadata.rightClipped(),
                                       (fragmentMetadata.position + fragmentMetadata.observedLength - currentContig->forward_.size()));
            fragmentMetadata.observedLength -= rightClip;
        }
        cigarBuffer_.addOperation(fragmentMetadata.observedLength, isaac::alignment::Cigar::ALIGN);
        ++fragmentMetadata.cigarLength;
        if (rightClip)
        {
            cigarBuffer_.addOperation(rightClip, isaac::alignment::Cigar::SOFT_CLIP);
            ++fragmentMetadata.cigarLength;
        }
        aligner.updateFragmentCigar(readMetadatList, contigList, fragmentMetadata, fragmentMetadata.contigId, fragmentMetadata.position, cigarBuffer_, fragmentMetadata.cigarOffset);
        if (contigList.end() != currentContig + 1)
        {
            ++currentContig;
        }
    }

    isaac::reference::ContigAnnotations contigAnnotations;
    std::transform(
        contigList.begin(), contigList.end(),
        std::back_inserter(contigAnnotations),
        boost::lambda::bind<isaac::reference::ContigAnnotation>(
            boost::lambda::constructor<isaac::reference::ContigAnnotation>(),
            boost::lambda::bind(&std::vector<char>::size, boost::lambda::bind(&isaac::reference::Contig::forward_, boost::lambda::_1)),
            isaac::reference::AnnotationValue(32)
        ));


    aligner.alignSimpleSv(cigarBuffer_, contigList, contigAnnotations, readMetadatList, isaac::alignment::TemplateLengthStatistics(), fragmentMetadataList);
}

void TestSplitReadAligner::testEverything()
{
    {

        isaac::alignment::FragmentMetadataList fragmentMetadataList;
        align("              AAAAAAAAAAAAAATTCTGTCTCTTATACAAC"
                                                             "AGAGACAGGTGCACCGCCTATACACATCTAGAATAAGAGACAGGTGCACCGCCTATACACATCTAGA",
              "TTTTTTTTTTTTTTAAAAAAAAAAAAAATTCTGTCTCTTATACAACAAGAGACAGGTGCACCGCCTATACACATCTAGAATAAGAGACAGGTGCACCGCCTATACACATCTAGA",
              fragmentMetadataList);

        CPPUNIT_ASSERT_EQUAL(std::size_t(3), fragmentMetadataList.size());
        CPPUNIT_ASSERT_EQUAL(std::string("32M1D67M"), fragmentMetadataList[2].getCigarString());
        CPPUNIT_ASSERT_EQUAL(14UL, fragmentMetadataList[2].getFStrandReferencePosition().getPosition());
        CPPUNIT_ASSERT_EQUAL(0U, fragmentMetadataList[2].getMismatchCount());
        CPPUNIT_ASSERT_EQUAL(1U, fragmentMetadataList[2].getEditDistance());
    }

    {

        isaac::alignment::FragmentMetadataList fragmentMetadataList;
        align("              AAAAAAAAAAAAAATTCTGTCTCTTATACAAC"
                                                           "AAGAGACAGGTGCACCGCCTATACACATCTAGAATAAGAGACAGGTGCACCGCCTATACACATCTAGA",
              "TTTTTTTTTTTTTTAAAAAAAAAAAAAATTCTGTCTCTTATACAACAGAGACAGGTGCACCGCCTATACACATCTAGAATAAGAGACAGGTGCACCGCCTATACACATCTAGA",
              fragmentMetadataList);

        CPPUNIT_ASSERT_EQUAL(std::size_t(3), fragmentMetadataList.size());
        CPPUNIT_ASSERT_EQUAL(std::string("32M1I67M"), fragmentMetadataList[2].getCigarString());
        CPPUNIT_ASSERT_EQUAL(14UL, fragmentMetadataList[2].getFStrandReferencePosition().getPosition());
        CPPUNIT_ASSERT_EQUAL(0U, fragmentMetadataList[2].getMismatchCount());
        CPPUNIT_ASSERT_EQUAL(1U, fragmentMetadataList[2].getEditDistance());
    }

    {
        isaac::alignment::FragmentMetadataList fragmentMetadataList(2);
        fragmentMetadataList[0].reverse = true;
        fragmentMetadataList[1].reverse = true;

        align(
              "TGTTACCGCAATGTTCTGTCTCTTATACAACATCTAGATGTGTATAAGAGACAGGTGCACCGCCTATACACATCTAGAATAAGAGACAGGTGCACCGCCTATA",
              "                                                                                                                                 TGTTACCGCAATGTTCTGTCTCTTATACAACATCTAGATGTGTATAAGAGACAGGTGCACCGCCTATACACATCTAGAATAAGAGACAGGTGCACCGCCTATA",
              "                                             AAGAGACAGGTGCACCGCCTATACACATCTAGAATAAGAGACAGGTGCACCGCCTATAAAAAAAAAAAAAAAAAAAAAAAAAAATGTTACCGCAATGTTCTGTCTCTTATACAACATCTAGATGTGTATAAGAGAC",
              fragmentMetadataList);

        CPPUNIT_ASSERT_EQUAL(std::size_t(2), fragmentMetadataList.size());
    }

    {
        isaac::alignment::FragmentMetadataList fragmentMetadataList(2);
        fragmentMetadataList[0].reverse = true;
        fragmentMetadataList[1].reverse = true;

        align(
              "TGTTACCGCAATGTTCTGTCTCTTATACAACATCTAGATGTGTATAAGAGACAGGTGCACCGCCTATACACATCTAGAATAAGAGACAGGTGCACCGCCTATA",
              "                                                                                                                                 TGTTACCGCAATGTTCTGTCTCTTATACAACATCTAGATGTGTATAAGAGACAGGTGCACCGCCTATACACATCTAGAATAAGAGACAGGTGCACCGCCTATA",
              "TGTTACCGCGACATATATCACATACATATATTTTTAGATGTGTATAAGAGACAGGTGCACCGCCTATACACATCTAGAATAAGAGACAGGTGCACCGCCTATAAAAAAAAAAAAAAAAAAAAAAAAAAATGTTACCGCAATGTTCTGTCTCTTATACAACATCTAGATGTGTATAAGAGAC",
              "TGTTACCGCGACATATATCACATACATATATTTTTAGATGTGTATAAGAGACAGGTGCACCGCCTATACACATCTAGAATAAGAGACAGGTGCACCGCCTATAAAAAAAAAAAAAAAAAAAAAAAAAAATGTTACCGCAATGTTCTGTCTCTTATACAACATCTAGATGTGTATAAGAGAC",
              fragmentMetadataList);

        CPPUNIT_ASSERT_EQUAL(std::size_t(2), fragmentMetadataList.size());
    }

    {
        isaac::alignment::FragmentMetadataList fragmentMetadataList(2);
        fragmentMetadataList[0].reverse = true;
        fragmentMetadataList[1].reverse = false;

        align(
              "TGTTACCGCAATGTTCTGTCTCTTATACAACATCTAGATGTGTATAAGAGACAGGTGCACCGCCTATACACATCTAGAATAAGAGACAGGTGCACCGCCTATA",
              "                                                                                                                                 TATAGGCGGTGCACCTGTCTCTTATTCTAGATGTGTATAGGCGGTGCACCTGTCTCTTATACACATCTAGATGTTGTATAAGAGACAGAACATTGCGGTAACA",
              "TGTTACCGCGACATATATCACATACATATATTTTTAGATGTGTATAAGAGACAGGTGCACCGCCTATACACATCTAGAATAAGAGACAGGTGCACCGCCTATAAAAAAAAAAAAAAAAAAAAAAAAAAATGTTACCGCAATGTTCTGTCTCTTATACAACATCTAGATGTGTATAAGAGAC",
              "TGTTACCGCGACATATATCACATACATATATTTTTAGATGTGTATAAGAGACAGGTGCACCGCCTATACACATCTAGAATAAGAGACAGGTGCACCGCCTATAAAAAAAAAAAAAAAAAAAAAAAAAAATGTTACCGCAATGTTCTGTCTCTTATACAACATCTAGATGTGTATAAGAGAC",
              fragmentMetadataList);

        CPPUNIT_ASSERT_EQUAL(std::size_t(2), fragmentMetadataList.size());
    }

    {
        isaac::alignment::FragmentMetadataList fragmentMetadataList(2);
        fragmentMetadataList[0].reverse = true;
        fragmentMetadataList[1].reverse = false;

        align("                                                                                                                                 TATAGGCGGTGCACCTGTCTCTTATTCTAGATGTGTATAGGCGGTGCACCTGTCTCTTATACACATCTAGATGTTGTATAAGAGACAGAACATTGCGGTAACA",
              "TGTTACCGCAATGTTCTGTCTCTTATACAACATCTAGATGTGTATAAGAGACAGGTGCACCGCCTATACACATCTAGAATAAGAGACAGGTGCACCGCCTATA",
              "TGTTACCGCAATGTTCTGTCTCTTATACAACATCTAGATGTGTATAAGAGACAGGAAACAAAAAACACACACAAACACCACCACACACAAAAACAAAAAACACAAAAAAAAAAAAAAAAAAAAAAAAAATATAGGCGGTGCACCTGTCTCTTATTCTAGATGTGTATAGGCGGTGCACCTGTCTCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAT",
              fragmentMetadataList);

        CPPUNIT_ASSERT_EQUAL(std::size_t(3), fragmentMetadataList.size());
        CPPUNIT_ASSERT_EQUAL(std::string("48M55F177B55M48S"), fragmentMetadataList[2].getCigarString());
        CPPUNIT_ASSERT_EQUAL(129L, fragmentMetadataList[2].getPosition());
        CPPUNIT_ASSERT_EQUAL(0U, fragmentMetadataList[2].getMismatchCount());
        CPPUNIT_ASSERT_EQUAL(177U, fragmentMetadataList[2].getEditDistance());
    }

    { //SAAC-604 crash due to inconsistent comparison for anchor overlap in translocations
        isaac::alignment::FragmentMetadataList fragmentMetadataList(2);
        fragmentMetadataList[0].reverse = true;
        fragmentMetadataList[1].reverse = true;

        align(
              "                                                                                                                                     GGCGGTGCACCTGTCTCTTATTCTAGATGTGTATAGGCGGTGCACCTGTCTCTTATACACATCTAGATGTTGTATAAGAGACAGAACATTGCGGTAACA",
              "                                                                                                                                                                                                                       GGCGGTGCACCTGTCTCTTATTCTAGATGTGTATAGGCGGTGCACCTGTCTCTTATACACATCTAGATGTTGTATAAGAGACAGAACATTGCGGTAACA",
              " GTTACCGCAATGTTCTGTCTCTTATACAACATCTAGATGTGTATAAGAGACAGGAAACAAAAAACACACACAAACACCACCACACACAAAAACAAAAAACACAAAAAAAAAAAAAAAAAAAAAAAAAATATAGGAGGAGTGCCCGTCTCATATTCTAGATGTGTATAGGCGGTGCACCTGTCTCTTATAAAAAAAAAAAAAAAAAATAAGAG",
              " ACAGAACATTGCGGTAACAAAAAAAAAAAAAAAAAAAAAAAAAAATATAGACACATCTAGATGACACATCTAGATGTTGTATAAGAGACAGAACATTGCGGTAACAAAAAAAAAAAAAAAAAAAAAAAAAAATATAGACACATCTAGATGTTGTATAAGAGACAGAACATTGCGGTAACAAAAAAAAAAAAAAAAAAAAAAAAAAATATAGATAGGCGGTGCACCTCTCTCTTATACTAGATGTGTATAGGCGGAGCATCTGACTCTTATACACATCTAGATGTTGTATAAGAGACTGAACATTGCGGTAACAAAAAAAACTGTCTCAAA",
              fragmentMetadataList);

        CPPUNIT_ASSERT_EQUAL(std::size_t(3), fragmentMetadataList.size());
        CPPUNIT_ASSERT_EQUAL(std::string("50M1C82D49M"), fragmentMetadataList[2].getCigarString());
        CPPUNIT_ASSERT_EQUAL(132L, fragmentMetadataList[2].getPosition());
        CPPUNIT_ASSERT_EQUAL(7U, fragmentMetadataList[2].getMismatchCount());
        // no particular reason to require 0 edit distance for inversion/translocation gap. Just so happens that it is at the moment...
        CPPUNIT_ASSERT_EQUAL(7U, fragmentMetadataList[2].getEditDistance());
    }

    {
        isaac::alignment::FragmentMetadataList fragmentMetadataList(2);
        fragmentMetadataList[0].reverse = true;
        fragmentMetadataList[1].reverse = false;

        align(
              "                                                                                                                                 TATAGGCGGTGCACCTGTCTCTTATTCTAGATGTGTATAGGCGGTGCACCTGTCTCTTATACACATCTAGATGTTGTATAAGAGACAGAACATTGCGGTAACA",
              "TGTTACCGCAATGTTCTGTCTCTTATACAACATCTAGATGTGTATAAGAGACAGGTGCACCGCCTATACACATCTAGAATAAGAGACAGGTGCACCGCCTATA",
              " GTTACCGCAATGTTCTGTCTCTTATACAACATCTAGATGTGTATAAGAGACAGGAAACAAAAAACACACACAAACACCACCACACACAAAAACAAAAAACACAAAAAAAAAAAAAAAAAAAAAAAAAATATAGGCGGTGCACCTGTCTCTTATTCTAGATGTGTATAGGCGGTGCACCTGTCTCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
              " GTTACCGCAATGTTCTGTCTCTTATACAACATCTAGATGTGTATAAGAGACAGGAAACAAAAAACACACACAAACACCACCACACACAAAAACAAAAAACACAAAAAAAAAAAAAAAAAAAAAAAAAATATAGGCGGTGCACCTGTCTCTTATTCTAGATGTGTATAGGCGGTGCACCTGTCTCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
              fragmentMetadataList);

        CPPUNIT_ASSERT_EQUAL(std::size_t(3), fragmentMetadataList.size());
        CPPUNIT_ASSERT_EQUAL(std::string("48M55F1C176B1S54M48S"), fragmentMetadataList[2].getCigarString());
        CPPUNIT_ASSERT_EQUAL(128L, fragmentMetadataList[2].getPosition());
        CPPUNIT_ASSERT_EQUAL(0U, fragmentMetadataList[2].getMismatchCount());
        // no particular reason to require 0 edit distance for inversion/translocation gap. Just so happens that it is at the moment...
        CPPUNIT_ASSERT_EQUAL(0U, fragmentMetadataList[2].getEditDistance());
    }

    {
        isaac::alignment::FragmentMetadataList fragmentMetadataList(2);
        fragmentMetadataList[0].reverse = true;
        fragmentMetadataList[1].reverse = false;

        align(
              "TGTTACCGCAATGTTCTGTCTCTTATACAACATCTAGATGTGTATAAGAGACAGGTGCACCGCCTATACACATCTAGAATAAGAGACAGGTGCACCGCCTATA",
              "                                                                                                                                 TATAGGCGGTGCACCTGTCTCTTATTCTAGATGTGTATAGGCGGTGCACCTGTCTCTTATACACATCTAGATGTTGTATAAGAGACAGAACATTGCGGTAACA",
              " GTTACCGCAATGTTCTGTCTCTTATACAACATCTAGATGTGTATAAGAGACAGGAAACAAAAAACACACACAAACACCACCACACACAAAAACAAAAAACACAAAAAAAAAAAAAAAAAAAAAAAAAATATAGGCGGTGCACCTGTCTCTTATTCTAGATGTGTATAGGCGGTGCACCTGTCTCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
              " GTTACCGCAATGTTCTGTCTCTTATACAACATCTAGATGTGTATAAGAGACAGGAAACAAAAAACACACACAAACACCACCACACACAAAAACAAAAAACACAAAAAAAAAAAAAAAAAAAAAAAAAATATAGGCGGTGCACCTGTCTCTTATTCTAGATGTGTATAGGCGGTGCACCTGTCTCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
              fragmentMetadataList);

        CPPUNIT_ASSERT_EQUAL(std::size_t(3), fragmentMetadataList.size());
        CPPUNIT_ASSERT_EQUAL(std::string("1S46M56F1C82D56M47S"), fragmentMetadataList[2].getCigarString());
        CPPUNIT_ASSERT_EQUAL(0L, fragmentMetadataList[2].getPosition());
        CPPUNIT_ASSERT_EQUAL(0U, fragmentMetadataList[2].getMismatchCount());
        // no particular reason to require 0 edit distance for inversion/translocation gap. Just so happens that it is at the moment...
        CPPUNIT_ASSERT_EQUAL(0U, fragmentMetadataList[2].getEditDistance());
    }

    { // same as above but making sure reference clipping does not cause accessing reference contig before the start
        isaac::alignment::FragmentMetadataList fragmentMetadataList(2);
        fragmentMetadataList[0].reverse = true;
        fragmentMetadataList[1].reverse = false;

        align("                                                                                                                                 TATAGGCGGTGCACCTGTCTCTTATTCTAGATGTGTATAGGCGGTGCACCTGTCTCTTATACACATCTAGATGTTGTATAAGAGACAGAACATTGCGGTAACA",
              "TGTTACCGCAATGTTCTGTCTCTTATACAACATCTAGATGTGTATAAGAGACAGGTGCACCGCCTATACACATCTAGAATAAGAGACAGGTGCACCGCCTATA",
              "    ACCGCAATGTTCTGTCTCTTATACAACATCTAGATGTGTATAAGAGACAGGAAACAAAAAACACACACAAACACCACCACACACAAAAACAAAAAACACAAAAAAAAAAAAAAAAAAAAAAAAAATATAGGCGGTGCACCTGTCTCTTATTCTAGATGTGTATAGGCGGTGCACCTGTCTCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAT",
              fragmentMetadataList);

        CPPUNIT_ASSERT_EQUAL(std::size_t(3), fragmentMetadataList.size());
        CPPUNIT_ASSERT_EQUAL(std::string("48M55F173B4S51M48S"), fragmentMetadataList[2].getCigarString());
        CPPUNIT_ASSERT_EQUAL(125L, fragmentMetadataList[2].getPosition());
        CPPUNIT_ASSERT_EQUAL(0U, fragmentMetadataList[2].getMismatchCount());
        CPPUNIT_ASSERT_EQUAL(173U, fragmentMetadataList[2].getEditDistance());
    }


    {
        isaac::alignment::FragmentMetadataList fragmentMetadataList(2);
        fragmentMetadataList[0].reverse = false;
        fragmentMetadataList[1].reverse = true;

        align("                                                                                                                                 TATAGGCGGTGCACCTGTCTCTTATTCTAGATGTGTATAGGCGGTGCACCTGTCTCTTATACACATCTAGATGTTGTATAAGAGACAGAACATTGCGGTAACA",
              "TGTTACCGCAATGTTCTGTCTCTTATACAACATCTAGATGTGTATAAGAGACAGGTGCACCGCCTATACACATCTAGAATAAGAGACAGGTGCACCGCCTATA",
              "TGTTACCGCAATGTTCTGTCTCTTATACAACATCTAGATGTGTATAAGAGACAGGAAACAAAAAACACACACAAACACCACCACACACAAAAACAAAAAACACAAAAAAAAAAAAAAAAAAAAAAAAAATATAGGCGGTGCACCTGTCTCTTATTCTAGATGTGTATAGGCGGTGCACCTGTCTCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAT",
              fragmentMetadataList);

        CPPUNIT_ASSERT_EQUAL(std::size_t(3), fragmentMetadataList.size());
        CPPUNIT_ASSERT_EQUAL(std::string("48M55F177B55M48S"), fragmentMetadataList[2].getCigarString());
        CPPUNIT_ASSERT_EQUAL(129L, fragmentMetadataList[2].getPosition());
        CPPUNIT_ASSERT_EQUAL(0U, fragmentMetadataList[2].getMismatchCount());
        CPPUNIT_ASSERT_EQUAL(177U, fragmentMetadataList[2].getEditDistance());
    }

    {
        isaac::alignment::FragmentMetadataList fragmentMetadataList(2);
        fragmentMetadataList[0].reverse = false;
        fragmentMetadataList[1].reverse = true;

        align("TGTTACCGCAATGTTCTGTCTCTTATACAACATCTAGATGTGTATAAGAGACAGGTGCACCGCCTATACACATCTAGAATAAGAGACAGGTGCACCGCCTATA",
              "                                                                                                                                 TATAGGCGGTGCACCTGTCTCTTATTCTAGATGTGTATAGGCGGTGCACCTGTCTCTTATACACATCTAGATGTTGTATAAGAGACAGAACATTGCGGTAACA",
              "TGTTACCGCAATGTTCTGTCTCTTATACAACATCTAGATGTGTATAAGAGACAGGAAACAAAAAACACACACAAACACCACCACACACAAAAACAAAAAACACAAAAAAAAAAAAAAAAAAAAAAAAAATATAGGCGGTGCACCTGTCTCTTATTCTAGATGTGTATAGGCGGTGCACCTGTCTCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAT",
              fragmentMetadataList);

        CPPUNIT_ASSERT_EQUAL(std::size_t(3), fragmentMetadataList.size());
        CPPUNIT_ASSERT_EQUAL(std::string("47M56F82D56M47S"), fragmentMetadataList[2].getCigarString());
        CPPUNIT_ASSERT_EQUAL(0L, fragmentMetadataList[2].getPosition());
        CPPUNIT_ASSERT_EQUAL(0U, fragmentMetadataList[2].getMismatchCount());
        CPPUNIT_ASSERT_EQUAL(82U, fragmentMetadataList[2].getEditDistance());
    }

    {
        isaac::alignment::FragmentMetadataList fragmentMetadataList(2);
        fragmentMetadataList[0].reverse = true;
        fragmentMetadataList[1].reverse = false;

        align("TGTTACCGCAATGTTCTGTCTCTTATACAACATCTAGATGTGTATAAGAGACAGGTGCACCGCCTATACACATCTAGAATAAGAGACAGGTGCACCGCCTATA",
              "                                                                                                                                 TATAGGCGGTGCACCTGTCTCTTATTCTAGATGTGTATAGGCGGTGCACCTGTCTCTTATACACATCTAGATGTTGTATAAGAGACAGAACATTGCGGTAACA",
              "TGTTACCGCAATGTTCTGTCTCTTATACAACATCTAGATGTGTATAAGAGACAGGAAACAAAAAACACACACAAACACCACCACACACAAAAACAAAAAACACAAAAAAAAAAAAAAAAAAAAAAAAAATATAGGCGGTGCACCTGTCTCTTATTCTAGATGTGTATAGGCGGTGCACCTGTCTCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAT",
              fragmentMetadataList);

        CPPUNIT_ASSERT_EQUAL(std::size_t(3), fragmentMetadataList.size());
        CPPUNIT_ASSERT_EQUAL(std::string("47M56F82D56M47S"), fragmentMetadataList[2].getCigarString());
        CPPUNIT_ASSERT_EQUAL(0L, fragmentMetadataList[2].getPosition());
        CPPUNIT_ASSERT_EQUAL(0U, fragmentMetadataList[2].getMismatchCount());
        CPPUNIT_ASSERT_EQUAL(82U, fragmentMetadataList[2].getEditDistance());
    }

    {
        isaac::alignment::FragmentMetadataList fragmentMetadataList(2);
        fragmentMetadataList[0].reverse = true;
        fragmentMetadataList[1].reverse = false;

        align("TGTTACCGCAATGTTCTGTCTCTTATACAACATCTAGATGTGTATAAGAGACAGGTGCACCGCCTATACACATCTAGAATAAGAGACAGGTGCACCGCCTATA",
              "                                                                                                                                 TATAGGCGGTGCACCTGTCTCTTATTCTAGATGTGTATAGGCGGTGCACCTGTCTCTTATACACATCTAGATGTTGTATAAGAGACAGAACATTGCGGTAACA",
              "TGTTAAAACGAAGCCCTGTCTCTATAAGGACATCTAGATGTGTATAAGAGACAGGTGCACCGCCTATACACATCTAGAATAAGAGACAGGTGCACCGCCTATAAAAAAAAAAAAAAAAAAAAAAAAAAATATATGAAAAATTTGGAAGGCCATAGACAGATGTGTATAGGCGGTGCACCTGTCTCTTATACACATCTAGATGTTGTATAAGAGACAGAACATTGCGGTAACA",
              fragmentMetadataList);

        CPPUNIT_ASSERT_EQUAL(std::size_t(3), fragmentMetadataList.size());
        CPPUNIT_ASSERT_EQUAL(std::string("32S71M0F97D71S32M"), fragmentMetadataList[2].getCigarString());
        CPPUNIT_ASSERT_EQUAL(32L, fragmentMetadataList[2].getPosition());
        CPPUNIT_ASSERT_EQUAL(0U, fragmentMetadataList[2].getMismatchCount());
        CPPUNIT_ASSERT_EQUAL(97U, fragmentMetadataList[2].getEditDistance());
    }

    {
        isaac::alignment::FragmentMetadataList fragmentMetadataList(2);
        fragmentMetadataList[0].reverse = false;
        fragmentMetadataList[1].reverse = true;

        align(      "TAGAGATGTTAAATATTAACTAAATTGTTACCACATCATATTGTTAGAGAATATTTAAAAGAAAAATTATTGAATAGTGCGAGAACCAAAAAGGCAAAA",
                    ///////////////////////////////////////////////|
              "TTTTGCCTTTTTGGTTCTCGCACTATTCAATAATTTTTCTTTTAAATATTCTCTAACAATATGATGTGGTAACAATTTAGTTAATATTTAACATCTCTA",

                    "TTTCATTTCTTTTCACCGTCTCGGGGTTTAGTCAGTTATTTTAACTTAGAATATTTAAAAGAAAAATTATTGAATAGTGCGAGAACCAAAAAGGCAAAA",
              "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACTAACAATATGATGTGGTAACAATTTAGTTAATATTTAACATCTCTA",
              fragmentMetadataList);

        CPPUNIT_ASSERT_EQUAL(std::size_t(3), fragmentMetadataList.size());
        CPPUNIT_ASSERT_EQUAL(std::string("47S52M0F1C47B52S47M"), fragmentMetadataList[2].getCigarString());
        CPPUNIT_ASSERT_EQUAL(47L, fragmentMetadataList[2].getPosition());
        CPPUNIT_ASSERT_EQUAL(0U, fragmentMetadataList[2].getMismatchCount());
        // Gap is not included in edit distance when contig changes
        CPPUNIT_ASSERT_EQUAL(0U, fragmentMetadataList[2].getEditDistance());
    }

    {
        isaac::alignment::FragmentMetadataList fragmentMetadataList(2);
        fragmentMetadataList[0].reverse = false;
        fragmentMetadataList[1].reverse = true;

        align(      "TAGAGATGTTAAATATTAACTAAATTGTTACCACATCATATTGTTAGAGAATATTTAAAAGAAAAATTATTGAATAGTGCGAGAACCAAAAAGGCAAAA",
                    ///////////////////////////////////////////////|
              "TTTTGCCTTTTTGGTTCTCGCACTATTCAATAATTTTTCTTTTAAATATTCTCTAACAATATGATGTGGTAACAATTTAGTTAATATTTAACATCTCTA",

                    "                                     CCCGGAAAACAGAATATTTAAAAGAAAAATTATTGAATAGTGCGAGAACCAAAAAGGCAAAA",
              "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACTAACAATATGATGTGGTAACAATTTAGTTAATATTTAACATCTCTA",
              fragmentMetadataList);

        CPPUNIT_ASSERT_EQUAL(std::size_t(3), fragmentMetadataList.size());
        CPPUNIT_ASSERT_EQUAL(std::string("47S52M0F1C10B52S47M"), fragmentMetadataList[2].getCigarString());
        CPPUNIT_ASSERT_EQUAL(10L, fragmentMetadataList[2].getPosition());
        CPPUNIT_ASSERT_EQUAL(0U, fragmentMetadataList[2].getMismatchCount());
        // Gap is not included in edit distance when contig changes
        CPPUNIT_ASSERT_EQUAL(0U, fragmentMetadataList[2].getEditDistance());
    }

    {
        isaac::alignment::FragmentMetadataList fragmentMetadataList(2);
        fragmentMetadataList[0].reverse = false;
        fragmentMetadataList[1].reverse = true;

        align(      ///////////////////////////////////////////////|
              "TTTTGCCTTTTTGGTTCTCGCACTATTCAATAATTTTTCTTTTAAATATTCTCTAACAATATGATGTGGTAACAATTTAGTTAATATTTAACATCTCTA",
                    "TAGAGATGTTAAATATTAACTAAATTGTTACCACATCATATTGTTAGAGAATATTTAAAAGAAAAATTATTGAATAGTGCGAGAACCAAAAAGGCAAAA",

              "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACTAACAATATGATGTGGTAACAATTTAGTTAATATTTAACATCTCTA",
                    "TAGGATAGCTGAATTTCAGCTAGAATGCTCCAAGTTCCCCGGAAAACAGAATATTTAAAAGAAAAATTATTGAATAGTGCGAGAACCAAAAAGGCAAA",
              fragmentMetadataList);

        CPPUNIT_ASSERT_EQUAL(std::size_t(3), fragmentMetadataList.size());
        CPPUNIT_ASSERT_EQUAL(std::string("52S47M0F1C52B47S51M1S"), fragmentMetadataList[2].getCigarString());
        CPPUNIT_ASSERT_EQUAL(52L, fragmentMetadataList[2].getPosition());
        CPPUNIT_ASSERT_EQUAL(0U, fragmentMetadataList[2].getMismatchCount());
        // Gap is not included in edit distance when contig changes
        CPPUNIT_ASSERT_EQUAL(0U, fragmentMetadataList[2].getEditDistance());
    }


    {
        isaac::alignment::FragmentMetadataList fragmentMetadataList(2);
        fragmentMetadataList[0].reverse = false;
        fragmentMetadataList[1].reverse = true;

        align(      "TAGAGATGTTAAATATTAACTAAATTGTTACCACATCATATTGTTAGAGAATATTTAAAAGAAAAATTATTGAATAGTGCGAGAACCAAAAAGGCAAAA",
                    ///////////////////////////////////////////////|
              "TTTTGCCTTTTTGGTTCTCGCACTATTCAATAATTTTTCTTTTAAATATTCTCTAACAATATGATGTGGTAACAATTTAGTTAATATTTAACATCTCTA",

                    "TTTCATTTCTTTTCACCGTCTCGGGGTTTAGTCAGTTATTTTAACTTAGAATATTTAAAAGAAAAATTATTGAATAGTGCGAGAACCAAAAAGGCA",
              "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACTAACAATATGATGTGGTAACAATTTAGTTAATATTTAACATCTCTA",
              fragmentMetadataList);

        CPPUNIT_ASSERT_EQUAL(std::size_t(3), fragmentMetadataList.size());
        CPPUNIT_ASSERT_EQUAL(std::string("47S49M3S0F1C44B52S47M"), fragmentMetadataList[2].getCigarString());
        CPPUNIT_ASSERT_EQUAL(47L, fragmentMetadataList[2].getPosition());
        CPPUNIT_ASSERT_EQUAL(0U, fragmentMetadataList[2].getMismatchCount());
        // Gap is not included in edit distance when contig changes
        CPPUNIT_ASSERT_EQUAL(0U, fragmentMetadataList[2].getEditDistance());
    }

    {
        isaac::alignment::FragmentMetadataList fragmentMetadataList(2);
        fragmentMetadataList[0].reverse = false;
        fragmentMetadataList[1].reverse = true;

        align(      "TAGAGATGTTAAATATTAACTAAATTGTTACCACATCATATTGTTAGAGAATATTTAAAAGAAAAATTATTGAATAGTGCGAGAACCAAAAAGGCAAAA",
                    ///////////////////////////////////////////////|
              "TTTTGCCTTTTTGGTTCTCGCACTATTCAATAATTTTTCTTTTAAATATTCTCTAACAATATGATGTGGTAACAATTTAGTTAATATTTAACATCTCTA",

                    "TTTCATTTCTTTTCACCGTCTCGGGGTTTAGTCAGTTATTTTAACTTAGAATATTTAAAAGAAAAATTATTGAATAGTGCGAGAACCAAAAAGGCA",
              "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACTAACAATATGATGTGGTAACAATTTAGTTAATATTTAACATCTCTA",
              fragmentMetadataList);

        CPPUNIT_ASSERT_EQUAL(std::size_t(3), fragmentMetadataList.size());
        CPPUNIT_ASSERT_EQUAL(std::string("47S49M3S0F1C44B52S47M"), fragmentMetadataList[2].getCigarString());
        CPPUNIT_ASSERT_EQUAL(47L, fragmentMetadataList[2].getPosition());
        CPPUNIT_ASSERT_EQUAL(0U, fragmentMetadataList[2].getMismatchCount());
        // Gap is not included in edit distance when contig changes
        CPPUNIT_ASSERT_EQUAL(0U, fragmentMetadataList[2].getEditDistance());
    }

    {
        isaac::alignment::FragmentMetadataList fragmentMetadataList(2);
        fragmentMetadataList[0].reverse = false;
        fragmentMetadataList[1].reverse = true;

        align(      "TAGAGATGTTAAATATTAACTAAATTGTTACCACATCATATTGTTAGAGAATATTTAAAAGAAAAATTATTGAATAGTGCGAGAACCAAAAAGGCAAAA",
                    ///////////////////////////////////////////////|
              "TTTTGCCTTTTTGGTTCTCGCACTATTCAATAATTTTTCTTTTAAATATTCTCTAACAATATGATGTGGTAACAATTTAGTTAATATTTAACATCTCTA",

                    "                                              GAGAATATTTAAAAGAAAAATTATTGAATAGTGCGAGAACCAAAAAGGCAAAA",
              "                                                    CTAACAATATGATGTGGTAACAATTTAGTTAATATTTAACATCTCTA",
              fragmentMetadataList);

        // although inversion is possible, we're not reducing any mismatches by it
        CPPUNIT_ASSERT_EQUAL(std::size_t(2), fragmentMetadataList.size());
    }

    {
        isaac::alignment::FragmentMetadataList fragmentMetadataList(2);
        fragmentMetadataList[0].reverse = false;
        fragmentMetadataList[1].reverse = true;

        align(      "TAGAGATGTTAAATATTAACTAAATTGTTACCACATCATATTGTTAGAGAATATTTAAAAGAAAAATTATTGAATAGTGCGAGAACCAAAAAGGCAAAA",
                    ///////////////////////////////////////////////|
              "TTTTGCCTTTTTGGTTCTCGCACTATTCAATAATTTTTCTTTTAAATATTCTCTAACAATATGATGTGGTAACAATTTAGTTAATATTTAACATCTCTA",

                    "                                              GAGAATATTTAAAAGAAAAATTATTGAATAGTGCGAGAACCAAAAAGGCAAAA",
              "                                                     TAACAATATGATGTGGTAACAATTTAGTTAATATTTAACATCTCTA",
              fragmentMetadataList);

        // although inversion is possible, we're not reducing any mismatches by it
        CPPUNIT_ASSERT_EQUAL(std::size_t(2), fragmentMetadataList.size());
    }

    {
        isaac::alignment::FragmentMetadataList fragmentMetadataList(2);
        fragmentMetadataList[0].reverse = false;
        fragmentMetadataList[1].reverse = true;

        align(
              "TTTTGCCTTTTTGGTTCTCGCACTATTCAATAATTTTTCTTTTAAATATTCTCTAACAATATGATGTGGTAACAATTTAGTTAATATTTAACATCTCTA",
                    "TAGAGATGTTAAATATTAACTAAATTGTTACCACATCATATTGTTAGAGAATATTTAAAAGAAAAATTATTGAATAGTGCGAGAACCAAAAAGGCAAAA",

              "                                  TTTTCTTTTAAATATTCTCTAACAATATGATGGTTATTGTTTTTAGTTAATATTTAACATCTCTA",
                    "TAGAGATGTTAAATATTAACTAAATTGTTACCACATCATATGGTTAGAGAATATTTAAAAGAAAA",
              fragmentMetadataList);

        CPPUNIT_ASSERT_EQUAL(std::size_t(3), fragmentMetadataList.size());
        CPPUNIT_ASSERT_EQUAL(std::string("34S32M33F1C32B33M66S"), fragmentMetadataList[2].getCigarString());
        CPPUNIT_ASSERT_EQUAL(0L, fragmentMetadataList[2].getPosition());
        CPPUNIT_ASSERT_EQUAL(0U, fragmentMetadataList[2].getMismatchCount());
        // Gap is not included in edit distance when contig changes
        CPPUNIT_ASSERT_EQUAL(0U, fragmentMetadataList[2].getEditDistance());
    }

    {
        isaac::alignment::FragmentMetadataList fragmentMetadataList(2);
        fragmentMetadataList[0].reverse = false;
        fragmentMetadataList[1].reverse = true;

        align(
              "TTTTGCCTTTTTGGTTCTCGCACTATTCAATAATTTTTCTTTTAAATATTCTCTAACAATATGATGTGGTAACAATTTAGTTAATATTTAACATCTCTA",
                    "TAGAGATGTTAAATATTAACTAAATTGTTACCACATCATATTGTTAGAGAATATTTAAAAGAAAAATTATTGAATAGTGCGAGAACCAAAAAGGCAAAA",

              "TTTTGCCTTTTTGGTTCTCGCACTATTCAATAATTTTTCTTTTAAATATTCTCTAACAATATGATGTGGTAACAATTTAGTTAATATTTAACATCTCTA",
                    "TAGAGATGTTAAATATTAACTAAATTGTTACCACATCATATGGTTAGAGAATATTTAAAAGAAAA",
              fragmentMetadataList);

        CPPUNIT_ASSERT_EQUAL(std::size_t(2), fragmentMetadataList.size());
    }

    {
        isaac::alignment::FragmentMetadataList fragmentMetadataList;

        align("                                                                                                                                 TGTTACCGCAATGTTCTGTCTCTTATACAACATCTAGATGTGTATAAGAGACAGGTGCACCGCCTATACACATCTAGAATAAGAGACAGGTGCACCGCCTATA",
              "TGTTACCGCAATGTTCTGTCTCTTATACAACATCTAGATGTGTATAAGAGACAGGTGCACCGCCTATACACATCTAGAATAAGAGACAGGTGCACCGCCTATA",
              "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCGAGACAGGTGCACCGCCTATACACATCTAGAATAAGAGACAGGTGCACCGCCTATAAAAAAAAAAAAAAAAAAAAAAAAGCGTGTTACCGCAATGTTCTGTCTCTTATACAACATCTAGATGTGTATAATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT",
              fragmentMetadataList);

        CPPUNIT_ASSERT_EQUAL(std::size_t(3), fragmentMetadataList.size());
        CPPUNIT_ASSERT_EQUAL(std::string("47M129B56M"), fragmentMetadataList[2].getCigarString());
        CPPUNIT_ASSERT_EQUAL(0U, fragmentMetadataList[2].getMismatchCount());
        CPPUNIT_ASSERT_EQUAL(129U, fragmentMetadataList[2].getEditDistance());
    }

    {
        isaac::alignment::FragmentMetadataList fragmentMetadataList;

        align("                                                                                                                                 TGTTACCGCAATGTTCTGTCTCTTATACAACATCTAGATGTGTATAAGAGACAGGTGCACCGCCTATACACATCTAGAATAAGAGACAGGTGCACCGCCTATA",
              "TGTTACCGCAATGTTCTGTCTCTTATACAACATCTAGATGTGTATAAGAGACAGGTGCACCGCCTATACACATCTAGAATAAGAGACAGGTGCACCGCCTATA",
              "                                               GAGACAGGTGCACCGTCTATACACATCTAGAATAAGAGACAGGTGCACCGCCTATAAAAAAAAAAAAAAAAAAAAAAAAGCGTGTTACCGCAATGTTCTGTCTCTTATACAACATCTAGATGTGTATAATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT",
              fragmentMetadataList);

        CPPUNIT_ASSERT_EQUAL(std::size_t(3), fragmentMetadataList.size());
        CPPUNIT_ASSERT_EQUAL(std::string("47M129B56M"), fragmentMetadataList[2].getCigarString());
        CPPUNIT_ASSERT_EQUAL(1U, fragmentMetadataList[2].getMismatchCount());
        CPPUNIT_ASSERT_EQUAL(130U, fragmentMetadataList[2].getEditDistance());
    }

    {
        isaac::alignment::FragmentMetadataList fragmentMetadataList;

        align("TGTTACCGCAATGTTCTGTCTCTTATACAACATCTAGATGTGTATAAGAGACAGGTGCACCGCCTATACACATCTAGAATAAGAGACAGGTGCACCGCCTATA",
              "                                                                                                                                 TGTTACCGCAATGTTCTGTCTCTTATACAACATCTAGATGTGTATAAGAGACAGGTGCACCGCCTATACACATCTAGAATAAGAGACAGGTGCACCGCCTATA",
              "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCGAGACAGGTGCACCGCCTATACACATCTAGAATAAGAGACAGGTGCACCGCCTATAAAAAAAAAAAAAAAAAAAAAAAAGCGTGTTACCGCAATGTTCTGTCTCTTATACAACATCTAGATGTGTATAATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT",
              fragmentMetadataList);

        CPPUNIT_ASSERT_EQUAL(std::size_t(3), fragmentMetadataList.size());
        CPPUNIT_ASSERT_EQUAL(std::string("47M129B56M"), fragmentMetadataList[2].getCigarString());
        CPPUNIT_ASSERT_EQUAL(0U, fragmentMetadataList[2].getMismatchCount());
        CPPUNIT_ASSERT_EQUAL(129U, fragmentMetadataList[2].getEditDistance());
    }

    {
        isaac::alignment::FragmentMetadataList fragmentMetadataList(2);
        fragmentMetadataList[0].reverse = true;
        fragmentMetadataList[1].reverse = true;

        align("                                                                                                                                 TGTTACCGCAATGTTCTGTCTCTTATACAACATCTAGATGTGTATAAGAGACAGGTGCACCGCCTATACACATCTAGAATAAGAGACAGGTGCACCGCCTATA",
              "TGTTACCGCAATGTTCTGTCTCTTATACAACATCTAGATGTGTATAAGAGACAGGTGCACCGCCTATACACATCTAGAATAAGAGACAGGTGCACCGCCTATA",
              "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCGAGACAGGTGCACCGCCTATACACATCTAGAATAAGAGACAGGTGCACCGCCTATAAAAAAAAAAAAAAAAAAAAAAAAGCGTGTTACCGCAATGTTCTGTCTCTTATACAACATCTAGATGTGTATAATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT",
              fragmentMetadataList);

        CPPUNIT_ASSERT_EQUAL(std::size_t(3), fragmentMetadataList.size());
        CPPUNIT_ASSERT_EQUAL(std::string("47M129B56M"), fragmentMetadataList[2].getCigarString());
        CPPUNIT_ASSERT_EQUAL(0U, fragmentMetadataList[2].getMismatchCount());
        CPPUNIT_ASSERT_EQUAL(129U, fragmentMetadataList[2].getEditDistance());
    }

    {
        isaac::alignment::FragmentMetadataList fragmentMetadataList;

        align(
              "                           GCCATCTGGGATTTTTTGTGCTTGTTAGCGCATTGTTCTGTCTCTTATACAACATCTTAGGTAAAAAAAAAAAAAAAAAAAAAAAAAATACACATCTAGAATAAGAGACAGGTGCACCGCCTATACATCTGGGATTTTTTGTGCTTTCTTCTGTCTTTGGGGTATGACGTCTCTAGCCT",
              "                            GCCATCTGGGATTTTTTGTGCTTGTTAGCGCATTGTTCTGTCTCTTATACAACATCTTAGGTAAAAAAAAAAAAAAAAAAAAAAAAAATACACATCTAGAATAAGAGACAGGTGCACCGCCTATACATCTGGGATTTTTTGTGCTTTCTTCTGTCTTTGGGGTATGACGTCTCTAGCCT",
              "GCCATCTGGGATTTTTTGTGCTTTCTTCTGTCTTTGGGGTATGACGTCTATGTTACCGCAATGTTCTGTCTCTTATACAACATCTTAGGTAAAAAAAAAAAAAAAAAAAAAAAAATACACATCTAGAATAAGAGACAGGTGCACCGCCTATACATCTGGGATTTTTTGTGCTTTCTTCTGTCTTTGGGGTATGACGTCTCTAGCCT",
              fragmentMetadataList);

        CPPUNIT_ASSERT_EQUAL(std::size_t(3), fragmentMetadataList.size());
        CPPUNIT_ASSERT_EQUAL(std::string("65M1I113M"), fragmentMetadataList[2].getCigarString());
        CPPUNIT_ASSERT_EQUAL(15U, fragmentMetadataList[2].getMismatchCount());
        CPPUNIT_ASSERT_EQUAL(16U, fragmentMetadataList[2].getEditDistance());
    }

    {
        isaac::alignment::FragmentMetadataList fragmentMetadataList(2);
        fragmentMetadataList[0].reverse = true;
        fragmentMetadataList[1].reverse = true;

        align("TGTTACCGCAATGTTCTGTCTCTTATACAACATCTAGATGTGTATAAGAGACAGGTGCACCGCCTATACACATCTAGAATAAGAGACAGGTGCACCGCCTATA",
              "                                                                                                                                 TGTTACCGCAATGTTCTGTCTCTTATACAACATCTAGATGTGTATAAGAGACAGGTGCACCGCCTATACACATCTAGAATAAGAGACAGGTGCACCGCCTATA",
              "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCGAGACAGGTGCACCGCCTATACACATCTAGAATAAGAGACAGGTGCACCGCCTATAAAAAAAAAAAAAAAAAAAAAAAAGCGTGTTACCGCAATGTTCTGTCTCTTATACAACATCTAGATGTGTATAATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT",
              fragmentMetadataList);

        CPPUNIT_ASSERT_EQUAL(std::size_t(3), fragmentMetadataList.size());
        CPPUNIT_ASSERT_EQUAL(std::string("47M129B56M"), fragmentMetadataList[2].getCigarString());
        CPPUNIT_ASSERT_EQUAL(0U, fragmentMetadataList[2].getMismatchCount());
        CPPUNIT_ASSERT_EQUAL(129U, fragmentMetadataList[2].getEditDistance());
    }

    {
        isaac::alignment::FragmentMetadataList fragmentMetadataList;

        align("ATTTGGTTAAGGTAGCGGTAAAAGCGTGTTACCGCAATGTTCTGTCTCTTATACAACATCTAGATGTGTAT"
                                                                                                   "AAGAGACAGGTGCACCGCCTATACACATCTAGAATAAGAGACAGGTGCACCGCCTATACACATCTAGA",
              "ATTTGGTTAAGGTAGCGGTAAAAGCGTGTTACCGCAATGTTCTGTCTCTTATACAACATCTAGATGTGTATAAAAAAAAAAAAAAAAGAGACAGGTGCACCGCCTATACACATCTAGAATAAGAGACAGGTGCACCGCCTATACACATCTAGA",
              fragmentMetadataList);

        CPPUNIT_ASSERT_EQUAL(std::size_t(3), fragmentMetadataList.size());
        CPPUNIT_ASSERT_EQUAL(std::string("71M14D68M"), fragmentMetadataList[2].getCigarString());
        CPPUNIT_ASSERT_EQUAL(0U, fragmentMetadataList[2].getMismatchCount());
        CPPUNIT_ASSERT_EQUAL(14U, fragmentMetadataList[2].getEditDistance());
    }

    { // verify proper preservation of alignment-independent clipping
        isaac::alignment::FragmentMetadataList fragmentMetadataList(2);
        fragmentMetadataList[0].leftClipped() = 8;
        fragmentMetadataList[1].rightClipped() = 7;

        align("ATTTGGTTAAGGTAGCGGTAAAAGCGTGTTACCGCAATGTTCTGTCTCTTATACAACATCTAGATGTGTAT"
                                                                                                   "AAGAGACAGGTGCACCGCCTATACACATCTAGAATAAGAGACAGGTGCACCGCCTATACACATCTAGA",
              "ATTTGGTTAAGGTAGCGGTAAAAGCGTGTTACCGCAATGTTCTGTCTCTTATACAACATCTAGATGTGTATAAAAAAAAAAAAAAAAGAGACAGGTGCACCGCCTATACACATCTAGAATAAGAGACAGGTGCACCGCCTATACACATCTAGA",
              fragmentMetadataList);

        CPPUNIT_ASSERT_EQUAL(std::size_t(3), fragmentMetadataList.size());
        CPPUNIT_ASSERT_EQUAL(std::string("8S63M14D61M7S"), fragmentMetadataList[2].getCigarString());
        CPPUNIT_ASSERT_EQUAL(0U, fragmentMetadataList[2].getMismatchCount());
        CPPUNIT_ASSERT_EQUAL(14U, fragmentMetadataList[2].getEditDistance());
        CPPUNIT_ASSERT_EQUAL(8U, unsigned(fragmentMetadataList[2].leftClipped()));
        CPPUNIT_ASSERT_EQUAL(7U, unsigned(fragmentMetadataList[2].rightClipped()));
    }

    {   // verify that despite the first 32 bases are occupied by an anchor, the earliest possible deletion
        // position is selected at offset 25
        isaac::alignment::FragmentMetadataList fragmentMetadataList;

        align("GGTGCAGACTAGTAACAGTTGGTGGGCCGGCA"
                                                                                 "CTGATCCATTAATATATATTGCCCAGGTGCCGTGGCTCACCTATAATCCCAGCACTTTGAGAGGCCAA",
              "GGTGCAGACTAGTAACAGTTGGTGGGCCGGCACTGATCCATTAATATATATTGCCCAGGTGCCGGCACTGATCCATTAATATATATTGCCCAGGTGCCGTGGCTCACCTATAATCCCAGCACTTTGAGAGGCCAA",
              fragmentMetadataList);

        CPPUNIT_ASSERT_EQUAL(std::size_t(3), fragmentMetadataList.size());
        CPPUNIT_ASSERT_EQUAL(std::string("25M35D75M"), fragmentMetadataList[2].getCigarString());
        CPPUNIT_ASSERT_EQUAL(0U, fragmentMetadataList[2].getMismatchCount());
        CPPUNIT_ASSERT_EQUAL(35U, fragmentMetadataList[2].getEditDistance());
    }

    {   // verify that despite the first 32 bases are occupied by an anchor, the earliest possible deletion
        // position is selected at offset 0 and causes alignment position change instead of cigar beginning with a deletion
        isaac::alignment::FragmentMetadataList fragmentMetadataList;

        align("GGTGCAGACTAGTAACAGTTGGTGGGCCGGCAC"
                                                                               "TGATCCATTAATATATATTGCCCAGGTGCCGTGGCTCACCTATAATCCCAGCACTTTGAGAGGCCAA",
              "GGTGCAGACTAGTAACAGTTGGTGGGCCGGCAGGTGCAGACTAGTAACAGTTGGTGGGCCGGCACTGATCTATTAATATATATTGCCCAGGTGCCGTGGCTCACCTATAATCCCAGCACTTTGAGAGGCCAA",
              fragmentMetadataList);

        CPPUNIT_ASSERT_EQUAL(std::size_t(3), fragmentMetadataList.size());
        CPPUNIT_ASSERT_EQUAL(std::string("100M"), fragmentMetadataList[2].getCigarString());
        CPPUNIT_ASSERT_EQUAL(32UL, fragmentMetadataList[2].getFStrandReferencePosition().getPosition());
        CPPUNIT_ASSERT_EQUAL(1U, fragmentMetadataList[2].getMismatchCount());
        CPPUNIT_ASSERT_EQUAL(1U, fragmentMetadataList[2].getEditDistance());
    }

    { //verifying that it picks the earliest possible position for the insertion gap
        isaac::alignment::FragmentMetadataList fragmentMetadataList;


        align("ATTTGGTTAAGGTAGCGGTAAAAGCGTGTTACCGCAATGTTCTGTCTCTTATACAACTTCTAGATGTGTAT"
                                                                       "AAGAGACAGGTGCACCGCCTATACACATCTAGAATAAGAGACAGGTGCACCGCCTATACACATCTAGA",
              "ATTTGGTTAAGGTAGCGGTAAAAGCGTGTTACCGCAATGTTCTGTCTCTTATACAACAAGAGACAGGTGCACCGCCTATACACATCTAGAATAAGAGACAGGTGCACCGCCTATACACATCTAGA",
              fragmentMetadataList);

        CPPUNIT_ASSERT_EQUAL(std::size_t(3), fragmentMetadataList.size());
        CPPUNIT_ASSERT_EQUAL(std::string("57M14I68M"), fragmentMetadataList[2].getCigarString());
        CPPUNIT_ASSERT_EQUAL(0U, fragmentMetadataList[2].getMismatchCount());
        CPPUNIT_ASSERT_EQUAL(14U, fragmentMetadataList[2].getEditDistance());
    }

    {
        isaac::alignment::FragmentMetadataList fragmentMetadataList;

        align("ATTTGGTTAAGGTAGCGGTAAAAGCGTGTTACCGCAATGTTCTGTCTCTTATACAACATCTAGATGTGTAT"
                                                                       "TAGAGACAGGTGCACCGCCTATACACATCTAGAATAAGAGACAGGTGCACCGCCTATACACATCTAGA",
              "ATTTGGTTAAGGTAGCGGTAAAAGCGTGTTACCGCAATGTTCTGTCTCTTATACAACTAGAGACAGGTGCACCGCCTATACACATCTAGAATAAGAGACAGGTGCACCGCCTATACACATCTAGA",
              fragmentMetadataList);

        CPPUNIT_ASSERT_EQUAL(std::size_t(3), fragmentMetadataList.size());
        CPPUNIT_ASSERT_EQUAL(std::string("57M14I68M"), fragmentMetadataList[2].getCigarString());
        CPPUNIT_ASSERT_EQUAL(0U, fragmentMetadataList[2].getMismatchCount());
        CPPUNIT_ASSERT_EQUAL(14U, fragmentMetadataList[2].getEditDistance());
    }

    { // verify proper preservation of alignment-independent clipping for insertions
        isaac::alignment::FragmentMetadataList fragmentMetadataList(2);
        fragmentMetadataList[0].leftClipped() = 8;
        fragmentMetadataList[1].rightClipped() = 7;

        align("ATTTGGTTAAGGTAGCGGTAAAAGCGTGTTACCGCAATGTTCTGTCTCTTATACAACATCTAGATGTGTAT"
                                                                       "TAGAGACAGGTGCACCGCCTATACACATCTAGAATAAGAGACAGGTGCACCGCCTATACACATCTAGA",
              "ATTTGGTTAAGGTAGCGGTAAAAGCGTGTTACCGCAATGTTCTGTCTCTTATACAACTAGAGACAGGTGCACCGCCTATACACATCTAGAATAAGAGACAGGTGCACCGCCTATACACATCTAGA",
              fragmentMetadataList);

        CPPUNIT_ASSERT_EQUAL(std::size_t(3), fragmentMetadataList.size());
        CPPUNIT_ASSERT_EQUAL(std::string("8S49M14I61M7S"), fragmentMetadataList[2].getCigarString());
        CPPUNIT_ASSERT_EQUAL(0U, fragmentMetadataList[2].getMismatchCount());
        CPPUNIT_ASSERT_EQUAL(14U, fragmentMetadataList[2].getEditDistance());
        CPPUNIT_ASSERT_EQUAL(8U, unsigned(fragmentMetadataList[2].leftClipped()));
        CPPUNIT_ASSERT_EQUAL(7U, unsigned(fragmentMetadataList[2].rightClipped()));
    }

    { //verifying that it picks the earliest possible position for the insertion but not within the anchoring seed
        isaac::alignment::FragmentMetadataList fragmentMetadataList;
        align("TTCTGTCTCTTATACAACAAGTGGATGTGTAA"
                                "AAGAGACAGGTGCTCCGCCTATACACATCTAGAATAAGAGACAGGTGCACCGCCTATACACATCTAGA",
              "TTCTGTCTCTTATACAACAAGTGGATGTGTAACCGCCTATACACATCTAGAATAAGAGACAGGTGCACCGCCTATACACATCTAGA",
              fragmentMetadataList);

        CPPUNIT_ASSERT_EQUAL(std::size_t(3), fragmentMetadataList.size());
        CPPUNIT_ASSERT_EQUAL(std::string("32M14I54M"), fragmentMetadataList[2].getCigarString());
        CPPUNIT_ASSERT_EQUAL(0U, fragmentMetadataList[2].getMismatchCount());
        CPPUNIT_ASSERT_EQUAL(14U, fragmentMetadataList[2].getEditDistance());

    }

    {

        isaac::alignment::FragmentMetadataList fragmentMetadataList;
        align("              AAAAAAAAAAAAAATTCTGTCTCTTATACAAC"
                                              "AAGAGACAGGTGCACCGCCTATACACATCTAGAATAAGAGACAGGTGCACCGCCTATACACATCTAGA",
              "TTTTTTTTTTTTTTAAAAAAAAAAAAAATTCTGTCTCTTATACAACCCGCCTATACACATCTAGAATAAGAGACAGGTGCACCGCCTATACACATCTAGA",
              fragmentMetadataList);

        CPPUNIT_ASSERT_EQUAL(std::size_t(3), fragmentMetadataList.size());
        CPPUNIT_ASSERT_EQUAL(std::string("32M14I54M"), fragmentMetadataList[2].getCigarString());
        CPPUNIT_ASSERT_EQUAL(14UL, fragmentMetadataList[2].getFStrandReferencePosition().getPosition());
        CPPUNIT_ASSERT_EQUAL(0U, fragmentMetadataList[2].getMismatchCount());
        CPPUNIT_ASSERT_EQUAL(14U, fragmentMetadataList[2].getEditDistance());
    }
//
//    { // ensure that the insertion is not placed before head seed or on tail seed
//        std::vector<isaac::alignment::SeedMetadata> seedMetadataList =
//            boost::assign::list_of
//            (isaac::alignment::SeedMetadata( 64, 32, 0, 0))
//            (isaac::alignment::SeedMetadata( 32, 32, 0, 1))
//            ;
//
//        isaac::alignment::FragmentMetadataList fragmentMetadataList;
//        align(                             "TGCCTCTCTGCGCCTGCGCCGGCGCTGTGGGCCTCTCTGCGCCTGCCCCGGCGCTGTGGGCCT"
//                                                        "CTCTGCGCCTGCCCCGGCGCTGTGGGCCTCTCTGCGCCTTTCGCCCGCGCTGTGCGCCTTTGCGA",
//              "GGTGTCTCACCTTCCCCTCATGGGCCTTCTGCCTCTCTGCGCCTGCGCCGGCGCTGTGGGCCTCTCTGCGCCTGCCCCGGCGCTGTGGGCCTCTCTGCGCCTGCCCCGGCGCTGTGGGCCTCTCTGCGCCTTTCGCCCGCGCTGTGCGCCTTTGCGA",
//              seedMetadataList,
//              fragmentMetadataList);
//
//        CPPUNIT_ASSERT_EQUAL(std::size_t(2), fragmentMetadataList.size());
//        CPPUNIT_ASSERT_EQUAL(std::string("128M"), fragmentMetadataList[0].getCigarString());
//        CPPUNIT_ASSERT_EQUAL(0UL, fragmentMetadataList[0].getFStrandReferencePosition().getPosition());
//        CPPUNIT_ASSERT_EQUAL(39U, fragmentMetadataList[0].getMismatchCount());
//        CPPUNIT_ASSERT_EQUAL(39U, fragmentMetadataList[0].getEditDistance());
//    }
//
//    { // ensure that the insertion is not placed before head seed
//        std::vector<isaac::alignment::SeedMetadata> seedMetadataList =
//            boost::assign::list_of
//            (isaac::alignment::SeedMetadata( 64, 32, 0, 0))
//            (isaac::alignment::SeedMetadata( 0, 32, 0, 1))
//            ;
//
//        isaac::alignment::FragmentMetadataList fragmentMetadataList;
//        align(                             "TGCCTCTCTGCGCCTGCGCCGGCGCTGTGGGCCTCTCTGCGCCTGCCCCGGCGCTGTGGGCCT"
//                                                        "CTCTGCGCCTGCCCCGGCGCTGTGGGCCTCTCTGCGCCTTTCGCCCGCGCTGTGCGCCTTTGCGA",
//              "GGTGTCTCACCTTCCCCTCATGGGCCTTCTGCCTCTCTGCGCCTGCGCCGGCGCTGTGGGCCTCTCTGCGCCTGCCCCGGCGCTGTGGGCCTCTCTGCGCCTGCCCCGGCGCTGTGGGCCTCTCTGCGCCTTTCGCCCGCGCTGTGCGCCTTTGCGA",
//              seedMetadataList,
//              fragmentMetadataList);
//
//        CPPUNIT_ASSERT_EQUAL(std::size_t(3), fragmentMetadataList.size());
//        CPPUNIT_ASSERT_EQUAL(std::string("32M29I67M"), fragmentMetadataList[2].getCigarString());
//        CPPUNIT_ASSERT_EQUAL(29UL, fragmentMetadataList[2].getFStrandReferencePosition().getPosition());
//        CPPUNIT_ASSERT_EQUAL(21U, fragmentMetadataList[2].getMismatchCount());
//        CPPUNIT_ASSERT_EQUAL(50U, fragmentMetadataList[2].getEditDistance());
//    }

    { // ensure that reference start soft clipping does not break
        std::vector<isaac::alignment::SeedMetadata> seedMetadataList =
            boost::assign::list_of
            (isaac::alignment::SeedMetadata( 32, 32, 0, 0))
            (isaac::alignment::SeedMetadata( 64, 32, 0, 1))
            ;

        isaac::alignment::FragmentMetadataList fragmentMetadataList;
        align("GGTGTCTCACCTTCCCCTCATGGGCCTTCTGCCTCTCTGCGCCTGCGCCGGCGCTGTGGGCCTC"
                                                                                 "TCTGCGCCTGCCCCGGCGCTGTGGGCCTCTCTGCGCCTGCCCCGGCGCTGTGGGCCTCTCTGCGCCTTTCGCCCGCGCTGTGCGCCTTTGCGA",
              "                     GGGCCTTCTGCCTCTCTGCGCCTGCGCCGGCGCTGTGGGCCTCAAATCTGCGCCTGCCCCGGCGCTGTGGGCCTCTCTGCGCCTGCCCCGGCGCTGTGGGCCTCTCTGCGCCTTTCGCCCGCGCTGTGCGCCTTTGCGA",
              fragmentMetadataList);

        CPPUNIT_ASSERT_EQUAL(std::size_t(3), fragmentMetadataList.size());
        CPPUNIT_ASSERT_EQUAL(std::string("21S43M3D93M"), fragmentMetadataList[2].getCigarString());
        CPPUNIT_ASSERT_EQUAL(0UL, fragmentMetadataList[2].getFStrandReferencePosition().getPosition());
        CPPUNIT_ASSERT_EQUAL(0U, fragmentMetadataList[2].getMismatchCount());
        CPPUNIT_ASSERT_EQUAL(3U, fragmentMetadataList[2].getEditDistance());
    }

    { // ensure that reference start soft clipping does not break
        isaac::alignment::FragmentMetadataList fragmentMetadataList;
        align("GGTGTCTCACCTTCCCCTCATGGGCCTTCTGCCTCTCTGCGCCTGCGCCGGCGCTGTGGGCCTC"
                                                                           "TCTGCGCCTGCCCCGGCGCTGTGGGCCTCTCTGCGCCTGCCCCGGCGCTGTGGGCCTCTCTGCGCCTTTCGCCCGCGCTGTGCGCCTTTGCGA",
              "                     GGGCCTTCTGCCTCTCTGCGCCTGCGCCGGCGCTGTGGGCCTCGCGCCTGCCCCGGCGCTGTGGGCCTCTCTGCGCCTGCCCCGGCGCTGTGGGCCTCTCTGCGCCTTTCGCCCGCGCTGTGCGCCTTTGCGA",
              fragmentMetadataList);

        CPPUNIT_ASSERT_EQUAL(std::size_t(3), fragmentMetadataList.size());
        CPPUNIT_ASSERT_EQUAL(std::string("21S43M3I90M"), fragmentMetadataList[2].getCigarString());
        CPPUNIT_ASSERT_EQUAL(0UL, fragmentMetadataList[2].getFStrandReferencePosition().getPosition());
        CPPUNIT_ASSERT_EQUAL(0U, fragmentMetadataList[2].getMismatchCount());
        CPPUNIT_ASSERT_EQUAL(3U, fragmentMetadataList[2].getEditDistance());
    }

    { // ensure that reference start soft clipping does not break
        isaac::alignment::FragmentMetadataList fragmentMetadataList;
        align("GGTGTCTCACCTTCCCCTCATGGGCCTTCTGCCTCTCTGCGCCTGCGCCGGCGCTGTGGGACTC"
                                                                           "TCTGCGCCTGCCCCGGCGCTGTGGGCCTCTCTGCGCCTGCCCCGGCGCTGTGGGCCTCTCTGCGCCTTTCGCCCGCGCTGTGCGCCTTTGCGA",
              "                 TCATGGGCCTTCTGCCTCTCTGCGCCTGCGCCGGCGCTGTGGGATCTGCGCCTGCCCCGGCGCTGTGGGCCTCTCTGCGCCTGCCCCGGCGCTGTGGGCCTCTCTGCGCCTTTCGCCCGCGCTGTGCGCCTTTGCGA",
              fragmentMetadataList);

        CPPUNIT_ASSERT_EQUAL(std::size_t(3), fragmentMetadataList.size());
        CPPUNIT_ASSERT_EQUAL(std::string("17S44M3I93M"), fragmentMetadataList[2].getCigarString());
        CPPUNIT_ASSERT_EQUAL(0UL, fragmentMetadataList[2].getFStrandReferencePosition().getPosition());
        CPPUNIT_ASSERT_EQUAL(0U, fragmentMetadataList[2].getMismatchCount());
        CPPUNIT_ASSERT_EQUAL(3U, fragmentMetadataList[2].getEditDistance());
    }


    { // avoid gaps with too many mismatches around them

        { // make sure gap gets accepted when there are less than 8 mismatches on each side
            isaac::alignment::FragmentMetadataList fragmentMetadataList;
            align("GGTGTCTCACTTTCCCCTCATGGGCCTTCTGCGAACGACCCCCTGCGCCGGCGCTGTGGGCCTCTCTGCGCCT"
                   /*                              xxxxxxx x                                */     "GCTGTGGGCCTCTCTGCGCCTGCCCCGGCGCTGTGGGCCTCTCTGCGCCTTTCGCCCGCGCTGTGCGCCTTTGCGA",
                  "GGTGTCTCACCTTCCCCTCATGGGCCTTCTGCCTCTCTGCGCCTGCGCCGGCGCTGTGGGCCTCTCTGCGCCTGCCCCGGCGCTGTGGGCCTCTCTGCGCCTGCCCCGGCGCTGTGGGCCTCTCTGCGCCTTTCGCCCGCGCTGTGCGCCTTTGCGA",
                  fragmentMetadataList);

            CPPUNIT_ASSERT_EQUAL(std::size_t(3), fragmentMetadataList.size());
            CPPUNIT_ASSERT_EQUAL(std::string("73M8D76M"), fragmentMetadataList[2].getCigarString());
            CPPUNIT_ASSERT_EQUAL(0UL, fragmentMetadataList[2].getFStrandReferencePosition().getPosition());
            CPPUNIT_ASSERT_EQUAL(9U, fragmentMetadataList[2].getMismatchCount());
            CPPUNIT_ASSERT_EQUAL(17U, fragmentMetadataList[2].getEditDistance());
        }

        { // too many mismatches in the right flank of the gap
            isaac::alignment::FragmentMetadataList fragmentMetadataList;
            align("GGTGTCTCACCTTCCCCTCATGGGCCTTCTGCGAACGACCCCCTGCGCCGGCGCTGTGGGCCTCTCTGCGCCT"     /*xxxxxxx                xx*/
                   /*                              xxxxxxx x                                */     "ATGTGATGCCTCTCTGCGCCTGCGTCGGCGCTGTGGGCCTCTCTGCGCCTTTCGCCCGCGCTGTGCGCCTTTGCGA",
                  "GGTGTCTCACCTTCCCCTCATGGGCCTTCTGCCTCTCTGCGCCTGCGCCGGCGCTGTGGGCCTCTCTGCGCCTGCCCCGGCGCTGTGGGCCTCTCTGCGCCTGCCCCGGCGCTGTGGGCCTCTCTGCGCCTTTCGCCCGCGCTGTGCGCCTTTGCGA",
                  fragmentMetadataList);

            CPPUNIT_ASSERT_EQUAL(std::string("149M"), fragmentMetadataList[0].getCigarString());
        }

        { // too many mismatches in the left flank of the gap
            isaac::alignment::FragmentMetadataList fragmentMetadataList;
            align("GGTGTCTCACCTTCCCCTCATGGGCCTTCTGCGAACGACCCCCTCCGCTGGGGCGGAGGTCCTCACCGCGACT"
                   /*                              xxxxxxx x   x   x  x  x x  x    x x   x  */     "TGTGTGGGCCTCTCTGCGCCTGCCCCGGCGCTGTGGGCCTCTCTGCGCCTTTCGCCCGCGCTGTGCGCCTTTGCGA",
                  "GGTGTCTCACCTTCCCCTCATGGGCCTTCTGCCTCTCTGCGCCTGCGCCGGCGCTGTGGGCCTCTCTGCGCCTGCCCCGGCGCTGTGGGCCTCTCTGCGCCTGCCCCGGCGCTGTGGGCCTCTCTGCGCCTTTCGCCCGCGCTGTGCGCCTTTGCGA",
                  fragmentMetadataList);

            CPPUNIT_ASSERT_EQUAL(std::string("149M"), fragmentMetadataList[0].getCigarString());
        }

        { // ok mismatches in the right flank of the gap
            isaac::alignment::FragmentMetadataList fragmentMetadataList;
            align("GGTGTCTCACCTTCCCCTCATGGGCCTTCTGCGAACGACCCCCTGCGCCGGCGCTGTGGGCCTCTCTGCGCCT"     /*xxx xxx                xx*/
                   /*                              xxxxxxx x                                */     "ATGGGATGCCTCTCTGCGCCTGCGTCGGCGCTGTGGGCCTCTCTGCGCCTTTCGCCCGCGCTGTGCGCCTTTGCGA",
                  "GGTGTCTCACCTTCCCCTCATGGGCCTTCTGCCTCTCTGCGCCTGCGCCGGCGCTGTGGGCCTCTCTGCGCCTGCCCCGGCGCTGTGGGCCTCTCTGCGCCTGCCCCGGCGCTGTGGGCCTCTCTGCGCCTTTCGCCCGCGCTGTGCGCCTTTGCGA",
                  fragmentMetadataList);

            CPPUNIT_ASSERT_EQUAL(std::size_t(3), fragmentMetadataList.size());
            CPPUNIT_ASSERT_EQUAL(std::string("73M8D76M"), fragmentMetadataList[2].getCigarString());
        }

        { // ok mismatches in the left flank of the gap
            isaac::alignment::FragmentMetadataList fragmentMetadataList;
            align("GGTGTCTCACCTTCCCCTCATGGGCCTTCTGCGAACGACCCCCTCCGCTGGGGCGGAGGTCCTCTCCGCGACT"
                   /*                              xxxxxxx x   x   x  x  x x  x      x   x  */     "TGTGTGGGCCTCTCTGCGCCTGCCCCGGCGCTGTGGGCCTCTCTGCGCCTTTCGCCCGCGCTGTGCGCCTTTGCGA",
                  "GGTGTCTCACCTTCCCCTCATGGGCCTTCTGCCTCTCTGCGCCTGCGCCGGCGCTGTGGGCCTCTCTGCGCCTGCCCCGGCGCTGTGGGCCTCTCTGCGCCTGCCCCGGCGCTGTGGGCCTCTCTGCGCCTTTCGCCCGCGCTGTGCGCCTTTGCGA",
                  fragmentMetadataList);

            CPPUNIT_ASSERT_EQUAL(std::size_t(3), fragmentMetadataList.size());
            CPPUNIT_ASSERT_EQUAL(std::string("73M8D76M"), fragmentMetadataList[2].getCigarString());
        }
    }

    {   // SAAC-692 Split Read aligner incorrectly fails assertion for lastBreakpointOffset
        isaac::alignment::FragmentMetadataList fragmentMetadataList;
        align(" CTGTGAGGAACTACTGTCTTCACGCAGAAAGCGTCTAGCCATGGCGTTAGTATGAGTGTCGTACAGCCTCCAGGACCCCCCTCCGGGAGAGCCATAGTGGTCTGCGGAACCGGTGAGTACACCGGAATTGCCAGGACGACCGGGTCCTTTCTTGGATCAACCCGCTCAATGCCTGGAGATTTGGGCGTGCCCCCGCAAGATCGCTAGCCGAGTAGTGTTGGGTCGCGAAAGGCCTTGTGGTACTGCCTGATAGCTGTCTCTTATACACATCTCCGAGCCCACGAGACCAGAGCGGATCT",
              "CTGTGAGGAACTACTGTCTTCACGCAGAAAGCGTCTAGCCATGGCGTTAGTATGAGTGTCGTACAGCCTCCAGGACCCCCCTCCGGGAGAGCCATAGTGGTCTGCGGAACCGGTGAGTACACCGGAATTGCCAGGACGACCGGGTCCTTTCTTGGATCAACCCGCTCAATGCCTGGAGATTTGGGCGTGCCCCCGCAAGATCGCTAGCCGAGTAGTGTTGGGTCGCGAAAGGCCTTGTGGTACTGCCTGATAGCTGTCTCTTATACACATCTCCGAGCCCACGAGACCAGAGCGGATCT",
              "                                              TTAGTATGAGTGTCGTACAGCCTCCAGGACCCCCCTCCCGGGAGAGCCATAGTGGTCTGCGGAACCGGTGAGTACACCGGAATTGCCAGGACGACCGGGTCCTTTCTTGGATCAACCCGCTCAATGCCTGGAGATTTGGGCGTGCCCCCGCAAGATCGCTAGCCGAGTAGTGTTGGGTCGCGAAAGGCCTTGTGGTACTGCCTGATAGCTGTCTCTTATACACATCTCCGAGCCCACGAGACCAGAGCGGATCT",
              fragmentMetadataList);

        CPPUNIT_ASSERT_EQUAL(std::size_t(3), fragmentMetadataList.size());
        CPPUNIT_ASSERT_EQUAL(std::string("46S36M1D217M"), fragmentMetadataList[2].getCigarString());
    }

    {   // SAAC-731 assertion failure in SplitReadAligner when the best insertion location is at the start of the read
        isaac::alignment::FragmentMetadataList fragmentMetadataList(2);
        fragmentMetadataList[1].rightClipped() = 6;
        fragmentMetadataList[1].firstAnchor_ = isaac::alignment::Anchor(0, 0, true);

        align(
              "GGCAGATGCTGGTGGCTGTCTTCCCTCCCAGGGGCAGATGCTGGTGGCTGTCTTCCCTCCCAGGGGCAGATGCTGGTGGCTGTCTTCCCTCCCAGGGGCAGGTGCTGGTGGATGTCTTCTCTCCCAGGAACAGTTCCTGGAGGCTGTCT",
              "                                                                GGCAGATGCTGGTGGCTGTCTTCCCTCCCAGGGGCAGATGCTGGTGGCTGTCTTCCCTCCCAGGGGCAGATGCTGGTGGCTGTCTTCCCTCCCAGGGGCAGGTGCTGGTGGATGTCTTCTCTCCCAGGAACAGTTCCTGGAGGCTGTCT",
              "AGCAGCTCCTGGAGGCTGTCTTCCCTCCCAGGGGCAGATGCTGGTGGCTGTCTTCCCTCCCAGAGGCAGATGCTGGTGGCTGTCTTCCCTCCCAGGGGCAGGTGCTGGTGGATGTCTTCTCTCCCAGGAACAGTTCCTGGAGGCTGTCTTCCCTCCCAGGGGCAGATGCAGATGCTGGTGGCTGGCTTTCTTCCCAGGGGCAGGTGCTGGAGGCTGTCTTCCCTCCCAGGGGCAGGTGCTGGTGGCTGTCTTCCCTCCCAGGGGCAGATGCTGGTGGCTGTCTTCCCTCCCAGGGGCAGG",
              fragmentMetadataList);


        CPPUNIT_ASSERT_EQUAL(std::size_t(2), fragmentMetadataList.size());
    }

    {   // SAAC-731 assertion failure in SplitReadAligner when the best insertion location is at the start of the read
        isaac::alignment::FragmentMetadataList fragmentMetadataList(2);
        fragmentMetadataList[1].rightClipped() = 6;
        fragmentMetadataList[1].firstAnchor_ = isaac::alignment::Anchor(0, 0, true);

        align(
              "GGCAGATGCTGGTGGCTGTCTTCCCTCCCAGGGGCAGATGCTGGTGGCTGTCTTCCCTCCCAGGGGCAGATGCTGGTGGCTGTCTTCCCTCCCAGGGGCAGGTGCTGGTGGATGTCTTCTCTCCCAGGAACAGTTCCTGGAGGCTGTCT",
              "                                                                GGCAGATGCTGGTGGCTGTCTTCCCTCCCAGGGGCAGATGCTGGTGGCTGTCTTCCCTCCCAGGGGCAGATGCTGGTGGCTGTCTTCCCTCCCAGGGGCAGGTGCTGGTGGATGTCTTCTCTCCCAGGAACAGTTCCTGGAGGCTGTCT",
              "AGCAGCTCCTGGAGGCTGTCTGGGCTCCCAGGGGCAGATGCTGGTGGCTGTCTTCCCTCCCAGAGGCAGATGCTGGTGGCTGTCTTCCCTCCCAGGGGCAGGTGCTGGTGGATGTCTTCTCTCCCAGGAACAGTTCCTGGAGGCTGTCTTCCCTCCCAGGGGCAGATGCAGATGCTGGTGGCTGGCTTTCTTCCCAGGGGCAGGTGCTGGAGGCTGTCTTCCCTCCCAGGGGCAGGTGCTGGTGGCTGTCTTCCCTCCCAGGGGCAGATGCTGGTGGCTGTCTTCCCTCCCAGGGGCAGG",
              fragmentMetadataList);


        CPPUNIT_ASSERT_EQUAL(std::size_t(3), fragmentMetadataList.size());
        CPPUNIT_ASSERT_EQUAL(std::string("64S85M"), fragmentMetadataList[2].getCigarString());
    }

}

