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
 ** \file testSaTagMaker.cpp
 **
 ** Test cases for SaTagMaker.
 **
 ** \author Roman Petrovski
 **/

#include <string>

#include <boost/algorithm/string.hpp>
#include <boost/foreach.hpp>

#include "build/SaTagMaker.hh"
#include "io/Fragment.hh"

using namespace isaac;


#include "RegistryName.hh"
#include "testSaTagMaker.hh"

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION( TestSaTagMaker, registryName("TestSaTagMaker"));

TestSaTagMaker::TestSaTagMaker()
{


}

void TestSaTagMaker::setUp()
{
}

void TestSaTagMaker::tearDown()
{
}

inline std::vector<char> vectorFromString(const std::string &str)
{
    return std::vector<char>(str.begin(), str.end());
}

struct TestFragmentAccessor : public io::FragmentAccessor
{
    static const unsigned READ_LENGTH_MAX_ = 1000;
    static const unsigned CIGAR_LENGTH_MAX_ = 1000;
    unsigned char buffer_[READ_LENGTH_MAX_ + CIGAR_LENGTH_MAX_ * sizeof(unsigned)];
    TestFragmentAccessor(
        const reference::ReferencePosition fStrandPosition,
        const bool reverse,
        const std::string &read,
        const alignment::Cigar &cigar,
        const unsigned observedLength,
        const unsigned short editDistance)
    {
        CPPUNIT_ASSERT(CIGAR_LENGTH_MAX_ > unsigned(cigar.size()));
        CPPUNIT_ASSERT(READ_LENGTH_MAX_ > unsigned(read.size()));

        fStrandPosition_ = fStrandPosition;
        flags_.reverse_ = reverse;

        unsigned char *b = const_cast<unsigned char*>(basesBegin());
        BOOST_FOREACH(unsigned char base,
                      std::make_pair(read.begin(), read.end()))
        {
            if (std::string::npos != std::string("ACGTN").find(base))
            {
                *b++ = (0x03 & oligo::getValue(base)) | ('N' == base ? 0 : 0x20);
                ++readLength_;
            }
        }
        std::copy(cigar.begin(), cigar.end(), const_cast<unsigned*>(cigarBegin()));
        cigarLength_ = cigar.size();
        observedLength_ = observedLength;
        editDistance_ = editDistance;
        mateFStrandPosition_ = fStrandPosition;
        alignmentScore_ = 1;
    }

    const char *begin() const {return reinterpret_cast<const char *>(this);}
    const char *end() const {return reinterpret_cast<const char *>(cigarEnd());}
};

isaac::alignment::Cigar::OpCode cigarOpCodeFromLetter(char l)
{
    return isaac::alignment::Cigar::OpCode(std::string("MIDNSHP=XCBF").find(l));
}

isaac::alignment::Cigar cigarFromString(const std::string cigarString)
{
    std::vector<std::string> numbers;
    std::vector<std::string> letters;
    boost::algorithm::split(numbers, cigarString, boost::algorithm::is_any_of("MIDNSHP=XCBF"));
    boost::algorithm::split(letters, cigarString, boost::algorithm::is_digit());
    letters.erase(std::remove_if(letters.begin(), letters.end(), boost::bind(&std::string::empty, _1)), letters.end());

    isaac::alignment::Cigar cigar;
    for (std::vector<std::string>::const_iterator n = numbers.begin(), l = letters.begin(); numbers.end() != n && letters.end() != l; ++n, ++l)
    {
        unsigned long len = boost::lexical_cast<unsigned long>(*n);
        do
        {
            cigar.push_back(isaac::alignment::Cigar::encode(len, cigarOpCodeFromLetter(*l->begin())));
        }
        while (len);
    }
    return cigar;
}

std::string makeSaTag(
    const isaac::reference::ReferencePosition &excludePos,
    const bool excludeReverse,
    const std::string &excludeCigarString,
    const unsigned splitGapLength,
    const std::string &cigarString,
    const std::string &read,
    const std::string &ref1,
    const std::string &ref2,
    unsigned &excludeNM)
{
    const isaac::alignment::Cigar excludeCigar = cigarFromString(excludeCigarString);

    const long long ref1Offset = ref1.find_first_not_of(' ');
    const std::string ref1WithoutSpaces = ref1.substr(ref1Offset);
    const size_t unclippedPos = std::distance(read.begin(),
                                        std::find_if(read.begin(), read.end(),
                                                     boost::bind(&boost::cref<char>, _1) != ' '));

    const isaac::alignment::Cigar cigar = cigarFromString(cigarString);
    isaac::alignment::Cigar::Component firstComponent = isaac::alignment::Cigar::decode(cigar.at(0));
    const reference::ReferencePosition fStrandPos(0, unclippedPos - ref1Offset +
                                                  (isaac::alignment::Cigar::SOFT_CLIP == firstComponent.second ? firstComponent.first : 0));

    const TestFragmentAccessor fragment(fStrandPos, false, read.substr(unclippedPos), cigar, 0, 0);
    reference::ContigList contigList(1, reference::Contig(0, "chr1"));
    contigList.at(0).forward_ = vectorFromString(ref1WithoutSpaces);
    if (!ref2.empty())
    {
        contigList.push_back(reference::Contig(1, "chr2"));
        const long long ref2Offset = ref2.find_first_not_of(' ');
        const std::string ref2WithoutSpaces = ref2.substr(ref2Offset);
        contigList.at(1).forward_ = vectorFromString(ref2WithoutSpaces);
    }

    std::string ret;
    isaac::build::makeSaTagString(
        fragment, excludePos, excludeReverse, excludeCigar.begin(), excludeCigar.end(),
        60, contigList, ret, splitGapLength, excludeNM);
    return ret;
}

std::string makeSaTag(
    const isaac::reference::ReferencePosition &excludePos,
    const bool excludeReverse,
    const std::string &excludeCigarString,
    const unsigned splitGapLength,
    const std::string &cigarString,
    const std::string &read,
    const std::string &ref1,
    unsigned &excludeNM)
{
    return makeSaTag(excludePos,
                     excludeReverse,
                     excludeCigarString,
                     splitGapLength,
                     cigarString,
                     read,
                     ref1, std::string(),
                     excludeNM);
}

std::string makeSaTag(
    const unsigned splitGapLength,
    const std::string &cigarString,
    const std::string &read,
    const std::string &ref1,
    const std::string &ref2,
    unsigned &excludeNM)
{
    return makeSaTag(
        isaac::reference::ReferencePosition(isaac::reference::ReferencePosition::NoMatch),
        false, "100M",
        splitGapLength, cigarString, read, ref1, ref2, excludeNM);

}

std::string makeSaTag(
    const unsigned splitGapLength,
    const std::string &cigarString,
    const std::string &read,
    const std::string &ref1,
    unsigned &excludeNM)
{
    return makeSaTag(
        splitGapLength, cigarString, read, ref1, std::string(), excludeNM);

}

void TestSaTagMaker::testFull()
{
    unsigned excludeNM = 0;
    {
        const std::string sa = makeSaTag(
            -1U,
            "56M26D47M",
            "                                               GAGACAGGTGCACCGCCTATACACATCTAGAATAAGAGACAGGTGCACCGCCTATATGTTACCGCAATGTTCTGTCTCTTATACAACATCTAGATGTGTATAA",
           //                                                                         GAGACAGGTGCACCGCCTATACACATCTAGAATAAGAGACAGGTGCACCGCCTATATGTTACCGCAATGTTCTGTCTCTTATACAACATCTAGATGTGTATAA"
            "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCGAGACAGGTGCACCGCCTATACACATCTAGAATAAGAGACAGGTGCACCGCCTATAAAAAAAAAAAAAAAAAAAAAAAAGCGTGTTACCGCAATGTTCTGTCTCTTATACAACATCTAGATGTGTATAATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT",
            excludeNM);

        CPPUNIT_ASSERT_EQUAL(std::string("chr1,48,+,56M26D47M,60,26;"), sa);
        CPPUNIT_ASSERT_EQUAL(-1U, excludeNM);
    }

    {
        const std::string sa = makeSaTag(
            10,
            "56M26D47M",
            "                                               GAGACAGGTGCACCGCCTATACACATCTAGAATAAGAGACAGGTGCACCGCCTATATGTTACCGCAATGTTCTGTCTCTTATACAACATCTAGATGTGTATAA",
           //                                                                         GAGACAGGTGCACCGCCTATACACATCTAGAATAAGAGACAGGTGCACCGCCTATATGTTACCGCAATGTTCTGTCTCTTATACAACATCTAGATGTGTATAA"
            "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCGAGACAGGTGCACCGCCTATACACATCTAGAATAAGAGACAGGTGCACCGCCTATAAAAAAAAAAAAAAAAAAAAAAAAGCGTGTTACCGCAATGTTCTGTCTCTTATACAACATCTAGATGTGTATAATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT",
            excludeNM);

        CPPUNIT_ASSERT_EQUAL(std::string("chr1,48,+,56M47S,60,0;chr1,130,+,56S47M,60,0;"), sa);
        CPPUNIT_ASSERT_EQUAL(-1U, excludeNM);
    }

    {
        const std::string sa = makeSaTag(
            isaac::reference::ReferencePosition(0, 47),
            false, "47S52M",
            -1U,
            "47S52M0F1C47B52S47M",
                  "TAGAGATGTTAAATATTAACTAAATTGTTACCACATCATATTGTTAGAGAATATTTAAAAGAAAAATTATTGAATAGTGCGAGAACCAAAAAGGCAAAA",
                  ///////////////////////////////////////////////|
//          "TTTTGCCTTTTTGGTTCTCGCACTATTCAATAATTTTTCTTTTAAATATTCTCTAACAATATGATGTGGTAACAATTTAGTTAATATTTAACATCTCTA",

                  "TTTCATTTCTTTTCACCGTCTCGGGGTTTAGTCAGTTATTTTAACTTAGAATATTTAAAAGAAAAATTATTGAATAGTGCGAGAACCAAAAAGGCAAAA",
            "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACTAACAATATGATGTGGTAACAATTTAGTTAATATTTAACATCTCTA",
            excludeNM
            );

        CPPUNIT_ASSERT_EQUAL(std::string("chr2,53,-,52S47M,60,0;"), sa);
        CPPUNIT_ASSERT_EQUAL(0U, excludeNM);
    }

    {
        const std::string sa = makeSaTag(
            -1U,
            "47S51M1S0F1C46B52S47M",
                  "TAGAGATGTTAAATATTAACTAAATTGTTACCACATCATATTGTTAGAGAATATTTAAAAGAAAAATTATTGAATAGTGCGAGAACCAAAAAGGCAAAA",
                  ///////////////////////////////////////////////|
//          "TTTTGCCTTTTTGGTTCTCGCACTATTCAATAATTTTTCTTTTAAATATTCTCTAACAATATGATGTGGTAACAATTTAGTTAATATTTAACATCTCTA",

                  "TTTCATTTCTTTTCACCGTCTCGGGGTTTAGTCAGTTATTTTAACTTAGAATATTTAAAAGAAAAATTATTGAATAGTGCGAGAACCAAAAAGGCAAA",
            "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACTAACAATATGATGTGGTAACAATTTAGTTAATATTTAACATCTCTA",
            excludeNM
            );

        CPPUNIT_ASSERT_EQUAL(std::string("chr1,48,+,47S51M1S,60,0;chr2,53,-,52S47M,60,0;"), sa);
        CPPUNIT_ASSERT_EQUAL(-1U, excludeNM);
    }


    {
        const std::string sa = makeSaTag(
            isaac::reference::ReferencePosition(0, 47),
            false, "47S51M1S",
            -1U,
            "47S51M1S0F1C46B52S47M",
                  "TAGAGATGTTAAATATTAACTAAATTGTTACCACATCATATTGTTAGAGAATATTTAAAAGAAAAATTATTGAATAGTGCGAGAACCAAAAAGGCAAAA",
                  ///////////////////////////////////////////////|
//          "TTTTGCCTTTTTGGTTCTCGCACTATTCAATAATTTTTCTTTTAAATATTCTCTAACAATATGATGTGGTAACAATTTAGTTAATATTTAACATCTCTA",

                  "TTTCATTTCTTTTCACCGTCTCGGGGTTTAGTCAGTTATTTTAACTTAGAATATTTAAAAGAAAAATTATTGAATAGTGCGAGAACCAAAAAGGCAAA",
            "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACTAACAATATGATGTGGTAACAATTTAGTTAATATTTAACATCTCTA",
            excludeNM
            );

        CPPUNIT_ASSERT_EQUAL(std::string("chr2,53,-,52S47M,60,0;"), sa);
        CPPUNIT_ASSERT_EQUAL(0U, excludeNM);
    }

    {
        const std::string sa = makeSaTag(
            -1U,
            "52S47M0F1C52B47S51M1S",
//                  "TAGAGATGTTAAATATTAACTAAATTGTTACCACATCATATTGTTAGAGAATATTTAAAAGAAAAATTATTGAATAGTGCGAGAACCAAAAAGGCAAAA",
                  ///////////////////////////////////////////////|
            "TTTTGCCTTTTTGGTTCTCGCACTATTCAATAATTTTTCTTTTAAATATTCTCTAACAATATGATGTGGTAACAATTTAGTTAATATTTAACATCTCTA",

            "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACTAACAATATGATGTGGTAACAATTTAGTTAATATTTAACATCTCTA",
                  "TTTCATTTCTTTTCACCGTCTCGGGGTTTAGTCAGTTATTTTAACTTAGAATATTTAAAAGAAAAATTATTGAATAGTGCGAGAACCAAAAAGGCAAA",
            excludeNM
            );

        CPPUNIT_ASSERT_EQUAL(std::string("chr1,53,+,52S47M,60,0;chr2,48,-,47S51M1S,60,0;"), sa);
        CPPUNIT_ASSERT_EQUAL(-1U, excludeNM);
    }


    {
        const std::string sa = makeSaTag(
            isaac::reference::ReferencePosition(0, 129),
            false,
            "56S47M",
            10,
            "56M26D47M",
            "                                               GAGACAGGTGCACCGCCTATACACATCTAGAATAAGAGACAGGTGCACCGCCTATATGTTACCGCAATGTTCTGTCTCTTATACAACATCTAGATGTGTATAA",
           //                                                                         GAGACAGGTGCACCGCCTATACACATCTAGAATAAGAGACAGGTGCACCGCCTATATGTTACCGCAATGTTCTGTCTCTTATACAACATCTAGATGTGTATAA"
            "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCGAGACAGGTGCACCGCCTATACACATCTAGAATAAGAGACAGGTGCACCGCCTATAAAAAAAAAAAAAAAAAAAAAAAAGCGTGTTACCGCAATGTTCTGTCTCTTATACAACATCTAGATGTGTATAATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT",
            excludeNM);

        CPPUNIT_ASSERT_EQUAL(std::string("chr1,48,+,56M47S,60,0;"/*"chr1,130,+,56S47M,60,0;"*/), sa);
        CPPUNIT_ASSERT_EQUAL(0U, excludeNM);
    }

    {
        const std::string sa = makeSaTag(
            isaac::reference::ReferencePosition(1, 52),
            true, "52S47M",
            -1U,
            "47S52M0F1C47B52S47M",
                  "TAGAGATGTTAAATATTAACTAAATTGTTACCACATCATATTGTTAGAGAATATTTAAAAGAAAAATTATTGAATAGTGCGAGAACCAAAAAGGCAAAA",
                  ///////////////////////////////////////////////|
//          "TTTTGCCTTTTTGGTTCTCGCACTATTCAATAATTTTTCTTTTAAATATTCTCTAACAATATGATGTGGTAACAATTTAGTTAATATTTAACATCTCTA",

                  "TTTCATTTCTTTTCACCGTCTCGGGGTTTAGTCAGTTATTTTAACTTAGAATATTTAAAAGAAAAATTATTGAATAGTGCGAGAACCAAAAAGGCAAAA",
            "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACTAACAATATGATGTGGTAACAATTTAGTTAATATTTAACATCTCTA",
            excludeNM
            );

        CPPUNIT_ASSERT_EQUAL(std::string("chr1,48,+,47S52M,60,0;"/*"chr2,53,-,52S47M,60,0;"*/), sa);
        CPPUNIT_ASSERT_EQUAL(0U, excludeNM);
    }

    {
        const std::string sa = makeSaTag(
            isaac::reference::ReferencePosition(0, 47),
            false, "47S52M",
            -1U,
            "47S52M0F1C47B52S47M",
                  "TAGAGATGTTAAATATTAACTAAATTGTTACCACATCATATTGTTAGAGAATATTTAAAAGAAAAATTATTGAATAGTGCGAGAACCAAAAAGGCAAAA",
                  ///////////////////////////////////////////////|
//          "TTTTGCCTTTTTGGTTCTCGCACTATTCAATAATTTTTCTTTTAAATATTCTCTAACAATATGATGTGGTAACAATTTAGTTAATATTTAACATCTCTA",

                  "TTTCATTTCTTTTCACCGTCTCGGGGTTTAGTCAGTTATTTTAACTTAGAATATTTAAAAGAAAAATTATTGAATAGTGCGAGAACCAAAAAGGCAAAA",
            "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACTAACAATATGATGTGGTAACAATTTAGTTAATATTTAACATCTCTA",
            excludeNM
            );

        CPPUNIT_ASSERT_EQUAL(std::string(/*"chr1,48,+,47S52M,60,0;"*/"chr2,53,-,52S47M,60,0;"), sa);
        CPPUNIT_ASSERT_EQUAL(0U, excludeNM);
    }

    {
        const std::string sa = makeSaTag(
            isaac::reference::ReferencePosition(0, 47),
            false, "47S52M",
            -1U,
            "47S52M0F1C47B52S47M",
                  "TAGAGATGTTAAATATTAACTAAATTGTTACCACATCATATTGTTAGAGAATATTTAAAAGAAAAATTATTGAATAGTGCGAGAACCAAAAAGGCAAAT",
                  ///////////////////////////////////////////////|
//          "ATTTGCCTTTTTGGTTCTCGCACTATTCAATAATTTTTCTTTTAAATATTCTCTAACAATATGATGTGGTAACAATTTAGTTAATATTTAACATCTCTA",

                  "TTTCATTTCTTTTCACCGTCTCGGGGTTTAGTCAGTTATTTTAACTTAGAATATTTAAAAGAAAAATTATTGAATAGTGCGAGAACCAAAAAGGCAAAA",
            "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACTAACAATATGATGTGGTAACAATTTAGTTAATATTTAACATCTCTA",
            excludeNM
            );

        CPPUNIT_ASSERT_EQUAL(std::string(/*"chr1,48,+,47S52M,60,0;"*/"chr2,53,-,52S47M,60,0;"), sa);
        CPPUNIT_ASSERT_EQUAL(1U, excludeNM);
    }

    {
        const std::string sa = makeSaTag(
            isaac::reference::ReferencePosition(0, 47),
            false,
            "56M47S",
            10,
            "56M26D47M",
            "                                               GAGACAGGTGCACCGCCTATACACATCTAGAATAAGAGACAGGTGCACCGCCTATATGTTACCGCAATGTTCTGTCTCTTATACAACATCTAGATGTGTATAA",
           //                                                                         GAGACAGGTGCACCGCCTATACACATCTAGAATAAGAGACAGGTGCACCGCCTATATGTTACCGCAATGTTCTGTCTCTTATACAACATCTAGATGTGTATAA"
            "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCGAGACAGGTGCACCGCCTATACACATCTAGAATAAGAGACAGGTGCACCGCCTATAAAAAAAAAAAAAAAAAAAAAAAAGCGTGTTACCGCAATGTTCTGTCTCTTATACAACATCTAGATGTGTATAATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT",
            excludeNM);

        CPPUNIT_ASSERT_EQUAL(std::string("chr1,130,+,56S47M,60,0;"), sa);
        CPPUNIT_ASSERT_EQUAL(0U, excludeNM);
    }

    {
        const std::string sa = makeSaTag(
            -1U,
            "47M56F82D56M47S",
            "TGTTACCGCAATGTTCTGTCTCTTATACAACATCTAGATGTGTATAAGAGACAGGTGCACCGCCTATACACATCTAGAATAAGAGACAGGTGCACCGCCTATA",
                                                                                                                                            //TATAGGCGGTGCACCTGTCTCTTATTCTAGATGTGTATAGGCGGTGCACCTGTCTCTTATACACATCTAGATGTTGTATAAGAGACAGAACATTGCGGTAACA
            "TGTTACCGCAATGTTCTGTCTCTTATACAACTTCTAGATGTGTATAAGAGACAGGAAACAAAAAACACACACAAACACCACCACACACAAAAACAAAAAACACAAAAAAAAAAAAAAAAAAAAAAAAAATATAGGCGGTGCACCTGTCTTCTATTCTAGATGTGTATAGGCGGTGCACCTGTCTCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAT",
            excludeNM);

        CPPUNIT_ASSERT_EQUAL(std::string("chr1,1,+,47M56S,60,1;chr1,130,-,56M47S,60,2;"), sa);
        CPPUNIT_ASSERT_EQUAL(-1U, excludeNM);
    }

// this one does not work as only one split per CIGAR is allowed
//    {
//        const std::string sa = makeSaTag(
//            -1U,
//            "47M56F82D36M1C20M",
//            "TGTTACCGCAATGTTCTGTCTCTTATACAACATCTAGATGTGTATAAGAGACAGGTGCACCGCCTATACACATCTAGAATAAGAGACAGGTGCACCGCCTATA",
//                                                                                                                                            //TATAGGCGGTGCACCTGTCTCTTATTCTAGATGTGTATAGGCGGTGCACCTGTCTCTTATACACATCTAGATGTTGTATAAGAGACAGAACATTGCGGTAACA
//            "TGTTACCGCAATGTTCTGTCTCTTATACAACTTCTAGATGTGTATAAGAGACAGGAAACAAAAAACACACACAAACACCACCACACACAAAAACAAAAAACACAAAAAAAAAAAAAAAAAAAAAAAAAATATAGGCGGTGCACCTGTCTTCTATTCTAGATGTGTATAGGCGGTGCACCTGTCTCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAT",
//            "TGTTACCGCAATGTTCTGTCTCTTATACAACTTCTAGATGTGTATAAGAGACAGGAAACAAAAAACACACACAAACACCACCACACACAAAAACAAAAAACACAAAAAAAAAAAAAAAAAAAAAAAAAATATAGGCGGTGCACCTGTCTTCTATTCTAGATGTGTATAGGCGGTGCACCTGTCTCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAT");
//
//        CPPUNIT_ASSERT_EQUAL(std::string("chr1,1,+,47M56S,60,1;chr1,130,-,36M67S,60,2;chr2,166,-,36S20M47S,60,0;"), sa);
//    }

    {
        const std::string sa = makeSaTag(
            -1U,
            "47M56F1C82D56M47S",
            "TGTTACCGCAATGTTCTGTCTCTTATACAACATCTAGATGTGTATAAGAGACAGGTGCACCGCCTATACACATCTAGAATAAGAGACAGGTGCACCGCCTATA",
                                                                                                                                            //TATAGGCGGTGCACCTGTCTCTTATTCTAGATGTGTATAGGCGGTGCACCTGTCTCTTATACACATCTAGATGTTGTATAAGAGACAGAACATTGCGGTAACA
            "TGTTACCGCAATGTTCTGTCTCTTATACAACTTCTAGATGTGTATAAGAGACAGGAAACAAAAAACACACACAAACACCACCACACACAAAAACAAAAAACACAAAAAAAAAAAAAAAAAAAAAAAAAATATAGGCGGTGCACCTGTCTTCTATTCTAGATGTGTATAGGCGGTGCACCTGTCTCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAT",
            "TGTTACCGCAATGTTCTGTCTCTTATACAACTTCTAGATGTGTATAAGAGACAGGAAACAAAAAACACACACAAACACCACCACACACAAAAACAAAAAACACAAAAAAAAAAAAAAAAAAAAAAAAAATATAGGCGGTGCACCTGTCTTCTATTCTAGATGTGTATAGGCGGTGCACCTGTCTCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAT",
            excludeNM);

        CPPUNIT_ASSERT_EQUAL(std::string("chr1,1,+,47M56S,60,1;chr2,130,-,56M47S,60,2;"), sa);
        CPPUNIT_ASSERT_EQUAL(-1U, excludeNM);
    }

    {
        const std::string sa = makeSaTag(
            -1U,
            "48M55F177B55M48S",
            "                                                                                                                                 TATAGGCGGTGCACCTGTCTCTTATTCTAGATGTGTATAGGCGGTGCACCTGTCTCTTATACACATCTAGATGTTGTATAAGAGACAGAACATTGCGGTAACA",
           //TGTTACCGCAATGTTCTGTCTCTTATACAACATCTAGATGTGTATAAGAGACAGGTGCACCGCCTATACACATCTAGAATAAGAGACAGGTGCACCGCCTATA"
            "TGTTACCGCAATGTTCTGTCTCTTATACAACTTCTAGATGTGTATAAGAGACAGGAAACAAAAAACACACACAAACACCACCACACACAAAAACAAAAAACACAAAAAAAAAAAAAAAAAAAAAAAAAATATAGGCGGTGCACCTGTCTTCTATTCTAGATGTGTATAGGCGGTGCACCTGTCTCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAT",
            excludeNM);

        CPPUNIT_ASSERT_EQUAL(std::string("chr1,130,+,48M55S,60,2;chr1,1,-,55M48S,60,1;"), sa);
        CPPUNIT_ASSERT_EQUAL(-1U, excludeNM);
    }

    {
        const std::string sa = makeSaTag(
            -1U,
            "48M55F173B4S51M48S",
            "                                                                                                                                 TATAGGCGGTGCACCTGTCTCTTATTCTAGATGTGTATAGGCGGTGCACCTGTCTCTTATACACATCTAGATGTTGTATAAGAGACAGAACATTGCGGTAACA",
           //TGTTACCGCAATGTTCTGTCTCTTATACAACATCTAGATGTGTATAAGAGACAGGTGCACCGCCTATACACATCTAGAATAAGAGACAGGTGCACCGCCTATA"                                             xx
           //                               x                       xxxxxxxxxxx
            "    ACCGCAATGTTCTGTCTCTTATACAACTTCTAGATGTGTATAAGAGACAGGAAACAAAAAACACACACAAACACCACCACACACAAAAACAAAAAACACAAAAAAAAAAAAAAAAAAAAAAAAAATATAGGCGGTGCACCTGTCTTCTATTCTAGATGTGTATAGGCGGTGCACCTGTCTCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAT",
            excludeNM);

        CPPUNIT_ASSERT_EQUAL(std::string("chr1,126,+,48M55S,60,2;chr1,1,-,4S51M48S,60,1;"), sa);
        CPPUNIT_ASSERT_EQUAL(-1U, excludeNM);
    }

    {
        const std::string sa = makeSaTag(
            -1U,
            "47M129B56M",
            "                                                                                                                                 TGTTACCGCAATGTTCTGTCTCTTATACAACATCTAGATGTGTATAAGAGACAGGTGCACCGCCTATACACATCTAGAATAAGAGACAGGTGCACCGCCTATA",
           //TGTTACCGCAATGTTCTGTCTCTTATACAACATCTAGATGTGTATAAGAGACAGGTGCACCGCCTATACACATCTAGAATAAGAGACAGGTGCACCGCCTATA"
            "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCGAGACAGGTGCACCGCCTATACACATCTAGAATAAGAGACAGGTGCACCGCCTATAAAAAAAAAAAAAAAAAAAAAAAAGCGTGTTACCGCAATGTTCTGTCTCTTATACAACATCTAGATGTGTATAATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT",
            excludeNM);

        CPPUNIT_ASSERT_EQUAL(std::string("chr1,130,+,47M56S,60,0;chr1,48,+,47S56M,60,0;"), sa);
        CPPUNIT_ASSERT_EQUAL(-1U, excludeNM);
    }

    {
        const std::string sa = makeSaTag(
            -1U,
            "47M1C129B56M",
            "                                                                                                                                 TGTTACCGCAATGTTCTGTCTCTTATACAACATCTAGATGTGTATAAGAGACAGGTGCACCGCCTATACACATCTAGAATAAGAGACAGGTGCACCGCCTATA",
           //TGTTACCGCAATGTTCTGTCTCTTATACAACATCTAGATGTGTATAAGAGACAGGTGCACCGCCTATACACATCTAGAATAAGAGACAGGTGCACCGCCTATA"
            "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCGAGACAGGTGCACCGCCTATACACATCTAGAATAAGAGACAGGTGCACCGCCTATAAAAAAAAAAAAAAAAAAAAAAAAGCGTGTTACCGCAATGTTCTGTCTCTTATACAACATCTAGATGTGTATAATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT",
            "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCGAGACAGGTGCACCGCCTATACACATCTAGAATAAGAGACAGGTGCACCGCCTATAAAAAAAAAAAAAAAAAAAAAAAAGCGTGTTACCGCAATGTTCTGTCTCTTATACAACATCTAGATGTGTATAATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT",
            excludeNM);

        CPPUNIT_ASSERT_EQUAL(std::string("chr1,130,+,47M56S,60,0;chr2,48,+,47S56M,60,0;"), sa);
        CPPUNIT_ASSERT_EQUAL(-1U, excludeNM);
    }

    {
        const std::string sa = makeSaTag(
            -1U,
            "47S52M0F1C47B52S47M",
                  "TAGAGATGTTAAATATTAACTAAATTGTTACCACATCATATTGTTAGAGAATATTTAAAAGAAAAATTATTGAATAGTGCGAGAACCAAAAAGGCAAAA",
                  ///////////////////////////////////////////////|
//          "TTTTGCCTTTTTGGTTCTCGCACTATTCAATAATTTTTCTTTTAAATATTCTCTAACAATATGATGTGGTAACAATTTAGTTAATATTTAACATCTCTA",

                  "TTTCATTTCTTTTCACCGTCTCGGGGTTTAGTCAGTTATTTTAACTTAGAATATTTAAAAGAAAAATTATTGAATAGTGCGAGAACCAAAAAGGCAAAA",
            "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACTAACAATATGATGTGGTAACAATTTAGTTAATATTTAACATCTCTA",
            excludeNM
            );

        CPPUNIT_ASSERT_EQUAL(std::string("chr1,48,+,47S52M,60,0;chr2,53,-,52S47M,60,0;"), sa);
        CPPUNIT_ASSERT_EQUAL(-1U, excludeNM);
    }


    {
        const std::string sa = makeSaTag(
            -1U,
            "35M1C35B65M",
            "GATTTGGTAAAGTTTCTAGTGTCATGGGTTACTCTGTCTGGGGGGTATGCACGCGATAGCATTGCGAGACGCTGGAGCCGGAGCACCCTATGTCGCAGTA",
            "GATTTGGTAAAGTTTCTAGTGTCATGGGTTACTCT",
                                               "GTCTGGGGGGTGTGCACGCGATAGCATTGCGAGACGCTGGAGCCGGAGCACCCTATGTCGCAGTA",
            excludeNM
            );

        CPPUNIT_ASSERT_EQUAL(std::string("chr1,1,+,35M65S,60,0;chr2,1,+,35S65M,60,1;"), sa);
        CPPUNIT_ASSERT_EQUAL(-1U, excludeNM);
    }

    {
        const std::string sa = makeSaTag(
            isaac::reference::ReferencePosition(1,0),
            false,
            "35S65M",
            -1U,
            "35M1C35B65M",
            "GATTTGGTAAAGTTTCTAGTGTCATGGGTTACTCTGTCTGGGGGGTATGCACGCGATAGCATTGCGAGACGCTGGAGCCGGAGCACCCTATGTCGCAGTA",
            "GATTTGGTAAAGTTTCTAGTGTCATGGGTTACTCT",
                                               "GTCTGGGGGGTGTGCACGCGATAGCATTGCGAGACGCTGGAGCCGGAGCACCCTATGTCGCAGTA",
            excludeNM
            );

        CPPUNIT_ASSERT_EQUAL(std::string("chr1,1,+,35M65S,60,0;"), sa);
        CPPUNIT_ASSERT_EQUAL(1U, excludeNM);
    }
}


