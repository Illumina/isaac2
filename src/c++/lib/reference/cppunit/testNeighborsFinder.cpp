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

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <boost/foreach.hpp>
#include <boost/assign.hpp>

#include <boost/assign/std/vector.hpp>
using namespace boost::assign;

using namespace std;

#include "RegistryName.hh"
#include "reference/Contig.hh"
#include "reference/NeighborsFinder.hh"
#include "reference/neighborsFinder/NeighborCounter.hh"

#include "testNeighborsFinder.hh"

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION( TestNeighborsFinder, registryName("NeighborsFinder"));

void TestNeighborsFinder::setUp()
{
}

void TestNeighborsFinder::tearDown()
{
}

// This is needed for cases when isaac::oligo::Kmer is defined as __uint128_t or else CPPUNIT_ASSERT_EQUAL fails to compile
inline std::ostream & operator <<(std::ostream &os, const __uint128_t &kmer)
{
    return isaac::common::traceHexValue(os, kmer);
}

inline std::vector<char> vectorFromString(const std::string &str)
{
    return std::vector<char>(str.begin(), str.end());
}

namespace isaac
{
namespace reference
{
template <typename KmerT>
class TestNeighborCounter: boost::noncopyable
{
public:
    TestNeighborCounter(
        const unsigned int maskWidth,
        const unsigned mask,
        const unsigned neighborhoodWidth,
        const unsigned repeatThreshold,
        const std::vector<char> &reference,
        const reference::SortedReferenceMetadata &sortedReferenceMetadata
        ) :
            neighborhoodWidth_(neighborhoodWidth),
            reference_(reference),
            sortedReferenceMetadata_(sortedReferenceMetadata),
            contigOffsets_(computeContigOffsets(sortedReferenceMetadata.getContigs())),
            // as we're not going to try all positions against all other positions, set 0 for pairs that will never match by prefix,
            annotation_(reference::genomeLength(sortedReferenceMetadata.getContigs()), 0)
    {
    }

    typedef neighborsFinder::AnnotatedKmer<KmerT, NeighborsCount> AnnotatedKmer;
    typedef std::vector<AnnotatedKmer > KmerList;
    void getKmers(const oligo::Permutate &permutate, KmerList &kmerList, common::ThreadVector &threads)
    {
        isaac::oligo::KmerGenerator<32, KmerT> kmerGenerator(reference_.begin(), reference_.end());
        KmerT kmer(0);
        while (reference_.end() != kmerGenerator.next(kmer))
        {
            kmerList.push_back(AnnotatedKmer(permutate(kmer), NeighborsCount(0)));
//            ISAAC_THREAD_CERR << "getKmers" << kmerList.back() << std::endl;
        }
    }

    bool update(
        const unsigned kmerMismatchCount,
        AnnotatedKmer &one,
        AnnotatedKmer &another)
    {
        ISAAC_ASSERT_MSG(!another.isReverseComplement(), "second of the two must be a forward-facing kmer: " << one << " " << another);

        // Exact equivalence is required as we produce each k annotation separately from others
        if (neighborhoodWidth_ == kmerMismatchCount)
        {
            if (!one.isReverseComplement())
            {
                one.annotation_ += another.repeats_;
            }

            another.annotation_ += one.repeats_;
            return true;
        }

//        ISAAC_THREAD_CERR << "update:" << one << "," << another << std::endl;
//        ISAAC_THREAD_CERR << "update kmerMismatchCount:" << kmerMismatchCount << std::endl;
//        ISAAC_THREAD_CERR << "update neighborhoodWidth_:" << neighborhoodWidth_ << std::endl;

        return false;
    }

    typedef std::vector<NeighborsCount> Annotation;
    template <typename ReferenceKmerT>
    void updateAnnotation(const KmerList &kmerList, const bool includeRepeats)
    {
        ISAAC_THREAD_CERR << "updateAnnotation start" << std::endl;
        KmerList sortedKmers(kmerList.begin(), kmerList.end());
        std::sort(sortedKmers.begin(), sortedKmers.end(), AnnotatedKmer::kmerLess);
        int i = 0;
        BOOST_FOREACH(const AnnotatedKmer &ak, sortedKmers)
        {
            annotation_.at(i++) += ak.annotation_;
//            ISAAC_THREAD_CERR << "updateAnnotation:" << ak << " " << annotation_.at(i-1) << std::endl;
        }
    }
    const Annotation &getAnnotation() const {return annotation_;}

private:
    const unsigned neighborhoodWidth_;
    const std::vector<char> &reference_;
    const reference::SortedReferenceMetadata &sortedReferenceMetadata_;
    const std::vector<unsigned long> contigOffsets_;

    Annotation annotation_;

    template <bool countUnique>
    void markRepeats(
        const std::vector<unsigned long> &contigOffsets,
        typename KmerList::iterator begin,
        typename KmerList::iterator end,
        Annotation &annotation) const
    {;}

};
}
}

void TestNeighborsFinder::testFindNeighbors()
{
    const std::vector<char> reference = vectorFromString(
        "AAAAAAAAAATAAAAAAAAAAAAAAAAAAAAA" //first 32-mer
        "CAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACC");

    isaac::reference::SortedReferenceMetadata sortedReference;
    sortedReference.putContig(0, "vassilij", "blah.txt",
                              0, reference.size(), reference.size(), reference.size(), 0, 0, "", "", "");

    isaac::common::ThreadVector threads(2);
    typedef isaac::reference::TestNeighborCounter<isaac::oligo::KmerType> NeighborCounter;
    const unsigned neighborhoodWidth = 1;
    NeighborCounter neighborCounter(0, 0, neighborhoodWidth, 1000, reference, sortedReference);
    isaac::reference::NeighborsFinder<isaac::oligo::KmerType, NeighborCounter> finder(neighborCounter, true, neighborhoodWidth, threads);
    const std::vector<isaac::reference::NeighborsCount> &annotation = finder.annotate<isaac::oligo::KmerType>();

    CPPUNIT_ASSERT_EQUAL(isaac::reference::NeighborsCount(2), annotation.at(0));   //AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAC
    CPPUNIT_ASSERT_EQUAL(isaac::reference::NeighborsCount(1), annotation.at(1));   //AAAAAAAAAAAAAAAAAAAAAAAAAAAAAACC
    CPPUNIT_ASSERT_EQUAL(isaac::reference::NeighborsCount(0), annotation.at(2));   //AAAAAAAAAAAAAAAAAAAAACAAAAAAAAAA
    CPPUNIT_ASSERT_EQUAL(isaac::reference::NeighborsCount(0), annotation.at(3));   //AAAAAAAAAAAAAAAAAAAACAAAAAAAAAAA
    CPPUNIT_ASSERT_EQUAL(isaac::reference::NeighborsCount(0), annotation.at(4));   //AAAAAAAAAAAAAAAAAAACAAAAAAAAAAAA
    CPPUNIT_ASSERT_EQUAL(isaac::reference::NeighborsCount(0), annotation.at(5));   //AAAAAAAAAAAAAAAAAACAAAAAAAAAAAAA
    CPPUNIT_ASSERT_EQUAL(isaac::reference::NeighborsCount(0), annotation.at(6));   //AAAAAAAAAAAAAAAAACAAAAAAAAAAAAAA
    CPPUNIT_ASSERT_EQUAL(isaac::reference::NeighborsCount(0), annotation.at(7));   //AAAAAAAAAAAAAAAACAAAAAAAAAAAAAAA
    CPPUNIT_ASSERT_EQUAL(isaac::reference::NeighborsCount(0), annotation.at(8));   //AAAAAAAAAAAAAAACAAAAAAAAAAAAAAAA
    CPPUNIT_ASSERT_EQUAL(isaac::reference::NeighborsCount(0), annotation.at(9));   //AAAAAAAAAAAAAACAAAAAAAAAAAAAAAAA
    CPPUNIT_ASSERT_EQUAL(isaac::reference::NeighborsCount(0), annotation.at(10));  //AAAAAAAAAAAAACAAAAAAAAAAAAAAAAAA
    CPPUNIT_ASSERT_EQUAL(isaac::reference::NeighborsCount(0), annotation.at(11));  //AAAAAAAAAAAACAAAAAAAAAAAAAAAAAAA
    CPPUNIT_ASSERT_EQUAL(isaac::reference::NeighborsCount(0), annotation.at(12));  //AAAAAAAAAAACAAAAAAAAAAAAAAAAAAAA
    CPPUNIT_ASSERT_EQUAL(isaac::reference::NeighborsCount(1), annotation.at(13));  //AAAAAAAAAACAAAAAAAAAAAAAAAAAAAAA
    CPPUNIT_ASSERT_EQUAL(isaac::reference::NeighborsCount(1), annotation.at(14));  //AAAAAAAAAATAAAAAAAAAAAAAAAAAAAAA
    CPPUNIT_ASSERT_EQUAL(isaac::reference::NeighborsCount(0), annotation.at(15));  //AAAAAAAAACAAAAAAAAAAAAAAAAAAAAAA
    CPPUNIT_ASSERT_EQUAL(isaac::reference::NeighborsCount(1), annotation.at(16));  //AAAAAAAAATAAAAAAAAAAAAAAAAAAAAAC
    CPPUNIT_ASSERT_EQUAL(isaac::reference::NeighborsCount(0), annotation.at(17));  //AAAAAAAACAAAAAAAAAAAAAAAAAAAAAAA
    CPPUNIT_ASSERT_EQUAL(isaac::reference::NeighborsCount(0), annotation.at(18));  //AAAAAAAATAAAAAAAAAAAAAAAAAAAAACA
    CPPUNIT_ASSERT_EQUAL(isaac::reference::NeighborsCount(0), annotation.at(19));  //AAAAAAACAAAAAAAAAAAAAAAAAAAAAAAA
    CPPUNIT_ASSERT_EQUAL(isaac::reference::NeighborsCount(0), annotation.at(20));  //AAAAAAATAAAAAAAAAAAAAAAAAAAAACAA
    CPPUNIT_ASSERT_EQUAL(isaac::reference::NeighborsCount(0), annotation.at(21));  //AAAAAACAAAAAAAAAAAAAAAAAAAAAAAAA
    CPPUNIT_ASSERT_EQUAL(isaac::reference::NeighborsCount(0), annotation.at(22));  //AAAAAATAAAAAAAAAAAAAAAAAAAAACAAA
    CPPUNIT_ASSERT_EQUAL(isaac::reference::NeighborsCount(0), annotation.at(23));  //AAAAACAAAAAAAAAAAAAAAAAAAAAAAAAA
    CPPUNIT_ASSERT_EQUAL(isaac::reference::NeighborsCount(0), annotation.at(24));  //AAAAATAAAAAAAAAAAAAAAAAAAAACAAAA
    CPPUNIT_ASSERT_EQUAL(isaac::reference::NeighborsCount(0), annotation.at(25));  //AAAACAAAAAAAAAAAAAAAAAAAAAAAAAAA
    CPPUNIT_ASSERT_EQUAL(isaac::reference::NeighborsCount(0), annotation.at(26));  //AAAATAAAAAAAAAAAAAAAAAAAAACAAAAA
    CPPUNIT_ASSERT_EQUAL(isaac::reference::NeighborsCount(0), annotation.at(27));  //AAACAAAAAAAAAAAAAAAAAAAAAAAAAAAA
    CPPUNIT_ASSERT_EQUAL(isaac::reference::NeighborsCount(0), annotation.at(28));  //AAATAAAAAAAAAAAAAAAAAAAAACAAAAAA
    CPPUNIT_ASSERT_EQUAL(isaac::reference::NeighborsCount(0), annotation.at(29));  //AACAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
    CPPUNIT_ASSERT_EQUAL(isaac::reference::NeighborsCount(0), annotation.at(30));  //AATAAAAAAAAAAAAAAAAAAAAACAAAAAAA
    CPPUNIT_ASSERT_EQUAL(isaac::reference::NeighborsCount(0), annotation.at(31));  //ACAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
    CPPUNIT_ASSERT_EQUAL(isaac::reference::NeighborsCount(0), annotation.at(32));  //ATAAAAAAAAAAAAAAAAAAAAACAAAAAAAA
    CPPUNIT_ASSERT_EQUAL(isaac::reference::NeighborsCount(0), annotation.at(33));  //CAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
    CPPUNIT_ASSERT_EQUAL(isaac::reference::NeighborsCount(0), annotation.at(34));  //TAAAAAAAAAAAAAAAAAAAAACAAAAAAAAA
}


