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
 **
 ** \file BandedSmithWaterman.cpp
 **
 ** \brief See BandedSmithWaterman.hh
 **
 ** \author Come Raczy
 **/

#include <iostream>
#include <cassert>
#include <cstdlib>
#include <algorithm>
#include <xmmintrin.h> 
#if __AVX2__
#include <immintrin.h>
#endif
#include <boost/format.hpp>

#include "alignment/BandedSmithWaterman.hh"

namespace isaac
{
namespace alignment
{

BandedSmithWaterman::BandedSmithWaterman(const int matchScore, const int mismatchScore,
                                         const int gapOpenScore, const int gapExtendScore,
                                         const int maxReadLength)
    : matchScore_(matchScore)
    , mismatchScore_(mismatchScore)
    , gapOpenScore_(gapOpenScore)
    , gapExtendScore_(gapExtendScore)
    , maxReadLength_(maxReadLength)
    , initialValue_(static_cast<int>(std::numeric_limits<short>::min()) + gapOpenScore_)
#if __AVX2__
    , T_((char *)_mm_malloc (maxReadLength_ * 3 * sizeof(__m256i), 32))
#else
    , T_((char *)_mm_malloc (maxReadLength_ * 3 * sizeof(__m128i), 16))
#endif
{
    // check that there won't be any overflows in the matrices
    const int maxScore = std::max(std::max(std::max(abs(matchScore_), abs(mismatchScore_)), abs(gapOpenScore_)), abs(gapExtendScore_));
    if ((maxReadLength_ * maxScore) >= abs(static_cast<int>(initialValue_)))
    {
        const std::string message = (boost::format("BandedSmithWaterman: unsupported read length (%i) for these scores (%i): use smaller scores or shorter reads") % maxReadLength_ % maxScore).str();
        BOOST_THROW_EXCEPTION(isaac::common::InvalidParameterException(message));
    }
}

BandedSmithWaterman::~BandedSmithWaterman()
{
    _mm_free(T_);
}

// insert in register 0 only -- workaround for missing sse4 instruction set
inline __m128i _mm_insert_epi8(__m128i v, char c, int)
{
    const unsigned int tmp = _mm_extract_epi16(v, 0) & 0xffffff00u;
    const unsigned s = c & 0xff;
    return _mm_insert_epi16(v, tmp | s, 0);
}

// helper functions for debugging
std::string epi8(__m128i v);
std::string epi8s(__m128i v);
std::string epi8c(__m128i v);
std::ostream &operator<<(std::ostream &os, const __m128i &H);

#if __AVX2__
inline void BandedSmithWaterman::max3_(__m256i E, __m256i G, __m256i F, __m256i *output, __m256i *outputType) const {
    __m256i cmpgtEg = _mm256_cmpgt_epi16(E, G);
    __m256i cmpgtEgMask = _mm256_srli_epi16(cmpgtEg, 15);

    __m256i maxEg = _mm256_max_epi16(G, E);

    __m256i cmpgtGf = _mm256_cmpgt_epi16(F, maxEg);
    __m256i cmpgtGfMask =_mm256_and_si256(cmpgtGf, _mm256_set1_epi16(2));

    __m256i maxGf = _mm256_max_epi16(maxEg, F);

    __m256i maxEgf = _mm256_max_epi16(cmpgtEgMask, cmpgtGfMask);
 
    _mm256_store_si256(outputType, maxEgf);
    _mm256_store_si256(output, maxGf);
}
  
inline void BandedSmithWaterman::max3_(__m256i E, __m256i G, __m256i F, __m256i gapOpen, __m256i gapExtend, __m256i *output, __m256i *outputType, int initialValue) const {
    __m256i cmpgtEg = _mm256_cmpgt_epi16(E, G);
    __m256i cmpgtEgMask = _mm256_and_si256(cmpgtEg, _mm256_set1_epi16(1));

    __m256i maxEg = _mm256_max_epi16(G, E);
    __m256i maxEgSubGapOpen = _mm256_sub_epi16(maxEg, gapOpen);

    __m256i cmpgtGf = _mm256_cmpgt_epi16(F, maxEgSubGapOpen);
    __m256i cmpgtGfMask =_mm256_and_si256(cmpgtGf, _mm256_set1_epi16(2));

    __m256i maxGf = _mm256_max_epi16(maxEgSubGapOpen, F);
    __m256i maxGfReset = _mm256_insert_epi16(maxGf, initialValue, 0);

    __m256i maxEgf = _mm256_max_epi16(cmpgtEgMask, cmpgtGfMask);
    __m256i maxEgfReset = _mm256_insert_epi16(maxEgf, 0, 0);

    _mm256_store_si256(outputType, maxEgfReset);
    _mm256_store_si256(output, maxGfReset);
}

inline void BandedSmithWaterman::max2_(__m256i G, __m256i F, int gapExtendScore, __m256i *output, __m256i *outputType, int initialValue) const {
    __m256i cmpgtFg = _mm256_cmpgt_epi16(F, G);
    __m256i cmpgtFgMask = _mm256_srli_epi16(cmpgtFg, 15);
    cmpgtFgMask = _mm256_slli_epi16(cmpgtFgMask, 1);

    __m256i cmpgtFgMaskSr = _mm256_alignr_epi8(_mm256_permute2x128_si256(cmpgtFgMask, cmpgtFgMask, _MM_SHUFFLE(2, 0, 0, 1)), cmpgtFgMask, 2);

    __m256i maxFg = _mm256_max_epi16(F, G);
    __m256i maxFgSr = _mm256_alignr_epi8(_mm256_permute2x128_si256(maxFg, maxFg, _MM_SHUFFLE(2, 0, 0, 1)), maxFg, 2);
    __m256i maxFgSrReset = _mm256_insert_epi16(maxFgSr, initialValue, 15);

    int16_t maxFgS[WIDEST_GAP_SIZE];
    int16_t maxFgSueFgS[WIDEST_GAP_SIZE];
    _mm256_store_si256((__m256i *) maxFgS, maxFg);

    short e = initialValue;
    short fg = initialValue;

    for (unsigned int j = WIDEST_GAP_SIZE; j > 0; --j)
    {
        short max = fg;
        if (e - gapExtendScore > fg)
        {
            max = e - gapExtendScore;
        }

        maxFgSueFgS[j - 1] = max;
        fg = maxFgS[j - 1];
    
        e = max;
    }

    __m256i maxFgSueFg = _mm256_load_si256((__m256i *) maxFgSueFgS);

    __m256i cmpgtFgSueFg = _mm256_cmpgt_epi16(maxFgSueFg, maxFgSrReset);
    __m256i cmpgtFgSueFgMask = _mm256_and_si256(cmpgtFgSueFg, _mm256_set1_epi16(5));

    maxFgSueFg = _mm256_max_epi16(maxFgSueFg, maxFgSrReset);
    __m256i cmpgtMerge = _mm256_and_si256(_mm256_max_epi16(cmpgtFgSueFgMask, cmpgtFgMaskSr), _mm256_set1_epi16(3));

    _mm256_store_si256(outputType, cmpgtMerge);
    _mm256_store_si256(output, maxFgSueFg);
}

#endif

unsigned BandedSmithWaterman::align(
    const std::vector<char> &query,
    const std::vector<char>::const_iterator databaseBegin,
    const std::vector<char>::const_iterator databaseEnd,
    Cigar &cigar) const
{
    return align(query.begin(), query.end(), databaseBegin, databaseEnd, cigar);
}

unsigned BandedSmithWaterman::trimTailIndels(Cigar& cigar, const size_t beginOffset) const
{
    unsigned ret = 0;
    unsigned long extend = 0;
    for (Cigar::Component component = Cigar::decode(cigar.back());
        cigar.size() != beginOffset; component = Cigar::decode(cigar.back()))
    {
        if (Cigar::DELETE == component.second)
        {
            //CASAVA does not like CIGAR beginning with a deletion in the data
            cigar.pop_back();
            ISAAC_ASSERT_MSG(
                Cigar::DELETE != Cigar::decode(cigar.back()).second,
                "two Cigar::DELETE cannot be next to each other");
            ret += component.first;
        }
        else if (Cigar::INSERT == component.second)
        {
            // tail and head insertions are biasing best alignment choice. remove them and extend the adjacent align operation
            extend += component.first;
            ret -= component.first;
            cigar.pop_back();
        }
        else
        {
            break;
        }
    }

    if (extend)
    {
        if(cigar.size() == beginOffset)
        {
            // this was a really bad s-w alignment. Something like this: 15D7I7D15I15D15I15D15I15D15I15D
            cigar.addOperation(extend, Cigar::ALIGN);
            ISAAC_ASSERT_MSG(!extend, "Unexpectedly long CIGAR operation");
        }
        else
        {
            const Cigar::Component component = Cigar::decode(cigar.back());
            ISAAC_ASSERT_MSG(Cigar::ALIGN == component.second, "Unexpected operation at the end of cigar " <<
                             Cigar::toString(cigar.begin() + beginOffset, cigar.end()) << " : "  << component.second);
            cigar.updateOperation(cigar.size() - 1, component.first + extend, Cigar::ALIGN);
        }
    }
    return ret;
}

void BandedSmithWaterman::removeAdjacentIndels(Cigar& cigar, const size_t beginOffset) const
{
//    ISAAC_THREAD_CERR << "BandedSmithWaterman::removeAdjacentIndels " << Cigar::toString(cigar.begin() + beginOffset, cigar.end()) << std::endl;
    const Cigar::iterator begin = cigar.begin() + beginOffset;
    for (Cigar::iterator it = begin + 1; cigar.end() != it;)
    {
        const Cigar::iterator last = it - 1;
        const Cigar::Component lastcomponent = Cigar::decode(*last);
        const Cigar::Component component = Cigar::decode(*it);
        if (Cigar::ALIGN == lastcomponent.second)
        {
            if (Cigar::ALIGN == component.second)
            {
                cigar.updateOperation(std::distance(cigar.begin(), last), lastcomponent.first + component.first, Cigar::ALIGN);
                it = cigar.erase(it);
            }
            else
            {
                ++it;
            }
        }
        else if (Cigar::DELETE == lastcomponent.second)
        {
            if (Cigar::DELETE == component.second)
            {
                cigar.updateOperation(std::distance(cigar.begin(), last), lastcomponent.first + component.first, Cigar::DELETE);
                it = cigar.erase(it);
            }
            else if (Cigar::INSERT == component.second)
            {
                if (component.first < lastcomponent.first)
                {
                    cigar.updateOperation(std::distance(cigar.begin(), last), lastcomponent.first - component.first, Cigar::DELETE);
                    cigar.updateOperation(std::distance(cigar.begin(), it), component.first, Cigar::ALIGN);
                }
                else if (lastcomponent.first < component.first)
                {
                    cigar.updateOperation(std::distance(cigar.begin(), last), component.first - lastcomponent.first, Cigar::INSERT);
                    cigar.updateOperation(std::distance(cigar.begin(), it), lastcomponent.first, Cigar::ALIGN);
                    if (begin != it - 1)
                    {
                        --it;
                    }
                }
                else
                {
                    cigar.updateOperation(std::distance(cigar.begin(), last), component.first, Cigar::ALIGN);
                    it = cigar.erase(it);
                    if (begin != it - 1)
                    {
                        --it;
                    }
                }
            }
            else
            {
                ++it;
            }
        }
        else if (Cigar::INSERT == lastcomponent.second)
        {
            if (Cigar::INSERT == component.second)
            {
                cigar.updateOperation(std::distance(cigar.begin(), last), lastcomponent.first + component.first, Cigar::INSERT);
                it = cigar.erase(it);
            }
            else if (Cigar::DELETE == component.second)
            {
                if (component.first < lastcomponent.first)
                {
                    cigar.updateOperation(std::distance(cigar.begin(), last), lastcomponent.first - component.first, Cigar::INSERT);
                    cigar.updateOperation(std::distance(cigar.begin(), it), component.first, Cigar::ALIGN);
                }
                else if (lastcomponent.first < component.first)
                {
                    cigar.updateOperation(std::distance(cigar.begin(), last), component.first - lastcomponent.first, Cigar::DELETE);
                    cigar.updateOperation(std::distance(cigar.begin(), it), lastcomponent.first, Cigar::ALIGN);
                    if (begin != it - 1)
                    {
                        --it;
                    }
                }
                else
                {
                    cigar.updateOperation(std::distance(cigar.begin(), last), component.first, Cigar::ALIGN);
                    it = cigar.erase(it);
                    if (begin != it - 1)
                    {
                        --it;
                    }
                }
            }
            else
            {
                ++it;
            }
        }
        else
        {
            ++it;
        }
    }
    //ISAAC_THREAD_CERR << "BandedSmithWaterman::removeAdjacentIndels done " << Cigar::toString(cigar.begin() + beginOffset, cigar.end()) << std::endl;
}


unsigned BandedSmithWaterman::align(
    const std::vector<char>::const_iterator queryBegin,
    const std::vector<char>::const_iterator queryEnd,
    const std::vector<char>::const_iterator databaseBegin,
    const std::vector<char>::const_iterator databaseEnd,
    Cigar &cigar) const
{
    assert(databaseEnd > databaseBegin);
    const size_t querySize = std::distance(queryBegin, queryEnd);
    ISAAC_ASSERT_MSG(querySize + WIDEST_GAP_SIZE - 1 == (unsigned long)(databaseEnd - databaseBegin), "q:" << std::string(queryBegin, queryEnd) << " db:" << std::string(databaseBegin, databaseEnd));
    assert(querySize <= size_t(maxReadLength_));
    const size_t originalCigarSize = cigar.size();
#if __AVX2__
    __m256i *t = (__m256i *)T_;
    const __m256i GapOpenScore = _mm256_set1_epi16(gapOpenScore_);
    const __m256i GapExtendScore = _mm256_set1_epi16(gapExtendScore_);
    const __m256i ONE = _mm256_set1_epi16(0xFF);
    // Initialize E, F and G
    __m256i E, F, G;
    E = _mm256_set1_epi16(initialValue_);
    F = _mm256_set1_epi16(0); // Should this be -10000?
    G = _mm256_set1_epi16(initialValue_);
    G = _mm256_insert_epi16(G, 0, 0);
    // initialize D -- Note that we must leave the leftmost position empty (else it will be discarded before use)
    __m256i D = _mm256_setzero_si256();
    int16_t dbBuffer[WIDEST_GAP_SIZE];
    std::reverse_copy(databaseBegin, databaseBegin + WIDEST_GAP_SIZE - 1, dbBuffer);
    dbBuffer[WIDEST_GAP_SIZE - 1] = dbBuffer[0];
    D = _mm256_load_si256((__m256i *)&dbBuffer);
#else

    __m128i *t = (__m128i *)T_;
    const __m128i GapOpenScore = _mm_set1_epi16(gapOpenScore_);
    const __m128i GapExtendScore = _mm_set1_epi16(gapExtendScore_);
    // Initialize E, F and G
    __m128i E[2], F[2], G[2];
    for (unsigned int i = 0; 2 > i; ++i)
    {
        E[i] = _mm_set1_epi16(initialValue_);
        F[i] = _mm_set1_epi16(0); // Should this be -10000?
        G[i] = _mm_set1_epi16(initialValue_);
    }
    G[0] = _mm_insert_epi16(G[0], 0, 0);
    // initialize D -- Note that we must leave the leftmost position empty (else it will be discarded before use)
    __m128i D = _mm_setzero_si128();
    for (unsigned int i = 0; WIDEST_GAP_SIZE > i + 1; ++i)
    {
        D = _mm_slli_si128(D, 1);
        D = _mm_insert_epi8(D, *(databaseBegin + i), 0);
    }
#endif
    // iterate over all bases in the query
    std::vector<char>::const_iterator queryCurrent = queryBegin;
    for (unsigned queryOffset = 0; queryEnd != queryCurrent; ++queryOffset, ++queryCurrent)
    {

#if __AVX2__
        __m256i G1, E1;
        // get F[i-1, j] - extend
        __m256i F1;
        __m256i TF;
        __m256i newF;

        __m256i TG;
        __m256i newG;

        G1 = _mm256_alignr_epi8(G, _mm256_permute2x128_si256(G, G, _MM_SHUFFLE(0, 0, 2, 0)), 14);
        E1 = _mm256_alignr_epi8(E, _mm256_permute2x128_si256(E, E, _MM_SHUFFLE(0, 0, 2, 0)), 14);
        F1 = _mm256_alignr_epi8(F, _mm256_permute2x128_si256(F, F, _MM_SHUFFLE(0, 0, 2, 0)), 14);
        F1 = _mm256_sub_epi16(F1, GapExtendScore);

        max3_(E1, G1, F1, GapOpenScore, GapExtendScore, &newF, &TF, initialValue_);
        max3_(E, G, F, &newG, &TG);

        __m256i TE = _mm256_setzero_si256();

        // add the match/mismatch score
        // load the query base in all 8 values of the register
        __m256i Q = _mm256_set1_epi16(*queryCurrent);
        // shift the database by 1 byte to the left and add the new base
        D = _mm256_alignr_epi8(D, _mm256_permute2x128_si256(D, D, _MM_SHUFFLE(0, 0, 2, 0)), 14);
        D = _mm256_insert_epi16(D, *(databaseBegin + queryOffset + WIDEST_GAP_SIZE - 1), 0);

        // compare query and database. 0xff if different (that also the sign bits)
        const __m256i B = _mm256_andnot_si256(_mm256_cmpeq_epi16(Q, D), ONE);

        // set match/mismatch scores, according to comparison
        const __m256i Match = _mm256_andnot_si256(B, _mm256_set1_epi16(matchScore_));
        const __m256i Mismatch = _mm256_and_si256(B, _mm256_set1_epi16(mismatchScore_));
        // add the match/mismatch scored to HH
        const __m256i W = _mm256_add_epi16(Match, Mismatch);

        newG = _mm256_add_epi16(newG, _mm256_slli_epi16(B, 8));
        newG = _mm256_add_epi16(newG, W);

        // E[i,j] = max(G[i, j-1] - open, E[i, j-1] - extend, F[i, j-1] - open)
        __m256i tmpNewG = _mm256_sub_epi16(newG, GapOpenScore);
        __m256i tmpNewF = _mm256_sub_epi16(newF, GapOpenScore);
        max2_(tmpNewG, tmpNewF, gapExtendScore_, &E, &TE, initialValue_);

        G = newG;
        F = newF;
        _mm256_store_si256(t++, TG);
        _mm256_store_si256(t++, TE);
        _mm256_store_si256(t++, TF);
#else
        __m128i tmp0[2], tmp1[2], tmp2[2];
        // F[i, j] = max(G[i-1, j] - open, E[i-1, j] - open, F[i-1, j] - extend)
        // G[i-1, j] and E[i-1, j]
        tmp0[1] = _mm_slli_si128(G[1], 2);
        tmp0[1] = _mm_insert_epi16(tmp0[1], _mm_extract_epi16(G[0], 7), 0);
        tmp0[0] = _mm_slli_si128(G[0], 2); // is the 0 initialisation alright?
        tmp1[1] = _mm_slli_si128(E[1], 2);
        tmp1[1] = _mm_insert_epi16(tmp1[1], _mm_extract_epi16(E[0], 7), 0);
        tmp1[0] = _mm_slli_si128(E[0], 2); // is the 0 initialisation alright?

        for (unsigned j = 0; 2 > j; ++j)
        {
            // identify which matrix provided the max (hack: true is 0xffff)
            tmp2[j] = _mm_cmplt_epi16(tmp0[j], tmp1[j]); // 0xffff if G[i-1, j] < E[i-1, j], 0 otherwise
            tmp2[j] = _mm_srli_epi16(tmp2[j], 15); // shift 15 bits to have only 1 for true
        }
        __m128i TF = _mm_packs_epi16(tmp2[0], tmp2[1]);
        // get max(G[i-1, j] - open, E[i-1, j] - open)
        __m128i newF[2];
        for (unsigned int j = 0; 2 > j; ++j)
        {
            newF[j] = _mm_max_epi16(tmp0[j], tmp1[j]);// newF = max(G[i-1, j], E[i-1, j])
            newF[j] = _mm_sub_epi16(newF[j], GapOpenScore); // newF = max(G[i-1, j] - open, E[i-1, j] - open)
        }
        // get F[i-1, j] - extend
        tmp0[1] = _mm_slli_si128(F[1], 2);
        tmp0[1] = _mm_insert_epi16(tmp0[1], _mm_extract_epi16(F[0], 7), 0);
        tmp0[1] = _mm_sub_epi16(tmp0[1], GapExtendScore);
        tmp0[0] = _mm_slli_si128(F[0], 2);
        tmp0[0] = _mm_sub_epi16(tmp0[0], GapExtendScore);
        // correct TF
        for (unsigned j = 0; 2 > j; ++j)
        {
            tmp2[j] = _mm_cmplt_epi16(newF[j], tmp0[j]); // 0xffff if max(G[i-1, j], E[i-1,j] - open < F[i-1, j] - extend
            //tmp2[j] = _mm_srli_epi16(tmp2[j], 14); // shift 14 bits to have 3 for true
            tmp2[j] = _mm_slli_epi16(_mm_srli_epi16(tmp2[j], 15), 1); // shift right 15 bits and left 1 bit to have 2 for true
        }
        TF = _mm_max_epu8(_mm_packs_epi16(tmp2[0], tmp2[1]), TF); // 0, 1, or 2 for G, E or F
        TF = _mm_insert_epi8(TF, 0, 0); // 0, 1, or 2 for G, E or F
        // correct F according to (F[i-1, j] - extend)
        for (unsigned int j = 0; 2 > j; ++j)
        {
            newF[j] = _mm_max_epi16(newF[j], tmp0[j]);
        }
        newF[0] = _mm_insert_epi16(newF[0], initialValue_, 0);
        // G[i, j] = max(G[i-1, j-1], E[i-1, j-1], F[i-1, j-1]
        __m128i newG[2];
        for (unsigned int j = 0; 2 > j; ++j)
        {
            tmp2[j] = _mm_cmplt_epi16(G[j], E[j]); // 0xffff if G[i-1, j-1] < E[i-1, j-1], 0 otherwise
            tmp2[j] = _mm_srli_epi16(tmp2[j], 15); // shift 15 bits to have only 1 for true
            newG[j] = _mm_max_epi16(G[j], E[j]);// newG = max(G[i-1, j-1], E[i-1, j-1])
        }
        __m128i TG = _mm_packs_epi16(tmp2[0], tmp2[1]);
        // correct G and TG
        //std::cerr << "newG: " << newG[1] << newG[0] << std::endl;
        for (unsigned int j = 0; 2 > j; ++j)
        {
            tmp2[j] =  _mm_cmplt_epi16(newG[j], F[j]); // 0xffff if max(G[i-1, j-1], E[i-1,j-1] < F[i-1, j-1]
            tmp2[j] = _mm_slli_epi16(_mm_srli_epi16(tmp2[j], 15), 1); // shift right 15 bits and left 1 bit to have 2 for true
            newG[j] = _mm_max_epi16(newG[j], F[j]);
        }
        //std::cerr << "newG: " << newG[1] << newG[0] << std::endl;
        //std::cerr << "   G: " << G[1] << G[0] << std::endl;
        //std::cerr << "   E: " << E[1] <<  E[0] << std::endl;
        //std::cerr << "   F: " << F[1] <<  F[0] << std::endl;
        //std::cerr << "tmp0: " << tmp0[1] << tmp0[0] << std::endl;
        //std::cerr << "tmp1: " << tmp1[1] << tmp1[0] << std::endl;
        TG = _mm_max_epu8(_mm_packs_epi16(tmp2[0], tmp2[1]), TG); // 0, 1, or 2 for G, E or F
        // add the match/mismatch score
        // load the query base in all 8 values of the register
        __m128i Q = _mm_set1_epi8(*queryCurrent);
        // shift the database by 1 byte to the left and add the new base
        D = _mm_slli_si128(D, 1);
        D = _mm_insert_epi8(D, *(databaseBegin + queryOffset + WIDEST_GAP_SIZE - 1), 0);
        // compare query and database. 0xff if different (that also the sign bits)
        const __m128i ONE = _mm_set1_epi8(0xFF);
        const __m128i B = _mm_andnot_si128(_mm_cmpeq_epi8(Q, D), ONE);

#if 0
        //std::cerr << std::endl << database << std::endl << query << std::endl;;
        std::cerr << (boost::format("%2d  D:") % (queryOffset+1)).str();
        for (unsigned j = 0; j < queryOffset; ++j)
        {
            std::cerr << " 0";
        }
        std::cerr << epi8c(D) << std::endl;

        std::cerr << (boost::format("%2d  Q:") % (queryOffset+1)).str();
        for (unsigned j = 0; j < queryOffset; ++j)
        {
            std::cerr << " 0";
        }
        std::cerr << epi8c(Q) << std::endl;
        //std::cerr << (boost::format("%2d  B:") % (i+1)).str();
        //for (unsigned j = 0; j < i; ++j)
        //{
        //    std::cerr << "   0";
        //}
        //std::cerr << epi8(B) << std::endl;
        //exit(1);
#endif
        // set match/mismatch scores, according to comparison
        const __m128i Match = _mm_andnot_si128(B, _mm_set1_epi8(matchScore_));
        const __m128i Mismatch = _mm_and_si128(B, _mm_set1_epi8(mismatchScore_));
        // add the match/mismatch scored to HH
        const __m128i W = _mm_add_epi8(Match, Mismatch);

#if 0
        std::cerr << (boost::format("%2d  W:") % (queryOffset+1)).str();
        for (unsigned j = 0; j < queryOffset; ++j)
        {
            std::cerr << " 0";
        }
        std::cerr << epi8s(W) << std::endl;
#endif
        newG[0] = _mm_add_epi16(newG[0], _mm_unpacklo_epi8(W, B));
        newG[1] = _mm_add_epi16(newG[1], _mm_unpackhi_epi8(W, B));
        // E[i,j] = max(G[i, j-1] - open, E[i, j-1] - extend, F[i, j-1] - open)
        __m128i TE = _mm_setzero_si128();
        // E should never be the maximum in the leftmost side of the window
        short g = initialValue_; // -gapOpenScore_;
        short e = initialValue_; // -gapExtendScore_;
        short f = initialValue_; // -gapOpenScore_;
        for (unsigned j = 0; 2 > j; ++j)
        {
            tmp0[j] = _mm_sub_epi16(newG[j], GapOpenScore);
            tmp1[j] = _mm_sub_epi16(newF[j], GapOpenScore);
        }
        //std::cerr << "newG: " << newG[1] << newG[0] << std::endl;
        //std::cerr << "tmp0: " << tmp0[1] << tmp0[0] << std::endl;
        //std::cerr << "newF: " << newF[1] << newF[0] << std::endl;
        //std::cerr << "tmp1: " << tmp1[1] << tmp1[0] << std::endl;
        //std::cerr <<  "WIDEST_GAP_SIZE = " << WIDEST_GAP_SIZE << std::endl;
        for (unsigned int j = 0; j < WIDEST_GAP_SIZE; ++j)
        {
            //std::cerr << "---- j = " << j << std::endl;
            if (8 == j)
            {
                //std::cerr << "   E: " << E[1] << E[0] << std::endl;
                //std::cerr << "tmp0: " << tmp0[1] << tmp0[0] << std::endl;
                //std::cerr << "tmp1: " << tmp1[1] << tmp1[0] << std::endl;
                tmp0[1] = tmp0[0];
                tmp1[1] = tmp1[0];
                E[1] = E[0];
            }
            short max = g;
            short tMax = 0;
            if (e > g && e > f)
            {
                max = e;
                tMax = 1;
            }
            else if (f > g)
            {
                max = f;
                tMax = 2;
            }
            //std::cerr << (boost::format("g = %d, e = %d, f = %d, max = %d, tMax = %d") % g % e % f % max % tMax ).str() << std::endl;
            TE = _mm_slli_si128(TE, 1);
            TE = _mm_insert_epi8(TE, tMax, 0);
            E[0] = _mm_slli_si128(E[0], 2);
            E[0] = _mm_insert_epi16(E[0], max, 0);
            g = _mm_extract_epi16(tmp0[1], 7);
            e = max - gapExtendScore_;
            f = _mm_extract_epi16(tmp1[1], 7);
            tmp0[1] = _mm_slli_si128(tmp0[1], 2);
            tmp1[1] = _mm_slli_si128(tmp1[1], 2);
            //std::cerr << "   E: " << E[1] << E[0] << std::endl;
            //std::cerr << "  TE: " << epi8(TE) << std::endl;
        }
        //E[0] = _mm_slli_si128(E[0], 2);
        for (unsigned j = 0; 2 > j; ++j)
        {
            G[j] = newG[j];
            F[j] = newF[j];
        }
        // TODO: add support for databases shorter than query + widestGapSize
        // store the matrix types
        _mm_store_si128(t++, TG);
        _mm_store_si128(t++, TE);
        _mm_store_si128(t++, TF);
#endif
    }
    // find the max of E, F and G at the end
#if __AVX2__
    short max = _mm256_extract_epi16(G, 15) - 1;
#else
    short max = _mm_extract_epi16(G[1], 7) - 1;
#endif

    int ii = querySize - 1;
    int jj = ii;
    unsigned maxType = 0;
#if __AVX2__
    __m256i TT, TMax;
    max3_(E, G, F, &TMax, &TT);

    int16_t TMaxS[16], TTS[16];
    _mm256_store_si256((__m256i *) TMaxS, TMax);
    _mm256_store_si256((__m256i *) TTS, TT);

    for(unsigned j = 16; 0 < j; --j) {
        const short value = TMaxS[j - 1];
        if(value > max)
        {
            max = value;
            jj = j - 1;
            maxType = (short) TTS[j - 1];
        }
    }
#else
    __m128i *TT[] = {G, E, F};
    for (unsigned int k = 0; 2 > k; ++k)
    {
        const unsigned int kk = (k + 1) % 2;

        for (unsigned j = 8; 0 < j; --j)
        {
            for (unsigned type = 0; 3 > type; ++type)
            {
                const short value = _mm_extract_epi16(TT[type][kk], 7);
                TT[type][kk] = _mm_slli_si128(TT[type][kk], 2);
                if (value > max)
                {
                    max = value;
                    jj = 8 * kk + j - 1;
                    maxType = type;
                }
            }
        }
    }
#endif
    const int jjIncrement[] = {0, 1, -1};
    const int iiIncrement[] = {-1, 0, -1};
    const Cigar::OpCode opCodes[] = {Cigar::ALIGN, Cigar::DELETE, Cigar::INSERT};
    unsigned opLength = 0;
    if (jj > 0)
    {
        cigar.addOperation(jj, Cigar::DELETE);
    }
    while(ii >= 0 && jj >= 0 && jj <= 15)
    {
        ++opLength;
#if __AVX2__
        const unsigned nextMaxType = T_[(ii * 3 + maxType) * sizeof(__m256i) + (jj * 2)];
#else
        const unsigned nextMaxType = T_[(ii * 3 + maxType) * sizeof(__m128i) + jj];
#endif
        if (nextMaxType != maxType)
        {
            cigar.addOperation(opLength, opCodes[maxType]);
            opLength = 0;
        }
        ii += iiIncrement[maxType];
        jj += jjIncrement[maxType];
        maxType = nextMaxType;
    }
    assert(-1 == ii);
    if (1 != maxType && opLength)
    {
        cigar.addOperation(opLength, opCodes[maxType]);
        opLength = 0;
    }
    if (15 > jj)
    {
        cigar.addOperation(opLength + 15 - jj, Cigar::DELETE);
        opLength = 0;
    }
    assert(0 == opLength);
    unsigned ret = trimTailIndels(cigar, originalCigarSize);
    std::reverse(cigar.begin() + originalCigarSize, cigar.end());
    trimTailIndels(cigar, originalCigarSize);
    removeAdjacentIndels(cigar, originalCigarSize);
    return ret;
}

std::string epi8(__m128i v)
{
    std::string result;
    for (unsigned int i = 0; 16 > i; ++i)
    {
        result += (boost::format("%3d") % (_mm_extract_epi16(v, 7) >> 8)).str();
        v = _mm_slli_si128(v, 1);
    }
    return result;
}

std::string epi8s(__m128i v)
{
    std::string result;
    for (unsigned int i = 0; 16 > i; ++i)
    {
        result += (boost::format("%3d") % ((int)(_mm_extract_epi16(v, 7)) >> 8)).str();
        v = _mm_slli_si128(v, 1);
    }
    return result;
}

std::string epi8c(__m128i v)
{
    std::string result;
    for (unsigned int i = 0; 16 > i; ++i)
    {
        result += (boost::format("  %c") % char(_mm_extract_epi16(v, 7) >> 8)).str();
        v = _mm_slli_si128(v, 1);
    }
    return result;
}

std::ostream &operator<<(std::ostream &os, const __m128i &H)
{   
    __m128i tmp0 = H;
    for (unsigned int i = 0; i < 8; ++i)
    {
        short v = _mm_extract_epi16(tmp0, 7);
        std::cerr << (boost::format("%3d") % std::max(short(-99), v)).str();
        tmp0 = _mm_slli_si128(tmp0, 2);
    }
    return os;
}

} // namespace alignment
} // namespace isaac
