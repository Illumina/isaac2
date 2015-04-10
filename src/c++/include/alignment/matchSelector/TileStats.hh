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
 ** \file MatchSelectorTileStats.hh
 **
 ** \brief Statistics aggregation helper.
 **
 ** \author Roman Petrovski
 **/

#ifndef ISAAC_ALIGNMENT_MATCH_SELECTOR_TILE_STATS_H
#define ISAAC_ALIGNMENT_MATCH_SELECTOR_TILE_STATS_H

#include "alignment/TemplateLengthStatistics.hh"

namespace isaac
{
namespace alignment
{
namespace matchSelector
{

struct TileStats
{
    explicit TileStats(): uniquelyAlignedFragmentCount_(0)
    {
        std::fill(alignmentScoreFragments_, alignmentScoreFragments_ + maxAlignmentScore_ + 1, 0);
        std::fill(alignmentScoreMismatches_, alignmentScoreMismatches_ + maxAlignmentScore_ + 1, 0);

        std::fill(alignmentScoreTemplates_, alignmentScoreTemplates_ + maxAlignmentScore_ + 1, 0);
        std::fill(alignmentScoreTemplateMismatches_, alignmentScoreTemplateMismatches_ + maxAlignmentScore_ + 1, 0);

        std::fill(cycleBlanks_, cycleBlanks_ + MAX_CYCLES, 0);
        std::fill(cycleMismatches_, cycleMismatches_ + MAX_CYCLES, 0);
        std::fill(cycle1MismatchFragments_, cycle1MismatchFragments_ + MAX_CYCLES, 0);
        std::fill(cycle2MismatchFragments_, cycle2MismatchFragments_ + MAX_CYCLES, 0);
        std::fill(cycle3MismatchFragments_, cycle3MismatchFragments_ + MAX_CYCLES, 0);
        std::fill(cycle4MismatchFragments_, cycle4MismatchFragments_ + MAX_CYCLES, 0);
        std::fill(cycleMoreMismatchFragments_, cycleMoreMismatchFragments_ + MAX_CYCLES, 0);

        std::fill(cycleUniquelyAlignedBlanks_, cycleUniquelyAlignedBlanks_ + MAX_CYCLES, 0);
        std::fill(cycleUniquelyAlignedMismatches_, cycleUniquelyAlignedMismatches_ + MAX_CYCLES, 0);
        std::fill(cycleUniquelyAligned1MismatchFragments_, cycleUniquelyAligned1MismatchFragments_ + MAX_CYCLES, 0);
        std::fill(cycleUniquelyAligned2MismatchFragments_, cycleUniquelyAligned2MismatchFragments_ + MAX_CYCLES, 0);
        std::fill(cycleUniquelyAligned3MismatchFragments_, cycleUniquelyAligned3MismatchFragments_ + MAX_CYCLES, 0);
        std::fill(cycleUniquelyAligned4MismatchFragments_, cycleUniquelyAligned4MismatchFragments_ + MAX_CYCLES, 0);
        std::fill(cycleUniquelyAlignedMoreMismatchFragments_, cycleUniquelyAlignedMoreMismatchFragments_ + MAX_CYCLES, 0);
    }

    void reset()
    {
        *this = TileStats();
    }

    void reserve(const unsigned reserveClusters)
    {
    }

    static const unsigned maxAlignmentScore_ = 0x1FFF;
    unsigned long alignmentScoreFragments_[maxAlignmentScore_ + 1];
    unsigned long alignmentScoreMismatches_[maxAlignmentScore_ + 1];

    unsigned long alignmentScoreTemplates_[maxAlignmentScore_ + 1];
    unsigned long alignmentScoreTemplateMismatches_[maxAlignmentScore_ + 1];

    static const unsigned MAX_CYCLES = 1024;
    unsigned long cycleBlanks_[MAX_CYCLES];
    unsigned long cycleUniquelyAlignedBlanks_[MAX_CYCLES];
    unsigned long cycleMismatches_[MAX_CYCLES];
    unsigned long cycleUniquelyAlignedMismatches_[MAX_CYCLES];

    long cycleUniquelyAligned1MismatchFragments_[MAX_CYCLES];
    long cycleUniquelyAligned2MismatchFragments_[MAX_CYCLES];
    long cycleUniquelyAligned3MismatchFragments_[MAX_CYCLES];
    long cycleUniquelyAligned4MismatchFragments_[MAX_CYCLES];
    long cycleUniquelyAlignedMoreMismatchFragments_[MAX_CYCLES];

    long cycle1MismatchFragments_[MAX_CYCLES];
    long cycle2MismatchFragments_[MAX_CYCLES];
    long cycle3MismatchFragments_[MAX_CYCLES];
    long cycle4MismatchFragments_[MAX_CYCLES];
    long cycleMoreMismatchFragments_[MAX_CYCLES];

    unsigned long uniquelyAlignedFragmentCount_;

    template<typename TemplateT>
    void recordTemplate(const TemplateT &templ)
    {
        if (templ.hasAlignmentScore())
        {
            ISAAC_ASSERT_MSG(maxAlignmentScore_ >= templ.getAlignmentScore(), "alignment score " << templ.getAlignmentScore() <<
                "is too big, if this is expected, change and recompile. " << templ);
            ++alignmentScoreTemplates_[templ.getAlignmentScore()];
            alignmentScoreTemplateMismatches_[templ.getAlignmentScore()] += templ.getMismatches();
        }
    }
    template<typename FragmentT>
    void recordFragment(const FragmentT &fragment, const flowcell::ReadMetadata &readMetadata)
    {

        if (fragment.hasAlignmentScore())
        {
            ISAAC_ASSERT_MSG(maxAlignmentScore_ >= fragment.getAlignmentScore(), "alignment score " << fragment.getAlignmentScore() <<
                "is too big, if this is expected, change and recompile");
            ++alignmentScoreFragments_[fragment.getAlignmentScore()];
            alignmentScoreMismatches_[fragment.getAlignmentScore()] += fragment.getMismatches();
        }

        ISAAC_ASSERT_MSG(readMetadata.getFirstCycle() + readMetadata.getLength() < MAX_CYCLES,
                         "Cycle number is too great, check the MAX_CYCLES constant.");

        std::transform(fragment.getForwardSequenceBegin(), fragment.getForwardSequenceEnd(),
                       cycleBlanks_ + readMetadata.getFirstCycle(),
                       cycleBlanks_ + readMetadata.getFirstCycle(),
                       boost::bind(std::plus<unsigned long>(), _2, boost::bind(&boost::cref<char>, _1) == oligo::SEQUENCE_OLIGO_N));
        std::for_each(fragment.mismatchCyclesBegin(), fragment.mismatchCyclesEnd(),
                      boost::bind(&TileStats::incrementCycleMismatches, this, _1));

        std::size_t mismatches = 0;
        BOOST_FOREACH(const unsigned short &cycle, std::make_pair(fragment.mismatchCyclesBegin(), fragment.mismatchCyclesEnd()))
        {
            incrementCycleXMismatchFragments(cycle, ++mismatches);
        }

        if (fragment.isUniquelyAligned())
        {
            std::transform(fragment.getForwardSequenceBegin(), fragment.getForwardSequenceEnd(),
                           cycleUniquelyAlignedBlanks_ + readMetadata.getFirstCycle(),
                           cycleUniquelyAlignedBlanks_ + readMetadata.getFirstCycle(),
                           boost::bind(std::plus<unsigned long>(), _2, boost::bind(&boost::cref<char>, _1) == oligo::SEQUENCE_OLIGO_N));
            std::for_each(fragment.mismatchCyclesBegin(), fragment.mismatchCyclesEnd(),
                          boost::bind(&TileStats::incrementCycleUniquelyAlignedMismatches, this, _1));

            mismatches = 0;
            BOOST_FOREACH(const unsigned short &cycle, std::make_pair(fragment.mismatchCyclesBegin(), fragment.mismatchCyclesEnd()))
            {
                incrementCycleUniquelyAlignedXMismatchFragments(cycle, ++mismatches);
            }
        }

        if (fragment.isUniquelyAligned())
        {
            ++uniquelyAlignedFragmentCount_;
        }
    }

    const TileStats &operator +=(const TileStats &right)
    {
        std::transform(alignmentScoreFragments_, alignmentScoreFragments_ + maxAlignmentScore_ + 1,
                       right.alignmentScoreFragments_, alignmentScoreFragments_,
                       std::plus<unsigned long>());

        std::transform(alignmentScoreMismatches_, alignmentScoreMismatches_ + maxAlignmentScore_ + 1,
                       right.alignmentScoreMismatches_, alignmentScoreMismatches_,
                       std::plus<unsigned long>());

        std::transform(alignmentScoreTemplates_, alignmentScoreTemplates_ + maxAlignmentScore_ + 1,
                       right.alignmentScoreTemplates_, alignmentScoreTemplates_,
                       std::plus<unsigned long>());

        std::transform(alignmentScoreTemplateMismatches_, alignmentScoreTemplateMismatches_ + maxAlignmentScore_ + 1,
                       right.alignmentScoreTemplateMismatches_, alignmentScoreTemplateMismatches_,
                       std::plus<unsigned long>());


        std::transform(cycleBlanks_, cycleBlanks_ + MAX_CYCLES,
                       right.cycleBlanks_, cycleBlanks_,
                       std::plus<unsigned long>());

        std::transform(cycleMismatches_, cycleMismatches_ + MAX_CYCLES,
                       right.cycleMismatches_, cycleMismatches_,
                       std::plus<unsigned long>());

        std::transform(cycle1MismatchFragments_, cycle1MismatchFragments_ + MAX_CYCLES,
                       right.cycle1MismatchFragments_, cycle1MismatchFragments_,
                       std::plus<unsigned long>());

        std::transform(cycle2MismatchFragments_, cycle2MismatchFragments_ + MAX_CYCLES,
                       right.cycle2MismatchFragments_, cycle2MismatchFragments_,
                       std::plus<unsigned long>());

        std::transform(cycle3MismatchFragments_, cycle3MismatchFragments_ + MAX_CYCLES,
                       right.cycle3MismatchFragments_, cycle3MismatchFragments_,
                       std::plus<unsigned long>());

        std::transform(cycle4MismatchFragments_, cycle4MismatchFragments_ + MAX_CYCLES,
                       right.cycle4MismatchFragments_, cycle4MismatchFragments_,
                       std::plus<unsigned long>());

        std::transform(cycleMoreMismatchFragments_, cycleMoreMismatchFragments_ + MAX_CYCLES,
                       right.cycleMoreMismatchFragments_, cycleMoreMismatchFragments_,
                       std::plus<unsigned long>());


        std::transform(cycleUniquelyAlignedBlanks_, cycleUniquelyAlignedBlanks_ + MAX_CYCLES,
                       right.cycleUniquelyAlignedBlanks_, cycleUniquelyAlignedBlanks_,
                       std::plus<unsigned long>());

        std::transform(cycleUniquelyAlignedMismatches_, cycleUniquelyAlignedMismatches_ + MAX_CYCLES,
                       right.cycleUniquelyAlignedMismatches_, cycleUniquelyAlignedMismatches_,
                       std::plus<unsigned long>());

        std::transform(cycleUniquelyAligned1MismatchFragments_, cycleUniquelyAligned1MismatchFragments_ + MAX_CYCLES,
                       right.cycleUniquelyAligned1MismatchFragments_, cycleUniquelyAligned1MismatchFragments_,
                       std::plus<unsigned long>());

        std::transform(cycleUniquelyAligned2MismatchFragments_, cycleUniquelyAligned2MismatchFragments_ + MAX_CYCLES,
                       right.cycleUniquelyAligned2MismatchFragments_, cycleUniquelyAligned2MismatchFragments_,
                       std::plus<unsigned long>());

        std::transform(cycleUniquelyAligned3MismatchFragments_, cycleUniquelyAligned3MismatchFragments_ + MAX_CYCLES,
                       right.cycleUniquelyAligned3MismatchFragments_, cycleUniquelyAligned3MismatchFragments_,
                       std::plus<unsigned long>());

        std::transform(cycleUniquelyAligned4MismatchFragments_, cycleUniquelyAligned4MismatchFragments_ + MAX_CYCLES,
                       right.cycleUniquelyAligned4MismatchFragments_, cycleUniquelyAligned4MismatchFragments_,
                       std::plus<unsigned long>());

        std::transform(cycleUniquelyAlignedMoreMismatchFragments_, cycleUniquelyAlignedMoreMismatchFragments_ + MAX_CYCLES,
                       right.cycleUniquelyAlignedMoreMismatchFragments_, cycleUniquelyAlignedMoreMismatchFragments_,
                       std::plus<unsigned long>());

        uniquelyAlignedFragmentCount_ += right.uniquelyAlignedFragmentCount_;

        return *this;
    }

    const TileStats operator +(const TileStats &right) const
    {
        TileStats ret(*this);
        ret += right;
        return ret;
    }

    void finalize()
    {
        // if at cycle 1 the fragment had 1 mismatch, this mismatch propagates to all subsequent cycles.
        std::transform(cycle1MismatchFragments_ + 1, cycle1MismatchFragments_ + MAX_CYCLES,
                       cycle1MismatchFragments_, cycle1MismatchFragments_ + 1,
                       std::plus<long>());

        std::transform(cycle2MismatchFragments_ + 1, cycle2MismatchFragments_ + MAX_CYCLES,
                       cycle2MismatchFragments_, cycle2MismatchFragments_ + 1,
                       std::plus<long>());

        std::transform(cycle3MismatchFragments_ + 1, cycle3MismatchFragments_ + MAX_CYCLES,
                       cycle3MismatchFragments_, cycle3MismatchFragments_ + 1,
                       std::plus<long>());

        std::transform(cycle4MismatchFragments_ + 1, cycle4MismatchFragments_ + MAX_CYCLES,
                       cycle4MismatchFragments_, cycle4MismatchFragments_ + 1,
                       std::plus<long>());

        std::transform(cycleMoreMismatchFragments_ + 1, cycleMoreMismatchFragments_ + MAX_CYCLES,
                       cycleMoreMismatchFragments_, cycleMoreMismatchFragments_ + 1,
                       std::plus<long>());

        // once a 1-mismatch fragments gets a second mismatch, it stops beinga 1-mismatch fragment. Substract.
        std::transform(cycle1MismatchFragments_, cycle1MismatchFragments_ + MAX_CYCLES,
                       cycle2MismatchFragments_, cycle1MismatchFragments_,
                       std::minus<long>());
        std::transform(cycle2MismatchFragments_, cycle2MismatchFragments_ + MAX_CYCLES,
                       cycle3MismatchFragments_, cycle2MismatchFragments_,
                       std::minus<long>());
        std::transform(cycle3MismatchFragments_, cycle3MismatchFragments_ + MAX_CYCLES,
                       cycle4MismatchFragments_, cycle3MismatchFragments_,
                       std::minus<long>());
        std::transform(cycle4MismatchFragments_, cycle4MismatchFragments_ + MAX_CYCLES,
                       cycleMoreMismatchFragments_, cycle4MismatchFragments_,
                       std::minus<long>());

        // two mismatch fragments include one mismatch fragments and so on, add sequentially
        std::transform(cycle2MismatchFragments_, cycle2MismatchFragments_ + MAX_CYCLES,
                       cycle1MismatchFragments_, cycle2MismatchFragments_,
                       std::plus<long>());
        std::transform(cycle3MismatchFragments_, cycle3MismatchFragments_ + MAX_CYCLES,
                       cycle2MismatchFragments_, cycle3MismatchFragments_,
                       std::plus<long>());
        std::transform(cycle4MismatchFragments_, cycle4MismatchFragments_ + MAX_CYCLES,
                       cycle3MismatchFragments_, cycle4MismatchFragments_,
                       std::plus<long>());
        std::transform(cycleMoreMismatchFragments_, cycleMoreMismatchFragments_ + MAX_CYCLES,
                       cycle4MismatchFragments_, cycleMoreMismatchFragments_,
                       std::plus<long>());

        // if at cycle 1 the fragment had 1 mismatch, this mismatch propagates to all subsequent cycles.
        std::transform(cycleUniquelyAligned1MismatchFragments_ + 1, cycleUniquelyAligned1MismatchFragments_ + MAX_CYCLES,
                       cycleUniquelyAligned1MismatchFragments_, cycleUniquelyAligned1MismatchFragments_ + 1,
                       std::plus<long>());

        std::transform(cycleUniquelyAligned2MismatchFragments_ + 1, cycleUniquelyAligned2MismatchFragments_ + MAX_CYCLES,
                       cycleUniquelyAligned2MismatchFragments_, cycleUniquelyAligned2MismatchFragments_ + 1,
                       std::plus<long>());

        std::transform(cycleUniquelyAligned3MismatchFragments_ + 1, cycleUniquelyAligned3MismatchFragments_ + MAX_CYCLES,
                       cycleUniquelyAligned3MismatchFragments_, cycleUniquelyAligned3MismatchFragments_ + 1,
                       std::plus<long>());

        std::transform(cycleUniquelyAligned4MismatchFragments_ + 1, cycleUniquelyAligned4MismatchFragments_ + MAX_CYCLES,
                       cycleUniquelyAligned4MismatchFragments_, cycleUniquelyAligned4MismatchFragments_ + 1,
                       std::plus<long>());

        std::transform(cycleUniquelyAlignedMoreMismatchFragments_ + 1, cycleUniquelyAlignedMoreMismatchFragments_ + MAX_CYCLES,
                       cycleUniquelyAlignedMoreMismatchFragments_, cycleUniquelyAlignedMoreMismatchFragments_ + 1,
                       std::plus<long>());

        // once a 1-mismatch fragments gets a second mismatch, it stops beinga 1-mismatch fragment. Substract.
        std::transform(cycleUniquelyAligned1MismatchFragments_, cycleUniquelyAligned1MismatchFragments_ + MAX_CYCLES,
                       cycleUniquelyAligned2MismatchFragments_, cycleUniquelyAligned1MismatchFragments_,
                       std::minus<long>());
        std::transform(cycleUniquelyAligned2MismatchFragments_, cycleUniquelyAligned2MismatchFragments_ + MAX_CYCLES,
                       cycleUniquelyAligned3MismatchFragments_, cycleUniquelyAligned2MismatchFragments_,
                       std::minus<long>());
        std::transform(cycleUniquelyAligned3MismatchFragments_, cycleUniquelyAligned3MismatchFragments_ + MAX_CYCLES,
                       cycleUniquelyAligned4MismatchFragments_, cycleUniquelyAligned3MismatchFragments_,
                       std::minus<long>());
        std::transform(cycleUniquelyAligned4MismatchFragments_, cycleUniquelyAligned4MismatchFragments_ + MAX_CYCLES,
                       cycleUniquelyAlignedMoreMismatchFragments_, cycleUniquelyAligned4MismatchFragments_,
                       std::minus<long>());

        // two mismatch fragments include one mismatch fragments and so on, add sequentially
        std::transform(cycleUniquelyAligned2MismatchFragments_, cycleUniquelyAligned2MismatchFragments_ + MAX_CYCLES,
                       cycleUniquelyAligned1MismatchFragments_, cycleUniquelyAligned2MismatchFragments_,
                       std::plus<long>());
        std::transform(cycleUniquelyAligned3MismatchFragments_, cycleUniquelyAligned3MismatchFragments_ + MAX_CYCLES,
                       cycleUniquelyAligned2MismatchFragments_, cycleUniquelyAligned3MismatchFragments_,
                       std::plus<long>());
        std::transform(cycleUniquelyAligned4MismatchFragments_, cycleUniquelyAligned4MismatchFragments_ + MAX_CYCLES,
                       cycleUniquelyAligned3MismatchFragments_, cycleUniquelyAligned4MismatchFragments_,
                       std::plus<long>());
        std::transform(cycleUniquelyAlignedMoreMismatchFragments_, cycleUniquelyAlignedMoreMismatchFragments_ + MAX_CYCLES,
                       cycleUniquelyAligned4MismatchFragments_, cycleUniquelyAlignedMoreMismatchFragments_,
                       std::plus<long>());
    }

private:
    void incrementCycleMismatches(const unsigned cycle)
    {
        ISAAC_ASSERT_MSG(cycle < MAX_CYCLES, "Cycle number too high.");
        ++cycleMismatches_[cycle];
    }

    void incrementCycleUniquelyAlignedMismatches(const unsigned cycle)
    {
        ISAAC_ASSERT_MSG(cycle < MAX_CYCLES, "Cycle number too high.");
        ++cycleUniquelyAlignedMismatches_[cycle];
    }

    void incrementCycleUniquelyAlignedXMismatchFragments(const unsigned short cycle, const unsigned mismatches)
    {
        ISAAC_ASSERT_MSG(cycle < MAX_CYCLES, "Cycle number too high.");
        switch (mismatches)
        {
        case 0:
            ISAAC_ASSERT_MSG(0 != mismatches, "incrementCycleUniquelyAlignedXMismatchFragments should not be called for 0 mismatches");
            break;
        case 1:
            ++cycleUniquelyAligned1MismatchFragments_[cycle];
            break;
        case 2:
            ++cycleUniquelyAligned2MismatchFragments_[cycle];
            break;
        case 3:
            ++cycleUniquelyAligned3MismatchFragments_[cycle];
            break;
        case 4:
            ++cycleUniquelyAligned4MismatchFragments_[cycle];
            break;
        case 5:
            ++cycleUniquelyAlignedMoreMismatchFragments_[cycle];
            break;
        default:
            // don't count fragments with more than 5 mismatches;
            break;
        }
    }

    void incrementCycleXMismatchFragments(const unsigned short cycle, const unsigned mismatches)
    {
        ISAAC_ASSERT_MSG(cycle < MAX_CYCLES, "Cycle number too high.");
        switch (mismatches)
        {
        case 0:
            ISAAC_ASSERT_MSG(0 != mismatches, "incrementCycleXMismatchFragments should not be called for 0 mismatches");
            break;
        case 1:
            ++cycle1MismatchFragments_[cycle];
            break;
        case 2:
            ++cycle2MismatchFragments_[cycle];
            break;
        case 3:
            ++cycle3MismatchFragments_[cycle];
            break;
        case 4:
            ++cycle4MismatchFragments_[cycle];
            break;
        case 5:
            ++cycleMoreMismatchFragments_[cycle];
            break;
        default:
            // don't count fragments with more than 5 mismatches;
            break;
        }
    }
};

} //namespace matchSelector
} //namespace alignment
} //namespace isaac

#endif //ISAAC_ALIGNMENT_MATCH_SELECTOR_TILE_STATS_H
