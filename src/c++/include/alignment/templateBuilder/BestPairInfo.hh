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
 ** \file BestPairInfo.hh
 **
 ** \brief Internal structures used by TemplateBuilder
 ** 
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_ALIGNMENT_TEMPLATE_BUILDER_BEST_PAIR_INFO_HH
#define iSAAC_ALIGNMENT_TEMPLATE_BUILDER_BEST_PAIR_INFO_HH

#include <boost/function_output_iterator.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>

#include "../../common/StaticVector.hh"
#include "common/Debug.hh"

namespace isaac
{
namespace alignment
{
namespace templateBuilder
{

static const unsigned TRACKED_REPEATS_MAX_ONE_READ = 1000;

/// arrays temporary used in buildDisjointTemplate and rescueShadow
class ShadowProbability
{
    reference::ReferencePosition pos_;
    double logProbability_;
    long observedLength_;
public:
    ShadowProbability(const FragmentMetadata &shadow) :
        pos_(shadow.getFStrandReferencePosition()), logProbability_(shadow.logProbability),
        observedLength_(shadow.getObservedLength())
    {
        // encode reverse in position to save space.
        // This structure can easily take 5 extra gigabytes if a boolean is introduced to store reverse (SAAC-478)
        pos_.setNeighbors(shadow.isReverse());
    }

    reference::ReferencePosition pos() const {return pos_;}
    double logProbability() const {return logProbability_;}
    double probability() const {return exp(logProbability_);}
    long observedLength() const {return observedLength_;}

    bool operator < (const ShadowProbability &that) const
    {
        //if we happen to have fragment and its inversion, let's have the higher probability on top so that it counts towards the total
        return pos_ < that.pos_ ||
            (pos_ == that.pos_ &&
                (logProbability_ < that.logProbability_ ||
                    (logProbability_ == that.logProbability_ &&
                        (observedLength_ < that.observedLength_)
                    )
                )
            );
    }

    bool operator == (const ShadowProbability &that) const
    {
        return pos_ == that.pos_ && logProbability_ == that.logProbability_ &&
            observedLength_ == that.observedLength_;
    }

    bool operator != (const ShadowProbability &that) const { return !(*this == that); }

    friend std::ostream & operator << (std::ostream &os, const ShadowProbability &sp)
    {
        return os << "sp(" << sp.pos_ << ":" << sp.logProbability_<< ":" << sp.observedLength_ << ")";
    }

};

struct PairProbability
{
    ShadowProbability r1_;
    ShadowProbability r2_;

    PairProbability(const FragmentMetadata &r1,
                    const FragmentMetadata &r2) : r1_(r1), r2_(r2){}

    bool operator < (const PairProbability &that) const
    {
        // since the sum of log probabilities of reads has to be considered we can't compare reads individually
        return
            (r1_.pos() < that.r1_.pos() || (r1_.pos() == that.r1_.pos() &&
                (r2_.pos() < that.r2_.pos() || (r2_.pos() == that.r2_.pos() &&
                    //if we happen to have pair and its inversion, let's have the higher probability on top so that it counts towards the total
                    (that.logProbability() < logProbability() || (logProbability() == that.logProbability() &&
                        (r1_.observedLength() < that.r1_.observedLength() || (r1_.observedLength() == that.r1_.observedLength() &&
                            r2_.observedLength() < that.r2_.observedLength()))))))));
    }

    bool operator == (const PairProbability &that) const
    {
        // since the sum of log probabilities of reads has to be considered we can't compare reads individually
        return r1_.pos() == that.r1_.pos() && r2_.pos() == that.r2_.pos() && logProbability() == that.logProbability() &&
            r1_.observedLength() == that.r1_.observedLength() && r2_.observedLength() == that.r2_.observedLength();
    }

    bool operator != (const PairProbability &that) const {return !(*this == that);}

    double logProbability() const {return r1_.logProbability() + r2_.logProbability();}
    double probability() const {return exp(logProbability());}

    friend std::ostream & operator << (std::ostream &os, const PairProbability &pp)
    {
        return os << "PairProbability(" << pp.r1_ << "," << pp.r2_<< ")";
    }
};

struct PairInfo
{
    PairInfo() : logProbability_(-std::numeric_limits<double>::max()), swScore_(-1), matchModel_(false) {}
    PairInfo(const FragmentMetadata &oneRead, const FragmentMetadata &anotherRead, bool matchModel):
        logProbability_(oneRead.logProbability + anotherRead.logProbability),
        swScore_(oneRead.smithWatermanScore + anotherRead.smithWatermanScore), matchModel_(matchModel) {}

    void clear()
    {
        *this = PairInfo();
    }

    double probability() const {return exp(logProbability_);}

    double logProbability_;
    unsigned long swScore_;
    bool matchModel_;

    friend std::ostream & operator << (std::ostream & os, const PairInfo& info)
    {
        return os << "PairInfo(" << info.logProbability_ << "lp " << info.swScore_ << "sws)";
    }
};

struct BestPairInfo
{
    typedef std::vector<FragmentMetadata>::const_iterator FragmentIterator;

    BestPairInfo(
        const unsigned repeatThreshold,
        const unsigned maxSeedsPerRead)
    {
        /// Max number of orphans times the max number of shadows that can be rescued for each orphan plus the max number of shadows that have been discovered during seed matching
        for (unsigned r = 0; r < READS_IN_A_PAIR; ++r)
        {
            // In the worst case each seed will generate --repeat-threshold alignments which then will generate TRACKED_REPEATS_MAX_ONE_READ shadows
            // So, max seeds per read times the max number of orphans times the max number of shadows that can be rescued for each orphan
            // plus the max number of shadows that have been discovered during seed matching...
            readProbabilities_[r].reserve(
                maxSeedsPerRead * repeatThreshold * TRACKED_REPEATS_MAX_ONE_READ + TRACKED_REPEATS_MAX_ONE_READ);
        }

        // In the worst case each seed will generate --repeat-threshold alignments which then will generate TRACKED_REPEATS_MAX_ONE_READ shadows
        // So, max seeds per read times max number of orphans times the max number of shadows that can be rescued for each orphan times the number of reads
        pairProbabilities_.reserve(maxSeedsPerRead * repeatThreshold * TRACKED_REPEATS_MAX_ONE_READ * READS_IN_A_PAIR);
        repeats_.reserve(TRACKED_REPEATS_MAX_ONE_READ * READS_IN_A_PAIR);
    }

    void clear()
    {
        info_.clear();

        repeats_.clear();

        readProbabilities_[0].clear();
        readProbabilities_[1].clear();

        pairProbabilities_.clear();
    }

    bool empty() const
    {
        return repeats_.empty();
    }

    void appendBest(const FragmentMetadata &oneRead, const FragmentMetadata &anotherRead)
    {
        if (oneRead.getReadIndex())
        {
            // make sure first is always r1
            appendBest(anotherRead, oneRead);
            return;
        }
        if (repeats_.capacity() != repeats_.size())
        {
            repeats_.push_back(BamTemplate(oneRead, anotherRead));
        }
        appendProbabilities(oneRead, anotherRead);
    }

    void appendProbabilities(const FragmentMetadata &oneRead, const FragmentMetadata &anotherRead)
    {
        if (oneRead.getReadIndex())
        {
            // make sure first is always r1
            appendProbabilities(anotherRead, oneRead);
            return;
        }
        ISAAC_ASSERT_MSG(oneRead.getReadIndex() != anotherRead.getReadIndex(), "Read indices match in the pair" << oneRead << " " << anotherRead);

        pairProbabilities_.push_back(templateBuilder::PairProbability(oneRead, anotherRead));
        readProbabilities_[0].push_back(oneRead);
        readProbabilities_[1].push_back(anotherRead);
    }

    void resetBest(const PairInfo &info, const FragmentMetadata &oneRead, const FragmentMetadata &anotherRead)
    {
        repeats_.clear();

        info_ = info;
        appendBest(oneRead, anotherRead);
    }

    inline bool isWorseThan(const PairInfo &that) const
    {
        return empty() ||
            ISAAC_LP_LESS(info_.logProbability_, that.logProbability_) ||
            (that.swScore_ < info_.swScore_ && ISAAC_LP_EQUALS(info_.logProbability_, that.logProbability_));
    }

    inline bool isAsGood(const PairInfo &that) const
    {
        return !empty() && that.swScore_ == info_.swScore_ && ISAAC_LP_EQUALS(that.logProbability_, info_.logProbability_);
    }

    double sumUniqueReadProbabilities(const std::size_t readIndex, const double ignoreLp)
    {
        bool ignoreIgnored = false;
        double ret = 0.0;
        std::sort(readProbabilities_[readIndex].begin(), readProbabilities_[readIndex].end());

        readProbabilities_[readIndex].erase(std::unique(readProbabilities_[readIndex].begin(), readProbabilities_[readIndex].end()), readProbabilities_[readIndex].end());
        BOOST_FOREACH(const ShadowProbability &sp, readProbabilities_[readIndex])
        {
            // individual probabilities are copies. simple comparison should work
            if (!ignoreIgnored && sp.logProbability() == ignoreLp)
            {
                ignoreIgnored = true;
            }
            else
            {
                ret += sp.probability();
            }
        }

        return ret;
    }

    double sumUniquePairProbabilities(const double ignoreLp)
    {
        bool ignoreIgnored = false;
        double ret = 0.0;
        std::sort(pairProbabilities_.begin(), pairProbabilities_.end());
        pairProbabilities_.erase(std::unique(pairProbabilities_.begin(), pairProbabilities_.end()), pairProbabilities_.end());

        BOOST_FOREACH(const PairProbability &pp, pairProbabilities_)
        {
            // individual probabilities are copies. simple comparison should work
            if (!ignoreIgnored && pp.logProbability() == ignoreLp)
            {
                ignoreIgnored = true;
                ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(repeats_.front().getFragmentMetadata(0).getCluster().getId(),
                                                       "sumUniquePairProbabilities: " << pp << " Ignored");
            }
            else
            {
                ret += pp.probability();
                ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(repeats_.front().getFragmentMetadata(0).getCluster().getId(),
                                                       "sumUniquePairProbabilities: " << pp << " total: " << ret);
            }
        }

        return ret;
    }

    /**
     * \return number of unique pair repeat alignments
     */
    std::size_t removeRepeatDuplicates()
    {
        ISAAC_ASSERT_MSG(!empty(), "Invalid method call on an empty BestPairInfo" << *this);
        if (repeats_.size() == 1)
        {
            return 1;
        }

        std::sort(repeats_.begin(), repeats_.end());

        repeats_.erase(std::unique(repeats_.begin(), repeats_.end()), repeats_.end());
        return repeats_.size();
    }

    unsigned bestPairMismatchCount() const
    {
        return repeats_.front().getFragmentMetadata(0).getMismatchCount() + repeats_.front().getFragmentMetadata(1).getMismatchCount();
    }

    bool isKUnique() const
    {
        return !repeats_.empty() && repeats_.front().isKUnique();
    }
    // potentially this needs to hold all the combinations of individual fragment alignments (TRACKED_REPEATS_MAX_ONE_READ * TRACKED_REPEATS_MAX_ONE_READ).
    // But repeat scattering hardly justifies to keeping this many equivalent pairs.
    typedef isaac::common::StaticVector<FragmentIterator, TRACKED_REPEATS_MAX_ONE_READ> FragmentIteratorVector;
    // corresponding entries represent read pairs. All read pairs in repeatAlignments_ represent alignments of the same bestTemplateScore_ to the repeat
    static const unsigned READS_IN_A_PAIR = 2;
    PairInfo info_;
    std::vector<templateBuilder::ShadowProbability> readProbabilities_[READS_IN_A_PAIR];
    std::vector<templateBuilder::PairProbability> pairProbabilities_;

    std::vector<BamTemplate> repeats_;

    friend std::ostream & operator << (std::ostream & os, const BestPairInfo& bestPairInfo)
    {
        os << "BestPairInfo(";
        if (bestPairInfo.empty())
        {
            return os << "empty)";
        }
        if (!bestPairInfo.repeats_.empty())
        {
            os << bestPairInfo.repeats_.front() << ",";
        }
        return os << bestPairInfo.info_ << " " << bestPairInfo.bestPairMismatchCount() << "bem)";
    }
};

} // namespace templateBuilder
} // namespace alignment
} // namespace isaac

#endif // #ifndef iSAAC_ALIGNMENT_TEMPLATE_BUILDER_BEST_PAIR_INFO_HH
