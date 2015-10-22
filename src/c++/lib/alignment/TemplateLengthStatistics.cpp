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
 ** \file TemplateLengthStatistics.cpp
 **
 ** \brief See TemplateLengthStatistics.hh
 ** 
 ** \author Come Raczy
 **/


#include "alignment/Quality.hh"
#include "alignment/TemplateLengthStatistics.hh"
#include "alignment/Cigar.hh"
#include "common/Debug.hh"

namespace isaac
{
namespace alignment
{

const double TemplateLengthDistribution::STANDARD_DEVIATIONS_MAX = 3.0;
const double TemplateLengthDistribution::FRAGMENT_LENGTH_CONFIDENCE_INTERVAL =
    boost::math::erf(STANDARD_DEVIATIONS_MAX/boost::math::constants::root_two<double>());
const double TemplateLengthDistribution::FRAGMENT_LENGTH_CONFIDENCE_INTERVAL_1Z = boost::math::erf(1.0/boost::math::constants::root_two<double>());
const double TemplateLengthDistribution::LOWER_PERCENT = (1.0 - FRAGMENT_LENGTH_CONFIDENCE_INTERVAL) / 2.0;
const double TemplateLengthDistribution::UPPER_PERCENT = (1.0 + FRAGMENT_LENGTH_CONFIDENCE_INTERVAL) / 2.0;
const double TemplateLengthDistribution::LOWER_PERCENT_1Z = (1.0 - FRAGMENT_LENGTH_CONFIDENCE_INTERVAL_1Z) / 2.0;
const double TemplateLengthDistribution::UPPER_PERCENT_1Z = (1.0 + FRAGMENT_LENGTH_CONFIDENCE_INTERVAL_1Z) / 2.0;

TemplateLengthStatistics::TemplateLengthStatistics()
    : min_(-1U)
    , max_(-1U)
    , median_(-1U)
    , lowStdDev_(-1U)
    , highStdDev_(-1U)
    , stable_(false)
    , mateMin_(-1U)
    , mateMax_(-1U)
{
    bestModels_[0] = InvalidAlignmentModel;
    bestModels_[1] = InvalidAlignmentModel;

}

void TemplateLengthStatistics::clear()
{
    min_ = -1U;
    max_ = -1U;
    median_ = -1U;
    lowStdDev_ = -1U;
    highStdDev_ = -1U;
    stable_ = false;
    bestModels_[0] = InvalidAlignmentModel;
    bestModels_[1] = InvalidAlignmentModel;
}

bool TemplateLengthStatistics::matchModel(const FragmentMetadata &f1, const FragmentMetadata &f2) const
{
    if (f1.getContigId() == f2.getContigId())
    {
        const unsigned long length = getLength(f1, f2);
        const AlignmentModel model = alignmentModel(f1, f2);
        const bool ret = (length <= max_ + TEMPLATE_LENGTH_THRESHOLD) && ((model == bestModels_[0]) || (model == bestModels_[1]));
        ISAAC_THREAD_CERR_DEV_TRACE((boost::format("TemplateLengthStatistics::matchModel %d: %s-%s (len: %u)") % ret %
                                     f1 % f2 % length).str());
        return ret;
    }
    return false;

}

TemplateLengthStatistics::AlignmentClass TemplateLengthStatistics::alignmentClass(const AlignmentModel model)
{
    assert(model < 8);
    return static_cast<AlignmentClass>((model < 4) ? model : ((~model) & 3));
}

const static std::vector<std::string> classNames = boost::assign::list_of("F+")("R+")("R-")("F-")("unknown");
const std::string &TemplateLengthStatistics::alignmentClassName(const AlignmentClass alignmentClass)
{
    return (classNames.size() > size_t(alignmentClass)) ? classNames.at(alignmentClass) : classNames.back();
}

const static std::vector<std::string> modelNames = boost::assign::list_of
        ("FF+")("FR+")("RF+")("RR+")("FF-")("FR-")("RF-")("RR-")("unknown");
const std::string &TemplateLengthStatistics::alignmentModelName(const AlignmentModel alignmentModel)
{
    return (modelNames.size() > size_t(alignmentModel)) ? modelNames.at(alignmentModel) : modelNames.back();
}

void TemplateLengthDistribution::reserve(const unsigned reserveClusters)
{
    ISAAC_ASSERT_MSG(reserveClusters, "reserveClusters must not be 0");
    std::for_each(histograms_.begin(), histograms_.end(), boost::bind(&std::vector<unsigned>::reserve, _1, reserveClusters));
    lengthList_.reserve(reserveClusters);
}

void TemplateLengthDistribution::updateStatistics()
{
    // save the old statistics
    const TemplateLengthStatistics oldStats = stats_;
    // find the two best alignment models
    stats_.setBestModel(histograms_[1].size() <= histograms_[0].size() ? TemplateLengthStatistics::FFp : TemplateLengthStatistics::FRp, 0);
    stats_.setBestModel(static_cast<TemplateLengthStatistics::AlignmentModel>((stats_.getBestModel(0) + 1) % 2), 1);
    for (size_t i = 2; histograms_.size() > i; ++i)
    {
        if (histograms_[i].size() > histograms_[stats_.getBestModel(0)].size())
        {
            stats_.setBestModel(stats_.getBestModel(0), 1);
            stats_.setBestModel(static_cast<TemplateLengthStatistics::AlignmentModel>(i), 0);
        }
        else if (histograms_[i].size() > histograms_[stats_.getBestModel(1)].size())
        {
            stats_.setBestModel(static_cast<TemplateLengthStatistics::AlignmentModel>(i), 1);
        }
    }

    // TODO: report incoherent models
    if (TemplateLengthStatistics::alignmentClass(stats_.getBestModel(0)) != TemplateLengthStatistics::alignmentClass(stats_.getBestModel(1)))
    {
        ISAAC_THREAD_CERR << "Incoherent alignment models:"
            << stats_.getBestModel(0) << " (" << TemplateLengthStatistics::alignmentModelName(stats_.getBestModel(0)) << "->"
                << TemplateLengthStatistics::alignmentClassName(TemplateLengthStatistics::alignmentClass(stats_.getBestModel(0))) << "):"
            << stats_.getBestModel(1) << " (" << TemplateLengthStatistics::alignmentModelName(stats_.getBestModel(1)) << "->"
                << TemplateLengthStatistics::alignmentClassName(TemplateLengthStatistics::alignmentClass(stats_.getBestModel(1))) << ")" << std::endl;
    }
    // sort all the lengths found for the two best models
    lengthList_.clear();
    lengthList_.insert(lengthList_.end(), histograms_[stats_.getBestModel(0)].begin(), histograms_[stats_.getBestModel(0)].end());
    lengthList_.insert(lengthList_.end(), histograms_[stats_.getBestModel(1)].begin(), histograms_[stats_.getBestModel(1)].end());
    std::sort(lengthList_.begin(), lengthList_.end());
    // compute the statistics
    setMin(lengthList_.empty() ? 0 : lengthList_[unsigned(lengthList_.size() * LOWER_PERCENT)]);
    setMedian(lengthList_.empty() ? TemplateLengthStatistics::TEMPLATE_LENGTH_THRESHOLD/2 : lengthList_[unsigned(lengthList_.size() * 0.5)]);
    setMax(lengthList_.empty() ? TemplateLengthStatistics::TEMPLATE_LENGTH_THRESHOLD : lengthList_[unsigned(lengthList_.size() * UPPER_PERCENT)]);
    stats_.setLowStdDev(lengthList_.empty() ?
        stats_.getMedian() : (stats_.getMedian() - lengthList_[unsigned(lengthList_.size() * LOWER_PERCENT_1Z)]));
    stats_.setHighStdDev(lengthList_.empty() ?
        stats_.getMedian() : (lengthList_[unsigned(lengthList_.size() * UPPER_PERCENT_1Z)] - stats_.getMedian()));

//    ISAAC_THREAD_CERR << "updateStatistics oldStats:" << oldStats <<  " stats_:" << stats_ << std::endl;
//    ISAAC_THREAD_CERR << "updateStatistics templateCount_:" << templateCount_ <<  " uniqueCount_:" << uniqueCount_ << std::endl;
    // check if we reached stability
    if (oldStats.getMin() == stats_.getMin() &&
        oldStats.getMedian() == stats_.getMedian() &&
        oldStats.getMax() == stats_.getMax() &&
        oldStats.getLowStdDev() == stats_.getLowStdDev() &&
        oldStats.getHighStdDev() == stats_.getHighStdDev() &&
        oldStats.getBestModel(0) == stats_.getBestModel(0) &&
        oldStats.getBestModel(1) == stats_.getBestModel(1))
    {
        stats_.setStable(true);
    }
}

/// Check that one of the two best models has the given orientation for the readIndex
bool TemplateLengthStatistics::isValidModel(const bool reverse, const unsigned readIndex) const
{
    ISAAC_ASSERT_MSG(READS_MAX > readIndex, "Invalid read index");
    // bit 0 is for read 2; bit 1 is for read 1
    const unsigned shift = (readIndex + 1) % 2;
    return (reverse == ((bestModels_[0] >> shift) & 1)) || (reverse == ((bestModels_[1] >> shift) & 1));
}

bool TemplateLengthStatistics::firstFragment(const bool reverse, const unsigned readIndex) const
{
    ISAAC_ASSERT_MSG(READS_MAX > readIndex, "Invalid read index");
    // bit 0 is for read 2; bit 1 is for read 1
    const unsigned shift = (readIndex + 1) % 2;
    for (unsigned i = 0; 2 > i; ++i)
    {
        if (reverse == ((bestModels_[i] >> shift) & 1))
        {
            return ((bestModels_[i] >> 2) & 1)== readIndex;
        }
    }
    assert(false); // shouldn't get there
    return false;
}

bool TemplateLengthStatistics::mateOrientation(unsigned readIndex, bool reverse) const
{
    // bit 0 is for readIndex 1 and vice-versa
    const unsigned shift = (readIndex + 1) % 2;
    for (unsigned i = 0; 2 > i; ++i)
    {
        if (reverse == ((bestModels_[i] >> shift) & 1))
        {
            return (bestModels_[i] >> readIndex) & 1;
        }
    }
    // for discorant models, return the expected orientation in the best model
    return (bestModels_[0] >> readIndex) & 1;
}

long TemplateLengthStatistics::mateMinPosition(
    const unsigned readIndex,
    const bool reverse,
    const long position,
    const unsigned *readLengths) const
{
    ISAAC_ASSERT_MSG(isValidModel(reverse, readIndex), "Invalid method call");
    if (firstFragment(reverse, readIndex))
    {
        return position + mateMin_/* / 2*/ - readLengths[(readIndex + 1) % 2];
    }
    else
    {
        return position - mateMax_/* * 2*/ + readLengths[readIndex];
    }
}

long TemplateLengthStatistics::mateMaxPosition(
    const unsigned readIndex,
    const bool reverse,
    const long position,
    const unsigned *readLengths) const
{
    ISAAC_ASSERT_MSG(isValidModel(reverse, readIndex), "Invalid method call");
    if (firstFragment(reverse, readIndex))
    {
        return position + mateMax_/* * 2*/ - readLengths[(readIndex + 1) % 2];
    }
    else
    {
        return position - mateMin_/* / 2*/ + readLengths[readIndex];
    }
}



TemplateLengthDistribution::TemplateLengthDistribution(int mateDriftRange)
    : mateDriftRange_(mateDriftRange)
    , templateCount_(0)
    , uniqueCount_(0)
    , count_(0)
    , histograms_(TemplateLengthStatistics::InvalidAlignmentModel)
{
}

void TemplateLengthDistribution::unreserve()
{
    std::vector<std::vector<unsigned> >().swap(histograms_);
    std::vector<unsigned>().swap(lengthList_);
}

void TemplateLengthDistribution::clear()
{
    stats_.clear();
    templateCount_ = 0;
    uniqueCount_ = 0;
    count_ = 0;
    std::for_each(histograms_.begin(), histograms_.end(), std::mem_fun_ref(&std::vector<unsigned>::clear));
    lengthList_.clear();
}


void TemplateLengthDistribution::reset(const std::vector<reference::Contig> &contigList,
                                     const flowcell::ReadMetadataList &readMetadataList)
{
    clear();
}


inline FragmentMetadataList::const_iterator findUniqueOrUniquePerfectAlignment(
    const FragmentMetadataList &list)
{
    // if the alignment is unique, assume it is good enough.
    if (1 == list.size())
    {
        return list.begin();
    }

    // if there are too many candidates, just skip it
    if (list.size() > 10)
    {
        return list.end();
    }

    // if there are multiple candidates, only look at the perfect ones and only if the perfect ones are unique.
    
    // Notice that there is still a problem

    const FragmentMetadataList::const_iterator ret = std::find_if(list.begin(), list.end(), !boost::bind(&FragmentMetadata::getEditDistance, _1));
    if (list.end() == ret)
    {
        return list.end();
    }

    const FragmentMetadataList::const_iterator it = std::find_if(ret + 1, list.end(), !boost::bind(&FragmentMetadata::getEditDistance, _1));
    if (list.end() == it)
    {
        return ret;
    }

    return list.end();
}

bool TemplateLengthDistribution::addTemplate(const std::vector<std::vector<FragmentMetadata> > &fragments)
{
    ISAAC_ASSERT_MSG(2 == fragments.size(), "Maximum of two fragments per template is supported");
    // discard templates where at least on fragment didn't align
    if (fragments[0].empty() || fragments[1].empty())
    {
        return stats_.isStable();
    }
    ++templateCount_;
    // discard templates where the alignment is not unique on both fragments
    const FragmentMetadataList::const_iterator perfect0 = findUniqueOrUniquePerfectAlignment(fragments[0]);
    const FragmentMetadataList::const_iterator perfect1 = findUniqueOrUniquePerfectAlignment(fragments[1]);
    if (fragments[0].end() == perfect0 || fragments[1].end() == perfect1)
    {
        return stats_.isStable();
    }
    ++uniqueCount_;
    // discard templates that span across several contigs
    if (perfect0->contigId != perfect1->contigId)
    {
        return stats_.isStable();
    }
    // discard fragments that ar not completely contained inside the contig
    // this is identified by the presence of leading or trailing inserts
    assert(perfect0->cigarBuffer == perfect1->cigarBuffer);

    // calculate the lenght of the template
    const unsigned long length = stats_.getLength(*perfect0, *perfect1);
    // discard excessively long templates
    if (length > TemplateLengthStatistics::TEMPLATE_LENGTH_THRESHOLD)
    {
        return stats_.isStable();
    }
    // update the histogram for the appropriate alignment model
    const TemplateLengthStatistics::AlignmentModel am = TemplateLengthStatistics::alignmentModel(*perfect0, *perfect1);
    if (TemplateLengthStatistics::InvalidAlignmentModel != am)
    {
        histograms_[am].push_back(length);
        // calculate the alignment statistics if appropriate
        ++count_;
        if (0 == (count_ % UPDATE_FREQUENCY))
        {
            // save the old statistics
            TemplateLengthStatistics oldStats = stats_;
            // get the new statistics
            updateStatistics();
            // check if we reached stability
            if (oldStats.getMin() == stats_.getMin() &&
                oldStats.getMedian() == stats_.getMedian() &&
                oldStats.getMax() == stats_.getMax() &&
                oldStats.getLowStdDev() == stats_.getLowStdDev() &&
                oldStats.getHighStdDev() == stats_.getHighStdDev())
            {
                stats_.setStable(true);
            }
        }
    }
    return stats_.isStable();
}

bool TemplateLengthDistribution::finalize()
{
    // save the old statistics
    TemplateLengthStatistics oldStats = stats_;
    // get the new statistics
    updateStatistics();
    // check if we reached stability
    if (oldStats.getMin() == stats_.getMin() &&
        oldStats.getMedian() == stats_.getMedian() &&
        oldStats.getMax() == stats_.getMax() &&
        oldStats.getLowStdDev() == stats_.getLowStdDev() &&
        oldStats.getHighStdDev() == stats_.getHighStdDev())
    {
        stats_.setStable(true);
    }
    return stats_.isStable();
}


} // namespace alignment
} // namespace isaac
