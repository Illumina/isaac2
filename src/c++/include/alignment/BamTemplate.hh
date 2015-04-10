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
 ** \file BamTemplate.hh
 **
 ** \brief DNA/RNA sequence composed of one or several Fragment(s), as defined
 ** by the SAM Format Specification
 ** 
 ** \author Come Raczy
 **/

#ifndef iSAAC_ALIGNMENT_TEMPLATE_HH
#define iSAAC_ALIGNMENT_TEMPLATE_HH

#include <iostream>
#include <numeric>
#include <vector>

#include "alignment/FragmentMetadata.hh"

namespace isaac
{
namespace alignment
{

/**
 ** \brief Container to encapsulate all the data and metadata associated to a DNA/RNA BamTemplate.
 **
 ** \sa Fragment
 ** \sa FragmentId
 ** \sa Cigar
 ** \sa FragmentMetadata
 **/
class BamTemplate
{
public:
    BamTemplate();
//    BamTemplate(const BamTemplate &that):
//        fragmentMetadataList_(that.fragmentMetadataList_),
//        alignmentScore_(that.alignmentScore_),
//        properPair_(that.properPair_)
//    {;}
    BamTemplate(const FragmentMetadata &read1, const FragmentMetadata &read2);
    /**
     ** \brief Initialization for a given cluster
     **
     ** Create the appropriate unaligned FragmentMetadata for each read in the
     ** cluster.
     **/
    void initialize(const flowcell::ReadMetadataList &tileReads, const Cluster &cluster);

    unsigned getMismatchCount() const
    {
        return std::accumulate(fragmentMetadataList_.begin(), fragmentMetadataList_.end(), 0,
                               boost::bind(std::plus<unsigned>(), _1,
                                           boost::bind(&FragmentMetadata::getMismatchCount, _2)));
    }

    unsigned getQuality() const
    {
        return std::accumulate(fragmentMetadataList_.begin(), fragmentMetadataList_.end(), 0U,
                               boost::bind(std::plus<unsigned>(), _1,
                                           boost::bind(&FragmentMetadata::getQuality, _2)));
    }

    unsigned getEditDistance() const
    {
        return std::accumulate(fragmentMetadataList_.begin(), fragmentMetadataList_.end(), 0,
                               boost::bind(std::plus<unsigned>(), _1,
                                           boost::bind(&FragmentMetadata::getEditDistance, _2)));
    }

    unsigned getTotalReadLength() const
    {
        return std::accumulate(fragmentMetadataList_.begin(), fragmentMetadataList_.end(), 0,
                               boost::bind(std::plus<unsigned>(), _1,
                                           boost::bind(&FragmentMetadata::getReadLength, _2)));
    }

    bool isUnanchored() const
    {
        return fragmentMetadataList_.end() ==
            std::find_if(fragmentMetadataList_.begin(), fragmentMetadataList_.end(),
                         boost::bind(&FragmentMetadata::getAlignmentScore, _1) != 0);
    }

    bool isKUnique() const
    {
        return std::accumulate(fragmentMetadataList_.begin(), fragmentMetadataList_.end(), true,
                               boost::bind(std::logical_and<bool>(), _1,
                                           boost::bind(&FragmentMetadata::isKUnique, _2)));
    }

    bool isRepeat() const
    {
        return std::accumulate(fragmentMetadataList_.begin(), fragmentMetadataList_.end(), true,
                               boost::bind(std::logical_and<bool>(), _1,
                                           boost::bind(&FragmentMetadata::isRepeat, _2)));
    }

    bool isHighRepeat() const
    {
        return std::accumulate(fragmentMetadataList_.begin(), fragmentMetadataList_.end(), true,
                               boost::bind(std::logical_and<bool>(), _1,
                                           boost::bind(&FragmentMetadata::isHighRepeat, _2)));
    }

    unsigned getFragmentCount() const {return fragmentMetadataList_.size();}
    const FragmentMetadata &getFragmentMetadata(unsigned fragmentIndex) const {return fragmentMetadataList_.at(fragmentIndex);}
    const FragmentMetadata &getMateFragmentMetadata(const FragmentMetadata &mate) const {return fragmentMetadataList_.at(getFragmentCount() - 1 - mate.getReadIndex());}
    FragmentMetadata &getFragmentMetadata(unsigned fragmentIndex) {return fragmentMetadataList_[fragmentIndex];}
    FragmentMetadata &getMateFragmentMetadata(const FragmentMetadata &mate) {return fragmentMetadataList_.at(getFragmentCount() - 1 - mate.getReadIndex());}
    unsigned getAlignmentScore() const {return alignmentScore_;}
    void resetAlignmentScore() {alignmentScore_ = -1U;}
    bool hasAlignmentScore() const {return -1U != alignmentScore_;}
    void setAlignmentScore(unsigned alignmentScore) {alignmentScore_ = alignmentScore;}
    bool getPassesFilter() const {return fragmentMetadataList_[0].getCluster().getPf();}
    void setProperPair(const bool properPair) {properPair_ = properPair;}
    bool isProperPair() const {return properPair_;}
    bool filterLowQualityFragments(const unsigned mapqThreshold);

    std::vector<char>::const_iterator nameBegin() const
    {
        return fragmentMetadataList_[0].getCluster().nameBegin();
    }

    std::vector<char>::const_iterator nameEnd() const
    {
        return fragmentMetadataList_[0].getCluster().nameEnd();
    }

    unsigned getNameLength() const {return std::distance(nameBegin(), nameEnd());}

    bool operator ==(const BamTemplate &right) const
    {
        ISAAC_ASSERT_MSG(getFragmentCount() == right.getFragmentCount(), "Incorrect template comparison");

        if (getFragmentMetadata(0) != right.getFragmentMetadata(0))
        {
            return false;
        }

        if (1 == getFragmentCount())
        {
            return true;
        }

        return getFragmentMetadata(1) == right.getFragmentMetadata(1);
    }

    bool operator <(const BamTemplate &right) const
    {
        ISAAC_ASSERT_MSG(getFragmentCount() == right.getFragmentCount(), "Incorrect template comparison");
        if (getFragmentMetadata(0) < right.getFragmentMetadata(0))
        {
            return true;
        }

        if (1 == getFragmentCount())
        {
            return false;
        }

        if (right.getFragmentMetadata(0) < getFragmentMetadata(0))
        {
            return false;
        }

        return getFragmentMetadata(1) < right.getFragmentMetadata(1);
    }
private:
    friend std::ostream &operator<<(std::ostream &os, const BamTemplate &bamTemplate);

    common::FiniteCapacityVector<FragmentMetadata, 2> fragmentMetadataList_;

    /**
     ** This depends on the all the pLogCorrect values for all the possible
     ** alignments for this template across the whole reference. It also takes into
     ** account the rest-of-genome correction. Value of -1U indicates unknown alignment score.
     **/
    unsigned alignmentScore_;
    bool properPair_;
};

inline std::ostream &operator<<(std::ostream &os, const BamTemplate &bamTemplate)
{
    return 2 == bamTemplate.getFragmentCount() ?
        os << "BamTemplate("
        << bamTemplate.fragmentMetadataList_.at(0) << "-" << bamTemplate.fragmentMetadataList_.at(1) << "," <<
        bamTemplate.alignmentScore_ << "as," << bamTemplate.properPair_<< ")" :
        os << "BamTemplate("
        << bamTemplate.fragmentMetadataList_.at(0) << "," <<
        bamTemplate.alignmentScore_ << "as," << bamTemplate.properPair_<< ")";
}

} // namespace alignment
} // namespace isaac

#endif // #ifndef iSAAC_ALIGNMENT_TEMPLATE_HH
