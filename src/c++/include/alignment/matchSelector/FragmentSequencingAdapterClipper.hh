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
 ** \file FragmentSequencingAdapterClipper.hh
 **
 ** \brief Utility classes for detecting and removing fragment parts that contain
 **        sequences of the adapters
 ** 
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_ALIGNMENT_FRAGMENT_SEQUENCING_ADAPTER_CLIPPER_HH
#define iSAAC_ALIGNMENT_FRAGMENT_SEQUENCING_ADAPTER_CLIPPER_HH

#include "alignment/matchSelector/SequencingAdapter.hh"
#include "reference/Contig.hh"

namespace isaac
{
namespace alignment
{
namespace matchSelector
{

class FragmentSequencingAdapterClipper: public boost::noncopyable
{
    // percent of mismatching bases below which the flank is assumed
    // to be too good for the real adapter-containing read.
    static const unsigned TOO_GOOD_READ_MISMATCH_PERCENT = 40;

    const matchSelector::SequencingAdapterList &sequencingAdapters_;
public:
    explicit FragmentSequencingAdapterClipper(
        const matchSelector::SequencingAdapterList &sequencingAdapters):
            sequencingAdapters_(sequencingAdapters)
    {
    }

    void checkInitStrand(
        const FragmentMetadata &fragmentMetadata,
        const reference::Contig &contig);

    void clip(
        const reference::Contig &contig,
        FragmentMetadata &fragment,
        std::vector<char>::const_iterator &sequenceBegin,
        std::vector<char>::const_iterator &sequenceEnd) const;
private:
    struct SequencingAdapterRange
    {
        SequencingAdapterRange() : initialized_(false), empty_(true), unbounded_(false){}
        bool initialized_;
        bool empty_;
        bool unbounded_;
        std::vector<char>::const_iterator adapterRangeBegin_;
        std::vector<char>::const_iterator adapterRangeEnd_;
    };

    struct StrandSequencingAdapterRange
    {
        // 0 - forward range, 1 - reverse range
        matchSelector::FragmentSequencingAdapterClipper::SequencingAdapterRange strandRange_[2];
    };

    StrandSequencingAdapterRange strandAdapters_;

    static bool decideWhichSideToClip(
        const reference::Contig &contig,
        const long contigPosition,
        const std::vector<char>::const_iterator sequenceBegin,
        const std::vector<char>::const_iterator sequenceEnd,
        const SequencingAdapterRange &adapterRange,
        bool &clipBackwards);

};

} // namespace matchSelector
} // namespace alignment
} // namespace isaac

#endif // #ifndef iSAAC_ALIGNMENT_FRAGMENT_SEQUENCING_ADAPTER_CLIPPER_HH
