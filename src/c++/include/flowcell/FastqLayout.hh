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
 ** \file BamLayout.hh
 **
 ** Specialization of Layout for bam flowcell.
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_FLOWCELL_FASTQ_LAYOUT_HH
#define iSAAC_FLOWCELL_FASTQ_LAYOUT_HH

#include "flowcell/Layout.hh"

namespace isaac
{
namespace flowcell
{

namespace fastq
{
    static const unsigned READ_NUMBER_MAX = 2;


    void getFastqFilePath(
        const boost::filesystem::path &baseCallsPath,
        const unsigned lane,
        const unsigned read,
        const bool compressed,
        boost::filesystem::path &result);
}

struct FastqFilePathAttributeTag
{
    typedef boost::filesystem::path value_type;
    friend std::ostream &operator << (std::ostream &os, const FastqFilePathAttributeTag &tag){return os << "FastqFilePathAttributeTag";}
};

struct FastqBaseQ0
{
    typedef char value_type;
    friend std::ostream &operator << (std::ostream &os, const FastqBaseQ0 &tag){return os << "FastqBaseQ0";}
};


template<>
const boost::filesystem::path & Layout::getLaneReadAttribute<Layout::Fastq, FastqFilePathAttributeTag>(
    const unsigned lane, const unsigned read, boost::filesystem::path &result) const;

template<>
const FastqBaseQ0::value_type & Layout::getAttribute<Layout::Fastq, FastqBaseQ0>(
    FastqBaseQ0::value_type &result) const;

template<>
inline boost::filesystem::path Layout::getLongestAttribute<Layout::Fastq, FastqFilePathAttributeTag>() const
{
    boost::filesystem::path filtersFilePath;
    getLaneReadAttribute<Layout::Fastq, FastqFilePathAttributeTag>(laneNumberMax_, fastq::READ_NUMBER_MAX, filtersFilePath);
    return filtersFilePath;
}



} // namespace flowcell
} // namespace isaac

#endif // #ifndef iSAAC_FLOWCELL_FASTQ_LAYOUT_HH
