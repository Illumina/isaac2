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
 ** \file Contig.hh
 **
 ** \brief Definition of a contig
 **
 ** \author Come Raczy
 **/

#ifndef iSAAC_REFERENCE_CONTIG_HH
#define iSAAC_REFERENCE_CONTIG_HH

#include <iostream>
#include <string>
#include <vector>

namespace isaac
{
namespace reference
{

struct Contig
{
    unsigned index_;
    std::string name_;
    std::vector<char> forward_;

    Contig(const unsigned index, const std::string &name) : index_(index), name_(name){;}
    size_t getLength() const {return forward_.size();}

    friend std::ostream &operator <<(std::ostream &os, const Contig &contig)
    {
        return os << "Contig(" << contig.index_ << "," << contig.name_ << "," << contig.forward_.size() << ")";
    }
};

typedef std::vector<reference::Contig> ContigList;
typedef std::vector<reference::ContigList> ContigLists;

/// Total length of all the contigs of a genome
size_t genomeLength(const std::vector<Contig> &contigList);

} // namespace reference
} // namespace isaac

#endif // #ifndef iSAAC_REFERENCE_CONTIG_HH
