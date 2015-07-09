/**
 ** Isaac Genome Alignment Software
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
 ** \file WigLoader.hh
 **
 ** \brief Loader for genomic position data stored in Wig format http://genome.ucsc.edu/goldenPath/help/wiggle.html
 **
 ** \author Roman Petrovski
 **/

#ifndef ISAAC_IO_WIG_LOADER_HH
#define ISAAC_IO_WIG_LOADER_HH

#include <fstream>
#include <vector>

#include <boost/filesystem.hpp>
#include <boost/noncopyable.hpp>

#include "reference/ReferencePosition.hh"
#include "reference/SortedReferenceMetadata.hh"

namespace isaac
{
namespace io
{

class WigReader: boost::noncopyable
{
    std::istream &is_;
public:
    WigReader(std::istream &is) : is_(is){}

    template <typename InsertF>
    void read(
        const reference::SortedReferenceMetadata::Contigs &contigs,
        InsertF insert);
};

class WigLoader: boost::noncopyable
{
    boost::filesystem::path filePath_;
    std::ifstream is_;
public:
    WigLoader(const boost::filesystem::path &filePath);

    template <typename InsertF>
    void load(
        const reference::SortedReferenceMetadata::Contigs &contigs,
        InsertF insert)
    {
        ISAAC_THREAD_CERR << "Loading wigh file: " << filePath_ << std::endl;
        WigReader reader(is_);
        reader.read(contigs, insert);
        ISAAC_THREAD_CERR << "Loading wigh file done: " << filePath_ << std::endl;
    }
};

inline std::pair<std::string , unsigned> referenceContigToNameId(const reference::SortedReferenceMetadata::Contig &contig)
{
    return std::make_pair(contig.name_, contig.karyotypeIndex_);
}

inline bool compareContigName(
    const std::pair<std::string , unsigned> &left,
    const std::pair<std::string , unsigned> &right)
{
    return left.first < right.first;
}

inline bool compareContigNameWithName(
    const std::pair<std::string , unsigned> &left,
    const std::string & name)
{
    return left.first < name;
}

template <typename InsertF>
void WigReader::read(
    const reference::SortedReferenceMetadata::Contigs &contigs,
    InsertF insert)
{
    typedef std::vector<std::pair<std::string , unsigned> > ContigNameMap;
    ContigNameMap contigNameMap;
    contigNameMap.reserve(contigs.size());
    std::transform(contigs.begin(), contigs.end(), std::back_inserter(contigNameMap), referenceContigToNameId);
    std::sort(contigNameMap.begin(), contigNameMap.end(), compareContigName);
    std::string line;
    unsigned contigId = -1U;
    while (std::getline(is_, line) || is_.eof())
    {
//        ISAAC_THREAD_CERR << line << std::endl;
        static const std::string NEW_CONTIG("variableStep chrom=");
        if (NEW_CONTIG == line.substr(0, NEW_CONTIG.length()))
        {
            const std::string name = line.substr(NEW_CONTIG.length(), line.find(NEW_CONTIG.length(), ' '));
            ContigNameMap::const_iterator contigNameIt = std::lower_bound(
                contigNameMap.begin(), contigNameMap.end(), name, compareContigNameWithName);
            if (contigNameMap.end() == contigNameIt)
            {
                const boost::format message = boost::format("Unexpected contig name found while reading wig file: name") % name;
                BOOST_THROW_EXCEPTION(common::IoException(EINVAL, message.str()));
            }
            contigId = contigNameIt->second;
            ISAAC_THREAD_CERR << line << " Name: " << contigNameIt->first << " Id: " << contigNameIt->second << std::endl;
        }
        else if (!line.empty())
        {
            const std::string positionStr = line.substr(0, line.find_first_of(" \t"));
//            ISAAC_THREAD_CERR << "positionStr:" << positionStr << std::endl;
            unsigned long pos = boost::lexical_cast<unsigned long>(positionStr) - 1; // wig is 1-based
            const std::string valueStr = line.substr(positionStr.length());
            long value = boost::lexical_cast<long>(valueStr.substr(valueStr.find_first_not_of(" \t")));
//            ISAAC_THREAD_CERR << "valueStr:" << valueStr << "value:" << value <<std::endl;
            insert(std::make_pair(reference::ReferencePosition(contigId, pos), value));
        }
        if (is_.eof())
        {
            break;
        }
    }

    if (!is_.eof())
    {
        const boost::format message = boost::format("Failed to reach the end of the file: ") % strerror(errno);
        BOOST_THROW_EXCEPTION(common::IoException(errno, message.str()));
    }
}


} // namespace reference
} //namespace isaac

#endif // #ifndef ISAAC_IO_WIG_LOADER_HH
