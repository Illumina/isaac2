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
 ** \file UseBasesMaskOption.cpp
 **
 ** Parsing of use-bases-mask
 **
 ** \author Roman Petrovski
 **/

#include <algorithm>
#include <string>
#include <vector>
#include <ostream>
#include <fstream>

#include <boost/algorithm/string.hpp>
#include <boost/assign.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/thread.hpp>
#include <boost/algorithm/string/regex.hpp>

#include "common/Debug.hh"
#include "common/Exceptions.hh"

#include "oligo/Kmer.hh"
#include "options/UseBasesMaskGrammar.hh"
#include "UseBasesMaskOption.hh"


namespace isaac
{
namespace options
{

static std::vector<unsigned int> figureReadFirstCycles(
    const std::vector<std::string > &readMasks)
{
    std::vector<unsigned int> readFirstCycles;
    readFirstCycles.push_back(1);
    BOOST_FOREACH(const std::string &readMask, readMasks)
    {
        readFirstCycles.push_back(readFirstCycles.back() + readMask.length());
    }
    // remove the last one that should not be there
    readFirstCycles.pop_back();
    return readFirstCycles;
}

static std::vector<std::string > expandUseBasesMask (
    const std::vector<unsigned int> &readLengths,
    const std::string &useBasesMask,
    const boost::filesystem::path &baseCallsDirectory)
{
    std::vector<std::string > result;
    std::string::const_iterator parseIt(useBasesMask.begin());
    const std::string::const_iterator parseEnd(useBasesMask.end());
    UseBasesMaskGrammar<std::string::const_iterator> parser(readLengths);
    if (!boost::spirit::qi::parse(parseIt, parseEnd, parser, result) ||
        parseEnd != parseIt)
    {
        const boost::format message = boost::format("\n   *** Could not parse the use-bases-mask '%s' for '%s' at: %s ***\n") %
                useBasesMask % baseCallsDirectory.string() % useBasesMask.substr(parseIt - useBasesMask.begin());
        BOOST_THROW_EXCEPTION(common::InvalidOptionException(message.str()));
    }

    ISAAC_THREAD_CERR << "use bases mask: " << boost::algorithm::join(result, ",") << "\n";
    ISAAC_THREAD_CERR << "reads parsed: " << parser.currentRead_ << "\n";

    if (result.size() != readLengths.size())
    {
        const boost::format message = boost::format("\n   *** use-bases-mask '%s' is incompatible with number of reads (%d) in %s ***\n") %
                useBasesMask % readLengths.size() % baseCallsDirectory.string();
        BOOST_THROW_EXCEPTION(common::InvalidOptionException(message.str()));
    }

    return result;
}

/**
 * \param cfgReadFirstCycles - first cycle number for each read. If empty, the cycle numbers are assigned based on
 *                             expansion of the useBasesMask
 */
ParsedUseBasesMask parseUseBasesMask (const std::vector<unsigned int> &cfgReadFirstCycles,
                                      const std::vector<unsigned int> &readLengths,
                                      const unsigned seedLength,
                                      const std::string &useBasesMask,
                                      const boost::filesystem::path &baseCallsDirectory)
{
    const std::vector<std::string > expandedUseBasesMasks = expandUseBasesMask(readLengths, useBasesMask, baseCallsDirectory);
    const std::vector<unsigned int> readFirstCycles = figureReadFirstCycles(expandedUseBasesMasks);

    ParsedUseBasesMask ret;
    unsigned dataReadOffset = 0;
    unsigned dataReadNumber = 1;
    BOOST_FOREACH(const std::string &readMask, expandedUseBasesMasks)
    {
        const unsigned readMaskIndex = &readMask - &expandedUseBasesMasks.front();
        const size_t currentReadFirstCycle(readFirstCycles.at(readMaskIndex));

        std::vector<unsigned> filteredIndexCycles;
        BOOST_FOREACH(const char &chref, readMask)
        {
            if ('i' == chref)
            {
                filteredIndexCycles.push_back(currentReadFirstCycle + &chref - &*readMask.begin());
            }
        }
        if (!filteredIndexCycles.empty())
        {
            const unsigned readIndex = ret.indexReads_.size();
            // at the moment index read numbers are not being used anywhere
            ret.indexReads_.push_back(flowcell::ReadMetadata(0, filteredIndexCycles, readIndex, -1U, currentReadFirstCycle));
            ISAAC_THREAD_CERR << "Discovered index read: " << ret.indexReads_.back() << std::endl;
        }

        std::vector<unsigned> filteredDataCycles;
        BOOST_FOREACH(const char &chref, readMask)
        {
            if ('y' == chref)
            {
                filteredDataCycles.push_back(currentReadFirstCycle + &chref - &*readMask.begin());
            }
        }
        if (!filteredDataCycles.empty())
        {
            const unsigned readIndex = ret.dataReads_.size();
            ret.dataReads_.push_back(flowcell::ReadMetadata(dataReadNumber, filteredDataCycles, readIndex, dataReadOffset, currentReadFirstCycle));
            ISAAC_THREAD_CERR << "Discovered data read: " << ret.dataReads_.back() << std::endl;
            dataReadOffset += ret.dataReads_.back().getLength();
            ++dataReadNumber;
        }
        else if (filteredIndexCycles.empty())
        {
            // this ensures read number preservation even if the data read is masked out
            ++dataReadNumber;
        }
    }

    BOOST_FOREACH(const flowcell::ReadMetadata &readMetadata, ret.dataReads_)
    {
        if (readMetadata.getLength() < seedLength)
        {
            const boost::format message = boost::format("\n   *** %s is too short: %d cycle < %d in %s ***\n") %
                readMetadata % readMetadata.getLength() % seedLength % baseCallsDirectory.string();
            BOOST_THROW_EXCEPTION(common::InvalidOptionException(message.str()));
        }
    }

    // Fix read numbers. Note: this is not possible to do right. in case of n*,y* we can't tell whether it was an index or data read masked out.
    // Just make sure if there are two data reads, they get numbered as 1 and 2. And if there is only one read, it is not numbered above 2
    if (2 == ret.dataReads_.size())
    {
        ret.dataReads_[0].setNumber(1);
        ret.dataReads_[1].setNumber(2);
    }
    else if (1 == ret.dataReads_.size())
    {
        ret.dataReads_[0].setNumber(std::min(2U, ret.dataReads_[0].getNumber()));
    }


    return ret;
}

} //namespace option
} // namespace isaac
