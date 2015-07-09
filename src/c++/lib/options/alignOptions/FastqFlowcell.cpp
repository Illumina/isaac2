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
 ** \file BclFlowcell.hh
 **
 ** Generate flowcell object out of BaseCalls/config.xml
 **
 ** \author Roman Petrovski
 **/

#include <boost/foreach.hpp>

#include <boost/algorithm/string/regex.hpp>

#include "alignOptions/UseBasesMaskOption.hh"
#include "flowcell/FastqLayout.hh"
#include "io/FastqReader.hh"

#include "FastqFlowcell.hh"
#include "SeedDescriptorOption.hh"


namespace isaac
{
namespace options
{
namespace alignOptions
{

class CasavaFastqParser : boost::noncopyable
{
    const io::FastqReader& fastq_;
public:
    CasavaFastqParser(const io::FastqReader& fastq) :
        fastq_(fastq){}
    std::string parseFlowcellId() const
    {
        if(!fastq_.hasData())
        {
            return "";
        }
        io::FastqReader::IteratorPair header = fastq_.getHeader();
        static const std::string fastqHeaderDelimiters(": ");

        if ('@' != *header.first)
        {
            BOOST_THROW_EXCEPTION(io::FastqFormatException((boost::format("Fastq header must begin with @: %s") %
                std::string(header.first, header.second)).str()));
        }
        // separator between instrument name and run number
        io::FastqReader::IteratorPair::first_type flowcellIdBegin =
            std::find_first_of(header.first + 1, header.second,
                               fastqHeaderDelimiters.begin(), fastqHeaderDelimiters.end());

        if (header.second != flowcellIdBegin)
        {
            // separator between run number and flowcell id
            flowcellIdBegin = std::find_first_of(flowcellIdBegin + 1, header.second,
                                    fastqHeaderDelimiters.begin(), fastqHeaderDelimiters.end());

            if (header.second != flowcellIdBegin)
            {
                // separator between flowcell id and lane number
                const io::FastqReader::IteratorPair::first_type end =
                    std::find_first_of(flowcellIdBegin + 1, header.second,
                                       fastqHeaderDelimiters.begin(), fastqHeaderDelimiters.end());

                return std::string(flowcellIdBegin + 1, end);
            }
        }
        return "";
    }
    unsigned parseReadLength() const
    {
        if(!fastq_.hasData())
        {
            return 0;
        }
        return fastq_.getReadLength();
    }

    unsigned parseReadNameLength() const
    {
        if(!fastq_.hasData())
        {
            return 0;
        }
        std::string name;
        // If your fastq names are longer than that, seek professional help.
        fastq_.extractReadName(1024, std::back_inserter(name));
        return std::distance(
            name.begin(),
            // we're interested in the part of the fastq header that ends at first space.
            // CASAVA names are pretty stable except for cluster position. That one can theoretically wiggle between 1
            // and 6 characters [0;10000] * 10 + 1000. So, assuming we're unlucky and the record we get has positions
            // in low 0's, then pad it to the maximum possible.
            // b.t.w. the getName will pad the name with 0 to 1024 length
            std::find(name.begin(), name.end(), 0)) + 10;
    }
};


FastqFlowcell::FastqPathPairList FastqFlowcell::findFastqPathPairs(
    const bool compressed,
    const unsigned laneNumberMax,
    const boost::filesystem::path &baseCallsDirectory)
{
    FastqPathPairList ret;

    for (unsigned lane = 1; laneNumberMax >= lane; ++lane)
    {
        boost::filesystem::path r1Path;
        flowcell::fastq::getFastqFilePath(baseCallsDirectory, lane, 1, compressed, r1Path);

        FastqPathPair p;
        if (boost::filesystem::exists(r1Path))
        {
            p.r1Path_ = r1Path;
        }

        boost::filesystem::path r2Path;
        flowcell::fastq::getFastqFilePath(baseCallsDirectory, lane, 2, compressed, r2Path);

        if (boost::filesystem::exists(r2Path))
        {
            p.r2Path_ = r2Path;
        }

        if (!p.r1Path_.empty() || !p.r2Path_.empty())
        {
            p.lane_ = lane;
            ret.push_back(p);
        }
    }

    return ret;
}

FastqFlowcellInfo FastqFlowcell::parseFastqFlowcellInfo(
    const FastqPathPair &laneFilePaths,
    const unsigned readNameLength)
{
    FastqFlowcellInfo ret;

    if (!laneFilePaths.r1Path_.empty())
    {
        io::FastqReader reader(false, 1);
        reader.open(laneFilePaths.r1Path_);
        CasavaFastqParser parser(reader);
        ret.readLengths_.first = parser.parseReadLength();

        ret.flowcellId_ = parser.parseFlowcellId();
        ret.readNameLength_ = readNameLength ? readNameLength : parser.parseReadNameLength();
    }
    if (!laneFilePaths.r2Path_.empty())
    {
        io::FastqReader reader(false, 1);
        reader.open(laneFilePaths.r2Path_);
        CasavaFastqParser parser(reader);
        ret.readLengths_.second = parser.parseReadLength();

        if (ret.flowcellId_.empty())
        {
            ret.flowcellId_ = parser.parseFlowcellId();
        }
        else if (ret.flowcellId_ != parser.parseFlowcellId())
        {
            BOOST_THROW_EXCEPTION(common::IoException(errno, (boost::format("Flowcell ID mismatch between fastq reads %s vs %s, %s, %s") %
                ret.flowcellId_ % parser.parseFlowcellId() % laneFilePaths.r1Path_ % laneFilePaths.r2Path_).str()));
        }

        if (!readNameLength && ret.readNameLength_ != parser.parseReadNameLength())
        {
            BOOST_THROW_EXCEPTION(common::IoException(errno, (boost::format("Read name length mismatch between two fastq reads %d vs %d, %s, %s") %
                ret.readNameLength_ % parser.parseReadNameLength() % laneFilePaths.r1Path_ % laneFilePaths.r2Path_).str()));
        }
    }

    ret.lanes_.push_back(laneFilePaths.lane_);

    return ret;
}

FastqFlowcellInfo FastqFlowcell::parseFastqFlowcellInfo(
    const FastqPathPairList &flowcellFilePaths,
    const bool allowVariableFastqLength,
    const unsigned readNameLength)
{
    FastqFlowcellInfo ret;
    bool flowcellInfoReady = false;
    for (FastqPathPairList::const_iterator it = flowcellFilePaths.begin();
        flowcellFilePaths.end() != it; ++it)
    {
        FastqFlowcellInfo anotherLane = parseFastqFlowcellInfo(*it, readNameLength);
        if (!flowcellInfoReady)
        {
            if (anotherLane.readLengths_.first || anotherLane.readLengths_.second)
            {
                ret = anotherLane;
                flowcellInfoReady = true;
            }
            else
            {
                ISAAC_THREAD_CERR << "WARNING: Skipping lane " << it->lane_ << " due to read length 0" << std::endl;
            }
        }
        else
        {
            if (!anotherLane.readLengths_.first && !anotherLane.readLengths_.second)
            {
                ISAAC_THREAD_CERR << "WARNING: Skipping lane " << it->lane_ << " due to read length 0" << std::endl;
            }
            else
            {
                // With allowVariableFastqLength the read lengths are normally forced by use-bases-mask
                // Ignore discrepancy here.
                if (!allowVariableFastqLength && anotherLane.readLengths_ != ret.readLengths_)
                {
                    BOOST_THROW_EXCEPTION(common::IoException(errno, (boost::format("Read lengths mismatch between lanes of the same flowcell %s vs %s") %
                        anotherLane % ret).str()));
                }

                if (anotherLane.flowcellId_ != ret.flowcellId_)
                {
                    ISAAC_THREAD_CERR << "WARNING: Flowcell id mismatch across the lanes of the same flowcell" <<
                        anotherLane << " vs " << ret << std::endl;
                }

                ret.lanes_.push_back(it->lane_);
            }
        }
    }
    ISAAC_THREAD_CERR << ret << std::endl;
    return ret;
}

flowcell::Layout FastqFlowcell::createFilteredFlowcell(
    const std::string &tilesFilter,
    const boost::filesystem::path &baseCallsDirectory,
    const bool compressed,
    const unsigned laneNumberMax,
    const unsigned readNameLength,
    std::string useBasesMask,
    const bool allowVariableFastqLength,
    const std::string &seedDescriptor,
    const unsigned seedLength,
    const reference::ReferenceMetadataList &referenceMetadataList,
    unsigned &firstPassSeeds)
{

    FastqPathPairList flowcellFilePaths = findFastqPathPairs(compressed, laneNumberMax, baseCallsDirectory);
    if (flowcellFilePaths.empty())
    {
        const boost::format message = boost::format("\n   *** Could not find any fastq lanes in: %s ***\n") %
            baseCallsDirectory;
        BOOST_THROW_EXCEPTION(common::InvalidOptionException(message.str()));
    }

    FastqFlowcellInfo flowcellInfo = parseFastqFlowcellInfo(flowcellFilePaths, allowVariableFastqLength, readNameLength);

    std::vector<unsigned int> readLengths;
    if (flowcellInfo.readLengths_.first)
    {
        readLengths.push_back(flowcellInfo.readLengths_.first);
    }
    if (flowcellInfo.readLengths_.second)
    {
        readLengths.push_back(flowcellInfo.readLengths_.second);
    }

    if ("default" == useBasesMask)
    {
        if (readLengths.size() == 1)
        {
            useBasesMask = "y*n";
        }
        else if (readLengths.size() == 2)
        {
            useBasesMask = "y*n,y*n";
        }
        else
        {
            const boost::format message =
                boost::format("\n   *** Could not guess the use-bases-mask for '%s', please supply the explicit value ***\n") %
                baseCallsDirectory.string();
            BOOST_THROW_EXCEPTION(common::InvalidOptionException(message.str()));
        }
    }

    std::vector<unsigned int> readFirstCycles;

    ParsedUseBasesMask parsedUseBasesMask;
    alignment::SeedMetadataList seedMetadataList;
    if (!readLengths.empty())
    {
        parsedUseBasesMask = parseUseBasesMask(readFirstCycles, readLengths, seedLength, useBasesMask, baseCallsDirectory);
        seedMetadataList = parseSeedDescriptor(parsedUseBasesMask.dataReads_, seedDescriptor, seedLength, firstPassSeeds);
    }

    flowcell::Layout fc(baseCallsDirectory,
                        flowcell::Layout::Fastq,
                        flowcell::FastqFlowcellData(compressed),
                        laneNumberMax,
                        flowcellInfo.readNameLength_,
                        std::vector<unsigned>(),
                        parsedUseBasesMask.dataReads_,
                        seedMetadataList, flowcellInfo.flowcellId_);

    std::string regexString(tilesFilter);
    std::replace(regexString.begin(), regexString.end(), ',', '|');
    boost::regex re(regexString);
    BOOST_FOREACH(const unsigned int lane, flowcellInfo.getLanes())
    {
        std::string laneString((boost::format("s_%d") % lane).str());
        if (boost::regex_search(laneString, re))
        {
            fc.addTile(lane, 1);
        }
    }

    return fc;
}


} // namespace alignOptions
} // namespace option
} // namespace isaac
