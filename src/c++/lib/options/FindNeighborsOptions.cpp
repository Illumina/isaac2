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
 ** \file FindNeighborsOptions.cpp
 **
 ** \brief See FindNeighborsOptions.hh.
 **
 ** \author Come Raczy
 **/

#include <string>
#include <vector>

#include <boost/assign.hpp>
#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <boost/thread.hpp>

#include "common/Exceptions.hh"
#include "oligo/Kmer.hh"
#include "options/FindNeighborsOptions.hh"
#include "reference/ReferenceKmer.hh"

namespace isaac
{
namespace options
{

namespace bpo = boost::program_options;
namespace bfs = boost::filesystem;

FindNeighborsOptions::FindNeighborsOptions()
    : seedLength(32)
    , maskWidth(0)
    , neighborhoodWidth(2)
    , mask(0)
    , referenceGenome("")
    , jobs(boost::thread::hardware_concurrency())
{
    namedOptions_.add_options()
        ("reference-genome,r",  bpo::value<bfs::path>(&referenceGenome),
                          "The input 'sorted-reference.xml' file")
        ("jobs,j", bpo::value<unsigned>(&jobs)->default_value(jobs),
                          "Maximum number of compute threads to run in parallel. Parallel sorting will use all cores regardless.")
        ("mask,m",         bpo::value<unsigned long>(&mask),
                          "mask used to filter the k-mers counted by this process (must be strictly less than 2^mask-width")
        ("mask-width,w",   bpo::value<unsigned int>(&maskWidth)->default_value(maskWidth),
                          "Width in bits of the mask used to split the sorted files")
        ("neighborhood-distance,d",   bpo::value<unsigned int>(&neighborhoodWidth)->default_value(neighborhoodWidth),
                          "Maximum hamming distance to consider the two k-mers as neighbors")
        ("output-file,o",  bpo::value<bfs::path>(&outputFile),
                          "The output file path")
        ("seed-length,s",  bpo::value<unsigned int>(&seedLength)->default_value(seedLength),
                          "Length of reference k-mer in bases. 64 or 32 is supported.")
        ;
}

common::Options::Action FindNeighborsOptions::parse(int argc, char *argv[])
{
    const std::vector<std::string> allOptions(argv, argv + argc);
    common::Options::Action ret = common::Options::parse(argc, argv);
    if (RUN == ret)
    {
        ISAAC_THREAD_CERR << "argc: " << argc << " argv: " << boost::join(allOptions, " ") << std::endl;
    }
    return ret;
}


void FindNeighborsOptions::postProcess(bpo::variables_map &vm)
{
    if(vm.count("help"))
    {
        return;
    }
    using isaac::common::InvalidOptionException;
    using boost::format;
    const std::vector<std::string> requiredOptions = boost::assign::list_of("reference-genome")("output-file");
    BOOST_FOREACH(const std::string &required, requiredOptions)
    {
        if(!vm.count(required))
        {
            const format message = format("\n   *** The '%s' option is required ***\n") % required;
            BOOST_THROW_EXCEPTION(InvalidOptionException(message.str()));
        }
    }
    typedef std::pair<bfs::path *, std::string> PathOption;
    const std::vector<PathOption> pathOptions = boost::assign::list_of
        (PathOption(&referenceGenome, "reference-genome"))
        ;
    BOOST_FOREACH(const PathOption &pathOption, pathOptions)
    {
        if(pathOption.first->empty())
        {
            const format message = format("\n   *** The '%s' can't be empty (use '.' for current directory) ***\n") % pathOption.second;
            BOOST_THROW_EXCEPTION(InvalidOptionException(message.str()));
        }
    }
    const std::vector<PathOption> existingPaths = boost::assign::list_of
        (PathOption(&referenceGenome, "reference-genome"))
        ;
    BOOST_FOREACH(const PathOption &pathOption, existingPaths)
    {
        if(!exists(*pathOption.first))
        {
            const format message = format("\n   *** The '%s' does not exist: %s ***\n") % pathOption.second % *pathOption.first;
            BOOST_THROW_EXCEPTION(InvalidOptionException(message.str()));
        }
    }
    const unsigned int maskCount = (1 << maskWidth);
    if(maskCount <= mask)
    {
        const format message = format("\n   *** The mask must be strictly less than %d: mask = %d ***\n") % maskCount % mask;
        BOOST_THROW_EXCEPTION(InvalidOptionException(message.str()));
    }

    if (!reference::isSupportedKmerLength(seedLength))
    {
        const format message = format("\n   *** The seed-length must be dividible by 4. Got: %d ***\n") % seedLength;
        BOOST_THROW_EXCEPTION(InvalidOptionException(message.str()));
    }

}

} // namespace options
} // namespace isaac
