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
 ** \file ExtractNeighborsOptions.cpp
 **
 ** Command line options for 'isaac-reorder-reference'
 **
 ** \author Roman Petrovski
 **/

#include <algorithm>
#include <string>
#include <vector>
#include <ostream>
#include <fstream>

#include <boost/algorithm/string/regex.hpp>
#include <boost/assign.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/foreach.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/lexical_cast.hpp>

#include "common/Debug.hh"
#include "common/Exceptions.hh"
#include "options/ExtractNeighborsOptions.hh"



namespace isaac
{
namespace options
{

namespace bpo = boost::program_options;
namespace bfs = boost::filesystem;
using common::InvalidOptionException;

ExtractNeighborsOptions::ExtractNeighborsOptions() : seedLength(32)
{
    namedOptions_.add_options()
        ("high-repeats-file,h"       , bpo::value<bfs::path>(&highRepeatsFilePath_),
                "Path for the output file where high repeat positions are flagged."
            )
        ("reference-genome,r"       , bpo::value<bfs::path>(&sortedReferenceMetadata_),
                "Full path to the reference genome XML descriptor."
            )
        ("output-file,o"       , bpo::value<bfs::path>(&outputFilePath_),
                "Path for the output file where neighbor positions are flagged."
            )
        ("seed-length,s",  bpo::value<unsigned int>(&seedLength)->default_value(seedLength),
                          "Length of reference k-mer in bases. 16,28,30,32,34,36,64 supported."
            )
            ;
}

common::Options::Action ExtractNeighborsOptions::parse(int argc, char *argv[])
{
    const std::vector<std::string> allOptions(argv, argv + argc);
    common::Options::Action ret = common::Options::parse(argc, argv);
    if (RUN == ret)
    {
        ISAAC_THREAD_CERR << "argc: " << argc << " argv: " << boost::join(allOptions, " ") << std::endl;
    }
    return ret;
}

void ExtractNeighborsOptions::postProcess(bpo::variables_map &vm)
{
    if(vm.count("help") ||  vm.count("version"))
    {
        return;
    }

    const std::vector<std::string> requiredOptions = boost::assign::list_of("output-file")("reference-genome");
    BOOST_FOREACH(const std::string &required, requiredOptions)
    {
        if(!vm.count(required))
        {
            const boost::format message = boost::format("\n   *** The '%s' option is required ***\n") % required;
            BOOST_THROW_EXCEPTION(InvalidOptionException(message.str()));
        }
    }

    outputFilePath_ = boost::filesystem::absolute(outputFilePath_);
    if (!highRepeatsFilePath_.empty())
    {
        highRepeatsFilePath_ = boost::filesystem::absolute(highRepeatsFilePath_);
    }

    if (boost::filesystem::exists(outputFilePath_))
    {
        BOOST_THROW_EXCEPTION(InvalidOptionException((boost::format("\n   *** Output file already exists. : %s ***\n") % outputFilePath_.string()).str()));
    }
    if (!highRepeatsFilePath_.empty() && boost::filesystem::exists(highRepeatsFilePath_))
    {
        BOOST_THROW_EXCEPTION(InvalidOptionException((boost::format("\n   *** Output file already exists. : %s ***\n") % highRepeatsFilePath_.string()).str()));
    }

    if (16 != seedLength &&
        28 != seedLength &&
        30 != seedLength &&
        32 != seedLength &&
        34 != seedLength &&
        36 != seedLength &&
        64 != seedLength)
    {
        BOOST_THROW_EXCEPTION(InvalidOptionException("\n   *** --seed-length other than 16,28,30,32,34,36 and 64 is not supported. ***\n"));
    }
}




} //namespace options
} // namespace isaac
