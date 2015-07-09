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
 ** \file ExtractNeighborsFromAnnotationOptions.cpp
 **
 ** Command line options for 'isaac-reorder-reference'
 **
 ** \author Roman Petrovski
 **/

#include <algorithm>
#include <string>
#include <vector>

#include <boost/assign.hpp>
#include <boost/foreach.hpp>

#include "common/Debug.hh"
#include "common/Exceptions.hh"
#include "options/ExtractNeighborsFromAnnotationOptions.hh"
#include "reference/ReferenceKmer.hh"



namespace isaac
{
namespace options
{

namespace bpo = boost::program_options;
namespace bfs = boost::filesystem;
using common::InvalidOptionException;

ExtractNeighborsFromAnnotationOptions::ExtractNeighborsFromAnnotationOptions()
{
    namedOptions_.add_options()
        ("reference-genome,r"       , bpo::value<bfs::path>(&contigsXmlPath_),
                "Xml file containing reference contig metadata."
        )
        ("input-file,i"       , bpo::value<bfs::path>(&highAnnotationFilePath_),
                "Path for the output file where high repeat positions are flagged."
        )
        ("bits-per-pos,b"    , bpo::value<unsigned>(&bitsPerValue_),
            "How many bits are stored in the bpb file for each genomic position."
        )
        ("output-file,o"       , bpo::value<bfs::path>(&outputFilePath_),
                "Path for the output file where neighbor positions are flagged."
        )
        ;
}

common::Options::Action ExtractNeighborsFromAnnotationOptions::parse(int argc, char *argv[])
{
    const std::vector<std::string> allOptions(argv, argv + argc);
    common::Options::Action ret = common::Options::parse(argc, argv);
    if (RUN == ret)
    {
        ISAAC_THREAD_CERR << "argc: " << argc << " argv: " << boost::join(allOptions, " ") << std::endl;
    }
    return ret;
}

void ExtractNeighborsFromAnnotationOptions::postProcess(bpo::variables_map &vm)
{
    if(vm.count("help") ||  vm.count("version"))
    {
        return;
    }

    const std::vector<std::string> requiredOptions = boost::assign::list_of("output-file")("bits-per-pos")("input-file")("reference-genome");
    BOOST_FOREACH(const std::string &required, requiredOptions)
    {
        if(!vm.count(required))
        {
            const boost::format message = boost::format("\n   *** The '%s' option is required ***\n") % required;
            BOOST_THROW_EXCEPTION(InvalidOptionException(message.str()));
        }
    }

    outputFilePath_ = boost::filesystem::absolute(outputFilePath_);
    highAnnotationFilePath_ = boost::filesystem::absolute(highAnnotationFilePath_);

    if (boost::filesystem::exists(outputFilePath_))
    {
        BOOST_THROW_EXCEPTION(InvalidOptionException((boost::format("\n   *** Output file already exists. : %s ***\n") % outputFilePath_.string()).str()));
    }
    if (!boost::filesystem::exists(highAnnotationFilePath_))
    {
        BOOST_THROW_EXCEPTION(InvalidOptionException((boost::format("\n   *** Input file is missing. : %s ***\n") % highAnnotationFilePath_.string()).str()));
    }
}




} //namespace options
} // namespace isaac
