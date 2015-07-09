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
 ** \file ReorderReferenceOptions.cpp
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
#include <boost/foreach.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/lexical_cast.hpp>

#include "common/Debug.hh"
#include "common/Exceptions.hh"
#include "options/ReorderReferenceOptions.hh"



namespace isaac
{
namespace options
{

namespace bpo = boost::program_options;
namespace bfs = boost::filesystem;
using common::InvalidOptionException;
using boost::format;

ReorderReferenceOptions::ReorderReferenceOptions()
    : basesPerLine_(70)
{
    namedOptions_.add_options()
        ("reference-genome,r"       , bpo::value<bfs::path>(&sortedReferenceMetadata_),
                "Full path to the reference genome XML descriptor."
            )
        ("order"       , bpo::value<std::string>(&newOrderString_),
                "Comma-separated list of contig names in the order in which they will appear in the new .fa file."
            )
        ("output-directory,d"       , bpo::value<bfs::path>(&newDataDirectory_),
                "Path for the reordered fasta and annotation files."
            )
        ("bases-per-line,b" , bpo::value<unsigned>(&basesPerLine_)->default_value(basesPerLine_),
                "Number of bases per line to print into .fa file."
            )
        ("output-xml,x"       , bpo::value<bfs::path>(&newXmlPath_),
                "Path for the new xml file."
            );
}

common::Options::Action ReorderReferenceOptions::parse(int argc, char *argv[])
{
    const std::vector<std::string> allOptions(argv, argv + argc);
    common::Options::Action ret = common::Options::parse(argc, argv);
    if (RUN == ret)
    {
        ISAAC_THREAD_CERR << "argc: " << argc << " argv: " << boost::join(allOptions, " ") << std::endl;
    }
    return ret;
}

void ReorderReferenceOptions::postProcess(bpo::variables_map &vm)
{
    if(vm.count("help") ||  vm.count("version"))
    {
        return;
    }

    const std::vector<std::string> requiredOptions = boost::assign::list_of("output-directory")("output-xml")("reference-genome");
    BOOST_FOREACH(const std::string &required, requiredOptions)
    {
        if(!vm.count(required))
        {
            const format message = format("\n   *** The '%s' option is required ***\n") % required;
            BOOST_THROW_EXCEPTION(InvalidOptionException(message.str()));
        }
    }

    newXmlPath_ = boost::filesystem::absolute(newXmlPath_);

    if (boost::filesystem::exists(newXmlPath_))
    {
        BOOST_THROW_EXCEPTION(InvalidOptionException((format("\n   *** Output file already exists. Please delete the existing file or change --new-name parameter.: %s ***\n") % newXmlPath_.string()).str()));
    }

    newDataDirectory_ = boost::filesystem::absolute(newDataDirectory_);

    if (!boost::filesystem::exists(newDataDirectory_))
    {
        BOOST_THROW_EXCEPTION(InvalidOptionException((format("\n   *** Output folder does not exist.: %s ***\n") % newDataDirectory_.string()).str()));
    }

    if (!newOrderString_.empty())
    {
        boost::split_regex(newOrder_, newOrderString_, boost::regex(","));
    }
}




} //namespace options
} // namespace isaac
