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
 ** \file BpbToWigOptions.cpp
 **
 ** Command line options for 'bpbToWig'
 **
 ** \author Roman Petrovski
 **/

#include <boost/assign.hpp>
#include <boost/foreach.hpp>

#include "common/Debug.hh"
#include "common/Exceptions.hh"
#include "options/BpbToWigOptions.hh"



namespace isaac
{
namespace options
{

namespace bpo = boost::program_options;
namespace bfs = boost::filesystem;
using common::InvalidOptionException;

BpbToWigOptions::BpbToWigOptions():
    wigDefaultValue_(32),
    outputFormatString_("wig"),
    bitsPerValue_(1),
    knownSitesOnly_(false),
    bedPrintAllPositions_(false)
{
    namedOptions_.add_options()
        ("bed-all-positions"   , bpo::value<bool>(&bedPrintAllPositions_)->default_value(bedPrintAllPositions_),
            "Print one line for each genomic position in bed file. When off, the range of positions with matching value are being collapsed."
            )
        ("bits-per-pos,b"    , bpo::value<unsigned>(&bitsPerValue_)->default_value(bitsPerValue_),
                "How many bits are stored in the bpb file for each genomic position."
            )
        ("input-file,i"       , bpo::value<bfs::path>(&inputFilePath_), "Path for the input file."
            )
        ("known-sites-only"   , bpo::value<bool>(&knownSitesOnly_)->default_value(knownSitesOnly_),
            "Don't print sites with non-ACGT nucleotide"
            )
        ("output-format,f"       , bpo::value<std::string>(&outputFormatString_),
                "Allowed options are:"
                "\n\tbed"
                "\n\twig"
            )
        ("reference-genome,r"    , bpo::value<bfs::path>(&sortedReferenceMetadata_),
                "Full path to the reference genome XML descriptor."
            )
        ("wig-default-value"    , bpo::value<unsigned>(&wigDefaultValue_)->default_value(wigDefaultValue_),
                "Default value to be excluded from wig file."
            )

            ;
}

common::Options::Action BpbToWigOptions::parse(int argc, char *argv[])
{
    const std::vector<std::string> allOptions(argv, argv + argc);
    common::Options::Action ret = common::Options::parse(argc, argv);
    if (RUN == ret)
    {
        ISAAC_THREAD_CERR << "argc: " << argc << " argv: " << boost::join(allOptions, " ") << std::endl;
    }
    return ret;
}

void BpbToWigOptions::postProcess(bpo::variables_map &vm)
{
    if(vm.count("help") ||  vm.count("version"))
    {
        return;
    }

    const std::vector<std::string> requiredOptions = boost::assign::list_of("input-file")("reference-genome");
    BOOST_FOREACH(const std::string &required, requiredOptions)
    {
        if(!vm.count(required))
        {
            const boost::format message = boost::format("\n   *** The '%s' option is required ***\n") % required;
            BOOST_THROW_EXCEPTION(InvalidOptionException(message.str()));
        }
    }

    inputFilePath_ = boost::filesystem::absolute(inputFilePath_);
}




} //namespace options
} // namespace isaac
