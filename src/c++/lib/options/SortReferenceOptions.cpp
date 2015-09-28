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
 ** \file SortReferenceOptions.cpp
 **
 ** Command line options for 'sortReference'
 **
 ** \author Come Raczy
 **/

#include <string>
#include <vector>
#include <boost/assign.hpp>
#include <boost/foreach.hpp>

#include "common/ParallelSort.hpp"
#include "reference/ReferenceKmer.hh"
#include "options/SortReferenceOptions.hh"

namespace isaac
{
namespace options
{

namespace bpo = boost::program_options;

SortReferenceOptions::SortReferenceOptions()
    : seedLength(32)
    , maskWidth(6)
    , mask(0)
    , repeatThreshold(1000)
{
     namedOptions_.add_options()
        ("reference-genome,r",  bpo::value<std::string>(&contigsXml),
                                "Full path to the reference genome XML descriptor. Only Contigs selection is required")
        ("genome-neighbors,n",  bpo::value<boost::filesystem::path>(&genomeNeighborsFile),
                                "Path to the file containing neighbor flags (one bit per genome file position)")
        ("mask,m",              bpo::value<unsigned long>(&mask),
                                "mask used to filter the k-mers counted by this process (must be strictly less than 2^mask-width")
        ("mask-width,w",        bpo::value<unsigned int>(&maskWidth)->default_value(maskWidth),
                                "Width in bits of the mask used to split the sorted files")
        ("repeat-threshold",    bpo::value<unsigned int>(&repeatThreshold)->default_value(repeatThreshold),
                                "k-mers occuring more than --repeat-threshold times in the genome are considered to be repeats. "
                                "Their positions are not stored. Special value of 0 forces inclusion the position of of every k-mer. "
                                "Value of 1 results in only the unique k-mer positions stored.")
        ("output-file,o",       bpo::value<boost::filesystem::path>(&outFile), "Output file path.")
        ("seed-length,s",       bpo::value<unsigned int>(&seedLength)->default_value(seedLength),
                                "Length of reference k-mer in bases. Lengths of 16,28,30,32,34,36 and 64 are supported.")
        ;
}

void SortReferenceOptions::postProcess(bpo::variables_map &vm)
{
    if(vm.count("help"))
    {
        return;
    }
    using isaac::common::InvalidOptionException;
    using boost::format;
    const std::vector<std::string> requiredOptions = boost::assign::list_of("mask")("reference-genome")("output-file");
    BOOST_FOREACH(const std::string &required, requiredOptions)
    {
        if(!vm.count(required))
        {
            const format message = format("\n   *** The '%s' option is required ***\n") % required;
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
        const format message = format("\n   *** The seed-length must be in range [%d,%d]. Got: %d ***\n") %
                                      reference::FIRST_SUPPORTED_KMER % reference::LAST_SUPPORTED_KMER % seedLength;
        BOOST_THROW_EXCEPTION(InvalidOptionException(message.str()));
    }
}

} //namespace option
} // namespace isaac
