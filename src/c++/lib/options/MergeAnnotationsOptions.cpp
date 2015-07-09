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
 ** \file MergeAnnotationsOptions.cpp
 **
 ** Command line options for 'mergeAnnotations'
 **
 ** \author Roman Petrovski
 **/

#include <boost/assign.hpp>
#include <boost/foreach.hpp>

#include "common/Debug.hh"
#include "common/Exceptions.hh"
#include "options/MergeAnnotationsOptions.hh"



namespace isaac
{
namespace options
{

namespace bpo = boost::program_options;
namespace bfs = boost::filesystem;
using common::InvalidOptionException;

MergeAnnotationsOptions::MergeAnnotationsOptions() :
    outputReadLengths_(false),
    aggregateFunction_()
{
    namedOptions_.add_options()
        ("aggregate-function"    , bpo::value<std::vector<std::string> >(&aggregateFunction_),
                "\nmin  - minimum of the values is produced."
                "\nmax  - maximum of the values is produced."
                "\nsum  - sum of the values is produced."
                "\nmask - Assuming first input is KUL annotation and the subsequent inputs are neighbor counts, the "
                         "result will the KUL annotation with positions where the neighbor counts are low downgraded "
                         "to 32."
            )
        ("input-file,i"    , bpo::value<std::vector<bfs::path> >(&filesToMerge_),
                "Paths of the files to be merged. First entry can be a initialization function name instead of file path"
            )
        ("output-file,o"       , bpo::value<bfs::path>(&outputFilePath_), "Path for the output file."
            )
        ("output-read-lengths"    , bpo::value<bool>(&outputReadLengths_)->default_value(outputReadLengths_),
                "For each annotated position, find minimum length of a k-unique read that would overlap it."
            )
        ("merged-type"    , bpo::value<std::string>(&mergedType_),
                "Type of the output data (kul or neighbor-counts)."
            )
        ("reference-genome,r",  bpo::value<bfs::path>(&referenceGenome),
                          "The input 'sorted-reference.xml' file")
            ;
}

common::Options::Action MergeAnnotationsOptions::parse(int argc, char *argv[])
{
    const std::vector<std::string> allOptions(argv, argv + argc);
    common::Options::Action ret = common::Options::parse(argc, argv);
    if (RUN == ret)
    {
        ISAAC_THREAD_CERR << "argc: " << argc << " argv: " << boost::join(allOptions, " ") << std::endl;
    }
    return ret;
}

void MergeAnnotationsOptions::postProcess(bpo::variables_map &vm)
{
    if(vm.count("help") ||  vm.count("version"))
    {
        return;
    }

    const std::vector<std::string> requiredOptions = boost::assign::list_of("input-file")("output-file")("reference-genome")("merged-type");
    BOOST_FOREACH(const std::string &required, requiredOptions)
    {
        if(!vm.count(required))
        {
            const boost::format message = boost::format("\n   *** The '%s' option is required ***\n") % required;
            BOOST_THROW_EXCEPTION(InvalidOptionException(message.str()));
        }
    }

    if (aggregateFunction_.empty() && filesToMerge_.size() > 1)
    {
        const boost::format message = boost::format("\n   *** --aggregate-function must be specified at least once. ***\n");
        BOOST_THROW_EXCEPTION(InvalidOptionException(message.str()));
    }
    if (aggregateFunction_.size() >= filesToMerge_.size())
    {
        const boost::format message = boost::format("\n   *** too many --aggregate-function parameters provided. ***\n");
        BOOST_THROW_EXCEPTION(InvalidOptionException(message.str()));
    }
    if (filesToMerge_.size() > 1)
    {
        aggregateFunction_.resize(std::max(filesToMerge_.size(), aggregateFunction_.size() - 1), aggregateFunction_.back());
    }

//    if ("max" == aggregateFunction_)
//    {
//    }
//    else if ("min" == aggregateFunction_)
//    {
//    }
//    else if ("sum" == aggregateFunction_)
//    {
//    }
//    else if ("zero-and-leq10" == aggregateFunction_);
//    else if ("zero-and-leq100" == aggregateFunction_);
//    else if ("fill-32" == aggregateFunction_);
//    else if ("fill-64" == aggregateFunction_);
//    else if ("fill-65535" == aggregateFunction_);
//    else if ("mask-16" == aggregateFunction_);
//    else if ("mask-20" == aggregateFunction_);
//    else if ("mask-24" == aggregateFunction_);
//    else if ("mask-28" == aggregateFunction_);
//    else if ("mask-32" == aggregateFunction_);
//    else if ("mask-36" == aggregateFunction_);
//    else if ("mask-40" == aggregateFunction_);
//    else if ("mask-44" == aggregateFunction_);
//    else if ("mask-48" == aggregateFunction_);
//    else if ("mask-52" == aggregateFunction_);
//    else if ("mask-56" == aggregateFunction_);
//    else if ("mask-60" == aggregateFunction_);
//    else if ("mask-64" == aggregateFunction_)
//    {
//    }
//    else
//    {
//        const boost::format message = boost::format("\n   *** --aggregate-function must be either 'min', 'max', 'sum' or 'mask'. Got: %s ***\n") % aggregateFunction_;
//        BOOST_THROW_EXCEPTION(InvalidOptionException(message.str()));
//    }
}




} //namespace options
} // namespace isaac
