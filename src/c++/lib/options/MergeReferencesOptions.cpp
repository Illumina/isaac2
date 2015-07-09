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
 ** \file MergeReferencesOptions.cpp
 **
 ** Command line options for 'mergeReferences'
 **
 ** \author Roman Petrovski
 **/

#include <boost/assign.hpp>
#include <boost/foreach.hpp>

#include "common/Debug.hh"
#include "common/Exceptions.hh"
#include "options/MergeReferencesOptions.hh"



namespace isaac
{
namespace options
{

namespace bpo = boost::program_options;
namespace bfs = boost::filesystem;
using common::InvalidOptionException;

MergeReferencesOptions::MergeReferencesOptions():
    mergeAnnotations_(false),
    makeAbsolutePaths_(true)
{
    namedOptions_.add_options()
        ("input-file,i"    , bpo::value<std::vector<bfs::path> >(&filesToMerge_), "Paths of the files to be merged.")
        ("output-file,o"       , bpo::value<bfs::path>(&outputFilePath_), "Path for the output file.")
        ("make-absolute-paths" , bpo::value<bool>(&makeAbsolutePaths_)->default_value(makeAbsolutePaths_), "Make sure paths in the output xml are absolute")
        ("merge-annotations"   , bpo::value<bool>(&mergeAnnotations_)->default_value(mergeAnnotations_),

            "If more than one reference contains annotations, merge annotations by selection minimum of the two values for each position"
            );
}

common::Options::Action MergeReferencesOptions::parse(int argc, char *argv[])
{
    const std::vector<std::string> allOptions(argv, argv + argc);
    common::Options::Action ret = common::Options::parse(argc, argv);
    if (RUN == ret)
    {
        ISAAC_THREAD_CERR << "argc: " << argc << " argv: " << boost::join(allOptions, " ") << std::endl;
    }
    return ret;
}

void MergeReferencesOptions::postProcess(bpo::variables_map &vm)
{
    if(vm.count("help") ||  vm.count("version"))
    {
        return;
    }

    const std::vector<std::string> requiredOptions = boost::assign::list_of("input-file")("output-file");
    BOOST_FOREACH(const std::string &required, requiredOptions)
    {
        if(!vm.count(required))
        {
            const boost::format message = boost::format("\n   *** The '%s' option is required ***\n") % required;
            BOOST_THROW_EXCEPTION(InvalidOptionException(message.str()));
        }
    }

    outputFilePath_ = boost::filesystem::absolute(outputFilePath_);
}




} //namespace options
} // namespace isaac
