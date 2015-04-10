/**
 ** Isaac Genome Alignment Software
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
 ** \file MergeAnnotationsOptions.cpp
 **
 ** Command line options for mergeAnnotations
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_OPTIONS_MERGE_ANNOTATIONS_OPTIONS_HH
#define iSAAC_OPTIONS_MERGE_ANNOTATIONS_OPTIONS_HH

#include <string>
#include <boost/filesystem.hpp>

#include "common/Program.hh"

namespace isaac
{
namespace options
{

class MergeAnnotationsOptions  : public common::Options
{
public:
    boost::filesystem::path referenceGenome;
    std::vector<boost::filesystem::path> filesToMerge_;
    boost::filesystem::path outputFilePath_;
    bool outputReadLengths_;
    std::vector<std::string> aggregateFunction_;
    std::string mergedType_;

public:
    MergeAnnotationsOptions();

    common::Options::Action parse(int argc, char *argv[]);

private:
    std::string usagePrefix() const {return "mergeReferences";}
    void postProcess(boost::program_options::variables_map &vm);
    void verifyMandatoryPaths(boost::program_options::variables_map &vm);
};

} // namespace options
} // namespace isaac

#endif // #ifndef iSAAC_OPTIONS_MERGE_ANNOTATIONS_OPTIONS_HH
