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
 ** \file extractNeighborsFromAnnotationOptions.hh
 **
 ** Command line options for extractNeighborsFromAnnotation
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_OPTIONS_EXTRACT_NEIGHBORS_FROM_ANNOTATION_OPTIONS_HH
#define iSAAC_OPTIONS_EXTRACT_NEIGHBORS_FROM_ANNOTATION_OPTIONS_HH

#include <string>
#include <boost/filesystem.hpp>

#include "common/Program.hh"

namespace isaac
{
namespace options
{

class ExtractNeighborsFromAnnotationOptions  : public common::Options
{
public:
    boost::filesystem::path contigsXmlPath_;
    unsigned bitsPerValue_;
    boost::filesystem::path highAnnotationFilePath_;
    boost::filesystem::path outputFilePath_;

public:
    ExtractNeighborsFromAnnotationOptions();

    common::Options::Action parse(int argc, char *argv[]);

private:
    std::string usagePrefix() const {return "extractNeighborsFromAnnotation";}
    void postProcess(boost::program_options::variables_map &vm);
    void verifyMandatoryPaths(boost::program_options::variables_map &vm);
};

} // namespace options
} // namespace isaac

#endif // #ifndef iSAAC_OPTIONS_EXTRACT_NEIGHBORS_FROM_ANNOTATION_OPTIONS_HH
