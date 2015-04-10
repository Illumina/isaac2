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
 ** \file SortReferenceOptions.hh
 **
 ** \brief Command line options for 'sortReference'
 **
 ** \author Come Raczy
 **/

#ifndef iSAAC_COMMON_SORT_REFERENCE_OPTIONS_HH
#define iSAAC_COMMON_SORT_REFERENCE_OPTIONS_HH

#include <string>
#include <boost/filesystem.hpp>

#include "common/Program.hh"

namespace isaac
{
namespace options
{

class SortReferenceOptions : public isaac::common::Options
{
public:
    SortReferenceOptions();
private:
    std::string usagePrefix() const {return "sortReference";}
    void postProcess(boost::program_options::variables_map &vm);
public:
    unsigned seedLength;
    unsigned maskWidth;
    unsigned long mask;
    std::string contigsXml;
    boost::filesystem::path genomeNeighborsFile;
    boost::filesystem::path outFile;
    unsigned int repeatThreshold;
    std::string parallelSort;
};

} // namespace options
} // namespace isaac

#endif // #ifndef iSAAC_COMMON_SORT_REFERENCE_OPTIONS_HH
