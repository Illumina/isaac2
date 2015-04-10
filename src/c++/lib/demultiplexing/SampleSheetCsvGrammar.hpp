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
 ** \file SampleSheetCsvGrammar.cpp
 **
 ** SampleSheet.csv grammar definition
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_DEMULTIPLEXING_SAMPLE_SHEET_CSV_GRAMMAR_HH
#define iSAAC_DEMULTIPLEXING_SAMPLE_SHEET_CSV_GRAMMAR_HH

#include "StandardSampleSheetCsvGrammar.hpp"

namespace isaac
{
namespace demultiplexing
{

namespace bs=boost::spirit;

template <typename Iterator>
struct SampleSheetCsvGrammar :
    bs::qi::grammar<Iterator, flowcell::BarcodeMetadataList()>
{
    SampleSheetCsvGrammar(const flowcell::SequencingAdapterMetadataList &defaultAdapters) :
        SampleSheetCsvGrammar::base_type(start_),
        standardGrammar_(defaultAdapters)
    {
        start_ = standardGrammar_.start_;
    }

    StandardSampleSheetCsvGrammar<Iterator> standardGrammar_;

    bs::qi::rule<Iterator, flowcell::BarcodeMetadataList()> start_;
};

} // namespace demultiplexing
} // namespace isaac

#endif //iSAAC_DEMULTIPLEXING_SAMPLE_SHEET_CSV_GRAMMAR_HH
