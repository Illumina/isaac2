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
 ** \file bpbToWig.cpp
 **
 ** \brief Merges multiple reference metadata files into one
 **
 ** \author Roman Petrovski
 **/

#include "common/Debug.hh"
#include "options/MergeAnnotationsOptions.hh"
#include "reference/AnnotationsMerger.hh"
#include "reference/KUniqueness.hh"
#include "workflow/MergeAnnotationsWorkflow.hh"

void mergeAnnotations(const isaac::options::MergeAnnotationsOptions &options)
{
    isaac::workflow::MergeAnnotationsWorkflow workflow(
        options.referenceGenome,
        options.filesToMerge_,
        options.outputFilePath_,
        options.aggregateFunction_,
        options.mergedType_);
    workflow.run();
}

int main(int argc, char *argv[])
{
    isaac::common::run(mergeAnnotations, argc, argv);
}

