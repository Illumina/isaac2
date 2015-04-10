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
 ** \file extractNeighbors.cpp
 **
 ** \brief Extract neighbor flags from kmers and store them in a file where
 **        each bit would correspond to a position in the reference
 **
 ** \author Roman Petrovski
 **/

#include "common/Debug.hh"
#include "options/ExtractNeighborsFromAnnotationOptions.hh"
#include "workflow/ExtractNeighborsFronAnnotationWorkflow.hh"

void extractNeighborsFromAnnotation(const isaac::options::ExtractNeighborsFromAnnotationOptions &options)
{
    isaac::workflow::ExtractNeighborsFromAnnotationWorkflow workflow(
        options.contigsXmlPath_,
        options.highAnnotationFilePath_,
        options.bitsPerValue_,
        options.outputFilePath_
        );

    workflow.run();
}

int main(int argc, char *argv[])
{
    isaac::common::run(extractNeighborsFromAnnotation, argc, argv);
}

