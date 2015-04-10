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
 ** \file findNeighbors.cpp
 **
 ** \brief The main for the identification of the neighbors (32-mers with 1 or 2 mismatches).
 **
 ** \author Come Raczy
 **/

#include "options/FindNeighborsOptions.hh"
#include "workflow/NeighborsFinderWorkflow.hh"


void findNeighbors(const isaac::options::FindNeighborsOptions &options)
{
    isaac::workflow::NeighborsFinderWorkflow workflow(
        options.seedLength,
        options.referenceGenome,
        options.outputFile,
        options.jobs);
    workflow.run(
        options.maskWidth,
        options.mask,
        options.neighborhoodWidth);
}

int main(int argc, char *argv[])
{
    isaac::common::configureMemoryManagement(true, true);
    isaac::common::run(findNeighbors, argc, argv);
}

