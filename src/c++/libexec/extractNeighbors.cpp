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
#include "options/ExtractNeighborsOptions.hh"
#include "workflow/ExtractNeighborsWorkflow.hh"

void extractNeighbors(const isaac::options::ExtractNeighborsOptions &options)
{
    isaac::workflow::ExtractNeighborsWorkflow workflow(
        options.sortedReferenceMetadata_,
        options.outputFilePath_,
        options.highRepeatsFilePath_
        );

    if (16 == options.seedLength)
    {
        workflow.run<isaac::oligo::ShortKmerType>();
    }
    else if (28 == options.seedLength)
    {
        workflow.run<isaac::oligo::BasicKmerType<28> >();
    }
    else if (30 == options.seedLength)
    {
        workflow.run<isaac::oligo::BasicKmerType<30> >();
    }
    else if (32 == options.seedLength)
    {
        workflow.run<isaac::oligo::KmerType>();
    }
    else if (34 == options.seedLength)
    {
        workflow.run<isaac::oligo::BasicKmerType<34> >();
    }
    else if (36 == options.seedLength)
    {
        workflow.run<isaac::oligo::BasicKmerType<36> >();
    }
    else if (64 == options.seedLength)
    {
        workflow.run<isaac::oligo::LongKmerType>();
    }
    else
    {
        ISAAC_ASSERT_MSG(false, "Unexpected seedLength " << options.seedLength)
    }
}

int main(int argc, char *argv[])
{
    isaac::common::run(extractNeighbors, argc, argv);
}

