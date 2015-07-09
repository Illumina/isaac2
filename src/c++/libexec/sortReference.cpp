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
 ** \file sortReference.cpp
 **
 ** The main for the reference sorter.
 **
 ** \author Come Raczy
 **/

#include "oligo/Kmer.hh"
#include "options/SortReferenceOptions.hh"
#include "reference/ReferenceSorter.hh"

template <typename KmerT>
void sortReferenceT(const isaac::options::SortReferenceOptions &options)
{
    isaac::reference::ReferenceSorter<KmerT> referenceSorter(
        options.maskWidth,
        options.mask,
        options.contigsXml,
        options.genomeNeighborsFile,
        options.outFile,
        options.repeatThreshold);
    referenceSorter.run();
}


template <unsigned seedLength>
void sortReference(const isaac::options::SortReferenceOptions &options)
{
    if (seedLength == options.seedLength)
    {
        sortReferenceT<isaac::oligo::BasicKmerType<seedLength> >(options);
    }
    else
    {
        sortReference<seedLength - 4>(options);
    }
}

template <>
void sortReference<4>(const isaac::options::SortReferenceOptions &options)
{
    ISAAC_ASSERT_MSG(false, "Unexpected seedLength requested: " << options.seedLength);
}

int main(int argc, char *argv[])
{
    isaac::common::run(&sortReference<80>, argc, argv);
}
