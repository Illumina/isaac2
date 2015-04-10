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
 ** \file BamLoader.cpp
 **
 ** Component to read Bam files.
 **
 ** \author Roman Petrovski
 **/

#include <boost/format.hpp>
#include <boost/tuple/tuple.hpp>

#include "common/Debug.hh"
#include "io/BamLoader.hh"
#include "oligo/Nucleotides.hh"

namespace isaac
{
namespace io
{

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

BamLoader::BamLoader(
    std::size_t maxPathLength,
    common::ThreadVector &threads,
    const unsigned coresMax) :
    // one of the threads is busy parsing
    bgzfReader_(std::max(1U, coresMax - 1), BUFFER_SIZE / std::max(1U, coresMax - 1) / bgzf::BgzfReader::UNCOMPRESSED_BGZF_BLOCK_SIZE / 10),
    lastUnparsedBytes_(0),
    // since the loading from bam has to happen sequentially, and we can have only one thread parsing,
    // the only parallelization can occur between a thread parsing and a bunch of threads loading and decompressing
    // so, 2 is a good ceiling here
    decompressParseParallelization_(std::min<unsigned>(coresMax, 2)),
    unparsedBytes_(decompressParseParallelization_),
    decompressionBuffers_(decompressParseParallelization_),
    decompressParseParallelizationThreads_(threads),
    nextDecompressorThread_(0),
    nextParserThread_(0)
{
    reserveBuffers(maxPathLength);
}

void BamLoader::reserveBuffers(std::size_t maxPathLength)
{
    lastPassBam_.reserve(BUFFER_SIZE);
    std::for_each(decompressionBuffers_.begin(), decompressionBuffers_.end(), boost::bind(&std::vector<char>::reserve, _1, BUFFER_SIZE));
}



} // namespace io
} // namespace isaac
