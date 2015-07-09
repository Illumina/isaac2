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
 ** \file Cluster.hh
 **
 ** Component containing the data associated to a cluster: sequence and quality
 ** strings for all the reads in the cluster.
 **
 ** \author Come Raczy
 **/

#include <numeric>
#include <boost/foreach.hpp>

#include "alignment/Cluster.hh"
#include "flowcell/ReadMetadata.hh"
#include "oligo/Nucleotides.hh"

namespace isaac
{
namespace alignment
{

static std::vector<char>::const_iterator uninitialized()
{
    static std::vector<char> blah;
    return blah.begin();
}

Cluster::Cluster(const unsigned maxReadLength)
    : tile_(0)
    , id_(0)
    , pf_(false)
    , barcodeLength_(0)
    , nonEmptyReads_(0)
    , bclData_(uninitialized())
    , readNameBegin_(uninitialized())
    , readNameEnd_(uninitialized())
{
    push_back(Read(maxReadLength, 0));
    push_back(Read(maxReadLength, 1));
}

void Cluster::init(
    const flowcell::ReadMetadataList &readMetadataList,
    std::vector<char>::const_iterator bclData,
    const unsigned tile,
    const unsigned long id,
    const ClusterXy &xy,
    const bool pf,
    const unsigned barcodeLength,
    const unsigned nameLengthMax)
{
    tile_ = tile;
    id_ = id;
    xy_ = xy;
    pf_ = pf;
    barcodeLength_ = barcodeLength;
    nonEmptyReads_ = 0;
    bclData_ = bclData;
    bclData += barcodeLength_;
    BOOST_FOREACH(const flowcell::ReadMetadata &readMetadata, readMetadataList)
    {
        //TODO: this check is to potentially allow 0-length reads (such as masked out) in readMetadataList
        // however, at the moment masked-out reads don't make into readMetadatList.
        if (readMetadata.getLength())
        {
            at(readMetadata.getIndex()).decodeBcl(
                bclData,
                bclData + readMetadata.getLength(),
                readMetadata.getIndex());
            bclData += readMetadata.getLength();
            ++nonEmptyReads_;
        }
    }
    readNameBegin_ = bclData;
    readNameEnd_ = std::find(readNameBegin_, readNameBegin_ + nameLengthMax, 0);
}

std::vector<char>::const_iterator Cluster::getBclData(const unsigned readIndex) const
{
    const unsigned offset =
        std::accumulate(begin(), begin() + readIndex, 0,
                        bind(std::plus<unsigned>(), _1, boost::bind(&Read::getLength, _2)));
    return bclData_ + barcodeLength_ + offset;
}

unsigned long Cluster::getBarcodeSequence() const
{
    return oligo::packBclBases(bclData_, bclData_ + barcodeLength_);
}

} // namespace alignment
} // namespace isaac
