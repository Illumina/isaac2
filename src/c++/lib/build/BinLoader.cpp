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
 ** \file BinLoader.cpp
 **
 ** Reorders aligned data and stores results in bam file.
 ** 
 ** \author Roman Petrovski
 **/

#include <boost/foreach.hpp>
#include <boost/function_output_iterator.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/ptr_container/ptr_vector.hpp>

#include "bam/Bam.hh"

#include "build/BinLoader.hh"
#include "common/Memory.hh"

namespace isaac
{
namespace build
{

void verifyFragmentIntegrity(const io::FragmentAccessor &fragment)
{
/*
    if (fragment.getTotalLength() != fragment.totalLength_)
    {
        ISAAC_THREAD_CERR << pFragment << std::endl;
        BOOST_THROW_EXCEPTION(common::IoException(
            errno, (boost::format("Corrupt fragment (fragment total length is broken at %ld) read from %s") %
                (reinterpret_cast<const char *>(pFragment) - &data_.front()) %
                bin.getPathString()).str()));
    }
    if (io::FragmentAccessor::magicValue_ != fragment.magic_)
    {
        ISAAC_THREAD_CERR << pFragment << std::endl;
        BOOST_THROW_EXCEPTION(common::IoException(
            errno, (boost::format("Corrupt fragment (magic is broken) read from %s") % bin.getPathString()).str()));
    }
*/
}

void BinLoader::loadData(BinData &binData)
{
    ISAAC_THREAD_CERR << "Loading unsorted data" << std::endl;
    const clock_t startLoad = clock();

    if(binData.isUnalignedBin())
    {
        loadUnalignedData(binData);
    }
    else
    {
        loadAlignedData(binData);
    }

    ISAAC_THREAD_CERR << "Loading unsorted data done in " << (clock() - startLoad) / 1000 << "ms" << std::endl;
}

void BinLoader::loadUnalignedData(BinData &binData)
{
    if(binData.bin_.getDataSize())
    {
        ISAAC_THREAD_CERR << "Reading unaligned records from " << binData.bin_ << std::endl;
        std::istream isData(binData.inputFileBuf_.get(binData.bin_.getPath()));
        if (!isData) {
            BOOST_THROW_EXCEPTION(common::IoException(errno, "Failed to open " + binData.bin_.getPathString()));
        }

        if (!isData.seekg(binData.bin_.getDataOffset(), std::ios_base::beg))
        {
            BOOST_THROW_EXCEPTION(common::IoException(
                errno, (boost::format("Failed to seek to position %d in %s") % binData.bin_.getDataOffset() % binData.bin_.getPathString()).str()));
        }
        // TODO: this takes time to fill it up with 0... 2 seconds per bin easily
        binData.data_.resize(binData.bin_);
        if (!isData.read(&binData.data_.front(), binData.bin_.getDataSize())) {
            BOOST_THROW_EXCEPTION(common::IoException(
                errno, (boost::format("Failed to read %d bytes from %s") % binData.bin_.getDataSize() % binData.bin_.getPathString()).str()));
        }

/*
        unsigned count = 0;
        unsigned long offset = 0;
        while(data_.size() != offset)
        {
            const io::FragmentAccessor &fragment = data_.getFragment(offset);
            offset += fragment.getTotalLength();
            ++count;
            ISAAC_THREAD_CERR << "offset " << offset << fragment << std::endl;
            ISAAC_ASSERT_MSG(offset <= data_.size(), "Fragment crosses bin boundary."
                " offset:" << offset << " data_.size():" << data_.size() << " " << bin_);
        }
*/
        ISAAC_THREAD_CERR << "Reading unaligned records done from " << binData.bin_ << std::endl;
    }
    else
    {
//        ISAAC_THREAD_CERR << "No unaligned records to read done from " << bin_ << std::endl;
    }
}


const io::FragmentAccessor &BinLoader::loadFragment(BinData &binData, std::istream &isData, unsigned long &offset)
{
    io::FragmentHeader header;
    if (!isData.read(reinterpret_cast<char*>(&header), sizeof(header))) {
        BOOST_THROW_EXCEPTION(common::IoException(
            errno, (boost::format("Failed to read FragmentHeader bytes from %s") % binData.bin_.getPathString()).str()));
    }

    const unsigned fragmentLength = header.getTotalLength();
    // fragments that don't belong to the bin are supposed to go into chunk 0
    offset = binData.dataDistribution_.addBytes(
        binData.bin_.hasPosition(header.fStrandPosition_) ? header.fStrandPosition_ - binData.bin_.getBinStart() : 0, fragmentLength);
//            ISAAC_THREAD_CERR << "offset:" << offset << " fragment: " << header << std::endl;

    // TODO: this takes time to fill it up with 0... 2 seconds per bin easily
    binData.data_.resize(std::max(binData.data_.size(), offset + fragmentLength));

    io::FragmentAccessor &fragment = binData.data_.getFragment(offset);
    io::FragmentHeader &headerRef = fragment;
    headerRef = header;
    if (!isData.read(reinterpret_cast<char*>(&fragment) + sizeof(header), fragmentLength - sizeof(header))) {
        BOOST_THROW_EXCEPTION(common::IoException(
            errno, (boost::format("Failed to read %d bytes from %s") % fragmentLength % binData.bin_.getPathString()).str()));
    }

//    ISAAC_THREAD_CERR << "LOADED: " << fragment << std::endl;

    return fragment;
}

void BinLoader::storeFragmentIndex(
    const io::FragmentAccessor& fragment,
    unsigned long offset,
    unsigned long mateOffset, BinData& binData)
{
    if (fragment.flags_.reverse_ || fragment.flags_.unmapped_)
    {
        RStrandOrShadowFragmentIndex rsIdx(
            fragment.fStrandPosition_, // shadows are stored at the position of their singletons,
            io::FragmentIndexAnchor(fragment),
            FragmentIndexMate(fragment.flags_.mateUnmapped_,
                              fragment.flags_.mateReverse_,
                              fragment.mateStorageBin_,
                              fragment.mateAnchor_),
            fragment.duplicateClusterRank_);
        rsIdx.dataOffset_ = offset;
        rsIdx.mateDataOffset_ = mateOffset;
        binData.rIdx_.push_back(rsIdx);
    }
    else
    {
        FStrandFragmentIndex fIdx(
            fragment.fStrandPosition_,
            FragmentIndexMate(fragment.flags_.mateUnmapped_,
                              fragment.flags_.mateReverse_,
                              fragment.mateStorageBin_,
                              fragment.mateAnchor_),
            fragment.duplicateClusterRank_);
        fIdx.dataOffset_ = offset;
        fIdx.mateDataOffset_ = mateOffset;
        binData.fIdx_.push_back(fIdx);
    }
}


void BinLoader::loadAlignedData(BinData &binData)
{
    if(binData.bin_.getDataSize())
    {
        ISAAC_THREAD_CERR << "Reading alignment records from " << binData.bin_ << std::endl;
        unsigned long dataSize = 0;
        std::istream isData(binData.inputFileBuf_.get(binData.bin_.getPath()));
        if (!isData) {
            BOOST_THROW_EXCEPTION(common::IoException(errno, "Failed to open " + binData.bin_.getPathString()));
        }
        if (!isData.seekg(binData.bin_.getDataOffset()))
        {
            BOOST_THROW_EXCEPTION(common::IoException(
                errno, (boost::format("Failed to seek to position %d in %s") % binData.bin_.getDataOffset() % binData.bin_.getPathString()).str()));
        }

        binData.rIdx_.clear();
        binData.fIdx_.clear();
        binData.seIdx_.clear();
        unsigned long binOffset = 0;

        while(isData && dataSize != binData.bin_.getDataSize())
        {
            unsigned long offset = 0;
            const io::FragmentAccessor &fragment = loadFragment(binData, isData, offset);
            const unsigned fragmentLength = fragment.getTotalLength();
            dataSize += fragmentLength;

            verifyFragmentIntegrity(fragment);

            binOffset += fragmentLength;
            if (!fragment.flags_.paired_)
            {
                SeFragmentIndex seIdx(fragment.fStrandPosition_);
                seIdx.dataOffset_ = offset;
                binData.seIdx_.push_back(seIdx);
            }
            else
            {
                unsigned long mateOffset = offset;

                // mates are present even if they belong to a different bin
                // if (bin_.coversPosition(fragment.mateFStrandPosition_))
                {
                    const io::FragmentAccessor &mateFragment = loadFragment(binData, isData, mateOffset);
                    ISAAC_ASSERT_MSG(mateFragment.clusterId_ == fragment.clusterId_, "mateFragment.clusterId_ != fragment.clusterId_");
                    ISAAC_ASSERT_MSG(mateFragment.flags_.unmapped_ == fragment.flags_.mateUnmapped_, "mateFragment.flags_.unmapped_ != fragment.flags_.mateUnmapped_");
                    ISAAC_ASSERT_MSG(mateFragment.flags_.reverse_ == fragment.flags_.mateReverse_,
                                     "mateFragment.flags_.reverse_ != fragment.flags_.mateReverse_" << fragment << " " << mateFragment);

                    const unsigned mateFragmentLength = mateFragment.getTotalLength();
                    dataSize += mateFragmentLength;

                    verifyFragmentIntegrity(mateFragment);

                    storeFragmentIndex(mateFragment, mateOffset, offset, binData);

                    binOffset += mateFragmentLength;
                }

                storeFragmentIndex(fragment, offset, mateOffset, binData);
            }
        }
        ISAAC_THREAD_CERR << "Reading alignment records done from " << binData.bin_ << std::endl;

        binData.finalize();
    }
}

} // namespace build
} // namespace isaac
