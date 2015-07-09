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
 ** \file BpbToWigWorkflow.cpp
 **
 ** \brief see BpbToWigWorkflow.hh
 **
 ** \author Roman Petrovski
 **/

#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>

#include "common/Debug.hh"
#include "common/Exceptions.hh"
#include "common/FastIo.hh"
#include "common/FileSystem.hh"
#include "io/BitsetSaver.hh"
#include "workflow/BpbToWigWorkflow.hh"

namespace isaac
{
namespace workflow
{

BpbToWigWorkflow::BpbToWigWorkflow(
    const unsigned wigDefaultValue,
    const bfs::path &sortedReferenceMetadata,
    const bfs::path &inputFilePath,
    const std::string &outputFormatString,
    const unsigned bitsPerValue,
    const bool knownSitesOnly,
    const bool bedPrintAllPositions
    )
    : wigDefaultValue_(wigDefaultValue),
      sortedReferenceMetadata_(sortedReferenceMetadata),
      inputFilePath_(inputFilePath),
      outputFormatString_(outputFormatString),
      bitsPerValue_(bitsPerValue),
      knownSitesOnly_(knownSitesOnly),
      bedPrintAllPositions_(bedPrintAllPositions),
      xml_(reference::loadSortedReferenceXml(sortedReferenceMetadata_)),
      threads_(knownSitesOnly_ ? boost::thread::hardware_concurrency() : 1),
      contigs_(knownSitesOnly_ ? reference::loadContigs(xml_.getContigs(), threads_) : reference::ContigList())
{

}

void BpbToWigWorkflow::run()
{
    boost::iostreams::filtering_istream is;
    if (common::isDotGzPath(inputFilePath_))
    {
        is.push(boost::iostreams::gzip_decompressor());
    }
    is.push(boost::iostreams::file_source(inputFilePath_.string()));

    if (!is)
    {
        const boost::format message = boost::format("Failed to open bitset file %s for reading: %s") %
            inputFilePath_ % strerror(errno);
        BOOST_THROW_EXCEPTION(common::IoException(errno, message.str()));
    }

    const reference::SortedReferenceMetadata::Contigs &contigs = xml_.getKaryotypeOrderedContigs();

    if ("wig" == outputFormatString_)
    {
        if (bitsPerValue_ / 8 > sizeof(unsigned char))
        {
            printWig<unsigned short>(is, bitsPerValue_, contigs);
        }
        else
        {
            printWig<unsigned char>(is, bitsPerValue_, contigs);
        }
    }
    else if ("bed" == outputFormatString_)
    {
        if (bitsPerValue_ / 8 > sizeof(unsigned char))
        {
            printBed<unsigned short>(is, bitsPerValue_, contigs);
        }
        else
        {
            printBed<unsigned char>(is, bitsPerValue_, contigs);
        }
    }
    else
    {
        ISAAC_ASSERT_MSG(false, "Unknown output format: " << outputFormatString_);
    }
}

template <typename ReadType>
void BpbToWigWorkflow::printWig(
    std::istream &bitsetFile,
    const unsigned bitsPerValue,
    const reference::SortedReferenceMetadata::Contigs &contigs )
{
    const ReadType valueMask = (sizeof(ReadType) * 8 == bitsPerValue) ? ~ReadType(0) :  ~(~ReadType(0) << (bitsPerValue));
    const std::size_t bytesPerValue = (bitsPerValue + 7) / 8;
    ISAAC_ASSERT_MSG(bytesPerValue <= sizeof(ReadType), "Single value requires more bytes than supported:" << sizeof(ReadType));
    ISAAC_ASSERT_MSG(!((sizeof(ReadType) * 8) % bitsPerValue), "bits must fully fit in " << sizeof(ReadType) << " bytes");
    ReadType bits = 0;

    ISAAC_THREAD_CERR << "bytesPerValue " << bytesPerValue << std::endl;
    ISAAC_THREAD_CERR << "bitsPerValue " << bitsPerValue << std::endl;
    ISAAC_THREAD_CERR << "valueMask " << std::hex << unsigned(valueMask) << std::endl;

    common::StaticVector<char, 1024> outBuf;

    bitsetFile.read(reinterpret_cast<char*>(&bits), bytesPerValue);
    std::size_t bitPos = bitsPerValue;
    BOOST_FOREACH(const reference::SortedReferenceMetadata::Contig &contig, contigs)
    {
        ISAAC_THREAD_CERR << contig << std::endl;
        std::cout << "variableStep chrom=" << contig.name_ << "\n";
        std::size_t contigBases = contig.totalBases_;
        while(bitsetFile && contigBases)
        {
            const std::size_t contigPosition = contig.totalBases_ - contigBases + 1;
            const unsigned value = (bits & valueMask);
            if (wigDefaultValue_ != value &&
                (!knownSitesOnly_ || oligo::REFERENCE_OLIGO_N != contigs_.at(contig.index_).forward_.at(contigPosition-1)))
            {
                outBuf.clear();
                common::appendUnsignedInteger(outBuf, contigPosition);
                outBuf.push_back('\t');
                common::appendUnsignedInteger(outBuf, value);
                outBuf.push_back('\r');
                outBuf.push_back('\n');
                std::cout.write(&outBuf.front(), outBuf.size());
            }
            if (!(bitPos % (sizeof(ReadType) * 8)))
            {
                bitsetFile.read(reinterpret_cast<char*>(&bits), bytesPerValue);
            }
            else
            {
                bits >>= bitsPerValue;
            }
            --contigBases;
            bitPos += bitsPerValue;
        }
        if (!std::cout)
        {
            const boost::format message = boost::format("Failed to write to cout") % strerror(errno);
            BOOST_THROW_EXCEPTION(common::IoException(errno, message.str()));
        }
        if (!bitsetFile && !bitsetFile.eof())
        {
            const boost::format message = boost::format("Failed to read bitset file %s: %s") % inputFilePath_ % strerror(errno);
            BOOST_THROW_EXCEPTION(common::IoException(errno, message.str()));
        }
    }
}

template <typename ReadType>
void BpbToWigWorkflow::printBed(
    std::istream &bitsetFile,
    const unsigned bitsPerValue,
    const reference::SortedReferenceMetadata::Contigs &contigs )
{
    std::cout << "track graphType=bar type=bedGraph" << std::endl;
/*
    "The count operand can be an immediate value or register CL. The count is masked to five bits, which limits the count range to 0 to 31."
    See http://www.intel.com/design/intarch/manuals/243191.htm
*/
    const ReadType valueMask = (sizeof(ReadType) * 8 == bitsPerValue) ? ~ReadType(0) :  ~(~ReadType(0) << (bitsPerValue));

    const std::size_t bytesPerValue = (bitsPerValue + 7) / 8;
    ISAAC_ASSERT_MSG(bytesPerValue <= sizeof(ReadType), "Single value requires more bytes than supported:" << sizeof(ReadType));
    ISAAC_ASSERT_MSG(!((sizeof(ReadType) * 8) % bitsPerValue), "bits must fully fit in " << sizeof(ReadType) << " bytes");
    ReadType bits = 0;

    ISAAC_THREAD_CERR << "bytesPerValue " << bytesPerValue << std::endl;
    ISAAC_THREAD_CERR << "bitsPerValue " << bitsPerValue << std::endl;
    ISAAC_THREAD_CERR << "valueMask " << std::hex << unsigned(valueMask) << std::endl;

    bitsetFile.read(reinterpret_cast<char*>(&bits), bytesPerValue);
    ReadType lastValue = 0;
    std::size_t lastContigPosition = 0;
    std::size_t bitPos = bitsPerValue;
    common::StaticVector<char, 1024> outBuf;
    BOOST_FOREACH(const reference::SortedReferenceMetadata::Contig &contig, contigs)
    {
        ISAAC_THREAD_CERR << contig << std::endl;
        std::size_t contigBases = contig.totalBases_;
        std::size_t acgtBases = 0;
        while(bitsetFile && contigBases)
        {
            const std::size_t contigPosition = contig.totalBases_ - contigBases;
            const ReadType value = (bits & valueMask);
            if (bedPrintAllPositions_ || value != lastValue)
            {
                if (!knownSitesOnly_ || oligo::REFERENCE_OLIGO_N != contigs_.at(contig.index_).forward_.at(contigPosition))
                {
                    outBuf.clear();
                    std::copy(contig.name_.begin(), contig.name_.end(), std::back_inserter(outBuf));
                    outBuf.push_back('\t');
                    common::appendUnsignedInteger(outBuf, lastContigPosition);
                    outBuf.push_back('\t');
                    common::appendUnsignedInteger(outBuf, contigPosition);
                    outBuf.push_back('\t');
                    common::appendUnsignedInteger(outBuf, lastValue);
                    outBuf.push_back('\r');
                    outBuf.push_back('\n');
                    std::cout.write(&outBuf.front(), outBuf.size());

//                    std::cout << contig.name_ << "\t" << lastContigPosition << "\t" << contigPosition << "\t" << unsigned(lastValue) << std::endl;
                }
                lastContigPosition = contigPosition;
                lastValue = value;
            }
            acgtBases += !knownSitesOnly_ || oligo::REFERENCE_OLIGO_N != contigs_.at(contig.index_).forward_.at(contigPosition);
            if (!(bitPos % (sizeof(ReadType) * 8)))
            {
                bitsetFile.read(reinterpret_cast<char*>(&bits), bytesPerValue);
            }
            else
            {
                bits >>= bitsPerValue;
            }
            --contigBases;
            bitPos += bitsPerValue;
        }
        ISAAC_ASSERT_MSG(!knownSitesOnly_ || contig.acgtBases_ == acgtBases, "Printed mismatching count of ACGT positions: " << acgtBases << " " << contig);
        if (!bitsetFile && !bitsetFile.eof())
        {
            const boost::format message = boost::format("Failed to read bitset file %s: %s") % inputFilePath_ % strerror(errno);
            BOOST_THROW_EXCEPTION(common::IoException(errno, message.str()));
        }
        const std::size_t contigPosition = contig.totalBases_ - contigBases;
        if (!bedPrintAllPositions_)
        {
            std::cout << contig.name_ << "\t" << lastContigPosition << "\t" << contigPosition << "\t" << unsigned(lastValue) << std::endl;
        }
        lastContigPosition = 0;
    }
}


} // namespace workflow
} // namespace isaac
