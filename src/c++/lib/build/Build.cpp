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
 ** \file Build.cpp
 **
 ** Reorders aligned data and stores results in bam file.
 **
 ** \author Roman Petrovski
 **/

#include "common/config.h"

#ifdef HAVE_NUMA
#include <numa.h>
#include <numaif.h>
#endif //HAVE_NUMA
 
#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/regex.hpp>

#include "bam/Bam.hh"
#include "bam/BamIndexer.hh"
#include "bgzf/BgzfCompressor.hh"
#include "build/Build.hh"
#include "common/Debug.hh"
#include "common/FileSystem.hh"
#include "common/Threads.hpp"
#include "io/Fragment.hh"
#include "reference/ContigLoader.hh"

#include "BuildStatsXml.hh"
#include "SortedReferenceXmlBamHeaderAdapter.hh"

namespace isaac
{
namespace build
{

const unsigned BuildContigMap::UNMAPPED_CONTIG;
/**
 * \return Returns the total memory in bytes required to load the bin data and indexes
 */
static unsigned long getBinTotalSize(const alignment::BinMetadata & binMetadata)
{
    return
        binMetadata.getDataSize() +
        binMetadata.getFIdxElements() * sizeof(FStrandFragmentIndex) +
        binMetadata.getRIdxElements() * sizeof(RStrandOrShadowFragmentIndex) +
        binMetadata.getSeIdxElements() * sizeof(SeFragmentIndex);
}

unsigned long Build::estimateBinCompressedDataRequirements(
    const alignment::BinMetadata & binMetadata,
    const unsigned outputFileIndex) const
{
    // TODO: put the real number in here.
    static const unsigned long EMPTY_BGZF_BLOCK_SIZE = 1234UL;
    if (!binMetadata.getTotalElements())
    {
        return EMPTY_BGZF_BLOCK_SIZE;
    }

    unsigned long thisOutputFileBarcodeElements = 0;
    // accumulate size required to store all barcodes that map to the same output file
    BOOST_FOREACH(const flowcell::BarcodeMetadata& barcode, barcodeMetadataList_)
    {
        const unsigned barcodeOutputFileIndex = barcodeBamMapping_.getSampleIndex(barcode.getIndex());
        if (outputFileIndex == barcodeOutputFileIndex)
        {
            const unsigned barcodeIndex = barcode.getIndex();
            ISAAC_ASSERT_MSG(0 != binMetadata.getTotalElements() || 0 == binMetadata.getBarcodeElements(barcodeIndex), "Can't have empty bin with non-empty bin barcode");

            thisOutputFileBarcodeElements += binMetadata.getBarcodeElements(barcodeIndex);
        }
    }

    // assume all data will take the same fraction or less than the number derived from demultiplexed fragments.
    return EMPTY_BGZF_BLOCK_SIZE +
        ((getBinTotalSize(binMetadata) * thisOutputFileBarcodeElements +
            binMetadata.getTotalElements() - 1) / binMetadata.getTotalElements()) * expectedBgzfCompressionRatio_;
}

inline boost::filesystem::path getSampleBamPath(
    const boost::filesystem::path &outputDirectory,
    const flowcell::BarcodeMetadata &barcode)
{
    return outputDirectory / barcode.getProject() / barcode.getSampleName() / "sorted.bam";
}
/**
 * \brief Produces mapping so that all barcodes having the same sample name go into the same output file.
 *
 * \return Returns the pair of a vector of that maps a barcode index to a unique output file index in the
 *         second vector so that two barcodes that are supposed to go into the same file will end up
 *         having the same mapping
 */
BarcodeBamMapping mapBarcodesToFiles(
    const boost::filesystem::path &outputDirectory,
    const flowcell::BarcodeMetadataList &barcodeMetadataList)
{
    // Map barcodes to projects
    std::vector<std::string> projects;
    std::transform(barcodeMetadataList.begin(), barcodeMetadataList.end(), std::back_inserter(projects),
                   boost::bind(&flowcell::BarcodeMetadata::getProject, _1));
    std::sort(projects.begin(), projects.end());
    projects.erase(std::unique(projects.begin(), projects.end()), projects.end());

    BarcodeBamMapping::BarcodeProjectIndexMap barcodeProject(barcodeMetadataList.size());
    BOOST_FOREACH(const flowcell::BarcodeMetadata &barcode, barcodeMetadataList)
    {
        barcodeProject.at(barcode.getIndex()) =
            std::distance(projects.begin(), std::lower_bound(projects.begin(), projects.end(), barcode.getProject()));
    }

    // Map barcodes to sample paths
    std::vector<boost::filesystem::path> samples;
    std::transform(barcodeMetadataList.begin(), barcodeMetadataList.end(), std::back_inserter(samples),
                   boost::bind(&getSampleBamPath, outputDirectory, _1));
    std::sort(samples.begin(), samples.end());
    samples.erase(std::unique(samples.begin(), samples.end()), samples.end());

    BarcodeBamMapping::BarcodeProjectIndexMap barcodeSample(barcodeMetadataList.size());
    BOOST_FOREACH(const flowcell::BarcodeMetadata &barcode, barcodeMetadataList)
    {
        barcodeSample.at(barcode.getIndex()) =
            std::distance(samples.begin(), std::lower_bound(samples.begin(), samples.end(),
                                                            getSampleBamPath(outputDirectory, barcode)));
    }

    return BarcodeBamMapping(barcodeProject, barcodeSample, samples);
}

inline bool orderBySampleIndex(
    const BarcodeBamMapping &barcodeBamMapping,
    const flowcell::BarcodeMetadata &left,
    const flowcell::BarcodeMetadata &right)
{
    return barcodeBamMapping.getSampleIndex(left.getIndex()) < barcodeBamMapping.getSampleIndex(right.getIndex());
}

std::vector<boost::shared_ptr<boost::iostreams::filtering_ostream> > Build::createOutputFileStreams(
    const flowcell::TileMetadataList &tileMetadataList,
    const flowcell::BarcodeMetadataList &barcodeMetadataList,
    boost::ptr_vector<bam::BamIndex> &bamIndexes) const
{
    unsigned sinkIndexToCreate = 0;
    std::vector<boost::shared_ptr<boost::iostreams::filtering_ostream> > ret;
    ret.reserve(barcodeBamMapping_.getTotalSamples());

    std::vector<boost::filesystem::path> directories;
    BOOST_FOREACH(const flowcell::BarcodeMetadata &barcode, barcodeMetadataList)
    {
        const boost::filesystem::path &bamPath = barcodeBamMapping_.getFilePath(barcode);
        directories.push_back(bamPath.parent_path());
        directories.push_back(directories.back().parent_path());
    }
    common::createDirectories(directories);

    flowcell::BarcodeMetadataList barcodesOrderedBySample(barcodeMetadataList);
    std::sort(barcodesOrderedBySample.begin(), barcodesOrderedBySample.end(),
              boost::bind(&orderBySampleIndex, boost::ref(barcodeBamMapping_), _1, _2));
    BOOST_FOREACH(const flowcell::BarcodeMetadata &barcode, barcodesOrderedBySample)
    {
        if (sinkIndexToCreate == barcodeBamMapping_.getSampleIndex(barcode.getIndex()))
        {
            const boost::filesystem::path &bamPath = barcodeBamMapping_.getFilePath(barcode);
            if (!barcode.isUnmappedReference())
            {
                ISAAC_THREAD_CERR << "Created BAM file: " << bamPath << std::endl;

                const reference::SortedReferenceMetadata &sampleReference =
                    sortedReferenceMetadataList_.at(barcode.getReferenceIndex());

                std::string compressedHeader;
                {
                    std::ostringstream oss(compressedHeader);
                    boost::iostreams::filtering_ostream bgzfStream;
                    bgzfStream.push(bgzf::BgzfCompressor(bamGzipLevel_));
                    bgzfStream.push(oss);
                    bam::serializeHeader(bgzfStream,
                                         argv_,
                                         description_,
                                         bamHeaderTags_,
                                         bamPuFormat_,
                                         makeSortedReferenceXmlBamHeaderAdapter(
                                             sampleReference,
                                             boost::bind(&BuildContigMap::isMapped, &contigMap_, barcode.getReferenceIndex(), _1),
                                             tileMetadataList, barcodeMetadataList,
                                             barcode.getSampleName()));
                    bgzfStream.strict_sync();
                    compressedHeader = oss.str();
                }

                ret.push_back(boost::shared_ptr<boost::iostreams::filtering_ostream>(new boost::iostreams::filtering_ostream()));
                boost::iostreams::filtering_ostream &bamStream = *ret.back();
                bamStream.push(io::FileSinkWithMd5(bamPath.c_str(), std::ios_base::binary));
                if (!bamStream) {
                    BOOST_THROW_EXCEPTION(common::IoException(errno, "Failed to open output BAM file " + bamPath.string()));
                }

                if (!bamStream.write(compressedHeader.c_str(), compressedHeader.size()))
                {
                    BOOST_THROW_EXCEPTION(
                        common::IoException(errno, (boost::format("Failed to write %d bytes into stream %s") %
                            compressedHeader.size() % bamPath.string()).str()));
                }

                // Create BAM Indexer
                unsigned headerCompressedLength = compressedHeader.size();
                unsigned contigCount = sampleReference.getContigsCount(
                    boost::bind(&BuildContigMap::isMapped, &contigMap_, barcode.getReferenceIndex(), _1));
                bamIndexes.push_back(new bam::BamIndex(bamPath, contigCount, headerCompressedLength));
            }
            else
            {
                ret.push_back(boost::shared_ptr<boost::iostreams::filtering_ostream>());
                bamIndexes.push_back(new bam::BamIndex());
                ISAAC_THREAD_CERR << "Skipped BAM file due to unmapped barcode reference: " << bamPath << " " << barcode << std::endl;
            }
            ++sinkIndexToCreate;
        }
    }
    ISAAC_ASSERT_MSG(barcodeBamMapping_.getTotalSamples() == sinkIndexToCreate, "must create all output file sinks");

    return ret;
}

const alignment::BinMetadataCRefList filterBins(
    const alignment::BinMetadataList& bins,
    const std::string &binRegexString)
{
    alignment::BinMetadataCRefList ret;

    if ("all" == binRegexString)
    {
        std::transform(bins.begin(), bins.end(), std::back_inserter(ret),
                       [](const alignment::BinMetadata &bm){return boost::cref(bm);});
    }
    else if ("skip-empty" == binRegexString)
    {
        BOOST_FOREACH(const alignment::BinMetadata &binMetadata, bins)
        {
            if (!binMetadata.isEmpty())
            {
                ret.push_back(boost::cref(binMetadata));
            }
        }
    }
    else // use regex to filter bins by name
    {
        std::string regexString(binRegexString);
        std::replace(regexString.begin(), regexString.end(), ',', '|');
        boost::regex re(regexString);
        BOOST_FOREACH(const alignment::BinMetadata &bin, bins)
        {
            if (!bin.isEmpty() && boost::regex_search(bin.getPath().filename().string(), re))
            {
                ret.push_back(boost::ref(bin));
            }
        }
        if (ret.empty())
        {
            ISAAC_THREAD_CERR << "WARNING: Bam files will be empty. No bins are left after applying the following regex filter: "
                << regexString << std::endl;
        }
    }
    return ret;
}

static void breakUpBin(
    const alignment::BinMetadata& bin,
    const unsigned long partsCount,
    alignment::BinMetadataList &ret)
{
    ISAAC_ASSERT_MSG(0 == bin.getIndex(), "At the moment only bin 0 is expected to be unaligned");
    const unsigned long newBinSize = bin.getDataSize() / partsCount;
    ISAAC_THREAD_CERR << "Breaking unaligned bin of " << bin.getDataSize()/1024/1024 << " megabytes into " <<
        partsCount << " bins of " << newBinSize /1024/1024 << " megabytes for better parallelization: " << bin << std::endl;

    for (unsigned long offset = 0; bin.getDataSize() > offset;)
    {
        alignment::BinMetadata part = bin.getChunks(offset, newBinSize);
        ISAAC_THREAD_CERR << " offset:" << offset << " " << part <<std::endl;
        offset += part.getDataSize();
        ret.push_back(part);
    }
}

/**
 * \brief breaks unaligned bin into about partsCount of roughly equivalent size bins
 */
static alignment::BinMetadataCRefList breakUpUnalignedBin(
    alignment::BinMetadataCRefList bins,
    const unsigned long partsCount,
    const bool keepUnaligned,
    const bool putUnalignedInTheBack,
    alignment::BinMetadataList &unalignedBinParts)
{
    if (!bins.empty())
    {
        // unaligned bins must occur at the start of the list
        const alignment::BinMetadata &bin = bins.front();
        if (bin.isUnalignedBin())
        {
            if (keepUnaligned)
            {
                breakUpBin(bin, partsCount, unalignedBinParts);
            }
            bins.erase(bins.begin());
        }

        if (putUnalignedInTheBack)
        {
            std::transform(unalignedBinParts.begin(), unalignedBinParts.end(), std::back_inserter(bins),
                           [](const alignment::BinMetadata &bm){return boost::cref(bm);});
        }
        else
        {
            std::transform(unalignedBinParts.begin(), unalignedBinParts.end(), std::inserter(bins, bins.begin()),
                           [](const alignment::BinMetadata &bm){return boost::cref(bm);});
        }
    }
    return bins;
}

Build::Build(const std::vector<std::string> &argv,
             const std::string &description,
             const flowcell::FlowcellLayoutList &flowcellLayoutList,
             const flowcell::TileMetadataList &tileMetadataList,
             const flowcell::BarcodeMetadataList &barcodeMetadataList,
             const alignment::BinMetadataList &bins,
             const reference::ReferenceMetadataList &referenceMetadataList,
             const std::vector<alignment::TemplateLengthStatistics> &barcodeTemplateLengthStatistics,
             const reference::SortedReferenceMetadataList &sortedReferenceMetadataList,
             const boost::filesystem::path outputDirectory,
             const unsigned maxLoaders,
             const unsigned maxComputers,
             const unsigned maxSavers,
             const build::GapRealignerMode realignGaps,
             const int bamGzipLevel,
             const std::string &bamPuFormat,
             const std::vector<std::string> &bamHeaderTags,
             const double expectedBgzfCompressionRatio,
             const bool singleLibrarySamples,
             const bool keepDuplicates,
             const bool markDuplicates,
             const bool anchorMate,
             const bool realignGapsVigorously,
             const bool realignDodgyFragments,
             const unsigned realignedGapsPerFragment,
             const bool clipSemialigned,
             const bool loadAllContigs,
             const std::string &binRegexString,
             const unsigned char forcedDodgyAlignmentScore,
             const bool keepUnaligned,
             const bool putUnalignedInTheBack,
             const IncludeTags includeTags,
             const bool pessimisticMapQ,
             const unsigned splitGapLength)
    :argv_(argv),
     description_(description),
     flowcellLayoutList_(flowcellLayoutList),
     tileMetadataList_(tileMetadataList),
     barcodeMetadataList_(barcodeMetadataList),
     unalignedBinParts_(),
     bins_(breakUpUnalignedBin(
         filterBins(bins, binRegexString), maxComputers, keepUnaligned, putUnalignedInTheBack, unalignedBinParts_)),
     barcodeTemplateLengthStatistics_(barcodeTemplateLengthStatistics),
     sortedReferenceMetadataList_(sortedReferenceMetadataList),
     contigMap_(barcodeMetadataList_, bins_, sortedReferenceMetadataList_, !loadAllContigs && "skip-empty" == binRegexString),
     outputDirectory_(outputDirectory),
     maxLoaders_(maxLoaders),
     maxComputers_(maxComputers),
     maxSavers_(maxSavers),
     bamGzipLevel_(bamGzipLevel),
     bamPuFormat_(bamPuFormat),
     bamHeaderTags_(bamHeaderTags),
     forcedDodgyAlignmentScore_(forcedDodgyAlignmentScore),
     singleLibrarySamples_(singleLibrarySamples),
     keepDuplicates_(keepDuplicates),
     markDuplicates_(markDuplicates),
     anchorMate_(anchorMate),
     realignGapsVigorously_(realignGapsVigorously),
     realignDodgyFragments_(realignDodgyFragments),
     realignedGapsPerFragment_(realignedGapsPerFragment),
     clipSemialigned_(clipSemialigned),
     realignGaps_(realignGaps),
     expectedBgzfCompressionRatio_(expectedBgzfCompressionRatio),
     maxReadLength_(getMaxReadLength(flowcellLayoutList_)),
     includeTags_(includeTags),
     pessimisticMapQ_(pessimisticMapQ),
     splitGapLength_(splitGapLength),
     forceTermination_(false),
     threads_(maxComputers_ + maxLoaders_ + maxSavers_),
     contigLists_(reference::loadContigs(sortedReferenceMetadataList_, contigMap_, threads_)),
     barcodeBamMapping_(mapBarcodesToFiles(outputDirectory_, barcodeMetadataList_)),
     bamIndexes_(),
     bamFileStreams_(createOutputFileStreams(tileMetadataList_, barcodeMetadataList_, bamIndexes_)),
     stats_(bins_, barcodeMetadataList_),
     threadBgzfBuffers_(threads_.size(), std::vector<std::vector<char> >(bamFileStreams_.size())),
     threadBgzfStreams_(threads_.size()),
     threadBamIndexParts_(threads_.size()),
     kUniquenessAnnotations_(reference::loadAnnotations(sortedReferenceMetadataList_)),
     gapRealigner_(maxComputers_,
         realignGapsVigorously, realignDodgyFragments, realignedGapsPerFragment, clipSemialigned,
         barcodeMetadataList, barcodeTemplateLengthStatistics, contigLists_),
     binSorter_(singleLibrarySamples_, keepDuplicates_, markDuplicates_, anchorMate_,
               barcodeBamMapping_, barcodeMetadataList_, contigLists_, splitGapLength_, kUniquenessAnnotations_)
{
    computeSlotWaitingBins_.reserve(threads_.size());
    while(threadBgzfStreams_.size() < threads_.size())
    {
        threadBgzfStreams_.push_back(new boost::ptr_vector<boost::iostreams::filtering_ostream>(bamFileStreams_.size()));
    }
    while(threadBamIndexParts_.size() < threads_.size())
    {
        threadBamIndexParts_.push_back(new boost::ptr_vector<bam::BamIndexPart>(bamFileStreams_.size()));
    }
    threads_.execute(boost::bind(&Build::allocateThreadData, this, _1));

    testBinsFitInRam();
}

/**
 * // test if all the bins will fit in remaining RAM
 */
void Build::testBinsFitInRam()
{
    ISAAC_THREAD_CERR << "Making sure all bins fit in memory" << std::endl;
    std::vector<std::vector<char> > testBgzfBuffers(bamFileStreams_.size());
    // no need to really block anything at the moment. just supply fake block to reserveBuffers so it compiles
    common::ScopedMallocBlock fakeMallocBlock(common::ScopedMallocBlock::Off);
    for(alignment::BinMetadataCRefList::const_iterator binIterator = bins_.begin(); bins_.end() != binIterator; ++binIterator)
    {
        {
            boost::unique_lock<boost::mutex> lock(stateMutex_);
            boost::shared_ptr<BinData> binDataPtr;
            const unsigned long bytesFailedToAllocate = reserveBuffers(lock, binIterator, bins_.end(), fakeMallocBlock, 0, binDataPtr);
            if (bytesFailedToAllocate)
            {
                BOOST_THROW_EXCEPTION(
                    common::MemoryException((
                        boost::format("%s requires %d bytes of memory for BAM generation.") %
                            *binIterator % bytesFailedToAllocate).str()));
            }
            BOOST_FOREACH(std::vector<char> &bgzfBuffer, testBgzfBuffers)
            {
                // release memory as next bin test might fail if we don't
                std::vector<char>().swap(bgzfBuffer);
            }
            binDataPtr.reset();
            threadBgzfStreams_.at(0).clear();
            threadBamIndexParts_.at(0).clear();
        }
    }
    ISAAC_THREAD_CERR << "Making sure all bins fit in memory done" << std::endl;
}

void Build::allocateThreadData(const std::size_t threadNumber)
{
#ifdef HAVE_NUMA
    if (-1 == numa_available())
    {
        ISAAC_THREAD_CERR << "WARNING: numa library is unavailable while the binary is compiled to use numa" << std::endl;
    }
    else if (-1 == numa_max_node())
    {
        ISAAC_THREAD_CERR << "WARNING: numa_max_node returned -1 while the binary is compiled to use numa" << std::endl;
    }
    else
    {
        int runOnNode = threadNumber % (numa_max_node() + 1);
        ISAAC_ASSERT_MSG(8 * sizeof(unsigned long) >= unsigned(runOnNode), "numa node is too high: " << runOnNode);
        ISAAC_ASSERT_MSG(-1 != numa_run_on_node(runOnNode), "numa_run_on_node " << runOnNode <<
            " failed, errno: " << errno  << ":" << strerror(errno));
        unsigned long nodemask = 1UL << runOnNode;
        ISAAC_ASSERT_MSG(-1 != set_mempolicy(MPOL_BIND/*|MPOL_F_STATIC_NODES*/, &nodemask, sizeof(nodemask) * 8),
                         "set_mempolicy for nodemask: " << nodemask <<
                         " failed, errno: " << errno << ":" << strerror(errno));
    }
#endif //HAVE_NUMA
}

void Build::run(common::ScopedMallocBlock &mallocBlock)
{
    alignment::BinMetadataCRefList::const_iterator nextUnprocessedBinIt(bins_.begin());
    alignment::BinMetadataCRefList::const_iterator nextUnallocatedBinIt(bins_.begin());
    alignment::BinMetadataCRefList::const_iterator nextUnloadedBinIt(bins_.begin());
    alignment::BinMetadataCRefList::const_iterator nextUncompressedBinIt(bins_.begin());
    alignment::BinMetadataCRefList::const_iterator nextUnsavedBinIt(bins_.begin());

    threads_.execute(boost::bind(&Build::sortBinParallel, this,
                                boost::ref(nextUnprocessedBinIt),
                                boost::ref(nextUnallocatedBinIt),
                                boost::ref(nextUnloadedBinIt),
                                boost::ref(nextUncompressedBinIt),
                                boost::ref(nextUnsavedBinIt),
                                boost::ref(mallocBlock),
                                _1));

    unsigned fileIndex = 0;
    BOOST_FOREACH(const boost::filesystem::path &bamFilePath, barcodeBamMapping_.getPaths())
    {
        // some of the streams are null_sink (that's when reference is unmapped for the sample).
        // this is the simplest way to ignore them...
        std::ostream *stm = bamFileStreams_.at(fileIndex).get();
        if (stm)
        {
            bam::serializeBgzfFooter(*stm);
            stm->flush();
            ISAAC_THREAD_CERR << "BAM file generated: " << bamFilePath << "\n";
            bamIndexes_.at(fileIndex).flush();
            ISAAC_THREAD_CERR << "BAM index generated for " << bamFilePath << "\n";
        }
        ++fileIndex;
    }
}

void Build::dumpStats(const boost::filesystem::path &statsXmlPath)
{
    BuildStatsXml statsXml(sortedReferenceMetadataList_, bins_, barcodeMetadataList_, stats_);
    std::ofstream os(statsXmlPath.string().c_str());
    if (!os) {
        BOOST_THROW_EXCEPTION(common::IoException(errno, "ERROR: Unable to open file for writing: " + statsXmlPath.string()));
    }
    statsXml.serialize(os);
}

unsigned long Build::estimateOptimumFragmentsPerBin(
    const unsigned int estimatedFragmentSize,
    const unsigned long availableMemory,
    const double expectedBgzfCompressionRatio,
    const unsigned computeThreads)
{
//    const size_t maxFragmentIndexBytes = std::max(sizeof(io::RStrandOrShadowFragmentIndex),
//                                                         sizeof(io::FStrandFragmentIndex));

    const std::size_t maxFragmentDedupedIndexBytes = sizeof(PackedFragmentBuffer::Index);
    const std::size_t maxFragmentCompressedBytes = estimatedFragmentSize * expectedBgzfCompressionRatio;

    const std::size_t fragmentMemoryRequirements =
        // assume the initial indexes don't stay in memory for too long //maxFragmentIndexBytes              //index containing duplicates
        + estimatedFragmentSize                 //data
        + maxFragmentDedupedIndexBytes     //deduplicated index
        + maxFragmentCompressedBytes       //bgzf chunk
        ;

    // reasonable amount of bins-in-progress to allow for no-delay input/compute/output overlap
//    const unsigned minOverlap = 3;;
    // try to increase granularity so that the CPU gets efficiently utilized.
    const unsigned minOverlap = computeThreads;
    return availableMemory / fragmentMemoryRequirements / minOverlap;
}

/**
 * \brief Attempts to reserve memory buffers required to process a bin.
 *
 * \return Non-zero size in bytes if the reservation failed.
 */
unsigned long Build::reserveBuffers(
    boost::unique_lock<boost::mutex> &lock,
    const alignment::BinMetadataCRefList::const_iterator thisThreadBinIt,
    const alignment::BinMetadataCRefList::const_iterator binsEnd,
    common::ScopedMallocBlock &mallocBlock,
    const std::size_t threadNumber,
    boost::shared_ptr<BinData> &binDataPtr)
{
    common::unlock_guard<boost::unique_lock<boost::mutex> > unlock(lock);
    boost::ptr_vector<boost::iostreams::filtering_ostream> &bgzfStreams = threadBgzfStreams_.at(threadNumber);
    boost::ptr_vector<bam::BamIndexPart> &bamIndexParts = threadBamIndexParts_.at(threadNumber);
    const alignment::BinMetadata &bin = *thisThreadBinIt;
    // bin stats have an entry per filtered bin reference.
    const unsigned binStatsIndex = std::distance(bins_.begin(), thisThreadBinIt);
    common::ScopedMallocBlockUnblock unblockMalloc(mallocBlock);
    try
    {
        binDataPtr = boost::shared_ptr<BinData>(
            new BinData(realignedGapsPerFragment_,
                        barcodeBamMapping_, barcodeMetadataList_,
                        realignGaps_, bin, binStatsIndex, tileMetadataList_, contigMap_, contigLists_, maxReadLength_,
                        forcedDodgyAlignmentScore_,  flowcellLayoutList_, includeTags_, pessimisticMapQ_, splitGapLength_));

        unsigned outputFileIndex = 0;
        BOOST_FOREACH(std::vector<char> &bgzfBuffer, threadBgzfBuffers_.at(threadNumber))
        {
            bgzfBuffer.reserve(estimateBinCompressedDataRequirements(bin, outputFileIndex++));
        }

        ISAAC_ASSERT_MSG(!bgzfStreams.size(), "Expecting empty pool of streams");
        while(bgzfStreams.size() < bamFileStreams_.size())
        {
            bgzfStreams.push_back(new boost::iostreams::filtering_ostream);
            bgzfStreams.back().push(bgzf::BgzfCompressor(bamGzipLevel_));
            bgzfStreams.back().push(
                boost::iostreams::back_insert_device<std::vector<char> >(
                    threadBgzfBuffers_.at(threadNumber).at(bgzfStreams.size()-1)));
        }

        ISAAC_ASSERT_MSG(!bamIndexParts.size(), "Expecting empty pool of bam index parts");
        while(bamIndexParts.size() < bamFileStreams_.size())
        {
            bamIndexParts.push_back(new bam::BamIndexPart);
        }
    }
    catch(std::bad_alloc &e)
    {
        errno = 0;
        bgzfStreams.clear();
        bamIndexParts.clear();
        // give a chance other threads to allocate what they need... TODO: this is not required anymore as allocation happens orderly
        binDataPtr.reset();
        unsigned long totalBuffersNeeded = 0UL;
        unsigned outputFileIndex = 0;
        BOOST_FOREACH(std::vector<char> &bgzfBuffer, threadBgzfBuffers_.at(threadNumber))
        {
            std::vector<char>().swap(bgzfBuffer);
            totalBuffersNeeded += estimateBinCompressedDataRequirements(bin, outputFileIndex++);
        }
        // reset errno, to prevent misleading error messages when failing code does not set errno
        errno = 0;
        return BinData::getMemoryRequirements(bin) + totalBuffersNeeded;
    }
    return 0;
}


boost::shared_ptr<BinData> Build::allocateBin(
    boost::unique_lock<boost::mutex> &lock,
    const alignment::BinMetadataCRefList::const_iterator thisThreadBinIt,
    const alignment::BinMetadataCRefList::const_iterator binsEnd,
    alignment::BinMetadataCRefList::const_iterator &nextUnallocatedBinIt,
    common::ScopedMallocBlock &mallocBlock,
    const std::size_t threadNumber)
{
    boost::shared_ptr<BinData> ret;
    unsigned long requiredMemory = 0;
    while(nextUnallocatedBinIt != thisThreadBinIt ||
          0 != (requiredMemory = reserveBuffers(lock, thisThreadBinIt, binsEnd, mallocBlock, threadNumber, ret)))
    {
        if (forceTermination_)
        {
            BOOST_THROW_EXCEPTION(common::ThreadingException("Terminating due to failures on other threads"));
        }

        if (nextUnallocatedBinIt == thisThreadBinIt)
        {
            ISAAC_THREAD_CERR << "WARNING: Holding up processing of bin: " <<
                thisThreadBinIt->get().getPath() << " until " << requiredMemory <<
                " bytes of allowed memory is available." << std::endl;
        }
        stateChangedCondition_.wait(lock);
    }
    ++nextUnallocatedBinIt;
    stateChangedCondition_.notify_all();
    return ret;
}

void Build::waitForLoadSlot(
    boost::unique_lock<boost::mutex> &lock,
    const alignment::BinMetadataCRefList::const_iterator thisThreadBinIt,
    const alignment::BinMetadataCRefList::const_iterator binsEnd,
    alignment::BinMetadataCRefList::const_iterator &nextUnloadedBinIt,
    common::ScopedMallocBlock &mallocBlock,
    const std::size_t threadNumber)
{
    while(nextUnloadedBinIt != thisThreadBinIt || !maxLoaders_)
    {
        if (forceTermination_)
        {
            BOOST_THROW_EXCEPTION(common::ThreadingException("Terminating due to failures on other threads"));
        }

        if (nextUnloadedBinIt == thisThreadBinIt)
        {
            ISAAC_THREAD_CERR << "WARNING: Holding up processing of bin: " <<
                thisThreadBinIt->get().getPath() << " until a load slot is available" << std::endl;
        }
        stateChangedCondition_.wait(lock);
    }

    ++nextUnloadedBinIt;
    --maxLoaders_;
    stateChangedCondition_.notify_all();
}

void Build::returnLoadSlot(const bool exceptionUnwinding)
{
    ++maxLoaders_;
    if (exceptionUnwinding)
    {
        forceTermination_ = true;
    }
    stateChangedCondition_.notify_all();
}

void Build::waitForComputeSlot(
    boost::unique_lock<boost::mutex> &lock,
    const alignment::BinMetadataCRefList::const_iterator thisThreadBinIt,
    alignment::BinMetadataCRefList::const_iterator &nextUncompressedBinIt)
{
    const unsigned binIndex = std::distance(bins_.begin(), thisThreadBinIt);
    computeSlotWaitingBins_.push_back(binIndex);

    while(!maxComputers_ ||
        binIndex != *std::min_element(computeSlotWaitingBins_.begin(), computeSlotWaitingBins_.end()))
    {
        if (forceTermination_)
        {
            BOOST_THROW_EXCEPTION(common::ThreadingException("Terminating due to failures on other threads"));
        }

        if (nextUncompressedBinIt == thisThreadBinIt)
        {
            ISAAC_THREAD_CERR << "WARNING: Holding up processing of bin: " <<
                thisThreadBinIt->get().getPath() << " until a compute slot is available." << std::endl;
        }

        stateChangedCondition_.wait(lock);
    }
    --maxComputers_;
    computeSlotWaitingBins_.erase(std::find(computeSlotWaitingBins_.begin() ,computeSlotWaitingBins_.end(), binIndex));
    ++nextUncompressedBinIt;
    stateChangedCondition_.notify_all();
}

void Build::returnComputeSlot(const bool exceptionUnwinding)
{
    ++maxComputers_;
    if (exceptionUnwinding)
    {
        forceTermination_ = true;
    }
    stateChangedCondition_.notify_all();
}

void Build::waitForSaveSlot(
    boost::unique_lock<boost::mutex> &lock,
    const alignment::BinMetadataCRefList::const_iterator thisThreadBinIt,
    alignment::BinMetadataCRefList::const_iterator &nextUnsavedBinIt)
{
    while(nextUnsavedBinIt != thisThreadBinIt)
    {
        if (forceTermination_)
        {
            BOOST_THROW_EXCEPTION(common::ThreadingException("Terminating due to failures on other threads"));
        }
        stateChangedCondition_.wait(lock);
    }
}

void Build::returnSaveSlot(
    alignment::BinMetadataCRefList::const_iterator &nextUnsavedBinIt,
    const bool exceptionUnwinding)
{
    ++nextUnsavedBinIt;
    if (exceptionUnwinding)
    {
        forceTermination_ = true;
    }
    stateChangedCondition_.notify_all();
}

void Build::sortBinParallel(alignment::BinMetadataCRefList::const_iterator &nextUnprocessedBinIt,
                            alignment::BinMetadataCRefList::const_iterator &nextUnallocatedBinIt,
                            alignment::BinMetadataCRefList::const_iterator &nextUnloadedBinIt,
                            alignment::BinMetadataCRefList::const_iterator &nextUncompressedBinIt,
                            alignment::BinMetadataCRefList::const_iterator &nextUnsavedBinIt,
                            common::ScopedMallocBlock &mallocBlock,
                            const std::size_t threadNumber)
{
    boost::unique_lock<boost::mutex> lock(stateMutex_);
    while(bins_.end() != nextUnprocessedBinIt)
    {
        alignment::BinMetadataCRefList::const_iterator thisThreadBinIt = nextUnprocessedBinIt++;

        // wait and allocate memory required for loading and compressing this bin
        boost::shared_ptr<BinData> binDataPtr =
            allocateBin(lock, thisThreadBinIt, bins_.end(), nextUnallocatedBinIt, mallocBlock, threadNumber);
        waitForLoadSlot(lock, thisThreadBinIt, bins_.end(), nextUnloadedBinIt, mallocBlock, threadNumber);
        ISAAC_BLOCK_WITH_CLENAUP(boost::bind(&Build::returnLoadSlot, this, _1))
        {
            {
                common::unlock_guard<boost::unique_lock<boost::mutex> > unlock(lock);
                BinLoader binLoader;
                binLoader.loadData(*binDataPtr);
            }
        }

        waitForComputeSlot(lock, thisThreadBinIt, nextUncompressedBinIt);
        ISAAC_BLOCK_WITH_CLENAUP(boost::bind(&Build::returnComputeSlot, this, _1))
        {
            {
                common::unlock_guard<boost::unique_lock<boost::mutex> > unlock(lock);
                binSorter_.resolveDuplicates(*binDataPtr, stats_);
                if (!binDataPtr->isUnalignedBin() && REALIGN_NONE != realignGaps_)
                {
                    boost::unique_lock<boost::mutex> lock(realignMutex_);
                    gapRealigner_.realignGaps(*binDataPtr);
                }

                binSorter_.serialize(
                    *binDataPtr, threadBgzfStreams_.at(threadNumber), threadBamIndexParts_.at(threadNumber));
                threadBgzfStreams_.at(threadNumber).clear();
            }
            // give back some memory to allow other threads to load
            // data while we're waiting for our turn to save
            binDataPtr.reset();
        }

        // wait for our turn to store bam data
        waitForSaveSlot(lock, thisThreadBinIt, nextUnsavedBinIt);
        ISAAC_BLOCK_WITH_CLENAUP(boost::bind(&Build::returnSaveSlot, this, boost::ref(nextUnsavedBinIt), _1))
        {
            saveAndReleaseBuffers(lock, thisThreadBinIt->get().getPath(), threadNumber);
        }
    }
}

/**
 * \brief Save bgzf compressed buffers into corresponding sample files and and release associated memory
 */
void Build::saveAndReleaseBuffers(
    boost::unique_lock<boost::mutex> &lock,
    const boost::filesystem::path &filePath,
    const std::size_t threadNumber)
{
    unsigned index = 0;
    BOOST_FOREACH(std::vector<char> &bgzfBuffer, threadBgzfBuffers_.at(threadNumber))
    {
        {
            common::unlock_guard<boost::unique_lock<boost::mutex> > unlock(lock);
            std::ostream *stm = bamFileStreams_.at(index).get();
            if (!stm)
            {
                ISAAC_ASSERT_MSG(bgzfBuffer.empty(), "Unexpected data for bam file belonging to a sample with unmapped reference");
            }
            else
            {
                saveBuffer(bgzfBuffer, *stm, threadBamIndexParts_.at(threadNumber).at(index), bamIndexes_.at(index), filePath);
            }
        }
        // release rest of the memory that was reserved for this bin
        std::vector<char>().swap(bgzfBuffer);
        ++index;
    }
    threadBamIndexParts_.at(threadNumber).clear();
}

void Build::saveBuffer(
    const std::vector<char> &bgzfBuffer,
    std::ostream &bamStream,
    const bam::BamIndexPart &bamIndexPart,
    bam::BamIndex &bamIndex,
    const boost::filesystem::path &filePath)
{
    ISAAC_THREAD_CERR << "Saving " << bgzfBuffer.size() << " bytes of sorted data for bin " << filePath << std::endl;
    const clock_t start = clock();
    if(!bgzfBuffer.empty() && !bamStream.write(&bgzfBuffer.front(), bgzfBuffer.size())/* ||
        !bamStream.strict_sync()*/){
        BOOST_THROW_EXCEPTION(common::IoException(
            errno, (boost::format("Failed to write bgzf block of %d bytes into bam stream") % bgzfBuffer.size()).str()));
    }
    bamIndex.processIndexPart( bamIndexPart, bgzfBuffer );

    ISAAC_THREAD_CERR << "Saving " << bgzfBuffer.size() << " bytes of sorted data for bin " << filePath << " done in " << (clock() - start) / 1000 << "ms\n";
}

} // namespace build
} // namespace isaac
