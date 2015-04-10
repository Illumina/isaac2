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
 ** \file FastqLoader.hh
 **
 ** Component to read FASTQ files.
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_IO_FASTQ_LOADER_HH
#define iSAAC_IO_FASTQ_LOADER_HH

#include "common/Threads.hpp"
#include "io/FastqReader.hh"

namespace isaac
{
namespace io
{

//#pragma GCC push_options
//#pragma GCC optimize ("0")

class FastqLoader
{
    const unsigned inputLoadersMax_;
    std::vector<boost::shared_ptr<FastqReader> >  readReaders_;
    bool paired_;
    common::ThreadVector &threads_;

public:
    /**
     * \brief initializes paired fastq loader.
     */
    /**
     * \brief Creates uninitialized fastq loader
     */
    FastqLoader(
        const bool allowVariableLength,
        const std::size_t maxPathLength,
        common::ThreadVector &threads,
        const unsigned inputLoadersMax) :
        inputLoadersMax_(inputLoadersMax),
        // notice that single-ended fastq will use only half the allowed threads to decompress bgzf.
        // That is not correct way to do it, but not particularly important. We mainly care here for
        // inputLoadersMax_=1 scenario which is important for debugging and such.
        paired_(false),
        threads_(threads)
    {
        readReaders_.resize(2);
        threads_.execute(boost::bind(&FastqLoader::initializeReaderThread, this, _1, allowVariableLength, maxPathLength), 2);
    }

    void open(
        const boost::filesystem::path &read1Path)
    {
        readReaders_[0]->open(read1Path);
        paired_ = false;
    }

    void open(
        const boost::filesystem::path &read1Path,
        const boost::filesystem::path &read2Path)
    {
        readReaders_[0]->open(read1Path);
        readReaders_[1]->open(read2Path);
        paired_ = true;
    }

    /**
     * \param clusterCount  Maximum number of clusters to load
     * \param clusterLength Expected cluster length. Will fail if the read cluster length does not match
     * \param it            Insert iterator for the buffer that is sufficient to load the clusterCount
     *                      clusters of clusterLength
     *
     * \return Actual number of loaded clusters
     */
    template <typename InsertIt>
    unsigned loadClusters(unsigned clusterCount, const unsigned nameLengthMax, const flowcell::ReadMetadataList &readMetadataList, InsertIt &it)
    {
        if (1 == readMetadataList.size())
        {
            return loadSingleRead(*readReaders_[0], clusterCount, readMetadataList.at(0), 0, nameLengthMax, it);
        }
        else
        {
            ISAAC_ASSERT_MSG(2 == readMetadataList.size(), "Only paired and single-ended data is supported");
            unsigned readClusters[2] = {0,0};
            InsertIt it1 = it;
            if (2 <= inputLoadersMax_)
            {
                it += readMetadataList.at(0).getLength();
                boost::reference_wrapper<InsertIt> insertIterators[] = {boost::ref(it1), boost::ref(it)};
                threads_.execute(boost::bind(&FastqLoader::threadLoadPairedReads<InsertIt>, this,
                                            clusterCount, boost::ref(readMetadataList), nameLengthMax,
                                            boost::ref(readClusters), insertIterators, _1),
                                 2);
            }
            else
            {
                ISAAC_ASSERT_MSG(1 == inputLoadersMax_, "At least one thread is expected for IO")
                readClusters[0] = loadSingleRead(*readReaders_[0], clusterCount, readMetadataList.at(0), readMetadataList.at(1).getLength(), 0, it1);
                it += readMetadataList.at(0).getLength();
                readClusters[1] = loadSingleRead(*readReaders_[1], clusterCount, readMetadataList.at(1), readMetadataList.at(0).getLength(), nameLengthMax, it);
            }

            if (readClusters[0] != readClusters[1])
            {
                BOOST_THROW_EXCEPTION(common::IoException(errno, (boost::format("Mismatching number of cluster read for r1/r2 = %d/%d, files: %s/%s") %
                    readClusters[0] % readClusters[1] % readReaders_[0]->getPath() % readReaders_[1]->getPath()).str()));
            }

            return readClusters[0];
        }

    }
private:
    template <typename InsertIt>
    static unsigned loadSingleRead(FastqReader &reader, unsigned clusterCount,
                            const flowcell::ReadMetadata &readMetadata,
                            const unsigned step, const unsigned nameLengthMax, InsertIt &it)
    {
        unsigned clustersToRead = clusterCount;
        for (;clustersToRead && reader.hasData();)
        {
            it = reader.extractBcl(readMetadata, it);
            if (nameLengthMax)
            {
                it = reader.extractReadName(nameLengthMax, it);
            }

            reader.next();
            // avoid debug glibc complaining about advancing iterator past the end of the container
            if (--clustersToRead)
            {
                std::advance(it, step);
            }
        }
        return clusterCount - clustersToRead;
    }

    template <typename InsertIt>
    void threadLoadPairedReads(unsigned clusterCount,
                              const flowcell::ReadMetadataList &readMetadataList,
                              const unsigned nameLengthMax,
                              unsigned readClusters[2],
                              boost::reference_wrapper<InsertIt> insertIterators[2],
                              const int threadNumber)
    {
        ISAAC_ASSERT_MSG(threadNumber < 2, "Only two threads are allowed to read paired data");
        if(!threadNumber)
        {
            readClusters[0] = loadSingleRead(
                *readReaders_[0],
                clusterCount, readMetadataList.at(0),
                readMetadataList.at(1).getLength() + nameLengthMax, 0,
                insertIterators[0].get());
        }
        else
        {
            readClusters[1] = loadSingleRead(
                *readReaders_[1],
                clusterCount, readMetadataList.at(1),
                readMetadataList.at(0).getLength(), nameLengthMax,
                insertIterators[1].get());
        }
    }

    void initializeReaderThread(const int threadNumber, const bool allowVariableLength, const std::size_t maxPathLength)
    {
        readReaders_.at(threadNumber).reset(new FastqReader(allowVariableLength, std::max(1U, inputLoadersMax_/2), maxPathLength));
    }
};

//#pragma GCC pop_options

} // namespace io
} // namespace isaac

#endif // #ifndef iSAAC_IO_FASTQ_LOADER_HH
