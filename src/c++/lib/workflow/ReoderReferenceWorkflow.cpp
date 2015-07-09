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
 ** \file ReoderReferenceWorkflow.cpp
 **
 ** \brief see ReoderReferenceWorkflow.hh
 **
 ** \author Roman Petrovski
 **/

#include <boost/algorithm/string/join.hpp>
#include <boost/foreach.hpp>
#include <boost/format.hpp>

#include "common/Debug.hh"
#include "common/Exceptions.hh"
#include "reference/AnnotationLoader.hh"
#include "reference/ContigLoader.hh"
#include "workflow/ReorderReferenceWorkflow.hh"

namespace isaac
{
namespace workflow
{

ReorderReferenceWorkflow::ReorderReferenceWorkflow(
    const bfs::path &sortedReferenceMetadata,
    const bfs::path &newXmlPath,
    const bfs::path &newDataFileDirectory,
    const std::vector<std::string> &newOrder,
    const unsigned basesPerLine
    )
    : sortedReferenceMetadata_(sortedReferenceMetadata),
      newXmlPath_(newXmlPath),
      newDataFileDirectory_(newDataFileDirectory),
      newOrder_(newOrder),
      basesPerLine_(basesPerLine),
      sameOrder_(false),
      threads_(boost::thread::hardware_concurrency()),
      xml_(reference::loadSortedReferenceXml(sortedReferenceMetadata_)),
      originalKaryotypeIndexes_(xml_.getContigs().size()),
      originalContigOffsets_(reference::computeContigOffsets(xml_.getKaryotypeOrderedContigs()))
{
    if (!newOrder_.empty())
    {
        {
            const reference::SortedReferenceMetadata::Contigs contigs = xml_.getKaryotypeOrderedContigs();
            std::vector<std::string> oldOrder;
            std::transform(contigs.begin(), contigs.end(), std::back_inserter(oldOrder),
                           boost::bind(&reference::SortedReferenceMetadata::Contig::name_, _1));

            ISAAC_THREAD_CERR << "old order: " << boost::algorithm::join(oldOrder, ",") << std::endl;
            ISAAC_THREAD_CERR << "new order: " << boost::algorithm::join(newOrder_, ",") << std::endl;

            sameOrder_ = (oldOrder == newOrder_);
        }

        std::vector<bool> present(newOrder_.size());
        BOOST_FOREACH(reference::SortedReferenceMetadata::Contig &xmlContig, xml_.getContigs())
        {
            const std::vector<std::string>::const_iterator newOrderIt = std::find(newOrder_.begin(), newOrder_.end(), xmlContig.name_);
            if (newOrder_.end() == newOrderIt)
            {
                BOOST_THROW_EXCEPTION(isaac::common::InvalidParameterException(
                    "Contig name not listed in the new order: " + xmlContig.name_));
            }
            const unsigned newKaryotypeIndex = newOrderIt - newOrder_.begin();
            present.at(newKaryotypeIndex) = true;
            originalKaryotypeIndexes_.at(newKaryotypeIndex) = xmlContig.karyotypeIndex_;
            xmlContig.karyotypeIndex_ = newKaryotypeIndex;
        }

        const std::vector<bool>::const_iterator firstUnused = std::find(present.begin(), present.end(), false);
        if(present.end() != firstUnused)
        {
            BOOST_THROW_EXCEPTION(isaac::common::InvalidParameterException(
                "Contig name listed in the new order not found in the reference: " + newOrder_.at(firstUnused - present.begin())));
        }

    }
    else
    {
        sameOrder_ = true;
        ISAAC_THREAD_CERR << "Preserving the existing order of contigs" << std::endl;
    }
}

bool ReorderReferenceWorkflow::orderByKaryotypeIndex(const reference::Contig& left, const reference::Contig& right)
{
    return xml_.getContigs().at(left.index_).karyotypeIndex_ < xml_.getContigs().at(right.index_).karyotypeIndex_;
}

void ReorderReferenceWorkflow::reorderContigs()
{
    if (!newOrder_.empty() && (!sameOrder_ || !xml_.singleFileReference()))
    {
        boost::filesystem::path newFastaPath = (newDataFileDirectory_ / "genome.fa");
        std::ofstream fastaOs(newFastaPath.c_str());
        if (!fastaOs)
        {
            BOOST_THROW_EXCEPTION(isaac::common::IoException(errno, "Failed to open output file: " + newFastaPath.string()));
        }

        std::vector<reference::Contig> contigs = reference::loadContigs(xml_.getContigs(), threads_);
        std::sort(contigs.begin(), contigs.end(),
                  boost::bind(&ReorderReferenceWorkflow::orderByKaryotypeIndex, this, _1, _2));

        BOOST_FOREACH(const reference::Contig &contig, contigs)
        {
            storeContig(fastaOs, contig, newFastaPath);
        }
    }
    else
    {
        const boost::filesystem::path &srcPath = xml_.getContigs().front().filePath_;
        const boost::filesystem::path targetPath = newDataFileDirectory_ / srcPath.filename();
        boost::filesystem::copy_file(srcPath, targetPath);
        ISAAC_THREAD_CERR << "Copied entire reference from " << srcPath << " to " << targetPath << std::endl;
        BOOST_FOREACH(reference::SortedReferenceMetadata::Contig &xmlContig, xml_.getContigs())
        {
            xmlContig.filePath_ = targetPath;
        }
    }
}

void ReorderReferenceWorkflow::reorderAnnotation()
{
    if (xml_.hasKUniquenessAnnotation())
    {
        const reference::SortedReferenceMetadata::AnnotationFile &annotationFile = xml_.getKUniquenessAnnotation();

        const boost::filesystem::path &srcPath = annotationFile.path_;
        const boost::filesystem::path targetPath = newDataFileDirectory_ / srcPath.filename();
        if (!newOrder_.empty() && !sameOrder_)
        {
            const reference::SortedReferenceMetadata::Contigs contigs = xml_.getKaryotypeOrderedContigs();
            boost::iostreams::filtering_ostream annotationOs;
            if (common::isDotGzPath(targetPath))
            {
                annotationOs.push(boost::iostreams::gzip_compressor());
            }
            annotationOs.push(boost::iostreams::file_sink(targetPath.string()));
            if (!annotationOs)
            {
                BOOST_THROW_EXCEPTION(isaac::common::IoException(errno, "Failed to open output file: " + targetPath.string()));
            }

            const reference::Annotation annotation =
                reference::loadAnnotationBlob<reference::Annotation>(annotationFile.path_, reference::genomeLength(contigs));

            BOOST_FOREACH(const reference::SortedReferenceMetadata::Contig &contig, contigs)
            {
                const unsigned long originalContigOffset = originalContigOffsets_.at(originalKaryotypeIndexes_.at(contig.karyotypeIndex_));
                ISAAC_THREAD_CERR << "Storing " << contig << " annotation from original offset: " << originalContigOffset << std::endl;
                if (!annotationOs.write(
                    reinterpret_cast<const char *>(&annotation.at(originalContigOffset)),
                    contig.totalBases_ * sizeof(reference::Annotation::value_type)))
                {
                    BOOST_THROW_EXCEPTION(isaac::common::IoException(errno, "Failed to write into output file: " + targetPath.string()));
                }
            }
        }
        else
        {
            boost::filesystem::copy_file(srcPath, targetPath);
            ISAAC_THREAD_CERR << "Copied entire annotation from " << srcPath << " to " << targetPath << std::endl;
        }
        xml_.setKUniquenessAnnotation(targetPath, annotationFile.k_);
    }
}


void ReorderReferenceWorkflow::run()
{
    std::ofstream xmlOs(newXmlPath_.c_str());
    if (!xmlOs)
    {
        BOOST_THROW_EXCEPTION(isaac::common::IoException(errno, "Failed to open output file: " + newXmlPath_.string()));
    }

    reorderContigs();
    reorderAnnotation();

    saveSortedReferenceXml(xmlOs, xml_);
}

void ReorderReferenceWorkflow::writeBase(std::ostream &os, const char base, const bool writeNewline)
{
    if (!(os << base) || (writeNewline && !(os << "\n")))
    {
        BOOST_THROW_EXCEPTION(
            isaac::common::IoException(errno, (boost::format("Failed to write data into output file %s: ") %
                newDataFileDirectory_.string()).str()));
    }
}

void ReorderReferenceWorkflow::storeContig(
    std::ostream &os,
    const reference::Contig &contig,
    const boost::filesystem::path &fastaPath)
{
    reference::SortedReferenceMetadata::Contig &xmlContig = xml_.getContigs().at(contig.index_);

    if (!(os << ">" << xmlContig.name_ << std::endl))
    {
        BOOST_THROW_EXCEPTION(
            isaac::common::IoException(errno, (boost::format("Failed to write contig name %s into output file: ") %
                xmlContig.name_ % fastaPath.string()).str()));
    }

    unsigned long startPos = os.tellp();
    std::vector<char>::const_iterator current = contig.forward_.begin();
    while(contig.forward_.end() != current)
    {
        writeBase(os, *current,
                  (!((current - contig.forward_.begin() + 1) % basesPerLine_)) ||
                  (contig.forward_.end() - 1 == current));
        ++current;
    }
    unsigned long endPos = os.tellp();

    xmlContig.filePath_ = fastaPath;
    xmlContig.offset_ = startPos;
    xmlContig.size_ = endPos - startPos;

    ISAAC_THREAD_CERR << "Stored contig: " << xmlContig.name_ << std::endl;
}

} // namespace workflow
} // namespace isaac
