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
 ** \file SortedReferenceMetadata.hh
 **
 ** Information about the pre-processed reference data files.
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_REFERENCE_SORTED_REFERENCE_METADATA_HH
#define iSAAC_REFERENCE_SORTED_REFERENCE_METADATA_HH

#include <numeric>

#include <boost/filesystem.hpp>
#include <boost/foreach.hpp>

#include "oligo/Kmer.hh"

namespace isaac
{
namespace reference
{

class SortedReferenceMetadata
{
public:
    static const unsigned OLDEST_SUPPORTED_REFERENCE_FORMAT_VERSION = 6;
    static const unsigned CURRENT_REFERENCE_FORMAT_VERSION = 6;

    struct Contig
    {
        Contig() : index_(0), karyotypeIndex_(0), offset_(0), size_(0), genomicPosition_(0), totalBases_(0), acgtBases_(0){}
        Contig(const unsigned int index, const unsigned karyotypeIndex, const std::string &name,
               const boost::filesystem::path &filePath, const unsigned long offset, const unsigned long size,
               const unsigned long genomicPosition, const unsigned long totalBases, const unsigned long acgtBases,
               const std::string &bamSqAs, const std::string &bamSqUr, const std::string &bamM5) :
                   index_(index), karyotypeIndex_(karyotypeIndex), name_(name), filePath_(filePath),
                   offset_(offset), size_(size),
                   genomicPosition_(genomicPosition), totalBases_(totalBases), acgtBases_(acgtBases),
                   bamSqAs_(bamSqAs), bamSqUr_(bamSqUr), bamM5_(bamM5) {}
        unsigned int index_;
        unsigned karyotypeIndex_;
        std::string name_;
        boost::filesystem::path filePath_;
        unsigned long offset_;
        unsigned long size_;
        unsigned long genomicPosition_;
        unsigned long totalBases_;
        unsigned long acgtBases_;
        std::string bamSqAs_;
        std::string bamSqUr_;
        std::string bamM5_;

        bool operator == (const Contig &that) const
        {
            return index_ == that.index_ && karyotypeIndex_ == that.karyotypeIndex_ &&
                name_ == that.name_ &&
                ((bamM5_.empty() && filePath_ == that.filePath_) || bamM5_ == that.bamM5_) &&
                offset_ == that.offset_ &&
                size_ == that.size_ && genomicPosition_ == that.genomicPosition_ && totalBases_ == that.totalBases_ &&
                acgtBases_ == that.acgtBases_ && bamSqAs_ == that.bamSqAs_ && bamSqUr_ == that.bamSqUr_;
        }
    };
    typedef std::vector<Contig> Contigs;

    struct MaskFile
    {
        MaskFile(): maskWidth(0), mask_(0), kmers(0){}
        MaskFile(
            const boost::filesystem::path &p,
            const unsigned mw,
            const unsigned m,
            const std::size_t km) : path(p), maskWidth(mw), mask_(m), kmers(km){}
        boost::filesystem::path path;
        unsigned maskWidth;
        unsigned mask_;
        std::size_t kmers;
        template<class Archive> friend void serialize(Archive & ar, MaskFile &, const unsigned int file_version);
    };
    typedef std::vector<MaskFile> MaskFiles;
    typedef std::map<unsigned, MaskFiles> AllMaskFiles;

    struct AnnotationFile
    {
        AnnotationFile(): k_(0){}
        AnnotationFile(
            const boost::filesystem::path &p,
            const unsigned k) : path_(p), k_(k){}
        boost::filesystem::path path_;
        unsigned k_;
        template<class Archive> friend void serialize(Archive & ar, AnnotationFile &, const unsigned int file_version);
        friend std::ostream& operator <<(std::ostream &os, const AnnotationFile& annotationFile)
        {
            return os << "AnnotationFile(" <<
                annotationFile.k_ << "," << annotationFile.path_ <<
                ")";
        }
    };
    typedef std::vector<AnnotationFile> AnnotationFiles;

private:
    AllMaskFiles maskFiles_;
    AnnotationFiles annotationFiles_;
    Contigs contigs_;
    unsigned formatVersion_;

public:
    SortedReferenceMetadata() :
        formatVersion_(CURRENT_REFERENCE_FORMAT_VERSION)
    {
    }

    void makeAbsolutePaths(const boost::filesystem::path &basePath);

    void putContig(const unsigned long genomicOffset,
                   const std::string& name,
                   const boost::filesystem::path &sequencePath,
                   const unsigned long byteOffset,
                   const unsigned long byteSize,
                   const unsigned long totalBases,
                   const unsigned long acgtBases,
                   const unsigned index,
                   const unsigned karyotypeIndex,
                   const std::string &bamSqAs,
                   const std::string &bamSqUr,
                   const std::string &bamM5);
    void addMaskFile(
        const unsigned seedLength,
        const unsigned int maskWidth,
        const unsigned mask, const boost::filesystem::path &filePath,
        const std::size_t kmers);

    /**
     ** Precondition: the contigs in the current instance are sequentially
     ** indexed from 0 and there are no duplicates.
     **/
    const Contigs &getContigs() const {return contigs_;}
    Contigs &getContigs() {return contigs_;}

    std::size_t getContigsCount() const {return contigs_.size();}

    /**
     ** \brief Return the number of contigs for which includeContig() returns true
     **/
    template <typename IncludeContigF>
    std::size_t getContigsCount(
        const IncludeContigF &includeContig) const
    {
        return std::count_if(contigs_.begin(), contigs_.end(),
                             boost::bind(includeContig, boost::bind(&Contig::index_, _1)));
    }

    /**
     ** \brief Return a list of contig where each contig is at the corresponding karyotype index.
     **
     ** Precondition: the contigs in the current instance are sequentially
     ** indexed from 0 and there are no duplicates.
     **/
    Contigs getKaryotypeOrderedContigs() const;

    /**
     ** \brief Return a list of contig where each contig is at the corresponding karyotype index.
     **        Only the contigs for which includeContig() returns true are included
     **
     ** Precondition: the contigs in the current instance are sequentially
     ** indexed from 0 and there are no duplicates.
     **/
    template <typename IncludeContigF>
    Contigs getKaryotypeOrderedContigs(
        IncludeContigF &includeContig) const
    {
        Contigs ret;
        ret.reserve(contigs_.size());
        std::remove_copy_if(contigs_.begin(), contigs_.end(), std::back_inserter(ret),
                            !boost::bind(includeContig, boost::bind(&Contig::karyotypeIndex_, _1)));
        std::sort(ret.begin(), ret.end(),
                  boost::bind(&Contig::karyotypeIndex_, _1) <
                  boost::bind(&Contig::karyotypeIndex_, _2));
        return ret;
    }


    /**
     ** \return Total number of kmers in all mask files
     **/
    unsigned long getTotalKmers(const unsigned seedLength) const;

    bool supportsSeedLength(const unsigned seedLength) const
    {
        return maskFiles_.end() != maskFiles_.find(seedLength);
    }
    const MaskFiles &getMaskFileList(const unsigned seedLength) const {return maskFiles_.at(seedLength);}
    MaskFiles &getMaskFileList(const unsigned seedLength) {return maskFiles_[seedLength];}

    bool hasKUniquenessAnnotation() const {return !annotationFiles_.empty();}
    const AnnotationFile& getKUniquenessAnnotation() const
        {ISAAC_ASSERT_MSG(hasKUniquenessAnnotation(), "Annotation requested for reference that does not have one"); return annotationFiles_.front();}
    void setKUniquenessAnnotation(const boost::filesystem::path &path, const unsigned k)
        {ISAAC_ASSERT_MSG(annotationFiles_.size() < 2, "Multiple annotations are not supported");
         annotationFiles_.resize(1); annotationFiles_.front() = AnnotationFile(path, k);}

    void clearAnnotations() {annotationFiles_.clear();}

    void clearMasks() {maskFiles_.clear();}

    void merge(SortedReferenceMetadata &that);

    bool singleFileReference() const;

    template<class Archive> friend void serialize(Archive & ar, SortedReferenceMetadata &, const unsigned int file_version);
    template<class Archive> friend void serialize(Archive & ar, const SortedReferenceMetadata &, const unsigned int file_version);

};

typedef std::vector<SortedReferenceMetadata> SortedReferenceMetadataList;

inline std::size_t genomeLength(const SortedReferenceMetadata::Contigs &contigList)
{
    return std::accumulate(
        contigList.begin(), contigList.end(),
        std::size_t(0), boost::bind<std::size_t>(std::plus<std::size_t>(), _1, boost::bind(&SortedReferenceMetadata::Contig::totalBases_, _2)));
}


inline unsigned long getLongestGenomeLength(const isaac::reference::SortedReferenceMetadataList &sortedReferenceMetadataList)
{
    std::size_t longestGenome = 0;
    BOOST_FOREACH(const isaac::reference::SortedReferenceMetadata &sortedReferenceMetadata, sortedReferenceMetadataList)
    {
        const std::vector<isaac::reference::SortedReferenceMetadata::Contig> xmlContigs = sortedReferenceMetadata.getContigs();
        const std::size_t genomeLength = reference::genomeLength(xmlContigs);
        longestGenome = std::max(longestGenome, genomeLength);
    }
    return longestGenome;
}

/**
 * \brief Builds a vector of global starts of contigs (all contig bases are considered)
 *        for the contigs in a given order
 */
inline std::vector<unsigned long> computeContigOffsets(const reference::SortedReferenceMetadata::Contigs &contigs)
{
    std::vector<unsigned long> ret(contigs.size());

    unsigned long lastOffset = 0UL;

    // the caller decides whether contigs are ordered by index_ or karyotypeIndex_
    unsigned orderIndex = 0;
    BOOST_FOREACH(const reference::SortedReferenceMetadata::Contig &contig, contigs)
    {
        ret.at(orderIndex) = lastOffset;
        lastOffset += contig.totalBases_;
        ++orderIndex;
    }

    return ret;
}


inline std::ostream &operator <<(std::ostream &os, const isaac::reference::SortedReferenceMetadata::Contig &xmlContig)
{
    return os << "SortedReferenceMetadata::Contig(" <<
        xmlContig.name_ << "," <<
        xmlContig.genomicPosition_ << "pos," <<
        xmlContig.totalBases_ << "tb," <<
        xmlContig.offset_ << "off" <<
        ")";
}

} // namespace reference
} // namespace isaac

#endif // #ifndef iSAAC_REFERENCE_SORTED_REFERENCE_METADATA_HH
