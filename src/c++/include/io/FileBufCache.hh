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
 ** \file FileBufCache.hh
 **
 ** Vector of file buffers which are kept open to reduce the cost of closing/opening files.
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_IO_FILE_BUF_CACHE_HH
#define iSAAC_IO_FILE_BUF_CACHE_HH

#include <boost/filesystem.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/ref.hpp>

#include "common/Exceptions.hh"
#include "io/FileBufWithReopen.hh"

namespace isaac
{
namespace io
{


template <typename FileBufT>
class FileBufHolder
{
    //TODO: either provide general-purpose implementation or fail at compile time with some explanations
};

template <>
class FileBufHolder<FileBufWithReopen>
{
public:
    typedef FileBufWithReopen FileBufType;
    typedef std::string::size_type path_size_type;
private:
    // strings are being used instead of boost::filesystem::path to have a better control over the memory
    // allocations.
    std::string filePath_;
    boost::scoped_ptr<FileBufType> fileBufPtr_;

public:

    FileBufHolder(const std::ios_base::openmode mode)
    : fileBufPtr_(new FileBufWithReopen(mode))
    {
    }

    FileBufHolder(const std::ios_base::openmode mode, const path_size_type reservePathLength)
    : fileBufPtr_(new FileBufWithReopen(mode))
    {
        reservePathBuffer(reservePathLength);
    }

    /**
     * \brief This copy constructor is only useful for reservation during initialization. It does not
     *        carry over the path, nor the handle
     */
    FileBufHolder(const FileBufHolder &that)
    : fileBufPtr_(new FileBufWithReopen(that.mode()))
    {
        filePath_.reserve(that.filePath_.capacity());
    }

    FileBufHolder& operator=(FileBufHolder that)
    {
        swap(that);
        return *this;
    }

    FileBufType *get() {return fileBufPtr_.get();}
    void reservePathBuffer(path_size_type reservePathLength)
    {
        filePath_.reserve(reservePathLength);
    }

    void reopen(const char *filePath, FileBufWithReopen::FadviseFlags fadvise)
    {
        if (!fileBufPtr_->reopen(filePath, fadvise)) {
            BOOST_THROW_EXCEPTION(common::IoException(errno, std::string("Failed to reopen a file handle for ") + filePath));
        }
        // ensure buffer does not become shared
        filePath_ = filePath;
    }

    std::ios_base::openmode mode() const {return fileBufPtr_->mode();}

    void flush()
    {
        if (!filePath_.empty())
        {
            fileBufPtr_->flush();
        }
    }

    const std::string &path() const {return filePath_;}

    void clear()
    {
        filePath_.clear();
    }


    void swap(FileBufHolder &that)
    {
        filePath_.swap(that.filePath_);
        fileBufPtr_.swap(that.fileBufPtr_);
    }
private:
    FileBufHolder();
};

/**
 * \brief The container is resizable. The file buffers are searchable by the file path.
 *        When requested file path is not found, a new entry is inserted into the sorted sequence.
 *        If there is no room for a new entry, the file buffer at the position that would hold the
 *        file path is reopened for the target file.
 *        Only one file buffer can be open for a given path at a time.
 *        The FileBufCache does not keep track of the file buffers in use by client code. It expects
 *        the client code to use one file buffer from the cache at a time.
 */
template<typename FileBufT>
class FileBufCache : public std::vector<FileBufHolder<FileBufT> >
{
public:
    typedef FileBufT FileBufType;
    typedef FileBufHolder<FileBufT> HolderType;
    typedef std::vector<HolderType > BaseType;
    typedef FileBufType *PointerType;
    typedef typename BaseType::iterator HolderIterator;

    const std::ios_base::openmode mode_;

    FileBufCache(const typename BaseType::size_type size,
                 const std::ios_base::openmode mode) : BaseType(size, HolderType(mode)), mode_(mode){
    }

    FileBufCache(const typename BaseType::size_type size,
                 const std::ios_base::openmode mode,
                 const typename HolderType::path_size_type reservePathLength) : BaseType(size, HolderType(mode)), mode_(mode)
    {
        reservePathBuffers(reservePathLength);
    }

    void reservePathBuffers(const typename HolderType::path_size_type reservePathLength)
    {
        std::for_each(this->begin(), this->end(), boost::bind(&HolderType::reservePathBuffer, _1, reservePathLength));
    }

    void unreserve()
    {
        std::vector<FileBufHolder<FileBufT> >().swap(*this);
    }

    /**
     * \brief Returns a cached handle for matching path. If missing, attempts to use an empty slot,
     *        if no slots are empty, evicts the handle at the position corresponding to the suppliled path.
     */
    PointerType get(const boost::filesystem::path &filePath, FileBufWithReopen::FadviseFlags fadvise = FileBufWithReopen::normal)
    {
        HolderIterator it = std::lower_bound(this->begin(), this->end(), filePath, compareFileBufHolderBoostPath);
        if (this->end() == it)
        {
            --it;
        }
        if (it->path() != filePath.c_str())
        {
            insertOrReopen(it, filePath, fadvise);
        }
        else
        {
            ISAAC_ASSERT_MSG(0 != it->get(), "Holders with non-empty file path must hold an open buffer.");
            // even thought the path is the same, reopen is needed to reset position and fadvise
            it->reopen(filePath.c_str(), fadvise);
        }
        return it->get();
    }

    /**
     * \brief forgets the paths of the cached file handles so that new gets are ensured to be cached as much as possible
     */
    void clear()
    {
        std::for_each(this->begin(), this->end(), boost::bind(&HolderType::clear, _1));
    }

    /**
     * \brief flushes output data if the stream is open for output
     */
    void flush()
    {
        std::for_each(this->begin(), this->end(), boost::bind(&HolderType::flush, _1));
    }

private:
    static bool compareFileBufHolderBoostPath(const HolderType &left, const boost::filesystem::path &right)
    {
        ISAAC_ASSERT_MSG(!right.empty(), "requesting empty path is not allowed");
        // all empty ones should go to the end
        // boost 1.46 filesystem::path > operator works on non-const iterators which causes the
        // underlying string to split shared buffers resulting in memory allocation. Direct string compare does not
        // have this issue.
        return left.path() > right.string();
    }

    static void swap(HolderType &left, HolderType &right)
    {
        ISAAC_ASSERT_MSG(left.mode() == right.mode(), "Access mode must be the same as it is set during the initialization");
        left.swap(right);
    }

    void insertOrReopen(HolderIterator it, const boost::filesystem::path &filePath, FileBufWithReopen::FadviseFlags fadvise)
    {
        // if the last element is not open yet, shift all to make room for new one
        if (!it->path().empty() && this->back().path().empty())
        {
            for (std::reverse_iterator<HolderIterator> itShift = this->rbegin();
                itShift < std::reverse_iterator<HolderIterator>(it+1); ++itShift)
            {
                swap(*itShift, *(itShift+1));
            }
        }
        // reopen element at the current position
        it->reopen(filePath.c_str(), fadvise);
    }

};

} // namespace io
} // namespace isaac

#endif // #ifndef iSAAC_IO_FILE_BUF_CACHE_HH
