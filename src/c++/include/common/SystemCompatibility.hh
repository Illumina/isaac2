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
 ** \file SystemCompatibility.hh
 **
 ** \brief Interface layer for system-dependent functionalities
 ** 
 ** \author Come Raczy
 **/

#ifndef iSAAC_COMMON_SYSTEM_COMPATIBILITY_HH
#define iSAAC_COMMON_SYSTEM_COMPATIBILITY_HH

#include <boost/filesystem.hpp>

namespace isaac
{
namespace common
{

/// Maximum number of files that a process can have opened at the same time
unsigned int getMaxOpenFiles();

/// File size in bytes as returned by stat
unsigned long getFileSize(const char *filePath);

/// Determine the processor time
long clock();

/// Check if the architecture is little endian
bool isLittleEndian();

/// limit virtual memory size available to process. (equivalent of ulimit -v)
bool ulimitV(const unsigned long availableMemory);
/// retrieves the current ulimit -v
bool ulimitV(unsigned long *pLimit);

/**
 * \brief Sets a hook that monitors memory allocations.
 *
 * \param hook  Hook to set. If hook returns false, allocation will attempt to terminateWithCoreDump and fail.
 *              Cals to hook happen under the boost::unique_lock<boost::mutex>
 */
void hookMalloc(bool (*hook)(size_t size, const void *caller));

/**
 * \brief Removes block set by blockMalloc. Returns the number of allocations made since last blockMalloc call.
 *
 * \param hook  Used to ensure that the hook removed is the one that was previously set.
 */
unsigned unhookMalloc(bool (*hook)(size_t size, const void *caller));

/**
 * \brief Generate a core dump with a meaningful backtrace
 */
void terminateWithCoreDump();

/**
 * \brief Disables memory management optimizations that are detrimental to the access pattern used in high-performance
 *        parts of the product
 */
void configureMemoryManagement(
    const bool disableMultipleArenas,
    const bool disableFastbins);

boost::filesystem::path getModuleFileName();

/**
 * \brief calls linux-specific fallocate with FALLOC_FL_KEEP_SIZE to pre-allocate file on disk in a way that does
 *        not mess up reopening for apend. Notice that posix_fallocate does not do the job.
 */
int linuxFallocate(int fd, std::size_t offset, std::size_t len);

} // namespace common
} // namespace isaac

#endif // #ifndef iSAAC_COMMON_SYSTEM_COMPATIBILITY_HH
