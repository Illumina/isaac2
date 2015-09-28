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
 ** \file Debug.cpp
 **
 ** \brief Various debugging-related helpers
 **
 ** \author Roman Petrovski
 **/

#include "common/Debug.hh"
#include "common/SystemCompatibility.hh"

namespace isaac
{
namespace common
{


namespace detail
{

// global variable that allows turning off cerr output for things such as unit tests.
std::atomic_int CerrBlocker::cerrBlocked_(0);
boost::recursive_mutex CerrLocker::cerrMutex_;


static bool malloc_warning_hook(size_t size, const void *caller)
{
    ISAAC_THREAD_CERR << "WARNING: blocked allocation of " << size << " bytes requested.\n";
    return true;
}

static bool malloc_strict_hook(size_t size, const void *caller)
{
    ISAAC_THREAD_CERR << "ERROR: blocked allocation of " << size << " bytes requested.\n";
    return false;
}

//static bool malloc_ignore_hook(size_t size, const void *caller)
//{
//    return true;
//}

} // namespace detail

ScopedMallocBlock::ScopedMallocBlock(const ScopedMallocBlock::Mode mode) :
    mode_(mode)
{
    block();
}

void ScopedMallocBlock::block()
{
    switch(mode_)
    {
    case Off:
//        hookMalloc(detail::malloc_ignore_hook);
        break;
    case Warning:
        hookMalloc(detail::malloc_warning_hook);
        break;
    case Strict:
        hookMalloc(detail::malloc_strict_hook);
        break;
    default:
        ISAAC_ASSERT_MSG(false, "invalid malloc block mode specified");
        break;
    }
}

ScopedMallocBlock::~ScopedMallocBlock()
{
    unblock();
}

void ScopedMallocBlock::unblock()
{
    switch(mode_)
    {
    case Off:
//        unhookMalloc(detail::malloc_ignore_hook);
        break;
    case Warning:
        unhookMalloc(detail::malloc_warning_hook);
        break;
    case Strict:
        unhookMalloc(detail::malloc_strict_hook);
        break;
    default:
        ISAAC_ASSERT_MSG(false, "invalid malloc block mode specified");
        break;
    }
}

ScopedMallocBlockUnblock::ScopedMallocBlockUnblock(ScopedMallocBlock &block) :
        block_(block)
{
    block_.unblock();
}

ScopedMallocBlockUnblock::~ScopedMallocBlockUnblock()
{
    block_.block();
}


} // namespace common
} // namespace isaac
