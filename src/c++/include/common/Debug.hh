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
 ** \file Debug.hh
 **
 ** \brief Various debugging-related helpers
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_LOG_THREAD_TIMESTAMP_HH
#define iSAAC_LOG_THREAD_TIMESTAMP_HH

#include <atomic>
#include <memory>
#include <typeinfo>

#include <boost/algorithm/string.hpp>
#include <boost/date_time.hpp>
#include <boost/io/ios_state.hpp>
#include <boost/iostreams/device/null.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/thread.hpp>

#include "common/SystemCompatibility.hh"

namespace isaac
{
namespace common
{

//TODO: check why simple CerrLocker(std::cerr) << ... << is not good enough
/**
 * \brief helper macro to simplify the thread-guarded logging. All elements on a single << line are serialized
 * under one CerrLocker
 */
#define ISAAC_THREAD_CERR \
    if(const ::isaac::common::detail::CerrLocker &isaac_cerr_lock = ::isaac::common::detail::CerrLocker()); \
    else (isaac_cerr_lock.cerrBlocked() ? \
        boost::iostreams::filtering_ostream(boost::iostreams::null_sink()) : \
            boost::iostreams::filtering_ostream(std::cerr)) << isaac::common::detail::ThreadTimestamp()

#define ISAAC_ASSERT_CERR \
    if(::isaac::common::detail::CerrLocker isaac_cerr_lock = ::isaac::common::detail::CerrLocker()); \
    else std::cerr << isaac::common::detail::ThreadTimestamp()


#define ISAAC_SCOPE_BLOCK_CERR if (const ::isaac::common::detail::CerrBlocker blocker = ::isaac::common::detail::CerrBlocker()); else

/**
 * \brief Evaluates expression always (even if NDEBUG is set and so on). Also uses ostream serialization which,
 *        unlike the standard assert, has shown not to allocate the dynamic memory at the time when you least
 *        expect this to happen.
 */
#define ISAAC_ASSERT_MSG(expr, msg) {if (expr) {} else \
{ ISAAC_ASSERT_CERR << "ERROR: ***** Internal Program Error - assertion (" << #expr << ") failed in " \
    << (BOOST_CURRENT_FUNCTION) << ":" << __FILE__ << '(' << __LINE__ << "): " << msg << std::endl; \
    ::isaac::common::terminateWithCoreDump();}}

inline std::string parseStat(const std::string &stat)
{
    std::vector<std::string> statFields;
    boost::algorithm::split(statFields, stat,  boost::algorithm::is_any_of(" "));
    return std::string(statFields.at(22) + "vm " + statFields.at(23) + "res");
}

#define ISAAC_TRACE_STAT(prefix) {\
    std::string statm; std::ifstream ifs("/proc/self/stat"); \
    std::getline(ifs, statm); \
    ISAAC_THREAD_CERR << "STAT: " << prefix << isaac::common::parseStat(statm) << std::endl;\
    }

class ScopedMallocBlock : boost::noncopyable
{
public:
    enum Mode
    {
        Invalid = 0,
        Off,
        Warning,
        Strict
    };

    ScopedMallocBlock(const Mode mode);
    ~ScopedMallocBlock();
private:
    const Mode mode_;

    friend class ScopedMallocBlockUnblock;
    void block();
    void unblock();
};

class ScopedMallocBlockUnblock : boost::noncopyable
{
    ScopedMallocBlock &block_;
public:
    ScopedMallocBlockUnblock(ScopedMallocBlock & block);
    ~ScopedMallocBlockUnblock();
};

namespace detail
{
class ThreadTimestamp
{
public:
};

/**
 * \brief formats time stamp and thread id to simplify threaded logging
 */
inline std::ostream & operator << (std::ostream &os, const ThreadTimestamp &) {

    // IMPORTANT: this is the way to serialize date without causing any dynamic memory operations to occur
    ::std::time_t t;
    ::std::time(&t);
    ::std::tm curr, *curr_ptr;
    curr_ptr = boost::date_time::c_time::localtime(&t, &curr);

    os << (curr_ptr->tm_year + 1900) << '-' <<
        std::setfill('0') << std::setw(2) << curr_ptr->tm_mon << '-'  <<
        std::setfill('0') << std::setw(2) << curr_ptr->tm_mday << ' '  <<

        std::setfill('0') << std::setw(2) << curr_ptr->tm_hour << ':' <<
        std::setfill('0') << std::setw(2) << curr_ptr->tm_min << ':' <<
        std::setfill('0') << std::setw(2) << curr_ptr->tm_sec << ' ' <<
        "\t[" << boost::this_thread::get_id() << "]\t";
    return os;
}

/**
 * \brief Blocks ISAAC_THREAD_CERR messages. Use for unit tests
 */
class CerrBlocker
{
    static std::atomic_int cerrBlocked_;
public:
    CerrBlocker()
    {
        ++cerrBlocked_;
    }

    ~CerrBlocker();

    operator bool () const {
        return false;
    }

    static bool blocked() {return cerrBlocked_;}
};

/**
 * \brief Guards std::cerr for the duration of CerrLocker existance
 *        Restores any changes made to ios::base
 */
class CerrLocker
{
    // some people allocate memory from under their trace code. For example by using boost::format.
    // if memory control is on, we don't want them to be dead-locked on their own thread cerrMutex_.
    static boost::recursive_mutex cerrMutex_;
    boost::lock_guard<boost::recursive_mutex> lock_;
    boost::io::ios_base_all_saver ias_;

public:

    CerrLocker(const CerrLocker &that) : lock_(cerrMutex_), ias_(std::cerr){
    }
    CerrLocker() : lock_(cerrMutex_), ias_(std::cerr) {
    }
    operator bool () const {
        return false;
    }

    bool cerrBlocked() const {return CerrBlocker::blocked();}
};

inline CerrBlocker::~CerrBlocker()
{
    ISAAC_ASSERT_MSG(cerrBlocked_, "Attempt to unblock more times than blocked. something is really wrong");
    --cerrBlocked_;
}


inline void assertion_failed_msg(char const * expr, char const * msg, char const * function,
                                 char const * file, long line)
{
    ISAAC_ASSERT_CERR
    << "ERROR: ***** Internal Program Error - assertion (" << expr << ") failed in "
    << function << ":" << file << '(' << line << "): " << msg << std::endl;

    common::terminateWithCoreDump();
}

} // namespace detail

/**
 ** \brief Provide a mechanism for detailed level of debugging
 **/
#ifdef ISAAC_THREAD_CERR_DEV_TRACE_ENABLED
    #define ISAAC_THREAD_CERR_DEV_TRACE(trace) {ISAAC_THREAD_CERR << trace << std::endl;}
    #define ISAAC_DEV_TRACE_BLOCK(block) block
    #define ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(clusterId, trace) ISAAC_THREAD_CERR_DEV_TRACE(trace)
#else
    #define ISAAC_THREAD_CERR_DEV_TRACE(blah)
    #ifdef ISAAC_THREAD_CERR_DEV_TRACE_ENABLED_CLUSTER_ID
        #define ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(clusterId, trace) {if(ISAAC_THREAD_CERR_DEV_TRACE_ENABLED_CLUSTER_ID == (clusterId)) {ISAAC_THREAD_CERR << trace << std::endl;}}
        #define ISAAC_DEV_TRACE_BLOCK(block) block
    #else
        #define ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(clusterId, blah)
        #define ISAAC_DEV_TRACE_BLOCK(block)
    #endif
#endif

inline std::time_t time()
{
    std::time_t ret;
    ISAAC_ASSERT_MSG(-1 != ::std::time(&ret), "std::time failed, errno: " << errno << strerror(errno));
    return ret;
}

} // namespace common
} // namespace isaac

#endif // #ifndef iSAAC_LOG_THREAD_TIMESTAMP_HH
