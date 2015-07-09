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
 ** \file FiniteCapacityVector.hh
 **
 ** Something that behaves more or less like std::vector but does not use dynamic memory
 ** and thus has fixed capacity().
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_COMMON_STATIC_VECTOR_HH
#define iSAAC_COMMON_STATIC_VECTOR_HH

#include <boost/numeric/ublas/storage.hpp>

namespace isaac
{
namespace common
{

// Bounded array - with allocator for size_type and difference_type
template<class T, std::size_t N, class ALLOC = std::allocator<T> >
class StaticVector: boost::numeric::ublas::bounded_array<T, N, ALLOC>
{
    typedef boost::numeric::ublas::bounded_array<T, N, ALLOC> BaseType;
public:
    using BaseType::begin;
    using BaseType::empty;
    using BaseType::iterator;
    using BaseType::const_iterator;
    using BaseType::const_reference;
    using BaseType::end;
    using BaseType::operator [];
    using BaseType::resize;
    using BaseType::size;
    using BaseType::value_type;

    StaticVector()
    {

    }

    StaticVector(std::size_t s, const T& v)
    {
        resize(s, v);
    }

    void clear()
    {
        BaseType::resize(0);
    }

    static std::size_t capacity() {return N;}

    bool full() const
    {
        return size() == capacity();
    }

    void push_back(const T& x)
    {
        BaseType::resize(BaseType::size() + 1, x);
    }

    T pop_back()
    {
        const T ret = back();
        BaseType::resize(BaseType::size() - 1);
        return ret;
    }

    T& front()
    {
        BOOST_UBLAS_CHECK (0 < size(), boost::numeric::ublas::bad_index ());
        return *begin();
    }

    const T& front() const
    {
        BOOST_UBLAS_CHECK (0 < BaseType::size(), boost::numeric::ublas::bad_index ());
        return *begin();
    }

    T& back()
    {
        BOOST_UBLAS_CHECK (0 < size(), boost::numeric::ublas::bad_index ());
        return *BaseType::rbegin();
    }

    const T& back() const
    {
        BOOST_UBLAS_CHECK (0 < BaseType::size(), boost::numeric::ublas::bad_index ());
        return *BaseType::rbegin();
    }

    void erase(typename BaseType::iterator b, typename BaseType::iterator e)
    {
        resize(std::distance(begin(), std::copy(e, end(), b)));
    }

    T& at(std::size_t index)
    {
        return this->operator[](index);
    }

    const T& at(std::size_t index)const
    {
        return this->operator[](index);
    }
};

} //namespace common
} //namespace isaac

#endif // #ifndef iSAAC_COMMON_STATIC_VECTOR_HH
