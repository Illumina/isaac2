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
 ** \file Kmer.hh
 **
 ** General definitions and tools for handling k-mers.
 **
 ** \author Come Raczy
 **/

#ifndef iSAAC_OLIGO_K_MER_HH
#define iSAAC_OLIGO_K_MER_HH

//#include <bitset>
#include <string>

#include <boost/foreach.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/mpl/front.hpp>
#include <boost/mpl/if.hpp>
#include <boost/mpl/modulus.hpp>
#include <boost/mpl/vector/vector10_c.hpp>


#include "common/SplitNumeric.hh"
#include "oligo/Nucleotides.hh"

namespace isaac
{
namespace oligo
{

static const unsigned BITS_PER_BASE = 2;
static const unsigned BITS_PER_BASE_MASK = 3;


template <unsigned K>
struct KmerBitsType
{
    typedef typename boost::uint_t<K*BITS_PER_BASE>::least BitsType;
};
typedef boost::mpl::vector1_c<int, 4> SUPPORTED_KMERS_4;
typedef boost::mpl::push_back<SUPPORTED_KMERS_4, boost::mpl::int_<6> >::type SUPPORTED_KMERS_6;
typedef boost::mpl::push_back<SUPPORTED_KMERS_6, boost::mpl::int_<8> >::type SUPPORTED_KMERS_8;
typedef boost::mpl::push_back<SUPPORTED_KMERS_8, boost::mpl::int_<10> >::type SUPPORTED_KMERS_10;
typedef boost::mpl::push_back<SUPPORTED_KMERS_10, boost::mpl::int_<12> >::type SUPPORTED_KMERS_12;
typedef boost::mpl::push_back<SUPPORTED_KMERS_12, boost::mpl::int_<14> >::type SUPPORTED_KMERS_14;
typedef boost::mpl::push_back<SUPPORTED_KMERS_14, boost::mpl::int_<16> >::type SUPPORTED_KMERS_16;
typedef boost::mpl::push_back<SUPPORTED_KMERS_16, boost::mpl::int_<18> >::type SUPPORTED_KMERS_18;
typedef boost::mpl::push_back<SUPPORTED_KMERS_18, boost::mpl::int_<20> >::type SUPPORTED_KMERS_20;
typedef boost::mpl::push_back<SUPPORTED_KMERS_20, boost::mpl::int_<22> >::type SUPPORTED_KMERS_22;
typedef boost::mpl::push_back<SUPPORTED_KMERS_22, boost::mpl::int_<24> >::type SUPPORTED_KMERS_24;
typedef boost::mpl::push_back<SUPPORTED_KMERS_24, boost::mpl::int_<26> >::type SUPPORTED_KMERS_26;
typedef boost::mpl::push_back<SUPPORTED_KMERS_26, boost::mpl::int_<28> >::type SUPPORTED_KMERS_28;
typedef boost::mpl::push_back<SUPPORTED_KMERS_28, boost::mpl::int_<30> >::type SUPPORTED_KMERS_30;
typedef boost::mpl::push_back<SUPPORTED_KMERS_30, boost::mpl::int_<32> >::type SUPPORTED_KMERS_32;

template <> struct KmerBitsType<34>{typedef __uint128_t BitsType;};
typedef boost::mpl::push_back<SUPPORTED_KMERS_32, boost::mpl::int_<34> >::type SUPPORTED_KMERS_34;
template <> struct KmerBitsType<36>{typedef __uint128_t BitsType;};
typedef boost::mpl::push_back<SUPPORTED_KMERS_34, boost::mpl::int_<36> >::type SUPPORTED_KMERS_36;
template <> struct KmerBitsType<38>{typedef __uint128_t BitsType;};
typedef boost::mpl::push_back<SUPPORTED_KMERS_36, boost::mpl::int_<38> >::type SUPPORTED_KMERS_38;
template <> struct KmerBitsType<40>{typedef __uint128_t BitsType;};
typedef boost::mpl::push_back<SUPPORTED_KMERS_38, boost::mpl::int_<40> >::type SUPPORTED_KMERS_40;
template <> struct KmerBitsType<44>{typedef __uint128_t BitsType;};
typedef boost::mpl::push_back<SUPPORTED_KMERS_40, boost::mpl::int_<44> >::type SUPPORTED_KMERS_44;
template <> struct KmerBitsType<48>{typedef __uint128_t BitsType;};
typedef boost::mpl::push_back<SUPPORTED_KMERS_44, boost::mpl::int_<48> >::type SUPPORTED_KMERS_48;
template <> struct KmerBitsType<52>{typedef __uint128_t BitsType;};
typedef boost::mpl::push_back<SUPPORTED_KMERS_48, boost::mpl::int_<52> >::type SUPPORTED_KMERS_52;
template <> struct KmerBitsType<56>{typedef __uint128_t BitsType;};
typedef boost::mpl::push_back<SUPPORTED_KMERS_52, boost::mpl::int_<56> >::type SUPPORTED_KMERS_56;
template <> struct KmerBitsType<60>{typedef __uint128_t BitsType;};
typedef boost::mpl::push_back<SUPPORTED_KMERS_56, boost::mpl::int_<60> >::type SUPPORTED_KMERS_60;
template <> struct KmerBitsType<64>{typedef __uint128_t BitsType;};
typedef boost::mpl::push_back<SUPPORTED_KMERS_60, boost::mpl::int_<64> >::type SUPPORTED_KMERS_64;
template <> struct KmerBitsType<68>{typedef common::Uint136 BitsType;};
typedef boost::mpl::push_back<SUPPORTED_KMERS_64, boost::mpl::int_<68> >::type SUPPORTED_KMERS_68;
template <> struct KmerBitsType<72>{typedef common::Uint144 BitsType;};
typedef boost::mpl::push_back<SUPPORTED_KMERS_68, boost::mpl::int_<72> >::type SUPPORTED_KMERS_72;
template <> struct KmerBitsType<76>{typedef common::Uint160 BitsType;};
typedef boost::mpl::push_back<SUPPORTED_KMERS_72, boost::mpl::int_<76> >::type SUPPORTED_KMERS_76;
template <> struct KmerBitsType<80>{typedef common::Uint160 BitsType;};
typedef boost::mpl::push_back<SUPPORTED_KMERS_76, boost::mpl::int_<80> >::type SUPPORTED_KMERS_80;
template <> struct KmerBitsType<128>{typedef common::Uint256 BitsType;};
typedef boost::mpl::push_back<SUPPORTED_KMERS_80, boost::mpl::int_<128> >::type SUPPORTED_KMERS_128;

typedef SUPPORTED_KMERS_128 SUPPORTED_KMERS;

template <typename BitsT, typename DerivedT>
class ArithmeticOperations
{
    DerivedT &derived() {return *static_cast<DerivedT*>(this);}
    const DerivedT &derived() const {return *static_cast<const DerivedT*>(this);}
public:
    typedef BitsT BitsType;

    template <typename ShiftT> DerivedT& operator <<=(const ShiftT &lshift) { derived().bits_ <<= lshift; derived().bits_ &= derived().BITS_MASK(); return derived();}
    template <typename ShiftT> DerivedT& operator >>=(const ShiftT &rshift) { derived().bits_ >>= rshift; return derived();}
    DerivedT& operator |=(const DerivedT &that) { derived().bits_ |= that.bits_; return derived();}
    DerivedT& operator &=(const DerivedT &that) { derived().bits_ &= that.bits_; return derived();}
    DerivedT& operator ^=(const DerivedT &that) { derived().bits_ ^= that.bits_; return derived();}

    template <typename ShiftT>
    DerivedT operator <<(const ShiftT &lshift) const {DerivedT ret(derived()); ret <<= lshift; return ret;}
    template <typename ShiftT>
    DerivedT operator >>(const ShiftT &rshift) const {DerivedT ret(derived()); ret >>= rshift; return ret;}
    DerivedT operator |(const DerivedT& right) const {DerivedT ret(derived()); ret |= right; return ret;}
    DerivedT operator &(const DerivedT& right) const {DerivedT ret(derived()); ret &= right; return ret;}

    DerivedT operator ^(const DerivedT& right) const {DerivedT ret(derived()); ret ^= right; return ret;}

    unsigned operator &(const unsigned& right) const {return derived().bits_ & right;}

    bool operator !() const { return !derived().bits_;}
    bool operator <(const DerivedT& right) const {return derived().bits_ < right.bits_;}
    bool operator >(const DerivedT& right) const {return right.bits_ < derived().bits_;}
    bool operator !=(const DerivedT& right) const {return derived().bits_ != right.bits_;}
    bool operator ==(const DerivedT& right) const {return !(derived().bits_ != right.bits_);}
    DerivedT operator ~() const { return DerivedT(~derived().bits_ & DerivedT::BITS_MASK());}
};


template <unsigned K>
struct BasicKmerType : public ArithmeticOperations<typename KmerBitsType<K>::BitsType, BasicKmerType<K> >
{
    typedef ArithmeticOperations<typename KmerBitsType<K>::BitsType, BasicKmerType<K> > BaseT;
    typedef typename BaseT::BitsType BitsType;
    BitsType bits_;

    static const unsigned BITS_PER_BYTE = 8;
    static const unsigned KMER_BASES = K;
    static const unsigned KMER_BITS = KMER_BASES * BITS_PER_BASE;

    static BitsType BITS_MASK(){return ~BitsType(0U) >> ((sizeof(BitsType) * BITS_PER_BYTE) - KMER_BITS);}

    explicit BasicKmerType(BitsType u) : bits_(u)
    {
        ISAAC_ASSERT_MSG((u & BITS_MASK()) == u, "Invalid initialization value " << *this << " supplied to " << KMER_BASES << "-mer");
    }

    template <unsigned N> BasicKmerType(const BasicKmerType<N>& that):
        bits_(that.bits_)
    {
    }

    //special-case operations allowed when K != N
    template <unsigned N> bool operator ==(const BasicKmerType<N>& right) const { return bits_ == right.bits_;}
    template <unsigned N> bool operator !=(const BasicKmerType<N>& right) const { return bits_ != right.bits_;}
    template <unsigned N> bool operator <(const BasicKmerType<N>& right) const { return bits_ < right.bits_;}
    template <unsigned N> bool operator >(const BasicKmerType<N>& right) const { return bits_ > right.bits_;}

}__attribute__ ((packed));

typedef BasicKmerType<64> LongKmerType;
BOOST_STATIC_ASSERT_MSG(16 == sizeof(LongKmerType), "Unexpected object type size");

typedef BasicKmerType<32> KmerType;
BOOST_STATIC_ASSERT_MSG(8 == sizeof(KmerType), "Unexpected object type size");

typedef BasicKmerType<16> ShortKmerType;
BOOST_STATIC_ASSERT_MSG(4 == sizeof(ShortKmerType), "Unexpected object type size");

typedef BasicKmerType<8> VeryShortKmerType;
BOOST_STATIC_ASSERT_MSG(2 == sizeof(VeryShortKmerType), "Unexpected object type size");

typedef BasicKmerType<4> TinyKmerType;
BOOST_STATIC_ASSERT_MSG(1 == sizeof(TinyKmerType), "Unexpected object type size");

template<typename KmerT, bool withSuffix>
struct KmerTraitsImpl
{
    static const unsigned BITS_PER_BYTE = 8;
    static const unsigned KMER_BASES = KmerT::KMER_BASES;
    static const unsigned KMER_BITS = KMER_BASES * BITS_PER_BASE;
    static const unsigned KMER_BYTES = (KMER_BITS + BITS_PER_BYTE - 1) / BITS_PER_BYTE;
    static const std::size_t UNUSED_BITS = KMER_BITS % BITS_PER_BYTE;

    BOOST_STATIC_ASSERT(!(KMER_BASES % 2));
    typedef BasicKmerType<KMER_BASES / 2> SuffixType;

    static const unsigned SUFFIX_BITS = KMER_BITS / 2;
};


template<typename KmerT>
struct KmerTraitsImpl<KmerT, false>
{
    static const unsigned BITS_PER_BYTE = 8;
    static const unsigned KMER_BASES = KmerT::KMER_BASES;
    static const unsigned KMER_BITS = KMER_BASES * BITS_PER_BASE;
    static const unsigned KMER_BYTES = (KMER_BITS + BITS_PER_BYTE - 1) / BITS_PER_BYTE;
    static const std::size_t UNUSED_BITS = KMER_BITS % BITS_PER_BYTE;
};

template<typename KmerT>
struct KmerTraits : public KmerTraitsImpl<KmerT, boost::mpl::equal_to<
    boost::mpl::modulus<boost::mpl::int_<KmerT::KMER_BASES>, boost::mpl::int_<2> >,
    boost::mpl::int_<0> >::value >
{
};


template <unsigned bytes>
struct IntegralKmerTraits
{
    static const unsigned BITS_PER_BYTE = 8;
    static const unsigned KMER_BITS = bytes * BITS_PER_BYTE;
    static const unsigned KMER_BYTES = bytes;
    static const unsigned KMER_BASES = KMER_BITS / BITS_PER_BASE;
    static const std::size_t UNUSED_BITS = 0;

    typedef typename boost::uint_t<KMER_BITS / 2>::exact SuffixType;
    static const unsigned SUFFIX_BITS = sizeof(SuffixType) * BITS_PER_BYTE;
};

template<>
struct KmerTraits<unsigned short> : public IntegralKmerTraits<sizeof(unsigned short)>
{
};

template<>
struct KmerTraits<unsigned int> : public IntegralKmerTraits<sizeof(unsigned int)>
{
};

/*
 * \brief shift left by specified amount of BASES. Also ensures that the result is not undefined in
 *        case the number of bases equals to KmerTraits<KmerT>::KMER_BASES. See below.
 *
 * "The count operand can be an immediate value or register CL.
 * The count is masked to five bits, which limits the count range to 0 to 31."
 * See http://www.intel.com/design/intarch/manuals/243191.htm
*/
template <typename KmerT>
KmerT shlBases(KmerT kmer, const std::size_t bases)
{
    ISAAC_ASSERT_MSG(bases <= KmerTraits<KmerT>::KMER_BASES, "Shifting more than max number of bases indicates an error in the code: " << bases);
    if (KmerTraits<KmerT>::KMER_BASES == bases)
    {
        return KmerT(0);
    }
    return kmer << (oligo::BITS_PER_BASE * bases);
}

/*
 * \brief shift left by specified amount of BASES. Also ensures that the result is not undefined in
 *        case the number of bases equals to KmerTraits<KmerT>::KMER_BASES. See below.
 *
 * "The count operand can be an immediate value or register CL.
 * The count is masked to five bits, which limits the count range to 0 to 31."
 * See http://www.intel.com/design/intarch/manuals/243191.htm
*/
template <typename KmerT>
KmerT safeShl(KmerT kmer, const std::size_t bits)
{
    ISAAC_ASSERT_MSG(bits <= KmerTraits<KmerT>::KMER_BITS, "Shifting more than max number of bits indicates an error in the code: " << bits << ">" << KmerTraits<KmerT>::KMER_BITS);
    if (KmerTraits<KmerT>::KMER_BITS == bits)
    {
        return KmerT(0);
    }
    return kmer << bits;
}


template <typename KmerT>
std::string bases(KmerT kmer)
{
    return bases<BITS_PER_BASE>(kmer, KmerTraits<KmerT>::KMER_BASES);
}

template <typename KmerT>
inline std::string reverseBases(KmerT kmer)
{
    std::string s = bases(~kmer);
    std::reverse(s.begin(), s.end());
    return s;
}

inline unsigned char rc(const unsigned char &forward)
{
    typedef unsigned char KmerT;
    KmerT kmer(forward);
    KmerT reversed(0);
    kmer = ~kmer; // complement all the bases
    for (unsigned i = 0; sizeof(KmerT) * 8 / BITS_PER_BASE > i; ++i)
    {
        reversed <<= BITS_PER_BASE;
        reversed |= (kmer & KmerT(BITS_PER_BASE_MASK));
        kmer >>= BITS_PER_BASE;
    }
    return reversed;
}

template <typename KmerT>
KmerT reverseComplement(KmerT kmer)
{
    static const unsigned char BYTE_REVERSE_COMPLEMENTS[] =
    {
     rc(0x00), rc(0x01), rc(0x02), rc(0x03), rc(0x04), rc(0x05), rc(0x06), rc(0x07), rc(0x08), rc(0x09), rc(0x0a), rc(0x0b), rc(0x0c), rc(0x0d), rc(0x0e), rc(0x0f),
     rc(0x10), rc(0x11), rc(0x12), rc(0x13), rc(0x14), rc(0x15), rc(0x16), rc(0x17), rc(0x18), rc(0x19), rc(0x1a), rc(0x1b), rc(0x1c), rc(0x1d), rc(0x1e), rc(0x1f),
     rc(0x20), rc(0x21), rc(0x22), rc(0x23), rc(0x24), rc(0x25), rc(0x26), rc(0x27), rc(0x28), rc(0x29), rc(0x2a), rc(0x2b), rc(0x2c), rc(0x2d), rc(0x2e), rc(0x2f),
     rc(0x30), rc(0x31), rc(0x32), rc(0x33), rc(0x34), rc(0x35), rc(0x36), rc(0x37), rc(0x38), rc(0x39), rc(0x3a), rc(0x3b), rc(0x3c), rc(0x3d), rc(0x3e), rc(0x3f),
     rc(0x40), rc(0x41), rc(0x42), rc(0x43), rc(0x44), rc(0x45), rc(0x46), rc(0x47), rc(0x48), rc(0x49), rc(0x4a), rc(0x4b), rc(0x4c), rc(0x4d), rc(0x4e), rc(0x4f),
     rc(0x50), rc(0x51), rc(0x52), rc(0x53), rc(0x54), rc(0x55), rc(0x56), rc(0x57), rc(0x58), rc(0x59), rc(0x5a), rc(0x5b), rc(0x5c), rc(0x5d), rc(0x5e), rc(0x5f),
     rc(0x60), rc(0x61), rc(0x62), rc(0x63), rc(0x64), rc(0x65), rc(0x66), rc(0x67), rc(0x68), rc(0x69), rc(0x6a), rc(0x6b), rc(0x6c), rc(0x6d), rc(0x6e), rc(0x6f),
     rc(0x70), rc(0x71), rc(0x72), rc(0x73), rc(0x74), rc(0x75), rc(0x76), rc(0x77), rc(0x78), rc(0x79), rc(0x7a), rc(0x7b), rc(0x7c), rc(0x7d), rc(0x7e), rc(0x7f),
     rc(0x80), rc(0x81), rc(0x82), rc(0x83), rc(0x84), rc(0x85), rc(0x86), rc(0x87), rc(0x88), rc(0x89), rc(0x8a), rc(0x8b), rc(0x8c), rc(0x8d), rc(0x8e), rc(0x8f),
     rc(0x90), rc(0x91), rc(0x92), rc(0x93), rc(0x94), rc(0x95), rc(0x96), rc(0x97), rc(0x98), rc(0x99), rc(0x9a), rc(0x9b), rc(0x9c), rc(0x9d), rc(0x9e), rc(0x9f),
     rc(0xa0), rc(0xa1), rc(0xa2), rc(0xa3), rc(0xa4), rc(0xa5), rc(0xa6), rc(0xa7), rc(0xa8), rc(0xa9), rc(0xaa), rc(0xab), rc(0xac), rc(0xad), rc(0xae), rc(0xaf),
     rc(0xb0), rc(0xb1), rc(0xb2), rc(0xb3), rc(0xb4), rc(0xb5), rc(0xb6), rc(0xb7), rc(0xb8), rc(0xb9), rc(0xba), rc(0xbb), rc(0xbc), rc(0xbd), rc(0xbe), rc(0xbf),
     rc(0xc0), rc(0xc1), rc(0xc2), rc(0xc3), rc(0xc4), rc(0xc5), rc(0xc6), rc(0xc7), rc(0xc8), rc(0xc9), rc(0xca), rc(0xcb), rc(0xcc), rc(0xcd), rc(0xce), rc(0xcf),
     rc(0xd0), rc(0xd1), rc(0xd2), rc(0xd3), rc(0xd4), rc(0xd5), rc(0xd6), rc(0xd7), rc(0xd8), rc(0xd9), rc(0xda), rc(0xdb), rc(0xdc), rc(0xdd), rc(0xde), rc(0xdf),
     rc(0xe0), rc(0xe1), rc(0xe2), rc(0xe3), rc(0xe4), rc(0xe5), rc(0xe6), rc(0xe7), rc(0xe8), rc(0xe9), rc(0xea), rc(0xeb), rc(0xec), rc(0xed), rc(0xee), rc(0xef),
     rc(0xf0), rc(0xf1), rc(0xf2), rc(0xf3), rc(0xf4), rc(0xf5), rc(0xf6), rc(0xf7), rc(0xf8), rc(0xf9), rc(0xfa), rc(0xfb), rc(0xfc), rc(0xfd), rc(0xfe), rc(0xff),
    };
    unsigned char * const begin = reinterpret_cast<unsigned char*>(&kmer);
    unsigned char * const end = begin + KmerTraits<KmerT>::KMER_BYTES;
    std::reverse(begin, end);

    BOOST_FOREACH(unsigned char &b, std::make_pair(begin, end))
    {
        b = BYTE_REVERSE_COMPLEMENTS[b];
    }

    static const std::size_t UNUSED_BITS = KmerTraits<KmerT>::UNUSED_BITS;
    if (UNUSED_BITS)
    {
        kmer >>= UNUSED_BITS;
    }

    return kmer;
}

template <typename KmerT>
inline std::ostream & operator <<(std::ostream &os, const oligo::Bases<BITS_PER_BASE, KmerT> &bases)
{
    return printBases(os, bases);
}

template <typename KmerT>
inline std::ostream & operator <<(std::ostream &os, const oligo::ReverseBases<BITS_PER_BASE, KmerT> &bases)
{
    return printReverseBases(os, bases);
}


template <unsigned K>
std::ostream & operator <<(std::ostream &os, const BasicKmerType<K>& kmer)
{
    return common::traceHexValue(os, kmer.bits_);
}

} // namespace oligo
} // namespace isaac

#endif // #ifndef iSAAC_OLIGO_K_MER_HH
