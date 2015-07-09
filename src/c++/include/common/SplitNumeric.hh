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
 ** \file Kmer.hh
 **
 ** General definitions and tools for handling k-mers.
 **
 ** \author Come Raczy
 **/

#ifndef iSAAC_COMMON_SPLIT_NUMERIC_HH
#define iSAAC_COMMON_SPLIT_NUMERIC_HH

#include <bitset>
#include <string>

#include "oligo/Nucleotides.hh"

namespace isaac
{
namespace common
{


template <typename T> std::ostream &traceHexValue(
    std::ostream &os, const T &kmer)
{
    const std::ios_base::fmtflags ff = os.flags();
    os << std::hex << std::setfill('0') << std::setw(sizeof(kmer) * 2) << (unsigned long)kmer;
    os.flags(ff);
    return os;
}

template <> inline std::ostream &traceHexValue<__uint128_t>(
    std::ostream &os, const __uint128_t &kmer)
{
    const unsigned long lobits(kmer);
    const unsigned long hibits(kmer >> 64);
    const std::ios_base::fmtflags ff = os.flags();
    os << std::hex << std::setw(16) << std::setfill('0') << hibits << "-" << std::setw(16) << std::setfill('0') << lobits;
    os.flags(ff);
    return os;
}

template <typename HiT>
struct SplitNumeric
{

    friend std::ostream &operator << (std::ostream& os, const SplitNumeric &split)
    {
        const std::ios_base::fmtflags ff = os.flags();
        traceHexValue(os, split.hi_);
        os << "-";
        traceHexValue(os, split.lo_);
        os.flags(ff);
        return os;
    }

    typedef __uint128_t LoT;
    static const unsigned BITS_PER_BYTE = 8;
    static const unsigned HI_BITS = sizeof(HiT) * BITS_PER_BYTE;
    static const unsigned LO_BITS = sizeof(LoT) * BITS_PER_BYTE;
    static const unsigned TOTAL_BITS = LO_BITS + HI_BITS;

    SplitNumeric(LoT lo) : lo_(lo), hi_(0U){;}
    SplitNumeric(HiT hi, LoT lo) : lo_(lo), hi_(hi){;}
    template <typename TargetHiT>
    explicit SplitNumeric(const SplitNumeric<TargetHiT> &that) : lo_(that.lo_), hi_(that.hi_){;}
    SplitNumeric& operator <<=(const std::size_t lshift)
    {
        ISAAC_ASSERT_MSG(lshift < TOTAL_BITS, "Refusing to shift " << TOTAL_BITS << "-bit data type by " << lshift << " bits");
        if (!lshift)
        {
            return *this;
        }

        hi_ = lshift >= HI_BITS ? HiT(0) : (hi_ << lshift);
        if (lshift >= LO_BITS)
        {
            hi_ |= HiT(lo_) << (lshift - LO_BITS);
            lo_ = LoT(0);
        }
        else
        {
            hi_ |= HiT(lo_ >> (LO_BITS - lshift));
            lo_ <<= lshift;
        }
        return *this;
    }

    SplitNumeric& operator >>=(const std::size_t rshift)
    {
        ISAAC_ASSERT_MSG(rshift < TOTAL_BITS, "Refusing to shift " << TOTAL_BITS << "-bit data type by " << rshift << " bits");

        if (!rshift)
        {
            return *this;
        }

        if (HI_BITS > rshift)
        {
//            const HiT carryMask = ~HiT(0) >> (HI_BITS - rshift);
            const HiT carry = (hi_ << (HI_BITS - rshift)) >> (HI_BITS - rshift);
            if (LO_BITS > rshift)
            {
                lo_ >>= rshift;
                lo_ |= LoT(carry) << (LO_BITS - rshift);
            }
            else
            {
                lo_ = LoT(carry >> (rshift - LO_BITS));
            }
            hi_ >>= rshift;
        }
        else if (LO_BITS > rshift)
        {
            lo_ >>= rshift;
            lo_ |= LoT(hi_) << (LO_BITS - rshift);
            hi_  = HiT(0);
        }
        else
        {
            lo_ = LoT(hi_ >> (rshift - LO_BITS));
            hi_  = HiT(0);
        }

        return *this;
    }

    template <typename OtherHiT>
    SplitNumeric& operator =(const SplitNumeric<OtherHiT> &that)
    {
        lo_ = that.lo_;
        hi_ = that.hi_;
        return *this;
    }

    SplitNumeric& operator |=(const SplitNumeric &that)
    {
        lo_ |= that.lo_;
        hi_ |= that.hi_;
        return *this;
    }

    SplitNumeric& operator &=(const SplitNumeric &that)
    {
        lo_ &= that.lo_;
        hi_ &= that.hi_;
        return *this;
    }

    SplitNumeric& operator ^=(const SplitNumeric &that)
    {
        lo_ ^= that.lo_;
        hi_ ^= that.hi_;
        return *this;
    }

    bool operator <(const SplitNumeric &that) const
    {
        return hi_ < that.hi_ || (hi_ == that.hi_ && lo_ < that.lo_);
    }

    bool operator !=(const SplitNumeric &that) const
    {
        return lo_ != that.lo_ || hi_ != that.hi_;
    }

    bool operator ==(const SplitNumeric &that) const
    {
        return lo_ == that.lo_ && hi_ == that.hi_;
    }

    unsigned operator &(const unsigned& right) const {return lo_ & right;}

    bool operator !() const {return !lo_ && !hi_;}

    SplitNumeric operator ~() const
    {
        SplitNumeric ret(~hi_, ~lo_);
        return ret;
    }

    SplitNumeric operator &(const SplitNumeric& right) const {SplitNumeric ret(*this); ret &= right; return ret;}
    template <typename ShifT> SplitNumeric operator >>(const ShifT &rshift) const {SplitNumeric ret(*this); ret >>= rshift; return ret;}
    template <typename ShifT> SplitNumeric operator <<(const ShifT &lshift) const {SplitNumeric ret(*this); ret <<= lshift; return ret;}

    operator LoT()const
    {
        return lo_;
    }

private:
    template <typename OtherHiT> friend class SplitNumeric;
    LoT lo_;
    HiT hi_;
}__attribute__ ((packed));

template <typename HiT>
bool operator == (const SplitNumeric<HiT> &left, const __uint128_t right)
{
    return left == SplitNumeric<HiT>(right);
}

template <typename HiT>
bool operator == (const __uint128_t left, const SplitNumeric<HiT> &right)
{
    return SplitNumeric<HiT>(left) == right;
}

template <typename HiT>
bool operator < (const SplitNumeric<HiT> &left, const __uint128_t right)
{
    return left < SplitNumeric<HiT>(right);
}

template <typename HiT>
bool operator < (const __uint128_t left, const SplitNumeric<HiT> &right)
{
    return SplitNumeric<HiT>(left) < right;
}

template <typename HiT>
inline std::ostream &traceHexValue(
    std::ostream &os, const SplitNumeric<HiT> &kmer)
{
    return os << kmer;
}


typedef SplitNumeric<unsigned char> Uint136;
typedef SplitNumeric<unsigned short> Uint144;
typedef SplitNumeric<unsigned> Uint160;
typedef SplitNumeric<unsigned long> Uint192;
typedef SplitNumeric<__uint128_t> Uint256;

} // namespace common
} // namespace isaac

#endif // #ifndef iSAAC_COMMON_SPLIT_NUMERIC_HH
