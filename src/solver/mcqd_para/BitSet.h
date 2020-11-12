/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 * 
 */

#ifndef BITSET_H
#define BITSET_H


// some algorithms are from: http://graphics.stanford.edu/~seander/bithacks.html#IntegerMinOrMax


#include <sstream>
#include <vector>
#include <bitset>
#include <algorithm>
#include <iostream>

template<unsigned int I>
struct static_log2 {
    enum {
        value = static_log2<(I >> 1)>::value+1
    };
};

template<>
struct static_log2<1> {
    enum {
        value = 0
    };
};


struct LogTableSetter {
    LogTableSetter(char* logTable256) {
        logTable256[0] = logTable256[1] = 0;
        for (int i = 2; i < 256; i++) 
            logTable256[i] = 1 + logTable256[i >> 1];
        logTable256[0] = -1; // if you want log(0) to return -1
    }
};

// log2 does not detect negative numbers and returns max(char) if zero is passed to it
template<class Int>
unsigned int log2(Int v) {
    static char logTable256[256];
    static LogTableSetter dummy(logTable256);

    register unsigned int shr = sizeof(Int) << 2;
    register unsigned int ofs = 0;
    register Int v2;

    while (shr >= 8) {
        v2 = v >> shr;
        if (v2) {
            ofs += shr;
            v = v2;
        }
        shr >>= 1;
    }
    return ofs + logTable256[v];
}

// reverses the order of bits
template<class Int>
Int reverseBits(Int v) {
    Int s = sizeof(v) << 3; // bit size
    Int mask = ~0;         
    while ((s >>= 1) > 0) {
        mask ^= (mask << s);
        v = ((v >> s) & mask) | ((v << s) & ~mask);
    }
    return v;
}

// find first bit that equals to 1 (only designed for 32 and 64 bit integers) that works on all machines and compilers
template<class Int>
unsigned int countTrailing0M(Int v) {
    if (sizeof(Int) < 8) {
        // smaller and possibly faster version for 32 and les bit numbers
        static const char multiplyDeBruijnBitPosition[32] = {
            0, 1, 28, 2, 29, 14, 24, 3, 30, 22, 20, 15, 25, 17, 4, 8, 
            31, 27, 13, 23, 21, 19, 16, 7, 26, 12, 18, 6, 11, 5, 10, 9
        };
        return multiplyDeBruijnBitPosition[((uint32_t)((v & -v) * 0x077CB531U)) >> 27];
    } else {
        static const char multiplyDeBruijnBitPosition[64] = {
            0, 1, 2, 56, 3, 32, 57, 46, 29, 4, 20, 33, 7, 58, 11, 47, 
            62, 30, 18, 5, 16, 21, 34, 23, 53, 8, 59, 36, 25, 12, 48, 39, 
            63, 55, 31, 45, 28, 19, 6, 10, 61, 17, 15, 22, 52, 35, 24, 38, 
            54, 44, 27, 9, 60, 14, 51, 37, 43, 26, 13, 50, 42, 49, 41, 40
        };
        return multiplyDeBruijnBitPosition[((uint64_t)((v & -v) * 0x26752B916FC7B0DULL)) >> 58];
    };
}

// CTZ and POPCNT instrinsics for various platforms

#if __cpp_lib_bitops >= 201907L
#include <bit>
#include <type_traits>

template<typename Int>
unsigned int countTrailing0(Int v) {
    typedef std::make_unsigned_t<Int> UInt;
    return std::countr_zero<UInt>(v);
}

template<class Int>
unsigned int countOnes(Int v) {
    typedef std::make_unsigned_t<Int> UInt;
    return std::popcount<UInt>(v);
}
#else

#ifdef _MSC_VER
#include <intrin.h>

template<class Int>
unsigned int countTrailing0(Int v) {
    unsigned long trailing_zero = 0;

    if (_BitScanForward64(&trailing_zero, value))
    {
        return trailing_zero;
    }
    else
    {
        return 64;
    }
}

#define __builtin_popcount __popcnt
#define __builtin_popcountl __popcnt64
#define __builtin_popcountll __popcnt64

#else
// find first bit that equals to 1 (GCC only, undefined return if the operand equals 0 - [often 0 or the maximum number of bits])
template<class Int>
unsigned int countTrailing0(Int v) {
    if (sizeof(Int) < sizeof(int)) { return countTrailing0<int>(v); }
    if (sizeof(Int) == sizeof(int)) return __builtin_ctz(v);
    if (sizeof(Int) == sizeof(long int)) return __builtin_ctzll(v);
    if (sizeof(Int) >= sizeof(long long int)) return __builtin_ctzll(v);
}
#endif

// count number of bits that equal 1 (GCC only)
template<class Int>
unsigned int countOnes(Int v) {
    if (sizeof(Int) < sizeof(int)) { return countOnes<int>((int)v); }
    if (sizeof(Int) == sizeof(int)) return __builtin_popcount(v);
    if (sizeof(Int) == sizeof(long int)) return __builtin_popcountl(v);
    if (sizeof(Int) >= sizeof(long long int)) return __builtin_popcountll(v);
}

#endif

class BitSet {
protected:
    // how large are bit blocks [number of bits]
    static const unsigned int res = sizeof(unsigned long) * 8;
    // factor for shr when converting offset in bits to offset in data
    static const unsigned int f_shr = static_log2<res>::value;
    // mask for converting absolute offset in bits to relative offset inside a data cell
    static const unsigned int shift_mask = res-1;
    
    typedef std::bitset<res> b64;
    std::vector<b64> data;
    size_t numUsed, numAllocated;
    
    class BoolProxy {
        BitSet* parent;
        size_t index;
    
    public:
        BoolProxy(BitSet* p, size_t i) : parent(p), index(i) {}
        bool operator= (bool b) {parent->setValue(index, b); return b;}
        operator bool() const {return parent->getValue(index);}
    };
        
public:
    //~BitSet() {std::cout << data.size() << " : " << numAllocated << " : " << numUsed << "\n";}
    BitSet() : numUsed(0), numAllocated(0) {}
    BitSet(size_t n) : numUsed(0), numAllocated(0) {resize(n);}
    BitSet(const BitSet& other) : numUsed(other.numUsed), numAllocated(other.numAllocated) {copy(other);}
    const BitSet& operator= (const BitSet& other) {
        numUsed = other.numUsed; 
        numAllocated = other.numAllocated;
        copy(other); 
        return *this;
    }
    
    void resize(size_t newSize) {data.resize((newSize+shift_mask) >> f_shr); numAllocated = data.size() << f_shr; numUsed = newSize;}
    void resize(size_t newSize, bool value) {data.resize((newSize+shift_mask) >> f_shr); numAllocated = data.size() << f_shr; numUsed = newSize;}
    void reserve(size_t newSize) {numAllocated = ((newSize+shift_mask) >> f_shr) << f_shr; data.reserve(numAllocated >> f_shr);}
    size_t size() const {return numUsed;}
    
    // size of the bit blocks
    static unsigned int resolution() {return res;}
    
    // individual bit set/get
    bool operator[](size_t index) const {return getValue(index);}
    BoolProxy operator[](size_t index)  {return BoolProxy(this, index);}
    bool operator= (bool b) {
        set(0, numUsed, b);
        return b;
    }
    
    BitSet operator~ () const {
        BitSet notB;
        notB.resize(size());
        for (size_t i = 0; i < data.size(); ++i)
            notB.data[i] = ~data[i];
        notB.data.back() &= ((1 << (numUsed & shift_mask)) - 1); // fix the last data item, which might not have all the bits used -> the unused bits must remain 0
        return notB;
    }
    
    bool isZero() const {
        unsigned long long int zero = 0;
        for (auto b : data)
            zero |= b.to_ullong();
        return zero == 0;
    }
    
    void clear() {reset(0, numUsed);}
    
    // multiple bit set
    void set(size_t index_from, size_t index_to, bool val) {
        val ? set(index_from, index_to) : reset(index_from, index_to);
    }
    
    // multiple bit set
    void set(size_t index_from, size_t index_to) {
        if (index_from == index_to) return;
        size_t data_i1 = index_from >> f_shr;
        size_t data_i2 = (index_to-1) >> f_shr;
        unsigned long int data_mask1 = -1l << (index_from & shift_mask);
        unsigned long int data_mask2 = (1l << ((index_to-1) & shift_mask)) - 1;
        if (data_i1 == data_i2) {
            data[data_i1] |= (data_mask1 & data_mask2);
        } else {
            data[data_i1] |= data_mask1;
            for (auto i = data_i1 + 1; i < data_i2; ++i)
                data[i] |= -1l;
            data[data_i2] |= data_mask2;
        }
    }
    
    // multiple bit reset
    void reset(size_t index_from, size_t index_to) {
        if (index_from == index_to) return;
        size_t data_i1 = index_from >> f_shr;
        size_t data_i2 = (index_to-1) >> f_shr;
        unsigned long int data_mask1 = -1l << (index_from & shift_mask);
        unsigned long int data_mask2 = (1l << ((index_to-1) & shift_mask)) - 1;
        if (data_i1 == data_i2) {
            data[data_i1] &= ~(data_mask1 & data_mask2);       
        } else {
            data[data_i1] &= ~data_mask1;
            for (auto i = data_i1 + 1; i < data_i2; ++i)
                data[i] = 0;
            data[data_i2] &= ~data_mask2;
        }
    }
    
    const BitSet& operator&= (const BitSet& other) {
        for (size_t i = 0; i < data.size(); ++i)
            data[i] &= other.data[i];
        return *this;
    }
    
    const BitSet& operator^= (const BitSet& other) {
        //std::cout << data.size() << ", " << other.data.size() << ";\n";
        for (size_t i = 0, im = data.size(); i < im; ++i)
            data[i] ^= other.data[i];
        return *this;
    }
    
    // number of bits set
    size_t count() const {
        size_t res = 0;
        for(const b64& b:data) res+=countOnes(b.to_ulong());
        return res;
        //return std::accumulate(data.begin(), data.end(), 0, [](size_t a, const b64& b){return a + countOnes(b.to_ulong());});
    }
    
    // position of the next set bit (=1)
    int nextSetBit() const {
        for (size_t i = 0; i < data.size(); ++i) {
            if (data[i] != 0) {
                auto ul = data[i].to_ulong();
                return i*res + countTrailing0(ul);
            }
        }
        return size();
    }
    
    std::string to_string() const {
        std::ostringstream s;
        for (size_t i = 0; i < size(); ++i)
            s << (getValue(i) ? '1' : '0');
        return s.str();
    }

protected:
    bool getValue(size_t index) const {
        if (index >= numUsed) {
        std::cout << "Error in BitSet.getValue:\n" <<
                    "requested " << index << ", holding only " << numUsed << std::endl;
        }
        return data[index >> f_shr][index & shift_mask];}
    void setValue(size_t index, bool b = true) {data[index >> f_shr][index & shift_mask] = b;}
    void copy(const BitSet& other) {
        if (&other != this) {
            resize(other.size());
            std::copy(other.data.cbegin(), other.data.cend(), data.begin());
        }
    }
};

template<class Out>
Out& operator<< (Out& out, const BitSet& b) {
    if (b.size() == 0) {
        out << "[/]";
    } else {
        out << "[";
        bool comma = false;
        for (size_t i = 0; i < b.size(); ++ i) {
            if (b[i]) {
                if (comma)
                    out << ",";
                comma = true;
                out << i;  
            }
        }
        out << "]";
    }
    return out;
}

#endif // BITSET_H
