#pragma once

// Arbitrary-precision unsigned integer used only for decimal <-> binary256
// conversions. Not constexpr; I/O is a runtime path.

#include "uint.hpp"

#include <algorithm>
#include <bit>
#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

namespace base256::detail {

class BigInt {
    std::vector<std::uint32_t> d; // base 2^32, little-endian, no leading zeros

    void normalize() {
        while (!d.empty() && d.back() == 0) d.pop_back();
    }

public:
    BigInt() = default;
    explicit BigInt(std::uint64_t x) {
        if (x == 0) return;
        d.push_back(static_cast<std::uint32_t>(x));
        if (x >> 32) d.push_back(static_cast<std::uint32_t>(x >> 32));
    }

    static BigInt from_u256(const u256& x) {
        BigInt r;
        for (int i = 3; i >= 0; --i) {
            r.shl(32);
            r.add_u32(static_cast<std::uint32_t>(x.d[static_cast<std::size_t>(i)] >> 32));
            r.shl(32);
            r.add_u32(static_cast<std::uint32_t>(x.d[static_cast<std::size_t>(i)]));
        }
        return r;
    }

    static BigInt from_dec(std::string_view digits) {
        BigInt r;
        for (char c : digits) {
            if (c < '0' || c > '9') continue;
            r.mul_u32(10);
            r.add_u32(static_cast<std::uint32_t>(c - '0'));
        }
        return r;
    }

    [[nodiscard]] bool is_zero() const noexcept { return d.empty(); }

    [[nodiscard]] int limb_count() const noexcept { return static_cast<int>(d.size()); }

    [[nodiscard]] int bit_length() const noexcept {
        if (d.empty()) return 0;
        const int last = static_cast<int>(d.size()) - 1;
        return last * 32 + (32 - std::countl_zero(d[static_cast<std::size_t>(last)]));
    }

    [[nodiscard]] bool bit(int i) const noexcept {
        if (i < 0) return false;
        const int limb = i / 32;
        if (limb >= static_cast<int>(d.size())) return false;
        return ((d[static_cast<std::size_t>(limb)] >> (i % 32)) & 1u) != 0;
    }

    void add_u32(std::uint32_t x) {
        if (x == 0) return;
        std::uint64_t c = x;
        for (std::size_t i = 0; i < d.size() && c; ++i) {
            c += d[i];
            d[i] = static_cast<std::uint32_t>(c);
            c >>= 32;
        }
        if (c) d.push_back(static_cast<std::uint32_t>(c));
    }

    void mul_u32(std::uint32_t x) {
        if (is_zero()) return;
        if (x == 0) {
            d.clear();
            return;
        }
        std::uint64_t c = 0;
        for (std::uint32_t& limb : d) {
            c += static_cast<std::uint64_t>(limb) * x;
            limb = static_cast<std::uint32_t>(c);
            c >>= 32;
        }
        if (c) d.push_back(static_cast<std::uint32_t>(c));
    }

    // Multiply by 5^n.
    void mul_pow5(unsigned n) {
        while (n >= 13) { // 5^13 = 1_220_703_125 < 2^32
            mul_u32(1220703125u);
            n -= 13;
        }
        std::uint32_t p = 1;
        while (n--) p *= 5;
        if (p != 1) mul_u32(p);
    }

    void shl(int bits) {
        if (bits <= 0 || is_zero()) return;
        const int limbs = bits / 32;
        const int b = bits % 32;
        if (limbs) d.insert(d.begin(), static_cast<std::size_t>(limbs), 0);
        if (b == 0) return;
        std::uint32_t c = 0;
        for (std::uint32_t& limb : d) {
            const std::uint64_t v = (static_cast<std::uint64_t>(limb) << b) | c;
            limb = static_cast<std::uint32_t>(v);
            c = static_cast<std::uint32_t>(v >> 32);
        }
        if (c) d.push_back(c);
    }

    void shr(int bits) {
        if (bits <= 0 || is_zero()) return;
        const int limbs = bits / 32;
        const int b = bits % 32;
        if (limbs >= static_cast<int>(d.size())) {
            d.clear();
            return;
        }
        d.erase(d.begin(), d.begin() + limbs);
        if (b == 0) {
            normalize();
            return;
        }
        std::uint32_t leftover = 0;
        for (int i = static_cast<int>(d.size()) - 1; i >= 0; --i) {
            const std::uint32_t cur = d[static_cast<std::size_t>(i)];
            d[static_cast<std::size_t>(i)] = (cur >> b) | leftover;
            leftover = cur << (32 - b);
        }
        normalize();
    }

    // Shift right by `bits`. Returns true if any 1-bit was discarded.
    bool shr_sticky(int bits) {
        if (bits <= 0 || is_zero()) return false;
        if (bits >= bit_length()) {
            d.clear();
            return true;
        }
        bool sticky = false;
        const int limbs = bits / 32;
        const int b = bits % 32;
        const int n = static_cast<int>(d.size());
        const int whole = limbs < n ? limbs : n;
        for (int i = 0; i < whole; ++i) {
            if (d[static_cast<std::size_t>(i)] != 0) {
                sticky = true;
                break;
            }
        }
        if (!sticky && b != 0 && limbs < n) {
            const std::uint32_t mask = (1u << b) - 1u;
            if ((d[static_cast<std::size_t>(limbs)] & mask) != 0) sticky = true;
        }
        shr(bits);
        return sticky;
    }

    // Divide by 10, return remainder 0-9.
    std::uint32_t divmod10() {
        if (is_zero()) return 0;
        std::uint64_t rem = 0;
        for (int i = static_cast<int>(d.size()) - 1; i >= 0; --i) {
            rem = (rem << 32) | d[static_cast<std::size_t>(i)];
            d[static_cast<std::size_t>(i)] = static_cast<std::uint32_t>(rem / 10);
            rem %= 10;
        }
        normalize();
        return static_cast<std::uint32_t>(rem);
    }

    // Divide by 5^n (n>=0). Returns true if the division was exact.
    bool div_pow5(unsigned n) {
        while (n--) {
            std::uint64_t rem = 0;
            for (int i = static_cast<int>(d.size()) - 1; i >= 0; --i) {
                rem = (rem << 32) | d[static_cast<std::size_t>(i)];
                d[static_cast<std::size_t>(i)] = static_cast<std::uint32_t>(rem / 5);
                rem %= 5;
            }
            normalize();
            if (rem != 0) return false;
        }
        return true;
    }

    // this := floor(this / divisor), returns remainder. Divisor is a BigInt.
    BigInt divmod(const BigInt& divisor) {
        if (divisor.is_zero()) return {};
        if (*this < divisor) {
            BigInt rem = *this;
            d.clear();
            return rem;
        }
        const int nbits = bit_length();
        BigInt quot;
        BigInt rem;
        for (int i = nbits - 1; i >= 0; --i) {
            rem.shl(1);
            if (bit(i)) rem.add_u32(1);
            if (!(rem < divisor)) {
                rem.sub(divisor);
                // set bit i of quot
                const int limb = i / 32;
                const int b = i % 32;
                if (limb >= static_cast<int>(quot.d.size())) {
                    quot.d.resize(static_cast<std::size_t>(limb) + 1, 0);
                }
                quot.d[static_cast<std::size_t>(limb)] |= (1u << b);
            }
        }
        quot.normalize();
        *this = std::move(quot);
        return rem;
    }

    void sub(const BigInt& b) {
        // assume *this >= b
        std::uint64_t br = 0;
        const std::size_t n = d.size();
        for (std::size_t i = 0; i < n; ++i) {
            const std::uint64_t bv = i < b.d.size() ? b.d[i] : 0;
            const std::int64_t t = static_cast<std::int64_t>(d[i]) - static_cast<std::int64_t>(bv) -
                                   static_cast<std::int64_t>(br);
            if (t < 0) {
                d[i] = static_cast<std::uint32_t>(t + (1LL << 32));
                br = 1;
            } else {
                d[i] = static_cast<std::uint32_t>(t);
                br = 0;
            }
        }
        normalize();
    }

    void add(const BigInt& b) {
        std::uint64_t c = 0;
        const std::size_t n = std::max(d.size(), b.d.size());
        d.resize(n, 0);
        for (std::size_t i = 0; i < n; ++i) {
            c += d[i];
            if (i < b.d.size()) c += b.d[i];
            d[i] = static_cast<std::uint32_t>(c);
            c >>= 32;
        }
        if (c) d.push_back(static_cast<std::uint32_t>(c));
    }

    [[nodiscard]] bool operator<(const BigInt& o) const noexcept {
        if (d.size() != o.d.size()) return d.size() < o.d.size();
        for (int i = static_cast<int>(d.size()) - 1; i >= 0; --i) {
            if (d[static_cast<std::size_t>(i)] != o.d[static_cast<std::size_t>(i)]) {
                return d[static_cast<std::size_t>(i)] < o.d[static_cast<std::size_t>(i)];
            }
        }
        return false;
    }

    [[nodiscard]] bool operator==(const BigInt& o) const noexcept { return d == o.d; }

    // Extract 256 bits starting at `hi` (inclusive) down to hi-255.
    // If hi < 255, the value is right-aligned (high bits zero).
    [[nodiscard]] u256 extract_u256(int hi) const noexcept {
        u256 r{};
        for (int i = 0; i < 256; ++i) {
            const int src = hi - 255 + i;
            if (src >= 0 && bit(src)) r.set_bit(i);
        }
        return r;
    }

    // Decimal string of the integer (no sign).
    [[nodiscard]] std::string to_dec() const {
        if (is_zero()) return "0";
        BigInt tmp = *this;
        std::string s;
        while (!tmp.is_zero()) {
            const std::uint32_t r = tmp.divmod10();
            s.push_back(static_cast<char>('0' + r));
        }
        std::reverse(s.begin(), s.end());
        return s;
    }
};

} // namespace base256::detail
