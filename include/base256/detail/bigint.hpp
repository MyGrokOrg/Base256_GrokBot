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

    // this := floor(this / divisor), returns remainder. Knuth, TAOCP 4.3.1
    // Algorithm D on 32-bit limbs: O(len(this) * len(divisor)).
    BigInt divmod(const BigInt& divisor) {
        if (divisor.is_zero()) return {};
        if (*this < divisor) {
            BigInt rem = *this;
            d.clear();
            return rem;
        }
        const std::size_t n = divisor.d.size();
        if (n == 1) {
            const std::uint64_t v = divisor.d[0];
            std::uint64_t rem = 0;
            for (std::size_t i = d.size(); i-- > 0;) {
                const std::uint64_t cur = (rem << 32) | d[i];
                d[i] = static_cast<std::uint32_t>(cur / v);
                rem = cur % v;
            }
            normalize();
            return BigInt(rem);
        }
        const std::size_t m = d.size() - n;
        const int sh = std::countl_zero(divisor.d.back());
        // Normalize: shift so the divisor's top limb has its high bit set.
        std::vector<std::uint32_t> v(n), u(d.size() + 1, 0);
        for (std::size_t i = n; i-- > 0;) {
            std::uint64_t x = static_cast<std::uint64_t>(divisor.d[i]) << sh;
            if (sh && i > 0) x |= static_cast<std::uint64_t>(divisor.d[i - 1]) >> (32 - sh);
            v[i] = static_cast<std::uint32_t>(x);
        }
        u[d.size()] = sh ? static_cast<std::uint32_t>(static_cast<std::uint64_t>(d.back()) >> (32 - sh)) : 0;
        for (std::size_t i = d.size(); i-- > 0;) {
            std::uint64_t x = static_cast<std::uint64_t>(d[i]) << sh;
            if (sh && i > 0) x |= static_cast<std::uint64_t>(d[i - 1]) >> (32 - sh);
            u[i] = static_cast<std::uint32_t>(x);
        }
        std::vector<std::uint32_t> q(m + 1, 0);
        const std::uint64_t B = 1ULL << 32;
        for (std::size_t j = m + 1; j-- > 0;) {
            const std::uint64_t num = (static_cast<std::uint64_t>(u[j + n]) << 32) | u[j + n - 1];
            std::uint64_t qhat = num / v[n - 1];
            std::uint64_t rhat = num % v[n - 1];
            while (qhat >= B || qhat * v[n - 2] > ((rhat << 32) | u[j + n - 2])) {
                --qhat;
                rhat += v[n - 1];
                if (rhat >= B) break;
            }
            // u[j..j+n] -= qhat * v
            std::int64_t borrow = 0;
            std::uint64_t carry = 0;
            for (std::size_t i = 0; i < n; ++i) {
                const std::uint64_t p = qhat * v[i] + carry;
                carry = p >> 32;
                const std::int64_t t = static_cast<std::int64_t>(u[i + j]) - borrow -
                                       static_cast<std::int64_t>(p & 0xFFFFFFFFULL);
                u[i + j] = static_cast<std::uint32_t>(t);
                borrow = t < 0 ? 1 : 0;
            }
            const std::int64_t t = static_cast<std::int64_t>(u[j + n]) - borrow -
                                   static_cast<std::int64_t>(carry);
            u[j + n] = static_cast<std::uint32_t>(t);
            if (t < 0) {
                // add back
                --qhat;
                std::uint64_t c = 0;
                for (std::size_t i = 0; i < n; ++i) {
                    const std::uint64_t s2 = static_cast<std::uint64_t>(u[i + j]) + v[i] + c;
                    u[i + j] = static_cast<std::uint32_t>(s2);
                    c = s2 >> 32;
                }
                u[j + n] = static_cast<std::uint32_t>(static_cast<std::uint64_t>(u[j + n]) + c);
            }
            q[j] = static_cast<std::uint32_t>(qhat);
        }
        // Remainder: u[0..n-1] >> sh
        BigInt rem;
        rem.d.resize(n);
        for (std::size_t i = 0; i < n; ++i) {
            std::uint64_t x = static_cast<std::uint64_t>(u[i]) >> sh;
            if (sh) x |= (static_cast<std::uint64_t>(u[i + 1]) << (32 - sh)) & 0xFFFFFFFFULL;
            rem.d[i] = static_cast<std::uint32_t>(x);
        }
        rem.normalize();
        d = std::move(q);
        normalize();
        return rem;
    }

    // this *= b (schoolbook; intended for a small b).
    void mul(const BigInt& b) {
        if (is_zero() || b.is_zero()) { d.clear(); return; }
        std::vector<std::uint32_t> r(d.size() + b.d.size(), 0);
        for (std::size_t i = 0; i < b.d.size(); ++i) {
            std::uint64_t c = 0;
            const std::uint64_t bi = b.d[i];
            for (std::size_t j = 0; j < d.size(); ++j) {
                const std::uint64_t t = static_cast<std::uint64_t>(r[i + j]) + bi * d[j] + c;
                r[i + j] = static_cast<std::uint32_t>(t);
                c = t >> 32;
            }
            std::size_t k = i + d.size();
            while (c) {
                const std::uint64_t t = static_cast<std::uint64_t>(r[k]) + c;
                r[k++] = static_cast<std::uint32_t>(t);
                c = t >> 32;
            }
        }
        d = std::move(r);
        normalize();
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
            // divide by 10^9 per step
            std::uint64_t rem = 0;
            for (std::size_t i = tmp.d.size(); i-- > 0;) {
                const std::uint64_t cur = (rem << 32) | tmp.d[i];
                tmp.d[i] = static_cast<std::uint32_t>(cur / 1000000000ULL);
                rem = cur % 1000000000ULL;
            }
            tmp.normalize();
            for (int k = 0; k < 9; ++k) {
                s.push_back(static_cast<char>('0' + rem % 10));
                rem /= 10;
                if (tmp.is_zero() && rem == 0) break;
            }
        }
        std::reverse(s.begin(), s.end());
        return s;
    }
};

} // namespace base256::detail
