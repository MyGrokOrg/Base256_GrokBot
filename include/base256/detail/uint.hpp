#pragma once

// 256-bit and 512-bit unsigned integers used by the binary256 soft-float core.
// Little-endian limbs: d[0] is the least-significant 64 bits.
// Portable: no dependence on __int128 (GCC/Clang use it as a fast path).
// Without it (MSVC, clang-cl, or -DBASE256_NO_INT128 to test that path with
// GCC/Clang) the same algorithms run on 64-bit words: MSVC's _umul128/__umulh
// and _udiv128 at run time, portable code in constant evaluation and elsewhere.
// Every routine returns the same (unique) integer result either way.

#include <array>
#include <bit>
#include <compare>
#include <cstdint>
#include <type_traits>
#include <utility>

#if defined(__SIZEOF_INT128__) && !defined(_MSC_VER) && !defined(BASE256_NO_INT128)
#define BASE256_HAS_INT128 1
#endif
#if defined(_MSC_VER) && !defined(__clang__) && (defined(_M_X64) || defined(_M_ARM64))
#include <intrin.h>
#define BASE256_MSVC_INTRIN 1
#endif

namespace base256::detail {

struct u256;
struct u512;

#if defined(BASE256_HAS_INT128)
    __extension__ using u128 = unsigned __int128;
#endif

// 64x64 -> 128 multiply, portable.
constexpr void mulu64(std::uint64_t a, std::uint64_t b,
                      std::uint64_t& lo, std::uint64_t& hi) noexcept {
#if defined(BASE256_HAS_INT128)
    const u128 p = static_cast<u128>(a) * b;
    lo = static_cast<std::uint64_t>(p);
    hi = static_cast<std::uint64_t>(p >> 64);
#else
#if defined(BASE256_MSVC_INTRIN)
    if (!std::is_constant_evaluated()) {
#if defined(_M_X64)
        lo = _umul128(a, b, &hi);
#else
        lo = a * b;
        hi = __umulh(a, b);
#endif
        return;
    }
#endif
    const std::uint64_t a0 = a & 0xffffffffULL, a1 = a >> 32;
    const std::uint64_t b0 = b & 0xffffffffULL, b1 = b >> 32;
    const std::uint64_t p0 = a0 * b0;
    const std::uint64_t p1 = a0 * b1;
    const std::uint64_t p2 = a1 * b0;
    const std::uint64_t p3 = a1 * b1;
    const std::uint64_t mid = (p0 >> 32) + (p1 & 0xffffffffULL) + (p2 & 0xffffffffULL);
    lo = (p0 & 0xffffffffULL) | (mid << 32);
    hi = p3 + (p1 >> 32) + (p2 >> 32) + (mid >> 32);
#endif
}

constexpr std::uint64_t addc64(std::uint64_t a, std::uint64_t b,
                               std::uint64_t& carry) noexcept {
#if defined(BASE256_HAS_INT128)
    const u128 s = static_cast<u128>(a) + b + carry;
    carry = static_cast<std::uint64_t>(s >> 64);
    return static_cast<std::uint64_t>(s);
#else
    const std::uint64_t s = a + b;
    const std::uint64_t c1 = s < a;
    const std::uint64_t s2 = s + carry;
    const std::uint64_t c2 = s2 < s;
    carry = c1 | c2;
    return s2;
#endif
}

constexpr std::uint64_t subb64(std::uint64_t a, std::uint64_t b,
                               std::uint64_t& borrow) noexcept {
#if defined(BASE256_HAS_INT128)
    const u128 d = static_cast<u128>(a) - b - borrow;
    borrow = static_cast<std::uint64_t>((d >> 64) & 1);
    return static_cast<std::uint64_t>(d);
#else
    const std::uint64_t t = a - borrow;
    const std::uint64_t b1 = a < borrow;
    const std::uint64_t d = t - b;
    const std::uint64_t b2 = t < b;
    borrow = b1 | b2;
    return d;
#endif
}

// (hi:lo) / d for a normalized d (top bit set) and hi < d: the 64-bit quotient,
// remainder in rem (Hacker's Delight divlu, two 32-bit digit steps).
constexpr std::uint64_t div128_64(std::uint64_t hi, std::uint64_t lo, std::uint64_t d,
                                  std::uint64_t& rem) noexcept {
#if defined(BASE256_MSVC_INTRIN) && defined(_M_X64)
    if (!std::is_constant_evaluated()) return _udiv128(hi, lo, d, &rem);
#endif
    constexpr std::uint64_t b = 1ULL << 32, mask = b - 1;
    const std::uint64_t d1 = d >> 32, d0 = d & mask;
    const std::uint64_t un1 = lo >> 32, un0 = lo & mask;
    std::uint64_t q1 = hi / d1, rhat = hi - q1 * d1;
    while (q1 >= b || q1 * d0 > ((rhat << 32) | un1)) {
        --q1;
        rhat += d1;
        if (rhat >= b) break;
    }
    const std::uint64_t un21 = (hi << 32) + un1 - q1 * d;
    std::uint64_t q0 = un21 / d1;
    rhat = un21 - q0 * d1;
    while (q0 >= b || q0 * d0 > ((rhat << 32) | un0)) {
        --q0;
        rhat += d1;
        if (rhat >= b) break;
    }
    rem = (un21 << 32) + un0 - q0 * d;
    return (q1 << 32) | q0;
}

struct u256 {
    std::uint64_t d[4]{}; // d[0] = LSB

    constexpr u256() noexcept = default;
    constexpr explicit u256(std::uint64_t x) noexcept : d{x, 0, 0, 0} {}

    static constexpr u256 from_words(std::uint64_t w0, std::uint64_t w1,
                                     std::uint64_t w2, std::uint64_t w3) noexcept {
        u256 r;
        r.d[0] = w0;
        r.d[1] = w1;
        r.d[2] = w2;
        r.d[3] = w3;
        return r;
    }

    [[nodiscard]] constexpr bool is_zero() const noexcept {
        return (d[0] | d[1] | d[2] | d[3]) == 0;
    }

    [[nodiscard]] constexpr std::strong_ordering operator<=>(const u256& o) const noexcept {
        for (int i = 3; i >= 0; --i) {
            if (d[i] < o.d[i]) return std::strong_ordering::less;
            if (d[i] > o.d[i]) return std::strong_ordering::greater;
        }
        return std::strong_ordering::equal;
    }
    [[nodiscard]] constexpr bool operator==(const u256& o) const noexcept = default;

    [[nodiscard]] constexpr int clz() const noexcept {
        for (int i = 3; i >= 0; --i) {
            if (d[i] != 0) {
                return (3 - i) * 64 + std::countl_zero(d[i]);
            }
        }
        return 256;
    }

    [[nodiscard]] constexpr int ctz() const noexcept {
        for (int i = 0; i < 4; ++i) {
            if (d[i] != 0) {
                return i * 64 + std::countr_zero(d[i]);
            }
        }
        return 256;
    }

    [[nodiscard]] constexpr int bit_width() const noexcept { return 256 - clz(); }

    [[nodiscard]] constexpr bool bit(int i) const noexcept {
        if (i < 0 || i >= 256) return false;
        return ((d[i >> 6] >> (i & 63)) & 1ULL) != 0;
    }

    constexpr void set_bit(int i) noexcept {
        if (i < 0 || i >= 256) return;
        d[i >> 6] |= 1ULL << (i & 63);
    }

    constexpr void clear_bit(int i) noexcept {
        if (i < 0 || i >= 256) return;
        d[i >> 6] &= ~(1ULL << (i & 63));
    }

    [[nodiscard]] constexpr u256 operator~() const noexcept {
        u256 r;
        r.d[0] = ~d[0];
        r.d[1] = ~d[1];
        r.d[2] = ~d[2];
        r.d[3] = ~d[3];
        return r;
    }
};

struct u512 {
    std::uint64_t d[8]{}; // d[0] = LSB

    constexpr u512() noexcept = default;
    constexpr explicit u512(const u256& lo) noexcept {
        d[0] = lo.d[0];
        d[1] = lo.d[1];
        d[2] = lo.d[2];
        d[3] = lo.d[3];
    }

    [[nodiscard]] constexpr bool is_zero() const noexcept {
        return (d[0] | d[1] | d[2] | d[3] | d[4] | d[5] | d[6] | d[7]) == 0;
    }

    [[nodiscard]] constexpr u256 lo256() const noexcept {
        return u256::from_words(d[0], d[1], d[2], d[3]);
    }
    [[nodiscard]] constexpr u256 hi256() const noexcept {
        return u256::from_words(d[4], d[5], d[6], d[7]);
    }

    [[nodiscard]] constexpr bool bit(int i) const noexcept {
        if (i < 0 || i >= 512) return false;
        return ((d[i >> 6] >> (i & 63)) & 1ULL) != 0;
    }

    constexpr void set_bit(int i) noexcept {
        if (i < 0 || i >= 512) return;
        d[i >> 6] |= 1ULL << (i & 63);
    }

    [[nodiscard]] constexpr int clz() const noexcept {
        for (int i = 7; i >= 0; --i) {
            if (d[i] != 0) {
                return (7 - i) * 64 + std::countl_zero(d[i]);
            }
        }
        return 512;
    }

    [[nodiscard]] constexpr std::strong_ordering operator<=>(const u512& o) const noexcept {
        for (int i = 7; i >= 0; --i) {
            if (d[i] < o.d[i]) return std::strong_ordering::less;
            if (d[i] > o.d[i]) return std::strong_ordering::greater;
        }
        return std::strong_ordering::equal;
    }
    [[nodiscard]] constexpr bool operator==(const u512& o) const noexcept = default;
};

[[nodiscard]] constexpr u256 add_u256(const u256& a, const u256& b,
                                      std::uint64_t& carry) noexcept {
    u256 r;
    carry = 0;
    r.d[0] = addc64(a.d[0], b.d[0], carry);
    r.d[1] = addc64(a.d[1], b.d[1], carry);
    r.d[2] = addc64(a.d[2], b.d[2], carry);
    r.d[3] = addc64(a.d[3], b.d[3], carry);
    return r;
}

[[nodiscard]] constexpr u256 operator+(const u256& a, const u256& b) noexcept {
    std::uint64_t c = 0;
    return add_u256(a, b, c);
}

[[nodiscard]] constexpr u256 sub_u256(const u256& a, const u256& b,
                                      std::uint64_t& borrow) noexcept {
    u256 r;
    borrow = 0;
    r.d[0] = subb64(a.d[0], b.d[0], borrow);
    r.d[1] = subb64(a.d[1], b.d[1], borrow);
    r.d[2] = subb64(a.d[2], b.d[2], borrow);
    r.d[3] = subb64(a.d[3], b.d[3], borrow);
    return r;
}

[[nodiscard]] constexpr u256 operator-(const u256& a, const u256& b) noexcept {
    std::uint64_t br = 0;
    return sub_u256(a, b, br);
}

[[nodiscard]] constexpr u256 operator<<(const u256& a, int n) noexcept {
    if (n <= 0) return a;
    if (n >= 256) return {};
    const int w = n >> 6;
    const int b = n & 63;
    u256 r{};
    if (b == 0) {
        for (int i = 3; i >= w; --i) r.d[i] = a.d[i - w];
    } else {
        for (int i = 3; i >= w; --i) {
            r.d[i] = a.d[i - w] << b;
            if (i - w - 1 >= 0) r.d[i] |= a.d[i - w - 1] >> (64 - b);
        }
    }
    return r;
}

[[nodiscard]] constexpr u256 operator>>(const u256& a, int n) noexcept {
    if (n <= 0) return a;
    if (n >= 256) return {};
    const int w = n >> 6;
    const int b = n & 63;
    u256 r{};
    if (b == 0) {
        for (int i = 0; i <= 3 - w; ++i) r.d[i] = a.d[i + w];
    } else {
        for (int i = 0; i <= 3 - w; ++i) {
            r.d[i] = a.d[i + w] >> b;
            if (i + w + 1 <= 3) r.d[i] |= a.d[i + w + 1] << (64 - b);
        }
    }
    return r;
}

// Shift right, OR 1 into the result if any 1-bits fall off (sticky).
[[nodiscard]] constexpr u256 shr_sticky(const u256& a, int n, bool& sticky) noexcept {
    if (n <= 0) return a;
    if (n >= 256) {
        sticky = sticky || !a.is_zero();
        return {};
    }
    // Bits that would be lost: the low n bits of a.
    if (n >= 1) {
        // Check low n bits.
        const int w = n >> 6;
        const int b = n & 63;
        for (int i = 0; i < w; ++i) {
            if (a.d[i] != 0) sticky = true;
        }
        if (b && ((a.d[w] & ((1ULL << b) - 1ULL)) != 0)) sticky = true;
    }
    return a >> n;
}

[[nodiscard]] constexpr u512 operator<<(const u512& a, int n) noexcept {
    if (n <= 0) return a;
    if (n >= 512) return {};
    const int w = n >> 6;
    const int b = n & 63;
    u512 r{};
    if (b == 0) {
        for (int i = 7; i >= w; --i) r.d[i] = a.d[i - w];
    } else {
        for (int i = 7; i >= w; --i) {
            r.d[i] = a.d[i - w] << b;
            if (i - w - 1 >= 0) r.d[i] |= a.d[i - w - 1] >> (64 - b);
        }
    }
    return r;
}

[[nodiscard]] constexpr u512 operator>>(const u512& a, int n) noexcept {
    if (n <= 0) return a;
    if (n >= 512) return {};
    const int w = n >> 6;
    const int b = n & 63;
    u512 r{};
    if (b == 0) {
        for (int i = 0; i <= 7 - w; ++i) r.d[i] = a.d[i + w];
    } else {
        for (int i = 0; i <= 7 - w; ++i) {
            r.d[i] = a.d[i + w] >> b;
            if (i + w + 1 <= 7) r.d[i] |= a.d[i + w + 1] << (64 - b);
        }
    }
    return r;
}

[[nodiscard]] constexpr u512 add_u512(const u512& a, const u512& b) noexcept {
    u512 r;
    std::uint64_t c = 0;
    for (int i = 0; i < 8; ++i) r.d[i] = addc64(a.d[i], b.d[i], c);
    return r;
}

[[nodiscard]] constexpr u512 sub_u512(const u512& a, const u512& b) noexcept {
    u512 r;
    std::uint64_t br = 0;
    for (int i = 0; i < 8; ++i) r.d[i] = subb64(a.d[i], b.d[i], br);
    return r;
}

[[nodiscard]] constexpr u512 mul_wide(const u256& a, const u256& b) noexcept {
    u512 r{};
#if defined(BASE256_HAS_INT128)
    u128 acc[8]{};
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            const u128 p = static_cast<u128>(a.d[i]) * b.d[j];
            acc[i + j] += static_cast<std::uint64_t>(p);
            acc[i + j + 1] += static_cast<std::uint64_t>(p >> 64);
        }
    }
    u128 c = 0;
    for (int i = 0; i < 8; ++i) {
        c += acc[i];
        r.d[i] = static_cast<std::uint64_t>(c);
        c >>= 64;
    }
#else
    for (int i = 0; i < 4; ++i) {
        std::uint64_t carry = 0;
        for (int j = 0; j < 4; ++j) {
            std::uint64_t lo, hi;
            mulu64(a.d[i], b.d[j], lo, hi);
            // r[i+j] + lo + carry, with the two carry-outs kept apart (the carry-in of the
            // second add must be 0, not the first carry-out): a*b + r + carry < 2^128, so
            // hi + c1 + c2 cannot overflow.
            std::uint64_t c1 = 0, c2 = 0;
            const std::uint64_t s0 = addc64(r.d[i + j], lo, c1);
            r.d[i + j] = addc64(s0, carry, c2);
            carry = hi + c1 + c2;
        }
        std::uint64_t c = 0;
        r.d[i + 4] = addc64(r.d[i + 4], carry, c);
    }
#endif
    return r;
}

// Restoring division of a 512-bit numerator by a 256-bit denominator.
// Quotient may be up to 257 bits; remainder < d.
struct Div512 {
    u512 quot;
    u256 rem;
};

[[nodiscard]] constexpr Div512 divmod_512_256_bitwise(u512 n, u256 d) noexcept {
    Div512 r{};
    if (d.is_zero()) return r; // caller must not divide by zero
    u512 rem{};
    for (int i = 511; i >= 0; --i) {
        // rem <<= 1
        rem = rem << 1;
        if (n.bit(i)) rem.d[0] |= 1ULL;
        // if rem >= d then rem -= d, quot bit = 1
        // rem is at most 257 bits here; compare low 256 + bit 256
        bool ge = false;
        if (rem.d[4] != 0 || rem.d[5] || rem.d[6] || rem.d[7]) {
            ge = true;
        } else {
            const u256 rlo = rem.lo256();
            ge = rlo >= d;
        }
        if (ge) {
            // rem -= d (d is 256-bit)
            std::uint64_t br = 0;
            rem.d[0] = subb64(rem.d[0], d.d[0], br);
            rem.d[1] = subb64(rem.d[1], d.d[1], br);
            rem.d[2] = subb64(rem.d[2], d.d[2], br);
            rem.d[3] = subb64(rem.d[3], d.d[3], br);
            rem.d[4] = subb64(rem.d[4], 0, br);
            r.quot.set_bit(i);
        }
    }
    r.rem = rem.lo256();
    return r;
}

// Knuth algorithm D on 64-bit limbs (TAOCP 4.3.1), with a 128-bit integer
// type or with 64-bit words (div128_64, mulu64).
[[nodiscard]] constexpr Div512 divmod_512_256(u512 n, u256 d) noexcept {
#if defined(BASE256_HAS_INT128)
    Div512 r{};
    int nd = 4;
    while (nd > 0 && d.d[nd - 1] == 0) --nd;
    if (nd == 0) return r;  // caller must not divide by zero
    const int s = std::countl_zero(d.d[nd - 1]);
    std::uint64_t v[4]{}, u[9]{};
    for (int i = nd - 1; i >= 0; --i)
        v[i] = (d.d[i] << s) | (s && i ? d.d[i - 1] >> (64 - s) : 0);
    u[8] = s ? n.d[7] >> (64 - s) : 0;
    for (int i = 7; i >= 0; --i) u[i] = (n.d[i] << s) | (s && i ? n.d[i - 1] >> (64 - s) : 0);
    const std::uint64_t vt = v[nd - 1], vs = nd > 1 ? v[nd - 2] : 0;
    for (int j = 8 - nd; j >= 0; --j) {
        const u128 num = (static_cast<u128>(u[j + nd]) << 64) | u[j + nd - 1];
        u128 qhat = num / vt, rhat = num % vt;
        while (qhat >> 64 || (nd > 1 && qhat * vs > ((rhat << 64) | u[j + nd - 2]))) {
            --qhat;
            rhat += vt;
            if (rhat >> 64) break;
        }
        // u[j .. j+nd] -= qhat * v
        std::uint64_t borrow = 0, carry = 0;
        for (int i = 0; i < nd; ++i) {
            const u128 p = qhat * v[i] + carry;
            carry = static_cast<std::uint64_t>(p >> 64);
            u[i + j] = subb64(u[i + j], static_cast<std::uint64_t>(p), borrow);
        }
        u[j + nd] = subb64(u[j + nd], carry, borrow);
        if (borrow) {  // qhat was one too large: add v back
            --qhat;
            std::uint64_t c = 0;
            for (int i = 0; i < nd; ++i) u[i + j] = addc64(u[i + j], v[i], c);
            u[j + nd] += c;
        }
        r.quot.d[j] = static_cast<std::uint64_t>(qhat);
    }
    for (int i = 0; i < nd; ++i) r.rem.d[i] = (u[i] >> s) | (s && i + 1 < nd ? u[i + 1] << (64 - s) : 0);
    return r;
#else
    Div512 r{};
    int nd = 4;
    while (nd > 0 && d.d[nd - 1] == 0) --nd;
    if (nd == 0) return r;  // caller must not divide by zero
    const int s = std::countl_zero(d.d[nd - 1]);
    std::uint64_t v[4]{}, u[9]{};
    for (int i = nd - 1; i >= 0; --i)
        v[i] = (d.d[i] << s) | (s && i ? d.d[i - 1] >> (64 - s) : 0);
    u[8] = s ? n.d[7] >> (64 - s) : 0;
    for (int i = 7; i >= 0; --i) u[i] = (n.d[i] << s) | (s && i ? n.d[i - 1] >> (64 - s) : 0);
    const std::uint64_t vt = v[nd - 1], vs = nd > 1 ? v[nd - 2] : 0;
    for (int j = 8 - nd; j >= 0; --j) {
        // qhat = min(floor((u[j+nd]:u[j+nd-1]) / vt), 2^64 - 1), rhat the matching remainder
        std::uint64_t qhat, rhat;
        bool rhat_big;  // rhat >= 2^64: the correction test below cannot succeed
        if (u[j + nd] >= vt) {  // (u[j+nd] == vt): the quotient digit is at most 2^64 - 1
            qhat = ~0ULL;
            rhat = u[j + nd - 1] + vt;
            rhat_big = rhat < vt;
        } else {
            qhat = div128_64(u[j + nd], u[j + nd - 1], vt, rhat);
            rhat_big = false;
        }
        while (!rhat_big && nd > 1) {  // qhat * vs > (rhat : u[j+nd-2]) ?
            std::uint64_t plo, phi;
            mulu64(qhat, vs, plo, phi);
            if (phi < rhat || (phi == rhat && plo <= u[j + nd - 2])) break;
            --qhat;
            rhat += vt;
            rhat_big = rhat < vt;
        }
        // u[j .. j+nd] -= qhat * v
        std::uint64_t borrow = 0, carry = 0;
        for (int i = 0; i < nd; ++i) {
            std::uint64_t plo, phi;
            mulu64(qhat, v[i], plo, phi);
            plo += carry;
            phi += plo < carry;
            carry = phi;
            u[i + j] = subb64(u[i + j], plo, borrow);
        }
        u[j + nd] = subb64(u[j + nd], carry, borrow);
        if (borrow) {  // qhat was one too large: add v back
            --qhat;
            std::uint64_t c = 0;
            for (int i = 0; i < nd; ++i) u[i + j] = addc64(u[i + j], v[i], c);
            u[j + nd] += c;
        }
        r.quot.d[j] = qhat;
    }
    for (int i = 0; i < nd; ++i) r.rem.d[i] = (u[i] >> s) | (s && i + 1 < nd ? u[i + 1] << (64 - s) : 0);
    return r;
#endif
}

[[nodiscard]] constexpr u512 operator|(const u512& a, const u512& b) noexcept {
    u512 r;
    for (int i = 0; i < 8; ++i) r.d[i] = a.d[i] | b.d[i];
    return r;
}

// Integer square root of a 512-bit value: floor(sqrt(n)) as a 256-bit integer.
// Hacker's Delight restoring algorithm over 256 even bit-places.
[[nodiscard]] constexpr u256 isqrt_512_bitwise(u512 x) noexcept {
    u512 m{};
    m.set_bit(510); // highest even bit of a 512-bit word
    u512 y{};
    while (!m.is_zero()) {
        const u512 b = y | m;
        y = y >> 1;
        if (x >= b) {
            x = sub_u512(x, b);
            y = y | m;
        }
        m = m >> 2;
    }
    return y.lo256();
}

// floor(sqrt(n)): a 64-bit seed from the top 128 bits, then Newton steps
// x <- (x + n/x)/2 from above (monotone decreasing to the floor). Used when
// the top limb is non-zero (always so for float256::sqrt); else bitwise.
[[nodiscard]] constexpr u256 isqrt_512(u512 n) noexcept {
#if defined(BASE256_HAS_INT128)
    if (n.d[7] < (1ULL << 62)) return isqrt_512_bitwise(n);
    const u128 t = (static_cast<u128>(n.d[7]) << 64) | n.d[6];
    // 64-bit bitwise isqrt of t
    u128 rem = t, root = 0, bit = static_cast<u128>(1) << 126;
    while (bit) {
        if (rem >= root + bit) {
            rem -= root + bit;
            root = (root >> 1) + bit;
        } else {
            root >>= 1;
        }
        bit >>= 2;
    }
    // sqrt(n) < (root + 1) * 2^192 <= 2^256
    u512 x{};
    const u128 r1 = root + 1;
    if (r1 >> 64) {  // root = 2^64 - 1: start just below 2^256
        for (int i = 0; i < 4; ++i) x.d[i] = ~0ULL;
    } else {
        x.d[3] = static_cast<std::uint64_t>(r1);
    }
    for (;;) {
        const Div512 q = divmod_512_256(n, x.lo256());
        const u512 y = add_u512(x, q.quot) >> 1;
        if (!(y < x)) break;
        x = y;
    }
    return x.lo256();
#else
    if (n.d[7] < (1ULL << 62)) return isqrt_512_bitwise(n);
    // 32-bit seed s = floor(sqrt(n.d[7])) >= 2^31: sqrt(n) < (s + 1) * 2^224.
    std::uint64_t rem = n.d[7], root = 0, bit = 1ULL << 62;
    while (bit) {
        if (rem >= root + bit) {
            rem -= root + bit;
            root = (root >> 1) + bit;
        } else {
            root >>= 1;
        }
        bit >>= 2;
    }
    u512 x{};
    if (root + 1 == (1ULL << 32)) {  // start just below 2^256
        for (int i = 0; i < 4; ++i) x.d[i] = ~0ULL;
    } else {
        x.d[3] = (root + 1) << 32;
    }
    for (;;) {
        const Div512 q = divmod_512_256(n, x.lo256());
        const u512 y = add_u512(x, q.quot) >> 1;
        if (!(y < x)) break;
        x = y;
    }
    return x.lo256();
#endif
}

} // namespace base256::detail
