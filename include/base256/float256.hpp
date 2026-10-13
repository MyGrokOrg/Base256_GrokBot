#pragma once

/// IEEE 754 binary256 (octuple-precision) floating-point type.
///
/// Layout (MSB first): 1 sign bit, 19 exponent bits (bias 262143),
/// 236 trailing significand bits. Precision p = 237 bits including the
/// implicit leading 1 of normals. Default rounding: round-to-nearest, ties
/// to even. Header-only, C++23.

#include "detail/uint.hpp"
#include "detail/bigint.hpp"

#include <array>
#include <bit>
#include <cctype>
#include <climits>
#include <cmath>
#include <compare>
#include <cstdint>
#include <cstring>
#include <limits>
#include <optional>
#include <ostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>

// glibc/libstdc++ define isnan/isinf/issubnormal as function-like macros.
#ifdef isnan
#undef isnan
#endif
#ifdef isinf
#undef isinf
#endif
#ifdef isfinite
#undef isfinite
#endif
#ifdef isnormal
#undef isnormal
#endif
#ifdef issubnormal
#undef issubnormal
#endif
#ifdef signbit
#undef signbit
#endif
#ifdef fpclassify
#undef fpclassify
#endif

namespace base256 {

class float256;

namespace detail {

inline constexpr int kExpBits      = 19;
inline constexpr int kFracBits     = 236;
inline constexpr int kPrecision    = 237;
inline constexpr int kBias         = 262143;
inline constexpr int kEmin         = -262142;
inline constexpr int kEmax         =  262143;
inline constexpr int kExpAllOnes   = 0x7FFFF;
inline constexpr int kLsbPos       = 19;            // 255 - 236
inline constexpr std::uint64_t kFracHiMask = 0xFFFFFFFFFFFULL; // 44 bits

enum class Kind : unsigned char { zero, subnormal, normal, inf, nan };

struct Unpacked {
    bool sign = false;
    Kind kind = Kind::zero;
    int  wexp = 0;          // value = ± (sig / 2^255) * 2^wexp
    u256 sig{};             // leading 1 at bit 255 for normals
    bool sticky = false;
    bool signaling = false;
};

} // namespace detail

class float256 {
public:
    using storage_type = std::array<std::uint64_t, 4>; // [0] = LSB word

    static constexpr int digits         = detail::kPrecision;
    static constexpr int digits10       = 71;
    static constexpr int max_digits10   = 73;
    static constexpr int max_exponent   = 262144;
    static constexpr int min_exponent   = -262141;
    static constexpr int max_exponent10 = 78912;
    static constexpr int min_exponent10 = -78913;
    static constexpr int radix          = 2;
    static constexpr int exponent_bias  = detail::kBias;

private:
    storage_type bits_{};

    friend struct std::numeric_limits<float256>;
    friend constexpr float256 copysign(float256, float256) noexcept;
    friend constexpr float256 nextafter(float256, float256) noexcept;
    friend constexpr float256 ldexp(float256, int) noexcept;
    friend constexpr float256 frexp(float256, int*) noexcept;
    friend constexpr float256 trunc(float256) noexcept;
    friend constexpr float256 floor(float256) noexcept;
    friend constexpr float256 ceil(float256) noexcept;
    friend constexpr float256 round(float256) noexcept;
    friend constexpr float256 rint(float256) noexcept;
    friend constexpr int ilogb(float256) noexcept;
    friend constexpr float256 logb(float256) noexcept;

    [[nodiscard]] constexpr bool sign_bit() const noexcept {
        return (bits_[3] >> 63) != 0;
    }
    [[nodiscard]] constexpr int biased_exp() const noexcept {
        return static_cast<int>((bits_[3] >> 44) & detail::kExpAllOnes);
    }
    [[nodiscard]] constexpr detail::u256 frac_bits() const noexcept {
        return detail::u256::from_words(
            bits_[0], bits_[1], bits_[2], bits_[3] & detail::kFracHiMask);
    }
    constexpr void set_sign_bit(bool s) noexcept {
        if (s) bits_[3] |= (1ULL << 63);
        else   bits_[3] &= ~(1ULL << 63);
    }

    static constexpr float256 from_encoding(std::uint64_t w0, std::uint64_t w1,
                                            std::uint64_t w2, std::uint64_t w3) noexcept {
        float256 r;
        r.bits_[0] = w0;
        r.bits_[1] = w1;
        r.bits_[2] = w2;
        r.bits_[3] = w3;
        return r;
    }

    [[nodiscard]] constexpr detail::Unpacked unpack() const noexcept {
        using namespace detail;
        Unpacked u;
        u.sign = sign_bit();
        const int be = biased_exp();
        const u256 frac = frac_bits();
        if (be == kExpAllOnes) {
            if (frac.is_zero()) {
                u.kind = Kind::inf;
            } else {
                u.kind = Kind::nan;
                u.sig = frac;
                u.signaling = (bits_[3] & (1ULL << 43)) == 0;
            }
            return u;
        }
        if (be == 0) {
            if (frac.is_zero()) {
                u.kind = Kind::zero;
                return u;
            }
            u.kind = Kind::subnormal;
            u.sig = frac;
            const int lz = u.sig.clz();
            u.sig = u.sig << lz;
            u.wexp = kEmin + kLsbPos - lz;
            return u;
        }
        u.kind = Kind::normal;
        u.sig = frac;
        u.sig.set_bit(kFracBits);
        u.sig = u.sig << kLsbPos;
        u.wexp = be - kBias;
        return u;
    }

    static constexpr bool round_ne(detail::u256& sig, int lsb_pos, bool sticky) noexcept {
        if (lsb_pos <= 0) return false;
        const bool g = sig.bit(lsb_pos - 1);
        bool rs = sticky;
        if (!rs && lsb_pos > 1) {
            if (lsb_pos - 1 <= 64) {
                rs = (sig.d[0] & ((1ULL << (lsb_pos - 1)) - 1ULL)) != 0;
            } else {
                for (int i = 0; i <= lsb_pos - 2; ++i) {
                    if (sig.bit(i)) { rs = true; break; }
                }
            }
        }
        const bool lsb = sig.bit(lsb_pos);
        const bool inc = g && (rs || lsb);
        if (lsb_pos >= 256) {
            sig = {};
            return inc;
        }
        sig = (sig >> lsb_pos) << lsb_pos;
        if (!inc) return false;
        detail::u256 one{};
        one.set_bit(lsb_pos);
        std::uint64_t carry = 0;
        sig = detail::add_u256(sig, one, carry);
        return carry != 0;
    }

    static constexpr float256 pack(detail::Unpacked u) noexcept {
        using namespace detail;
        if (u.kind == Kind::nan) {
            std::uint64_t w3 = (1ULL << 63) * (u.sign ? 1ULL : 0ULL);
            w3 |= static_cast<std::uint64_t>(kExpAllOnes) << 44;
            w3 |= (1ULL << 43);
            w3 |= u.sig.d[3] & kFracHiMask;
            return from_encoding(u.sig.d[0], u.sig.d[1], u.sig.d[2], w3);
        }
        if (u.kind == Kind::inf) return inf(u.sign);
        if (u.kind == Kind::zero || (u.sig.is_zero() && !u.sticky)) return zero(u.sign);

        const std::int64_t biased = static_cast<std::int64_t>(u.wexp) + kBias;
        if (biased >= kExpAllOnes) return inf(u.sign);

        u256 sig = u.sig;
        bool sticky = u.sticky;
        int wexp = u.wexp;

        if (biased <= 0) {
            const int rshift = kLsbPos - wexp + kEmin;
            if (rshift >= 256 + 2) return zero(u.sign);
            if (rshift > 0) {
                const bool ovf = round_ne(sig, rshift, sticky);
                if (ovf) {
                    return pack_normal(u.sign, kEmin, u256::from_words(0, 0, 0, 1ULL << 63));
                }
                sig = sig >> rshift;
            }
            if (sig.bit(kFracBits)) {
                return pack_normal(u.sign, kEmin, sig << kLsbPos);
            }
            return from_encoding(sig.d[0], sig.d[1], sig.d[2],
                                 (sig.d[3] & kFracHiMask) | (u.sign ? (1ULL << 63) : 0ULL));
        }

        const bool ovf = round_ne(sig, kLsbPos, sticky);
        if (ovf) {
            ++wexp;
            sig = {};
            sig.set_bit(255);
        }
        if (wexp > kEmax) return inf(u.sign);
        if (wexp < kEmin) {
            Unpacked v;
            v.sign = u.sign;
            v.kind = Kind::subnormal;
            v.wexp = wexp;
            v.sig = sig;
            return pack(v);
        }
        return pack_normal(u.sign, wexp, sig);
    }

    static constexpr float256 pack_normal(bool sign, int wexp, detail::u256 sig_left) noexcept {
        using namespace detail;
        const u256 f = sig_left >> kLsbPos;
        const int be = wexp + kBias;
        std::uint64_t w3 = (sign ? (1ULL << 63) : 0ULL)
                         | (static_cast<std::uint64_t>(be) << 44)
                         | (f.d[3] & kFracHiMask);
        return from_encoding(f.d[0], f.d[1], f.d[2], w3);
    }

    static constexpr float256 from_unpacked(detail::Unpacked u) noexcept { return pack(u); }

    static constexpr detail::Unpacked make_nan(bool signaling = false,
                                               bool sign = false) noexcept {
        detail::Unpacked u;
        u.sign = sign;
        u.kind = detail::Kind::nan;
        u.signaling = signaling;
        u.sig.set_bit(235);
        if (signaling) u.sig.clear_bit(235);
        u.sig.set_bit(0);
        return u;
    }

public:
    constexpr float256() noexcept = default;
    constexpr float256(const float256&) noexcept = default;
    constexpr float256(float256&&) noexcept = default;
    constexpr float256& operator=(const float256&) noexcept = default;
    constexpr float256& operator=(float256&&) noexcept = default;

    constexpr explicit float256(double v) noexcept { *this = from_double(v); }
    constexpr explicit float256(float v) noexcept { *this = from_double(static_cast<double>(v)); }
    explicit float256(long double v) noexcept { *this = from_double(static_cast<double>(v)); }

    constexpr explicit float256(std::int64_t v) noexcept { *this = from_int(v); }
    constexpr explicit float256(std::int32_t v) noexcept { *this = from_int(static_cast<std::int64_t>(v)); }
    constexpr explicit float256(std::uint64_t v) noexcept { *this = from_uint(v); }
    constexpr explicit float256(std::uint32_t v) noexcept { *this = from_uint(static_cast<std::uint64_t>(v)); }

    explicit float256(std::string_view s) {
        auto r = parse(s);
        if (!r) throw std::invalid_argument("float256: invalid numeric string");
        *this = *r;
    }

    [[nodiscard]] static constexpr float256 zero(bool negative = false) noexcept {
        return from_encoding(0, 0, 0, negative ? (1ULL << 63) : 0ULL);
    }
    [[nodiscard]] static constexpr float256 inf(bool negative = false) noexcept {
        const std::uint64_t w3 = (negative ? (1ULL << 63) : 0ULL)
                               | (static_cast<std::uint64_t>(detail::kExpAllOnes) << 44);
        return from_encoding(0, 0, 0, w3);
    }
    [[nodiscard]] static constexpr float256 nan(bool signaling = false) noexcept {
        const std::uint64_t quiet = signaling ? 0ULL : (1ULL << 43);
        const std::uint64_t w3 = (static_cast<std::uint64_t>(detail::kExpAllOnes) << 44)
                               | quiet | 1ULL;
        return from_encoding(1, 0, 0, w3);
    }
    [[nodiscard]] static constexpr float256 qnan() noexcept { return nan(false); }
    [[nodiscard]] static constexpr float256 snan() noexcept { return nan(true); }

    [[nodiscard]] static constexpr float256 from_bits(storage_type b) noexcept {
        float256 r;
        r.bits_ = b;
        return r;
    }
    [[nodiscard]] static constexpr float256 from_bits(std::uint64_t w0, std::uint64_t w1,
                                                      std::uint64_t w2, std::uint64_t w3) noexcept {
        return from_encoding(w0, w1, w2, w3);
    }

    [[nodiscard]] constexpr bool signbit()     const noexcept { return sign_bit(); }
    [[nodiscard]] constexpr bool isnan()       const noexcept { return biased_exp() == detail::kExpAllOnes && !frac_bits().is_zero(); }
    [[nodiscard]] constexpr bool isinf()       const noexcept { return biased_exp() == detail::kExpAllOnes &&  frac_bits().is_zero(); }
    [[nodiscard]] constexpr bool isfinite()    const noexcept { return biased_exp() != detail::kExpAllOnes; }
    [[nodiscard]] constexpr bool isnormal()    const noexcept {
        const int be = biased_exp();
        return be != 0 && be != detail::kExpAllOnes;
    }
    [[nodiscard]] constexpr bool is_subnormal() const noexcept {
        return biased_exp() == 0 && !frac_bits().is_zero();
    }
    [[nodiscard]] constexpr bool iszero() const noexcept {
        return biased_exp() == 0 && frac_bits().is_zero();
    }
    [[nodiscard]] constexpr int fpclassify() const noexcept {
        const int be = biased_exp();
        if (be == detail::kExpAllOnes) return frac_bits().is_zero() ? FP_INFINITE : FP_NAN;
        if (be == 0) return frac_bits().is_zero() ? FP_ZERO : FP_SUBNORMAL;
        return FP_NORMAL;
    }

    [[nodiscard]] constexpr storage_type raw_bits() const noexcept { return bits_; }

    [[nodiscard]] static constexpr float256 from_double(double v) noexcept {
        using namespace detail;
        const std::uint64_t u = std::bit_cast<std::uint64_t>(v);
        const bool sign = (u >> 63) != 0;
        const int de = static_cast<int>((u >> 52) & 0x7FF);
        const std::uint64_t dm = u & 0x000FFFFFFFFFFFFFULL;
        if (de == 0x7FF) {
            if (dm == 0) return inf(sign);
            Unpacked n = make_nan(false, sign);
            n.sig.d[0] = dm;
            return pack(n);
        }
        if (de == 0) {
            if (dm == 0) return zero(sign);
            Unpacked uo;
            uo.sign = sign;
            uo.kind = Kind::normal;
            uo.sig = u256{dm};
            const int lz = uo.sig.clz();
            uo.sig = uo.sig << lz;
            uo.wexp = 255 - lz - 1074;
            return pack(uo);
        }
        Unpacked uo;
        uo.sign = sign;
        uo.kind = Kind::normal;
        uo.wexp = de - 1023;
        const std::uint64_t s53 = (1ULL << 52) | dm;
        uo.sig = u256{s53} << (255 - 52);
        return pack(uo);
    }

    [[nodiscard]] static constexpr float256 from_uint(std::uint64_t v) noexcept {
        using namespace detail;
        if (v == 0) return zero(false);
        Unpacked u;
        u.kind = Kind::normal;
        u.sig = u256{v};
        const int lz = u.sig.clz();
        u.sig = u.sig << lz;
        u.wexp = 255 - lz;
        return pack(u);
    }

    [[nodiscard]] static constexpr float256 from_int(std::int64_t v) noexcept {
        if (v < 0) {
            const auto mag = static_cast<std::uint64_t>(-(v + 1)) + 1ULL;
            float256 r = from_uint(mag);
            r.set_sign_bit(true);
            return r;
        }
        return from_uint(static_cast<std::uint64_t>(v));
    }

    [[nodiscard]] constexpr double to_double() const noexcept {
        using namespace detail;
        const Unpacked u = unpack();
        if (u.kind == Kind::nan) return std::numeric_limits<double>::quiet_NaN();
        if (u.kind == Kind::inf) {
            return u.sign ? -std::numeric_limits<double>::infinity()
                          :  std::numeric_limits<double>::infinity();
        }
        if (u.kind == Kind::zero) return u.sign ? -0.0 : 0.0;

        u256 sig = u.sig;
        int wexp = u.wexp;
        const int dbl_lsb = 255 - 52;
        const bool ovf = round_ne(sig, dbl_lsb, u.sticky);
        if (ovf) {
            ++wexp;
            sig = {};
            sig.set_bit(255);
        }
        if (wexp > 1023) {
            return u.sign ? -std::numeric_limits<double>::infinity()
                          :  std::numeric_limits<double>::infinity();
        }
        if (wexp < -1022) {
            const int rshift = (-1022 - wexp) + dbl_lsb;
            bool st = false;
            u256 s = sig;
            const bool ov = round_ne(s, rshift, st);
            s = s >> rshift;
            if (ov || s.bit(52)) {
                const std::uint64_t bits = (u.sign ? (1ULL << 63) : 0ULL) | (1ULL << 52);
                return std::bit_cast<double>(bits);
            }
            const std::uint64_t mant = s.d[0] & 0x000FFFFFFFFFFFFFULL;
            if (mant == 0) return u.sign ? -0.0 : 0.0;
            const std::uint64_t bits = (u.sign ? (1ULL << 63) : 0ULL) | mant;
            return std::bit_cast<double>(bits);
        }
        const std::uint64_t mant = (sig.d[3] >> 11) & 0x000FFFFFFFFFFFFFULL;
        const std::uint64_t be = static_cast<std::uint64_t>(wexp + 1023);
        const std::uint64_t bits = (u.sign ? (1ULL << 63) : 0ULL) | (be << 52) | mant;
        return std::bit_cast<double>(bits);
    }

    [[nodiscard]] constexpr float to_float() const noexcept {
        return static_cast<float>(to_double());
    }

    [[nodiscard]] constexpr std::partial_ordering operator<=>(const float256& o) const noexcept {
        if (isnan() || o.isnan()) return std::partial_ordering::unordered;
        if (iszero() && o.iszero()) return std::partial_ordering::equivalent;
        const bool sa = sign_bit();
        const bool sb = o.sign_bit();
        if (sa != sb) return sa ? std::partial_ordering::less : std::partial_ordering::greater;
        for (int i = 3; i >= 0; --i) {
            if (bits_[i] < o.bits_[i]) return sa ? std::partial_ordering::greater : std::partial_ordering::less;
            if (bits_[i] > o.bits_[i]) return sa ? std::partial_ordering::less    : std::partial_ordering::greater;
        }
        return std::partial_ordering::equivalent;
    }

    [[nodiscard]] constexpr bool operator==(const float256& o) const noexcept {
        if (isnan() || o.isnan()) return false;
        if (iszero() && o.iszero()) return true;
        return bits_ == o.bits_;
    }
    [[nodiscard]] constexpr bool operator!=(const float256& o) const noexcept { return !(*this == o); }
    [[nodiscard]] constexpr bool operator<(const float256& o)  const noexcept { return (*this <=> o) == std::partial_ordering::less; }
    [[nodiscard]] constexpr bool operator>(const float256& o)  const noexcept { return (*this <=> o) == std::partial_ordering::greater; }
    [[nodiscard]] constexpr bool operator<=(const float256& o) const noexcept {
        auto c = *this <=> o;
        return c == std::partial_ordering::less || c == std::partial_ordering::equivalent;
    }
    [[nodiscard]] constexpr bool operator>=(const float256& o) const noexcept {
        auto c = *this <=> o;
        return c == std::partial_ordering::greater || c == std::partial_ordering::equivalent;
    }

    [[nodiscard]] constexpr float256 operator+() const noexcept { return *this; }
    [[nodiscard]] constexpr float256 operator-() const noexcept {
        float256 r = *this;
        r.set_sign_bit(!r.sign_bit());
        return r;
    }

    [[nodiscard]] friend constexpr float256 operator+(float256 a, float256 b) noexcept {
        return add_sub(a, b, false);
    }
    [[nodiscard]] friend constexpr float256 operator-(float256 a, float256 b) noexcept {
        return add_sub(a, b, true);
    }
    [[nodiscard]] friend constexpr float256 operator*(float256 a, float256 b) noexcept {
        return mul(a, b);
    }
    [[nodiscard]] friend constexpr float256 operator/(float256 a, float256 b) noexcept {
        return div(a, b);
    }

    constexpr float256& operator+=(float256 o) noexcept { return *this = *this + o; }
    constexpr float256& operator-=(float256 o) noexcept { return *this = *this - o; }
    constexpr float256& operator*=(float256 o) noexcept { return *this = *this * o; }
    constexpr float256& operator/=(float256 o) noexcept { return *this = *this / o; }

    [[nodiscard]] static constexpr float256 fma(float256 a, float256 b, float256 c) noexcept;
    [[nodiscard]] static constexpr float256 sqrt(float256 x) noexcept;

    [[nodiscard]] std::string to_string(int precision = 0) const;
    [[nodiscard]] std::string to_hex_string() const;
    [[nodiscard]] std::string to_bit_string() const;
    [[nodiscard]] static std::optional<float256> parse(std::string_view s);

    friend std::ostream& operator<<(std::ostream& os, const float256& v) {
        return os << v.to_string();
    }

    [[nodiscard]] static constexpr float256 pi() noexcept {
        return from_bits(0xf98e804177d4c762ULL, 0x839a252049c1114cULL,
                         0x18469898cc51701bULL, 0x40000921fb54442dULL);
    }
    [[nodiscard]] static constexpr float256 e() noexcept {
        return from_bits(0xa6d2b53c26c8228dULL, 0xa79e3b1738b079c5ULL,
                         0x695355fb8ac404e7ULL, 0x400005bf0a8b1457ULL);
    }
    [[nodiscard]] static constexpr float256 ln2() noexcept {
        return from_bits(0x16c5b141a2eb7175ULL, 0x5ed5e81e6864ce53ULL,
                         0xef35793c7673007eULL, 0x3fffe62e42fefa39ULL);
    }
    [[nodiscard]] static constexpr float256 ln10() noexcept {
        return from_bits(0xbf21d078c3d0403eULL, 0x61451c51fd9f3b4bULL,
                         0x1582dd4adac5705aULL, 0x4000026bb1bbb555ULL);
    }
    [[nodiscard]] static constexpr float256 log2e() noexcept {
        return from_bits(0x2b4b1164a2cd9a34ULL, 0xa7d11d6aef551badULL,
                         0xfe1777d0ffda0d23ULL, 0x3ffff71547652b82ULL);
    }
    [[nodiscard]] static constexpr float256 log10e() noexcept {
        return from_bits(0x34404747e5a89ef2ULL, 0x7b8647dc68c048b9ULL,
                         0x0e32a6ab7555f5a6ULL, 0x3fffdbcb7b1526e5ULL);
    }
    [[nodiscard]] static constexpr float256 sqrt2() noexcept {
        return from_bits(0x75099da2f590b066ULL, 0x57d3e3adec175127ULL,
                         0xcc908b2fb1366ea9ULL, 0x3ffff6a09e667f3bULL);
    }
    [[nodiscard]] static constexpr float256 one()  noexcept { return from_bits(0, 0, 0, 0x3ffff00000000000ULL); }
    [[nodiscard]] static constexpr float256 two()  noexcept { return from_bits(0, 0, 0, 0x4000000000000000ULL); }
    [[nodiscard]] static constexpr float256 half() noexcept { return from_bits(0, 0, 0, 0x3fffe00000000000ULL); }

    [[nodiscard]] static constexpr float256 min() noexcept {
        return from_bits(0, 0, 0, 0x0000100000000000ULL);
    }
    [[nodiscard]] static constexpr float256 denorm_min() noexcept {
        return from_bits(1, 0, 0, 0);
    }
    [[nodiscard]] static constexpr float256 max() noexcept {
        return from_bits(~0ULL, ~0ULL, ~0ULL,
                         (static_cast<std::uint64_t>(0x7FFFE) << 44) | detail::kFracHiMask);
    }
    [[nodiscard]] static constexpr float256 epsilon() noexcept {
        return from_bits(0, 0, 0, static_cast<std::uint64_t>(detail::kBias - (detail::kPrecision - 1)) << 44);
    }

private:
    static constexpr float256 add_sub(float256 a, float256 b, bool subtract) noexcept;
    static constexpr float256 mul(float256 a, float256 b) noexcept;
    static constexpr float256 div(float256 a, float256 b) noexcept;
};

constexpr float256 float256::add_sub(float256 a, float256 b, bool subtract) noexcept {
    using namespace detail;
    if (subtract) b.set_sign_bit(!b.sign_bit());
    Unpacked ua = a.unpack();
    Unpacked ub = b.unpack();

    if (ua.kind == Kind::nan) return pack(ua);
    if (ub.kind == Kind::nan) return pack(ub);
    if (ua.kind == Kind::inf && ub.kind == Kind::inf) {
        if (ua.sign != ub.sign) return pack(make_nan());
        return pack(ua);
    }
    if (ua.kind == Kind::inf) return pack(ua);
    if (ub.kind == Kind::inf) return pack(ub);
    if (ua.kind == Kind::zero && ub.kind == Kind::zero) {
        if (ua.sign == ub.sign) return zero(ua.sign);
        return zero(false);
    }
    if (ua.kind == Kind::zero) return pack(ub);
    if (ub.kind == Kind::zero) return pack(ua);

    if (ua.wexp < ub.wexp || (ua.wexp == ub.wexp && ua.sig < ub.sig)) {
        std::swap(ua, ub);
    }
    const int diff = ua.wexp - ub.wexp;
    bool sticky = ub.sticky;
    u256 sb = shr_sticky(ub.sig, diff, sticky);

    Unpacked r;
    if (ua.sign == ub.sign) {
        std::uint64_t carry = 0;
        r.sig = add_u256(ua.sig, sb, carry);
        r.wexp = ua.wexp;
        r.sign = ua.sign;
        r.sticky = sticky;
        r.kind = Kind::normal;
        if (carry) {
            sticky = sticky || r.sig.bit(0);
            r.sig = (r.sig >> 1);
            r.sig.set_bit(255);
            r.sticky = sticky;
            ++r.wexp;
        }
    } else {
        std::uint64_t br = 0;
        r.sig = sub_u256(ua.sig, sb, br);
        r.sign = ua.sign;
        r.wexp = ua.wexp;
        r.sticky = sticky;
        r.kind = Kind::normal;
        if (r.sig.is_zero()) return zero(false);
        const int lz = r.sig.clz();
        r.sig = r.sig << lz;
        r.wexp -= lz;
    }
    return pack(r);
}

constexpr float256 float256::mul(float256 a, float256 b) noexcept {
    using namespace detail;
    Unpacked ua = a.unpack();
    Unpacked ub = b.unpack();
    const bool sign = ua.sign ^ ub.sign;

    if (ua.kind == Kind::nan) return pack(ua);
    if (ub.kind == Kind::nan) return pack(ub);
    const bool a_inf = ua.kind == Kind::inf, b_inf = ub.kind == Kind::inf;
    const bool a_z   = ua.kind == Kind::zero, b_z   = ub.kind == Kind::zero;
    if ((a_inf && b_z) || (b_inf && a_z)) return pack(make_nan());
    if (a_inf || b_inf) return inf(sign);
    if (a_z || b_z) return zero(sign);

    const u512 p = mul_wide(ua.sig, ub.sig);
    Unpacked r;
    r.sign = sign;
    r.kind = Kind::normal;
    r.wexp = ua.wexp + ub.wexp;
    if (p.bit(511)) {
        r.sig = p.hi256();
        r.sticky = !p.lo256().is_zero();
        ++r.wexp;
    } else {
        const u512 s = p << 1;
        r.sig = s.hi256();
        r.sticky = !s.lo256().is_zero();
    }
    return pack(r);
}

constexpr float256 float256::div(float256 a, float256 b) noexcept {
    using namespace detail;
    Unpacked ua = a.unpack();
    Unpacked ub = b.unpack();
    const bool sign = ua.sign ^ ub.sign;

    if (ua.kind == Kind::nan) return pack(ua);
    if (ub.kind == Kind::nan) return pack(ub);
    if (ua.kind == Kind::inf && ub.kind == Kind::inf) return pack(make_nan());
    if (ub.kind == Kind::zero) {
        if (ua.kind == Kind::zero) return pack(make_nan());
        return inf(sign);
    }
    if (ua.kind == Kind::inf) return inf(sign);
    if (ua.kind == Kind::zero) return zero(sign);
    if (ub.kind == Kind::inf) return zero(sign);

    u512 num{};
    num.d[4] = ua.sig.d[0];
    num.d[5] = ua.sig.d[1];
    num.d[6] = ua.sig.d[2];
    num.d[7] = ua.sig.d[3];
    const Div512 dr = divmod_512_256(num, ub.sig);

    Unpacked r;
    r.sign = sign;
    r.kind = Kind::normal;
    r.wexp = ua.wexp - ub.wexp;
    r.sticky = !dr.rem.is_zero();
    if (dr.quot.bit(256)) {
        r.sticky = r.sticky || dr.quot.bit(0);
        r.sig = (dr.quot >> 1).lo256();
    } else {
        --r.wexp;
        r.sig = dr.quot.lo256();
        if (r.sig.bit(255) == false) {
            const int lz = r.sig.clz();
            r.sig = r.sig << lz;
            r.wexp -= lz;
        }
    }
    return pack(r);
}

constexpr float256 float256::fma(float256 a, float256 b, float256 c) noexcept {
    using namespace detail;
    Unpacked ua = a.unpack();
    Unpacked ub = b.unpack();
    Unpacked uc = c.unpack();
    const bool psign = ua.sign ^ ub.sign;

    if (ua.kind == Kind::nan) return pack(ua);
    if (ub.kind == Kind::nan) return pack(ub);
    if (uc.kind == Kind::nan) return pack(uc);

    const bool a_inf = ua.kind == Kind::inf, b_inf = ub.kind == Kind::inf;
    const bool a_z = ua.kind == Kind::zero, b_z = ub.kind == Kind::zero;
    if ((a_inf && b_z) || (b_inf && a_z)) return pack(make_nan());

    if (a_inf || b_inf) {
        if (uc.kind == Kind::inf && uc.sign != psign) return pack(make_nan());
        return inf(psign);
    }
    if (a_z || b_z) return pack(uc);

    const u512 p = mul_wide(ua.sig, ub.sig);
    int pexp = ua.wexp + ub.wexp;
    u512 psig = p;
    if (!psig.bit(511)) {
        psig = psig << 1;
    } else {
        ++pexp;
    }

    if (uc.kind == Kind::zero) {
        Unpacked r;
        r.sign = psign;
        r.kind = Kind::normal;
        r.wexp = pexp;
        r.sig = psig.hi256();
        r.sticky = !psig.lo256().is_zero();
        return pack(r);
    }
    if (uc.kind == Kind::inf) return pack(uc);

    u512 cs{};
    cs.d[4] = uc.sig.d[0];
    cs.d[5] = uc.sig.d[1];
    cs.d[6] = uc.sig.d[2];
    cs.d[7] = uc.sig.d[3];
    int cexp = uc.wexp;

    int rexp = pexp;
    u512 as = psig;
    u512 bs = cs;
    bool bsign = uc.sign;
    bool asticky = false, bsticky = false;
    if (cexp > pexp) {
        const int d = cexp - pexp;
        if (d >= 512) {
            asticky = !as.is_zero();
            as = {};
        } else {
            for (int i = 0; i < d; ++i) if (as.bit(i)) asticky = true;
            as = as >> d;
        }
        rexp = cexp;
    } else if (pexp > cexp) {
        const int d = pexp - cexp;
        if (d >= 512) {
            bsticky = !bs.is_zero();
            bs = {};
        } else {
            for (int i = 0; i < d; ++i) if (bs.bit(i)) bsticky = true;
            bs = bs >> d;
        }
    }

    Unpacked r;
    r.kind = Kind::normal;
    r.wexp = rexp;
    if (psign == bsign) {
        r.sign = psign;
        u512 s = add_u512(as, bs);
        r.sticky = asticky || bsticky;
        if (as.bit(511) && bs.bit(511)) {
            r.sticky = r.sticky || s.bit(0);
            s = s >> 1;
            s.set_bit(511);
            ++r.wexp;
        } else if (as.bit(511) != bs.bit(511) && !s.bit(511)) {
            r.sticky = r.sticky || s.bit(0);
            s = s >> 1;
            s.set_bit(511);
            ++r.wexp;
        }
        r.sig = s.hi256();
        r.sticky = r.sticky || !s.lo256().is_zero();
        return pack(r);
    }
    bool prod_ge = !(as < bs);
    u512 s = prod_ge ? sub_u512(as, bs) : sub_u512(bs, as);
    r.sign = prod_ge ? psign : bsign;
    r.sticky = asticky || bsticky;
    if (s.is_zero() && !r.sticky) return zero(false);
    const int lz = s.clz();
    s = s << lz;
    r.wexp -= lz;
    r.sig = s.hi256();
    r.sticky = r.sticky || !s.lo256().is_zero();
    return pack(r);
}

constexpr float256 float256::sqrt(float256 x) noexcept {
    using namespace detail;
    Unpacked u = x.unpack();
    if (u.kind == Kind::nan) return pack(u);
    if (u.kind == Kind::zero) return x;
    if (u.sign) return pack(make_nan());
    if (u.kind == Kind::inf) return x;

    int exp = u.wexp;
    u512 n{};
    n.d[0] = u.sig.d[0];
    n.d[1] = u.sig.d[1];
    n.d[2] = u.sig.d[2];
    n.d[3] = u.sig.d[3];
    n = n << 255;
    if (exp & 1) {
        n = n << 1;
        --exp;
    }
    const u256 s = isqrt_512(n);
    const u512 s2 = mul_wide(s, s);
    const bool inexact = !(s2 == n);

    Unpacked r;
    r.sign = false;
    r.kind = Kind::normal;
    r.wexp = exp / 2;
    r.sticky = inexact;
    const int lz = s.clz();
    r.sig = s << lz;
    r.wexp -= lz;
    return pack(r);
}

inline std::optional<float256> float256::parse(std::string_view s) {
    using namespace detail;
    std::size_t i = 0;
    while (i < s.size() && std::isspace(static_cast<unsigned char>(s[i]))) ++i;
    if (i >= s.size()) return std::nullopt;

    bool sign = false;
    if (s[i] == '+' || s[i] == '-') {
        sign = s[i] == '-';
        ++i;
    }
    if (i >= s.size()) return std::nullopt;

    auto ieq = [](char a, char b) {
        return std::tolower(static_cast<unsigned char>(a)) ==
               std::tolower(static_cast<unsigned char>(b));
    };

    if (i + 3 <= s.size() && ieq(s[i], 'i') && ieq(s[i + 1], 'n') && ieq(s[i + 2], 'f')) {
        return inf(sign);
    }
    if (i + 3 <= s.size() && ieq(s[i], 'n') && ieq(s[i + 1], 'a') && ieq(s[i + 2], 'n')) {
        return qnan();
    }

    if (i + 2 <= s.size() && s[i] == '0' && (s[i + 1] == 'x' || s[i + 1] == 'X')) {
        i += 2;
        u256 sig{};
        int sig_bits = 0;
        bool seen_digit = false;
        bool seen_dot = false;
        int frac_bits = 0;
        auto push_nibble = [&](int v) {
            if (sig_bits >= 256) return;
            sig = sig << 4;
            sig.d[0] |= static_cast<std::uint64_t>(v);
            sig_bits += 4;
            if (seen_dot) frac_bits += 4;
            seen_digit = true;
        };
        for (; i < s.size(); ++i) {
            const char c = s[i];
            int v = -1;
            if (c >= '0' && c <= '9') v = c - '0';
            else if (c >= 'a' && c <= 'f') v = 10 + c - 'a';
            else if (c >= 'A' && c <= 'F') v = 10 + c - 'A';
            else if (c == '.' && !seen_dot) { seen_dot = true; continue; }
            else break;
            if (v >= 0) push_nibble(v);
        }
        if (!seen_digit) return std::nullopt;
        if (i >= s.size() || (s[i] != 'p' && s[i] != 'P')) return std::nullopt;
        ++i;
        bool esign = false;
        if (i < s.size() && (s[i] == '+' || s[i] == '-')) { esign = s[i] == '-'; ++i; }
        int eacc = 0;
        bool eany = false;
        while (i < s.size() && s[i] >= '0' && s[i] <= '9') {
            eany = true;
            if (eacc < 10'000'000) eacc = eacc * 10 + (s[i] - '0');
            ++i;
        }
        if (!eany) return std::nullopt;
        const int exp2 = (esign ? -eacc : eacc) - frac_bits;
        if (sig.is_zero()) return zero(sign);
        const int lz = sig.clz();
        Unpacked u;
        u.sign = sign;
        u.kind = Kind::normal;
        u.sig = sig << lz;
        u.wexp = exp2 + (255 - lz);
        return pack(u);
    }

    std::string digits;
    int frac_count = 0;
    bool seen_dot = false;
    bool seen_digit = false;
    for (; i < s.size(); ++i) {
        const char c = s[i];
        if (c >= '0' && c <= '9') {
            digits.push_back(c);
            if (seen_dot) ++frac_count;
            seen_digit = true;
        } else if (c == '.' && !seen_dot) {
            seen_dot = true;
        } else {
            break;
        }
    }
    if (!seen_digit) return std::nullopt;

    int exp10 = 0;
    if (i < s.size() && (s[i] == 'e' || s[i] == 'E')) {
        ++i;
        bool esign = false;
        if (i < s.size() && (s[i] == '+' || s[i] == '-')) { esign = s[i] == '-'; ++i; }
        int acc = 0;
        bool any = false;
        while (i < s.size() && s[i] >= '0' && s[i] <= '9') {
            any = true;
            if (acc < 10'000'000) acc = acc * 10 + (s[i] - '0');
            ++i;
        }
        if (!any) return std::nullopt;
        exp10 = esign ? -acc : acc;
    }
    while (i < s.size() && std::isspace(static_cast<unsigned char>(s[i]))) ++i;
    if (i != s.size()) return std::nullopt;

    exp10 -= frac_count;
    std::size_t z = 0;
    while (z < digits.size() && digits[z] == '0') ++z;
    if (z == digits.size()) return zero(sign);
    digits = digits.substr(z);

    BigInt m = BigInt::from_dec(digits);
    constexpr int extra = 80;

    Unpacked u;
    u.sign = sign;
    u.kind = Kind::normal;

    if (exp10 >= 0) {
        m.mul_pow5(static_cast<unsigned>(exp10));
        const int bl = m.bit_length();
        if (bl == 0) return zero(sign);
        const int hi = bl - 1;
        u.sig = m.extract_u256(hi);
        bool sticky = false;
        for (int b = hi - 256; b >= 0; --b) {
            if (m.bit(b)) { sticky = true; break; }
        }
        u.wexp = exp10 + hi;
        u.sticky = sticky;
        return pack(u);
    }

    const int q = -exp10;
    BigInt fiveq{1};
    fiveq.mul_pow5(static_cast<unsigned>(q));
    const int mbl = m.bit_length();
    const int fbl = fiveq.bit_length();
    int shift = (256 + extra) + fbl - mbl;
    if (shift < extra) shift = extra;
    m.shl(shift);
    const BigInt rem = m.divmod(fiveq);
    bool sticky = !rem.is_zero();
    const int bl = m.bit_length();
    if (bl == 0) return zero(sign);
    const int hi = bl - 1;
    u.sig = m.extract_u256(hi);
    for (int b = hi - 256; b >= 0; --b) {
        if (m.bit(b)) { sticky = true; break; }
    }
    u.wexp = hi - shift - q;
    u.sticky = sticky;
    return pack(u);
}

inline std::string float256::to_string(int precision) const {
    using namespace detail;
    if (isnan()) return "nan";
    if (isinf()) return signbit() ? "-inf" : "inf";
    if (iszero()) return signbit() ? "-0" : "0";
    if (precision <= 0) precision = max_digits10;
    if (precision > 80) precision = 80;

    const Unpacked u = unpack();
    const std::uint64_t top = u.sig.d[3];
    const std::uint64_t frac53 = (top >> 11) & 0x000FFFFFFFFFFFFFULL;
    const std::uint64_t dbits = (1023ULL << 52) | frac53;
    const double mant = std::bit_cast<double>(dbits);
    const double log10v =
        static_cast<double>(u.wexp) * 0.30102999566398119521373889472449 +
        std::log10(mant);
    int n = static_cast<int>(std::floor(log10v));

    auto digits_at = [&](int k) -> std::pair<std::string, bool> {
        BigInt num = BigInt::from_u256(u.sig);
        const int two_exp = u.wexp - 255 - k;
        bool sticky = false;
        if (k >= 0) {
            BigInt den{1};
            den.mul_pow5(static_cast<unsigned>(k));
            if (two_exp >= 0) {
                num.shl(two_exp);
            } else {
                den.shl(-two_exp);
            }
            const BigInt rem = num.divmod(den);
            sticky = !rem.is_zero();
        } else {
            const int pk = -k;
            num.mul_pow5(static_cast<unsigned>(pk));
            const int te2 = u.wexp - 255 + pk;
            if (te2 >= 0) {
                num.shl(te2);
            } else {
                sticky = num.shr_sticky(-te2);
            }
        }
        return {num.to_dec(), sticky};
    };

    // The double log10 estimate is off by one near powers of ten (the top 53
    // significand bits can sit on the wrong side of an integer boundary).
    // Regenerate until floor(|v| / 10^(n-precision)) has exactly precision+1
    // digits, then round ties-to-even using the extra digit and any sticky bit.
    std::string core;
    for (int attempt = 0; attempt < 8; ++attempt) {
        auto [digs, sticky1] = digits_at(n - precision);
        if (digs == "0") {
            n -= precision > 4 ? precision : 4;
            continue;
        }
        const int want = precision + 1;
        const int got = static_cast<int>(digs.size());
        if (got != want) {
            n += got - want;
            continue;
        }
        const int last = digs.back() - '0';
        digs.pop_back();
        const bool odd = !digs.empty() && ((digs.back() - '0') & 1);
        const bool inc = last > 5 || (last == 5 && (sticky1 || odd));
        if (inc) {
            int i = static_cast<int>(digs.size()) - 1;
            while (i >= 0 && digs[static_cast<std::size_t>(i)] == '9') {
                digs[static_cast<std::size_t>(i)] = '0';
                --i;
            }
            if (i < 0) {
                digs.assign(static_cast<std::size_t>(precision), '0');
                digs[0] = '1';
                ++n;
            } else {
                ++digs[static_cast<std::size_t>(i)];
            }
        }
        core = std::move(digs);
        break;
    }
    if (core.empty()) core.assign(static_cast<std::size_t>(precision), '0');
    if (static_cast<int>(core.size()) > precision)
        core.resize(static_cast<std::size_t>(precision));
    while (static_cast<int>(core.size()) < precision) core.push_back('0');

    std::string out;
    if (signbit()) out.push_back('-');
    out.push_back(core[0]);
    if (precision > 1) {
        out.push_back('.');
        out.append(core, 1, std::string::npos);
    }
    if (precision == max_digits10 && out.find('.') != std::string::npos) {
        while (out.size() > 2 && out.back() == '0') out.pop_back();
        if (out.back() == '.') out.pop_back();
    }
    out.push_back('e');
    if (n >= 0) out.push_back('+');
    out += std::to_string(n);
    return out;
}

inline std::string float256::to_hex_string() const {
    using namespace detail;
    if (isnan()) return "nan";
    if (isinf()) return signbit() ? "-inf" : "inf";
    const Unpacked u = unpack();
    std::string out;
    if (u.sign) out.push_back('-');
    if (u.kind == Kind::zero) {
        out += "0x0p+0";
        return out;
    }
    out += "0x1.";
    const u256 f = u.sig << 1;
    char buf[64];
    for (int nibble = 0; nibble < 59; ++nibble) {
        int v = 0;
        for (int b = 0; b < 4; ++b) {
            const int bit = 255 - (nibble * 4 + b);
            if (f.bit(bit)) v |= 1 << (3 - b);
        }
        buf[nibble] = "0123456789abcdef"[v];
    }
    buf[59] = 0;
    out += buf;
    out.push_back('p');
    const int e = u.wexp;
    if (e >= 0) out.push_back('+');
    out += std::to_string(e);
    return out;
}

inline std::string float256::to_bit_string() const {
    std::string s(256, '0');
    for (int i = 0; i < 256; ++i) {
        const int word = i / 64;
        const int b = i % 64;
        if ((bits_[word] >> b) & 1ULL) s[255 - i] = '1';
    }
    return s;
}

[[nodiscard]] constexpr bool isnan(float256 x) noexcept { return x.isnan(); }
[[nodiscard]] constexpr bool isinf(float256 x) noexcept { return x.isinf(); }
[[nodiscard]] constexpr bool isfinite(float256 x) noexcept { return x.isfinite(); }
[[nodiscard]] constexpr bool isnormal(float256 x) noexcept { return x.isnormal(); }
[[nodiscard]] constexpr bool signbit(float256 x) noexcept { return x.signbit(); }
[[nodiscard]] constexpr int  fpclassify(float256 x) noexcept { return x.fpclassify(); }

[[nodiscard]] constexpr float256 abs(float256 x) noexcept {
    auto b = x.raw_bits();
    b[3] &= ~(1ULL << 63);
    return float256::from_bits(b);
}
[[nodiscard]] constexpr float256 fabs(float256 x) noexcept { return abs(x); }

[[nodiscard]] constexpr float256 copysign(float256 mag, float256 sgn) noexcept {
    auto b = mag.raw_bits();
    b[3] &= ~(1ULL << 63);
    if (sgn.signbit()) b[3] |= (1ULL << 63);
    return float256::from_bits(b);
}

[[nodiscard]] constexpr float256 fmin(float256 a, float256 b) noexcept {
    if (a.isnan()) return b;
    if (b.isnan()) return a;
    // IEEE 754: minNum of opposite-signed zeros is −0.
    if (a.iszero() && b.iszero())
        return float256::zero(a.signbit() || b.signbit());
    return a < b ? a : b;
}
[[nodiscard]] constexpr float256 fmax(float256 a, float256 b) noexcept {
    if (a.isnan()) return b;
    if (b.isnan()) return a;
    // IEEE 754: maxNum of opposite-signed zeros is +0.
    if (a.iszero() && b.iszero())
        return float256::zero(a.signbit() && b.signbit());
    return a > b ? a : b;
}

[[nodiscard]] constexpr float256 fma(float256 a, float256 b, float256 c) noexcept {
    return float256::fma(a, b, c);
}
[[nodiscard]] constexpr float256 sqrt(float256 x) noexcept { return float256::sqrt(x); }

[[nodiscard]] constexpr float256 ldexp(float256 x, int n) noexcept {
    using namespace detail;
    Unpacked u = x.unpack();
    if (u.kind == Kind::nan || u.kind == Kind::inf || u.kind == Kind::zero) return x;
    const std::int64_t ne = static_cast<std::int64_t>(u.wexp) + n;
    if (ne > kEmax + 16) return float256::inf(u.sign);
    if (ne < static_cast<std::int64_t>(kEmin) - 256) return float256::zero(u.sign);
    u.wexp = static_cast<int>(ne);
    return float256::from_unpacked(u);
}

[[nodiscard]] constexpr float256 scalbn(float256 x, int n) noexcept { return ldexp(x, n); }

[[nodiscard]] constexpr float256 frexp(float256 x, int* exp) noexcept {
    using namespace detail;
    Unpacked u = x.unpack();
    if (u.kind == Kind::nan || u.kind == Kind::inf) {
        if (exp) *exp = 0;
        return x;
    }
    if (u.kind == Kind::zero) {
        if (exp) *exp = 0;
        return x;
    }
    if (exp) *exp = u.wexp + 1;
    u.wexp = -1;
    return float256::from_unpacked(u);
}

[[nodiscard]] constexpr int ilogb(float256 x) noexcept {
    using namespace detail;
    Unpacked u = x.unpack();
    if (u.kind == Kind::nan) return FP_ILOGBNAN;
    if (u.kind == Kind::inf) return INT_MAX;
    if (u.kind == Kind::zero) return FP_ILOGB0;
    return u.wexp;
}

[[nodiscard]] constexpr float256 logb(float256 x) noexcept {
    using namespace detail;
    Unpacked u = x.unpack();
    if (u.kind == Kind::nan) return x;
    if (u.kind == Kind::inf) return float256::inf(false);
    if (u.kind == Kind::zero) return float256::inf(true);
    return float256::from_int(static_cast<std::int64_t>(u.wexp));
}

[[nodiscard]] constexpr float256 nextafter(float256 from, float256 to) noexcept {
    if (from.isnan()) return from;
    if (to.isnan()) return to;
    if (from == to) return to;
    const bool up = from < to;
    float256::storage_type b = from.raw_bits();
    if (from.iszero()) {
        return to.signbit() ? float256::from_bits(1, 0, 0, 1ULL << 63)
                            : float256::from_bits(1, 0, 0, 0);
    }
    const bool neg = from.signbit();
    const bool inc = (up && !neg) || (!up && neg);
    if (inc) {
        std::uint64_t c = 1;
        b[3] &= ~(1ULL << 63);
        for (int i = 0; i < 4; ++i) {
            const std::uint64_t s = b[i] + c;
            c = s < b[i] ? 1 : 0;
            b[i] = s;
            if (!c) break;
        }
        if (neg) b[3] |= (1ULL << 63);
    } else {
        b[3] &= ~(1ULL << 63);
        std::uint64_t br = 1;
        for (int i = 0; i < 4; ++i) {
            const std::uint64_t s = b[i] - br;
            br = b[i] < br ? 1 : 0;
            b[i] = s;
            if (!br) break;
        }
        if (neg) b[3] |= (1ULL << 63);
    }
    return float256::from_bits(b);
}

[[nodiscard]] constexpr float256 trunc(float256 x) noexcept {
    using namespace detail;
    Unpacked u = x.unpack();
    if (u.kind != Kind::normal && u.kind != Kind::subnormal) return x;
    if (u.wexp >= kPrecision - 1) return x;
    if (u.wexp < 0) return float256::zero(u.sign);
    const int frac_lsb = 255 - u.wexp;
    u.sig = (u.sig >> frac_lsb) << frac_lsb;
    u.sticky = false;
    return float256::from_unpacked(u);
}

[[nodiscard]] constexpr float256 floor(float256 x) noexcept {
    const float256 t = trunc(x);
    if (x.signbit() && t != x) return t - float256::one();
    return t;
}
[[nodiscard]] constexpr float256 ceil(float256 x) noexcept {
    const float256 t = trunc(x);
    if (!x.signbit() && t != x) return t + float256::one();
    return t;
}
[[nodiscard]] constexpr float256 round(float256 x) noexcept {
    using namespace detail;
    Unpacked u = x.unpack();
    if (u.kind != Kind::normal && u.kind != Kind::subnormal) return x;
    if (u.wexp >= kPrecision - 1) return x;
    if (u.wexp < -1) return float256::zero(u.sign);
    const float256 t = trunc(x);
    const float256 frac = abs(x - t);
    if (frac > float256::half() || (frac == float256::half())) {
        return x.signbit() ? t - float256::one() : t + float256::one();
    }
    return t;
}

// Round to nearest integer, ties to even (the default binary256 mode).
[[nodiscard]] constexpr float256 rint(float256 x) noexcept {
    using namespace detail;
    Unpacked u = x.unpack();
    if (u.kind != Kind::normal && u.kind != Kind::subnormal) return x;
    if (u.wexp >= kPrecision - 1) return x;
    if (u.wexp < -1) return float256::zero(u.sign);
    if (u.wexp < 0) {
        // |x| ∈ [0.5, 1). Exact ±0.5 ties to even, which is 0.
        const bool tie = u.wexp == -1 && (u.sig << 1).is_zero();
        if (tie) return float256::zero(u.sign);
        return u.sign ? -float256::one() : float256::one();
    }
    const int unit = 255 - u.wexp;
    const bool guard = u.sig.bit(unit - 1);
    bool sticky = false;
    for (int i = 0; i < unit - 1; ++i) {
        if (u.sig.bit(i)) { sticky = true; break; }
    }
    const bool odd = u.sig.bit(unit);
    u.sig = (u.sig >> unit) << unit;
    u.sticky = false;
    u.kind = Kind::normal;
    if (guard && (sticky || odd)) {
        u256 one{};
        one.set_bit(unit);
        std::uint64_t carry = 0;
        u.sig = add_u256(u.sig, one, carry);
        if (carry) {
            u.sig = {};
            u.sig.set_bit(255);
            ++u.wexp;
        }
    }
    return float256::from_unpacked(u);
}

// Integer part in *iptr, signed fraction returned. modf(±∞) yields ±0.
[[nodiscard]] constexpr float256 modf(float256 x, float256* iptr) noexcept {
    if (x.isnan()) {
        if (iptr) *iptr = x;
        return x;
    }
    if (x.isinf()) {
        if (iptr) *iptr = x;
        return float256::zero(x.signbit());
    }
    const float256 i = trunc(x);
    if (iptr) *iptr = i;
    if (x == i) return float256::zero(x.signbit());
    return x - i;
}

[[nodiscard]] constexpr float256 nextup(float256 x) noexcept {
    return nextafter(x, float256::inf(false));
}
[[nodiscard]] constexpr float256 nextdown(float256 x) noexcept {
    return nextafter(x, float256::inf(true));
}

[[nodiscard]] inline float256 hypot(float256 x, float256 y) noexcept {
    // IEEE 754: hypot(∞, qNaN) = hypot(qNaN, ∞) = +∞.
    if (x.isnan()) return y.isinf() ? float256::inf(false) : x;
    if (y.isnan()) return x.isinf() ? float256::inf(false) : y;
    if (x.isinf() || y.isinf()) return float256::inf(false);
    x = abs(x);
    y = abs(y);
    if (y > x) {
        const float256 t = x;
        x = y;
        y = t;
    }
    if (x.iszero()) return float256::zero();
    // Scale into [0.5, 1) so x² + y² cannot overflow or underflow.
    int ex = 0;
    const float256 xm = frexp(x, &ex);
    const float256 ym = ldexp(y, -ex);
    return ldexp(sqrt(xm * xm + ym * ym), ex);
}

} // namespace base256

template <>
class std::numeric_limits<base256::float256> {
public:
    static constexpr bool is_specialized = true;
    static constexpr bool is_signed = true;
    static constexpr bool is_integer = false;
    static constexpr bool is_exact = false;
    static constexpr bool has_infinity = true;
    static constexpr bool has_quiet_NaN = true;
    static constexpr bool has_signaling_NaN = true;
    static constexpr std::float_denorm_style has_denorm = std::denorm_present;
    static constexpr bool has_denorm_loss = false;
    static constexpr bool is_iec559 = true;
    static constexpr bool is_bounded = true;
    static constexpr bool is_modulo = false;
    static constexpr bool traps = false;
    static constexpr bool tinyness_before = false;
    static constexpr int digits = 237;
    static constexpr int digits10 = 71;
    static constexpr int max_digits10 = 73;
    static constexpr int radix = 2;
    static constexpr int min_exponent = -262141;
    static constexpr int min_exponent10 = -78913;
    static constexpr int max_exponent = 262144;
    static constexpr int max_exponent10 = 78912;
    static constexpr std::float_round_style round_style = std::round_to_nearest;

    static constexpr base256::float256 min() noexcept { return base256::float256::min(); }
    static constexpr base256::float256 max() noexcept { return base256::float256::max(); }
    static constexpr base256::float256 lowest() noexcept { return -base256::float256::max(); }
    static constexpr base256::float256 epsilon() noexcept { return base256::float256::epsilon(); }
    static constexpr base256::float256 round_error() noexcept { return base256::float256::half(); }
    static constexpr base256::float256 infinity() noexcept { return base256::float256::inf(); }
    static constexpr base256::float256 quiet_NaN() noexcept { return base256::float256::qnan(); }
    static constexpr base256::float256 signaling_NaN() noexcept { return base256::float256::snan(); }
    static constexpr base256::float256 denorm_min() noexcept { return base256::float256::denorm_min(); }
};

namespace base256::literals {
[[nodiscard]] inline float256 operator""_f256(const char* s) {
    auto r = float256::parse(s ? std::string_view(s) : std::string_view{});
    if (!r) throw std::invalid_argument("float256 literal");
    return *r;
}
[[nodiscard]] inline float256 operator""_f256(const char* s, std::size_t n) {
    auto r = float256::parse(std::string_view(s, n));
    if (!r) throw std::invalid_argument("float256 literal");
    return *r;
}
[[nodiscard]] constexpr float256 operator""_f256(long double x) {
    return float256(static_cast<double>(x));
}
} // namespace base256::literals
