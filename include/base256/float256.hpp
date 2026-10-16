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
#include <format>
#include <system_error>
#include <charconv>
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
#include <type_traits>
#include <utility>

// <cmath> (included above) guarantees isnan/isinf/signbit/... are functions,
// not macros, so the member and free functions below need no #undef tricks
// and user code keeps working whether it includes <math.h> or <cmath>.

// The rounding kernels (add_sub_m, mul_m, ... below) are templates on the mode
// source. GCC 14 kept add_sub_m<dyn_mode> out of line, so the operators lost
// the inlining that the 4.25.0 functions got (a dependent chain s = s + x * y
// became about 1.5x slower; CHANGELOG 4.26.0); forcing the kernels inline into
// their non-template callers (the operators and the per-call overloads)
// restores it. Clang 19 inlines as before without the attribute and was
// slightly slower with it, so it is GCC only.
#if defined(__GNUC__) && !defined(__clang__)
#define BASE256_KERNEL_INLINE [[gnu::always_inline]]
#else
#define BASE256_KERNEL_INLINE
#endif

namespace base256 {

class float256;

/// IEEE 754 rounding-direction attributes.
enum class round_mode : unsigned char {
    to_nearest_even,   ///< roundTiesToEven (default)
    to_nearest_away,   ///< roundTiesToAway
    toward_zero,       ///< roundTowardZero
    upward,            ///< roundTowardPositive
    downward,          ///< roundTowardNegative
};

/// IEEE 754 exception flags (bit mask).
enum fe_flag : unsigned {
    fe_invalid   = 1u << 0,
    fe_divbyzero = 1u << 1,
    fe_overflow  = 1u << 2,
    fe_underflow = 1u << 3,
    fe_inexact   = 1u << 4,
    fe_all       = 0x1Fu,
};

namespace detail {
struct Env {
    round_mode mode = round_mode::to_nearest_even;
    unsigned flags = 0;
};
inline thread_local Env tls_env{};

// During constant evaluation the environment is fixed: round-to-nearest-even,
// flags are not recorded.
[[nodiscard]] constexpr round_mode current_mode() noexcept {
    if consteval {
        return round_mode::to_nearest_even;
    } else {
        return tls_env.mode;
    }
}
constexpr void raise_flag(unsigned f) noexcept {
    if !consteval {
        tls_env.flags |= f;
    }
}
// The rounding mode of an operation: either an explicit round_mode (per-call
// overloads) or dyn_mode, which reads the thread's current mode at the point of
// rounding, exactly where the operators always read it.
struct dyn_mode {
    [[nodiscard]] constexpr operator round_mode() const noexcept { return current_mode(); }
};
} // namespace detail

/// Per-thread floating-point environment for float256 (like <cfenv>).
namespace fenv {
[[nodiscard]] inline round_mode get_round() noexcept { return detail::tls_env.mode; }
inline void set_round(round_mode m) noexcept { detail::tls_env.mode = m; }
/// Flags among `mask` that are currently raised.
[[nodiscard]] inline unsigned test_flags(unsigned mask = fe_all) noexcept { return detail::tls_env.flags & mask; }
inline void clear_flags(unsigned mask = fe_all) noexcept { detail::tls_env.flags &= ~mask; }
inline void raise_flags(unsigned mask) noexcept { detail::tls_env.flags |= (mask & fe_all); }
[[nodiscard]] inline unsigned get_flags() noexcept { return detail::tls_env.flags; }
inline void set_flags(unsigned f) noexcept { detail::tls_env.flags = f & fe_all; }

/// RAII: switch the rounding mode for a scope, restore it afterwards.
class scoped_round {
public:
    explicit scoped_round(round_mode m) noexcept : saved_(get_round()) { set_round(m); }
    ~scoped_round() { set_round(saved_); }
    scoped_round(const scoped_round&) = delete;
    scoped_round& operator=(const scoped_round&) = delete;
private:
    round_mode saved_;
};

/// RAII: run a computation in round-to-nearest-even with the caller's flags
/// saved and restored afterwards (used by the elementary functions).
class scoped_default_env {
public:
    scoped_default_env() noexcept : saved_(detail::tls_env) { detail::tls_env = detail::Env{}; }
    ~scoped_default_env() { detail::tls_env = saved_; }
    scoped_default_env(const scoped_default_env&) = delete;
    scoped_default_env& operator=(const scoped_default_env&) = delete;
private:
    detail::Env saved_;
};
} // namespace fenv

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

#if defined(__SIZEOF_INT128__)
__extension__ typedef unsigned __int128 u128_t;
__extension__ typedef __int128 i128_t;
template <class T> inline constexpr bool is_signed128_v = std::is_same_v<std::remove_cv_t<T>, i128_t>;
template <class T> inline constexpr bool is_int_like_v =
    (std::is_integral_v<T> && !std::is_same_v<std::remove_cv_t<T>, bool>) ||
    std::is_same_v<std::remove_cv_t<T>, i128_t> || std::is_same_v<std::remove_cv_t<T>, u128_t>;
#else
// No 128-bit integer type (e.g. MSVC): integer conversions top out at 64 bits.
typedef std::uint64_t u128_t;
typedef std::int64_t i128_t;
template <class T> inline constexpr bool is_signed128_v = false;
template <class T> inline constexpr bool is_int_like_v =
    std::is_integral_v<T> && !std::is_same_v<std::remove_cv_t<T>, bool>;
#endif
inline constexpr int kU128Bits = static_cast<int>(sizeof(u128_t) * 8);

struct DecFloor;

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
    friend constexpr float256 ldexp(float256, int, round_mode) noexcept;
    friend constexpr float256 frexp(float256, int*) noexcept;
    friend constexpr float256 trunc(float256) noexcept;
    friend constexpr float256 floor(float256) noexcept;
    friend constexpr float256 ceil(float256) noexcept;
    friend constexpr float256 round(float256) noexcept;
    friend constexpr float256 rint(float256) noexcept;
    friend constexpr float256 rint(float256, round_mode) noexcept;
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

    /// Round `sig` to a multiple of 2^lsb_pos (bits below are discarded,
    /// `sticky` stands for nonzero bits below bit 0). Returns true if the
    /// increment carried out of bit 255. Sets `inexact` when anything was
    /// discarded.
    static constexpr bool round_bits(detail::u256& sig, int lsb_pos, bool sticky, bool sign,
                                     round_mode mode, bool& inexact) noexcept {
        if (lsb_pos <= 0) { inexact = sticky; return false; }
        if (lsb_pos >= 256) {
            const bool g = sig.bit(lsb_pos - 1);
            bool rs = sticky;
            for (int i = 0; i < 4 && !rs; ++i) {
                std::uint64_t w = sig.d[i];
                if (lsb_pos == 256 && i == 3) w &= ~(1ULL << 63);  // exclude the guard bit
                rs = w != 0;
            }
            if (lsb_pos > 256) rs = rs || g;  // all of sig is below the guard
            const bool gg = lsb_pos == 256 ? g : false;
            inexact = gg || rs;
            bool inc = false;
            switch (mode) {
                case round_mode::to_nearest_even: inc = gg && rs; break;
                case round_mode::to_nearest_away: inc = gg; break;
                case round_mode::toward_zero:     inc = false; break;
                case round_mode::upward:          inc = inexact && !sign; break;
                case round_mode::downward:        inc = inexact && sign; break;
            }
            sig = {};
            return inc;
        }
        // Word/bit position of the result lsb and of the guard bit.
        const int w = lsb_pos >> 6, b = lsb_pos & 63;
        const int gp = lsb_pos - 1, gw = gp >> 6, gb = gp & 63;
        const bool g = ((sig.d[gw] >> gb) & 1ULL) != 0;
        bool rs = sticky;
        if (!rs) {
            for (int i = 0; i < gw; ++i) rs = rs || sig.d[i] != 0;
            if (gb) rs = rs || (sig.d[gw] & ((1ULL << gb) - 1ULL)) != 0;
        }
        inexact = g || rs;
        const bool lsb = ((sig.d[w] >> b) & 1ULL) != 0;
        bool inc = false;
        switch (mode) {
            case round_mode::to_nearest_even: inc = g && (rs || lsb); break;
            case round_mode::to_nearest_away: inc = g; break;
            case round_mode::toward_zero:     inc = false; break;
            case round_mode::upward:          inc = inexact && !sign; break;
            case round_mode::downward:        inc = inexact && sign; break;
        }
        // Clear the bits below lsb_pos.
        for (int i = 0; i < w; ++i) sig.d[i] = 0;
        if (b) sig.d[w] &= ~((1ULL << b) - 1ULL);
        if (!inc) return false;
        std::uint64_t add = 1ULL << b;
        for (int i = w; i < 4; ++i) {
            const std::uint64_t t = sig.d[i] + add;
            add = t < sig.d[i] ? 1ULL : 0ULL;
            sig.d[i] = t;
            if (!add) break;
        }
        return add != 0;
    }

    static constexpr bool round_ne(detail::u256& sig, int lsb_pos, bool sticky) noexcept {
        bool inexact = false;
        return round_bits(sig, lsb_pos, sticky, false, round_mode::to_nearest_even, inexact);
    }

    /// Overflow result for `mode`: infinity, or the largest finite value when
    /// the rounding direction points toward zero.
    static constexpr float256 overflow_result(bool sign, round_mode mode) noexcept {
        detail::raise_flag(fe_overflow | fe_inexact);
        const bool to_max = mode == round_mode::toward_zero ||
                            (mode == round_mode::upward && sign) ||
                            (mode == round_mode::downward && !sign);
        if (to_max) {
            float256 m = max();
            m.set_sign_bit(sign);
            return m;
        }
        return inf(sign);
    }

    /// Exact zero from cancelling operands: +0, or -0 when rounding downward.
    static constexpr float256 exact_zero_sum() noexcept {
        return zero(detail::current_mode() == round_mode::downward);
    }
    template <class M>
    static constexpr float256 exact_zero_sum(M m) noexcept {
        return zero(static_cast<round_mode>(m) == round_mode::downward);
    }

    /// Round and encode. Tininess is detected after rounding (as if the
    /// exponent range were unbounded); underflow is raised only when the
    /// result is tiny and inexact.
    static constexpr float256 pack(detail::Unpacked u) noexcept {
        return pack(u, detail::current_mode());
    }

    static constexpr float256 pack(detail::Unpacked u, round_mode mode) noexcept {
        using namespace detail;
#ifndef BASE256_NO_FAST_PATHS
        if (u.kind == Kind::normal && normal_range(u.wexp) && u.sig.bit(255))
            return pack_normal_range(u.sign, u.wexp, u.sig, u.sticky, mode);  // common case first
#endif
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
        if (biased >= kExpAllOnes) return overflow_result(u.sign, mode);

        u256 sig = u.sig;
        bool sticky = u.sticky;
        int wexp = u.wexp;
        bool inexact = false;

        if (biased <= 0) {
            // Tininess after rounding: round to full precision with an
            // unbounded exponent and see whether it is still below 2^Emin.
            bool tiny = true;
            {
                u256 t = sig;
                bool ix = false;
                const bool c = round_bits(t, kLsbPos, sticky, u.sign, mode, ix);
                tiny = !(c && wexp + 1 >= kEmin);
            }
            const std::int64_t rshift64 = static_cast<std::int64_t>(kLsbPos) - wexp + kEmin;
            if (rshift64 > 256) {
                // |x| < 2^-1 subnormal units: rounds to 0 or the smallest subnormal.
                const bool inc = (mode == round_mode::upward && !u.sign) ||
                                 (mode == round_mode::downward && u.sign);
                detail::raise_flag(fe_underflow | fe_inexact);
                return inc ? from_encoding(1, 0, 0, u.sign ? (1ULL << 63) : 0ULL) : zero(u.sign);
            }
            const int rshift = static_cast<int>(rshift64);
            if (rshift > 0) {
                const bool ovf = round_bits(sig, rshift, sticky, u.sign, mode, inexact);
                if (ovf) {
                    // Rounding carried out of bit 255: the result is exactly
                    // 2^(256 - rshift) subnormal units (min normal when that is
                    // bit 236, handled below).
                    sig = {};
                    sig.set_bit(256 - rshift);
                } else {
                    sig = sig >> rshift;
                }
            } else {
                inexact = sticky;
            }
            if (inexact) detail::raise_flag(fe_inexact | (tiny ? fe_underflow : 0u));
            if (sig.bit(kFracBits)) {
                return pack_normal(u.sign, kEmin, sig << kLsbPos);
            }
            return from_encoding(sig.d[0], sig.d[1], sig.d[2],
                                 (sig.d[3] & kFracHiMask) | (u.sign ? (1ULL << 63) : 0ULL));
        }

        const bool ovf = round_bits(sig, kLsbPos, sticky, u.sign, mode, inexact);
        if (ovf) {
            ++wexp;
            sig = {};
            sig.set_bit(255);
        }
        if (wexp > kEmax) return overflow_result(u.sign, mode);
        if (inexact) detail::raise_flag(fe_inexact);
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
    static constexpr float256 from_unpacked(detail::Unpacked u, round_mode m) noexcept { return pack(u, m); }

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
    // --- error-free transformations on the integer significands ------------------
    // Same results as the classic float formulas (s = RNE(a op b), e = RNE(exact - s)),
    // computed in one 512-bit pass. Return false (caller falls back to the float
    // formula) for non-finite/zero operands or when s overflows.
private:
    // RNE-round the signed 512-bit integer (neg, m) · 2^scale.
    static constexpr float256 round_u512(bool neg, const detail::u512& m, int scale) noexcept {
        using namespace detail;
        if (m.is_zero()) return zero(false);
        const int lz = m.clz();
        const u512 n = m << lz;
        Unpacked u;
        u.kind = Kind::normal;
        u.sign = neg;
        u.wexp = scale + 511 - lz;
        u.sig = n.hi256();
        u.sticky = !n.lo256().is_zero();
        return pack(u, round_mode::to_nearest_even);
    }
    // e = exact - s, with exact = (neg, m) · 2^scale and s its rounding (finite).
    static constexpr float256 eft_error(bool neg, const detail::u512& m, int scale, float256 s) noexcept {
        using namespace detail;
        if (s.iszero()) return round_u512(neg, m, scale);  // underflowed: e = RNE(exact)
        u512 si{};
        {
            const Unpacked us = s.unpack();
            const int sh = us.wexp - 255 - scale;  // s = sig · 2^(wexp-255)
            if (sh < 0) {
                // only possible when s is far below m's lowest bit: cannot happen for
                // a correct rounding of m; fall back to the conservative path
                si = u512(us.sig) >> (-sh);
            } else {
                si = u512(us.sig) << sh;  // mod 2^512: the error is < 2^511
            }
        }
        const u512 d = sub_u512(m, si);  // m - s_int (two's complement, |d| < 2^511)
        if (d.d[7] >> 63) {
            const u512 nd = sub_u512(u512{}, d);
            return round_u512(!neg, nd, scale);
        }
        return round_u512(neg, d, scale);
    }

    // RNE of a normalized significand (leading bit 255) with sticky, for a result
    // known to be normal: wexp in [kEmin, kEmax].
    static constexpr float256 pack_rne_normal(bool neg, int wexp, detail::u256 sig, bool sticky) noexcept {
        bool ix = false;
        if (round_bits(sig, detail::kLsbPos, sticky, neg, round_mode::to_nearest_even, ix)) {
            ++wexp;
            sig = {};
            sig.set_bit(255);
        }
        return pack_normal(neg, wexp, sig);
    }
    // s = RNE(m · 2^scale), e = RNE(m · 2^scale - s) directly on the 512-bit
    // significand. False when s or e is not comfortably normal (the caller then
    // takes the general path).
    static constexpr bool eft_round(bool neg, const detail::u512& m, int scale, float256& s, float256& e) noexcept {
        using namespace detail;
        const int lz = m.clz();
        if (lz >= 512) return false;
        const u512 n = m << lz;  // leading bit 511; value = n · 2^(scale - lz)
        const int wexp = scale + 511 - lz;
        if (wexp < kEmin + 300 || wexp >= kEmax) return false;
        constexpr int kCut = 256 + kLsbPos;  // bits below kCut are rounded off
        // L = low kCut bits of n.
        u512 low = n;
        for (int i = 7; i > kCut / 64; --i) low.d[i] = 0;
        low.d[kCut / 64] &= (1ULL << (kCut % 64)) - 1ULL;
        const bool g = ((n.d[(kCut - 1) / 64] >> ((kCut - 1) % 64)) & 1ULL) != 0;
        const bool lsb = ((n.d[kCut / 64] >> (kCut % 64)) & 1ULL) != 0;
        bool rest = (n.d[(kCut - 1) / 64] & ((1ULL << ((kCut - 1) % 64)) - 1ULL)) != 0;
        for (int i = 0; i < (kCut - 1) / 64 && !rest; ++i) rest = n.d[i] != 0;
        const bool inc = g && (rest || lsb);
        u256 top = n.hi256();
        top.d[0] &= ~((1ULL << kLsbPos) - 1ULL);
        int ws = wexp;
        if (inc) {
            std::uint64_t c = 0;
            top = add_u256(top, u256{1ULL << kLsbPos}, c);
            if (c) { ++ws; top = {}; top.set_bit(255); }
        }
        if (ws > kEmax) return false;
        s = pack_normal(neg, ws, top);
        // e = L - inc · 2^kCut, at scale (scale - lz).
        bool eneg = neg;
        u512 d = low;
        if (inc) {
            u512 big{};
            big.set_bit(kCut);
            d = sub_u512(big, low);
            eneg = !neg;
        }
        if (d.is_zero()) { e = zero(false); return true; }
        const int dz = d.clz();
        const int we = scale - lz + 511 - dz;
        if (we < kEmin + 1) return false;
        const u512 dn = d << dz;
        e = pack_rne_normal(eneg, we, dn.hi256(), !dn.lo256().is_zero());
        return true;
    }

    // --- fast paths for normal operands (bit-identical to the general path) -----
    [[nodiscard]] constexpr bool is_normal_enc() const noexcept {
        const int be = biased_exp();
        return be != 0 && be != detail::kExpAllOnes;
    }
    // Significand of a normal number with the leading bit at 255 (as unpack()).
    [[nodiscard]] constexpr detail::u256 normal_sig() const noexcept {
        const std::uint64_t w3 = (bits_[3] & detail::kFracHiMask) | (1ULL << 44);
        constexpr int s = detail::kLsbPos;
        return detail::u256::from_words(bits_[0] << s, (bits_[1] << s) | (bits_[0] >> (64 - s)),
                                        (bits_[2] << s) | (bits_[1] >> (64 - s)),
                                        (w3 << s) | (bits_[2] >> (64 - s)));
    }
    // pack() restricted to results with kEmin <= wexp < kEmax (normal before and
    // after rounding): the normal branch of pack(), without the special cases.
    static constexpr float256 pack_normal_range(bool sign, int wexp, detail::u256 sig, bool sticky,
                                                round_mode mode) noexcept {
        bool inexact = false;
        if (round_bits(sig, detail::kLsbPos, sticky, sign, mode, inexact)) {
            ++wexp;
            sig = {};
            sig.set_bit(255);
        }
        if (inexact) detail::raise_flag(fe_inexact);
        return pack_normal(sign, wexp, sig);
    }
    [[nodiscard]] static constexpr bool normal_range(int wexp) noexcept {
        return wexp >= detail::kEmin && wexp < detail::kEmax;
    }

    static constexpr bool finite_nonzero(const detail::Unpacked& u) noexcept {
        return u.kind == detail::Kind::normal || u.kind == detail::Kind::subnormal;
    }
    // Add sig · 2^(w - 255) (sign `neg`) into the two's-complement window W whose
    // bit 508 is 2^maxw (so the sum of four terms below 2^(maxw+1) cannot overflow).
    static constexpr void acc_term(detail::u512& W, const detail::Unpacked& u, int maxw) noexcept {
        using namespace detail;
        const int shift = 253 - (maxw - u.wexp);
        if (shift <= -256) return;
        const u512 t = shift >= 0 ? (u512(u.sig) << shift) : (u512(u.sig) >> (-shift));
        W = u.sign ? sub_u512(W, t) : add_u512(W, t);
    }

public:
    /// (hi, lo) = normalized dd sum of up to four finite terms, accurate to about
    /// 2^-505 relative (one truncation of the 512-bit window). False when an
    /// operand is not finite or the result is outside the comfortable range.
    static constexpr bool dd_add_parts(float256 a0, float256 a1, float256 b0, float256 b1, float256& hi,
                                       float256& lo) noexcept {
        using namespace detail;
        const Unpacked t[4] = {a0.unpack(), a1.unpack(), b0.unpack(), b1.unpack()};
        int maxw = 0;
        bool any = false;
        for (const Unpacked& u : t) {
            if (u.kind == Kind::zero) continue;
            if (!finite_nonzero(u)) return false;
            if (!any || u.wexp > maxw) maxw = u.wexp;
            any = true;
        }
        if (!any) return false;
        u512 W{};
        for (const Unpacked& u : t)
            if (u.kind != Kind::zero) acc_term(W, u, maxw);
        bool neg = false;
        if (W.d[7] >> 63) { W = sub_u512(u512{}, W); neg = true; }
        if (W.is_zero()) { hi = zero(false); lo = zero(false); return true; }
        return eft_round(neg, W, maxw - 508, hi, lo);
    }
    /// (hi, lo) = normalized dd product (a0 + a1)(b0 + b1), omitting a1·b1
    /// (about 2^-474 relative), as the classic dd product does.
    static constexpr bool dd_mul_parts(float256 a0, float256 a1, float256 b0, float256 b1, float256& hi,
                                       float256& lo) noexcept {
        using namespace detail;
        const Unpacked ua = a0.unpack(), ub = b0.unpack(), ual = a1.unpack(), ubl = b1.unpack();
        if (!finite_nonzero(ua) || !finite_nonzero(ub)) return false;
        if (!(finite_nonzero(ual) || ual.kind == Kind::zero) || !(finite_nonzero(ubl) || ubl.kind == Kind::zero))
            return false;
        // W in units of 2^(wa + wb - 508): the main product shifted right by 2.
        u512 W = mul_wide(ua.sig, ub.sig) >> 2;
        // big · small, where small (a low part) stands in for the high part `replaced`.
        auto cross = [&](const Unpacked& big, const Unpacked& small, const Unpacked& replaced) {
            if (small.kind == Kind::zero) return;
            const int sh = (replaced.wexp - small.wexp) + 2;
            if (sh >= 512) return;
            const u512 c = mul_wide(big.sig, small.sig) >> sh;
            W = small.sign != replaced.sign ? sub_u512(W, c) : add_u512(W, c);
        };
        cross(ua, ubl, ub);
        cross(ub, ual, ua);
        if (W.d[7] >> 63) return false;  // cannot happen for normalized dd inputs
        return eft_round(ua.sign != ub.sign, W, ua.wexp + ub.wexp - 508, hi, lo);
    }

    static constexpr bool eft_two_sum(float256 a, float256 b, float256& s, float256& e) noexcept {
        using namespace detail;
        Unpacked ua = a.unpack(), ub = b.unpack();
        if ((ua.kind != Kind::normal && ua.kind != Kind::subnormal) ||
            (ub.kind != Kind::normal && ub.kind != Kind::subnormal))
            return false;
        if (ua.wexp < ub.wexp || (ua.wexp == ub.wexp && ua.sig < ub.sig)) std::swap(ua, ub);
        const int diff = ua.wexp - ub.wexp;
        if (diff > 255) return false;
        const int scale = ua.wexp - 510;
        const u512 A = u512(ua.sig) << 255;
        const u512 B = u512(ub.sig) << (255 - diff);
        const u512 m = ua.sign == ub.sign ? add_u512(A, B) : sub_u512(A, B);
        if (m.is_zero()) { s = zero(false); e = zero(false); return true; }
        if (eft_round(ua.sign, m, scale, s, e)) return true;
        s = round_u512(ua.sign, m, scale);
        if (!s.isfinite()) return false;
        e = eft_error(ua.sign, m, scale, s);
        return true;
    }
    static constexpr bool eft_two_prod(float256 a, float256 b, float256& p, float256& e) noexcept {
        using namespace detail;
        const Unpacked ua = a.unpack(), ub = b.unpack();
        if ((ua.kind != Kind::normal && ua.kind != Kind::subnormal) ||
            (ub.kind != Kind::normal && ub.kind != Kind::subnormal))
            return false;
        const bool neg = ua.sign != ub.sign;
        const int scale = ua.wexp + ub.wexp - 510;
        const u512 m = mul_wide(ua.sig, ub.sig);
        if (eft_round(neg, m, scale, p, e)) return true;
        p = round_u512(neg, m, scale);
        if (!p.isfinite()) return false;
        e = eft_error(neg, m, scale, p);
        return true;
    }

    /// Round (hi + lo) · 2^k to binary256 in `mode` with one rounding, raising
    /// the flags as a single operation would. hi and lo are finite and form a
    /// normalized double-float256 (|lo| <= ulp(hi)/2). `tail` describes the
    /// exact value relative to hi + lo: 0 = exact, +1 = slightly above,
    /// -1 = slightly below, 2 = unknown but inexact (taken in the direction of
    /// lo, or away from zero when lo is 0). The tail only matters when hi + lo
    /// is itself representable at the rounding position.
    static constexpr float256 round_dd(float256 hi, float256 lo, int k, round_mode mode, int tail) noexcept {
        using namespace detail;
        Unpacked ua = hi.unpack();
        Unpacked ub = lo.unpack();
        if (ua.kind == Kind::zero) {
            ua = ub;
            ub = Unpacked{};
        }
        if (ua.kind == Kind::zero) return hi;  // exact zero
        const bool rsign = ua.sign;
        bool mag_up = true;  // the tail increases |value|
        if (tail == 2) mag_up = ub.kind == Kind::zero || ub.sign == rsign;
        else if (tail != 0) mag_up = (tail > 0) != rsign;
        ua.wexp += k;
        Unpacked r;
        r.kind = Kind::normal;
        r.sign = rsign;
        bool sticky = false;
        if (ub.kind == Kind::zero) {
            r.sig = ua.sig;
            r.wexp = ua.wexp;
        } else {
            ub.wexp += k;
            if (ua.wexp < ub.wexp || (ua.wexp == ub.wexp && ua.sig < ub.sig)) std::swap(ua, ub);
            const long long diff = static_cast<long long>(ua.wexp) - ub.wexp;
            u256 sb = shr_sticky(ub.sig, diff > 300 ? 300 : static_cast<int>(diff), sticky);
            r.sign = ua.sign;
            r.wexp = ua.wexp;
            if (ua.sign == ub.sign) {
                std::uint64_t carry = 0;
                r.sig = add_u256(ua.sig, sb, carry);
                if (carry) {
                    sticky = sticky || r.sig.bit(0);
                    r.sig = r.sig >> 1;
                    r.sig.set_bit(255);
                    ++r.wexp;
                }
            } else {
                std::uint64_t br = 0;
                r.sig = sub_u256(ua.sig, sb, br);
                if (sticky) {  // see add(): the truncated subtrahend borrows one unit
                    std::uint64_t br2 = 0;
                    r.sig = sub_u256(r.sig, u256{1}, br2);
                }
                if (r.sig.is_zero() && !sticky) {
                    if (tail == 0) return zero(mode == round_mode::downward);
                    r.sig.set_bit(0);  // degenerate: keep a nonzero trace of the tail
                    r.sign = tail > 0 ? false : true;
                }
                const int lz = r.sig.clz();
                r.sig = r.sig << lz;
                r.wexp -= lz;
            }
        }
        if (tail != 0 && !sticky) {
            if (!mag_up) {
                std::uint64_t br = 0;
                r.sig = sub_u256(r.sig, u256{1}, br);
                const int lz = r.sig.clz();
                r.sig = r.sig << lz;
                r.wexp -= lz;
            }
            sticky = true;
        }
        r.sticky = sticky;
        return pack(r, mode);
    }

    static constexpr float256 round_integral(float256 x, round_mode mode, bool signal) noexcept {
        using namespace detail;
        Unpacked u = x.unpack();
        if (u.kind == Kind::nan) {
            if (u.signaling) raise_flag(fe_invalid);
            return pack(u);
        }
        if (u.kind != Kind::normal && u.kind != Kind::subnormal) return x;
        if (u.wexp >= kPrecision - 1) return x;
        const int unit = 255 - u.wexp;            // bit position of 2^0
        u256 sig = u.sig;
        bool inexact = false;
        const bool carry = round_bits(sig, unit > 300 ? 300 : unit, false, u.sign, mode, inexact);
        if (inexact && signal) raise_flag(fe_inexact);
        if (unit >= 256) {
            if (!carry) return zero(u.sign);
            return u.sign ? -one() : one();
        }
        if (carry) {
            sig = {};
            sig.set_bit(255);
            ++u.wexp;
        }
        u.sig = sig;
        u.sticky = false;
        u.kind = Kind::normal;
        return pack(u, round_mode::to_nearest_even); // exact
    }
private:
public:
    /// Exact fmod (ieee == false, quotient truncated) or IEEE remainder
    /// (ieee == true, quotient rounded to nearest, ties to even). The result is
    /// always exact. If quo is non-null it receives the sign of x/y and the low
    /// 31 bits of the integral quotient (C remquo).
    static constexpr float256 remainder_impl(float256 x, float256 y, bool ieee, int* quo) noexcept {
        using namespace detail;
        if (quo) *quo = 0;
        Unpacked ux = x.unpack();
        Unpacked uy = y.unpack();
        if (ux.signaling || uy.signaling) raise_flag(fe_invalid);
        if (ux.kind == Kind::nan) return pack(ux);
        if (uy.kind == Kind::nan) return pack(uy);
        if (ux.kind == Kind::inf || uy.kind == Kind::zero) return invalid_nan();
        if (uy.kind == Kind::inf || ux.kind == Kind::zero) return x;
        const bool qneg = ux.sign != uy.sign;
        // Integer significands: value = m * 2^(wexp - 236), m < 2^237.
        u256 mx = ux.sig >> kLsbPos;
        const u256 my = uy.sig >> kLsbPos;
        const int ex = ux.wexp;
        const int ey = uy.wexp;
        u256 r{};
        std::uint64_t qlow = 0;
        if (ex < ey) {
            r = mx; // |x| < |y|: quotient 0
        } else {
            r = mx;
            std::uint64_t br = 0;
            if (!(r < my)) { r = sub_u256(r, my, br); qlow = 1; }
            for (int i = ex - ey; i > 0; --i) {
                r = r << 1;
                qlow <<= 1;
                if (!(r < my)) { r = sub_u256(r, my, br); qlow |= 1; }
            }
        }
        // Now |x| = q*|y| + r*2^(ey-236), 0 <= r < my.
        bool flip = false;
        if (ieee) {
            const int c = ex < ey - 1 ? -1 : 0; // |x| < |y|/2 when ex < ey - 1
            bool round_up = false;
            if (c == 0) {
                // compare 2r with my (r may be scaled by 2^(ex-ey) if ex < ey)
                if (ex >= ey) {
                    const u256 r2 = r << 1;
                    round_up = my < r2 || (r2 == my && (qlow & 1));
                } else { // ex == ey - 1: |x| = mx * 2^(ex-236), |y|/2 = my * 2^(ey-237)
                    round_up = my < mx; // ties (equal) go to even quotient 0
                }
            }
            if (round_up) {
                ++qlow;
                flip = true;
            }
        }
        if (quo) {
            const int q31 = static_cast<int>(qlow & 0x7FFFFFFFULL);
            *quo = qneg ? -q31 : q31;
        }
        Unpacked u;
        u.kind = Kind::normal;
        if (!flip) {
            if (r.is_zero()) return zero(ux.sign);
            u.sign = ux.sign;
            u.sig = r;
            u.wexp = (ex < ey ? ex : ey) + (255 - 236) - 19;
        } else {
            // remainder = r - |y| (negative relative to x's sign)
            std::uint64_t br = 0;
            if (ex >= ey) {
                u.sig = sub_u256(my, r, br);
                u.wexp = ey;
            } else {
                // |x| - |y| with ex == ey - 1: my*2 - mx in units of 2^(ex-236)
                u.sig = sub_u256(my << 1, mx, br);
                u.wexp = ex;
            }
            u.sign = !ux.sign;
            if (u.sig.is_zero()) return zero(ux.sign);
        }
        // u.sig is an integer in units of 2^(wexp - 236); normalize to bit 255.
        const int lz = u.sig.clz();
        u.wexp = u.wexp + (255 - lz) - 236;
        u.sig = u.sig << lz;
        return pack(u, round_mode::to_nearest_even); // exact
    }
private:
    static constexpr float256 invalid_nan() noexcept {
        detail::raise_flag(fe_invalid);
        return pack(make_nan());
    }

public:
    constexpr float256() noexcept = default;
    constexpr float256(const float256&) noexcept = default;
    constexpr float256(float256&&) noexcept = default;
    constexpr float256& operator=(const float256&) noexcept = default;
    constexpr float256& operator=(float256&&) noexcept = default;

    constexpr explicit float256(double v) noexcept { *this = from_double(v); }
    constexpr explicit float256(float v) noexcept { *this = from_double(static_cast<double>(v)); }
    explicit float256(long double v) noexcept { *this = from_long_double(v); }

    /// Any integer type (including __int128 where available). Exact up to
    /// 237 bits; wider values are rounded in the current mode.
    template <class T>
        requires(detail::is_int_like_v<T>)
    constexpr explicit float256(T v) noexcept {
        if constexpr (std::is_signed_v<T> || detail::is_signed128_v<T>) {
            const bool neg = v < 0;
            detail::u128_t mag = neg ? detail::u128_t(0) - static_cast<detail::u128_t>(v)
                                     : static_cast<detail::u128_t>(v);
            *this = from_u128(mag, neg);
        } else {
            *this = from_u128(static_cast<detail::u128_t>(v), false);
        }
    }

    [[nodiscard]] static constexpr float256 from_u128(detail::u128_t v, bool negative) noexcept {
        using namespace detail;
        if (v == 0) return zero(negative);
        Unpacked u;
        u.sign = negative;
        u.kind = Kind::normal;
        if constexpr (kU128Bits > 64) {
            u.sig = u256::from_words(static_cast<std::uint64_t>(v), static_cast<std::uint64_t>(v >> (kU128Bits / 2)), 0, 0);
        } else {
            u.sig = u256::from_words(static_cast<std::uint64_t>(v), 0, 0, 0);
        }
        const int lz = u.sig.clz();
        u.sig = u.sig << lz;
        u.wexp = 255 - lz;
        return pack(u);  // exact: at most 128 significant bits
    }

    /// IEEE convertToInteger: round in `mode` (default toward zero, like a C
    /// cast). NaN or out-of-range results raise invalid and saturate (NaN -> 0).
    /// With exact = true, a value change also raises inexact.
    template <class T>
        requires(detail::is_int_like_v<T>)
    [[nodiscard]] constexpr T to_integer(round_mode mode = round_mode::toward_zero,
                                         bool exact = false) const noexcept {
        using namespace detail;
        constexpr bool is_signed = std::is_signed_v<T> || is_signed128_v<T>;
        constexpr int bits = static_cast<int>(sizeof(T) * 8) - (is_signed ? 1 : 0);
        constexpr u128_t maxmag = (bits >= kU128Bits) ? ~u128_t(0) : ((u128_t(1) << bits) - 1);
        const T tmax = static_cast<T>(maxmag);
        const T tmin = is_signed ? static_cast<T>(-static_cast<T>(maxmag) - 1) : T(0);
        if (isnan()) {
            raise_flag(fe_invalid);
            return T(0);
        }
        if (isinf()) {
            raise_flag(fe_invalid);
            return signbit() ? tmin : tmax;
        }
        const float256 r = round_integral(*this, mode, exact);
        if (r.iszero()) return T(0);
        const Unpacked u = r.unpack();
        const bool neg = u.sign;
        // |r| = sig / 2^255 * 2^wexp, an integer
        if (u.wexp >= kU128Bits) {
            raise_flag(fe_invalid);
            return neg ? tmin : tmax;
        }
        const u256 iv = u.sig >> (255 - u.wexp);
        u128_t mag = static_cast<u128_t>(iv.d[0]);
        if constexpr (kU128Bits > 64) mag |= static_cast<u128_t>(iv.d[1]) << (kU128Bits / 2);
        if (!neg) {
            if (mag > maxmag) { raise_flag(fe_invalid); return tmax; }
            return static_cast<T>(mag);
        }
        if constexpr (!is_signed) {
            raise_flag(fe_invalid);
            return T(0);
        } else {
            if (mag > maxmag + 1) { raise_flag(fe_invalid); return tmin; }
            return static_cast<T>(u128_t(0) - mag);
        }
    }

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
    /// True for a signaling NaN (quiet bit, fraction bit 235, clear).
    [[nodiscard]] constexpr bool issignaling() const noexcept {
        return isnan() && (bits_[3] & (1ULL << 43)) == 0;
    }
    /// The same NaN with the quiet bit set (payload kept); other values unchanged.
    [[nodiscard]] constexpr float256 quieted() const noexcept {
        float256 r = *this;
        if (isnan()) r.bits_[3] |= (1ULL << 43);
        return r;
    }
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

    /// Exact conversion from long double (x87 80-bit, binary128, or binary64,
    /// whatever the platform uses): every finite long double is representable.
    /// NaN payloads are not carried over; the result is the default quiet NaN
    /// with the input's sign.
    [[nodiscard]] static float256 from_long_double(long double v) noexcept {
        using namespace detail;
        const bool sign = std::signbit(v);
        if (std::isnan(v)) {
            float256 r = qnan();
            r.set_sign_bit(sign);
            return r;
        }
        if (std::isinf(v)) return inf(sign);
        if (v == 0.0L) return zero(sign);
        int e = 0;
        long double m = std::frexp(std::fabs(v), &e); // m in [0.5, 1), exact
        Unpacked u;
        u.sign = sign;
        u.kind = Kind::normal;
        // Peel off 64 significand bits at a time; ldexp and the subtraction
        // are exact, and long double has at most 113 significand bits.
        for (int w = 3; w >= 0 && m != 0.0L; --w) {
            m = std::ldexp(m, 64);
            const auto chunk = static_cast<std::uint64_t>(m);
            m -= static_cast<long double>(chunk);
            u.sig.d[w] = chunk;
        }
        u.wexp = e - 1;
        return pack(u);
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

    struct Narrowed {
        int cls = 0;          // 0 zero, 1 finite nonzero, 2 infinity, 3 nan, 4 largest finite
        bool sign = false;
        detail::u256 sig{};   // rounded, leading 1 at bit 255
        int wexp = 0;
    };
    /// Round a finite/special value to a binary format with precision P and
    /// normal exponent range [emin, emax], honouring the rounding mode and
    /// raising inexact/underflow/overflow.
    [[nodiscard]] constexpr Narrowed narrow_to(int P, int emin, int emax) const noexcept {
        return narrow_to(P, emin, emax, detail::dyn_mode{});
    }
    template <class M>
    [[nodiscard]] constexpr Narrowed narrow_to(int P, int emin, int emax, M m) const noexcept {
        using namespace detail;
        const Unpacked u = unpack();
        Narrowed r;
        r.sign = u.sign;
        if (u.kind == Kind::nan) {
            if (u.signaling) raise_flag(fe_invalid);
            r.cls = 3;
            return r;
        }
        if (u.kind == Kind::inf) { r.cls = 2; return r; }
        if (u.kind == Kind::zero) { r.cls = 0; return r; }
        round_mode mode = m; // not const: a const integral/enum local would be
        // manifestly constant-evaluated, making `if consteval` take the compile-time branch
        int keep = P;
        if (u.wexp < emin) keep = P - (emin - u.wexp);
        if (keep < 0) {
            raise_flag(fe_underflow | fe_inexact);
            const bool inc = (mode == round_mode::upward && !u.sign) ||
                             (mode == round_mode::downward && u.sign);
            if (!inc) { r.cls = 0; return r; }
            r.cls = 1;
            r.sig.set_bit(255);
            r.wexp = emin - P + 1;
            return r;
        }
        bool tiny = false;
        if (u.wexp < emin) {
            u256 t = u.sig;
            bool ix = false;
            const bool c = round_bits(t, 256 - P, u.sticky, u.sign, mode, ix);
            tiny = !(c && u.wexp + 1 >= emin);
        }
        u256 sig = u.sig;
        int wexp = u.wexp;
        bool inexact = false;
        if (round_bits(sig, 256 - keep, u.sticky, u.sign, mode, inexact)) {
            sig = {};
            sig.set_bit(255);
            ++wexp;
        }
        if (wexp > emax) {
            raise_flag(fe_overflow | fe_inexact);
            const bool to_max = mode == round_mode::toward_zero ||
                                (mode == round_mode::upward && u.sign) ||
                                (mode == round_mode::downward && !u.sign);
            r.cls = to_max ? 4 : 2;
            return r;
        }
        if (inexact) raise_flag(fe_inexact | (tiny ? fe_underflow : 0u));
        if (sig.is_zero()) { r.cls = 0; return r; }
        r.cls = 1;
        r.sig = sig;
        r.wexp = wexp;
        return r;
    }

    /// Correctly rounded conversion (current rounding mode, flags raised).
    [[nodiscard]] long double to_long_double() const noexcept { return to_long_double_m(detail::dyn_mode{}); }
    /// The same, rounded in `m` instead of the current mode.
    [[nodiscard]] long double to_long_double(round_mode m) const noexcept { return to_long_double_m(m); }

private:
    template <class M>
    [[nodiscard]] long double to_long_double_m(M m) const noexcept {
        using L = std::numeric_limits<long double>;
        const Narrowed n = narrow_to(L::digits, L::min_exponent - 1, L::max_exponent - 1, m);
        switch (n.cls) {
            case 0: return n.sign ? -0.0L : 0.0L;
            case 2: return n.sign ? -L::infinity() : L::infinity();
            case 3: return L::quiet_NaN();
            case 4: return n.sign ? -L::max() : L::max();
            default: break;
        }
        // At most L::digits significant bits remain, so every partial sum is exact.
        long double r = 0.0L;
        for (int w = 0; w < 4; ++w)
            r += std::ldexp(static_cast<long double>(n.sig.d[w]), 64 * w - 255 + n.wexp);
        return n.sign ? -r : r;
    }

    template <class M>
    [[nodiscard]] constexpr double to_double_m(M m) const noexcept {
        const Narrowed n = narrow_to(53, -1022, 1023, m);
        const std::uint64_t sb = n.sign ? (1ULL << 63) : 0ULL;
        switch (n.cls) {
            case 0: return std::bit_cast<double>(sb);
            case 2: return std::bit_cast<double>(sb | 0x7FF0000000000000ULL);
            case 3: return std::numeric_limits<double>::quiet_NaN();
            case 4: return std::bit_cast<double>(sb | 0x7FEFFFFFFFFFFFFFULL);
            default: break;
        }
        if (n.wexp >= -1022) {
            const std::uint64_t mant = (n.sig.d[3] >> 11) & 0x000FFFFFFFFFFFFFULL;
            const auto be = static_cast<std::uint64_t>(n.wexp + 1023);
            return std::bit_cast<double>(sb | (be << 52) | mant);
        }
        const int sh = 11 + (-1022 - n.wexp); // <= 63, exact after narrowing
        return std::bit_cast<double>(sb | (n.sig.d[3] >> sh));
    }

    template <class M>
    [[nodiscard]] constexpr float to_float_m(M m) const noexcept {
        const Narrowed n = narrow_to(24, -126, 127, m);
        const std::uint32_t sb = n.sign ? (1U << 31) : 0U;
        switch (n.cls) {
            case 0: return std::bit_cast<float>(sb);
            case 2: return std::bit_cast<float>(sb | 0x7F800000U);
            case 3: return std::numeric_limits<float>::quiet_NaN();
            case 4: return std::bit_cast<float>(sb | 0x7F7FFFFFU);
            default: break;
        }
        const auto top = static_cast<std::uint32_t>(n.sig.d[3] >> 32);
        if (n.wexp >= -126) {
            const auto be = static_cast<std::uint32_t>(n.wexp + 127);
            return std::bit_cast<float>(sb | (be << 23) | ((top >> 8) & 0x7FFFFFU));
        }
        const int sh = 8 + (-126 - n.wexp);
        return std::bit_cast<float>(sb | (top >> sh));
    }

public:
    /// Correctly rounded conversion (current rounding mode, flags raised).
    [[nodiscard]] constexpr double to_double() const noexcept { return to_double_m(detail::dyn_mode{}); }
    /// The same, rounded in `m` instead of the current mode (also in constant expressions).
    [[nodiscard]] constexpr double to_double(round_mode m) const noexcept { return to_double_m(m); }
    /// Correctly rounded (single rounding, current mode).
    [[nodiscard]] constexpr float to_float() const noexcept { return to_float_m(detail::dyn_mode{}); }
    /// The same, rounded in `m` instead of the current mode (also in constant expressions).
    [[nodiscard]] constexpr float to_float(round_mode m) const noexcept { return to_float_m(m); }

    constexpr void signal_unordered(const float256& o) const noexcept {
        if (isnan() || o.isnan()) detail::raise_flag(fe_invalid);
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

    // == and != are quiet predicates (invalid only for a signaling NaN);
    // <, <=, >, >= signal invalid for any NaN operand (IEEE 754 5.11).
    [[nodiscard]] constexpr bool operator==(const float256& o) const noexcept {
        if (issignaling() || o.issignaling()) detail::raise_flag(fe_invalid);
        if (isnan() || o.isnan()) return false;
        if (iszero() && o.iszero()) return true;
        return bits_ == o.bits_;
    }
    [[nodiscard]] constexpr bool operator!=(const float256& o) const noexcept { return !(*this == o); }
    [[nodiscard]] constexpr bool operator<(const float256& o)  const noexcept { signal_unordered(o); return (*this <=> o) == std::partial_ordering::less; }
    [[nodiscard]] constexpr bool operator>(const float256& o)  const noexcept { signal_unordered(o); return (*this <=> o) == std::partial_ordering::greater; }
    [[nodiscard]] constexpr bool operator<=(const float256& o) const noexcept {
        signal_unordered(o);
        auto c = *this <=> o;
        return c == std::partial_ordering::less || c == std::partial_ordering::equivalent;
    }
    [[nodiscard]] constexpr bool operator>=(const float256& o) const noexcept {
        signal_unordered(o);
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

    /// Scientific notation. precision > 0: that many significant digits,
    /// rounded in the current mode. precision == 0: the shortest string that
    /// parses back to the same value (round-to-nearest-even parsing).
    [[nodiscard]] std::string to_string(int precision = 0) const;
    /// to_string with the decimal rounding in `m` instead of the current mode.
    [[nodiscard]] std::string to_string(int precision, round_mode m) const;
    /// C99 hex-float, full 236-bit fraction (0x1.<59 hex digits>p<exp>).
    [[nodiscard]] std::string to_hex_string() const;
    /// std::to_chars-style formatting. precision < 0 means shortest round-trip.
    ///  scientific: d.ddde+XX, `precision` digits after the point
    ///  fixed:      ddd.ddd, `precision` digits after the point
    ///  general:    like %g with `precision` significant digits (0 -> 1); for
    ///              shortest, whichever of fixed/scientific is shorter
    ///  hex:        h.hhhp+X (no 0x prefix), `precision` hex digits
    [[nodiscard]] std::string format(std::chars_format fmt, int precision = -1) const;
    /// format with the decimal rounding in `m` instead of the current mode.
    [[nodiscard]] std::string format(std::chars_format fmt, int precision, round_mode m) const;

    // Decimal digit generation (exact, big-integer based).
    struct DecimalDigits {
        std::string digits;   // significant digits, no leading zeros
        int exp10 = 0;        // value = d1.d2d3... * 10^exp10
    };
    /// Exactly `ndig` significant digits rounded in `mode` (flags: inexact).
    [[nodiscard]] DecimalDigits decimal_significant(int ndig, round_mode mode) const;
    /// floor(|x|) to exactly `ndig` significant digits, plus a sticky flag.
    [[nodiscard]] detail::DecFloor decimal_floor(int ndig) const;
    /// The shortest digit string that round-trips under ties-to-even parsing.
    [[nodiscard]] DecimalDigits decimal_shortest() const;
    /// Digits of round(|x| * 10^frac) in `mode` (no leading zeros; "0" if zero).
    [[nodiscard]] std::string decimal_fixed(int frac, round_mode mode) const;
    [[nodiscard]] std::string to_bit_string() const;
    [[nodiscard]] static std::optional<float256> parse(std::string_view s);
    /// parse with the result rounded in `m` instead of the current mode.
    [[nodiscard]] static std::optional<float256> parse(std::string_view s, round_mode m);

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

public:
    // The operations with their rounding mode as a parameter: M is round_mode
    // (the per-call overloads add(a, b, m), ...) or detail::dyn_mode (the
    // operators, which use the thread's current mode). Not for direct use.
    template <class M> BASE256_KERNEL_INLINE static constexpr float256 add_sub_m(float256 a, float256 b, bool subtract, M m) noexcept;
    template <class M> BASE256_KERNEL_INLINE static constexpr float256 mul_m(float256 a, float256 b, M m) noexcept;
    template <class M> BASE256_KERNEL_INLINE static constexpr float256 div_m(float256 a, float256 b, M m) noexcept;
    template <class M> BASE256_KERNEL_INLINE static constexpr float256 fma_m(float256 a, float256 b, float256 c, M m) noexcept;
    template <class M> BASE256_KERNEL_INLINE static constexpr float256 sqrt_m(float256 x, M m) noexcept;
};

template <class M>
constexpr float256 float256::add_sub_m(float256 a, float256 b, bool subtract, M m) noexcept {
    using namespace detail;
    if (subtract) b.set_sign_bit(!b.sign_bit());
#ifndef BASE256_NO_FAST_PATHS
    if (a.is_normal_enc() && b.is_normal_enc()) {
        // Fast path: both normal. Same arithmetic as below; falls back to pack()
        // only when the result may be subnormal, zero or overflow.
        int ea = a.biased_exp(), eb = b.biased_exp();
        u256 sa = a.normal_sig(), sbig = b.normal_sig();
        bool na = a.sign_bit(), nb = b.sign_bit();
        if (ea < eb || (ea == eb && sa < sbig)) {
            std::swap(ea, eb);
            std::swap(sa, sbig);
            std::swap(na, nb);
        }
        bool sticky = false;
        const u256 sb = shr_sticky(sbig, ea - eb, sticky);
        int wexp = ea - kBias;
        u256 sig;
        if (na == nb) {
            std::uint64_t carry = 0;
            sig = add_u256(sa, sb, carry);
            if (carry) {
                sticky = sticky || sig.bit(0);
                sig = sig >> 1;
                sig.set_bit(255);
                ++wexp;
            }
        } else {
            std::uint64_t br = 0;
            sig = sub_u256(sa, sb, br);
            if (sticky) {
                std::uint64_t br2 = 0;
                sig = sub_u256(sig, u256{1}, br2);
            }
            if (sig.is_zero() && !sticky) return exact_zero_sum(m);
            const int lz = sig.clz();
            sig = sig << lz;
            wexp -= lz;
        }
        if (normal_range(wexp)) return pack_normal_range(na, wexp, sig, sticky, m);
        Unpacked r;
        r.sign = na;
        r.kind = Kind::normal;
        r.wexp = wexp;
        r.sig = sig;
        r.sticky = sticky;
        return pack(r, m);
    }
#endif
    Unpacked ua = a.unpack();
    Unpacked ub = b.unpack();

    if (ua.signaling || ub.signaling) detail::raise_flag(fe_invalid);
    if (ua.kind == Kind::nan) return pack(ua, m);
    if (ub.kind == Kind::nan) return pack(ub, m);
    if (ua.kind == Kind::inf && ub.kind == Kind::inf) {
        if (ua.sign != ub.sign) return invalid_nan();
        return pack(ua, m);
    }
    if (ua.kind == Kind::inf) return pack(ua, m);
    if (ub.kind == Kind::inf) return pack(ub, m);
    if (ua.kind == Kind::zero && ub.kind == Kind::zero) {
        if (ua.sign == ub.sign) return zero(ua.sign);
        return exact_zero_sum(m);
    }
    if (ua.kind == Kind::zero) return pack(ub, m);
    if (ub.kind == Kind::zero) return pack(ua, m);

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
        // The aligned subtrahend was truncated: its true value is sb + eps with
        // 0 < eps < 1 (in units of bit 0). The true difference is therefore
        // (r.sig - 1) + (1 - eps), so borrow one unit and keep sticky set.
        // Rounding the un-borrowed value would treat results just below a
        // rounding midpoint as being above it.
        if (sticky) {
            std::uint64_t br2 = 0;
            r.sig = sub_u256(r.sig, u256{1}, br2);
        }
        r.sign = ua.sign;
        r.wexp = ua.wexp;
        r.sticky = sticky;
        r.kind = Kind::normal;
        if (r.sig.is_zero() && !sticky) return exact_zero_sum(m);
        const int lz = r.sig.clz();
        r.sig = r.sig << lz;
        r.wexp -= lz;
    }
    return pack(r, m);
}

template <class M>
constexpr float256 float256::mul_m(float256 a, float256 b, M m) noexcept {
    using namespace detail;
#ifndef BASE256_NO_FAST_PATHS
    if (a.is_normal_enc() && b.is_normal_enc()) {
        const bool sign = a.sign_bit() != b.sign_bit();
        const u512 p = mul_wide(a.normal_sig(), b.normal_sig());
        int wexp = (a.biased_exp() - kBias) + (b.biased_exp() - kBias);
        u256 sig;
        bool sticky;
        if (p.bit(511)) {
            sig = p.hi256();
            sticky = !p.lo256().is_zero();
            ++wexp;
        } else {
            const u512 q = p << 1;
            sig = q.hi256();
            sticky = !q.lo256().is_zero();
        }
        if (normal_range(wexp)) return pack_normal_range(sign, wexp, sig, sticky, m);
        Unpacked r;
        r.sign = sign;
        r.kind = Kind::normal;
        r.wexp = wexp;
        r.sig = sig;
        r.sticky = sticky;
        return pack(r, m);
    }
#endif
    Unpacked ua = a.unpack();
    Unpacked ub = b.unpack();
    const bool sign = ua.sign ^ ub.sign;

    if (ua.signaling || ub.signaling) detail::raise_flag(fe_invalid);
    if (ua.kind == Kind::nan) return pack(ua, m);
    if (ub.kind == Kind::nan) return pack(ub, m);
    const bool a_inf = ua.kind == Kind::inf, b_inf = ub.kind == Kind::inf;
    const bool a_z   = ua.kind == Kind::zero, b_z   = ub.kind == Kind::zero;
    if ((a_inf && b_z) || (b_inf && a_z)) return invalid_nan();
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
    return pack(r, m);
}

template <class M>
constexpr float256 float256::div_m(float256 a, float256 b, M m) noexcept {
    using namespace detail;
    Unpacked ua = a.unpack();
    Unpacked ub = b.unpack();
    const bool sign = ua.sign ^ ub.sign;

    if (ua.signaling || ub.signaling) detail::raise_flag(fe_invalid);
    if (ua.kind == Kind::nan) return pack(ua, m);
    if (ub.kind == Kind::nan) return pack(ub, m);
    if (ua.kind == Kind::inf && ub.kind == Kind::inf) return invalid_nan();
    if (ub.kind == Kind::zero) {
        if (ua.kind == Kind::zero) return invalid_nan();
        if (ua.kind != Kind::inf) detail::raise_flag(fe_divbyzero);
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
    return pack(r, m);
}

template <class M>
constexpr float256 float256::fma_m(float256 a, float256 b, float256 c, M m) noexcept {
    using namespace detail;
    Unpacked ua = a.unpack();
    Unpacked ub = b.unpack();
    Unpacked uc = c.unpack();
    const bool psign = ua.sign ^ ub.sign;

    if (ua.signaling || ub.signaling || uc.signaling) detail::raise_flag(fe_invalid);
    if (ua.kind == Kind::nan) return pack(ua, m);
    if (ub.kind == Kind::nan) return pack(ub, m);
    if (uc.kind == Kind::nan) return pack(uc, m);

    const bool a_inf = ua.kind == Kind::inf, b_inf = ub.kind == Kind::inf;
    const bool a_z = ua.kind == Kind::zero, b_z = ub.kind == Kind::zero;
    if ((a_inf && b_z) || (b_inf && a_z)) return invalid_nan();

    if (a_inf || b_inf) {
        if (uc.kind == Kind::inf && uc.sign != psign) return invalid_nan();
        return inf(psign);
    }
    if (a_z || b_z) {
        // (±0) + c: exact. Two zeros of opposite sign sum to +0.
        if (uc.kind == Kind::zero) return psign == uc.sign ? zero(psign) : exact_zero_sum(m);
        return pack(uc, m);
    }

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
        return pack(r, m);
    }
    if (uc.kind == Kind::inf) return pack(uc, m);

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
        return pack(r, m);
    }
    bool prod_ge = !(as < bs);
    u512 s = prod_ge ? sub_u512(as, bs) : sub_u512(bs, as);
    r.sign = prod_ge ? psign : bsign;
    r.sticky = asticky || bsticky;
    // Only the smaller (shifted) operand can have lost bits; its true value
    // exceeds the truncated one, so borrow one unit (see add_sub).
    if (r.sticky) {
        u512 one{};
        one.d[0] = 1;
        s = sub_u512(s, one);
    }
    if (s.is_zero() && !r.sticky) return exact_zero_sum(m);
    const int lz = s.clz();
    s = s << lz;
    r.wexp -= lz;
    r.sig = s.hi256();
    r.sticky = r.sticky || !s.lo256().is_zero();
    return pack(r, m);
}

template <class M>
constexpr float256 float256::sqrt_m(float256 x, M m) noexcept {
    using namespace detail;
    Unpacked u = x.unpack();
    if (u.signaling) detail::raise_flag(fe_invalid);
    if (u.kind == Kind::nan) return pack(u, m);
    if (u.kind == Kind::zero) return x;
    if (u.sign) return invalid_nan();
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
    return pack(r, m);
}

// The operators: the thread's current mode, read where the result is rounded.
constexpr float256 float256::add_sub(float256 a, float256 b, bool subtract) noexcept {
    return add_sub_m(a, b, subtract, detail::dyn_mode{});
}
constexpr float256 float256::mul(float256 a, float256 b) noexcept { return mul_m(a, b, detail::dyn_mode{}); }
constexpr float256 float256::div(float256 a, float256 b) noexcept { return div_m(a, b, detail::dyn_mode{}); }
constexpr float256 float256::fma(float256 a, float256 b, float256 c) noexcept {
    return fma_m(a, b, c, detail::dyn_mode{});
}
constexpr float256 float256::sqrt(float256 x) noexcept { return sqrt_m(x, detail::dyn_mode{}); }

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

    // Far outside the binary256 range: avoid huge big-integer powers of ten.
    // |x| >= 10^(exp10 + nd - 1) and < 10^(exp10 + nd).
    {
        const long long mag = static_cast<long long>(exp10) + static_cast<long long>(digits.size());
        if (mag > 78915) {           // > 10^78914 > max()
            Unpacked big;
            big.sign = sign;
            big.kind = Kind::normal;
            big.sig.set_bit(255);
            big.wexp = kEmax + 8;
            return pack(big);
        }
        if (mag < -78985) {          // < 10^-78985 < denorm_min()/2
            Unpacked tiny;
            tiny.sign = sign;
            tiny.kind = Kind::normal;
            tiny.sig.set_bit(255);
            tiny.wexp = kEmin - 1000;
            tiny.sticky = true;
            return pack(tiny);
        }
    }

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

namespace detail {

// Decimal conversion at large decimal exponents. The exact path multiplies or
// divides by 5^|k| as a BigInt (|k| up to ~79000, i.e. 180000-bit numbers, built
// limb by limb: ~10 ms per power). For |k| >= kDecFastK the paths below use a
// 640-bit lower bound of 5^|k| instead and only answer when that decides the
// result exactly; otherwise (a value within ~2^-470 relative of a decimal
// boundary, or more than ~140 digits) they report failure and the caller takes
// the exact path. Exact decimal ties cannot occur there: c·10^q = B·2^E with
// c < 10^140, B < 2^256 would need 5^|q| to divide B (q > 0) or c (q < 0).
inline constexpr int kDecFastK = 160;
inline constexpr int kPow5Bits = 640;
struct Pow5Lower {
    BigInt m;     // 5^n ∈ [m·2^e, m·2^e·(1 + 2^-600)] for n < 2^24
    int e = 0;
};
// Binary powering, truncating every product to kPow5Bits bits (always down,
// so a lower bound). Squaring doubles the relative error of the base, so after
// j <= 24 squarings it is below 2^(j+2−kPow5Bits); the products add up to
// < 2^(28−kPow5Bits) < 2^-600 in total.
inline Pow5Lower pow5_lower(unsigned n) {
    Pow5Lower r;
    r.m = BigInt{1};
    BigInt b{5};
    int be = 0;
    auto trunc = [](BigInt& x, int& e) {
        const int bl = x.bit_length();
        if (bl > kPow5Bits) { x.shr(bl - kPow5Bits); e += bl - kPow5Bits; }
    };
    for (unsigned k = n;;) {
        if (k & 1) { r.m.mul(b); r.e += be; trunc(r.m, r.e); }
        k >>= 1;
        if (k == 0) break;
        const BigInt t = b;
        b.mul(t);
        be *= 2;
        trunc(b, be);
    }
    return r;
}
// X << sh for a signed shift (right shifts truncate).
inline BigInt shift_signed(BigInt x, long sh) {
    if (sh > 0) x.shl(static_cast<int>(sh));
    else if (sh < 0) x.shr(static_cast<int>(-sh));
    return x;
}
// floor(|v| / 10^k), v = sig·2^(wexp−255), |k| >= kDecFastK, with the sticky bit.
inline bool dec_floor_scaled_fast(const u256& sig, int wexp, int k, BigInt& q, bool& sticky) {
    constexpr int G = 128;  // guard bits below the integer part
    if (k < kDecFastK && k > -kDecFastK) return false;
    if (k >= (1 << 24) || k <= -(1 << 24)) return false;
    const BigInt S = BigInt::from_u256(sig);
    if (k > 0) {
        // v = S·2^t / 5^k, t = wexp − 255 − k; 5^k ∈ [m2^e, m2^e(1+δ)]
        // U = S·2^(t+G−e)/m >= v·2^G >= U(1 − δ); Z = floor(U)
        const Pow5Lower p = pow5_lower(static_cast<unsigned>(k));
        const long sh = static_cast<long>(wexp) - 255 - k + G - p.e;
        BigInt num = S, den = p.m;
        if (sh >= 0) num.shl(static_cast<int>(sh)); else den.shl(static_cast<int>(-sh));
        BigInt Z = num;
        (void)Z.divmod(den);
        if (Z.bit_length() > 600) return false;
        // v·2^G ∈ [Z − 1 − Z·2^-600, Z + 1): floor(v) is Z >> G if the low end agrees
        BigInt lo = Z;
        BigInt d = Z; d.shr(600); d.add_u32(1);
        if (lo < d) return false;
        lo.sub(d);
        BigInt zq = Z; zq.shr(G);
        lo.shr(G);
        if (!(lo == zq)) return false;
        q = std::move(zq);
        sticky = true;  // 5^k > 2^256 > sig: v is never an integer
        return true;
    }
    const int pk = -k;
    const long te2 = static_cast<long>(wexp) - 255 + pk;
    if (te2 >= 0) return false;  // v an integer with ~te2 bits more than the digits asked: exact path
    // v = S·5^pk / 2^sh: an integer iff 2^sh divides S (5^pk is odd)
    const long sh = -te2;
    sticky = sh > static_cast<long>(sig.ctz());
    if (!sticky) return false;  // exact multiple: the bound below would straddle it
    const Pow5Lower p = pow5_lower(static_cast<unsigned>(pk));
    BigInt X = S;
    X.mul(p.m);
    // v·2^G ∈ [X·2^(e−sh+G), X·2^(e−sh+G)·(1+δ)]; Z = floor of the low end
    const BigInt Z = shift_signed(X, p.e - sh + G);
    if (Z.bit_length() > 600) return false;
    BigInt hi = Z; BigInt d = Z; d.shr(600); hi.add(d); hi.add_u32(2);
    BigInt zq = Z; zq.shr(G);
    hi.shr(G);
    if (!(hi == zq)) return false;
    q = std::move(zq);
    return true;
}
// sign(c·10^q − B·2^E) for |q| >= kDecFastK, or 2 if the 640-bit bound does not decide it.
inline int cmp_dec_bin_fast(const BigInt& c, int q, const BigInt& B, int E, const Pow5Lower& p) {
    // q > 0: c·5^q·2^q vs B·2^E: L = c·m·2^(e+q) <= lhs <= L(1+δ).
    // q < 0: c·2^q vs B·5^|q|·2^E (both sides times 5^|q|): the bound is on the right.
    BigInt X = q > 0 ? c : B;
    X.mul(p.m);
    const long xe = static_cast<long>(p.e) + (q > 0 ? q : E);  // the bounded side: X·2^xe
    const BigInt& Y = q > 0 ? B : c;
    const long ye = q > 0 ? E : q;                               // the exact side: Y·2^ye
    // align to integers: X·2^(xe−mn), Y·2^(ye−mn)
    const long mn = xe < ye ? xe : ye;
    const BigInt Xs = shift_signed(X, xe - mn), Ys = shift_signed(Y, ye - mn);
    BigInt Xu = Xs; BigInt d = Xs; d.shr(600); Xu.add(d); Xu.add_u32(1);
    int s;
    if (Ys < Xs) s = 1;          // bounded side > exact side
    else if (Xu < Ys) s = -1;
    else return 2;
    return q > 0 ? s : -s;       // s compares the bounded side with the exact side
}

// floor(|v| / 10^k) for v = sig * 2^(wexp - 255), plus whether it was inexact,
// with exact BigInt arithmetic.
inline std::pair<BigInt, bool> dec_floor_scaled_exact(const u256& sig, int wexp, int k) {
    BigInt num = BigInt::from_u256(sig);
    bool sticky = false;
    if (k >= 0) {
        BigInt den{1};
        den.mul_pow5(static_cast<unsigned>(k));
        const int two_exp = wexp - 255 - k;
        if (two_exp >= 0) num.shl(two_exp);
        else den.shl(-two_exp);
        const BigInt rem = num.divmod(den);
        sticky = !rem.is_zero();
    } else {
        const int pk = -k;
        num.mul_pow5(static_cast<unsigned>(pk));
        const int te2 = wexp - 255 + pk;
        if (te2 >= 0) num.shl(te2);
        else sticky = num.shr_sticky(-te2);
    }
    return {std::move(num), sticky};
}
inline std::pair<BigInt, bool> dec_floor_scaled(const u256& sig, int wexp, int k) {
    BigInt q;
    bool st = false;
    if (dec_floor_scaled_fast(sig, wexp, k, q, st)) return {std::move(q), st};
    return dec_floor_scaled_exact(sig, wexp, k);
}

// Increment decision for decimal rounding given the first discarded digit.
inline bool dec_round_up(int last, bool sticky, bool odd, bool sign, round_mode mode) {
    const bool inexact = last != 0 || sticky;
    switch (mode) {
        case round_mode::to_nearest_even: return last > 5 || (last == 5 && (sticky || odd));
        case round_mode::to_nearest_away: return last >= 5;
        case round_mode::toward_zero:     return false;
        case round_mode::upward:          return inexact && !sign;
        case round_mode::downward:        return inexact && sign;
    }
    return false;
}

// Add one to a decimal digit string; returns true if it grew a digit.
inline bool dec_increment(std::string& d) {
    int i = static_cast<int>(d.size()) - 1;
    while (i >= 0 && d[static_cast<std::size_t>(i)] == '9') {
        d[static_cast<std::size_t>(i)] = '0';
        --i;
    }
    if (i < 0) {
        d.insert(d.begin(), '1');
        return true;
    }
    ++d[static_cast<std::size_t>(i)];
    return false;
}

// Compare c * 10^q with b * 2^e exactly (all integers nonnegative).
inline int cmp_dec_bin(const BigInt& c, int q, const BigInt& b, int e) {
    BigInt l = c, r = b;
    const int a = q < 0 ? -q : 0;  // multiply both sides by 10^a
    const int z = e < 0 ? -e : 0;  // and by 2^z
    const int ql = q + a;
    l.mul_pow5(static_cast<unsigned>(ql));
    l.shl(ql + z);
    r.mul_pow5(static_cast<unsigned>(a));
    r.shl(e + z + a);
    if (l < r) return -1;
    if (r < l) return 1;
    return 0;
}

} // namespace detail

namespace detail {
struct DecFloor {
    std::string digits;  // exactly the requested count, floor of |x|
    int exp10 = 0;       // first digit at 10^exp10
    bool sticky = false; // nonzero digits below
};
// Round a floor digit string to k <= digits.size() significant digits.
inline std::pair<std::string, int> round_floor_digits(const DecFloor& f, int k, bool sign,
                                                      round_mode mode, bool& inexact) {
    std::string d = f.digits.substr(0, static_cast<std::size_t>(k));
    int n = f.exp10;
    const int last = static_cast<std::size_t>(k) < f.digits.size() ? f.digits[static_cast<std::size_t>(k)] - '0' : 0;
    bool sticky = f.sticky;
    for (std::size_t i = static_cast<std::size_t>(k) + 1; !sticky && i < f.digits.size(); ++i)
        sticky = f.digits[i] != '0';
    inexact = last != 0 || sticky;
    const bool odd = ((d.back() - '0') & 1) != 0;
    if (dec_round_up(last, sticky, odd, sign, mode)) {
        if (dec_increment(d)) {
            d.pop_back();
            ++n;
        }
    }
    return {std::move(d), n};
}
} // namespace detail

inline detail::DecFloor float256::decimal_floor(int ndig) const {
    using namespace detail;
    DecFloor out;
    const Unpacked u = unpack();
    // Estimate the decimal exponent from the top 53 bits, then correct it.
    const std::uint64_t frac53 = (u.sig.d[3] >> 11) & 0x000FFFFFFFFFFFFFULL;
    const double mant = std::bit_cast<double>((1023ULL << 52) | frac53);
    int n = static_cast<int>(std::floor(static_cast<double>(u.wexp) * 0.30102999566398119521373889472449 +
                                        std::log10(mant)));
    for (int attempt = 0; attempt < 8; ++attempt) {
        auto [q, sticky] = dec_floor_scaled(u.sig, u.wexp, n - ndig + 1);
        std::string digs = q.to_dec();
        const int got = digs == "0" ? 0 : static_cast<int>(digs.size());
        if (got != ndig) {
            n += got - ndig;
            continue;
        }
        out.digits = std::move(digs);
        out.exp10 = n;
        out.sticky = sticky;
        return out;
    }
    out.digits.assign(static_cast<std::size_t>(ndig), '0');  // not reached
    out.exp10 = n;
    return out;
}

inline float256::DecimalDigits float256::decimal_significant(int ndig, round_mode mode) const {
    using namespace detail;
    if (ndig < 1) ndig = 1;
    const DecFloor f = decimal_floor(ndig + 1);
    bool inexact = false;
    auto [d, n] = round_floor_digits(f, ndig, signbit(), mode, inexact);
    if (inexact) raise_flag(fe_inexact);
    return DecimalDigits{std::move(d), n};
}

inline float256::DecimalDigits float256::decimal_shortest() const {
    using namespace detail;
    const Unpacked u = unpack();
    // Exact integer significand on the format's grid: |x| = M * 2^qe.
    const int qe = (u.wexp > kEmin ? u.wexp : kEmin) - 236;
    const int shift = (qe + 236) - u.wexp;  // >= 0
    const u256 M = (u.sig >> kLsbPos) >> shift;
    const bool even = !M.bit(0);
    const bool pow2_boundary = u.wexp > kEmin && M == (u256{1} << 236);
    // Rounding interval [lo, hi] (inclusive iff M even), as B * 2^E.
    BigInt hiB = BigInt::from_u256(M);
    hiB.shl(1);
    hiB.add_u32(1);
    const int hiE = qe - 1;
    BigInt loB = BigInt::from_u256(M);
    int loE;
    if (pow2_boundary) {
        loB.shl(2);
        loB.sub(BigInt{1});
        loE = qe - 2;
    } else {
        loB.shl(1);
        loB.sub(BigInt{1});
        loE = qe - 1;
    }
    // Candidates are roundings of the exact value to k significant digits:
    // c * 10^q with q = n - k + 1 >= q0 := n - (max_digits10 - 1). Scale the
    // bounds once by 10^-q0 so each check is a small-by-big multiply:
    //   c * 10^(q - q0) * Bd  vs  Bn,  where bound = Bn / Bd = B * 2^E / 10^q0.
    const unsigned saved = fenv::get_flags();
    const DecFloor fl = decimal_floor(max_digits10 + 1);
    auto cand = [&](int k) {
        bool ix = false;
        auto [d, n] = round_floor_digits(fl, k, false, round_mode::to_nearest_even, ix);
        return DecimalDigits{std::move(d), n};
    };
    DecimalDigits best = cand(max_digits10);
    const int q0 = best.exp10 - (max_digits10 - 1);
    auto scaled = [&](const BigInt& B, int E, BigInt& num, BigInt& den) {
        num = B;
        den = BigInt{1};
        if (E >= 0) num.shl(E); else den.shl(-E);
        if (q0 >= 0) { den.mul_pow5(static_cast<unsigned>(q0)); den.shl(q0); }
        else { num.mul_pow5(static_cast<unsigned>(-q0)); num.shl(-q0); }
    };
    // The exact bounds are built only if a check needs them (at large |q| the
    // 640-bit test in cmp_dec_bin_fast decides, see kDecFastK).
    BigInt loN, loD, hiN, hiD;
    bool have_exact = false;
    std::vector<std::pair<int, Pow5Lower>> p5;  // the reference returned is used before the next call
    auto pow5 = [&](int q) -> const Pow5Lower& {
        for (const auto& e : p5) if (e.first == q) return e.second;
        p5.emplace_back(q, pow5_lower(static_cast<unsigned>(q < 0 ? -q : q)));
        return p5.back().second;
    };
    auto ok = [&](const DecimalDigits& d) {
        const int q = d.exp10 - static_cast<int>(d.digits.size()) + 1;
        if (q < q0) return false;
        if ((q >= kDecFastK || q <= -kDecFastK) && q < (1 << 24) && q > -(1 << 24)) {
            const BigInt c = BigInt::from_dec(d.digits);
            const Pow5Lower& p = pow5(q);
            const int cl = cmp_dec_bin_fast(c, q, loB, loE, p);
            const int ch = cl == 2 ? 2 : cmp_dec_bin_fast(c, q, hiB, hiE, p);
            if (cl != 2 && ch != 2) {
                const bool above_lo = even ? cl >= 0 : cl > 0;
                const bool below_hi = even ? ch <= 0 : ch < 0;
                return above_lo && below_hi;
            }
        }
        if (!have_exact) {
            scaled(loB, loE, loN, loD);
            scaled(hiB, hiE, hiN, hiD);
            have_exact = true;
        }
        const BigInt c = BigInt::from_dec(d.digits + std::string(static_cast<std::size_t>(q - q0), '0'));
        BigInt l = c, h = c;
        l.mul(loD);
        h.mul(hiD);
        const bool above_lo = even ? !(l < loN) : (loN < l);
        const bool below_hi = even ? !(hiN < h) : (h < hiN);
        return above_lo && below_hi;
    };
    int lo = 1, hi = max_digits10;
    while (lo < hi) {
        const int mid = (lo + hi) / 2;
        DecimalDigits d = cand(mid);
        if (ok(d)) {
            best = std::move(d);
            hi = mid;
        } else {
            lo = mid + 1;
        }
    }
    if (static_cast<int>(best.digits.size()) != lo) best = cand(lo);
    // Trim trailing zeros (possible after a rounding carry).
    while (best.digits.size() > 1 && best.digits.back() == '0') best.digits.pop_back();
    fenv::set_flags(saved);  // shortest output is a representation, not a rounding
    return best;
}

inline std::string float256::decimal_fixed(int frac, round_mode mode) const {
    using namespace detail;
    const Unpacked u = unpack();
    auto [q, sticky] = dec_floor_scaled(u.sig, u.wexp, -frac - 1);
    std::string digs = q.to_dec();
    const int last = digs.back() - '0';
    digs.pop_back();
    if (digs.empty()) digs = "0";
    const bool odd = ((digs.back() - '0') & 1) != 0;
    if (last != 0 || sticky) raise_flag(fe_inexact);
    if (dec_round_up(last, sticky, odd, u.sign, mode)) dec_increment(digs);
    return digs;
}

namespace detail {
inline std::string exp_suffix(char e, int n, int min_digits) {
    std::string out(1, e);
    out.push_back(n < 0 ? '-' : '+');
    std::string num = std::to_string(n < 0 ? -n : n);
    while (static_cast<int>(num.size()) < min_digits) num.insert(num.begin(), '0');
    return out + num;
}
// Render digits d1d2... * 10^(exp10) (d1 at 10^exp10) in fixed notation.
inline std::string render_fixed(const std::string& digits, int exp10, int min_frac) {
    std::string out;
    const int nd = static_cast<int>(digits.size());
    if (exp10 >= 0) {
        if (exp10 + 1 >= nd) {
            out = digits + std::string(static_cast<std::size_t>(exp10 + 1 - nd), '0');
        } else {
            out = digits.substr(0, static_cast<std::size_t>(exp10 + 1)) + "." +
                  digits.substr(static_cast<std::size_t>(exp10 + 1));
        }
    } else {
        out = "0." + std::string(static_cast<std::size_t>(-exp10 - 1), '0') + digits;
    }
    if (min_frac > 0) {
        std::size_t dot = out.find('.');
        if (dot == std::string::npos) { out.push_back('.'); dot = out.size() - 1; }
        const int have = static_cast<int>(out.size() - dot - 1);
        if (have < min_frac) out.append(static_cast<std::size_t>(min_frac - have), '0');
    }
    return out;
}
inline std::string render_sci(const std::string& digits, int exp10, int min_frac) {
    std::string out(1, digits[0]);
    std::string rest = digits.substr(1);
    while (static_cast<int>(rest.size()) < min_frac) rest.push_back('0');
    if (!rest.empty()) out += "." + rest;
    return out + exp_suffix('e', exp10, 2);
}
} // namespace detail

inline std::string float256::to_string(int precision) const {
    using namespace detail;
    if (isnan()) return signbit() ? "-nan" : "nan";
    if (isinf()) return signbit() ? "-inf" : "inf";
    if (iszero()) return signbit() ? "-0" : "0";
    const DecimalDigits d = precision > 0 ? decimal_significant(precision, current_mode())
                                          : decimal_shortest();
    std::string out;
    if (signbit()) out.push_back('-');
    out.push_back(d.digits[0]);
    if (d.digits.size() > 1) {
        out.push_back('.');
        out.append(d.digits, 1, std::string::npos);
    }
    out.push_back('e');
    if (d.exp10 >= 0) out.push_back('+');
    out += std::to_string(d.exp10);
    return out;
}

inline std::string float256::format(std::chars_format fmt, int precision) const {
    using namespace detail;
    std::string sign = signbit() ? "-" : "";
    if (isnan()) return sign + "nan";
    if (isinf()) return sign + "inf";
    round_mode mode = current_mode();  // not const (see narrow_to)
    if (fmt == std::chars_format::hex) {
        if (iszero()) {
            std::string z = sign + "0";
            if (precision > 0) z += "." + std::string(static_cast<std::size_t>(precision), '0');
            return z + "p+0";
        }
        Unpacked u = unpack();
        u256 sig = u.sig;
        int e = u.wexp;
        int nhex = 59;  // 236 fraction bits
        if (precision >= 0 && precision < 59) {
            bool ix = false;
            if (round_bits(sig, 255 - 4 * precision, false, u.sign, mode, ix)) {
                sig = {};
                sig.set_bit(255);
                ++e;
            }
            if (ix) raise_flag(fe_inexact);
            nhex = precision;
        }
        std::string h;
        const u256 f = sig << 1;
        for (int i = 0; i < nhex; ++i) {
            int v = 0;
            for (int b2 = 0; b2 < 4; ++b2)
                if (f.bit(255 - (i * 4 + b2))) v |= 1 << (3 - b2);
            h.push_back("0123456789abcdef"[v]);
        }
        if (precision < 0)
            while (!h.empty() && h.back() == '0') h.pop_back();
        else
            while (static_cast<int>(h.size()) < precision) h.push_back('0');
        std::string out = sign + "1";
        if (!h.empty()) out += "." + h;
        return out + "p" + (e >= 0 ? "+" : "-") + std::to_string(e < 0 ? -e : e);
    }
    if (iszero()) {
        if (fmt == std::chars_format::scientific) {
            std::string z = sign + "0";
            if (precision > 0) z += "." + std::string(static_cast<std::size_t>(precision), '0');
            return z + "e+00";
        }
        if (fmt == std::chars_format::fixed && precision > 0)
            return sign + "0." + std::string(static_cast<std::size_t>(precision), '0');
        return sign + "0";
    }
    if (fmt == std::chars_format::fixed) {
        if (precision < 0) {
            const DecimalDigits d = decimal_shortest();
            return sign + render_fixed(d.digits, d.exp10, 0);
        }
        std::string digs = decimal_fixed(precision, mode);
        if (precision == 0) return sign + digs;
        while (static_cast<int>(digs.size()) <= precision) digs.insert(digs.begin(), '0');
        const std::size_t ip = digs.size() - static_cast<std::size_t>(precision);
        return sign + digs.substr(0, ip) + "." + digs.substr(ip);
    }
    if (fmt == std::chars_format::scientific) {
        const DecimalDigits d = precision < 0 ? decimal_shortest() : decimal_significant(precision + 1, mode);
        return sign + render_sci(d.digits, d.exp10, precision < 0 ? 0 : precision);
    }
    // general
    if (precision < 0) {
        const DecimalDigits d = decimal_shortest();
        const std::string e = render_sci(d.digits, d.exp10, 0);
        // length of render_fixed(d.digits, d.exp10, 0), without building it (78914
        // characters at the top of the range)
        const long nd = static_cast<long>(d.digits.size()), x = d.exp10;
        const long flen = x >= 0 ? (x + 1 >= nd ? x + 1 : nd + 1) : 2 + (-x - 1) + nd;
        if (flen <= static_cast<long>(e.size())) return sign + render_fixed(d.digits, d.exp10, 0);
        return sign + e;
    }
    const int P = precision == 0 ? 1 : precision;
    DecimalDigits d = decimal_significant(P, mode);
    std::string digs = d.digits;
    while (digs.size() > 1 && digs.back() == '0') digs.pop_back();  // %g trims zeros
    if (d.exp10 < -4 || d.exp10 >= P) return sign + render_sci(digs, d.exp10, 0);
    return sign + render_fixed(digs, d.exp10, 0);
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

[[nodiscard]] constexpr bool issignaling(float256 x) noexcept { return x.issignaling(); }

// fmin/fmax (C23 semantics, IEEE 754-2008 minNum/maxNum): a quiet NaN operand
// is treated as missing data and the other operand is returned; a signaling
// NaN operand yields a quiet NaN (its payload, quieted). Opposite-signed zeros:
// fmin gives -0, fmax gives +0.
[[nodiscard]] constexpr float256 fmin(float256 a, float256 b) noexcept {
    if (a.issignaling()) return a.quieted();
    if (b.issignaling()) return b.quieted();
    if (a.isnan()) return b;
    if (b.isnan()) return a;
    if (a.iszero() && b.iszero())
        return float256::zero(a.signbit() || b.signbit());
    return a < b ? a : b;
}
[[nodiscard]] constexpr float256 fmax(float256 a, float256 b) noexcept {
    if (a.issignaling()) return a.quieted();
    if (b.issignaling()) return b.quieted();
    if (a.isnan()) return b;
    if (b.isnan()) return a;
    if (a.iszero() && b.iszero())
        return float256::zero(a.signbit() && b.signbit());
    return a > b ? a : b;
}

// IEEE 754-2019 minimum/maximum: any NaN operand propagates (quieted);
// -0 < +0.
[[nodiscard]] constexpr float256 minimum(float256 a, float256 b) noexcept {
    if (a.isnan()) return a.quieted();
    if (b.isnan()) return b.quieted();
    if (a.iszero() && b.iszero()) return float256::zero(a.signbit() || b.signbit());
    return a < b ? a : b;
}
[[nodiscard]] constexpr float256 maximum(float256 a, float256 b) noexcept {
    if (a.isnan()) return a.quieted();
    if (b.isnan()) return b.quieted();
    if (a.iszero() && b.iszero()) return float256::zero(a.signbit() && b.signbit());
    return a > b ? a : b;
}
// IEEE 754-2019 minimumNumber/maximumNumber: a NaN operand (quiet or
// signaling) is ignored in favour of a number; NaN only if both are NaN.
[[nodiscard]] constexpr float256 minimum_number(float256 a, float256 b) noexcept {
    if (a.isnan() && b.isnan()) return a.quieted();
    if (a.isnan()) return b;
    if (b.isnan()) return a;
    return minimum(a, b);
}
[[nodiscard]] constexpr float256 maximum_number(float256 a, float256 b) noexcept {
    if (a.isnan() && b.isnan()) return a.quieted();
    if (a.isnan()) return b;
    if (b.isnan()) return a;
    return maximum(a, b);
}

[[nodiscard]] constexpr float256 fma(float256 a, float256 b, float256 c) noexcept {
    return float256::fma(a, b, c);
}
[[nodiscard]] constexpr float256 sqrt(float256 x) noexcept { return float256::sqrt(x); }

// Per-call rounding: the same operations rounded in `m` instead of the thread's
// current mode (README "Per-call rounding"). The result and the flags raised are
// those of the operator evaluated with fenv::set_round(m); the thread's mode is
// neither read nor changed, so these also work in constant expressions.
[[nodiscard]] constexpr float256 add(float256 a, float256 b, round_mode m) noexcept {
    return float256::add_sub_m(a, b, false, m);
}
[[nodiscard]] constexpr float256 sub(float256 a, float256 b, round_mode m) noexcept {
    return float256::add_sub_m(a, b, true, m);
}
[[nodiscard]] constexpr float256 mul(float256 a, float256 b, round_mode m) noexcept { return float256::mul_m(a, b, m); }
[[nodiscard]] constexpr float256 div(float256 a, float256 b, round_mode m) noexcept { return float256::div_m(a, b, m); }
[[nodiscard]] constexpr float256 fma(float256 a, float256 b, float256 c, round_mode m) noexcept {
    return float256::fma_m(a, b, c, m);
}
[[nodiscard]] constexpr float256 sqrt(float256 x, round_mode m) noexcept { return float256::sqrt_m(x, m); }

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
/// ldexp rounded in `m` (only a subnormal result rounds).
[[nodiscard]] constexpr float256 ldexp(float256 x, int n, round_mode m) noexcept {
    using namespace detail;
    Unpacked u = x.unpack();
    if (u.kind == Kind::nan || u.kind == Kind::inf || u.kind == Kind::zero) return x;
    const std::int64_t ne = static_cast<std::int64_t>(u.wexp) + n;
    if (ne > kEmax + 16) return float256::inf(u.sign);  // as ldexp(x, n)
    if (ne < static_cast<std::int64_t>(kEmin) - 256) return float256::zero(u.sign);
    u.wexp = static_cast<int>(ne);
    return float256::from_unpacked(u, m);
}

[[nodiscard]] constexpr float256 scalbn(float256 x, int n) noexcept { return ldexp(x, n); }
[[nodiscard]] constexpr float256 scalbn(float256 x, int n, round_mode m) noexcept { return ldexp(x, n, m); }

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

/// Round to an integral value in the current rounding mode; raises inexact
/// when the value changes (C `rint`, IEEE roundToIntegralExact).
[[nodiscard]] constexpr float256 rint(float256 x) noexcept {
    return float256::round_integral(x, detail::current_mode(), true);
}
/// As rint, but never raises inexact (C `nearbyint`).
[[nodiscard]] constexpr float256 nearbyint(float256 x) noexcept {
    return float256::round_integral(x, detail::current_mode(), false);
}
/// IEEE roundToIntegral with an explicit direction (no flags).
[[nodiscard]] constexpr float256 round_to_integral(float256 x, round_mode m) noexcept {
    return float256::round_integral(x, m, false);
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


/// IEEE 754 totalOrder: -NaN < -inf < ... < -0 < +0 < ... < +inf < +NaN,
/// NaNs ordered by payload (signaling before quiet for positive NaNs).
[[nodiscard]] constexpr bool total_order(float256 a, float256 b) noexcept {
    auto key = [](float256 x) {
        auto w = x.raw_bits();
        if (w[3] >> 63) {
            for (auto& v : w) v = ~v;
        } else {
            w[3] |= 1ULL << 63;
        }
        return w;
    };
    const auto ka = key(a), kb = key(b);
    for (int i = 3; i >= 0; --i) {
        if (ka[static_cast<std::size_t>(i)] != kb[static_cast<std::size_t>(i)])
            return ka[static_cast<std::size_t>(i)] < kb[static_cast<std::size_t>(i)];
    }
    return true;
}
/// totalOrder on absolute values.
[[nodiscard]] constexpr bool total_order_mag(float256 a, float256 b) noexcept {
    return total_order(abs(a), abs(b));
}

/// std::to_chars-style output. precision < 0 (the default) gives the
/// shortest round-trip representation.
inline std::to_chars_result to_chars(char* first, char* last, float256 v,
                                     std::chars_format fmt = std::chars_format::general,
                                     int precision = -1) {
    const std::string s = v.format(fmt, precision);
    if (static_cast<std::size_t>(last - first) < s.size()) return {last, std::errc::value_too_large};
    std::memcpy(first, s.data(), s.size());
    return {first + s.size(), std::errc{}};
}

/// std::from_chars-style parsing (no leading '+' or whitespace; hex input
/// has no "0x" prefix). The value is rounded in the current mode. Overflow
/// to infinity or underflow of a nonzero input to zero reports
/// result_out_of_range and leaves `value` unchanged.
inline std::from_chars_result from_chars(const char* first, const char* last, float256& value,
                                         std::chars_format fmt = std::chars_format::general) {
    const char* p = first;
    bool neg = false;
    if (p < last && *p == '-') { neg = true; ++p; }
    auto ieq = [&](const char* q, const char* word) {
        for (; *word; ++word, ++q) {
            if (q >= last || (static_cast<unsigned char>(*q) | 0x20) != static_cast<unsigned char>(*word)) return false;
        }
        return true;
    };
    if (ieq(p, "inf")) {
        const char* e = p + 3;
        if (ieq(p, "infinity")) e = p + 8;
        value = float256::inf(neg);
        return {e, std::errc{}};
    }
    if (ieq(p, "nan")) {
        const char* e = p + 3;
        if (e < last && *e == '(') {
            const char* q = e + 1;
            while (q < last && (std::isalnum(static_cast<unsigned char>(*q)) || *q == '_')) ++q;
            if (q < last && *q == ')') e = q + 1;
        }
        value = neg ? -float256::qnan() : float256::qnan();
        return {e, std::errc{}};
    }
    const bool hex = fmt == std::chars_format::hex;
    auto is_dig = [&](char c) {
        return hex ? std::isxdigit(static_cast<unsigned char>(c)) != 0 : (c >= '0' && c <= '9');
    };
    const char* q = p;
    bool any = false;
    while (q < last && is_dig(*q)) { ++q; any = true; }
    if (q < last && *q == '.') {
        ++q;
        while (q < last && is_dig(*q)) { ++q; any = true; }
    }
    if (!any) return {first, std::errc::invalid_argument};
    const char* mant_end = q;
    const char ech = hex ? 'p' : 'e';
    bool has_exp = false;
    if (q < last && (*q | 0x20) == ech && fmt != std::chars_format::fixed) {
        const char* r = q + 1;
        if (r < last && (*r == '+' || *r == '-')) ++r;
        const char* ds = r;
        while (r < last && *r >= '0' && *r <= '9') ++r;
        if (r > ds) { q = r; has_exp = true; }
    }
    if (fmt == std::chars_format::scientific && !has_exp) return {first, std::errc::invalid_argument};
    std::string buf;
    if (neg) buf.push_back('-');
    if (hex) buf += "0x";
    buf.append(p, static_cast<std::size_t>(q - p));
    if (hex && !has_exp) buf += "p0";
    (void)mant_end;
    const auto r = float256::parse(buf);
    if (!r) return {first, std::errc::invalid_argument};
    bool nonzero_input = false;
    for (const char* c = p; c < mant_end; ++c) {
        if (is_dig(*c) && *c != '0') { nonzero_input = true; break; }
    }
    if (r->isinf() || (r->iszero() && nonzero_input)) return {q, std::errc::result_out_of_range};
    value = *r;
    return {q, std::errc{}};
}


// ---- Per-call rounding for the remaining operations (README "Per-call rounding").
namespace detail {
/// Run f() with the thread's rounding mode set to m, then restore the caller's
/// mode (also if f throws). Flags raised by f() stay raised, as with set_round.
template <class Fn>
inline auto with_round(round_mode m, Fn&& f) -> decltype(f()) {
    const fenv::scoped_round g(m);
    return f();
}
} // namespace detail

inline std::string float256::to_string(int precision, round_mode m) const {
    return detail::with_round(m, [&] { return to_string(precision); });
}
inline std::string float256::format(std::chars_format fmt, int precision, round_mode m) const {
    return detail::with_round(m, [&] { return format(fmt, precision); });
}
inline std::optional<float256> float256::parse(std::string_view s, round_mode m) {
    return detail::with_round(m, [&] { return parse(s); });
}

/// rint / nearbyint rounded in `m` (rint raises inexact when the value changes).
[[nodiscard]] constexpr float256 rint(float256 x, round_mode m) noexcept { return float256::round_integral(x, m, true); }
[[nodiscard]] constexpr float256 nearbyint(float256 x, round_mode m) noexcept { return float256::round_integral(x, m, false); }
[[nodiscard]] inline float256 hypot(float256 x, float256 y, round_mode m) { return detail::with_round(m, [&] { return hypot(x, y); }); }
/// to_chars / from_chars with the decimal rounding in `m`.
inline std::to_chars_result to_chars(char* first, char* last, float256 v, std::chars_format fmt, int precision,
                                     round_mode m) {
    return detail::with_round(m, [&] { return to_chars(first, last, v, fmt, precision); });
}
inline std::from_chars_result from_chars(const char* first, const char* last, float256& value, std::chars_format fmt,
                                         round_mode m) {
    return detail::with_round(m, [&] { return from_chars(first, last, value, fmt); });
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
    // Not claimed (yet): the core operations, flags and all five rounding
    // modes are verified against MPFR, but decimal output is not
    // mode-aware, and totalOrder / convertToInteger are incomplete. See the
    // README "IEEE 754 conformance" table.
    static constexpr bool is_iec559 = false;
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
    // Default mode; the dynamic mode is base256::fenv::get_round().
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
// Numeric literals (3.14_f256, 1e-80_f256, 0x1p-3_f256, 42_f256) use the raw
// literal form: the compiler passes the exact source spelling, which is
// parsed with correct rounding. (A cooked `long double` overload would take
// precedence and round the literal to long double/double first.)
// Digit separators (1'000.5_f256) are accepted.
[[nodiscard]] inline float256 operator""_f256(const char* s) {
    std::string buf;
    for (const char* p = s; p && *p; ++p) {
        if (*p != '\'') buf.push_back(*p);
    }
    auto r = float256::parse(buf);
    if (!r) throw std::invalid_argument("float256 literal");
    return *r;
}
// String literals: "0.1"_f256, "-inf"_f256, "0x1.8p+1"_f256.
[[nodiscard]] inline float256 operator""_f256(const char* s, std::size_t n) {
    auto r = float256::parse(std::string_view(s, n));
    if (!r) throw std::invalid_argument("float256 literal");
    return *r;
}
} // namespace base256::literals

/// std::format support. Spec: [[fill]align][sign][#][0][width][.precision][type]
/// with type one of a A e E f F g G (none = shortest round-trip, general).
template <>
struct std::formatter<base256::float256, char> {
    char fill = ' ';
    char align = 0;      // '<', '>', '^' or 0
    char sign = '-';     // '+', '-', ' '
    bool alt = false;
    bool zero_pad = false;
    int width = 0;
    int precision = -1;
    char type = 0;

    constexpr auto parse(std::format_parse_context& ctx) {
        auto it = ctx.begin();
        const auto end = ctx.end();
        auto is_align = [](char c) { return c == '<' || c == '>' || c == '^'; };
        if (it != end && it + 1 != end && is_align(*(it + 1)) && *it != '}') {
            fill = *it;
            align = *(it + 1);
            it += 2;
        } else if (it != end && is_align(*it)) {
            align = *it++;
        }
        if (it != end && (*it == '+' || *it == '-' || *it == ' ')) sign = *it++;
        if (it != end && *it == '#') { alt = true; ++it; }
        if (it != end && *it == '0') { zero_pad = true; ++it; }
        while (it != end && *it >= '0' && *it <= '9') width = width * 10 + (*it++ - '0');
        if (it != end && *it == '.') {
            ++it;
            precision = 0;
            while (it != end && *it >= '0' && *it <= '9') precision = precision * 10 + (*it++ - '0');
        }
        if (it != end && *it != '}') {
            switch (*it) {
                case 'a': case 'A': case 'e': case 'E': case 'f': case 'F': case 'g': case 'G':
                    type = *it++;
                    break;
                default:
                    throw std::format_error("invalid float256 format type");
            }
        }
        if (it != end && *it != '}') throw std::format_error("invalid float256 format spec");
        return it;
    }

    template <class Ctx>
    auto format(const base256::float256& v, Ctx& ctx) const {
        std::chars_format fmt = std::chars_format::general;
        int prec = precision;
        switch (type) {
            case 'a': case 'A': fmt = std::chars_format::hex; break;
            case 'e': case 'E': fmt = std::chars_format::scientific; if (prec < 0) prec = 6; break;
            case 'f': case 'F': fmt = std::chars_format::fixed; if (prec < 0) prec = 6; break;
            case 'g': case 'G': fmt = std::chars_format::general; if (prec < 0) prec = 6; break;
            default: break;
        }
        std::string body = base256::float256(v).format(fmt, prec);
        bool neg = !body.empty() && body[0] == '-';
        if (neg) body.erase(0, 1);
        if (type == 'A' || type == 'E' || type == 'F' || type == 'G')
            for (char& c : body) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
        if (type == 'a' || type == 'A') body.insert(0, type == 'A' ? "0X" : "0x");
        if (alt && body.find('.') == std::string::npos && v.isfinite()) {
            const auto epos = body.find_first_of("eEpP");
            body.insert(epos == std::string::npos ? body.size() : epos, ".");
        }
        std::string sgn = neg ? "-" : (sign == '+' ? "+" : (sign == ' ' ? " " : ""));
        std::string out;
        const int len = static_cast<int>(sgn.size() + body.size());
        if (width > len && zero_pad && align == 0 && v.isfinite()) {
            out = sgn + std::string(static_cast<std::size_t>(width - len), '0') + body;
        } else if (width > len) {
            const int padn = width - len;
            const char a = align ? align : '>';
            const int left = a == '<' ? 0 : (a == '^' ? padn / 2 : padn);
            out = std::string(static_cast<std::size_t>(left), fill) + sgn + body +
                  std::string(static_cast<std::size_t>(padn - left), fill);
        } else {
            out = sgn + body;
        }
        return std::copy(out.begin(), out.end(), ctx.out());
    }
};
