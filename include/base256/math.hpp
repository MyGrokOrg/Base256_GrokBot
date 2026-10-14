#pragma once

/// Elementary functions for binary256.
/// Accuracy: < 1 ulp (≈ ½ ulp observed against MPFR) over the full range;
/// every function evaluates in double-float256 (detail/dd.hpp) and rounds
/// once. Trig reduction is Cody–Waite below 2^100 and Payne–Hanek (full
/// 2/π table, detail/two_over_pi.hpp) above, accurate for every finite x.

#include "float256.hpp"
#include "detail/dd.hpp"

namespace base256 {
namespace detail_math {

[[nodiscard]] inline bool is_odd_integer(float256 y) noexcept {
    if (!y.isfinite()) return false;
    const float256 yi = trunc(y);
    if (yi != y) return false;
    // |y| >= 2^237 is a multiple of a large power of two, hence even.
    const float256 half_i = yi * float256::half();
    return trunc(half_i) != half_i;
}

// base^n for n >= 0 by squaring. Exact when every product is a dyadic rational
// inside the 237-bit significand (powers of two, small factorials, 10^k, …).
[[nodiscard]] inline float256 integer_pow(float256 base, std::uint64_t n) {
    float256 result = float256::one();
    while (n > 0) {
        if ((n & 1ULL) != 0) result = result * base;
        n >>= 1ULL;
        if (n > 0) base = base * base;
    }
    return result;
}

[[nodiscard]] inline float256 taylor_exp(float256 r) {
    // exp(r) for |r| <= ln2/2 ≈ 0.347
    float256 term = float256::one();
    float256 sum  = float256::one();
    for (int i = 1; i <= 80; ++i) {
        term = term * r / float256::from_int(i);
        const float256 next = sum + term;
        if (next == sum) break;
        sum = next;
    }
    return sum;
}

[[nodiscard]] inline float256 taylor_sin(float256 r) {
    // sin(r) = r - r^3/3! + r^5/5! - ...   |r| <= π/4
    const float256 r2 = r * r;
    float256 term = r;
    float256 sum  = r;
    for (int n = 1; n <= 60; ++n) {
        term = term * r2 / float256::from_int(static_cast<std::int64_t>((2 * n) * (2 * n + 1)));
        term = -term;
        const float256 next = sum + term;
        if (next == sum) break;
        sum = next;
    }
    return sum;
}

[[nodiscard]] inline float256 taylor_cos(float256 r) {
    const float256 r2 = r * r;
    float256 term = float256::one();
    float256 sum  = float256::one();
    for (int n = 1; n <= 60; ++n) {
        term = term * r2 / float256::from_int(static_cast<std::int64_t>((2 * n - 1) * (2 * n)));
        term = -term;
        const float256 next = sum + term;
        if (next == sum) break;
        sum = next;
    }
    return sum;
}

} // namespace detail_math

namespace detail_math {
[[nodiscard]] inline int sgn(float256 v) noexcept { return v.signbit() ? -1 : 1; }
inline constexpr int kTiny = -120;  // below 2^-120 the leading Taylor terms suffice
[[nodiscard]] inline bool is_tiny(float256 a) { return ilogb(a) < kTiny; }

// x·(1 + c) with the correction c = x^2·coef (coef > 0 or < 0) for tiny x, as
// an exact-direction dd: when x^2·coef underflows to zero the tail still
// records its sign. `next_sign` is the sign of the first omitted term
// relative to x (+1 same sign as x, -1 opposite).
[[nodiscard]] inline float256 tiny_odd(detail_dd::FnEnv& env, float256 x, int coef_num, int coef_den, int next_sign) {
    using detail_dd::dd;
    const float256 c = x * x * x * float256::from_int(coef_num) / float256::from_int(coef_den);
    if (c.iszero()) return env.round(x, float256::zero(), 0, sgn(x) * (coef_num > 0 ? 1 : -1));
    const dd v = dd(x) + dd(c);
    return env.round(v.hi, v.lo, 0, sgn(x) * next_sign);
}
} // namespace detail_math

[[nodiscard]] inline float256 exp(float256 x) {
    detail_dd::FnEnv env;
    using detail_dd::dd;
    if (x.isnan()) return env.nan(x);
    if (x.isinf()) return env.exact(x.signbit() ? float256::zero() : x);
    if (x.iszero()) return env.exact(float256::one());
    if (x > float256(181800)) return env.overflow(false);
    if (x < float256(-182000)) return env.underflow(false);
    if (detail_math::is_tiny(abs(x))) {
        // 1 + x + x^2/2; the omitted x^3/6 is far below the rounding position
        const dd v = detail_dd::two_sum(float256::one(), x) + dd(ldexp(x * x, -1));
        return env.round(v.hi, v.lo, 0, +1);
    }
    const detail_dd::Scaled e = detail_dd::exp_scaled(dd(x));
    return env.round(e.m.hi, e.m.lo, e.k);
}

[[nodiscard]] inline float256 exp2(float256 x) {
    detail_dd::FnEnv env;
    using detail_dd::dd;
    if (x.isnan()) return env.nan(x);
    if (x.isinf()) return env.exact(x.signbit() ? float256::zero() : x);
    if (x.iszero()) return env.exact(float256::one());
    if (x > float256(262200)) return env.overflow(false);
    if (x < float256(-262500)) return env.underflow(false);
    // 2^x = 2^k * exp(f ln2) with k = round(x) and f = x - k exact.
    const float256 kf = round(x);
    const float256 f = x - kf;
    const int k = static_cast<int>(kf.to_double());
    if (f.iszero()) return env.round(float256::one(), float256::zero(), k, 0);  // exact unless out of range
    const dd t = detail_dd::ln2() * f;
    if (detail_math::is_tiny(abs(t.hi))) {
        const dd v = dd(float256::one()) + t + dd(ldexp(t.hi * t.hi, -1));
        return env.round(v.hi, v.lo, k, +1);
    }
    const detail_dd::Scaled e = detail_dd::exp_scaled(t);
    return env.round(e.m.hi, e.m.lo, e.k + k);
}

[[nodiscard]] inline float256 expm1(float256 x) {
    detail_dd::FnEnv env;
    using detail_dd::dd;
    if (x.isnan()) return env.nan(x);
    if (x.iszero()) return env.exact(x);
    if (x.isinf()) return env.exact(x.signbit() ? -float256::one() : x);
    if (x > float256(181800)) return env.overflow(false);
    if (x < float256(-200)) return env.round(-float256::one(), float256::zero(), 0, +1);  // -1 + e^x
    const float256 ax = abs(x);
    if (ilogb(ax) < -200) {
        // x + x^2/2 (+ x^3/6 of the sign of x)
        const float256 c = ldexp(x * x, -1);
        if (c.iszero()) return env.round(x, float256::zero(), 0, +1);
        const dd v = dd(x) + dd(c);
        return env.round(v.hi, v.lo, 0, detail_math::sgn(x));
    }
    if (ax <= float256(0.34)) return detail_dd::expm1_reduced(dd(x)).out(env);
    const detail_dd::Scaled e = detail_dd::exp_scaled(dd(x));
    if (e.k > 600) return env.round(e.m.hi, e.m.lo, e.k);
    return (detail_dd::ldexp(e.m, e.k) - dd(float256::one())).out(env);
}

[[nodiscard]] inline float256 log(float256 x) {
    detail_dd::FnEnv env;
    if (x.isnan()) return env.nan(x);
    if (x.iszero()) return env.pole(true);
    if (x.signbit()) return env.invalid();
    if (x.isinf()) return env.exact(x);
    if (x == float256::one()) return env.exact(float256::zero());
    return detail_dd::log(detail_dd::dd(x)).out(env);
}

[[nodiscard]] inline float256 log2(float256 x) {
    detail_dd::FnEnv env;
    if (x.isnan()) return env.nan(x);
    if (x.iszero()) return env.pole(true);
    if (x.signbit()) return env.invalid();
    if (x.isinf()) return env.exact(x);
    // frexp returns m in [0.5, 1) with x = m * 2^e, so a pure power of two
    // is m == 0.5 and log2(x) = e - 1 exactly.
    int e = 0;
    const float256 m = frexp(x, &e);
    if (m == float256::half()) return env.exact(float256::from_int(static_cast<std::int64_t>(e) - 1));
    return (detail_dd::log(detail_dd::dd(x)) / detail_dd::ln2()).out(env);
}

[[nodiscard]] inline float256 log10(float256 x) {
    detail_dd::FnEnv env;
    if (x.isnan()) return env.nan(x);
    if (x.iszero()) return env.pole(true);
    if (x.signbit()) return env.invalid();
    if (x.isinf()) return env.exact(x);
    const detail_dd::dd r = detail_dd::log(detail_dd::dd(x)) / detail_dd::ln10();
    // log10 is exact only at x = 10^n, n = 0..102 (10^n has an odd part 5^n < 2^237).
    const float256 n = round(r.hi);
    if (n >= float256::zero() && n <= float256(102) && abs(r.hi - n) < ldexp(float256::one(), -100)) {
        detail_dd::dd p{float256::one()}, b{float256(10)};
        for (unsigned u = static_cast<unsigned>(n.to_double()); u; u >>= 1) {
            if (u & 1U) p = p * b;
            if (u > 1) b = b * b;
        }
        if (p.lo.iszero() && p.hi == x) return env.exact(n);
    }
    return r.out(env);
}

[[nodiscard]] inline float256 log1p(float256 x) {
    detail_dd::FnEnv env;
    using detail_dd::dd;
    if (x.isnan()) return env.nan(x);
    if (x.iszero()) return env.exact(x); // log1p(±0) = ±0
    if (x == -float256::one()) return env.pole(true);
    if (x < -float256::one()) return env.invalid();
    if (x.isinf()) return env.exact(x); // +∞ only; −∞ was rejected above
    if (detail_math::is_tiny(abs(x))) {
        // x - x^2/2 + x^3/3; the omitted -x^4/4 is negative.
        const dd x2 = detail_dd::two_prod(x, x);
        if (x2.hi.iszero()) return env.round(x, float256::zero(), 0, -1);
        const dd v = dd(x) - detail_dd::ldexp(x2, -1) + dd(x2.hi * x / float256::from_int(3));
        return env.round(v.hi, v.lo, 0, -1);
    }
    return detail_dd::log(detail_dd::two_sum(float256::one(), x)).out(env);
}

namespace detail_math {
// Exact product a·b (true when the rounding error is zero).
[[nodiscard]] inline bool exact_mul(float256 a, float256 b, float256& out) {
    const detail_dd::dd p = detail_dd::two_prod(a, b);
    out = p.hi;
    return p.lo.iszero() && p.hi.isfinite() && !(p.hi.iszero() && !a.iszero() && !b.iszero());
}
// base^n exactly (n >= 1), or false.
[[nodiscard]] inline bool exact_ipow(float256 base, std::uint64_t n, float256& out) {
    float256 r = float256::one(), b = base;
    while (n) {
        if (n & 1ULL) { if (!exact_mul(r, b, r)) return false; }
        n >>= 1ULL;
        if (n) { if (!exact_mul(b, b, b)) return false; }
    }
    out = r;
    return true;
}
// Is c = x^y exactly? For y = m / 2^j (j <= 16, |m| <= 2^20): c^(2^j) == x^m.
[[nodiscard]] inline bool pow_is_exact(float256 x, float256 y, float256 c) {
    if (!c.isfinite() || c.iszero()) return false;
    int j = 0;
    float256 m = y;
    while (trunc(m) != m && j < 16) { m = m + m; ++j; }
    if (trunc(m) != m || abs(m) > float256(1 << 20)) return false;
    const std::int64_t mi = static_cast<std::int64_t>(m.to_double());
    float256 lhs;  // c^(2^j)
    if (!exact_ipow(c, 1ULL << j, lhs)) return false;
    float256 rhs;
    if (!exact_ipow(x, static_cast<std::uint64_t>(mi < 0 ? -mi : mi), rhs)) return false;
    if (mi < 0) {  // c^(2^j) · x^|m| == 1
        float256 prod;
        return exact_mul(lhs, rhs, prod) && prod == float256::one();
    }
    return lhs == rhs;
}
} // namespace detail_math

[[nodiscard]] inline float256 pow(float256 x, float256 y) {
    detail_dd::FnEnv env;
    // C11 F.10.4.4 / IEEE 754-2019: pow(x, ±0) and pow(+1, y) are +1 for every
    // quiet input, NaN included; pow(−1, ±∞) is +1. A signaling NaN is invalid.
    if (x.issignaling() || y.issignaling()) return env.nan2(x, y);
    if (y.iszero()) return env.exact(float256::one());
    if (x == float256::one()) return env.exact(float256::one());
    if (x.isnan() || y.isnan()) return env.nan2(x, y);
    if (y.isinf()) {
        if (x == -float256::one()) return env.exact(float256::one());
        const float256 ax = abs(x);
        if (ax > float256::one()) return env.exact(y.signbit() ? float256::zero() : float256::inf(false));
        return env.exact(y.signbit() ? float256::inf(false) : float256::zero());
    }
    if (x.iszero()) {
        const bool odd = detail_math::is_odd_integer(y);
        // pow(±0, −odd) = ±∞, pow(±0, −other) = +∞: divide-by-zero.
        if (y.signbit()) return env.pole(odd && x.signbit());
        return env.exact(odd ? float256::zero(x.signbit()) : float256::zero(false));
    }
    if (x.isinf()) {
        const bool yneg = y.signbit();
        if (!x.signbit()) return env.exact(yneg ? float256::zero() : float256::inf(false));
        // −∞. Non-integers lose the sign; odd integers keep it.
        const bool odd = detail_math::is_odd_integer(y);
        if (yneg) return env.exact(float256::zero(odd));
        return env.exact(float256::inf(odd));
    }
    if (x.signbit() && trunc(y) != y) return env.invalid();
    const bool neg = x.signbit() && detail_math::is_odd_integer(y);
    const float256 ax = abs(x);
    // Integer exponents with |y| <= 4096, when the result stays well inside
    // the exponent range: dd repeated squaring. Exact whenever the true
    // power is representable (every partial power then fits in 474 bits).
    if (trunc(y) == y && abs(y) <= float256(4096)) {
        const std::int64_t n = static_cast<std::int64_t>(y.to_double());
        const double mag = static_cast<double>(ilogb(ax)) * static_cast<double>(n < 0 ? -n : n);
        if (mag < 261000.0 && mag > -261000.0) {
            detail_dd::dd r{float256::one()}, b{ax};
            for (std::uint64_t u = static_cast<std::uint64_t>(n < 0 ? -n : n); u; u >>= 1) {
                if (u & 1ULL) r = r * b;
                if (u > 1) b = b * b;
            }
            bool exact = r.lo.iszero();
            if (n < 0) {
                int e = 0;
                exact = exact && frexp(r.hi, &e) == float256::half();
                r = detail_dd::dd(float256::one()) / r;
            }
            if (neg) r = -r;
            return env.round(r.hi, r.lo, 0, exact ? 0 : 2);
        }
    }
    const detail_dd::dd L = detail_dd::log(detail_dd::dd(ax)) * y;
    if (L.hi > float256(181800)) return env.overflow(neg);
    if (L.hi < float256(-182000)) return env.underflow(neg);
    detail_dd::Scaled e = detail_dd::exp_scaled(L);
    if (neg) e.m = -e.m;
    // Exact non-integer cases such as pow(4, 0.5) or pow(8, 1/2^k·m).
    {
        const float256 c = abs(detail_dd::round_scaled(e.m, e.k));
        if (abs(e.m.lo) <= ldexp(abs(e.m.hi), -400) && detail_math::pow_is_exact(ax, y, c))
            return env.exact(neg ? -c : c);
    }
    return env.round(e.m.hi, e.m.lo, e.k);
}

[[nodiscard]] inline float256 cbrt(float256 x) {
    detail_dd::FnEnv env;
    if (x.isnan()) return env.nan(x);
    if (x.isinf() || x.iszero()) return env.exact(x);
    const bool s = x.signbit();
    // Scale by 2^(3j) to [1, 8): the cube root scales by 2^j exactly.
    const int j = ilogb(x) >= 0 ? ilogb(x) / 3 : -((-ilogb(x) + 2) / 3);
    const float256 ax = ldexp(abs(x), -3 * j);
    float256 y = float256(std::cbrt(ax.to_double()));
    // Halley: y <- y * (y^3 + 2x) / (2 y^3 + x)
    for (int i = 0; i < 6; ++i) {
        const float256 y3 = y * y * y;
        const float256 nxt = y * ((y3 + ax * float256::two()) / (y3 * float256::two() + ax));
        if (nxt == y) break;
        y = nxt;
    }
    // One Newton correction from the exact dd residual x - y^3.
    const detail_dd::dd y3 = detail_dd::two_prod(y, y) * y;
    const detail_dd::dd r = detail_dd::dd(ax) - y3;
    const float256 corr = r.hi / (float256::from_int(3) * y * y);
    detail_dd::dd res = detail_dd::dd(y) + detail_dd::dd(corr);
    // Perfect cubes are exact.
    {
        const float256 c = res.round();
        float256 c2, c3;
        if (detail_math::exact_mul(c, c, c2) && detail_math::exact_mul(c2, c, c3) && c3 == ax)
            return env.exact(ldexp(s ? -c : c, j));
    }
    if (s) res = -res;
    return env.round(res.hi, res.lo, j);
}

[[nodiscard]] inline float256 two_pi() {
    return float256::pi() * float256::two();
}

/// Legacy float256 reduction to [0, π/2) with quadrant 0..3 (256-bit π; loses
/// accuracy for large |x|). Kept for compatibility; sin/cos/tan no longer use it.
[[nodiscard]] inline float256 reduce_trig(float256 x, int& quadrant) {
    const float256 tp = two_pi();
    const float256 half_pi = float256::pi() * float256::half();
    float256 a = abs(x);
    // n = floor(a / (π/2))
    const float256 qf = trunc(a / half_pi);
    double qd = qf.to_double();
    if (!std::isfinite(qd) || std::fabs(qd) > 1.0e16) {
        // last-ditch: modulo via remainder of huge division
        const float256 n = trunc(a / tp);
        a = a - n * tp;
        if (a < float256::zero()) a = a + tp;
        qd = trunc(a / half_pi).to_double();
    }
    long qn = static_cast<long>(qd);
    a = a - qf * half_pi;
    if (a.signbit()) a = -a;
    if (a >= half_pi) {
        a = a - half_pi;
        ++qn;
    }
    quadrant = static_cast<int>(qn & 3);
    if (x.signbit()) {
        // odd function for sin; handled by caller
    }
    return a;
}

namespace detail_math {
// sin/cos/tan of a finite |x| through the dd path: Cody–Waite below 2^100,
// Payne–Hanek above.
inline void trig_dd(float256 ax, detail_dd::dd& s, detail_dd::dd& c, int& quad) {
    const detail_dd::dd r = ilogb(ax) >= detail_dd::kTrigExactExp ? detail_dd::payne_hanek(ax, quad)
                                                                   : detail_dd::reduce_pio2(ax, quad);
    detail_dd::sincos_reduced(r, s, c);
}
} // namespace detail_math

[[nodiscard]] inline float256 sin(float256 x) {
    detail_dd::FnEnv env;
    if (x.isnan()) return env.nan(x);
    if (x.iszero()) return env.exact(x);
    if (x.isinf()) return env.invalid();
    if (detail_math::is_tiny(abs(x))) return detail_math::tiny_odd(env, x, -1, 6, +1);  // x - x^3/6 (+x^5/120)
    int quad = 0;
    detail_dd::dd sd, cd;
    detail_math::trig_dd(abs(x), sd, cd, quad);
    detail_dd::dd y = (quad & 1) ? cd : sd;
    if (quad >= 2) y = -y;
    if (x.signbit()) y = -y;
    return y.out(env);
}

[[nodiscard]] inline float256 cos(float256 x) {
    detail_dd::FnEnv env;
    if (x.isnan()) return env.nan(x);
    if (x.isinf()) return env.invalid();
    if (x.iszero()) return env.exact(float256::one());
    if (detail_math::is_tiny(abs(x))) {
        // 1 - x^2/2 (+ x^4/24)
        const detail_dd::dd v = detail_dd::dd(float256::one()) - detail_dd::dd(ldexp(x * x, -1));
        return env.round(v.hi, v.lo, 0, v.lo.iszero() ? -1 : +1);
    }
    int quad = 0;
    detail_dd::dd sd, cd;
    detail_math::trig_dd(abs(x), sd, cd, quad);
    switch (quad) {
        case 0: return cd.out(env);
        case 1: return (-sd).out(env);
        case 2: return (-cd).out(env);
        default: return sd.out(env);
    }
}

[[nodiscard]] inline float256 tan(float256 x) {
    detail_dd::FnEnv env;
    if (x.isnan()) return env.nan(x);
    if (x.isinf()) return env.invalid();
    if (x.iszero()) return env.exact(x);
    if (detail_math::is_tiny(abs(x))) return detail_math::tiny_odd(env, x, 1, 3, +1);  // x + x^3/3 (+2x^5/15)
    int quad = 0;
    detail_dd::dd sd, cd;
    detail_math::trig_dd(abs(x), sd, cd, quad);
    // tan = s/c in even quadrants, -c/s in odd ones.
    detail_dd::dd t = (quad & 1) ? -(cd / sd) : sd / cd;
    if (x.signbit()) t = -t;
    return t.out(env);
}

namespace detail_math {
// float256 arctangent, a few ulps; seed for the dd correction below.
[[nodiscard]] inline float256 atan_seed(float256 x) {
    if (x.isnan() || x.isinf()) {
        if (x.isinf()) return copysign(float256::pi() * float256::half(), x);
        return x;
    }
    if (x.iszero()) return x;
    bool sign = x.signbit();
    float256 a = abs(x);
    bool recip = false;
    if (a > float256::one()) {
        a = float256::one() / a;
        recip = true;
    }
    // atan(a) for a in [0,1]: series in terms of z, or Taylor of atan
    // atan(a) = a - a^3/3 + a^5/5 - ...  slow near 1, use identity
    // atan(a) = atan(c) + atan((a-c)/(1+a c)) with c = 2-sqrt(3) ≈ tan(π/12)
    const float256 tan15 = float256::parse("0.267949192431122706472553658494127633").value_or(float256(0.2679491924311227));
    int k = 0;
    while (a > tan15 && k < 8) {
        // atan(a) = π/6 + atan((a*√3 - 1)/(√3 + a))  -- reduce toward 0
        const float256 s3 = float256::sqrt2() * float256::parse("1.224744871391589049098642037352945696").value_or(float256(1.224744871391589));
        (void)s3;
        // simpler: half-angle: atan(a) = 2 atan(a/(1+sqrt(1+a^2)))
        a = a / (float256::one() + sqrt(float256::one() + a * a));
        ++k;
    }
    const float256 a2 = a * a;
    float256 term = a;
    float256 sum  = a;
    for (int n = 1; n <= 80; ++n) {
        term = term * a2;
        float256 add = term / float256::from_int(static_cast<std::int64_t>(2 * n + 1));
        if (n & 1) add = -add;
        const float256 next = sum + add;
        if (next == sum) break;
        sum = next;
    }
    for (int i = 0; i < k; ++i) sum = sum + sum;
    if (recip) sum = float256::pi() * float256::half() - sum;
    return sign ? -sum : sum;
}


// atan2 seed for finite (y, x) not both zero.
[[nodiscard]] inline float256 atan2_seed(float256 y, float256 x) {
    if (x.iszero()) return copysign(float256::pi() * float256::half(), y);
    const float256 a = atan_seed(y / x);
    if (x.signbit()) return y.signbit() ? a - float256::pi() : a + float256::pi();
    return a;
}

// sin and cos of a float256 |a| <= 4 as dd.
inline void sincos_small(float256 a, detail_dd::dd& s, detail_dd::dd& c) {
    int q = 0;
    detail_dd::dd sd, cd;
    const detail_dd::dd r = detail_dd::reduce_pio2(abs(a), q);
    detail_dd::sincos_reduced(r, sd, cd);
    s = (q & 1) ? cd : sd;
    c = (q & 1) ? -sd : cd;
    if (q >= 2) { s = -s; c = -c; }
    if (a.signbit()) s = -s;
}

// atan2 of dd (Y, X), finite, not both zero: θ0 from the seed, then
// θ = θ0 + atan(t), t = (Y cos θ0 − X sin θ0)/(X cos θ0 + Y sin θ0), |t| ~ 2^-230.
[[nodiscard]] inline detail_dd::dd atan2_dd(detail_dd::dd Y, detail_dd::dd X) {
    using detail_dd::dd;
    // Bring both to a moderate exponent so the products neither overflow nor go subnormal.
    const int ey = Y.hi.iszero() ? -1000000 : ilogb(Y.hi);
    const int ex = X.hi.iszero() ? -1000000 : ilogb(X.hi);
    const int e = ey > ex ? ey : ex;
    Y = detail_dd::ldexp(Y, -e);
    X = detail_dd::ldexp(X, -e);
    const float256 t0 = atan2_seed(Y.round(), X.round());
    dd s, c;
    sincos_small(t0, s, c);
    const dd num = Y * c - X * s;
    const dd den = X * c + Y * s;
    const dd t = num / den;
    return dd(t0) + (t - t * t * t / dd(float256::from_int(3)));
}
} // namespace detail_math

[[nodiscard]] inline float256 atan(float256 x) {
    detail_dd::FnEnv env;
    if (x.isnan()) return env.nan(x);
    if (x.isinf()) {
        const detail_dd::dd h = x.signbit() ? -detail_dd::pio2() : detail_dd::pio2();
        return h.out(env);
    }
    if (x.iszero()) return env.exact(x);
    if (detail_math::is_tiny(abs(x))) return detail_math::tiny_odd(env, x, -1, 3, +1);  // x - x^3/3 (+x^5/5)
    return detail_math::atan2_dd(detail_dd::dd(x), detail_dd::dd(float256::one())).out(env);
}

[[nodiscard]] inline float256 atan2(float256 y, float256 x) {
    detail_dd::FnEnv env;
    using detail_dd::dd;
    if (x.isnan() || y.isnan()) return env.nan2(y, x);
    const bool ys = y.signbit();
    auto signed_out = [&](dd v) { return (ys ? -v : v).out(env); };
    if (y.isinf() && x.isinf()) {
        // atan2(±∞, +∞) = ±π/4, atan2(±∞, −∞) = ±3π/4
        const dd pi4 = detail_dd::ldexp(detail_dd::pi(), -2);
        return signed_out(x.signbit() ? pi4 * float256::from_int(3) : pi4);
    }
    if (x.isinf()) {
        if (x.signbit()) return signed_out(detail_dd::pi());
        return env.exact(float256::zero(ys));
    }
    if (y.isinf()) return signed_out(detail_dd::pio2());
    if (y.iszero()) {
        if (x.signbit()) return signed_out(detail_dd::pi());
        return env.exact(y);  // ±0
    }
    if (x.iszero()) return signed_out(detail_dd::pio2());
    // Tiny ratio with x > 0: atan(q) = q - q^3/3, q = y/x as dd.
    if (!x.signbit() && ilogb(y) - ilogb(x) < detail_math::kTiny - 2) {
        const dd q = dd(y) / dd(x);
        const float256 c = q.hi * q.hi * q.hi / float256::from_int(-3);
        if (c.iszero()) {
            // q itself may be exact; the omitted -q^3/3 points toward zero.
            return env.round(q.hi, q.lo, 0, q.hi.signbit() ? +1 : -1);
        }
        const dd v = q + dd(c);
        return env.round(v.hi, v.lo, 0, q.hi.signbit() ? -1 : +1);
    }
    return detail_math::atan2_dd(dd(y), dd(x)).out(env);
}

namespace detail_math {
// sqrt((1 - x)(1 + x)) in dd; 1 ± x are exact as dd.
[[nodiscard]] inline detail_dd::dd sqrt_one_minus_sq(float256 x) {
    const detail_dd::dd a = detail_dd::two_sum(float256::one(), -x);
    const detail_dd::dd b = detail_dd::two_sum(float256::one(), x);
    return detail_dd::sqrt(a * b);
}
} // namespace detail_math

[[nodiscard]] inline float256 asin(float256 x) {
    detail_dd::FnEnv env;
    if (x.isnan()) return env.nan(x);
    if (x.iszero()) return env.exact(x);
    const float256 a = abs(x);
    if (a > float256::one()) return env.invalid();
    if (detail_math::is_tiny(a)) return detail_math::tiny_odd(env, x, 1, 6, +1);  // x + x^3/6 (+3x^5/40)
    return detail_math::atan2_dd(detail_dd::dd(x), detail_math::sqrt_one_minus_sq(x)).out(env);
}

[[nodiscard]] inline float256 acos(float256 x) {
    detail_dd::FnEnv env;
    if (x.isnan()) return env.nan(x);
    const float256 a = abs(x);
    if (a > float256::one()) return env.invalid();
    if (x == float256::one()) return env.exact(float256::zero());
    return detail_math::atan2_dd(detail_math::sqrt_one_minus_sq(x), detail_dd::dd(x)).out(env);
}

[[nodiscard]] inline float256 sinh(float256 x) {
    detail_dd::FnEnv env;
    using detail_dd::dd;
    if (x.isnan()) return env.nan(x);
    if (x.isinf() || x.iszero()) return env.exact(x);
    const float256 a = abs(x);
    if (a > float256(181800)) return env.overflow(x.signbit());
    if (detail_math::is_tiny(a)) return detail_math::tiny_odd(env, x, 1, 6, +1);  // x + x^3/6 (+x^5/120)
    dd y;
    int k = 0;
    if (a <= float256(0.34)) {
        // (em + em/(em+1))/2 with em = expm1(a): no cancellation.
        const dd em = detail_dd::expm1_reduced(dd(a));
        y = detail_dd::ldexp(em + em / (em + dd(float256::one())), -1);
    } else {
        const detail_dd::Scaled e = detail_dd::exp_scaled(dd(a));
        if (e.k > 300) {
            y = e.m;
            k = e.k - 1;
        } else {
            const dd E = detail_dd::ldexp(e.m, e.k);
            y = detail_dd::ldexp(E - dd(float256::one()) / E, -1);
        }
    }
    if (x.signbit()) y = -y;
    return y.out(env, k);
}

[[nodiscard]] inline float256 cosh(float256 x) {
    detail_dd::FnEnv env;
    using detail_dd::dd;
    if (x.isnan()) return env.nan(x);
    if (x.isinf()) return env.exact(float256::inf(false));
    if (x.iszero()) return env.exact(float256::one());
    const float256 a = abs(x);
    if (a > float256(181800)) return env.overflow(false);
    if (detail_math::is_tiny(a)) {
        const dd v = dd(float256::one()) + dd(ldexp(a * a, -1));  // 1 + a^2/2 (+a^4/24)
        return env.round(v.hi, v.lo, 0, +1);
    }
    if (a <= float256(0.34)) {
        // 1 + em^2 / (2 (1 + em))
        const dd em = detail_dd::expm1_reduced(dd(a));
        return (dd(float256::one()) + detail_dd::ldexp(em * em / (em + dd(float256::one())), -1)).out(env);
    }
    const detail_dd::Scaled e = detail_dd::exp_scaled(dd(a));
    if (e.k > 300) return e.m.out(env, e.k - 1);
    const dd E = detail_dd::ldexp(e.m, e.k);
    return detail_dd::ldexp(E + dd(float256::one()) / E, -1).out(env);
}

[[nodiscard]] inline float256 tanh(float256 x) {
    detail_dd::FnEnv env;
    using detail_dd::dd;
    if (x.isnan()) return env.nan(x);
    if (x.isinf()) return env.exact(copysign(float256::one(), x));
    if (x.iszero()) return env.exact(x);
    const float256 a = abs(x);
    // 1 - tanh(x) ≈ 2 exp(-2|x|) is below ulp(1)/4 once |x| >= 83.
    if (a >= float256::from_int(83))
        return env.round(copysign(float256::one(), x), float256::zero(), 0, x.signbit() ? +1 : -1);
    if (detail_math::is_tiny(a)) return detail_math::tiny_odd(env, x, -1, 3, +1);  // x - x^3/3 (+2x^5/15)
    const dd em = detail_dd::expm1(dd(a + a));
    dd y = em / (em + dd(float256::two()));
    if (x.signbit()) y = -y;
    return y.out(env);
}

/// Exact: x - trunc(x/y)*y, sign of x. invalid for y == 0 or x infinite.
[[nodiscard]] constexpr float256 fmod(float256 x, float256 y) noexcept {
    return float256::remainder_impl(x, y, false, nullptr);
}

/// Exact IEEE remainder: x - n*y with n = x/y rounded to nearest, ties to even.
[[nodiscard]] constexpr float256 remainder(float256 x, float256 y) noexcept {
    return float256::remainder_impl(x, y, true, nullptr);
}

} // namespace base256
