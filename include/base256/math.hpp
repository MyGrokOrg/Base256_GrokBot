#pragma once

/// Elementary functions for binary256.
/// Target accuracy: a few ulps on the primary domain. Argument reduction for
/// trig uses a 256-bit 2π, so huge arguments lose the last bits of the
/// reduced angle (documented IEEE-754-style limitation without extra-precise π).

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

[[nodiscard]] inline float256 exp(float256 x) {
    const fenv::scoped_default_env env_guard_; // RNE inside, caller's flags kept
    if (x.isnan()) return x;
    if (x.isinf()) return x.signbit() ? float256::zero() : x;
    if (x.iszero()) return float256::one();
    if (x > float256(181800)) return float256::inf(false);
    if (x < float256(-182000)) return float256::zero();
    const detail_dd::Scaled e = detail_dd::exp_scaled(detail_dd::dd(x));
    return detail_dd::round_scaled(e.m, e.k);
}

[[nodiscard]] inline float256 exp2(float256 x) {
    const fenv::scoped_default_env env_guard_; // RNE inside, caller's flags kept
    if (x.isnan()) return x;
    if (x.isinf()) return x.signbit() ? float256::zero() : x;
    if (x.iszero()) return float256::one();
    if (x > float256(262200)) return float256::inf(false);
    if (x < float256(-262500)) return float256::zero();
    // 2^x = 2^k * exp(f ln2) with k = round(x) and f = x - k exact.
    const float256 kf = round(x);
    const float256 f = x - kf;
    const int k = static_cast<int>(kf.to_double());
    if (f.iszero()) return ldexp(float256::one(), k);
    const detail_dd::Scaled e = detail_dd::exp_scaled(detail_dd::ln2() * f);
    return detail_dd::round_scaled(e.m, e.k + k);
}

[[nodiscard]] inline float256 expm1(float256 x) {
    const fenv::scoped_default_env env_guard_; // RNE inside, caller's flags kept
    using detail_dd::dd;
    if (x.isnan()) return x;
    if (x.iszero()) return x;
    if (x.isinf()) return x.signbit() ? -float256::one() : x;
    if (x > float256(181800)) return float256::inf(false);
    if (x < float256(-200)) return -float256::one();  // e^x < 2^-288 << ulp(1)/2
    const float256 ax = abs(x);
    if (ax < ldexp(float256::one(), -200)) {
        // x + x^2/2; x^3/6 is below 2^-400 relative.
        return (dd(x) + dd(ldexp(x * x, -1))).round();
    }
    if (ax <= float256(0.34)) return detail_dd::expm1_reduced(dd(x)).round();
    const detail_dd::Scaled e = detail_dd::exp_scaled(dd(x));
    if (e.k > 600) return detail_dd::round_scaled(e.m, e.k);
    return (detail_dd::ldexp(e.m, e.k) - dd(float256::one())).round();
}

[[nodiscard]] inline float256 log(float256 x) {
    const fenv::scoped_default_env env_guard_; // RNE inside, caller's flags kept
    if (x.isnan()) return x;
    if (x.signbit() && !x.iszero()) return float256::qnan();
    if (x.iszero()) return float256::inf(true);
    if (x.isinf()) return x;
    return detail_dd::log(detail_dd::dd(x)).round();
}

[[nodiscard]] inline float256 log2(float256 x) {
    const fenv::scoped_default_env env_guard_; // RNE inside, caller's flags kept
    if (x.isnan()) return x;
    if (x.signbit() && !x.iszero()) return float256::qnan();
    if (x.iszero()) return float256::inf(true);
    if (x.isinf()) return x;
    // frexp returns m in [0.5, 1) with x = m * 2^e, so a pure power of two
    // is m == 0.5 and log2(x) = e - 1 exactly.
    int e = 0;
    const float256 m = frexp(x, &e);
    if (m == float256::half())
        return float256::from_int(static_cast<std::int64_t>(e) - 1);
    return (detail_dd::log(detail_dd::dd(x)) / detail_dd::ln2()).round();
}
[[nodiscard]] inline float256 log10(float256 x) {
    const fenv::scoped_default_env env_guard_; // RNE inside, caller's flags kept
    if (x.isnan()) return x;
    if (x.signbit() && !x.iszero()) return float256::qnan();
    if (x.iszero()) return float256::inf(true);
    if (x.isinf()) return x;
    return (detail_dd::log(detail_dd::dd(x)) / detail_dd::ln10()).round();
}

[[nodiscard]] inline float256 log1p(float256 x) {
    const fenv::scoped_default_env env_guard_; // RNE inside, caller's flags kept
    using detail_dd::dd;
    if (x.isnan()) return x;
    if (x.iszero()) return x; // log1p(±0) = ±0
    if (x == float256::from_int(-1)) return float256::inf(true);
    if (x < float256::from_int(-1)) return float256::qnan();
    if (x.isinf()) return x; // +∞ only; −∞ was rejected above
    if (abs(x) < ldexp(float256::one(), -120)) {
        // x - x^2/2 + x^3/3; the next term is below 2^-360 relative.
        const dd x2 = detail_dd::two_prod(x, x);
        return (dd(x) - detail_dd::ldexp(x2, -1) + dd(x2.hi * x / float256::from_int(3))).round();
    }
    return detail_dd::log(detail_dd::two_sum(float256::one(), x)).round();
}

[[nodiscard]] inline float256 pow(float256 x, float256 y) {
    const fenv::scoped_default_env env_guard_; // RNE inside, caller's flags kept
    // C11 F.10.4.4 / IEEE 754-2019 recommended specials. These fire before
    // the NaN rule: pow(x, ±0) and pow(+1, y) are +1 for every input, NaN
    // included. pow(−1, ±∞) is +1 as well.
    if (y.iszero()) return float256::one();
    if (x == float256::one()) return float256::one();
    if (x.isnan() || y.isnan()) return float256::qnan();
    if (y.isinf()) {
        if (x == -float256::one()) return float256::one();
        const float256 ax = abs(x);
        if (ax > float256::one())
            return y.signbit() ? float256::zero() : float256::inf(false);
        if (ax < float256::one())
            return y.signbit() ? float256::inf(false) : float256::zero();
        return float256::qnan();
    }
    if (x.iszero()) {
        const bool odd = detail_math::is_odd_integer(y);
        if (y.signbit()) {
            // pow(±0, −odd) = ±∞; otherwise +∞.
            return odd ? float256::inf(x.signbit()) : float256::inf(false);
        }
        return odd ? float256::zero(x.signbit()) : float256::zero(false);
    }
    if (x.isinf()) {
        const bool yneg = y.signbit();
        if (!x.signbit())
            return yneg ? float256::zero() : float256::inf(false);
        // −∞. Non-integers lose the sign; odd integers keep it.
        if (trunc(y) != y)
            return yneg ? float256::zero() : float256::inf(false);
        const bool odd = detail_math::is_odd_integer(y);
        if (yneg) return odd ? float256::zero(true) : float256::zero(false);
        return odd ? float256::inf(true) : float256::inf(false);
    }
    // Integer exponents in this range are a short product, and exact whenever
    // the mathematical power is a binary256 number. exp(y log x) is not.
    if (x.signbit() && trunc(y) != y) return float256::qnan();
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
            if (n < 0) r = detail_dd::dd(float256::one()) / r;
            const float256 p = r.round();
            return neg ? -p : p;
        }
    }
    const detail_dd::dd L = detail_dd::log(detail_dd::dd(ax)) * y;
    float256 p;
    if (L.hi > float256(181800)) p = float256::inf(false);
    else if (L.hi < float256(-182000)) p = float256::zero();
    else {
        const detail_dd::Scaled e = detail_dd::exp_scaled(L);
        p = detail_dd::round_scaled(e.m, e.k);
    }
    return neg ? -p : p;
}

[[nodiscard]] inline float256 cbrt(float256 x) {
    const fenv::scoped_default_env env_guard_; // RNE inside, caller's flags kept
    if (x.isnan() || x.isinf() || x.iszero()) return x;
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
    const float256 res = ldexp((detail_dd::dd(y) + detail_dd::dd(corr)).round(), j);
    return s ? -res : res;
}

[[nodiscard]] inline float256 two_pi() {
    return float256::pi() * float256::two();
}

[[nodiscard]] inline float256 reduce_trig(float256 x, int& quadrant) {
    // Reduce to [0, π/2) and report quadrant 0..3. Uses 256-bit π.
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
// sin/cos/tan of |x| through the dd path. Returns false for |x| >= 2^100
// (caller falls back to the float256 reduction).
inline bool trig_dd(float256 ax, detail_dd::dd& s, detail_dd::dd& c, int& quad) {
    if (ilogb(ax) >= detail_dd::kTrigExactExp) return false;
    const detail_dd::dd r = detail_dd::reduce_pio2(ax, quad);
    detail_dd::sincos_reduced(r, s, c);
    return true;
}
} // namespace detail_math

[[nodiscard]] inline float256 sin(float256 x) {
    const fenv::scoped_default_env env_guard_; // RNE inside, caller's flags kept
    if (x.isnan() || x.iszero()) return x;
    if (x.isinf()) return float256::qnan();
    int quad = 0;
    const bool neg = x.signbit();
    detail_dd::dd sd, cd;
    if (detail_math::trig_dd(abs(x), sd, cd, quad)) {
        detail_dd::dd y = (quad & 1) ? cd : sd;
        if (quad >= 2) y = -y;
        const float256 out = y.round();
        return neg ? -out : out;
    }
    float256 r = reduce_trig(x, quad);
    if (r.iszero()) {
        switch (quad) {
            case 0: return copysign(float256::zero(), x);
            case 1: return neg ? -float256::one() : float256::one();
            case 2: return copysign(float256::zero(), x);
            default: return neg ? float256::one() : -float256::one();
        }
    }
    float256 s, c;
    s = detail_math::taylor_sin(r);
    c = detail_math::taylor_cos(r);
    float256 y;
    switch (quad) {
        case 0: y = s; break;
        case 1: y = c; break;
        case 2: y = -s; break;
        default: y = -c; break;
    }
    return neg ? -y : y;
}

[[nodiscard]] inline float256 cos(float256 x) {
    const fenv::scoped_default_env env_guard_; // RNE inside, caller's flags kept
    if (x.isnan()) return x;
    if (x.isinf()) return float256::qnan();
    if (x.iszero()) return float256::one();
    int quad = 0;
    detail_dd::dd sd, cd;
    if (detail_math::trig_dd(abs(x), sd, cd, quad)) {
        switch (quad) {
            case 0: return cd.round();
            case 1: return (-sd).round();
            case 2: return (-cd).round();
            default: return sd.round();
        }
    }
    float256 r = reduce_trig(x, quad);
    if (r.iszero()) {
        switch (quad) {
            case 0: return float256::one();
            case 1: return float256::zero();
            case 2: return -float256::one();
            default: return float256::zero();
        }
    }
    const float256 s = detail_math::taylor_sin(r);
    const float256 c = detail_math::taylor_cos(r);
    switch (quad) {
        case 0: return c;
        case 1: return -s;
        case 2: return -c;
        default: return s;
    }
}

[[nodiscard]] inline float256 tan(float256 x) {
    const fenv::scoped_default_env env_guard_; // RNE inside, caller's flags kept
    if (x.isnan()) return x;
    if (x.isinf()) return float256::qnan();
    if (x.iszero()) return x;
    int quad = 0;
    detail_dd::dd sd, cd;
    if (detail_math::trig_dd(abs(x), sd, cd, quad)) {
        // tan = s/c in even quadrants, -c/s in odd ones.
        const float256 t = ((quad & 1) ? -(cd / sd) : sd / cd).round();
        return x.signbit() ? -t : t;
    }
    const float256 r = reduce_trig(x, quad);
    if (r.iszero()) {
        if ((quad & 1) == 0) return copysign(float256::zero(), x);
        return float256::inf(x.signbit());
    }
    return sin(x) / cos(x);
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
    const fenv::scoped_default_env env_guard_; // RNE inside, caller's flags kept
    if (x.isnan()) return x;
    if (x.isinf()) return copysign(float256::pi() * float256::half(), x);
    if (x.iszero()) return x;
    if (abs(x) < ldexp(float256::one(), -120)) {
        // x - x^3/3; next term is below 2^-480 relative.
        return (detail_dd::dd(x) - detail_dd::dd(x * x * x / float256::from_int(3))).round();
    }
    return detail_math::atan2_dd(detail_dd::dd(x), detail_dd::dd(float256::one())).round();
}

[[nodiscard]] inline float256 atan2(float256 y, float256 x) {
    const fenv::scoped_default_env env_guard_; // RNE inside, caller's flags kept
    if (x.isnan() || y.isnan()) return float256::qnan();
    if (y.isinf() && x.isinf()) {
        // atan2(±∞, +∞) = ±π/4, atan2(±∞, −∞) = ±3π/4
        const float256 pi4 = float256::pi() * float256::half() * float256::half();
        const float256 mag = x.signbit() ? pi4 * float256::from_int(3) : pi4;
        return copysign(mag, y);
    }
    if (x.isinf()) {
        if (x.signbit()) return copysign(float256::pi(), y);
        return copysign(float256::zero(), y);
    }
    if (y.isinf()) return copysign(float256::pi() * float256::half(), y);
    if (x.iszero() && y.iszero()) {
        if (x.signbit()) return copysign(float256::pi(), y);
        return y; // ±0
    }
    if (x.signbit() && y.iszero()) return copysign(float256::pi(), y);
    if (y.iszero()) return y;
    return detail_math::atan2_dd(detail_dd::dd(y), detail_dd::dd(x)).round();
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
    const fenv::scoped_default_env env_guard_; // RNE inside, caller's flags kept
    if (x.isnan() || x.iszero()) return x;
    const float256 a = abs(x);
    if (a > float256::one()) return float256::qnan();
    if (a < ldexp(float256::one(), -120)) {
        // x + x^3/6
        return (detail_dd::dd(x) + detail_dd::dd(x * x * x / float256::from_int(6))).round();
    }
    return detail_math::atan2_dd(detail_dd::dd(x), detail_math::sqrt_one_minus_sq(x)).round();
}

[[nodiscard]] inline float256 acos(float256 x) {
    const fenv::scoped_default_env env_guard_; // RNE inside, caller's flags kept
    if (x.isnan()) return x;
    const float256 a = abs(x);
    if (a > float256::one()) return float256::qnan();
    if (x == float256::one()) return float256::zero();
    return detail_math::atan2_dd(detail_math::sqrt_one_minus_sq(x), detail_dd::dd(x)).round();
}

[[nodiscard]] inline float256 sinh(float256 x) {
    const fenv::scoped_default_env env_guard_; // RNE inside, caller's flags kept
    using detail_dd::dd;
    if (x.isnan() || x.isinf() || x.iszero()) return x;
    const float256 a = abs(x);
    float256 y;
    if (a > float256(181800)) {
        y = float256::inf(false);
    } else if (a < ldexp(float256::one(), -120)) {
        y = (dd(a) + dd(a * a * a / float256::from_int(6))).round();
    } else if (a <= float256(0.34)) {
        // (em + em/(em+1))/2 with em = expm1(a): no cancellation.
        const dd em = detail_dd::expm1_reduced(dd(a));
        y = detail_dd::ldexp(em + em / (em + dd(float256::one())), -1).round();
    } else {
        const detail_dd::Scaled e = detail_dd::exp_scaled(dd(a));
        if (e.k > 300) {
            y = detail_dd::round_scaled(e.m, e.k - 1);
        } else {
            const dd E = detail_dd::ldexp(e.m, e.k);
            y = detail_dd::ldexp(E - dd(float256::one()) / E, -1).round();
        }
    }
    return x.signbit() ? -y : y;
}
[[nodiscard]] inline float256 cosh(float256 x) {
    const fenv::scoped_default_env env_guard_; // RNE inside, caller's flags kept
    using detail_dd::dd;
    if (x.isnan()) return x;
    if (x.isinf()) return float256::inf(false);
    const float256 a = abs(x);
    if (a > float256(181800)) return float256::inf(false);
    if (a < ldexp(float256::one(), -120)) return (dd(float256::one()) + dd(ldexp(a * a, -1))).round();
    if (a <= float256(0.34)) {
        // 1 + em^2 / (2 (1 + em))
        const dd em = detail_dd::expm1_reduced(dd(a));
        return (dd(float256::one()) + detail_dd::ldexp(em * em / (em + dd(float256::one())), -1)).round();
    }
    const detail_dd::Scaled e = detail_dd::exp_scaled(dd(a));
    if (e.k > 300) return detail_dd::round_scaled(e.m, e.k - 1);
    const dd E = detail_dd::ldexp(e.m, e.k);
    return detail_dd::ldexp(E + dd(float256::one()) / E, -1).round();
}
[[nodiscard]] inline float256 tanh(float256 x) {
    const fenv::scoped_default_env env_guard_; // RNE inside, caller's flags kept
    using detail_dd::dd;
    if (x.isnan()) return x;
    if (x.isinf()) return copysign(float256::one(), x);
    if (x.iszero()) return x;
    const float256 a = abs(x);
    // 1 - tanh(x) ≈ 2 exp(-2|x|) is below ulp(1)/4 once |x| >= 83.
    if (a >= float256::from_int(83)) return copysign(float256::one(), x);
    float256 y;
    if (a < ldexp(float256::one(), -120)) {
        y = (dd(a) - dd(a * a * a / float256::from_int(3))).round();
    } else {
        const dd em = detail_dd::expm1(dd(a + a));
        y = (em / (em + dd(float256::two()))).round();
    }
    return x.signbit() ? -y : y;
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
