#pragma once

/// Elementary functions for binary256.
/// Target accuracy: a few ulps on the primary domain. Argument reduction for
/// trig uses a 256-bit 2π, so huge arguments lose the last bits of the
/// reduced angle (documented IEEE-754-style limitation without extra-precise π).

#include "float256.hpp"

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

    // exp(x) = 2^{x/ln2} = 2^{n+f} = 2^n * exp(f ln2), |f|<=0.5
    const float256 y = x / float256::ln2();
    const float256 nf = round(y);
    // n fits in int for any finite binary256 exp input that doesn't overflow
    int n = 0;
    if (nf.isfinite() && !nf.isnan()) {
        const double nd = nf.to_double();
        if (nd > 300000.0) return float256::inf(false);
        if (nd < -300000.0) return float256::zero(false);
        n = static_cast<int>(nd);
    }
    const float256 r = x - nf * float256::ln2();
    return ldexp(detail_math::taylor_exp(r), n);
}

[[nodiscard]] inline float256 exp2(float256 x) {
    const fenv::scoped_default_env env_guard_; // RNE inside, caller's flags kept
    if (x.isnan()) return x;
    if (x.isinf()) return x.signbit() ? float256::zero() : x;
    if (x.iszero()) return float256::one();
    // 2^n is exact for every integer n in the binary256 exponent range.
    if (trunc(x) == x) {
        const double nd = x.to_double();
        if (!(nd <= 300000.0)) return float256::inf(false);
        if (nd < -300000.0) return float256::zero(false);
        return ldexp(float256::one(), static_cast<int>(nd));
    }
    return exp(x * float256::ln2());
}

[[nodiscard]] inline float256 expm1(float256 x) {
    const fenv::scoped_default_env env_guard_; // RNE inside, caller's flags kept
    if (x.isnan()) return x;
    if (x.iszero()) return x;
    // For small x the Taylor of exp(x)-1 is more accurate
    const float256 ax = abs(x);
    if (ax < float256::half()) {
        float256 term = x;
        float256 sum  = x;
        for (int i = 2; i <= 80; ++i) {
            term = term * x / float256::from_int(i);
            const float256 next = sum + term;
            if (next == sum) break;
            sum = next;
        }
        return sum;
    }
    return exp(x) - float256::one();
}

[[nodiscard]] inline float256 log(float256 x) {
    const fenv::scoped_default_env env_guard_; // RNE inside, caller's flags kept
    if (x.isnan()) return x;
    if (x.signbit() && !x.iszero()) return float256::qnan();
    if (x.iszero()) return float256::inf(true);
    if (x.isinf()) return x;

    int e = 0;
    const float256 m = frexp(x, &e); // m in [0.5, 1)
    // Write x = 2^{e-1} * (1+f) with 1+f in [1, 2)
    const float256 one = float256::one();
    float256 y = m;
    int adj = e;
    if (y < one) {
        y = y + y; // [1, 2)
        --adj;
    }
    const float256 f = y - one;                  // [0, 1)
    const float256 z = f / (y + one);            // f/(2+f) = atanh argument
    const float256 z2 = z * z;
    float256 term = z;
    float256 sum  = z;
    for (int n = 1; n <= 120; ++n) {
        term = term * z2;
        const float256 add = term / float256::from_int(static_cast<std::int64_t>(2 * n + 1));
        const float256 next = sum + add;
        if (next == sum) break;
        sum = next;
    }
    const float256 log_m = sum + sum; // 2 * atanh(z) = log(1+f)
    return float256::from_int(static_cast<std::int64_t>(adj)) * float256::ln2() + log_m;
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
    return log(x) * float256::log2e();
}
[[nodiscard]] inline float256 log10(float256 x) {
    const fenv::scoped_default_env env_guard_; // RNE inside, caller's flags kept
    return log(x) * float256::log10e();
}

[[nodiscard]] inline float256 log1p(float256 x) {
    const fenv::scoped_default_env env_guard_; // RNE inside, caller's flags kept
    if (x.isnan()) return x;
    if (x.iszero()) return x; // log1p(±0) = ±0
    if (x == float256::from_int(-1)) return float256::inf(true);
    if (x < float256::from_int(-1)) return float256::qnan();
    if (x.isinf()) return x; // +∞ only; −∞ was rejected above
    // Mercator series. |x| < 1/16 converges in well under 80 terms at
    // 237-bit precision. Below 2^{-237}, 1+x rounds to 1, so log(1+x)
    // would flush the whole result to zero.
    const float256 ax = abs(x);
    if (ax < ldexp(float256::one(), -4) || float256::one() + x == float256::one()) {
        float256 s = x;
        float256 sum = x;
        for (int n = 2; n <= 80; ++n) {
            s = s * x;
            float256 term = s / float256::from_int(n);
            if ((n & 1) == 0) term = -term;
            const float256 next = sum + term;
            if (next == sum) break;
            sum = next;
        }
        return sum;
    }
    return log(float256::one() + x);
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
    if (trunc(y) == y && abs(y) <= float256(4096)) {
        const std::int64_t n = static_cast<std::int64_t>(y.to_double());
        if (n < 0) {
            const float256 mag = detail_math::integer_pow(x, static_cast<std::uint64_t>(-n));
            return float256::one() / mag;
        }
        return detail_math::integer_pow(x, static_cast<std::uint64_t>(n));
    }
    if (x.signbit()) {
        if (trunc(y) != y) return float256::qnan();
        const float256 p = exp(y * log(abs(x)));
        if (detail_math::is_odd_integer(y)) return -p;
        return p;
    }
    return exp(y * log(x));
}

[[nodiscard]] inline float256 cbrt(float256 x) {
    const fenv::scoped_default_env env_guard_; // RNE inside, caller's flags kept
    if (x.isnan() || x.isinf() || x.iszero()) return x;
    const bool s = x.signbit();
    const float256 ax = abs(x);
    // initial guess from double
    float256 y = float256(std::cbrt(ax.to_double()));
    if (y.iszero() || y.isinf()) {
        // scale
        int e = 0;
        const float256 m = frexp(ax, &e);
        y = float256(std::cbrt(m.to_double()));
        y = ldexp(y, e / 3);
    }
    // Halley: y <- y * (y^3 + 2x) / (2 y^3 + x)
    for (int i = 0; i < 8; ++i) {
        const float256 y3 = y * y * y;
        const float256 num = y3 + ax * float256::two();
        const float256 den = y3 * float256::two() + ax;
        const float256 nxt = y * (num / den);
        if (nxt == y) break;
        y = nxt;
    }
    return s ? -y : y;
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

[[nodiscard]] inline float256 sin(float256 x) {
    const fenv::scoped_default_env env_guard_; // RNE inside, caller's flags kept
    if (x.isnan() || x.iszero()) return x;
    if (x.isinf()) return float256::qnan();
    int quad = 0;
    const bool neg = x.signbit();
    float256 r = reduce_trig(x, quad);
    // Exact multiples of π/2 reduce to 0. Taylor then yields a spurious −0
    // in the odd quadrants (sin(π), cos(π/2)).
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
    const float256 r = reduce_trig(x, quad);
    if (r.iszero()) {
        // Exact multiple of π/2. Even quadrants are kπ (signed zero);
        // odd quadrants are the poles.
        if ((quad & 1) == 0) return copysign(float256::zero(), x);
        return float256::inf(x.signbit());
    }
    return sin(x) / cos(x);
}

[[nodiscard]] inline float256 atan(float256 x) {
    const fenv::scoped_default_env env_guard_; // RNE inside, caller's flags kept
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
    const float256 a = atan(y / x);
    if (x.signbit()) {
        return y.signbit() ? a - float256::pi() : a + float256::pi();
    }
    return a;
}

[[nodiscard]] inline float256 asin(float256 x) {
    const fenv::scoped_default_env env_guard_; // RNE inside, caller's flags kept
    if (x.isnan()) return x;
    const float256 a = abs(x);
    if (a > float256::one()) return float256::qnan();
    if (a == float256::one()) return copysign(float256::pi() * float256::half(), x);
    // asin(x) = atan(x / sqrt(1-x^2))
    return atan(x / sqrt(float256::one() - x * x));
}

[[nodiscard]] inline float256 acos(float256 x) {
    const fenv::scoped_default_env env_guard_; // RNE inside, caller's flags kept
    if (x.isnan()) return x;
    const float256 a = abs(x);
    if (a > float256::one()) return float256::qnan();
    return float256::pi() * float256::half() - asin(x);
}

[[nodiscard]] inline float256 sinh(float256 x) {
    const fenv::scoped_default_env env_guard_; // RNE inside, caller's flags kept
    if (x.isnan() || x.isinf() || x.iszero()) return x;
    const float256 a = abs(x);
    float256 y;
    if (a < float256::one()) {
        // (e^a - e^{-a})/2 = (em + em/(em+1))/2, exact in expm1.
        // exp(a) - exp(-a) cancels the whole result when |a| << 1.
        const float256 em = expm1(a);
        y = (em + em / (em + float256::one())) * float256::half();
    } else {
        const float256 e = exp(a);
        y = (e - float256::one() / e) * float256::half();
    }
    return x.signbit() ? -y : y;
}
[[nodiscard]] inline float256 cosh(float256 x) {
    const fenv::scoped_default_env env_guard_; // RNE inside, caller's flags kept
    if (x.isnan()) return x;
    if (x.isinf()) return float256::inf(false);
    const float256 a = abs(x);
    if (a < float256::one()) {
        // cosh(x) = 1 + 2 sinh(x/2)^2 keeps the departure from 1.
        const float256 s = sinh(a * float256::half());
        return float256::one() + (s + s) * s;
    }
    const float256 e = exp(a);
    return (e + float256::one() / e) * float256::half();
}
[[nodiscard]] inline float256 tanh(float256 x) {
    const fenv::scoped_default_env env_guard_; // RNE inside, caller's flags kept
    if (x.isnan()) return x;
    if (x.isinf()) return copysign(float256::one(), x);
    if (x.iszero()) return x;
    const float256 a = abs(x);
    // 1 - tanh(x) ≈ 2 exp(-2|x|) drops below ½ ulp(1) once |x| >= 83,
    // so the result is exactly ±1. Below that, exp(2x) is still finite;
    // past it, exp(2x) overflows and (inf-1)/(inf+1) would be NaN.
    if (a >= float256::from_int(83)) return copysign(float256::one(), x);
    const float256 e = expm1(x + x);
    return e / (e + float256::two());
}

[[nodiscard]] inline float256 fmod(float256 x, float256 y) {
    const fenv::scoped_default_env env_guard_; // RNE inside, caller's flags kept
    if (x.isnan() || y.isnan()) return float256::qnan();
    if (y.iszero() || x.isinf()) return float256::qnan();
    if (y.isinf()) return x;
    const float256 q = trunc(x / y);
    return x - q * y;
}

[[nodiscard]] inline float256 remainder(float256 x, float256 y) {
    const fenv::scoped_default_env env_guard_; // RNE inside, caller's flags kept
    if (x.isnan() || y.isnan()) return float256::qnan();
    if (y.iszero() || x.isinf()) return float256::qnan();
    if (y.isinf()) return x;
    const float256 q = rint(x / y); // remainder: ties to even
    return x - q * y;
}

} // namespace base256
