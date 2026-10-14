#pragma once

/// Extra elementary functions layered on math.hpp.
/// Gamma, joint sine/cosine, and IEEE remquo live here.
#include <mutex>
#include <vector>
#include "math.hpp"

namespace base256 {

// nearbyint lives in float256.hpp (current rounding mode, no inexact).

[[nodiscard]] inline float256 fdim(float256 x, float256 y) {
    if (x.isnan() || y.isnan()) {
        if (x.issignaling() || y.issignaling()) detail::raise_flag(fe_invalid);
        return x.isnan() ? x.quieted() : y.quieted();
    }
    if (x > y) return x - y;
    return float256::zero(false);
}

// logb lives in float256.hpp (constexpr). A second definition here made
// <base256/base256.hpp> ill-formed.


namespace detail_spec {
// c_n = 2^n / (1·3·…·(2n+1)) for the erf series, as dd (computed once).
inline constexpr int kErfTerms = 400;
[[nodiscard]] inline const std::array<detail_dd::dd, kErfTerms>& erf_coefs() {
    static const std::array<detail_dd::dd, kErfTerms> t = [] {
        const fenv::scoped_default_env g;
        std::array<detail_dd::dd, kErfTerms> c{};
        c[0] = detail_dd::dd(float256::one());
        for (int n = 1; n < kErfTerms; ++n)
            c[static_cast<std::size_t>(n)] =
                detail_dd::ldexp(c[static_cast<std::size_t>(n - 1)], 1) / float256::from_int(2 * n + 1);
        return c;
    }();
    return t;
}

// Taylor expansions of erfc at c_j = j/16, j = 1..256 (MPFR-generated erfc(c_j)
// and K_j = 2/√π e^{-c_j²}, see tools/gen_tables.c):
//   erfc(c + h) = erfc(c) - K Σ_{n>=0} b_n h^{n+1},  b_n = a_n / (n+1),
// where a_n are the Taylor coefficients of e^{-2ch-h²}: a_0 = 1, a_1 = -2c,
// a_{n+1} = (-2c a_n - 2 a_{n-1}) / (n+1). Coefficients are built on first use
// of each center (|h| <= 1/32); the evaluation truncates by the actual |h|.
struct ErfCenter {
    std::once_flag once;
    std::vector<detail_dd::dd> b;
};
[[nodiscard]] inline const std::vector<detail_dd::dd>& erf_center_coefs(int j) {
    static ErfCenter centers[detail_dd::tables::kErfCenters];
    ErfCenter& ce = centers[j];
    std::call_once(ce.once, [&] {
        using detail_dd::dd;
        const fenv::scoped_default_env g;
        const float256 m2c = -ldexp(float256::from_int(j), -3);  // -2c, exact
        dd am1{float256::zero()}, a0{float256::one()};
        const float256 lim = ldexp(float256::one(), -512);  // |b_n| 2^{-5n} below this: done
        int small = 0;
        for (int n = 0; n < 400; ++n) {
            const dd bn = a0 / float256::from_int(n + 1);
            ce.b.push_back(bn);
            const bool tiny = bn.hi.iszero() || abs(ldexp(bn.hi, -5 * n)) < lim;
            small = tiny ? small + 1 : 0;
            if (small >= 2) break;
            const dd an1 = (a0 * m2c - detail_dd::ldexp(am1, 1)) / float256::from_int(n + 1);
            am1 = a0;
            a0 = an1;
        }
    });
    return ce.b;
}

// erfc(a) for 1/32 <= a < 16 in dd: the expansion at the nearest c_j = j/16
// (h = a - c_j is exact, |h| <= 1/32; K e^{...}/erfc stays O(1), so at most a
// couple of bits cancel).
[[nodiscard]] inline detail_dd::dd erfc_center_dd(float256 a) {
    using detail_dd::dd;
    const int j = static_cast<int>(std::lround(a.to_double() * 16.0));
    const float256 h = a - ldexp(float256::from_int(j), -4);
    namespace T = detail_dd::tables;
    const dd ec = detail_dd::dd_from_words(T::kErfcVal[j]);
    if (h.iszero()) return ec;
    const auto& b = erf_center_coefs(j);
    // Terms: stop once |b_n h^n| < 2^-510 for two consecutive n.
    std::size_t N = b.size();
    const int eh = ilogb(h) + 1;  // |h| < 2^eh
    int small = 0;
    for (std::size_t n = 0; n < b.size(); ++n) {
        const bool tiny = b[n].hi.iszero() || ilogb(b[n].hi) + eh * static_cast<int>(n) < -510;
        small = tiny ? small + 1 : 0;
        if (small >= 2) { N = n + 1; break; }
    }
    dd p = b[N - 1];
    for (std::size_t n = N - 1; n-- > 0;) p = p * h + b[n];
    return ec - detail_dd::dd_from_words(T::kErfK[j]) * (p * h);
}

// erf(a) for 0 < a < 16 in dd. Below 1/32: 2/√π e^{-a²} a Σ c_n a^{2n} by
// Horner on the precomputed c_n (all terms positive); above, 1 - erfc(a).
[[nodiscard]] inline detail_dd::dd erf_series_dd(float256 a) {
    if (a >= ldexp(float256::one(), -5)) return detail_dd::dd(float256::one()) - erfc_center_dd(a);
    using detail_dd::dd;
    const dd x2 = detail_dd::two_prod(a, a);
    // Number of terms: c_n a^{2n} < 2^-496 (the sum is >= 1).
    const double l2 = std::log2(x2.hi.to_double());
    double lc = 0;
    int n = 1;
    for (; n < kErfTerms - 1; ++n) {
        lc += 1.0 - std::log2(2.0 * n + 1.0);
        if (lc + n * l2 < -496.0) break;
    }
    const auto& c = erf_coefs();
    dd p = c[static_cast<std::size_t>(n)];
    for (int k = n - 1; k >= 0; --k) p = p * x2 + c[static_cast<std::size_t>(k)];
    const auto& cs = detail_dd::consts();
    const dd two_over_sqrtpi = dd(float256::two()) / dd(cs.sqrtpi[0], cs.sqrtpi[1]);
    return p * a * detail_dd::exp(-x2) * two_over_sqrtpi;
}

// erfc(a) for a >= 16 (any a >= 4 works) as m·2^k: e^{-a²} / (√π · K(a)), with the Laplace
// continued fraction K = a + (1/2)/(a + 1/(a + (3/2)/(a + …))) truncated at a
// depth chosen from a, evaluated by the forward
// convergent recurrence P_n = a P_{n-1} + (n/2) P_{n-2} (likewise Q): all
// terms are positive, and only one division is needed.
[[nodiscard]] inline detail_dd::Scaled erfc_cf_dd(float256 a) {
    using detail_dd::dd;
    // depth ~ 6400/a² + 50 keeps the truncation error below 2^-300 (measured
    // against depth 3000 for a = 4 … 40).
    const double ad = a.to_double();
    const int depth = static_cast<int>(6400.0 / (ad * ad)) + 50;
    dd p2{float256::one()}, p1{a}, q2{float256::zero()}, q1{float256::one()};
    for (int n = 1; n <= depth; ++n) {
        const float256 bn = float256::from_int(n) * float256::half();
        const dd p = p1 * a + p2 * bn;
        const dd q = q1 * a + (q2.hi.iszero() ? q2 : q2 * bn);
        p2 = p1; p1 = p; q2 = q1; q1 = q;
        if ((n & 31) == 0) {  // keep the convergents in range
            const int e = ilogb(q1.hi);
            p1 = detail_dd::ldexp(p1, -e); p2 = detail_dd::ldexp(p2, -e);
            q1 = detail_dd::ldexp(q1, -e); q2 = detail_dd::ldexp(q2, -e);
        }
    }
    const dd f = p1 / q1;
    const auto& c = detail_dd::consts();
    detail_dd::Scaled e = detail_dd::exp_scaled(-detail_dd::two_prod(a, a));
    e.m = e.m / (f * dd(c.sqrtpi[0], c.sqrtpi[1]));
    return e;
}
} // namespace detail_spec

[[nodiscard]] inline float256 erf(float256 x) {
    detail_dd::FnEnv env;
    using detail_dd::dd;
    if (x.isnan()) return env.nan(x);
    if (x.isinf()) return env.exact(copysign(float256::one(), x));
    if (x.iszero()) return env.exact(x);
    const float256 a = abs(x);
    const int s = x.signbit() ? -1 : 1;
    dd y;
    int tail = 2;
    if (detail_math::is_tiny(a)) {
        // 2a/√π (1 - a²/3): the omitted a³ term points toward zero.
        const auto& c = detail_dd::consts();
        y = dd(a + a) / dd(c.sqrtpi[0], c.sqrtpi[1]);
        tail = -s;
    } else if (a < float256(16)) {
        y = detail_spec::erf_series_dd(a);
    } else {
        y = dd(float256::one());  // erfc(16) < 2^-370
        tail = -s;
    }
    if (s < 0) y = -y;
    return env.round(y.hi, y.lo, 0, tail);
}

[[nodiscard]] inline float256 erfc(float256 x) {
    detail_dd::FnEnv env;
    using detail_dd::dd;
    if (x.isnan()) return env.nan(x);
    if (x.isinf()) return env.exact(x.signbit() ? float256::two() : float256::zero());
    if (x.iszero()) return env.exact(float256::one());
    const float256 a = abs(x);
    const auto& c = detail_dd::consts();
    if (x.signbit()) {
        // 1 + erf(a)
        if (a > float256(16)) return env.round(float256::two(), float256::zero(), 0, -1);
        if (detail_math::is_tiny(a))
            return (dd(float256::one()) + dd(a + a) / dd(c.sqrtpi[0], c.sqrtpi[1])).out(env, 0, -1);
        if (a < float256(16)) return (dd(float256::one()) + detail_spec::erf_series_dd(a)).out(env);
        const detail_dd::Scaled e = detail_spec::erfc_cf_dd(a);
        return (dd(float256::two()) - detail_dd::ldexp(e.m, e.k)).out(env);
    }
    if (detail_math::is_tiny(a))
        return (dd(float256::one()) - dd(a + a) / dd(c.sqrtpi[0], c.sqrtpi[1])).out(env, 0, +1);
    if (a < ldexp(float256::one(), -5)) return (dd(float256::one()) - detail_spec::erf_series_dd(a)).out(env);
    if (a < float256(16)) return detail_spec::erfc_center_dd(a).out(env);
    if (a > float256(520)) return env.underflow(false);  // erfc(520) < 2^-390000
    const detail_dd::Scaled e = detail_spec::erfc_cf_dd(a);
    return e.m.out(env, e.k);
}

[[nodiscard]] inline float256 asinh(float256 x) {
    detail_dd::FnEnv env;
    using detail_dd::dd;
    if (x.isnan()) return env.nan(x);
    if (x.isinf() || x.iszero()) return env.exact(x);
    const float256 a = abs(x);
    if (detail_math::is_tiny(a)) return detail_math::tiny_odd(env, x, -1, 6, +1);  // x - x^3/6 (+3x^5/40)
    dd y;
    if (a > ldexp(float256::one(), 300)) {
        y = detail_dd::log(dd(a)) + detail_dd::ln2();  // log(2a); 1/(4a^2) is negligible
    } else {
        const dd r = detail_dd::sqrt(detail_dd::two_prod(a, a) + dd(float256::one()));
        y = detail_dd::log(dd(a) + r);
    }
    if (x.signbit()) y = -y;
    return y.out(env);
}

[[nodiscard]] inline float256 acosh(float256 x) {
    detail_dd::FnEnv env;
    using detail_dd::dd;
    if (x.isnan()) return env.nan(x);
    if (x < float256::one()) return env.invalid();
    if (x == float256::one()) return env.exact(float256::zero());
    if (x.isinf()) return env.exact(x);
    if (x > ldexp(float256::one(), 300)) return (detail_dd::log(dd(x)) + detail_dd::ln2()).out(env);
    const dd xm1 = detail_dd::two_sum(x, -float256::one());
    const dd xp1 = detail_dd::two_sum(x, float256::one());
    return detail_dd::log(dd(x) + detail_dd::sqrt(xm1 * xp1)).out(env);
}

[[nodiscard]] inline float256 atanh(float256 x) {
    detail_dd::FnEnv env;
    using detail_dd::dd;
    if (x.isnan()) return env.nan(x);
    if (x.iszero()) return env.exact(x);
    const float256 a = abs(x);
    if (a > float256::one()) return env.invalid();
    if (a == float256::one()) return env.pole(x.signbit());
    if (detail_math::is_tiny(a)) return detail_math::tiny_odd(env, x, 1, 3, +1);  // x + x^3/3 (+x^5/5)
    // (log(1 + a) - log(1 - a)) / 2, with 1 ± a exact as dd.
    const dd lp = detail_dd::log(detail_dd::two_sum(float256::one(), a));
    const dd lm = detail_dd::log(detail_dd::two_sum(float256::one(), -a));
    dd y = detail_dd::ldexp(lp - lm, -1);
    if (x.signbit()) y = -y;
    return y.out(env);
}


namespace detail_gamma {

// Stirling series in dd for ln Γ(w), w >= kStirlingShift:
//   (w-1/2) ln w - w + ln(2π)/2 + Σ_k B_{2k} / (2k(2k-1) w^{2k-1})
// Coefficients c_k = B_{2k}/(2k(2k-1)), k = 1..140, are MPFR dd values in
// detail/tables.hpp. The asymptotic series reaches 2^-505 relative to the
// result for w >= kStirlingShift before its terms start to grow.
#ifndef BASE256_STIRLING_SHIFT
#define BASE256_STIRLING_SHIFT 64
#endif
inline constexpr int kStirlingShift = BASE256_STIRLING_SHIFT;
inline constexpr int kStirlingTerms = detail_dd::tables::kStirlingN;
inline constexpr int kStirlingTailBits = 270;

[[nodiscard]] inline detail_dd::dd stirling_coef(int k) {
    return detail_dd::dd_from_words(detail_dd::tables::kStirling[k]);
}

// Σ_k c_k / w^(2k+1) (the Stirling correction), by Horner in 1/w² with the
// number of terms chosen from log2(w) (each omitted term below 2^-505 of `ref`).
[[nodiscard]] inline detail_dd::dd stirling_series(const detail_dd::dd& w, int ref_exp) {
    using detail_dd::dd;
    const double lw = std::log2(w.hi.to_double()) * (1.0 - 1e-12);  // 1/w <= 2^-lw
    // K terms in all; terms from H on are below 2^-270 of `ref` and are
    // summed in float256 (their rounding errors stay below 2^-500 of `ref`).
    int K = kStirlingTerms, H = -1;
    for (int k = 0; k < kStirlingTerms; ++k) {
        const double e = static_cast<double>(ilogb(stirling_coef(k).hi) + 1) - (2 * k + 1) * lw;
        if (H < 0 && e < ref_exp - kStirlingTailBits) H = k;
        if (e < ref_exp - 505) { K = k; break; }
    }
    if (H < 0) H = K;
    if (K == 0) return dd(float256::zero());
    const dd inv = dd(float256::one()) / w;
    const dd inv2 = inv * inv;
    dd p;
    int k = K - 1;
    if (H < K) {
        float256 t = stirling_coef(k).hi;
        for (--k; k >= H; --k) t = t * inv2.hi + stirling_coef(k).hi;
        p = k < 0 ? dd(t) : dd(t) * inv2 + stirling_coef(k);
    } else {
        p = stirling_coef(k);
    }
    for (--k; k >= 0; --k) p = p * inv2 + stirling_coef(k);
    return p * inv;
}

// (w - 1/2) ln w - w + ln √(2π) + Σ c_k / w^(2k+1)  (ln Γ(w), w >= 128).
[[nodiscard]] inline detail_dd::dd stirling_dd(const detail_dd::dd& w) {
    using detail_dd::dd;
    const auto& cs = detail_dd::consts();
    const dd lw = detail_dd::log(w);
    const dd sum = (w - dd(float256::half())) * lw - w + dd(cs.half_ln_2pi[0], cs.half_ln_2pi[1]);
    return sum + stirling_series(w, ilogb(sum.hi));
}

// z (z+1) … (z+n-1) as dd, pairing (z+i)(z+n-1-i) = z(z+n-1) + i(n-1-i).
[[nodiscard]] inline detail_dd::dd rising(const detail_dd::dd& z, int n) {
    using detail_dd::dd;
    if (n <= 0) return dd(float256::one());
    if (n == 1) return z;
    const dd q = z * (z + dd(float256::from_int(n - 1)));
    dd prod = q;
    int i = 1;
    for (; i < n - 1 - i; ++i) prod = prod * (q + dd(float256::from_int(static_cast<std::int64_t>(i) * (n - 1 - i))));
    if (i == n - 1 - i) prod = prod * (z + dd(float256::from_int(i)));
    return prod;
}

// ln Γ(z) for a positive dd z: shift to z + n >= 128 and divide out the
// rising product z (z+1) … (z+n-1), itself kept as dd.
[[nodiscard]] inline detail_dd::dd lgamma_pos_dd(const detail_dd::dd& z) {
    using detail_dd::dd;
    if (z.hi >= float256(kStirlingShift)) return stirling_dd(z);
    const int n = kStirlingShift - static_cast<int>(floor(z.hi).to_double());
    return stirling_dd(z + dd(float256::from_int(n))) - detail_dd::log(rising(z, n));
}

// Γ(z) as m·2^k for a positive dd z.
[[nodiscard]] inline detail_dd::Scaled gamma_abs_scaled_pos(const detail_dd::dd& z) {
    if (z.hi >= float256(kStirlingShift)) return detail_dd::exp_scaled(stirling_dd(z));
    const int n = kStirlingShift - static_cast<int>(floor(z.hi).to_double());
    detail_dd::Scaled e = detail_dd::exp_scaled(stirling_dd(z + detail_dd::dd(float256::from_int(n))));
    e.m = e.m / rising(z, n);
    return e;
}

// sin(π g) for g in (0, 1/2], as dd. Positive.
[[nodiscard]] inline detail_dd::dd sin_pi_dd(float256 g) {
    detail_dd::dd s, c;
    if (g <= float256(0.25)) {
        detail_dd::sincos_reduced(detail_dd::pi() * g, s, c);
        return s;
    }
    detail_dd::sincos_reduced(detail_dd::pi() * (float256::half() - g), s, c);
    return c;
}

// ln|Γ(x)| for finite x that is not a pole, as dd.
[[nodiscard]] inline detail_dd::dd lgamma_dd(float256 x) {
    using detail_dd::dd;
    if (!x.signbit()) return lgamma_pos_dd(dd(x));
    // ln|Γ(x)| = ln π − ln|sin(πx)| − ln Γ(1−x); 1 − x is exact as dd.
    const float256 f = x - round(x);
    const dd sp = sin_pi_dd(abs(f));
    return detail_dd::log(detail_dd::pi()) - detail_dd::log(sp) -
           lgamma_pos_dd(detail_dd::two_sum(float256::one(), -x));
}

// ln Γ(1 + e) and ln Γ(2 + e) for |e| < 2^-100 (three Taylor terms; the
// next one is below 2^-300 relative).
[[nodiscard]] inline detail_dd::dd lgamma_near_1_2(float256 e, bool at_two) {
    using detail_dd::dd;
    using detail_dd::hexc;
    const dd a1 = at_two
        ? dd(hexc("0x6.c3b981c824f385b2e41c07efead34a95e3133c509a33fe6f3fc20cb8f64p-4"),
             hexc("0x1.004271b4a605fc560f112f9b63349defa82eefa9516ecdeca5f71bc4b99p-240"))
        : -dd(hexc("0x9.3c467e37db0c7a4d1be3f810152cb56a1cecc3af65cc0190c03df347098p-4"),
              hexc("0x2.ffbd8e4b59fa03a9f0eed0649ccb621057d11056ae9132135a08e43b468p-240"));
    const float256 a2 = at_two ? hexc("0x5.28d3312983e991873d8912200bace5eb5c88daa81162d8b33582c031a6p-4")
                               : hexc("0xd.28d3312983e991873d8912200bace5eb5c88daa81162d8b33582c031a6p-4");
    const float256 a3 = at_two ? hexc("0x1.13e001a5576068125d0741df7857fb3b2ae7f3b26dbc01f45c3428e5b9dp-4")
                               : hexc("0x6.693556faacb5bd67b25c9734cdad5090803d4907c3115749b1897e3b0f4p-4");
    const float256 e2 = e * e;
    return a1 * e + dd(a2 * e2) - dd(a3 * e2 * e);
}

// |Γ(x)| as m·2^k for finite x that is not a pole, without going through ln:
// Γ(z) = exp(Stirling(z + n)) / (z (z+1) … (z+n-1)) for z > 0, and the
// reflection |Γ(x)| = π / (|sin(πx)| Γ(1 - x)) for x < 0.
[[nodiscard]] inline detail_dd::Scaled gamma_abs_scaled(float256 x) {
    using detail_dd::dd;
    if (x.signbit()) {
        const float256 f = x - round(x);
        const dd sp = sin_pi_dd(abs(f));
        detail_dd::Scaled e = gamma_abs_scaled_pos(detail_dd::two_sum(float256::one(), -x));
        e.m = detail_dd::pi() / (sp * e.m);
        e.k = -e.k;
        return e;
    }
    return gamma_abs_scaled_pos(dd(x));
}

} // namespace detail_gamma

// Natural log of |Γ(x)|. Poles (±0 and non-positive integers) return +∞ with divide-by-zero.
[[nodiscard]] inline float256 lgamma(float256 x) {
    detail_dd::FnEnv env;
    if (x.isnan()) return env.nan(x);
    if (x.isinf()) return env.exact(float256::inf(false));
    if (x.iszero()) return env.pole(false);
    if (x.signbit() && trunc(x) == x) return env.pole(false); // non-positive integer
    if (x == float256::one() || x == float256::two()) return env.exact(float256::zero());
    const float256 tiny = ldexp(float256::one(), -100);
    // Taylor at 1 and 2; the omitted ε^4 term is positive.
    if (abs(x - float256::one()) < tiny) return detail_gamma::lgamma_near_1_2(x - float256::one(), false).out(env, 0, +1);
    if (abs(x - float256::two()) < tiny) return detail_gamma::lgamma_near_1_2(x - float256::two(), true).out(env, 0, +1);
    if (x > ldexp(float256::one(), 262100)) return env.overflow(false);  // x ln x > max
    return detail_gamma::lgamma_dd(x).out(env);
}

// Γ(x). Negative integers and −∞ are invalid; tgamma(±0) is ±∞ (divide-by-zero).
[[nodiscard]] inline float256 tgamma(float256 x) {
    detail_dd::FnEnv env;
    if (x.isnan()) return env.nan(x);
    if (x.iszero()) return env.pole(x.signbit());
    if (x.isinf()) return x.signbit() ? env.invalid() : env.exact(x);
    if (x.signbit() && trunc(x) == x) return env.invalid();
    if (x > float256(20000)) return env.overflow(false);  // Γ(20000) > 2^262144
    // Γ(n) = (n-1)! is exact when representable: dd product, exact iff lo == 0.
    if (!x.signbit() && trunc(x) == x && x <= float256(300)) {
        detail_dd::dd p{float256::one()};
        const int n = static_cast<int>(x.to_double());
        for (int i = 2; i < n; ++i) p = p * float256::from_int(i);
        return env.round(p.hi, p.lo, 0, p.lo.iszero() ? 0 : 2);
    }
    bool neg = false;
    if (x.signbit()) {
        // sign of Γ(x) = sign of sin(πx) = (-1)^k sign(x - k), k = round(x)
        const float256 k = round(x);
        const bool odd = trunc(k * float256::half()) != k * float256::half();
        neg = odd != (x - k).signbit();
        if (x < float256(-20000)) return env.underflow(neg);  // |Γ(x)| < 2^-262378
    }
    detail_dd::Scaled e = detail_gamma::gamma_abs_scaled(x);
    if (neg) e.m = -e.m;
    return e.m.out(env, e.k);
}

// sin and cos of the same argument, one reduction.
inline void sincos(float256 x, float256* sin_out, float256* cos_out) {
    if (sin_out == nullptr && cos_out == nullptr) return;
    if (x.isnan() || x.isinf() || x.iszero() || detail_math::is_tiny(abs(x))) {
        // Same results and flags as the separate calls.
        if (sin_out) *sin_out = sin(x);
        if (cos_out) *cos_out = cos(x);
        return;
    }
    detail_dd::FnEnv env;
    int quad = 0;
    detail_dd::dd sd, cd;
    detail_math::trig_dd(abs(x), sd, cd, quad);
    detail_dd::dd sv = (quad & 1) ? cd : sd;
    detail_dd::dd cv = (quad & 1) ? -sd : cd;
    if (quad >= 2) { sv = -sv; cv = -cv; }
    if (x.signbit()) sv = -sv;
    const float256 so = sv.out(env);  // restores the caller's environment
    if (sin_out) *sin_out = so;
    if (cos_out) *cos_out = float256::round_dd(cv.hi, cv.lo, 0, env.mode(), 2);
}

// IEEE remainder, plus the low 31 bits of the integer quotient (sign of x/y).
/// IEEE remainder plus the sign and low 31 bits of the rounded quotient.
[[nodiscard]] constexpr float256 remquo(float256 x, float256 y, int* quo) noexcept {
    return float256::remainder_impl(x, y, true, quo);
}

} // namespace base256
