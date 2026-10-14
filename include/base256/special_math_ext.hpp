#pragma once
/// More special functions (included by base256.hpp): the C++17 orthogonal
/// polynomials (legendre, assoc_legendre, sph_legendre, laguerre,
/// assoc_laguerre, hermite) and the exponential integrals expint (Ei) and
/// expint(n, x) (E_n), the sine/cosine integrals sinint, cosint, sinhint,
/// coshint and the logarithmic integral logint. Each evaluates in double-float256 (~474 bits)
/// and rounds once in the caller's mode with IEEE flags. See README "Special
/// functions" for measured accuracy, domains and costs.

#include <algorithm>
#include <array>

#include "special_math.hpp"

namespace base256 {

namespace detail_sf {
namespace poly {

[[nodiscard]] inline float256 U(std::uint64_t v) { return float256::from_uint(v); }

// Keep a recurrence pair near 1 in magnitude; the scale goes into k.
inline void rescale(dd& p0, dd& p1, int& k) {
    const float256 big = abs(p1.hi) > abs(p0.hi) ? p1.hi : p0.hi;
    if (big.iszero() || !big.isfinite()) return;
    const int e = ilogb(big);
    if (e > 1000 || e < -1000) {
        p0 = detail_dd::ldexp(p0, -e);
        p1 = detail_dd::ldexp(p1, -e);
        k += e;
    }
}

// Round a polynomial value. Polynomials at binary256 points often have exact
// dyadic values, and some of them are binary256 numbers or midpoints between
// two. A dd value within 2^-440 (relative) of such a number is taken as exactly
// that number (the dd error is ~n·2^-474), so exact results raise no inexact
// flag and exact ties round as ties.
[[nodiscard]] inline float256 out_snap(detail_dd::FnEnv& env, const SD& v, bool neg_zero = false) {
    if (v.m.hi.iszero()) return env.exact(neg_zero ? -float256::zero() : float256::zero());
    if (v.k > -260000 && v.k < 260000) {
        const float256 r = v.m.hi + v.m.lo;
        const dd d = v.m - dd(r);
        const int lim = ilogb(r) - 440;
        if (d.hi.iszero() || ilogb(d.hi) < lim) return env.round(r, float256::zero(), v.k, 0);
        // half an ulp on the side of d (a quarter below a power of two)
        float256 h = ldexp(float256::one(), ilogb(r) - 237);
        if (d.hi.signbit() != r.signbit() && abs(r) == ldexp(float256::one(), ilogb(r))) h = ldexp(h, -1);
        const dd e = detail_sf::dabs(d) - dd(h);
        if (e.hi.iszero() || ilogb(e.hi) < lim) return env.round(r, copysign(h, d.hi), v.k, 0);
    }
    return out(env, v);
}

// P_l(x), any real x (three-term recurrence; stable on [-1, 1], dominant outside).
[[nodiscard]] inline SD legendre_sd(unsigned l, const dd& x) {
    if (l == 0) return {dd(float256::one()), 0};
    dd p0{float256::one()}, p1 = x;
    int k = 0;
    for (unsigned j = 1; j < l; ++j) {
        const dd p2 = (x * p1 * U(2ull * j + 1) - p0 * U(j)) / U(j + 1ull);
        p0 = p1; p1 = p2;
        rescale(p0, p1, k);
    }
    return norm({p1, k});
}

// H_n(x) (physicists'): H_{j+1} = 2x H_j − 2j H_{j−1}.
[[nodiscard]] inline SD hermite_sd(unsigned n, float256 x) {
    if (n == 0) return {dd(float256::one()), 0};
    const dd x2 = detail_dd::ldexp(dd(x), 1);
    dd p0{float256::one()}, p1 = x2;
    int k = 0;
    for (unsigned j = 1; j < n; ++j) {
        const dd p2 = x2 * p1 - detail_dd::ldexp(p0 * U(j), 1);
        p0 = p1; p1 = p2;
        rescale(p0, p1, k);
    }
    return norm({p1, k});
}

// L_n^(m)(x): (j+1) L_{j+1} = (2j+1+m−x) L_j − (j+m) L_{j−1}.
[[nodiscard]] inline SD laguerre_sd(unsigned n, unsigned m, float256 x) {
    if (n == 0) return {dd(float256::one()), 0};
    dd p0{float256::one()}, p1 = detail_dd::two_sum(U(1ull + m), -x);
    int k = 0;
    for (unsigned j = 1; j < n; ++j) {
        const dd c = detail_dd::two_sum(U(2ull * j + 1 + m), -x);
        const dd p2 = (c * p1 - p0 * U(static_cast<std::uint64_t>(j) + m)) / U(j + 1ull);
        p0 = p1; p1 = p2;
        rescale(p0, p1, k);
    }
    return norm({p1, k});
}

// P_l^m(x) without the Condon–Shortley phase, given x and s = sqrt(1 − x²) >= 0:
// P_m^m = (2m−1)!! s^m, then the recurrence in l.
[[nodiscard]] inline SD assoc_legendre_sd(unsigned l, unsigned m, const dd& x, const dd& s) {
    if (m > l) return {dd(float256::zero()), 0};
    if (m == 0) return legendre_sd(l, x);
    dd p0{float256::one()};
    int k = 0;
    for (unsigned i = 1; i <= m; ++i) {
        p0 = s * p0 * U(2ull * i - 1);
        dd z{float256::zero()};
        rescale(z, p0, k);
    }
    if (l == m) return norm({p0, k});
    dd p1 = x * p0 * U(2ull * m + 1);
    for (unsigned j = m + 1; j < l; ++j) {
        const dd p2 = (x * p1 * U(2ull * j + 1) - p0 * U(static_cast<std::uint64_t>(j) + m)) / U(j - m + 1ull);
        p0 = p1; p1 = p2;
        rescale(p0, p1, k);
    }
    return norm({p1, k});
}

} // namespace poly
} // namespace detail_sf

/// Legendre polynomial P_l(x) for any real x (C++17 requires |x| <= 1; this
/// extends it by the same polynomial).
[[nodiscard]] inline float256 legendre(unsigned l, float256 x) {
    detail_dd::FnEnv env;
    if (x.isnan()) return env.nan(x);
    if (l == 0) return env.exact(float256::one());
    if (x.isinf()) return env.exact((l & 1) ? x : float256::inf(false));
    return detail_sf::poly::out_snap(env, detail_sf::poly::legendre_sd(l, detail_dd::dd(x)), (l & 1) && x.signbit());
}

/// Associated Legendre function P_l^m(x) = (1−x²)^{m/2} d^m/dx^m P_l(x) for
/// |x| <= 1, without the Condon–Shortley phase (C++17). m > l gives 0;
/// |x| > 1 is a domain error.
[[nodiscard]] inline float256 assoc_legendre(unsigned l, unsigned m, float256 x) {
    detail_dd::FnEnv env;
    if (x.isnan()) return env.nan(x);
    if (abs(x) > float256::one()) return env.invalid();
    if (m > l) return env.exact(float256::zero());
    if (m == 0) return detail_sf::poly::out_snap(env, detail_sf::poly::legendre_sd(l, detail_dd::dd(x)), (l & 1) && x.signbit());
    using detail_dd::dd;
    const dd s2 = detail_dd::two_sum(float256::one(), -x) * detail_dd::two_sum(float256::one(), x);
    const dd s = detail_dd::sqrt(s2);
    return detail_sf::poly::out_snap(env, detail_sf::poly::assoc_legendre_sd(l, m, dd(x), s), ((l - m) & 1) && x.signbit());
}

/// Spherical harmonic Y_l^m(θ, 0) = (−1)^m sqrt((2l+1)/(4π) (l−m)!/(l+m)!) P_l^m(cos θ)
/// (C++17 sph_legendre), any finite θ; it depends on θ through cos θ and |sin θ|.
/// m > l gives 0.
[[nodiscard]] inline float256 sph_legendre(unsigned l, unsigned m, float256 theta) {
    detail_dd::FnEnv env;
    if (theta.isnan()) return env.nan(theta);
    if (theta.isinf()) return env.invalid();
    if (m > l) return env.exact(float256::zero());
    using detail_dd::dd;
    using detail_sf::SD;
    dd s, c;
    detail_sf::sincos_x(theta, s, c);
    s = detail_sf::dabs(s);
    SD p = detail_sf::poly::assoc_legendre_sd(l, m, c, s);
    if (p.m.hi.iszero()) return env.exact((m & 1) ? -float256::zero() : float256::zero());
    // N² = (2l+1) / (4π · (l+m)!/(l−m)!)
    SD r{dd(float256::one()), 0};
    for (unsigned j = l - m + 1; j <= l + m && j > l - m; ++j) r = detail_sf::mul(r, dd(detail_sf::poly::U(j)));
    SD n2 = detail_sf::norm({dd(detail_sf::poly::U(2ull * l + 1)) / (detail_dd::ldexp(detail_dd::pi(), 2) * r.m), -r.k});
    if (n2.k & 1) { n2.m = detail_dd::ldexp(n2.m, 1); --n2.k; }
    p = detail_sf::mul(p, SD{detail_dd::sqrt(n2.m), n2.k / 2});
    if (m & 1) p.m = -p.m;
    return detail_sf::out(env, p);
}

/// Laguerre polynomial L_n(x), any real x (C++17 requires x >= 0).
[[nodiscard]] inline float256 laguerre(unsigned n, float256 x) {
    detail_dd::FnEnv env;
    if (x.isnan()) return env.nan(x);
    if (n == 0) return env.exact(float256::one());
    if (x.isinf()) return env.exact((x.signbit() || !(n & 1)) ? float256::inf(false) : float256::inf(true));
    return detail_sf::poly::out_snap(env, detail_sf::poly::laguerre_sd(n, 0, x));
}

/// Associated Laguerre polynomial L_n^(m)(x) = (−1)^m d^m/dx^m L_{n+m}(x), any real x.
[[nodiscard]] inline float256 assoc_laguerre(unsigned n, unsigned m, float256 x) {
    detail_dd::FnEnv env;
    if (x.isnan()) return env.nan(x);
    if (n == 0) return env.exact(float256::one());
    if (x.isinf()) return env.exact((x.signbit() || !(n & 1)) ? float256::inf(false) : float256::inf(true));
    return detail_sf::poly::out_snap(env, detail_sf::poly::laguerre_sd(n, m, x));
}

/// Hermite polynomial H_n(x) (physicists'), any real x.
[[nodiscard]] inline float256 hermite(unsigned n, float256 x) {
    detail_dd::FnEnv env;
    if (x.isnan()) return env.nan(x);
    if (n == 0) return env.exact(float256::one());
    if (x.isinf()) return env.exact((n & 1) ? x : float256::inf(false));
    return detail_sf::poly::out_snap(env, detail_sf::poly::hermite_sd(n, x), (n & 1) && x.signbit());
}

namespace detail_sf {
namespace eint {

// The series and fractions below take a float256 or (for li) a dd argument.
[[nodiscard]] inline float256 H(float256 v) { return v; }
[[nodiscard]] inline float256 H(const dd& v) { return v.hi; }
[[nodiscard]] inline dd add(float256 a, float256 b) { return detail_dd::two_sum(a, b); }
[[nodiscard]] inline dd add(const dd& a, float256 b) { return a + dd(b); }

// S(x) = Σ_{k>=1} x^k/(k·k!), so Ei(x) = γ + ln|x| + S(x) and E1(t) = −γ − ln t − S(−t).
// Once the terms decrease (k > |x|) and fall below 2^-245 of the sum, the rest
// is summed in float256.
template <class T>
[[nodiscard]] dd s_series(const T& x) {
    dd term{float256::one()}, sum{float256::zero()};
    const float256 ax = abs(H(x)), xh = H(x);
    std::uint64_t k = 1;
    for (;; ++k) {
        const float256 K = float256::from_uint(k);
        term = term * x / K;
        const dd ad = term / K;
        sum = sum + ad;
        if (K > ax && f256_tail(ad, sum)) break;
    }
    float256 tf = term.hi, st{};
    for (++k;; ++k) {
        const float256 K = float256::from_uint(k);
        tf = tf * xh / K;
        const float256 ad = tf / K;
        st = st + ad;
        if (small_rel(dd(ad), sum)) break;
    }
    return sum + dd(st);
}

inline constexpr int kEiAsym = 360;  // Ei asymptotic series from here (min term < 2^-500)

// Ei(x) for x >= kEiAsym: e^x/x Σ k!/x^k.
template <class T>
[[nodiscard]] SD ei_asym(const T& x) {
    const float256 xh = H(x);
    dd term{float256::one()}, sum{float256::one()};
    std::uint64_t k = 1;
    for (;; ++k) {
        term = term * float256::from_uint(k) / x;
        sum = sum + term;
        if (f256_tail(term, sum)) break;
    }
    float256 tf = term.hi, st{};
    for (++k;; ++k) {
        tf = tf * float256::from_uint(k) / xh;
        st = st + tf;
        if (small_rel(dd(tf), sum)) break;
    }
    sum = sum + dd(st);
    SD e = exp_sd(dd(x));
    return mul(e, sum / x);
}

// The continued fraction E_n(x) = e^{-x} / (b0 + a1/(b1 + a2/(b2 + ...))),
// b_k = x + n + 2k, a_k = −k(n − 1 + k) (division-free Wallis form, rescaled).
template <class T>
[[nodiscard]] dd en_cf(std::uint64_t n, const T& x) {
    const float256 nn = float256::from_uint(n);
    dd A0{float256::one()}, B0{float256::zero()};
    dd A1 = add(x, nn), B1{float256::one()};
    dd h = A1, hprev = h;
    for (std::uint64_t k = 1;; ++k) {
        const float256 K = float256::from_uint(k);
        const dd b = add(x, nn + K + K);
        const dd a = -(detail_dd::two_prod(K, nn - float256::one() + K));
        const dd A2 = b * A1 + a * A0, B2 = b * B1 + a * B0;
        A0 = A1; B0 = B1; A1 = A2; B1 = B2;
        if ((k & 7) == 0) {
            const int e = ilogb(A1.hi);
            if (e > 1000 || e < -1000) {
                A0 = detail_dd::ldexp(A0, -e); A1 = detail_dd::ldexp(A1, -e);
                B0 = detail_dd::ldexp(B0, -e); B1 = detail_dd::ldexp(B1, -e);
            }
            h = A1 / B1;
            if (near_one(h / hprev)) break;
            hprev = h;
        }
    }
    return h;
}

// E_n(x) for x > 0 by the series
//   E_n(x) = (−x)^{n−1}/(n−1)! (ψ(n) − ln x) − Σ_{k≠n−1} (−x)^k / (k! (k−n+1)),
// ψ(n) = −γ + H_{n−1}. Used for x < 40 (it cancels about 2.9x bits).
template <class T>
[[nodiscard]] dd en_series(std::uint64_t n, const T& x) {
    const float256 xh = H(x);
    const dd one{float256::one()};
    dd sum{float256::zero()}, term = one;  // term = (−x)^k/k!
    const T nx = -x;
    const float256 nxh = -xh;
    std::uint64_t k = 0;
    for (;; ++k) {
        const float256 K = float256::from_uint(k);
        if (k > 0) term = term * nx / K;
        if (k + 1 != n) {
            const dd ad = term / (K - float256::from_uint(n - 1));
            sum = sum - ad;
            if (K > xh && k + 1 > n && f256_tail(ad, sum)) break;
        } else {
            dd psi = -euler_gamma();
            for (std::uint64_t j = 1; j < n; ++j) psi = psi + one / float256::from_uint(j);
            sum = sum + term * (psi - detail_dd::log(dd(x)));
        }
    }
    // float256 tail (k + 1 > n here, so k − n + 1 > 0)
    float256 tf = term.hi, st{};
    for (++k;; ++k) {
        const float256 K = float256::from_uint(k);
        tf = tf * nxh / K;
        const float256 ad = tf / (K - float256::from_uint(n - 1));
        st = st - ad;
        if (small_rel(dd(ad), sum)) break;
    }
    return sum + dd(st);
}

// Measured costs (µs): the fraction about 20000/(x + n/8), the series about
// 120 + 15x + 3n. The series is only used below x = 40.
[[nodiscard]] inline bool en_use_cf(std::uint64_t n, float256 x) {
    if (x >= float256(40)) return true;
    const double xd = x.to_double(), nd = static_cast<double>(n);
    return 20000.0 / (xd + nd / 8) < 120 + 15 * xd + 3 * nd;
}

// E_n(x), x > 0, as m·2^k.
template <class T>
[[nodiscard]] SD en_sd(std::uint64_t n, const T& x) {
    if (n == 0) return mul(exp_sd(-dd(x)), dd(float256::one()) / x);  // e^{-x}/x
    if (en_use_cf(n, H(x))) {
        const dd f = en_cf(n, x);
        return mul(exp_sd(-dd(x)), dd(float256::one()) / f);
    }
    return norm({en_series(n, x), 0});
}

// Ei(x), x != 0 finite, as m·2^k.
template <class T>
[[nodiscard]] SD ei_sd(const T& x) {
    if (H(x).signbit()) {
        SD e = en_sd(1, -x);
        e.m = -e.m;
        return e;
    }
    if (H(x) >= float256(kEiAsym)) return ei_asym(x);
    return norm({euler_gamma() + detail_dd::log(dd(x)) + s_series(x), 0});
}

} // namespace eint
} // namespace detail_sf

/// Exponential integral Ei(x) = −PV∫_{−x}^∞ e^{−t}/t dt (C++17 expint), any
/// real x: Ei(0) = −∞ (pole), Ei(−∞) = −0, Ei(+∞) = +∞. Overflows above x ≈ 181704.
[[nodiscard]] inline float256 expint(float256 x) {
    detail_dd::FnEnv env;
    if (x.isnan()) return env.nan(x);
    if (x.iszero()) return env.pole(true);
    if (x.isinf()) return env.exact(x.signbit() ? -float256::zero() : x);
    if (x > float256(190000)) return env.overflow(false);
    if (x < float256(-190000)) return env.underflow(true);
    return detail_sf::out(env, detail_sf::eint::ei_sd(x));
}

/// Generalized exponential integral E_n(x) = ∫_1^∞ e^{−xt}/t^n dt (Boost's
/// expint(n, x)) for x >= 0; E_n(0) = 1/(n − 1) for n >= 2 and a pole for
/// n <= 1. E_1(x) for x < 0 is −Ei(−x) (the real part, as in GSL); other
/// n with x < 0 are a domain error.
[[nodiscard]] inline float256 expint(unsigned n, float256 x) {
    detail_dd::FnEnv env;
    if (x.isnan()) return env.nan(x);
    if (x.signbit() && !x.iszero()) {
        if (n != 1) return env.invalid();
        if (x.isinf()) return env.exact(-float256::inf(false));
        if (x < float256(-190000)) return env.overflow(true);
        const detail_sf::SD e = detail_sf::eint::ei_sd(-x);
        return detail_sf::out(env, {-e.m, e.k});
    }
    if (x.iszero()) {
        if (n <= 1) return env.pole(false);
        const detail_dd::dd q = detail_dd::dd(float256::one()) / float256::from_uint(n - 1ull);
        return env.round(q.hi, q.lo, 0, q.lo.iszero() ? 0 : 2);
    }
    if (x.isinf()) return env.exact(float256::zero());
    if (x > float256(190000)) return env.underflow(false);
    return detail_sf::out(env, detail_sf::eint::en_sd(n, x));
}

namespace detail_sf {
namespace trigint {

// Odd and even parts of Σ_{k>=1} c_k x^k/(k·k!), c_k = 1 (hyperbolic) or the
// signs of i^{k−1} (odd, Si) and i^k (even, Ci): odd = Shi or Si,
// even + γ + ln x = Chi or Ci.
inline void series(float256 x, bool trig, dd& odd, dd& even) {
    dd term{float256::one()};
    odd = even = dd(float256::zero());
    std::uint64_t k = 1;
    bool odone = false, edone = false;
    for (;; ++k) {
        const float256 K = float256::from_uint(k);
        term = term * x / K;
        dd ad = term / K;
        if (trig && ((k & 3) == 3 || (k & 3) == 2)) ad = -ad;
        if (k & 1) { odd = odd + ad; odone = K > x && f256_tail(ad, odd); }
        else { even = even + ad; edone = K > x && f256_tail(ad, even); }
        if (odone && edone) break;
    }
    float256 tf = term.hi, so{}, se{};
    odone = edone = false;
    for (++k;; ++k) {
        const float256 K = float256::from_uint(k);
        tf = tf * x / K;
        float256 ad = tf / K;
        if (trig && ((k & 3) == 3 || (k & 3) == 2)) ad = -ad;
        if (k & 1) { so = so + ad; odone = small_rel(dd(ad), odd); }
        else { se = se + ad; edone = small_rel(dd(ad), even); }
        if (odone && edone) break;
    }
    odd = odd + dd(so);
    even = even + dd(se);
}

inline constexpr int kSeriesX = 80;   // Si/Ci series below this (cancels about 1.44x bits)
inline constexpr int kAsymX = 360;    // f, g asymptotic series from here

// E1(ix) by the complex Wallis form of E1(z) = e^{-z}/(z+1 − 1/(z+3 − 4/(z+5 − ...)))
// gives Ci(x) = −Re E1(ix), Si(x) = π/2 + Im E1(ix).
inline void cf(float256 x, dd& si, dd& ci) {
    const dd zero{float256::zero()}, one{float256::one()};
    dd Ar0 = one, Ai0 = zero, Br0 = zero, Bi0 = zero;
    dd Ar1{float256::one()}, Ai1{x}, Br1 = one, Bi1 = zero;
    dd hr = one, hi = zero, pr = one, pi = zero;  // h = B/A
    for (std::uint64_t k = 1;; ++k) {
        const float256 K = float256::from_uint(k);
        const dd br{K + K + float256::one()};
        const dd a = -detail_dd::two_prod(K, K);
        // A2 = (br + i x) A1 + a A0
        const dd Ar2 = br * Ar1 - Ai1 * x + a * Ar0, Ai2 = br * Ai1 + Ar1 * x + a * Ai0;
        const dd Br2 = br * Br1 - Bi1 * x + a * Br0, Bi2 = br * Bi1 + Br1 * x + a * Bi0;
        Ar0 = Ar1; Ai0 = Ai1; Br0 = Br1; Bi0 = Bi1;
        Ar1 = Ar2; Ai1 = Ai2; Br1 = Br2; Bi1 = Bi2;
        if ((k & 7) == 0) {
            const float256 big = abs(Ar1.hi) > abs(Ai1.hi) ? abs(Ar1.hi) : abs(Ai1.hi);
            const int e = ilogb(big);
            if (e > 1000 || e < -1000) {
                for (dd* v : {&Ar0, &Ai0, &Br0, &Bi0, &Ar1, &Ai1, &Br1, &Bi1}) *v = detail_dd::ldexp(*v, -e);
            }
            const dd d = Ar1 * Ar1 + Ai1 * Ai1;
            hr = (Br1 * Ar1 + Bi1 * Ai1) / d;
            hi = (Bi1 * Ar1 - Br1 * Ai1) / d;
            const dd dr = hr - pr, di = hi - pi;
            const float256 m = abs(dr.hi) > abs(di.hi) ? abs(dr.hi) : abs(di.hi);
            const float256 sz = abs(hr.hi) > abs(hi.hi) ? abs(hr.hi) : abs(hi.hi);
            if (m.iszero() || ilogb(m) < ilogb(sz) - 468) break;
            pr = hr; pi = hi;
        }
    }
    // E1(ix) = e^{-ix} h = (cos x − i sin x)(hr + i hi)
    dd s, c;
    sincos_x(x, s, c);
    const dd er = c * hr + s * hi, ei = c * hi - s * hr;
    ci = -er;
    si = detail_dd::pio2() + ei;
}

// x >= kAsymX: Si = π/2 − f cos x − g sin x, Ci = f sin x − g cos x with
// f = F/x, g = G/x², F = Σ (−1)^k (2k)!/x^{2k}, G = Σ (−1)^k (2k+1)!/x^{2k}.
inline void asym(float256 x, dd& si, dd& ci, int& ck) {
    const dd one{float256::one()};
    dd F = one, G = one;
    const bool huge = ilogb(x) > 1000;
    if (!huge) {
        const dd ix2 = one / (dd(x) * x);
        dd tF = one, tG = one;
        for (std::uint64_t k = 1;; ++k) {
            const float256 a = float256::from_uint(2 * k - 1), b = float256::from_uint(2 * k), c = float256::from_uint(2 * k + 1);
            tF = -(tF * ix2 * a * b);
            tG = -(tG * ix2 * b * c);
            F = F + tF; G = G + tG;
            if (small_rel(tF, F) && small_rel(tG, G)) break;
        }
    }
    dd s, c;
    sincos_x(x, s, c);
    // Ci = (F sin x − (G/x) cos x)/x as m·2^ck (x may be near the top of the range).
    const dd gx = huge ? dd(float256::zero()) : G / x;
    const int ex = ilogb(x);
    const dd xm{ldexp(x, -ex)};
    const SD civ = norm({(F * s - gx * c) / xm, -ex});
    ci = civ.m; ck = civ.k;
    const SD sv = norm({(F * c + gx * s) / xm, -ex});
    si = detail_dd::pio2() - (sv.k < -1200 ? dd(float256::zero()) : detail_dd::ldexp(sv.m, sv.k));
}

} // namespace trigint
} // namespace detail_sf

/// Sine integral Si(x) = ∫_0^x sin t/t dt, any real x (odd; Si(±∞) = ±π/2).
[[nodiscard]] inline float256 sinint(float256 x) {
    detail_dd::FnEnv env;
    if (x.isnan()) return env.nan(x);
    if (x.iszero()) return env.exact(x);
    using detail_dd::dd;
    const bool neg = x.signbit();
    const float256 ax = abs(x);
    dd si, ci;
    if (ax.isinf()) si = detail_dd::pio2();
    else if (ax < float256(detail_sf::trigint::kSeriesX)) detail_sf::trigint::series(ax, true, si, ci);
    else if (ax < float256(detail_sf::trigint::kAsymX)) detail_sf::trigint::cf(ax, si, ci);
    else { int k; detail_sf::trigint::asym(ax, si, ci, k); }
    if (neg) si = -si;
    return env.round(si.hi, si.lo);
}

/// Cosine integral Ci(x) = γ + ln x + ∫_0^x (cos t − 1)/t dt for x >= 0:
/// Ci(0) = −∞ (pole), Ci(+∞) = +0; x < 0 is a domain error.
[[nodiscard]] inline float256 cosint(float256 x) {
    detail_dd::FnEnv env;
    if (x.isnan()) return env.nan(x);
    if (x.iszero()) return env.pole(true);
    if (x.signbit()) return env.invalid();
    if (x.isinf()) return env.exact(float256::zero());
    using detail_dd::dd;
    dd si, ci;
    int k = 0;
    if (x < float256(detail_sf::trigint::kSeriesX)) {
        detail_sf::trigint::series(x, true, si, ci);
        ci = ci + detail_sf::euler_gamma() + detail_dd::log(dd(x));
    } else if (x < float256(detail_sf::trigint::kAsymX)) detail_sf::trigint::cf(x, si, ci);
    else detail_sf::trigint::asym(x, si, ci, k);
    return detail_sf::out(env, detail_sf::norm({ci, k}));
}

/// Hyperbolic sine integral Shi(x) = ∫_0^x sinh t/t dt, any real x (odd).
/// Overflows above |x| ≈ 181705.
[[nodiscard]] inline float256 sinhint(float256 x) {
    detail_dd::FnEnv env;
    if (x.isnan()) return env.nan(x);
    if (x.iszero() || x.isinf()) return env.exact(x);
    using detail_dd::dd;
    const bool neg = x.signbit();
    const float256 ax = abs(x);
    if (ax > float256(190000)) return env.overflow(neg);
    detail_sf::SD v;
    if (ax >= float256(detail_sf::eint::kEiAsym)) {
        v = detail_sf::eint::ei_asym(ax);  // Shi = (Ei(x) + E1(x))/2, E1 below 2^-1000 relative
        --v.k;
    } else {
        dd o, e;
        detail_sf::trigint::series(ax, false, o, e);
        v = detail_sf::norm({o, 0});
    }
    if (neg) v.m = -v.m;
    return detail_sf::out(env, v);
}

/// Hyperbolic cosine integral Chi(x) = γ + ln x + ∫_0^x (cosh t − 1)/t dt for
/// x >= 0: Chi(0) = −∞ (pole); x < 0 is a domain error. Overflows above x ≈ 181705.
[[nodiscard]] inline float256 coshint(float256 x) {
    detail_dd::FnEnv env;
    if (x.isnan()) return env.nan(x);
    if (x.iszero()) return env.pole(true);
    if (x.signbit()) return env.invalid();
    if (x.isinf()) return env.exact(x);
    if (x > float256(190000)) return env.overflow(false);
    using detail_dd::dd;
    detail_sf::SD v;
    if (x >= float256(detail_sf::eint::kEiAsym)) {
        v = detail_sf::eint::ei_asym(x);  // Chi = (Ei(x) − E1(x))/2
        --v.k;
    } else {
        dd o, e;
        detail_sf::trigint::series(x, false, o, e);
        v = detail_sf::norm({e + detail_sf::euler_gamma() + detail_dd::log(dd(x)), 0});
    }
    return detail_sf::out(env, v);
}

/// Logarithmic integral li(x) = Ei(ln x) for x >= 0: li(0) = 0, li(1) = −∞
/// (pole), li(+∞) = +∞; x < 0 is a domain error.
[[nodiscard]] inline float256 logint(float256 x) {
    detail_dd::FnEnv env;
    if (x.isnan()) return env.nan(x);
    if (x.iszero()) return env.exact(x);
    if (x.signbit()) return env.invalid();
    if (x == float256::one()) return env.pole(true);
    if (x.isinf()) return env.exact(x);
    return detail_sf::out(env, detail_sf::eint::ei_sd(detail_dd::log(detail_dd::dd(x))));
}

namespace detail_sf {
namespace psi {

inline constexpr int kAsymW = 60;  // asymptotic series at w >= 60: the k = 140 term is below 2^-520

// ψ(w) = ln w − 1/(2w) − Σ_{k>=1} B_2k/(2k) w^{−2k}, w >= kAsymW (B_2k/(2k) = bern_fact(k)·(2k−1)!).
[[nodiscard]] inline dd asym(const dd& w) {
    const dd one{float256::one()};
    const dd iw = one / w, iw2 = iw * iw;
    dd sum{float256::zero()}, p = iw2, fact = one;  // fact = (2k−1)!
    for (int k = 1; k <= detail_dd::tables::kBernFactN; ++k) {
        if (k > 1) fact = fact * float256::from_int((2 * k - 2) * (2 * k - 1));
        const dd t = zeta::bern_fact(k) * fact * p;
        sum = sum + t;
        if (small_rel(t, sum) && k > 2) break;
        p = p * iw2;
    }
    return detail_dd::log(w) - detail_dd::ldexp(iw, -1) - sum;
}

// ψ(x) for x > 0: ψ(x) = ψ(x + N) − Σ_{j<N} 1/(x + j), the sum kept as one fraction P/Q.
[[nodiscard]] inline dd pos(const dd& x) {
    if (x.hi >= float256(kAsymW)) return asym(x);
    const int N = kAsymW - static_cast<int>(floor(x.hi).to_double());
    dd P{float256::zero()}, Q{float256::one()};
    for (int j = 0; j < N; ++j) {
        const dd t = x + dd(float256::from_int(j));
        P = P * t + Q;
        Q = Q * t;
        const int e = ilogb(Q.hi);
        if (e > 1000) { P = detail_dd::ldexp(P, -e); Q = detail_dd::ldexp(Q, -e); }
    }
    return asym(x + dd(float256::from_int(N))) - P / Q;
}

// Γ(a)/Γ(b), a, b > 0: gamma_sd below 2^20, else exp(lnΓ(a) − lnΓ(b)) (absolute
// error ~2^-474·lnΓ there, so correctly rounded only while lnΓ < 2^200 or so).
[[nodiscard]] inline SD gamma_ratio(const dd& a, const dd& b) {
    const float256 lim = ldexp(float256::one(), 20);
    if (a.hi < lim && b.hi < lim) {
        const SD ga = gamma_sd(a), gb = gamma_sd(b);
        return norm({ga.m / gb.m, ga.k - gb.k});
    }
    return exp_sd(detail_gamma::lgamma_pos_dd(a) - detail_gamma::lgamma_pos_dd(b));
}

} // namespace psi
} // namespace detail_sf

/// Digamma ψ(x) = Γ'(x)/Γ(x), any real x: ψ(±0) = ∓∞ (pole), negative
/// integers are a domain error, ψ(+∞) = +∞, ψ(−∞) is a domain error.
[[nodiscard]] inline float256 digamma(float256 x) {
    detail_dd::FnEnv env;
    using detail_dd::dd;
    if (x.isnan()) return env.nan(x);
    if (x.iszero()) return env.pole(!x.signbit());
    if (x.isinf()) return x.signbit() ? env.invalid() : env.exact(x);
    if (!x.signbit()) {
        if (x > ldexp(float256::one(), 2000)) return detail_dd::log(dd(x)).out(env);  // 1/(2x) below 2^-1990 relative
        return detail_sf::psi::pos(dd(x)).out(env);
    }
    if (x == floor(x)) return env.invalid();
    // ψ(x) = ψ(1 − x) − π cot(πx)
    dd s, c;
    detail_sf::sincos_pi(x, s, c);
    const dd y = detail_dd::two_sum(float256::one(), -x);
    const dd r = detail_sf::psi::pos(y) - detail_dd::pi() * c / s;
    return r.out(env);
}

/// Polygamma ψ^(n)(x) = (−1)^{n+1} n! ζ(n+1, x) for n >= 1 and x > 0
/// (n = 0 is digamma). x <= 0 is a domain error for n >= 1 (x = 0: pole).
[[nodiscard]] inline float256 polygamma(unsigned n, float256 x) {
    if (n == 0) return digamma(x);
    detail_dd::FnEnv env;
    using detail_dd::dd;
    if (x.isnan()) return env.nan(x);
    if (x.iszero()) return env.pole((n & 1) == 0);
    if (x.signbit()) return env.invalid();
    if (x.isinf()) return env.exact((n & 1) ? float256::zero() : -float256::zero());
    const float256 s = float256::from_uint(n + 1ull);
    detail_sf::SD z = detail_sf::zeta::hurwitz(dd(s), dd(float256::from_uint(n)), x);
    z = detail_sf::mul(z, detail_sf::gamma_sd(dd(s)));
    if ((n & 1) == 0) z.m = -z.m;
    return detail_sf::out(env, z);
}

/// Trigamma ψ'(x) = ζ(2, x) for x > 0.
[[nodiscard]] inline float256 trigamma(float256 x) { return polygamma(1, x); }

/// ln B(a, b) for a, b > 0 (zero argument: +∞ pole; negative: domain error).
[[nodiscard]] inline float256 lbeta(float256 a, float256 b) {
    detail_dd::FnEnv env;
    using detail_dd::dd;
    if (a.isnan() || b.isnan()) return env.nan2(a, b);
    if ((a.signbit() && !a.iszero()) || (b.signbit() && !b.iszero())) return env.invalid();
    if (a.iszero() || b.iszero()) return env.pole(false);
    if (a.isinf() || b.isinf()) return env.exact(-float256::inf(false));
    const float256 one = float256::one();
    if (a == one || b == one) {  // ln(1/b)
        const float256 o = a == one ? b : a;
        if (o == one) return env.exact(float256::zero());
        return (-detail_dd::log(dd(o))).out(env);
    }
    if (a > float256(65536) || b > float256(65536))  // B may leave the SD range: lnΓ(a) + lnΓ(b) − lnΓ(a+b)
        return (detail_gamma::lgamma_pos_dd(dd(a)) + detail_gamma::lgamma_pos_dd(dd(b)) -
                detail_gamma::lgamma_pos_dd(detail_dd::two_sum(a, b))).out(env);
    const detail_sf::SD v = detail_sf::beta_sd(a, b);
    const dd r = detail_dd::log(v.m) + detail_dd::ln2() * float256::from_int(v.k);
    return r.out(env);
}

/// Γ(a)/Γ(b) for a, b > 0 (Boost's tgamma_ratio).
[[nodiscard]] inline float256 tgamma_ratio(float256 a, float256 b) {
    detail_dd::FnEnv env;
    using detail_dd::dd;
    if (a.isnan() || b.isnan()) return env.nan2(a, b);
    if (a.signbit() || b.signbit() || (a.isinf() && b.isinf())) return env.invalid();
    if (b.iszero()) return a.iszero() ? env.invalid() : env.exact(float256::zero());
    if (a.iszero()) return env.pole(false);
    if (a.isinf()) return env.exact(a);
    if (b.isinf()) return env.exact(float256::zero());
    if (a == b) return env.exact(float256::one());
    return detail_sf::out(env, detail_sf::psi::gamma_ratio(dd(a), dd(b)));
}

/// Γ(a)/Γ(a + δ) for a > 0, a + δ > 0 (Boost's tgamma_delta_ratio); a + δ is
/// formed exactly.
[[nodiscard]] inline float256 tgamma_delta_ratio(float256 a, float256 delta) {
    detail_dd::FnEnv env;
    using detail_dd::dd;
    if (a.isnan() || delta.isnan()) return env.nan2(a, delta);
    if (a.signbit() || a.isinf() || delta.isinf()) return env.invalid();
    if (a.iszero()) return env.pole(false);
    if (delta.iszero()) return env.exact(float256::one());
    const dd b = detail_dd::two_sum(a, delta);
    if (b.hi.signbit() || b.hi.iszero()) return env.invalid();
    return detail_sf::out(env, detail_sf::psi::gamma_ratio(dd(a), b));
}

namespace detail_sf {
namespace airy {

struct Consts { dd c1, c2, sqrt3, rpi, rsqrtpi; };
// c1 = Ai(0) = 3^{-2/3}/Γ(2/3), c2 = −Ai'(0) = 3^{-1/3}/Γ(1/3).
inline const Consts& consts() {
    static const Consts c = [] {
        Consts r;
        const dd one{float256::one()}, three{float256(3)};
        const dd third = one / float256(3);
        const dd l3 = detail_dd::log(three);
        const SD g1 = gamma_sd(third), g2 = gamma_sd(one - third);
        r.c2 = detail_dd::exp(-(l3 * third)) / detail_dd::ldexp(g1.m, g1.k);
        r.c1 = detail_dd::exp(-(l3 * (one - third))) / detail_dd::ldexp(g2.m, g2.k);
        r.sqrt3 = detail_dd::sqrt(three);
        r.rpi = one / detail_dd::pi();
        r.rsqrtpi = detail_dd::sqrt(r.rpi);
        return r;
    }();
    return c;
}

// Maclaurin series: f = Σ x^{3k} 3^k (1/3)_k/(3k)!, g = Σ x^{3k+1} 3^k (2/3)_k/(3k+1)! and
// their derivatives; Ai = c1 f − c2 g, Bi = √3 (c1 f + c2 g).
inline void series(float256 x, dd& f, dd& g, dd& fp, dd& gp) {
    const dd x3 = detail_dd::two_prod(x, x) * x;
    const dd one{float256::one()};
    dd a = one, b{x}, ap = detail_dd::ldexp(detail_dd::two_prod(x, x), -1), bp = one;
    f = a; g = b; fp = dd(float256::zero()); gp = bp;  // the f' series starts at k = 1
    for (int k = 1; k < 100000; ++k) {
        const float256 K3 = float256::from_int(3 * k);
        if (k > 1) ap = ap * x3 / ((K3 - float256(3)) * (K3 - float256::one()));
        a = a * x3 / ((K3 - float256::one()) * K3);
        b = b * x3 / (K3 * (K3 + float256::one()));
        bp = bp * x3 / ((K3 - float256::two()) * K3);
        f = f + a; g = g + b; fp = fp + ap; gp = gp + bp;
        if (k > 2 && small_rel(a, f) && small_rel(b, g) && small_rel(ap, fp) && small_rel(bp, gp)) break;
    }
}

// K_{1/3}(z), K_{2/3}(z) for a dd z >= 2 (Steed's CF2 with μ = −1/3, as in
// bessel::steed_k), as m·2^kexp with a shared exponent.
inline void steed_k13(const dd& X, dd& k13, dd& k23, int& kexp) {
    const dd one{float256::one()};
    const dd mu = -(one / float256(3));
    const dd mu2 = mu * mu;
    dd b = (one + X) * float256::two();
    dd d = one / b;
    dd h = d, delh = d;
    dd q1{float256::zero()}, q2 = one;
    const dd a1 = dd(float256(0.25)) - mu2;
    dd q = a1, c = a1, a = -a1;
    dd s = one + q * delh;
    for (int i = 1; i < 1000000; ++i) {
        a = a - D(2 * i);
        c = -a * c / float256::from_int(i + 1);
        const dd qnew = (q1 - b * q2) / a;
        q1 = q2; q2 = qnew;
        q = q + c * qnew;
        b = b + D(2);
        d = one / (b + a * d);
        delh = (b * d - one) * delh;
        h = h + delh;
        const dd dels = q * delh;
        s = s + dels;
        if (small_rel(dels, s) && small_rel(delh, h)) break;
    }
    h = a1 * h;
    const SD ex = exp_sd(-X);
    const dd pre = detail_dd::sqrt(detail_dd::pi() / (X * float256::two())) * ex.m / s;
    k13 = pre;
    k23 = pre * (mu + X + dd(float256::half()) - h) / X;
    kexp = ex.k;
}

// Asymptotic sums Σ (±1)^k w_k ζ^{-k} split into even and odd k, w = u (deriv
// false) or v; u_k = u_{k−1}(6k−5)(6k−3)(6k−1)/((2k−1)·216k), v_k = −u_k(6k+1)/(6k−1).
// sign: 0 all +, 1 (−1)^k, 2 (−1)^{⌊k/2⌋}. Stops at 2^-490 or at the smallest term.
inline void asym_sums(const dd& zeta, bool deriv, int sign, dd& ev, dd& od) {
    const dd one{float256::one()};
    const dd iz = one / zeta;
    dd u = one, p = one;
    ev = one; od = dd(float256::zero());
    float256 prev = float256::inf(false);
    for (int k = 1; k < 100000; ++k) {
        const float256 a = float256::from_int(6 * k - 5) * float256::from_int(6 * k - 3) * float256::from_int(6 * k - 1);
        u = u * a / (float256::from_int(2 * k - 1) * float256::from_int(216 * k));
        p = p * iz;
        dd t = u * p;
        if (deriv) t = -(t * float256::from_int(6 * k + 1) / float256::from_int(6 * k - 1));
        if ((sign == 1 && (k & 1)) || (sign == 2 && (k & 2))) t = -t;
        const float256 at = abs(t.hi);
        if (at > prev) break;  // past the smallest term
        prev = at;
        if (k & 1) od = od + t; else ev = ev + t;
        if (small_rel(t, ev)) break;
    }
}

inline constexpr int kSerPos = 60;    // x > 0: Ai series while ζ < 60 (cancels ~2.9ζ bits)
inline constexpr int kAsym = 104;     // asymptotic for ζ >= 104 (smallest term ~e^{-2ζ} < 2^-300)
inline constexpr int kSerNeg = 110;   // x < 0: series while ζ < 110 (cancels ~1.44ζ bits)

// sin ζ, cos ζ for ζ = (2/3)t^{3/2} past dd range (t >= 2^140): ζ·(2/π) mod 4 by an
// exact multiprecision reduction with the 2/π table (395264 bits, enough for every
// finite binary256 x). False if ζ ever needed more of the table than it has.
[[nodiscard]] inline bool big_phase(float256 t, dd& sz, dd& cz) {
    namespace M = detail_mp;
    const int Q = 3 * ilogb(t) / 2 + 400;
    if (Q + 64 + 127 > 64 * detail_mp::kTwoOverPiAllWords) return false;
    const M::MP T = M::from_f256(t);
    const M::MP z = M::div_u(M::mul_u(M::mul(T, M::sqrt(T, Q), Q), 2, Q), 3, Q);
    const M::MP y = M::mul(z, M::two_over_pi(Q + 64), Q);
    const M::MP f = M::ldexp(zeta::zmp::mod2(M::ldexp(y, -1)), 1);  // in [−2, 2]
    const long n = std::lround(M::to_dd(f).hi.to_double());
    const M::MP r = M::sub(f, M::from_u64(static_cast<std::uint64_t>(n < 0 ? -n : n), n < 0), Q);
    const dd th = M::to_dd(M::mul(r, M::ldexp(M::pi(600), -1), 600));
    dd s1, c1, s2, c2;
    sincos_x(th.hi, s1, c1);
    sincos_x(th.lo, s2, c2);
    const dd s = s1 * c2 + c1 * s2, c = c1 * c2 - s1 * s2;
    switch (((n % 4) + 4) % 4) {
        case 0: sz = s; cz = c; break;
        case 1: sz = c; cz = -s; break;
        case 2: sz = -s; cz = -c; break;
        default: sz = -c; cz = s; break;
    }
    return true;
}

// which: 0 Ai, 1 Ai', 2 Bi, 3 Bi'. Finite x.
[[nodiscard]] inline SD eval(float256 x, int which) {
    const Consts& C = consts();
    const dd one{float256::one()};
    const bool deriv = which & 1, bi = which >= 2;
    const float256 t = abs(x);
    const dd sq = detail_dd::sqrt(dd(t));
    const dd zeta = sq * t * float256(2) / float256(3);
    const float256 zh = zeta.hi;
    const bool neg = x.signbit();
    const bool use_series = neg ? zh < float256(kSerNeg) : (bi ? zh < float256(kAsym) : zh < float256(kSerPos));
    if (use_series) {
        dd f, g, fp, gp;
        series(x, f, g, fp, gp);
        dd r;
        if (!deriv) r = bi ? C.sqrt3 * (C.c1 * f + C.c2 * g) : C.c1 * f - C.c2 * g;
        else r = bi ? C.sqrt3 * (C.c1 * fp + C.c2 * gp) : C.c1 * fp - C.c2 * gp;
        return norm({r, 0});
    }
    const dd q4 = detail_dd::sqrt(sq);  // |x|^{1/4}
    if (!neg) {
        if (!bi && zh < float256(kAsym)) {
            // Ai = (1/π)√(x/3) K_{1/3}(ζ), Ai' = −(x/(π√3)) K_{2/3}(ζ)
            dd k13, k23;
            int ke;
            steed_k13(zeta, k13, k23, ke);
            const dd r = deriv ? -(dd(x) * C.rpi / C.sqrt3 * k23) : C.rpi * sq / C.sqrt3 * k13;
            return norm({r, ke});
        }
        dd ev, od;
        asym_sums(zeta, deriv, bi ? 0 : 1, ev, od);
        const dd s = ev + od;
        const SD e = exp_sd(bi ? zeta : -zeta);
        dd pre = deriv ? q4 * C.rsqrtpi : C.rsqrtpi / q4;
        if (!bi) pre = detail_dd::ldexp(pre, -1);
        if (deriv && !bi) pre = -pre;
        return mul(e, pre * s);
    }
    // x < 0, oscillatory: θ = ζ − π/4.
    dd ev, od;
    asym_sums(zeta, deriv, 2, ev, od);  // ev = Σ(−1)^k w_{2k}ζ^{-2k}, od = Σ(−1)^k w_{2k+1}ζ^{-2k−1}
    dd sz, cz;
    if (zh >= ldexp(float256::one(), 140)) {
        if (!big_phase(t, sz, cz)) return SD{dd(float256::nan()), 0};  // past the table: NaN + invalid (unreachable for finite x)
    } else {
        dd s1, c1, s2, c2;
        sincos_x(zeta.hi, s1, c1);
        sincos_x(zeta.lo, s2, c2);
        sz = s1 * c2 + c1 * s2; cz = c1 * c2 - s1 * s2;
    }
    const dd r2 = detail_dd::sqrt(dd(float256::half()));
    const dd st = (sz - cz) * r2, ct = (cz + sz) * r2;  // sin θ, cos θ
    dd r;
    switch (which) {
        case 0: r = (ct * ev + st * od) * C.rsqrtpi / q4; break;
        case 2: r = (ct * od - st * ev) * C.rsqrtpi / q4; break;
        case 1: r = (st * ev - ct * od) * C.rsqrtpi * q4; break;
        default: r = (ct * ev + st * od) * C.rsqrtpi * q4; break;
    }
    return norm({r, 0});
}

[[nodiscard]] inline float256 pub(float256 x, int which) {
    detail_dd::FnEnv env;
    if (x.isnan()) return env.nan(x);
    if (x.isinf()) {
        if (!x.signbit()) {
            if (which == 0) return env.exact(float256::zero());
            if (which == 1) return env.exact(-float256::zero());
            return env.exact(x);
        }
        return (which & 1) ? env.invalid() : env.exact(float256::zero());  // Ai', Bi' grow like |x|^{1/4}
    }
    if (!x.signbit() && x > float256(6000)) {  // ζ > 3e5: beyond the range
        return which >= 2 ? env.overflow(false) : env.underflow(which == 1);
    }
    // the phase reduction needs 1.5·ilogb|x| + 600 bits of 2/π: NaN + invalid past the table
    // (the table covers every finite x, so this is a guard only)
    if (x.signbit() && 3 * ilogb(x) / 2 + 400 + 64 + 127 > 64 * detail_mp::kTwoOverPiAllWords) return env.invalid();
    return out(env, eval(x, which));
}

} // namespace airy
} // namespace detail_sf

/// Airy functions Ai(x), Bi(x) and their derivatives, any real x. For x < 0
/// they oscillate; for x below −2^93 the phase ζ = (2/3)|x|^{3/2} is reduced
/// mod π/2 exactly in multiprecision with a 395264-bit 2/π table (cost grows like
/// log2|x|²), so the results stay correctly rounded down to the most negative
/// finite value. Ai, Ai' underflow and Bi,
/// Bi' overflow above x ≈ 4205. Ai(±∞) = +0, Bi(+∞) = +∞, Bi(−∞) = +0;
/// Ai'(−∞), Bi'(−∞) are domain errors.
[[nodiscard]] inline float256 airy_ai(float256 x) { return detail_sf::airy::pub(x, 0); }
[[nodiscard]] inline float256 airy_ai_prime(float256 x) { return detail_sf::airy::pub(x, 1); }
[[nodiscard]] inline float256 airy_bi(float256 x) { return detail_sf::airy::pub(x, 2); }
[[nodiscard]] inline float256 airy_bi_prime(float256 x) { return detail_sf::airy::pub(x, 3); }

namespace detail_sf {
namespace inv {

// erfc(a) for finite a >= 0 as m·2^k (dd accuracy), from the elementary-function internals.
[[nodiscard]] inline SD erfc_sd(float256 a) {
    if (a < ldexp(float256::one(), -5)) return norm({dd(float256::one()) - detail_spec::erf_series_dd(a), 0});
    if (a < float256(16)) return norm({detail_spec::erfc_center_dd(a), 0});
    const detail_dd::Scaled e = detail_spec::erfc_cf_dd(a);
    return norm({e.m, e.k});
}
[[nodiscard]] inline dd two_over_sqrtpi() {
    const auto& c = detail_dd::consts();
    return dd(float256::two()) / dd(c.sqrtpi[0], c.sqrtpi[1]);
}

// erf^{-1}(x) for 0 < x <= 1/2: second-order Newton steps from a double guess,
// each from a float256 point w (the dd erf is exact-argument there):
//   r = erf(w) − x, δ1 = −r/erf'(w), δ = δ1 + w δ1² (erf''/erf' = −2w).
[[nodiscard]] inline dd erfinv_small(float256 x) {
    const double xd = x.to_double();
    double w = 0.886226925452758 * xd;
    for (int i = 0; i < 6; ++i) w -= (std::erf(w) - xd) / (1.1283791670955126 * std::exp(-w * w));
    dd W{float256(w)};
    for (int it = 0; it < 2; ++it) {
        const float256 w0 = W.round();
        const dd e = w0 < ldexp(float256::one(), -5) ? detail_spec::erf_series_dd(w0)
                                                     : dd(float256::one()) - detail_spec::erfc_center_dd(w0);
        const dd fp = two_over_sqrtpi() * detail_dd::exp(-detail_dd::two_prod(w0, w0));
        const dd d1 = -(e - dd(x)) / fp;
        W = dd(w0) + d1 + d1 * d1 * w0;
    }
    return W;
}

// erfc^{-1}(y) for 0 < y < 1/2 (w > 0.47) in the log form g(w) = ln(erfc(w)/y):
//   R = erf'(w)/erfc(w), δ1 = g/R, δ = δ1 − (R − 2w) δ1²/2.
[[nodiscard]] inline dd erfcinv_small(float256 y) {
    double w;
    const double ly = static_cast<double>(log(y).to_double());  // ln y, > −181800
    {
        // std::erfc below w = 25, else ln erfc(w) ≈ −w² − ln(w√π) + ln(1 − 1/(2w²) + 3/(4w⁴))
        // (relative error ~2^-50 there).
        w = std::sqrt(std::max(-ly - 0.5 * std::log(-3.14159 * ly), 0.25));
        for (int i = 0; i < 40; ++i) {
            double g, gp;
            if (w < 25) {
                const double e = std::erfc(w);
                g = std::log(e) - ly;
                gp = -1.1283791670955126 * std::exp(-w * w) / e;
            } else {
                const double w2 = w * w, s = 1 - 0.5 / w2 + 0.75 / (w2 * w2);
                g = -w2 - std::log(w * 1.772453850905516) + std::log(s) - ly;
                gp = -2 * w - 1 / w + (1 / (w2 * w) - 3 / (w2 * w2 * w)) / s;
            }
            const double d = g / gp;
            w -= d;
            if (std::abs(d) < 1e-17 * w) break;
        }
    }
    dd W{float256(w)};
    for (int it = 0; it < 3; ++it) {
        const float256 w0 = W.round();
        const SD E = erfc_sd(w0);
        // E/y and R = (2/√π) e^{-w²} / E without leaving the exponent range
        const int ey = ilogb(y);
        const dd q = detail_dd::ldexp(E.m / ldexp(y, -ey), E.k - ey);
        const SD ex = exp_sd(-detail_dd::two_prod(w0, w0));
        const dd R = detail_dd::ldexp(two_over_sqrtpi() * ex.m / E.m, ex.k - E.k);
        const dd g = detail_dd::log(q);
        const dd d1 = g / R;
        W = dd(w0) + d1 - (R - dd(w0 + w0)) * d1 * d1 * float256::half();
    }
    return W;
}

// W(x) by Halley steps on f = w e^w − x (near w = −1 and for |w| <~ 2) or
// Newton steps on w + ln|w| = ln|x| (large |w|), from a double or float256 guess.
[[nodiscard]] inline dd lambert(float256 x, bool minus1) {
    const dd one{float256::one()};
    const dd e1 = detail_dd::exp(one);
    const dd ep = e1 * x + one;  // e·x + 1 >= 0 (checked by the caller)
    double w;
    const double xd = x.to_double();
    const bool huge = abs(x) > float256(1e300) || (abs(x) < float256(1e-300));
    if (ep.hi < float256(0.25)) {
        const double p = std::sqrt(2 * ep.hi.to_double()) * (minus1 ? -1 : 1);
        w = -1 + p - p * p / 3 + 11.0 / 72 * p * p * p - 43.0 / 540 * p * p * p * p;
    } else if (huge) {
        const double L1 = log(abs(x)).to_double(), L2 = std::log(std::abs(L1));
        w = L1 - L2 + L2 / L1;
        if (!minus1 && x.signbit()) w = xd;  // unreachable (tiny negative x handled by the caller)
    } else if (!minus1) {
        w = xd < 3 ? std::log1p(xd) * (1 - std::log1p(std::log1p(xd)) / (2 + std::log1p(xd))) : std::log(xd) - std::log(std::log(xd));
        if (xd < 0.1 && xd > -0.1) w = xd - xd * xd;
    } else {
        const double L1 = std::log(-xd), L2 = std::log(-L1);
        w = L1 - L2 + L2 / L1;
    }
    if (!huge && ep.hi >= float256(0.25)) {
        for (int i = 0; i < 30; ++i) {  // Halley in double
            const double ew = std::exp(w), f = w * ew - xd;
            const double d = f / (ew * (w + 1) - (w + 2) * f / (2 * w + 2));
            w -= d;
            if (std::abs(d) <= 1e-16 * std::abs(w)) break;
        }
    }
    dd W{float256(w)};
    const bool logform = std::abs(w) > 3;
    const dd L = logform ? detail_dd::log(dd(abs(x))) : dd(float256::zero());
    for (int it = 0; it < 12; ++it) {
        dd d;
        if (logform) {  // w ← w (1 + L − ln|w|)/(1 + w)
            const dd lw = detail_dd::log(detail_sf::dabs(W));
            d = W * (W + lw - L) / (W + one);
        } else {
            const dd ew = detail_dd::exp(W), f = W * ew - dd(x), w1 = W + one;
            d = f / (ew * w1 - (W + dd(float256::two())) * f / (w1 + w1));
        }
        W = W - d;
        if (d.hi.iszero() || ilogb(d.hi) < ilogb(W.hi) - 470) break;
    }
    return W;
}

} // namespace inv
} // namespace detail_sf

/// Inverse error function: erf_inv(x) for −1 <= x <= 1 (±1 give ±∞ with a
/// pole; |x| > 1 is a domain error).
[[nodiscard]] inline float256 erf_inv(float256 x) {
    detail_dd::FnEnv env;
    using detail_dd::dd;
    if (x.isnan()) return env.nan(x);
    if (x.iszero()) return env.exact(x);
    const float256 a = abs(x), one = float256::one();
    if (a > one) return env.invalid();
    if (a == one) return env.pole(x.signbit());
    dd r;
    if (a < ldexp(one, -300)) {  // √π/2·x (1 + πx²/12): the omitted term is away from zero
        const auto& c = detail_dd::consts();
        r = dd(c.sqrtpi[0], c.sqrtpi[1]) * float256::half() * a;
        if (x.signbit()) r = -r;
        return env.round(r.hi, r.lo, 0, x.signbit() ? -1 : 1);
    }
    r = a <= float256::half() ? detail_sf::inv::erfinv_small(a) : detail_sf::inv::erfcinv_small(one - a);
    if (x.signbit()) r = -r;
    return r.out(env);
}

/// Inverse complementary error function erfc_inv(y) for 0 <= y <= 2
/// (0 → +∞ and 2 → −∞ with a pole; outside is a domain error).
[[nodiscard]] inline float256 erfc_inv(float256 y) {
    detail_dd::FnEnv env;
    using detail_dd::dd;
    if (y.isnan()) return env.nan(y);
    const float256 one = float256::one(), two = float256::two();
    if (y.signbit() && !y.iszero()) return env.invalid();
    if (y > two) return env.invalid();
    if (y.iszero()) return env.pole(false);
    if (y == two) return env.pole(true);
    if (y == one) return env.exact(float256::zero());
    const bool neg = y > one;
    const float256 t = neg ? two - y : y;  // exact (Sterbenz)
    dd r;
    if (t < float256::half()) r = detail_sf::inv::erfcinv_small(t);
    else {
        const float256 x = one - t;  // exact
        if (x < ldexp(one, -300)) {
            const auto& c = detail_dd::consts();
            r = dd(c.sqrtpi[0], c.sqrtpi[1]) * float256::half() * x;
            return env.round(neg ? -r.hi : r.hi, neg ? -r.lo : r.lo, 0, neg ? -1 : 1);
        }
        r = detail_sf::inv::erfinv_small(x);
    }
    if (neg) r = -r;
    return r.out(env);
}

/// Lambert W, principal branch W0(x) for x >= −1/e (W0(−1/e) = −1 is not a
/// binary256 argument; x < −1/e is a domain error), W0(+∞) = +∞.
[[nodiscard]] inline float256 lambert_w0(float256 x) {
    detail_dd::FnEnv env;
    using detail_dd::dd;
    if (x.isnan()) return env.nan(x);
    if (x.iszero() || (x.isinf() && !x.signbit())) return env.exact(x);
    if (x.isinf()) return env.invalid();
    const dd ep = detail_dd::exp(dd(float256::one())) * x + dd(float256::one());
    if (ep.hi.signbit()) return env.invalid();
    if (abs(x) < ldexp(float256::one(), -300)) return env.round(x, float256::zero(), 0, -1);  // x − x² + …
    return detail_sf::inv::lambert(x, false).out(env);
}

/// Lambert W, lower branch W−1(x) for −1/e <= x < 0: W−1(−0) = −∞ (pole);
/// x >= 0 (other than −0) and x < −1/e are domain errors.
[[nodiscard]] inline float256 lambert_wm1(float256 x) {
    detail_dd::FnEnv env;
    using detail_dd::dd;
    if (x.isnan()) return env.nan(x);
    if (x.iszero()) return env.pole(true);
    if (!x.signbit() || x.isinf()) return env.invalid();
    const dd ep = detail_dd::exp(dd(float256::one())) * x + dd(float256::one());
    if (ep.hi.signbit()) return env.invalid();
    return detail_sf::inv::lambert(x, true).out(env);
}

namespace detail_sf {
namespace fres {

// sin(πt), cos(πt) for t = x²/2 given exactly as a dd (exact reduction of both parts).
inline void sincos_half_x2(float256 x, dd& s, dd& c) {
    const dd t = detail_dd::ldexp(detail_dd::two_prod(x, x), -1);
    dd s1, c1, s2, c2;
    sincos_pi(t.hi, s1, c1);
    sincos_pi(t.lo, s2, c2);
    s = s1 * c2 + c1 * s2;
    c = c1 * c2 - s1 * s2;
}

// Series: with z = πx²/2, T_n = x z^n/n!: C = Σ_{n even} (−1)^{n/2} T_n/(2n+1),
// S = Σ_{n odd} (−1)^{(n−1)/2} T_n/(2n+1). Cancels about 2.27x² bits.
inline void series(float256 x, dd& C, dd& S) {
    const dd z = detail_dd::two_prod(x, x) * detail_dd::pio2();
    dd T{x};
    C = T; S = dd(float256::zero());
    for (int n = 1; n < 100000; ++n) {
        T = T * z / float256::from_int(n);
        dd t = T / float256::from_int(2 * n + 1);
        if (n & 2) t = -t;
        if (n & 1) S = S + t; else C = C + t;
        if (n > 4 && small_rel(t, (n & 1) ? S : C) && float256::from_int(n) > z.hi) {
            // one more term of the other sum
            T = T * z / float256::from_int(n + 1);
            dd u = T / float256::from_int(2 * n + 3);
            if ((n + 1) & 2) u = -u;
            if ((n + 1) & 1) S = S + u; else C = C + u;
            break;
        }
    }
}

// Complex Laplace continued fraction of erfc(w), w = √π (1 − i) x/2:
//   C + iS = (1+i)/2 (1 − e^{−w²}/(√π K)), K = w + (1/2)/(w + 1/(w + (3/2)/(w + …))),
// e^{−w²} = e^{iπx²/2}. Wallis recurrences A_n = w A_{n−1} + (n/2) A_{n−2}.
inline void cf(float256 x, dd& C, dd& S) {
    const dd zero{float256::zero()}, one{float256::one()};
    const auto& cs = detail_dd::consts();
    const dd sqpi{cs.sqrtpi[0], cs.sqrtpi[1]};
    const dd wr = sqpi * x * float256::half(), wi = -wr;
    dd Ar0 = one, Ai0 = zero, Br0 = zero, Bi0 = zero;  // A_{-1} = 1, B_{-1} = 0
    dd Ar1 = wr, Ai1 = wi, Br1 = one, Bi1 = zero;      // A_0 = w, B_0 = 1
    dd hr = zero, hi = zero, pr = one, pi = zero;      // h = B/A
    for (int n = 1; n < 10000000; ++n) {
        const float256 a = float256::from_int(n) * float256::half();
        const dd Ar2 = wr * Ar1 - wi * Ai1 + Ar0 * a, Ai2 = wr * Ai1 + wi * Ar1 + Ai0 * a;
        const dd Br2 = wr * Br1 - wi * Bi1 + Br0 * a, Bi2 = wr * Bi1 + wi * Br1 + Bi0 * a;
        Ar0 = Ar1; Ai0 = Ai1; Br0 = Br1; Bi0 = Bi1;
        Ar1 = Ar2; Ai1 = Ai2; Br1 = Br2; Bi1 = Bi2;
        if ((n & 7) == 0) {
            const float256 big = abs(Ar1.hi) > abs(Ai1.hi) ? abs(Ar1.hi) : abs(Ai1.hi);
            const int e = ilogb(big);
            if (e > 1000 || e < -1000)
                for (dd* v : {&Ar0, &Ai0, &Br0, &Bi0, &Ar1, &Ai1, &Br1, &Bi1}) *v = detail_dd::ldexp(*v, -e);
            const dd d = Ar1 * Ar1 + Ai1 * Ai1;
            hr = (Br1 * Ar1 + Bi1 * Ai1) / d;
            hi = (Bi1 * Ar1 - Br1 * Ai1) / d;
            const dd dr = hr - pr, di = hi - pi;
            const float256 m = abs(dr.hi) > abs(di.hi) ? abs(dr.hi) : abs(di.hi);
            const float256 sz = abs(hr.hi) > abs(hi.hi) ? abs(hr.hi) : abs(hi.hi);
            if (m.iszero() || ilogb(m) < ilogb(sz) - 468) break;
            pr = hr; pi = hi;
        }
    }
    dd s, c;
    sincos_half_x2(x, s, c);
    // Q = e^{iπx²/2} h / √π
    const dd qr = (c * hr - s * hi) / sqpi, qi = (c * hi + s * hr) / sqpi;
    const dd u = one - qr, v = -qi;
    C = detail_dd::ldexp(u - v, -1);
    S = detail_dd::ldexp(u + v, -1);
}

// x >= kAsym: C = 1/2 + f sin(πx²/2) − g cos, S = 1/2 − f cos − g sin,
// f = (1/(πx)) Σ (−1)^k (4k−1)!!/(πx²)^{2k}, g = (1/(π²x³)) Σ (−1)^k (4k+1)!!/(πx²)^{2k}.
inline void asym(float256 x, dd& C, dd& S) {
    const dd one{float256::one()};
    const dd pix = detail_dd::pi() * x;
    const dd y = one / (pix * x);   // 1/(πx²)
    const dd y2 = y * y;
    dd F = one, G = one, tf = one, tg = one;
    float256 prev = float256::inf(false);
    for (int k = 1; k < 100000 && ilogb(x) < 1000; ++k) {
        tf = -(tf * y2 * float256::from_int(4 * k - 3) * float256::from_int(4 * k - 1));
        tg = -(tg * y2 * float256::from_int(4 * k - 1) * float256::from_int(4 * k + 1));
        if (abs(tf.hi) > prev) break;
        prev = abs(tf.hi);
        F = F + tf; G = G + tg;
        if (small_rel(tf, F) && small_rel(tg, G)) break;
    }
    const dd f = F / pix, g = G * y / pix;
    dd s, c;
    sincos_half_x2(x, s, c);
    const dd h{float256::half()};
    C = h + f * s - g * c;
    S = h - f * c - g * s;
}

inline constexpr double kSeriesX = 4.5, kAsymX = 12;

inline float256 pub(float256 x, bool sine) {
    detail_dd::FnEnv env;
    if (x.isnan()) return env.nan(x);
    if (x.iszero()) return env.exact(x);
    if (x.isinf()) return env.exact(copysign(float256::half(), x));
    const bool neg = x.signbit();
    const float256 a = abs(x);
    if (a < ldexp(float256::one(), -100)) {
        if (!sine) return env.round(x, float256::zero(), 0, neg ? 1 : -1);  // x − π²x⁵/40
        const int e = ilogb(a);
        const dd m{ldexp(a, -e)};
        SD v = norm({m * m * m * detail_dd::pi() / float256(6), 3 * e});  // πx³/6 − …
        if (neg) v.m = -v.m;
        return out(env, v, neg ? 1 : -1);
    }
    dd C, S;
    if (a < float256(kSeriesX)) series(a, C, S);
    else if (a < float256(kAsymX)) cf(a, C, S);
    else asym(a, C, S);
    dd r = sine ? S : C;
    if (neg) r = -r;
    return r.out(env);
}

// Dawson F(x) = e^{−x²} Σ x^{2k+1}/(k!(2k+1)) (positive terms) for |x| < 16;
// asymptotically (1/(2x)) Σ (2k−1)!!/(2x²)^k.
[[nodiscard]] inline SD dawson_sd(float256 a) {
    const dd one{float256::one()};
    if (a < float256(16)) {
        const dd x2 = detail_dd::two_prod(a, a);
        dd t{a}, sum{a};
        int k = 1;
        for (;; ++k) {
            t = t * x2 / float256::from_int(k);
            const dd ad = t / float256::from_int(2 * k + 1);
            sum = sum + ad;
            if (float256::from_int(k) > x2.hi && f256_tail(ad, sum)) break;
        }
        float256 tf = t.hi, st{};
        for (++k;; ++k) {
            tf = tf * x2.hi / float256::from_int(k);
            const float256 ad = tf / float256::from_int(2 * k + 1);
            st = st + ad;
            if (small_rel(dd(ad), sum)) break;
        }
        sum = sum + dd(st);
        return mul(exp_sd(-x2), sum);
    }
    const int e = ilogb(a);
    const dd m{ldexp(a, -e)};
    dd sum = one;
    if (e < 1000) {
        const dd y = one / (detail_dd::two_prod(a, a) * float256::two());
        dd t = one;
        float256 prev = float256::inf(false);
        for (int k = 1; k < 100000; ++k) {
            t = t * y * float256::from_int(2 * k - 1);
            if (abs(t.hi) > prev) break;
            prev = abs(t.hi);
            sum = sum + t;
            if (small_rel(t, sum)) break;
        }
    }
    return norm({sum / m, -e - 1});
}

} // namespace fres
} // namespace detail_sf

/// Fresnel integrals S(x) = ∫_0^x sin(πt²/2) dt and C(x) = ∫_0^x cos(πt²/2) dt,
/// any real x (odd; S(±∞) = C(±∞) = ±1/2).
[[nodiscard]] inline float256 fresnel_s(float256 x) { return detail_sf::fres::pub(x, true); }
[[nodiscard]] inline float256 fresnel_c(float256 x) { return detail_sf::fres::pub(x, false); }

/// Dawson's integral F(x) = e^{−x²} ∫_0^x e^{t²} dt, any real x (odd; F(±∞) = ±0).
[[nodiscard]] inline float256 dawson(float256 x) {
    detail_dd::FnEnv env;
    if (x.isnan()) return env.nan(x);
    if (x.iszero()) return env.exact(x);
    if (x.isinf()) return env.exact(copysign(float256::zero(), x));
    const bool neg = x.signbit();
    const float256 a = abs(x);
    if (a < ldexp(float256::one(), -300)) return env.round(x, float256::zero(), 0, neg ? 1 : -1);  // x − 2x³/3
    detail_sf::SD v = detail_sf::fres::dawson_sd(a);
    if (neg) v.m = -v.m;
    // above 2^1000 the sum is 1 + (below 2^-2000): the true value is just above 1/(2x) in magnitude
    return detail_sf::out(env, v, ilogb(a) >= 1000 ? (neg ? -1 : 1) : 2);
}

namespace detail_sf {
namespace plog {

inline constexpr int kZetaN = 64;
// ζ(j), j = 0..kZetaN (ζ(1) unused), computed once.
inline const std::array<dd, kZetaN + 1>& zeta_int() {
    static const std::array<dd, kZetaN + 1> z = [] {
        const fenv::scoped_default_env g;
        std::array<dd, kZetaN + 1> r{};
        r[0] = dd(-float256::half());
        for (int j = 2; j <= kZetaN; ++j) r[static_cast<std::size_t>(j)] = zeta::pos(D(j), D(j - 1));
        return r;
    }();
    return z;
}
// ζ(m) for integer m <= kZetaN, m != 1: ζ(−2k) = 0, ζ(1−2k) = −B_2k/(2k) = −bern_fact(k)(2k−1)!.
[[nodiscard]] inline dd zeta_at(int m) {
    if (m >= 0) return zeta_int()[static_cast<std::size_t>(m)];
    if ((-m) % 2 == 0) return dd(float256::zero());
    const int k = (1 - m) / 2;
    dd f{float256::one()};
    for (int i = 2; i <= 2 * k - 1; ++i) f = f * float256::from_int(i);
    return -(zeta::bern_fact(k) * f);
}

// Σ_{k>=1} x^k / k^n, |x| <= 1/2 or n large.
[[nodiscard]] inline dd series(unsigned n, const dd& x) {
    dd p = x, sum = x;
    for (std::uint64_t k = 2; k < 100000; ++k) {
        p = p * x;
        const float256 K = float256::from_uint(k);
        dd t = p;
        for (unsigned i = 0; i < n; ++i) t = t / K;
        sum = sum + t;
        if (small_rel(t, sum) || t.hi.iszero()) break;
    }
    return sum;
}

// Li_n(e^μ), −ln 2 <= μ <= 0 (x in [1/2, 1]), n >= 2, n <= kZetaN:
//   Σ_{k≠n−1} ζ(n−k) μ^k/k! + μ^{n−1}/(n−1)! (H_{n−1} − ln(−μ)).
[[nodiscard]] inline dd mu_exp(unsigned n, const dd& mu) {
    const dd one{float256::one()};
    if (mu.hi.iszero()) return zeta_at(static_cast<int>(n));
    dd sum{float256::zero()}, p = one;  // p = μ^k/k!
    for (int k = 0; k < 2000; ++k) {
        if (k > 0) p = p * mu / float256::from_int(k);
        if (p.hi.iszero() || k - static_cast<int>(n) > 2 * detail_dd::tables::kBernFactN - 2) break;  // |μ/2π|^k < 2^-750 by then
        dd t;
        if (k == static_cast<int>(n) - 1) {
            if (mu.hi.iszero()) continue;
            dd H{float256::zero()};
            for (unsigned j = 1; j < n; ++j) H = H + one / float256::from_uint(j);
            t = p * (H - detail_dd::log(-mu));
        } else {
            t = zeta_at(static_cast<int>(n) - k) * p;
            if (t.hi.iszero()) continue;  // trivial zeros ζ(−2j)
        }
        sum = sum + t;
        if (k > static_cast<int>(n) && small_rel(t, sum) && !p.hi.iszero()) break;
        if (p.hi.iszero()) break;
    }
    return sum;
}

// Li2(X) for −1 <= X <= 1/2 by the Bernoulli series in u = −ln(1 − X):
//   Li2 = u − u²/4 + Σ_k B_2k u^{2k+1}/(2k+1)!, |u| <= ln 2.
[[nodiscard]] inline dd li2_core(const dd& X) {
    if (X.hi.iszero()) return X;
    if (ilogb(X.hi) < -20) return series(2, X);
    const dd one{float256::one()};
    const dd u = -detail_dd::log(one - X);
    const dd u2 = u * u;
    dd sum = u - u2 * float256(0.25), p = u;
    for (int k = 1; k <= detail_dd::tables::kBernFactN; ++k) {
        p = p * u2;
        const dd t = zeta::bern_fact(k) * p / float256::from_int(2 * k + 1);
        sum = sum + t;
        if (small_rel(t, sum)) break;
    }
    return sum;
}

// Li_n(X) for −1 <= X <= 1, n >= 2.
[[nodiscard]] inline dd lin(unsigned n, const dd& X) {
    const dd one{float256::one()};
    const float256 half = float256::half();
    if (n == 2 && X.hi <= half) return li2_core(X);
    if (abs(X.hi) <= half || n > static_cast<unsigned>(kZetaN)) return series(n, X);
    if (X.hi.signbit()) {  // Li_n(x) = 2^{1−n} Li_n(x²) − Li_n(−x)
        return detail_dd::ldexp(lin(n, X * X), 1 - static_cast<int>(n)) - lin(n, -X);
    }
    if (n == 2) {  // Li2(x) = π²/6 − ln x ln(1−x) − Li2(1−x)
        const dd y = one - X;
        const dd pi2_6 = detail_dd::pi() * detail_dd::pi() / float256(6);
        if (y.hi.iszero()) return pi2_6;
        return pi2_6 - detail_dd::log(X) * detail_dd::log(y) - li2_core(y);
    }
    return mu_exp(n, detail_dd::log(X));
}

} // namespace plog
} // namespace detail_sf

/// Dilogarithm Li₂(x) = −∫_0^x ln(1−t)/t dt for real x <= 1 (x > 1 is complex:
/// domain error); Li₂(1) = π²/6, Li₂(−∞) = −∞.
[[nodiscard]] inline float256 dilog(float256 x) {
    detail_dd::FnEnv env;
    using detail_dd::dd;
    if (x.isnan()) return env.nan(x);
    if (x.iszero()) return env.exact(x);
    const float256 one = float256::one();
    if (x > one) return env.invalid();
    if (x.isinf()) return env.exact(x);
    if (abs(x) < ldexp(one, -300)) return env.round(x, float256::zero(), 0, 1);  // x + x²/4
    if (x >= -one) return detail_sf::plog::lin(2, dd(x)).out(env);
    // x < −1: Li2(x) = −π²/6 − ln²(−x)/2 − Li2(1/x)
    const dd l = detail_dd::log(dd(-x));
    const dd r = -(detail_dd::pi() * detail_dd::pi() / float256(6)) - l * l * float256::half() -
                 detail_sf::plog::li2_core(dd(one) / x);
    return r.out(env);
}

/// Polylogarithm Li_n(x) = Σ_{k>=1} x^k/k^n for integer n >= 0: n = 0, 1 for
/// any x < 1 (x/(1−x), −ln(1−x)), n = 2 as dilog, n >= 3 for −1 <= x <= 1
/// (|x| > 1 is not supported: domain error). Li_n(1) = ζ(n) (n >= 2; n <= 1 pole).
[[nodiscard]] inline float256 polylog(unsigned n, float256 x) {
    if (n == 2) return dilog(x);  // before FnEnv: dilog keeps its own mode and flags
    detail_dd::FnEnv env;
    using detail_dd::dd;
    if (x.isnan()) return env.nan(x);
    if (x.iszero()) return env.exact(x);
    const float256 one = float256::one();
    if (n <= 1) {
        if (x == one) return env.pole(false);
        if (x > one) return env.invalid();
        if (n == 0) {
            if (x.isinf()) return env.exact(-one);  // x/(1−x) → −1 from above
            return (dd(x) / detail_dd::two_sum(one, -x)).out(env);
        }
        if (x.isinf()) return env.exact(-x);
        if (abs(x) < ldexp(one, -300)) return env.round(x, float256::zero(), 0, 1);  // x + x²/2
        if (abs(x) <= float256::half()) return detail_sf::plog::series(1, dd(x)).out(env);
        return (-detail_dd::log(detail_dd::two_sum(one, -x))).out(env);
    }
    if (abs(x) > one) return env.invalid();
    if (abs(x) < ldexp(one, -300)) return env.round(x, float256::zero(), 0, x.signbit() ? -1 : 1);
    if (x == one && n > static_cast<unsigned>(detail_sf::plog::kZetaN))
        return env.round(one, float256::zero(), 0, 1);
    return detail_sf::plog::lin(n, dd(x)).out(env);
}

/// Dirichlet eta η(s) = Σ (−1)^{n−1} n^{−s} = (1 − 2^{1−s}) ζ(s), all real s
/// (η(1) = ln 2, η(−2k) = 0, η(0) = 1/2); overflows below s ≈ −30000.
[[nodiscard]] inline float256 dirichlet_eta(float256 s) {
    detail_dd::FnEnv env;
    using detail_dd::dd;
    const float256 one = float256::one();
    if (s.isnan()) return env.nan(s);
    if (s.isinf()) return s.signbit() ? env.invalid() : env.exact(one);
    if (s.iszero()) return env.exact(float256::half());
    if (s == one) return detail_dd::ln2().out(env);
    const dd sm1 = detail_dd::two_sum(s, -one);
    if (!s.signbit()) {
        if (s > float256(600)) return env.round(one, float256::zero(), 0, -1);  // 1 − 2^-s
        const dd f = -detail_dd::expm1(-(sm1 * detail_dd::ln2()));  // 1 − 2^{1−s}
        return (f * detail_sf::zeta::pos(dd(s), sm1)).out(env);
    }
    const float256 h = ldexp(s, -1);
    if (trunc(h) == h) return env.exact(float256::zero());
    dd sn, cn;
    detail_sf::sincos_pi(h, sn, cn);
    if (s < float256(-30000)) return env.overflow(!sn.hi.signbit());  // (1 − 2^{1−s}) < 0
    // ζ(s) = 2^s π^{s−1} sin(πs/2) Γ(1−s) ζ(1−s); 1 − 2^{1−s} = −2^{1−s}(1 − 2^{s−1})
    const dd oms = detail_dd::two_sum(one, -s);
    const dd z = detail_sf::zeta::pos(oms, dd(-s));
    const detail_sf::SD g = detail_sf::gamma_sd(oms);
    // 2^s · 2^{1−s} = 2: η = −2 π^{s−1} (1 − 2^{s−1}) sin(πs/2) Γ(1−s) ζ(1−s)
    const dd L = (dd(s) - dd(one)) * detail_dd::log(detail_dd::pi());
    const dd f = -detail_dd::expm1(sm1 * detail_dd::ln2());  // 1 − 2^{s−1}
    detail_sf::SD r = detail_sf::mul(detail_sf::exp_sd(L), g);
    r = detail_sf::mul(r, -(sn * z * f));
    ++r.k;
    return detail_sf::out(env, r);
}

/// Bernoulli number B_2n (B_0 = 1, B_2 = 1/6, ...) rounded to binary256, from
/// B_2n = (−1)^{n+1} 2 (2n)! ζ(2n)/(2π)^{2n}; |B_2n| overflows from n = 13704 on.
[[nodiscard]] inline float256 bernoulli_b2n(unsigned n) {
    detail_dd::FnEnv env;
    using detail_dd::dd;
    if (n == 0) return env.exact(float256::one());
    if (n <= static_cast<unsigned>(detail_dd::tables::kBernFactN)) {
        detail_sf::SD f{detail_sf::zeta::bern_fact(static_cast<int>(n)), 0};
        f = detail_sf::mul(f, detail_sf::gamma_sd(dd(float256::from_uint(2ull * n + 1))));
        return detail_sf::out(env, f);
    }
    if (n > 20000) return env.overflow(!(n & 1));
    const float256 N2 = float256::from_uint(2ull * n);
    // ζ(2n) = 1 + 2^{-2n} + … (n > 140: 1 + below 2^-280 rel.; include 2^-2n + 3^-2n)
    const dd z = dd(float256::one()) + dd(ldexp(float256::one(), -static_cast<int>(2 * n))) +
                 detail_dd::exp(-(dd(N2) * detail_dd::log(dd(float256(3)))));
    detail_sf::SD r = detail_sf::gamma_sd(dd(N2 + float256::one()));
    r = detail_sf::mul(r, detail_sf::exp_sd(-(dd(N2) * detail_dd::log(detail_dd::pi() * float256::two()))));
    r = detail_sf::mul(r, z);
    ++r.k;
    if (!(n & 1)) r.m = -r.m;
    return detail_sf::out(env, r);
}

} // namespace base256
