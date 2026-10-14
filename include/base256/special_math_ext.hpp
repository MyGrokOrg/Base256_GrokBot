#pragma once
/// More special functions (included by base256.hpp): the C++17 orthogonal
/// polynomials (legendre, assoc_legendre, sph_legendre, laguerre,
/// assoc_laguerre, hermite) and the exponential integrals expint (Ei) and
/// expint(n, x) (E_n). Each evaluates in double-float256 (~474 bits)
/// and rounds once in the caller's mode with IEEE flags. See README "Special
/// functions" for measured accuracy, domains and costs.

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

// S(x) = Σ_{k>=1} x^k/(k·k!), so Ei(x) = γ + ln|x| + S(x) and E1(t) = −γ − ln t − S(−t).
// Once the terms decrease (k > |x|) and fall below 2^-245 of the sum, the rest
// is summed in float256.
[[nodiscard]] inline dd s_series(float256 x) {
    dd term{float256::one()}, sum{float256::zero()};
    const float256 ax = abs(x);
    std::uint64_t k = 1;
    for (;; ++k) {
        const float256 K = float256::from_uint(k);
        term = term * x / K;
        const dd add = term / K;
        sum = sum + add;
        if (K > ax && f256_tail(add, sum)) break;
    }
    float256 tf = term.hi, st{};
    for (++k;; ++k) {
        const float256 K = float256::from_uint(k);
        tf = tf * x / K;
        const float256 add = tf / K;
        st = st + add;
        if (small_rel(dd(add), sum)) break;
    }
    return sum + dd(st);
}

inline constexpr int kEiAsym = 360;  // Ei asymptotic series from here (min term < 2^-500)

// Ei(x) for x >= kEiAsym: e^x/x Σ k!/x^k.
[[nodiscard]] inline SD ei_asym(float256 x) {
    dd term{float256::one()}, sum{float256::one()};
    std::uint64_t k = 1;
    for (;; ++k) {
        term = term * float256::from_uint(k) / x;
        sum = sum + term;
        if (f256_tail(term, sum)) break;
    }
    float256 tf = term.hi, st{};
    for (++k;; ++k) {
        tf = tf * float256::from_uint(k) / x;
        st = st + tf;
        if (small_rel(dd(tf), sum)) break;
    }
    sum = sum + dd(st);
    SD e = exp_sd(dd(x));
    return mul(e, sum / x);
}

// The continued fraction E_n(x) = e^{-x} / (b0 + a1/(b1 + a2/(b2 + ...))),
// b_k = x + n + 2k, a_k = −k(n − 1 + k) (division-free Wallis form, rescaled).
[[nodiscard]] inline dd en_cf(std::uint64_t n, float256 x) {
    const float256 nn = float256::from_uint(n);
    dd A0{float256::one()}, B0{float256::zero()};
    dd A1 = detail_dd::two_sum(x, nn), B1{float256::one()};
    dd h = A1, hprev = h;
    for (std::uint64_t k = 1;; ++k) {
        const float256 K = float256::from_uint(k);
        const dd b = detail_dd::two_sum(x, nn + K + K);
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
[[nodiscard]] inline dd en_series(std::uint64_t n, float256 x) {
    const dd one{float256::one()};
    dd sum{float256::zero()}, term = one;  // term = (−x)^k/k!
    const float256 nx = -x;
    std::uint64_t k = 0;
    for (;; ++k) {
        const float256 K = float256::from_uint(k);
        if (k > 0) term = term * nx / K;
        if (k + 1 != n) {
            const dd add = term / (K - float256::from_uint(n - 1));
            sum = sum - add;
            if (K > x && k + 1 > n && f256_tail(add, sum)) break;
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
        tf = tf * nx / K;
        const float256 add = tf / (K - float256::from_uint(n - 1));
        st = st - add;
        if (small_rel(dd(add), sum)) break;
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
[[nodiscard]] inline SD en_sd(std::uint64_t n, float256 x) {
    if (n == 0) return mul(exp_sd(-dd(x)), dd(float256::one()) / x);  // e^{-x}/x
    if (en_use_cf(n, x)) {
        const dd f = en_cf(n, x);
        return mul(exp_sd(-dd(x)), dd(float256::one()) / f);
    }
    return norm({en_series(n, x), 0});
}

// Ei(x), x != 0 finite, as m·2^k.
[[nodiscard]] inline SD ei_sd(float256 x) {
    if (x.signbit()) {
        SD e = en_sd(1, -x);
        e.m = -e.m;
        return e;
    }
    if (x >= float256(kEiAsym)) return ei_asym(x);
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

} // namespace base256
