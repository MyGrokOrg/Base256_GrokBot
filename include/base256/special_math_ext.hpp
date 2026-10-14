#pragma once
/// More special functions (included by base256.hpp): the C++17 orthogonal
/// polynomials (legendre, assoc_legendre, sph_legendre, laguerre,
/// assoc_laguerre, hermite). Each evaluates in double-float256 (~474 bits)
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

} // namespace base256
