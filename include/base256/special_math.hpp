#pragma once
/// Special functions (C++17 <cmath> special-math names): Bessel functions,
/// beta and the regularized incomplete beta, Riemann zeta, incomplete gamma,
/// complete and incomplete elliptic integrals (Carlson R_F, R_D, R_J; AGM).
/// Every function evaluates in double-float256 (~474 bits) and rounds once in
/// the caller's rounding mode, raising IEEE flags like the elementary
/// functions. See README "Special functions" for measured accuracy and limits.

#include <cmath>
#include <cstdint>
#include <vector>

#include "detail/dd.hpp"
#include "extra_math.hpp"
#include "float256.hpp"
#include "math.hpp"

namespace base256 {

namespace detail_sf {

using detail_dd::dd;

[[nodiscard]] inline dd D(int i) { return dd(float256::from_int(i)); }
[[nodiscard]] inline dd dabs(const dd& a) { return a.hi.signbit() ? -a : a; }
[[nodiscard]] inline dd euler_gamma() { return detail_dd::dd_from_words(detail_dd::tables::kEulerGamma); }
// Relative convergence threshold for the series and continued fractions.
inline const float256& eps() { static const float256 e = ldexp(float256::one(), -490); return e; }
[[nodiscard]] inline bool small_rel(const dd& t, const dd& s) {
    if (t.hi.iszero()) return true;
    if (s.hi.iszero()) return false;
    return ilogb(t.hi) < ilogb(s.hi) - 490;
}

// Continued-fraction convergence: a factor within ~2^-468 of 1 (the dd
// working precision is ~2^-474, so tighter tests may never be met).
[[nodiscard]] inline bool near_one(const dd& del) {
    const dd e = del - dd(float256::one());
    return e.hi.iszero() || ilogb(e.hi) < -468;
}
// Term t below 2^-245 of the sum s: the rest of a decreasing series only needs float256.
[[nodiscard]] inline bool f256_tail(const dd& t, const dd& s) {
    return t.hi.iszero() || (!s.hi.iszero() && ilogb(t.hi) < ilogb(s.hi) - 245);
}
[[nodiscard]] inline bool near_zero(const dd& v) { return v.hi.iszero() || ilogb(v.hi) < -468; }

// A value m·2^k with a dd significand; keeps huge and tiny intermediates in range.
struct SD {
    dd m;
    int k = 0;
};
[[nodiscard]] inline SD norm(SD a) {
    if (a.m.hi.iszero() || !a.m.hi.isfinite()) return a;
    const int e = ilogb(a.m.hi);
    return {detail_dd::ldexp(a.m, -e), a.k + e};
}
[[nodiscard]] inline SD mul(const SD& a, const SD& b) { return norm({a.m * b.m, a.k + b.k}); }
[[nodiscard]] inline SD mul(const SD& a, const dd& b) { return norm({a.m * b, a.k}); }
// exp(x) as m·2^k; |x| beyond the binary256 range is clamped (the result
// then overflows or underflows when rounded).
[[nodiscard]] inline SD exp_sd(const dd& x) {
    if (x.hi > float256(400000)) return {dd(float256::one()), 600000};
    if (x.hi < float256(-400000)) return {dd(float256::one()), -600000};
    const detail_dd::Scaled s = detail_dd::exp_scaled(x);
    return norm({s.m, s.k});
}
[[nodiscard]] inline float256 out(detail_dd::FnEnv& env, const SD& v, int tail = 2) {
    if (v.k > 300000) return env.overflow(v.m.hi.signbit());
    if (v.k < -300000) return env.underflow(v.m.hi.signbit());
    return env.round(v.m.hi, v.m.lo, v.k, tail);
}
// 1/Γ(1+μ) = E(μ²) + μ·O(μ²) for |μ| <= 1/2 from the MPFR Taylor table
// (detail/tables.hpp); trailing terms below 2^-492 are dropped.
inline void rgamma1p_eo(const dd& mu, dd& E, dd& O) {
    namespace T = detail_dd::tables;
    E = O = dd(float256::zero());
    if (mu.hi.iszero()) { E = dd(float256::one()); return; }
    const int lm = ilogb(mu.hi) + 1;  // |μ| < 2^lm
    int n = T::kRGamma1pN - 1;
    while (n > 2) {
        const float256 c = detail_dd::dd_from_words(T::kRGamma1p[n]).hi;
        if (c.iszero() || ilogb(c) + n * lm >= -492) break;
        --n;
    }
    // Terms below 2^-245 of the sum (|a_k μ^k| < 2^-245 and a_0 = 1) only need
    // float256 precision: run that part of Horner in float256.
    int nf = n | 1;
    while (nf > 1) {
        const float256 c = detail_dd::dd_from_words(T::kRGamma1p[nf - 1]).hi;
        if (!c.iszero() && ilogb(c) + (nf - 1) * lm >= -245) break;
        nf -= 2;
    }
    const dd m2 = mu * mu;
    float256 Of{}, Ef{};
    for (int k = n | 1; k > nf; k -= 2) {
        Of = Of * m2.hi + detail_dd::dd_from_words(T::kRGamma1p[k]).hi;
        Ef = Ef * m2.hi + detail_dd::dd_from_words(T::kRGamma1p[k - 1]).hi;
    }
    O = dd(Of); E = dd(Ef);
    for (int k = nf; k >= 1; k -= 2) {
        O = O * m2 + detail_dd::dd_from_words(T::kRGamma1p[k]);
        E = E * m2 + detail_dd::dd_from_words(T::kRGamma1p[k - 1]);
    }
}
inline constexpr int kGammaTableMax = 40;  // table + product for z < this
// Γ(z) as m·2^k for a dd z > 0.
[[nodiscard]] inline SD gamma_sd(const dd& z) {
    if (z.hi >= float256::half() && z.hi < float256(kGammaTableMax)) {
        // z = 1 + μ + n, |μ| <= 1/2: Γ(z) = (1+μ)…(n+μ) / (1/Γ(1+μ)).
        const float256 r = floor(z.hi + float256::half());
        const int n = static_cast<int>(r.to_double()) - 1;
        const dd mu = z - dd(r);
        if (n >= 0 && abs(mu.hi) <= float256::half()) {
            dd E, O;
            rgamma1p_eo(mu, E, O);
            dd p{float256::one()};
            for (int j = 1; j <= n; ++j) p = p * (dd(float256::from_int(j)) + mu);
            return norm({p / (E + O * mu), 0});
        }
    }
    const detail_dd::Scaled s = detail_gamma::gamma_abs_scaled_pos(z);
    return norm({s.m, s.k});
}
// 1/Γ(z) for a dd z > 0 as m·2^k.
[[nodiscard]] inline SD rgamma_sd(const dd& z) {
    const SD g = gamma_sd(z);
    return norm({dd(float256::one()) / g.m, -g.k});
}
// sin(π t) and cos(π t) for any finite float256 t (exact reduction mod 2).
inline void sincos_pi(float256 t, dd& s, dd& c) {
    const bool neg = t.signbit();
    const float256 r = fmod(abs(t), float256::two());  // exact, r in [0, 2)
    // quadrant q = floor(2r), f = r − q/2 in [0, 1/2) (exact)
    const float256 q = floor(r + r);
    const float256 f = r - q * float256::half();
    dd sf, cf;
    detail_dd::sincos_reduced(detail_dd::pi() * f, sf, cf);
    switch (static_cast<int>(q.to_double()) & 3) {
        case 0: s = sf; c = cf; break;
        case 1: s = cf; c = -sf; break;
        case 2: s = -sf; c = -cf; break;
        default: s = -cf; c = sf; break;
    }
    if (neg) s = -s;
}
[[nodiscard]] inline dd sinpi(float256 t) { dd s, c; sincos_pi(t, s, c); return s; }

// sin x, cos x for a finite float256 x (any size).
inline void sincos_x(float256 x, dd& s, dd& c) {
    int quad = 0;
    dd sd, cd;
    detail_math::trig_dd(abs(x), sd, cd, quad);
    s = (quad & 1) ? cd : sd;
    c = (quad & 1) ? -sd : cd;
    if (quad >= 2) { s = -s; c = -c; }
    if (x.signbit()) s = -s;
}

// cosh(e) and sinh(e)/e for a dd e (|e| < 4e5).
inline void cosh_sinhc(const dd& e, dd& ch, dd& shc) {
    const dd one{float256::one()};
    if (abs(e.hi) < float256(0.5)) {
        // Σ e^{2k}/(2k)!, Σ e^{2k}/(2k+1)!
        // Terms with e^{2k}/(2k)! < 2^-492 are dropped; those below 2^-245 run in float256.
        const auto& f = detail_dd::inv_factorials();
        const dd e2 = e * e;
        if (e2.hi.iszero()) { ch = shc = one; return; }
        const int le = ilogb(e2.hi) + 1;  // e² < 2^le
        int n = 50, nf = 50;
        while (n > 1 && ilogb(f[static_cast<std::size_t>(2 * n)].hi) + n * le < -492) --n;
        nf = n;
        while (nf > 0 && ilogb(f[static_cast<std::size_t>(2 * nf)].hi) + nf * le < -245) --nf;
        float256 af{}, bf{};
        for (int k = n; k > nf; --k) {
            af = af * e2.hi + f[static_cast<std::size_t>(2 * k)].hi;
            bf = bf * e2.hi + f[static_cast<std::size_t>(2 * k + 1)].hi;
        }
        dd a{af}, b{bf};
        for (int k = nf; k >= 0; --k) {
            a = a * e2 + f[static_cast<std::size_t>(2 * k)];
            b = b * e2 + f[static_cast<std::size_t>(2 * k + 1)];
        }
        ch = a; shc = b;
        return;
    }
    const dd E = detail_dd::exp(e);
    const dd iE = one / E;
    ch = detail_dd::ldexp(E + iE, -1);
    shc = (E - iE) / (e + e);
}

namespace bessel {

// Hankel asymptotic sums at order m: t_k = a_k(m)/x^k, a_k = Π_{j<=k} (4m² − (2j−1)²) / (k! 8^k).
// kind 0 (J/Y): P = Σ (−1)^j t_{2j}, Q = Σ (−1)^j t_{2j+1};
// kind 1 (I):   S = Σ (−1)^k t_k;  kind 2 (K): S = Σ t_k.
// Returns false if the terms reach their minimum above 2^-330 relative.
// (The expansion is exact for half-integer m; Hankel's P, Q are used from x >= 120,
// where the smallest term is e^{-2x} < 2^-346.)
[[nodiscard]] inline bool asym(float256 m, float256 x, int kind, dd& P, dd& Q) {
    const dd m4 = detail_dd::two_prod(m, m) * float256(4);
    dd t{float256::one()};
    P = t; Q = dd(float256::zero());
    float256 prev = float256::inf(false);
    const dd x8 = dd(x) * float256(8);
    for (int k = 1; k < 100000; ++k) {
        const std::int64_t o = 2 * static_cast<std::int64_t>(k) - 1;
        const dd num = m4 - dd(float256::from_int(o * o));
        if (num.hi.iszero()) return true;  // terminating series (half-integer order)
        t = t * num / (x8 * float256::from_int(k));
        const float256 at = abs(t.hi);
        const dd& ref = P;
        // Asymptotic: stop at the smallest term; accept if it is below 2^-330 relative.
        if (at > prev) return ilogb(prev) < ilogb(ref.hi) - 330;
        prev = at;
        if (kind == 0) {
            const bool odd = k & 1;
            const int j = k / 2;
            const dd term = (j & 1) ? -t : t;
            if (odd) Q = Q + term; else P = P + term;
        } else if (kind == 1) {
            P = P + ((k & 1) ? -t : t);
        } else {
            P = P + t;
        }
        if (ilogb(at) < ilogb(ref.hi) - 495) return true;
        if (ilogb(at) < -495 && kind == 0) return true;
    }
    return false;
}

// Hankel: J_m(x), Y_m(x) for large x. Returns false if not converged.
[[nodiscard]] inline bool hankel_jy(float256 m, float256 x, dd& J, dd& Y) {
    dd P, Q;
    if (!asym(m, x, 0, P, Q)) return false;
    // χ = x − (m/2 + 1/4)π = x − φ, φ = (m + 1/2)·π/2
    const dd t = detail_dd::two_sum(m, float256::half());
    const float256 q = floor(t.hi);
    dd f = t - dd(q);
    dd sp, cp;
    detail_dd::sincos_reduced(detail_dd::pio2() * f, sp, cp);
    const int qi = static_cast<int>(fmod(q, float256(4)).to_double());
    for (int i = 0; i < qi; ++i) { const dd tmp = sp; sp = cp; cp = -tmp; }  // φ += π/2
    dd sx, cx;
    sincos_x(x, sx, cx);
    const dd cchi = cx * cp + sx * sp;
    const dd schi = sx * cp - cx * sp;
    const dd amp = detail_dd::sqrt(dd(float256::two()) / (detail_dd::pi() * x));
    J = amp * (P * cchi - Q * schi);
    Y = amp * (P * schi + Q * cchi);
    return true;
}

// Power series: (x/2)^ν/Γ(ν+1) · Σ (s x²/4)^k/(k!(ν+1)_k), s = −1 (J) or +1 (I); ν > −1.
[[nodiscard]] inline SD series(const dd& nu, float256 x, int s) {
    const dd one{float256::one()};
    const dd np1 = nu + one;
    SD pre;
    if (nu.hi.iszero()) pre = {one, 0};
    else pre = mul(exp_sd(nu * detail_dd::log(dd(ldexp(x, -1)))), rgamma_sd(np1));
    const dd q = detail_dd::two_prod(ldexp(x, -1), ldexp(x, -1));
    const dd qs = s < 0 ? -q : q;
    dd term = one, sum = one;
    int k = 1;
    for (; k < 2000000; ++k) {
        term = term * qs / (dd(float256::from_int(k)) * (np1 + D(k - 1)));
        sum = sum + term;
        if (float256::from_int(k) * float256::from_int(k) > abs(q.hi)) {
            if (small_rel(term, sum)) return mul(pre, sum);
            // Past the peak, terms below 2^-245 of the sum only need float256.
            if (term.hi.iszero() || ilogb(term.hi) < ilogb(sum.hi) - 245) break;
        }
    }
    float256 t = term.hi, tail{};
    const float256 np1f = np1.hi;
    const int lim = ilogb(sum.hi) - 492;
    for (++k; k < 2000000; ++k) {
        t = t * qs.hi / (float256::from_int(k) * (np1f + float256::from_int(k - 1)));
        tail = tail + t;
        if (t.iszero() || ilogb(t) < lim) break;
    }
    return mul(pre, sum + dd(tail));
}

// Forward recurrence v_{j+1} = (2(μ+j)/x) v_j + s v_{j−1}, n steps from (v0, v1) = (v_μ, v_{μ+1}).
// Returns v_{μ+n} as m·2^k (rescaling as needed).
[[nodiscard]] inline SD recur(dd v0, dd v1, float256 mu, float256 x, int n, int s) {
    if (n == 0) return norm({v0, 0});
    int k = 0;
    const dd xi2 = dd(float256::two()) / dd(x);
    for (int j = 1; j < n; ++j) {
        const dd v2 = (dd(mu) + D(j)) * xi2 * v1 + (s < 0 ? -v0 : v0);
        v0 = v1; v1 = v2;
        const int e = ilogb(v1.hi);
        if (e > 4000 || e < -4000) { v0 = detail_dd::ldexp(v0, -e); v1 = detail_dd::ldexp(v1, -e); k += e; }
    }
    return norm({v1, k});
}

// 1/Γ(1+μ) = E(μ²) + μ O(μ²) from the MPFR Taylor table (|μ| <= 1/2), so
// gam1 = (1/Γ(1−μ) − 1/Γ(1+μ))/(2μ) = −O and gam2 = (1/Γ(1−μ) + 1/Γ(1+μ))/2 = E.
inline void temme_gammas(float256 mu, dd& gam1, dd& gam2, dd& gampl, dd& gammi) {
    if (mu.iszero()) { gam1 = -euler_gamma(); gam2 = gampl = gammi = dd(float256::one()); return; }
    dd E, O;
    rgamma1p_eo(dd(mu), E, O);
    const dd mo = O * mu;
    gampl = E + mo;
    gammi = E - mo;
    gam2 = E;
    gam1 = -O;
}

// Shared setup for Temme's series at |μ| <= 1/2, 0 < x: fact·(gam1 cosh e + gam2 d sinh(e)/e)
// with d = −ln(x/2), e = μd, and E = e^e; also sin(πμ/2)/(πμ/2).
struct TemmeSetup {
    dd gam1, gam2, gampl, gammi, ff0, E, fact3;
};
inline TemmeSetup temme_setup(float256 mu, float256 x) {
    const dd one{float256::one()};
    const dd pi = detail_dd::pi();
    TemmeSetup t;
    temme_gammas(mu, t.gam1, t.gam2, t.gampl, t.gammi);
    const dd d = -detail_dd::log(dd(ldexp(x, -1)));
    const dd e = d * mu;
    dd ch, shc;
    cosh_sinhc(e, ch, shc);
    t.E = abs(e.hi) < float256(0.5) ? ch + e * shc : detail_dd::exp(e);
    dd fact = one;
    t.fact3 = one;
    if (!mu.iszero()) {
        // sin(π|μ|) = 2 sin(π|μ|/2) cos(π|μ|/2): one reduced sincos for both factors.
        const float256 hm = ldexp(abs(mu), -1);
        const dd a = pi * hm;
        dd sh, chh;
        detail_dd::sincos_reduced(a, sh, chh);
        fact = a / (sh * chh);   // π|μ| / sin(π|μ|)
        t.fact3 = sh / a;        // sin(π|μ|/2)/(π|μ|/2)
    }
    t.ff0 = fact * (t.gam1 * ch + t.gam2 * shc * d);
    return t;
}

// Y_μ(x), Y_{μ+1}(x) for 0 < x, |μ| <= 1/2 (Temme's series; cancellation ~e^x).
inline void temme_y(float256 mu, float256 x, dd& ymu, dd& y1) {
    const dd one{float256::one()};
    const dd pi = detail_dd::pi();
    const TemmeSetup t = temme_setup(mu, x);
    const float256 x2 = ldexp(x, -1);
    dd ff = t.ff0 * float256::two() / pi;
    dd p = t.E / (t.gampl * pi);
    dd q = one / (t.E * pi * t.gammi);
    const dd r_s = pi * pi * ldexp(mu, -1) * t.fact3 * t.fact3;  // π·(πμ/2)·(sin(πμ/2)/(πμ/2))²
    dd c = one;
    const dd dq = -detail_dd::two_prod(x2, x2);
    dd sum = ff + r_s * q, sum1 = p;
    const dd mu2 = detail_dd::two_prod(mu, mu);
    const dd M{mu};
    for (int i = 1; i < 100000; ++i) {
        const dd di = D(i);
        const dd r = one / (di * di - mu2);  // 1/(i−μ) = (i+μ)r, 1/(i+μ) = (i−μ)r
        ff = (di * ff + p + q) * r;
        c = c * dq / float256::from_int(i);
        p = p * (di + M) * r;
        q = q * (di - M) * r;
        const dd del = c * (ff + r_s * q);
        sum = sum + del;
        const dd del1 = c * p - di * del;
        sum1 = sum1 + del1;
        if (small_rel(del, sum + one) && small_rel(del1, sum1)) break;
        if (float256::from_int(i) * float256::from_int(i) > abs(dq.hi) && f256_tail(del, sum + one) && f256_tail(c * p, sum1) && f256_tail(di * del, sum1)) {
            // Remaining terms are below 2^-245 of the sums: finish in float256.
            float256 F = ff.hi, P = p.hi, Q = q.hi, C = c.hi, ts{}, ts1{};
            const float256 m = mu, rs = r_s.hi, dqf = dq.hi;
            const int l0 = ilogb((sum + one).hi) - 492, l1 = ilogb(sum1.hi) - 492;
            for (int j = i + 1; j < 100000; ++j) {
                const float256 dj = float256::from_int(j);
                const float256 rinv = float256::one() / (dj * dj - m * m);
                F = (dj * F + P + Q) * rinv;
                C = C * dqf / dj;
                P = P * (dj + m) * rinv;
                Q = Q * (dj - m) * rinv;
                const float256 d0 = C * (F + rs * Q), d1 = C * P - dj * d0;
                ts = ts + d0; ts1 = ts1 + d1;
                if ((d0.iszero() || ilogb(d0) < l0) && (d1.iszero() || ilogb(d1) < l1)) break;
            }
            sum = sum + dd(ts); sum1 = sum1 + dd(ts1);
            break;
        }
    }
    ymu = -sum;
    y1 = -sum1 * float256::two() / dd(x);
}

// K_μ(x), K_{μ+1}(x) for 0 < x, |μ| <= 1/2 (Temme's series; cancellation ~e^{2x}).
inline void temme_k(float256 mu, float256 x, dd& kmu, dd& k1) {
    const dd one{float256::one()};
    const TemmeSetup t = temme_setup(mu, x);
    const float256 x2 = ldexp(x, -1);
    dd ff = t.ff0;
    dd sum = ff;
    dd p = detail_dd::ldexp(t.E / t.gampl, -1);
    dd q = detail_dd::ldexp(one / (t.E * t.gammi), -1);
    dd c = one;
    const dd dq = detail_dd::two_prod(x2, x2);
    dd sum1 = p;
    const dd mu2 = detail_dd::two_prod(mu, mu);
    const dd M{mu};
    for (int i = 1; i < 100000; ++i) {
        const dd di = D(i);
        const dd r = one / (di * di - mu2);
        ff = (di * ff + p + q) * r;
        c = c * dq / float256::from_int(i);
        p = p * (di + M) * r;
        q = q * (di - M) * r;
        const dd del = c * ff;
        sum = sum + del;
        const dd del1 = c * (p - di * ff);
        sum1 = sum1 + del1;
        if (small_rel(del, sum) && small_rel(del1, sum1)) break;
        if (float256::from_int(i) * float256::from_int(i) > abs(dq.hi) && f256_tail(del, sum) && f256_tail(c * p, sum1) && f256_tail(c * di * ff, sum1)) {
            float256 F = ff.hi, P = p.hi, Q = q.hi, C = c.hi, ts{}, ts1{};
            const float256 m = mu, dqf = dq.hi;
            const int l0 = ilogb(sum.hi) - 492, l1 = ilogb(sum1.hi) - 492;
            for (int j = i + 1; j < 100000; ++j) {
                const float256 dj = float256::from_int(j);
                const float256 rinv = float256::one() / (dj * dj - m * m);
                F = (dj * F + P + Q) * rinv;
                C = C * dqf / dj;
                P = P * (dj + m) * rinv;
                Q = Q * (dj - m) * rinv;
                const float256 d0 = C * F, d1 = C * (P - dj * F);
                ts = ts + d0; ts1 = ts1 + d1;
                if ((d0.iszero() || ilogb(d0) < l0) && (d1.iszero() || ilogb(d1) < l1)) break;
            }
            sum = sum + dd(ts); sum1 = sum1 + dd(ts1);
            break;
        }
    }
    kmu = sum;
    k1 = sum1 * float256::two() / dd(x);
}

// K_μ, K_{μ+1} for x >= 2, |μ| <= 1/2 (Steed's continued fraction CF2), as
// m·2^k with a shared exponent (e^{−x} factor).
inline void steed_k(float256 mu, float256 x, dd& kmu, dd& k1, int& kexp) {
    const dd one{float256::one()};
    const dd X{x};
    const dd mu2 = detail_dd::two_prod(mu, mu);
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
    kmu = pre;
    k1 = pre * (dd(mu) + X + dd(float256::half()) - h) / X;
    kexp = ex.k;
}

// K_μ, K_{μ+1} for moderate x (|μ| <= 1/2) without Steed's slowly converging
// normalization sum: from the Wronskian I_μ K_{μ+1} + I_{μ+1} K_μ = 1/x,
//   K_μ = 1 / (x I_μ (r + ρ)),  r = K_{μ+1}/K_μ,  ρ = I_{μ+1}/I_μ,
// with I_μ from its (cancellation-free) power series and r, ρ from continued
// fractions evaluated by division-free Wallis recurrences:
//   r = (μ + x + 1/2 − a_1 T)/x, T = 1/(b_1 + a_2/(b_2 + …)), a_n = −(a_1 + n(n−1)),
//   a_1 = 1/4 − μ², b_n = 2(x + n);   ρ = 1/(c_1 + 1/(c_2 + …)), c_n = 2(μ+n)/x.
// Returns K_μ, K_{μ+1} as m·2^kexp (shared exponent).
[[nodiscard]] inline dd wallis_ratio(int kind, const dd& mu, float256 x) {
    // kind 0: T (K ratio), kind 1: ρ·x (I ratio, scaled: 1/(c'_1 + x²/(c'_2 + …)), c'_n = 2(μ+n)).
    const dd one{float256::one()};
    const dd a1 = dd(float256(0.25)) - mu * mu;
    const dd x2 = detail_dd::two_prod(x, x);
    dd A0 = one, A1{float256::zero()}, B0{float256::zero()}, B1 = one;
    dd r{float256::zero()}, rprev{float256::one()};
    for (int n = 1; n < 10000000; ++n) {
        dd an, bn;
        if (kind == 0) { an = n == 1 ? one : -(a1 + D(n * (n - 1))); bn = detail_dd::two_sum(float256::from_int(2 * n), x + x); }
        else { an = n == 1 ? one : x2; bn = (mu + D(n)) * float256::two(); }
        const dd A2 = bn * A1 + an * A0, B2 = bn * B1 + an * B0;
        A0 = A1; A1 = A2; B0 = B1; B1 = B2;
        const int e = ilogb(B1.hi);
        if (e > 2000 || e < -2000) {
            A0 = detail_dd::ldexp(A0, -e); A1 = detail_dd::ldexp(A1, -e);
            B0 = detail_dd::ldexp(B0, -e); B1 = detail_dd::ldexp(B1, -e);
        }
        if ((n & 7) == 0) {
            r = A1 / B1;
            if (near_one(r / rprev)) break;
            rprev = r;
        }
    }
    return r;
}
// I_μ(x) = pre·S0 and ρ = I_{μ+1}/I_μ = (x/2) S1/S0 from one series loop:
// t_k = q^k/(k!(μ+1)_k), q = x²/4, S0 = Σ t_k, S1 = Σ u_k, u_k = t_k/(μ+1+k),
// t_{k+1} = u_k q/(k+1): one dd division per term; terms below 2^-245 of the
// sums finish in float256.
inline void i_pair(float256 mu, float256 x, SD& I, dd& rho) {
    const dd one{float256::one()};
    const dd M{mu};
    const dd q = detail_dd::two_prod(ldexp(x, -1), ldexp(x, -1));
    dd t = one, S0 = one, S1{float256::zero()};
    int k = 0;
    for (; k < 2000000; ++k) {
        const dd u = t / (M + D(k + 1));
        S1 = S1 + u;
        t = u * q / float256::from_int(k + 1);
        S0 = S0 + t;
        if (float256::from_int(k) * float256::from_int(k) > abs(q.hi) && f256_tail(t, S0) && f256_tail(u, S1)) break;
    }
    float256 tf = t.hi, s0{}, s1{};
    const int l0 = ilogb(S0.hi) - 492, l1 = ilogb(S1.hi) - 492;
    for (++k; k < 2000000; ++k) {
        const float256 u = tf / (mu + float256::from_int(k + 1));
        s1 = s1 + u;
        tf = u * q.hi / float256::from_int(k + 1);
        s0 = s0 + tf;
        if ((tf.iszero() || ilogb(tf) < l0) && (u.iszero() || ilogb(u) < l1)) break;
    }
    S0 = S0 + dd(s0);
    S1 = S1 + dd(s1);
    const SD pre = mul(exp_sd(M * detail_dd::log(dd(ldexp(x, -1)))), rgamma_sd(M + one));
    I = mul(pre, S0);
    rho = S1 * ldexp(x, -1) / S0;
}
inline void k_wronskian(float256 mu, float256 x, dd& kmu, dd& k1, int& kexp) {
    const dd one{float256::one()};
    const dd M{mu}, X{x};
    const dd a1 = dd(float256(0.25)) - detail_dd::two_prod(mu, mu);
    const dd T = wallis_ratio(0, M, x);
    const dd r = (M + X + dd(float256::half()) - a1 * T) / X;
    SD I;
    dd rh;
    i_pair(mu, x, I, rh);
    kmu = one / (X * I.m * (r + rh));
    k1 = kmu * r;
    kexp = -I.k;
}

// Y_μ, Y_{μ+1} for x >= 2, |μ| <= 1/2 (Steed: CF1 for J'/J, CF2 for (J'+iY')/(J+iY)).
inline void steed_y(float256 mu, float256 x, dd& ymu, dd& y1) {
    const dd one{float256::one()};
    const dd X{x};
    const dd xi = one / X;
    const dd xi2 = xi * float256::two();
    // CF1 for f = J'_μ/J_μ = μ/x − J_{μ+1}/J_μ, J_{μ+1}/J_μ = 1/(b_1 − 1/(b_2 − …)),
    // b_n = 2(μ+n)/x: Wallis recurrence (no divisions), power-of-two rescaling,
    // tested every 8 steps. isign = sign of B_n (the unnormalized J_μ of the
    // downward recurrence), as the product of Lentz's D_n = B_{n−1}/B_n.
    int isign = 1;
    dd h;
    {
        // A_{−1} = 1, A_0 = 0, B_{−1} = 0, B_0 = 1; the ratio is A_n/B_n.
        dd A0 = one, A1{float256::zero()}, B0{float256::zero()}, B1 = one;
        dd bn = xi2 * mu;
        dd r{float256::zero()}, rprev{float256::one()};
        for (int n = 1; n < 10000000; ++n) {
            bn = bn + xi2;
            const dd A2 = n == 1 ? A0 : bn * A1 - A0, B2 = n == 1 ? bn : bn * B1 - B0;  // a_1 = +1, a_n = −1
            A0 = A1; A1 = A2; B0 = B1; B1 = B2;
            const int e = ilogb(B1.hi);
            if (e > 2000 || e < -2000) {
                A0 = detail_dd::ldexp(A0, -e); A1 = detail_dd::ldexp(A1, -e);
                B0 = detail_dd::ldexp(B0, -e); B1 = detail_dd::ldexp(B1, -e);
            }
            if ((n & 7) == 0) {
                r = A1 / B1;
                if (near_one(r / rprev)) break;
                rprev = r;
            }
        }
        if (B1.hi.signbit()) isign = -1;
        h = dd(mu) * xi - r;
    }
    const dd f = h;
    // CF2: p + iq = (J'+iY')/(J+iY) = −1/(2x) + i + (i/x) T, T = a_1/(b_1 + a_2/(b_2 + …)),
    // a_k = ((2k−1)/2)² − μ², b_k = 2x + 2k i. Complex Wallis recurrence (no
    // divisions; b_k has float256 parts), rescaled by powers of two, tested every 8 steps.
    const dd mu2 = detail_dd::two_prod(mu, mu);
    const float256 x2 = x + x;
    dd Ar0 = one, Ai0{float256::zero()}, Ar1{float256::zero()}, Ai1{float256::zero()};  // A_{k−2}, A_{k−1}
    dd Br0{float256::zero()}, Bi0{float256::zero()}, Br1 = one, Bi1{float256::zero()};
    dd Tr{float256::zero()}, Ti{float256::zero()}, Pr{float256::zero()}, Pi{float256::zero()};
    for (int k = 1; k < 1000000; ++k) {
        const dd ak = dd(ldexp(float256::from_int((2 * k - 1) * (2 * k - 1)), -2)) - mu2;
        const float256 bk = float256::from_int(2 * k);
        const dd Ar2 = Ar1 * x2 - Ai1 * bk + ak * Ar0, Ai2 = Ai1 * x2 + Ar1 * bk + ak * Ai0;
        const dd Br2 = Br1 * x2 - Bi1 * bk + ak * Br0, Bi2 = Bi1 * x2 + Br1 * bk + ak * Bi0;
        Ar0 = Ar1; Ai0 = Ai1; Ar1 = Ar2; Ai1 = Ai2;
        Br0 = Br1; Bi0 = Bi1; Br1 = Br2; Bi1 = Bi2;
        const int e = ilogb(abs(Br1.hi) > abs(Bi1.hi) ? Br1.hi : Bi1.hi);
        if (e > 2000 || e < -2000) {
            for (dd* v : {&Ar0, &Ai0, &Ar1, &Ai1, &Br0, &Bi0, &Br1, &Bi1}) *v = detail_dd::ldexp(*v, -e);
        }
        if ((k & 7) == 0) {
            const dd iden = one / (Br1 * Br1 + Bi1 * Bi1);
            Tr = (Ar1 * Br1 + Ai1 * Bi1) * iden;
            Ti = (Ai1 * Br1 - Ar1 * Bi1) * iden;
            const dd dr = Tr - Pr, di = Ti - Pi;
            const int lt = ilogb(abs(Tr.hi) > abs(Ti.hi) ? Tr.hi : Ti.hi);
            if ((dr.hi.iszero() || ilogb(dr.hi) < lt - 468) && (di.hi.iszero() || ilogb(di.hi) < lt - 468)) break;
            Pr = Tr; Pi = Ti;
        }
    }
    const dd p = -detail_dd::ldexp(xi, -1) - Ti * xi;
    const dd q = one + Tr * xi;
    const dd w = xi2 / detail_dd::pi();
    const dd gam = (p - f) / q;
    dd rjmu = detail_dd::sqrt(w / ((p - f) * gam + q));
    if (isign < 0) rjmu = -rjmu;
    ymu = rjmu * gam;
    const dd ymup = ymu * (p + q / gam);
    y1 = dd(mu) * xi * ymu - ymup;
}

// CF1 (modified Lentz) for f = J'_ν/J_ν at x > 0; isign = sign of the
// unnormalized J_ν that the downward recurrence starts from.
[[nodiscard]] inline dd cf1_j(float256 nu, float256 x, int& isign) {
    const dd one{float256::one()};
    const dd xi = one / dd(x);
    const dd xi2 = xi * float256::two();
    const float256 tiny = ldexp(float256::one(), -20000);
    isign = 1;
    dd h = dd(nu) * xi;
    if (abs(h.hi) < tiny) h = dd(tiny);
    dd b = xi2 * nu;
    dd d{float256::zero()};
    dd c = h;
    for (int i = 1; i < 100000000; ++i) {
        b = b + xi2;
        d = b - d;
        if (abs(d.hi) < tiny) d = dd(tiny);
        c = b - one / c;
        if (abs(c.hi) < tiny) c = dd(tiny);
        d = one / d;
        const dd del = c * d;
        h = del * h;
        if (d.hi.signbit()) isign = -isign;
        if (near_one(del)) break;
    }
    return h;
}

// J_ν(x) = J_μ(x) · (J_ν/J_μ), the ratio from CF1 and the (stable) downward
// recurrence ν → μ; J_μ, J_μ+1 from Hankel (x >= 120).
[[nodiscard]] inline bool j_down(float256 nu, float256 x, int n, float256 mu, SD& out_v) {
    dd Jm, Ym;
    if (!hankel_jy(mu, x, Jm, Ym)) return false;
    int isign = 1;
    const dd f = cf1_j(nu, x, isign);
    const dd xi = dd(float256::one()) / dd(x);
    dd rjl{float256::from_int(isign)};
    dd rjpl = f * rjl;
    const dd rjl1 = rjl;
    int k = 0;
    dd fact = dd(nu) * xi;
    for (int l = n; l >= 1; --l) {
        const dd rjtemp = fact * rjl + rjpl;
        fact = fact - xi;
        rjpl = fact * rjtemp - rjl;
        rjl = rjtemp;
        const int e = ilogb(rjl.hi);
        if (e > 4000) { rjl = detail_dd::ldexp(rjl, -e); rjpl = detail_dd::ldexp(rjpl, -e); k += e; }
    }
    out_v = norm({Jm * rjl1 / rjl, -k});
    return true;
}

// Splits ν >= 0 into n + μ, n = round(ν), μ in [−1/2, 1/2).
inline void split(float256 nu, int& n, float256& mu) {
    float256 r = floor(nu + float256::half());
    if (r - nu > float256::half()) r = r - float256::one();
    n = static_cast<int>(r.to_double());
    mu = nu - r;  // exact
    if (mu >= float256::half()) { mu = mu - float256::one(); ++n; }
}

inline constexpr int kAsymX = 120;  // Hankel / large-x expansions from here
// Temme's series below these x (it cancels ~e^x for Y and ~e^{2x} for K). Above:
// Y by Steed (division-free CF1/CF2), K by the I-series Wronskian (no cancellation),
// whichever was measured fastest; Steed's K remains the fallback from kAsymX.
inline constexpr int kTemmeYX = 56;
inline constexpr int kTemmeKX = 30;  // then the I-series Wronskian (k_wronskian) up to kAsymX

// J_ν(x), ν >= 0, x > 0 finite.
[[nodiscard]] inline SD j(float256 nu, float256 x) {
    if (x < float256(kAsymX)) return series(dd(nu), x, -1);
    dd J, Y;
    if (hankel_jy(nu, x, J, Y)) return norm({J, 0});
    int n; float256 mu;
    split(nu, n, mu);
    if (nu >= x - float256(1)) {
        SD v;
        if (j_down(nu, x, n, mu, v)) return v;
        return series(dd(nu), x, -1);
    }
    dd J0, Y0, J1, Y1;
    const bool ok0 = hankel_jy(mu, x, J0, Y0);
    const bool ok1 = hankel_jy(mu + float256::one(), x, J1, Y1);
    if (!ok0 || !ok1) return series(dd(nu), x, -1);
    return recur(J0, J1, mu, x, n, -1);
}

// Y_ν(x), ν >= 0, x > 0 finite.
[[nodiscard]] inline SD y(float256 nu, float256 x) {
    dd J, Y;
    if (x >= float256(kAsymX) && hankel_jy(nu, x, J, Y)) return norm({Y, 0});
    int n; float256 mu;
    split(nu, n, mu);
    dd y0, y1;
    if (x < float256(kTemmeYX)) temme_y(mu, x, y0, y1);
    else if (x < float256(kAsymX)) steed_y(mu, x, y0, y1);
    else {
        dd J0, J1;
        if (!hankel_jy(mu, x, J0, y0) || !hankel_jy(mu + float256::one(), x, J1, y1)) steed_y(mu, x, y0, y1);
    }
    return recur(y0, y1, mu, x, n, -1);
}

// I_ν(x), ν >= 0, x > 0 finite.
[[nodiscard]] inline SD i(float256 nu, float256 x) {
    if (x >= float256(kAsymX)) {
        dd S, Q;
        if (asym(nu, x, 1, S, Q)) {
            const SD e = exp_sd(dd(x));
            const dd amp = detail_dd::sqrt(dd(float256::one()) / (detail_dd::pi() * x * float256::two()));
            return mul(e, amp * S);
        }
    }
    return series(dd(nu), x, +1);
}

// K_ν(x), ν >= 0, x > 0 finite.
[[nodiscard]] inline SD k(float256 nu, float256 x) {
    const dd X{x};
    if (x >= float256(kAsymX)) {
        dd S, Q;
        if (asym(nu, x, 2, S, Q)) {
            const SD e = exp_sd(-X);
            return mul(e, detail_dd::sqrt(detail_dd::pi() / (X * float256::two())) * S);
        }
    }
    int n; float256 mu;
    split(nu, n, mu);
    dd k0, k1;
    int ke = 0;
    if (x < float256(kTemmeKX)) temme_k(mu, x, k0, k1);
    else if (x < float256(kAsymX)) k_wronskian(mu, x, k0, k1, ke);
    else steed_k(mu, x, k0, k1, ke);
    SD r = recur(k0, k1, mu, x, n, +1);
    r.k += ke;
    return r;
}

} // namespace bessel
} // namespace detail_sf

namespace detail_sf {
// Common argument handling for the cylinder functions. Returns true and sets
// `r` when the result is decided without evaluation.
enum class BesselKind { J, Y, I, K };
[[nodiscard]] inline bool bessel_special(detail_dd::FnEnv& env, BesselKind kind, float256 nu, float256 x, float256& r) {
    if (nu.isnan() || x.isnan()) { r = env.nan2(nu, x); return true; }
    if (nu.signbit() && !nu.iszero()) { r = env.invalid(); return true; }
    if (x.signbit() && !x.iszero()) { r = env.invalid(); return true; }
    if (nu.isinf()) { r = env.invalid(); return true; }
    if (x.isinf()) {
        r = env.exact(kind == BesselKind::I ? float256::inf(false) : float256::zero());
        return true;
    }
    if (x.iszero()) {
        switch (kind) {
            case BesselKind::J:
            case BesselKind::I: r = env.exact(nu.iszero() ? float256::one() : float256::zero()); return true;
            case BesselKind::Y: r = env.pole(true); return true;
            case BesselKind::K: r = env.pole(false); return true;
        }
    }
    return false;
}
// |x|^n·2^s/D, when it lies exactly on the 238-bit grid (a binary256 number,
// or a midpoint between two), as hi + lo (lo = 0 or half an ulp of hi) times 2^k.
// The first series term of J/I and j_n for tiny x is such a monomial; there the
// dd evaluation (error ~2^-470) can sit on the wrong side of the grid point,
// while the true value is off it by only ~x^2 (~2^-600), so the tail decides
// the rounding and the point itself must be exact.
[[nodiscard]] inline bool exact_monomial(float256 x, unsigned n, int s, const detail::BigInt& D, float256& hi,
                                         float256& lo, int& k) {
    if (n == 0 || n > 60 || x.iszero() || !x.isfinite()) return false;
    const auto w = abs(x).raw_bits();
    const long be = static_cast<long>((w[3] >> 44) & 0x7FFFF);
    std::vector<std::uint32_t> limbs;
    for (int i = 0; i < 4; ++i) {
        std::uint64_t v = w[i];
        if (i == 3) v = (v & 0xFFFFFFFFFFFULL) | (be ? (1ULL << 44) : 0);
        limbs.push_back(static_cast<std::uint32_t>(v));
        limbs.push_back(static_cast<std::uint32_t>(v >> 32));
    }
    detail::BigInt M = detail::BigInt::from_limbs(limbs);
    long E = (be ? be : 1) - 262143 - 236;
    while (!M.is_zero() && !M.bit(0)) { M.shr(1); ++E; }
    // A 238-bit quotient needs M^n <= 2^238·D: skip long significands early.
    if (static_cast<long>(M.bit_length() - 1) * n > static_cast<long>(238 + D.bit_length())) return false;
    detail::BigInt N(1);
    for (unsigned i = 0; i < n; ++i) N.mul(M);
    detail::BigInt Dd = D;
    long Et = E * static_cast<long>(n) + s;
    while (!Dd.is_zero() && !Dd.bit(0)) { Dd.shr(1); --Et; }
    const detail::BigInt rem = N.divmod(Dd);
    if (!rem.is_zero() || N.bit_length() > 238) return false;
    while (!N.bit(0)) { N.shr(1); ++Et; }
    const bool mid = N.bit_length() == 238;
    if (mid) { N.shr(1); ++Et; }  // (Q − 1)/2, Q odd: 237 bits
    float256 h{};
    for (int i = N.limb_count() - 1; i >= 0; --i) h = h + ldexp(float256::from_uint(N.limb(i)), 32 * i);
    if (Et < -400000 || Et > 400000) return false;
    k = static_cast<int>(Et);
    if (mid) { hi = ldexp(h, 1); lo = float256::one(); --k; }  // (2·(Q−1)/2 + 1)·2^(k−1)
    else { hi = h; lo = float256::zero(); }
    return true;
}
// J or I for x so small that the series is its first term: (x/2)^ν/Γ(ν+1)·(1 ∓ tiny).
[[nodiscard]] inline bool bessel_tiny(detail_dd::FnEnv& env, BesselKind kind, float256 nu, float256 x, float256& r) {
    if ((kind != BesselKind::J && kind != BesselKind::I) || ilogb(x) >= -120) return false;
    const int tail = kind == BesselKind::J ? -1 : +1;  // x > 0 here (x < 0 is reflected)
    // Below 2^-120 the correction x²/(4(ν+1)) is under a quarter ulp, so an exact
    // grid-point first term rounds by the tail alone.
    if (nu == trunc(nu) && nu >= float256::one() && nu <= float256(60)) {
        // integer ν: (x/2)^ν/ν! can be exactly representable (e.g. J_1 = x/2·(1 − x²/8))
        const unsigned n = static_cast<unsigned>(nu.to_double());
        detail::BigInt fact(1);
        for (unsigned i = 2; i <= n; ++i) fact.mul_u32(i);
        float256 hi, lo;
        int k = 0;
        if (exact_monomial(x, n, -static_cast<int>(n), fact, hi, lo, k)) {
            r = env.round(hi, lo, k, tail);
            return true;
        }
    }
    if (ilogb(x) >= -300) return false;
    SD v = bessel::series(dd(nu), x, kind == BesselKind::J ? -1 : 1);
    r = out(env, v, tail);
    return true;
}
} // namespace detail_sf

namespace detail_sf {
// −f() rounded in the current mode: f evaluated in the mirrored mode, then negated
// (round_up(−v) = −round_down(v)); flags are the same. f is a public function.
template <class F>
[[nodiscard]] inline float256 negated(F f) {
    const round_mode m = fenv::get_round();
    if (m == round_mode::upward || m == round_mode::downward)
        fenv::set_round(m == round_mode::upward ? round_mode::downward : round_mode::upward);
    const float256 v = f();
    fenv::set_round(m);
    return -v;
}
// x < 0 (or −0) and integer ν >= 0, where C_ν(−x) = (−1)^ν C_ν(x) for C = J, I.
[[nodiscard]] inline bool bessel_reflect(float256 nu, float256 x) {
    return x.signbit() && !x.isnan() && !nu.isnan() && !nu.isinf() && !nu.signbit() && nu == trunc(nu);
}
[[nodiscard]] inline bool odd_int(float256 nu) { return ldexp(nu, -1) != trunc(ldexp(nu, -1)); }
} // namespace detail_sf

/// Bessel function of the first kind J_ν(x), ν >= 0. x < 0: for integer ν,
/// J_ν(−x) = (−1)^ν J_ν(x) (so J_ν(−0) = −0 for odd ν); for non-integer ν the
/// value is complex, so NaN with `invalid` (C++17 leaves x < 0 to the implementation).
[[nodiscard]] inline float256 cyl_bessel_j(float256 nu, float256 x) {
    if (detail_sf::bessel_reflect(nu, x)) {
        if (detail_sf::odd_int(nu)) return detail_sf::negated([&] { return cyl_bessel_j(nu, -x); });
        return cyl_bessel_j(nu, -x);
    }
    detail_dd::FnEnv env;
    float256 r;
    if (detail_sf::bessel_special(env, detail_sf::BesselKind::J, nu, x, r)) return r;
    if (detail_sf::bessel_tiny(env, detail_sf::BesselKind::J, nu, x, r)) return r;
    return detail_sf::out(env, detail_sf::bessel::j(nu, x));
}
/// Bessel function of the second kind (Neumann) Y_ν(x), ν >= 0, x >= 0; x < 0
/// is complex for every ν, so NaN with `invalid`.
[[nodiscard]] inline float256 cyl_neumann(float256 nu, float256 x) {
    detail_dd::FnEnv env;
    float256 r;
    if (detail_sf::bessel_special(env, detail_sf::BesselKind::Y, nu, x, r)) return r;
    return detail_sf::out(env, detail_sf::bessel::y(nu, x));
}
/// Modified Bessel function of the first kind I_ν(x), ν >= 0. x < 0: for integer
/// ν, I_ν(−x) = (−1)^ν I_ν(x); for non-integer ν NaN with `invalid` (complex).
[[nodiscard]] inline float256 cyl_bessel_i(float256 nu, float256 x) {
    if (detail_sf::bessel_reflect(nu, x)) {
        if (detail_sf::odd_int(nu)) return detail_sf::negated([&] { return cyl_bessel_i(nu, -x); });
        return cyl_bessel_i(nu, -x);
    }
    detail_dd::FnEnv env;
    float256 r;
    if (detail_sf::bessel_special(env, detail_sf::BesselKind::I, nu, x, r)) return r;
    if (detail_sf::bessel_tiny(env, detail_sf::BesselKind::I, nu, x, r)) return r;
    return detail_sf::out(env, detail_sf::bessel::i(nu, x));
}
/// Modified Bessel function of the second kind K_ν(x), ν >= 0, x >= 0; x < 0
/// is complex for every ν, so NaN with `invalid`.
[[nodiscard]] inline float256 cyl_bessel_k(float256 nu, float256 x) {
    detail_dd::FnEnv env;
    float256 r;
    if (detail_sf::bessel_special(env, detail_sf::BesselKind::K, nu, x, r)) return r;
    return detail_sf::out(env, detail_sf::bessel::k(nu, x));
}

namespace detail_sf {
// sqrt(π/(2x)) · C_{n+1/2}(x)
[[nodiscard]] inline float256 sph(unsigned n, float256 x, bool second) {
    detail_dd::FnEnv env;
    const float256 nu = float256::from_uint(n) + float256::half();
    if (x.isnan()) return env.nan(x);
    if (x.signbit() && !x.iszero()) return env.invalid();
    if (x.isinf()) return env.exact(float256::zero());
    if (x.iszero()) {
        if (second) return env.pole(true);
        return env.exact(n == 0 ? float256::one() : float256::zero());
    }
    if (!second && ilogb(x) < -120 && n >= 1) {
        // j_n(x) = x^n/(2n+1)!!·(1 − x²/(2(2n+3)) + …): exact first term, if on the grid
        detail::BigInt df(1);
        for (unsigned i = 3; i <= 2 * n + 1; i += 2) df.mul_u32(i);
        float256 hi, lo;
        int k = 0;
        if (exact_monomial(x, n, 0, df, hi, lo, k)) return env.round(hi, lo, k, -1);
    }
    SD v = second ? bessel::y(nu, x) : (ilogb(x) < -300 ? bessel::series(dd(nu), x, -1) : bessel::j(nu, x));
    v = mul(v, detail_dd::sqrt(detail_dd::pi() / (dd(x) * float256::two())));
    return out(env, v, (!second && ilogb(x) < -300) ? -1 : 2);
}
} // namespace detail_sf

/// Spherical Bessel function j_n(x), any real x: j_n(−x) = (−1)^n j_n(x).
[[nodiscard]] inline float256 sph_bessel(unsigned n, float256 x) {
    if (x.signbit() && !x.isnan()) {
        if (n & 1) return detail_sf::negated([&] { return detail_sf::sph(n, -x, false); });
        return detail_sf::sph(n, -x, false);
    }
    return detail_sf::sph(n, x, false);
}
/// Spherical Neumann function y_n(x), any real x: y_n(−x) = (−1)^{n+1} y_n(x)
/// (y_n(±0) is a pole, −∞ at +0, (−1)^n·∞ at −0).
[[nodiscard]] inline float256 sph_neumann(unsigned n, float256 x) {
    if (x.signbit() && !x.isnan()) {
        if (!(n & 1)) return detail_sf::negated([&] { return detail_sf::sph(n, -x, true); });
        return detail_sf::sph(n, -x, true);
    }
    return detail_sf::sph(n, x, true);
}

/// J_n(x) for any integer n and real x (J_{−n} = (−1)^n J_n, J_n(−x) = (−1)^n J_n(x)).
[[nodiscard]] inline float256 jn(int n, float256 x) {
    const bool odd = (n & 1) != 0;
    const unsigned an = n < 0 ? 0u - static_cast<unsigned>(n) : static_cast<unsigned>(n);
    const bool neg = odd && ((n < 0) != x.signbit());
    if (x.isnan()) return cyl_bessel_j(float256::zero(), x);
    const float256 r = cyl_bessel_j(float256::from_uint(an), abs(x));
    if (!neg) return r;
    // Rounding is symmetric except in the directed modes: evaluate −J in the flipped mode.
    const round_mode m = fenv::get_round();
    if (m == round_mode::upward || m == round_mode::downward) {
        fenv::set_round(m == round_mode::upward ? round_mode::downward : round_mode::upward);
        const float256 v = cyl_bessel_j(float256::from_uint(an), abs(x));
        fenv::set_round(m);
        return -v;
    }
    return -r;
}
/// Y_n(x) for any integer n and x >= 0 (Y_{−n} = (−1)^n Y_n).
[[nodiscard]] inline float256 yn(int n, float256 x) {
    const unsigned an = n < 0 ? 0u - static_cast<unsigned>(n) : static_cast<unsigned>(n);
    const bool neg = n < 0 && (n & 1) != 0;
    if (!neg) return cyl_neumann(float256::from_uint(an), x);
    const round_mode m = fenv::get_round();
    if (m == round_mode::upward || m == round_mode::downward)
        fenv::set_round(m == round_mode::upward ? round_mode::downward : round_mode::upward);
    const float256 v = cyl_neumann(float256::from_uint(an), x);
    fenv::set_round(m);
    return -v;
}

namespace detail_sf {
// ln Γ(z) for a dd z > 0.
[[nodiscard]] inline dd lgam(const dd& z) { return detail_gamma::lgamma_pos_dd(z); }
// B(a, b) = Γ(a)Γ(b)/Γ(a+b) as m·2^k for a, b > 0 finite (a + b exact as dd).
[[nodiscard]] inline SD beta_sd(float256 a, float256 b) {
    const dd ab = detail_dd::two_sum(a, b);
    const float256 big = ldexp(float256::one(), 20);
    if (a < big && b < big) {
        const SD ga = gamma_sd(dd(a)), gb = gamma_sd(dd(b)), gab = gamma_sd(ab);
        return norm({ga.m * gb.m / gab.m, ga.k + gb.k - gab.k});
    }
    return exp_sd(lgam(dd(a)) + lgam(dd(b)) - lgam(ab));
}

// Continued fraction for I_x(a, b): h = 1/(1 + e_1/(1 + e_2/(1 + …))) with
// e_n = N_n/Q_n, Q_n = (a+n−1)(a+n), N_{2m} = m(b−m)x, N_{2m+1} = −(a+m)(a+b+m)x.
// The equivalence transform c_n = Q_n clears the denominators (partial
// numerators Q_{n−1}N_n, denominators Q_n), so the Wallis recurrence for the
// convergents A_n/B_n needs no divisions; A, B are rescaled by powers of two,
// and h = B_n/A_n is tested for convergence every 8 steps.
[[nodiscard]] inline dd betacf(const dd& a, const dd& b, const dd& x) {
    const dd one{float256::one()};
    const dd ab = a + b;
    dd A0 = one, A1 = one, B0{float256::zero()}, B1 = one;  // A_{n−2}, A_{n−1}, B_{n−2}, B_{n−1}
    dd Qprev = one;  // Q_{n−1} (c_0 = 1)
    dd h = one, hprev{float256::zero()};
    for (int n = 1; n < 100000000; ++n) {
        const int m = n >> 1;
        const dd Q = (a + D(n - 1)) * (a + D(n));
        const dd N = (n & 1) ? -((a + D(m)) * (ab + D(m)) * x) : D(m) * (b - D(m)) * x;
        const dd an = Qprev * N;
        const dd A2 = Q * A1 + an * A0, B2 = Q * B1 + an * B0;
        A0 = A1; A1 = A2; B0 = B1; B1 = B2;
        Qprev = Q;
        const int e = ilogb(A1.hi);
        if (e > 2000 || e < -2000) {
            A0 = detail_dd::ldexp(A0, -e); A1 = detail_dd::ldexp(A1, -e);
            B0 = detail_dd::ldexp(B0, -e); B1 = detail_dd::ldexp(B1, -e);
        }
        if ((n & 7) == 0) {
            h = B1 / A1;
            if (near_one(h / hprev)) break;
            hprev = h;
        }
    }
    return h;
}
} // namespace detail_sf

/// Beta function B(x, y) = Γ(x)Γ(y)/Γ(x+y) for x, y > 0. A zero argument is a
/// pole (+∞, divide-by-zero); negative arguments are a domain error.
[[nodiscard]] inline float256 beta(float256 x, float256 y) {
    detail_dd::FnEnv env;
    if (x.isnan() || y.isnan()) return env.nan2(x, y);
    if ((x.signbit() && !x.iszero()) || (y.signbit() && !y.iszero())) return env.invalid();
    if (x.iszero() || y.iszero()) return env.pole(false);
    if (x.isinf() || y.isinf()) return env.exact(float256::zero());
    if (x == float256::one() || y == float256::one()) {  // B(1, y) = 1/y
        const detail_dd::dd q = detail_dd::dd(float256::one()) / (x == float256::one() ? y : x);
        return env.round(q.hi, q.lo, 0, q.lo.iszero() ? 0 : 2);
    }
    return detail_sf::out(env, detail_sf::beta_sd(x, y));
}

/// Regularized incomplete beta I_x(a, b) = B(x; a, b)/B(a, b) for a, b > 0 and
/// 0 <= x <= 1 (Boost's ibeta); ibetac = 1 − I_x(a, b) computed directly.
namespace detail_sf {
[[nodiscard]] inline float256 ibeta_impl(float256 a, float256 b, float256 x, bool complement) {
    detail_dd::FnEnv env;
    if (a.isnan() || b.isnan() || x.isnan()) {
        if (a.isnan() || b.isnan()) return env.nan2(a, b);
        return env.nan(x);
    }
    if (a.signbit() || b.signbit() || a.iszero() || b.iszero() || a.isinf() || b.isinf()) return env.invalid();
    if (x.signbit() && !x.iszero()) return env.invalid();
    if (x > float256::one()) return env.invalid();
    const float256 one = float256::one();
    if (x.iszero()) return env.exact(complement ? one : float256::zero());
    if (x == one) return env.exact(complement ? float256::zero() : one);
    const dd X{x};
    const dd Xc = detail_dd::two_sum(one, -x);
    // x^a (1−x)^b / B(a, b)
    const dd L = dd(a) * detail_dd::log(X) + dd(b) * detail_dd::log(Xc);
    const SD bt = beta_sd(a, b);
    const SD front = mul(exp_sd(L), SD{dd(one) / bt.m, -bt.k});
    // Use the fraction that converges fast: x < (a+1)/(a+b+2).
    const dd lhs = X * (dd(a) + dd(b) + D(2)) - (dd(a) + dd(one));
    const bool swap = !lhs.hi.signbit() && !lhs.hi.iszero();
    SD part;
    if (!swap) part = mul(front, betacf(dd(a), dd(b), X) / dd(a));
    else part = mul(front, betacf(dd(b), dd(a), Xc) / dd(b));
    // part = I (no swap) or 1 − I (swap)
    if (swap == complement) return out(env, part);
    if (part.k < -400) return env.round(one, float256::zero(), 0, -1);  // 1 − tiny
    return (dd(one) - detail_dd::ldexp(part.m, part.k)).out(env);
}
} // namespace detail_sf
[[nodiscard]] inline float256 ibeta(float256 a, float256 b, float256 x) { return detail_sf::ibeta_impl(a, b, x, false); }
[[nodiscard]] inline float256 ibetac(float256 a, float256 b, float256 x) { return detail_sf::ibeta_impl(a, b, x, true); }

namespace detail_sf {
namespace zeta {

// B_2k/(2k)!, k >= 1 (MPFR table, k <= 140).
[[nodiscard]] inline dd bern_fact(int k) { return detail_dd::dd_from_words(detail_dd::tables::kBernFact[k - 1]); }

// Plans Euler–Maclaurin at split point w: the first k <= 140 whose term
// B_2k/(2k)! (s)_{2k−1} w^{−s−2k+1} is below 2^-506 of the reference
// 2^-ref_lg · w^{1−s}, estimated in double; −1 if none is.
[[nodiscard]] inline int em_plan(double s, double w, double ref_lg) {
    // log2 |(s)_{2k−1}| accumulated factor by factor (any real s; a zero factor
    // at a non-positive integer s makes every later term vanish).
    const double l2pw = std::log2(6.283185307179586 * w);
    double lp = 0;
    for (int k = 1; k <= detail_dd::tables::kBernFactN; ++k) {
        for (int i = (k == 1 ? 0 : 2 * k - 3); i <= 2 * k - 2; ++i) {
            const double f = std::fabs(s + i);
            if (f == 0) return k;
            lp += std::log2(f);
        }
        if (1 + lp - 2 * k * l2pw + ref_lg < -506) return k;
    }
    return -1;
}

// Σ_{k=1}^{K} B_2k/(2k)! (s)_{2k−1} w^{−2k} (the Euler–Maclaurin corrections
// divided by w^{1−s}), stopping early once a term is below 2^-490 of `ref`.
[[nodiscard]] inline dd em_corr(const dd& s, const dd& w, int K, const dd& ref) {
    const dd iw = dd(float256::one()) / w;
    const dd iw2 = iw * iw;
    dd t = s * iw2, sum{float256::zero()};
    for (int k = 1; k <= K; ++k) {
        const dd term = bern_fact(k) * t;
        sum = sum + term;
        if (small_rel(term, ref + sum)) break;
        t = t * (s + D(2 * k - 1)) * (s + D(2 * k)) * iw2;
    }
    return sum;
}

// ζ(s) for a dd s > 0, s ≠ 1, given sm1 = s − 1 exactly. Direct sum to N − 1
// (n^{−s} from prime powers, one exp per prime, ln p from the table) plus
// Euler–Maclaurin at the smallest N in {24, 32, …, 64} the plan accepts.
[[nodiscard]] inline dd pos(const dd& s, const dd& sm1) {
    const dd one{float256::one()};
    if (s.hi >= float256(240)) {
        // 1 + 2^-s + 3^-s + …, terms below 2^-500 dropped
        dd sum = one;
        for (int n = 2; n < 1000; ++n) {
            const dd t = detail_dd::exp(-(s * detail_dd::log(D(n))));
            if (ilogb(t.hi) < -500) break;
            sum = sum + t;
        }
        return sum;
    }
    const double sd = s.hi.to_double();
    int N = 64, K = 140;
    for (int c = 24; c <= 64; c += 8) {
        const int k = em_plan(sd, c, 1 + (1 - sd) * std::log2(double(c)));
        if (k > 0) { N = c; K = k; break; }
    }
    dd pw[64];
    pw[1] = one;
    int pi = 0;
    for (int n = 2; n < N; ++n) {
        int p = 2;
        while (p * p <= n && n % p) ++p;
        if (p * p > n) pw[n] = detail_dd::exp(-(s * detail_dd::dd_from_words(detail_dd::tables::kLnPrime64[pi++])));
        else pw[n] = pw[p] * pw[n / p];
    }
    dd sum{float256::zero()};
    for (int n = N - 1; n >= 1; --n) sum = sum + pw[n];
    // N^{1−s} (1/(s−1) + 1/(2N) + Σ_k B_2k/(2k)! (s)_{2k−1} N^{−2k})
    const dd W = D(N);
    const dd N1s = detail_dd::exp(-(sm1 * detail_dd::log(W)));
    const dd br = one / sm1 + one / (W * float256::two());
    return sum + N1s * (br + em_corr(s, W, K, br + sum / N1s));
}

// Hurwitz ζ(s, a) = Σ_{n>=0} (n+a)^{−s}, s > 1, a > 0 finite, as m·2^k.
[[nodiscard]] inline SD hurwitz(const dd& s, const dd& sm1, float256 a, double* cancel = nullptr) {
    if (cancel) *cancel = 0;
    const dd one{float256::one()};
    const double sd = s.hi.to_double(), ad = a.to_double();
    const double lsm1 = std::log2(std::fabs(sm1.hi.to_double()));
    const bool neg = s.hi.signbit();  // s < 0: terms grow with n; scale by w^{-s} instead of a^{-s}
    // Option 1: direct sum of nd terms (rest < 2^-506 relative): s log2(1 + nd/a) >= 506 + log2(1 + (a+nd)/(s−1)).
    long long nd = -1;
    if (sm1.hi > float256::zero()) {
        long long lo = 1;
        while (lo < (1LL << 40)) {
            const double n = double(lo);
            const double lg = sd * std::log1p(n / ad) * 1.4426950408889634 - std::log2(1 + (ad + n) / sm1.hi.to_double());
            if (lg >= 506) break;
            lo *= 2;
        }
        if (lo < (1LL << 40)) {
            long long hi = lo; lo = lo / 2;  // refine
            while (hi - lo > 1) {
                const long long mid = (lo + hi) / 2;
                const double n = double(mid);
                const double lg = sd * std::log1p(n / ad) * 1.4426950408889634 - std::log2(1 + (ad + n) / sm1.hi.to_double());
                if (lg >= 506) hi = mid; else lo = mid;
            }
            nd = hi;
        }
    }
    // Option 2: N direct terms + Euler–Maclaurin at w = N + a. An exp costs ~20 correction terms.
    int bestN = -1, bestK = 0;
    double bestCost = nd > 0 ? double(nd) : 1e300;
    for (int N : {0, 4, 8, 16, 24, 32, 48, 64, 96, 128, 192, 256, 384, 512, 768, 1024}) {
        const int k = em_plan(sd, double(N) + ad, -lsm1);
        if (k < 0) continue;
        const double cost = N + k / 20.0;
        if (cost < bestCost) { bestCost = cost; bestN = N; bestK = k; }
    }
    if (bestN < 0 && nd < 0) { bestN = 1024; bestK = detail_dd::tables::kBernFactN; }
    const dd A{a};
    const dd la = detail_dd::log(A);
    if (bestN == 0) {
        // a^{1−s} (1/(s−1) + 1/(2a) + Σ_k B_2k/(2k)! (s)_{2k−1} a^{−2k})
        const dd br = one / sm1 + one / (A * float256::two());
        const dd v = br + em_corr(s, A, bestK, br);
        if (cancel && !v.hi.iszero()) {
            const float256 big = abs(one.hi / sm1.hi) > abs(one.hi / (a + a)) ? abs(one.hi / sm1.hi) : abs(one.hi / (a + a));
            *cancel = ilogb(big) - ilogb(v.hi);
        }
        return mul(exp_sd(-(sm1 * la)), v);
    }
    // Relative to c^{−s}, c = a (s > 0) or c = w = N + a (s < 0), so every
    // term is <= 1: Σ_{n<N} ((n+a)/c)^{−s} [+ EM tail at w].
    const long long Nt = bestN > 0 ? bestN : nd;
    const dd W = D(static_cast<int>(bestN > 0 ? bestN : 0)) + A;
    const dd C = neg ? W : A;
    const dd lc = neg ? detail_dd::log(C) : la;
    dd sum{float256::zero()};
    // integer s in [2, 512] (polygamma, trigamma, ζ(n, a)): (c/(n+a))^s by binary
    // powering, about (s + 2·log2 s)·2^-476 relative, against a log and an exp
    // (error ~2^-474·s·ln((n+a)/c)) that cost ~20 dd multiplications each
    const bool sint = !neg && s.lo.iszero() && s.hi == trunc(s.hi) && s.hi <= float256(512);
    const unsigned si = sint ? static_cast<unsigned>(s.hi.to_double()) : 0;
    for (long long n = Nt - 1; n >= 0; --n) {
        if (n == 0 && !neg) { sum = sum + one; break; }
        const dd na = dd(float256::from_int(static_cast<int>(n))) + A;
        if (sint) {
            dd b = C / na, p = one;
            for (unsigned e = si;;) {
                if (e & 1) p = p * b;
                e >>= 1;
                if (e == 0) break;
                b = b * b;
            }
            sum = sum + p;
            continue;
        }
        const dd r = na / C;
        sum = sum + detail_dd::exp(-(s * detail_dd::log(r)));
    }
    if (bestN > 0) {
        // (w/c)^{−s} · w · (1/(s−1) + 1/(2w) + corrections)
        const dd f = neg ? W : detail_dd::exp(-(s * detail_dd::log(W / A))) * W;
        const dd br = one / sm1 + one / (W * float256::two());
        const dd tail = f * (br + em_corr(s, W, bestK, br));
        sum = sum + tail;
        if (cancel && !sum.hi.iszero()) *cancel = (ilogb(tail.hi) > 0 ? ilogb(tail.hi) : 0) - ilogb(sum.hi);
    }
    return mul(exp_sd(-(s * lc)), sum);
}

// ζ(s, a) for s <= −kReflS (σ = 1 − s >= 1 + kReflS), a > 0, where the
// Euler–Maclaurin sum would cancel more than ~e^{|s|}. Hurwitz's formula for
// 0 < b <= 1:
//   ζ(s, b) = 2 Γ(σ) (2π)^{−σ} Σ_{n>=1} cos(πσ/2 − 2πnb) n^{−σ}
// (n^{−σ} < 2^-506 beyond n = 49; cos/sin(2πnb) by rotation), and for a > 1
// the shift ζ(s, a) = ζ(s, b) − Σ_{j<m} (b+j)^{−s}, b = a − m ∈ (0, 1], all
// combined in the log domain. The result has an absolute error ~2^-470 of
// its largest part, so values near a zero in a lose relative accuracy.
inline constexpr int kReflS = 89;
[[nodiscard]] inline SD hurwitz_refl(float256 s, float256 a, double* cancel = nullptr) {
    const dd one{float256::one()};
    if (cancel) *cancel = 0;
    // σ = 1 − s as an exact dd: in float256 it needs 238 bits when |s| + 1 crosses
    // a power of two (s ∈ (−128, −127], (−256, −255], … with the last bit set),
    // and the rounded σ cost up to ~650 ulps (hz_vneg, s = −127.2338…).
    const dd Sg = detail_dd::two_sum(float256::one(), -s);
    float256 m = ceil(a) - float256::one();
    if (m.signbit()) m = float256::zero();
    const float256 b = a - m;  // exact, in (0, 1]
    // S = Σ cos(πσ/2 − 2πnb) n^{−σ} = Σ [cos(πσ/2) cos(2πnb) + sin(πσ/2) sin(2πnb)] n^{−σ},
    // with sin(πσ/2) = cos(πs/2) and cos(πσ/2) = sin(πs/2) (s/2 is exact)
    dd sa, ca, st, ct;
    sincos_pi(ldexp(s, -1), ca, sa);
    sincos_pi(b + b, st, ct);
    dd cn = ct, sn = st, S{float256::zero()};
    dd pw[50];
    int pi = 0;
    for (int n = 1; n < 50; ++n) {
        if (n == 1) pw[1] = one;
        else {
            int p = 2;
            while (p * p <= n && n % p) ++p;
            if (p * p > n) pw[n] = detail_dd::exp(-(Sg * detail_dd::dd_from_words(detail_dd::tables::kLnPrime64[pi++])));
            else pw[n] = pw[p] * pw[n / p];
        }
        S = S + (ca * cn + sa * sn) * pw[n];
        if (ilogb(pw[n].hi) < -506) break;
        const dd c2 = cn * ct - sn * st, s2 = sn * ct + cn * st;
        cn = c2; sn = s2;
    }
    // log magnitude of the reflection part: ln 2 + lnΓ(σ) − σ ln 2π + ln|S|
    const dd lg = detail_dd::ln2() + lgam(Sg) - Sg * detail_dd::log(detail_dd::pi() * float256::two());
    if (S.hi.iszero()) { if (cancel) *cancel = 1e9; return SD{dd(float256::zero()), 0}; }
    if (cancel) *cancel = -ilogb(S.hi);
    const long mi = static_cast<long>(m.to_double());
    if (mi == 0) return mul(exp_sd(lg), S);
    // minus Σ_{j<m} (b+j)^{|s|}: log terms −s ln(b+j); the largest is j = m−1.
    const dd neg_s{-s};
    const dd lmax = neg_s * detail_dd::log(dd(b) + D(static_cast<int>(mi - 1)));
    const dd lsS = lg + detail_dd::log(dabs(S));
    const dd M = lmax.hi > lsS.hi ? lmax : lsS;
    dd tot = (S.hi.signbit() ? -one : one) * detail_dd::exp(lsS - M);
    for (long j = mi - 1; j >= 0; --j) {
        const dd e = neg_s * detail_dd::log(dd(b) + D(static_cast<int>(j))) - M;
        if (e.hi < float256(-400)) break;
        tot = tot - detail_dd::exp(e);
    }
    if (cancel) {
        // errors: ~2^-470 of each unit-size part, and of e^{lsS−M}/|S| for the S part
        const double sw = (lsS.hi - M.hi).to_double() / 0.6931471805599453 - static_cast<double>(ilogb(S.hi));
        *cancel = tot.hi.iszero() ? 1e9 : (sw > 0 ? sw : 0) - ilogb(tot.hi);
    }
    return mul(exp_sd(M), tot);
}

// Multiprecision fallbacks for ζ(s, a), s < 1, next to a zero in a, where the
// dd result above cancels below its error. Same formulas at P bits; `big`
// returns the magnitude of the largest part (to measure the cancellation).
namespace zmp {
namespace M = detail_mp;
[[nodiscard]] inline M::MP pow_neg(const M::MP& x, const M::MP& s, int Q) { return M::exp(M::neg(M::mul(s, M::log(x, Q), Q)), Q); }
// s as an integer in [1, 2^20] if it is one exactly, else 0.
[[nodiscard]] inline std::uint32_t small_uint(const M::MP& s) {
    if (s.zero() || s.neg || s.top() > 21) return 0;
    auto m = s.m;
    if (s.e < 0) {
        if (m.shr_sticky(static_cast<int>(-s.e))) return 0;  // fraction bits
    } else {
        m.shl(static_cast<int>(s.e));
    }
    return m.limb(0);
}
// x^{−s} for an integer s >= 1 by binary powering and one division: about
// 2·log2(s) roundings of 2^-(Q+32) and one of 2^-Q, against an exp and a log
// (each several multiplications at Q bits) in pow_neg.
[[nodiscard]] inline M::MP pow_neg_int(const M::MP& x, std::uint32_t s, int Q) {
    const int Q2 = Q + 32;
    M::MP r = M::from_u64(1), b = x;
    for (std::uint32_t e = s;;) {
        if (e & 1) r = M::mul(r, b, Q2);
        e >>= 1;
        if (e == 0) break;
        b = M::mul(b, b, Q2);
    }
    return M::div(M::from_u64(1), r, Q);
}
[[nodiscard]] inline M::MP absmax(const M::MP& a, const M::MP& b) { return M::cmp_abs(a, b) >= 0 ? a : b; }
// t − 2·round(t/2) ∈ [−1, 1] (exact).
[[nodiscard]] inline M::MP mod2(M::MP t) {
    if (t.zero() || t.e >= 1) return M::MP{};  // e >= 1: an even integer
    // integer part bits: m >> −e; keep the low bit of the integer part and the fraction
    const int fb = static_cast<int>(-t.e);
    M::MP r = t;
    // reduce r.m modulo 2^{fb+1} exactly
    auto q = r.m; q.shr(fb + 1); q.shl(fb + 1); r.m.sub(q);
    // r in [0, 2) for t >= 0 (sign kept): map to [−1, 1)
    M::MP one = M::from_u64(1);
    if (M::cmp_abs(r, one) >= 0) r = M::sub(r, M::from_u64(2, r.neg), 1 << 20);
    return r;
}
inline void sincos_pi(const M::MP& t, int Q, M::MP& s, M::MP& c) {
    M::sincos(M::mul(M::pi(Q), mod2(t), Q), Q, s, c);
}

// ζ(S, A) for an MP s and a > 0 (sd, ad: their values as doubles). Large s: the
// terms fall so fast that the direct sum is used.
[[nodiscard]] inline M::MP em_mp(const M::MP& S, const M::MP& A, double sd, double ad, int P, M::MP& big) {
    const int Q = P + 64;
    const M::MP one = M::from_u64(1);
    const std::uint32_t si = small_uint(S);  // integer s (polygamma): powers without exp/log
    auto pow_neg = [&](const M::MP& x, const M::MP& s, int q) { return si ? pow_neg_int(x, si, q) : zmp::pow_neg(x, s, q); };
    if (sd > 0.5 * P + 50) {
        M::MP sum;
        big = M::MP{};
        for (std::uint64_t n = 0;; ++n) {
            const M::MP t = pow_neg(M::add(A, M::from_u64(n), Q), S, Q);
            sum = M::add(sum, t, Q);
            big = absmax(big, t);
            if (t.zero() || t.top() < big.top() - Q) break;
        }
        big.neg = false;
        return sum;
    }
    const double w0 = 0.12 * P + std::fabs(sd) / 2 + 10;
    const long N = ad < w0 ? static_cast<long>(std::ceil(w0 - ad)) : 0;
    M::MP sum, W = A;
    big = M::MP{};
    for (long n = 0; n < N; ++n) {
        const M::MP t = pow_neg(M::add(A, M::from_u64(static_cast<std::uint64_t>(n)), Q), S, Q);
        sum = M::add(sum, t, Q);
        big = absmax(big, t);
    }
    W = M::add(A, M::from_u64(static_cast<std::uint64_t>(N)), Q);
    const M::MP Wms = pow_neg(W, S, Q), W1ms = M::mul(Wms, W, Q);
    const M::MP t1 = M::div(W1ms, M::sub(S, one, Q), Q), t2 = M::ldexp(Wms, -1);
    big = absmax(absmax(big, t1), t2);
    M::MP corr, poch = S, pw = M::div(one, W, Q);
    const M::MP iw2 = M::mul(pw, pw, Q);
    M::MP prev;
    for (int k = 1; k < 4000; ++k) {
        const M::MP term = M::mul(M::mul(M::bern_fact(k, Q), poch, Q), M::mul(pw, Wms, Q), Q);
        if (term.zero()) break;  // (s)_{2k−1} = 0 at a non-positive integer s
        if (k > 2 && M::cmp_abs(term, prev) > 0) break;
        corr = M::add(corr, term, Q);
        big = absmax(big, term);
        if (term.top() < big.top() - Q) break;
        prev = term;
        poch = M::mul(M::mul(poch, M::add(S, M::from_u64(static_cast<std::uint64_t>(2 * k - 1)), Q), Q),
                      M::add(S, M::from_u64(static_cast<std::uint64_t>(2 * k)), Q), Q);
        pw = M::mul(pw, iw2, Q);
    }
    big.neg = false;
    return M::add(M::add(M::add(sum, t1, Q), t2, Q), corr, Q);
}
[[nodiscard]] inline M::MP em(float256 s, float256 a, int P, M::MP& big) {
    return em_mp(M::from_f256(s), M::from_f256(a), s.to_double(), a.to_double(), P, big);
}

[[nodiscard]] inline M::MP refl(float256 s, float256 a, int P, M::MP& big) {
    const int Q = P + 64;
    const M::MP one = M::from_u64(1);
    const M::MP sig = M::sub(one, M::from_f256(s), 1 << 20);  // exact
    float256 m = ceil(a) - float256::one();
    if (m.signbit()) m = float256::zero();
    const float256 b = a - m;
    M::MP sa, ca, st, ct;
    sincos_pi(M::ldexp(sig, -1), Q, sa, ca);
    sincos_pi(M::ldexp(M::from_f256(b), 1), Q, st, ct);
    M::MP cn = ct, sn = st, S;
    std::vector<M::MP> pw(2);
    pw[1] = one;
    for (std::uint32_t n = 1;; ++n) {
        if (n >= 2) {
            std::uint32_t p = 2;
            while (p * p <= n && n % p) ++p;
            pw.push_back(p * p > n ? pow_neg(M::from_u64(n), sig, Q) : M::mul(pw[p], pw[n / p], Q));
        }
        S = M::add(S, M::mul(M::add(M::mul(ca, cn, Q), M::mul(sa, sn, Q), Q), pw[n], Q), Q);
        if (pw[n].zero() || pw[n].top() < -Q) break;
        const M::MP c2 = M::sub(M::mul(cn, ct, Q), M::mul(sn, st, Q), Q), s2 = M::add(M::mul(sn, ct, Q), M::mul(cn, st, Q), Q);
        cn = c2; sn = s2;
    }
    // 2 Γ(σ) (2π)^{−σ}
    const M::MP lg = M::sub(M::lgamma_pos(sig, Q), M::mul(sig, M::log(M::ldexp(M::pi(Q), 1), Q), Q), Q);
    const M::MP amp = M::ldexp(M::exp(lg, Q), 1);
    M::MP r = M::mul(amp, S, Q);
    big = amp;
    const long mi = static_cast<long>(m.to_double());
    const M::MP B = M::from_f256(b), Sm = M::from_f256(s);
    for (long j = mi - 1; j >= 0; --j) {
        const M::MP t = pow_neg(M::add(B, M::from_u64(static_cast<std::uint64_t>(j)), Q), Sm, Q);
        r = M::sub(r, t, Q);
        big = absmax(big, t);
        if (t.top() < big.top() - Q) break;
    }
    big.neg = false;
    return r;
}

// ζ(s, a) to ~330 bits after the cancellation, as m·2^k.
[[nodiscard]] inline SD hurwitz_mp(float256 s, float256 a, double cancel) {
    int P = 330 + 64 + static_cast<int>(cancel < 2000 ? cancel : 2000);
    for (;;) {
        M::MP big;
        const M::MP v = s <= float256(-kReflS) ? refl(s, a, P, big) : em(s, a, P, big);
        const long c = v.zero() ? 100000 : big.top() - v.top();
        if (P - c >= 330 || P > 6000) {
            if (v.zero()) return SD{dd(float256::zero()), 0};
            const long k = v.top();
            return norm({M::to_dd(M::ldexp(v, -k)), static_cast<int>(k)});
        }
        P = static_cast<int>(c) + 400 > 2 * P ? static_cast<int>(c) + 400 : 2 * P;
    }
}
} // namespace zmp

} // namespace zeta
} // namespace detail_sf

/// Riemann zeta function ζ(x) for all real x: Euler–Maclaurin for x >= 0 and
/// the functional equation ζ(x) = 2^x π^{x−1} sin(πx/2) Γ(1−x) ζ(1−x) for x < 0.
/// ζ(1) is a pole (+∞, divide-by-zero); ζ(−2n) = +0 exactly; ζ(0) = −1/2.
[[nodiscard]] inline float256 riemann_zeta(float256 x) {
    detail_dd::FnEnv env;
    using detail_dd::dd;
    const float256 one = float256::one();
    if (x.isnan()) return env.nan(x);
    if (x.isinf()) return x.signbit() ? env.invalid() : env.exact(one);
    if (x == one) return env.pole(false);
    if (x.iszero()) return env.exact(-float256::half());
    if (!x.signbit()) {
        if (x > float256(600)) return env.round(one, float256::zero(), 0, +1);  // 1 + 2^-x, x > 600
        const dd v = detail_sf::zeta::pos(dd(x), detail_dd::two_sum(x, -one));
        return v.out(env);
    }
    // x < 0
    const float256 h = ldexp(x, -1);
    if (trunc(h) == h) return env.exact(float256::zero());  // trivial zero
    if (x < float256(-30000)) {
        // |ζ(x)| ~ 2Γ(1−x)/(2π)^{1−x} > 2^262144; sign of sin(πx/2)
        dd s, c;
        detail_sf::sincos_pi(h, s, c);
        return env.overflow(s.hi.signbit());
    }
    const dd oms = detail_dd::two_sum(one, -x);  // 1 − x, exact
    dd s, c;
    detail_sf::sincos_pi(h, s, c);
    const dd z = detail_sf::zeta::pos(oms, dd(-x));
    const detail_sf::SD g = detail_sf::gamma_sd(oms);
    const dd L = dd(x) * detail_dd::ln2() + (dd(x) - dd(one)) * detail_dd::log(detail_dd::pi());
    detail_sf::SD r = detail_sf::mul(detail_sf::exp_sd(L), g);
    r = detail_sf::mul(r, s * z);
    return detail_sf::out(env, r);
}

/// Hurwitz zeta function ζ(s, a) = Σ_{n>=0} (n + a)^{−s} for real s ≠ 1 and
/// a > 0 (analytic continuation for s < 1). s = 1 is a pole (+∞,
/// divide-by-zero); a <= 0 and s = −∞ are domain errors. ζ(s, 1) = ζ(s),
/// ζ(0, a) = 1/2 − a. For s < 1 the error is absolute relative to the
/// largest part of the sum, so results near a zero in a (e.g. ζ(0, 1/2))
/// lose relative accuracy; see README.
[[nodiscard]] inline float256 hurwitz_zeta(float256 s, float256 a) {
    if (a == float256::one()) return riemann_zeta(s);  // before FnEnv: it keeps its own mode and flags
    detail_dd::FnEnv env;
    using detail_dd::dd;
    const float256 one = float256::one();
    if (s.isnan() || a.isnan()) return env.nan2(s, a);
    if (a.signbit() || a.iszero()) return env.invalid();
    if (s == one) return env.pole(false);
    if (s.isinf()) {
        if (s.signbit()) return env.invalid();
        if (a < one) return env.exact(float256::inf(false));  // a^{−s} dominates
        return env.exact(a == one ? one : float256::zero());
    }
    if (a.isinf()) return env.exact(s > one ? float256::zero() : -float256::inf(false));  // s < 1: ~ a^{1−s}/(s−1) → −∞
    if (s.iszero()) {  // ζ(0, a) = 1/2 − a, rounded once (exact when representable)
        const dd q = detail_dd::two_sum(float256::half(), -a);
        return env.round(q.hi, q.lo, 0, q.lo.iszero() ? 0 : 2);
    }
    // ζ(−2n, 1/2) = (2^{−2n} − 1) ζ(−2n) = 0: the only zeros in a at a binary256 point
    // (the rational roots of Bernoulli polynomials are 0, 1/2 and 1).
    if (s.signbit() && a == float256::half() && trunc(ldexp(s, -1)) == ldexp(s, -1)) return env.exact(float256::zero());
    double cancel = 0;
    detail_sf::SD r;
    if (s <= float256(-detail_sf::zeta::kReflS) &&
        !(a > one && detail_sf::zeta::em_plan(s.to_double(), a.to_double(), -std::log2(1 - s.to_double())) > 0))  // EM at N = 0 (no cancellation) when it converges
        r = detail_sf::zeta::hurwitz_refl(s, a, &cancel);
    else
        r = detail_sf::zeta::hurwitz(dd(s), detail_dd::two_sum(s, -one), a, &cancel);
    // next to a zero in a (s < 1) the dd sum cancels: redo it in multiprecision
    if (s < one && cancel > 140) r = detail_sf::zeta::zmp::hurwitz_mp(s, a, cancel);
    return detail_sf::out(env, r);
}

namespace detail_sf {
namespace igamma {

// Σ_{n>=0} x^n / (a(a+1)…(a+n)) (γ(a,x) = x^a e^{-x} · this), for x < a + 1.
// Positive terms; once they decrease and fall below 2^-245 of the sum the rest
// runs in float256.
[[nodiscard]] inline dd series(const dd& a, const dd& x) {
    const dd one{float256::one()};
    dd term = one / a, sum = term;
    int n = 1;
    for (; n < 100000000; ++n) {
        term = term * x / (a + D(n));
        sum = sum + term;
        if (small_rel(term, sum)) return sum;
        if ((a + D(n)).hi > x.hi && f256_tail(term, sum)) break;
    }
    float256 t = term.hi, tail{};
    const int lim = ilogb(sum.hi) - 492;
    for (++n; n < 100000000; ++n) {
        t = t * x.hi / (a.hi + float256::from_int(n));
        tail = tail + t;
        if (t.iszero() || ilogb(t) < lim) break;
    }
    return sum + dd(tail);
}
// Continued fraction h = 1/(b_1 + a_2/(b_2 + …)), b_i = x + 2i − 1 − a,
// a_i = −(i−1)(i−1−a), with Γ(a,x) = x^a e^{-x} · h, for x >= a + 1. Wallis
// recurrence (no divisions) with power-of-two rescaling, tested every 8 steps.
[[nodiscard]] inline dd cf(const dd& a, const dd& x) {
    const dd one{float256::one()};
    dd b = x + one - a;
    dd A0 = one, A1 = b, B0{float256::zero()}, B1 = one;
    dd h = one / b, hprev{float256::zero()};
    for (int i = 1; i < 100000000; ++i) {
        const dd an = -(D(i) * (D(i) - a));
        b = b + D(2);
        const dd A2 = b * A1 + an * A0, B2 = b * B1 + an * B0;
        A0 = A1; A1 = A2; B0 = B1; B1 = B2;
        const int e = ilogb(A1.hi);
        if (e > 2000 || e < -2000) {
            A0 = detail_dd::ldexp(A0, -e); A1 = detail_dd::ldexp(A1, -e);
            B0 = detail_dd::ldexp(B0, -e); B1 = detail_dd::ldexp(B1, -e);
        }
        if ((i & 7) == 0) {
            h = B1 / A1;
            if (near_one(h / hprev)) break;
            hprev = h;
        }
    }
    return h;
}
// ln Γ(1 + a) for 0 < a < 1, accurate relative to its size near a = 0.
[[nodiscard]] inline dd lgamma1p(float256 a) {
    if (ilogb(a) < -100) return detail_gamma::lgamma_near_1_2(a, false);
    return detail_gamma::lgamma_pos_dd(detail_dd::two_sum(float256::one(), a));
}
// Γ(a) as m·2^k, a > 0.
[[nodiscard]] inline SD gam(float256 a) {
    if (a < ldexp(float256::one(), 20)) return gamma_sd(dd(a));
    return exp_sd(detail_gamma::lgamma_pos_dd(dd(a)));
}

// a − b for m·2^k values (the smaller one shifted; far below the larger it is dropped).
[[nodiscard]] inline SD sd_sub(const SD& a, const SD& b) {
    if (b.m.hi.iszero()) return a;
    if (a.m.hi.iszero()) return {-b.m, b.k};
    if (a.k >= b.k) {
        if (a.k - b.k > 1200) return a;
        return norm({a.m - detail_dd::ldexp(b.m, b.k - a.k), a.k});
    }
    if (b.k - a.k > 1200) return {-b.m, b.k};
    return norm({detail_dd::ldexp(a.m, a.k - b.k) - b.m, b.k});
}
// Γ(a, x) for a <= 0, 0 < x < 36, a >= −1000 (as m·2^k). With N = round(−a),
// e = a + N ∈ [−1/2, 1/2] (exact):
//   Γ(a, x) = (−1)^N x^e expm1(L)/(e N!) − x^a Σ_{n>=0, n≠N} (−x)^n/(n!(a+n)),
//   L = lnΓ(1+e) − Σ_{k=1}^N ln(1 − e/k) − e ln x,
// which is Γ(a) − γ(a, x) with the two poles at a = −N cancelled analytically
// (e = 0, a negative integer or zero: expm1(L)/e → ψ(N+1) − ln x). The rest
// cancels at most ~x^N e^x (N+x)/N! < 2^110 of the ~2^-474 working precision.
[[nodiscard]] inline SD upper_nonpos_series(float256 a, const dd& X, const dd& lx) {
    const float256 one = float256::one();
    const float256 nf = floor(-a + float256::half());
    const int N = static_cast<int>(nf.to_double());
    const float256 e = a + nf;  // exact (Sterbenz)
    dd fact{one};
    for (int k = 2; k <= N; ++k) fact = fact * D(k);
    dd q;  // expm1(L)/e
    if (e.iszero()) {
        dd h{float256::zero()};
        for (int k = N; k >= 1; --k) h = h + dd(one) / D(k);
        q = h - euler_gamma() - lx;
    } else {
        const dd E{e};
        const dd lg = ilogb(e) < -100 ? detail_gamma::lgamma_near_1_2(e, false)
                                      : detail_gamma::lgamma_pos_dd(detail_dd::two_sum(one, e));
        dd sl{float256::zero()};  // Σ ln(1 − e/k)
        if (ilogb(e) < -60) {
            // −Σ_j e^j H_N^(j)/j, j <= 9 (e^10 < 2^-600)
            dd hj[10]{};
            for (int k = N; k >= 1; --k) {
                const dd r = dd(one) / D(k);
                dd rp = r;
                for (int j = 1; j <= 9; ++j) { hj[j] = hj[j] + rp; rp = rp * r; }
            }
            dd ej = E;
            for (int j = 1; j <= 9 && N > 0; ++j) {
                sl = sl - ej * hj[j] / D(j);
                ej = ej * E;
            }
        } else if (N > 0) {
            dd p{one};
            for (int k = 1; k <= N; ++k) p = p * (dd(one) - E / D(k));
            sl = detail_dd::log(p);
        }
        const dd L = lg - sl - E * lx;
        q = detail_dd::expm1(L) / E * detail_dd::exp(E * lx);
    }
    dd t1 = q / fact;
    if (N & 1) t1 = -t1;
    // S = Σ_{n≠N} (−x)^n/(n!(a+n)); stop once n > N, n > 2x and the term is
    // 2^-620 below the largest partial term.
    const dd A{a};
    dd term{one}, S{float256::zero()};
    int big = -1000000;
    for (int n = 0; n < 100000; ++n) {
        if (n > 0) term = -(term * X) / D(n);
        if (n != N) {
            const dd tn = term / (A + D(n));
            S = S + tn;
            if (!tn.hi.iszero()) {
                const int lt = ilogb(tn.hi);
                if (lt > big) big = lt;
                if (n > N && dd(D(n)).hi > X.hi + X.hi && lt < big - 620) break;
            } else if (n > N) break;
        }
    }
    const SD xa = exp_sd(A * lx);
    return sd_sub(norm({t1, 0}), mul(xa, S));
}
// Upper incomplete gamma Γ(a, x) for a <= 0 (finite or −0), x not NaN.
[[nodiscard]] inline float256 upper_nonpos(float256 a, float256 x) {
    detail_dd::FnEnv env;
    if (a.isinf()) return env.invalid();
    if (x.signbit() && !x.iszero()) return env.invalid();
    if (x.iszero()) return env.pole(false);      // ∫_0 t^{a−1} diverges
    if (x.isinf()) return env.exact(float256::zero());
    const float256 one = float256::one();
    const float256 two1000 = ldexp(one, 1000);
    if (x > two1000) return env.underflow(false);  // < e^{−x}
    const dd X{x};
    if (-a > two1000) {
        // x^a overflows (x < 1) or underflows (x > 1) by more than 2^(2^763);
        // at x = 1, Γ(a, 1) = e^{−1}/(1 − a)·(1 + O(a^-2)).
        if (x < one) return env.overflow(false);
        if (x > one) return env.underflow(false);
        return out(env, mul(exp_sd(-X), dd(one) / (X - dd(a))));
    }
    const dd lx = detail_dd::log(X);
    if (x < float256(36) && -a <= float256(1000)) return out(env, upper_nonpos_series(a, X, lx));
    // Legendre continued fraction; for a <= 0 every partial numerator is
    // negative and it converges fast once x >= 36 or −a > 1000.
    const dd A{a};
    return out(env, mul(exp_sd(A * lx - X), cf(A, X)));
}

// kind: 0 lower γ, 1 upper Γ, 2 P, 3 Q.
[[nodiscard]] inline float256 eval(float256 a, float256 x, int kind) {
    // γ(a, ∞) = Γ(a, 0) = Γ(a): tgamma's own rounding (exact for small integers).
    if (kind < 2 && !a.isnan() && !a.signbit() && !a.iszero() && !a.isinf() &&
        ((kind == 0 && x.isinf() && !x.signbit()) || (kind == 1 && x.iszero())))
        return base256::tgamma(a);
    if (kind == 1 && !a.isnan() && !x.isnan() && (a.signbit() || a.iszero())) return upper_nonpos(a, x);
    detail_dd::FnEnv env;
    const float256 one = float256::one();
    if (a.isnan() || x.isnan()) return env.nan2(a, x);
    if (a.signbit() || a.iszero() || a.isinf()) return env.invalid();
    if (x.signbit() && !x.iszero()) return env.invalid();
    const bool regularized = kind >= 2;
    const bool upper = kind == 1 || kind == 3;
    if (x.iszero() || x.isinf()) {
        const bool full = upper == x.iszero();  // the integral covers (0, ∞)
        if (!full) return env.exact(float256::zero());
        return env.exact(one);  // regularized (the unregularized full case returned above)
    }
    const dd A{a}, X{x};
    const dd lx = detail_dd::log(X);
    const SD pre = exp_sd(A * lx - X);  // x^a e^{-x}
    // The series for γ also runs for a + 1 <= x < 36, where the fraction for Γ
    // converges slowly at this precision; the complement 1 − P then cancels at
    // most ~e^36 (< 2^52) of the ~2^-474 working precision.
    const bool small_x = x < a + one;
    const bool use_series = small_x || x < float256(36);
    // Small a, x < a + 1, upper part: Γ(a,x) = Γ(a)·(−expm1(a ln x − lnΓ(1+a))) + x^a Σ_{n>=1} (−1)^{n+1} x^n/(n!(a+n)).
    if (upper && small_x && a < one) {
        const dd e = A * lx - lgamma1p(a);
        const dd t1 = -detail_dd::expm1(e);  // Q part
        dd term = dd(one), sum{float256::zero()};
        for (int n = 1; n < 100000000; ++n) {
            term = -(term * X) / D(n);
            const dd tn = term / (A + D(n));
            sum = sum - tn;
            if (small_rel(tn, sum)) break;
        }
        // x^a Σ / Γ(a) = sum · exp(a ln x − lnΓ(a)), lnΓ(a) = lnΓ(1+a) − ln a
        const SD r2 = exp_sd(A * lx - lgamma1p(a) + detail_dd::log(A));
        const dd q = t1 + detail_dd::ldexp(r2.m * sum, r2.k);
        if (regularized) return q.out(env);
        return out(env, mul(gam(a), q));
    }
    SD part;  // lower (series) or upper (fraction), unregularized
    if (use_series) part = mul(pre, series(A, X));
    else part = mul(pre, cf(A, X));
    const bool direct = use_series != upper;  // the part is what was asked for
    const SD g = gam(a);
    if (regularized) {
        const SD r = norm({part.m / g.m, part.k - g.k});
        if (direct) return out(env, r);
        if (r.k < -400) return env.round(one, float256::zero(), 0, -1);
        return (dd(one) - detail_dd::ldexp(r.m, r.k)).out(env);
    }
    if (direct) return out(env, part);
    // Γ(a) − part, as m·2^k with Γ(a)'s exponent
    const int dk = part.k - g.k;
    if (dk < -400) return out(env, g, -1);
    return out(env, norm({g.m - detail_dd::ldexp(part.m, dk), g.k}));
}

} // namespace igamma
} // namespace detail_sf

/// Lower incomplete gamma γ(a, x) = ∫_0^x t^{a−1} e^{−t} dt, a > 0, x >= 0.
[[nodiscard]] inline float256 tgamma_lower(float256 a, float256 x) { return detail_sf::igamma::eval(a, x, 0); }
/// Upper incomplete gamma Γ(a, x) = ∫_x^∞ t^{a−1} e^{−t} dt, x >= 0; any a > 0, and
/// a <= 0 (including 0 and negative integers) for x > 0. For a <= 0, Γ(a, 0) = +∞
/// (divide-by-zero) and a = −∞ is invalid.
[[nodiscard]] inline float256 tgamma_upper(float256 a, float256 x) { return detail_sf::igamma::eval(a, x, 1); }
/// Regularized lower incomplete gamma P(a, x) = γ(a, x)/Γ(a).
[[nodiscard]] inline float256 gamma_p(float256 a, float256 x) { return detail_sf::igamma::eval(a, x, 2); }
/// Regularized upper incomplete gamma Q(a, x) = Γ(a, x)/Γ(a) = 1 − P(a, x).
[[nodiscard]] inline float256 gamma_q(float256 a, float256 x) { return detail_sf::igamma::eval(a, x, 3); }

namespace detail_sf {
namespace ell {

using detail_dd::sqrt;
// max |A − v|/|A| below 2^-41 (|Z_i| < 2^-40): the degree-12 tails of
// carlson_tail are then < 2^-500.
[[nodiscard]] inline bool close(const dd& A, const dd& a, const dd& b, const dd& c, const dd& d) {
    const int ea = ilogb(A.hi) - 41;
    for (const dd* v : {&a, &b, &c, &d}) {
        const dd e = A - *v;
        if (!e.hi.iszero() && ilogb(e.hi) >= ea) return false;
    }
    return true;
}

// Carlson's R_{−a}(b; z) after duplication, from its Taylor series about the
// weighted mean A (DLMF 19.19.7, 19.36): A^{−a} Σ_N (a)_N/(c)_N T_N, with
// T_N = [t^N] Π_i (1 − Z_i t)^{−b_i}, Z_i = 1 − z_i/A, c = Σ b_i. This returns
// the sum (without A^{−a}) through degree kTailDeg; with max |Z_i| < 2^-40
// (`close`) the dropped terms are below 2^-500. Degrees >= 7 (< 2^-280) run in
// float256. Arguments are doubled: b2[i] = 2 b_i, a2 = 2a, c2 = 2c.
inline constexpr int kTailDeg = 12;
inline constexpr int kTailDD = 6;
[[nodiscard]] inline dd carlson_tail(const dd* Z, const int* b2, int nv, int a2, int c2) {
    constexpr int D = kTailDeg;
    dd prod[D + 1];
    float256 prodf[D + 1];
    for (int i = 0; i < nv; ++i) {
        // (1 − Z t)^{−b}: s_k = s_{k−1} Z (2b + 2k − 2)/(2k)
        dd ser[D + 1];
        ser[0] = dd(float256::one());
        for (int k = 1; k <= D; ++k) {
            if (k <= kTailDD) ser[k] = ser[k - 1] * Z[i] * float256::from_int(b2[i] + 2 * k - 2) / float256::from_int(2 * k);
            else ser[k] = dd(ser[k - 1].hi * Z[i].hi * float256::from_int(b2[i] + 2 * k - 2) / float256::from_int(2 * k));
        }
        if (i == 0) {
            for (int k = 0; k <= D; ++k) { prod[k] = ser[k]; prodf[k] = ser[k].hi; }
            continue;
        }
        for (int n = D; n >= 0; --n) {
            if (n <= kTailDD) {
                dd acc = prod[n];
                for (int k = 1; k <= n; ++k) acc = acc + prod[n - k] * ser[k];
                prod[n] = acc;
                prodf[n] = acc.hi;
            } else {
                float256 acc = prodf[n];
                for (int k = 1; k <= n; ++k) acc = acc + prodf[n - k] * ser[k].hi;
                prodf[n] = acc;
            }
        }
    }
    // Σ r_N T_N, r_N = (a)_N/(c)_N
    float256 rf = float256::one(), tailf{};
    dd r{float256::one()}, sum = prod[0];
    for (int n = 1; n <= D; ++n) {
        if (n <= kTailDD) {
            r = r * float256::from_int(a2 + 2 * n - 2) / float256::from_int(c2 + 2 * n - 2);
            rf = r.hi;
            if (n >= 2) sum = sum + r * prod[n];  // T_1 = Σ b_i Z_i = 0
        } else {
            rf = rf * float256::from_int(a2 + 2 * n - 2) / float256::from_int(c2 + 2 * n - 2);
            tailf = tailf + rf * prodf[n];
        }
    }
    return sum + dd(tailf);
}

// R_F(x, y, z), x, y, z >= 0, at most one zero.
[[nodiscard]] inline dd rf(dd x, dd y, dd z) {
    const dd one{float256::one()};
    dd A = (x + y + z) / float256(3);
    for (int it = 0; it < 400 && !close(A, x, y, z, z); ++it) {
        const dd sx = sqrt(x), sy = sqrt(y), sz = sqrt(z);
        const dd lam = sx * sy + sy * sz + sz * sx;
        x = detail_dd::ldexp(x + lam, -2); y = detail_dd::ldexp(y + lam, -2); z = detail_dd::ldexp(z + lam, -2);
        A = detail_dd::ldexp(A + lam, -2);
    }
    const dd iA = one / A;
    const dd Zs[3] = {one - x * iA, one - y * iA, one - z * iA};
    const int b2[3] = {1, 1, 1};
    return carlson_tail(Zs, b2, 3, 1, 3) / sqrt(A);
}

// R_D(x, y, z), x, y >= 0 (at most one zero), z > 0.
[[nodiscard]] inline dd rd(dd x, dd y, dd z) {
    const dd one{float256::one()};
    dd A = (x + y + z * float256(3)) / float256(5);
    dd sum{float256::zero()};
    int m = 0;
    for (int it = 0; it < 400 && !close(A, x, y, z, z); ++it, ++m) {
        const dd sx = sqrt(x), sy = sqrt(y), sz = sqrt(z);
        const dd lam = sx * sy + sy * sz + sz * sx;
        sum = sum + detail_dd::ldexp(one / (sz * (z + lam)), -2 * m);
        x = detail_dd::ldexp(x + lam, -2); y = detail_dd::ldexp(y + lam, -2); z = detail_dd::ldexp(z + lam, -2);
        A = detail_dd::ldexp(A + lam, -2);
    }
    const dd iA = one / A;
    const dd Zs[3] = {one - x * iA, one - y * iA, one - z * iA};
    const int b2[3] = {1, 1, 3};
    const dd poly = carlson_tail(Zs, b2, 3, 3, 5);
    return sum * float256(3) + detail_dd::ldexp(poly / (A * sqrt(A)), -2 * m);
}

// R_C(x, y) = ∫ dt / (2 √(t+x) (t+y)), x >= 0, y > 0.
[[nodiscard]] inline dd rc(const dd& x, const dd& y) {
    const dd one{float256::one()};
    if (x.hi.iszero()) return detail_dd::pi() / (sqrt(y) * float256::two());
    const dd t = (y - x) / x;  // R_C = atan(√t)/√t / √x  (atanh form for t < 0)
    dd g;
    if (t.hi.iszero() || ilogb(t.hi) < -8) {
        // Σ (−t)^n/(2n+1): terms below 2^-492 dropped, those below 2^-245 in float256.
        const auto& io = detail_dd::inv_odd();
        if (t.hi.iszero()) {
            g = one;
        } else {
            const int lt = -(ilogb(t.hi) + 1);  // |t| < 2^-lt, lt >= 8
            int n = 492 / lt + 1, nf = 245 / lt + 1;
            if (n > detail_dd::kMaxOddInv) n = detail_dd::kMaxOddInv;
            if (nf > n) nf = n;
            float256 pf{};
            for (int j = n; j > nf; --j) pf = io[static_cast<std::size_t>(j)].hi - t.hi * pf;
            dd p{pf};
            for (int j = nf; j >= 0; --j) p = io[static_cast<std::size_t>(j)] - t * p;
            g = p;
        }
    } else if (!t.hi.signbit()) {
        const dd st = sqrt(t);
        g = detail_math::atan_dd(st) / st;
    } else {
        const dd u = sqrt(-t);  // atanh(u)/u, 0 < u < 1
        g = detail_dd::log((one + u) / (one - u)) / (u * float256::two());
    }
    return g / sqrt(x);
}

// R_J(x, y, z, p), x, y, z >= 0 (at most one zero), p > 0.
[[nodiscard]] inline dd rj(dd x, dd y, dd z, dd p) {
    const dd one{float256::one()};
    dd A = (x + y + z + p * float256::two()) / float256(5);
    dd sum{float256::zero()};
    int m = 0;
    for (int it = 0; it < 400 && !close(A, x, y, z, p); ++it, ++m) {
        const dd sx = sqrt(x), sy = sqrt(y), sz = sqrt(z);
        const dd lam = sx * sy + sy * sz + sz * sx;
        const dd al0 = p * (sx + sy + sz) + sx * sy * sz;
        const dd pl = p + lam;
        sum = sum + detail_dd::ldexp(rc(al0 * al0, p * pl * pl), -2 * m);
        x = detail_dd::ldexp(x + lam, -2); y = detail_dd::ldexp(y + lam, -2);
        z = detail_dd::ldexp(z + lam, -2); p = detail_dd::ldexp(p + lam, -2);
        A = detail_dd::ldexp(A + lam, -2);
    }
    const dd iA = one / A;
    const dd Zs[4] = {one - x * iA, one - y * iA, one - z * iA, one - p * iA};
    const int b2[4] = {1, 1, 1, 2};
    const dd poly = carlson_tail(Zs, b2, 4, 3, 5);
    return sum * float256(3) + detail_dd::ldexp(poly / (A * sqrt(A)), -2 * m);
}

// 1 − k² = (1 − k)(1 + k), exact inputs, as dd.
[[nodiscard]] inline dd kc2(float256 k) {
    const dd one{float256::one()};
    return (one - dd(k)) * (one + dd(k));
}
// K(k) and E(k) by the AGM (|k| < 1).
inline void agm_ke(float256 k, dd& K, dd* E) {
    dd a{float256::one()}, b = sqrt(kc2(k));
    dd csum = detail_dd::ldexp(detail_dd::two_prod(k, k), -1);  // 2^{-1} c_0²
    for (int n = 1; n < 100; ++n) {
        const dd c = detail_dd::ldexp(a - b, -1);
        const dd an = detail_dd::ldexp(a + b, -1);
        b = sqrt(a * b);
        a = an;
        csum = csum + detail_dd::ldexp(c * c, n - 1);
        if (c.hi.iszero() || ilogb(c.hi) < ilogb(a.hi) - 250) break;  // c_{n+1} ~ c_n² : next term < 2^-500
    }
    K = detail_dd::pio2() / a;
    if (E) *E = K * (dd(float256::one()) - csum);
}

// Bulirsch's cel(k_c, p, 1, 1) = Π(1 − p, k) for p > 0 (quadratically
// convergent Gauss transformation; NR's cel with a = b = 1). Stops when
// |g − q_c| <= 2^-240 g, after which the error is ~2^-480.
[[nodiscard]] inline dd cel(const dd& kc2v, const dd& p0) {
    const dd one{float256::one()};
    dd qc = sqrt(kc2v);
    dd a = one, p = sqrt(p0);
    dd b = one / p;
    dd e = qc, em = one;
    for (int it = 0; it < 100; ++it) {
        const dd f = a;
        a = a + b / p;
        dd g = e / p;
        b = b + f * g;
        b = b + b;
        p = g + p;
        g = em;
        em = em + qc;
        const dd d = g - qc;
        if (d.hi.iszero() || ilogb(d.hi) < ilogb(g.hi) - 240) break;
        qc = sqrt(e);
        qc = qc + qc;
        e = qc * em;
    }
    return detail_dd::pio2() * (b + a * em) / (em * (em + p));
}

} // namespace ell
} // namespace detail_sf

/// Complete elliptic integral of the first kind K(k), |k| < 1 (|k| = 1: +∞, divide-by-zero).
[[nodiscard]] inline float256 comp_ellint_1(float256 k) {
    detail_dd::FnEnv env;
    if (k.isnan()) return env.nan(k);
    if (abs(k) > float256::one()) return env.invalid();
    if (abs(k) == float256::one()) return env.pole(false);
    detail_dd::dd K;
    detail_sf::ell::agm_ke(k, K, nullptr);
    return K.out(env);
}
/// Complete elliptic integral of the second kind E(k), |k| <= 1.
[[nodiscard]] inline float256 comp_ellint_2(float256 k) {
    detail_dd::FnEnv env;
    if (k.isnan()) return env.nan(k);
    if (abs(k) > float256::one()) return env.invalid();
    if (abs(k) == float256::one()) return env.exact(float256::one());
    detail_dd::dd K, E;
    detail_sf::ell::agm_ke(k, K, &E);
    return E.out(env);
}
/// Complete elliptic integral of the third kind Π(ν, k) = ∫_0^{π/2} dθ/((1 − ν sin²θ)√(1 − k² sin²θ)),
/// |k| < 1, ν < 1 (ν = 1 or |k| = 1: +∞, divide-by-zero; ν > 1: domain error).
[[nodiscard]] inline float256 comp_ellint_3(float256 k, float256 nu) {
    detail_dd::FnEnv env;
    using detail_dd::dd;
    if (k.isnan() || nu.isnan()) return env.nan2(k, nu);
    if (abs(k) > float256::one() || nu > float256::one()) return env.invalid();
    if (nu.isinf()) return env.exact(float256::zero());  // ν → −∞: Π → 0
    if (abs(k) == float256::one() || nu == float256::one()) return env.pole(false);
    if (nu.iszero()) { dd K; detail_sf::ell::agm_ke(k, K, nullptr); return K.out(env); }
    return detail_sf::ell::cel(detail_sf::ell::kc2(k), detail_dd::two_sum(float256::one(), -nu)).out(env);
}

namespace detail_sf {
namespace ell {
// kind 1: F, 2: E, 3: Π. φ = mπ + r, |r| <= π/2: X(φ) = 2m X_complete + X(r).
[[nodiscard]] inline float256 incomplete(int kind, float256 k, float256 nu, float256 phi) {
    detail_dd::FnEnv env;
    const float256 one = float256::one();
    if (k.isnan() || phi.isnan() || nu.isnan()) {
        if (nu.isnan()) return env.nan(nu);
        return env.nan2(k, phi);
    }
    if (abs(k) > one) return env.invalid();
    if (phi.isinf()) return (kind == 3 && nu.isinf()) ? env.invalid() : env.exact(phi);  // grows like φ
    if (nu.isinf()) return env.invalid();
    if (phi.iszero()) return env.exact(phi);
    if (k.iszero() && (kind == 1 || kind == 2 || nu.iszero())) return env.exact(phi);  // F(φ,0) = E(φ,0) = φ
    // tiny φ: X ≈ φ + c φ³, c = k²/6 (F), −k²/6 (E), (2ν + k²)/6 (Π)
    if (ilogb(phi) < -250) {
        int c = kind == 2 ? -1 : 1;
        if (kind == 3) {
            const dd t = dd(nu) * float256::two() + detail_dd::two_prod(k, k);
            c = t.hi.iszero() ? 0 : (t.hi.signbit() ? -1 : 1);
        }
        const int tail = c == 0 ? 2 : (phi.signbit() ? -c : c);
        return env.round(phi, float256::zero(), 0, tail);
    }
    const dd pi = detail_dd::pi();
    const float256 mf = round(phi / pi.hi);
    const dd r = dd(phi) - pi * mf;
    dd s, c;
    detail_dd::sincos_reduced(r, s, c);
    const dd s2 = s * s, c2 = c * c;
    const dd k2 = detail_dd::two_prod(k, k);
    const dd delta = c2 + kc2(k) * s2;  // 1 − k² sin² r
    const dd one_d{one};
    dd part;
    if (kind == 3) {
        const dd p = one_d - dd(nu) * s2;
        if (p.hi.signbit() || p.hi.iszero()) return env.invalid();  // 1 − ν sin²φ <= 0: not supported
        part = s * rf(c2, delta, one_d) + dd(nu) * s * s2 * rj(c2, delta, one_d, p) / float256(3);
    } else if (abs(k) == one && mf.iszero() && c2.hi.iszero()) {
        return env.pole(phi.signbit());
    } else {
        part = s * rf(c2, delta, one_d);
        if (kind == 2) part = part - k2 * s * s2 * rd(c2, delta, one_d) / float256(3);
    }
    if (!mf.iszero()) {
        dd full;
        if (kind == 1) {
            if (abs(k) == one) return env.pole(phi.signbit());
            agm_ke(k, full, nullptr);
        } else if (kind == 2) {
            if (abs(k) == one) full = one_d;
            else { dd K; agm_ke(k, K, &full); }
        } else {
            if (abs(k) == one || nu >= one) return nu > one ? env.invalid() : env.pole(phi.signbit());
            full = cel(kc2(k), detail_dd::two_sum(one, -nu));
        }
        part = part + full * (mf + mf);
    }
    return part.out(env);
}
} // namespace ell
} // namespace detail_sf

/// Incomplete elliptic integral of the first kind F(k, φ) = ∫_0^φ dθ/√(1 − k² sin²θ), |k| <= 1.
[[nodiscard]] inline float256 ellint_1(float256 k, float256 phi) { return detail_sf::ell::incomplete(1, k, float256::zero(), phi); }
/// Incomplete elliptic integral of the second kind E(k, φ) = ∫_0^φ √(1 − k² sin²θ) dθ, |k| <= 1.
[[nodiscard]] inline float256 ellint_2(float256 k, float256 phi) { return detail_sf::ell::incomplete(2, k, float256::zero(), phi); }
/// Incomplete elliptic integral of the third kind Π(ν; φ, k) = ∫_0^φ dθ/((1 − ν sin²θ)√(1 − k² sin²θ)),
/// |k| <= 1, 1 − ν sin²φ > 0 (and ν < 1 when |φ| > π/2); otherwise a domain error.
[[nodiscard]] inline float256 ellint_3(float256 k, float256 nu, float256 phi) { return detail_sf::ell::incomplete(3, k, nu, phi); }


// ---- Per-call rounding (README "Per-call rounding"): f(args..., m) is f(args...)
// evaluated with the rounding mode m; the thread's mode is restored afterwards.
[[nodiscard]] inline float256 cyl_bessel_j(float256 nu, float256 x, round_mode m) { return detail::with_round(m, [&] { return cyl_bessel_j(nu, x); }); }
[[nodiscard]] inline float256 cyl_neumann(float256 nu, float256 x, round_mode m) { return detail::with_round(m, [&] { return cyl_neumann(nu, x); }); }
[[nodiscard]] inline float256 cyl_bessel_i(float256 nu, float256 x, round_mode m) { return detail::with_round(m, [&] { return cyl_bessel_i(nu, x); }); }
[[nodiscard]] inline float256 cyl_bessel_k(float256 nu, float256 x, round_mode m) { return detail::with_round(m, [&] { return cyl_bessel_k(nu, x); }); }
[[nodiscard]] inline float256 sph_bessel(unsigned n, float256 x, round_mode m) { return detail::with_round(m, [&] { return sph_bessel(n, x); }); }
[[nodiscard]] inline float256 sph_neumann(unsigned n, float256 x, round_mode m) { return detail::with_round(m, [&] { return sph_neumann(n, x); }); }
[[nodiscard]] inline float256 jn(int n, float256 x, round_mode m) { return detail::with_round(m, [&] { return jn(n, x); }); }
[[nodiscard]] inline float256 yn(int n, float256 x, round_mode m) { return detail::with_round(m, [&] { return yn(n, x); }); }
[[nodiscard]] inline float256 beta(float256 x, float256 y, round_mode m) { return detail::with_round(m, [&] { return beta(x, y); }); }
[[nodiscard]] inline float256 ibeta(float256 a, float256 b, float256 x, round_mode m) { return detail::with_round(m, [&] { return ibeta(a, b, x); }); }
[[nodiscard]] inline float256 ibetac(float256 a, float256 b, float256 x, round_mode m) { return detail::with_round(m, [&] { return ibetac(a, b, x); }); }
[[nodiscard]] inline float256 riemann_zeta(float256 x, round_mode m) { return detail::with_round(m, [&] { return riemann_zeta(x); }); }
[[nodiscard]] inline float256 hurwitz_zeta(float256 s, float256 a, round_mode m) { return detail::with_round(m, [&] { return hurwitz_zeta(s, a); }); }
[[nodiscard]] inline float256 tgamma_lower(float256 a, float256 x, round_mode m) { return detail::with_round(m, [&] { return tgamma_lower(a, x); }); }
[[nodiscard]] inline float256 tgamma_upper(float256 a, float256 x, round_mode m) { return detail::with_round(m, [&] { return tgamma_upper(a, x); }); }
[[nodiscard]] inline float256 gamma_p(float256 a, float256 x, round_mode m) { return detail::with_round(m, [&] { return gamma_p(a, x); }); }
[[nodiscard]] inline float256 gamma_q(float256 a, float256 x, round_mode m) { return detail::with_round(m, [&] { return gamma_q(a, x); }); }
[[nodiscard]] inline float256 comp_ellint_1(float256 k, round_mode m) { return detail::with_round(m, [&] { return comp_ellint_1(k); }); }
[[nodiscard]] inline float256 comp_ellint_2(float256 k, round_mode m) { return detail::with_round(m, [&] { return comp_ellint_2(k); }); }
[[nodiscard]] inline float256 comp_ellint_3(float256 k, float256 nu, round_mode m) { return detail::with_round(m, [&] { return comp_ellint_3(k, nu); }); }
[[nodiscard]] inline float256 ellint_1(float256 k, float256 phi, round_mode m) { return detail::with_round(m, [&] { return ellint_1(k, phi); }); }
[[nodiscard]] inline float256 ellint_2(float256 k, float256 phi, round_mode m) { return detail::with_round(m, [&] { return ellint_2(k, phi); }); }
[[nodiscard]] inline float256 ellint_3(float256 k, float256 nu, float256 phi, round_mode m) { return detail::with_round(m, [&] { return ellint_3(k, nu, phi); }); }
} // namespace base256
