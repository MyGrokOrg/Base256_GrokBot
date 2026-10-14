#pragma once
/// Special functions (C++17 <cmath> special-math names): Bessel functions,
/// beta and the regularized incomplete beta.
/// Every function evaluates in double-float256 (~474 bits) and rounds once in
/// the caller's rounding mode, raising IEEE flags like the elementary
/// functions. See README "Special functions" for measured accuracy and limits.

#include <cmath>
#include <cstdint>

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
// Γ(z) as m·2^k for a dd z > 0.
[[nodiscard]] inline SD gamma_sd(const dd& z) {
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
    float256 r = fmod(t, float256::two());  // exact, |r| < 2
    if (r.signbit()) r = r + float256::two();  // exact for |r| >= 2^-236·...; tiny r handled by sign
    // r in [0, 2): quadrant q = floor(2r), f = r - q/2 in [0, 1/2)
    const float256 q = floor(r + r);
    const float256 f = r - q * float256::half();
    dd sf, cf;
    detail_dd::sincos_reduced(detail_dd::pi() * f, sf, cf);
    const int qi = static_cast<int>(q.to_double()) & 3;
    switch (qi) {
        case 0: s = sf; c = cf; break;
        case 1: s = cf; c = -sf; break;
        case 2: s = -sf; c = -cf; break;
        default: s = -cf; c = sf; break;
    }
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
        const auto& f = detail_dd::inv_factorials();
        const dd e2 = e * e;
        dd a = f[100], b = f[101];
        for (int k = 49; k >= 0; --k) {
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
    for (int k = 1; k < 2000000; ++k) {
        term = term * qs / (dd(float256::from_int(k)) * (np1 + D(k - 1)));
        sum = sum + term;
        if (small_rel(term, sum) && float256::from_int(k) * float256::from_int(k) > abs(q.hi)) break;
    }
    return mul(pre, sum);
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

// gam1 = (1/Γ(1−μ) − 1/Γ(1+μ))/(2μ), gam2 = (1/Γ(1−μ) + 1/Γ(1+μ))/2, plus 1/Γ(1±μ); |μ| <= 1/2.
inline void temme_gammas(float256 mu, dd& gam1, dd& gam2, dd& gampl, dd& gammi) {
    const dd one{float256::one()};
    if (mu.iszero()) { gam1 = -euler_gamma(); gam2 = one; gampl = gammi = one; return; }
    const SD a = rgamma_sd(detail_dd::two_sum(float256::one(), mu));
    const SD b = rgamma_sd(detail_dd::two_sum(float256::one(), -mu));
    gampl = detail_dd::ldexp(a.m, a.k);
    gammi = detail_dd::ldexp(b.m, b.k);
    gam2 = detail_dd::ldexp(gammi + gampl, -1);
    if (ilogb(mu) < -100) {
        const dd c3 = detail_dd::dd_from_words(detail_dd::tables::kInvGammaC3);
        gam1 = -(euler_gamma() + c3 * detail_dd::two_prod(mu, mu));
    } else {
        gam1 = (gammi - gampl) / (dd(mu) * float256::two());
    }
}

// Y_μ(x), Y_{μ+1}(x) for 0 < x < 2, |μ| <= 1/2 (Temme's series).
inline void temme_y(float256 mu, float256 x, dd& ymu, dd& y1) {
    const dd one{float256::one()};
    const dd pi = detail_dd::pi();
    dd gam1, gam2, gampl, gammi;
    temme_gammas(mu, gam1, gam2, gampl, gammi);
    const float256 x2 = ldexp(x, -1);
    const dd fact = mu.iszero() ? one : pi * abs(mu) / detail_gamma::sin_pi_dd(abs(mu));
    const dd d = -detail_dd::log(dd(x2));
    const dd e = d * mu;
    dd ch, shc;
    cosh_sinhc(e, ch, shc);
    dd ff = fact * (gam1 * ch + gam2 * shc * d) * float256::two() / pi;
    const dd E = detail_dd::exp(e);
    dd p = E / (gampl * pi);
    dd q = one / (E * pi * gammi);
    const float256 hm = ldexp(abs(mu), -1);
    const dd fact3 = mu.iszero() ? one : detail_gamma::sin_pi_dd(hm) / (pi * hm);
    const dd r_s = pi * pi * ldexp(mu, -1) * fact3 * fact3;  // π·(πμ/2)·(sin(πμ/2)/(πμ/2))²
    dd c = one;
    const dd dq = -detail_dd::two_prod(x2, x2);
    dd sum = ff + r_s * q, sum1 = p;
    const dd mu2 = detail_dd::two_prod(mu, mu);
    for (int i = 1; i < 100000; ++i) {
        const dd di = D(i);
        ff = (di * ff + p + q) / (D(i) * D(i) - mu2);
        c = c * dq / di;
        p = p / (di - dd(mu));
        q = q / (di + dd(mu));
        const dd del = c * (ff + r_s * q);
        sum = sum + del;
        const dd del1 = c * p - di * del;
        sum1 = sum1 + del1;
        if (small_rel(del, sum + one) && small_rel(del1, sum1)) break;
    }
    ymu = -sum;
    y1 = -sum1 * float256::two() / dd(x);
}

// K_μ(x), K_{μ+1}(x) for 0 < x < 2, |μ| <= 1/2 (Temme's series).
inline void temme_k(float256 mu, float256 x, dd& kmu, dd& k1) {
    const dd one{float256::one()};
    const dd pi = detail_dd::pi();
    dd gam1, gam2, gampl, gammi;
    temme_gammas(mu, gam1, gam2, gampl, gammi);
    const float256 x2 = ldexp(x, -1);
    const dd fact = mu.iszero() ? one : pi * abs(mu) / detail_gamma::sin_pi_dd(abs(mu));
    const dd d = -detail_dd::log(dd(x2));
    const dd e = d * mu;
    dd ch, shc;
    cosh_sinhc(e, ch, shc);
    dd ff = fact * (gam1 * ch + gam2 * shc * d);
    dd sum = ff;
    const dd E = detail_dd::exp(e);
    dd p = detail_dd::ldexp(E / gampl, -1);
    dd q = detail_dd::ldexp(one / (E * gammi), -1);
    dd c = one;
    const dd dq = detail_dd::two_prod(x2, x2);
    dd sum1 = p;
    const dd mu2 = detail_dd::two_prod(mu, mu);
    for (int i = 1; i < 100000; ++i) {
        const dd di = D(i);
        ff = (di * ff + p + q) / (di * di - mu2);
        c = c * dq / di;
        p = p / (di - dd(mu));
        q = q / (di + dd(mu));
        const dd del = c * ff;
        sum = sum + del;
        const dd del1 = c * (p - di * ff);
        sum1 = sum1 + del1;
        if (small_rel(del, sum) && small_rel(del1, sum1)) break;
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
        c = -a * c / D(i + 1);
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

// Y_μ, Y_{μ+1} for x >= 2, |μ| <= 1/2 (Steed: CF1 for J'/J, CF2 for (J'+iY')/(J+iY)).
inline void steed_y(float256 mu, float256 x, dd& ymu, dd& y1) {
    const dd one{float256::one()};
    const dd X{x};
    const dd xi = one / X;
    const dd xi2 = xi * float256::two();
    const float256 tiny = ldexp(float256::one(), -20000);
    // CF1 (modified Lentz) for f = J'_μ/J_μ.
    int isign = 1;
    dd h = dd(mu) * xi;
    if (abs(h.hi) < tiny) h = dd(tiny);
    dd b = xi2 * mu;
    dd dd_{float256::zero()};
    dd c = h;
    for (int i = 1; i < 10000000; ++i) {
        b = b + xi2;
        dd_ = b - dd_;
        if (abs(dd_.hi) < tiny) dd_ = dd(tiny);
        c = b - one / c;
        if (abs(c.hi) < tiny) c = dd(tiny);
        dd_ = one / dd_;
        const dd del = c * dd_;
        h = del * h;
        if (dd_.hi.signbit()) isign = -isign;
        if (near_one(del)) break;
    }
    const dd f = h;
    // CF2.
    const dd mu2 = detail_dd::two_prod(mu, mu);
    dd a = dd(float256(0.25)) - mu2;
    dd p = -detail_dd::ldexp(xi, -1);
    dd q = one;
    const dd br = X * float256::two();
    dd bi = D(2);
    dd fact = a * xi / (p * p + q * q);
    dd cr = br + q * fact;
    dd ci = bi + p * fact;
    dd den = br * br + bi * bi;
    dd dr = br / den;
    dd di = -bi / den;
    dd dlr = cr * dr - ci * di;
    dd dli = cr * di + ci * dr;
    dd temp = p * dlr - q * dli;
    q = p * dli + q * dlr;
    p = temp;
    for (int i = 2; i < 1000000; ++i) {
        a = a + D(2 * (i - 1));
        bi = bi + D(2);
        dr = a * dr + br;
        di = a * di + bi;
        fact = a / (cr * cr + ci * ci);
        cr = br + cr * fact;
        ci = bi - ci * fact;
        den = dr * dr + di * di;
        dr = dr / den;
        di = -di / den;
        dlr = cr * dr - ci * di;
        dli = cr * di + ci * dr;
        temp = p * dlr - q * dli;
        q = p * dli + q * dlr;
        p = temp;
        if (near_one(dlr) && near_zero(dli)) break;
    }
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
inline constexpr int kTemmeX = 40;  // Temme series below (cancellation < 2^116), Steed CF2 above

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
    if (x < float256(kTemmeX)) temme_y(mu, x, y0, y1);
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
    if (x < float256(kTemmeX)) temme_k(mu, x, k0, k1);
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
// J or I for x so small that the series is its first term: (x/2)^ν/Γ(ν+1)·(1 ∓ tiny).
[[nodiscard]] inline bool bessel_tiny(detail_dd::FnEnv& env, BesselKind kind, float256 nu, float256 x, float256& r) {
    if ((kind != BesselKind::J && kind != BesselKind::I) || ilogb(x) >= -300) return false;
    SD v = bessel::series(dd(nu), x, kind == BesselKind::J ? -1 : 1);
    r = out(env, v, kind == BesselKind::J ? -1 : +1);
    return true;
}
} // namespace detail_sf

/// Bessel function of the first kind J_ν(x), ν >= 0, x >= 0.
[[nodiscard]] inline float256 cyl_bessel_j(float256 nu, float256 x) {
    detail_dd::FnEnv env;
    float256 r;
    if (detail_sf::bessel_special(env, detail_sf::BesselKind::J, nu, x, r)) return r;
    if (detail_sf::bessel_tiny(env, detail_sf::BesselKind::J, nu, x, r)) return r;
    return detail_sf::out(env, detail_sf::bessel::j(nu, x));
}
/// Bessel function of the second kind (Neumann) Y_ν(x), ν >= 0, x >= 0.
[[nodiscard]] inline float256 cyl_neumann(float256 nu, float256 x) {
    detail_dd::FnEnv env;
    float256 r;
    if (detail_sf::bessel_special(env, detail_sf::BesselKind::Y, nu, x, r)) return r;
    return detail_sf::out(env, detail_sf::bessel::y(nu, x));
}
/// Modified Bessel function of the first kind I_ν(x), ν >= 0, x >= 0.
[[nodiscard]] inline float256 cyl_bessel_i(float256 nu, float256 x) {
    detail_dd::FnEnv env;
    float256 r;
    if (detail_sf::bessel_special(env, detail_sf::BesselKind::I, nu, x, r)) return r;
    if (detail_sf::bessel_tiny(env, detail_sf::BesselKind::I, nu, x, r)) return r;
    return detail_sf::out(env, detail_sf::bessel::i(nu, x));
}
/// Modified Bessel function of the second kind K_ν(x), ν >= 0, x >= 0.
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
    SD v = second ? bessel::y(nu, x) : (ilogb(x) < -300 ? bessel::series(dd(nu), x, -1) : bessel::j(nu, x));
    v = mul(v, detail_dd::sqrt(detail_dd::pi() / (dd(x) * float256::two())));
    return out(env, v, (!second && ilogb(x) < -300) ? -1 : 2);
}
} // namespace detail_sf

/// Spherical Bessel function j_n(x), x >= 0.
[[nodiscard]] inline float256 sph_bessel(unsigned n, float256 x) { return detail_sf::sph(n, x, false); }
/// Spherical Neumann function y_n(x), x >= 0.
[[nodiscard]] inline float256 sph_neumann(unsigned n, float256 x) { return detail_sf::sph(n, x, true); }

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

// Continued fraction for I_x(a, b) (modified Lentz), x and 1 − x given as dd.
[[nodiscard]] inline dd betacf(const dd& a, const dd& b, const dd& x) {
    const dd one{float256::one()};
    const float256 tiny = ldexp(float256::one(), -20000);
    const dd qab = a + b, qap = a + one, qam = a - one;
    dd c = one;
    dd d = one - qab * x / qap;
    if (abs(d.hi) < tiny) d = dd(tiny);
    d = one / d;
    dd h = d;
    for (int m = 1; m < 100000000; ++m) {
        const dd dm = D(m), m2 = D(2 * m);
        dd aa = dm * (b - dm) * x / ((qam + m2) * (a + m2));
        d = one + aa * d;
        if (abs(d.hi) < tiny) d = dd(tiny);
        c = one + aa / c;
        if (abs(c.hi) < tiny) c = dd(tiny);
        d = one / d;
        h = h * d * c;
        aa = -((a + dm) * (qab + dm) * x) / ((a + m2) * (qap + m2));
        d = one + aa * d;
        if (abs(d.hi) < tiny) d = dd(tiny);
        c = one + aa / c;
        if (abs(c.hi) < tiny) c = dd(tiny);
        d = one / d;
        const dd del = d * c;
        h = h * del;
        if (near_one(del)) break;
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
    if (part.k < -600) return env.round(one, float256::zero(), 0, -1);  // 1 − tiny
    return (dd(one) - detail_dd::ldexp(part.m, part.k)).out(env);
}
} // namespace detail_sf
[[nodiscard]] inline float256 ibeta(float256 a, float256 b, float256 x) { return detail_sf::ibeta_impl(a, b, x, false); }
[[nodiscard]] inline float256 ibetac(float256 a, float256 b, float256 x) { return detail_sf::ibeta_impl(a, b, x, true); }

} // namespace base256
