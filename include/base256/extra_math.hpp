#pragma once

/// Extra elementary functions layered on math.hpp.
/// Gamma, joint sine/cosine, and IEEE remquo live here.
#include "math.hpp"

namespace base256 {

// nearbyint lives in float256.hpp (current rounding mode, no inexact).

[[nodiscard]] inline float256 fdim(float256 x, float256 y) {
    if (x.isnan() || y.isnan()) return float256::qnan();
    if (x > y) return x - y;
    return float256::zero(false);
}

// logb lives in float256.hpp (constexpr). A second definition here made
// <base256/base256.hpp> ill-formed.


namespace detail_spec {
// erf(a) for 0 < a < 4 in dd: 2/√π e^{-a²} Σ 2^n a^{2n+1} / (1·3·…·(2n+1)).
// All terms are positive, so the sum has no cancellation.
[[nodiscard]] inline detail_dd::dd erf_series_dd(float256 a) {
    using detail_dd::dd;
    const dd x2 = detail_dd::two_prod(a, a);
    const dd two_x2 = detail_dd::ldexp(x2, 1);
    dd term{a}, sum{a};
    const float256 rel = ldexp(float256::one(), -490);
    for (int n = 1; n <= 2000; ++n) {
        term = term * two_x2 / dd(float256::from_int(2 * n + 1));
        sum = sum + term;
        if (term.hi < abs(sum.hi) * rel) break;
    }
    const auto& c = detail_dd::consts();
    const dd two_over_sqrtpi = dd(float256::two()) / dd(c.sqrtpi[0], c.sqrtpi[1]);
    return sum * detail_dd::exp(-x2) * two_over_sqrtpi;
}

// erfc(a) for a >= 4 as m·2^k: e^{-a²} / (√π · K(a)), with the Laplace
// continued fraction K = a + (1/2)/(a + 1/(a + (3/2)/(a + …))) evaluated
// backward from a fixed depth (chosen so the truncation error is < 2^-260).
[[nodiscard]] inline detail_dd::Scaled erfc_cf_dd(float256 a) {
    using detail_dd::dd;
    const int depth = a < float256(6) ? 400 : a < float256(8) ? 250 : a < float256(10) ? 200
                    : a < float256(20) ? 120 : 80;
    dd f{a};
    for (int n = depth; n >= 1; --n) f = dd(a) + dd(float256::from_int(n) * float256::half()) / f;
    const auto& c = detail_dd::consts();
    detail_dd::Scaled e = detail_dd::exp_scaled(-detail_dd::two_prod(a, a));
    e.m = e.m / (f * dd(c.sqrtpi[0], c.sqrtpi[1]));
    return e;
}
} // namespace detail_spec

[[nodiscard]] inline float256 erf(float256 x) {
    const fenv::scoped_default_env env_guard_; // RNE inside, caller's flags kept
    using detail_dd::dd;
    if (x.isnan()) return x;
    if (x.isinf()) return copysign(float256::one(), x);
    if (x.iszero()) return x;
    const float256 a = abs(x);
    float256 y;
    if (a < ldexp(float256::one(), -120)) {
        // 2a/√π (1 - a²/3): the a³ term is below 2^-240 relative.
        const auto& c = detail_dd::consts();
        y = (dd(a + a) / dd(c.sqrtpi[0], c.sqrtpi[1])).round();
    } else if (a < float256(4)) {
        y = detail_spec::erf_series_dd(a).round();
    } else if (a > float256(16)) {
        y = float256::one();  // erfc(16) < 2^-370
    } else {
        const detail_dd::Scaled e = detail_spec::erfc_cf_dd(a);
        y = (dd(float256::one()) - detail_dd::ldexp(e.m, e.k)).round();
    }
    return x.signbit() ? -y : y;
}

[[nodiscard]] inline float256 erfc(float256 x) {
    const fenv::scoped_default_env env_guard_; // RNE inside, caller's flags kept
    using detail_dd::dd;
    if (x.isnan()) return x;
    if (x.isinf()) return x.signbit() ? float256::two() : float256::zero();
    if (x.iszero()) return float256::one();
    const float256 a = abs(x);
    if (x.signbit()) {
        // 1 + erf(a)
        if (a > float256(16)) return float256::two();
        if (a < ldexp(float256::one(), -120)) return float256::one();
        if (a < float256(4)) return (dd(float256::one()) + detail_spec::erf_series_dd(a)).round();
        const detail_dd::Scaled e = detail_spec::erfc_cf_dd(a);
        return (dd(float256::two()) - detail_dd::ldexp(e.m, e.k)).round();
    }
    if (a < ldexp(float256::one(), -120)) return float256::one();
    if (a < float256(4)) return (dd(float256::one()) - detail_spec::erf_series_dd(a)).round();
    if (a > float256(520)) return float256::zero();  // erfc(520) < 2^-390000
    const detail_dd::Scaled e = detail_spec::erfc_cf_dd(a);
    return detail_dd::round_scaled(e.m, e.k);
}

[[nodiscard]] inline float256 asinh(float256 x) {
    const fenv::scoped_default_env env_guard_; // RNE inside, caller's flags kept
    using detail_dd::dd;
    if (x.isnan() || x.isinf() || x.iszero()) return x;
    const float256 a = abs(x);
    float256 y;
    if (a < ldexp(float256::one(), -120)) {
        y = (dd(a) - dd(a * a * a / float256::from_int(6))).round();
    } else if (a > ldexp(float256::one(), 300)) {
        y = (detail_dd::log(dd(a)) + detail_dd::ln2()).round();  // log(2a); 1/(4a^2) is negligible
    } else {
        const dd r = detail_dd::sqrt(detail_dd::two_prod(a, a) + dd(float256::one()));
        y = detail_dd::log(dd(a) + r).round();
    }
    return x.signbit() ? -y : y;
}

[[nodiscard]] inline float256 acosh(float256 x) {
    const fenv::scoped_default_env env_guard_; // RNE inside, caller's flags kept
    using detail_dd::dd;
    if (x.isnan()) return x;
    if (x < float256::one()) return float256::qnan();
    if (x == float256::one()) return float256::zero();
    if (x.isinf()) return x;
    if (x > ldexp(float256::one(), 300)) return (detail_dd::log(dd(x)) + detail_dd::ln2()).round();
    const dd xm1 = detail_dd::two_sum(x, -float256::one());
    const dd xp1 = detail_dd::two_sum(x, float256::one());
    return detail_dd::log(dd(x) + detail_dd::sqrt(xm1 * xp1)).round();
}

[[nodiscard]] inline float256 atanh(float256 x) {
    const fenv::scoped_default_env env_guard_; // RNE inside, caller's flags kept
    using detail_dd::dd;
    if (x.isnan()) return x;
    if (x.iszero()) return x;
    const float256 a = abs(x);
    if (a > float256::one()) return float256::qnan();
    if (a == float256::one()) return float256::inf(x.signbit());
    float256 y;
    if (a < ldexp(float256::one(), -120)) {
        y = (dd(a) + dd(a * a * a / float256::from_int(3))).round();
    } else {
        // (log(1 + a) - log(1 - a)) / 2, with 1 ± a exact as dd.
        const dd lp = detail_dd::log(detail_dd::two_sum(float256::one(), a));
        const dd lm = detail_dd::log(detail_dd::two_sum(float256::one(), -a));
        y = detail_dd::ldexp(lp - lm, -1).round();
    }
    return x.signbit() ? -y : y;
}


namespace detail_gamma {

// Stirling series for ln Γ(z), z >= 48:
//   (z-1/2) ln z - z + ln(2π)/2 + Σ_{k=1}^{40} B_{2k} / (2k(2k-1) z^{2k-1})
// Forty terms leave a truncation under 2^{-250} for every z >= 48.
inline constexpr int kStirlingShift = 48;
inline constexpr int kStirlingTerms = 40;

[[nodiscard]] inline const float256* stirling_coeff() {
    static const float256 c[kStirlingTerms] = {
        float256("8.333333333333333333333333333333333333333333333333333333333333333333333333333333e-2"),
        float256("-2.777777777777777777777777777777777777777777777777777777777777777777777777777778e-3"),
        float256("7.936507936507936507936507936507936507936507936507936507936507936507936507936508e-4"),
        float256("-5.952380952380952380952380952380952380952380952380952380952380952380952380952381e-4"),
        float256("8.417508417508417508417508417508417508417508417508417508417508417508417508417508e-4"),
        float256("-1.917526917526917526917526917526917526917526917526917526917526917526917526917527e-3"),
        float256("6.410256410256410256410256410256410256410256410256410256410256410256410256410256e-3"),
        float256("-2.955065359477124183006535947712418300653594771241830065359477124183006535947712e-2"),
        float256("1.796443723688305731649384900158893966943502547217717496355267253100070437531738e-1"),
        float256("-1.392432216905901116427432216905901116427432216905901116427432216905901116427432e+0"),
        float256("1.340286404416839199447895100069013112491373360938578329882677708764665286404417e+1"),
        float256("-1.568482846260020173063651324520889738281042628868715825237564367999150607846260e+2"),
        float256("2.193103333333333333333333333333333333333333333333333333333333333333333333333333e+3"),
        float256("-3.610877125372498935717326521924223073648361004682843763303533418475947211579395e+4"),
        float256("6.914722688513130671083952507756734675533340716877980504231894665710016099337568e+5"),
        float256("-1.523822153940741619228336495888678051865907653383934218848829854522454142947502e+7"),
        float256("3.829007513914141414141414141414141414141414141414141414141414141414141414141414e+8"),
        float256("-1.088226603578439108901514916552510537472943487981081966044372059409653394615801e+10"),
        float256("3.473202837650022522522522522522522522522522522522522522522522522522522522522523e+11"),
        float256("-1.236960214226927445425171034927132488108097864195425171034927132488108097864195e+13"),
        float256("4.887880647930793350758151625180229021084705389056738218070362953273576399741216e+14"),
        float256("-2.132033396091937389697505898213683855746545331985170205594876980114593865857769e+16"),
        float256("1.021775296525700077565287628053585500394011032308904649330181245074862096138692e+18"),
        float256("-5.357547217330020361082770919196920448484904054365881649986781401049235842727708e+19"),
        float256("3.061578263704883415043151051329622758194186765615337043908472479901051329622758e+21"),
        float256("-1.899991742639920405029371429306942902947342458996177087187076088296954001726477e+23"),
        float256("1.276337403382883414923495137769782597654163360882990144823974681637707125956050e+25"),
        float256("-9.252847176120416307230242348347622779519331243469174503657262277951933124346917e+26"),
        float256("7.218822595185610297836050187301637922489840420259688769947467538903756566557637e+28"),
        float256("-6.045183405995856967743148238754547286066144395967196207406301608096013351953622e+30"),
        float256("5.420670471570094545193477814826100013661202185792349726775956284153005464480874e+32"),
        float256("-5.192957815314081946700194764391857684699706271397447868036103330220977279800809e+34"),
        float256("5.303658855119700596654839243069758643699292635405549097956625536119371291337547e+36"),
        float256("-5.763325348164964013894435850780992555190737562189054726368159203980099502487562e+38"),
        float256("6.651155714848453937516520145810555951039739359454928958909363073478199014022651e+40"),
        float256("-8.137378358136680538716172632093575691840689164973879262367945004499733015020849e+42"),
        float256("1.053696695335714180375480492764181018964837337501141552511415525114155251141553e+45"),
        float256("-1.441818059996220626180537780151181280957033206366421111111111111111111111111111e+47"),
        float256("2.081735652208956546242480824126356231131734326414997918933588011345696698318760e+49"),
        float256("-3.167022663488666182741349556774256134291806983042075303039154461638706957705196e+51"),
    };
    return c;
}

[[nodiscard]] inline float256 stirling_lgamma(float256 z) {
    const float256 ln = log(z);
    float256 sum = (z - float256::half()) * ln - z;
    sum = sum + log(float256::pi() * float256::two()) * float256::half();
    const float256* c = stirling_coeff();
    const float256 z2 = z * z;
    float256 p = float256::one() / z;
    for (int k = 0; k < kStirlingTerms; ++k) {
        const float256 term = c[k] * p;
        const float256 next = sum + term;
        if (next == sum) break;
        sum = next;
        p = p / z2;
    }
    return sum;
}

// Rising product 1*2*...*n for 0 <= n <= 53, where n! is an integer < 2^237
// and the product is therefore exact.
[[nodiscard]] inline float256 factorial(int n) {
    float256 p = float256::one();
    for (int k = 2; k <= n; ++k) p = p * float256(k);
    return p;
}

[[nodiscard]] inline float256 lgamma_pos(float256 x) {
    // Γ(n) = (n-1)! exactly for n = 1..54.
    if (trunc(x) == x && x <= float256(54)) {
        const int n = static_cast<int>(x.to_double());
        if (n <= 2) return float256::zero();
        return log(factorial(n - 1));
    }
    if (x >= float256(kStirlingShift)) return stirling_lgamma(x);
    // ln Γ(x) = ln Γ(x+s) - ln(x(x+1)...(x+s-1)), with x+s >= 48.
    // Step an integer counter: x+1 would stick at 1 for a subnormal x.
    const int steps = kStirlingShift - static_cast<int>(floor(x).to_double());
    float256 acc = float256::zero();
    for (int i = 0; i < steps; ++i) acc = acc + log(x + float256(i));
    return stirling_lgamma(x + float256(steps)) - acc;
}

// sin(π f) for f in (0, 1), reduced into [0, π/4] so the Taylor argument
// stays where it is accurate. Positive.
[[nodiscard]] inline float256 sin_pi_unit(float256 f) {
    float256 w = f;
    if (w > float256::half()) w = float256::one() - w;
    if (w > float256(0.25))
        return cos(float256::pi() * (float256::half() - w));
    return sin(float256::pi() * w);
}

} // namespace detail_gamma

// Natural log of |Γ(x)|. Poles (±0 and non-positive integers) return +∞.
[[nodiscard]] inline float256 lgamma(float256 x) {
    const fenv::scoped_default_env env_guard_; // RNE inside, caller's flags kept
    if (x.isnan()) return x;
    if (x.isinf()) return float256::inf(false);
    if (x.iszero()) return float256::inf(false);
    if (x > float256::zero()) return detail_gamma::lgamma_pos(x);
    const float256 fl = floor(x);
    if (fl == x) return float256::inf(false); // non-positive integer
    const float256 f = x - fl;                // (0, 1)
    const float256 s = detail_gamma::sin_pi_unit(f);
    if (s.iszero()) return float256::inf(false);
    // ln|Γ(x)| = ln π − ln|sin(πx)| − ln Γ(1−x), and 1−x > 1.
    const float256 lg = detail_gamma::lgamma_pos(float256::one() - x);
    if (lg.isinf()) return float256::inf(true);
    return log(float256::pi()) - log(s) - lg;
}

// Γ(x). Negative integers and −∞ are NaN; tgamma(±0) is ±∞.
[[nodiscard]] inline float256 tgamma(float256 x) {
    const fenv::scoped_default_env env_guard_; // RNE inside, caller's flags kept
    if (x.isnan()) return x;
    if (x.iszero()) return float256::inf(x.signbit());
    if (x.isinf()) return x.signbit() ? float256::qnan() : x;
    if (!x.signbit() && trunc(x) == x && x <= float256(54)) {
        const int n = static_cast<int>(x.to_double());
        return detail_gamma::factorial(n - 1);
    }
    if (x.signbit() && trunc(x) == x) return float256::qnan();
    const float256 mag = exp(lgamma(x));
    if (!x.signbit()) return mag;
    const float256 half_i = floor(x) * float256::half();
    const bool odd = trunc(half_i) != half_i;
    return odd ? -mag : mag;
}

// sin and cos of the same argument, one reduction.
inline void sincos(float256 x, float256* sin_out, float256* cos_out) {
    const fenv::scoped_default_env env_guard_; // RNE inside, caller's flags kept
    if (sin_out == nullptr && cos_out == nullptr) return;
    if (x.isnan() || x.isinf()) {
        const float256 nan = x.isnan() ? x : float256::qnan();
        if (sin_out) *sin_out = nan;
        if (cos_out) *cos_out = nan;
        return;
    }
    if (x.iszero()) {
        if (sin_out) *sin_out = x;
        if (cos_out) *cos_out = float256::one();
        return;
    }
    int quad = 0;
    const bool neg = x.signbit();
    {
        detail_dd::dd sd, cd;
        if (detail_math::trig_dd(abs(x), sd, cd, quad)) {
            detail_dd::dd sv = (quad & 1) ? cd : sd;
            detail_dd::dd cv = (quad & 1) ? -sd : cd;
            if (quad >= 2) { sv = -sv; cv = -cv; }
            const float256 so = sv.round();
            if (sin_out) *sin_out = neg ? -so : so;
            if (cos_out) *cos_out = cv.round();
            return;
        }
    }
    const float256 r = reduce_trig(x, quad);
    float256 ss;
    float256 cc;
    if (r.iszero()) {
        switch (quad) {
            case 0:
                ss = copysign(float256::zero(), x);
                cc = float256::one();
                break;
            case 1:
                ss = neg ? -float256::one() : float256::one();
                cc = float256::zero();
                break;
            case 2:
                ss = copysign(float256::zero(), x);
                cc = -float256::one();
                break;
            default:
                ss = neg ? float256::one() : -float256::one();
                cc = float256::zero();
                break;
        }
    } else {
        const float256 s = detail_math::taylor_sin(r);
        const float256 c = detail_math::taylor_cos(r);
        switch (quad) {
            case 0: ss = s; cc = c; break;
            case 1: ss = c; cc = -s; break;
            case 2: ss = -s; cc = -c; break;
            default: ss = -c; cc = s; break;
        }
        if (neg) ss = -ss;
    }
    if (sin_out) *sin_out = ss;
    if (cos_out) *cos_out = cc;
}

// IEEE remainder, plus the low 31 bits of the integer quotient (sign of x/y).
/// IEEE remainder plus the sign and low 31 bits of the rounded quotient.
[[nodiscard]] constexpr float256 remquo(float256 x, float256 y, int* quo) noexcept {
    return float256::remainder_impl(x, y, true, quo);
}

} // namespace base256
