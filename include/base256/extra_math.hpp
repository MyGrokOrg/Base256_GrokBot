#pragma once

/// Extra elementary functions layered on math.hpp.
/// Gamma, joint sine/cosine, and IEEE remquo live here.
#include "math.hpp"

namespace base256 {

[[nodiscard]] inline float256 nearbyint(float256 x) {
    // Default rounding mode is roundTiesToEven, same as rint.
    return rint(x);
}

[[nodiscard]] inline float256 fdim(float256 x, float256 y) {
    if (x.isnan() || y.isnan()) return float256::qnan();
    if (x > y) return x - y;
    return float256::zero(false);
}

// logb lives in float256.hpp (constexpr). A second definition here made
// <base256/base256.hpp> ill-formed.

namespace detail_spec {

// erf series on |x| <= 2. Odd, so a signed argument is fine.
[[nodiscard]] inline float256 erf_series(float256 x) {
    const float256 x2 = x * x;
    float256 term = x;
    float256 sum  = x;
    for (int n = 1; n <= 250; ++n) {
        term = -term * x2 * float256::from_int(static_cast<std::int64_t>(2 * n - 1))
               / (float256::from_int(static_cast<std::int64_t>(n)) *
                  float256::from_int(static_cast<std::int64_t>(2 * n + 1)));
        const float256 next = sum + term;
        if (next == sum) break;
        sum = next;
    }
    return sum * (float256::two() / sqrt(float256::pi()));
}

// erfc(x) for x > 2 via the continued fraction
//   exp(x^2) sqrt(pi) erfc(x) = 1 / (x + (1/2)/(x + (2/2)/(x + (3/2)/(x + ...))))
[[nodiscard]] inline float256 erfc_cf(float256 x) {
    const float256 x2 = x * x;
    if (x2.isinf() || x > float256::from_int(600)) return float256::zero();
    const float256 tiny = ldexp(float256::one(), -100);
    const float256 tol  = ldexp(float256::one(), -200);
    float256 f = x;
    float256 C = f;
    float256 D = float256::zero();
    for (int j = 1; j <= 500; ++j) {
        const float256 a = float256::from_int(static_cast<std::int64_t>(j)) * float256::half();
        D = x + a * D;
        if (D.iszero()) D = tiny;
        C = x + a / C;
        if (C.iszero()) C = tiny;
        D = float256::one() / D;
        const float256 delta = C * D;
        f = f * delta;
        if (abs(delta - float256::one()) <= tol) break;
    }
    const float256 num = exp(-x2);
    if (num.iszero()) return float256::zero();
    return num / (sqrt(float256::pi()) * f);
}

} // namespace detail_spec

[[nodiscard]] inline float256 erf(float256 x) {
    if (x.isnan()) return x;
    if (x.isinf()) return copysign(float256::one(), x);
    if (x.iszero()) return x;
    if (abs(x) <= float256::two()) return detail_spec::erf_series(x);
    const float256 y = float256::one() - detail_spec::erfc_cf(abs(x));
    return copysign(y, x);
}

[[nodiscard]] inline float256 erfc(float256 x) {
    if (x.isnan()) return x;
    if (x.isinf()) return x.signbit() ? float256::two() : float256::zero();
    if (abs(x) <= float256::two()) return float256::one() - detail_spec::erf_series(x);
    if (x.signbit()) return float256::two() - detail_spec::erfc_cf(abs(x));
    return detail_spec::erfc_cf(x);
}

[[nodiscard]] inline float256 asinh(float256 x) {
    if (x.isnan() || x.isinf() || x.iszero()) return x;
    const bool s = x.signbit();
    const float256 a = abs(x);
    float256 y;
    // asinh(x) = log(x + sqrt(x^2+1)); for large x, log(2x).
    const float256 thresh = ldexp(float256::one(), 64);
    if (a > thresh) {
        y = log(a) + float256::ln2();
    } else {
        y = log1p(a + a * a / (float256::one() + sqrt(a * a + float256::one())));
    }
    return s ? -y : y;
}

[[nodiscard]] inline float256 acosh(float256 x) {
    if (x.isnan()) return x;
    if (x < float256::one()) return float256::qnan();
    if (x == float256::one()) return float256::zero();
    if (x.isinf()) return x;
    const float256 thresh = ldexp(float256::one(), 64);
    if (x > thresh) return log(x) + float256::ln2();
    // acosh(x) = log(x + sqrt(x^2-1)) = log1p(x-1 + sqrt((x-1)(x+1)))
    const float256 xm1 = x - float256::one();
    return log1p(xm1 + sqrt(xm1 * (x + float256::one())));
}

[[nodiscard]] inline float256 atanh(float256 x) {
    if (x.isnan()) return x;
    if (x.iszero()) return x;
    const float256 a = abs(x);
    if (a > float256::one()) return float256::qnan();
    if (a == float256::one()) return float256::inf(x.signbit());
    // atanh(x) = 1/2 log((1+x)/(1-x)) = 1/2 log1p(2x/(1-x))
    const float256 y = log1p((a + a) / (float256::one() - a)) * float256::half();
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
[[nodiscard]] inline float256 remquo(float256 x, float256 y, int* quo) {
    if (quo) *quo = 0;
    if (x.isnan() || y.isnan() || y.iszero() || x.isinf()) return float256::qnan();
    if (y.isinf()) return x;
    const float256 q = rint(x / y);
    const float256 r = x - q * y;
    if (quo && q.isfinite() && !q.iszero()) {
        const float256 aq = abs(q);
        const float256 high = trunc(ldexp(aq, -31));
        const float256 low = aq - ldexp(high, 31);
        int bits = static_cast<int>(low.to_double());
        if (q.signbit()) bits = -bits;
        *quo = bits;
    }
    return r;
}

} // namespace base256
