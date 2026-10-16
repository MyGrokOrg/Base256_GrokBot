#pragma once

/// Complex numbers with binary256 components: base256::complex.
///
///     #include <base256/complex.hpp>
///     using namespace base256::literals;
///     base256::complex z(1.5_f256, -2_f256);
///     auto w = base256::exp(z) * 3_i256;
///
/// The interface follows std::complex (constructors, real/imag, arithmetic
/// with scalars and complex values, ==, abs, arg, norm, conj, proj, polar,
/// exp, log, log10, pow, sqrt, the trigonometric, hyperbolic and inverse
/// functions, stream I/O) plus std::format. Each component of each result is
/// rounded once in the current float256 rounding mode with IEEE flags; the
/// special values (infinities, NaNs, signed zeros, branch cuts) follow C
/// Annex G. See the README section "Complex numbers" for the exact guarantees.
///
/// A separate header (like <complex>): <base256/base256.hpp> does not include it.

#include <format>
#include <complex>
#include <istream>
#include <ostream>
#include <sstream>
#include <string>
#include <string_view>
#include <type_traits>

#include "float256.hpp"
#include "math.hpp"
#include "extra_math.hpp"
#include "detail/cx_engine.hpp"

namespace base256 {

class complex {
public:
    using value_type = float256;

    constexpr complex(float256 re = float256(), float256 im = float256()) noexcept : re_(re), im_(im) {}
    constexpr complex(const complex&) noexcept = default;
    constexpr complex& operator=(const complex&) noexcept = default;
    /// From std::complex<float/double/long double> (exact).
    template <class T>
        requires std::is_floating_point_v<T>
    explicit complex(const std::complex<T>& z) noexcept : re_(float256(z.real())), im_(float256(z.imag())) {}
    /// To std::complex<float/double/long double> (each component rounded in the current mode).
    template <class T>
        requires std::is_floating_point_v<T>
    explicit operator std::complex<T>() const noexcept {
        if constexpr (std::is_same_v<T, float>) return {re_.to_float(), im_.to_float()};
        else if constexpr (std::is_same_v<T, double>) return {re_.to_double(), im_.to_double()};
        else return {static_cast<T>(re_.to_long_double()), static_cast<T>(im_.to_long_double())};
    }

    [[nodiscard]] constexpr float256 real() const noexcept { return re_; }
    [[nodiscard]] constexpr float256 imag() const noexcept { return im_; }
    constexpr void real(float256 v) noexcept { re_ = v; }
    constexpr void imag(float256 v) noexcept { im_ = v; }

    constexpr complex& operator=(float256 v) noexcept { re_ = v; im_ = float256(); return *this; }
    constexpr complex& operator+=(float256 v) noexcept { re_ += v; return *this; }
    constexpr complex& operator-=(float256 v) noexcept { re_ -= v; return *this; }
    constexpr complex& operator*=(float256 v) noexcept { re_ *= v; im_ *= v; return *this; }
    constexpr complex& operator/=(float256 v) noexcept { re_ /= v; im_ /= v; return *this; }
    constexpr complex& operator+=(const complex& w) noexcept { re_ += w.re_; im_ += w.im_; return *this; }
    constexpr complex& operator-=(const complex& w) noexcept { re_ -= w.re_; im_ -= w.im_; return *this; }
    complex& operator*=(const complex& w);
    complex& operator/=(const complex& w);

private:
    float256 re_, im_;
};

// ------------------------------------------------------------ basic operations

[[nodiscard]] constexpr float256 real(const complex& z) noexcept { return z.real(); }
[[nodiscard]] constexpr float256 imag(const complex& z) noexcept { return z.imag(); }
[[nodiscard]] constexpr complex conj(const complex& z) noexcept { return {z.real(), -z.imag()}; }
/// Riemann-sphere projection (Annex G cproj): any infinite component gives +∞ + i·copysign(0, imag).
[[nodiscard]] constexpr complex proj(const complex& z) noexcept {
    if (z.real().isinf() || z.imag().isinf()) return {float256::inf(false), float256::zero(z.imag().signbit())};
    return z;
}

[[nodiscard]] constexpr complex operator+(const complex& z) noexcept { return z; }
[[nodiscard]] constexpr complex operator-(const complex& z) noexcept { return {-z.real(), -z.imag()}; }
[[nodiscard]] constexpr complex operator+(const complex& z, const complex& w) noexcept {
    return {z.real() + w.real(), z.imag() + w.imag()};
}
[[nodiscard]] constexpr complex operator+(const complex& z, float256 v) noexcept { return {z.real() + v, z.imag()}; }
[[nodiscard]] constexpr complex operator+(float256 v, const complex& z) noexcept { return {v + z.real(), z.imag()}; }
[[nodiscard]] constexpr complex operator-(const complex& z, const complex& w) noexcept {
    return {z.real() - w.real(), z.imag() - w.imag()};
}
[[nodiscard]] constexpr complex operator-(const complex& z, float256 v) noexcept { return {z.real() - v, z.imag()}; }
/// Annex G: x − (u + iv) = (x − u) − iv (the imaginary part is negated, not 0 − v).
[[nodiscard]] constexpr complex operator-(float256 v, const complex& z) noexcept { return {v - z.real(), -z.imag()}; }
[[nodiscard]] constexpr complex operator*(const complex& z, float256 v) noexcept { return {z.real() * v, z.imag() * v}; }
[[nodiscard]] constexpr complex operator*(float256 v, const complex& z) noexcept { return {v * z.real(), v * z.imag()}; }
[[nodiscard]] constexpr complex operator/(const complex& z, float256 v) noexcept { return {z.real() / v, z.imag() / v}; }
[[nodiscard]] complex operator*(const complex& z, const complex& w);
[[nodiscard]] complex operator/(const complex& z, const complex& w);
[[nodiscard]] inline complex operator/(float256 v, const complex& w) { return complex(v, float256()) / w; }

[[nodiscard]] constexpr bool operator==(const complex& z, const complex& w) noexcept {
    return z.real() == w.real() && z.imag() == w.imag();
}
[[nodiscard]] constexpr bool operator==(const complex& z, float256 v) noexcept {
    return z.real() == v && z.imag() == float256();
}

inline complex& complex::operator*=(const complex& w) { return *this = *this * w; }
inline complex& complex::operator/=(const complex& w) { return *this = *this / w; }

namespace detail_cx {

using F = float256;

inline void quiet(F& v) {
    if (v.issignaling()) { detail::raise_flag(fe_invalid); v = v.quieted(); }
}
[[nodiscard]] inline bool sb(F v) { return v.signbit(); }
[[nodiscard]] inline F pinf() { return F::inf(false); }
[[nodiscard]] inline F qnan() { return F::qnan(); }

// real functions as plain function pointers (evaluated after the environment is restored)
inline F r_exp(F t) { return base256::exp(t); }
inline F r_log(F t) { return base256::log(t); }
inline F r_log10(F t) { return base256::log10(t); }
inline F r_sqrt(F t) { return base256::sqrt(t); }
inline F r_sin(F t) { return base256::sin(t); }
inline F r_cos(F t) { return base256::cos(t); }
inline F r_tan(F t) { return base256::tan(t); }
inline F r_asin(F t) { return base256::asin(t); }
inline F r_acos(F t) { return base256::acos(t); }
inline F r_atan(F t) { return base256::atan(t); }
inline F r_sinh(F t) { return base256::sinh(t); }
inline F r_cosh(F t) { return base256::cosh(t); }
inline F r_tanh(F t) { return base256::tanh(t); }
inline F r_asinh(F t) { return base256::asinh(t); }
inline F r_acosh(F t) { return base256::acosh(t); }
inline F r_atanh(F t) { return base256::atanh(t); }
inline F r_atan2(F a, F b) { return base256::atan2(a, b); }
inline F r_pow(F a, F b) { return base256::pow(a, b); }

// ±π, ±π/2, ±π/4, ±3π/4 rounded in the caller's mode (atan2 is correctly rounded)
[[nodiscard]] inline Comp c_pi(bool n) { return Comp::fn2(r_atan2, F::zero(n), -F::one()); }
[[nodiscard]] inline Comp c_pio2(bool n) { return Comp::fn2(r_atan2, n ? -F::one() : F::one(), F::zero()); }
[[nodiscard]] inline Comp c_pio4(bool n) { return Comp::fn2(r_atan2, n ? -F::one() : F::one(), F::one()); }
[[nodiscard]] inline Comp c_3pio4(bool n) { return Comp::fn2(r_atan2, n ? -F::one() : F::one(), -F::one()); }
[[nodiscard]] inline Comp c_ex(F v) { return Comp::ex(v); }
[[nodiscard]] inline Comp c_zero(bool n) { return Comp::ex(F::zero(n)); }
[[nodiscard]] inline Comp c_inf(bool n) { return Comp::ex(F::inf(n)); }
[[nodiscard]] inline Comp c_nan() { return Comp::ex(F::qnan()); }

[[nodiscard]] inline bool cos_neg(F y) { fenv::scoped_default_env g; return base256::cos(y).signbit(); }
[[nodiscard]] inline bool sin_neg(F y) { fenv::scoped_default_env g; return base256::sin(y).signbit(); }
inline void raise(unsigned f) { detail::raise_flag(f); }

enum class Out : unsigned char { id, neg_i, pos_i };
inline void apply(Out t, Comp& re, Comp& im) {
    if (t == Out::neg_i) { std::swap(re, im); im.negate(); }       // −i(u + iv) = v − iu
    else if (t == Out::pos_i) { std::swap(re, im); re.negate(); }  // i(u + iv) = −v + iu
}
[[nodiscard]] inline F value(const Comp& c) { return c.k == Comp::exact ? c.val : eval_realfn(c); }

/// special(x, y, re, im) fills both components and returns true when no
/// multiprecision evaluation is needed; otherwise kernel(x, y, P, re, im)
/// runs in the Ziv loop. `t` maps the result (derived functions).
template <class Sp, class Kr>
[[nodiscard]] complex run(F x, F y, Sp&& special, Kr&& kernel, Out t = Out::id) {
    Comp re, im;
    if (special(x, y, re, im)) {
        apply(t, re, im);
        const F a = value(re);
        const F b = value(im);
        return {a, b};
    }
    const auto r = ziv([&](int P, Comp& a, Comp& b) {
        kernel(x, y, P, a, b);
        apply(t, a, b);
    });
    return {r[0], r[1]};
}

[[nodiscard]] inline bool gt(F a, double b) { return a > F(b); }
[[nodiscard]] inline MP half(const MP& a) { return M::ldexp(a, -1); }
[[nodiscard]] inline MP sqr(const MP& a, int Q) { return M::mul(a, a, Q); }
[[nodiscard]] inline MP signed_(MP a, bool n) { if (n) a = M::neg(a); return a; }
[[nodiscard]] inline int sgnf(F v) { return v.signbit() ? -1 : 1; }


// ------------------------------------------------------------------ tiny arguments
// For |x|, |y| < 2^-64 the near-identity components (sinh z ≈ z, cosh z ≈ 1 + z²/2,
// e^z ≈ 1 + z, tanh, asin, acos, atanh ≈ z) are evaluated as an exact leading term
// (the Comp base) plus a power-series rest known to full relative precision, so the
// rounding direction is known even when the lead is itself a rounding boundary.
struct CM {
    MP r, i;
};
[[nodiscard]] inline CM cm_mul(const CM& a, const CM& b, int Q) {
    return {M::sub(M::mul(a.r, b.r, Q), M::mul(a.i, b.i, Q), Q), M::add(M::mul(a.r, b.i, Q), M::mul(a.i, b.r, Q), Q)};
}
[[nodiscard]] inline CM cm_add(const CM& a, const CM& b, int Q) { return {M::add(a.r, b.r, Q), M::add(a.i, b.i, Q)}; }
[[nodiscard]] inline CM cm_scale(const CM& a, std::uint32_t num, std::uint32_t den, int Q) {
    CM r = a;
    if (num != 1) { r.r = M::mul_u(r.r, num, Q); r.i = M::mul_u(r.i, num, Q); }
    if (den != 1) { r.r = M::div_u(r.r, den, Q); r.i = M::div_u(r.i, den, Q); }
    return r;
}
[[nodiscard]] inline long cm_top(const CM& a) {
    const long tr = a.r.zero() ? LONG_MIN / 4 : a.r.top(), ti = a.i.zero() ? LONG_MIN / 4 : a.i.top();
    return std::max(tr, ti);
}
/// Series stop rule, per component: the next term t is below the partial sum by
/// 2^-(Q+8) (a component whose sum is still 0 waits for a nonzero term, up to 2^-8Q
/// below the first term `scale`).
[[nodiscard]] inline bool cm_converged(const CM& sum, const CM& t, long scale, int Q) {
    if (cm_top(t) < scale - 8L * Q) return true;
    auto done = [&](const MP& sc, const MP& tc) { return tc.zero() || (!sc.zero() && tc.top() < sc.top() - Q - 8); };
    return done(sum.r, t.r) && done(sum.i, t.i);
}
[[nodiscard]] inline bool tiny_arg(F x, F y) {
    return (x.iszero() || ilogb(x) < -64) && (y.iszero() || ilogb(y) < -64);
}
/// Σ_{k>=k0} c_k·w^k with t_{k0} = first, t_{k+1} = t_k·step·num(k)/den(k), until the
/// terms fall 2^-(Q+8) below the first one (step = z or z², |step| < 2^-128).
template <class Ratio>
[[nodiscard]] inline CM cm_series(CM t, const CM& step, long k, Ratio&& ratio, int Q) {
    CM sum{MP{}, MP{}};
    const long scale = cm_top(t);
    while (!(t.r.zero() && t.i.zero())) {
        sum = cm_add(sum, t, Q);
        std::uint32_t num = 1, den = 1;
        ratio(k, num, den);
        t = cm_scale(cm_mul(t, step, Q), num, den, Q);
        ++k;
        if (cm_converged(sum, t, scale, Q)) break;
    }
    return sum;
}
/// Σ_{k>=1} z^(2k+1)/(2k+1)!  (sinh z − z)
[[nodiscard]] inline CM tiny_sinh_rest(const CM& z, int Q) {
    const CM z2 = cm_mul(z, z, Q);
    const CM t = cm_scale(cm_mul(z2, z, Q), 1, 6, Q);
    return cm_series(t, z2, 1, [](long k, std::uint32_t& n, std::uint32_t& d) {
        n = 1; d = static_cast<std::uint32_t>((2 * k + 2) * (2 * k + 3));
    }, Q);
}
/// Σ_{k>=2} z^(2k)/(2k)!  (cosh z − 1 − z²/2)
[[nodiscard]] inline CM tiny_cosh_rest(const CM& z, int Q) {
    const CM z2 = cm_mul(z, z, Q);
    const CM t = cm_scale(cm_mul(z2, z2, Q), 1, 24, Q);
    return cm_series(t, z2, 2, [](long k, std::uint32_t& n, std::uint32_t& d) {
        n = 1; d = static_cast<std::uint32_t>((2 * k + 1) * (2 * k + 2));
    }, Q);
}
/// Σ_{k>=2} z^k/k!  (e^z − 1 − z)
[[nodiscard]] inline CM tiny_exp_rest(const CM& z, int Q) {
    const CM t = cm_scale(cm_mul(z, z, Q), 1, 2, Q);
    return cm_series(t, z, 2, [](long k, std::uint32_t& n, std::uint32_t& d) {
        n = 1; d = static_cast<std::uint32_t>(k + 1);
    }, Q);
}
/// Σ_{k>=1} z^(2k+1)/(2k+1)  (atanh z − z)
[[nodiscard]] inline CM tiny_atanh_rest(const CM& z, int Q) {
    const CM z2 = cm_mul(z, z, Q);
    const CM p = cm_mul(z2, z, Q);  // z^(2k+1), divided by 2k+1 per term below
    CM sum{MP{}, MP{}}, t = p;
    const long scale = cm_top(p);
    for (long k = 1; !(t.r.zero() && t.i.zero()); ++k) {
        sum = cm_add(sum, cm_scale(t, 1, static_cast<std::uint32_t>(2 * k + 1), Q), Q);
        t = cm_mul(t, z2, Q);
        if (cm_converged(sum, t, scale, Q)) break;
    }
    return sum;
}
/// Σ_{k>=1} a_k z^(2k+1)/(2k+1), a_k = Π_{j<=k} (2j−1)/(2j)  (asin z − z)
[[nodiscard]] inline CM tiny_asin_rest(const CM& z, int Q) {
    const CM z2 = cm_mul(z, z, Q);
    CM sum{MP{}, MP{}}, t = cm_scale(cm_mul(z2, z, Q), 1, 2, Q);  // a_1 z³
    const long scale = cm_top(t);
    for (long k = 1; !(t.r.zero() && t.i.zero()); ++k) {
        sum = cm_add(sum, cm_scale(t, 1, static_cast<std::uint32_t>(2 * k + 1), Q), Q);
        t = cm_scale(cm_mul(t, z2, Q), static_cast<std::uint32_t>(2 * k + 1), static_cast<std::uint32_t>(2 * k + 2), Q);
        if (cm_converged(sum, t, scale, Q)) break;
    }
    return sum;
}
/// tanh z − z = (sinh z − z cosh z)/cosh z, numerator −Σ_{k>=1} 2k·z^(2k+1)/(2k+1)!
[[nodiscard]] inline CM tiny_tanh_rest(const CM& z, int Q) {
    const CM z2 = cm_mul(z, z, Q);
    CM num{MP{}, MP{}}, t = cm_scale(cm_mul(z2, z, Q), 1, 6, Q);  // z^(2k+1)/(2k+1)!
    const long scale = cm_top(t);
    for (long k = 1; !(t.r.zero() && t.i.zero()); ++k) {
        num = cm_add(num, cm_scale(t, static_cast<std::uint32_t>(2 * k), 1, Q), Q);
        t = cm_scale(cm_mul(t, z2, Q), 1, static_cast<std::uint32_t>((2 * k + 2) * (2 * k + 3)), Q);
        if (cm_converged(num, t, scale, Q)) break;
    }
    num = {M::neg(num.r), M::neg(num.i)};
    const CM c = tiny_cosh_rest(z, Q);
    const CM D{M::add(M::add(one(), half(z2.r), Q), c.r, Q), M::add(half(z2.i), c.i, Q)};
    const MP d2 = M::add(sqr(D.r, Q), sqr(D.i, Q), Q);
    const CM n = cm_mul(num, CM{D.r, M::neg(D.i)}, Q);
    return {M::div(n.r, d2, Q), M::div(n.i, d2, Q)};
}

// ------------------------------------------------------------------ exp
inline bool exp_special(F x, F y, Comp& re, Comp& im) {
    if (x.isnan()) {
        if (y.iszero()) { re = c_ex(x); im = c_ex(y); return true; }
        re = c_nan(); im = c_nan(); return true;
    }
    if (x.isinf()) {
        if (!sb(x)) {
            if (y.iszero()) { re = c_ex(x); im = c_ex(y); return true; }
            if (y.isinf()) { raise(fe_invalid); re = c_inf(false); im = c_nan(); return true; }
            if (y.isnan()) { re = c_inf(false); im = c_nan(); return true; }
            re = c_inf(cos_neg(y)); im = c_inf(sin_neg(y)); return true;
        }
        if (!y.isfinite()) { re = c_zero(false); im = c_zero(sb(y)); return true; }
        re = c_zero(cos_neg(y)); im = c_zero(sin_neg(y)); return true;
    }
    if (y.isinf()) { raise(fe_invalid); re = c_nan(); im = c_nan(); return true; }
    if (y.isnan()) { re = c_nan(); im = c_nan(); return true; }
    if (y.iszero()) { re = Comp::fn(r_exp, x); im = c_ex(y); return true; }
    return false;
}
inline void exp_kernel(F x, F y, int P, Comp& re, Comp& im) {
    const int Q = P + 32;
    const MP X = mpf(x), Y = mpf(y);
    if (tiny_arg(x, y)) {  // 1 + z exact, rest Σ z^k/k!
        const CM r = tiny_exp_rest(CM{X, Y}, Q);
        re = Comp::ap(trunc(r.r, P), xadd(one(), X));
        im = Comp::ap(trunc(r.i, P), Y);
        return;
    }
    MP s, c;
    msincos(Y, Q, s, c);
    if (gt(x, 190000)) { re = Comp::ap(huge(sgn(c), true)); im = Comp::ap(huge(sgn(s), true)); return; }
    if (gt(-x, 190000)) { re = Comp::ap(huge(sgn(c), false)); im = Comp::ap(huge(sgn(s), false)); return; }
    const MP e = mexp(X, Q);
    im = Comp::ap(trunc(M::mul(e, s, Q), P));
    if (!gt(abs(x), 1) && !gt(abs(y), 1)) {
        MP sh, ch;
        M::sincos(half(Y), Q, sh, ch);
        const MP rest = M::sub(M::mul(mexpm1(X, Q), c, Q), M::ldexp(sqr(sh, Q), 1), Q);
        re = Comp::ap(trunc(rest, P), F::one());
    } else {
        re = Comp::ap(trunc(M::mul(e, c, Q), P));
    }
}

// ------------------------------------------------------------------ sinh / cosh / tanh
inline bool sinh_special(F x, F y, Comp& re, Comp& im) {
    if (x.isnan()) {
        if (y.iszero()) { re = c_ex(x); im = c_ex(y); return true; }
        re = c_nan(); im = c_nan(); return true;
    }
    if (x.isinf()) {
        if (y.iszero()) { re = c_ex(x); im = c_ex(y); return true; }
        if (y.isinf()) { raise(fe_invalid); re = c_ex(x); im = c_nan(); return true; }
        if (y.isnan()) { re = c_ex(x); im = c_nan(); return true; }
        re = c_inf(sb(x) != cos_neg(y)); im = c_inf(sin_neg(y)); return true;
    }
    if (!y.isfinite()) {
        if (y.isinf()) raise(fe_invalid);
        if (x.iszero()) { re = c_ex(x); im = c_nan(); return true; }
        re = c_nan(); im = c_nan(); return true;
    }
    if (y.iszero()) { re = Comp::fn(r_sinh, x); im = c_ex(y); return true; }
    if (x.iszero()) { re = c_zero(sb(x) != cos_neg(y)); im = Comp::fn(r_sin, y); return true; }
    return false;
}
inline void sinh_kernel(F x, F y, int P, Comp& re, Comp& im) {
    const int Q = P + 32;
    const MP X = mpf(x), Y = mpf(y);
    if (tiny_arg(x, y)) {
        const CM r = tiny_sinh_rest(CM{X, Y}, Q);
        re = Comp::ap(trunc(r.r, P), X);
        im = Comp::ap(trunc(r.i, P), Y);
        return;
    }
    MP s, c;
    msincos(Y, Q, s, c);
    if (gt(abs(x), 190000)) {
        re = Comp::ap(huge(sgnf(x) * sgn(c), true)); im = Comp::ap(huge(sgn(s), true)); return;
    }
    re = Comp::ap(trunc(M::mul(msinh(X, Q), c, Q), P));
    im = Comp::ap(trunc(M::mul(mcosh(X, Q), s, Q), P));
}

inline bool cosh_special(F x, F y, Comp& re, Comp& im) {
    if (x.isnan()) {
        if (y.iszero()) { re = c_nan(); im = c_zero(false); return true; }
        re = c_nan(); im = c_nan(); return true;
    }
    if (x.isinf()) {
        if (y.iszero()) { re = c_inf(false); im = c_zero(sb(x) != sb(y)); return true; }
        if (y.isinf()) { raise(fe_invalid); re = c_inf(false); im = c_nan(); return true; }
        if (y.isnan()) { re = c_inf(false); im = c_nan(); return true; }
        re = c_inf(cos_neg(y)); im = c_inf(sb(x) != sin_neg(y)); return true;
    }
    if (!y.isfinite()) {
        if (y.isinf()) raise(fe_invalid);
        if (x.iszero()) { re = c_nan(); im = c_zero(false); return true; }
        re = c_nan(); im = c_nan(); return true;
    }
    if (y.iszero()) { re = Comp::fn(r_cosh, x); im = c_zero(sb(x) != sb(y)); return true; }
    if (x.iszero()) { re = Comp::fn(r_cos, y); im = c_zero(sb(x) != sin_neg(y)); return true; }
    return false;
}
inline void cosh_kernel(F x, F y, int P, Comp& re, Comp& im) {
    const int Q = P + 32;
    const MP X = mpf(x), Y = mpf(y);
    if (tiny_arg(x, y)) {  // 1 + z²/2 exact: 1 + (x² − y²)/2 and x·y
        const CM r = tiny_cosh_rest(CM{X, Y}, Q);
        re = Comp::ap(trunc(r.r, P), xadd(one(), half(xsub(xmul(X, X), xmul(Y, Y)))));
        im = Comp::ap(trunc(r.i, P), xmul(X, Y));
        return;
    }
    MP s, c;
    msincos(Y, Q, s, c);
    if (gt(abs(x), 190000)) {
        re = Comp::ap(huge(sgn(c), true)); im = Comp::ap(huge(sgnf(x) * sgn(s), true)); return;
    }
    im = Comp::ap(trunc(M::mul(msinh(X, Q), s, Q), P));
    if (!gt(abs(x), 1) && !gt(abs(y), 1)) {
        MP sh, ch;
        M::sincos(half(Y), Q, sh, ch);
        const MP rest = M::sub(M::mul(mcosh_m1(X, Q), c, Q), M::ldexp(sqr(sh, Q), 1), Q);
        re = Comp::ap(trunc(rest, P), F::one());
    } else {
        re = Comp::ap(trunc(M::mul(mcosh(X, Q), c, Q), P));
    }
}

inline bool tanh_special(F x, F y, Comp& re, Comp& im) {
    if (x.isnan()) {
        if (y.iszero()) { re = c_ex(x); im = c_ex(y); return true; }
        re = c_nan(); im = c_nan(); return true;
    }
    if (x.isinf()) {
        const F one = sb(x) ? -F::one() : F::one();
        if (!y.isfinite()) { re = c_ex(one); im = c_zero(sb(y)); return true; }
        re = c_ex(one); im = c_zero(y.iszero() ? sb(y) : (sin_neg(y) != cos_neg(y))); return true;
    }
    if (!y.isfinite()) {
        if (y.isinf()) raise(fe_invalid);
        if (x.iszero()) { re = c_ex(x); im = c_nan(); return true; }
        re = c_nan(); im = c_nan(); return true;
    }
    if (y.iszero()) { re = Comp::fn(r_tanh, x); im = c_ex(y); return true; }
    if (x.iszero()) { re = c_ex(x); im = Comp::fn(r_tan, y); return true; }
    return false;
}
inline void tanh_kernel(F x, F y, int P, Comp& re, Comp& im) {
    const int Q = P + 32;
    const MP X = mpf(x), Y = mpf(y);
    if (tiny_arg(x, y)) {
        const CM r = tiny_tanh_rest(CM{X, Y}, Q);
        re = Comp::ap(trunc(r.r, P), X);
        im = Comp::ap(trunc(r.i, P), Y);
        return;
    }
    const int sx = sgnf(x);
    const F one = sx < 0 ? -F::one() : F::one();
    MP s, c;
    msincos(Y, Q, s, c);
    const int c2 = M::cmp_abs(c, s);  // sign of cos 2y = c² − s²
    if (gt(abs(x), 200000)) {
        re = Comp::ap(tiny(-sx * c2), one);
        im = Comp::ap(huge(sgn(s) * sgn(c), false));
        return;
    }
    const MP sh = msinh(X, Q);
    const MP D = M::add(sqr(sh, Q), sqr(c, Q), Q);
    im = Comp::ap(trunc(M::div(M::mul(s, c, Q), D, Q), P));
    if (!(abs(x) < F::one())) {
        if (gt(abs(x), 4000)) { re = Comp::ap(tiny(-sx * c2), one); return; }
        // Re = sign(x)·(1 − δ), δ = (cos 2y + e^{−2|x|}) / (2D)
        const MP cos2y = M::mul(M::sub(c, s, Q), M::add(c, s, Q), Q);
        const MP e2 = mexp(M::neg(M::ldexp(absv(X), 1)), Q);
        MP rest = M::div(M::add(cos2y, e2, Q), M::ldexp(D, 1), Q);
        if (sx > 0) rest = M::neg(rest);
        re = Comp::ap(trunc(rest, P), one);
        return;
    }
    re = Comp::ap(trunc(M::div(M::mul(sh, mcosh(X, Q), Q), D, Q), P));
}

// ------------------------------------------------------------------ asin / acos (Hull, Fairgrieve, Tang)
inline void asin_acos_kernel(F x, F y, int P, bool want_acos, Comp& re, Comp& im) {
    const int Q = P + 64;
    const MP X = mpf(x), Y = mpf(y), ax = absv(X), ay = absv(Y);
    if (tiny_arg(x, y)) {  // asin z = z + rest; acos z = π/2 − asin z
        const CM r = tiny_asin_rest(CM{X, Y}, Q);
        if (!want_acos) {
            re = Comp::ap(trunc(r.r, P), X);
            im = Comp::ap(trunc(r.i, P), Y);
        } else {
            re = Comp::ap(trunc(M::sub(half(M::pi(Q)), M::add(X, r.r, Q), Q), P));
            im = Comp::ap(trunc(M::neg(r.i), P), M::neg(Y));
        }
        return;
    }
    const MP ay2 = trunc(xmul(ay, ay), Q);
    const MP xp1 = xadd(ax, one()), xm1 = xsub(ax, one());
    const MP r = M::sqrt(trunc(xadd(xmul(xp1, xp1), xmul(ay, ay)), Q), Q);
    const MP s = M::sqrt(trunc(xadd(xmul(xm1, xm1), xmul(ay, ay)), Q), Q);
    const MP A = half(M::add(r, s, Q));
    const MP rx1 = M::add(r, trunc(xp1, Q), Q);
    const MP t1 = M::div(ay2, rx1, Q);
    MP Am1, D;
    const MP Apx = M::add(A, ax, Q);
    if (M::cmp_abs(ax, one()) < 0) {
        const MP t2 = M::add(s, trunc(xsub(one(), ax), Q), Q);
        Am1 = half(M::add(t1, M::div(ay2, t2, Q), Q));
        D = M::sqrt(half(M::mul(Apx, M::add(t1, t2, Q), Q)), Q);
    } else {
        const MP t2 = M::add(s, trunc(xm1, Q), Q);
        Am1 = half(M::add(t1, t2, Q));
        D = M::mul(ay, M::sqrt(half(M::add(M::div(Apx, rx1, Q), M::div(Apx, t2, Q), Q)), Q), Q);
    }
    const MP L = mlog1p(M::add(Am1, M::sqrt(M::mul(Am1, M::add(A, one(), Q), Q), Q), Q), Q);
    if (!want_acos) {
        re = Comp::ap(trunc(signed_(matan2(ax, D, Q), sb(x)), P));
        im = Comp::ap(trunc(signed_(L, sb(y)), P));
    } else {
        re = Comp::ap(trunc(matan2(D, X, Q), P));
        im = Comp::ap(trunc(signed_(L, !sb(y)), P));
    }
}

inline bool acos_special(F x, F y, Comp& re, Comp& im) {
    if (x.isnan()) {
        if (y.isinf()) { re = c_nan(); im = c_inf(!sb(y)); return true; }
        re = c_nan(); im = c_nan(); return true;
    }
    if (y.isnan()) {
        if (x.isinf()) { re = c_nan(); im = c_inf(true); return true; }
        if (x.iszero()) { re = c_pio2(false); im = c_nan(); return true; }
        re = c_nan(); im = c_nan(); return true;
    }
    if (x.isinf() || y.isinf()) {
        im = c_inf(!sb(y));
        if (y.isinf()) re = x.isinf() ? (sb(x) ? c_3pio4(false) : c_pio4(false)) : c_pio2(false);
        else re = sb(x) ? c_pi(false) : c_zero(false);
        return true;
    }
    if (y.iszero()) {
        if (!gt(abs(x), 1)) { re = Comp::fn(r_acos, x); im = c_zero(!sb(y)); return true; }
        re = sb(x) ? c_pi(false) : c_zero(false);
        im = Comp::fn(r_acosh, abs(x), !sb(y));
        return true;
    }
    if (x.iszero()) { re = c_pio2(false); im = Comp::fn(r_asinh, y, true); return true; }
    return false;
}
inline void acos_kernel(F x, F y, int P, Comp& re, Comp& im) { asin_acos_kernel(x, y, P, true, re, im); }

inline bool acosh_special(F x, F y, Comp& re, Comp& im) {
    if (x.isnan()) {
        if (y.isinf()) { re = c_inf(false); im = c_nan(); return true; }
        re = c_nan(); im = c_nan(); return true;
    }
    if (y.isnan()) {
        if (x.isinf()) { re = c_inf(false); im = c_nan(); return true; }
        re = c_nan(); im = c_nan(); return true;
    }
    if (x.isinf() || y.isinf()) {
        re = c_inf(false);
        if (y.isinf()) im = x.isinf() ? (sb(x) ? c_3pio4(sb(y)) : c_pio4(sb(y))) : c_pio2(sb(y));
        else im = sb(x) ? c_pi(sb(y)) : c_zero(sb(y));
        return true;
    }
    if (y.iszero()) {
        if (x.iszero()) { re = c_zero(false); im = c_pio2(sb(y)); return true; }
        if (!(x < F::one())) { re = Comp::fn(r_acosh, x); im = c_ex(y); return true; }
        if (x > -F::one()) { re = c_zero(false); im = Comp::fn(r_acos, x, sb(y)); return true; }
        re = Comp::fn(r_acosh, -x); im = c_pi(sb(y)); return true;
    }
    return false;
}
inline void acosh_kernel(F x, F y, int P, Comp& re, Comp& im) {
    asin_acos_kernel(x, y, P, true, re, im);
    apply(sb(y) ? Out::neg_i : Out::pos_i, re, im);
}

inline bool asinh_special(F x, F y, Comp& re, Comp& im) {
    if (x.isnan()) {
        if (y.iszero()) { re = c_ex(x); im = c_ex(y); return true; }
        if (y.isinf()) { re = c_inf(false); im = c_nan(); return true; }
        re = c_nan(); im = c_nan(); return true;
    }
    if (y.isnan()) {
        if (x.isinf()) { re = c_ex(x); im = c_nan(); return true; }
        re = c_nan(); im = c_nan(); return true;
    }
    if (x.isinf() || y.isinf()) {
        re = c_inf(sb(x));
        if (y.isinf()) im = x.isinf() ? c_pio4(sb(y)) : c_pio2(sb(y));
        else im = c_zero(sb(y));
        return true;
    }
    if (y.iszero()) { re = Comp::fn(r_asinh, x); im = c_ex(y); return true; }
    if (x.iszero()) {
        if (!gt(abs(y), 1)) { re = c_ex(x); im = Comp::fn(r_asin, y); return true; }
        re = Comp::fn(r_acosh, abs(y), sb(x)); im = c_pio2(sb(y)); return true;
    }
    return false;
}
/// asinh(z) = −i asin(iz)
inline void asinh_kernel(F x, F y, int P, Comp& re, Comp& im) {
    asin_acos_kernel(-y, x, P, false, re, im);
    apply(Out::neg_i, re, im);
}

inline bool atanh_special(F x, F y, Comp& re, Comp& im) {
    if (x.isnan()) {
        if (y.isinf()) { re = c_zero(false); im = c_pio2(sb(y)); return true; }
        re = c_nan(); im = c_nan(); return true;
    }
    if (y.isnan()) {
        if (x.isinf()) { re = c_zero(sb(x)); im = c_nan(); return true; }
        if (x.iszero()) { re = c_ex(x); im = c_nan(); return true; }
        re = c_nan(); im = c_nan(); return true;
    }
    if (x.isinf() || y.isinf()) { re = c_zero(sb(x)); im = c_pio2(sb(y)); return true; }
    if (y.iszero()) {
        const F ax = abs(x);
        if (ax < F::one()) { re = Comp::fn(r_atanh, x); im = c_ex(y); return true; }
        if (ax == F::one()) { raise(fe_divbyzero); re = c_inf(sb(x)); im = c_ex(y); return true; }
        return false;
    }
    if (x.iszero()) { re = c_ex(x); im = Comp::fn(r_atan, y); return true; }
    return false;
}
inline void atanh_kernel(F x, F y, int P, Comp& re, Comp& im) {
    const int Q = P + 64;
    const MP X = mpf(x), Y = mpf(y);
    if (tiny_arg(x, y)) {
        const CM r = tiny_atanh_rest(CM{X, Y}, Q);
        re = Comp::ap(trunc(r.r, P), X);
        im = Comp::ap(trunc(r.i, P), Y);
        return;
    }
    const MP omx = xsub(one(), X), opx = xadd(one(), X), y2 = xmul(Y, Y);
    const MP den = xadd(xmul(omx, omx), y2);   // (1 − x)² + y²
    const MP num = xadd(xmul(opx, opx), y2);   // (1 + x)² + y² = den + 4x
    const MP q = M::div(M::ldexp(X, 2), trunc(den, Q), Q);
    MP l;
    if (q.top() <= -1) l = mlog1p(q, Q);
    else l = M::sub(mlog(num, Q), mlog(den, Q), Q);
    re = Comp::ap(trunc(M::ldexp(l, -2), P));
    if (y.iszero()) { im = c_pio2(sb(y)); return; }
    const MP n1 = xsub(one(), xadd(xmul(X, X), y2));  // 1 − x² − y²
    im = Comp::ap(trunc(half(matan2(M::ldexp(Y, 1), trunc(n1, Q), Q)), P));
}

// ------------------------------------------------------------------ log
inline bool log_special(F x, F y, Comp& re, Comp& im) {
    if (x.isnan() || y.isnan()) {
        if (x.isinf() || y.isinf()) { re = c_inf(false); im = c_nan(); return true; }
        re = c_nan(); im = c_nan(); return true;
    }
    if (x.isinf() || y.isinf()) { re = c_inf(false); im = Comp::fn2(r_atan2, y, x); return true; }
    if (x.iszero() && y.iszero()) {
        raise(fe_divbyzero);
        re = c_inf(true); im = Comp::fn2(r_atan2, y, x); return true;
    }
    if (y.iszero()) { re = Comp::fn(r_log, abs(x)); im = Comp::fn2(r_atan2, y, x); return true; }
    if (x.iszero()) { re = Comp::fn(r_log, abs(y)); im = Comp::fn2(r_atan2, y, x); return true; }
    return false;
}
/// ln |z| = ½ ln(x² + y²), through log1p of the exact x² + y² − 1 near |z| = 1.
[[nodiscard]] inline MP log_abs(const MP& N, int Q) {
    const MP t = xsub(N, one());
    if (t.top() <= -1) return half(mlog1p(trunc(t, Q), Q));
    return half(mlog(N, Q));
}

/// x² + y² = 10^n exactly (n >= 0)? (Then log10|z| = n/2.)
[[nodiscard]] inline bool is_pow10(const MP& N, long& n) {
    if (N.zero() || N.neg) return false;
    int tz = 0;
    while (!N.m.bit(tz)) ++tz;
    n = N.e + tz;
    if (n < 0) return false;
    BigInt o = N.m;
    o.shr(tz);
    if (static_cast<double>(n) * 2.3219280948873622 > o.bit_length() + 1.0) return false;
    BigInt p = M::from_u64(1).m;
    p.mul_pow5(static_cast<unsigned>(n));
    return p == o;
}

} // namespace detail_cx

// ============================================================ public functions

namespace detail_cx {

[[nodiscard]] inline complex mul_special(F a, F b, F c, F d) {
    const F inf = pinf();
    F ac = a * c, bd = b * d, ad = a * d, bc = b * c;
    F x = ac - bd, y = ad + bc;
    if (x.isnan() && y.isnan()) {
        bool recalc = false;
        auto box = [](F v) { return copysign(v.isinf() ? F::one() : F::zero(), v); };
        auto unnan = [](F v) { return v.isnan() ? copysign(F::zero(), v) : v; };
        if (a.isinf() || b.isinf()) { a = box(a); b = box(b); c = unnan(c); d = unnan(d); recalc = true; }
        if (c.isinf() || d.isinf()) { c = box(c); d = box(d); a = unnan(a); b = unnan(b); recalc = true; }
        if (!recalc && (ac.isinf() || bd.isinf() || ad.isinf() || bc.isinf())) {
            a = unnan(a); b = unnan(b); c = unnan(c); d = unnan(d); recalc = true;
        }
        if (recalc) { x = inf * (a * c - b * d); y = inf * (a * d + b * c); }
    }
    return {x, y};
}

[[nodiscard]] inline complex div_special(F a, F b, F c, F d) {
    const F inf = pinf();
    const F logbw = logb(fmax(abs(c), abs(d)));
    int ilogbw = 0;
    if (logbw.isfinite()) {
        ilogbw = static_cast<int>(logbw.to_double());
        c = ldexp(c, -ilogbw);
        d = ldexp(d, -ilogbw);
    }
    const F denom = c * c + d * d;
    F x = ldexp((a * c + b * d) / denom, -ilogbw);
    F y = ldexp((b * c - a * d) / denom, -ilogbw);
    if (x.isnan() && y.isnan()) {
        auto box = [](F v) { return copysign(v.isinf() ? F::one() : F::zero(), v); };
        if (denom.iszero() && (!a.isnan() || !b.isnan())) {
            x = copysign(inf, c) * a;
            y = copysign(inf, c) * b;
        } else if ((a.isinf() || b.isinf()) && c.isfinite() && d.isfinite()) {
            a = box(a); b = box(b);
            x = inf * (a * c + b * d);
            y = inf * (b * c - a * d);
        } else if (logbw.isinf() && a.isfinite() && b.isfinite()) {
            c = box(c); d = box(d);
            x = F::zero() * (a * c + b * d);
            y = F::zero() * (b * c - a * d);
        }
    }
    return {x, y};
}

/// Plan for the exact quotient n/d (d > 0, n != 0).
[[nodiscard]] inline Plan quotient_plan(const MP& n, const MP& d) {
    Plan p;
    for (int P = 320;; P *= 2) {
        const MP q = M::div(trunc(n, P + 32), trunc(d, P + 32), P);
        auto oracle = [&](const MP& B) { return sgn(xsub(n, xmul(B, d))); };
        if (plan_algebraic(q, q.top() - P + 6, oracle, p)) return p;
    }
}

[[nodiscard]] inline bool prod_neg(F a, F b) { return a.signbit() != b.signbit(); }
[[nodiscard]] inline bool prod_zero(F a, F b) { return a.iszero() || b.iszero(); }

} // namespace detail_cx

/// Componentwise correctly rounded: Re = round(ac − bd), Im = round(ad + bc)
/// from the exact products (finite operands); otherwise Annex G G.5.1.
[[nodiscard]] inline complex operator*(const complex& z, const complex& w) {
    using namespace detail_cx;
    const F a = z.real(), b = z.imag(), c = w.real(), d = w.imag();
    if (!(a.isfinite() && b.isfinite() && c.isfinite() && d.isfinite())) return mul_special(a, b, c, d);
    CxEnv env;
    const MP A = mpf(a), B = mpf(b), C = mpf(c), D = mpf(d);
    const MP re = xsub(xmul(A, C), xmul(B, D)), im = xadd(xmul(A, D), xmul(B, C));
    Plan pr, pi;
    if (!re.zero()) pr = plan_from(re, 0);
    if (!im.zero()) pi = plan_from(im, 0);
    env.restore();
    const round_mode m = env.mode();
    const F x = re.zero() ? F::zero(zero_sum_negative(prod_zero(a, c), prod_neg(a, c), prod_zero(b, d), !prod_neg(b, d), m))
                          : pr.round(m);
    const F y = im.zero() ? F::zero(zero_sum_negative(prod_zero(a, d), prod_neg(a, d), prod_zero(b, c), prod_neg(b, c), m))
                          : pi.round(m);
    return {x, y};
}

/// Componentwise correctly rounded: Re = round((ac + bd)/(c² + d²)),
/// Im = round((bc − ad)/(c² + d²)) exactly (finite operands, nonzero divisor);
/// otherwise Annex G G.5.1 (division by zero also raises divide-by-zero).
[[nodiscard]] inline complex operator/(const complex& z, const complex& w) {
    using namespace detail_cx;
    const F a = z.real(), b = z.imag(), c = w.real(), d = w.imag();
    const bool fin = a.isfinite() && b.isfinite() && c.isfinite() && d.isfinite();
    if (fin && c.iszero() && d.iszero()) {
        if (a.iszero() && b.iszero()) { raise(fe_invalid); return {qnan(), qnan()}; }
        raise(fe_divbyzero);
        const F inf = copysign(pinf(), c);
        return {a.iszero() ? qnan() : inf * a, b.iszero() ? qnan() : inf * b};
    }
    if (!fin) return div_special(a, b, c, d);
    CxEnv env;
    const MP A = mpf(a), B = mpf(b), C = mpf(c), D = mpf(d);
    const MP nr = xadd(xmul(A, C), xmul(B, D)), ni = xsub(xmul(B, C), xmul(A, D));
    const MP den = xadd(xmul(C, C), xmul(D, D));
    Plan pr, pi;
    if (!nr.zero()) pr = quotient_plan(nr, den);
    if (!ni.zero()) pi = quotient_plan(ni, den);
    env.restore();
    const round_mode m = env.mode();
    const F x = nr.zero() ? F::zero(zero_sum_negative(prod_zero(a, c), prod_neg(a, c), prod_zero(b, d), prod_neg(b, d), m))
                          : pr.round(m);
    const F y = ni.zero() ? F::zero(zero_sum_negative(prod_zero(b, c), prod_neg(b, c), prod_zero(a, d), !prod_neg(a, d), m))
                          : pi.round(m);
    return {x, y};
}

/// |z| = hypot(re, im), correctly rounded; +∞ if either part is infinite (even with a NaN).
[[nodiscard]] inline float256 abs(const complex& z) {
    using namespace detail_cx;
    F x = z.real(), y = z.imag();
    if (x.isinf() || y.isinf()) { quiet(x); quiet(y); return pinf(); }
    if (x.isnan() || y.isnan()) { quiet(x); quiet(y); return qnan(); }
    if (x.iszero()) return abs(y);
    if (y.iszero()) return abs(x);
    CxEnv env;
    const MP X = mpf(x), Y = mpf(y);
    const MP N = xadd(xmul(X, X), xmul(Y, Y));
    Plan p;
    for (int P = 320;; P *= 2) {
        const MP v = M::sqrt(trunc(N, P + 32), P);
        auto oracle = [&](const MP& B) { return sgn(xsub(N, xmul(B, B))); };
        if (plan_algebraic(v, v.top() - P + 6, oracle, p)) break;
    }
    env.restore();
    return p.round(env.mode());
}

/// re² + im², correctly rounded; +∞ if either part is infinite.
[[nodiscard]] inline float256 norm(const complex& z) {
    using namespace detail_cx;
    F x = z.real(), y = z.imag();
    if (x.isinf() || y.isinf()) { quiet(x); quiet(y); return pinf(); }
    if (x.isnan() || y.isnan()) { quiet(x); quiet(y); return qnan(); }
    CxEnv env;
    const MP X = mpf(x), Y = mpf(y);
    const MP N = xadd(xmul(X, X), xmul(Y, Y));
    if (N.zero()) return F::zero();
    const Plan p = plan_from(N, 0);
    env.restore();
    return p.round(env.mode());
}

/// atan2(im, re), correctly rounded.
[[nodiscard]] inline float256 arg(const complex& z) { return atan2(z.imag(), z.real()); }

/// rho·cos θ + i·rho·sin θ, componentwise correctly rounded.
/// polar(rho, ±0) = rho ± i0·sign(rho); NaN or infinite θ gives NaN + iNaN.
[[nodiscard]] inline complex polar(float256 rho, float256 theta = float256()) {
    using namespace detail_cx;
    quiet(rho);
    quiet(theta);
    if (rho.isnan() || theta.isnan()) return {qnan(), qnan()};
    if (theta.isinf()) { raise(fe_invalid); return {qnan(), qnan()}; }
    if (theta.iszero()) return {rho, F::zero(sb(rho) != sb(theta))};
    if (rho.isinf()) return {F::inf(sb(rho) != cos_neg(theta)), F::inf(sb(rho) != sin_neg(theta))};
    if (rho.iszero()) return {F::zero(sb(rho) != cos_neg(theta)), F::zero(sb(rho) != sin_neg(theta))};
    const bool tiny_th = ilogb(theta) < -64;
    const auto r = ziv([&](int P, Comp& re, Comp& im) {
        const int Q = P + 32;
        const MP R = mpf(rho), T = mpf(theta);
        if (tiny_th) {  // ρ − 2ρ·sin²(θ/2) and ρθ + ρ(sin θ − θ)
            MP sh, ch;
            M::sincos(half(T), Q, sh, ch);
            re = Comp::ap(trunc(M::neg(M::ldexp(M::mul(R, sqr(sh, Q), Q), 1)), P), R);
            const CM r = tiny_sinh_rest(CM{MP{}, T}, Q);  // sinh(iθ) − iθ = i(sin θ − θ)
            im = Comp::ap(trunc(M::mul(R, r.i, Q), P), xmul(R, T));
            return;
        }
        MP s, c;
        msincos(T, Q, s, c);
        re = Comp::ap(trunc(M::mul(R, c, Q), P));
        im = Comp::ap(trunc(M::mul(R, s, Q), P));
    });
    return {r[0], r[1]};
}

/// e^z (Annex G cexp).
[[nodiscard]] inline complex exp(const complex& z) {
    using namespace detail_cx;
    F x = z.real(), y = z.imag();
    quiet(x); quiet(y);
    return run(x, y, exp_special, exp_kernel);
}

/// Natural logarithm, branch cut along the negative real axis, Im in [−π, π] (Annex G clog).
[[nodiscard]] inline complex log(const complex& z) {
    using namespace detail_cx;
    F x = z.real(), y = z.imag();
    quiet(x); quiet(y);
    Comp re, im;
    if (log_special(x, y, re, im)) {
        const F a = value(re);
        const F b = value(im);
        return {a, b};
    }
    const MP X = mpf(x), Y = mpf(y);
    const MP N = xadd(xmul(X, X), xmul(Y, Y));
    const MP t = xsub(N, one());
    const bool tiny_t = !t.zero() && t.top() < -64;
    const auto r = ziv([&](int P, Comp& a, Comp& b) {
        if (tiny_t) {  // ½ log1p(t) = t/2 − t²/4 (exact) + ½ Σ_{k>=3} (−1)^(k+1) t^k/k
            const int Q = P + 32;
            const MP tq = trunc(t, Q);
            MP sum, pw = M::mul(M::mul(tq, tq, Q), tq, Q);
            for (long k = 3; !pw.zero() && pw.top() > 3 * t.top() - Q - 8; ++k) {
                const MP term = M::div_u(pw, static_cast<std::uint32_t>(k), Q);
                sum = (k & 1) ? M::add(sum, term, Q) : M::sub(sum, term, Q);
                pw = M::mul(pw, tq, Q);
            }
            a = Comp::ap(trunc(half(sum), P), xsub(half(t), M::ldexp(xmul(t, t), -2)));
        } else {
            a = Comp::ap(trunc(log_abs(N, P + 32), P));
        }
        b = Comp::fn2(r_atan2, y, x);
    });
    return {r[0], r[1]};
}

/// log(z)/ln 10 (C++ std::log10), the same branch cut as log.
[[nodiscard]] inline complex log10(const complex& z) {
    using namespace detail_cx;
    F x = z.real(), y = z.imag();
    quiet(x); quiet(y);
    // θ = arg z as a multiple k·π/4 when it is one (special values, axes)
    int k = 99;
    Comp re;
    bool re_set = false;
    if (x.isnan() || y.isnan()) {
        if (x.isinf() || y.isinf()) return {pinf(), qnan()};
        return {qnan(), qnan()};
    }
    if (x.isinf() || y.isinf()) {
        re = c_inf(false); re_set = true;
        if (x.isinf() && y.isinf()) k = sb(x) ? 3 : 1;
        else if (y.isinf()) k = 2;
        else k = sb(x) ? 4 : 0;
        if (sb(y)) k = -k;
    } else if (x.iszero() && y.iszero()) {
        raise(fe_divbyzero);
        re = c_inf(true); re_set = true;
        k = sb(x) ? (sb(y) ? -4 : 4) : 0;
    } else if (y.iszero()) {
        re = Comp::fn(r_log10, abs(x)); re_set = true;
        k = sb(x) ? (sb(y) ? -4 : 4) : 0;
    } else if (x.iszero()) {
        re = Comp::fn(r_log10, abs(y)); re_set = true;
        k = sb(y) ? -2 : 2;
    }
    if (k == 0 && re_set) {
        const F a = value(re);
        return {a, F::zero(sb(y))};
    }
    MP N;
    long n10 = 0;
    bool exact_re = false;
    if (!re_set) {
        const MP X = mpf(x), Y = mpf(y);
        N = xadd(xmul(X, X), xmul(Y, Y));
        exact_re = is_pow10(N, n10);
    }
    const auto r = ziv([&](int P, Comp& a, Comp& b) {
        const int Q = P + 32;
        const MP l10 = ln10(Q);
        if (re_set) a = re;
        else if (exact_re) a = Comp::ex(ldexp(F::from_int(n10), -1));
        else a = Comp::ap(trunc(M::div(log_abs(N, Q), l10, Q), P));
        if (k != 99) {
            MP th = M::div_u(M::mul_u(M::pi(Q), static_cast<std::uint32_t>(k < 0 ? -k : k), Q), 4, Q);
            b = Comp::ap(trunc(signed_(M::div(th, l10, Q), k < 0), P));
        } else {
            b = Comp::ap(trunc(M::div(matan2(mpf(y), mpf(x), Q), l10, Q), P));
        }
    });
    return {r[0], r[1]};
}

namespace detail_cx {
/// Componentwise correctly rounded square root (algebraic: exact boundary test).
[[nodiscard]] inline complex sqrt_finite(F x, F y) {
    CxEnv env;
    const MP X = mpf(x), Y = mpf(y);
    const MP N = xadd(xmul(X, X), xmul(Y, Y));
    // c = sqrt((|z| + s·x)/2) vs B >= 0: |z| vs K = 2B² − s·x
    auto oracle_abs = [&](const MP& Bm, int s) {
        const MP B = absv(Bm);
        const MP K = xsub(M::ldexp(xmul(B, B), 1), s > 0 ? X : M::neg(X));
        if (K.neg) return 1;
        return sgn(xsub(N, xmul(K, K)));
    };
    Plan pu, pv;
    bool du = false, dv = false;
    for (int P = 320; !(du && dv); P *= 2) {
        const int Q = P + 32;
        const MP az = M::sqrt(trunc(N, Q), Q);
        MP u, v;
        if (!sb(x)) {
            u = M::sqrt(half(M::add(az, X, Q)), Q);
            v = M::div(absv(Y), M::ldexp(u, 1), Q);
        } else {
            v = M::sqrt(half(M::sub(az, X, Q)), Q);
            u = M::div(absv(Y), M::ldexp(v, 1), Q);
        }
        u = trunc(u, P);
        v = trunc(signed_(v, sb(y)), P);
        if (!du) du = plan_algebraic(u, u.top() - P + 8, [&](const MP& B) { return oracle_abs(B, 1); }, pu);
        if (!dv) dv = plan_algebraic(v, v.top() - P + 8, [&](const MP& B) {
            const int s = oracle_abs(B, -1);
            return sb(y) ? -s : s;
        }, pv);
    }
    env.restore();
    const round_mode m = env.mode();
    const F a = pu.round(m);
    const F b = pv.round(m);
    return {a, b};
}
} // namespace detail_cx

/// Principal square root, Re >= 0, branch cut along the negative real axis (Annex G csqrt).
[[nodiscard]] inline complex sqrt(const complex& z) {
    using namespace detail_cx;
    F x = z.real(), y = z.imag();
    quiet(x); quiet(y);
    if (y.isinf()) return {pinf(), y};
    if (x.isnan()) return {qnan(), qnan()};
    if (x.isinf()) {
        if (!sb(x)) return {pinf(), y.isnan() ? qnan() : F::zero(sb(y))};
        if (y.isnan()) return {qnan(), pinf()};
        return {F::zero(), F::inf(sb(y))};
    }
    if (y.isnan()) return {qnan(), qnan()};
    if (y.iszero()) {
        if (x.iszero()) return {F::zero(), y};
        if (!sb(x)) { const F a = sqrt(x); return {a, y}; }
        return {F::zero(), eval_realfn(Comp::fn(r_sqrt, -x, sb(y)))};
    }
    return sqrt_finite(x, y);
}

/// sinh z (Annex G csinh).
[[nodiscard]] inline complex sinh(const complex& z) {
    using namespace detail_cx;
    F x = z.real(), y = z.imag();
    quiet(x); quiet(y);
    return run(x, y, sinh_special, sinh_kernel);
}
/// cosh z (Annex G ccosh).
[[nodiscard]] inline complex cosh(const complex& z) {
    using namespace detail_cx;
    F x = z.real(), y = z.imag();
    quiet(x); quiet(y);
    return run(x, y, cosh_special, cosh_kernel);
}
/// tanh z (Annex G ctanh, C17 for ±0 + i∞ and ±0 + iNaN).
[[nodiscard]] inline complex tanh(const complex& z) {
    using namespace detail_cx;
    F x = z.real(), y = z.imag();
    quiet(x); quiet(y);
    return run(x, y, tanh_special, tanh_kernel);
}
/// sin z = −i sinh(iz).
[[nodiscard]] inline complex sin(const complex& z) {
    using namespace detail_cx;
    F x = z.real(), y = z.imag();
    quiet(x); quiet(y);
    return run(-y, x, sinh_special, sinh_kernel, Out::neg_i);
}
/// cos z = cosh(iz).
[[nodiscard]] inline complex cos(const complex& z) {
    using namespace detail_cx;
    F x = z.real(), y = z.imag();
    quiet(x); quiet(y);
    return run(-y, x, cosh_special, cosh_kernel);
}
/// tan z = −i tanh(iz).
[[nodiscard]] inline complex tan(const complex& z) {
    using namespace detail_cx;
    F x = z.real(), y = z.imag();
    quiet(x); quiet(y);
    return run(-y, x, tanh_special, tanh_kernel, Out::neg_i);
}
/// acos z, branch cuts (−∞, −1) and (1, ∞) on the real axis, Re in [0, π] (Annex G cacos).
[[nodiscard]] inline complex acos(const complex& z) {
    using namespace detail_cx;
    F x = z.real(), y = z.imag();
    quiet(x); quiet(y);
    return run(x, y, acos_special, acos_kernel);
}
/// asin z = −i asinh(iz), Re in [−π/2, π/2].
[[nodiscard]] inline complex asin(const complex& z) {
    using namespace detail_cx;
    F x = z.real(), y = z.imag();
    quiet(x); quiet(y);
    // the special values of casinh(iz), mapped; the finite values of asin directly
    Comp re, im;
    if (asinh_special(-y, x, re, im)) {
        apply(Out::neg_i, re, im);
        const F a = value(re);
        const F b = value(im);
        return {a, b};
    }
    const auto r = ziv([&](int P, Comp& a, Comp& b) { asin_acos_kernel(x, y, P, false, a, b); });
    return {r[0], r[1]};
}
/// atan z = −i atanh(iz), branch cuts (−i∞, −i) and (i, i∞).
[[nodiscard]] inline complex atan(const complex& z) {
    using namespace detail_cx;
    F x = z.real(), y = z.imag();
    quiet(x); quiet(y);
    return run(-y, x, atanh_special, atanh_kernel, Out::neg_i);
}
/// asinh z, branch cuts (−i∞, −i) and (i, i∞) (Annex G casinh).
[[nodiscard]] inline complex asinh(const complex& z) {
    using namespace detail_cx;
    F x = z.real(), y = z.imag();
    quiet(x); quiet(y);
    return run(x, y, asinh_special, asinh_kernel);
}
/// acosh z, branch cut (−∞, 1), Re >= 0 (Annex G cacosh).
[[nodiscard]] inline complex acosh(const complex& z) {
    using namespace detail_cx;
    F x = z.real(), y = z.imag();
    quiet(x); quiet(y);
    return run(x, y, acosh_special, acosh_kernel);
}
/// atanh z, branch cuts (−∞, −1] and [1, ∞) (Annex G catanh).
[[nodiscard]] inline complex atanh(const complex& z) {
    using namespace detail_cx;
    F x = z.real(), y = z.imag();
    quiet(x); quiet(y);
    return run(x, y, atanh_special, atanh_kernel);
}

namespace detail_cx {

/// z^n for |n| <= 64 from the exact Gaussian power (componentwise correctly rounded).
[[nodiscard]] inline complex pow_int_exact(F x, F y, long n) {
    CxEnv env;
    MP br = mpf(x), bi = mpf(y), rr = one(), ri;
    for (long e = n < 0 ? -n : n; e; e >>= 1) {
        if (e & 1) {
            const MP t = xsub(xmul(rr, br), xmul(ri, bi));
            ri = xadd(xmul(rr, bi), xmul(ri, br));
            rr = t;
        }
        if (e > 1) {
            const MP t = xsub(xmul(br, br), xmul(bi, bi));
            bi = M::ldexp(xmul(br, bi), 1);
            br = t;
        }
    }
    Plan pr, pi;
    MP nr = rr, ni = ri;
    if (n > 0) {
        if (!nr.zero()) pr = plan_from(nr, 0);
        if (!ni.zero()) pi = plan_from(ni, 0);
    } else {
        const MP den = xadd(xmul(rr, rr), xmul(ri, ri));
        ni = M::neg(ri);
        if (!nr.zero()) pr = quotient_plan(nr, den);
        if (!ni.zero()) pi = quotient_plan(ni, den);
    }
    env.restore();
    const round_mode m = env.mode();
    const F a = nr.zero() ? F::zero(m == round_mode::downward) : pr.round(m);
    const F b = ni.zero() ? F::zero(m == round_mode::downward) : pi.round(m);
    return {a, b};
}

inline constexpr int kPowExtraCap = 4096;

/// exp(w·log z) for finite z != 0 and finite w (the Ziv loop; a value still on a
/// rounding boundary at the last precision is taken as exact).
[[nodiscard]] inline complex pow_general(F x, F y, F wr, F wi) {
    const MP X = mpf(x), Y = mpf(y), WR = mpf(wr), WI = mpf(wi);
    const MP N = xadd(xmul(X, X), xmul(Y, Y));
    long wt = std::max(WR.zero() ? LONG_MIN / 2 : WR.top(), WI.zero() ? LONG_MIN / 2 : WI.top());
    const int extra = static_cast<int>(std::clamp(wt + 24, 0L, static_cast<long>(kPowExtraCap)));
    const auto r = ziv([&](int P, Comp& re, Comp& im) {
        const int Q = P + 64 + extra;
        const MP Lr = log_abs(N, Q), Li = matan2(Y, X, Q);
        const MP ur = M::sub(M::mul(WR, Lr, Q), M::mul(WI, Li, Q), Q);
        const MP ui = M::add(M::mul(WR, Li, Q), M::mul(WI, Lr, Q), Q);
        MP s, c;
        msincos(ui, P + 32, s, c);
        const double urd = ur.zero() ? 0.0 : (ur.top() > 40 ? (ur.neg ? -1e30 : 1e30) : detail_mp::to_dd(ur).hi.to_double());
        if (urd > 190000) { re = Comp::ap(huge(sgn(c), true)); im = Comp::ap(huge(sgn(s), true)); return; }
        if (urd < -190000) { re = Comp::ap(huge(sgn(c), false)); im = Comp::ap(huge(sgn(s), false)); return; }
        const int R = P + 32;
        const MP e = mexp(ur, R);
        im = Comp::ap(trunc(M::mul(e, s, R), P));
        if (ur.top() <= 0 && ui.top() <= 0) {
            MP sh, ch;
            M::sincos(half(ui), R, sh, ch);
            const MP rest = M::sub(M::mul(mexpm1(ur, R), c, R), M::ldexp(sqr(sh, R), 1), R);
            re = Comp::ap(trunc(rest, P), F::one());
        } else {
            re = Comp::ap(trunc(M::mul(e, c, R), P));
        }
    });
    return {r[0], r[1]};
}

[[nodiscard]] inline complex pow_impl(complex z, F wr, F wi, bool real_w) {
    F x = z.real(), y = z.imag();
    quiet(x); quiet(y); quiet(wr); quiet(wi);
    if (wr.iszero() && wi.iszero()) return {F::one(), F::zero()};
    if (x.isnan() || y.isnan() || wr.isnan() || wi.isnan()) return {qnan(), qnan()};
    if (x.iszero() && y.iszero()) {
        if (wr > F::zero()) return {F::zero(), F::zero()};
        if (wi.iszero()) { raise(fe_divbyzero); return {pinf(), F::zero()}; }
        raise(fe_invalid);
        return {qnan(), qnan()};
    }
    if (!(x.isfinite() && y.isfinite() && wr.isfinite() && wi.isfinite())) {
        const complex l = log(complex(x, y));
        return exp(real_w ? l * wr : complex(wr, wi) * l);
    }
    if (wi.iszero()) {
        const F t = wr;
        if (t == F::half()) return sqrt(complex(x, y));
        if (y.iszero()) {
            if (!sb(x)) {
                const F a = pow(x, t);
                return {a, F::zero(sb(t) != sb(y))};
            }
            if (t == trunc(t)) {
                const bool odd = detail_math::is_odd_integer(t);
                const F a = eval_realfn(Comp::fn2(r_pow, -x, t, odd));
                return {a, F::zero((sb(t) != !odd) != sb(y))};
            }
            const F t2 = ldexp(t, 1);
            if (t2 == trunc(t2)) {
                // t = k + 1/2: |x|^t · (±i)·(−1)^k
                const bool kodd = detail_math::is_odd_integer(t - F::half());
                const F b = eval_realfn(Comp::fn2(r_pow, -x, t, sb(y) != kodd));
                return {F::zero(), b};
            }
        }
        if (t == trunc(t) && abs(t) <= F(64)) return pow_int_exact(x, y, static_cast<long>(t.to_double()));
    }
    return pow_general(x, y, wr, wi);
}

} // namespace detail_cx

/// z^w = exp(w·log z) (principal branch); pow(z, 0) = 1 for every z (NaN too);
/// pow(0, w) = 0 for Re w > 0, +∞ (divide-by-zero) for real w < 0, else NaN (invalid).
[[nodiscard]] inline complex pow(const complex& z, const complex& w) {
    return detail_cx::pow_impl(z, w.real(), w.imag(), false);
}
[[nodiscard]] inline complex pow(const complex& z, float256 t) { return detail_cx::pow_impl(z, t, float256(), true); }
[[nodiscard]] inline complex pow(float256 x, const complex& w) {
    return detail_cx::pow_impl(complex(x, float256()), w.real(), w.imag(), false);
}

// ------------------------------------------------------------ stream I/O

/// "(re,im)" with each part as float256's operator<< writes it (shortest round-trip).
inline std::ostream& operator<<(std::ostream& os, const complex& z) {
    std::ostringstream s;
    s << '(' << z.real() << ',' << z.imag() << ')';
    return os << s.str();
}

/// Reads "re", "(re)" or "(re,im)" (each part a float256 literal: decimal,
/// hex, inf, nan); sets failbit on malformed input.
inline std::istream& operator>>(std::istream& is, complex& z) {
    auto token = [&](float256& out) -> bool {
        is >> std::ws;
        std::string t;
        while (true) {
            const int c = is.peek();
            if (c == std::char_traits<char>::eof() || c == ',' || c == ')' || std::isspace(c)) break;
            t.push_back(static_cast<char>(is.get()));
        }
        const auto v = float256::parse(t);
        if (!v) return false;
        out = *v;
        return true;
    };
    float256 re, im;
    is >> std::ws;
    if (is.peek() == '(') {
        is.get();
        if (!token(re)) { is.setstate(std::ios::failbit); return is; }
        is >> std::ws;
        const int c = is.get();
        if (c == ')') { z = complex(re, float256()); return is; }
        if (c != ',' || !token(im)) { is.setstate(std::ios::failbit); return is; }
        is >> std::ws;
        if (is.get() != ')') { is.setstate(std::ios::failbit); return is; }
        z = complex(re, im);
        return is;
    }
    if (!token(re)) { is.setstate(std::ios::failbit); return is; }
    z = complex(re, float256());
    return is;
}

} // namespace base256

namespace base256::literals {
/// Imaginary literals: 2.5_i256 is complex(0, 2.5) (parsed with correct rounding).
[[nodiscard]] inline complex operator""_i256(const char* s) { return {float256(), operator""_f256(s)}; }
[[nodiscard]] inline complex operator""_i256(const char* s, std::size_t n) { return {float256(), operator""_f256(s, n)}; }
} // namespace base256::literals

/// std::format: "(re,im)". The spec is float256's ([[fill]align][sign][#][width][.precision][type]);
/// sign, #, precision and type apply to each part, fill/align/width to the whole; '0' is not accepted.
template <>
struct std::formatter<base256::complex, char> {
    std::formatter<base256::float256, char> part;
    char fill = ' ';
    char align = 0;
    int width = 0;

    constexpr auto parse(std::format_parse_context& ctx) {
        auto it = part.parse(ctx);
        if (part.zero_pad) throw std::format_error("'0' is not supported for base256::complex");
        fill = part.fill;
        align = part.align;
        width = part.width;
        part.width = 0;
        part.align = 0;
        part.fill = ' ';
        return it;
    }

    template <class Ctx>
    auto format(const base256::complex& z, Ctx& ctx) const {
        std::string spec = "{:";
        if (part.sign != '-') spec += part.sign;
        if (part.alt) spec += '#';
        if (part.precision >= 0) spec += "." + std::to_string(part.precision);
        if (part.type) spec += part.type;
        spec += "}";
        const base256::float256 re = z.real(), im = z.imag();
        std::string body = "(" + std::vformat(spec, std::make_format_args(re)) + "," +
                           std::vformat(spec, std::make_format_args(im)) + ")";
        const int len = static_cast<int>(body.size());
        if (width > len) {
            const int padn = width - len;
            const char a = align ? align : '>';
            const int left = a == '<' ? 0 : (a == '^' ? padn / 2 : padn);
            body = std::string(static_cast<std::size_t>(left), fill) + body +
                   std::string(static_cast<std::size_t>(padn - left), fill);
        }
        return std::copy(body.begin(), body.end(), ctx.out());
    }
};
