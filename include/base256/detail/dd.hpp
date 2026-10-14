#pragma once

/// Double-float256 ("dd") arithmetic used internally by the elementary
/// functions: a value is an unevaluated sum hi + lo with |lo| <= ulp(hi)/2,
/// about 474 bits of precision. The error-free transformations rely on the
/// correctly rounded float256 +, × and fma (round-to-nearest-even), so callers
/// must run under fenv::scoped_default_env. Results are rounded to float256
/// once, by the final hi + lo addition.

#include <algorithm>
#include <array>
#include <cmath>
#include "../float256.hpp"
#include "bigint.hpp"
#include "two_over_pi.hpp"

namespace base256::detail_dd {

struct dd {
    float256 hi;
    float256 lo;
    dd() = default;
    dd(float256 h) : hi(h), lo(float256::zero()) {}  // NOLINT(google-explicit-constructor)
    dd(float256 h, float256 l) : hi(h), lo(l) {}
    [[nodiscard]] float256 round() const { return hi + lo; }
    /// Round once in the caller's mode through `env` (inexact result).
    [[nodiscard]] float256 out(class FnEnv& env, int k = 0, int tail = 2) const;
};

/// Environment guard for the public elementary functions. The body runs in
/// round-to-nearest-even with cleared flags (internal operations never touch
/// the caller's environment). The result is then produced by one of the
/// member functions, which restore the caller's environment first and round
/// or signal in the caller's mode, so only the final operation raises flags.
class FnEnv {
public:
    FnEnv() noexcept : saved_(base256::detail::tls_env) { base256::detail::tls_env = base256::detail::Env{}; }
    ~FnEnv() { if (!done_) base256::detail::tls_env = saved_; }
    FnEnv(const FnEnv&) = delete;
    FnEnv& operator=(const FnEnv&) = delete;

    [[nodiscard]] round_mode mode() const noexcept { return saved_.mode; }
    void restore() noexcept { base256::detail::tls_env = saved_; done_ = true; }

    /// An exact result: no flags.
    float256 exact(float256 x) noexcept { restore(); return x; }
    /// (hi + lo)·2^k rounded once in the caller's mode. `tail`: 0 exact,
    /// ±1 the exact value is slightly above/below, 2 inexact (direction from lo).
    float256 round(float256 hi, float256 lo, int k = 0, int tail = 2) noexcept {
        restore();
        return float256::round_dd(hi, lo, k, saved_.mode, tail);
    }
    /// NaN operand: invalid for a signaling NaN, result quieted.
    float256 nan(float256 x) noexcept {
        restore();
        if (x.issignaling()) base256::detail::raise_flag(fe_invalid);
        return x.quieted();
    }
    float256 nan2(float256 x, float256 y) noexcept {
        restore();
        if (x.issignaling() || y.issignaling()) base256::detail::raise_flag(fe_invalid);
        return x.isnan() ? x.quieted() : y.quieted();
    }
    /// Domain error.
    float256 invalid() noexcept { restore(); base256::detail::raise_flag(fe_invalid); return float256::qnan(); }
    /// Exact infinite result from finite operands (pole).
    float256 pole(bool negative) noexcept { restore(); base256::detail::raise_flag(fe_divbyzero); return float256::inf(negative); }
    /// Result far beyond max() / far below the subnormal range, rounded per mode.
    float256 overflow(bool negative) noexcept { return round(negative ? -float256::one() : float256::one(), float256::zero(), 300000); }
    float256 underflow(bool negative) noexcept { return round(negative ? -float256::one() : float256::one(), float256::zero(), -300000); }

private:
    base256::detail::Env saved_;
    bool done_ = false;
};

inline float256 dd::out(FnEnv& env, int k, int tail) const { return env.round(hi, lo, k, tail); }

[[nodiscard]] inline dd two_sum(float256 a, float256 b) {
    {
        float256 s, e;
        if (float256::eft_two_sum(a, b, s, e)) return {s, e};
    }
    const float256 s = a + b;
    const float256 bb = s - a;
    const float256 e = (a - (s - bb)) + (b - bb);
    return {s, e};
}
[[nodiscard]] inline dd fast_two_sum(float256 a, float256 b) {  // |a| >= |b|
    {
        float256 s, e;
        if (float256::eft_two_sum(a, b, s, e)) return {s, e};
    }
    const float256 s = a + b;
    return {s, b - (s - a)};
}
[[nodiscard]] inline dd two_prod(float256 a, float256 b) {
    {
        float256 p, e;
        if (float256::eft_two_prod(a, b, p, e)) return {p, e};
    }
    const float256 p = a * b;
    return {p, fma(a, b, -p)};
}

[[nodiscard]] inline dd operator-(const dd& a) { return {-a.hi, -a.lo}; }

[[nodiscard]] inline dd operator+(const dd& a, const dd& b) {
    {
        dd r;
        if (float256::dd_add_parts(a.hi, a.lo, b.hi, b.lo, r.hi, r.lo)) return r;
    }
    if (!a.hi.isfinite() || !b.hi.isfinite()) return dd(a.hi + b.hi);
    dd s = two_sum(a.hi, b.hi);
    const dd t = two_sum(a.lo, b.lo);
    s.lo = s.lo + t.hi;
    s = fast_two_sum(s.hi, s.lo);
    s.lo = s.lo + t.lo;
    return fast_two_sum(s.hi, s.lo);
}
[[nodiscard]] inline dd operator-(const dd& a, const dd& b) { return a + (-b); }

[[nodiscard]] inline dd operator*(const dd& a, const dd& b) {
    {
        dd r;
        if (float256::dd_mul_parts(a.hi, a.lo, b.hi, b.lo, r.hi, r.lo)) return r;
    }
    dd p = two_prod(a.hi, b.hi);
    if (!p.hi.isfinite()) return dd(p.hi);
    p.lo = p.lo + (a.hi * b.lo + a.lo * b.hi);
    return fast_two_sum(p.hi, p.lo);
}
[[nodiscard]] inline dd operator*(const dd& a, float256 b) {
    {
        dd r;
        if (float256::dd_mul_parts(a.hi, a.lo, b, float256{}, r.hi, r.lo)) return r;
    }
    dd p = two_prod(a.hi, b);
    if (!p.hi.isfinite()) return dd(p.hi);
    p.lo = p.lo + a.lo * b;
    return fast_two_sum(p.hi, p.lo);
}

[[nodiscard]] inline dd operator/(const dd& a, const dd& b) {
    const float256 q1 = a.hi / b.hi;
    if (!q1.isfinite() || q1.iszero()) return dd(q1);
    dd r = a - b * q1;
    const float256 q2 = r.hi / b.hi;
    r = r - b * q2;
    const float256 q3 = r.hi / b.hi;
    return fast_two_sum(q1, q2) + dd(q3);
}

/// dd / float256 (used for division by small integers): one fewer correction
/// step than dd / dd, error about 2^-470 relative.
[[nodiscard]] inline dd operator/(const dd& a, float256 b) {
    const float256 q1 = a.hi / b;
    if (!q1.isfinite() || q1.iszero()) return dd(q1);
    const dd p = two_prod(q1, b);
    const float256 r = ((a.hi - p.hi) - p.lo) + a.lo;
    return fast_two_sum(q1, r / b);
}

[[nodiscard]] inline dd ldexp(const dd& a, int k) { return {base256::ldexp(a.hi, k), base256::ldexp(a.lo, k)}; }

[[nodiscard]] inline dd sqrt(const dd& a) {
    if (a.hi.iszero()) return a;
    const float256 s = base256::sqrt(a.hi);
    const dd r = a - two_prod(s, s);
    return fast_two_sum(s, r.hi / (s + s));
}

[[nodiscard]] inline bool lt(const dd& a, float256 b) { return a.hi < b || (a.hi == b && a.lo.signbit() && !a.lo.iszero()); }

// Constants split into float256 parts (generated with MPFR).
[[nodiscard]] inline float256 hexc(const char* s) { return float256::parse(s).value(); }

struct Consts {
    float256 ln2[3];
    float256 pio2[5];
    float256 ln10[3];
    float256 pi[3];
    float256 sqrtpi[3];
    float256 half_ln_2pi[3];
    Consts() {
        const char* l2[3] = {"0xb.17217f7d1cf79abc9e3b39803f2f6af40f343267298b62d8a0d175b8ba8p-4",
                             "0x2.fa2be7b876206debac98559552fb4afa1b10ed2eae35c138214427573b2p-240",
                             "0x9.1169b8253e96ca16224ae8c51acbda11317c387eb9ea9bc3b136603b258p-480"};
        const char* p2[5] = {"0x1.921fb54442d18469898cc51701b839a252049c1114cf98e804177d4c762p+0",
                             "0x7.3644a29410f31c6809bbdf2a33679a748636605614dbe4be286e9fc26acp-240",
                             "0x1.adaa3848bc90b6aecc4bcfd8de89885d34c6fdad617feb96de80d6fdbdcp-476",
                             "0x7.0d7f6b5133f4b5d3e4822f8963fcc9250cca3d9c8b67b8400f97142c78p-716",
                             "-0x1.f4ce4b6f93c75458cb2dd380ae05b66140f93545b846b8a4d3c73a1953cp-952"};
        const char* l10[3] = {"0x2.4d763776aaa2b05ba95b58ae0b4c28a38a3fb3e76977e43a0f187a0807cp+0",
                              "0xb.5ca58bc0b5ec6a0417331c32f00b17c35a0b1889061042f8b6bee3de21p-244",
                              "0xb.945b59e0b3e28a2a324479d96a9b0ec360c7efbd9b3ac12acf1be94587p-488"};
        const char* pp[3] = {"0x3.243f6a8885a308d313198a2e03707344a4093822299f31d0082efa98ec4p+0",
                             "0xe.6c89452821e638d01377be5466cf34e90c6cc0ac29b7c97c50dd3f84d58p-240",
                             "0x3.5b54709179216d5d98979fb1bd1310ba698dfb5ac2ffd72dbd01adfb7b8p-476"};
        const char* sp[3] = {"0x1.c5bf891b4ef6aa79c3b0520d5db9383fe3921546f63b252dca100bd3ea1p+0",
                             "0x4.746ed76ffd6f941f1dbacd8e76141e0c77476c521446ddc35851524490cp-240",
                             "0x3.4a42174dcbafdf109d1b3acfb553ade465bea6c5032aa8a111a626e2c7p-480"};
        const char* hl[3] = {"0xe.b3f8e4325f5a53494bc900144192023cfb08f8d13458b4ddec6a3133da8p-4",
                             "0x2.155d212f9d7fe00e86bf93eabf905c5569bbb05cab571b4cda5bcf66ecep-240",
                             "0x1.772badf2926770eca54487a7307d867635c116969254d1304a59fb7e271p-480"};
        for (int i = 0; i < 3; ++i) {
            ln2[i] = hexc(l2[i]); ln10[i] = hexc(l10[i]); pi[i] = hexc(pp[i]);
            sqrtpi[i] = hexc(sp[i]); half_ln_2pi[i] = hexc(hl[i]);
        }
        for (int i = 0; i < 5; ++i) pio2[i] = hexc(p2[i]);
    }
};
[[nodiscard]] inline const Consts& consts() {
    static const Consts c = [] { const fenv::scoped_default_env g; return Consts{}; }();
    return c;
}
[[nodiscard]] inline dd ln2() { const auto& c = consts(); return {c.ln2[0], c.ln2[1]}; }
[[nodiscard]] inline dd ln10() { const auto& c = consts(); return {c.ln10[0], c.ln10[1]}; }
[[nodiscard]] inline dd pi() { const auto& c = consts(); return {c.pi[0], c.pi[1]}; }
[[nodiscard]] inline dd pio2() { const auto& c = consts(); return {c.pio2[0], c.pio2[1]}; }

// 1/n! for n = 0..kMaxInvFact as dd (computed once, ~2^-465 relative).
inline constexpr int kMaxInvFact = 120;
[[nodiscard]] inline const std::array<dd, kMaxInvFact + 1>& inv_factorials() {
    static const std::array<dd, kMaxInvFact + 1> t = [] {
        const fenv::scoped_default_env g;
        std::array<dd, kMaxInvFact + 1> a{};
        a[0] = dd(float256::one());
        for (int n = 1; n <= kMaxInvFact; ++n)
            a[static_cast<std::size_t>(n)] = a[static_cast<std::size_t>(n - 1)] / float256::from_int(n);
        return a;
    }();
    return t;
}
// Smallest N with |x|^(N-1) / N! < 2^-bits for |x| < 2^(e+1): the Taylor term
// of index N relative to the term of index 1.
[[nodiscard]] inline int series_terms(int e, int bits) {
    double lf = 0;  // log2(N!)
    const double m = -1.0 - e;  // |x| < 2^-m
    for (int n = 1; n < kMaxInvFact; ++n) {
        lf += std::log2(static_cast<double>(n));
        if (m * (n - 1) + lf >= bits) return n;
    }
    return kMaxInvFact;
}

// k * ln2 to ~710 bits for an integer k with |k| < 2^30.
[[nodiscard]] inline dd k_ln2(float256 k) {
    const auto& c = consts();
    return two_prod(k, c.ln2[0]) + (two_prod(k, c.ln2[1]) + dd(k * c.ln2[2]));
}

/// expm1 of a reduced dd argument |r| <= 0.35: Taylor on r/2^12, then
/// 12 squarings in the form em <- em (em + 2).
[[nodiscard]] inline dd expm1_reduced(const dd& r) {
    constexpr int kSq = 12;
    const dd x = ldexp(r, -kSq);
    if (x.hi.iszero()) return x;
    // Horner on 1/n! (precomputed): expm1(x) = x (1/1! + x (1/2! + x (...))).
    const int n = series_terms(base256::ilogb(x.hi), 500);
    const auto& f = inv_factorials();
    dd p = f[static_cast<std::size_t>(n)];
    for (int k = n - 1; k >= 1; --k) p = p * x + f[static_cast<std::size_t>(k)];
    dd sum = p * x;
    for (int i = 0; i < kSq; ++i) sum = sum * (sum + dd(float256::two()));
    return sum;
}

/// exp(x) = m * 2^k with m a dd near [0.7, 1.42]; x finite, |x| < 2^20.
struct Scaled { dd m; int k; };
[[nodiscard]] inline Scaled exp_scaled(const dd& x) {
    if (base256::ilogb(x.hi) < -200) {
        // 1 + x + x^2/2: x/2^12 below would underflow for very small x.
        return {dd(float256::one()) + x + dd(base256::ldexp(x.hi * x.hi, -1)), 0};
    }
    const float256 kf = base256::round(x.hi / consts().ln2[0]);
    const dd r = x - k_ln2(kf);
    const dd m = expm1_reduced(r) + dd(float256::one());
    return {m, static_cast<int>(kf.to_double())};
}
[[nodiscard]] inline float256 round_scaled(const dd& m, int k) {
    // Split the scaling so neither part overflows prematurely.
    const float256 hi = base256::ldexp(m.hi, k);
    if (!hi.isfinite()) return hi;
    // lo below the subnormal range is irrelevant except for the last bit.
    return hi + base256::ldexp(m.lo, k);
}
[[nodiscard]] inline dd exp(const dd& x) {
    const Scaled s = exp_scaled(x);
    return ldexp(s.m, s.k);
}
[[nodiscard]] inline dd expm1(const dd& x) {  // any finite x with |x| < 2^20
    if (base256::abs(x.hi) <= float256(0.34)) return expm1_reduced(x);
    return exp(x) - dd(float256::one());
}

/// Natural log of a positive finite dd, via one Newton step on exp from a
/// float256 starting value: y = y0 + log1p(m e^{-y0} - 1).
// float256 log of m in [0.7, 1.42] (atanh series), a few ulps.
[[nodiscard]] inline float256 log_seed(float256 m) {
    const float256 f = m - float256::one();
    if (f.iszero()) return f;
    const float256 z = f / (m + float256::one());
    const float256 z2 = z * z;
    float256 term = z, sum = z;
    for (int n = 1; n <= 200; ++n) {
        term = term * z2;
        const float256 next = sum + term / float256::from_int(2 * n + 1);
        if (next == sum) break;
        sum = next;
    }
    return sum + sum;
}
[[nodiscard]] inline dd log(const dd& a) {
    int e = 0;
    (void)base256::frexp(a.hi, &e);
    dd m = ldexp(a, -e);  // [0.5, 1)
    if (m.hi < float256(0.70710678118654752440)) { m = m + m; --e; }
    const float256 y0 = log_seed(m.round());
    const dd t = m * exp(dd(-y0)) - dd(float256::one());
    // log1p(t) = t - t^2/2 + t^3/3, |t| ~ 2^-230
    const dd t2 = t * t;
    const dd l1p = t - ldexp(t2, -1) + t2 * t / float256::from_int(3);
    dd y = dd(y0) + l1p;
    if (e != 0) y = k_ln2(float256::from_int(e)) + y;
    return y;
}

/// Cody–Waite reduction by π/2 split into five float256 parts (~1190 bits):
/// x = k·π/2 + r, |r| <= π/4, quadrant = k mod 4. Valid for |x| < 2^100
/// (k is an exact float256 integer and the reduced angle keeps > 474 bits).
inline constexpr int kTrigExactExp = 100;
[[nodiscard]] inline dd reduce_pio2(float256 x, int& quadrant) {
    const auto& c = consts();
    const float256 kf = base256::round(x / c.pio2[0]);
    dd r = dd(x) - two_prod(kf, c.pio2[0]);
    r = r - two_prod(kf, c.pio2[1]);
    r = r - two_prod(kf, c.pio2[2]);
    r = r - two_prod(kf, c.pio2[3]);
    r = r - dd(kf * c.pio2[4]);
    const float256 q = float256::remainder_impl(kf, float256::from_int(4), false, nullptr);
    int qi = static_cast<int>(q.to_double());
    if (qi < 0) qi += 4;
    quadrant = qi;
    return r;
}

/// Payne–Hanek reduction for any finite |x| >= 2^kTrigExactExp (x > 0):
/// x = M·2^E with M a 237-bit integer. Bits b_i of 2/π with i <= E-2 add
/// multiples of 4 to x·2/π and are skipped; bits E-1 .. E+kPhFrac are
/// multiplied by M exactly, which leaves kPhFrac fraction bits (error below
/// 2^(237-kPhFrac) ulps of the fraction). Returns r = (x mod π/2) with
/// |r| <= π/4 as dd and the quadrant.
inline constexpr int kPhFrac = 1216;  // multiple of 32
static_assert(kPhFrac % 32 == 0);
[[nodiscard]] inline dd payne_hanek(float256 x, int& quadrant) {
    using base256::detail::BigInt;
    const auto w = x.raw_bits();
    const int be = static_cast<int>((w[3] >> 44) & 0x7FFFF);
    const int E = be - 262143 - 236;
    // M as little-endian 32-bit limbs.
    const std::uint64_t top = (w[3] & 0xFFFFFFFFFFFULL) | (1ULL << 44);
    const std::uint64_t mw[4] = {w[0], w[1], w[2], top};
    std::vector<std::uint32_t> ml;
    for (std::uint64_t v : mw) {
        ml.push_back(static_cast<std::uint32_t>(v));
        ml.push_back(static_cast<std::uint32_t>(v >> 32));
    }
    BigInt M = BigInt::from_limbs(std::move(ml));
    // Window of 2/π bits i0..i1 (1-based fraction positions) as an integer.
    const int i0 = E - 1 > 1 ? E - 1 : 1;
    const int i1 = E + kPhFrac;
    const int nbits = i1 - i0 + 1;
    std::vector<std::uint32_t> wl(static_cast<std::size_t>((nbits + 31) / 32), 0);
    for (int i = i0; i <= i1; ++i) {
        const int k = (i - 1) >> 6, b = 63 - ((i - 1) & 63);
        if ((kTwoOverPi[k] >> b) & 1ULL) {
            const int pos = i1 - i;  // weight 2^pos in W
            wl[static_cast<std::size_t>(pos >> 5)] |= 1U << (pos & 31);
        }
    }
    const BigInt W = BigInt::from_limbs(std::move(wl));
    M.mul(W);  // x·2/π mod 4 = M / 2^kPhFrac (mod 4)
    const int q = (M.bit(kPhFrac) ? 1 : 0) | (M.bit(kPhFrac + 1) ? 2 : 0);
    // Fraction f in [0, 1) with kPhFrac bits.
    BigInt f = M;
    f.keep_low_limbs(kPhFrac / 32);  // kPhFrac is a multiple of 32
    bool neg = false;
    int quad = q;
    if (f.bit(kPhFrac - 1)) {  // f >= 1/2: use f - 1 and the next quadrant
        std::vector<std::uint32_t> ol(static_cast<std::size_t>(kPhFrac / 32) + 1, 0);
        ol.back() = 1;
        BigInt one = BigInt::from_limbs(std::move(ol));
        one.sub(f);
        f = std::move(one);
        neg = true;
        quad = (q + 1) & 3;
    }
    // f to dd: sum its top 32-bit limbs (exact float256 values) in dd.
    dd r{float256::zero()};
    const int n = f.limb_count();
    for (int j = n - 1; j >= 0 && j >= n - 20; --j) {
        r = r + dd(base256::ldexp(float256::from_int(static_cast<std::int64_t>(f.limb(j))),
                                  32 * j - kPhFrac));
    }
    const auto& c = consts();
    r = r * dd(c.pio2[0], c.pio2[1]);
    quadrant = quad;
    return neg ? -r : r;
}

/// sin and cos of a dd |r| <= ~0.8: Taylor at r/2^8, then 8 double-angle steps.
inline void sincos_reduced(const dd& r, dd& s, dd& c) {
    const dd one{float256::one()};
    if (r.hi.iszero()) { s = r; c = one; return; }
    if (base256::abs(r.hi) < base256::ldexp(float256::one(), -120)) {
        const dd r2 = r * r;
        s = r - r2 * r / float256::from_int(6);
        c = one - ldexp(r2, -1);
        return;
    }
    constexpr int kHalf = 8;
    const dd x = ldexp(r, -kHalf);
    const dd x2 = x * x;
    // Horner in x2 on ±1/(2k+1)! and ±1/(2k)! (precomputed).
    const auto& f = inv_factorials();
    int n = series_terms(base256::ilogb(x.hi), 495);  // terms of the full series x^n/n!
    n = std::min(n + 2, kMaxInvFact - 1);
    const int ks = (n - 1) / 2, kc = n / 2;  // highest k used in sin / cos
    auto coef = [&](int idx, int k) { const dd c = f[static_cast<std::size_t>(idx)]; return (k & 1) ? -c : c; };
    dd ps = coef(2 * ks + 1, ks), pc = coef(2 * kc, kc);
    for (int k = ks - 1; k >= 0; --k) ps = ps * x2 + coef(2 * k + 1, k);
    for (int k = kc - 1; k >= 0; --k) pc = pc * x2 + coef(2 * k, k);
    dd ss = ps * x, sc = pc;
    for (int i = 0; i < kHalf; ++i) {
        const dd s2 = ldexp(ss * sc, 1);
        sc = one - ldexp(ss * ss, 1);
        ss = s2;
    }
    s = ss;
    c = sc;
}

} // namespace base256::detail_dd
