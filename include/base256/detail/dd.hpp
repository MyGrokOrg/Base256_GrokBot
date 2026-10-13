#pragma once

/// Double-float256 ("dd") arithmetic used internally by the elementary
/// functions: a value is an unevaluated sum hi + lo with |lo| <= ulp(hi)/2,
/// about 474 bits of precision. The error-free transformations rely on the
/// correctly rounded float256 +, × and fma (round-to-nearest-even), so callers
/// must run under fenv::scoped_default_env. Results are rounded to float256
/// once, by the final hi + lo addition.

#include "../float256.hpp"

namespace base256::detail_dd {

struct dd {
    float256 hi;
    float256 lo;
    dd() = default;
    dd(float256 h) : hi(h), lo(float256::zero()) {}  // NOLINT(google-explicit-constructor)
    dd(float256 h, float256 l) : hi(h), lo(l) {}
    [[nodiscard]] float256 round() const { return hi + lo; }
};

[[nodiscard]] inline dd two_sum(float256 a, float256 b) {
    const float256 s = a + b;
    const float256 bb = s - a;
    const float256 e = (a - (s - bb)) + (b - bb);
    return {s, e};
}
[[nodiscard]] inline dd fast_two_sum(float256 a, float256 b) {  // |a| >= |b|
    const float256 s = a + b;
    return {s, b - (s - a)};
}
[[nodiscard]] inline dd two_prod(float256 a, float256 b) {
    const float256 p = a * b;
    return {p, fma(a, b, -p)};
}

[[nodiscard]] inline dd operator-(const dd& a) { return {-a.hi, -a.lo}; }

[[nodiscard]] inline dd operator+(const dd& a, const dd& b) {
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
    dd p = two_prod(a.hi, b.hi);
    if (!p.hi.isfinite()) return dd(p.hi);
    p.lo = p.lo + (a.hi * b.lo + a.lo * b.hi);
    return fast_two_sum(p.hi, p.lo);
}
[[nodiscard]] inline dd operator*(const dd& a, float256 b) {
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
    dd term = x, sum = x;
    for (int n = 2; n <= 60; ++n) {
        term = term * x / dd(float256::from_int(n));
        sum = sum + term;
        if (term.hi.iszero() || base256::abs(term.hi) < base256::ldexp(base256::abs(sum.hi), -500)) break;
    }
    for (int i = 0; i < kSq; ++i) sum = sum * (sum + dd(float256::two()));
    return sum;
}

/// exp(x) = m * 2^k with m a dd near [0.7, 1.42]; x finite, |x| < 2^20.
struct Scaled { dd m; int k; };
[[nodiscard]] inline Scaled exp_scaled(const dd& x) {
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
    const dd l1p = t - ldexp(t2, -1) + t2 * t / dd(float256::from_int(3));
    dd y = dd(y0) + l1p;
    if (e != 0) y = k_ln2(float256::from_int(e)) + y;
    return y;
}

} // namespace base256::detail_dd
