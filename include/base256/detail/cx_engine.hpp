#pragma once
// Internal engine of base256::complex (complex.hpp): exact arithmetic on
// binary floats with BigInt mantissas, correct rounding of exact values,
// rounding decisions for approximations (with an exact boundary oracle for the
// algebraic operations), and the multiprecision real kernels (sin/cos with
// Payne–Hanek reduction, exp/expm1, log/log1p, atan/atan2, sinh/cosh) that the
// transcendental complex functions evaluate in a Ziv loop.
// Runtime only (not constexpr).

#include <algorithm>
#include <array>
#include <climits>
#include <cstdint>
#include <vector>

#include "../extra_math.hpp"
#include "mp.hpp"

namespace base256::detail_cx {

using detail::BigInt;
using detail_mp::MP;
namespace M = detail_mp;

// ---------------------------------------------------------------- conversions

/// Exact value of a finite float256 (zero for ±0, ±∞, NaN).
[[nodiscard]] inline MP mpf(float256 x) {
    MP r;
    if (!x.isfinite() || x.iszero()) return r;
    const auto w = x.raw_bits();
    const int be = static_cast<int>((w[3] >> 44) & detail::kExpAllOnes);
    std::uint64_t top = w[3] & detail::kFracHiMask;
    if (be != 0) top |= (1ULL << 44);
    std::vector<std::uint32_t> limbs;
    limbs.reserve(8);
    for (int i = 0; i < 3; ++i) {
        limbs.push_back(static_cast<std::uint32_t>(w[static_cast<std::size_t>(i)]));
        limbs.push_back(static_cast<std::uint32_t>(w[static_cast<std::size_t>(i)] >> 32));
    }
    limbs.push_back(static_cast<std::uint32_t>(top));
    limbs.push_back(static_cast<std::uint32_t>(top >> 32));
    r.m = BigInt::from_limbs(std::move(limbs));
    r.e = (be == 0 ? detail::kEmin : be - detail::kBias) - detail::kFracBits;
    r.neg = x.signbit();
    return r;
}

/// m·2^e as a float256 for 1 <= m < 2^237 and a result in the normal range.
[[nodiscard]] inline float256 f256_exact(const BigInt& m, long e) {
    const int bl = m.bit_length();
    BigInt s = m;
    s.shl(237 - bl);  // leading bit at 236
    const long be = e - (237 - bl) + detail::kFracBits + detail::kBias;
    std::uint64_t w[4];
    for (int i = 0; i < 4; ++i)
        w[i] = static_cast<std::uint64_t>(s.limb(2 * i)) | (static_cast<std::uint64_t>(s.limb(2 * i + 1)) << 32);
    w[3] = (w[3] & detail::kFracHiMask) | (static_cast<std::uint64_t>(be) << 44);
    return float256::from_bits(w[0], w[1], w[2], w[3]);
}

[[nodiscard]] inline MP one() { return M::from_u64(1); }
[[nodiscard]] inline MP pow2(long e, bool neg = false) { MP r = M::from_u64(1, neg); r.e = e; return r; }
[[nodiscard]] inline int sgn(const MP& a) { return a.zero() ? 0 : (a.neg ? -1 : 1); }
[[nodiscard]] inline MP absv(MP a) { a.neg = false; return a; }
[[nodiscard]] inline MP trunc(MP a, int P) { M::norm(a, P); return a; }

// ---------------------------------------------------------- exact arithmetic

[[nodiscard]] inline MP xmul(const MP& a, const MP& b) {
    MP r;
    if (a.zero() || b.zero()) return r;
    r.m = a.m;
    r.m.mul(b.m);
    r.e = a.e + b.e;
    r.neg = a.neg != b.neg;
    return r;
}
[[nodiscard]] inline MP xadd(const MP& x, const MP& y) {
    if (x.zero()) return y;
    if (y.zero()) return x;
    MP a = x, b = y;
    const long emin = a.e < b.e ? a.e : b.e;
    a.m.shl(static_cast<int>(a.e - emin));
    a.e = emin;
    b.m.shl(static_cast<int>(b.e - emin));
    b.e = emin;
    if (a.neg == b.neg) a.m.add(b.m);
    else if (b.m < a.m) a.m.sub(b.m);
    else { BigInt t = b.m; t.sub(a.m); a.m = t; a.neg = b.neg; }
    if (a.m.is_zero()) { a.neg = false; a.e = 0; }
    return a;
}
[[nodiscard]] inline MP xsub(const MP& a, const MP& b) { return xadd(a, M::neg(b)); }

// ------------------------------------------------------------------ rounding

/// A value prepared for float256::round_dd: (hi + lo)·2^k, `tail` the sign of
/// the exact value minus that (0: exact).
struct Plan {
    float256 hi, lo;
    int k = 0;
    int tail = 0;
    [[nodiscard]] float256 round(round_mode m) const { return float256::round_dd(hi, lo, k, m, tail); }
};

/// Plan for a nonzero MP value v: the top 474 bits as hi + lo, the bits below
/// as the tail; `tail_if_exact` is the tail when no bits are discarded.
[[nodiscard]] inline Plan plan_from(const MP& v, int tail_if_exact) {
    BigInt m = v.m;
    long e = v.e;
    const int bl = m.bit_length();
    bool sticky = false;
    if (bl > 474) { sticky = m.shr_sticky(bl - 474); e += bl - 474; }
    else { m.shl(474 - bl); e -= 474 - bl; }
    BigInt h = m;
    h.shr(237);
    BigInt hs = h;
    hs.shl(237);
    BigInt l = m;
    l.sub(hs);
    Plan p;
    p.hi = f256_exact(h, -237);
    p.lo = l.is_zero() ? float256::zero() : f256_exact(l, -474);
    if (v.neg) { p.hi = -p.hi; p.lo = -p.lo; }
    long k = e + 474;
    k = std::clamp(k, -(1L << 28), 1L << 28);  // far outside the range either way
    p.k = static_cast<int>(k);
    p.tail = sticky ? (v.neg ? -1 : 1) : tail_if_exact;
    return p;
}

/// Plan for an exact float256 plus a nonzero amount far below its last bit.
[[nodiscard]] inline Plan plan_tiny(float256 base, int tail) {
    Plan p;
    p.hi = base;
    p.lo = float256::zero();
    p.tail = tail;
    return p;
}

enum class Dec { decided, boundary, unknown };
inline constexpr long kExact = LONG_MIN;

/// low t bits of m
[[nodiscard]] inline BigInt low_bits(const BigInt& m, long t) {
    if (t >= m.bit_length()) return m;
    BigInt h = m;
    h.shr(static_cast<int>(t));
    h.shl(static_cast<int>(t));
    BigInt l = m;
    l.sub(h);
    return l;
}

/// Where the value t, known as v with |t − v| <= 2^err_e, rounds: `decided`
/// when [v − err, v + err] holds no multiple of the half-unit grid 2^(E−237)
/// (E = exponent of v, unbounded: every rounding boundary of every mode,
/// subnormal ones included, lies on it); `boundary` with the one grid point B
/// in it when err < grid/4; `unknown` otherwise.
[[nodiscard]] inline Dec decide(const MP& v, long err_e, MP& B) {
    if (v.zero()) return Dec::unknown;
    const long g = v.top() - 238;
    if (err_e != kExact && err_e >= g - 2) return Dec::unknown;
    const long t = g - v.e;
    if (t <= 0) { B = v; return Dec::boundary; }
    const BigInt low = low_bits(v.m, t);
    bool dec;
    if (err_e == kExact || err_e < v.e) {
        dec = !low.is_zero();
    } else {
        BigInt E = M::from_u64(1).m;
        E.shl(static_cast<int>(err_e - v.e));
        BigInt full = M::from_u64(1).m;
        full.shl(static_cast<int>(t));
        BigInt high = full;
        high.sub(low);
        dec = E < low && E < high;
    }
    if (dec) return Dec::decided;
    BigInt q = v.m;
    q.shr(static_cast<int>(t));
    if (low.bit(static_cast<int>(t - 1))) q.add_u32(1);
    B.m = q;
    B.e = g;
    B.neg = v.neg;
    return Dec::boundary;
}

/// Correct rounding of an algebraic value approximated by v (|t − v| <= 2^err_e)
/// with `oracle(B)` = sign(t − B) exact. False: the error is too large.
template <class Oracle>
[[nodiscard]] bool plan_algebraic(const MP& v, long err_e, Oracle&& oracle, Plan& out) {
    MP B;
    switch (decide(v, err_e, B)) {
        case Dec::decided: out = plan_from(v, v.neg ? -1 : 1); return true;
        case Dec::boundary: {
            const int s = oracle(B);
            out = plan_from(B, s);
            return true;
        }
        default: return false;
    }
}

/// sign of a + b + ... for exact terms
[[nodiscard]] inline int sign_of(std::initializer_list<MP> terms) {
    MP s;
    for (const MP& t : terms) s = xadd(s, t);
    return sgn(s);
}

/// Sign of an exact zero sum p1 + p2 of two (possibly zero) exact products in
/// IEEE 754 arithmetic: same-signed zeros keep their sign; otherwise +0, or −0
/// when rounding downward.
[[nodiscard]] inline bool zero_sum_negative(bool z1, bool n1, bool z2, bool n2, round_mode m) {
    if (z1 && z2 && n1 == n2) return n1;
    return m == round_mode::downward;
}

// ------------------------------------------------------------- real kernels
// All take and return MP values with a relative error of a few units of 2^-P
// (the Ziv driver checks two precisions against each other).

[[nodiscard]] inline int isqrt_i(int P) { return static_cast<int>(std::sqrt(static_cast<double>(P))); }

/// e^x − 1, relative accuracy also for tiny x.
[[nodiscard]] inline MP mexpm1(const MP& x, int P) {
    if (x.zero()) return x;
    const int Q = P + 32;
    if (x.top() <= -1) {  // |x| < 1/2
        MP sum = x, t = x;
        for (std::uint32_t n = 2;; ++n) {
            t = M::div_u(M::mul(t, x, Q), n, Q);
            sum = M::add(sum, t, Q);
            if (t.zero() || t.top() < sum.top() - Q) break;
        }
        return trunc(sum, P);
    }
    return trunc(M::sub(M::exp(x, Q), one(), Q), P);
}

/// e^x for |x| < 2^30.
[[nodiscard]] inline MP mexp(const MP& x, int P) { return M::exp(x, P + 16); }

/// sinh, cosh, cosh − 1 = 2 sinh²(x/2)
[[nodiscard]] inline MP msinh(const MP& x, int P) {
    const int Q = P + 32;
    if (x.zero()) return x;
    if (x.top() <= -1) {
        const MP em = mexpm1(x, Q);
        return trunc(M::ldexp(M::add(em, M::div(em, M::add(em, one(), Q), Q), Q), -1), P);
    }
    const MP e = mexp(x, Q);
    return trunc(M::ldexp(M::sub(e, M::div(one(), e, Q), Q), -1), P);
}
[[nodiscard]] inline MP mcosh(const MP& x, int P) {
    const int Q = P + 32;
    const MP e = mexp(x, Q);
    return trunc(M::ldexp(M::add(e, M::div(one(), e, Q), Q), -1), P);
}
[[nodiscard]] inline MP mcosh_m1(const MP& x, int P) {
    const int Q = P + 32;
    const MP s = msinh(M::ldexp(x, -1), Q);
    return trunc(M::ldexp(M::mul(s, s, Q), 1), P);
}

/// log x, x > 0
[[nodiscard]] inline MP mlog(const MP& x, int P) { return M::log(trunc(x, P + 64), P + 16); }

/// log(1 + x), x > −1, relative accuracy also for tiny x: 2 atanh(x/(2 + x)).
[[nodiscard]] inline MP mlog1p(const MP& x, int P) {
    if (x.zero()) return x;
    const int Q = P + 32;
    if (x.top() <= -1) {
        const MP u = M::div(x, M::add(M::from_u64(2), x, Q), Q);
        const MP u2 = M::mul(u, u, Q);
        MP sum = u, t = u;
        for (std::uint32_t k = 1;; ++k) {
            t = M::mul(t, u2, Q);
            const MP s = M::div_u(t, 2 * k + 1, Q);
            sum = M::add(sum, s, Q);
            if (s.zero() || s.top() < sum.top() - Q) break;
        }
        return trunc(M::ldexp(sum, 1), P);
    }
    return mlog(M::add(one(), x, Q), P);
}

/// atan t for |t| <= 1: halvings t → t/(1 + sqrt(1 + t²)), then the series.
[[nodiscard]] inline MP atan_small(const MP& t0, int P) {
    if (t0.zero()) return t0;
    const int h = isqrt_i(P) / 3 + 2;
    const int Q = P + 2 * h + 32;
    MP t = t0;
    for (int i = 0; i < h; ++i)
        t = M::div(t, M::add(one(), M::sqrt(M::add(one(), M::mul(t, t, Q), Q), Q), Q), Q);
    const MP t2 = M::mul(t, t, Q);
    MP sum = t, p = t;
    for (std::uint32_t n = 1;; ++n) {
        p = M::neg(M::mul(p, t2, Q));
        const MP s = M::div_u(p, 2 * n + 1, Q);
        sum = M::add(sum, s, Q);
        if (s.zero() || s.top() < sum.top() - Q) break;
    }
    return trunc(M::ldexp(sum, h), P);
}

/// atan2(y, x) (not both zero)
[[nodiscard]] inline MP matan2(const MP& y, const MP& x, int P) {
    const int Q = P + 32;
    const MP pi = M::pi(Q);
    if (y.zero()) return x.neg ? trunc(pi, P) : MP{};
    if (x.zero()) { MP r = trunc(M::ldexp(pi, -1), P); r.neg = y.neg; return r; }
    const MP ax = absv(x), ay = absv(y);
    MP th;
    if (M::cmp_abs(ay, ax) <= 0) th = atan_small(M::div(ay, ax, Q), Q);
    else th = M::sub(M::ldexp(pi, -1), atan_small(M::div(ax, ay, Q), Q), Q);
    if (x.neg) th = M::sub(pi, th, Q);
    th = trunc(th, P);
    th.neg = y.neg;
    return th;
}

/// sin and cos of an exact value y (any size): y·(2/π) = n + f with |f| <= 1/2
/// from the 2/π table (Payne–Hanek, with P + 128 bits below the binary point),
/// r = f·π/2, then the quadrant.
inline void msincos(const MP& y, int P, MP& s, MP& c) {
    const int Q = P + 32;
    if (y.zero()) { s = MP{}; c = one(); return; }
    if (y.top() <= 0) { M::sincos(y, Q, s, c); s = trunc(s, P); c = trunc(c, P); return; }
    const long W = std::max(0L, y.top()) + Q + 128;
    const MP T = M::two_over_pi(static_cast<int>(W));
    BigInt prod = y.m;
    prod.mul(T.m);
    const long F = -(y.e + T.e);  // fraction bits (> 0: T has more bits than y's exponent)
    BigInt I = prod;
    I.shr(static_cast<int>(F));
    BigInt frac = low_bits(prod, F);
    unsigned q = I.limb(0) & 3u;
    bool fneg = false;
    if (frac.bit(static_cast<int>(F - 1))) {  // round to nearest: f = frac − 1
        q = (q + 1) & 3u;
        BigInt full = M::from_u64(1).m;
        full.shl(static_cast<int>(F));
        full.sub(frac);
        frac = full;
        fneg = true;
    }
    MP f;
    f.m = frac;
    f.e = -F;
    f.neg = fneg;
    if (f.m.is_zero()) f = MP{};
    M::norm(f, Q + 32);
    const MP r = M::mul(f, M::ldexp(M::pi(Q + 32), -1), Q + 32);
    MP sr, cr;
    if (r.zero()) { sr = MP{}; cr = one(); }
    else M::sincos(r, Q, sr, cr);
    switch (q) {
        case 0: s = sr; c = cr; break;
        case 1: s = cr; c = M::neg(sr); break;
        case 2: s = M::neg(sr); c = M::neg(cr); break;
        default: s = M::neg(cr); c = sr; break;
    }
    if (y.neg) s = M::neg(s);
    s = trunc(s, P);
    c = trunc(c, P);
}

/// sincos of an approximate MP argument: carries top(y) extra bits so the
/// absolute error of the reduction stays at 2^-P (the caller supplies y with
/// that accuracy).


/// ln 10
[[nodiscard]] inline MP ln10(int P) { return M::log(M::from_u64(10), P + 16); }

// ---------------------------------------------------------- the Ziv driver

/// One component of a complex result as produced by a kernel at precision P.
struct Comp {
    enum Kind : unsigned char { none, exact, realfn, approx };
    Kind k = none;
    float256 val;                        // exact
    float256 (*f1)(float256) = nullptr;  // realfn: f1(a) or f2(a, b), negated if `neg`
    float256 (*f2)(float256, float256) = nullptr;
    float256 a, b;
    bool neg = false;
    MP base;                             // approx: value = base + v (base exact)
    MP v;

    static Comp ex(float256 x) { Comp c; c.k = exact; c.val = x; return c; }
    static Comp fn(float256 (*f)(float256), float256 x, bool n = false) {
        Comp c; c.k = realfn; c.f1 = f; c.a = x; c.neg = n; return c;
    }
    static Comp fn2(float256 (*f)(float256, float256), float256 x, float256 y, bool n = false) {
        Comp c; c.k = realfn; c.f2 = f; c.a = x; c.b = y; c.neg = n; return c;
    }
    static Comp ap(MP x) { Comp c; c.k = approx; c.v = std::move(x); return c; }
    /// base + x with an exact base (a float256 or an exact MP leading term).
    static Comp ap(MP x, float256 bse) { return ap(std::move(x), mpf(bse)); }
    static Comp ap(MP x, MP bse) { Comp c; c.k = approx; c.v = std::move(x); c.base = std::move(bse); return c; }
    void negate() {
        switch (k) {
            case exact: val = -val; break;
            case realfn: neg = !neg; break;
            case approx: base = M::neg(base); v = M::neg(v); break;
            default: break;
        }
    }
};

/// Rounding mode with the directions of upward/downward exchanged (rounding
/// −t in mode m is −(rounding t in flip(m))).
[[nodiscard]] constexpr round_mode flip(round_mode m) noexcept {
    return m == round_mode::upward ? round_mode::downward
         : m == round_mode::downward ? round_mode::upward : m;
}

/// The caller's environment, saved; the work runs in round-to-nearest with no
/// flags recorded; results are rounded after restore().
class CxEnv {
public:
    CxEnv() noexcept : saved_(detail::tls_env) { detail::tls_env = detail::Env{}; }
    ~CxEnv() { if (!done_) detail::tls_env = saved_; }
    CxEnv(const CxEnv&) = delete;
    CxEnv& operator=(const CxEnv&) = delete;
    void restore() noexcept { if (!done_) { detail::tls_env = saved_; done_ = true; } }
    [[nodiscard]] round_mode mode() const noexcept { return saved_.mode; }
private:
    detail::Env saved_;
    bool done_ = false;
};

inline float256 eval_realfn(const Comp& c) {
    auto call = [&] { return c.f1 ? c.f1(c.a) : c.f2(c.a, c.b); };
    if (!c.neg) return call();
    const round_mode m = fenv::get_round();
    fenv::set_round(flip(m));
    const float256 r = call();
    fenv::set_round(m);
    return -r;
}

inline constexpr int kZivPrec[] = {300, 360, 720, 1440, 2880};
inline constexpr int kZivLevels = 5;

/// Decision for an approx component from the evaluations at the previous and
/// the current precision. True and `out` set when decided; at the last level
/// (`last`) a value on a boundary is taken as exact and an undecidable value
/// near zero as zero.
[[nodiscard]] inline bool ziv_step(const Comp& prev, const Comp& cur, int Pcur, bool last, Plan& out, bool& is_zero) {
    is_zero = false;
    const MP& v2 = cur.v;
    // An exact base of at most 474 bits is a Plan of its own; a tiny rest only adds its sign as the tail.
    const bool short_base = !cur.base.zero() && cur.base.m.bit_length() <= 474;
    if (v2.zero()) {
        if (!last) return false;
        if (cur.base.zero()) { is_zero = true; return true; }
        out = plan_from(cur.base, 0);
        return true;
    }
    const MP d = xsub(v2, prev.v);
    long err = v2.top() - Pcur + 24;
    if (!d.zero()) err = std::max(err, d.top());
    err += 1;
    MP v = v2;
    if (!cur.base.zero()) {
        if (short_base && v2.top() < cur.base.top() - 600) {
            if (err < v2.top() - 2) { out = plan_from(cur.base, sgn(v2)); return true; }
            if (!last) return false;
            out = plan_from(cur.base, 0);
            return true;
        }
        v = xadd(cur.base, v2);
        if (v.zero()) {
            if (!last) return false;
            is_zero = true;
            return true;
        }
    }
    MP B;
    switch (decide(v, err, B)) {
        case Dec::decided: out = plan_from(v, v.neg ? -1 : 1); return true;
        case Dec::boundary:
            if (!last) return false;
            out = plan_from(B, 0);
            return true;
        default:
            if (!last) return false;
            if (err >= v.top()) { is_zero = true; return true; }
            out = plan_from(v, 0);
            return true;
    }
}

/// Run `kernel(P, re, im)` at increasing precision until both approx
/// components are decided, then round in the caller's mode (flags raised there).
template <class Kernel>
[[nodiscard]] std::array<float256, 2> ziv(Kernel&& kernel) {
    CxEnv env;
    std::array<Comp, 2> prev, cur;
    kernel(kZivPrec[0], prev[0], prev[1]);
    std::array<bool, 2> done{prev[0].k != Comp::approx, prev[1].k != Comp::approx};
    std::array<Plan, 2> plan;
    std::array<bool, 2> zero{false, false};
    for (int lv = 1; lv < kZivLevels && !(done[0] && done[1]); ++lv) {
        kernel(kZivPrec[lv], cur[0], cur[1]);
        for (int i = 0; i < 2; ++i) {
            if (done[i]) continue;
            bool z = false;
            if (ziv_step(prev[i], cur[i], kZivPrec[lv], lv == kZivLevels - 1, plan[i], z)) {
                done[i] = true;
                zero[i] = z;
            }
        }
        prev = cur;
    }
    env.restore();
    const round_mode m = env.mode();
    std::array<float256, 2> r;
    for (int i = 0; i < 2; ++i) {
        const Comp& c = prev[i];
        if (c.k == Comp::exact) r[i] = c.val;
        else if (c.k == Comp::realfn) r[i] = eval_realfn(c);
        else if (zero[i]) r[i] = float256::zero(m == round_mode::downward);
        else r[i] = plan[i].round(m);
    }
    return r;
}

/// A value far outside the binary256 range (overflows or underflows in every
/// mode) with the given sign: 2^±300000.
[[nodiscard]] inline MP huge(int sign, bool big) {
    MP r;  // (2^400 − 1)·2^(±300000 − 400): many bits, so never on a rounding boundary
    r.m = M::from_u64(1).m;
    r.m.shl(400);
    r.m.sub(M::from_u64(1).m);
    r.e = (big ? 300000L : -300000L) - 400;
    r.neg = sign < 0;
    return r;
}
/// A nonzero value far below any rounding position of a base of order 1.
[[nodiscard]] inline MP tiny(int sign) { return huge(sign, false); }

} // namespace base256::detail_cx
