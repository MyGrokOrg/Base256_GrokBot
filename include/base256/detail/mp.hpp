#pragma once
// A small arbitrary-precision binary float (BigInt mantissa) for the few
// special-function paths whose cancellation exceeds what double-float256 can
// carry (lgamma next to its negative zeros, Hurwitz zeta next to zeros in a).
// Values are ±m·2^e, truncated to a working precision of P bits; every
// operation's error is a few units of 2^-P relative. Runtime only, not
// constexpr, and slow (tens of µs to ms): use only as a fallback.

#include <cmath>
#include <cstdint>
#include <mutex>
#include <vector>

#include "dd.hpp"
#include "two_over_pi_ext.hpp"

namespace base256::detail_mp {

using detail::BigInt;
using detail_dd::dd;

struct MP {
    bool neg = false;
    BigInt m;        // 0 means zero
    long e = 0;      // value = (−1)^neg · m · 2^e
    [[nodiscard]] bool zero() const { return m.is_zero(); }
    [[nodiscard]] long top() const { return e + m.bit_length(); }  // |v| < 2^top
};

inline void norm(MP& a, int P) {
    const int bl = a.m.bit_length();
    if (bl > P) { a.m.shr(bl - P); a.e += bl - P; }
    if (a.m.is_zero()) { a.neg = false; a.e = 0; }
}

[[nodiscard]] inline MP from_u64(std::uint64_t v, bool neg = false) {
    MP r; r.m = BigInt(v); r.neg = neg && v != 0; return r;
}

[[nodiscard]] inline MP from_f256(float256 x) {
    MP r;
    if (x.iszero() || !x.isfinite()) return r;
    int ex = 0;
    float256 f = frexp(abs(x), &ex);  // [1/2, 1)
    f = ldexp(f, 237);                 // integer < 2^237 (subnormals: fewer bits, still exact)
    const float256 two32 = ldexp(float256::one(), 32);
    std::vector<std::uint32_t> limbs;
    for (int i = 0; i < 8; ++i) {
        const float256 q = floor(f / two32);
        const float256 rem = f - q * two32;  // exact
        limbs.push_back(static_cast<std::uint32_t>(rem.to_double()));
        f = q;
    }
    r.m = BigInt::from_limbs(limbs);
    r.e = static_cast<long>(ex) - 237;
    r.neg = x.signbit();
    return r;
}

[[nodiscard]] inline MP from_dd(const dd& v) {
    MP a = from_f256(v.hi), b = from_f256(v.lo);
    // exact sum (the parts may be far apart; P large enough to keep both)
    if (b.zero()) return a;
    if (a.zero()) return b;
    const long emin = a.e < b.e ? a.e : b.e;
    a.m.shl(static_cast<int>(a.e - emin)); a.e = emin;
    b.m.shl(static_cast<int>(b.e - emin)); b.e = emin;
    if (a.neg == b.neg) a.m.add(b.m);
    else if (b.m < a.m) a.m.sub(b.m);
    else { BigInt t = b.m; t.sub(a.m); a.m = t; a.neg = b.neg; }
    if (a.m.is_zero()) a.neg = false;
    return a;
}

// The top ~480 bits as a dd (relative error ~2^-470).
[[nodiscard]] inline dd to_dd(const MP& a) {
    if (a.zero()) return dd(float256::zero());
    BigInt m = a.m;
    long e = a.e;
    const int bl = m.bit_length();
    if (bl > 480) { m.shr(bl - 480); e += bl - 480; }
    dd acc{float256::zero()};
    for (int i = m.limb_count() - 1; i >= 0; --i)
        acc = detail_dd::ldexp(acc, 32) + dd(float256::from_uint(m.limb(i)));
    // apply 2^e in two steps (e may be far outside one ldexp's range of the parts)
    acc = detail_dd::ldexp(acc, static_cast<int>(e));
    return a.neg ? -acc : acc;
}

[[nodiscard]] inline MP neg(MP a) { if (!a.zero()) a.neg = !a.neg; return a; }
[[nodiscard]] inline MP ldexp(MP a, long k) { if (!a.zero()) a.e += k; return a; }

[[nodiscard]] inline MP add(const MP& x, const MP& y, int P) {
    if (x.zero()) { MP r = y; norm(r, P); return r; }
    if (y.zero()) { MP r = x; norm(r, P); return r; }
    // a value far below the other's precision only matters as a sticky bit
    if (x.top() - y.top() > P + 64) { MP r = x; norm(r, P); return r; }
    if (y.top() - x.top() > P + 64) { MP r = y; norm(r, P); return r; }
    MP a = x, b = y;
    // drop bits far below the working precision before aligning
    const long lo = (a.top() > b.top() ? a.top() : b.top()) - P - 96;
    for (MP* v : {&a, &b})
        if (v->e < lo) { v->m.shr(static_cast<int>(lo - v->e)); v->e = lo; }
    const long emin = a.e < b.e ? a.e : b.e;
    a.m.shl(static_cast<int>(a.e - emin)); a.e = emin;
    b.m.shl(static_cast<int>(b.e - emin)); b.e = emin;
    if (a.neg == b.neg) a.m.add(b.m);
    else if (b.m < a.m) a.m.sub(b.m);
    else { BigInt t = b.m; t.sub(a.m); a.m = t; a.neg = b.neg; }
    norm(a, P);
    return a;
}
[[nodiscard]] inline MP sub(const MP& x, const MP& y, int P) { return add(x, neg(y), P); }

[[nodiscard]] inline MP mul(const MP& x, const MP& y, int P) {
    MP r;
    if (x.zero() || y.zero()) return r;
    r.m = x.m; r.m.mul(y.m);
    r.e = x.e + y.e;
    r.neg = x.neg != y.neg;
    norm(r, P);
    return r;
}
[[nodiscard]] inline MP mul_u(const MP& x, std::uint32_t k, int P) {
    MP r = x; r.m.mul_u32(k); if (k == 0) r = MP{}; norm(r, P); return r;
}

[[nodiscard]] inline MP div(const MP& x, const MP& y, int P) {
    MP r;
    if (x.zero()) return r;
    BigInt a = x.m;
    const int sh = P + 2 + y.m.bit_length() - x.m.bit_length();
    long e = x.e - y.e;
    if (sh > 0) { a.shl(sh); e -= sh; }
    (void)a.divmod(y.m);
    r.m = a; r.e = e; r.neg = x.neg != y.neg;
    norm(r, P);
    return r;
}
[[nodiscard]] inline MP div_u(const MP& x, std::uint32_t k, int P) { return div(x, from_u64(k), P); }

[[nodiscard]] inline int cmp_abs(const MP& a, const MP& b) {
    if (a.zero() || b.zero()) return a.zero() ? (b.zero() ? 0 : -1) : 1;
    if (a.top() != b.top()) return a.top() < b.top() ? -1 : 1;
    MP d = sub(MP{false, a.m, a.e}, MP{false, b.m, b.e}, 1 << 20);
    return d.zero() ? 0 : (d.neg ? -1 : 1);
}

// Newton steps from the dd value: y ← (y + x/y)/2.
[[nodiscard]] inline MP sqrt(const MP& x, int P) {
    if (x.zero()) return x;
    // scale x to near 1 so the dd seed does not overflow
    long k = x.top() / 2 * 2;
    MP xs = ldexp(x, -k);
    MP y = from_dd(detail_dd::sqrt(to_dd(xs)));
    for (int prec = 470; prec < P + 64; prec *= 2) {
        const int Q = (2 * prec < P + 64 ? 2 * prec : P + 64) + 32;
        y = ldexp(add(y, div(xs, y, Q), Q), -1);
    }
    y = ldexp(y, k / 2);
    norm(y, P);
    return y;
}

// Fixed-point helpers for constants: Σ of a series of positive terms.
struct Consts {
    std::mutex mu;
    int P_ln2 = 0, P_pi = 0;
    MP ln2, pi;
};
inline Consts& consts() { static Consts c; return c; }

// ln 2 = 2 atanh(1/3) = 2 Σ 3^{−(2k+1)}/(2k+1).
[[nodiscard]] inline MP ln2(int P) {
    Consts& c = consts();
    std::lock_guard<std::mutex> g(c.mu);
    if (c.P_ln2 >= P) { MP r = c.ln2; norm(r, P); return r; }
    const int Q = P + 64;
    MP sum, p = div_u(from_u64(1), 3, Q);
    for (std::uint32_t k = 0;; ++k) {
        const MP t = div_u(p, 2 * k + 1, Q);
        sum = add(sum, t, Q);
        if (t.zero() || t.top() < sum.top() - Q) break;
        p = div_u(p, 9, Q);
    }
    c.ln2 = ldexp(sum, 1); c.P_ln2 = P;
    MP r = c.ln2; norm(r, P); return r;
}

// Word k of the 2/π fraction bits (two_over_pi.hpp, then two_over_pi_ext.hpp).
inline constexpr int kTwoOverPiAllWords = detail_dd::kTwoOverPiWords + detail_dd::kTwoOverPiExtWords;
[[nodiscard]] inline std::uint64_t two_over_pi_word(int k) {
    return k < detail_dd::kTwoOverPiWords ? detail_dd::kTwoOverPi[k] : detail_dd::kTwoOverPiExt[k - detail_dd::kTwoOverPiWords];
}

// π = 2 / (2/π) from the 2/π table (395264 bits).
[[nodiscard]] inline MP pi(int P) {
    Consts& c = consts();
    std::lock_guard<std::mutex> g(c.mu);
    if (c.P_pi >= P) { MP r = c.pi; norm(r, P); return r; }
    const int words = (P + 127) / 64;
    std::vector<std::uint32_t> limbs;
    limbs.reserve(static_cast<std::size_t>(2 * words));
    for (int i = words - 1; i >= 0; --i) {
        limbs.push_back(static_cast<std::uint32_t>(two_over_pi_word(i)));
        limbs.push_back(static_cast<std::uint32_t>(two_over_pi_word(i) >> 32));
    }
    MP t; t.m = BigInt::from_limbs(limbs); t.e = -64L * words;  // 2/π
    c.pi = div(from_u64(2), t, P + 32); c.P_pi = P;
    MP r = c.pi; norm(r, P); return r;
}

// 2/π to P bits straight from the table (P <= 64·kTwoOverPiAllWords − 127).
[[nodiscard]] inline MP two_over_pi(int P) {
    const int words = (P + 127) / 64;
    std::vector<std::uint32_t> limbs;
    limbs.reserve(static_cast<std::size_t>(2 * words));
    for (int i = words - 1; i >= 0; --i) {
        limbs.push_back(static_cast<std::uint32_t>(two_over_pi_word(i)));
        limbs.push_back(static_cast<std::uint32_t>(two_over_pi_word(i) >> 32));
    }
    MP t; t.m = BigInt::from_limbs(limbs); t.e = -64L * words;
    norm(t, P); return t;
}

// e^x for |x| < 2^30: x = k ln2 + r, r/2^s by Taylor, then s squarings.
[[nodiscard]] inline MP exp(const MP& x, int P) {
    if (x.zero()) return from_u64(1);
    const double xd = to_dd(x).hi.to_double();
    const long k = std::lround(xd / 0.6931471805599453);
    const int s = static_cast<int>(std::sqrt(static_cast<double>(P))) / 2 + 4;
    const int Q = P + 2 * s + 64 + static_cast<int>(std::log2(std::fabs(static_cast<double>(k)) + 2));
    MP r = sub(x, mul(from_u64(static_cast<std::uint64_t>(k < 0 ? -k : k), k < 0), ln2(Q), Q), Q);
    r = ldexp(r, -s);
    MP sum = from_u64(1), t = from_u64(1);
    for (std::uint32_t n = 1;; ++n) {
        t = div_u(mul(t, r, Q), n, Q);
        sum = add(sum, t, Q);
        if (t.zero() || t.top() < sum.top() - Q) break;
    }
    for (int i = 0; i < s; ++i) sum = mul(sum, sum, Q);
    sum = ldexp(sum, k);
    norm(sum, P);
    return sum;
}

// ln x, x > 0: Newton y ← y + x e^{−y} − 1 from the dd log.
[[nodiscard]] inline MP log(const MP& x, int P) {
    const long k = x.top();
    const MP xs = ldexp(x, -k);  // in [1/2, 1)
    MP y = from_dd(detail_dd::log(to_dd(xs)));
    for (int prec = 460; prec < P + 64; prec *= 2) {
        const int Q = (2 * prec < P + 64 ? 2 * prec : P + 64) + 32;
        y = add(y, sub(mul(xs, exp(neg(y), Q), Q), from_u64(1), Q), Q);
    }
    if (k != 0) y = add(y, mul(from_u64(static_cast<std::uint64_t>(k < 0 ? -k : k), k < 0), ln2(P + 64), P + 64), P + 64);
    norm(y, P);
    return y;
}

// sin x, cos x for |x| <= 1 (callers reduce): Taylor at x/2^s, then doubling.
inline void sincos(const MP& x, int P, MP& s, MP& c) {
    const int h = static_cast<int>(std::sqrt(static_cast<double>(P))) / 2 + 2;
    const int Q = P + 2 * h + 64;
    const MP r = ldexp(x, -h), r2 = mul(r, r, Q);
    MP ss = r, t = r;
    for (std::uint32_t n = 1;; ++n) {  // sin
        t = neg(div_u(mul(t, r2, Q), (2 * n) * (2 * n + 1), Q));
        ss = add(ss, t, Q);
        if (t.zero() || t.top() < ss.top() - Q) break;
    }
    // cos = 1 − 2 sin²(r/2)… use cos r = sqrt(1 − sin² r) (r small, cos > 0)
    MP cc = sqrt(sub(from_u64(1), mul(ss, ss, Q), Q), Q);
    for (int i = 0; i < h; ++i) {
        const MP s2 = ldexp(mul(ss, cc, Q), 1);
        cc = sub(mul(cc, cc, Q), mul(ss, ss, Q), Q);
        ss = s2;
    }
    s = ss; c = cc; norm(s, P); norm(c, P);
}

// B_2k/(2k)! for k = 1..n from the tangent numbers T_k (Brent–Harvey):
// B_2k = (−1)^{k−1} 2k T_k / (2^{2k}(2^{2k} − 1)).
struct Bern {
    std::mutex mu;
    std::vector<BigInt> T;  // T[k], k >= 1
};
inline Bern& bern() { static Bern b; return b; }
[[nodiscard]] inline MP bern_fact(int k, int P) {
    Bern& B = bern();
    BigInt Tk;
    {
        std::lock_guard<std::mutex> g(B.mu);
        if (static_cast<int>(B.T.size()) <= k) {
            // recompute up to n = 2k (amortized)
            const int n = 2 * k + 8;
            std::vector<BigInt> T(static_cast<std::size_t>(n + 1));
            T[1] = BigInt(1);
            for (int i = 2; i <= n; ++i) { T[static_cast<std::size_t>(i)] = T[static_cast<std::size_t>(i - 1)]; T[static_cast<std::size_t>(i)].mul_u32(static_cast<std::uint32_t>(i - 1)); }
            for (int i = 2; i <= n; ++i)
                for (int j = i; j <= n; ++j) {
                    BigInt a = T[static_cast<std::size_t>(j - 1)];
                    a.mul_u32(static_cast<std::uint32_t>(j - i));
                    BigInt b = T[static_cast<std::size_t>(j)];
                    b.mul_u32(static_cast<std::uint32_t>(j - i + 2));
                    a.add(b);
                    T[static_cast<std::size_t>(j)] = a;
                }
            B.T = std::move(T);
        }
        Tk = B.T[static_cast<std::size_t>(k)];
    }
    // denominator 2^{2k}(2^{2k} − 1)(2k)! / (2k)
    BigInt den(1);
    for (int i = 1; i < 2 * k; ++i) den.mul_u32(static_cast<std::uint32_t>(i));  // (2k−1)!
    BigInt p2(1); p2.shl(2 * k); p2.sub(BigInt(1));
    den.mul(p2);
    MP num; num.m = Tk; num.neg = (k % 2) == 0;
    MP d; d.m = den; d.e = 2 * k;
    return div(num, d, P);
}

// ln Γ(x), x > 0: shift to y >= 0.12 P + 16, then Stirling.
[[nodiscard]] inline MP lgamma_pos(const MP& x, int P) {
    const int Q = P + 64;
    const double yd = 0.12 * P + 16;
    MP y = x, prod = from_u64(1);
    while (to_dd(y).hi.to_double() < yd) {
        prod = mul(prod, y, Q);
        y = add(y, from_u64(1), Q);
    }
    const MP half = ldexp(from_u64(1), -1);
    const MP ly = log(y, Q);
    MP s = sub(mul(sub(y, half, Q), ly, Q), y, Q);
    s = add(s, ldexp(log(ldexp(pi(Q), 1), Q), -1), Q);
    const MP iy = div(from_u64(1), y, Q), iy2 = mul(iy, iy, Q);
    MP p = iy;  // y^{1−2k}
    for (int k = 1;; ++k) {
        // B_2k/(2k(2k−1)) y^{1−2k} = bern_fact(k)·(2k−2)! y^{1−2k}
        MP c = bern_fact(k, Q);
        for (int i = 2; i <= 2 * k - 2; ++i) c = mul_u(c, static_cast<std::uint32_t>(i), Q);
        const MP t = mul(c, p, Q);
        s = add(s, t, Q);
        if (t.zero() || t.top() < s.top() - Q) break;
        p = mul(p, iy2, Q);
    }
    if (!(prod.m.bit_length() == 1 && prod.e == 0)) s = sub(s, log(prod, Q), Q);
    norm(s, P);
    return s;
}

} // namespace base256::detail_mp
