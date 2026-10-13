// Differential test of base256::float256 against GNU MPFR.
// MPFR is configured to emulate binary256 exactly: precision 237,
// emin/emax matching binary256, and mpfr_subnormalize after each op.
//
// Usage: mpfr_diff [cases_per_op] [seed]
#include <base256/base256.hpp>
#include <mpfr.h>

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <cstring>
#include <random>
#include <string>

using base256::float256;
using W = float256::storage_type;

namespace {

constexpr int kP = 237;
constexpr long kEminN = -262142;                // min normal exponent
constexpr long kMpfrEmin = kEminN - kP + 2;     // -262377
constexpr long kMpfrEmax = 262144;

std::mt19937_64 rng;

std::uint64_t r64() { return rng(); }

// Decode float256 -> mpfr (exact).
void to_mpfr(mpfr_t out, float256 x) {
    const W w = x.raw_bits();
    const bool s = (w[3] >> 63) != 0;
    const long be = static_cast<long>((w[3] >> 44) & 0x7FFFF);
    const std::uint64_t f3 = w[3] & 0xFFFFFFFFFFFULL;
    if (be == 0x7FFFF) {
        if (f3 == 0 && w[2] == 0 && w[1] == 0 && w[0] == 0) mpfr_set_inf(out, s ? -1 : 1);
        else mpfr_set_nan(out);
        return;
    }
    mpz_t m;
    mpz_init(m);
    mpz_set_ui(m, f3);
    for (int i = 2; i >= 0; --i) {
        mpz_mul_2exp(m, m, 64);
        mpz_add_ui(m, m, w[i]);
    }
    long e;
    if (be == 0) {
        e = kEminN - 236;
    } else {
        mpz_setbit(m, 236);
        e = be - 262143 - 236;
    }
    mpfr_set_z_2exp(out, m, e, MPFR_RNDN); // exact: m < 2^237
    if (s) mpfr_neg(out, out, MPFR_RNDN);
    if (mpfr_zero_p(out) && s) mpfr_setsign(out, out, 1, MPFR_RNDN);
    mpz_clear(m);
}

// Encode an mpfr value already rounded to binary256 -> float256 bits.
float256 from_mpfr(mpfr_t x) {
    if (mpfr_nan_p(x)) return float256::qnan();
    const bool s = mpfr_signbit(x) != 0;
    if (mpfr_inf_p(x)) return float256::inf(s);
    if (mpfr_zero_p(x)) return float256::zero(s);
    mpz_t m;
    mpz_init(m);
    long e = mpfr_get_z_2exp(m, x);  // x = m * 2^e
    mpz_abs(m, m);
    const long L = static_cast<long>(mpz_sizeinbase(m, 2));
    const long E = e + L - 1;
    long biased;
    if (E >= kEminN) {
        if (L < 237) mpz_mul_2exp(m, m, static_cast<unsigned long>(237 - L));
        else if (L > 237) mpz_fdiv_q_2exp(m, m, static_cast<unsigned long>(L - 237));
        mpz_clrbit(m, 236);
        biased = E + 262143;
    } else {
        const long sh = e - (kEminN - 236);
        if (sh >= 0) mpz_mul_2exp(m, m, static_cast<unsigned long>(sh));
        else mpz_fdiv_q_2exp(m, m, static_cast<unsigned long>(-sh));
        biased = 0;
    }
    W w{};
    for (int i = 0; i < 4; ++i) {
        mpz_t t; mpz_init(t);
        mpz_fdiv_r_2exp(t, m, 64);
        w[i] = mpz_get_ui(t);
        mpz_fdiv_q_2exp(m, m, 64);
        mpz_clear(t);
    }
    w[3] |= static_cast<std::uint64_t>(biased) << 44;
    if (s) w[3] |= 1ULL << 63;
    mpz_clear(m);
    return float256::from_bits(w);
}

float256 from_mpfr_ld(mpfr_t tmp, long double v) {
    mpfr_set_ld(tmp, v, MPFR_RNDN); // exact at 237 bits
    return from_mpfr(tmp);
}

// Random operand generator with emphasis on hard cases.
std::uint64_t frac_word() {
    switch (r64() % 6) {
        case 0: return 0;
        case 1: return ~0ULL;
        case 2: return 1ULL << (r64() % 64);
        case 3: return ~(1ULL << (r64() % 64));
        default: return r64();
    }
}

long rand_biased_exp() {
    switch (r64() % 10) {
        case 0: return 0;                                  // subnormal / zero
        case 1: return 1 + static_cast<long>(r64() % 300);  // near min normal
        case 2: return 0x7FFFE - static_cast<long>(r64() % 300); // near max
        case 3: return 0x7FFFF;                            // inf / nan (rare-ish)
        case 4: case 5: case 6:
            return 262143 - 300 + static_cast<long>(r64() % 600); // around 1
        default: return 1 + static_cast<long>(r64() % 0x7FFFE);
    }
}

float256 rand_f256() {
    W w{frac_word(), frac_word(), frac_word(), frac_word() & 0xFFFFFFFFFFFULL};
    long be = rand_biased_exp();
    if (be == 0x7FFFF && (r64() % 4)) be = 0x7FFFE; // keep specials rarer
    w[3] |= static_cast<std::uint64_t>(be) << 44;
    if (r64() & 1) w[3] |= 1ULL << 63;
    return float256::from_bits(w);
}

// A second operand correlated with the first: nearby exponent, often a
// nearly identical significand (catastrophic cancellation, sticky cases).
float256 rand_near(float256 a) {
    if (r64() % 3 == 0) return rand_f256();
    W w = a.raw_bits();
    long be = static_cast<long>((w[3] >> 44) & 0x7FFFF);
    if (be == 0x7FFFF) return rand_f256();
    long d = static_cast<long>(r64() % 300) - 20;
    if (r64() % 2) d = static_cast<long>(r64() % 4);
    long nb = be - d;
    if (nb < 0) nb = 0;
    if (nb > 0x7FFFE) nb = 0x7FFFE;
    w[3] = (w[3] & ~(0x7FFFFULL << 44)) | (static_cast<std::uint64_t>(nb) << 44);
    int k = static_cast<int>(r64() % 4);
    switch (r64() % 4) {
        case 0: w[k] ^= 1ULL << (r64() % 64); break;
        case 1: w[k] = frac_word(); break;
        case 2: w[0] += 1; break;
        default: break;
    }
    w[3] = (w[3] & ~((0x7FFFFULL << 44) | (1ULL << 63))) |
           (static_cast<std::uint64_t>(nb) << 44);
    if (r64() & 1) w[3] |= 1ULL << 63;
    return float256::from_bits(w);
}

bool same(float256 a, float256 b) {
    if (a.isnan() && b.isnan()) return true;
    return a.raw_bits() == b.raw_bits();
}

std::string hex(float256 x) {
    char buf[96];
    const W w = x.raw_bits();
    std::snprintf(buf, sizeof buf, "%016llx_%016llx_%016llx_%016llx",
                  (unsigned long long)w[3], (unsigned long long)w[2],
                  (unsigned long long)w[1], (unsigned long long)w[0]);
    return buf;
}

struct Ctx {
    mpfr_t a, b, c, r;
    Ctx() { mpfr_inits2(kP, a, b, c, r, (mpfr_ptr)0); }
    ~Ctx() { mpfr_clears(a, b, c, r, (mpfr_ptr)0); }
};

long g_fail = 0;
long g_total = 0;

const char* kModeName[5] = {"RNE", "RNA", "RTZ", "RUP", "RDN"};
const char* kOpName[6] = {"add", "sub", "mul", "div", "fma", "sqrt"};

int mp_op(int op, mpfr_ptr r, Ctx& k, mpfr_rnd_t rnd) {
    switch (op) {
        case 0: return mpfr_add(r, k.a, k.b, rnd);
        case 1: return mpfr_sub(r, k.a, k.b, rnd);
        case 2: return mpfr_mul(r, k.a, k.b, rnd);
        case 3: return mpfr_div(r, k.a, k.b, rnd);
        case 4: return mpfr_fma(r, k.a, k.b, k.c, rnd);
        default: return mpfr_sqrt(r, k.a, rnd);
    }
}

float256 f_op(int op, float256 a, float256 b, float256 c) {
    switch (op) {
        case 0: return a + b;
        case 1: return a - b;
        case 2: return a * b;
        case 3: return a / b;
        case 4: return base256::fma(a, b, c);
        default: return base256::sqrt(a);
    }
}

void full_range() { mpfr_set_emin(mpfr_get_emin_min()); mpfr_set_emax(mpfr_get_emax_max()); }
void b256_range() { mpfr_set_emin(kMpfrEmin); mpfr_set_emax(kMpfrEmax); }

struct Ref { float256 v; unsigned flags; };

// Is a nonzero value, rounded to 237 bits with unbounded exponent, below 2^Emin?
bool tiny_237(int op, Ctx& k, mpfr_rnd_t rnd) {
    full_range();
    mpfr_t t;
    mpfr_init2(t, kP);
    mp_op(op, t, k, rnd);
    const bool tiny = mpfr_regular_p(t) && mpfr_get_exp(t) <= kEminN;
    mpfr_clear(t);
    b256_range();
    return tiny;
}

// Round-to-nearest-away reference: MPFR has no subnormal-aware RNDNA, so take
// a 4000-bit truncated result and round it to the binary256 grid by hand.
Ref ref_rna(int op, Ctx& k, bool& handled_nan) {
    full_range();
    mpfr_t big, sc;
    mpfr_inits2(4000, big, sc, (mpfr_ptr)0);
    mpfr_clear_flags();
    const int t = mp_op(op, big, k, MPFR_RNDZ);
    Ref r{float256::zero(), 0};
    handled_nan = false;
    if (mpfr_nan_p(big) || mpfr_inf_p(big) || mpfr_zero_p(big)) {
        handled_nan = mpfr_nan_p(big);
        r.v = from_mpfr(big);
        if (mpfr_divby0_p()) r.flags |= base256::fe_divbyzero;
        mpfr_clears(big, sc, (mpfr_ptr)0);
        b256_range();
        return r;
    }
    const bool neg = mpfr_signbit(big) != 0;
    mpfr_abs(big, big, MPFR_RNDN);
    const long E = mpfr_get_exp(big) - 1;
    auto round_at = [&](long q, mpz_t n) -> bool {  // returns inexact
        mpfr_mul_2si(sc, big, -q, MPFR_RNDN);         // exact
        mpfr_get_z(n, sc, MPFR_RNDZ);
        mpfr_sub_z(sc, sc, n, MPFR_RNDN);             // fraction, exact
        const bool ix = !mpfr_zero_p(sc) || t != 0;
        if (mpfr_cmp_d(sc, 0.5) >= 0) mpz_add_ui(n, n, 1);
        return ix;
    };
    mpz_t n;
    mpz_init(n);
    // tininess after rounding (unbounded exponent)
    const bool ix_unb = round_at(E - 236, n);
    mpfr_t v;
    mpfr_init2(v, 400);
    mpfr_set_z_2exp(v, n, E - 236, MPFR_RNDN);
    const bool tiny = mpfr_get_exp(v) <= kEminN;
    const long q = std::max(E - 236, kEminN - 236);
    const bool ix = round_at(q, n);
    (void)ix_unb;
    mpfr_set_z_2exp(v, n, q, MPFR_RNDN);
    if (ix) r.flags |= base256::fe_inexact;
    if (ix && tiny) r.flags |= base256::fe_underflow;
    // max finite = (2^237 - 1) * 2^(262143 - 236)
    if (mpfr_get_exp(v) > kMpfrEmax) {
        r.flags |= base256::fe_overflow | base256::fe_inexact;
        r.v = float256::inf(neg);
    } else {
        mpfr_t w;
        mpfr_init2(w, kP);
        mpfr_set(w, v, MPFR_RNDN);  // exact (<= 237 significant bits)
        if (neg) mpfr_neg(w, w, MPFR_RNDN);
        r.v = from_mpfr(w);
        mpfr_clear(w);
    }
    mpfr_clear(v);
    mpz_clear(n);
    mpfr_clears(big, sc, (mpfr_ptr)0);
    b256_range();
    return r;
}

Ref ref_mpfr(int op, Ctx& k, mpfr_rnd_t rnd, bool& is_nan) {
    b256_range();
    mpfr_clear_flags();
    int t = mp_op(op, k.r, k, rnd);
    t = mpfr_check_range(k.r, t, rnd);
    t = mpfr_subnormalize(k.r, t, rnd);
    Ref r{from_mpfr(k.r), 0};
    is_nan = mpfr_nan_p(k.r);
    const bool inexact = t != 0 || mpfr_overflow_p();
    if (inexact) r.flags |= base256::fe_inexact;
    if (mpfr_overflow_p()) r.flags |= base256::fe_overflow;
    if (mpfr_divby0_p()) r.flags |= base256::fe_divbyzero;
    if (inexact && !mpfr_overflow_p() && tiny_237(op, k, rnd)) r.flags |= base256::fe_underflow;
    return r;
}

void report(const char* op, const char* mode, float256 a, float256 b, float256 c, float256 got,
            float256 want, unsigned gf, unsigned wf) {
    if (++g_fail <= 20) {
        std::printf("MISMATCH %s %s\n  a=%s\n  b=%s\n  c=%s\n  got =%s flags=%02x\n  want=%s flags=%02x\n", op,
                    mode, hex(a).c_str(), hex(b).c_str(), hex(c).c_str(), hex(got).c_str(), gf,
                    hex(want).c_str(), wf);
    }
}

} // namespace

int main(int argc, char** argv) {
    const long n = argc > 1 ? std::atol(argv[1]) : 200000;
    const unsigned long long seed = argc > 2 ? std::strtoull(argv[2], nullptr, 0) : 12345;
    // optional third argument: comma-free list of modes, e.g. "0" or "01234"
    const std::string modes = argc > 3 ? argv[3] : "01234";
    rng.seed(seed);
    b256_range();
    Ctx k;

    for (long i = 0; i < 1000; ++i) {
        float256 a = rand_f256();
        if (a.isnan()) continue;
        to_mpfr(k.a, a);
        if (!same(from_mpfr(k.a), a)) { std::printf("codec failure %s\n", hex(a).c_str()); return 2; }
    }

    const base256::round_mode bm[5] = {base256::round_mode::to_nearest_even,
                                       base256::round_mode::to_nearest_away,
                                       base256::round_mode::toward_zero,
                                       base256::round_mode::upward,
                                       base256::round_mode::downward};
    const mpfr_rnd_t mm[5] = {MPFR_RNDN, MPFR_RNDNA, MPFR_RNDZ, MPFR_RNDU, MPFR_RNDD};

    for (long i = 0; i < n; ++i) {
        const float256 a = rand_f256();
        const float256 b = rand_near(a);
        const float256 c = rand_near(r64() & 1 ? a : b);
        to_mpfr(k.a, a);
        to_mpfr(k.b, b);
        to_mpfr(k.c, c);
        const bool any_snan = a.issignaling() || b.issignaling() || c.issignaling();
        for (char mc : modes) {
            const int m = mc - '0';
            for (int op = 0; op < 6; ++op) {
                const bool in_nan = op == 5 ? a.isnan()
                                  : op == 4 ? (a.isnan() || b.isnan() || c.isnan())
                                            : (a.isnan() || b.isnan());
                const bool in_snan = op == 5 ? a.issignaling()
                                   : op == 4 ? any_snan
                                             : (a.issignaling() || b.issignaling());
                bool rnan = false;
                Ref want = m == 1 ? ref_rna(op, k, rnan) : ref_mpfr(op, k, mm[m], rnan);
                if (in_nan) want.flags = 0;  // NaN in, NaN out: no other flags
                if (in_snan || (rnan && !in_nan)) want.flags |= base256::fe_invalid;
                base256::fenv::set_round(bm[m]);
                base256::fenv::clear_flags();
                const float256 got = f_op(op, a, b, c);
                const unsigned gf = base256::fenv::test_flags();
                base256::fenv::set_round(base256::round_mode::to_nearest_even);
                ++g_total;
                if (!same(got, want.v) || gf != want.flags)
                    report(kOpName[op], kModeName[m], a, b, c, got, want.v, gf, want.flags);
            }
            // to_double / to_long_double in the same mode (MPFR has no RNDNA here)
            if (m != 1 && !a.isnan()) {
                base256::fenv::set_round(bm[m]);
                const double gd = a.to_double();
                const long double gl = a.to_long_double();
                base256::fenv::set_round(base256::round_mode::to_nearest_even);
                const double wd = mpfr_get_d(k.a, mm[m]);
                const long double wl = mpfr_get_ld(k.a, mm[m]);
                g_total += 2;
                if (std::memcmp(&gd, &wd, sizeof gd) != 0 && ++g_fail <= 20)
                    std::printf("MISMATCH to_double %s a=%s got=%a want=%a\n", kModeName[m], hex(a).c_str(), gd, wd);
                if (!(gl == wl && std::signbit(gl) == std::signbit(wl)) && ++g_fail <= 20)
                    std::printf("MISMATCH to_long_double %s a=%s got=%La want=%La\n", kModeName[m], hex(a).c_str(), gl, wl);
            }
        }
    }
    // exact from_long_double round trip
    for (long i = 0; i < n / 4; ++i) {
        const float256 a = rand_f256();
        if (a.isnan()) continue;
        to_mpfr(k.a, a);
        const long double ld = mpfr_get_ld(k.a, MPFR_RNDN);
        ++g_total;
        if (!same(float256(ld), from_mpfr_ld(k.c, ld)) && ++g_fail <= 20)
            std::printf("MISMATCH from_long_double %La\n", ld);
    }
    std::printf("mpfr_diff seed=%llu cases=%ld modes=%s checks=%ld failures=%ld\n", seed, n,
                modes.c_str(), g_total, g_fail);
    return g_fail ? 1 : 0;
}
