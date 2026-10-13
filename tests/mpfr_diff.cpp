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

using MpfrRnd = mpfr_rnd_t;

float256 finish(Ctx& k, int t, MpfrRnd rnd) {
    t = mpfr_check_range(k.r, t, rnd);
    mpfr_subnormalize(k.r, t, rnd);
    return from_mpfr(k.r);
}

long g_fail = 0;
long g_total = 0;

void report(const char* op, float256 a, float256 b, float256 c, float256 got, float256 want) {
    if (++g_fail <= 20) {
        std::printf("MISMATCH %s\n  a=%s\n  b=%s\n  c=%s\n  got =%s\n  want=%s\n", op,
                    hex(a).c_str(), hex(b).c_str(), hex(c).c_str(), hex(got).c_str(),
                    hex(want).c_str());
    }
}

} // namespace

int main(int argc, char** argv) {
    const long n = argc > 1 ? std::atol(argv[1]) : 200000;
    const unsigned long long seed = argc > 2 ? std::strtoull(argv[2], nullptr, 0) : 12345;
    rng.seed(seed);
    mpfr_set_emin(kMpfrEmin);
    mpfr_set_emax(kMpfrEmax);
    Ctx k;
    const MpfrRnd rnd = MPFR_RNDN;

    // Self-check of the codec.
    for (long i = 0; i < 1000; ++i) {
        float256 a = rand_f256();
        if (a.isnan()) continue;
        to_mpfr(k.a, a);
        if (!same(from_mpfr(k.a), a)) { std::printf("codec failure %s\n", hex(a).c_str()); return 2; }
    }

    long per[7] = {};
    for (long i = 0; i < n; ++i) {
        const float256 a = rand_f256();
        const float256 b = rand_near(a);
        const float256 c = rand_near(r64() & 1 ? a : b);
        to_mpfr(k.a, a);
        to_mpfr(k.b, b);
        to_mpfr(k.c, c);
        const auto chk = [&](const char* op, int idx, float256 got, float256 want) {
            ++g_total;
            ++per[idx];
            if (!same(got, want)) report(op, a, b, c, got, want);
        };
        chk("add", 0, a + b, finish(k, mpfr_add(k.r, k.a, k.b, rnd), rnd));
        chk("sub", 1, a - b, finish(k, mpfr_sub(k.r, k.a, k.b, rnd), rnd));
        chk("mul", 2, a * b, finish(k, mpfr_mul(k.r, k.a, k.b, rnd), rnd));
        chk("div", 3, a / b, finish(k, mpfr_div(k.r, k.a, k.b, rnd), rnd));
        chk("fma", 4, base256::fma(a, b, c), finish(k, mpfr_fma(k.r, k.a, k.b, k.c, rnd), rnd));
        chk("sqrt", 5, base256::sqrt(a), finish(k, mpfr_sqrt(k.r, k.a, rnd), rnd));
        // to_double vs mpfr_get_d (binary64 rounding incl. subnormals)
        if (!a.isnan()) {
            ++g_total; ++per[6];
            const double got = a.to_double();
            const double want = mpfr_get_d(k.a, rnd);
            if (std::memcmp(&got, &want, sizeof got) != 0) {
                if (++g_fail <= 20) std::printf("MISMATCH to_double a=%s got=%a want=%a\n", hex(a).c_str(), got, want);
            }
        }
    }
    // to_long_double vs mpfr_get_ld, and the exact from_long_double round-trip.
    for (long i = 0; i < n / 4; ++i) {
        const float256 a = rand_f256();
        if (a.isnan()) continue;
        to_mpfr(k.a, a);
        ++g_total;
        const long double got = a.to_long_double();
        const long double want = mpfr_get_ld(k.a, rnd);
        if (!(got == want || (got != got && want != want)) || std::signbit(got) != std::signbit(want)) {
            if (++g_fail <= 20) std::printf("MISMATCH to_long_double a=%s got=%La want=%La\n", hex(a).c_str(), got, want);
        }
        ++g_total;
        if (!same(float256(want), from_mpfr_ld(k.c, want))) {
            if (++g_fail <= 20) std::printf("MISMATCH from_long_double %La\n", want);
        }
    }
    std::printf("mpfr_diff seed=%llu cases/op=%ld total=%ld failures=%ld\n", seed, n, g_total, g_fail);
    return g_fail ? 1 : 0;
}
