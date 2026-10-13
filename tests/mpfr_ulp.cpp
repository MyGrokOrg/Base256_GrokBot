// Measures elementary-function error in ulps against MPFR (exact reference
// at 1024 bits). Usage: mpfr_ulp [cases_per_fn] [seed] [bound_check]
// With bound_check=1 the program exits nonzero if any function exceeds its
// documented bound (see README "Elementary function accuracy").
#include <base256/base256.hpp>
#include <mpfr.h>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <random>
#include <string>
#include <vector>

using base256::float256;
using W = float256::storage_type;
namespace {
constexpr long kEminN = -262142;
std::mt19937_64 rng;
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


// Error of y (a float256) relative to exact ref, in ulps of the binary256 format at ref.
double ulp_err(float256 y, mpfr_t ref) {
    mpfr_t a, u;
    mpfr_inits2(1100, a, u, (mpfr_ptr)0);
    if (mpfr_nan_p(ref) || mpfr_inf_p(ref)) { mpfr_clears(a, u, (mpfr_ptr)0); return (y.isnan() == (mpfr_nan_p(ref) != 0)) ? 0 : 1e9; }
    to_mpfr(a, y);
    if (mpfr_inf_p(a) || mpfr_nan_p(a)) { mpfr_clears(a, u, (mpfr_ptr)0); return 1e9; }
    mpfr_sub(a, a, ref, MPFR_RNDN);
    long e = mpfr_zero_p(ref) ? kEminN : mpfr_get_exp(ref) - 1;  // ref in [2^e, 2^(e+1))
    if (e < kEminN) e = kEminN;
    mpfr_set_ui_2exp(u, 1, e - 236, MPFR_RNDN);
    mpfr_div(a, a, u, MPFR_RNDN);
    const double r = std::fabs(mpfr_get_d(a, MPFR_RNDN));
    mpfr_clears(a, u, (mpfr_ptr)0);
    return r;
}

float256 rnd_in(double lo, double hi, bool logscale) {
    // random 237-bit value in [lo, hi] (linear or log-uniform magnitude, same sign)
    std::uniform_real_distribution<double> d(0, 1);
    float256 t = float256(static_cast<double>(rng() >> 11) * 0x1p-53) +
                 float256(static_cast<double>(rng() >> 11) * 0x1p-106) +
                 float256(static_cast<double>(rng() >> 11) * 0x1p-159) +
                 float256(static_cast<double>(rng() >> 11) * 0x1p-212) +
                 float256(static_cast<double>(rng() >> 11) * 0x1p-265);
    if (!logscale) return float256(lo) + (float256(hi) - float256(lo)) * t;
    // log scale: lo/hi are binary exponents
    const int ex = static_cast<int>(std::floor(lo + (hi - lo) * d(rng)));
    return base256::ldexp(float256::one() + t, ex);
}

// x = M·2^E (M a 237-bit integer) that is unusually close to a multiple of
// π/2: M is a multiple of the largest continued-fraction denominator q <= 2^237
// of α = frac(2^E · 2/π), so M·α is within about M/q² of an integer.
float256 hard_trig(long E) {
    const long prec = (E > 0 ? E : 0) + 1000;
    mpfr_t c; mpfr_init2(c, prec);
    mpfr_const_pi(c, MPFR_RNDN); mpfr_ui_div(c, 2, c, MPFR_RNDN);
    if (E >= 0) mpfr_mul_2ui(c, c, static_cast<unsigned long>(E), MPFR_RNDN);
    else mpfr_div_2ui(c, c, static_cast<unsigned long>(-E), MPFR_RNDN);
    mpfr_frac(c, c, MPFR_RNDN);
    mpfr_mul_2ui(c, c, 700, MPFR_RNDN);
    mpz_t a, b, q0, q1, t, lim, qq;
    mpz_inits(a, b, q0, q1, t, lim, qq, (mpz_ptr)0);
    mpfr_get_z(a, c, MPFR_RNDZ);
    mpz_set_ui(b, 1); mpz_mul_2exp(b, b, 700);  // α = a/b
    mpz_set_ui(lim, 1); mpz_mul_2exp(lim, lim, 237);
    mpz_set_ui(q0, 0); mpz_set_ui(q1, 1);  // q_{-1}, q_0 for α < 1: CF [0; a1, a2, ...]
    // Euclid on b/a gives partial quotients of α = a/b.
    mpz_t num, den; mpz_inits(num, den, (mpz_ptr)0);
    mpz_set(num, b); mpz_set(den, a);
    mpz_set_ui(q0, 1); mpz_set_ui(q1, 0);  // q_{-1}=0? use standard: q_{-2}=1? keep q_prev=0, q_cur=1
    mpz_set_ui(q0, 0); mpz_set_ui(q1, 1);
    mpz_set(qq, q1);
    while (mpz_sgn(den) != 0) {
        mpz_fdiv_qr(t, num, num, den);       // t = partial quotient
        mpz_swap(num, den);
        mpz_mul(t, t, q1); mpz_add(t, t, q0);  // q_new = t q1 + q0
        if (mpz_cmp(t, lim) >= 0) break;
        mpz_set(q0, q1); mpz_set(q1, t);
        mpz_set(qq, q1);
    }
    // M = m·q in [2^236, 2^237)
    mpz_t M; mpz_init(M);
    mpz_sub_ui(t, lim, 1);
    mpz_fdiv_q(t, t, qq);
    mpz_mul(M, t, qq);
    mpfr_t xm; mpfr_init2(xm, 237);
    mpfr_set_z_2exp(xm, M, E, MPFR_RNDN);
    char* str = nullptr; mpfr_asprintf(&str, "%Ra", xm);
    const float256 x = float256::parse(str).value();
    mpfr_free_str(str);
    mpfr_clears(c, xm, (mpfr_ptr)0);
    mpz_clears(a, b, q0, q1, t, lim, qq, num, den, M, (mpz_ptr)0);
    return x;
}

int mpfr_lngamma_abs(mpfr_ptr r, mpfr_srcptr x, mpfr_rnd_t rnd) { int sg; return mpfr_lgamma(r, &sg, x, rnd); }

struct Fn {
    const char* name;
    std::function<float256(float256)> f;
    int (*m)(mpfr_ptr, mpfr_srcptr, mpfr_rnd_t);
    double lo, hi; bool logscale; bool neg_too;
    double bound;
    double maxabs;
    int gen = 0;  // 1: near multiples of pi/2, 2: 1 - 2^[lo,hi), 3: 1 + 2^[lo,hi), 4: negated, 5: 2 ± 2^[lo,hi), 6: hard_trig(E in [lo,hi))
};
} // namespace

int main(int argc, char** argv) {
    const long n = argc > 1 ? std::atol(argv[1]) : 2000;
    rng.seed(argc > 2 ? std::strtoull(argv[2], nullptr, 10) : 1);
    const bool check = argc > 3 && std::atoi(argv[3]) != 0;
    using namespace base256;
    std::vector<Fn> fns = {
        {"exp",   [](float256 x){return exp(x);},   mpfr_exp,   -120, 18, true, true, 1, 181000},
        {"exp2",  [](float256 x){return exp2(x);},  mpfr_exp2,  -120, 18, true, true, 1, 262000},
        {"expm1", [](float256 x){return expm1(x);}, mpfr_expm1, -300, 18, true, true, 1, 181000},
        {"log",   [](float256 x){return log(x);},   mpfr_log,   -262000, 262000, true, false, 1, 0},
        {"log2",  [](float256 x){return log2(x);},  mpfr_log2,  -262000, 262000, true, false, 1, 0},
        {"log10", [](float256 x){return log10(x);}, mpfr_log10, -262000, 262000, true, false, 1, 0},
        {"logn1", [](float256 x){return log(x);},   mpfr_log,   -2, 1, true, false, 1, 0},
        {"log1p", [](float256 x){return log1p(x);}, mpfr_log1p, -300, 300, true, false, 1, 0},
        {"log1pn",[](float256 x){return log1p(x);}, mpfr_log1p, -300, -1, true, false, 1, 0},
        {"sin",   [](float256 x){return sin(x);},   mpfr_sin,   -120, 100, true, true, 1, 0},
        {"cos",   [](float256 x){return cos(x);},   mpfr_cos,   -120, 100, true, true, 1, 0},
        {"tan",   [](float256 x){return tan(x);},   mpfr_tan,   -120, 100, true, true, 1, 0},
        {"sin_kpi",[](float256 x){return sin(x);},  mpfr_sin,   0, 99, false, true, 1, 0, 1},
        {"cos_kpi",[](float256 x){return cos(x);},  mpfr_cos,   0, 99, false, true, 1, 0, 1},
        {"tan_kpi",[](float256 x){return tan(x);},  mpfr_tan,   0, 99, false, true, 1, 0, 1},
        {"asin1", [](float256 x){return asin(x);},  mpfr_asin,  -236, -1, true, true, 1, 0, 2},
        {"acos1", [](float256 x){return acos(x);},  mpfr_acos,  -236, -1, true, true, 1, 0, 2},
        {"atanh1",[](float256 x){return atanh(x);}, mpfr_atanh, -236, -1, true, true, 1, 0, 2},
        {"acosh1",[](float256 x){return acosh(x);}, mpfr_acosh, -236, -1, true, false, 1, 0, 3},
        {"sin_all",[](float256 x){return sin(x);},  mpfr_sin,   -120, 262143, true, true, 1, 0},
        {"cos_all",[](float256 x){return cos(x);},  mpfr_cos,   -120, 262143, true, true, 1, 0},
        {"tan_all",[](float256 x){return tan(x);},  mpfr_tan,   -120, 262143, true, true, 1, 0},
        {"sin_kbig",[](float256 x){return sin(x);}, mpfr_sin,   100, 262140, false, true, 1, 0, 1},
        {"sin_hard",[](float256 x){return sin(x);}, mpfr_sin,   -236, 261907, false, true, 1, 0, 6},
        {"cos_hard",[](float256 x){return cos(x);}, mpfr_cos,   -236, 261907, false, true, 1, 0, 6},
        {"tan_hard",[](float256 x){return tan(x);}, mpfr_tan,   -236, 261907, false, true, 1, 0, 6},
        {"sin_hardcw",[](float256 x){return sin(x);}, mpfr_sin, -236, -137, false, true, 1, 0, 6},
        {"atan",  [](float256 x){return atan(x);},  mpfr_atan,  -120, 300, true, true, 1, 0},
        {"asin",  [](float256 x){return asin(x);},  mpfr_asin,  -120, 0, true, true, 1, 1},
        {"acos",  [](float256 x){return acos(x);},  mpfr_acos,  -120, 0, true, true, 1, 1},
        {"sinh",  [](float256 x){return sinh(x);},  mpfr_sinh,  -120, 18, true, true, 1, 181000},
        {"cosh",  [](float256 x){return cosh(x);},  mpfr_cosh,  -120, 18, true, true, 1, 181000},
        {"tanh",  [](float256 x){return tanh(x);},  mpfr_tanh,  -120, 8, true, true, 1, 0},
        {"asinh", [](float256 x){return asinh(x);}, mpfr_asinh, -120, 3000, true, true, 1, 0},
        {"acosh", [](float256 x){return acosh(x);}, mpfr_acosh, 0, 3000, true, false, 1, 0},
        {"atanh", [](float256 x){return atanh(x);}, mpfr_atanh, -120, 0, true, true, 1, 1},
        {"cbrt",  [](float256 x){return cbrt(x);},  mpfr_cbrt,  -262000, 262000, true, true, 1, 0},
        {"erf",   [](float256 x){return erf(x);},   mpfr_erf,   -120, 5, true, true, 1, 0},
        {"erfc",  [](float256 x){return erfc(x);},  mpfr_erfc,  -120, 9, true, false, 1, 500},
        {"erfcn", [](float256 x){return erfc(x);},  mpfr_erfc,  -120, 5, true, false, 1, 0, 4},
        {"tgamma",[](float256 x){return tgamma(x);},mpfr_gamma, -10, 10, true, false, 1, 1000},
        {"tgamman",[](float256 x){return tgamma(x);},mpfr_gamma, -10, 10, true, false, 1, 1000, 4},
        {"lgamman",[](float256 x){return lgamma(x);},mpfr_lngamma_abs, -10, 20, true, false, 1, 0, 4},
        {"lgamma1",[](float256 x){return lgamma(x);},mpfr_lngamma, -236, -1, true, false, 1, 0, 3},
        {"lgamma2",[](float256 x){return lgamma(x);},mpfr_lngamma_abs, -236, -1, true, true, 1, 0, 5},
        {"lgamma",[](float256 x){return lgamma(x);},mpfr_lngamma,-10, 20, true, false, 1, 0},
    };    if (argc > 4) { std::string only = argv[4]; std::erase_if(fns, [&](const Fn& f){ return only.find(f.name) == std::string::npos; }); }
    mpfr_t xm, r;
    mpfr_inits2(1100, xm, r, (mpfr_ptr)0);
    mpfr_set_emin(mpfr_get_emin_min()); mpfr_set_emax(mpfr_get_emax_max());
    int bad = 0;
    for (const Fn& f : fns) {
        double worst = 0; float256 wx;
        for (long i = 0; i < n; ++i) {
            float256 x = rnd_in(f.lo, f.hi, f.logscale);
            if (f.gen == 4) x = -x;
            if (f.gen == 5) x = float256::two() + x;
            if (f.gen == 3) x = float256::one() + rnd_in(f.lo, f.hi, true);
            if (f.gen == 2) x = copysign(float256::one() - rnd_in(f.lo, f.hi, true), x);
            if (f.gen == 6) x = hard_trig(static_cast<long>(f.lo + (f.hi - f.lo) * std::uniform_real_distribution<double>(0, 1)(rng)));
            if (f.gen == 1) {
                // the float256 nearest to k*pi/2 for k up to 2^(lo..hi)
                const int bits = static_cast<int>(f.lo + (f.hi - f.lo) * std::uniform_real_distribution<double>(0, 1)(rng));
                float256 k = trunc(ldexp(rnd_in(0, 1, false), bits)) + float256::one();
                mpfr_t t; mpfr_init2(t, bits + 1400);
                mpfr_const_pi(t, MPFR_RNDN); mpfr_div_2ui(t, t, 1, MPFR_RNDN);
                mpfr_t km; mpfr_init2(km, 300 + 0 * bits); to_mpfr(km, k);
                mpfr_mul(t, t, km, MPFR_RNDN);
                mpfr_t q; mpfr_init2(q, 237); mpfr_set(q, t, MPFR_RNDN);
                // read back via hex string
                char* str = nullptr; mpfr_asprintf(&str, "%Ra", q);
                x = float256::parse(str).value();
                mpfr_free_str(str); mpfr_clears(t, km, q, (mpfr_ptr)0);
            }
            if (f.neg_too && (rng() & 1)) x = -x;
            if (f.maxabs > 0 && abs(x) > float256(f.maxabs)) { --i; continue; }
            to_mpfr(xm, x);
            f.m(r, xm, MPFR_RNDN);
            const double e = ulp_err(f.f(x), r);
            if (e > worst) { worst = e; wx = x; }
        }
        const bool over = worst > f.bound;
        bad += over;
        std::printf("%-7s max %.3g ulp at %s%s\n", f.name, worst, wx.to_string().c_str(), over ? "  OVER BOUND" : "");
    }
    if (argc <= 4 || std::string(argv[4]).find("pow") != std::string::npos) {
        // pow(x, y): x log-uniform in [2^-64, 2^64), y in [-2000, 2000] (half integers, half not).
        mpfr_t ym;
        mpfr_init2(ym, 1100);
        double worst = 0; float256 wx, wy;
        for (long i = 0; i < n; ++i) {
            const float256 x = rnd_in(-64, 64, true);
            float256 y = rnd_in(-2000, 2000, false);
            if (rng() & 1) y = trunc(y);
            if (abs(y * log2(x)) > float256(250000)) { --i; continue; }
            to_mpfr(xm, x); to_mpfr(ym, y);
            mpfr_pow(r, xm, ym, MPFR_RNDN);
            const double e = ulp_err(pow(x, y), r);
            if (e > worst) { worst = e; wx = x; wy = y; }
        }
        const bool over = worst > 1.0;
        bad += over;
        std::printf("%-7s max %.3g ulp at %s, %s%s\n", "pow", worst, wx.to_string(20).c_str(), wy.to_string(20).c_str(), over ? "  OVER BOUND" : "");
        mpfr_clear(ym);
    }
    if (argc <= 4 || std::string(argv[4]).find("atan2") != std::string::npos) {
        mpfr_t ym;
        mpfr_init2(ym, 1100);
        double worst = 0; float256 wx, wy;
        for (long i = 0; i < n; ++i) {
            float256 y = rnd_in(-300, 300, true), x = rnd_in(-300, 300, true);
            if (rng() & 1) y = -y;
            if (rng() & 1) x = -x;
            to_mpfr(xm, x); to_mpfr(ym, y);
            mpfr_atan2(r, ym, xm, MPFR_RNDN);
            const double e = ulp_err(atan2(y, x), r);
            if (e > worst) { worst = e; wx = x; wy = y; }
        }
        const bool over = worst > 1.0;
        bad += over;
        std::printf("%-7s max %.3g ulp at %s, %s%s\n", "atan2", worst, wy.to_string(20).c_str(), wx.to_string().c_str(), over ? "  OVER BOUND" : "");
        mpfr_clear(ym);
    }
    mpfr_clears(xm, r, (mpfr_ptr)0);
    return check && bad ? 1 : 0;
}
