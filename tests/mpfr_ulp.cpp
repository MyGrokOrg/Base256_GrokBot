// Measures elementary-function error in ulps against MPFR (exact reference
// at 1100 bits). Usage: mpfr_ulp [cases_per_fn] [seed] [bound_check] ["fn list"] [modes 1..5]
// Every case is evaluated in each rounding mode (RNE, RNA, RZ, RU, RD) and
// checked for value (vs. MPFR rounded in that mode), direction and flags.
// With bound_check=1 the program exits nonzero if any function exceeds its
// documented bound (see README "Elementary function accuracy").
#include <base256/base256.hpp>
#include <mpfr.h>
#include <cmath>
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

// --- rounding modes and flags ------------------------------------------------
constexpr base256::round_mode kModes[5] = {base256::round_mode::to_nearest_even, base256::round_mode::to_nearest_away,
                                           base256::round_mode::toward_zero, base256::round_mode::upward,
                                           base256::round_mode::downward};
// MPFR_RNDNA (ties away) is emulated in round_b256.
constexpr mpfr_rnd_t kMpfrModes[5] = {MPFR_RNDN, MPFR_RNDNA, MPFR_RNDZ, MPFR_RNDU, MPFR_RNDD};
const char* const kModeNames[5] = {"RNE", "RNA", "RZ", "RU", "RD"};
int g_modes = 5;  // how many of kModes to test

struct Expected {
    bool inexact, overflow, tiny;
};
// Round the (1100-bit) reference to binary256 in rnd, emulating its exponent range.
Expected round_b256(mpfr_t out, mpfr_srcptr ref, mpfr_rnd_t rnd) {
    if (rnd == MPFR_RNDNA) {
        // MPFR's set() has no ties-away mode: round away when ref is exactly
        // halfway between its toward-zero and away-from-zero roundings.
        mpfr_t z, s;
        mpfr_init2(z, 237);
        mpfr_init2(s, 2400);
        round_b256(z, ref, MPFR_RNDZ);
        Expected e = round_b256(out, ref, MPFR_RNDA);
        bool tie = false;
        if (mpfr_number_p(z) && mpfr_number_p(out) && !mpfr_equal_p(z, out)) {
            mpfr_add(s, z, out, MPFR_RNDN);
            mpfr_div_2ui(s, s, 1, MPFR_RNDN);
            tie = mpfr_equal_p(s, ref) != 0;
        }
        if (!tie) e = round_b256(out, ref, MPFR_RNDN);
        mpfr_clears(z, s, (mpfr_ptr)0);
        return e;
    }
    const mpfr_exp_t emin = mpfr_get_emin(), emax = mpfr_get_emax();
    mpfr_set_emin(kEminN - 237 + 2);  // -262377: smallest subnormal 2^-262378
    mpfr_set_emax(262144);
    mpfr_clear_flags();
    int t = mpfr_set(out, ref, rnd);
    t = mpfr_check_range(out, t, rnd);
    t = mpfr_subnormalize(out, t, rnd);
    Expected e{};
    e.inexact = t != 0;
    e.overflow = mpfr_overflow_p() != 0;
    mpfr_set_emin(emin);
    mpfr_set_emax(emax);
    // Tininess after rounding (unbounded exponent): |round_237(ref)| < 2^-262142.
    mpfr_t u;
    mpfr_init2(u, 237);
    mpfr_set(u, ref, rnd);
    e.tiny = !mpfr_zero_p(u) && !mpfr_inf_p(u) && !mpfr_nan_p(u) && mpfr_get_exp(u) <= kEminN;
    mpfr_clear(u);
    return e;
}

// If the 1100-bit reference was rounded (ternary != 0), move it one 1100-bit ulp
// toward the true value, so that a true value just off a binary256 number (e.g.
// expm1(-80000) = -1 + 2^-115000) is still seen as inexact and on the right side.
void nudge(mpfr_t r, int ternary) {
    if (ternary > 0) mpfr_nextbelow(r);
    else if (ternary < 0) mpfr_nextabove(r);
}

struct Stat {
    double worst[5] = {0, 0, 0, 0, 0};
    long misround[5] = {0, 0, 0, 0, 0};
    long flag_bad = 0, dir_bad = 0, cases = 0;
    std::string first_bad;
    bool failed(double bound) const {
        for (int m = 0; m < 5; ++m) if (worst[m] > (m < 2 ? bound : 1.0)) return true;  // directed: at most 1 ulp (=1 only at a binade edge)
        return flag_bad > 0 || dir_bad > 0;
    }
};

// Evaluate `call` in every mode and check value, direction and flags against `ref`.
void check_case(Stat& st, const std::function<float256()>& call, mpfr_srcptr ref, const std::string& where) {
    mpfr_t o, ym;
    mpfr_init2(o, 237);
    mpfr_init2(ym, 1100);
    ++st.cases;
    for (int m = 0; m < g_modes; ++m) {
        base256::fenv::set_round(kModes[m]);
        base256::fenv::clear_flags();
        const float256 y = call();
        const unsigned fl = base256::fenv::get_flags();
        base256::fenv::set_round(base256::round_mode::to_nearest_even);
        base256::fenv::clear_flags();
        const Expected ex = round_b256(o, ref, kMpfrModes[m]);
        to_mpfr(ym, y);
        const bool same = (mpfr_nan_p(ym) && mpfr_nan_p(o)) || (mpfr_equal_p(ym, o) && mpfr_signbit(ym) == mpfr_signbit(o));
        if (!same) ++st.misround[m];
        double err = 0;
        // (an overflow saturated to the largest finite value in a directed mode is not an ulp error)
        if (!ex.overflow && !mpfr_inf_p(ym) && !mpfr_nan_p(ym) && !mpfr_inf_p(o)) err = ulp_err(y, const_cast<mpfr_ptr>(ref));
        if (!same) {
            if (mpfr_inf_p(ym) || mpfr_inf_p(o) || mpfr_nan_p(ym) || mpfr_nan_p(o)) err = 1e9;
            // direction
            const int c = mpfr_cmp(ym, ref);
            bool bad_dir = false;
            if (kMpfrModes[m] == MPFR_RNDU) bad_dir = c < 0;
            if (kMpfrModes[m] == MPFR_RNDD) bad_dir = c > 0;
            if (kMpfrModes[m] == MPFR_RNDZ) bad_dir = mpfr_cmpabs(ym, ref) > 0;
            if (bad_dir) {
                ++st.dir_bad;
                if (st.first_bad.empty()) st.first_bad = std::string("direction ") + kModeNames[m] + " at " + where;
            }
        }
        if (err > st.worst[m]) st.worst[m] = err;
        unsigned want = 0;
        if (ex.inexact) want |= base256::fe_inexact;
        if (ex.overflow) want |= base256::fe_overflow;
        if (ex.tiny && ex.inexact) want |= base256::fe_underflow;
        if (fl != want) {
            ++st.flag_bad;
            if (st.first_bad.empty()) {
                char buf[96];
                std::snprintf(buf, sizeof buf, "flags %s got %#x want %#x at ", kModeNames[m], fl, want);
                st.first_bad = buf + where;
            }
        }
    }
    mpfr_clears(o, ym, (mpfr_ptr)0);
}

void report(const char* name, const Stat& st, double bound, int& bad) {
    const bool over = st.failed(bound);
    bad += over;
    std::printf("%-10s n=%ld max ulp", name, st.cases);
    for (int m = 0; m < g_modes; ++m) std::printf(" %s %.3g", kModeNames[m], st.worst[m]);
    long mis = 0;
    for (int m = 0; m < g_modes; ++m) mis += st.misround[m];
    std::printf(" | misrounded %ld, flag errors %ld, direction errors %ld%s%s\n", mis, st.flag_bad, st.dir_bad,
                over ? "  OVER BOUND " : "", st.first_bad.empty() ? "" : (" [" + st.first_bad + "]").c_str());
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
    if (argc > 5) g_modes = std::atoi(argv[5]) < 1 ? 1 : std::atoi(argv[5]) > 5 ? 5 : std::atoi(argv[5]);
    using namespace base256;
    std::vector<Fn> fns = {
        {"exp",   [](float256 x){return exp(x);},   mpfr_exp,   -120, 18, true, true, 1, 181000},
        {"exp2",  [](float256 x){return exp2(x);},  mpfr_exp2,  -120, 18, true, true, 1, 262000},
        {"expm1", [](float256 x){return expm1(x);}, mpfr_expm1, -300, 18, true, true, 1, 181000},
        // Overflow / gradual-underflow / underflow-to-zero regions (flags and directed modes).
        {"exp_big", [](float256 x){return exp(x);},   mpfr_exp,   17, 18, true, true, 1, 0},
        {"exp2_big",[](float256 x){return exp2(x);},  mpfr_exp2,  17, 19, true, true, 1, 0},
        {"expm1_big",[](float256 x){return expm1(x);},mpfr_expm1, 17, 18, true, true, 1, 0},
        {"sinh_big",[](float256 x){return sinh(x);},  mpfr_sinh,  17, 18, true, true, 1, 0},
        {"cosh_big",[](float256 x){return cosh(x);},  mpfr_cosh,  17, 18, true, true, 1, 0},
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
        {"erf_lin", [](float256 x){return erf(x);},  mpfr_erf,  -17, 17, false, false, 1, 0},
        {"erfc_lin",[](float256 x){return erfc(x);}, mpfr_erfc, -17, 30, false, false, 1, 0},
        {"erfcn", [](float256 x){return erfc(x);},  mpfr_erfc,  -120, 5, true, false, 1, 0, 4},
        {"tgamma",[](float256 x){return tgamma(x);},mpfr_gamma, -10, 10, true, false, 1, 1000},
        {"tgamman",[](float256 x){return tgamma(x);},mpfr_gamma, -10, 10, true, false, 1, 1000, 4},
        {"lgamman",[](float256 x){return lgamma(x);},mpfr_lngamma_abs, -10, 20, true, false, 1, 0, 4},
        {"lgamma1",[](float256 x){return lgamma(x);},mpfr_lngamma, -236, -1, true, false, 1, 0, 3},
        {"zeta",   [](float256 x){return riemann_zeta(x);}, mpfr_zeta, -120, 8, true, true, 1, 0},
        {"zeta1",  [](float256 x){return riemann_zeta(x);}, mpfr_zeta, -236, -1, true, false, 1, 0, 3},
        {"zetab1", [](float256 x){return riemann_zeta(x);}, mpfr_zeta, -236, -1, true, false, 1, 0, 2},
        {"zetan",  [](float256 x){return riemann_zeta(x);}, mpfr_zeta, 8, 14, true, false, 1, 0, 4},
        {"lgamma2",[](float256 x){return lgamma(x);},mpfr_lngamma_abs, -236, -1, true, true, 1, 0, 5},
        {"lgamma",[](float256 x){return lgamma(x);},mpfr_lngamma,-10, 20, true, false, 1, 0},
    };
    const bool all_fns = argc <= 4 || std::string(argv[4]).empty() || std::string(argv[4]) == "all";
    if (!all_fns) { std::string only = argv[4]; std::erase_if(fns, [&](const Fn& f){ return only.find(f.name) == std::string::npos; }); }
    mpfr_t xm, r;
    mpfr_inits2(1100, xm, r, (mpfr_ptr)0);
    mpfr_set_emin(mpfr_get_emin_min()); mpfr_set_emax(mpfr_get_emax_max());
    int bad = 0;
    for (const Fn& f : fns) {
        Stat st;
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
            nudge(r, f.m(r, xm, MPFR_RNDN));
            check_case(st, [&] { return f.f(x); }, r, x.to_hex_string());
        }
        report(f.name, st, f.bound, bad);
    }
    if (all_fns || std::string(argv[4]).find("pow") != std::string::npos) {
        // pow(x, y): x log-uniform in [2^-64, 2^64), y in [-2000, 2000] (half integers, half not).
        mpfr_t ym;
        mpfr_init2(ym, 1100);
        Stat st;
        for (long i = 0; i < n; ++i) {
            const float256 x = rnd_in(-64, 64, true);
            float256 y = rnd_in(-2000, 2000, false);
            if (rng() & 1) y = trunc(y);
            if (abs(y * log2(x)) > float256(250000)) { --i; continue; }
            to_mpfr(xm, x); to_mpfr(ym, y);
            nudge(r, mpfr_pow(r, xm, ym, MPFR_RNDN));
            check_case(st, [&] { return pow(x, y); }, r, x.to_hex_string() + ", " + y.to_hex_string());
        }
        report("pow", st, 1.0, bad);
        mpfr_clear(ym);
    }
    if (all_fns || std::string(argv[4]).find("atan2") != std::string::npos) {
        mpfr_t ym;
        mpfr_init2(ym, 1100);
        Stat st;
        for (long i = 0; i < n; ++i) {
            float256 y = rnd_in(-300, 300, true), x = rnd_in(-300, 300, true);
            if (rng() & 1) y = -y;
            if (rng() & 1) x = -x;
            to_mpfr(xm, x); to_mpfr(ym, y);
            nudge(r, mpfr_atan2(r, ym, xm, MPFR_RNDN));
            check_case(st, [&] { return atan2(y, x); }, r, y.to_hex_string() + ", " + x.to_hex_string());
        }
        report("atan2", st, 1.0, bad);
        mpfr_clear(ym);
    }
    // Two-argument special functions: gen() draws (a, b); ref(r, a, b) is the MPFR reference.
    auto two_arg = [&](const char* name, const std::function<void(float256&, float256&)>& gen,
                       const std::function<float256(float256, float256)>& ours,
                       const std::function<int(mpfr_ptr, mpfr_srcptr, mpfr_srcptr)>& ref) {
        if (!all_fns && std::string(argv[4]).find(name) == std::string::npos) return;
        mpfr_t am, bm;
        mpfr_inits2(1100, am, bm, (mpfr_ptr)0);
        Stat st;
        for (long i = 0; i < n; ++i) {
            float256 a, b;
            gen(a, b);
            to_mpfr(am, a); to_mpfr(bm, b);
            nudge(r, ref(r, am, bm));
            check_case(st, [&] { return ours(a, b); }, r, a.to_hex_string() + ", " + b.to_hex_string());
        }
        report(name, st, 1.0, bad);
        mpfr_clears(am, bm, (mpfr_ptr)0);
    };
    auto rnd_int = [&](int lo, int hi) { return lo + static_cast<int>(rng() % static_cast<unsigned>(hi - lo + 1)); };
    // jn / yn (MPFR has integer orders): x log-uniform, jn also for negative n and x.
    two_arg("jn", [&](float256& a, float256& b) {
                a = float256::from_int(rnd_int(-40, 40));
                b = rnd_in(-12, 9, true);
                if (rng() & 1) b = -b;
            },
            [](float256 a, float256 b) { return base256::jn(static_cast<int>(a.to_double()), b); },
            [](mpfr_ptr r, mpfr_srcptr a, mpfr_srcptr b) { return mpfr_jn(r, mpfr_get_si(a, MPFR_RNDN), b, MPFR_RNDN); });
    two_arg("yn", [&](float256& a, float256& b) {
                a = float256::from_int(rnd_int(-40, 40));
                b = rnd_in(-12, 9, true);
            },
            [](float256 a, float256 b) { return base256::yn(static_cast<int>(a.to_double()), b); },
            [](mpfr_ptr r, mpfr_srcptr a, mpfr_srcptr b) { return mpfr_yn(r, mpfr_get_si(a, MPFR_RNDN), b, MPFR_RNDN); });
    // cyl_bessel_j / cyl_neumann at integer orders (the general-ν code paths).
    two_arg("cyl_j", [&](float256& a, float256& b) { a = float256::from_int(rnd_int(0, 30)); b = rnd_in(-4, 11, true); },
            [](float256 a, float256 b) { return base256::cyl_bessel_j(a, b); },
            [](mpfr_ptr r, mpfr_srcptr a, mpfr_srcptr b) { return mpfr_jn(r, mpfr_get_si(a, MPFR_RNDN), b, MPFR_RNDN); });
    two_arg("cyl_y", [&](float256& a, float256& b) { a = float256::from_int(rnd_int(0, 30)); b = rnd_in(-4, 11, true); },
            [](float256 a, float256 b) { return base256::cyl_neumann(a, b); },
            [](mpfr_ptr r, mpfr_srcptr a, mpfr_srcptr b) { return mpfr_yn(r, mpfr_get_si(a, MPFR_RNDN), b, MPFR_RNDN); });
    // beta(x, y), x and y log-uniform in [2^-30, 2^30).
    two_arg("beta", [&](float256& a, float256& b) { a = rnd_in(-30, 30, true); b = rnd_in(-30, 30, true); },
            [](float256 a, float256 b) { return base256::beta(a, b); },
            [](mpfr_ptr r, mpfr_srcptr a, mpfr_srcptr b) { return mpfr_beta(r, a, b, MPFR_RNDN); });
    // Incomplete gamma: a in [2^-8, 2^8), x = a·2^[-4, 4) or log-uniform; references
    // from mpfr_gamma_inc (upper) and Γ(a) at 3000 bits for the others.
    auto ig_gen = [&](float256& a, float256& b) {
        a = rnd_in(-8, 8, true);
        b = (rng() & 1) ? a * rnd_in(-4, 4, true) : rnd_in(-12, 9, true);
    };
    auto ig_ref = [](int kind) {
        return [kind](mpfr_ptr r, mpfr_srcptr a, mpfr_srcptr x) {
            // Γ(a) − Γ(a,x) cancels about log2(Γ(a)/γ(a,x)) bits; γ ≈ x^a e^{-x}/a for small x.
            const double ad = mpfr_get_d(a, MPFR_RNDN), xd = mpfr_get_d(x, MPFR_RNDN);
            const double cancel = (std::lgamma(ad) - ad * std::log(xd) + xd + std::log(ad)) / std::log(2.0);
            const long prec = 3000 + (cancel > 0 ? static_cast<long>(cancel) : 0);
            mpfr_t u, g;
            mpfr_inits2(prec, u, g, (mpfr_ptr)0);
            mpfr_gamma_inc(u, a, x, MPFR_RNDN);
            mpfr_gamma(g, a, MPFR_RNDN);
            if (kind == 0 || kind == 2) mpfr_sub(u, g, u, MPFR_RNDN);
            if (kind >= 2) mpfr_div(u, u, g, MPFR_RNDN);
            const int t = mpfr_set(r, u, MPFR_RNDN);
            mpfr_clears(u, g, (mpfr_ptr)0);
            return t;
        };
    };
    two_arg("tgamma_lower", ig_gen, [](float256 a, float256 b) { return base256::tgamma_lower(a, b); }, ig_ref(0));
    two_arg("tgamma_upper", ig_gen, [](float256 a, float256 b) { return base256::tgamma_upper(a, b); }, ig_ref(1));
    two_arg("gamma_p", ig_gen, [](float256 a, float256 b) { return base256::gamma_p(a, b); }, ig_ref(2));
    two_arg("gamma_q", ig_gen, [](float256 a, float256 b) { return base256::gamma_q(a, b); }, ig_ref(3));
    // Orthogonal polynomials of integer degree n in [0, 40]: exact references from
    // the integer-coefficient forms at a precision that holds every term
    // exactly (one final rounding), so exact results and flags are checked too.
    // x is log-uniform in [2^-8, 2^6) with either sign, or a multiple of 1/64.
    auto poly_gen = [&](float256& a, float256& b) {
        a = float256::from_int(rnd_int(0, 40));
        b = (rng() & 1) ? rnd_in(-8, 6, true) : ldexp(float256::from_int(rnd_int(-256, 256)), -6);
        if (rng() & 1) b = -b;
    };
    auto poly_prec = [](long n) { return 260 * (n + 1) + 64 * (n + 1) + 256; };
    // P_l(x) = 2^-l Σ_k (−1)^k C(l,k) C(2l−2k,l) x^{l−2k}
    auto leg_ref = [&](mpfr_ptr r, mpfr_srcptr a, mpfr_srcptr x) {
        const long l = mpfr_get_si(a, MPFR_RNDN);
        mpfr_t s, p; mpz_t c, c2;
        mpfr_inits2(poly_prec(l), s, p, (mpfr_ptr)0); mpz_inits(c, c2, nullptr);
        mpfr_set_zero(s, 1);
        for (long k = 0; 2 * k <= l; ++k) {
            mpz_bin_uiui(c, static_cast<unsigned long>(l), static_cast<unsigned long>(k));
            mpz_bin_uiui(c2, static_cast<unsigned long>(2 * l - 2 * k), static_cast<unsigned long>(l));
            mpz_mul(c, c, c2);
            if (k & 1) mpz_neg(c, c);
            mpfr_pow_ui(p, x, static_cast<unsigned long>(l - 2 * k), MPFR_RNDN);
            mpfr_mul_z(p, p, c, MPFR_RNDN);
            mpfr_add(s, s, p, MPFR_RNDN);
        }
        mpfr_div_2ui(s, s, static_cast<unsigned long>(l), MPFR_RNDN);
        if (mpfr_zero_p(s)) mpfr_setsign(s, s, (l & 1) && mpfr_signbit(x), MPFR_RNDN);  // odd: P_l(±0) = ±0
        const int tr = mpfr_set(r, s, MPFR_RNDN);
        mpfr_clears(s, p, (mpfr_ptr)0); mpz_clears(c, c2, nullptr);
        return tr;
    };
    two_arg("legendre", poly_gen, [](float256 a, float256 b) { return base256::legendre(static_cast<unsigned>(a.to_double()), b); }, leg_ref);
    // H_n by its recurrence (integer arithmetic on x, exact at this precision).
    two_arg("hermite", poly_gen, [](float256 a, float256 b) { return base256::hermite(static_cast<unsigned>(a.to_double()), b); },
            [&](mpfr_ptr r, mpfr_srcptr a, mpfr_srcptr x) {
                const long n = mpfr_get_si(a, MPFR_RNDN);
                mpfr_t p0, p1, p2, t;
                mpfr_inits2(poly_prec(n) + 8 * n, p0, p1, p2, t, (mpfr_ptr)0);
                mpfr_set_ui(p0, 1, MPFR_RNDN);
                mpfr_mul_2ui(p1, x, 1, MPFR_RNDN);
                if (n == 0) mpfr_set(p1, p0, MPFR_RNDN);
                for (long j = 1; j < n; ++j) {
                    mpfr_mul(p2, x, p1, MPFR_RNDN);
                    mpfr_mul_ui(t, p0, static_cast<unsigned long>(j), MPFR_RNDN);
                    mpfr_sub(p2, p2, t, MPFR_RNDN);
                    mpfr_mul_2ui(p2, p2, 1, MPFR_RNDN);
                    mpfr_set(p0, p1, MPFR_RNDN); mpfr_set(p1, p2, MPFR_RNDN);
                }
                if (mpfr_zero_p(p1)) mpfr_setsign(p1, p1, (n & 1) && mpfr_signbit(x), MPFR_RNDN);
                const int tr = mpfr_set(r, p1, MPFR_RNDN);
                mpfr_clears(p0, p1, p2, t, (mpfr_ptr)0);
                return tr;
            });
    // L_n^(m)(x) = (1/n!) Σ_k (−1)^k C(n+m, n−k) (n!/k!) x^k, m in [0, 10] packed as a = n + 64·m.
    auto lag_gen = [&](float256& a, float256& b) {
        poly_gen(a, b);
        a = a + float256::from_int(64 * ((rng() & 1) ? 0 : rnd_int(0, 10)));
    };
    two_arg("laguerre", lag_gen,
            [](float256 a, float256 b) {
                const unsigned v = static_cast<unsigned>(a.to_double());
                return v < 64 ? base256::laguerre(v, b) : base256::assoc_laguerre(v % 64, v / 64, b);
            },
            [&](mpfr_ptr r, mpfr_srcptr a, mpfr_srcptr x) {
                const long v = mpfr_get_si(a, MPFR_RNDN), n = v % 64, m = v / 64;
                mpfr_t s, p; mpz_t c, f;
                mpfr_inits2(poly_prec(n) + 16 * n, s, p, (mpfr_ptr)0); mpz_inits(c, f, nullptr);
                mpfr_set_zero(s, 1);
                for (long k = 0; k <= n; ++k) {
                    mpz_bin_uiui(c, static_cast<unsigned long>(n + m), static_cast<unsigned long>(n - k));
                    mpz_fac_ui(f, static_cast<unsigned long>(n));
                    mpz_mul(c, c, f);
                    mpz_fac_ui(f, static_cast<unsigned long>(k));
                    mpz_divexact(c, c, f);
                    if (k & 1) mpz_neg(c, c);
                    mpfr_pow_ui(p, x, static_cast<unsigned long>(k), MPFR_RNDN);
                    mpfr_mul_z(p, p, c, MPFR_RNDN);
                    mpfr_add(s, s, p, MPFR_RNDN);
                }
                mpz_fac_ui(f, static_cast<unsigned long>(n));
                const int tr = mpfr_div_z(r, s, f, MPFR_RNDN);
                mpfr_clears(s, p, (mpfr_ptr)0); mpz_clears(c, f, nullptr);
                return tr;
            });
    mpfr_clears(xm, r, (mpfr_ptr)0);
    return check && bad ? 1 : 0;
}
