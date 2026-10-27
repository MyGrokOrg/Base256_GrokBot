// Elementary functions: rounding modes and IEEE exception flags.
#include <base256/base256.hpp>
#include <functional>
#include <iostream>
#include <string>

using base256::float256;
using base256::round_mode;
namespace fenv = base256::fenv;
using F = std::function<float256(float256)>;

static int g_fails = 0, g_pass = 0;
#define CHECK(cond)                                                              \
    do {                                                                         \
        if (cond) ++g_pass;                                                      \
        else { ++g_fails; std::cerr << "FAIL " << __LINE__ << "  " #cond "\n"; } \
    } while (0)

constexpr unsigned kInv = base256::fe_invalid, kDiv = base256::fe_divbyzero, kOvf = base256::fe_overflow,
                   kUnf = base256::fe_underflow, kInx = base256::fe_inexact;
const round_mode kModes[5] = {round_mode::to_nearest_even, round_mode::to_nearest_away, round_mode::toward_zero,
                              round_mode::upward, round_mode::downward};

// Evaluate in every mode; the flags raised must be exactly `want` and the value `expect`
// (bitwise, so NaN compares by isnan and zeros by sign).
static bool same(float256 a, float256 b) {
    if (a.isnan() || b.isnan()) return a.isnan() && b.isnan();
    return a == b && a.signbit() == b.signbit();
}
static void expect(const char* what, int line, const std::function<float256()>& f, float256 value, unsigned want) {
    for (round_mode m : kModes) {
        fenv::scoped_round g(m);
        fenv::clear_flags();
        const float256 y = f();
        const unsigned got = fenv::test_flags();
        fenv::clear_flags();
        if (same(y, value) && got == want) { ++g_pass; continue; }
        ++g_fails;
        std::cerr << "FAIL " << line << " " << what << " mode " << int(m) << ": got " << y.to_string(10) << " flags "
                  << got << ", want " << value.to_string(10) << " flags " << want << "\n";
    }
}
#define EXPECT(expr, value, flags) expect(#expr, __LINE__, [&] { return (expr); }, (value), (flags))

// As expect, with one expected value per mode (RNE RNA RZ RU RD), given as hex literals.
static void expect5(const char* what, int line, const std::function<float256()>& f, const char* const (&values)[5],
                    unsigned want) {
    for (int i = 0; i < 5; ++i) {
        fenv::scoped_round g(kModes[i]);
        const float256 value = float256::parse(values[i]).value();
        fenv::clear_flags();
        const float256 y = f();
        const unsigned got = fenv::test_flags();
        fenv::clear_flags();
        if (same(y, value) && got == want) { ++g_pass; continue; }
        ++g_fails;
        std::cerr << "FAIL " << line << " " << what << " mode " << i << ": got " << y.to_hex_string() << " flags " << got
                  << ", want " << value.to_hex_string() << " flags " << want << "\n";
    }
}
#define EXPECT5(expr, values, flags) expect5(#expr, __LINE__, [&] { return (expr); }, (values), (flags))

int main() {
    using namespace base256;
    const float256 one = float256::one(), zero{}, nzero = -zero, two(2.0), half(0.5), inf = float256::inf(),
                   qnan = float256::qnan(), snan = float256::snan();

    // --- exact cases: no flags in any mode ---------------------------------------
    EXPECT(exp(zero), one, 0);
    EXPECT(exp(nzero), one, 0);
    EXPECT(exp2(float256(10)), float256(1024), 0);
    EXPECT(exp2(float256(-3)), float256(0.125), 0);
    EXPECT(expm1(zero), zero, 0);
    EXPECT(expm1(nzero), nzero, 0);
    EXPECT(log(one), zero, 0);
    EXPECT(log2(float256(1024)), float256(10), 0);
    EXPECT(log2(float256(0.25)), float256(-2), 0);
    EXPECT(log10(float256(1000)), float256(3), 0);
    EXPECT(log1p(zero), zero, 0);
    EXPECT(sin(zero), zero, 0);
    EXPECT(sin(nzero), nzero, 0);
    EXPECT(cos(zero), one, 0);
    EXPECT(tan(nzero), nzero, 0);
    EXPECT(asin(zero), zero, 0);
    EXPECT(acos(one), zero, 0);
    EXPECT(atan(nzero), nzero, 0);
    EXPECT(sinh(zero), zero, 0);
    EXPECT(cosh(zero), one, 0);
    EXPECT(tanh(nzero), nzero, 0);
    EXPECT(asinh(zero), zero, 0);
    EXPECT(acosh(one), zero, 0);
    EXPECT(atanh(zero), zero, 0);
    EXPECT(erf(zero), zero, 0);
    EXPECT(erfc(zero), one, 0);
    EXPECT(cbrt(float256(-27)), float256(-3), 0);
    EXPECT(cbrt(float256(0.125)), half, 0);
    EXPECT(tgamma(float256(5)), float256(24), 0);
    EXPECT(tgamma(one), one, 0);
    EXPECT(lgamma(one), zero, 0);
    EXPECT(lgamma(two), zero, 0);
    EXPECT(pow(two, float256(10)), float256(1024), 0);
    EXPECT(pow(float256(3), float256(4)), float256(81), 0);
    EXPECT(pow(float256(4), half), two, 0);
    EXPECT(pow(float256(8), -one), float256(0.125), 0);
    EXPECT(pow(float256(16), float256(0.75)), float256(8), 0);
    EXPECT(pow(float256(7), zero), one, 0);
    EXPECT(pow(one, qnan), one, 0);
    EXPECT(atan2(zero, one), zero, 0);
    EXPECT(exp(inf), inf, 0);
    EXPECT(exp(-inf), zero, 0);
    EXPECT(log(inf), inf, 0);
    EXPECT(atan(qnan), qnan, 0);

    // --- invalid: domain errors and signaling NaNs --------------------------------
    EXPECT(log(-one), qnan, kInv);
    EXPECT(log2(-inf), qnan, kInv);
    EXPECT(log1p(float256(-2)), qnan, kInv);
    EXPECT(sqrt(-one), qnan, kInv);
    EXPECT(sin(inf), qnan, kInv);
    EXPECT(cos(-inf), qnan, kInv);
    EXPECT(tan(inf), qnan, kInv);
    EXPECT(asin(two), qnan, kInv);
    EXPECT(acos(-two), qnan, kInv);
    EXPECT(acosh(half), qnan, kInv);
    EXPECT(atanh(two), qnan, kInv);
    EXPECT(pow(-two, half), qnan, kInv);
    EXPECT(tgamma(float256(-3)), qnan, kInv);
    EXPECT(tgamma(-inf), qnan, kInv);
    for (const F& f : {F([](float256 x) { return exp(x); }), F([](float256 x) { return log(x); }),
                       F([](float256 x) { return sin(x); }), F([](float256 x) { return atanh(x); }),
                       F([](float256 x) { return erfc(x); }), F([](float256 x) { return lgamma(x); }),
                       F([](float256 x) { return cbrt(x); }), F([](float256 x) { return tanh(x); })}) {
        EXPECT(f(snan), qnan, kInv);
        EXPECT(f(qnan), qnan, 0);
    }
    EXPECT(pow(snan, one), qnan, kInv);
    EXPECT(atan2(one, snan), qnan, kInv);

    // --- divbyzero: poles --------------------------------------------------------
    EXPECT(log(zero), -inf, kDiv);
    EXPECT(log10(nzero), -inf, kDiv);
    EXPECT(log1p(-one), -inf, kDiv);
    EXPECT(atanh(one), inf, kDiv);
    EXPECT(atanh(-one), -inf, kDiv);
    EXPECT(lgamma(zero), inf, kDiv);
    EXPECT(lgamma(float256(-4)), inf, kDiv);
    EXPECT(tgamma(nzero), -inf, kDiv);
    EXPECT(pow(zero, -one), inf, kDiv);
    EXPECT(pow(nzero, float256(-3)), -inf, kDiv);
    EXPECT(pow(nzero, float256(-2)), inf, kDiv);

    // --- overflow / underflow follow the mode -----------------------------------
    const float256 big(200000), max = float256::max(), dmin = float256::denorm_min();
    for (round_mode m : kModes) {
        fenv::scoped_round g(m);
        fenv::clear_flags();
        const float256 y = exp(big);
        CHECK(fenv::test_flags() == (kOvf | kInx));
        const bool to_max = m == round_mode::toward_zero || m == round_mode::downward;
        CHECK(same(y, to_max ? max : inf));
        fenv::clear_flags();
        const float256 z = exp(-big);
        CHECK(fenv::test_flags() == (kUnf | kInx));
        CHECK(same(z, m == round_mode::upward ? dmin : zero));
        fenv::clear_flags();
    }
    // Exact subnormal result: no underflow, no inexact.
    EXPECT(exp2(float256(-262300)), ldexp(one, -262300), 0);

    // --- directed modes bracket the true value ------------------------------------
    for (const F& f : {F([](float256 x) { return exp(x); }), F([](float256 x) { return log(x); }),
                       F([](float256 x) { return sin(x); }), F([](float256 x) { return cos(x); }),
                       F([](float256 x) { return atan(x); }), F([](float256 x) { return erf(x); }),
                       F([](float256 x) { return tgamma(x); }), F([](float256 x) { return cbrt(x); })}) {
        for (float256 x : {float256(3.3), float256(-0.7), float256(1e10)}) {
            float256 r[5];
            for (int i = 0; i < 5; ++i) {
                fenv::scoped_round g(kModes[i]);
                r[i] = f(x);
            }
            fenv::clear_flags();
            if (r[0].isnan()) continue;  // log(-0.7)
            CHECK(r[3] == nextup(r[4]));                                     // RU = RD + 1 ulp
            CHECK(r[0] == r[3] || r[0] == r[4]);                             // RNE is one of them
            CHECK(r[2] == (r[4].signbit() ? r[3] : r[4]));                   // RZ = toward zero
        }
    }


    // --- tiny arguments in the directed modes: f(t) = lead + a correction far below an
    // ulp, so only its sign decides (|t| < 2^-237: the correction is below the
    // double-float256 precision and must still set the direction) ---------------------
    {
        auto near_val = [](float256 v, int s, round_mode m) {
            if (m == round_mode::upward && s > 0) return nextup(v);
            if (m == round_mode::downward && s < 0) return nextdown(v);
            if (m == round_mode::toward_zero && ((s > 0) == v.signbit())) return nextafter(v, float256());
            return v;
        };
        auto tiny_case = [&](const char* name, int line, const std::function<float256()>& f, float256 v, int s) {
            for (round_mode m : kModes) {
                fenv::scoped_round g(m);
                fenv::clear_flags();
                const float256 y = f();
                const unsigned got = fenv::test_flags();
                fenv::clear_flags();
                if (same(y, near_val(v, s, m)) && got == kInx) { ++g_pass; continue; }
                ++g_fails;
                std::cerr << "FAIL " << line << " " << name << " mode " << int(m) << "\n";
            }
        };
        for (int e : {-240, -300, -2000, -100000}) {
            for (int sg : {1, -1}) {
                const float256 t = ldexp(float256(1.5) * float256(sg), e);
                tiny_case("sin", __LINE__, [&] { return sin(t); }, t, -sg);
                tiny_case("tan", __LINE__, [&] { return tan(t); }, t, sg);
                tiny_case("atan", __LINE__, [&] { return atan(t); }, t, -sg);
                tiny_case("asin", __LINE__, [&] { return asin(t); }, t, sg);
                tiny_case("asinh", __LINE__, [&] { return asinh(t); }, t, -sg);
                tiny_case("atanh", __LINE__, [&] { return atanh(t); }, t, sg);
                tiny_case("sinh", __LINE__, [&] { return sinh(t); }, t, sg);
                tiny_case("tanh", __LINE__, [&] { return tanh(t); }, t, -sg);
                tiny_case("expm1", __LINE__, [&] { return expm1(t); }, t, 1);
                tiny_case("log1p", __LINE__, [&] { return log1p(t); }, t, -1);
                tiny_case("atan2", __LINE__, [&] { return atan2(t, one); }, t, -sg);
                if (e < -240) {
                    tiny_case("exp", __LINE__, [&] { return exp(t); }, one, sg);
                    tiny_case("exp2", __LINE__, [&] { return exp2(t); }, one, sg);
                    tiny_case("cos", __LINE__, [&] { return cos(t); }, one, -1);
                    tiny_case("cosh", __LINE__, [&] { return cosh(t); }, one, 1);
                    tiny_case("erfc", __LINE__, [&] { return erfc(t); }, one, -sg);
                }
            }
        }
    }

    // --- no internal flags leak, caller's flags kept --------------------------------
    fenv::clear_flags();
    fenv::raise_flags(kDiv);
    (void)exp(zero);
    (void)sin(float256(1e100));
    (void)lgamma(float256(-2.5));
    (void)pow(float256(3), float256(2.5));
    CHECK(fenv::test_flags() == (kDiv | kInx));
    fenv::clear_flags();
    CHECK(fenv::get_round() == round_mode::to_nearest_even);

    // --- atan2 with x > 0 and y/x in or below the subnormal range (4.26.0 lost underflow/inexact
    // and returned 0 in every mode below the range, and rounded twice inside it) ------------
    {
        auto H = [](const char* t) { return float256::parse(t).value(); };
        const float256 p1 = H("0x1p-200000"), p2 = H("0x1p+100000");  // y/x = 2^-300000
        const char* const below_pos[5] = {"0x0p+0", "0x0p+0", "0x0p+0", "0x1p-262378", "0x0p+0"};
        const char* const below_neg[5] = {"-0x0p+0", "-0x0p+0", "-0x0p+0", "-0x0p+0", "-0x1p-262378"};
        EXPECT5(atan2(p1, p2), below_pos, kUnf | kInx);
        EXPECT5(atan2(-p1, p2), below_neg, kUnf | kInx);
        // The nightly-fuzz point (cx_check log_wide, seed 20261019): y/x about 2^-275000.
        const float256 fx = H("0x14187d14ec0fd805938558f1444a2f62aba612033aaa6f0c1b092035dff8p+120770"),
                       fy = H("0x1ad5b9f3f5eea603e991c4e874c1d6c91a7d4003416b3c4dee0388e6e160p-154322");
        EXPECT5(atan2(fy, fx), below_pos, kUnf | kInx);
        // y/x inside the subnormal range (references: MPFR atan2 at 3000 bits, rounded per mode).
        const char* const sub_pos[5] = {
            "0x112257693dc5eb2b0894053922a1496aa800000000000000000000000000p-262485",
            "0x112257693dc5eb2b0894053922a1496aa800000000000000000000000000p-262485",
            "0x112257693dc5eb2b0894053922a1496aa800000000000000000000000000p-262485",
            "0x112257693dc5eb2b0894053922a1496ab000000000000000000000000000p-262485",
            "0x112257693dc5eb2b0894053922a1496aa800000000000000000000000000p-262485"};
        EXPECT5(atan2(H("0x1ff202c7bff261b339ff248174e5598b88dbaa99e07987751d4ca8501e2cp-227121"),
                      H("0x1dd4b0caae1c75d0dd66cf72f858a4b66f8c462804db7b87a9e25fefe911p+35128")),
                sub_pos, kUnf | kInx);
        const char* const sub_neg[5] = {
            "-0x14d06aaba62f576720e5bf0a6f175288d13b029762bfda205d76f2000000p-262403",
            "-0x14d06aaba62f576720e5bf0a6f175288d13b029762bfda205d76f2000000p-262403",
            "-0x14d06aaba62f576720e5bf0a6f175288d13b029762bfda205d76f2000000p-262403",
            "-0x14d06aaba62f576720e5bf0a6f175288d13b029762bfda205d76f2000000p-262403",
            "-0x14d06aaba62f576720e5bf0a6f175288d13b029762bfda205d76f4000000p-262403"};
        EXPECT5(atan2(H("-0x1e3ecf44dd3f89e7d15f17362f25244caf9c4dabb4817253edc618187993p-237976"),
                      H("0x173ff6fa5db8656abd72fb710734986e86cb0ab8ab67a26b7f62b1852f27p+24191")),
                sub_neg, kUnf | kInx);
        const char* const sub_neg2[5] = {
            "-0x16b85642da51f49e7af3bd45e2ca512858a49e4eb6af230a9b6d12e08000p-262391",
            "-0x16b85642da51f49e7af3bd45e2ca512858a49e4eb6af230a9b6d12e08000p-262391",
            "-0x16b85642da51f49e7af3bd45e2ca512858a49e4eb6af230a9b6d12e08000p-262391",
            "-0x16b85642da51f49e7af3bd45e2ca512858a49e4eb6af230a9b6d12e08000p-262391",
            "-0x16b85642da51f49e7af3bd45e2ca512858a49e4eb6af230a9b6d12e0a000p-262391"};
        EXPECT5(atan2(H("-0x136a0023b682af5570eed8e94b150452ef05f542441d111b8aaa62f28d1ap-224249"),
                      H("0x1b57faf8cda9601e5b45785116080d650372e90794dfed52a24135b00a54p+37905")),
                sub_neg2, kUnf | kInx);
    }

    std::cout << "passed=" << g_pass << " failed=" << g_fails << "\n";
    return g_fails ? 1 : 0;
}
