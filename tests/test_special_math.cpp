// Special functions: reference values (mpmath, see tools/sf_refgen.py),
// special arguments, domain errors, IEEE flags and rounding modes.
#include <base256/base256.hpp>
#include <functional>
#include <iostream>
#include <string>

using base256::float256;
using base256::round_mode;
namespace fenv = base256::fenv;

static int g_fails = 0, g_pass = 0;
#define CHECK(cond)                                                              \
    do {                                                                         \
        if (cond) ++g_pass;                                                      \
        else { ++g_fails; std::cerr << "FAIL " << __LINE__ << "  " #cond "\n"; } \
    } while (0)

constexpr unsigned kInv = base256::fe_invalid, kDiv = base256::fe_divbyzero, kInx = base256::fe_inexact;
const round_mode kModes[5] = {round_mode::to_nearest_even, round_mode::to_nearest_away, round_mode::toward_zero,
                              round_mode::upward, round_mode::downward};

// Literals are parsed in round-to-nearest (parsing follows the current mode).
static float256 P(const std::string& s) {
    fenv::scoped_round g(round_mode::to_nearest_even);
    return float256::parse(s).value();
}
static bool same(float256 a, float256 b) {
    if (a.isnan() || b.isnan()) return a.isnan() && b.isnan();
    return a == b && a.signbit() == b.signbit();
}

static float256 eval(const std::string& op, const std::string& a, const std::string& b, const std::string& c) {
    using namespace base256;
    if (op == "j") return cyl_bessel_j(P(a), P(b));
    if (op == "y") return cyl_neumann(P(a), P(b));
    if (op == "i") return cyl_bessel_i(P(a), P(b));
    if (op == "k") return cyl_bessel_k(P(a), P(b));
    if (op == "sj") return sph_bessel(static_cast<unsigned>(std::stoul(a)), P(b));
    if (op == "sy") return sph_neumann(static_cast<unsigned>(std::stoul(a)), P(b));
    if (op == "jn") return jn(std::stoi(a), P(b));
    if (op == "yn") return yn(std::stoi(a), P(b));
    if (op == "beta") return beta(P(a), P(b));
    if (op == "ibeta") return ibeta(P(a), P(b), P(c));
    if (op == "ibetac") return ibetac(P(a), P(b), P(c));
    std::cerr << "unknown op " << op << "\n";
    return float256::qnan();
}

struct Ref { const char *op, *a, *b, *c, *want; };
const Ref kRefs[] = {
    {"j", "0", "1", "", "0x187c7fdbd7b8ef839e3530c22b34cedc3cbaf7ccbfde50538b5581f88245p-237"},
    {"j", "0.5", "3", "", "0x10a460536e4369cdeea10a105e02c4d339d175bb19b68ab3bad4748baf03p-240"},
    {"j", "2.25", "17.5", "", "0x11e90ccdd217718c7f58e6dc5dcebfdb9d0294d3d9976efae66fade762bfp-239"},
    {"j", "10", "1000", "", "-0x191bef1b8729ba9c86f779422721b421d95e3a6403194a55c53953cf2e35p-242"},
    {"j", "150", "100", "", "0x139ede31fcbf6c53b8b620481a3dce2a3666524dd1ceef28d9ae84b0a7d5p-288"},
    {"y", "0", "1", "", "0x16980226f358df29efcc5d79265550eb6a378dcf2886fc1c0a09de8b0620p-240"},
    {"y", "1.75", "0.125", "", "-0x12d27997639523a381cbea362edb935e6ec16ed7340a3198ef7752f0c014p-231"},
    {"y", "3", "60", "", "-0x18464d396d8dd927d80a1f109da42be740af51fb4587a56bd01a4e9ad035p-240"},
    {"y", "0.5", "500", "", "0x1025bd880e0ed349a5ef870416c3903d80906fd332b6b8c5831a98dd6c53p-241"},
    {"i", "0", "1", "", "0x1441ce4b386c2caea109def999685cf4119e467a3fc15a41f6af801536f7p-236"},
    {"i", "2.5", "30", "", "0x1476aea133e68109e4e300d13ba2df7a602516abdda2c78cf1d829447fa8p-197"},
    {"i", "1", "1000", "", "0x14e7fd269bb9b3bdf4586f634914be05b7b5fadaaaa789ac46abcf599d90p+1200"},
    {"k", "0", "1", "", "0x1af2107c43e119d373d510b59feb07626add9584e1b15087a4be61723466p-238"},
    {"k", "0.25", "0.0625", "", "0x1a49a6f8c89a47c067e1da8295e0cb26a06b5d9e81132a11590d599a6b44p-235"},
    {"k", "7", "70", "", "0x1b4d1167a72bfd23711b3d2e3313efbbae0d85314dc2a2c12a9c459bb130p-340"},
    {"k", "1", "1000", "", "0x1913fdf0aa7e0c1fa817ba3ae25f882fe74225a321d00661627c54c55fcfp-1684"},
    {"sj", "3", "2", "", "0x1f16f77f368fbe312b6268f8d9581a93c7f689ff039ef7ffa9ac70ae35b5p-241"},
    {"sy", "2", "1.5", "", "-0x15880a088e972991acc4f3441c8ddb3839cb1bc4a72d1fa60afa386c7fa4p-236"},
    {"jn", "-3", "-2.5", "", "0x1bb98fc5e82abb20d377d0fadde40255df35beecc62ddbed9ed35ab495d0p-239"},
    {"yn", "-3", "2.5", "", "0x18319b4c1093dfe62d671ffdb32be6e7c68c6efecb13bda0c2cd721f8a40p-237"},
    {"beta", "0.5", "0.5", "", "0x1921fb54442d18469898cc51701b839a252049c1114cf98e804177d4c762p-235"},
    {"beta", "3.25", "1e-5", "", "0x1869e67e354ec82c32f8c9d4e6353bad1d5bf1e7dd44298bf37a6edafb34p-220"},
    {"beta", "2000.5", "3000.25", "", "0x1879218301e67217c8570d69842cf4e7b940ac91d918a1dfac085628536bp-5096"},
    {"ibeta", "2.5", "3.5", "0.25", "0x1ac9d438759618f7f2aa3fcae85c861403ee6aefeafcd97b6bcfc6daa460p-239"},
    {"ibeta", "0.001", "50", "0.5", "0x1fffffffffffffffeb58a617f0efef9affeca5e38a288b507323fe50d2ebp-237"},
    {"ibeta", "400", "600", "0.375", "0x1adc195d15bd1fcd36a25c2eb71fa70b106eece9dee55e563613a3bbf209p-241"},
    {"ibetac", "2.5", "3.5", "0.875", "0x12d53224bc5bec7b2c13201d7fcf46f98bb7e15a41318e24e01cc7d51bfdp-244"},
    {"ibetac", "30", "2", "0.0625", "0x1ffffffffffffffffffffffffffffc5c0000000000000000000000000000p-237"},
};

// Value and flags in every rounding mode.
static void expect(int line, const std::function<float256()>& f, float256 value, unsigned want) {
    for (round_mode m : kModes) {
        fenv::scoped_round g(m);
        fenv::clear_flags();
        const float256 y = f();
        const unsigned got = fenv::test_flags();
        fenv::clear_flags();
        if (same(y, value) && got == want) { ++g_pass; continue; }
        ++g_fails;
        std::cerr << "FAIL " << line << " mode " << int(m) << ": got " << y.to_string(10) << " flags " << got
                  << ", want " << value.to_string(10) << " flags " << want << "\n";
    }
}
#define EXPECT(expr, value, flags) expect(__LINE__, [&] { return (expr); }, (value), (flags))
// Expected value depending on the rounding mode.
static void expect_m(int line, const std::function<float256()>& f, const std::function<float256(round_mode)>& value,
                     unsigned want) {
    for (round_mode m : kModes) {
        fenv::scoped_round g(m);
        fenv::clear_flags();
        const float256 y = f();
        const unsigned got = fenv::test_flags();
        fenv::clear_flags();
        const float256 v = value(m);
        if (same(y, v) && got == want) { ++g_pass; continue; }
        ++g_fails;
        std::cerr << "FAIL " << line << " mode " << int(m) << ": got " << y.to_string(10) << " flags " << got
                  << ", want " << v.to_string(10) << " flags " << want << "\n";
    }
}
#define EXPECT_M(expr, valuefn, flags) expect_m(__LINE__, [&] { return (expr); }, (valuefn), (flags))
static bool down(round_mode m) { return m == round_mode::toward_zero || m == round_mode::downward; }

int main() {
    using namespace base256;
    // Reference values: round-to-nearest result within 1/2 ulp (bitwise equal),
    // the directed results bracket it and only inexact is raised.
    for (const Ref& r : kRefs) {
        const float256 want = P(r.want);
        float256 got[5];
        unsigned fl[5];
        for (int m = 0; m < 5; ++m) {
            fenv::scoped_round g(kModes[m]);
            fenv::clear_flags();
            got[m] = eval(r.op, r.a, r.b, r.c);
            fl[m] = fenv::test_flags();
            fenv::clear_flags();
        }
        const bool ok = same(got[0], want) && got[3] >= want && got[4] <= want &&
                        (got[3] == want || got[3] == nextafter(want, float256::inf(false))) &&
                        (got[4] == want || got[4] == nextafter(want, float256::inf(true)));
        bool flags_ok = true;
        for (unsigned f : fl) flags_ok = flags_ok && f == kInx;
        if (ok && flags_ok) { ++g_pass; continue; }
        ++g_fails;
        std::cerr << "FAIL ref " << r.op << "(" << r.a << ", " << r.b << (r.c[0] ? ", " : "") << r.c << ") = "
                  << got[0].to_hex_string() << " want " << want.to_hex_string() << " flags " << fl[0] << "\n";
    }

    const float256 one = float256::one(), zero = float256::zero(), inf = float256::inf(false);
    const float256 nan = float256::qnan();
    // Bessel: special arguments.
    EXPECT(cyl_bessel_j(zero, zero), one, 0);
    EXPECT(cyl_bessel_j(float256(2.5), zero), zero, 0);
    EXPECT(cyl_bessel_i(zero, zero), one, 0);
    EXPECT(cyl_bessel_i(float256(3), zero), zero, 0);
    EXPECT(cyl_neumann(float256(1), zero), -inf, kDiv);
    EXPECT(cyl_bessel_k(float256(0.5), zero), inf, kDiv);
    EXPECT(cyl_bessel_j(float256(1), inf), zero, 0);
    EXPECT(cyl_neumann(float256(1), inf), zero, 0);
    EXPECT(cyl_bessel_i(float256(1), inf), inf, 0);
    EXPECT(cyl_bessel_k(float256(1), inf), zero, 0);
    EXPECT(cyl_bessel_j(float256(-1), one), nan, kInv);
    EXPECT(cyl_bessel_j(one, float256(-1)), nan, kInv);
    EXPECT(cyl_bessel_k(one, float256(-1)), nan, kInv);
    EXPECT(cyl_bessel_j(inf, one), nan, kInv);
    EXPECT(cyl_bessel_j(nan, one), nan, 0);
    EXPECT(cyl_bessel_j(one, float256::snan()), nan, kInv);
    EXPECT(sph_bessel(0, zero), one, 0);
    EXPECT(sph_bessel(4, zero), zero, 0);
    EXPECT(sph_neumann(2, zero), -inf, kDiv);
    EXPECT(yn(3, zero), -inf, kDiv);
    EXPECT(yn(-3, zero), inf, kDiv);
    EXPECT(yn(2, float256(-1)), nan, kInv);
    EXPECT(jn(-2, zero), zero, 0);
    // Parity: J_n(-x) = (-1)^n J_n(x), J_{-n} = (-1)^n J_n in every mode.
    for (round_mode m : kModes) {
        fenv::scoped_round g(m);
        const float256 x = float256(7.25);
        CHECK(same(jn(3, -x), -[&] { fenv::scoped_round h(m == round_mode::upward ? round_mode::downward
                                                          : m == round_mode::downward ? round_mode::upward : m);
                                       return jn(3, x); }()));
        CHECK(same(jn(4, -x), jn(4, x)));
        CHECK(same(jn(-4, x), jn(4, x)));
        CHECK(same(yn(-2, x), yn(2, x)));
    }
    // Tiny argument: J_0(x) = 1 - x^2/4 is just below 1.
    EXPECT_M(cyl_bessel_j(zero, ldexp(one, -400)),
             [](round_mode m) { return down(m) ? nextafter(float256::one(), float256::zero()) : float256::one(); }, kInx);
    EXPECT_M(cyl_bessel_i(zero, ldexp(one, -400)),
             [](round_mode m) { return m == round_mode::upward ? nextafter(float256::one(), float256::inf(false)) : float256::one(); }, kInx);
    // Overflow / underflow.
    EXPECT_M(cyl_bessel_i(zero, float256(200000)),
             [](round_mode m) { return down(m) ? float256::max() : float256::inf(false); }, base256::fe_overflow | kInx);
    EXPECT_M(cyl_bessel_k(zero, float256(200000)),
             [](round_mode m) { return m == round_mode::upward ? float256::denorm_min() : float256::zero(); },
             base256::fe_underflow | kInx);

    // Beta and incomplete beta: special arguments.
    EXPECT(beta(zero, one), inf, kDiv);
    EXPECT(beta(float256(-1), one), nan, kInv);
    EXPECT(beta(inf, one), zero, 0);
    EXPECT(beta(one, float256(4)), float256(0.25), 0);
    EXPECT(beta(one, one), one, 0);
    EXPECT(beta(nan, one), nan, 0);
    EXPECT(ibeta(float256(2), float256(3), zero), zero, 0);
    EXPECT(ibeta(float256(2), float256(3), one), one, 0);
    EXPECT(ibetac(float256(2), float256(3), zero), one, 0);
    EXPECT(ibeta(zero, float256(3), float256(0.5)), nan, kInv);
    EXPECT(ibeta(one, float256(3), float256(1.5)), nan, kInv);
    EXPECT(ibeta(one, float256(3), float256(-0.5)), nan, kInv);
    // I_x(a, b) + I_{1−x}(b, a) = 1
    for (double xv : {0.1, 0.4, 0.7, 0.95}) {
        const float256 x(xv), a(2.75), b(4.5);
        CHECK(abs(ibeta(a, b, x) - ibetac(b, a, float256::one() - x)) <= ldexp(one, -236));
    }

    std::cout << g_pass << " passed, " << g_fails << " failed\n";
    return g_fails ? 1 : 0;
}
