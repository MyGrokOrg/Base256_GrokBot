// base256::complex: C Annex G special values, branch cuts, signed zeros, IEEE flags
// in all five rounding modes, exact cases, and stream / std::format I/O.
#include <base256/complex.hpp>
#include <functional>
#include <iostream>
#include <numbers>
#include <sstream>
#include <string>

using base256::complex;
using base256::float256;
using base256::round_mode;
namespace fenv = base256::fenv;

static int g_fails = 0, g_pass = 0;
#define CHECK(cond)                                                              \
    do {                                                                         \
        if (cond) ++g_pass;                                                      \
        else { ++g_fails; std::cerr << "FAIL " << __LINE__ << "  " #cond "\n"; } \
    } while (0)

constexpr unsigned kInv = base256::fe_invalid, kDiv = base256::fe_divbyzero, kOvf = base256::fe_overflow,
                   kUnf = base256::fe_underflow, kInx = base256::fe_inexact, kAll = kInv | kDiv | kOvf | kUnf | kInx;
const round_mode kModes[5] = {round_mode::to_nearest_even, round_mode::to_nearest_away, round_mode::toward_zero,
                              round_mode::upward, round_mode::downward};

// Expected component: a value compared bitwise (NaN by isnan, zero and infinity by sign),
// or the same with the sign left open where Annex G says "the sign is unspecified".
struct V {
    float256 v;
    bool any_sign = false;
    std::function<float256()> f{};  // when set: the reference, evaluated in the mode under test
    V(float256 x, bool a = false) : v(x), any_sign(a) {}
    V(std::function<float256()> g, bool a = false) : v(), any_sign(a), f(std::move(g)) {}
    V at_mode() const { return f ? V(f(), any_sign) : *this; }
};
static V anys(float256 v) { return {v, true}; }
// −f() correctly rounded in the current mode: f evaluated in the mirrored mode, then negated.
static float256 neg_rounded(const std::function<float256()>& f) {
    const round_mode m = fenv::get_round();
    const round_mode mm = m == round_mode::upward ? round_mode::downward
                          : m == round_mode::downward ? round_mode::upward : m;
    fenv::scoped_round g(mm);
    return -f();
}
static bool match(float256 got, const V& w) {
    if (w.v.isnan()) return got.isnan();
    if (got.isnan()) return false;
    if (w.any_sign) return abs(got) == abs(w.v);
    return got == w.v && got.signbit() == w.v.signbit();
}

// Evaluate f in all five modes; components must match and (flags & mask) must equal want.
static void expect(const char* what, int line, const std::function<complex()>& f, V re, V im, unsigned want,
                   unsigned mask) {
    const V re0 = re, im0 = im;
    for (round_mode m : kModes) {
        fenv::scoped_round g(m);
        re = re0.at_mode();
        im = im0.at_mode();
        fenv::clear_flags();
        const complex z = f();
        const unsigned got = fenv::test_flags() & mask;
        fenv::clear_flags();
        if (match(z.real(), re) && match(z.imag(), im) && got == want) { ++g_pass; continue; }
        ++g_fails;
        std::cerr << "FAIL " << line << " " << what << " mode " << int(m) << ": got " << z << " flags " << got
                  << ", want (" << re.v << (re.any_sign ? "[±]" : "") << "," << im.v << (im.any_sign ? "[±]" : "")
                  << ") flags " << want << "\n";
    }
}
// EXPECT: exact flags. EXPECTM: only the flags in mask are checked (π-valued results are inexact).
#define EXPECT(expr, re, im, flags) expect(#expr, __LINE__, [&] { return (expr); }, V(re), V(im), (flags), kAll)
#define EXPECTA(expr, re, im, flags) expect(#expr, __LINE__, [&] { return (expr); }, (re), (im), (flags), kAll)
#define EXPECTM(expr, re, im, flags, mask) \
    expect(#expr, __LINE__, [&] { return (expr); }, V(re), V(im), (flags), (mask))
#define EXPECTMA(expr, re, im, flags, mask) \
    expect(#expr, __LINE__, [&] { return (expr); }, (re), (im), (flags), (mask))

// Per-mode check against a reference float256 expression evaluated in the same mode.
static void expect_modes(const char* what, int line, const std::function<complex()>& f,
                         const std::function<complex()>& ref, unsigned want) {
    for (round_mode m : kModes) {
        fenv::scoped_round g(m);
        const complex r = ref();
        fenv::clear_flags();
        const complex z = f();
        const unsigned got = fenv::test_flags();
        fenv::clear_flags();
        if (match(z.real(), V{r.real()}) && match(z.imag(), V{r.imag()}) && got == want) { ++g_pass; continue; }
        ++g_fails;
        std::cerr << "FAIL " << line << " " << what << " mode " << int(m) << ": got " << z << " flags " << got
                  << ", want " << r << " flags " << want << "\n";
    }
}
#define EXPECT_MODES(expr, ref, flags) \
    expect_modes(#expr, __LINE__, [&] { return (expr); }, [&] { return (ref); }, (flags))


// v + s·ε (ε > 0 far below an ulp of v) rounded in the current mode.
static V near(float256 v, int s) {
    return V(std::function<float256()>([v, s] {
        switch (fenv::get_round()) {
            case round_mode::upward: return s > 0 ? nextup(v) : v;
            case round_mode::downward: return s < 0 ? nextdown(v) : v;
            case round_mode::toward_zero: return ((s > 0) == v.signbit()) ? nextafter(v, float256()) : v;
            default: return v;
        }
    }));
}

int main() {
    using namespace base256;
    using namespace base256::literals;
    const float256 one = float256::one(), z0{}, n0 = -z0, two(2.0), inf = float256::inf(), ninf = -inf,
                   nan = float256::qnan(), snan = float256::snan();
    // π multiples correctly rounded in the mode under test (real atan2/atan are correctly rounded)
    const V pi([=] { return atan2(z0, -one); }), mpi([=] { return atan2(n0, -one); }),
        pi2([=] { return atan2(one, z0); }), mpi2([=] { return atan2(-one, z0); }),
        pi4([=] { return atan(one); }), pi34([=] { return atan2(one, -one); });
    auto RPI = [=] { return atan2(z0, -one); };
    auto RPI2 = [=] { return atan2(one, z0); };
    auto cx = [](float256 a, float256 b) { return complex(a, b); };
    const unsigned ID = kInv | kDiv;  // mask for results containing a rounded π multiple

    // ---- construction, access, conversions, literals --------------------------------
    {
        constexpr complex c(float256(1.5), float256(-2.0));
        static_assert(c.real() == float256(1.5) && c.imag() == float256(-2.0));
        complex d;
        CHECK(d.real() == z0 && d.imag() == z0 && !d.real().signbit());
        d.real(two);
        d.imag(-one);
        CHECK(d == cx(two, -one));
        CHECK(real(d) == two && imag(d) == -one);
        CHECK(2.5_i256 == cx(z0, float256(2.5)));
        const complex fromstd(std::complex<double>(0.5, -3.0));
        CHECK(fromstd == cx(float256(0.5), float256(-3.0)));
        const auto back = static_cast<std::complex<double>>(fromstd);
        CHECK(back == std::complex<double>(0.5, -3.0));
        CHECK(complex(one) == one);
        CHECK(!(complex(one, one) == one));
        CHECK(conj(cx(one, z0)).imag().signbit());
        CHECK(proj(cx(one, two)) == cx(one, two));
        const complex p = proj(cx(nan, ninf));
        CHECK(p.real() == inf && p.imag().iszero() && p.imag().signbit());
        CHECK(proj(cx(ninf, nan)).imag().iszero() && !proj(cx(ninf, nan)).imag().signbit());
    }

    // ---- arithmetic: exact cases and Annex G mixed-operand rules ----------------------
    EXPECT(cx(one, two) * cx(float256(3), float256(-1)), float256(5), float256(5), 0u);
    EXPECT(cx(float256(5), float256(5)) / cx(float256(3), float256(-1)), one, two, 0u);
    // x − (u + iv) = (x − u) − iv: the imaginary part keeps a negated zero sign
    EXPECT_MODES(one - cx(one, z0), cx(one - one, n0), 0u);
    // exact cancellation in the real part: +0, except −0 rounding downward (checked per mode)
    for (round_mode m : kModes) {
        fenv::scoped_round g(m);
        const complex r = cx(one, one) * cx(one, -one);  // 2 + 0i
        CHECK(r.real() == two && r.imag().iszero() && r.imag().signbit() == (m == round_mode::downward));
    }
    // a single rounding of a·c − b·d: (1 + 2^-200 i)² = (1 − 2^-400) + 2^-199 i
    {
        const float256 t = ldexp(one, -200);
        EXPECT_MODES(cx(one, t) * cx(one, t), cx(one - ldexp(one, -400), ldexp(one, -199)), kInx);
        // naive evaluation would give exactly 1 in every mode; the downward/toward-zero results differ
        fenv::scoped_round g(round_mode::downward);
        CHECK((cx(one, t) * cx(one, t)).real() < one);
    }
    // division: a single rounding per component, 1/(3) both parts in each mode
    EXPECT_MODES(cx(one, one) / cx(float256(3), z0), cx(one / float256(3), one / float256(3)), kInx);
    EXPECT_MODES(norm(cx(one, ldexp(one, -200))) * one, cx(one + ldexp(one, -400), z0), kInx);
    // Annex G G.5.1 infinities and zeros
    {
        const complex r = cx(inf, inf) * cx(one, z0);
        CHECK(r.real().isinf() || r.imag().isinf());
        const complex q = cx(one, one) / cx(z0, z0);
        CHECK(q.real().isinf() && q.imag().isinf());
        const complex s = cx(one, one) / cx(inf, inf);
        CHECK(s.real().iszero() && s.imag().iszero());
        const complex u = cx(inf, nan) * cx(two, z0);
        CHECK(u.real().isinf() || u.imag().isinf());
    }
    EXPECTA(cx(one, one) / cx(z0, z0), V{inf}, V{inf}, kDiv);
    EXPECTA(cx(z0, z0) / cx(z0, z0), V{nan}, V{nan}, kInv);
    // overflow and underflow of the exact product
    EXPECT_MODES(cx(ldexp(one, 200000), z0) * cx(ldexp(one, 100000), ldexp(one, 100000)),
                 cx(ldexp(one, 200000) * ldexp(one, 100000), ldexp(one, 200000) * ldexp(one, 100000)),
                 kOvf | kInx);
    // signalling NaN inputs: invalid, quiet NaN results
    {
        fenv::clear_flags();
        const complex r = exp(cx(snan, one));
        CHECK(r.real().isnan() && !r.real().issignaling() && (fenv::test_flags() & kInv));
        fenv::clear_flags();
    }

    // ---- abs / norm / arg / polar --------------------------------------------------
    EXPECT(complex(abs(cx(float256(3), float256(-4))), z0), float256(5), z0, 0u);
    EXPECT(complex(norm(cx(float256(3), float256(-4))), z0), float256(25), z0, 0u);
    EXPECT(complex(abs(cx(ninf, nan)), z0), inf, z0, 0u);
    EXPECT(complex(abs(cx(nan, inf)), z0), inf, z0, 0u);
    EXPECT(complex(arg(cx(one, z0)), z0), z0, z0, 0u);
    EXPECT(complex(arg(cx(one, n0)), z0), n0, z0, 0u);
    EXPECTM(complex(arg(cx(n0, n0)), z0), mpi, z0, 0u, ID);
    EXPECT_MODES(complex(arg(cx(-one, z0)), z0), cx(RPI(), z0), kInx);
    EXPECT_MODES(complex(abs(cx(one, one)), z0), cx(sqrt(two), z0), kInx);
    EXPECT(polar(two, z0), two, z0, 0u);
    EXPECT_MODES(polar(two, one), cx(two * cos(one), two * sin(one)), kInx);

    // ---- cexp (G.6.3.1) --------------------------------------------------------------
    EXPECT(exp(cx(z0, z0)), one, z0, 0u);
    EXPECT(exp(cx(n0, n0)), one, n0, 0u);
    EXPECT(exp(cx(one, inf)), nan, nan, kInv);
    EXPECT(exp(cx(one, nan)), nan, nan, 0u);
    EXPECT(exp(cx(inf, z0)), inf, z0, 0u);
    EXPECT(exp(cx(inf, n0)), inf, n0, 0u);
    EXPECT_MODES(exp(cx(ninf, one)), cx(z0 * cos(one), z0 * sin(one)), 0u);
    EXPECT_MODES(exp(cx(inf, two)), cx(inf * cos(two), inf * sin(two)), 0u);
    EXPECTA(exp(cx(ninf, inf)), anys(z0), anys(z0), 0u);
    EXPECTA(exp(cx(inf, inf)), anys(inf), V{nan}, kInv);
    EXPECTA(exp(cx(ninf, nan)), anys(z0), anys(z0), 0u);
    EXPECTA(exp(cx(inf, nan)), anys(inf), V{nan}, 0u);
    EXPECT(exp(cx(nan, z0)), nan, z0, 0u);
    EXPECT(exp(cx(nan, n0)), nan, n0, 0u);
    EXPECT(exp(cx(nan, one)), nan, nan, 0u);
    EXPECT_MODES(exp(cx(one, z0)), cx(exp(one), z0), kInx);
    EXPECT_MODES(exp(cx(z0, one)), cx(cos(one), sin(one)), kInx);
    EXPECT_MODES(exp(cx(float256(200000), one)), cx(exp(float256(200000)), exp(float256(200000))), kOvf | kInx);

    // ---- clog (G.6.3.2) ----------------------------------------------------------------
    EXPECTM(log(cx(n0, z0)), ninf, pi, kDiv, ID);
    EXPECTM(log(cx(n0, n0)), ninf, mpi, kDiv, ID);
    EXPECT(log(cx(z0, z0)), ninf, z0, kDiv);
    EXPECT(log(cx(z0, n0)), ninf, n0, kDiv);
    EXPECTM(log(cx(one, inf)), inf, pi2, 0u, ID);
    EXPECTM(log(cx(ninf, one)), inf, pi, 0u, ID);
    EXPECT(log(cx(inf, one)), inf, z0, 0u);
    EXPECTM(log(cx(ninf, inf)), inf, pi34, 0u, ID);
    EXPECTM(log(cx(inf, inf)), inf, pi4, 0u, ID);
    EXPECT(log(cx(inf, nan)), inf, nan, 0u);
    EXPECT(log(cx(ninf, nan)), inf, nan, 0u);
    EXPECT(log(cx(nan, inf)), inf, nan, 0u);
    EXPECT(log(cx(nan, one)), nan, nan, 0u);
    EXPECT(log(cx(nan, nan)), nan, nan, 0u);
    EXPECT(log(cx(one, z0)), z0, z0, 0u);
    EXPECT(log(cx(one, n0)), z0, n0, 0u);
    EXPECT_MODES(log(cx(-one, z0)), cx(z0, RPI()), kInx);   // branch cut, upper side
    EXPECT_MODES(log(cx(-one, n0)), cx(z0, atan2(n0, -one)), kInx);  // lower side
    EXPECT_MODES(log(cx(two, z0)), cx(log(two), z0), kInx);
    EXPECT(log10(cx(float256(100), z0)), two, z0, 0u);
    for (round_mode m : kModes) {  // |z| = 10: the real part is exactly 1 in every mode
        fenv::scoped_round g(m);
        CHECK(log10(cx(float256(-6), float256(8))).real() == one);
    }
    EXPECT_MODES(log10(cx(two, z0)), cx(log10(two), z0), kInx);

    // ---- csqrt (G.6.4.2) -----------------------------------------------------------------
    EXPECT(sqrt(cx(z0, z0)), z0, z0, 0u);
    EXPECT(sqrt(cx(n0, z0)), z0, z0, 0u);
    EXPECT(sqrt(cx(n0, n0)), z0, n0, 0u);
    EXPECT(sqrt(cx(one, inf)), inf, inf, 0u);
    EXPECT(sqrt(cx(nan, inf)), inf, inf, 0u);
    EXPECT(sqrt(cx(nan, ninf)), inf, ninf, 0u);
    EXPECT(sqrt(cx(ninf, one)), z0, inf, 0u);
    EXPECT(sqrt(cx(ninf, -one)), z0, ninf, 0u);
    EXPECT(sqrt(cx(inf, one)), inf, z0, 0u);
    EXPECTA(sqrt(cx(ninf, nan)), V{nan}, anys(inf), 0u);
    EXPECT(sqrt(cx(inf, nan)), inf, nan, 0u);
    EXPECT(sqrt(cx(nan, one)), nan, nan, 0u);
    EXPECT(sqrt(cx(float256(-4), z0)), z0, two, 0u);  // branch cut, exact
    EXPECT(sqrt(cx(float256(-4), n0)), z0, -two, 0u);
    EXPECT(sqrt(cx(float256(-5), float256(12))), two, float256(3), 0u);
    EXPECT(sqrt(cx(float256(4), z0)), two, z0, 0u);
    EXPECT_MODES(sqrt(cx(two, z0)), cx(sqrt(two), z0), kInx);
    EXPECT_MODES(sqrt(cx(-two, z0)), cx(z0, sqrt(two)), kInx);

    // ---- csinh (G.6.2.5), ccosh (G.6.2.4), ctanh (G.6.2.6) -----------------------------
    EXPECT(sinh(cx(z0, z0)), z0, z0, 0u);
    EXPECT(sinh(cx(n0, n0)), n0, n0, 0u);
    EXPECTA(sinh(cx(z0, inf)), anys(z0), V{nan}, kInv);
    EXPECTA(sinh(cx(z0, nan)), anys(z0), V{nan}, 0u);
    EXPECT(sinh(cx(one, inf)), nan, nan, kInv);
    EXPECT(sinh(cx(one, nan)), nan, nan, 0u);
    EXPECT(sinh(cx(inf, z0)), inf, z0, 0u);
    EXPECT(sinh(cx(ninf, n0)), ninf, n0, 0u);
    EXPECT_MODES(sinh(cx(inf, two)), cx(inf * cos(two), inf * sin(two)), 0u);
    EXPECTA(sinh(cx(inf, inf)), anys(inf), V{nan}, kInv);
    EXPECTA(sinh(cx(inf, nan)), anys(inf), V{nan}, 0u);
    EXPECT(sinh(cx(nan, z0)), nan, z0, 0u);
    EXPECT(sinh(cx(nan, one)), nan, nan, 0u);
    EXPECT_MODES(sinh(cx(one, z0)), cx(sinh(one), z0), kInx);
    EXPECT_MODES(sinh(cx(z0, one)), cx(z0, sin(one)), kInx);

    EXPECT(cosh(cx(z0, z0)), one, z0, 0u);
    EXPECT(cosh(cx(n0, z0)), one, n0, 0u);
    EXPECTA(cosh(cx(z0, inf)), V{nan}, anys(z0), kInv);
    EXPECTA(cosh(cx(z0, nan)), V{nan}, anys(z0), 0u);
    EXPECT(cosh(cx(one, inf)), nan, nan, kInv);
    EXPECT(cosh(cx(one, nan)), nan, nan, 0u);
    EXPECT(cosh(cx(inf, z0)), inf, z0, 0u);
    EXPECT(cosh(cx(ninf, z0)), inf, n0, 0u);
    EXPECT_MODES(cosh(cx(inf, two)), cx(inf * cos(two), inf * sin(two)), 0u);
    EXPECTA(cosh(cx(inf, inf)), anys(inf), V{nan}, kInv);
    EXPECT(cosh(cx(inf, nan)), inf, nan, 0u);
    EXPECTA(cosh(cx(nan, z0)), V{nan}, anys(z0), 0u);
    EXPECT(cosh(cx(nan, one)), nan, nan, 0u);
    EXPECT_MODES(cosh(cx(one, z0)), cx(cosh(one), z0), kInx);
    EXPECT_MODES(cosh(cx(z0, one)), cx(cos(one), z0), kInx);

    EXPECT(tanh(cx(z0, z0)), z0, z0, 0u);
    EXPECT(tanh(cx(n0, n0)), n0, n0, 0u);
    EXPECT(tanh(cx(z0, inf)), z0, nan, kInv);  // C17/C23 (DR 471): +0 + i∞ gives +0 + iNaN
    EXPECT(tanh(cx(one, inf)), nan, nan, kInv);
    EXPECT(tanh(cx(z0, nan)), z0, nan, 0u);
    EXPECT(tanh(cx(one, nan)), nan, nan, 0u);
    EXPECT(tanh(cx(inf, one)), one, z0, 0u);   // sin(2) > 0
    EXPECT(tanh(cx(inf, two)), one, n0, 0u);   // sin(4) < 0
    EXPECTA(tanh(cx(inf, inf)), V{one}, anys(z0), 0u);
    EXPECTA(tanh(cx(inf, nan)), V{one}, anys(z0), 0u);
    EXPECT(tanh(cx(ninf, one)), -one, z0, 0u);
    EXPECT(tanh(cx(nan, z0)), nan, z0, 0u);
    EXPECT(tanh(cx(nan, one)), nan, nan, 0u);
    EXPECT_MODES(tanh(cx(one, z0)), cx(tanh(one), z0), kInx);
    EXPECT_MODES(tanh(cx(z0, one)), cx(z0, tan(one)), kInx);

    // ---- sin / cos / tan via the hyperbolic identities: signed zeros on the real axis ---
    EXPECT_MODES(sin(cx(one, z0)), cx(sin(one), z0), kInx);
    EXPECT_MODES(sin(cx(-two, z0)), cx(sin(-two), n0), kInx);  // cos(−2) < 0
    EXPECT_MODES(cos(cx(one, z0)), cx(cos(one), n0), kInx);    // −sin(1)·sinh(+0)
    EXPECT_MODES(cos(cx(-one, z0)), cx(cos(one), z0), kInx);
    EXPECT_MODES(tan(cx(one, z0)), cx(tan(one), z0), kInx);
    EXPECT_MODES(sin(cx(z0, one)), cx(z0, sinh(one)), kInx);
    EXPECT_MODES(cos(cx(z0, one)), cx(cosh(one), n0), kInx);
    EXPECT_MODES(tan(cx(z0, one)), cx(z0, tanh(one)), kInx);
    EXPECT(sin(cx(z0, z0)), z0, z0, 0u);
    EXPECT(cos(cx(z0, z0)), one, n0, 0u);
    EXPECT(tan(cx(n0, z0)), n0, z0, 0u);
    EXPECT(sin(cx(inf, z0)), nan, z0, kInv);  // −i·csinh(i·z) with csinh(+0 + i∞) = ±0 + iNaN
    EXPECT(cos(cx(nan, z0)), nan, z0, 0u);

    // ---- cacos (G.6.1.1) and casin --------------------------------------------------------
    EXPECTM(acos(cx(z0, z0)), pi2, n0, 0u, ID);
    EXPECTM(acos(cx(n0, n0)), pi2, z0, 0u, ID);
    EXPECTM(acos(cx(z0, nan)), pi2, nan, 0u, ID);
    EXPECTM(acos(cx(one, inf)), pi2, ninf, 0u, ID);
    EXPECT(acos(cx(one, nan)), nan, nan, 0u);
    EXPECTM(acos(cx(ninf, one)), pi, ninf, 0u, ID);
    EXPECT(acos(cx(inf, one)), z0, ninf, 0u);
    EXPECT(acos(cx(inf, -one)), z0, inf, 0u);
    EXPECTM(acos(cx(ninf, inf)), pi34, ninf, 0u, ID);
    EXPECTM(acos(cx(inf, inf)), pi4, ninf, 0u, ID);
    EXPECTA(acos(cx(inf, nan)), V{nan}, anys(inf), 0u);
    EXPECT(acos(cx(nan, one)), nan, nan, 0u);
    EXPECT(acos(cx(nan, inf)), nan, ninf, 0u);
    EXPECT(acos(cx(nan, nan)), nan, nan, 0u);
    EXPECT(acos(cx(one, z0)), z0, n0, 0u);
    EXPECT_MODES(acos(cx(two, z0)), cx(z0, neg_rounded([&] { return acosh(two); })), kInx);  // branch cut x > 1
    EXPECT_MODES(acos(cx(two, n0)), cx(z0, acosh(two)), kInx);
    EXPECT_MODES(acos(cx(float256(0.5), z0)), cx(acos(float256(0.5)), n0), kInx);
    EXPECT(asin(cx(z0, z0)), z0, z0, 0u);
    EXPECT(asin(cx(n0, n0)), n0, n0, 0u);
    EXPECTM(asin(cx(two, inf)), z0, inf, 0u, ID);
    EXPECTM(asin(cx(inf, one)), pi2, inf, 0u, ID);
    EXPECT_MODES(asin(cx(two, z0)), cx(RPI2(), acosh(two)), kInx);
    EXPECT_MODES(asin(cx(float256(0.5), z0)), cx(asin(float256(0.5)), z0), kInx);

    // ---- casinh (G.6.2.2) ------------------------------------------------------------------
    EXPECT(asinh(cx(z0, z0)), z0, z0, 0u);
    EXPECT(asinh(cx(n0, n0)), n0, n0, 0u);
    EXPECTM(asinh(cx(one, inf)), inf, pi2, 0u, ID);
    EXPECTM(asinh(cx(-one, ninf)), ninf, mpi2, 0u, ID);
    EXPECT(asinh(cx(one, nan)), nan, nan, 0u);
    EXPECT(asinh(cx(inf, one)), inf, z0, 0u);
    EXPECTM(asinh(cx(inf, inf)), inf, pi4, 0u, ID);
    EXPECT(asinh(cx(inf, nan)), inf, nan, 0u);
    EXPECT(asinh(cx(nan, z0)), nan, z0, 0u);
    EXPECT(asinh(cx(nan, n0)), nan, n0, 0u);
    EXPECT(asinh(cx(nan, one)), nan, nan, 0u);
    EXPECTA(asinh(cx(nan, inf)), anys(inf), V{nan}, 0u);
    EXPECT(asinh(cx(nan, nan)), nan, nan, 0u);
    EXPECT_MODES(asinh(cx(z0, two)), cx(acosh(two), RPI2()), kInx);   // cut, right side
    EXPECT_MODES(asinh(cx(n0, two)), cx(neg_rounded([&] { return acosh(two); }), RPI2()), kInx);  // left side
    EXPECT_MODES(asinh(cx(one, z0)), cx(asinh(one), z0), kInx);

    // ---- cacosh (G.6.2.1) -------------------------------------------------------------------
    EXPECTM(acosh(cx(z0, z0)), z0, pi2, 0u, ID);
    EXPECTM(acosh(cx(n0, z0)), z0, pi2, 0u, ID);
    EXPECTM(acosh(cx(z0, n0)), z0, mpi2, 0u, ID);
    EXPECTM(acosh(cx(one, inf)), inf, pi2, 0u, ID);
    EXPECT(acosh(cx(one, nan)), nan, nan, 0u);
    EXPECTM(acosh(cx(ninf, one)), inf, pi, 0u, ID);
    EXPECT(acosh(cx(inf, one)), inf, z0, 0u);
    EXPECTM(acosh(cx(ninf, inf)), inf, pi34, 0u, ID);
    EXPECTM(acosh(cx(inf, inf)), inf, pi4, 0u, ID);
    EXPECT(acosh(cx(inf, nan)), inf, nan, 0u);
    EXPECT(acosh(cx(nan, one)), nan, nan, 0u);
    EXPECT(acosh(cx(nan, inf)), inf, nan, 0u);
    EXPECT(acosh(cx(nan, nan)), nan, nan, 0u);
    EXPECT(acosh(cx(one, z0)), z0, z0, 0u);
    EXPECT_MODES(acosh(cx(two, z0)), cx(acosh(two), z0), kInx);
    EXPECT_MODES(acosh(cx(float256(-2), z0)), cx(acosh(two), RPI()), kInx);  // cut x < 1
    EXPECT_MODES(acosh(cx(float256(-2), n0)), cx(acosh(two), atan2(n0, -one)), kInx);

    // ---- catanh (G.6.2.3) and catan --------------------------------------------------------
    EXPECT(atanh(cx(z0, z0)), z0, z0, 0u);
    EXPECT(atanh(cx(n0, n0)), n0, n0, 0u);
    EXPECT(atanh(cx(z0, nan)), z0, nan, 0u);
    EXPECT(atanh(cx(one, z0)), inf, z0, kDiv);
    EXPECT(atanh(cx(-one, n0)), ninf, n0, kDiv);
    EXPECTM(atanh(cx(one, inf)), z0, pi2, 0u, ID);
    EXPECT(atanh(cx(one, nan)), nan, nan, 0u);
    EXPECTM(atanh(cx(inf, one)), z0, pi2, 0u, ID);
    EXPECTM(atanh(cx(inf, inf)), z0, pi2, 0u, ID);
    EXPECT(atanh(cx(inf, nan)), z0, nan, 0u);
    EXPECT(atanh(cx(nan, one)), nan, nan, 0u);
    EXPECTMA(atanh(cx(nan, inf)), anys(z0), pi2, 0u, ID);
    EXPECT(atanh(cx(nan, nan)), nan, nan, 0u);
    EXPECT_MODES(atanh(cx(two, z0)), cx(atanh(one / two), RPI2()), kInx);
    EXPECT_MODES(atanh(cx(float256(0.5), z0)), cx(atanh(float256(0.5)), z0), kInx);
    EXPECT_MODES(atan(cx(one, z0)), cx(atan(one), z0), kInx);
    EXPECT(atan(cx(z0, one)), z0, inf, kDiv);
    EXPECTM(atan(cx(inf, one)), pi2, z0, 0u, ID);


    // ---- tiny arguments: the result is the exact lead plus a rest far below an ulp,
    // so only the rounding direction shows (z = t + it, t = 2^-2000: z³ = t³(−2 + 2i))
    {
        const float256 t = ldexp(one, -2000), t2 = ldexp(one, -4000);
        EXPECTA(sinh(cx(t, t)), near(t, -1), near(t, 1), kInx);
        EXPECTA(sin(cx(t, t)), near(t, 1), near(t, -1), kInx);
        EXPECTA(tan(cx(t, t)), near(t, -1), near(t, 1), kInx);
        EXPECTA(tanh(cx(t, t)), near(t, 1), near(t, -1), kInx);
        EXPECTA(atanh(cx(t, t)), near(t, -1), near(t, 1), kInx);
        EXPECTA(atan(cx(t, t)), near(t, 1), near(t, -1), kInx);
        EXPECTA(asin(cx(t, t)), near(t, -1), near(t, 1), kInx);
        EXPECTA(asinh(cx(t, t)), near(t, 1), near(t, -1), kInx);
        EXPECTA(complex(z0, acos(cx(t, t)).imag()), V(z0), near(-t, -1), kInx);
        EXPECTA(exp(cx(t, t)), near(one, 1), near(t, 1), kInx);
        EXPECTA(cosh(cx(t, t)), near(one, -1), near(t2, -1), kInx);
        EXPECTA(cos(cx(t, t)), near(one, -1), near(-t2, 1), kInx);  // 1 − t⁴/6 and −t² + t⁶/90
        EXPECTA(log(cx(one, t)), near(ldexp(t2, -1), -1), near(t, -1), kInx);
        EXPECTA(polar(two, t), near(two, -1), near(ldexp(t, 1), -1), kInx);
    }

    // ---- pow --------------------------------------------------------------------------------
    // integer powers are exact Gaussian powers: an exactly cancelling part is +0 (−0 rounding downward)
    EXPECT_MODES(pow(cx(one, one), complex(two)), cx(one - one, two), 0u);
    EXPECT_MODES(pow(cx(one, one), float256(4)), cx(float256(-4), one - one), 0u);
    // real z > 0, real w: the imaginary part is the zero w·Im z (exp(w·log z))
    EXPECT(pow(cx(two, z0), complex(float256(-2))), float256(0.25), n0, 0u);
    EXPECT(pow(cx(two, z0), complex(float256(3))), float256(8), z0, 0u);
    EXPECT(pow(cx(float256(-2), z0), complex(float256(3))), float256(-8), z0, 0u);
    EXPECT(pow(cx(float256(-2), n0), complex(float256(3))), float256(-8), n0, 0u);
    EXPECT(pow(cx(float256(3), float256(7)), complex(z0)), one, z0, 0u);
    EXPECT(pow(cx(z0, z0), complex(two, one)), z0, z0, 0u);
    EXPECT(pow(cx(z0, z0), complex(-two)), inf, z0, kDiv);
    EXPECT(pow(cx(float256(-4), z0), complex(float256(0.5))), z0, two, 0u);
    EXPECT(pow(float256(-8), complex(float256(3))), float256(-512), z0, 0u);
    EXPECT(pow(cx(float256(-4), z0), complex(float256(1.5))), z0, float256(-8), 0u);
    EXPECT_MODES(pow(cx(two, z0), complex(float256(0.5))), cx(sqrt(two), z0), kInx);
    EXPECT_MODES(pow(cx(two, z0), complex(one / float256(3))), cx(pow(two, one / float256(3)), z0), kInx);
    EXPECT_MODES(pow(cx(one, one), complex(float256(3))), cx(float256(-2), two), 0u);

    // ---- stream I/O and std::format ----------------------------------------------------------
    {
        std::ostringstream os;
        os << cx(float256(1.5), float256(-2));
        std::ostringstream want;
        want << '(' << float256(1.5) << ',' << float256(-2) << ')';
        CHECK(os.str() == want.str());
        complex r;
        std::istringstream is("(0.25, -3)  7  (inf,nan) (1.5)");
        is >> r;
        CHECK(r == cx(float256(0.25), float256(-3)));
        is >> r;
        CHECK(r == cx(float256(7), z0));
        is >> r;
        CHECK(r.real().isinf() && r.imag().isnan());
        is >> r;
        CHECK(r == cx(float256(1.5), z0));
        std::istringstream bad("(1;2)");
        bad >> r;
        CHECK(bad.fail());
        // round trip of non-representable decimal parts
        const complex t = cx(one / float256(3), -two / float256(7));
        std::ostringstream o2;
        o2 << t;
        std::istringstream i2(o2.str());
        complex u;
        i2 >> u;
        CHECK(u == t && u.imag().signbit());
        CHECK(std::format("{}", cx(float256(1.5), float256(-2))) == "(1.5,-2)");
        CHECK(std::format("{:+}", cx(one, two)) == "(+1,+2)");
        CHECK(std::format("{:.3f}", cx(one / float256(3), two)) == "(0.333,2.000)");
        CHECK(std::format("{:*^11}", cx(one, two)) == "***(1,2)***");
        CHECK(std::format("{:<7}|", cx(one, two)) == "(1,2)  |");
    }

    std::cout << "test_complex: " << g_pass << " passed, " << g_fails << " failed\n";
    return g_fails ? 1 : 0;
}
