#include <base256/base256.hpp>

#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>

using base256::float256;
using namespace base256::literals;

static int g_fails = 0;
static int g_pass  = 0;

#define CHECK(cond)                                                            \
    do {                                                                       \
        if (cond) {                                                            \
            ++g_pass;                                                          \
        } else {                                                               \
            ++g_fails;                                                         \
            std::cerr << "FAIL " << __FILE__ << ":" << __LINE__ << "  "        \
                      << #cond << "\n";                                        \
        }                                                                      \
    } while (0)

#define CHECK_EQ(a, b)                                                         \
    do {                                                                       \
        const auto _va = (a);                                                  \
        const auto _vb = (b);                                                  \
        if (_va == _vb) {                                                      \
            ++g_pass;                                                          \
        } else {                                                               \
            ++g_fails;                                                         \
            std::cerr << "FAIL " << __FILE__ << ":" << __LINE__ << "  "        \
                      << #a << " == " << #b << "\n    lhs=" << _va             \
                      << "\n    rhs=" << _vb << "\n";                          \
        }                                                                      \
    } while (0)

int main() {
    using L = std::numeric_limits<float256>;

    // --- specials ----------------------------------------------------------
    CHECK(float256::zero().iszero());
    CHECK((-float256::zero()).iszero());
    CHECK(float256::zero().signbit() == false);
    CHECK((-float256::zero()).signbit() == true);
    CHECK(float256::inf().isinf());
    CHECK((-float256::inf()).isinf());
    CHECK(float256::qnan().isnan());
    CHECK(!float256::qnan().isfinite());
    CHECK(L::has_infinity);
    CHECK(L::digits == 237);
    CHECK(L::digits10 == 71);
    CHECK(L::is_iec559);

    // --- integer / double round-trip ---------------------------------------
    CHECK_EQ(float256(0).to_double(), 0.0);
    CHECK_EQ(float256(1).to_double(), 1.0);
    CHECK_EQ(float256(2).to_double(), 2.0);
    CHECK_EQ(float256(-3).to_double(), -3.0);
    CHECK_EQ(float256(42.0).to_double(), 42.0);
    CHECK_EQ(float256(-0.5).to_double(), -0.5);
    CHECK_EQ(float256(3.141592653589793).to_double(), 3.141592653589793);

    // --- constants vs double -----------------------------------------------
    CHECK(std::fabs(float256::pi().to_double() - 3.141592653589793) < 1e-15);
    CHECK(std::fabs(float256::e().to_double()  - 2.718281828459045) < 1e-15);
    CHECK(std::fabs(float256::ln2().to_double() - 0.6931471805599453) < 1e-15);

    // --- 1.0 encoding ------------------------------------------------------
    {
        auto b = float256::one().raw_bits();
        CHECK(b[0] == 0 && b[1] == 0 && b[2] == 0);
        CHECK(b[3] == 0x3ffff00000000000ULL);
    }
    {
        auto b = float256::two().raw_bits();
        CHECK(b[3] == 0x4000000000000000ULL);
    }

    // --- arithmetic: small integers ----------------------------------------
    CHECK_EQ((float256(1) + float256(1)).to_double(), 2.0);
    CHECK_EQ((float256(10) - float256(3)).to_double(), 7.0);
    CHECK_EQ((float256(6) * float256(7)).to_double(), 42.0);
    CHECK_EQ((float256(84) / float256(2)).to_double(), 42.0);
    CHECK_EQ((float256(1) / float256(2)).to_double(), 0.5);

    // --- signed zero -------------------------------------------------------
    CHECK((float256::zero() + float256::zero()).signbit() == false);
    CHECK(((-float256::zero()) + (-float256::zero())).signbit() == true);
    CHECK(float256::zero() == -float256::zero());

    // --- inf / nan propagation ---------------------------------------------
    CHECK((float256::inf() + float256(1)).isinf());
    CHECK((float256::inf() - float256::inf()).isnan());
    CHECK((float256::inf() * float256::zero()).isnan());
    CHECK((float256(1) / float256::zero()).isinf());
    CHECK((float256::zero() / float256::zero()).isnan());
    CHECK(!(float256::qnan() == float256::qnan()));
    CHECK(!(float256::qnan() < float256(1)));

    // --- precision beyond double -------------------------------------------
    // ulp(1.0) for double is 2^-52; for binary256 it is 2^-236.
    {
        const float256 one = float256::one();
        const float256 delta = base256::ldexp(float256::one(), -80);
        const float256 sum = one + delta;
        CHECK(sum != one);
        CHECK((sum - one) == delta);
        // double cannot represent 1 + 2^-80
        CHECK(1.0 + std::ldexp(1.0, -80) == 1.0);
    }

    // --- 1/3 * 3 recovers 1 ------------------------------------------------
    {
        const float256 third = float256::one() / float256(3);
        const float256 back  = third * float256(3);
        CHECK(back == float256::one());
    }

    // --- parse / format ----------------------------------------------------
    {
        auto p = float256::parse("1.0");
        CHECK(p.has_value());
        CHECK(*p == float256::one());
        auto p2 = float256::parse("2");
        CHECK(p2.has_value() && *p2 == float256::two());
        auto pi = float256::parse("3.14159265358979323846");
        CHECK(pi.has_value());
        CHECK(std::fabs(pi->to_double() - 3.141592653589793) < 1e-15);
        auto infs = float256::parse("-inf");
        CHECK(infs.has_value() && infs->isinf() && infs->signbit());
    }

    // 0.1 + 0.2 is closer to 0.3 than IEEE-754 binary64
    {
        auto a = float256::parse("0.1");
        auto b = float256::parse("0.2");
        auto c = float256::parse("0.3");
        CHECK(a && b && c);
        const float256 s = *a + *b;
        const float256 err256 = base256::abs(s - *c);
        const double err64 = std::fabs((0.1 + 0.2) - 0.3);
        CHECK(err256.to_double() < err64);
        CHECK(err256.to_double() < 1e-20);
    }

    // user-defined literal
    {
        auto x = 2.5_f256;
        CHECK(x.to_double() == 2.5);
        auto y = "0.125"_f256;
        CHECK(y.to_double() == 0.125);
    }

    // --- sqrt --------------------------------------------------------------
    {
        const float256 s = base256::sqrt(float256(4));
        CHECK(s == float256(2));
        const float256 s2 = base256::sqrt(float256(9));
        CHECK(s2 == float256(3));
        const float256 s0 = base256::sqrt(float256::zero());
        CHECK(s0.iszero());
        CHECK(base256::sqrt(-float256::one()).isnan());
        const float256 r2 = base256::sqrt(float256(2));
        const float256 back = r2 * r2;
        CHECK(base256::abs(back - float256(2)).to_double() < 1e-15);
        // high-precision: |sqrt(2)^2 - 2| should be at most a few ulps of 2
        const float256 ulp2 = base256::ldexp(float256::one(), 1 - 236);
        CHECK(base256::abs(back - float256(2)) <= ulp2 * float256(8));
    }

    // --- comparisons / nextafter -------------------------------------------
    CHECK(float256(1) < float256(2));
    CHECK(float256(-1) < float256(1));
    CHECK(float256::max() > float256::one());
    {
        const float256 n = base256::nextafter(float256::one(), float256(2));
        CHECK(n > float256::one());
        CHECK(base256::nextafter(n, float256::one()) == float256::one());
    }

    // --- ldexp / frexp -----------------------------------------------------
    {
        int e = 0;
        const float256 m = base256::frexp(float256(8), &e);
        CHECK(e == 4);
        CHECK(m.to_double() == 0.5);
        CHECK(base256::ldexp(m, e) == float256(8));
    }

    // --- floor / ceil / trunc ----------------------------------------------
    CHECK(base256::floor(float256(1.5)) == float256(1));
    CHECK(base256::ceil(float256(1.5)) == float256(2));
    CHECK(base256::trunc(float256(-1.5)) == float256(-1));
    CHECK(base256::floor(float256(-1.5)) == float256(-2));

    // --- exp / log ---------------------------------------------------------
    {
        const float256 e1 = base256::exp(float256::one());
        CHECK(std::fabs(e1.to_double() - 2.718281828459045) < 1e-12);
        const float256 back = base256::log(e1);
        CHECK(std::fabs(back.to_double() - 1.0) < 1e-12);
        const float256 l2 = base256::log(float256(2));
        CHECK(std::fabs((l2 - float256::ln2()).to_double()) < 1e-12);
    }

    // --- sin / cos ---------------------------------------------------------
    {
        CHECK(std::fabs(base256::sin(float256::zero()).to_double()) < 1e-30);
        CHECK(std::fabs(base256::cos(float256::zero()).to_double() - 1.0) < 1e-15);
        const float256 s = base256::sin(float256::pi() / float256(2));
        CHECK(std::fabs(s.to_double() - 1.0) < 1e-10);
        const float256 c = base256::cos(float256::pi());
        CHECK(std::fabs(c.to_double() + 1.0) < 1e-10);
    }

    // --- subnormals exist --------------------------------------------------
    {
        const float256 dmin = float256::denorm_min();
        CHECK(dmin.is_subnormal());
        CHECK(dmin > float256::zero());
        CHECK(base256::nextafter(float256::zero(), float256::one()) == dmin);
    }

    // --- max / min / epsilon -----------------------------------------------
    CHECK(float256::max().isnormal());
    CHECK(float256::min().isnormal());
    CHECK((float256::one() + float256::epsilon()) != float256::one());
    CHECK((float256::one() + float256::epsilon() * float256::half()) == float256::one() ||
          (float256::one() + float256::epsilon() * float256::half()) != float256::one());
    // 1 + eps/2 ties-to-even should stay 1 (lsb of 1.0 is 0)
    CHECK((float256::one() + float256::epsilon() * float256::half()) == float256::one());

    // --- to_string ---------------------------------------------------------
    {
        CHECK(float256(1).to_string(1).rfind("1", 0) == 0);
        const std::string s2 = float256(2).to_string(6);
        CHECK(s2.find("2.00000") != std::string::npos || s2.find("2e+") != std::string::npos);
        const std::string spi = float256::pi().to_string(20);
        CHECK(spi.find("3.141592653589793238") != std::string::npos);
        const std::string s01 = float256::parse("0.1")->to_string(20);
        CHECK(s01.find("1.0000000000000000000e-1") != std::string::npos);
    }

    // --- hex round-trip ----------------------------------------------------
    {
        const auto h = float256::pi().to_hex_string();
        auto p = float256::parse(h);
        CHECK(p.has_value());
        CHECK(*p == float256::pi());
    }

    // --- decimal round-trip of powers of ten (to_string exponent) ---------
    {
        int pow_fail = 0;
        for (int e = -40; e <= 40; ++e) {
            const std::string src = "1e" + std::to_string(e);
            auto p = float256::parse(src);
            if (!p) { ++pow_fail; continue; }
            const std::string printed = p->to_string(73);
            auto q = float256::parse(printed);
            if (!q || *q != *p) ++pow_fail;
            // The printed scientific exponent must be the true order of magnitude.
            const std::string shortp = p->to_string(8);
            const auto ep = shortp.rfind('e');
            if (ep == std::string::npos || std::stoi(shortp.substr(ep + 1)) != e)
                ++pow_fail;
        }
        CHECK(pow_fail == 0);
        const auto hundredth = float256::parse("0.01");
        CHECK(hundredth.has_value());
        CHECK(hundredth->to_string(8).find("1.0000000e-2") != std::string::npos);
        const auto tiny = float256::parse("1e-28");
        CHECK(tiny.has_value());
        CHECK(tiny->to_string(12).find("e-28") != std::string::npos);
    }

    // --- fmin / fmax signed zero (IEEE 754 minNum / maxNum) ----------------
    {
        const float256 pz = float256::zero(false);
        const float256 nz = float256::zero(true);
        CHECK(base256::fmin(pz, nz).signbit());
        CHECK(base256::fmin(nz, pz).signbit());
        CHECK(!base256::fmax(pz, nz).signbit());
        CHECK(!base256::fmax(nz, pz).signbit());
        CHECK(base256::fmin(float256(2), float256::qnan()) == float256(2));
        CHECK(base256::fmax(float256::qnan(), float256(-3)) == float256(-3));
    }

    // --- hypot does not overflow or underflow on a pure scale -------------
    {
        const float256 big = base256::ldexp(float256::one(), 200000);
        const float256 h = base256::hypot(big, big);
        CHECK(h.isfinite());
        const float256 expect = big * float256::sqrt2();
        const float256 ulp = base256::ldexp(float256::one(), 200000 - 236);
        CHECK(base256::abs(h - expect) <= ulp * float256(4));

        const float256 tiny = base256::ldexp(float256::one(), -200000);
        const float256 ht = base256::hypot(tiny, tiny);
        CHECK(ht.isfinite());
        CHECK(!ht.iszero());
        const float256 expect_t = tiny * float256::sqrt2();
        CHECK(base256::abs(ht - expect_t) <=
              base256::ldexp(float256::one(), -200000 - 236) * float256(4));

        CHECK(base256::hypot(float256::inf(), float256(1)).isinf());
        CHECK(!base256::hypot(float256::inf(), float256(1)).signbit());
        CHECK(base256::hypot(float256::qnan(), float256(1)).isnan());
        CHECK(base256::hypot(float256::inf(), float256::qnan()).isinf());
        CHECK(base256::hypot(float256::qnan(), float256::inf()).isinf());
    }

    // --- pow of signed zero and negative bases -----------------------------
    {
        const float256 nz = float256::zero(true);
        const float256 p3 = base256::pow(nz, float256(3));
        CHECK(p3.iszero() && p3.signbit());
        const float256 p2 = base256::pow(nz, float256(2));
        CHECK(p2.iszero() && !p2.signbit());
        const float256 n1 = base256::pow(nz, -float256::one());
        CHECK(n1.isinf() && n1.signbit());
        const float256 n2 = base256::pow(nz, -float256::two());
        CHECK(n2.isinf() && !n2.signbit());
        CHECK(base256::pow(-float256::two(), float256(3)) == -float256(8));
        CHECK(base256::pow(-float256::two(), float256::half()).isnan());
    }

    // --- atan2 on infinities -----------------------------------------------
    {
        const float256 pinf = float256::inf();
        const float256 ninf = float256::inf(true);
        const float256 pi4 = float256::pi() * float256::half() * float256::half();
        CHECK(base256::atan2(pinf, pinf) == pi4);
        CHECK(base256::atan2(ninf, pinf) == -pi4);
        CHECK(base256::atan2(pinf, ninf) == pi4 * float256(3));
        CHECK(base256::atan2(ninf, ninf) == -(pi4 * float256(3)));
        CHECK(base256::atan2(float256::one(), pinf).iszero());
        CHECK(!base256::atan2(float256::one(), pinf).signbit());
        CHECK(base256::atan2(float256::one(), ninf) == float256::pi());
    }

    // --- trig at exact multiples of the library π --------------------------
    {
        CHECK(base256::sin(float256::pi()).iszero() && !base256::sin(float256::pi()).signbit());
        CHECK(base256::sin(-float256::pi()).iszero() && base256::sin(-float256::pi()).signbit());
        CHECK(base256::cos(float256::pi() * float256::half()).iszero());
        CHECK(!base256::cos(float256::pi() * float256::half()).signbit());
        CHECK(base256::cos(float256::pi()) == -float256::one());
        CHECK(base256::tan(float256::pi()).iszero() && !base256::tan(float256::pi()).signbit());
        CHECK(base256::tan(float256::pi() * float256::half()).isinf());
        CHECK(!base256::tan(float256::pi() * float256::half()).signbit());
    }

    // --- rint (ties to even), modf, nextup ---------------------------------
    {
        CHECK(base256::rint(float256(1.5)) == float256(2));
        CHECK(base256::rint(float256(2.5)) == float256(2));
        CHECK(base256::rint(float256(3.5)) == float256(4));
        CHECK(base256::rint(float256(-1.5)) == float256(-2));
        CHECK(base256::rint(float256(-2.5)) == float256(-2));
        CHECK(base256::rint(float256::half()).iszero() && !base256::rint(float256::half()).signbit());
        CHECK(base256::rint(-float256::half()).iszero() && base256::rint(-float256::half()).signbit());
        CHECK(base256::round(float256(2.5)) == float256(3));

        float256 ip;
        CHECK(base256::modf(float256(3.25), &ip) == float256(0.25));
        CHECK(ip == float256(3));
        CHECK(base256::modf(float256(-3.25), &ip) == float256(-0.25));
        CHECK(ip == float256(-3));
        const float256 zf = base256::modf(float256::inf(true), &ip);
        CHECK(ip.isinf() && ip.signbit());
        CHECK(zf.iszero() && zf.signbit());

        CHECK(base256::nextup(float256::one()) == base256::nextafter(float256::one(), float256::inf()));
        CHECK(base256::nextdown(float256::one()) == base256::nextafter(float256::one(), -float256::inf()));
    }

    // --- log1p keeps arguments that 1+x flushes away -----------------
    {
        CHECK(base256::log1p(float256::zero()).iszero());
        CHECK(!base256::log1p(float256::zero()).signbit());
        CHECK(base256::log1p(-float256::zero()).iszero());
        CHECK(base256::log1p(-float256::zero()).signbit());
        const float256 tiny = base256::ldexp(float256::one(), -300);
        CHECK(base256::log1p(tiny) == tiny);
        CHECK(base256::log1p(-tiny) == -tiny);
        const float256 x80 = base256::ldexp(float256::one(), -80);
        const float256 merc = x80 - (x80 * x80) * float256::half();
        CHECK(base256::abs(base256::log1p(x80) - merc) <=
              base256::ldexp(float256::one(), -220));
        CHECK(std::fabs(base256::log1p(float256(-0.125)).to_double() -
                        std::log1p(-0.125)) < 1e-15);
        CHECK(base256::log1p(-float256::one()).isinf());
        CHECK(base256::log1p(float256(-1.5)).isnan());
    }

    // --- pow specials and exact integer powers ----------------------------
    {
        CHECK(base256::pow(float256::qnan(), float256::zero()) == float256::one());
        CHECK(base256::pow(float256::one(), float256::qnan()) == float256::one());
        CHECK(base256::pow(-float256::one(), float256::inf()) == float256::one());
        CHECK(base256::pow(-float256::one(), float256::inf(true)) == float256::one());
        CHECK(base256::pow(float256::two(), float256::inf()).isinf());
        CHECK(base256::pow(float256::two(), float256::inf(true)).iszero());
        CHECK(base256::pow(float256::half(), float256::inf()).iszero());
        CHECK(base256::pow(float256::half(), float256::inf(true)).isinf());
        CHECK(base256::pow(float256::inf(), float256(3)).isinf());
        CHECK(base256::pow(float256::inf(), -float256::one()).iszero());
        CHECK(base256::pow(float256::inf(true), float256(3)).isinf());
        CHECK(base256::pow(float256::inf(true), float256(3)).signbit());
        CHECK(base256::pow(float256::inf(true), float256(2)).isinf());
        CHECK(!base256::pow(float256::inf(true), float256(2)).signbit());
        CHECK(base256::pow(float256::inf(true), float256::half()).isinf());
        CHECK(!base256::pow(float256::inf(true), float256::half()).signbit());

        CHECK(base256::pow(float256::two(), float256(100)) ==
              base256::ldexp(float256::one(), 100));
        CHECK(base256::pow(float256::two(), float256(-80)) ==
              base256::ldexp(float256::one(), -80));
        const float256 weird = float256::one() + base256::ldexp(float256::one(), -80);
        CHECK(base256::pow(weird, float256::one()) == weird);
        CHECK(base256::pow(weird, float256::two()) == weird * weird);
        auto ten20 = float256::parse("1e20");
        CHECK(ten20.has_value());
        CHECK(base256::pow(float256(10), float256(20)) == *ten20);
        CHECK(base256::pow(-float256::two(), float256(-3)) == -float256::one() / float256(8));
    }

    // --- regression: effective subtraction with lost (sticky) bits --------
    {
        const float256 one = float256::one();
        // exact 1 - 2^-238 - 2^-400 lies just below the midpoint 1 - 2^-238
        const float256 r = one - (base256::ldexp(one, -238) + base256::ldexp(one, -400));
        CHECK(r == base256::nextdown(one));
        const float256 r2 = one + (-(base256::ldexp(one, -238) + base256::ldexp(one, -400)));
        CHECK(r2 == base256::nextdown(one));
        // same through fma: 1*1 - (2^-238 + 2^-400)
        const float256 c = -(base256::ldexp(one, -238) + base256::ldexp(one, -400));
        CHECK(base256::fma(one, one, c) == base256::nextdown(one));
        // exactly at the midpoint still ties to even (up to 1)
        CHECK(one - base256::ldexp(one, -238) == one);
    }

    // --- regression: subnormal results that round up ---------------------
    {
        const float256 dm = float256::denorm_min();
        CHECK(dm * float256(1.5) == dm + dm);            // tie -> even (2 ulp)
        CHECK(dm * float256(2.5) == dm + dm);            // tie -> even (2 ulp)
        CHECK(dm * float256(0.75) == dm);
        CHECK((dm * float256(0.5)).iszero());            // tie -> even (0)
        // guard bit exactly at bit 64 of the working significand (rshift 65)
        const float256 x = base256::ldexp(float256::one() + base256::ldexp(float256::one(), -191)
                                          + base256::ldexp(float256::one(), -200), -262188);
        CHECK(x == base256::ldexp(float256::one() + base256::ldexp(float256::one(), -190), -262188));
        // binary64 subnormal conversion rounds once
        const double dmin = std::numeric_limits<double>::denorm_min();
        CHECK(float256(dmin).to_double() == dmin);
        CHECK((float256(dmin) * float256(1.5)).to_double() == 2 * dmin);
        const float256 y = base256::ldexp(float256::one() - base256::ldexp(float256::one(), -200), -1022 - 43);
        CHECK(y.to_double() == std::ldexp(1.0, -1065));
    }

    // --- regression: fma with an exact zero product ----------------------
    {
        const float256 z = float256::zero(), nz = -float256::zero();
        CHECK(!base256::fma(z, z, nz).signbit());   // +0 + -0 = +0
        CHECK(!base256::fma(nz, z, z).signbit());   // -0 + +0 = +0
        CHECK(base256::fma(nz, z, nz).signbit());   // -0 + -0 = -0
        CHECK(base256::fma(z, float256::one(), float256(3.0)) == float256(3.0));
    }

    // --- regression: numeric _f256 literals are parsed exactly -----------
    {
        CHECK(0.1_f256 == "0.1"_f256);
        CHECK(0.1_f256 != float256(0.1));
        CHECK(3.14159265358979323846264338327950288419716939937510582097494459_f256 ==
              "3.14159265358979323846264338327950288419716939937510582097494459"_f256);
        CHECK(1e-100_f256 == "1e-100"_f256);
        CHECK(0x1.8p+1_f256 == float256(3.0));
        CHECK(42_f256 == float256(42.0));
        CHECK(1'000.5_f256 == float256(1000.5));
    }

    // --- regression: exact long double construction ----------------------
    {
        const long double tenth = 0.1L;
        const float256 f(tenth);
        CHECK(f.to_long_double() == tenth);
        if constexpr (std::numeric_limits<long double>::digits > 53) {
            CHECK(f != float256(0.1));
            CHECK(f == float256::parse(std::numeric_limits<long double>::digits == 64
                                           ? "0x1.999999999999999ap-4"
                                           : "0x1.999999999999999999999999999ap-4").value());
        }
        const long double lmin = std::numeric_limits<long double>::denorm_min();
        CHECK(float256(lmin).to_long_double() == lmin);
        CHECK(float256(-lmin).signbit());
        const long double lmax = std::numeric_limits<long double>::max();
        CHECK(float256(lmax).to_long_double() == lmax);
        CHECK(float256(std::numeric_limits<long double>::infinity()).isinf());
        CHECK(float256(std::numeric_limits<long double>::quiet_NaN()).isnan());
        CHECK(float256(-0.0L).signbit());
        CHECK(float256(1.0L / 3.0L).to_long_double() == 1.0L / 3.0L);
    }

    // --- signaling NaN handling and payload propagation ------------------
    {
        const float256 one = float256::one();
        const float256 sn = float256::snan();
        CHECK(sn.issignaling());
        CHECK(!float256::qnan().issignaling());
        CHECK(!sn.quieted().issignaling() && sn.quieted().isnan());
        CHECK(base256::fmin(sn, one).isnan());
        CHECK(base256::fmax(one, sn).isnan());
        CHECK(!base256::fmin(sn, one).issignaling());
        CHECK(base256::fmin(float256::qnan(), one) == one);
        CHECK(base256::minimum(float256::qnan(), one).isnan());
        CHECK(base256::minimum_number(sn, one) == one);
        CHECK(base256::minimum(-float256::zero(), float256::zero()).signbit());
        CHECK(!base256::maximum(-float256::zero(), float256::zero()).signbit());
        // arithmetic: result is the first NaN operand's payload, quieted
        auto bits = float256::snan().raw_bits();
        bits[0] = 0xABCDULL;
        const float256 payload = float256::from_bits(bits);
        CHECK(payload.issignaling());
        const float256 r = payload + one;
        CHECK(r.isnan() && !r.issignaling());
        CHECK(r.raw_bits()[0] == 0xABCDULL);
        CHECK((one * payload).raw_bits()[0] == 0xABCDULL);
        CHECK((-payload).issignaling());          // sign ops are quiet-computational
        CHECK(base256::abs(payload).issignaling());
    }

    std::cout << "passed=" << g_pass << " failed=" << g_fails << "\n";
    return g_fails ? 1 : 0;
}
