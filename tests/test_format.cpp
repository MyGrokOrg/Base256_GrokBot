// Decimal/hex output, std::format, to_chars/from_chars, integer conversions,
// totalOrder.
#include <base256/base256.hpp>

#include <chrono>
#include <cstdint>
#include <format>
#include <iostream>
#include <limits>
#include <random>
#include <string>

using base256::float256;
using base256::round_mode;
namespace fenv = base256::fenv;
using namespace base256::literals;

static int g_fails = 0, g_pass = 0;
#define CHECK(cond)                                                              \
    do {                                                                         \
        if (cond) ++g_pass;                                                      \
        else { ++g_fails; std::cerr << "FAIL " << __LINE__ << "  " #cond "\n"; } \
    } while (0)
#define CHECK_EQ(a, b)                                                                     \
    do {                                                                                   \
        const auto va_ = (a);                                                              \
        const auto vb_ = (b);                                                              \
        if (va_ == vb_) ++g_pass;                                                          \
        else { ++g_fails; std::cerr << "FAIL " << __LINE__ << "  " #a " == " #b "\n    lhs=" \
                                    << va_ << "\n    rhs=" << vb_ << "\n"; }               \
    } while (0)

static float256 rnd(std::mt19937_64& g) {
    std::uint64_t w3 = g() & ~(0x7FFFFULL << 44);
    std::uint64_t be;
    switch (g() % 6) {
        case 0: be = 0; break;                                  // subnormal
        case 1: be = 0x7FFFE - g() % 64; break;                 // near max
        case 2: be = 1 + g() % 64; break;                       // near min normal
        default: be = 262143 - 2000 + g() % 4000; break;        // moderate
    }
    w3 |= be << 44;
    return float256::from_bits(g(), g(), g(), w3);
}

int main() {
    // --- shortest round-trip ------------------------------------------------
    CHECK_EQ(float256(0.5).to_string(), std::string("5e-1"));
    CHECK_EQ("0.1"_f256.to_string(), std::string("1e-1"));
    CHECK_EQ(float256(123.456).to_string().substr(0, 12), std::string("1.2345600000"));  // the double's value
    CHECK_EQ("123.456"_f256.to_string(), std::string("1.23456e+2"));
    CHECK_EQ(float256::one().to_string(), std::string("1e+0"));
    CHECK_EQ((-"2.5e-300"_f256).to_string(), std::string("-2.5e-300"));
    CHECK_EQ(float256::denorm_min().to_string(), std::string("2e-78984"));
    {
        std::mt19937_64 g(2024);
        int bad = 0, nonminimal = 0;
        for (int i = 0; i < 400; ++i) {
            const float256 x = rnd(g);
            const std::string s = x.to_string();
            const auto back = float256::parse(s);
            if (!back || !(back->raw_bits() == x.raw_bits())) { ++bad; std::cerr << "rt " << s << "\n"; }
            // minimality: one digit fewer must not round-trip
            const std::size_t e = s.find('e');
            std::string mant = s.substr(0, e);
            std::string digs;
            for (char c : mant) if (c >= '0' && c <= '9') digs.push_back(c);
            if (digs.size() > 1) {
                const auto shorter = x.decimal_significant(static_cast<int>(digs.size()) - 1,
                                                           round_mode::to_nearest_even);
                std::string t = shorter.digits.substr(0, 1);
                if (shorter.digits.size() > 1) t += "." + shorter.digits.substr(1);
                t += "e" + std::to_string(shorter.exp10);
                if (float256::parse(t)->raw_bits() == x.raw_bits()) ++nonminimal;
            }
        }
        CHECK(bad == 0);
        CHECK(nonminimal == 0);
    }

    // --- fixed precision, rounding modes in decimal output -----------------
    const float256 third = float256::one() / float256(3.0);
    CHECK_EQ(third.to_string(5), std::string("3.3333e-1"));
    CHECK_EQ((float256(2.0) / float256(3.0)).to_string(5), std::string("6.6667e-1"));
    {
        fenv::scoped_round g(round_mode::upward);
        CHECK_EQ(third.to_string(5), std::string("3.3334e-1"));
        CHECK_EQ((-third).to_string(5), std::string("-3.3333e-1"));
    }
    {
        fenv::scoped_round g(round_mode::toward_zero);
        CHECK_EQ((float256(2.0) / float256(3.0)).to_string(5), std::string("6.6666e-1"));
    }
    fenv::clear_flags();
    (void)float256(0.5).to_string(3);
    CHECK(fenv::test_flags() == 0);
    (void)third.to_string(3);
    CHECK(fenv::test_flags() == base256::fe_inexact);
    fenv::clear_flags();

    // --- format() / to_chars ----------------------------------------------
    using cf = std::chars_format;
    CHECK_EQ(float256(1234.5).format(cf::fixed, 2), std::string("1234.50"));
    CHECK_EQ(float256(1234.5).format(cf::fixed), std::string("1234.5"));
    CHECK_EQ(float256(0.00125).format(cf::fixed, 3), std::string("0.001"));   // tie -> even
    CHECK_EQ(float256(0.00375).format(cf::fixed, 4), std::string("0.0037"));  // double 0.00375 < tie
    CHECK_EQ(float256(2.5).format(cf::fixed, 0), std::string("2"));
    CHECK_EQ(float256(1234.5).format(cf::scientific, 3), std::string("1.234e+03"));
    CHECK_EQ(float256(1234.5).format(cf::scientific), std::string("1.2345e+03"));
    CHECK_EQ(float256(1e-5).format(cf::general, 3), std::string("1e-05"));
    CHECK_EQ(float256(123456.0).format(cf::general, 3), std::string("1.23e+05"));
    CHECK_EQ(float256(0.0001).format(cf::general, 3), std::string("0.0001"));
    CHECK_EQ(float256(100.0).format(cf::general), std::string("100"));
    CHECK_EQ("1e30"_f256.format(cf::general), std::string("1e+30"));
    CHECK_EQ(float256(3.0).format(cf::hex), std::string("1.8p+1"));
    CHECK_EQ(float256(3.0).format(cf::hex, 3), std::string("1.800p+1"));
    CHECK_EQ(float256::zero(true).format(cf::fixed, 2), std::string("-0.00"));
    CHECK_EQ(float256::inf(true).format(cf::general), std::string("-inf"));
    CHECK_EQ(float256::qnan().format(cf::general), std::string("nan"));
    {
        char buf[16];
        auto r = base256::to_chars(buf, buf + sizeof buf, third);
        CHECK(r.ec == std::errc::value_too_large);
        char big[128];
        r = base256::to_chars(big, big + sizeof big, float256(0.25));
        CHECK(r.ec == std::errc{} && std::string(big, r.ptr) == "0.25");
    }

    // --- std::format --------------------------------------------------------
    CHECK_EQ(std::format("{}", float256(0.25)), std::string("0.25"));
    CHECK_EQ(std::format("{:.3f}", third), std::string("0.333"));
    CHECK_EQ(std::format("{:e}", float256(1.5)), std::string("1.500000e+00"));
    CHECK_EQ(std::format("{:E}", float256(1.5)), std::string("1.500000E+00"));
    CHECK_EQ(std::format("{:g}", float256(1.5)), std::string("1.5"));
    CHECK_EQ(std::format("{:a}", float256(3.0)), std::string("0x1.8p+1"));
    CHECK_EQ(std::format("{:>8.2f}", float256(3.14159)), std::string("    3.14"));
    CHECK_EQ(std::format("{:*^9.1f}", float256(2.0)), std::string("***2.0***"));
    CHECK_EQ(std::format("{:+.1f}", float256(2.0)), std::string("+2.0"));
    CHECK_EQ(std::format("{:08.2f}", float256(-3.5)), std::string("-0003.50"));
    CHECK_EQ(std::format("{}", float256::inf()), std::string("inf"));

    // --- from_chars -----------------------------------------------------------
    {
        float256 v;
        const std::string s1 = "1.5e3xyz";
        auto r = base256::from_chars(s1.data(), s1.data() + s1.size(), v);
        CHECK(r.ec == std::errc{} && r.ptr == s1.data() + 5 && v == float256(1500.0));
        const std::string s2 = "abc";
        r = base256::from_chars(s2.data(), s2.data() + s2.size(), v);
        CHECK(r.ec == std::errc::invalid_argument && r.ptr == s2.data());
        const std::string s3 = "1.8p1";
        r = base256::from_chars(s3.data(), s3.data() + s3.size(), v, std::chars_format::hex);
        CHECK(r.ec == std::errc{} && v == float256(3.0));
        const std::string s4 = "1e99999";
        v = float256(7.0);
        r = base256::from_chars(s4.data(), s4.data() + s4.size(), v);
        CHECK(r.ec == std::errc::result_out_of_range && v == float256(7.0));
        const std::string s5 = "1e-99999";
        r = base256::from_chars(s5.data(), s5.data() + s5.size(), v);
        CHECK(r.ec == std::errc::result_out_of_range);
        const std::string s6 = "-infinity";
        r = base256::from_chars(s6.data(), s6.data() + s6.size(), v);
        CHECK(r.ec == std::errc{} && v.isinf() && v.signbit() && r.ptr == s6.data() + 9);
        const std::string s7 = "12e";
        r = base256::from_chars(s7.data(), s7.data() + s7.size(), v);
        CHECK(r.ec == std::errc{} && r.ptr == s7.data() + 2 && v == float256(12.0));
        const std::string s8 = "1.25";
        r = base256::from_chars(s8.data(), s8.data() + s8.size(), v, std::chars_format::scientific);
        CHECK(r.ec == std::errc::invalid_argument);
        // to_chars -> from_chars round trip
        std::mt19937_64 g(7);
        int bad = 0;
        for (int i = 0; i < 200; ++i) {
            const float256 x = rnd(g);
            char buf[512];
            for (auto f : {cf::general, cf::scientific, cf::hex}) {
                auto w = base256::to_chars(buf, buf + sizeof buf, x, f);
                float256 y;
                auto rr = base256::from_chars(buf, w.ptr, y, f);
                if (rr.ec != std::errc{} || !(y.raw_bits() == x.raw_bits())) ++bad;
            }
        }
        CHECK(bad == 0);
    }

    // --- parse of extreme exponents is fast and correctly rounded -------------
    {
        const auto t0 = std::chrono::steady_clock::now();
        CHECK(float256::parse("1e-1000000000").value().iszero());
        CHECK(float256::parse("1e1000000000").value().isinf());
        CHECK(float256::parse("-1e999999999999").value().isinf());
        {
            fenv::scoped_round g(round_mode::upward);
            CHECK(float256::parse("1e-1000000000").value() == float256::denorm_min());
        }
        {
            fenv::scoped_round g(round_mode::toward_zero);
            CHECK(float256::parse("1e1000000000").value() == float256::max());
        }
        const auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - t0).count();
        CHECK(ms < 2000);
        CHECK(float256::parse("1.61132571748576047361957211845200501064402387454966951747637125049607183e78913").value() == float256::max());
    }

    // --- integer conversions --------------------------------------------------
    CHECK(float256(static_cast<long long>(-5)) == float256(-5.0));
    CHECK(float256(static_cast<unsigned char>(200)) == float256(200.0));
    CHECK(float256(static_cast<short>(-7)) == float256(-7.0));
    CHECK(float256(std::numeric_limits<std::int64_t>::min()) == -base256::ldexp(float256::one(), 63));
    CHECK(float256(std::numeric_limits<std::uint64_t>::max()) ==
          base256::ldexp(float256::one(), 64) - float256::one());
#if defined(__SIZEOF_INT128__)
    {
        const base256::detail::u128_t u = ~static_cast<base256::detail::u128_t>(0);
        CHECK(float256(u) == base256::ldexp(float256::one(), 128) - float256::one());
        const auto m = static_cast<base256::detail::i128_t>(static_cast<base256::detail::u128_t>(1) << 127);
        CHECK(float256(m) == -base256::ldexp(float256::one(), 127));
        CHECK(float256(m).to_integer<base256::detail::i128_t>() == m);
        CHECK(float256(u).to_integer<base256::detail::u128_t>() == u);
    }
#endif
    CHECK(float256(2.7).to_integer<int>() == 2);
    CHECK(float256(-2.7).to_integer<int>() == -2);
    CHECK(float256(2.5).to_integer<int>(round_mode::to_nearest_even) == 2);
    CHECK(float256(2.5).to_integer<int>(round_mode::to_nearest_away) == 3);
    CHECK(float256(-2.5).to_integer<int>(round_mode::downward) == -3);
    CHECK(float256(-2.1).to_integer<int>(round_mode::upward) == -2);
    fenv::clear_flags();
    CHECK(float256(1e30).to_integer<std::int64_t>() == std::numeric_limits<std::int64_t>::max());
    CHECK(fenv::test_flags(base256::fe_invalid));
    fenv::clear_flags();
    CHECK(float256(-1.0).to_integer<unsigned>() == 0u);
    CHECK(fenv::test_flags(base256::fe_invalid));
    fenv::clear_flags();
    CHECK(float256::qnan().to_integer<int>() == 0);
    CHECK(fenv::test_flags(base256::fe_invalid));
    fenv::clear_flags();
    CHECK(float256(-128.0).to_integer<signed char>() == -128);
    CHECK(fenv::test_flags() == 0);
    (void)float256(2.5).to_integer<int>(round_mode::toward_zero, true);
    CHECK(fenv::test_flags() == base256::fe_inexact);
    fenv::clear_flags();
    (void)float256(2.5).to_integer<int>();
    CHECK(fenv::test_flags() == 0);

    // --- totalOrder -------------------------------------------------------------
    CHECK(base256::total_order(-float256::zero(), float256::zero()));
    CHECK(!base256::total_order(float256::zero(), -float256::zero()));
    CHECK(base256::total_order(float256::inf(), float256::qnan()));
    CHECK(base256::total_order(-float256::qnan(), -float256::inf()));
    CHECK(base256::total_order(float256::snan(), float256::qnan()));
    CHECK(base256::total_order(float256(1.0), float256(2.0)));
    CHECK(!base256::total_order(float256(2.0), float256(1.0)));
    CHECK(base256::total_order(float256(-2.0), float256(-1.0)));
    CHECK(base256::total_order(float256(1.0), float256(1.0)));
    CHECK(base256::total_order_mag(float256(-1.0), float256(2.0)));

    std::cout << "passed=" << g_pass << " failed=" << g_fails << "\n";
    return g_fails ? 1 : 0;
}
