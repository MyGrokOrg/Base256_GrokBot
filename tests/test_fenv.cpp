// Rounding modes and exception flags.
#include <base256/base256.hpp>
#include <iostream>

using base256::float256;
using base256::round_mode;
namespace fenv = base256::fenv;

static int g_fails = 0, g_pass = 0;
#define CHECK(cond)                                                              \
    do {                                                                         \
        if (cond) ++g_pass;                                                      \
        else { ++g_fails; std::cerr << "FAIL " << __LINE__ << "  " #cond "\n"; } \
    } while (0)

int main() {
    const float256 one = float256::one(), three(3.0);
    const float256 ulp1 = base256::ldexp(one, -236);   // ulp of 1

    // Default environment.
    CHECK(fenv::get_round() == round_mode::to_nearest_even);
    fenv::clear_flags();
    CHECK(fenv::test_flags() == 0);

    // 1/3 in every mode: RD/RZ agree, RU = RD + ulp.
    float256 q[5];
    const round_mode modes[5] = {round_mode::to_nearest_even, round_mode::to_nearest_away,
                                 round_mode::toward_zero, round_mode::upward, round_mode::downward};
    for (int i = 0; i < 5; ++i) {
        fenv::scoped_round g(modes[i]);
        q[i] = one / three;
    }
    CHECK(fenv::get_round() == round_mode::to_nearest_even);  // restored
    CHECK(q[2] == q[4]);
    CHECK(q[3] == base256::nextup(q[4]));
    CHECK(q[0] == q[4] || q[0] == q[3]);
    CHECK((-one / three) == -q[0]);
    {
        fenv::scoped_round g(round_mode::upward);
        CHECK((-one) / three == -q[4]);   // toward +inf on a negative value = smaller magnitude
    }

    // Ties: 1 + ulp/2 is a tie between 1 and 1+ulp.
    const float256 half_ulp = base256::ldexp(one, -237);
    {
        fenv::scoped_round g(round_mode::to_nearest_even);
        CHECK(one + half_ulp == one);
        CHECK(one + ulp1 + half_ulp == one + ulp1 + ulp1);
    }
    {
        fenv::scoped_round g(round_mode::to_nearest_away);
        CHECK(one + half_ulp == one + ulp1);
        CHECK(-one - half_ulp == -(one + ulp1));
    }

    // Exact zero sums: +0 except in downward mode.
    CHECK(!(one - one).signbit());
    {
        fenv::scoped_round g(round_mode::downward);
        CHECK((one - one).signbit());
        CHECK(base256::fma(one, one, -one).signbit());
    }

    // Flags.
    fenv::clear_flags();
    (void)(one / three);
    CHECK(fenv::test_flags() == base256::fe_inexact);
    fenv::clear_flags();
    (void)(one + one);
    CHECK(fenv::test_flags() == 0);
    fenv::clear_flags();
    (void)(one / float256::zero());
    CHECK(fenv::test_flags() == base256::fe_divbyzero);
    fenv::clear_flags();
    (void)(float256::zero() / float256::zero());
    CHECK(fenv::test_flags() == base256::fe_invalid);
    fenv::clear_flags();
    (void)base256::sqrt(-one);
    CHECK(fenv::test_flags(base256::fe_invalid));
    fenv::clear_flags();
    (void)(float256::inf() - float256::inf());
    CHECK(fenv::test_flags(base256::fe_invalid));
    fenv::clear_flags();
    (void)(float256::qnan() + one);
    CHECK(fenv::test_flags() == 0);                      // quiet NaN: no flag
    (void)(float256::snan() + one);
    CHECK(fenv::test_flags() == base256::fe_invalid);    // signaling NaN: invalid
    fenv::clear_flags();
    (void)(float256::qnan() < one);
    CHECK(fenv::test_flags() == base256::fe_invalid);    // ordered compare signals
    fenv::clear_flags();
    (void)(float256::qnan() == one);
    CHECK(fenv::test_flags() == 0);                      // == is quiet

    // Overflow: inf, or max in the direction toward zero.
    const float256 big = float256::max();
    fenv::clear_flags();
    CHECK((big * float256(2.0)).isinf());
    CHECK(fenv::test_flags() == (base256::fe_overflow | base256::fe_inexact));
    {
        fenv::scoped_round g(round_mode::toward_zero);
        CHECK(big * float256(2.0) == big);
        CHECK(-big * float256(2.0) == -big);
    }
    {
        fenv::scoped_round g(round_mode::downward);
        CHECK(big * float256(2.0) == big);
        CHECK((-big * float256(2.0)).isinf());
    }

    // Underflow: tiny and inexact. An exact subnormal raises nothing.
    const float256 dm = float256::denorm_min();
    fenv::clear_flags();
    (void)(dm + dm);
    CHECK(fenv::test_flags() == 0);
    fenv::clear_flags();
    CHECK((dm * float256(0.5)).iszero());
    CHECK(fenv::test_flags() == (base256::fe_underflow | base256::fe_inexact));
    {
        fenv::scoped_round g(round_mode::upward);
        CHECK(dm * float256(0.25) == dm);
        CHECK((-dm * float256(0.25)).iszero() && (-dm * float256(0.25)).signbit());
    }
    // Tininess after rounding: a result that rounds up to min normal is not tiny.
    fenv::clear_flags();
    const float256 almost = base256::nextdown(float256::min());
    const float256 r = almost * (one + ulp1);  // = (1 - 2^-472) * 2^Emin
    CHECK(r == float256::min());
    CHECK(fenv::test_flags() == base256::fe_inexact);

    // Conversions honour the mode.
    const float256 x = one - base256::ldexp(one, -100);
    {
        fenv::scoped_round g(round_mode::toward_zero);
        CHECK(x.to_double() == 1.0 - 0x1p-53);
        CHECK(float256::parse("0.1").value() < float256::parse("0.1").value() + ulp1); // sanity
    }
    {
        fenv::scoped_round g(round_mode::upward);
        CHECK(x.to_double() == 1.0);
        CHECK(float256::parse("1e-80000").value() == dm);   // parse rounds in mode
    }
    CHECK(x.to_double() == 1.0);

    // rint / nearbyint follow the mode; only rint raises inexact.
    const float256 h(2.5);
    CHECK(base256::rint(h) == float256(2.0));
    {
        fenv::scoped_round g(round_mode::upward);
        CHECK(base256::rint(h) == float256(3.0));
        CHECK(base256::rint(float256(-2.5)) == float256(-2.0));
        CHECK(base256::rint(float256(0.25)) == one);
    }
    {
        fenv::scoped_round g(round_mode::to_nearest_away);
        CHECK(base256::rint(h) == float256(3.0));
        CHECK(base256::rint(float256(0.5)) == one);
    }
    fenv::clear_flags();
    (void)base256::nearbyint(h);
    CHECK(fenv::test_flags() == 0);
    (void)base256::rint(h);
    CHECK(fenv::test_flags() == base256::fe_inexact);

    // Elementary functions run in round-to-nearest and keep the caller's flags.
    fenv::clear_flags();
    fenv::raise_flags(base256::fe_overflow);
    {
        fenv::scoped_round g(round_mode::upward);
        (void)base256::exp(one);
    }
    CHECK(fenv::test_flags() == base256::fe_overflow);
    fenv::clear_flags();

    // Constant evaluation uses round-to-nearest-even and records nothing.
    constexpr float256 c3 = float256::one() / float256(3.0);
    CHECK(c3 == q[0]);

    std::cout << "passed=" << g_pass << " failed=" << g_fails << "\n";
    return g_fails ? 1 : 0;
}
