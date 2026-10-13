#include <base256/base256.hpp>

#include <cmath>
#include <iomanip>
#include <iostream>

int main() {
    using namespace base256;

    // 1 + 2^-80 is lost in IEEE-754 binary64, kept in binary256.
    const float256 one = float256::one();
    const float256 delta = ldexp(float256::one(), -80);
    const float256 sum = one + delta;

    std::cout << std::boolalpha;
    std::cout << "binary64:  1 + 2^-80 == 1  -> " << (1.0 + std::ldexp(1.0, -80) == 1.0) << "\n";
    std::cout << "binary256: 1 + 2^-80 == 1  -> " << (sum == one) << "\n";
    std::cout << "binary256: (1 + 2^-80) - 1 = " << (sum - one).to_hex_string() << "\n";

    const auto a = *float256::parse("0.1");
    const auto b = *float256::parse("0.2");
    const auto c = *float256::parse("0.3");
    std::cout << "\n0.1 + 0.2  (binary256) = " << (a + b).to_string(40) << "\n";
    std::cout << "0.3          (binary256) = " << c.to_string(40) << "\n";
    std::cout << "| (0.1+0.2) - 0.3 |      = " << abs((a + b) - c).to_string(10) << "\n";
    std::cout << "| (0.1+0.2) - 0.3 | f64  = " << std::fabs((0.1 + 0.2) - 0.3) << "\n";

    const float256 third = float256::one() / float256(3);
    std::cout << "\n1/3          = " << third.to_string(50) << "\n";
    std::cout << "(1/3)*3 == 1 = " << (third * float256(3) == float256::one()) << "\n";

    const float256 g = tgamma(float256::half());
    std::cout << "\ngamma(1/2) = " << g.to_string(40) << "\n";
    std::cout << "sqrt(pi)    = " << sqrt(float256::pi()).to_string(40) << "\n";
    return 0;
}
