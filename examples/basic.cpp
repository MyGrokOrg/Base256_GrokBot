#include <base256/base256.hpp>

#include <iostream>

int main() {
    using namespace base256;
    using namespace base256::literals;

    const auto pi = float256::pi();
    const auto two = 2.0_f256;
    const auto area = pi * two * two;

    std::cout << "pi       = " << pi.to_string(40) << "\n";
    std::cout << "pi hex   = " << pi.to_hex_string() << "\n";
    std::cout << "area 4pi = " << area.to_string(40) << "\n";
    std::cout << "sqrt(2)  = " << sqrt(two).to_string(40) << "\n";
    return 0;
}
