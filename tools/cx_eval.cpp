// Reads "<op> <args...>" lines (float256 literals, usually hex) and prints the
// base256::complex result as hex floats. Used by tools/cx_check.py.
//   unary complex ops   exp log log10 sqrt sin cos tan asin acos atan
//                       sinh cosh tanh asinh acosh atanh       args: re im
//   real-valued         abs norm arg                           args: re im
//   binary              mul div pow                            args: re im re im
//   powr (z^t, real t)  args: re im t      polar   args: rho theta
// Output: "re im" (or one value for abs/norm/arg). With CX_MODES set, every line
// is evaluated in all five rounding modes (RNE RNA RZ RU RD) and printed as
// "re im flags" groups (or "v flags").
#include <base256/complex.hpp>
#include <cstdio>
#include <cstdlib>
#include <format>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
using namespace base256;
static float256 P(const std::string& s) { return float256::parse(s).value(); }

static int eval(const std::string& f, const std::vector<std::string>& a, complex& r, float256& v) {
    auto z = [&](int i) { return complex(P(a[i]), P(a[i + 1])); };
    if (f == "abs") { v = abs(z(0)); return 1; }
    if (f == "norm") { v = norm(z(0)); return 1; }
    if (f == "arg") { v = arg(z(0)); return 1; }
    if (f == "mul") r = z(0) * z(2);
    else if (f == "div") r = z(0) / z(2);
    else if (f == "pow") r = pow(z(0), z(2));
    else if (f == "powr") r = pow(z(0), P(a[2]));
    else if (f == "polar") r = polar(P(a[0]), P(a[1]));
    else if (f == "exp") r = exp(z(0));
    else if (f == "log") r = log(z(0));
    else if (f == "log10") r = log10(z(0));
    else if (f == "sqrt") r = sqrt(z(0));
    else if (f == "sin") r = sin(z(0));
    else if (f == "cos") r = cos(z(0));
    else if (f == "tan") r = tan(z(0));
    else if (f == "asin") r = asin(z(0));
    else if (f == "acos") r = acos(z(0));
    else if (f == "atan") r = atan(z(0));
    else if (f == "sinh") r = sinh(z(0));
    else if (f == "cosh") r = cosh(z(0));
    else if (f == "tanh") r = tanh(z(0));
    else if (f == "asinh") r = asinh(z(0));
    else if (f == "acosh") r = acosh(z(0));
    else if (f == "atanh") r = atanh(z(0));
    else return 0;
    return 2;
}

int main() {
    const bool modes = std::getenv("CX_MODES") != nullptr;
    std::string line;
    while (std::getline(std::cin, line)) {
        std::istringstream is(line);
        std::string f, t;
        std::vector<std::string> a;
        is >> f;
        while (is >> t) a.push_back(t);
        complex r;
        float256 v;
        if (!modes) {
            const int k = eval(f, a, r, v);
            if (k == 0) std::printf("?\n");
            else if (k == 1) std::cout << std::format("{:a}", v) << "\n";
            else std::cout << std::format("{:a} {:a}", r.real(), r.imag()) << "\n";
            continue;
        }
        std::string out;
        bool ok = true;
        for (round_mode m : {round_mode::to_nearest_even, round_mode::to_nearest_away, round_mode::toward_zero,
                             round_mode::upward, round_mode::downward}) {
            fenv::set_round(m);
            fenv::clear_flags();
            const int k = eval(f, a, r, v);
            const unsigned fl = fenv::test_flags();
            fenv::set_round(round_mode::to_nearest_even);
            if (k == 0) { ok = false; break; }
            if (k == 1) out += std::format("{:a} {} ", v, fl);
            else out += std::format("{:a} {:a} {} ", r.real(), r.imag(), fl);
        }
        std::cout << (ok ? out : std::string("?")) << "\n";
    }
}
