// Reads "<op> <args...>" lines (float256 literals, usually hex) and prints the
// base256::complex result as hex floats. Used by tools/cx_check.py.
//   unary complex ops   exp log log10 sqrt sin cos tan asin acos atan
//                       sinh cosh tanh asinh acosh atanh       args: re im
//   real-valued         abs norm arg                           args: re im
//   binary              mul div pow                            args: re im re im
//   powr (z^t, real t)  args: re im t      polar   args: rho theta
// Output: "re im" (or one value for abs/norm/arg). With CX_MODES set, every line
// is evaluated in all five rounding modes (RNE RNA RZ RU RD) and printed as
// "re im flags" groups (or "v flags"). With CX_PERCALL set, every operation goes
// through its per-call rounding overload (mul(z, w, m), exp(z, m), ...) while the
// thread's mode is set to a different (decoy) mode; a line whose call changes the
// thread's mode prints "?MODE".
#include <base256/complex.hpp>
#include <cstdio>
#include <cstdlib>
#include <format>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
using namespace base256;
static bool g_pc = false;                         // per-call overloads (CX_PERCALL)
static round_mode g_m = round_mode::to_nearest_even;  // their mode
static float256 P(const std::string& s) { return (g_pc ? float256::parse(s, g_m) : float256::parse(s)).value(); }
#define C(fn, ...) (g_pc ? fn(__VA_ARGS__, g_m) : fn(__VA_ARGS__))

static int eval(const std::string& f, const std::vector<std::string>& a, complex& r, float256& v) {
    auto z = [&](int i) { return complex(P(a[i]), P(a[i + 1])); };
    if (f == "abs") { v = C(abs, z(0)); return 1; }
    if (f == "norm") { v = C(norm, z(0)); return 1; }
    if (f == "arg") { v = C(arg, z(0)); return 1; }
    if (f == "mul") r = g_pc ? mul(z(0), z(2), g_m) : z(0) * z(2);
    else if (f == "div") r = g_pc ? div(z(0), z(2), g_m) : z(0) / z(2);
    else if (f == "pow") r = C(pow, z(0), z(2));
    else if (f == "powr") r = C(pow, z(0), P(a[2]));
    else if (f == "polar") r = C(polar, P(a[0]), P(a[1]));
    else if (f == "exp") r = C(exp, z(0));
    else if (f == "log") r = C(log, z(0));
    else if (f == "log10") r = C(log10, z(0));
    else if (f == "sqrt") r = C(sqrt, z(0));
    else if (f == "sin") r = C(sin, z(0));
    else if (f == "cos") r = C(cos, z(0));
    else if (f == "tan") r = C(tan, z(0));
    else if (f == "asin") r = C(asin, z(0));
    else if (f == "acos") r = C(acos, z(0));
    else if (f == "atan") r = C(atan, z(0));
    else if (f == "sinh") r = C(sinh, z(0));
    else if (f == "cosh") r = C(cosh, z(0));
    else if (f == "tanh") r = C(tanh, z(0));
    else if (f == "asinh") r = C(asinh, z(0));
    else if (f == "acosh") r = C(acosh, z(0));
    else if (f == "atanh") r = C(atanh, z(0));
    else return 0;
    return 2;
}

static const round_mode kModes[5] = {round_mode::to_nearest_even, round_mode::to_nearest_away,
                                     round_mode::toward_zero, round_mode::upward, round_mode::downward};
static round_mode decoy(round_mode m) { return kModes[(static_cast<int>(m) + 2) % 5]; }
int main() {
    const bool modes = std::getenv("CX_MODES") != nullptr;
    g_pc = std::getenv("CX_PERCALL") != nullptr;
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
            g_m = round_mode::to_nearest_even;
            if (g_pc) fenv::set_round(decoy(g_m));
            const int k = eval(f, a, r, v);
            const bool kept = !g_pc || fenv::get_round() == decoy(g_m);
            fenv::set_round(round_mode::to_nearest_even);
            if (k == 0) std::printf("?\n");
            else if (!kept) std::printf("?MODE\n");
            else if (k == 1) std::cout << std::format("{:a}", v) << "\n";
            else std::cout << std::format("{:a} {:a}", r.real(), r.imag()) << "\n";
            continue;
        }
        std::string out;
        bool ok = true, kept = true;
        for (round_mode m : kModes) {
            g_m = m;
            fenv::set_round(g_pc ? decoy(m) : m);
            fenv::clear_flags();
            const int k = eval(f, a, r, v);
            const unsigned fl = fenv::test_flags();
            kept = kept && (!g_pc || fenv::get_round() == decoy(m));
            fenv::set_round(round_mode::to_nearest_even);
            if (k == 0) { ok = false; break; }
            if (k == 1) out += std::format("{:a} {} ", v, fl);
            else out += std::format("{:a} {:a} {} ", r.real(), r.imag(), fl);
        }
        std::cout << (!ok ? std::string("?") : !kept ? std::string("?MODE") : out) << "\n";
    }
}
