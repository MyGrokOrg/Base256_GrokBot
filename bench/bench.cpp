// Micro-benchmark: average time per call for basic operations and elementary functions.
// Usage: base256_bench [min_seconds_per_item=0.2] [filter]
#include <base256/base256.hpp>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <string>
#include <vector>

using base256::float256;

namespace {
volatile unsigned g_sink;
void sink(float256 x) { g_sink = g_sink + static_cast<unsigned>(x.signbit()) + static_cast<unsigned>(x.isnan()); }

struct Item {
    const char* name;
    std::function<float256(float256, float256)> f;
    double lo, hi;  // argument range
};

std::vector<float256> args(double lo, double hi, int n, unsigned seed) {
    std::vector<float256> v;
    unsigned s = seed;
    for (int i = 0; i < n; ++i) {
        s = s * 1664525u + 1013904223u;
        const float256 t = float256(static_cast<double>(s >> 8) / 16777216.0);
        v.push_back(float256(lo) + (float256(hi) - float256(lo)) * (t + float256(1e-3) * float256(i)) / float256(1.0 + 1e-3 * n));
    }
    return v;
}
}  // namespace

int main(int argc, char** argv) {
    using namespace base256;
    const double min_s = argc > 1 ? std::atof(argv[1]) : 0.2;
    const std::string filter = argc > 2 ? argv[2] : "";
    std::vector<Item> items = {
        {"add", [](float256 x, float256 y) { return x + y; }, 0.1, 100},
        {"mul", [](float256 x, float256 y) { return x * y; }, 0.1, 100},
        {"div", [](float256 x, float256 y) { return x / y; }, 0.1, 100},
        {"fma", [](float256 x, float256 y) { return fma(x, y, x); }, 0.1, 100},
        {"sqrt", [](float256 x, float256) { return sqrt(x); }, 0.1, 100},
        {"exp", [](float256 x, float256) { return exp(x); }, -50, 50},
        {"exp2", [](float256 x, float256) { return exp2(x); }, -50, 50},
        {"expm1", [](float256 x, float256) { return expm1(x); }, -2, 2},
        {"log", [](float256 x, float256) { return log(x); }, 0.01, 1e6},
        {"log2", [](float256 x, float256) { return log2(x); }, 0.01, 1e6},
        {"log10", [](float256 x, float256) { return log10(x); }, 0.01, 1e6},
        {"log1p", [](float256 x, float256) { return log1p(x); }, -0.5, 10},
        {"pow", [](float256 x, float256 y) { return pow(x, y); }, 0.1, 30},
        {"cbrt", [](float256 x, float256) { return cbrt(x); }, 0.1, 1e6},
        {"sin", [](float256 x, float256) { return sin(x); }, -10, 10},
        {"cos", [](float256 x, float256) { return cos(x); }, -10, 10},
        {"tan", [](float256 x, float256) { return tan(x); }, -10, 10},
        {"sin_1e30", [](float256 x, float256) { return sin(x * float256(1e30)); }, 1, 10},
        {"sin_1e9999", [](float256 x, float256) { return sin(ldexp(x, 33000)); }, 1, 2},
        {"atan", [](float256 x, float256) { return atan(x); }, -10, 10},
        {"asin", [](float256 x, float256) { return asin(x); }, -1, 1},
        {"acos", [](float256 x, float256) { return acos(x); }, -1, 1},
        {"atan2", [](float256 x, float256 y) { return atan2(x, y - float256(50)); }, -10, 10},
        {"sinh", [](float256 x, float256) { return sinh(x); }, -10, 10},
        {"cosh", [](float256 x, float256) { return cosh(x); }, -10, 10},
        {"tanh", [](float256 x, float256) { return tanh(x); }, -10, 10},
        {"asinh", [](float256 x, float256) { return asinh(x); }, -10, 10},
        {"acosh", [](float256 x, float256) { return acosh(x); }, 1, 100},
        {"atanh", [](float256 x, float256) { return atanh(x); }, -0.99, 0.99},
        {"erf", [](float256 x, float256) { return erf(x); }, -3, 3},
        {"erfc", [](float256 x, float256) { return erfc(x); }, -1, 20},
        {"tgamma", [](float256 x, float256) { return tgamma(x); }, 0.1, 50},
        {"lgamma", [](float256 x, float256) { return lgamma(x); }, 0.1, 1000},
        // Special functions (special_math.hpp): order/parameter from x, argument from y.
        {"cyl_bessel_j", [](float256 x, float256 y) { return cyl_bessel_j(x * float256(0.2), y); }, 0.1, 50},
        {"cyl_neumann", [](float256 x, float256 y) { return cyl_neumann(x * float256(0.2), y); }, 0.1, 50},
        {"cyl_bessel_i", [](float256 x, float256 y) { return cyl_bessel_i(x * float256(0.2), y); }, 0.1, 50},
        {"cyl_bessel_k", [](float256 x, float256 y) { return cyl_bessel_k(x * float256(0.2), y); }, 0.1, 50},
        {"jn", [](float256 x, float256 y) { return jn(static_cast<int>(x.to_double() * 0.2), y); }, 0.1, 50},
        {"sph_bessel", [](float256 x, float256 y) { return sph_bessel(static_cast<unsigned>(x.to_double() * 0.2), y); }, 0.1, 50},
        {"beta", [](float256 x, float256 y) { return beta(x, y); }, 0.1, 50},
        {"ibeta", [](float256 x, float256 y) { return ibeta(x, y, float256(0.4)); }, 0.1, 50},
        {"riemann_zeta", [](float256 x, float256) { return riemann_zeta(x); }, -20, 40},
        {"gamma_p", [](float256 x, float256 y) { return gamma_p(x, y); }, 0.1, 50},
        {"gamma_q", [](float256 x, float256 y) { return gamma_q(x, y); }, 0.1, 50},
        {"comp_ellint_1", [](float256 x, float256) { return comp_ellint_1(x); }, -0.99, 0.99},
        {"comp_ellint_2", [](float256 x, float256) { return comp_ellint_2(x); }, -0.99, 0.99},
        {"comp_ellint_3", [](float256 x, float256 y) { return comp_ellint_3(x, y * float256(0.9)); }, -0.99, 0.99},
        {"ellint_1", [](float256 x, float256 y) { return ellint_1(x, y * float256(3)); }, -0.99, 0.99},
        {"ellint_2", [](float256 x, float256 y) { return ellint_2(x, y * float256(3)); }, -0.99, 0.99},
        {"ellint_3", [](float256 x, float256 y) { return ellint_3(x, y * float256(0.5), y * float256(1.5)); }, -0.99, 0.99},
    };
    std::printf("%-12s %14s\n", "function", "ns/call");
    for (const Item& it : items) {
        if (!filter.empty() && filter.find(it.name) == std::string::npos) continue;
        const auto xs = args(it.lo, it.hi, 64, 1), ys = args(it.lo, it.hi, 64, 2);
        long calls = 0;
        const auto t0 = std::chrono::steady_clock::now();
        double el = 0;
        do {
            for (std::size_t i = 0; i < xs.size(); ++i) sink(it.f(xs[i], ys[i]));
            calls += static_cast<long>(xs.size());
            el = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        } while (el < min_s);
        std::printf("%-12s %14.1f\n", it.name, el / static_cast<double>(calls) * 1e9);
    }
    return 0;
}
