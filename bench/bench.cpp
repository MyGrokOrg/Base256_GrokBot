// Micro-benchmark: average time per call for basic operations and elementary functions.
// Usage: base256_bench [min_seconds_per_item=0.2] [filter]
#include <base256/base256.hpp>
#include <base256/complex.hpp>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <string>
#include <vector>

using base256::float256;

namespace {
volatile unsigned g_sink;
float256 cx_out(const base256::complex& z) { return z.real() + z.imag(); }
float256 out_len(const std::string& s) { return float256(static_cast<double>(s.size() + static_cast<unsigned char>(s[s.size() / 2]))); }
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
        {"hurwitz_zeta", [](float256 x, float256 y) { return hurwitz_zeta(x + float256(1), y); }, 0.1, 50},
        {"hurwitz_zeta(3,a)", [](float256, float256 y) { return hurwitz_zeta(float256(3), y); }, 0.1, 50},
        {"hurwitz_zeta(s<1)", [](float256 x, float256 y) { return hurwitz_zeta(x, y + float256(20.5)); }, -20, 0.9},
        {"hurwitz_zeta(2..4,a<1)", [](float256 x, float256 y) { return hurwitz_zeta(x, y * float256(0.2)); }, 2, 4},
        {"legendre(20,x)", [](float256 x, float256) { return legendre(20, x); }, -1, 1},
        {"hermite(30,x)", [](float256 x, float256) { return hermite(30, x); }, -10, 10},
        {"laguerre(30,x)", [](float256 x, float256) { return laguerre(30, x); }, 0, 50},
        {"sph_legendre(10,3,t)", [](float256 x, float256) { return sph_legendre(10, 3, x); }, 0, 3},
        {"expint (Ei)", [](float256 x, float256) { return expint(x); }, -50, 50},
        {"expint(3,x) (E_n)", [](float256 x, float256) { return expint(3, x); }, 0.1, 50},
        {"sinint", [](float256 x, float256) { return sinint(x); }, 0.1, 50},
        {"cosint", [](float256 x, float256) { return cosint(x); }, 0.1, 50},
        {"sinhint", [](float256 x, float256) { return sinhint(x); }, 0.1, 50},
        {"coshint", [](float256 x, float256) { return coshint(x); }, 0.1, 50},
        {"logint", [](float256 x, float256) { return logint(x); }, 0.1, 50},
        {"digamma", [](float256 x, float256) { return digamma(x); }, -20, 50},
        {"trigamma", [](float256 x, float256) { return trigamma(x); }, 0.1, 50},
        {"polygamma(3,x)", [](float256 x, float256) { return polygamma(3, x); }, 0.1, 50},
        {"polygamma(1,x<0)", [](float256 x, float256) { return polygamma(1, x); }, -50, -0.1},
        {"polygamma(2,x<0)", [](float256 x, float256) { return polygamma(2, x); }, -50, -0.1},
        {"polygamma(3,x<0)", [](float256 x, float256) { return polygamma(3, x); }, -50, -0.1},
        {"polygamma(10,x<0)", [](float256 x, float256) { return polygamma(10, x); }, -20, -0.1},
        {"lbeta", [](float256 x, float256 y) { return lbeta(x, y); }, 0.1, 50},
        {"tgamma_ratio", [](float256 x, float256 y) { return tgamma_ratio(x, y); }, 0.1, 50},
        {"airy_ai", [](float256 x, float256) { return airy_ai(x); }, -50, 50},
        {"airy_bi", [](float256 x, float256) { return airy_bi(x); }, -50, 50},
        {"airy_ai_prime", [](float256 x, float256) { return airy_ai_prime(x); }, -50, 50},
        {"lambert_w0", [](float256 x, float256) { return lambert_w0(x); }, -0.36, 50},
        {"lambert_wm1", [](float256 x, float256) { return lambert_wm1(x); }, -0.36, -0.001},
        {"erf_inv", [](float256 x, float256) { return erf_inv(x); }, -0.999, 0.999},
        {"erfc_inv", [](float256 x, float256) { return erfc_inv(x); }, 0.001, 1.999},
        {"dawson", [](float256 x, float256) { return dawson(x); }, -50, 50},
        {"fresnel_s", [](float256 x, float256) { return fresnel_s(x); }, -50, 50},
        {"fresnel_c", [](float256 x, float256) { return fresnel_c(x); }, -50, 50},
        {"dilog", [](float256 x, float256) { return dilog(x); }, -50, 1},
        {"polylog(3,x)", [](float256 x, float256) { return polylog(3, x); }, -1, 1},
        {"dirichlet_eta", [](float256 x, float256) { return dirichlet_eta(x); }, -20, 40},
        {"gamma_p", [](float256 x, float256 y) { return gamma_p(x, y); }, 0.1, 50},
        {"gamma_q", [](float256 x, float256 y) { return gamma_q(x, y); }, 0.1, 50},
        {"comp_ellint_1", [](float256 x, float256) { return comp_ellint_1(x); }, -0.99, 0.99},
        {"comp_ellint_2", [](float256 x, float256) { return comp_ellint_2(x); }, -0.99, 0.99},
        {"comp_ellint_3", [](float256 x, float256 y) { return comp_ellint_3(x, y * float256(0.9)); }, -0.99, 0.99},
        {"ellint_1", [](float256 x, float256 y) { return ellint_1(x, y * float256(3)); }, -0.99, 0.99},
        {"ellint_2", [](float256 x, float256 y) { return ellint_2(x, y * float256(3)); }, -0.99, 0.99},
        {"ellint_3", [](float256 x, float256 y) { return ellint_3(x, y * float256(0.5), y * float256(1.5)); }, -0.99, 0.99},
        // base256::complex (complex.hpp): z = x + iy, both parts from the range.
        {"cx mul", [](float256 x, float256 y) { return cx_out(complex(x, y) * complex(y, x)); }, -10, 10},
        {"cx div", [](float256 x, float256 y) { return cx_out(complex(x, y) / complex(y, x)); }, -10, 10},
        {"cx abs", [](float256 x, float256 y) { return abs(complex(x, y)); }, -10, 10},
        {"cx sqrt", [](float256 x, float256 y) { return cx_out(sqrt(complex(x, y))); }, -10, 10},
        {"cx exp", [](float256 x, float256 y) { return cx_out(exp(complex(x, y))); }, -10, 10},
        {"cx log", [](float256 x, float256 y) { return cx_out(log(complex(x, y))); }, -10, 10},
        {"cx sin", [](float256 x, float256 y) { return cx_out(sin(complex(x, y))); }, -10, 10},
        {"cx tan", [](float256 x, float256 y) { return cx_out(tan(complex(x, y))); }, -10, 10},
        {"cx asin", [](float256 x, float256 y) { return cx_out(asin(complex(x, y))); }, -10, 10},
        {"cx atanh", [](float256 x, float256 y) { return cx_out(atanh(complex(x, y))); }, -10, 10},
        {"cx pow", [](float256 x, float256 y) { return cx_out(pow(complex(x, y), complex(y * float256(0.3), x * float256(0.1)))); }, -10, 10},
        {"cx pow(z,5)", [](float256 x, float256 y) { return cx_out(pow(complex(x, y), float256(5))); }, -10, 10},
        // Parts 2^400000 apart (x, y in [1, 2)): the exact operations work on ~400000-bit BigInts.
        {"cx mul (gap)", [](float256 x, float256 y) { return cx_out(complex(ldexp(x, 200000), ldexp(y, -200000)) * complex(ldexp(y, -200000), ldexp(x, 200000))); }, 1, 2},
        {"cx div (gap)", [](float256 x, float256 y) { return cx_out(complex(ldexp(x, 200000), ldexp(y, -200000)) / complex(ldexp(y, -200000), ldexp(x, 200000))); }, 1, 2},
        {"cx abs (gap)", [](float256 x, float256 y) { return abs(complex(ldexp(x, 200000), ldexp(y, -200000))); }, 1, 2},
        {"cx sqrt (gap)", [](float256 x, float256 y) { return cx_out(sqrt(complex(ldexp(x, 200000), ldexp(y, -200000)))); }, 1, 2},
        // Decimal output at the extreme exponents (x in [1, 2) scaled to near max / deep subnormal).
        {"to_string(~2^262140)", [](float256 x, float256) { return out_len(ldexp(x, 262140).to_string()); }, 1, 2},
        {"to_string(~2^-262370)", [](float256 x, float256) { return out_len(ldexp(x, -262370).to_string()); }, 1, 2},
        {"to_string(~2^-262000)", [](float256 x, float256) { return out_len(ldexp(x, -262000).to_string()); }, 1, 2},
        {"format(e,75)(~2^262140)", [](float256 x, float256) { return out_len(ldexp(x, 262140).format(std::chars_format::scientific, 75)); }, 1, 2},
        {"format(e,75)(~2^-262370)", [](float256 x, float256) { return out_len(ldexp(x, -262370).format(std::chars_format::scientific, 75)); }, 1, 2},
        {"format(f,0)(~2^262140)", [](float256 x, float256) { return out_len(ldexp(x, 262140).format(std::chars_format::fixed, 0)); }, 1, 2},
        {"to_chars(~2^262140)", [](float256 x, float256) { char b[100]; const auto r = to_chars(b, b + 100, ldexp(x, 262140)); return float256(static_cast<double>(r.ptr - b)); }, 1, 2},
        {"to_string(1..2)", [](float256 x, float256) { return out_len(x.to_string()); }, 1, 2},
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
