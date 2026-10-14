// Reads "<func> <a> [<b> [<c>]]" lines (decimal or hex float256 literals) and
// prints the binary256 result as a hex float. Used by tools/sf_check.py.
#include <base256/base256.hpp>
#include <cstdio>
#include <format>
#include <iostream>
#include <sstream>
#include <string>
using namespace base256;
static float256 P(const std::string& s) { return float256::parse(s).value(); }
int main() {
    std::string line;
    while (std::getline(std::cin, line)) {
        std::istringstream is(line);
        std::string f, a, b, c;
        is >> f >> a >> b >> c;
        float256 r;
        if (f == "j") r = cyl_bessel_j(P(a), P(b));
        else if (f == "y") r = cyl_neumann(P(a), P(b));
        else if (f == "i") r = cyl_bessel_i(P(a), P(b));
        else if (f == "k") r = cyl_bessel_k(P(a), P(b));
        else if (f == "sj") r = sph_bessel(static_cast<unsigned>(std::stoul(a)), P(b));
        else if (f == "sy") r = sph_neumann(static_cast<unsigned>(std::stoul(a)), P(b));
        else if (f == "jn") r = jn(std::stoi(a), P(b));
        else if (f == "yn") r = yn(std::stoi(a), P(b));
        else if (f == "beta") r = beta(P(a), P(b));
        else if (f == "ibeta") r = ibeta(P(a), P(b), P(c));
        else if (f == "ibetac") r = ibetac(P(a), P(b), P(c));
        else if (f == "zeta") r = riemann_zeta(P(a));
        else if (f == "gl") r = tgamma_lower(P(a), P(b));
        else if (f == "gu") r = tgamma_upper(P(a), P(b));
        else if (f == "gp") r = gamma_p(P(a), P(b));
        else if (f == "gq") r = gamma_q(P(a), P(b));
#ifdef SF_MORE
        SF_MORE
#endif
        else { std::printf("?\n"); continue; }
        std::cout << std::format("{:a}", r) << "\n";
    }
}
