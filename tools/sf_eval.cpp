// Reads "<func> <a> [<b> [<c>]]" lines (decimal or hex float256 literals) and
// prints the binary256 result as a hex float. Used by tools/sf_check.py.
// With SF_MODES set in the environment, each line is evaluated in all five
// rounding modes (RNE, RNA, RZ, RU, RD) and printed as "hex flags" pairs.
#include <base256/base256.hpp>
#include <cstdio>
#include <cstdlib>
#include <format>
#include <iostream>
#include <sstream>
#include <string>
using namespace base256;
static float256 P(const std::string& s) { return float256::parse(s).value(); }
static bool eval(const std::string& f, const std::string& a, const std::string& b, const std::string& c, float256& r) {
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
        else if (f == "hz") r = hurwitz_zeta(P(a), P(b));
        else if (f == "gl") r = tgamma_lower(P(a), P(b));
        else if (f == "gu") r = tgamma_upper(P(a), P(b));
        else if (f == "gp") r = gamma_p(P(a), P(b));
        else if (f == "gq") r = gamma_q(P(a), P(b));
        else if (f == "K") r = comp_ellint_1(P(a));
        else if (f == "E") r = comp_ellint_2(P(a));
        else if (f == "Pi") r = comp_ellint_3(P(a), P(b));
        else if (f == "F1") r = ellint_1(P(a), P(b));
        else if (f == "E2") r = ellint_2(P(a), P(b));
        else if (f == "P3") r = ellint_3(P(a), P(b), P(c));
        else if (f == "leg") r = legendre(static_cast<unsigned>(std::stoul(a)), P(b));
        else if (f == "aleg") r = assoc_legendre(static_cast<unsigned>(std::stoul(a)), static_cast<unsigned>(std::stoul(b)), P(c));
        else if (f == "sleg") r = sph_legendre(static_cast<unsigned>(std::stoul(a)), static_cast<unsigned>(std::stoul(b)), P(c));
        else if (f == "lag") r = laguerre(static_cast<unsigned>(std::stoul(a)), P(b));
        else if (f == "alag") r = assoc_laguerre(static_cast<unsigned>(std::stoul(a)), static_cast<unsigned>(std::stoul(b)), P(c));
        else if (f == "her") r = hermite(static_cast<unsigned>(std::stoul(a)), P(b));
        else if (f == "ei") r = expint(P(a));
        else if (f == "en") r = expint(static_cast<unsigned>(std::stoul(a)), P(b));
        else if (f == "si") r = sinint(P(a));
        else if (f == "ci") r = cosint(P(a));
        else if (f == "shi") r = sinhint(P(a));
        else if (f == "chi") r = coshint(P(a));
        else if (f == "li") r = logint(P(a));
        else if (f == "psi") r = digamma(P(a));
        else if (f == "pg") r = polygamma(static_cast<unsigned>(std::stoul(a)), P(b));
        else if (f == "lbeta") r = lbeta(P(a), P(b));
        else if (f == "gr") r = tgamma_ratio(P(a), P(b));
        else if (f == "gdr") r = tgamma_delta_ratio(P(a), P(b));
        else if (f == "ai") r = airy_ai(P(a));
        else if (f == "aip") r = airy_ai_prime(P(a));
        else if (f == "bi") r = airy_bi(P(a));
        else if (f == "bip") r = airy_bi_prime(P(a));
        else if (f == "erfinv") r = erf_inv(P(a));
        else if (f == "erfcinv") r = erfc_inv(P(a));
        else if (f == "w0") r = lambert_w0(P(a));
        else if (f == "wm1") r = lambert_wm1(P(a));
        else if (f == "fs") r = fresnel_s(P(a));
        else if (f == "fc") r = fresnel_c(P(a));
        else if (f == "daw") r = dawson(P(a));
        else if (f == "li2") r = dilog(P(a));
        else if (f == "plog") r = polylog(static_cast<unsigned>(std::stoul(a)), P(b));
        else if (f == "eta") r = dirichlet_eta(P(a));
        else if (f == "bern") r = bernoulli_b2n(static_cast<unsigned>(std::stoul(a)));
#ifdef SF_MORE
        SF_MORE
#endif
        else return false;
        return true;
}
int main() {
    const bool modes = std::getenv("SF_MODES") != nullptr;
    std::string line;
    while (std::getline(std::cin, line)) {
        std::istringstream is(line);
        std::string f, a, b, c;
        is >> f >> a >> b >> c;
        float256 r;
        if (!modes) {
            if (!eval(f, a, b, c, r)) { std::printf("?\n"); continue; }
            std::cout << std::format("{:a}", r) << "\n";
            continue;
        }
        std::string out;
        bool ok = true;
        for (round_mode m : {round_mode::to_nearest_even, round_mode::to_nearest_away, round_mode::toward_zero,
                             round_mode::upward, round_mode::downward}) {
            fenv::set_round(m);
            fenv::clear_flags();
            ok = eval(f, a, b, c, r);
            const unsigned fl = fenv::test_flags();
            fenv::set_round(round_mode::to_nearest_even);
            if (!ok) break;
            out += std::format("{:a} {} ", r, fl);
        }
        std::cout << (ok ? out : std::string("?")) << "\n";
    }
}
