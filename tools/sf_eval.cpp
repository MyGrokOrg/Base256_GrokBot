// Reads "<func> <a> [<b> [<c>]]" lines (decimal or hex float256 literals) and
// prints the binary256 result as a hex float. Used by tools/sf_check.py.
// With SF_MODES set in the environment, each line is evaluated in all five
// rounding modes (RNE, RNA, RZ, RU, RD) and printed as "hex flags" pairs.
// With SF_PERCALL set, every function is called through its per-call rounding
// overload f(args..., m) while the thread's mode is set to a different (decoy)
// mode; a line whose call changes the thread's mode prints "?MODE".
#include <base256/base256.hpp>
#include <cstdio>
#include <cstdlib>
#include <format>
#include <iostream>
#include <sstream>
#include <string>
using namespace base256;
static bool g_pc = false;                         // per-call overloads (SF_PERCALL)
static round_mode g_m = round_mode::to_nearest_even;  // their mode
static float256 P(const std::string& s) { return (g_pc ? float256::parse(s, g_m) : float256::parse(s)).value(); }
#define C(fn, ...) (g_pc ? fn(__VA_ARGS__, g_m) : fn(__VA_ARGS__))
static bool eval(const std::string& f, const std::string& a, const std::string& b, const std::string& c, float256& r) {
        if (f == "j") r = C(cyl_bessel_j, P(a), P(b));
        else if (f == "y") r = C(cyl_neumann, P(a), P(b));
        else if (f == "i") r = C(cyl_bessel_i, P(a), P(b));
        else if (f == "k") r = C(cyl_bessel_k, P(a), P(b));
        else if (f == "sj") r = C(sph_bessel, static_cast<unsigned>(std::stoul(a)), P(b));
        else if (f == "sy") r = C(sph_neumann, static_cast<unsigned>(std::stoul(a)), P(b));
        else if (f == "jn") r = C(jn, std::stoi(a), P(b));
        else if (f == "yn") r = C(yn, std::stoi(a), P(b));
        else if (f == "beta") r = C(beta, P(a), P(b));
        else if (f == "ibeta") r = C(ibeta, P(a), P(b), P(c));
        else if (f == "ibetac") r = C(ibetac, P(a), P(b), P(c));
        else if (f == "zeta") r = C(riemann_zeta, P(a));
        else if (f == "hz") r = C(hurwitz_zeta, P(a), P(b));
        else if (f == "gl") r = C(tgamma_lower, P(a), P(b));
        else if (f == "gu") r = C(tgamma_upper, P(a), P(b));
        else if (f == "gp") r = C(gamma_p, P(a), P(b));
        else if (f == "gq") r = C(gamma_q, P(a), P(b));
        else if (f == "K") r = C(comp_ellint_1, P(a));
        else if (f == "E") r = C(comp_ellint_2, P(a));
        else if (f == "Pi") r = C(comp_ellint_3, P(a), P(b));
        else if (f == "F1") r = C(ellint_1, P(a), P(b));
        else if (f == "E2") r = C(ellint_2, P(a), P(b));
        else if (f == "P3") r = C(ellint_3, P(a), P(b), P(c));
        else if (f == "leg") r = C(legendre, static_cast<unsigned>(std::stoul(a)), P(b));
        else if (f == "aleg") r = C(assoc_legendre, static_cast<unsigned>(std::stoul(a)), static_cast<unsigned>(std::stoul(b)), P(c));
        else if (f == "sleg") r = C(sph_legendre, static_cast<unsigned>(std::stoul(a)), static_cast<unsigned>(std::stoul(b)), P(c));
        else if (f == "lag") r = C(laguerre, static_cast<unsigned>(std::stoul(a)), P(b));
        else if (f == "alag") r = C(assoc_laguerre, static_cast<unsigned>(std::stoul(a)), static_cast<unsigned>(std::stoul(b)), P(c));
        else if (f == "her") r = C(hermite, static_cast<unsigned>(std::stoul(a)), P(b));
        else if (f == "ei") r = C(expint, P(a));
        else if (f == "en") r = C(expint, static_cast<unsigned>(std::stoul(a)), P(b));
        else if (f == "si") r = C(sinint, P(a));
        else if (f == "ci") r = C(cosint, P(a));
        else if (f == "shi") r = C(sinhint, P(a));
        else if (f == "chi") r = C(coshint, P(a));
        else if (f == "li") r = C(logint, P(a));
        else if (f == "psi") r = C(digamma, P(a));
        else if (f == "pg") r = C(polygamma, static_cast<unsigned>(std::stoul(a)), P(b));
        else if (f == "lbeta") r = C(lbeta, P(a), P(b));
        else if (f == "gr") r = C(tgamma_ratio, P(a), P(b));
        else if (f == "gdr") r = C(tgamma_delta_ratio, P(a), P(b));
        else if (f == "ai") r = C(airy_ai, P(a));
        else if (f == "aip") r = C(airy_ai_prime, P(a));
        else if (f == "bi") r = C(airy_bi, P(a));
        else if (f == "bip") r = C(airy_bi_prime, P(a));
        else if (f == "erfinv") r = C(erf_inv, P(a));
        else if (f == "erfcinv") r = C(erfc_inv, P(a));
        else if (f == "w0") r = C(lambert_w0, P(a));
        else if (f == "wm1") r = C(lambert_wm1, P(a));
        else if (f == "fs") r = C(fresnel_s, P(a));
        else if (f == "fc") r = C(fresnel_c, P(a));
        else if (f == "daw") r = C(dawson, P(a));
        else if (f == "lg") r = C(lgamma, P(a));
        else if (f == "li2") r = C(dilog, P(a));
        else if (f == "plog") r = C(polylog, static_cast<unsigned>(std::stoul(a)), P(b));
        else if (f == "eta") r = C(dirichlet_eta, P(a));
        else if (f == "bern") r = C(bernoulli_b2n, static_cast<unsigned>(std::stoul(a)));
#ifdef SF_MORE
        SF_MORE
#endif
        else return false;
        return true;
}
static const round_mode kModes[5] = {round_mode::to_nearest_even, round_mode::to_nearest_away,
                                     round_mode::toward_zero, round_mode::upward, round_mode::downward};
static round_mode decoy(round_mode m) { return kModes[(static_cast<int>(m) + 2) % 5]; }
int main() {
    const bool modes = std::getenv("SF_MODES") != nullptr;
    g_pc = std::getenv("SF_PERCALL") != nullptr;
    std::string line;
    while (std::getline(std::cin, line)) {
        std::istringstream is(line);
        std::string f, a, b, c;
        is >> f >> a >> b >> c;
        float256 r;
        if (!modes) {
            g_m = round_mode::to_nearest_even;
            if (g_pc) fenv::set_round(decoy(g_m));
            const bool ok = eval(f, a, b, c, r);
            const bool kept = !g_pc || fenv::get_round() == decoy(g_m);
            fenv::set_round(round_mode::to_nearest_even);
            if (!ok) { std::printf("?\n"); continue; }
            if (!kept) { std::printf("?MODE\n"); continue; }
            std::cout << std::format("{:a}", r) << "\n";
            continue;
        }
        std::string out;
        bool ok = true, kept = true;
        for (round_mode m : kModes) {
            g_m = m;
            fenv::set_round(g_pc ? decoy(m) : m);
            fenv::clear_flags();
            ok = eval(f, a, b, c, r);
            const unsigned fl = fenv::test_flags();
            kept = kept && (!g_pc || fenv::get_round() == decoy(m));
            fenv::set_round(round_mode::to_nearest_even);
            if (!ok) break;
            out += std::format("{:a} {} ", r, fl);
        }
        std::cout << (!ok ? std::string("?") : !kept ? std::string("?MODE") : out) << "\n";
    }
}
