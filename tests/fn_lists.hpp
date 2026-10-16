// X-macro lists of the public functions with per-call rounding overloads, by
// signature. Shared by tests/test_percall.cpp (per-call vs global mode) and
// tools/mode_digest.cpp (digest of the global-mode results, built against two
// revisions to show the default path is unchanged).
#pragma once

// float256 f(float256)
#define B256_FN_X(X) X(exp) X(exp2) X(expm1) X(log) X(log2) X(log10) X(log1p) X(cbrt) X(sin) X(cos) X(tan) \
    X(atan) X(asin) X(acos) X(sinh) X(cosh) X(tanh) X(erf) X(erfc) X(asinh) X(acosh) X(atanh) X(lgamma)   \
    X(tgamma) X(riemann_zeta) X(comp_ellint_1) X(comp_ellint_2) X(expint) X(sinint) X(cosint) X(sinhint)  \
    X(coshint) X(logint) X(digamma) X(trigamma) X(airy_ai) X(airy_ai_prime) X(airy_bi) X(airy_bi_prime)   \
    X(erf_inv) X(erfc_inv) X(lambert_w0) X(lambert_wm1) X(fresnel_s) X(fresnel_c) X(dawson) X(dilog)      \
    X(dirichlet_eta)
// float256 f(float256, float256)
#define B256_FN_XY(X) X(pow) X(atan2) X(fdim) X(hypot) X(cyl_bessel_j) X(cyl_neumann) X(cyl_bessel_i)      \
    X(cyl_bessel_k) X(beta) X(hurwitz_zeta) X(tgamma_lower) X(tgamma_upper) X(gamma_p) X(gamma_q)          \
    X(comp_ellint_3) X(ellint_1) X(ellint_2) X(lbeta) X(tgamma_ratio) X(tgamma_delta_ratio)
// float256 f(float256, float256, float256)
#define B256_FN_XYZ(X) X(ibeta) X(ibetac) X(ellint_3)
// float256 f(unsigned, float256)
#define B256_FN_NX(X) X(sph_bessel) X(sph_neumann) X(legendre) X(laguerre) X(hermite) X(expint) X(polygamma) \
    X(polylog)
// float256 f(int, float256)
#define B256_FN_IX(X) X(jn) X(yn)
// float256 f(unsigned, unsigned, float256)
#define B256_FN_NMX(X) X(assoc_legendre) X(sph_legendre) X(assoc_laguerre)
// complex f(complex)
#define B256_CX_Z(X) X(exp) X(log) X(log10) X(sqrt) X(sinh) X(cosh) X(tanh) X(sin) X(cos) X(tan) X(acos) \
    X(asin) X(atan) X(asinh) X(acosh) X(atanh)
// float256 f(complex)
#define B256_CX_R(X) X(abs) X(norm) X(arg)
