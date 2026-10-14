#!/usr/bin/env python3
"""Accuracy check of the special functions against mpmath (independent reference).

usage: sf_check.py EVAL_BINARY [cases_per_range] [seed] [filter]

Each range draws random binary256 arguments, evaluates them with EVAL_BINARY
(tools/sf_eval.cpp, round-to-nearest) and compares with mpmath at 700 bits.
Prints per range the max error in ulps and the number of results that are not
the correctly rounded value (error > 0.5 ulp).
"""
import random, subprocess, sys
import mpmath as mp

mp.mp.prec = 700
P = 237
EMIN = -262142

def rnd256(lo, hi, logscale, rng):
    """random binary256 number in [lo, hi] (or 2^[lo,hi] if logscale), as (hex, mpf)."""
    if logscale:
        e = rng.randint(int(lo), int(hi) - 1)
        m = rng.getrandbits(P - 1) | (1 << (P - 1))
        v = mp.ldexp(m, e - (P - 1))
    else:
        v = mp.mpf(lo) + (mp.mpf(hi) - mp.mpf(lo)) * mp.mpf(rng.getrandbits(P)) / (1 << P)
        if v != 0:
            e = int(mp.floor(mp.log(abs(v), 2)))
            m = int(mp.nint(mp.ldexp(v, P - 1 - e)))
            v = mp.ldexp(m, e - (P - 1))
    return to_hex(v), v

def to_hex(v):
    if v == 0:
        return "0x0p+0"
    s = "-" if v < 0 else ""
    v = abs(v)
    e = int(mp.floor(mp.log(v, 2)))
    m = int(mp.ldexp(v, P - 1 - e))
    if m >= (1 << P):
        e += 1; m = int(mp.ldexp(v, P - 1 - e))
    if m < (1 << (P - 1)):
        e -= 1; m = int(mp.ldexp(v, P - 1 - e))
    return f"{s}0x{m:x}p{e - (P - 1):+d}"

def parse_hex(s):
    s = s.strip().replace("0xinf", "inf")
    if s in ("inf", "+inf"): return mp.inf
    if s == "-inf": return -mp.inf
    if "nan" in s: return mp.nan
    neg = s.startswith("-")
    s = s.lstrip("+-")[2:]
    mant, ex = s.split("p")
    if "." in mant:
        a, b = mant.split(".")
    else:
        a, b = mant, ""
    v = mp.ldexp(int(a + b, 16), int(ex) - 4 * len(b))
    return -v if neg else v

def ulp_err(got, ref):
    if mp.isnan(got) or mp.isinf(got):
        return float("inf") if not (mp.isinf(ref) and got == ref) else 0.0
    if ref == 0:
        return 0.0 if got == 0 else float("inf")
    e = max(int(mp.floor(mp.log(abs(ref), 2))), EMIN)
    u = mp.ldexp(1, e - (P - 1))
    return float(abs(got - ref) / u)

# name, eval-op, reference, list of (arg ranges): each arg (lo, hi, logscale) or ("int", lo, hi)
def J(nu, x): return mp.besselj(nu, x)
def Y(nu, x): return mp.bessely(nu, x)
def I(nu, x): return mp.besseli(nu, x)
def K(nu, x): return mp.besselk(nu, x)
def SJ(n, x): return mp.sqrt(mp.pi / (2 * x)) * mp.besselj(n + mp.mpf(1) / 2, x)
def SY(n, x): return mp.sqrt(mp.pi / (2 * x)) * mp.bessely(n + mp.mpf(1) / 2, x)

def BETA(a, b): return mp.beta(a, b)
# mpmath's betainc can lose bits to cancellation; evaluate it at 2000 bits.
def IB(a, b, x):
    with mp.workprec(2000): return mp.betainc(a, b, 0, x, regularized=True)
def IBC(a, b, x):
    with mp.workprec(2000): return mp.betainc(a, b, x, 1, regularized=True)

def ZETA(x): return mp.zeta(x)
# mpmath's zeta(s, a) loses hundreds of bits at large a; evaluate it at 3000 bits.
def HZ(s, a):
    with mp.workprec(3000): return mp.zeta(s, a)

def GL(a, x): return mp.gammainc(a, 0, x)
def GU(a, x): return mp.gammainc(a, x, mp.inf)
def GP(a, x): return mp.gammainc(a, 0, x, regularized=True)
def GQ(a, x): return mp.gammainc(a, x, mp.inf, regularized=True)

def EK(k): return mp.ellipk(k * k)
def EE(k): return mp.ellipe(k * k)
def EPI(k, n): return mp.ellippi(n, k * k)
def EF(k, phi): return mp.ellipf(phi, k * k)
def EE2(k, phi): return mp.ellipe(phi, k * k)
def EP3(k, n, phi): return mp.ellippi(n, phi, k * k)

# Orthogonal polynomials: three-term recurrences at 2000 bits (mpmath's
# hypergeometric evaluation is far slower at large degree; the two agree, see
# --selftest).
def _rec(n, p0, p1, step):
    if n == 0: return p0
    for j in range(1, n):
        p0, p1 = p1, step(j, p0, p1)
    return p1
def LEG(n, x):
    n = int(n)
    with mp.workprec(2000): return _rec(n, mp.mpf(1), x, lambda j, a, b: ((2*j+1)*x*b - j*a) / (j+1))
def ALEG(n, m, x):
    n, m = int(n), int(m)
    if m > n: return mp.mpf(0)
    with mp.workprec(2000):
        pmm = mp.fac2(2*m-1) * mp.sqrt(1 - x*x) ** m if m else mp.mpf(1)
        if n == m: return pmm
        return _rec(n - m, pmm, x*(2*m+1)*pmm, lambda j, a, b: ((2*(j+m)+1)*x*b - (j+2*m)*a) / (j+1))
def SLEG(n, m, t):
    n, m = int(n), int(m)
    if m > n: return mp.mpf(0)
    with mp.workprec(2000):
        x, s = mp.cos(t), abs(mp.sin(t))
        pmm = mp.fac2(2*m-1) * s ** m if m else mp.mpf(1)
        p = pmm if n == m else _rec(n - m, pmm, x*(2*m+1)*pmm, lambda j, a, b: ((2*(j+m)+1)*x*b - (j+2*m)*a) / (j+1))
        return (-1)**m * mp.sqrt((2*n+1) / (4*mp.pi) * mp.factorial(n-m) / mp.factorial(n+m)) * p
def ALAG(n, m, x):
    n, m = int(n), int(m)
    with mp.workprec(2000): return _rec(n, mp.mpf(1), 1 + m - x, lambda j, a, b: ((2*j+1+m-x)*b - (j+m)*a) / (j+1))
def LAG(n, x): return ALAG(n, 0, x)
def HER(n, x):
    n = int(n)
    with mp.workprec(2000): return _rec(n, mp.mpf(1), 2*x, lambda j, a, b: 2*x*b - 2*j*a)

def selftest():
    with mp.workprec(600):
        for n, m, x in [(7, 3, mp.mpf('0.3')), (20, 0, mp.mpf('-0.7')), (13, 13, mp.mpf('0.99'))]:
            assert abs(ALEG(n, m, x) - mp.legenp(n, m, x, type=2) * (-1)**m) < mp.mpf(2)**-500 * max(1, abs(ALEG(n, m, x)))
            assert abs(SLEG(n, m, x) - mp.re(mp.spherharm(n, m, x, 0))) < mp.mpf(2)**-500
            assert abs(ALAG(n, m, 5*x) - mp.laguerre(n, m, 5*x)) < mp.mpf(2)**-500 * max(1, abs(ALAG(n, m, 5*x)))
            assert abs(HER(n, 3*x) - mp.hermite(n, 3*x)) < mp.mpf(2)**-500 * max(1, abs(HER(n, 3*x)))
            assert abs(LEG(n, x) - mp.legendre(n, x)) < mp.mpf(2)**-500
    print("selftest ok")

def EI(x):
    with mp.workprec(1500): return mp.ei(x)
def EN(n, x):
    with mp.workprec(1500): return mp.expint(int(n), x)

def W1500(f):
    def g(*a):
        with mp.workprec(1500): return f(*a)
    return g
SI, CI, SHI, CHI = W1500(mp.si), W1500(mp.ci), W1500(mp.shi), W1500(mp.chi)
LI = W1500(lambda x: mp.li(x))

PSI = W1500(mp.digamma)
PG = W1500(lambda n, x: mp.polygamma(int(n), x))
LBETA = W1500(lambda a, b: mp.log(mp.beta(a, b)))
GR = W1500(lambda a, b: mp.exp(mp.loggamma(a) - mp.loggamma(b)))
GDR = W1500(lambda a, d: mp.exp(mp.loggamma(a) - mp.loggamma(a + d)))

AI = W1500(lambda x: mp.airyai(x)); AIP = W1500(lambda x: mp.airyai(x, 1))
BI = W1500(lambda x: mp.airybi(x)); BIP = W1500(lambda x: mp.airybi(x, 1))

RANGES = [
    ("ai",       "ai", AI, [(-20, 5, True)]),
    ("ai_pos",   "ai", AI, [(0, 60, False)]),
    ("ai_neg",   "ai", AI, [(-60, 0, False)]),
    ("ai_far",   "ai", AI, [(6, 100, True)]),
    ("ai_big",   "ai", AI, [(60, 4900, False)]),
    ("aip_pos",  "aip", AIP, [(0, 60, False)]),
    ("aip_neg",  "aip", AIP, [(-60, 0, False)]),
    ("aip_far",  "aip", AIP, [(6, 100, True)]),
    ("bi_pos",   "bi", BI, [(0, 60, False)]),
    ("bi_neg",   "bi", BI, [(-60, 0, False)]),
    ("bi_far",   "bi", BI, [(6, 12, True)]),
    ("bip_pos",  "bip", BIP, [(0, 60, False)]),
    ("bip_neg",  "bip", BIP, [(-60, 0, False)]),
    ("bip_far",  "bip", BIP, [(6, 12, True)]),
    ("ai_nfar",  "ai", AI, [(-10000, -60, False)]),
    ("bi_nfar",  "bi", BI, [(-10000, -60, False)]),
    ("aip_nfar", "aip", AIP, [(-10000, -60, False)]),
    ("bip_nfar", "bip", BIP, [(-1e30, -60, False)]),
    ("psi",      "psi", PSI, [(-20, 10, True)]),
    ("psi_mid",  "psi", PSI, [(0, 100, False)]),
    ("psi_neg",  "psi", PSI, [(-60, 0, False)]),
    ("psi_big",  "psi", PSI, [(10, 3000, True)]),
    ("psi_root", "psi", PSI, [(1.4616, 1.4617, False)]),
    ("pg",       "pg", PG, [("int", 1, 10), (-10, 8, True)]),
    ("pg_bign",  "pg", PG, [("int", 10, 300), (-4, 8, True)]),
    ("lbeta",    "lbeta", LBETA, [(-20, 30, True), (-20, 30, True)]),
    ("lbeta_1",  "lbeta", LBETA, [(0.999, 1.001, False), (0.5, 2, False)]),
    ("gr",       "gr", GR, [(-10, 12, True), (-10, 12, True)]),
    ("gdr_big",  "gdr", GDR, [(20, 100, True), (-10, 5, True)]),
    ("gdr",      "gdr", GDR, [(-10, 15, True), (-60, -1, True)]),
    ("si",       "si", SI, [(-10, 6, True)]),
    ("si_mid",   "si", SI, [(-400, 400, False)]),
    ("si_big",   "si", SI, [(8, 1000, True)]),
    ("ci",       "ci", CI, [(-10, 6, True)]),
    ("ci_mid",   "ci", CI, [(0, 400, False)]),
    ("ci_big",   "ci", CI, [(8, 1000, True)]),
    ("si_huge",  "si", SI, [(900, 262000, True)]),
    ("ci_huge",  "ci", CI, [(900, 262000, True)]),
    ("ci_root",  "ci", CI, [(0.6165, 0.6166, False)]),
    ("shi",      "shi", SHI, [(-20, 9, True)]),
    ("shi_big",  "shi", SHI, [(9, 17.47, True)]),
    ("chi",      "chi", CHI, [(-20, 9, True)]),
    ("chi_big",  "chi", CHI, [(9, 17.47, True)]),
    ("chi_root", "chi", CHI, [(0.5238, 0.5239, False)]),
    ("li",       "li", LI, [(-20, 20, True)]),
    ("li_1",     "li", LI, [(0.999, 1.001, False)]),
    ("li_big",   "li", LI, [(20, 3000, True)]),
    ("li_tiny",  "li", LI, [(-262000, -20, True)]),
    ("ei",       "ei", EI, [(-10, 6, True)]),
    ("ei_neg",   "ei", EI, [(-400, 0, False)]),
    ("ei_mid",   "ei", EI, [(0, 500, False)]),
    ("ei_big",   "ei", EI, [(9, 17.4, True)]),
    ("ei_tiny",  "ei", EI, [(-300, -10, True)]),
    ("ei_root",  "ei", EI, [(0.37, 0.375, False)]),
    ("e1",       "en", EN, [("int", 1, 1), (-10, 8, True)]),
    ("en",       "en", EN, [("int", 0, 30), (-10, 8, True)]),
    ("en_mid",   "en", EN, [("int", 2, 200), (0, 80, False)]),
    ("en_bign",  "en", EN, [("int", 200, 5000), (-20, 10, True)]),
    ("en_tiny",  "en", EN, [("int", 1, 30), (-300, -10, True)]),
    ("leg",      "leg", LEG, [("int", 0, 60), (-1, 1, False)]),
    ("leg_big",  "leg", LEG, [("int", 100, 1000), (-1, 1, False)]),
    ("leg_out",  "leg", LEG, [("int", 0, 60), (0, 40, True)]),
    ("aleg",     "aleg", ALEG, [("int", 0, 40), ("int", 0, 40), (-1, 1, False)]),
    ("aleg_big", "aleg", ALEG, [("int", 100, 400), ("int", 0, 100), (-1, 1, False)]),
    ("aleg_1",   "aleg", ALEG, [("int", 0, 30), ("int", 0, 30), (0.999, 1, False)]),
    ("sleg",     "sleg", SLEG, [("int", 0, 40), ("int", 0, 40), (-4, 4, False)]),
    ("sleg_big", "sleg", SLEG, [("int", 100, 300), ("int", 0, 100), (0, 3.2, False)]),
    ("sleg_0",   "sleg", SLEG, [("int", 0, 30), ("int", 0, 30), (-60, -2, True)]),
    ("lag",      "lag", LAG, [("int", 0, 60), (0, 200, False)]),
    ("lag_small","lag", LAG, [("int", 0, 60), (-30, 3, True)]),
    ("lag_big",  "lag", LAG, [("int", 100, 1000), (0, 4000, False)]),
    ("lag_neg",  "lag", LAG, [("int", 0, 60), (-100, 0, False)]),
    ("alag",     "alag", ALAG, [("int", 0, 60), ("int", 0, 60), (0, 300, False)]),
    ("her",      "her", HER, [("int", 0, 60), (-12, 12, False)]),
    ("her_big",  "her", HER, [("int", 100, 1000), (-50, 50, False)]),
    ("her_huge", "her", HER, [("int", 0, 60), (3, 100, True)]),
    ("ell_K",    "K", EK, [(-1, 1, False)]),
    ("ell_K1",   "K", EK, [(0.999, 1, False)]),
    ("ell_E",    "E", EE, [(-1, 1, False)]),
    ("ell_Pi",   "Pi", EPI, [(-1, 1, False), (-10, 1, False)]),
    ("ell_Pi1",  "Pi", EPI, [(0.999, 1, False), (0.999, 1, False)]),
    ("ell_Pibig","Pi", EPI, [(-1, 1, False), (-1e6, -10, False)]),
    ("ell_F",    "F1", EF, [(-1, 1, False), (-1.57, 1.57, False)]),
    ("ell_Fbig", "F1", EF, [(-1, 1, False), (-100, 100, False)]),
    ("ell_E2",   "E2", EE2, [(-1, 1, False), (-1.57, 1.57, False)]),
    ("ell_E2big","E2", EE2, [(-1, 1, False), (-100, 100, False)]),
    ("ell_P3",   "P3", EP3, [(-1, 1, False), (-10, 0.9, False), (-1.57, 1.57, False)]),
    ("ell_P3big","P3", EP3, [(-0.99, 0.99, False), (-10, 0.9, False), (-20, 20, False)]),
    ("gl",       "gl", GL, [(-10, 8, True), (-10, 9, True)]),
    ("gu",       "gu", GU, [(-10, 8, True), (-10, 9, True)]),
    ("gp",       "gp", GP, [(-10, 8, True), (-10, 9, True)]),
    ("gq",       "gq", GQ, [(-10, 8, True), (-10, 9, True)]),
    ("gq_mid",   "gq", GQ, [(0, 40, False), (0, 40, False)]),
    ("gp_mid",   "gp", GP, [(0, 40, False), (0, 40, False)]),
    ("gu_mid",   "gu", GU, [(0, 3, False), (1, 40, False)]),
    ("gq_tinya", "gq", GQ, [(-200, -10, True), (-10, 3, True)]),
    ("gu_tinya", "gu", GU, [(-200, -10, True), (-10, 3, True)]),
    ("gp_near",  "gp", GP, [(100, 3000, False), (100, 3000, False)]),
    ("gq_near",  "gq", GQ, [(100, 3000, False), (100, 3000, False)]),
    ("zeta",     "zeta", ZETA, [(-60, 300, False)]),
    ("hz",       "hz", HZ, [(1, 60, False), (-20, 20, True)]),
    ("hz_01",    "hz", HZ, [(0, 1, False), (-10, 10, True)]),
    ("hz_neg",   "hz", HZ, [(-88, 0, False), (-10, 10, True)]),
    ("hz_vneg",  "hz", HZ, [(-600, -89, False), (0, 4, False)]),
    ("hz_negba", "hz", HZ, [(-600, -89, False), (100, 2000, False)]),
    ("hz_s1",    "hz", HZ, [(1, 1.01, False), (-5, 5, True)]),
    ("hz_bigs",  "hz", HZ, [(60, 3000, False), (-3, 10, True)]),
    ("hz_biga",  "hz", HZ, [(1, 30, False), (20, 300, True)]),
    ("zeta_neg", "zeta", ZETA, [(-3000, -60, False)]),
    ("zeta_01",  "zeta", ZETA, [(-1, 1, False)]),
    ("zeta_tz",  "zeta", ZETA, [("nearnegeven", 1, 500)]),
    ("zeta_tiny","zeta", ZETA, [(-300, -2, True)]),
    ("beta",     "beta", BETA, [(-20, 12, True), (-20, 12, True)]),
    ("beta_mix", "beta", BETA, [(-200, -20, True), (0, 40, True)]),
    ("beta_big", "beta", BETA, [(20, 60, True), (-10, 60, True)]),
    ("ibeta_skw","ibeta", IB, [(-60, -5, True), (5, 15, True), (-30, 0, True)]),
    ("ibeta_x1", "ibeta", IB, [(-3, 4, True), (-3, 4, True), (0.999, 1, False)]),
    ("ibeta",    "ibeta", IB, [(-6, 7, True), (-6, 7, True), (0, 1, False)]),
    ("ibeta_sm", "ibeta", IB, [(0, 10, False), (0, 10, False), (-60, -1, True)]),
    ("ibetac",   "ibetac", IBC, [(-6, 7, True), (-6, 7, True), (0, 1, False)]),
    ("ibeta_big","ibeta", IB, [(1000, 5000, False), (1000, 5000, False), (0.3, 0.7, False)]),
    ("j_small", "j", J, [(0, 5, False), (-8, 1, True)]),
    ("j_mid",   "j", J, [(0, 20, False), (0, 100, False)]),
    ("j_big",   "j", J, [(0, 20, False), (7, 12, True)]),
    ("j_hinu",  "j", J, [(20, 200, False), (0, 200, False)]),
    ("y_small", "y", Y, [(0, 5, False), (-8, 1, True)]),
    ("y_mid",   "y", Y, [(0, 20, False), (0, 100, False)]),
    ("y_big",   "y", Y, [(0, 20, False), (7, 12, True)]),
    ("y_temme", "y", Y, [(0, 5, False), (40, 120, False)]),
    ("k_temme", "k", K, [(0, 5, False), (30, 90, False)]),
    ("k_wr",    "k", K, [(0, 3, False), (28, 125, False)]),
    ("y_steed", "y", Y, [(0, 3, False), (50, 125, False)]),
    ("y_int",   "y", Y, [("int", 0, 30), (0, 60, False)]),
    ("i_small", "i", I, [(0, 5, False), (-8, 1, True)]),
    ("i_mid",   "i", I, [(0, 20, False), (0, 150, False)]),
    ("i_big",   "i", I, [(0, 20, False), (7, 14, True)]),
    ("k_small", "k", K, [(0, 5, False), (-8, 1, True)]),
    ("k_mid",   "k", K, [(0, 20, False), (0, 150, False)]),
    ("k_big",   "k", K, [(0, 20, False), (7, 14, True)]),
    ("k_int",   "k", K, [("int", 0, 30), (0, 60, False)]),
    ("j_hugex",  "j", J, [(0, 10, False), (14, 200, True)]),
    ("y_hugex",  "y", Y, [(0, 10, False), (14, 200, True)]),
    ("j_tinyx",  "j", J, [(0, 10, False), (-400, -8, True)]),
    ("y_tinyx",  "y", Y, [(0, 10, False), (-400, -8, True)]),
    ("k_tinyx",  "k", K, [(0, 10, False), (-400, -8, True)]),
    ("j_nugtx",  "j", J, [(300, 1000, False), (100, 900, False)]),
    ("y_nugtx",  "y", Y, [(300, 1000, False), (100, 900, False)]),
    ("i_nugtx",  "i", I, [(100, 400, False), (100, 900, False)]),
    ("k_nugtx",  "k", K, [(100, 400, False), (100, 900, False)]),
    ("y_nearint","y", Y, [("nearint", 0, 5), (0, 200, False)]),
    ("k_nearint","k", K, [("nearint", 0, 5), (0, 200, False)]),
    ("sph_j",   "sj", SJ, [("int", 0, 20), (0, 50, False)]),
    ("sph_y",   "sy", SY, [("int", 0, 20), (0, 50, False)]),
]

def main():
    if sys.argv[1] == "--selftest":
        selftest(); return
    exe = sys.argv[1]
    n = int(sys.argv[2]) if len(sys.argv) > 2 else 200
    rng = random.Random(int(sys.argv[3]) if len(sys.argv) > 3 else 1)
    flt = sys.argv[4] if len(sys.argv) > 4 else ""
    bad_total = 0
    for name, op, ref, args in RANGES:
        if flt and flt not in name:
            continue
        cases = []
        for _ in range(n):
            hs, vs = [], []
            for a in args:
                if a[0] == "nearnegeven":
                    v = -2 * mp.mpf(rng.randint(a[1], a[2])) + mp.ldexp(rng.choice([-1, 1]), -rng.randint(20, 200))
                    h = to_hex(v); v = parse_hex(h); hs.append(h); vs.append(v)
                elif a[0] == "nearint":
                    v = mp.mpf(rng.randint(a[1], a[2])) + mp.ldexp(rng.choice([-1, 1]), -rng.randint(20, 200))
                    if v < 0: v = -v
                    h = to_hex(v); v = parse_hex(h); hs.append(h); vs.append(v)
                elif a[0] == "int":
                    k = rng.randint(a[1], a[2]); hs.append(str(k)); vs.append(mp.mpf(k))
                else:
                    h, v = rnd256(a[0], a[1], a[2], rng); hs.append(h); vs.append(v)
            cases.append((hs, vs))
        inp = "".join(f"{op} {' '.join(h)}\n" for h, _ in cases)
        res = subprocess.run([exe], input=inp, capture_output=True, text=True, check=True).stdout.split("\n")
        worst, wcase, mis = 0.0, None, 0
        for (hs, vs), r in zip(cases, res):
            got = parse_hex(r)
            want = ref(*vs)
            err = ulp_err(got, want)
            if err > 0.5: mis += 1
            if err > worst: worst, wcase = err, hs
        bad_total += worst > 1
        print(f"{name:10s} n={n} max ulp {worst:.3g}  not correctly rounded {mis}" + (f"  worst at {wcase}" if worst > 0.5 else ""), flush=True)
    sys.exit(1 if bad_total else 0)

if __name__ == "__main__":
    main()
