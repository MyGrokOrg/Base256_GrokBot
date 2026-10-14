#!/usr/bin/env python3
"""Accuracy check of the special functions against mpmath (independent reference).

usage: sf_check.py EVAL_BINARY [cases_per_range] [seed] [filter] [modes]

With "modes" as the fifth argument, every case is evaluated in all five
rounding modes (SF_MODES=1) and checked against the correctly rounded value
of a 2000-bit reference in that mode, including the IEEE flags
(inexact, overflow, underflow with tininess after rounding).

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
# x < 0, integer order: real by the reflection identities
def JN(n, x): return (-1) ** int(n) * mp.besselj(n, -x)
def IN(n, x): return (-1) ** int(n) * mp.besseli(n, -x)
def SJN(n, x): return (-1) ** int(n) * SJ(n, -x)
def SYN(n, x): return (-1) ** (int(n) + 1) * SY(n, -x)

def BETA(a, b): return mp.beta(a, b)
# mpmath's betainc can lose bits to cancellation; evaluate it at 2000 bits.
def IB(a, b, x):
    with mp.workprec(2000): return mp.betainc(a, b, 0, x, regularized=True)
def IBC(a, b, x):
    with mp.workprec(2000): return mp.betainc(a, b, x, 1, regularized=True)

def ZETA(x): return mp.zeta(x)
# mpmath's zeta(s, a) loses hundreds of bits at large a; evaluate it at 3000 bits.
def _hz_hurwitz(s, a, bits):
    """ζ(s, a), s <= −89, 0 < a <= 4, by Hurwitz's formula at `bits` (mpmath's zeta
    can stall here, e.g. s ≈ −519.84, a ≈ 0.541 at 3000 bits)."""
    with mp.workprec(bits + 64):
        sig = 1 - s
        m = int(mp.ceil(a)) - 1
        b = a - m
        S = mp.mpf(0); n = 1
        while True:
            t = mp.power(n, -sig)
            S += mp.cos(mp.pi * sig / 2 - 2 * mp.pi * n * b) * t
            if t < mp.mpf(2)**(-bits - 32): break
            n += 1
        r = 2 * mp.gamma(sig) / (2 * mp.pi)**sig * S
        for j in range(m): r -= mp.power(b + j, -s)
        return r

def HZ(s, a):
    if s <= -89 and a <= 4:
        # two precisions must agree (guards cancellation near zeros)
        r1, r2 = _hz_hurwitz(s, a, 1200), _hz_hurwitz(s, a, 1800)
        if r2 == 0 or abs(r1 - r2) <= mp.mpf(2)**-700 * abs(r2): return r2
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
    with mp.workprec(REFPREC[0]): return mp.ei(x)
def EN(n, x):
    with mp.workprec(REFPREC[0]): return mp.expint(int(n), x)

REFPREC = [1500]
def W1500(f):
    def g(*a):
        with mp.workprec(REFPREC[0]): return f(*a)
    return g
SI, CI, SHI, CHI = W1500(mp.si), W1500(mp.ci), W1500(mp.shi), W1500(mp.chi)
LI = W1500(lambda x: mp.li(x))

PSI = W1500(mp.digamma)
_PG = W1500(lambda n, x: mp.polygamma(int(n), x))
_COTP = {}
def _cot_poly(n):
    """integer coefficients of P_n with d^n/dθ^n cot θ = P_n(cot θ)."""
    if n not in _COTP:
        p = [0, 1]
        for _ in range(n):
            d = [i * p[i] for i in range(1, len(p))] + [0, 0]  # P'
            q = [0] * (len(d) + 2)
            for i, c in enumerate(d): q[i] -= c; q[i + 2] -= c   # × (−1 − c²)
            while len(q) > 1 and q[-1] == 0: q.pop()
            p = q
        _COTP[n] = p
    return _COTP[n]
def PG(n, x):
    """mp.polygamma, except x < −64: reflection ψ^(n)(x) = (−1)^n ψ^(n)(1−x) − π^{n+1} cot^{(n)}(πx),
    with x mod 1 taken exactly (mp.polygamma is very slow there)."""
    n = int(n)
    # next to a pole the δ^{-s} term exceeds the rest by s·log2(1/δ) bits: carry that many more
    extra = 0
    if x < 0:
        with mp.workprec(REFPREC[0] + 400):
            d = abs(x - mp.nint(x))
            if d != 0: extra = int((n + 1) * max(0, -mp.log(d, 2))) + 64
    if extra:
        with mp.workprec(REFPREC[0] + extra):
            if x >= -64: return mp.polygamma(n, x)
    if x >= -64: return _PG(n, x)
    with mp.workprec(REFPREC[0] + 64 + extra):
        fr = x - mp.floor(x)  # exact: x is a binary256 value
        c = mp.cot(mp.pi * fr)
        cp = sum(k * c**i for i, k in enumerate(_cot_poly(n)))
        return (-1)**n * mp.polygamma(n, 1 - x) - mp.pi**(n + 1) * cp
LBETA = W1500(lambda a, b: mp.log(mp.beta(a, b)))
GR = W1500(lambda a, b: mp.exp(mp.loggamma(a) - mp.loggamma(b)))
GDR = W1500(lambda a, d: mp.exp(mp.loggamma(a) - mp.loggamma(a + d)))

def _airy_huge(k):
    def f(x):
        e = int(mp.log(abs(x), 2)) if x != 0 else 0
        with mp.workprec(max(REFPREC[0], 3 * e // 2 + 900)):
            return mp.airyai(x, k % 2) if k < 2 else mp.airybi(x, k % 2)
    return f
AIH, AIPH, BIH, BIPH = (_airy_huge(k) for k in range(4))

def _airy_asym(k):
    """Ai, Ai', Bi, Bi' at x <= −2^90 from the DLMF 9.7.9–9.7.12 asymptotic series
    (mpmath's own airyai is too slow past −2^12000). ζ is reduced mod 2π at
    1.5·log2|x| + 900 bits; the series then runs at 1200 bits (ζ^-1 < 2^-135, so a
    few terms reach 2^-1000). Agrees with mp.airyai to ~2^-897 relative at −2^95…−2^3000 (16 spot checks)."""
    def f(x):
        z = -x
        e = int(mp.log(z, 2))
        with mp.workprec(3 * e // 2 + 900):
            zeta = 2 * z * mp.sqrt(z) / 3
            th = zeta - mp.pi / 4
            th = th - 2 * mp.pi * mp.floor(th / (2 * mp.pi))
        with mp.workprec(1200):
            zeta = mp.mpf(2) * z * mp.sqrt(z) / 3
            st, ct = mp.sin(th), mp.cos(th)
            ev = od = mp.mpf(0); u = mp.mpf(1); kk = 0
            while True:
                c = u if k % 2 == 0 else -(6 * kk + 1) * u / (6 * kk - 1)
                term = c / zeta**kk
                sgn = (-1) ** (kk // 2)
                if kk % 2 == 0: ev += sgn * term
                else: od += sgn * term
                if kk > 2 and abs(term) < mp.mpf(2)**-1150: break
                kk += 1
                u = u * (6 * kk - 5) * (6 * kk - 3) * (6 * kk - 1) / ((2 * kk - 1) * 216 * kk)
            pre = 1 / mp.sqrt(mp.pi) * (z ** (mp.mpf(1) / 4) if k % 2 else z ** (-mp.mpf(1) / 4))
            if k == 0: return pre * (ct * ev + st * od)
            if k == 2: return pre * (ct * od - st * ev)
            if k == 1: return pre * (st * ev - ct * od)
            return pre * (ct * ev + st * od)
    return f
AIA, AIPA, BIA, BIPA = (_airy_asym(k) for k in range(4))
AI = W1500(lambda x: mp.airyai(x)); AIP = W1500(lambda x: mp.airyai(x, 1))
BI = W1500(lambda x: mp.airybi(x)); BIP = W1500(lambda x: mp.airybi(x, 1))

def ERFINV(x):
    with mp.workprec(REFPREC[0]): return mp.erfinv(x)
def ERFCINV(y):
    with mp.workprec(REFPREC[0]):
        if y < 1: return mp.findroot(lambda w: mp.log(mp.erfc(w)) - mp.log(y), mp.sqrt(-mp.log(y)) if y < 0.1 else mp.erfinv(1 - y))
        return -mp.erfinv(y - 1) if y != 1 else mp.mpf(0)
W0 = W1500(lambda x: mp.lambertw(x).real)
WM1 = W1500(lambda x: mp.lambertw(x, -1).real)

FS = W1500(mp.fresnels); FC = W1500(mp.fresnelc)
def DAW(x):
    with mp.workprec(REFPREC[0]):
        if abs(x) < 40: return mp.sqrt(mp.pi) / 2 * mp.exp(-x * x) * mp.erfi(x)
        s, t, k = mp.mpf(1), mp.mpf(1), 1
        while abs(t) > mp.mpf(2)**-1400:
            t *= (2 * k - 1) / (2 * x * x); s += t; k += 1
        return s / (2 * x)

LI2 = W1500(lambda x: mp.polylog(2, x).real)
PLOG = W1500(lambda n, x: mp.polylog(int(n), x).real)
ETA = W1500(lambda s: mp.altzeta(s))
BERN = W1500(lambda n: mp.bernoulli(2 * int(n)))

def LG(x):
    with mp.workprec(max(REFPREC[0], 3000)): return mp.log(abs(mp.gamma(x)))
_LGZ = []
def lg_zeros():
    """zeros of ln|Γ(x)| for x in (−41, −2) (two per unit interval)."""
    if not _LGZ:
        with mp.workprec(3000):
            f = lambda x: mp.log(abs(mp.gamma(x)))
            for n in range(2, 41):
                if n == 2: guesses = [mp.mpf(-2.457), mp.mpf(-2.747)]
                elif n == 3: guesses = [mp.mpf(-3.1436), mp.mpf(-3.9553)]
                else: guesses = [-n - 1 / mp.factorial(n), -n - 1 + 1 / mp.factorial(n + 1)]
                for g in guesses:
                    x = g
                    for _ in range(200):  # Newton, (ln|Γ|)' = ψ
                        d = f(x) / mp.digamma(x)
                        x -= d
                        if abs(d) < mp.mpf(2)**-2800 * abs(x): break
                    assert abs(f(x)) < mp.mpf(2)**-2500 * abs(mp.digamma(x)), (n, g)
                    _LGZ.append(x)
    return _LGZ

_HZZ = []
def hz_zeros():
    """(s, a) with ζ(s, a) = 0, a in (0, 1), for a few s < 1 (both EM and reflection paths)."""
    if not _HZZ:
        for s in [mp.mpf(-1), mp.mpf(-3), mp.mpf(-4), mp.mpf(-0.5), mp.mpf(0.5), mp.mpf(-2.5), mp.mpf(-7.25),
                  mp.mpf(-30.75), mp.mpf(-100.5)]:
            with mp.workprec(400):
                g = [mp.mpf(i) / 100 for i in range(1, 100)]
                v = [mp.zeta(s, x) for x in g]
            for i in range(len(g) - 1):
                if v[i] == 0 or v[i] * v[i + 1] < 0:
                    with mp.workprec(1200):
                        x = mp.findroot(lambda a: mp.zeta(s, a), (mp.mpf(g[i]), mp.mpf(g[i + 1])), solver="anderson", tol=mp.mpf(2)**-1100)
                        if abs(x - mp.mpf(1) / 2) < mp.mpf(2)**-1000: continue  # exact zero, tested separately
                        _HZZ.append((s, x))
    return _HZZ

_PGZ = []
def pg_zeros():
    """(n, x) with ψ^(n)(x) = 0, even n, one per interval (−k, −k+1)."""
    if not _PGZ:
        for n in (2, 4, 6):
            for k in range(1, 13):
                with mp.workprec(1200):
                    lo, hi = mp.mpf(-k) + mp.mpf(2)**-20, mp.mpf(-k + 1) - mp.mpf(2)**-20
                    if k == 1: hi = -mp.mpf(2)**-20
                    x = mp.findroot(lambda y: mp.polygamma(n, y), (lo, hi), solver="anderson", tol=mp.mpf(2)**-1100)
                _PGZ.append((n, x))
    return _PGZ

RANGES = [
    ("pg_negm",  "pg", PG, [("int", 1, 10), (-50, 0, False)]),
    ("pg_negs",  "pg", PG, [("int", 1, 10), (-1, 0, False)]),
    ("pg_negbig","pg", PG, [("int", 1, 10), ("negexp", 1, 230)]),
    ("pg_negtiny","pg", PG, [("int", 1, 10), ("negexp", -300, -1)]),
    ("pg_negn",  "pg", PG, [("int", 10, 300), (-20, 0, False)]),
    ("pg_negpole","pg", PG, [("int", 1, 10), ("negnearint", 1, 40)]),
    ("pg_negz",  "pg", PG, [("pairzero", pg_zeros, 1, 236), None]),
    ("hz_nz",    "hz", HZ, [("pairzero", hz_zeros, 1, 236), None]),
    ("lg_negz",  "lg", LG, [("nearzero", lg_zeros, 1, 236)]),
    ("lg_neg",   "lg", LG, [(-60, -2, False)]),
    ("li2",      "li2", LI2, [(-3, 1, False)]),
    ("li2_s",    "li2", LI2, [(-300, -1, True)]),
    ("li2_neg",  "li2", LI2, [(-1e6, -1, False)]),
    ("li2_1",    "li2", LI2, [(0.999, 1, False)]),
    ("plog",     "plog", PLOG, [("int", 0, 12), (-1, 1, False)]),
    ("plog_1",   "plog", PLOG, [("int", 3, 70), (0.99, 1, False)]),
    ("plog_n",   "plog", PLOG, [("int", 3, 70), (-1, -0.5, False)]),
    ("eta",      "eta", ETA, [(-30, 40, False)]),
    ("eta_1",    "eta", ETA, [(0.999, 1.001, False)]),
    ("eta_neg",  "eta", ETA, [(-3000, -30, False)]),
    ("bern",     "bern", BERN, [("int", 1, 400)]),
    ("fs",       "fs", FS, [(-20, 4, True)]),
    ("fs_mid",   "fs", FS, [(-20, 20, False)]),
    ("fs_big",   "fs", FS, [(4, 200, True)]),
    ("fs_tiny",  "fs", FS, [(-3000, -20, True)]),
    ("fc",       "fc", FC, [(-20, 4, True)]),
    ("fc_mid",   "fc", FC, [(-20, 20, False)]),
    ("fc_big",   "fc", FC, [(4, 200, True)]),
    ("daw",      "daw", DAW, [(-20, 5, True)]),
    ("daw_mid",  "daw", DAW, [(-30, 30, False)]),
    ("daw_big",  "daw", DAW, [(5, 262000, True)]),
    ("erfinv",   "erfinv", ERFINV, [(-1, 1, False)]),
    ("erfinv_s", "erfinv", ERFINV, [(-300, -2, True)]),
    ("erfinv_1", "erfinv", ERFINV, [(0.999999, 1, False)]),
    ("erfcinv",  "erfcinv", ERFCINV, [(0, 2, False)]),
    ("erfcinv_t","erfcinv", ERFCINV, [(-262000, -2, True)]),
    ("erfcinv_m","erfcinv", ERFCINV, [(-60, -1, True)]),
    ("w0",       "w0", W0, [(-0.3678794411714423, 3, False)]),
    ("w0_log",   "w0", W0, [(-300, 300, True)]),
    ("w0_big",   "w0", W0, [(300, 262000, True)]),
    ("w0_br",    "w0", W0, [(-0.36787944117144233, -0.3678794, False)]),
    ("wm1",      "wm1", WM1, [(-0.3678794411714423, 0, False)]),
    ("wm1_br",   "wm1", WM1, [(-0.36787944117144233, -0.3678794, False)]),
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
    ("ai_huge",  "ai", AIH, [("negexp", 94, 3000)]),
    ("aip_huge", "aip", AIPH, [("negexp", 94, 3000)]),
    ("bi_huge",  "bi", BIH, [("negexp", 94, 3000)]),
    ("bip_huge", "bip", BIPH, [("negexp", 94, 3000)]),
    ("ai_vhuge", "ai", AIH, [("negexp", 3000, 12000)]),
    ("ai_deep",  "ai", AIA, [("negexp", 12000, 175000)]),
    ("aip_deep", "aip", AIPA, [("negexp", 12000, 175000)]),
    ("bi_deep",  "bi", BIA, [("negexp", 12000, 175000)]),
    ("bip_deep", "bip", BIPA, [("negexp", 12000, 175000)]),
    ("ai_full",  "ai", AIA, [("negexp", 175000, 262144)]),
    ("aip_full", "aip", AIPA, [("negexp", 175000, 262144)]),
    ("bi_full",  "bi", BIA, [("negexp", 175000, 262144)]),
    ("bip_full", "bip", BIPA, [("negexp", 175000, 262144)]),
    ("ai_nfar",  "ai", AI, [(-10000, -60, False)]),
    ("ai_unf",   "ai", AI, [(4180, 4230, False)]),
    ("aip_unf",  "aip", AIP, [(4180, 4230, False)]),
    ("bi_ovf",   "bi", BI, [(4180, 4230, False)]),
    ("bip_ovf",  "bip", BIP, [(4180, 4230, False)]),
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
    ("j_negx",   "j", JN, [("int", 0, 40), (-100, 0, False)]),
    ("j_negbig", "j", JN, [("int", 0, 20), ("negexp", 7, 200)]),
    ("j_negtiny","j", JN, [("int", 0, 10), ("negexp", -400, -8)]),
    ("i_negx",   "i", IN, [("int", 0, 30), (-150, 0, False)]),
    ("i_negtiny","i", IN, [("int", 0, 10), ("negexp", -400, -8)]),
    ("sph_jneg", "sj", SJN, [("int", 0, 20), (-50, 0, False)]),
    ("sph_yneg", "sy", SYN, [("int", 0, 20), (-50, 0, False)]),
    ("sph_j",   "sj", SJ, [("int", 0, 20), (0, 50, False)]),
    ("sph_y",   "sy", SY, [("int", 0, 20), (0, 50, False)]),
]

EMAX = 262143
MODES = ("RNE", "RNA", "RZ", "RU", "RD")
F_INV, F_DIV, F_OVF, F_UNF, F_INX = 1, 2, 4, 8, 16

def _rnd(a, s, e, mode):
    """|a| rounded to 237 bits at exponent floor e (as integer q·2^(e−236)); sign s."""
    q = mp.ldexp(a, P - 1 - e)
    fl = int(mp.floor(q)); frac = q - fl
    if frac == 0: return fl, True
    up = {"RNE": frac > 0.5 or (frac == 0.5 and fl & 1), "RNA": frac >= 0.5, "RZ": False,
          "RU": s > 0, "RD": s < 0}[mode]
    return fl + (1 if up else 0), False

def round_b256(v, mode):
    """Correctly rounded binary256 value of v in mode, and the flags it raises."""
    if isinstance(v, mp.mpf) and v._mpf_[3] + 64 > mp.mp.prec:  # a reference carried more bits
        with mp.workprec(v._mpf_[3] + 64): return _round_b256(v, mode)
    return _round_b256(v, mode)

def _round_b256(v, mode):
    if mp.isnan(v): return mp.nan, 0
    if mp.isinf(v): return v, 0
    if v == 0: return v, 0
    s = 1 if v > 0 else -1
    a = abs(v)
    e = int(mp.frexp(a)[1]) - 1
    # tininess after rounding with an unbounded exponent
    qu, _ = _rnd(a, s, e, mode)
    tiny = mp.ldexp(qu, e - (P - 1)) < mp.ldexp(1, EMIN)
    e2 = max(e, EMIN)
    q, exact = _rnd(a, s, e2, mode)
    r = mp.ldexp(q, e2 - (P - 1))
    flags = 0 if exact else F_INX
    if tiny and not exact: flags |= F_UNF
    if r >= mp.ldexp(1, EMAX + 1):
        flags |= F_OVF | F_INX
        to_inf = mode in ("RNE", "RNA") or (mode == "RU" and s > 0) or (mode == "RD" and s < 0)
        r = mp.inf if to_inf else mp.ldexp((1 << P) - 1, EMAX - (P - 1))
    return s * r, flags

def check_modes(exe, rng, n, flt):
    import os
    REFPREC[0] = 2000
    mp.mp.prec = 2000
    env = dict(os.environ, SF_MODES="1")
    bad_total = 0
    for name, op, ref, args in RANGES:
        if flt and flt not in name.split(",") and not any(name.startswith(f) for f in flt.split(",")):
            continue
        cases = gen_cases(args, rng, n)
        inp = "".join(f"{op} {' '.join(h)}\n" for h, _ in cases)
        res = subprocess.run([exe], input=inp, capture_output=True, text=True, check=True, env=env).stdout.split("\n")
        worst = [0.0] * 5
        mis = flagerr = 0
        first = None
        for (hs, vs), r in zip(cases, res):
            parts = r.split()
            want = ref(*vs)
            for i, m in enumerate(MODES):
                got = parse_hex(parts[2 * i]); gf = int(parts[2 * i + 1])
                exp_v, exp_f = round_b256(want, m)
                err = ulp_err(got, want)
                if not (mp.isinf(got) and mp.isinf(exp_v)) and err > worst[i]: worst[i] = err
                same = (got == exp_v) or (mp.isnan(got) and mp.isnan(exp_v))
                if not same: mis += 1
                if (gf & 0x1F) != exp_f and same: flagerr += 1
                if (not same or (gf & 0x1F) != exp_f) and first is None:
                    first = (m, hs, parts[2 * i], gf, to_hex(exp_v) if not mp.isinf(exp_v) else str(exp_v), exp_f)
        bad = mis + flagerr
        bad_total += bad
        print(f"{name:10s} n={n} max ulp " + " ".join(f"{m} {w:.3g}" for m, w in zip(MODES, worst)) +
              f" | misrounded {mis}, flag errors {flagerr}" + (f"  first {first}" if bad else ""), flush=True)
    sys.exit(1 if bad_total else 0)

def gen_cases(args, rng, n):
    cases = []
    for _ in range(n):
        hs, vs = [], []
        for a in args:
            if a is None: continue  # filled by the previous pair generator
            if a[0] == "negexp":
                h, v = rnd256(a[1], a[2], True, rng); v = -v; h = to_hex(v); hs.append(h); vs.append(v)
            elif a[0] == "pairzero":
                s, z = rng.choice(a[1]())
                k = rng.randint(a[2], a[3])
                v = z * (1 + rng.choice([-1, 1]) * mp.ldexp(rng.random(), -k))
                hs.append(str(s) if isinstance(s, int) else to_hex(s)); vs.append(mp.mpf(s))
                h = to_hex(v); v = parse_hex(h); hs.append(h); vs.append(v)
            elif a[0] == "negnearint":
                v = -mp.mpf(rng.randint(a[1], a[2])) + mp.ldexp(rng.choice([-1, 1]), -rng.randint(1, 200))
                h = to_hex(v); v = parse_hex(h); hs.append(h); vs.append(v)
            elif a[0] == "nearnegeven":
                v = -2 * mp.mpf(rng.randint(a[1], a[2])) + mp.ldexp(rng.choice([-1, 1]), -rng.randint(20, 200))
                h = to_hex(v); v = parse_hex(h); hs.append(h); vs.append(v)
            elif a[0] == "nearint":
                v = mp.mpf(rng.randint(a[1], a[2])) + mp.ldexp(rng.choice([-1, 1]), -rng.randint(20, 200))
                if v < 0: v = -v
                h = to_hex(v); v = parse_hex(h); hs.append(h); vs.append(v)
            elif a[0] == "int":
                k = rng.randint(a[1], a[2]); hs.append(str(k)); vs.append(mp.mpf(k))
            elif a[0] == "nearzero":
                # a binary256 point within 2^-k (relative) of a listed zero, k up to 236 (the nearest points)
                z = rng.choice(a[1]())
                k = rng.randint(a[2], a[3])
                v = z * (1 + rng.choice([-1, 1]) * mp.ldexp(rng.random(), -k))
                h = to_hex(v); v = parse_hex(h); hs.append(h); vs.append(v)
            else:
                h, v = rnd256(a[0], a[1], a[2], rng); hs.append(h); vs.append(v)
        cases.append((hs, vs))
    return cases

def main():
    if sys.argv[1] == "--selftest":
        selftest(); return
    exe = sys.argv[1]
    n = int(sys.argv[2]) if len(sys.argv) > 2 else 200
    rng = random.Random(int(sys.argv[3]) if len(sys.argv) > 3 else 1)
    flt = sys.argv[4] if len(sys.argv) > 4 else ""
    if len(sys.argv) > 5 and sys.argv[5] == "modes":
        check_modes(exe, rng, n, flt); return
    bad_total = 0
    for name, op, ref, args in RANGES:
        if flt and flt not in name:
            continue
        cases = gen_cases(args, rng, n)
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
