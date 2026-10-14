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
    s = s.strip()
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

RANGES = [
    ("ell_K",    "K", EK, [(-1, 1, False)]),
    ("ell_K1",   "K", EK, [(0.999, 1, False)]),
    ("ell_E",    "E", EE, [(-1, 1, False)]),
    ("ell_Pi",   "Pi", EPI, [(-1, 1, False), (-10, 1, False)]),
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
