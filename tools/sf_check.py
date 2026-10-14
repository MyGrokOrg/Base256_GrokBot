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

RANGES = [
    ("j_small", "j", J, [(0, 5, False), (-8, 1, True)]),
    ("j_mid",   "j", J, [(0, 20, False), (0, 100, False)]),
    ("j_big",   "j", J, [(0, 20, False), (7, 12, True)]),
    ("j_hinu",  "j", J, [(20, 200, False), (0, 200, False)]),
    ("y_small", "y", Y, [(0, 5, False), (-8, 1, True)]),
    ("y_mid",   "y", Y, [(0, 20, False), (0, 100, False)]),
    ("y_big",   "y", Y, [(0, 20, False), (7, 12, True)]),
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
                if a[0] == "nearint":
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
