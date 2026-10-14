#!/usr/bin/env python3
"""Prints reference values (mpmath, 700 bits, rounded to binary256 nearest) for
tests/test_special_math.cpp: {"op", "a", "b", "c", "result"} lines."""
import sys
sys.path.insert(0, __file__.rsplit("/", 1)[0])
import mpmath as mp
from sf_check import to_hex
mp.mp.prec = 700

def r256(v):
    # round to 237 bits (nearest even) in a normal range
    e = int(mp.floor(mp.log(abs(v), 2)))
    m = mp.nint(mp.ldexp(v, 236 - e))
    return to_hex(mp.ldexp(m, e - 236))

H = mp.mpf(1) / 2

def rd(s):
    """the binary256 value of a decimal literal (round to nearest)"""
    v = mp.mpf(s)
    if v == 0 or (v == int(v) and abs(v) < 2**200):
        return v
    e = int(mp.floor(mp.log(abs(v), 2)))
    return mp.ldexp(mp.nint(mp.ldexp(v, 236 - e)), e - 236)

def sph(f, n, x): return mp.sqrt(mp.pi / (2 * x)) * f(n + H, x)
F = {
    "j": mp.besselj, "y": mp.bessely, "i": mp.besseli, "k": mp.besselk,
    "sj": lambda n, x: sph(mp.besselj, n, x), "sy": lambda n, x: sph(mp.bessely, n, x),
    "jn": mp.besselj, "yn": mp.bessely,
    "beta": mp.beta,
    "ibeta": lambda a, b, x: mp.betainc(a, b, 0, x, regularized=True),
    "ibetac": lambda a, b, x: mp.betainc(a, b, x, 1, regularized=True),
}
CASES = [
    ("j", "0", "1"),
    ("j", "0.5", "3"),
    ("j", "2.25", "17.5"),
    ("j", "10", "1000"),
    ("j", "150", "100"),
    ("y", "0", "1"),
    ("y", "1.75", "0.125"),
    ("y", "3", "60"),
    ("y", "0.5", "500"),
    ("i", "0", "1"),
    ("i", "2.5", "30"),
    ("i", "1", "1000"),
    ("k", "0", "1"),
    ("k", "0.25", "0.0625"),
    ("k", "7", "70"),
    ("k", "1", "1000"),
    ("sj", "3", "2"),
    ("sy", "2", "1.5"),
    ("jn", "-3", "-2.5"),
    ("yn", "-3", "2.5"),
    ("beta", "0.5", "0.5"),
    ("beta", "3.25", "1e-5"),
    ("beta", "2000.5", "3000.25"),
    ("ibeta", "2.5", "3.5", "0.25"),
    ("ibeta", "0.001", "50", "0.5"),
    ("ibeta", "400", "600", "0.375"),
    ("ibetac", "2.5", "3.5", "0.875"),
    ("ibetac", "30", "2", "0.0625"),
]

def main():
    which = set(sys.argv[1:])
    for c in CASES:
        op, args = c[0], c[1:]
        if which and op not in which:
            continue
        a = list(args) + [""] * (3 - len(args))
        print(f'    {{"{op}", "{a[0]}", "{a[1]}", "{a[2]}", "{r256(F[op](*[rd(x) for x in args]))}"}},')

if __name__ == "__main__":
    main()
