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
CASES = [
    ("j", "0", "1", lambda: mp.besselj(0, 1)),
    ("j", "0.5", "3", lambda: mp.besselj(H, 3)),
    ("j", "2.25", "17.5", lambda: mp.besselj(mp.mpf("2.25"), mp.mpf("17.5"))),
    ("j", "10", "1000", lambda: mp.besselj(10, 1000)),
    ("j", "150", "100", lambda: mp.besselj(150, 100)),
    ("y", "0", "1", lambda: mp.bessely(0, 1)),
    ("y", "1.75", "0.125", lambda: mp.bessely(mp.mpf("1.75"), mp.mpf("0.125"))),
    ("y", "3", "60", lambda: mp.bessely(3, 60)),
    ("y", "0.5", "500", lambda: mp.bessely(H, 500)),
    ("i", "0", "1", lambda: mp.besseli(0, 1)),
    ("i", "2.5", "30", lambda: mp.besseli(mp.mpf("2.5"), 30)),
    ("i", "1", "1000", lambda: mp.besseli(1, 1000)),
    ("k", "0", "1", lambda: mp.besselk(0, 1)),
    ("k", "0.25", "0.0625", lambda: mp.besselk(mp.mpf("0.25"), mp.mpf("0.0625"))),
    ("k", "7", "70", lambda: mp.besselk(7, 70)),
    ("k", "1", "1000", lambda: mp.besselk(1, 1000)),
    ("sj", "3", "2", lambda: mp.sqrt(mp.pi / 4) * mp.besselj(3 + H, 2)),
    ("sy", "2", "1.5", lambda: mp.sqrt(mp.pi / 3) * mp.bessely(2 + H, mp.mpf("1.5"))),
    ("jn", "-3", "-2.5", lambda: mp.besselj(-3, mp.mpf("-2.5"))),
    ("yn", "-3", "2.5", lambda: mp.bessely(-3, mp.mpf("2.5"))),
]

def main():
    which = set(sys.argv[1:])
    for c in CASES:
        op, args, f = c[0], c[1:-1], c[-1]
        if which and op not in which:
            continue
        a = list(args) + [""] * (3 - len(args))
        print(f'    {{"{op}", "{a[0]}", "{a[1]}", "{a[2]}", "{r256(f())}"}},')

if __name__ == "__main__":
    main()
