#!/usr/bin/env python3
"""Accuracy check of base256::complex against GNU MPC (through gmpy2).

usage: cx_check.py EVAL_BINARY [cases_per_range] [seed] [filter] [modes]
       cx_check.py --list                 (range names, one per line)

EVAL_BINARY is tools/cx_eval.cpp. Each range draws random binary256 operands,
evaluates them with EVAL_BINARY and compares every component with the
correctly rounded value of a reference:
  * mul, div, norm: exact rational arithmetic (Python integers);
  * everything else: GNU MPC/MPFR at 2000 bits; a component within 2^200 units
    of that precision of a rounding boundary is recomputed at 40000 bits, and
    one still that close is taken as exactly on the boundary (the exact cases:
    perfect squares, integer powers, ...);
  * |x|, |y| < 2^-64 for exp, sin ... atanh (and θ for polar): the exact leading
    Taylor term plus the rest of the series summed with MPC numbers at 3000 bits
    (MPC itself would need about 2·|log2 z| bits to round those correctly).
Prints per range the max error in ulps (both components) and the number of
components that are not the correctly rounded value.

With "modes" as the fifth argument every case runs in all five rounding modes
(RNE RNA RZ RU RD) and is compared with the reference rounded in that mode,
including the IEEE flags (inexact, overflow, underflow with tininess after
rounding; invalid and divide-by-zero must stay clear on these finite cases),
and the sign of exactly-zero components. Exit status 1 on any difference.

Environment (as tools/sf_check.py, used by the nightly fuzz workflow):
  SF_RANGE_SEED=1  reseed per range from (seed, range name).
  SF_SHARD=i/K     run only the ranges whose index in RANGES is i mod K.
  SF_MAXFAIL=n     print at most n FAIL lines per range (default 50).
  SF_TIMING=1      print the seconds per range.
"""
import os, random, subprocess, sys, time
from fractions import Fraction
import gmpy2
from gmpy2 import mpfr, mpc

P = 237
EMIN = -262142
EMAX = 262143
MODES = ("RNE", "RNA", "RZ", "RU", "RD")
F_INV, F_DIV, F_OVF, F_UNF, F_INX = 1, 2, 4, 8, 16

def set_prec(p):
    c = gmpy2.get_context()
    c.precision = p
    c.emin = gmpy2.get_emin_min()
    c.emax = gmpy2.get_emax_max()
    c.round = gmpy2.RoundToNearest

# ---------------------------------------------------------------- values
# An exact binary256 operand is (m, e): m·2^e, m an int with |m| < 2^237.
# A reference value is a Fraction (exact) or ("zero", negative).

def from_frac(x):
    """x (Fraction, nonzero) truncated to 237 bits, as (m, e)."""
    s = -1 if x < 0 else 1
    a = abs(x)
    e = a.numerator.bit_length() - a.denominator.bit_length()
    if Fraction(2) ** e > a: e -= 1
    m = int(a * Fraction(2) ** (P - 1 - e))
    return s * m, e - (P - 1)

def to_hex(op):
    m, e = op
    if m == 0: return "-0x0p+0" if e == "-" else "0x0p+0"
    s = "-" if m < 0 else ""
    return f"{s}0x{abs(m):x}p{e:+d}"

def to_mpfr(op):
    m, e = op
    return gmpy2.mul_2exp(mpfr(m, max(P, abs(m).bit_length() + 1)), e)

def to_frac(op):
    m, e = op
    return Fraction(m) * Fraction(2) ** e

def parse_hex(s):
    """C++ output token → ("nan",), ("inf", neg), ("zero", neg) or Fraction."""
    s = s.strip()
    if "nan" in s: return ("nan",)
    if "inf" in s: return ("inf", s.startswith("-"))
    neg = s.startswith("-")
    t = s.lstrip("+-")[2:]
    mant, ex = t.split("p")
    a, b = (mant.split(".") + [""])[:2]
    v = Fraction(int(a + b, 16)) * Fraction(2) ** (int(ex) - 4 * len(b))
    if v == 0: return ("zero", neg)
    return -v if neg else v

def mpfr_to_ref(v):
    if gmpy2.is_nan(v): return ("nan",)
    if gmpy2.is_infinite(v): return ("inf", v < 0)
    if v == 0: return ("zero", gmpy2.is_signed(v))
    m, e = v.as_mantissa_exp()
    return Fraction(int(m)) * Fraction(2) ** int(e)

def exp2_floor(a):
    """floor(log2 a) for a positive Fraction."""
    e = a.numerator.bit_length() - a.denominator.bit_length()
    if Fraction(2) ** e > a: e -= 1
    return e

def round_ref(v, mode):
    """Correctly rounded binary256 value of the reference v in mode, and flags."""
    if not isinstance(v, Fraction): return v, 0
    s = -1 if v < 0 else 1
    a = abs(v)
    E = exp2_floor(a)
    def rnd(Eg):
        q, r = divmod(a.numerator * 2 ** max(0, P - 1 - Eg), a.denominator * 2 ** max(0, Eg - (P - 1)))
        if r == 0: return q, True
        twice = 2 * r
        den = a.denominator * 2 ** max(0, Eg - (P - 1))
        up = {"RNE": twice > den or (twice == den and q & 1), "RNA": twice >= den, "RZ": False,
              "RU": s > 0, "RD": s < 0}[mode]
        return q + (1 if up else 0), False
    qu, _ = rnd(E)
    tiny = E < EMIN and not (qu == 1 << P and E + 1 >= EMIN)
    E2 = max(E, EMIN)
    q, exact = rnd(E2)
    flags = 0 if exact else F_INX
    if tiny and not exact: flags |= F_UNF
    if q == 0: return ("zero", s < 0), flags
    r = Fraction(q) * Fraction(2) ** (E2 - (P - 1))
    if r >= Fraction(2) ** (EMAX + 1):
        flags |= F_OVF | F_INX
        to_inf = mode in ("RNE", "RNA") or (mode == "RU" and s > 0) or (mode == "RD" and s < 0)
        if to_inf: return ("inf", s < 0), flags
        r = Fraction((1 << P) - 1) * Fraction(2) ** (EMAX - (P - 1))
    return s * r, flags

def same(got, want):
    """(equal value, equal sign of zero)"""
    if isinstance(got, Fraction) or isinstance(want, Fraction): return got == want, True
    if got[0] != want[0]: return False, True
    if got[0] == "zero": return True, got[1] == want[1]
    if got[0] == "inf": return got[1] == want[1], True
    return True, True

MAXF = Fraction((1 << P) - 1) * Fraction(2) ** (EMAX - (P - 1))  # largest finite

def ulp_err(got, want):
    """|got − want| in ulps of want (want the unrounded reference)."""
    if not isinstance(want, Fraction):
        if want[0] == "zero":
            return 0.0 if (not isinstance(got, Fraction) and got[0] == "zero") else float("inf")
        return 0.0 if same(got, want)[0] else float("inf")
    if not isinstance(got, Fraction):
        return float("inf") if got[0] != "zero" else float(abs(want) / Fraction(2) ** (max(exp2_floor(abs(want)), EMIN) - (P - 1)))
    e = max(exp2_floor(abs(want)), EMIN)
    return float(abs(got - want) / Fraction(2) ** (e - (P - 1)))

def near_grid(v, prec):
    """Is the mpfr value v within 2^200 units of 2^-prec (relative) of a multiple
    of the half-unit grid 2^(E-237)?"""
    if gmpy2.is_nan(v) or gmpy2.is_infinite(v) or v == 0: return False
    m, e = v.as_mantissa_exp()
    m, e = abs(int(m)), int(e)
    E = m.bit_length() + e - 1
    g = E - P
    if e >= g: return True
    t = g - e
    low = m & ((1 << t) - 1)
    d = min(low, (1 << t) - low)
    return d == 0 or (d.bit_length() + e) < E - prec + 200

def snap(v):
    """the grid point nearest to v (taken as the exact value)"""
    m, e = v.as_mantissa_exp()
    s = -1 if m < 0 else 1
    m, e = abs(int(m)), int(e)
    E = m.bit_length() + e - 1
    g = E - P
    if e >= g: return Fraction(s * m) * Fraction(2) ** e
    t = g - e
    q = (m + (1 << (t - 1))) >> t
    return Fraction(s * q) * Fraction(2) ** g

# ---------------------------------------------------------------- references

def Z(ops, i, p):
    set_prec(p)
    return mpc(to_mpfr(ops[i]), to_mpfr(ops[i + 1]))

def mp_ref(op, ops, p):
    """list of mpfr results at precision p"""
    set_prec(p)
    if op in ("exp", "log", "log10", "sqrt", "sin", "cos", "tan", "asin", "acos", "atan",
              "sinh", "cosh", "tanh", "asinh", "acosh", "atanh"):
        z = Z(ops, 0, p)
        r = getattr(gmpy2, op)(z)
        return [r.real, r.imag]
    if op == "abs": return [abs(Z(ops, 0, p))]
    if op == "arg": return [gmpy2.phase(Z(ops, 0, p))]
    if op == "pow":
        r = Z(ops, 0, p) ** Z(ops, 2, p)
        return [r.real, r.imag]
    if op == "powr":
        r = Z(ops, 0, p) ** to_mpfr(ops[2])
        return [r.real, r.imag]
    if op == "polar":
        set_prec(p)
        r = gmpy2.rect(to_mpfr(ops[0]), to_mpfr(ops[1]))
        return [r.real, r.imag]
    raise ValueError(op)

def exact_ref(op, ops):
    if op in ("mul", "div", "norm"):
        a, b = to_frac(ops[0]), to_frac(ops[1])
        if op == "norm": return [a * a + b * b]
        c, d = to_frac(ops[2]), to_frac(ops[3])
        if op == "mul": re, im = a * c - b * d, a * d + b * c
        else:
            den = c * c + d * d
            re, im = (a * c + b * d) / den, (b * c - a * d) / den
        return [re, im]
    return None


# ---------------------------------------------------------------- tiny arguments
# MPC rounds correctly at the requested precision, so for |z| < 2^-64 (where
# f(z) = lead + O(|z|^3) lies within 2^-(2·|log2 z|) of the exactly
# representable lead) it would need ~2·|log2 z| bits internally. There the
# reference is the exact lead (a Fraction) plus the Taylor rest, summed in MPC
# arithmetic at 3000 bits (relative to the rest) from rational coefficients.

TINY_OPS = ("exp", "sin", "cos", "tan", "sinh", "cosh", "tanh", "asin", "acos", "atan",
            "asinh", "acosh", "atanh")

def _series(z, first, ratio, stop_rel=3100):
    """Σ t_k with t_0 = first and t_{k+1} = ratio(k, t_k), until each component of
    the next term is 2^-stop_rel below that component of the sum (a component still
    0 waits for a nonzero term, up to 2^-(8·stop_rel) below the first term)."""
    total, t, k = first * 0, first, 0
    two = gmpy2.mpfr(2)
    floor = abs(first) * two ** (-8 * stop_rel)
    def done(sc, tc):
        return tc == 0 or (sc != 0 and abs(tc) < abs(sc) * two ** (-stop_rel))
    while t != 0:
        total += t
        t = ratio(k, t)
        k += 1
        if abs(t) < floor: break
        if isinstance(t, type(first)) and hasattr(t, "real") and hasattr(t, "imag") and not isinstance(t, type(two)):
            if done(total.real, t.real) and done(total.imag, t.imag): break
        elif done(total, t): break
    return total

def _frac(v):  # mpfr → Fraction (0 for zero)
    if v == 0: return Fraction(0)
    m, e = v.as_mantissa_exp()
    return Fraction(int(m)) * Fraction(2) ** int(e)

def tiny_ref(op, ops):
    if op == "polar":
        th = to_frac(ops[1])
        if ops[1][0] == 0 or exp2_floor(abs(th)) >= -64: return None
        set_prec(3000)
        rho, t = to_mpfr(ops[0]), to_mpfr(ops[1])
        cos_rest = _series(t, -t * t / 2, lambda k, u: -u * t * t / ((2 * k + 3) * (2 * k + 4)))
        sin_rest = _series(t, -t * t * t / 6, lambda k, u: -u * t * t / ((2 * k + 4) * (2 * k + 5)))
        r = to_frac(ops[0])
        return [r + _frac(rho * cos_rest), r * th + _frac(rho * sin_rest)]
    if op not in TINY_OPS: return None
    x, y = to_frac(ops[0]), to_frac(ops[1])
    for v, o in ((x, ops[0]), (y, ops[1])):
        if o[0] != 0 and exp2_floor(abs(v)) >= -64: return None
    if x == 0 or y == 0: return None  # the axes are special-cased (real functions)
    set_prec(3000)
    z = mpc(to_mpfr(ops[0]), to_mpfr(ops[1]))
    z2 = z * z
    alt = op in ("sin", "cos", "tan", "atan", "asinh")  # (−1)^k coefficients
    sg = -1 if alt else 1
    lead = (x, y)
    if op in ("sin", "sinh"):
        rest = _series(z, sg * z2 * z / 6, lambda k, u: sg * u * z2 / ((2 * k + 4) * (2 * k + 5)))
    elif op == "exp":
        lead = (1 + x, y)
        rest = _series(z, z2 / 2, lambda k, u: u * z / (k + 3))
    elif op in ("cos", "cosh"):
        zr2, zi2 = x * x - y * y, 2 * x * y  # z² exact
        lead = (1 + sg * zr2 / 2, sg * zi2 / 2)
        rest = _series(z, z2 * z2 / 24, lambda k, u: sg * u * z2 / ((2 * k + 5) * (2 * k + 6)))
    elif op in ("tan", "tanh"):
        # tan z − z = (sin z − z cos z)/cos z (and the hyperbolic analogue)
        num = _series(z, sg * z2 * z * (mpfr(1) / 6 - mpfr(1) / 2),
                      lambda k, u: sg * u * z2 * (mpfr(1) / gmpy2.fac(2 * k + 5) - mpfr(1) / gmpy2.fac(2 * k + 4))
                      / (mpfr(1) / gmpy2.fac(2 * k + 3) - mpfr(1) / gmpy2.fac(2 * k + 2)))
        den = 1 + _series(z, sg * z2 / 2, lambda k, u: sg * u * z2 / ((2 * k + 3) * (2 * k + 4)))
        rest = num / den
    elif op in ("atan", "atanh"):
        p3 = z2 * z
        rest = _series(z, sg * p3 / 3, lambda k, u: sg * u * z2 * (2 * k + 3) / (2 * k + 5))
    else:  # asin, asinh, acos, acosh: asin z − z = Σ a_k z^(2k+1)/(2k+1), a_k = Π (2j−1)/(2j)
        sa = -1 if op == "asinh" else 1
        rest = _series(z, sa * z2 * z / 6,
                       lambda k, u: sa * u * z2 * (2 * k + 3) * (2 * k + 3) / ((2 * k + 4) * (2 * k + 5)))
        if op in ("acos", "acosh"):
            half_pi = gmpy2.const_pi() / 2
            acos_re = mpfr_to_ref(half_pi - to_mpfr(ops[0]) - rest.real)  # π/2 − Re asin z: not near a lead
            acos_im = (-y, -rest.imag)
            if op == "acos": return [acos_re, acos_im[0] + _frac(acos_im[1])]
            # acosh z = ±i·acos z, Re >= 0: i·acos for Im z >= 0, −i·acos otherwise
            neg_y = ops[1][0] < 0
            re_lead, re_rest = (-acos_im[0], -acos_im[1]) if not neg_y else acos_im
            im = acos_re if not neg_y else -acos_re
            return [re_lead + _frac(re_rest), im]
    return [lead[0] + _frac(rest.real), lead[1] + _frac(rest.imag)]

def reference(op, ops):
    ex = exact_ref(op, ops)
    if ex is not None:
        # an exactly zero component: its sign is checked only through the
        # IEEE rule in the unit tests (Fraction 0 compares equal to ±0 here)
        return [v if v != 0 else ("zero", None) for v in ex]
    tr = tiny_ref(op, ops)
    if tr is not None: return tr
    vals = mp_ref(op, ops, 2000)
    if any(near_grid(v, 2000) for v in vals):
        vals = mp_ref(op, ops, 40000)
        return [snap(v) if near_grid(v, 40000) else mpfr_to_ref(v) for v in vals]
    return [mpfr_to_ref(v) for v in vals]

# ---------------------------------------------------------------- operand generators

def L(elo, ehi, signed=True):
    """±m·2^k, |value| in [2^elo, 2^(ehi+1)), 237-bit random m"""
    def g(rng):
        e = rng.randint(elo, ehi)
        m = rng.getrandbits(P - 1) | (1 << (P - 1))
        if signed and rng.random() < 0.5: m = -m
        return (m, e - (P - 1))
    return g

def U(lo, hi):
    def g(rng):
        x = Fraction(lo) + (Fraction(hi) - Fraction(lo)) * Fraction(rng.getrandbits(P), 1 << P)
        return from_frac(x) if x != 0 else (0, 0)
    return g

def small_int_bits(rng, bits):
    m = rng.getrandbits(bits) | 1
    return -m if rng.random() < 0.5 else m

def g_each(*gs):
    return lambda rng: [g(rng) for g in gs]

def g_mul_cancel(rng):
    """a·c − b·d nearly cancels: d = a·c/b truncated."""
    a, b, c = U(1, 4)(rng), U(1, 4)(rng), U(1, 4)(rng)
    d = from_frac(to_frac(a) * to_frac(c) / to_frac(b))
    return [a, b, c, d]

def g_div_cancel(rng):
    """b·c − a·d nearly cancels (Im of the quotient)."""
    a, b, c = U(1, 4)(rng), U(1, 4)(rng), U(1, 4)(rng)
    d = from_frac(to_frac(b) * to_frac(c) / to_frac(a))
    return [a, b, c, d]

def g_div_exact(rng):
    """z = w·q exactly (components of w, q: <= 100-bit integers times a common
    power of two each): the quotient q is exact."""
    kw, kq = rng.randint(-150, 50), rng.randint(-150, 50)
    w = [(small_int_bits(rng, 100), kw) for _ in range(2)]
    q = [(small_int_bits(rng, 100), kq) for _ in range(2)]
    wr, wi, qr, qi = map(to_frac, w + q)
    zr, zi = wr * qr - wi * qi, wr * qi + wi * qr
    if zr == 0 or zi == 0: return g_div_exact(rng)
    return [from_frac(zr), from_frac(zi), w[0], w[1]]

def g_abs_exact(rng):
    """(m(p² − q²), 2mpq)·2^k: |z| = m(p² + q²)·2^k exactly."""
    p, q = rng.getrandbits(50) | 1, rng.getrandbits(49) | 1
    m = rng.getrandbits(120) | 1
    k = rng.randint(-200, 200)
    a, b = m * (p * p - q * q), 2 * m * p * q
    if a == 0: return g_abs_exact(rng)
    return [(a, k), (b, k)]

def g_sqrt_exact(rng):
    """z = w² with w = (p, q)·2^k, p, q odd integers of <= 117 bits: sqrt(z) = w exactly."""
    k = rng.randint(-120, 80)
    a = Fraction(rng.getrandbits(117) | 1) * Fraction(2) ** k
    b = Fraction(small_int_bits(rng, 117)) * Fraction(2) ** k
    if rng.random() < 0.5: a = -a  # w = sqrt(z) has Re >= 0: z = (±a + ib)² covers both signs of x
    zr, zi = a * a - b * b, 2 * a * b
    if zr == 0: return g_sqrt_exact(rng)
    return [from_frac(zr), from_frac(zi)]

def g_unit_near(rng):
    """|z| = 1 + O(2^-237): e^{iθ}(1 + δ), components rounded"""
    set_prec(600)
    th = mpfr(rng.getrandbits(64)) / mpfr(2) ** 64 * 2 * gmpy2.const_pi()
    d = mpfr(rng.choice([-1, 1])) * gmpy2.mul_2exp(mpfr(rng.random()), -rng.randint(30, 400))
    x, y = gmpy2.cos(th) * (1 + d), gmpy2.sin(th) * (1 + d)
    return [from_frac(mpfr_to_ref(x)), from_frac(mpfr_to_ref(y))]

def g_branch(c):
    """c + small complex (c = ±1 or ±i): near the branch points"""
    def g(rng):
        cr, ci = c
        dr = L(-600, -2)(rng)
        di = L(-600, -2)(rng)
        x = Fraction(cr) + to_frac(dr)
        y = Fraction(ci) + to_frac(di)
        return [from_frac(x) if x else (0, 0), from_frac(y) if y else (0, 0)]
    return g

def g_cut(xg):
    """just above/below a branch cut on the real axis: tiny imaginary part"""
    return lambda rng: [xg(rng), L(-3000, -240)(rng)]

def g_int_w(lo, hi):
    def g(rng):
        n = 0
        while n == 0: n = rng.randint(lo, hi)
        return [U(-3, 3)(rng), U(-3, 3)(rng), (n, 0), (0, 0)]
    return g

def g_int_big(rng):
    """z near the unit circle, integer w in [65, 3000] (the exp·log path)"""
    z = g_unit_near(rng)
    n = rng.randint(65, 3000) * rng.choice([-1, 1])
    return z + [(n, 0), (0, 0)]

def g_powr(rng):
    return [U(-10, 10)(rng), U(-10, 10)(rng), U(-10, 10)(rng)]

A = U(-10, 10)
RANGES = [
    ("mul_gen", "mul", g_each(A, A, A, A)),
    ("mul_wide", "mul", g_each(*[L(-200000, 200000)] * 4)),
    ("mul_sub", "mul", g_each(*[L(-131180, -131060)] * 4)),
    ("mul_cancel", "mul", g_mul_cancel),
    ("div_gen", "div", g_each(A, A, A, A)),
    ("div_wide", "div", g_each(*[L(-100000, 100000)] * 4)),
    ("div_cancel", "div", g_div_cancel),
    ("div_exact", "div", g_div_exact),
    ("norm_gen", "norm", g_each(A, A)),
    ("norm_wide", "norm", g_each(L(-140000, 140000), L(-140000, 140000))),
    ("abs_gen", "abs", g_each(A, A)),
    ("abs_wide", "abs", g_each(L(-5000, 5000), L(-5000, 5000))),
    ("abs_exact", "abs", g_abs_exact),
    ("arg_gen", "arg", g_each(A, A)),
    ("sqrt_gen", "sqrt", g_each(A, A)),
    ("sqrt_wide", "sqrt", g_each(L(-5000, 5000), L(-5000, 5000))),
    ("sqrt_exact", "sqrt", g_sqrt_exact),
    ("sqrt_cut", "sqrt", g_cut(U(-10, -0.1))),
    ("exp_gen", "exp", g_each(A, A)),
    ("exp_small", "exp", g_each(L(-3000, -1), L(-3000, -1))),
    ("exp_big", "exp", g_each(U(-180000, 180000), A)),
    ("exp_ybig", "exp", g_each(U(-1, 1), L(10, 1000))),
    ("log_gen", "log", g_each(A, A)),
    ("log_unit", "log", g_unit_near),
    ("log_wide", "log", g_each(L(-200000, 200000), L(-200000, 200000))),
    ("log_cut", "log", g_cut(U(-10, -0.1))),
    ("log10_gen", "log10", g_each(A, A)),
    ("log10_unit", "log10", g_unit_near),
    ("pow_gen", "pow", g_each(A, A, U(-4, 4), U(-4, 4))),
    ("pow_int", "pow", g_int_w(-64, 64)),
    ("pow_intbig", "pow", g_int_big),
    ("pow_real", "powr", g_powr),
    ("pow_wide", "pow", g_each(L(-1000, 1000), L(-1000, 1000), U(-3, 3), U(-3, 3))),
    ("sin_gen", "sin", g_each(A, A)),
    ("sin_small", "sin", g_each(L(-3000, -1), L(-3000, -1))),
    ("sin_bigre", "sin", g_each(L(0, 1000), U(-1, 1))),
    ("cos_gen", "cos", g_each(A, A)),
    ("cos_small", "cos", g_each(L(-3000, -1), L(-3000, -1))),
    ("cos_bigre", "cos", g_each(L(0, 1000), U(-1, 1))),
    ("tan_gen", "tan", g_each(A, A)),
    ("tan_small", "tan", g_each(L(-3000, -1), L(-3000, -1))),
    ("tan_bigim", "tan", g_each(A, U(-3000, 3000))),
    ("sinh_gen", "sinh", g_each(A, A)),
    ("sinh_big", "sinh", g_each(U(-180000, 180000), A)),
    ("cosh_gen", "cosh", g_each(A, A)),
    ("cosh_small", "cosh", g_each(L(-3000, -1), L(-3000, -1))),
    ("cosh_big", "cosh", g_each(U(-180000, 180000), A)),
    ("tanh_gen", "tanh", g_each(A, A)),
    ("tanh_big", "tanh", g_each(L(0, 12), A)),
    ("asin_gen", "asin", g_each(A, A)),
    ("asin_branch", "asin", g_branch((1, 0))),
    ("asin_cut", "asin", g_cut(U(1.01, 10))),
    ("asin_wide", "asin", g_each(L(-100000, 100000), L(-100000, 100000))),
    ("acos_gen", "acos", g_each(A, A)),
    ("acos_branch", "acos", g_branch((-1, 0))),
    ("acos_cut", "acos", g_cut(U(-10, -1.01))),
    ("atan_gen", "atan", g_each(A, A)),
    ("atan_branch", "atan", g_branch((0, 1))),
    ("atan_wide", "atan", g_each(L(-100000, 100000), L(-100000, 100000))),
    ("asinh_gen", "asinh", g_each(A, A)),
    ("asinh_branch", "asinh", g_branch((0, -1))),
    ("acosh_gen", "acosh", g_each(A, A)),
    ("acosh_branch", "acosh", g_branch((1, 0))),
    ("acosh_cut", "acosh", g_cut(U(-10, 0.99))),
    ("atanh_gen", "atanh", g_each(A, A)),
    ("atanh_branch", "atanh", g_branch((1, 0))),
    ("atanh_cut", "atanh", g_cut(U(1.01, 10))),
    ("polar_gen", "polar", g_each(L(-1000, 1000, False), A)),
    ("polar_bigth", "polar", g_each(U(0.5, 2), L(0, 2000))),
]

# ---------------------------------------------------------------- driver

def selected(flt):
    shard = os.environ.get("SF_SHARD", "")
    si, sk = (int(t) for t in shard.split("/")) if shard else (0, 1)
    out = []
    for idx, rg in enumerate(RANGES):
        if idx % sk != si: continue
        if flt and not any(rg[0] == f or rg[0].startswith(f) for f in flt.split(",")): continue
        out.append(rg)
    return out

def range_rng(rng, seed, name):
    return random.Random(f"{seed}:{name}") if os.environ.get("SF_RANGE_SEED") == "1" else rng

def fmt(v):
    if isinstance(v, Fraction):
        return to_hex(from_frac(v)) if v.denominator & (v.denominator - 1) == 0 and abs(v.numerator).bit_length() <= P + 1 else f"{float(v):.6g}"
    if v[0] == "zero": return "-0" if v[1] else "0"
    if v[0] == "inf": return "-inf" if v[1] else "inf"
    return "nan"

def run(exe, n, seed, flt, modes):
    rng = random.Random(seed)
    env = dict(os.environ)
    if modes: env["CX_MODES"] = "1"
    maxfail = int(os.environ.get("SF_MAXFAIL", "50"))
    total_bad = 0
    crashed = False
    for name, op, gen in selected(flt):
        t0 = time.time()
        r = range_rng(rng, seed, name)
        cases = [gen(r) for _ in range(n)]
        inp = "".join(f"{op} {' '.join(to_hex(o) for o in ops)}\n" for ops in cases)
        pr = subprocess.run([exe], input=inp, capture_output=True, text=True, env=env)
        if pr.returncode != 0:
            print(f"CRASH range={name} op={op} rc={pr.returncode} seed={seed} stderr={pr.stderr.strip()[:300]!r}", flush=True)
            crashed = True
            continue
        lines = pr.stdout.split("\n")
        nmodes = 5 if modes else 1
        worst = [0.0] * nmodes
        mis = flagerr = zsign = shown = 0
        ovf = 0  # components whose exact value exceeds the largest finite: excluded from max ulp
        for ops, line in zip(cases, lines):
            refs = reference(op, ops)
            k = len(refs)
            tok = line.split()
            stride = k + (1 if modes else 0)
            for mi in range(nmodes):
                m = MODES[mi] if modes else "RNE"
                got = [parse_hex(tok[mi * stride + j]) for j in range(k)]
                gflags = int(tok[mi * stride + k]) & 0x1F if modes else 0
                wflags = 0
                bad = []
                for j in range(k):
                    w = refs[j]
                    if not isinstance(w, Fraction) and w[0] == "zero" and w[1] is None:
                        want = ("zero", m == "RD") if op in ("mul", "div") else ("zero", False)
                        fl = 0
                        eq, sg = isinstance(got[j], tuple) and got[j][0] == "zero", True
                    else:
                        want, fl = round_ref(w, m)
                        eq, sg = same(got[j], want)
                    wflags |= fl
                    if isinstance(w, Fraction) and abs(w) > MAXF:  # overflow: ±inf or ±max finite
                        ovf += 1
                        e = 0.0 if eq else float("inf")
                    else:
                        e = ulp_err(got[j], w if not (isinstance(w, tuple) and w[1] is None) else ("zero", False))
                    if e > worst[mi]: worst[mi] = e
                    if not eq: mis += 1; bad.append(f"c{j} got={fmt(got[j])} want={fmt(want)} ulp={e:.3g}")
                    elif not sg: zsign += 1; bad.append(f"c{j} zero sign got={fmt(got[j])} want={fmt(want)}")
                if modes and gflags != wflags:
                    flagerr += 1
                    bad.append(f"flags got={gflags:#x} want={wflags:#x}")
                if bad and shown < maxfail:
                    shown += 1
                    print(f"FAIL range={name} op={op} in={','.join(to_hex(o) for o in ops)} mode={m} " + "; ".join(bad) + f" seed={seed}", flush=True)
        total_bad += mis + flagerr + zsign
        tm = f"  [{time.time() - t0:.0f} s]" if os.environ.get("SF_TIMING") == "1" else ""
        if modes:
            print(f"{name:13s} n={n} max ulp " + " ".join(f"{MODES[i]} {worst[i]:.3g}" for i in range(5)) +
                  f" | misrounded {mis}, flag errors {flagerr}, zero-sign errors {zsign}" +
                  (f", overflowed {ovf}" if ovf else "") + tm, flush=True)
        else:
            print(f"{name:13s} n={n} max ulp {worst[0]:.3g}  misrounded {mis}" + (f", zero-sign errors {zsign}" if zsign else "") +
                  (f", overflowed {ovf}" if ovf else "") + tm, flush=True)
    print(f"SUMMARY seed={seed} cases={n} {'modes' if modes else 'RNE'} errors={total_bad}" + (" CRASH" if crashed else ""), flush=True)
    sys.exit(2 if crashed else (1 if total_bad else 0))

def main():
    if len(sys.argv) > 1 and sys.argv[1] == "--list":
        print("\n".join(r[0] for r in RANGES)); return
    exe = sys.argv[1]
    n = int(sys.argv[2]) if len(sys.argv) > 2 else 40
    seed = int(sys.argv[3]) if len(sys.argv) > 3 else 1
    flt = sys.argv[4] if len(sys.argv) > 4 else ""
    modes = len(sys.argv) > 5 and sys.argv[5] == "modes"
    run(exe, n, seed, flt, modes)

if __name__ == "__main__":
    main()
