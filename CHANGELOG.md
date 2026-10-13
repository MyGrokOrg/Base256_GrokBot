# Changelog

## 4.6.5

- **`log1p`** of a tiny argument is the argument, not zero. The old cutoff was `2^{-245760}` instead of `2^{-16}`, so `log1p(2^{-300})` went through `log(1+x)` and flushed. `log1p(±0)` keeps its sign.
- **`pow`** follows the IEEE special cases: `pow(x, ±0) = 1` and `pow(+1, y) = 1` even when the other argument is NaN, and `pow(−1, ±∞) = 1`. Infinities as a base pick up the odd-integer sign. Exponents that are integers with absolute value at most 4096 are computed by squaring, so `pow(2, 100)` and `pow(10, 20)` are exact.
- **`tgamma` / `lgamma`**. Positive integers through 54 use an exact factorial (`53!` still fits in 237 bits). Everything else is a 40-term Stirling series shifted up to 48, with the reflection formula on the negative axis. `tgamma(±0)` is `±∞`; negative integers are NaN for `tgamma` and `+∞` for `lgamma`.
- **`sincos`** evaluates sine and cosine from one argument reduction. **`remquo`** is `remainder` plus the low 31 bits of the quotient.

## 4.6.4

- **`<base256/base256.hpp>` compiles again.** `logb` was defined in both `float256.hpp` and `extra_math.hpp`, so the umbrella header was ill-formed. The constexpr `float256.hpp` version is the one that remains.
- **`tanh`** of a large finite argument is `±1`, not NaN. `exp(2x)` overflows past `|x| ≈ 9·10^4`, and the old `(e^{2x}-1)/(e^{2x}+1)` became `inf/inf`.
- **`sinh` / `cosh`** near zero go through `expm1`, so `sinh(2^{-200})` stays `2^{-200}` instead of flushing to zero.
- **`exp2`** of an integer is the exact power of two (`exp2(10) == 1024`). **`log2`** of a power of two is the exact integer (`log2(2^{-80}) == -80`).
- **`erf` / `erfc`**, odd/even specials included (`erf(±0)`, `erf(±∞)`, `erfc(-∞) = 2`). They agree with `erf(x) + erfc(x) = 1`.

## 4.6.3

- **`asinh` / `acosh` / `atanh`** with large-argument `log(2|x|)` reduction so they stay finite past `2^64`.
- **`fdim(x, y)`** returns `max(x-y, +0)`. **`logb`** is the floating `ilogb`. **`nearbyint`** matches `rint` (ties to even).

## 4.6.2

- **`pow(±0, y)`** keeps the sign of zero when `y` is an odd integer. `pow(−0, 3)` is `−0`; `pow(−0, −1)` is `−∞`.
- **`atan2`** on infinities returns `±π/4` or `±3π/4`. `atan2(y, +∞)` is `±0` and `atan2(y, −∞)` is `±π`.
- **`sin` / `cos` / `tan`** at exact multiples of this library's π no longer produce a spurious `−0`. `sin(π)` is `+0`, `cos(π/2)` is `+0`, `tan(π/2)` is `+∞`.
- **`rint`** rounds to nearest, ties to even. **`modf`** splits a value into a signed integer and a signed fraction. **`nextup`** / **`nextdown`** wrap `nextafter`. **`remainder`** uses ties-to-even.
- **`hypot(∞, NaN)`** returns `+∞`, matching IEEE 754.

## 4.6.1

- **`to_string`** no longer trusts a double `log10` estimate that lands on the wrong side of a power of ten (`0.01` used to print as `1e-1`).
- **`hypot`** scales into `[0.5, 1)` before squaring, so large finite inputs stay finite.
- **`fmin` / `fmax`** of two zeros follow IEEE 754.

## 4.6.0

Initial release: header-only IEEE 754 binary256 with correctly rounded `+ − × ÷`, FMA, `sqrt`, decimal and hex I/O, and elementary functions.
