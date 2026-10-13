# Changelog

## 4.9.0

- **Shortest round-trip `to_string()`**: with no precision it prints the fewest digits that parse back to the same value. Candidates are checked against the exact rounding interval with big-integer arithmetic. `to_string(n)` rounds to `n` significant digits in the current mode.
- **`format(std::chars_format, precision)`** supports fixed, scientific, general (`%g`-style) and hex output. Every format is exact at any precision and honours the rounding mode.
- **`base256::to_chars` / `base256::from_chars`** follow the `<charconv>` contract: `errc::value_too_large` for a short buffer, `errc::result_out_of_range` on overflow or underflow, and the format restrictions. Parsing is correctly rounded.
- **`std::formatter<float256>`** handles fill, alignment, sign, `#`, `0`, width, precision and the types `a A e E f F g G`.
- **Integer conversions**: every integer type, including `__int128`, converts in. `to_integer<T>(mode, exact)` is IEEE `convertToInteger`; it raises invalid and saturates when the value is out of range. `total_order` and `total_order_mag` are added.
- `BigInt` has Knuth Algorithm D division and a decimal conversion that works in chunks of 10^9. Decimal output at the extreme exponents takes about 35 ms.
- Parsing clamps absurd decimal exponents instead of overflowing `int`.

## 4.8.0

- **Rounding modes**: all five IEEE 754 directions (`base256::round_mode`), per thread, via `base256::fenv::{get_round, set_round, scoped_round}`. They apply to arithmetic, `fma`, `sqrt`, parsing, the narrowing conversions and `rint`/`nearbyint`. New `round_to_integral(x, mode)`.
- **Exception flags**: invalid, divide-by-zero, overflow, underflow (tininess after rounding) and inexact, through `fenv::{test_flags, clear_flags, raise_flags, get_flags, set_flags}`. Ordered comparisons signal invalid on NaN, while `==` signals only for a signaling NaN.
- An exact zero sum is `−0` when rounding downward. Overflow honours the direction (`max()` toward zero).
- **`to_float`** rounds once instead of going through `double`. `to_double`, `to_float` and `to_long_double` share one mode-aware narrowing routine.
- **`fmod`, `remainder`, `remquo`** are exact (shift-and-subtract on the significands) and `constexpr`. They were wrong once the quotient exceeded 237 bits.
- Elementary functions evaluate in round-to-nearest and leave the caller's flags unchanged.
- `numeric_limits<float256>::is_iec559` is now `false`, with a conformance table in the README. The previous `true` overstated conformance.
- The MPFR differential test covers all modes and flags, plus `fmod` and `remainder`.

## 4.7.0 (Base256_GrokBot)

Correctness fixes found with the new MPFR differential test.

- **Effective subtraction** (`a − b` with opposite-sign `+`, and `fma`) misrounded when bits of the smaller operand were shifted out: results just below a rounding midpoint were rounded up. Example: `1 − (2^-238 + 2^-400)` returned `1` instead of `nextdown(1)`.
- **Subnormal results** that rounded up across a power of two came back as the minimum normal (`denorm_min × 1.5` was `min()`), and a guard bit at position 64 evaluated `1 << 64` (undefined behaviour). `to_double` rounded binary64 subnormals twice.
- **`fma((±0), y, (∓0))`** is `+0`.
- **Numeric literals** (`0.1_f256`) are parsed exactly. The `long double` literal overload is gone, so numeric literals are no longer `constexpr`. Digit separators are accepted.
- **`float256(long double)`** is exact (x87 80-bit, binary128). New `to_long_double()` rounds once.
- The header no longer `#undef`s `isnan`, `signbit`, `fpclassify` and the rest.
- **Signaling NaNs**: `issignaling()` and `quieted()`. `fmin`/`fmax` return a quiet NaN for an sNaN operand. New IEEE 754-2019 `minimum`, `maximum`, `minimum_number`, `maximum_number`. The NaN payload rules are documented in the README.
- **Testing**: `tests/mpfr_diff.cpp` (`-DBASE256_MPFR_TESTS=ON`) compares `+ − × ÷ fma sqrt`, `to_double` and `to_long_double` bit-for-bit with MPFR. CI adds MPFR and ASan/UBSan jobs.

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
