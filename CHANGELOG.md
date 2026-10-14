# Changelog

## 4.13.0

- **Speed: elementary functions 10–21× faster, `/` 14×, `sqrt` 5×** with no accuracy change: the full `mpfr_ulp` suite in all five modes still shows 0 misrounded results and 0 flag errors.
  - `divmod_512_256`: Knuth algorithm D on 64-bit limbs; `isqrt_512`: a 64-bit seed plus Newton steps. The bitwise versions remain as the fallback when `__int128` is missing.
  - `round_bits` works on whole limbs.
  - dd `two_sum`/`two_prod` work directly on the 512-bit significands and give results identical to the float formulas. dd `+`/`*` accumulate in a 512-bit window and round once, which is also slightly more accurate (about 2^-505).
  - The Taylor series for `expm1` and `sin`/`cos` use Horner on a precomputed 1/n! table, with the number of terms chosen from the exponent. Division by small integers uses `dd / float256`.
- New `base256_bench` target. See the README section "Performance" for the table (e.g. `sin` 450 µs → 23 µs, `exp` 328 µs → 16 µs).

## 4.12.0

- **Elementary functions honor the rounding mode and raise IEEE flags.** All of `math.hpp` and `extra_math.hpp` round their ~474-bit dd result once in the current mode (new `float256::round_dd(hi, lo, k, mode, tail)`, with a tail hint for results such as `tanh → ±1` or `expm1 → -1` where the correction is below the dd precision). Before, they always rounded to nearest and raised no flags.
  - `fe_invalid`: domain errors (`log(-1)`, `asin(2)`, `sin(inf)`, `pow(-2, 0.5)`, `tgamma(-3)`, …) and signaling-NaN inputs (the result is the quieted NaN).
  - `fe_divbyzero`: poles — `log*(±0)`, `log1p(-1)`, `atanh(±1)`, `lgamma(0, -1, -2, …)`, `tgamma(±0)`, `pow(±0, y<0)`.
  - `fe_overflow`/`fe_underflow`: library convention (tininess after rounding; underflow only with inexact); directed modes saturate to `max()`/`denorm_min()` as IEEE requires.
  - `fe_inexact`: on every inexact result, never on exact ones (`exp(0)`, `exp2(n)`, `log(1)`, `log2(2^n)`, `log10(10^n)`, `sin(±0)`, exact `pow`, `cbrt` of perfect cubes, `tgamma` of small integers, `lgamma(1|2)`).
  - Internal steps run in a private environment (`detail::FnEnv`), so no intermediate flag leaks; the caller's flags are kept.
- `mpfr_ulp` evaluates every case in all five modes and compares with MPFR rounded in that mode (binary256 range emulated), checking direction, the 1-ulp bound and the exact flag set. New overflow/underflow ranges (`exp_big`, `exp2_big`, `expm1_big`, `sinh_big`, `cosh_big`). Result: 0 misrounded values and 0 flag errors over 265k mode-evaluations.
- New test `tests/test_math_flags.cpp` (ctest `math_flags`).
- `test_fenv` now expects elementary functions to follow the mode.

## 4.11.0

- **`sin`, `cos`, `tan`, `sincos` are accurate for every finite argument.** For |x| ≥ 2^100 they now use Payne–Hanek reduction: the significand is multiplied exactly by a window of 2/π bits, leaving 1216 guard bits. The bits come from a 264,192-bit table generated with MPFR (`detail/two_over_pi.hpp`, 4128 `uint64_t` words, about 100 KB). Before this, results above 2^100 could be wrong in every digit. Now they are under 1 ulp up to `max()`. The old float256 path is gone; `reduce_trig` is kept only for compatibility.
- `mpfr_ulp` adds tests across the whole exponent range, kπ/2 for k up to 2^262140, and continued-fraction hard cases: arguments within about 2^-238 of a multiple of π/2. Regression tests cover `sin`/`cos`/`tan` of `max()` and two hard cases near the top of the range.
- `BigInt` gains `from_limbs`, `limb` and `keep_low_limbs`.

## 4.10.0

- **Elementary functions to < 1 ulp.** Every function in `math.hpp` and `extra_math.hpp` now evaluates in double-float256 (`detail/dd.hpp`, about 474 bits) and rounds once. Before, errors ranged from a few ulps to completely wrong results.
  - `exp`: up to 6e4 ulp for large |x|.
  - `erfc`: about 4e10 ulp near x = 3.6.
  - `asinh`/`acosh` above 2^64: about 1e30 ulp.
  - `sin`/`cos`/`tan` above 1e17: `inf` ulp (wrong in every digit).
  - `lgamma`: about 1e5 ulp near 1.
- Changes by function:
  - `exp`, `exp2`, `expm1`: reduction by a 3-part ln 2, then a Taylor series on r/2^12 with 12 squarings of `expm1`.
  - `log` family: one Newton step on `exp` from a float256 seed.
  - `pow`: `exp(y·log x)` in dd. Integer exponents up to 4096 use dd squaring (exact when representable).
  - `sin`, `cos`, `tan`, `sincos`: Cody–Waite reduction with a 1190-bit π/2 for |x| < 2^100.
  - `atan`, `atan2`, `asin`, `acos`: a dd correction step `θ0 + atan(t)`.
  - Hyperbolic functions and their inverses: built on dd `exp`, `expm1` and `log`.
  - `erf`, `erfc`: a positive-term series below 4 and the Laplace continued fraction above, including the subnormal tail.
  - `lgamma`, `tgamma`: a 60-term dd Stirling series (shift to 128), the reflection formula, and Taylor series at 1 and 2.
  - `cbrt`: a dd residual correction.
- `sin(pi())` is now the true residual (≈ 1.2e-72) rather than `+0`. `pi()` is not π.
- **New test** `tests/mpfr_ulp.cpp` (ctest `mpfr_ulp`, with `-DBASE256_MPFR_TESTS=ON`) measures ulp error against MPFR and fails on any result over 1 ulp. The README has an "Elementary function accuracy" section.

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
