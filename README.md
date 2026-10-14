# Base256

**IEEE 754 binary256** (octuple-precision) floating-point for C++23.

Header-only. Real 256-bit software arithmetic — not a `double` wrapper.

| Field | Width |
|---|---|
| Sign | 1 bit |
| Exponent | 19 bits (bias `262143`) |
| Trailing significand | 236 bits |
| Precision *p* | **237 bits** (~71 decimal digits) |
| *E*<sub>min</sub> … *E*<sub>max</sub> | −262142 … +262143 |
| Max finite | ≈ 1.61 × 10<sup>78913</sup> |
| Min positive subnormal | ≈ 2.25 × 10<sup>−78984</sup> |

```
  255  254            236  235                                                    0
   [s][     exponent     ][                  trailing significand                 ]
    1          19                                   236
```

## Features

- Correctly rounded **`+ − × ÷`**, **FMA**, **`sqrt`**, exact **`fmod`/`remainder`**, all checked bit-for-bit against MPFR
- All five IEEE 754 rounding modes and the five exception flags (`base256::fenv`)
- IEEE 754 specials: `±0`, `±∞`, quiet/signaling NaN, subnormals
- Unordered NaN comparisons (`std::partial_ordering`) and `+0 == −0`
- Conversions to/from `float`, `double`, integers, decimal, and C99 hex-float
- `std::numeric_limits<float256>` (`is_iec559` is `false` for now; see [IEEE 754 conformance](#ieee-754-conformance))
- Elementary functions: `exp`, `log`, `pow`, `sin`/`cos`/`tan`, `atan`/`asin`/`acos`, hyperbolics and inverse hyperbolics, `cbrt`, `erf`/`erfc`, `tgamma`/`lgamma`
- Special functions (`special_math.hpp`, C++17 names): `cyl_bessel_j/i/k`, `cyl_neumann`, `sph_bessel`, `sph_neumann`, `jn`/`yn`, `beta`, `ibeta`/`ibetac`, `riemann_zeta`, `hurwitz_zeta`, `tgamma_lower`/`tgamma_upper`, `gamma_p`/`gamma_q`, `comp_ellint_1/2/3`, `ellint_1/2/3`
- `sincos` and IEEE `remquo`
- `ldexp`, `frexp`, `ilogb` / `logb`, `nextafter` / `nextup` / `nextdown`, `floor`/`ceil`/`trunc`/`round` (ties away), `rint` / `nearbyint` (ties to even), `modf`, `hypot`, `fmod`, `remainder`, `fdim`
- User-defined literals `3.1415_f256` and `"0.1"_f256`, both parsed exactly (correctly rounded)
- Exact `float256(long double)` (x87 80-bit / binary128) and correctly rounded `to_long_double()`

## Install

Header-only. Point your include path at `include/`:

```bash
git clone https://github.com/MyGrokOrg/Base256_GrokBot.git
cmake -B build -S Base256_GrokBot
cmake --build build
ctest --test-dir build --output-on-failure
```

```cmake
add_subdirectory(Base256_GrokBot)
target_link_libraries(your_app PRIVATE Base256::base256)
```

Requires a C++23 compiler (GCC 12+, Clang 16+, MSVC 19.4+).

## Usage

```cpp
#include <base256/base256.hpp>
#include <iostream>

int main() {
    using namespace base256;
    using namespace base256::literals;

    // 1 + 2^-80 vanishes in binary64; binary256 keeps it.
    auto one   = float256::one();
    auto delta = ldexp(one, -80);
    std::cout << std::boolalpha
              << ((one + delta) == one) << "\n";          // false

    auto third = 1.0_f256 / 3.0_f256;
    std::cout << (third * 3.0_f256 == one) << "\n";       // true

    auto a = "0.1"_f256, b = "0.2"_f256, c = "0.3"_f256;
    std::cout << abs((a + b) - c).to_string(10) << "\n";  // ~1e-72, not ~1e-17

    std::cout << pi().to_string(50) << "\n";
    std::cout << sqrt(2.0_f256).to_hex_string() << "\n";
}
```

Include `<base256/float256.hpp>` for the type alone, or `<base256/base256.hpp>` for the type plus elementary functions.

## Representation

Words are little-endian `std::array<std::uint64_t,4>` (`[0]` is the least-significant 64 bits of the 256-bit encoding):

- Word 3: `[sign:1][exponent:19][fraction 235:192]`
- Word 2: fraction bits 191:128
- Word 1: fraction bits 127:64
- Word 0: fraction bits 63:0

`raw_bits()` / `from_bits()` expose this packing. Hex-float I/O uses C99 `0x1.ffffp+e` with a 236-bit fraction.

## NaN rules

- The quiet bit is the most significant fraction bit (bit 235). `issignaling()` tests it; `quieted()` sets it.
- Arithmetic (`+ − × ÷`, `fma`, `sqrt`, conversions) on NaN operands returns the **first NaN operand** (in argument order) with its payload kept and the quiet bit set.
- Invalid operations (`∞ − ∞`, `0 × ∞`, `0/0`, `sqrt(−x)`, …) return the default NaN `float256::qnan()`: positive, quiet, payload `1`.
- `-x`, `abs`, `copysign` only touch the sign bit and never quiet a signaling NaN.
- `fmin`/`fmax` (C23 / IEEE 754-2008 minNum/maxNum) ignore a quiet NaN but return a quiet NaN for a signaling NaN. `minimum`/`maximum` (IEEE 754-2019) propagate any NaN; `minimum_number`/`maximum_number` ignore NaNs.
- Converting a NaN to or from `double`/`long double` does not keep the payload.

## Testing

Unit tests run with `ctest`. A differential test against GNU MPFR (237-bit precision, binary256 exponent range, `mpfr_subnormalize`) with randomized operands is available:

```bash
cmake -B build -DBASE256_MPFR_TESTS=ON -DBASE256_MPFR_CASES=1000000
cmake --build build && ctest --test-dir build --output-on-failure
./build/mpfr_diff 5000000 42   # cases per op, seed
./build/mpfr_ulp 5000 42 1     # elementary-function ulp check: cases per function, seed, fail over bound
./build/mpfr_ulp 5000 42 0 "sin cos"   # report only, selected functions
python3 tools/sf_check.py build/sf_eval 200 1   # special functions vs mpmath (pip install mpmath)
```

## Rounding modes and exception flags

```cpp
#include <base256/float256.hpp>
using namespace base256;

fenv::set_round(round_mode::upward);          // or scoped: fenv::scoped_round g(round_mode::downward);
fenv::clear_flags();
float256 q = float256::one() / float256(3.0);
if (fenv::test_flags(fe_inexact)) { /* ... */ }
fenv::set_round(round_mode::to_nearest_even);
```

- Modes: `to_nearest_even` (default), `to_nearest_away`, `toward_zero`, `upward`, `downward`. The environment is per thread.
- Flags: `fe_invalid`, `fe_divbyzero`, `fe_overflow`, `fe_underflow`, `fe_inexact` (`fe_all`). Underflow uses tininess **after** rounding and is raised only when the result is also inexact.
- The mode applies to `+ − × ÷`, `fma`, `sqrt`, `parse`/literals, `to_double`/`to_float`/`to_long_double`, and `rint`/`nearbyint`. `round_to_integral(x, mode)` takes an explicit direction.
- Overflow gives `±∞`, or `±max()` when the direction is toward zero. An exact zero sum is `+0`, or `−0` in `downward` mode.
- **Constant evaluation** (`constexpr`) always uses `to_nearest_even` and records no flags.
- Elementary functions (`exp`, `sin`, `pow`, …) honor the current rounding mode: the ~474-bit internal result is rounded once, in that mode. They raise `fe_invalid` for domain errors and signaling NaNs, `fe_divbyzero` at poles (`log(0)`, `atanh(±1)`, `lgamma` at non-positive integers, `pow(0, y<0)`, `tgamma(±0)`), `fe_overflow`/`fe_underflow` by the convention above, and `fe_inexact` on every inexact result (never on exact ones such as `exp(0)`, `log(1)`, `sin(±0)`, exact `pow`, `cbrt(27)`, `tgamma(5)`). Internal steps leave no flags behind (see [Elementary function accuracy](#elementary-function-accuracy)).

## IEEE 754 conformance

| Area | Status |
|---|---|
| `+ − × ÷`, `fma`, `sqrt`, `fmod`, `remainder` | Correctly rounded or exact in all 5 modes, flags included. `+ − ×` take a fast path for normal operands and results; it is bit-identical to the general path (ctest `fastpath`; define `BASE256_NO_FAST_PATHS` to disable it). Verified bit-for-bit against MPFR (237-bit, subnormalized) on randomized and edge-case operands |
| Subnormals, signed zeros, ±∞, NaN payloads | Supported (see NaN rules) |
| Comparisons | `==`/`!=` quiet, `<`/`<=`/`>`/`>=` signal invalid on NaN; `std::partial_ordering` |
| Conversions from/to `double`, `float`, `long double`, integers | Exact or correctly rounded in the current mode |
| Decimal/hex parsing, literals | Correctly rounded in the current mode |
| Decimal output (`to_string`, `format`, `to_chars`, `std::format`) | Shortest round-trip, or exact at any precision and correctly rounded in the current mode |
| `totalOrder`, `convertToInteger` | `total_order`, `total_order_mag`, `to_integer<T>(mode, exact)` |
| Elementary functions | All 5 modes, flags included. Under 1 ulp in every mode, correct direction in directed modes; no misrounding observed against MPFR (correct rounding not guaranteed) |

Elementary functions are not guaranteed correctly rounded, so `std::numeric_limits<float256>::is_iec559` stays `false`.

## Elementary function accuracy

Every elementary function evaluates in double-float256 (`detail/dd.hpp`, an unevaluated sum `hi + lo` of about 474 bits). The error-free sums and products come from the correctly rounded float256 `+` and `fma`, and the result is rounded to binary256 once. Trigonometric reduction is accurate for every finite argument. Below 2^100 it uses Cody–Waite with a 1190-bit π/2. From 2^100 up to `max()` it uses Payne–Hanek with a 264,192-bit table of 2/π (`detail/two_over_pi.hpp`, 4128 `uint64_t` words): the 237-bit significand is multiplied exactly by the needed 2/π window, keeping 1216 fraction bits. Constants (ln 2, π/2 in five parts, ln 10, √π, ½ ln 2π, 140 Stirling coefficients, γ, ζ(2), ζ(3)) were generated with MPFR.

**Documented bound: < 1 ulp** (round-to-nearest) for all of `exp exp2 expm1 log log2 log10 log1p pow cbrt sin cos tan sincos atan atan2 asin acos sinh cosh tanh asinh acosh atanh erf erfc tgamma lgamma`, over the full binary256 range, with the exceptions below. The evaluation error is about 2^-460 relative, so in practice results are correctly rounded except when the exact value lies within roughly 2^-220 ulp of a rounding midpoint. `tests/mpfr_ulp.cpp` checks this against a 1100-bit MPFR reference on random, log-uniform arguments plus targeted hard regions: near kπ/2 for k up to 2^262140, sin/cos/tan over the whole exponent range, continued-fraction worst cases within ~2^-238 of a multiple of π/2 for every exponent, near ±1 for `asin`/`acos`/`atanh`/`acosh`, near 1 and 2 for `lgamma`, and negative arguments for the gamma functions. With 4 seeds × 5000 cases per function (≈ 840k evaluations), the worst observed error is 0.500 ulp for every function and none exceeded the bound. CI runs the check on every push.

Exceptions and limits:
- `lgamma` within ~2^-220 of its non-trivial zeros on the negative axis (and `tgamma` there, where it is ±1) has an absolute error of about 2^-460, not a relative bound.
- `float256::pi()` is π rounded, so `sin(pi())` is the residual `π − pi()` (≈ 1.2e-72), not zero. The same goes for `cos(pi()/2)` and `tan(pi())`.
- Since 4.12.0 every case is checked in all five rounding modes against MPFR rounded in that mode (237-bit, subnormalized, binary256 exponent range), with the flags: inexact iff the result is inexact, overflow iff MPFR overflows, underflow iff tiny after rounding and inexact, never invalid/divbyzero for these arguments. Overflow and gradual-underflow ranges of `exp`, `exp2`, `expm1`, `sinh`, `cosh` are included. Over 53 checks × 1000 cases × 5 modes (265k evaluations), 0 results differed from MPFR's correctly rounded value and 0 flag mismatches. `tests/test_math_flags.cpp` covers exact cases, poles, domain errors, sNaN, overflow/underflow per mode and flag hygiene.


## Special functions

`#include <base256/special_math.hpp>`, which `base256.hpp` also includes. These functions use the C++17 `<cmath>` special-math names and argument order, such as `cyl_bessel_j(ν, x)`, `ellint_3(k, ν, φ)` and `sph_bessel(n, x)`. Boost names are used where `<cmath>` has none: `ibeta`, `ibetac`, `tgamma_lower`, `tgamma_upper`, `gamma_p`, `gamma_q`.

How they behave:
- Like the elementary functions, each one is evaluated in double-float256 and rounded once in the current mode.
- They raise `inexact`, `overflow` and `underflow` like those functions do.
- A domain error returns NaN with `invalid`. A pole returns ±∞ with `divbyzero`. A NaN argument propagates, and a signaling NaN also raises `invalid`.

**Measured accuracy.** In every case below, the round-to-nearest result was the correctly rounded value (at most 0.5 ulp). The directed modes were checked against MPFR wherever MPFR has the function: the error is at most 1 ulp, the direction is correct and the flags are exact.

| Function | Reference | Ranges checked (cases) | Max error |
|---|---|---|---|
| `jn`, `yn`, `cyl_bessel_j`, `cyl_neumann` (integer ν) | MPFR `jn`/`yn`, all 5 modes | n ∈ [−40, 40], \|x\| ∈ [2^-12, 2^9]; ν ∈ [0, 30], x ∈ [2^-4, 2^11] (1500 each) | 0.5 ulp RN, 1 ulp directed |
| `cyl_bessel_j/i/k`, `cyl_neumann` (real ν) | mpmath, 700 bits | ν ∈ [0, 1000], x ∈ [2^-400, 2^200], incl. ν ≥ x, ν within 2^-20…2^-200 of an integer (200 per range, 26 ranges) | 0.5 ulp |
| `sph_bessel`, `sph_neumann` | mpmath | n ≤ 20, x ≤ 50 | 0.5 ulp |
| `beta` | MPFR `beta`, all modes; mpmath | x, y ∈ [2^-200, 2^60] | 0.5 ulp RN, 1 ulp directed |
| `ibeta`, `ibetac` | mpmath at 2000 bits | a, b ∈ [2^-60, 5000], x ∈ [2^-60, 1] | 0.5 ulp |
| `riemann_zeta` | MPFR `zeta`, all modes; mpmath | \|x\| ∈ [2^-300, 2^14], x within 2^-236 of 1, near the trivial zeros | 0.5 ulp RN, 1 ulp directed |
| `hurwitz_zeta` | mpmath at 3000 bits | s ∈ [−2000, 3000], a ∈ [2^-20, 2^300], s within 0.01 of 1 (100 per range, 8 ranges) | 0.5 ulp |
| `tgamma_lower/upper`, `gamma_p/q` | MPFR `gamma_inc` (Γ(a) at up to ~5000 bits), all modes; mpmath | a ∈ [2^-200, 3000], x ∈ [2^-12, 2^9] | 0.5 ulp RN, 1 ulp directed |
| `comp_ellint_1/2/3`, `ellint_1/2/3` | mpmath | \|k\| ≤ 1 (incl. 0.999…1), ν ∈ [−10, 1), \|φ\| ≤ 100 | 0.5 ulp |

Over 1500 cases × 5 modes per MPFR-checked function, no result was misrounded. mpmath has no binary256 rounding modes, so the mpmath-only functions were checked in round-to-nearest. Their directed modes are covered by the reference values in `tests/test_special_math.cpp`. As with the elementary functions, these results are measured, not guaranteed.

Methods:
- **Bessel:** power series for x < 120. For Y and K, Temme's series below x = 56 (Y) and 30 (K). Above that, up to 120, Y uses division-free Steed CF1/CF2, and K uses the I series with the Wronskian. The Temme gamma factors come from an MPFR Taylor table of 1/Γ(1+μ). Terms below 2^-245 of a sum are finished in float256. Hankel and large-x expansions from 120 up, with forward recurrence (Y, K, and J for ν < x) or CF1 plus downward recurrence (J for ν ≥ x).
- **Zeta:** Euler–Maclaurin for x ≥ 0 at the smallest N in 24…64 whose remainder estimate is below 2^-506, with B_2k/(2k)! and ln p from MPFR tables. The functional equation handles x < 0. `hurwitz_zeta` picks by estimated cost between a direct sum and N terms plus Euler–Maclaurin at N + a. For s < 1 it uses the analytic continuation of the same expansion, and Hurwitz's formula for s ≤ −89.
- **Incomplete gamma:** the series for x < max(a + 1, 36), otherwise a division-free (Wallis) continued fraction. A small-a form handles the upper part.
- **Incomplete beta:** division-free (Wallis) continued fraction.
- **Elliptic:** AGM for K and E, and Bulirsch's `cel` for complete Π. Carlson R_F, R_D and R_J handle the rest, stopping duplication at 2^-40 and finishing with a degree-12 Taylor tail.

Domain and limits:
- **Bessel:** ν ≥ 0 and x ≥ 0 for `cyl_*`, `sph_*` and `yn`. Negative values are a domain error. `jn` accepts any n and any x.
  - J and Y have *absolute* error ~2^-470 relative to their envelope, so a result within ~2^-200 of a zero can lose relative accuracy.
  - Cost grows with ν (recurrences take ν steps) and with x for ν > x (CF1 takes ~x steps).
- **`beta`:** x, y > 0. Zero is a pole and negative values are a domain error. Above about 2^200 the ln-gamma difference loses accuracy (2^-474 · x ln x); this was checked only up to 2^60.
- **`ibeta`/`ibetac`:** a, b > 0 and 0 ≤ x ≤ 1. The continued fraction needs ~√max(a, b) iterations, so very large a or b is slow.
- **`riemann_zeta`:** below −30000 it overflows. ζ(−2n) is exactly +0, ζ(1) is a pole and ζ(−∞) is a domain error. `hurwitz_zeta(s, a)` takes all real s ≠ 1 and a > 0. s = 1 is a pole, a ≤ 0 or s = −∞ is a domain error, ζ(0, a) = 1/2 − a exactly, and a = +∞ gives 0 (s > 1) or −∞ (s < 1). For s < 1 the error is absolute next to zeros in a, about 2^-470 of the largest term. It takes about 2–4 ms for s > 1 with a ≲ 50.
- **Incomplete gamma:** a > 0 and x ≥ 0 (Γ(a, x) for a ≤ 0 is not provided). Large a costs ~√a iterations.
- **Elliptic:** |k| ≤ 1. `comp_ellint_1(±1)`, `comp_ellint_3(k, 1)` and F(±1, φ) for |φ| ≥ π/2 are poles.
  - `ellint_3` needs 1 − ν sin²φ > 0, and ν < 1 when |φ| > π/2. Cauchy principal values (ν > 1) are not supported.
  - `ellint_3` uses the C++17 sign convention (1 − ν sin²θ).

**Speed.** These functions cost 18 µs to 0.7 ms per call (`hurwitz_zeta` about 1.9 ms) (see the Performance table), well above the elementary functions. Most of the time goes into the dd recurrences and continued fractions.

## Performance

`base256_bench [seconds] [filter]` (built with the tests) prints ns/call. GCC 13 `-O2`, one x86-64 core:

| Function | 4.12.0 | 4.13.0 | 4.14.0 | 4.15.0 | 4.16.0 | Total speedup |
|---|---:|---:|---:|---:|---:|---:|
| `add` | 72 ns | 62 ns | 42 ns | 39 ns | 36 ns | 2.0× |
| `mul` | 66 ns | 55 ns | 38 ns | 43 ns | 38 ns | 1.7× |
| `div` | 2.0 µs | 140 ns | 144 ns | 159 ns | 142 ns | 14.1× |
| `fma` | 118 ns | 107 ns | 117 ns | 121 ns | 104 ns | 1.1× |
| `sqrt` | 2.1 µs | 397 ns | 410 ns | 418 ns | 386 ns | 5.4× |
| `exp` | 328.5 µs | 16.0 µs | 16.4 µs | 16.9 µs | 15.7 µs | 21.0× |
| `exp2` | 325.3 µs | 16.5 µs | 16.4 µs | 16.9 µs | 16.3 µs | 19.9× |
| `expm1` | 326.1 µs | 15.5 µs | 15.6 µs | 18.1 µs | 15.4 µs | 21.1× |
| `log` | 436.5 µs | 27.4 µs | 27.1 µs | 11.4 µs | 8.8 µs | 49.4× |
| `log2` | 447.6 µs | 29.1 µs | 28.4 µs | 12.2 µs | 10.1 µs | 44.3× |
| `log10` | 446.9 µs | 28.4 µs | 29.1 µs | 12.7 µs | 10.5 µs | 42.5× |
| `log1p` | 434.5 µs | 27.1 µs | 26.9 µs | 10.0 µs | 8.9 µs | 48.9× |
| `pow` | 751.0 µs | 43.8 µs | 43.7 µs | 25.1 µs | 25.2 µs | 29.8× |
| `cbrt` | 13.6 µs | 3.5 µs | 3.0 µs | 3.1 µs | 3.2 µs | 4.2× |
| `sin` | 450.4 µs | 22.7 µs | 22.4 µs | 19.6 µs | 20.0 µs | 22.5× |
| `cos` | 447.0 µs | 24.8 µs | 22.2 µs | 19.4 µs | 19.8 µs | 22.6× |
| `tan` | 463.6 µs | 25.2 µs | 23.8 µs | 21.0 µs | 20.9 µs | 22.2× |
| `sin_1e30` | 466.8 µs | 27.7 µs | 26.4 µs | 24.6 µs | 24.4 µs | 19.1× |
| `sin_1e9999` | 469.0 µs | 28.3 µs | 26.9 µs | 24.5 µs | 24.3 µs | 19.3× |
| `atan` | 600.3 µs | 45.6 µs | 43.1 µs | 16.0 µs | 17.2 µs | 34.9× |
| `asin` | 631.9 µs | 47.6 µs | 45.5 µs | 17.8 µs | 17.1 µs | 36.9× |
| `acos` | 639.2 µs | 47.4 µs | 44.8 µs | 17.8 µs | 17.9 µs | 35.8× |
| `atan2` | 521.9 µs | 38.4 µs | 36.6 µs | 16.6 µs | 16.3 µs | 32.0× |
| `sinh` | 338.0 µs | 18.7 µs | 17.9 µs | 17.6 µs | 17.6 µs | 19.3× |
| `cosh` | 334.3 µs | 18.8 µs | 17.2 µs | 17.7 µs | 17.7 µs | 18.9× |
| `tanh` | 335.0 µs | 18.6 µs | 17.5 µs | 17.6 µs | 17.2 µs | 19.5× |
| `asinh` | 427.5 µs | 30.2 µs | 27.3 µs | 10.8 µs | 10.6 µs | 40.2× |
| `acosh` | 445.8 µs | 30.2 µs | 29.6 µs | 10.8 µs | 11.0 µs | 40.4× |
| `atanh` | 856.8 µs | 55.2 µs | 53.1 µs | 18.2 µs | 17.7 µs | 48.3× |
| `erf` | 1616.6 µs | 143.8 µs | 57.2 µs | 25.0 µs | 24.1 µs | 67.0× |
| `erfc` | 2485.1 µs | 265.6 µs | 113.5 µs | 39.1 µs | 38.7 µs | 64.2× |
| `tgamma` | 1427.9 µs | 128.6 µs | 79.4 µs | 71.6 µs | 54.5 µs | 26.2× |
| `lgamma` | 602.7 µs | 56.3 µs | 48.7 µs | 32.3 µs | 20.7 µs | 29.1× |

Special functions (same machine):

| Function | Arguments | 4.17.0 | 4.18.0 |
|---|---|---:|---:|
| `cyl_bessel_j` | ν ∈ [0.02, 10], x ∈ [0.1, 50] | 351 µs | 255 µs |
| `cyl_neumann` | same | 1.29 ms | 666 µs |
| `cyl_bessel_i` | same | 326 µs | 238 µs |
| `cyl_bessel_k` | same | 1.34 ms | 627 µs |
| `jn` | n ∈ [0, 9], x ∈ [0.1, 50] | 281 µs | 203 µs |
| `sph_bessel` | same | 295 µs | 228 µs |
| `beta` | x, y ∈ [0.1, 50] | 157 µs | 108 µs |
| `ibeta` | a, b ∈ [0.1, 50], x = 0.4 | 1.09 ms | 552 µs |
| `riemann_zeta` | x ∈ [−20, 40] | 835 µs | 471 µs |
| `hurwitz_zeta` | s ∈ [1.1, 51], a ∈ [0.1, 50] | — | 1.74 ms |
| `gamma_p` | a, x ∈ [0.1, 50] | 691 µs | 333 µs |
| `gamma_q` | same | 702 µs | 344 µs |
| `comp_ellint_1` | \|k\| < 0.99 | 19 µs | 18 µs |
| `comp_ellint_2` | same | 19 µs | 18 µs |
| `comp_ellint_3` | ν ∈ (−0.9, 0.9) | 1.81 ms | 39 µs |
| `ellint_1` | \|φ\| < 3 | 276 µs | 204 µs |
| `ellint_2` | same | 612 µs | 421 µs |
| `ellint_3` | ν ∈ (−0.5, 0.5), \|φ\| < 1.5 | 1.75 ms | 590 µs |

Run-to-run noise is about ±10% for the basic operations. `fma` has not changed since 4.13.0: its column (117, 121, 104 ns) only shows that noise. Callgrind counts 1427 → 1406 instructions per call from 4.13 to 4.14, and five repeated runs in 4.16.0 give 107–110 ns. Since 4.15.0, `log`, `atan`, `sin`/`cos` and `erf`/`erfc` use MPFR-generated tables (`detail/tables.hpp`, generated by `tools/gen_tables.c`). The `erf`/`erfc` Taylor coefficients at each center are built lazily on first use (about 100 µs per center, thread-safe).

## Notes on specific functions

`round` is half away from zero, matching C `round`. `to_string()` prints the shortest scientific string that round-trips. `to_string(n)` gives `n` significant digits. Use `format`, `to_chars` or `std::format("{:.40f}", x)` for fixed, general or hex output. `hypot` is scaled, so `hypot(2^200000, 2^200000)` stays finite, and `hypot(∞, NaN)` is `+∞`. `fmin` / `fmax` follow IEEE 754 signed-zero rules (`fmin(+0, −0) = −0`). `pow(−0, odd)` keeps the sign of zero. `exp2` on integers and `log2` on powers of two are exact. `tanh` of a huge finite value is `±1` (it does not overflow into NaN). `log1p` of tiny arguments is not flushed. The IEEE special cases `pow(NaN, 0) = 1`, `pow(1, NaN) = 1`, and `pow(−1, ±∞) = 1` hold. Integer `pow` uses dd squaring, so it is exact whenever the power is representable. For accuracy, see [Elementary function accuracy](#elementary-function-accuracy).

## License

Dual-licensed: you may use this software under either the **MIT License** ([LICENSE](LICENSE)) or **THE BEER-WARE LICENSE** (Revision 42, [LICENSE-BEERWARE](LICENSE-BEERWARE)), at your option.

**4.19.0**: `hurwitz_zeta` for all real s ≠ 1, faster Bessel Y/K for x in 30–120. **4.18.0**: faster special functions, `hurwitz_zeta`. **4.17.0**: special functions (`special_math.hpp`). **4.16.0**: `tgamma`/`lgamma` with 140 MPFR Stirling coefficients and shift 64. **4.15.0**: MPFR table-driven `log`, `atan`, `sin`/`cos`, `erf`/`erfc`. **4.8.0**: rounding modes and exception flags, exact `fmod`/`remainder`. **4.7.0**: correctness fixes found by MPFR differential testing. See [CHANGELOG](CHANGELOG.md).

Created by Craig D. Mansfield, PhD, EI (MyGrokOrg).
