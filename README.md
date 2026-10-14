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
| `+ − × ÷`, `fma`, `sqrt`, `fmod`, `remainder` | Correctly rounded or exact in all 5 modes, flags included. Verified bit-for-bit against MPFR (237-bit, subnormalized) on randomized and edge-case operands |
| Subnormals, signed zeros, ±∞, NaN payloads | Supported (see NaN rules) |
| Comparisons | `==`/`!=` quiet, `<`/`<=`/`>`/`>=` signal invalid on NaN; `std::partial_ordering` |
| Conversions from/to `double`, `float`, `long double`, integers | Exact or correctly rounded in the current mode |
| Decimal/hex parsing, literals | Correctly rounded in the current mode |
| Decimal output (`to_string`, `format`, `to_chars`, `std::format`) | Shortest round-trip, or exact at any precision and correctly rounded in the current mode |
| `totalOrder`, `convertToInteger` | `total_order`, `total_order_mag`, `to_integer<T>(mode, exact)` |
| Elementary functions | All 5 modes, flags included. Under 1 ulp in every mode, correct direction in directed modes; no misrounding observed against MPFR (correct rounding not guaranteed) |

Elementary functions are not guaranteed correctly rounded, so `std::numeric_limits<float256>::is_iec559` stays `false`.

## Elementary function accuracy

Every elementary function evaluates in double-float256 (`detail/dd.hpp`, an unevaluated sum `hi + lo` of about 474 bits). The error-free sums and products come from the correctly rounded float256 `+` and `fma`, and the result is rounded to binary256 once. Trigonometric reduction is accurate for every finite argument. Below 2^100 it uses Cody–Waite with a 1190-bit π/2. From 2^100 up to `max()` it uses Payne–Hanek with a 264,192-bit table of 2/π (`detail/two_over_pi.hpp`, 4128 `uint64_t` words): the 237-bit significand is multiplied exactly by the needed 2/π window, keeping 1216 fraction bits. Constants (ln 2, π/2 in five parts, ln 10, √π, ½ ln 2π, 60 Stirling coefficients, γ, ζ(2), ζ(3)) were generated with MPFR.

**Documented bound: < 1 ulp** (round-to-nearest) for all of `exp exp2 expm1 log log2 log10 log1p pow cbrt sin cos tan sincos atan atan2 asin acos sinh cosh tanh asinh acosh atanh erf erfc tgamma lgamma`, over the full binary256 range, with the exceptions below. The evaluation error is about 2^-460 relative, so in practice results are correctly rounded except when the exact value lies within roughly 2^-220 ulp of a rounding midpoint. `tests/mpfr_ulp.cpp` checks this against a 1100-bit MPFR reference on random, log-uniform arguments plus targeted hard regions: near kπ/2 for k up to 2^262140, sin/cos/tan over the whole exponent range, continued-fraction worst cases within ~2^-238 of a multiple of π/2 for every exponent, near ±1 for `asin`/`acos`/`atanh`/`acosh`, near 1 and 2 for `lgamma`, and negative arguments for the gamma functions. With 4 seeds × 5000 cases per function (≈ 840k evaluations), the worst observed error is 0.500 ulp for every function and none exceeded the bound. CI runs the check on every push.

Exceptions and limits:
- `lgamma` within ~2^-220 of its non-trivial zeros on the negative axis (and `tgamma` there, where it is ±1) has an absolute error of about 2^-460, not a relative bound.
- `float256::pi()` is π rounded, so `sin(pi())` is the residual `π − pi()` (≈ 1.2e-72), not zero. The same goes for `cos(pi()/2)` and `tan(pi())`.
- Since 4.12.0 every case is checked in all five rounding modes against MPFR rounded in that mode (237-bit, subnormalized, binary256 exponent range), with the flags: inexact iff the result is inexact, overflow iff MPFR overflows, underflow iff tiny after rounding and inexact, never invalid/divbyzero for these arguments. Overflow and gradual-underflow ranges of `exp`, `exp2`, `expm1`, `sinh`, `cosh` are included. Over 53 checks × 1000 cases × 5 modes (265k evaluations), 0 results differed from MPFR's correctly rounded value and 0 flag mismatches. `tests/test_math_flags.cpp` covers exact cases, poles, domain errors, sNaN, overflow/underflow per mode and flag hygiene.


## Notes on specific functions

`round` is half away from zero, matching C `round`. `to_string()` prints the shortest scientific string that round-trips. `to_string(n)` gives `n` significant digits. Use `format`, `to_chars` or `std::format("{:.40f}", x)` for fixed, general or hex output. `hypot` is scaled, so `hypot(2^200000, 2^200000)` stays finite, and `hypot(∞, NaN)` is `+∞`. `fmin` / `fmax` follow IEEE 754 signed-zero rules (`fmin(+0, −0) = −0`). `pow(−0, odd)` keeps the sign of zero. `exp2` on integers and `log2` on powers of two are exact. `tanh` of a huge finite value is `±1` (it does not overflow into NaN). `log1p` of tiny arguments is not flushed. The IEEE special cases `pow(NaN, 0) = 1`, `pow(1, NaN) = 1`, and `pow(−1, ±∞) = 1` hold. Integer `pow` uses dd squaring, so it is exact whenever the power is representable. For accuracy, see [Elementary function accuracy](#elementary-function-accuracy).

## License

Dual-licensed: you may use this software under either the **MIT License** ([LICENSE](LICENSE)) or **THE BEER-WARE LICENSE** (Revision 42, [LICENSE-BEERWARE](LICENSE-BEERWARE)), at your option.

**4.8.0**: rounding modes and exception flags, exact `fmod`/`remainder`. **4.7.0**: correctness fixes found by MPFR differential testing. See [CHANGELOG](CHANGELOG.md).

Created by Craig D. Mansfield, PhD, EI (MyGrokOrg).
