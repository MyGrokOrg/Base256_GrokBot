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

- Correctly-rounded **`+ − × ÷`**, **FMA**, and **`sqrt`** (round-to-nearest, ties to even)
- IEEE 754 specials: `±0`, `±∞`, quiet/signaling NaN, subnormals
- Unordered NaN comparisons (`std::partial_ordering`) and `+0 == −0`
- Conversions to/from `float`, `double`, integers, decimal, and C99 hex-float
- `std::numeric_limits<float256>` (`is_iec559 == true`)
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
```

## Rounding and compliance

Default rounding is IEEE 754 **roundTiesToEven**. Add, subtract, multiply, divide, fused multiply-add, square root, and conversions from integers/decimal apply that mode. `rint` follows that mode; `round` is half away from zero, matching C `round`. `to_string` prints a correctly scaled scientific form (the exponent is corrected when a value sits just below a power of ten). `hypot` is scaled, so `hypot(2^200000, 2^200000)` stays finite, and `hypot(∞, NaN)` is `+∞`. `fmin` / `fmax` follow IEEE 754 signed-zero rules (`fmin(+0, −0) = −0`). `pow(−0, odd)` keeps the sign of zero. `exp2` on integers and `log2` on powers of two are exact. `tanh` of a huge finite value is `±1` (it does not overflow into NaN). `log1p` uses its series for `|x| < 1/16`, so values smaller than `2^{-237}` are not flushed. Integer `pow` with `|n| ≤ 4096` is a product (exact when the power is representable), and the IEEE special cases `pow(NaN, 0) = 1`, `pow(1, NaN) = 1`, and `pow(−1, ±∞) = 1` hold. `tgamma` / `lgamma` use a 40-term Stirling series. Elementary transcendental functions, including `erf` / `erfc`, target a few ulps on the primary domain; trigonometric range reduction uses a 256-bit 2π (very large arguments lose low bits of the reduced angle). Exact multiples of this library's π reduce cleanly: `sin(π)` is `+0` and `cos(π/2)` is `+0`.

## License

Dual-licensed: you may use this software under either the **MIT License** ([LICENSE](LICENSE)) or **THE BEER-WARE LICENSE** (Revision 42, [LICENSE-BEERWARE](LICENSE-BEERWARE)), at your option.

**4.6.5** — `tgamma`/`lgamma`, `sincos`, `remquo`, exact integer `pow`, and a `log1p` fix for tiny arguments. See [CHANGELOG](CHANGELOG.md).

Created by Craig D. Mansfield, PhD, EI (MyGrokOrg).
