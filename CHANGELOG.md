# Changelog

## Unreleased

- **Portability: MSVC and Apple Clang are tested.** New CI job `portability` (PR #2): MSVC 19.51 (`windows-latest`, `/W4 /WX /permissive-`) and Apple Clang 21 (`macos-latest`, arm64, `-Werror`) build and run all 13 ctests. README "Supported toolchains" lists what is tested and what is not.
  - **`detail/uint.hpp` without a 128-bit integer type** (MSVC, or `-DBASE256_NO_INT128`): `mulu64`, `addc64`, `subb64`, `div128_64` and the 256/512-bit routines run on 64-bit words (`_umul128`, `__umulh`, `_udiv128` at run time on MSVC). `test_uint` (ctest `uint`, and `uint_noint128` for the word path) checks 64x64, 256x256 multiplication, 512/256 division and the 512-bit square root against independent references.
  - **Fix: `mul_wide` in the 64-bit-word path** added the first carry-out into the same limb instead of the next one, so products were wrong whenever that carry occurred (e.g. `mul(1/3 rounded up, 3, upward)` with `-DBASE256_NO_INT128`). Found by `test_percall`'s `static_assert`s. The new multiplication tests: 94556 of 700004 checks fail with the old code and `-DBASE256_NO_INT128`, 0 of 700004 with the fix (g++-14).
  - **Fix: the Hurwitz zeta kernel's direct-sum term count (`zeta`, `hurwitz_zeta`)** used `long` and `1L << 40`, which is 32 bits on Windows (C4293); it is now `long long`.
  - **Golden default-path digests:** ctest `digest_golden` compares `tools/mode_digest 20000 4` (133 digests) with `tests/digest_golden.txt` (g++-14, x86-64) on every platform; locally it also matches clang++-19 and g++-14 `-DBASE256_NO_INT128`. `tools/mode_digest.cpp` drew several random numbers inside one expression (unspecified evaluation order), so its arguments, and so its digests, differed between GCC and Clang (115 of 133 keys); every draw is now its own statement. The identity check against 98a3143 (`tools/mode_digest.sh`) is unaffected: it builds the same file against both header sets.
  - **MSVC `/W4 /WX` and libc++ warnings fixed in the headers:** local `digits` hiding a member (C4458), local `r` hiding another (C4456), signed/unsigned comparison (C4018), unary minus on an unsigned type in `to_integer` (C4146, same value), `numeric_limits` declared `struct` vs MSVC's `class` (C4099, pragma), `numeric_limits::has_denorm` deprecated in C++23 by libc++ (pragma; the members stay).
- **Fix: `atan2(y, x)` (and so complex `arg` and the imaginary part of `log`) when y/x is in or below the subnormal range.** The tiny-ratio path (x > 0, y/x < 2^-122) divided y by x unscaled inside the internal round-to-nearest environment, so a quotient in the subnormal range was rounded twice (wrong in 3 of 5 modes on affected inputs) and one below 2^-262378 was flushed to 0 with underflow and inexact lost. The quotient is now formed at exponent 0 and rounded once, with flags, in the caller's mode. Found by the nightly `cx_check` `log_wide` range (seed 20261019, failing since Oct 23). Regression cases in `test_math_flags` (atan2 below and inside the subnormal range, both signs, all 5 modes, exact flags) and `test_complex` (`log` and `arg` at the fuzz point). Results (g++-14, `SF_RANGE_SEED=1`):
  - `cx_check.py cx_eval 20 20261019 log_wide modes`: 6 errors before, 0 after. `cx_check.py cx_eval 200 S log_wide modes`: seeds 1, 2, 3 gave 36, 18, 30 errors before; seeds 20261019, 1, 2, 3 give 0 after.
  - atan2 sweep against the exact cx_check reference (gmpy2), seed 7, 200 inputs per class, all 5 modes with flags: y/x below the subnormal range 200/200 inputs wrong before (1,000 mode results), 0/200 after; y/x in the subnormal range 110/200 wrong before (332 mode results), 0/200 after; tiny normal y/x 0/200 before and after.
  - The default path differs from 4.25.0 only for these inputs, where 4.25.0 was wrong (the exception is noted in the `default-path` job in `ci.yml`). `tools/mode_digest.sh 98a3143 <cxx> 200000 24` still gives `identical: 133 digests` with g++-14 and clang++-19, since its arguments never reach this path.
- **Fix (test checker): `mpfr_ulp` misjudged results within one 1100-bit ulp of a binary256 number.** A rounded 1100-bit reference was moved a full 1100-bit ulp toward the true value, which could land exactly on a binary256 number: `gamma_q(a, x)` with 1 − Q ≈ 2^-1100.94 was taken as exactly 1, so the correct RZ/RD result 1 − 2^-237 counted as misrounded and the inexact flag as wrong (nightly seed 20261020). The reference is now moved half an ulp, with one exact extra bit, and reset to 1100 bits before each reference call. The failing point is now the first `gamma_q` case of every run. The library was correct. `mpfr_ulp 2000 20261020 2 all` (the CI reproduction): before, `gamma_q` failed with 5 FAIL lines (misrounded 2, flag errors 5); after, exit 0, `gamma_q` misrounded 0, flag errors 0, and no FAIL in any function.

## 4.26.0

- **Release CI.** Run 38986138171 (push of the release commit 50b0e3f): success, all jobs green (`test` g++ and clang++, `mpfr` g++ and clang++, `sanitizers`, `default-path` g++-14 and clang++-19, `compile-time`).
- **New: per-call rounding.** Every operation that rounds has an overload that takes the rounding mode as its last argument: `add(a, b, m)`, `sub`, `mul`, `div`, `fma(a, b, c, m)`, `sqrt(x, m)`, `ldexp`/`scalbn(x, n, m)`, `rint`/`nearbyint(x, m)`, `x.to_double(m)`, `to_float(m)`, `to_long_double(m)`, `float256::parse(s, m)`, `to_string(precision, m)`, `format(fmt, precision, m)`, `to_chars`/`from_chars(..., m)`, all 86 real-function overloads (`exp(x, m)`, `pow(x, y, m)`, `hypot`, `fdim`, `sincos(x, &s, &c, m)`, every special function; `tests/fn_lists.hpp` lists them), and for `base256::complex` `add`/`sub`/`mul`/`div` (complex and real operands), `abs`, `norm`, `arg`, `polar(ρ, θ, m)`, the 16 functions, the three `pow` and `to_std<T>(m)`. The result and the flags are those of the same operation under `fenv::set_round(m)`; the caller's thread mode is unchanged afterwards. README "Per-call rounding" has the table, the semantics, the reasons for the API and what has no per-call form.
  - `constexpr` (the mode is passed down to the rounding kernels and the thread mode is never read or written; usable in constant expressions, where the global mode is fixed to round-to-nearest-even): the arithmetic, `fma`, `sqrt`, `ldexp`/`scalbn`, `rint`/`nearbyint`, `to_double`/`to_float`, complex `add`/`sub` and multiplication and division by a real.
  - The others run the existing code under `detail::with_round(m, f)`, an RAII guard that sets the thread mode for the call and restores it on return and when the call throws.
  - The kernels (`add_sub_m`, `mul_m`, `div_m`, `fma_m`, `sqrt_m`, packing and narrowing) are templates on the mode source: `round_mode` for the per-call forms, `detail::dyn_mode` for the operators, which reads `current_mode()` at the same rounding point as before.
- **The default (thread-mode) path is bit-identical to 4.25.0.** `tools/mode_digest.sh 98a3143` (new) builds `tools/mode_digest.cpp` against the 4.25.0 headers and the new ones and compares FNV digests of the result bits and flags of the global path in all 5 modes: arithmetic, `fma`, `sqrt` (200,000 cases per operation, biased to hard cases), `ldexp`, `rint`, the conversions, decimal and hex text in both directions, every real function and every complex operation (24 arguments each). Result: `identical: 133 digests, 98a3143 vs working tree`, with g++-14 and with clang++-19. All existing tests pass with both compilers (`-Werror`).
- **The identity is checked in CI.** New CI job `default-path` (g++-14 and clang++-19, ubuntu-24.04) runs `tools/mode_digest.sh 98a3143` and `test_percall`, and fails on any difference. Run 38983733663 (push of 79870bb): success, `identical: 133 digests, 98a3143 vs working tree, arith_cases=200000 fn_args=24` with both compilers, and `test_percall: 2533361 passed, 0 failed` with both.
- **Per-call == global, flags included.** New `tests/test_percall.cpp` (ctest `percall`) evaluates every per-call overload in all 5 modes with the thread mode set to a different mode (a decoy), and compares bits and flags with the global evaluation in that mode. It also checks that the decoy mode is left unchanged, that `with_round` restores the mode when the call throws, and the `constexpr` forms with `static_assert`. Default run: 2,533,361 checks passed, 0 failed (g++-14 and clang++-19).
- **Per-call accuracy against MPFR, MPC and mpmath.** The checkers have a per-call mode in which the evaluator calls `f(args..., m)` with the thread mode set to a decoy, and a changed thread mode counts as a failure: `SF_PERCALL=1` (sf_eval / sf_check.py), `CX_PERCALL=1` (cx_eval / cx_check.py), `MPFR_PERCALL=1` (mpfr_diff), `ULP_PERCALL=1` (mpfr_ulp). New ctests `mpfr_diff_percall` (seed 54321) and `mpfr_ulp_percall` (seed 778, strict, all modes). Results (g++-14 build, `SF_RANGE_SEED=1`):
  - `SF_PERCALL=1 sf_check.py sf_eval 40 1` (round-to-nearest, 40 cases per range, seed 1): 233 ranges, 0 not correctly rounded, max 0.5 ulp (as printed with 3 digits). In `bi_ovf` and `bip_ovf`, 20 and 25 results overflow to +∞ (correct); they are counted as overflowed (see the sf_check fix below).
  - `SF_PERCALL=1 sf_check.py sf_eval 20 2 "" modes` (all 5 modes, exact flags, 20 cases per range, seed 2): 233 ranges, `SUMMARY seed=2 cases=20 misrounded+flag errors=0 reference errors=0`. Max ulp over the ranges: RNE and RNA 0.5, RZ, RU and RD 1, except the RZ/RD results of `bi_ovf` and `bip_ovf` that are the largest finite value for an exact value above it (correct; sf_check counts them in max ulp as about 2.1e+71).
  - `CX_PERCALL=1 cx_check.py cx_eval 40 1`: 69 ranges, `SUMMARY seed=1 cases=40 RNE errors=0`, 0 misrounded, max 0.5 ulp, 15 overflowed components.
  - `CX_PERCALL=1 cx_check.py cx_eval 20 2 "" modes`: 69 ranges, `SUMMARY seed=2 cases=20 modes errors=0` (0 misrounded, 0 flag errors, 0 zero-sign errors), max ulp RNE/RNA 0.5, RZ/RU/RD 1, 40 overflowed components (`mul_wide`).
  - `MPFR_PERCALL=1 mpfr_diff 20000 12345`: 960,909 checks, 0 failures (arithmetic, `fma`, `sqrt` and the conversions to `double` and `long double` through the per-call forms, all 5 modes).
  - `ULP_PERCALL=1 mpfr_ulp 30 777 2` (strict, all 5 modes): 82 functions, 0 misrounded, 0 flag errors, 0 direction errors, max ulp RNE/RNA 0.5, directed 1; the same with the clang++-19 build, seed 779.
- **Default-path speed.** `base256_bench 1`, g++-14 14.2.0 `-O3 -DNDEBUG -std=c++23`, shared 8-core box ("Intel(R) Xeon(R) Processor"), pinned with `taskset -c 5`, ≥ 1 s per item over 64 arguments, 98a3143 headers and the new headers built with the same `bench.cpp`, 3 runs interleaved, median (ns/call): `add` 44.3 → 43.8, `mul` 37.9 → 37.0, `div` 139.9 → 129.9, `fma` 104.8 → 98.7, `sqrt` 385.9 → 356.9, `exp` 14,323 → 14,117, `log` 8,098 → 7,945, `sin` 17,452 → 17,525, `pow` 22,374 → 22,634, complex `z * w` 700.7 → 674.6, `z / w` 2,342 → 2,211, `exp` 92,344 → 87,312, `to_string(1..2)` 4,221 → 4,553. Over all 111 items the medians differ by −11.2 % to +7.9 %, with overlapping ranges. The largest increase (`to_string(1..2)`, +7.9 %) is noise: 7 interleaved runs of the to_string, to_chars and complex `sqrt`/`abs` items give 4,192 → 4,184 ns (−0.2 %), and −1.6 % to +1.2 % for the others.
  - **A regression found and fixed before release.** `bench.cpp` times independent calls. A dependent chain `s = s + x[i] * y[i]` (new `bench/chain.cpp`; 6.4 M iterations, wall clock, median of 7, pinned) took 0.324 s with 4.25.0 and 0.480 s with the first draft under g++-14: GCC kept `add_sub_m<dyn_mode>` out of line where it had inlined the old `add_sub`. `BASE256_KERNEL_INLINE` (`[[gnu::always_inline]]` on the kernels, GCC only) restores it: 0.314 s. Clang 19 inlined as before (0.239 s for 4.25.0 and the draft) and was slower with the attribute (0.249 s), so it is not applied there.
  - With the release headers, callgrind instruction counts of `bench/chain.cpp` (128,000 iterations; 25,600 for `div` and `sqrt`) vs 4.25.0: g++-14 `muladd` −0.53 %, `add` −0.26 %, `mul` −0.57 %, `div` +1.10 %, `fma` −0.59 %, `sqrt` +1.65 %; clang++-19 0.00 % for all except `fma` −0.11 %. Timings (median of 7, pinned, ns per iteration, 4.25.0 → 4.26.0): g++-14 `muladd` 70.07 → 73.22, `add` 42.30 → 43.91, `mul` 248.20 → 254.65, `div` 322.86 → 321.62, `fma` 107.15 → 112.22, `sqrt` 392.16 → 380.32; clang++-19 37.46 → 38.41, 18.09 → 18.09, 190.78 → 191.08, 236.73 → 236.26, 60.39 → 61.07, 425.06 → 438.69. The results are bit-identical between the builds.
- **Compile time (4.26.0).** The per-call overloads add 20,592 bytes to the core headers (812,822 → 833,414, +2.5 %) and 5,639 to the complex group (78,679 → 84,318, +7.2 %). CI run 38982941051 (workflow_dispatch on e5369c1, whose headers equal 4.26.0, `compile_refs=98a3143`; compiler CPU s, mean of 5, 4.25.0 → 4.26.0):

  | TU | g++-14 | clang++-19 |
  |---|---|---|
  | core `include_only` | 1.277 → 1.303 (+2.0 %) | 1.795 → 1.840 (+2.5 %) |
  | core `arith` | 1.550 → 1.560 (+0.6 %) | 1.751 → 1.779 (+1.6 %) |
  | core `math` | 3.294 → 3.310 (+0.5 %) | 2.970 → 2.984 (+0.5 %) |
  | core `special` | 6.296 → 6.448 (+2.4 %) | 5.457 → 5.468 (+0.2 %) |
  | core `igamma` | 3.548 → 3.570 (+0.6 %) | 3.197 → 3.253 (+1.8 %) |
  | complex `complex` | 6.152 → 6.206 (+0.9 %) | 5.705 → 5.766 (+1.1 %) |

  Against the group baselines in the same run, 0 warnings: core vs ad368d3 (4.23.0) +1.1 % to +7.3 % per TU (largest g++-14 `special` 6.008 → 6.448), header bytes 797,147 → 833,414 (+4.5 %); complex vs 98a3143 as above. The push runs 38980718638 (e5369c1) and 38983733663 (79870bb) also report 0 warnings.
- sf_check: the round-to-nearest path now counts a +∞/−∞ result whose exact value rounds to ±∞ as "overflowed" (excluded from max ulp and from "not correctly rounded"), as `cx_check` already did. Before, the overflow ranges `bi_ovf` and `bip_ovf` showed max ulp inf and 20 / 25 "not correctly rounded" (40 cases, seed 1, `SF_RANGE_SEED=1`, the same with and without `SF_PERCALL`); now 0, with 20 / 25 overflowed. The `modes` path is unchanged.
- **Known, not changed (to keep the default path bit-identical):** `ldexp`/`scalbn` take a shortcut when the result exponent is more than 16 above E_max or more than 256 below E_min: they return ±∞ or ±0 without flags, in every mode. The correct results are ±∞ (±max in toward-zero mode and when rounding toward the opposite infinity) with `overflow` and `inexact`, and ±0 (or the smallest subnormal in the direction of the sign) with `underflow` and `inexact`. Closer to the range the result is rounded with the right flags (e.g. toward-zero 2^262144 → max with overflow + inexact, but 2^262200 → +∞ without flags; 1.5·2^-262380 → 0 with underflow + inexact, but 1.5·2^-262500 → 0 without flags). The per-call `ldexp(x, n, m)` matches the global path. To be fixed in a later release, together with the global path.
- **Compile time of code that does not use `complex` is unchanged by 4.25.0.** 4.25.0 changed none of the headers that `base256.hpp` includes: `git diff cb33f18 98a3143 -- include` only adds `complex.hpp` and `detail/cx_engine.hpp`, and the 13 headers of the `base256.hpp` include closure are byte-identical (812,822 bytes). 4.24.0 (9f9d15a) → cb33f18 changed `math.hpp` and `extra_math.hpp` (+626 bytes, the fix of directed rounding for |x| < 2^-237). CI run 38974440910 (workflow_dispatch on 91d040f, whose headers equal 98a3143; ubuntu-24.04 runner, `-std=c++23 -O2 -c`, compiler CPU s, mean of 5, the three revisions measured in the same job):

  | TU | 9f9d15a | cb33f18 | 98a3143 | vs 9f9d15a | vs cb33f18 |
  |---|---|---|---|---|---|
  | g++-14 `include_only` | 1.275 | 1.268 | 1.267 | −0.6 % | −0.1 % |
  | g++-14 `arith` | 1.540 | 1.530 | 1.533 | −0.5 % | +0.2 % |
  | g++-14 `math` | 3.292 | 3.256 | 3.284 | −0.2 % | +0.9 % |
  | g++-14 `special` | 6.312 | 6.330 | 6.309 | −0.0 % | −0.3 % |
  | g++-14 `igamma` | 3.576 | 3.578 | 3.563 | −0.4 % | −0.4 % |
  | clang++-19 `include_only` | 1.807 | 1.801 | 1.820 | +0.7 % | +1.1 % |
  | clang++-19 `arith` | 1.775 | 1.765 | 1.771 | −0.2 % | +0.3 % |
  | clang++-19 `math` | 2.969 | 2.971 | 2.963 | −0.2 % | −0.3 % |
  | clang++-19 `special` | 5.427 | 5.409 | 5.439 | +0.2 % | +0.6 % |
  | clang++-19 `igamma` | 3.206 | 3.207 | 3.196 | −0.3 % | −0.3 % |

  The standard deviations are 0.02–0.12 s, so every difference is within about one sd. The 4.25.0 compile-time warning ("headers 797,147 → 891,501 bytes, +11.8 %") came from counting every file in `include/` (812,196 / 812,822 / 891,501 bytes at the three revisions), including `complex.hpp`, which `base256.hpp` does not include.
- **`tools/compile_time.py` measures two header groups, each against its own baseline** (`tools/compile_time_baseline.json`, `"groups"`): `core`, the include closure of `base256/base256.hpp` (13 files), and `complex`, what `complex.hpp` adds to it (`complex.hpp` + `detail/cx_engine.hpp`, 78,679 bytes). The `complex` TU belongs to the `complex` group, all others to `core`. A group's baseline revision is built only for that group's TUs, and the +10 % warning applies to each group's TUs and header bytes on its own. The CI job and `--ref` (also the new workflow_dispatch input `compile_refs`) work as before.
  - `complex` baseline: 98a3143 (4.25.0).
  - `core` baseline: **kept at ad368d3 (4.23.0), not re-baselined.** The core headers did not change in 4.25.0 apart from the cb33f18 fix, so there is nothing new to accept; moving the baseline would only hide the drift accumulated since 4.23.0 (+0.9 % to +6.0 % per TU, +2.0 % bytes, all below the 10 % threshold) from the next comparisons.
  - CI run 38974915314 (push of fb61069, mean of 5, CPU s), **0 warnings**:

    | group, TU | g++-14 | clang++-19 |
    |---|---|---|
    | core `include_only` | 1.235 → 1.270 (+2.8 %) | 1.730 → 1.781 (+3.0 %) |
    | core `arith` | 1.516 → 1.536 (+1.3 %) | 1.719 → 1.745 (+1.5 %) |
    | core `math` | 3.250 → 3.279 (+0.9 %) | 2.907 → 2.943 (+1.2 %) |
    | core `special` | 5.973 → 6.333 (+6.0 %) | 5.188 → 5.407 (+4.2 %) |
    | core `igamma` | 3.509 → 3.551 (+1.2 %) | 3.168 → 3.196 (+0.9 %) |
    | complex `complex` | 6.163 → 6.175 (+0.2 %) | 5.701 → 5.721 (+0.3 %) |

    Header bytes: core 797,147 → 812,822 (+2.0 %), complex 78,679 → 78,679.

## 4.25.0

- **New: `base256::complex`** (`#include <base256/complex.hpp>`, not included by `base256.hpp`). It is its own class with std::complex's interface: `real`/`imag`, `+ − × ÷` with complex and real operands, `==`, `abs`, `norm`, `arg`, `conj`, `proj`, `polar`, `exp`, `log`, `log10`, `sqrt`, `pow` (complex/real × complex/real), `sin`, `cos`, `tan`, `asin`, `acos`, `atan`, `sinh`, `cosh`, `tanh`, `asinh`, `acosh`, `atanh`, `operator<<`/`>>`, `std::format`, the `_i256` literal, and explicit conversions to and from `std::complex<float/double/long double>`. The standard leaves `std::complex<T>` unspecified for other `T` and does not allow overloads in `std`, hence the own class (README "Complex numbers", which has the interface table and the per-function guarantees).
  - Correctly rounded parts, all 5 modes, exact flags: `z * w`, `z / w`, `x / w`, `norm`, `abs`, `sqrt` (exact BigInt arithmetic, one rounding, exact boundary tests), and `pow(z, n)` with Im z ≠ 0 and integral |n| ≤ 64 (exact Gaussian power).
  - Correctly rounded unless the exact part lies within about 2^-2850 relative of a rounding boundary without being on it: `exp`, Re `log`, `log10`, `sin`, `cos`, `tan`, `sinh`, `cosh`, `tanh`, `asin`, `acos`, `atan`, `asinh`, `acosh`, `atanh`, `polar`, other `pow`. These use a Ziv loop at 300, 360, 720, 1440 and 2880 bits; a part still undecided at 2880 bits is taken to be on the boundary. Exact parts, Annex G special values, and arguments with |x|, |y| < 2^-64 (exact leading Taylor term plus the series rest) are decided exactly.
  - **Not proven correctly rounded** (real-library accuracy, < 1 ulp, no misrounding observed against MPFR): `arg` and Im `log` (`atan2`), parts on the real or imaginary axis that are computed by the real functions (so that e.g. `exp(complex(x, 0)).real() == exp(x)`), the π multiples of special values and cuts, and `pow` with real z (Im z = ±0).
  - Not guaranteed: `pow` with |w| ≥ 2^4072 (the extra precision for w·log z is capped at 4096 bits).
  - Special values follow C17 Annex G (G.5.1 `_Cmultd`/`_Cdivd` for `*` and `/`, G.6 for the functions). The README lists the choices where Annex G leaves a sign or flag open (no "may raise invalid", `tanh(+0 + i∞)` = +0 + iNaN per DR 471, `pow` conventions). Known `pow` zero signs: for real z > 0 and real w the imaginary part is the zero w·Im z; an exactly cancelling part of an integer power is +0 (−0 in downward).
  - Unit tests: new `tests/test_complex.cpp` (ctest `complex`), 1235 checks passing with g++-14 and clang++-19. It covers construction, the G.5.1 arithmetic and zero signs, single rounding of `*` and `/`, `abs`/`norm`/`arg`/`polar`, the Annex G tables of every function, branch cuts and signed zeros, exact cases, tiny-argument rounding directions, and stream and `std::format` I/O, in all 5 modes with flags.
- **New checker `tools/cx_check.py`** (69 ranges, `--list`) with the evaluator `cx_eval`. The reference is GNU MPC 1.3.1 (MPFR 4.2.1) through gmpy2 2.2.1 at 2000 bits, redone at 40000 bits when the value is near the rounding grid. Exact rationals are used for `mul`, `div` and `norm`. For tiny arguments (both parts nonzero with exponent < −64) it uses the exact leading term plus a 3000-bit series. `modes` checks all 5 modes with exact flags and zero signs. Components whose exact value exceeds the largest finite are excluded from max ulp and counted as "overflowed" (summed over the modes in `modes`).
- Accuracy, `cx_eval` built from 65dd666's headers (g++-14 `-O2`), `SF_RANGE_SEED=1`. Round-to-nearest, 40 cases per range, seed 1: `SUMMARY … errors=0`, 0 misrounded in every range. `modes` (all 5 modes, exact flags), 20 cases per range, seed 2: `SUMMARY … errors=0`, with 0 misrounded, 0 flag errors and 0 zero-sign errors in every range. The max ulp is the largest error over both parts. "U(a, b)" means uniform with 237-bit resolution, and "±2^[a, b]" means a random 237-bit significand with an exponent drawn uniformly from [a, b]. The directed-mode maxima are below 1 ulp (a "1" is < 1 rounded to 3 digits).

  | range | inputs | RN 40, seed 1: max ulp | `modes` 20, seed 2: max ulp RNE / RNA / RZ / RU / RD |
  |---|---|---|---|
  | `mul_gen` | z, w parts U(−10,10) | 0.498 | 0.494 / 0.494 / 0.973 / 0.946 / 0.996 |
  | `mul_wide` | parts ±2^[−200000,200000] | 0.497 (12 overflowed) | 0.491 / 0.491 / 0.996 / 0.971 / 0.996 (40 overflowed) |
  | `mul_sub` | parts ±2^[−131180,−131060] (subnormal results) | 0.495 | 0.486 / 0.486 / 0.943 / 0.96 / 0.995 |
  | `mul_cancel` | ac − bd nearly cancels | 0.477 | 0.484 / 0.484 / 0.975 / 0.919 / 0.975 |
  | `div_gen` | parts U(−10,10) | 0.5 | 0.476 / 0.476 / 0.94 / 0.999 / 0.975 |
  | `div_wide` | parts ±2^[−100000,100000] | 0.497 | 0.495 / 0.495 / 0.991 / 0.999 / 0.991 |
  | `div_cancel` | bc − ad nearly cancels | 0.5 | 0.479 / 0.479 / 0.946 / 0.906 / 0.946 |
  | `div_exact` | exact quotients | 0 | 0 / 0 / 0 / 0 / 0 |
  | `norm_gen` | parts U(−10,10) | 0.496 | 0.471 / 0.471 / 0.98 / 0.984 / 0.98 |
  | `norm_wide` | parts ±2^[−140000,140000] | 0.483 (3 overflowed) | 0.496 / 0.496 / 0.911 / 0.987 / 0.911 |
  | `abs_gen` | parts U(−10,10) | 0.488 | 0.46 / 0.46 / 0.956 / 0.965 / 0.956 |
  | `abs_wide` | parts ±2^[−5000,5000] | 8.23e-142 | 1.88e-79 / 1.88e-79 / 1.88e-79 / 1 / 1.88e-79 |
  | `abs_exact` | exact |z| (Pythagorean) | 0 | 0 / 0 / 0 / 0 / 0 |
  | `arg_gen` | parts U(−10,10) | 0.498 | 0.422 / 0.422 / 0.929 / 0.929 / 0.988 |
  | `sqrt_gen` | parts U(−10,10) | 0.499 | 0.498 / 0.498 / 0.996 / 0.97 / 0.996 |
  | `sqrt_wide` | parts ±2^[−5000,5000] | 0.497 | 0.441 / 0.441 / 0.938 / 0.988 / 0.968 |
  | `sqrt_exact` | exact roots | 0 | 0 / 0 / 0 / 0 / 0 |
  | `sqrt_cut` | x ∈ (−10,−0.1), |y| ∈ 2^[−3000,−240] | 0.489 | 0.486 / 0.486 / 0.962 / 0.997 / 0.962 |
  | `exp_gen` | parts U(−10,10) | 0.496 | 0.498 / 0.498 / 0.986 / 0.986 / 0.961 |
  | `exp_small` | parts ±2^[−3000,−1] | 0.34 | 0.377 / 0.377 / 1 / 1 / 1 |
  | `exp_big` | x ∈ (−180000,180000), y U(−10,10) | 0.49 | 0.494 / 0.494 / 0.974 / 0.974 / 0.925 |
  | `exp_ybig` | x ∈ (−1,1), y ±2^[10,1000] | 0.498 | 0.496 / 0.496 / 1 / 0.938 / 1 |
  | `log_gen` | parts U(−10,10) | 0.484 | 0.494 / 0.494 / 0.996 / 0.974 / 0.996 |
  | `log_unit` | |z| = 1 + O(2^-30…2^-400) | 0.493 | 0.498 / 0.498 / 0.992 / 0.989 / 0.992 |
  | `log_wide` | parts ±2^[−200000,200000] | 0.495 | 0.495 / 0.495 / 0.986 / 0.986 / 0.798 |
  | `log_cut` | x ∈ (−10,−0.1), |y| ∈ 2^[−3000,−240] | 0.483 | 0.486 / 0.486 / 0.847 / 0.989 / 0.847 |
  | `log10_gen` | parts U(−10,10) | 0.471 | 0.498 / 0.498 / 0.985 / 0.919 / 0.993 |
  | `log10_unit` | |z| near 1 | 0.496 | 0.455 / 0.455 / 0.986 / 0.986 / 0.934 |
  | `pow_gen` | z U(−10,10), w U(−4,4) | 0.476 | 0.481 / 0.481 / 0.997 / 0.997 / 0.969 |
  | `pow_int` | z U(−3,3), integer w ∈ [−64,64] | 0.497 | 0.494 / 0.494 / 1 / 1 / 0.989 |
  | `pow_intbig` | |z| near 1, integer |w| ∈ [65,3000] | 0.487 | 0.5 / 0.5 / 0.983 / 0.998 / 0.983 |
  | `pow_real` | pow(z, real w), U(−10,10) | 0.498 | 0.464 / 0.464 / 0.985 / 0.985 / 0.988 |
  | `pow_wide` | z parts ±2^[−1000,1000], w U(−3,3) | 0.5 | 0.499 / 0.499 / 0.989 / 0.947 / 0.998 |
  | `sin_gen` | parts U(−10,10) | 0.496 | 0.478 / 0.478 / 0.983 / 0.983 / 0.961 |
  | `sin_small` | parts ±2^[−3000,−1] | 0.472 | 0.207 / 0.207 / 1 / 1 / 1 |
  | `sin_bigre` | x ±2^[0,1000], y U(−1,1) | 0.496 | 0.495 / 0.495 / 0.976 / 0.976 / 0.976 |
  | `cos_gen` | parts U(−10,10) | 0.478 | 0.489 / 0.489 / 0.98 / 0.98 / 0.936 |
  | `cos_small` | parts ±2^[−3000,−1] | 0.497 | 0.479 / 0.479 / 1 / 1 / 1 |
  | `cos_bigre` | x ±2^[0,1000], y U(−1,1) | 0.467 | 0.488 / 0.488 / 0.978 / 0.992 / 0.971 |
  | `tan_gen` | parts U(−10,10) | 0.484 | 0.453 / 0.453 / 0.96 / 0.967 / 0.96 |
  | `tan_small` | parts ±2^[−3000,−1] | 0.17 | 0.472 / 0.472 / 1 / 1 / 1 |
  | `tan_bigim` | x U(−10,10), y U(−3000,3000) | 0.494 | 0.487 / 0.487 / 1 / 1 / 1 |
  | `sinh_gen` | parts U(−10,10) | 0.496 | 0.491 / 0.491 / 0.966 / 0.97 / 0.983 |
  | `sinh_big` | x ∈ (−180000,180000), y U(−10,10) | 0.496 | 0.48 / 0.48 / 0.985 / 0.988 / 0.997 |
  | `cosh_gen` | parts U(−10,10) | 0.498 | 0.496 / 0.496 / 0.956 / 0.972 / 0.995 |
  | `cosh_small` | parts ±2^[−3000,−1] | 0.486 | 0.489 / 0.489 / 1 / 1 / 1 |
  | `cosh_big` | x ∈ (−180000,180000), y U(−10,10) | 0.486 | 0.486 / 0.486 / 0.989 / 0.998 / 0.989 |
  | `tanh_gen` | parts U(−10,10) | 0.499 | 0.499 / 0.499 / 0.976 / 0.931 / 0.976 |
  | `tanh_big` | x ±2^[0,12], y U(−10,10) | 0.499 | 0.483 / 0.483 / 1 / 1 / 1 |
  | `asin_gen` | parts U(−10,10) | 0.495 | 0.496 / 0.496 / 0.937 / 0.996 / 0.976 |
  | `asin_branch` | 1 + small | 0.5 | 0.494 / 0.494 / 0.977 / 0.986 / 0.959 |
  | `asin_cut` | x ∈ (1.01,10), |y| ∈ 2^[−3000,−240] | 0.484 | 0.475 / 0.475 / 0.948 / 0.901 / 0.948 |
  | `asin_wide` | parts ±2^[−100000,100000] | 0.495 | 0.493 / 0.493 / 0.938 / 1 / 1 |
  | `acos_gen` | parts U(−10,10) | 0.497 | 0.494 / 0.494 / 0.887 / 0.936 / 0.881 |
  | `acos_branch` | −1 + small | 0.496 | 0.462 / 0.462 / 0.979 / 0.979 / 0.972 |
  | `acos_cut` | x ∈ (−10,−1.01), |y| ∈ 2^[−3000,−240] | 0.492 | 0.497 / 0.497 / 0.985 / 0.872 / 0.985 |
  | `atan_gen` | parts U(−10,10) | 0.495 | 0.475 / 0.475 / 0.982 / 0.985 / 0.982 |
  | `atan_branch` | i + small | 0.496 | 0.483 / 0.483 / 0.989 / 0.914 / 0.989 |
  | `atan_wide` | parts ±2^[−100000,100000] | 0.474 | 0.485 / 0.485 / 1 / 1 / 1 |
  | `asinh_gen` | parts U(−10,10) | 0.497 | 0.485 / 0.485 / 0.981 / 0.981 / 0.986 |
  | `asinh_branch` | −i + small | 0.493 | 0.484 / 0.484 / 0.991 / 0.969 / 0.991 |
  | `acosh_gen` | parts U(−10,10) | 0.487 | 0.483 / 0.483 / 0.996 / 0.993 / 0.996 |
  | `acosh_branch` | 1 + small | 0.493 | 0.496 / 0.496 / 0.973 / 0.958 / 0.977 |
  | `acosh_cut` | x ∈ (−10,0.99), |y| ∈ 2^[−3000,−240] | 0.499 | 0.497 / 0.497 / 0.961 / 0.992 / 0.961 |
  | `atanh_gen` | parts U(−10,10) | 0.497 | 0.484 / 0.484 / 0.97 / 0.994 / 0.97 |
  | `atanh_branch` | 1 + small | 0.5 | 0.5 / 0.5 / 0.992 / 0.99 / 0.992 |
  | `atanh_cut` | x ∈ (1.01,10), |y| ∈ 2^[−3000,−240] | 0.499 | 0.49 / 0.49 / 0.972 / 0.984 / 0.972 |
  | `polar_gen` | ρ 2^[−1000,1000], θ U(−10,10) | 0.489 | 0.495 / 0.495 / 0.936 / 0.936 / 0.918 |
  | `polar_bigth` | ρ ∈ (0.5,2), θ ±2^[0,2000] | 0.478 | 0.487 / 0.487 / 1 / 0.973 / 1 |

- Cost of `base256::complex`, `base256_bench 1` (new `cx …` items in `bench/bench.cpp`), built at b66cfed, g++-14 14.2.0 `-O3 -DNDEBUG -std=c++23` (CMake Release), shared 8-core box ("Intel(R) Xeon(R) Processor"), pinned with `taskset -c 5`, ≥ 1 s per item over 64 arguments, median of 3 runs. ns/call, z = x + iy and w with parts in [−10, 10]:

  | item | ns/call | item | ns/call |
  |---|---|---|---|
  | `cx mul` | 683 | `cx sin` | 120,831 |
  | `cx div` | 2,161 | `cx tan` | 122,652 |
  | `cx abs` | 9,300 | `cx asin` | 363,738 |
  | `cx sqrt` | 19,162 | `cx atanh` | 304,273 |
  | `cx exp` | 85,013 | `cx pow` (z^w) | 428,686 |
  | `cx log` | 69,001 | `cx pow(z,5)` | 2,741 |

  The real functions in the same runs (each on its own argument range, README "Performance") were `exp` 14,081, `log` 7,929, `sin` 17,449, `tan` 18,798, `asin` 15,427, `atanh` 15,223, `pow` 22,307 and `sqrt` 361 ns. So the complex transcendental functions take 6–24× the time of the real function of the same name: `exp` 6.0×, `log` 8.7×, `sin` 6.9×, `tan` 6.5×, `asin` 23.6×, `atanh` 20.0×, `pow` 19.2×. The exact operations grow with the exponent gap between the parts. With parts x·2^200000 and y·2^-200000 (x, y ∈ [1, 2), new "(gap)" items), `cx mul` took 49.1 µs, `cx div` 417 µs, `cx abs` 99.4 µs and `cx sqrt` 37.2 µs.
- Fix (cb33f18): the real functions `sin`, `atan`, `asinh`, `tanh`, `atan2`, `exp`, `exp2`, `expm1` and `erfc` misrounded in the directed modes for |x| < 2^-237. When the first Taylor correction was more than ~600 bits below the leading term, the double-float256 sum dropped it, and the tail passed to the rounding had the wrong sign. Round-to-nearest was not affected. Now only the correction's sign is passed when it is below the dd precision. The complex unit tests found this: `log(1 + i·2^-2000)` takes its imaginary part from `atan2`. A new `test_math_flags` block checks 15 functions at t = ±1.5·2^e for e ∈ {−240, −300, −2000, −100000}, in all 5 modes, value and flags. It had 90 failures before the fix and 0 after (passed=1121, failed=0).
- Nightly fuzz: a new `cx_check` job runs all 69 complex ranges in `modes` (4 shards, workflow input `cx_cases`, default 20, same seed as the other jobs, gmpy2 2.2.1 from pip).
- `tools/compile_time.py`: new TU `complex` (`exp`, `log`, `sin`, `sqrt`, `pow`, `atanh`, `*`, `/` of `base256::complex`). It is reported as n/a for revisions without `complex.hpp`, and as "new TU" in `--compare`.
- Compile time, `tools/compile_time.py --compare` against the 4.23.0 baseline ad368d3, threshold +10 % (warn only). CPU s, ad368d3 → 52f6d24. The CI figures are from the `compile-time` job of run 38970314664 (mean of 5). The dev-box figures are a mean of 10 on the shared box, sd 0.04–0.18 s (another job was compiling at the same time).

  | TU | CI job of run 38970314664 (mean of 5) | dev box (mean of 10) |
  |---|---|---|
  | g++-14 include_only | 1.205 → 1.240 (+2.9 %) | 0.967 → 1.010 (+4.4 %) |
  | g++-14 arith | 1.484 → 1.506 (+1.5 %) | 1.213 → 1.215 (+0.2 %) |
  | g++-14 math | 3.186 → 3.212 (+0.8 %) | 2.623 → 2.663 (+1.5 %) |
  | g++-14 special | 5.910 → 6.224 (+5.3 %) | 4.850 → 5.121 (+5.6 %) |
  | g++-14 igamma | 3.444 → 3.480 (+1.0 %) | 2.810 → 2.822 (+0.4 %) |
  | g++-14 complex (new) | 6.037 | 4.918 |
  | clang++-19 include_only | 1.700 → 1.762 (+3.7 %) | 1.101 → 1.166 (+5.9 %) |
  | clang++-19 arith | 1.691 → 1.717 (+1.5 %) | 1.126 → 1.151 (+2.2 %) |
  | clang++-19 math | 2.870 → 2.908 (+1.3 %) | 2.048 → 2.088 (+1.9 %) |
  | clang++-19 special | 5.080 → 5.295 (+4.2 %) | 3.785 → 3.939 (+4.1 %) |
  | clang++-19 igamma | 3.110 → 3.134 (+0.8 %) | 2.244 → 2.276 (+1.4 %) |
  | clang++-19 complex (new) | 5.602 | 4.268 |

  `complex.hpp` is not included by `base256.hpp`, so the existing TUs do not compile it. Their differences from ad368d3 are about the size measured for 4.24.0 (run 38688772366: include_only +3.0 % / +3.1 %, special +4.7 % / +5.3 %). Only the `complex` TU pays for the new header. On the CI runner it costs 6.04 s (g++-14) and 5.60 s (clang++-19), about 1.9× the `math` TU. Header bytes went from 797,147 to 891,501 (+11.8 %), above the +10 % threshold, so the job printed its one `::warning::`. Of that, 78,679 bytes are the new `complex.hpp` (58,520) and `detail/cx_engine.hpp` (20,159); the rest is 4.24.0 and cb33f18 (812,822 bytes at cb33f18).

## 4.24.0

- Per-call timings of the speedups above (README "Performance", 4.24.0 table): `base256_bench 1.0`, g++-14 14.2.0, `-O3 -DNDEBUG -std=c++23` (CMake Release), on the shared 8-core box ("Intel(R) Xeon(R) Processor"), pinned with `taskset -c 5`, ≥ 1 s per item over 64 arguments, median of 3 runs, 7da2567 and 76f07bb interleaved with the same `bench.cpp`.
  - `polygamma(1/2/3, x ∈ [−50, −0.1])`: 11.4/11.0/10.6 ms → 465/486/501 µs. `polygamma(10, x ∈ [−20, −0.1])`: 12.7 ms → 615 µs. `polygamma(3, x ∈ [0.1, 50])`: 1.20 ms → 270 µs. `trigamma`: 1.24 ms → 255 µs. `hurwitz_zeta(3, a)`: 1.12 ms → 266 µs.
  - Unchanged (non-integer s): `hurwitz_zeta` 1.63 → 1.67 ms (s ∈ [1.1, 51]), 1.74 → 1.70 ms (s ∈ [2, 4], a ∈ [0.4, 0.8]), 1.21 → 1.21 ms (s < 1). `digamma` 145 → 144 µs, `riemann_zeta` 443 → 429 µs.
  - Decimal at the extreme exponents: `to_string` 31.4 ms → 34.3 µs (x·2^262140), 31.8 ms → 30.5 µs (x·2^-262370), 31.6 ms → 31.5 µs (x·2^-262000). `format(scientific, 75)` 10.2 ms → 5.37 µs and 10.1 ms → 4.41 µs. `to_chars` 31.3 ms → 34.8 µs. `format(fixed, 0)` of x·2^262140 (78,913 digits) 119 → 40.2 ms.
  - `to_string` of x ∈ [1, 2) (exact path in both): 4.21 → 4.39 µs (+4 %).
    - Investigated afterwards: it is noise, not a slowdown. Callgrind (2000 calls, g++-14 `-O3`) counts 132,814,701 instructions for 7da2567 and 131,141,095 for 4.24.0 (−1.3 %), with identical malloc/free costs. The only new work on this path is one |k| < 160 early-out per conversion step and an empty vector.
    - Re-measured with the same method: full item list, median of 3 interleaved pinned runs, 4.62 µs (4.51–4.71) vs 5.00 µs (4.52–5.08), overlapping. Ten interleaved runs of the item alone: medians 4.51 µs (4.44–5.10) vs 4.55 µs (4.42–5.19), +1.0 %. Standalone loop: 4.00–4.33 vs 4.06–4.26 µs.
    - No code change. Output re-checked on 9f9d15a: `test_format` passes, and the 800-value side-by-side against 7da2567 is still identical.
- Compile time with these changes, `tools/compile_time.py` against the 4.23.0 baseline ad368d3, threshold +10 % (warn only), 0 warnings in both runs. Headers 797,147 → 812,196 bytes (+1.9 %). CPU s, ad368d3 → 76f07bb:

  | TU | CI job of run 38688772366 (mean of 5) | dev box (mean of 10) |
  |---|---|---|
  | g++-14 include_only | 1.246 → 1.284 (+3.0 %) | 0.976 → 1.014 (+3.8 %) |
  | g++-14 arith | 1.541 → 1.559 (+1.2 %) | 1.181 → 1.228 (+4.0 %) |
  | g++-14 math | 3.278 → 3.320 (+1.3 %) | 2.605 → 2.633 (+1.1 %) |
  | g++-14 special | 6.067 → 6.355 (+4.7 %) | 4.912 → 5.236 (+6.6 %) |
  | g++-14 igamma | 3.528 → 3.567 (+1.1 %) | 2.846 → 2.846 (+0.0 %) |
  | clang++-19 include_only | 1.756 → 1.810 (+3.1 %) | 1.119 → 1.174 (+4.9 %) |
  | clang++-19 arith | 1.736 → 1.766 (+1.7 %) | 1.160 → 1.189 (+2.4 %) |
  | clang++-19 math | 2.957 → 2.993 (+1.2 %) | 2.086 → 2.099 (+0.6 %) |
  | clang++-19 special | 5.249 → 5.527 (+5.3 %) | 3.798 → 3.905 (+2.8 %) |
  | clang++-19 igamma | 3.188 → 3.238 (+1.5 %) | 2.265 → 2.279 (+0.6 %) |

  The box sd is 0.04–0.16 s, so its include-only and arith changes are about one sd. On the CI runner (sd ≤ 0.05 s) include-only is +3 % for both compilers (the cause was not isolated; `float256.hpp`, `bigint.hpp`, `mp.hpp` and both special-function headers changed), and `special` is +4.7 % (g++) and +5.3 % (clang). The `special` figure also includes the Bessel fix of 92396d0, which was +2.4 % / +2.5 % on its own.
- **Speed: decimal output at the extreme exponents, about 1000× faster (shortest and `e` formats), 3× for long `fixed` output.** Before, every conversion multiplied or divided by an exact 5^|k| BigInt (|k| up to ~79000, a 180000-bit number built 13 powers of 5 at a time), and `to_string` did that 3–4 times per call.
  - For |k| ≥ 160, `dec_floor_scaled` (`to_string`, `format`, `to_chars`, `std::format`, `operator<<`) and the interval checks of the shortest search now use a 640-bit lower bound of 5^|k| (binary powering, truncated down, relative error below 2^-600). They answer only when that bound decides the digits, or the comparison, exactly. Otherwise, or with more than ~140 digits, they fall back to the exact BigInt path, which is unchanged (`dec_floor_scaled_exact`). Exact ties cannot occur there (they would need 5^|k| to divide a 256-bit significand or a < 10^140 digit string), so the fallback is only for values within ~2^-600 relative of a boundary.
  - `BigInt::to_dec` splits large integers by 10^(9·2^i) (Knuth division) down to 48 limbs instead of dividing the whole number by 10^9 per 9 digits (`format(fixed, 0)` of 2^262140: 78,913 digits).
  - `format(general)` without precision no longer builds the fixed string just to compare lengths.
  - Output is unchanged. New `test_format` checks: the fast floor answered 296 of 300 random cases (the other 4 fell back) and agreed with the exact floor on all 296. The fast comparison answered 200 of 200 (bounds within a few units of c·10^q), 0 disagree. `format(e, 75)`, `to_string` and `to_chars` round-trip 300 random values with extreme exponents. A side-by-side run against 7da2567 printed identical strings and flags for 800 random binary256 values, uniform over all exponents, in all 5 modes: `to_string`, `format(e, 75)`, `format(e, 3)`, `format(g, 200)` and `format(f, 0…399)` (seeds 1 and 2, 400 values each).
- **Speed: `polygamma` 4.5–24× faster, `trigamma` 4.9×, `hurwitz_zeta` 4.2× at integer s.**
  - x < 0 (`poly_neg`, multiprecision Euler–Maclaurin): s = n + 1 is an integer, so (a + k)^{−s} is now binary powering plus one division (`zmp::pow_neg_int`) instead of exp(−s·log) at P + 64 bits. That was 90 % of the time (callgrind). The B_2k/(2k)! values of `detail_mp::bern_fact` are now cached at the highest precision asked; a cached value truncated to P bits equals the value divided at P bits (both are the exact quotient truncated), so they are bit-identical.
  - x > 0 and `hurwitz_zeta` with an integer s in [2, 512]: the direct terms of the dd Euler–Maclaurin sum are (a/(n + a))^s by binary powering (about (s + 2 log2 s)·2^-476 relative) instead of a dd log and exp per term.
  - Unchanged: `hurwitz_zeta` for non-integer s (about 56 − a direct terms, each a dd log and exp) and `digamma`.
  - Checks, `sf_eval` built from the working tree with these changes, `SF_RANGE_SEED=1`, round-to-nearest 40 cases seed 1 and `modes` (all 5 modes, exact flags) 20 cases seed 2, every range 0 not correctly rounded / 0 misrounded / 0 flag errors (max ulp RNE; `modes` RNE/RZ/RU/RD): `pg_negm` 0.486; 0.488/0.926/0.981/0.926 · `pg_negs` 0.496; 0.478/0.965/0.99/0.965 · `pg_negbig` 0.499; 0.457/0.968/0.993/0.968 · `pg_negtiny` 0.438; 0.49/0.891/0.995/0.891 · `pg_negn` 0.499; 0.447/0.981/0.962/0.981 · `pg_negpole` 0.442; 0.401/1/1/1 · `pg_negz` 0.498; 0.5/0.994/0.983/0.994 · `pg` 0.498; 0.49/0.979/0.918/0.997 · `pg_bign` 0.495; 0.476/0.902/0.996/0.913 · `hz_nz` 0.478; 0.395/0.975/0.965/0.975 · `hz` 0.495; 0.492/0.914/0.925/0.914 · `hz_01` 0.494; 0.494/0.991/0.991/0.763 · `hz_neg` 0.492; 0.48/0.969/0.969/0.761 · `hz_vneg` 0.483; 0.495/0.999/0.958/0.999 · `hz_negba` 0.441; 0.495/0.955/0.955/0.906 · `hz_s1` 0.49; 0.482/0.912/0.991/0.912 · `hz_bigs` 0.492; 0.497/0.94/0.977/0.94 · `hz_biga` 0.439; 0.476/0.903/0.893/0.903 · `lg_negz` 0.486; 0.462/0.996/0.996/0.948 · `lg_neg` 0.498; 0.465/0.983/0.983/0.956 (the `lgamma` multiprecision path shares `bern_fact`).
- Fix: in the reflection path of `hurwitz_zeta` (s ≤ −89, a in (0, 1] or a > 1 where Euler–Maclaurin does not converge at a), σ = 1 − s was formed in float256 and described as exact. It needs 238 bits when |s| + 1 crosses a power of two and the last bit of s is set (s in (−128, −127], (−256, −255], …), so the rounded σ shifted every term. σ is now an exact dd, and sin/cos(πσ/2) are taken as cos/sin(πs/2). Found by the task-5 checks: `modes`, 20 cases, seed 2, `SF_RANGE_SEED=1`, range `hz_vneg`, s = −0x1fcef6d02a34436af90379951a2759aa38e9c36f67a82a869ff5c4d876c3p-230 (−127.2338…), a = 0x1dc0356c32c6efaba38152cf6f3fa98463506e5792824c29b4fbb44f03bcp-235 (3.7189…). On 7da2567 it was misrounded in all 5 modes (648 ulp RNE, `misrounded 5`); the nightly fuzz seeds had not drawn it. That value is now a unit-test reference (correctly rounded −0x1d729c2580b7c9860abaf8049e1b4b082d11f275940e5e80dee1f8d6666ep+134).
- Fuzz results for the workflow above, all from the run logs and artifacts:
  - Run 38673348535 (`workflow_dispatch` on 295e751, seed 20261015 from the UTC date): success. `sf_check` `modes` ran 12 cases × 230 ranges × 5 rounding modes, with all 230 range lines showing 0 misrounded and 0 flag errors (16 shards, `SUMMARY … misrounded+flag errors=0 reference errors=0`). `mpfr_ulp` strict ran 2000 cases × 82 functions with 0 misrounded, 0 flag errors and 0 direction errors. `mpfr_diff` ran 20000 cases, 961,099 checks, 0 failures.
  - Run 38674313074 (on 1196954, seed 4231, 60 cases, `mpfr_cases` 10000): failure. Shards 9 and 14 failed on `ibeta_skw` and `ibeta_big` (reference errors, fixed in sf_check above). Shards 15 and 1 failed on `j_negtiny` and `i_negtiny` (a library bug, fixed below). `mpfr_ulp` showed 82 functions with 0 misrounded, flag or direction errors. `mpfr_diff` ran 4,803,846 checks with 0 failures.
  - Run 38677614231 (on 92396d0, after both fixes, seed 4231): success. `sf_check` `modes` ran 60 cases × 233 ranges × 5 rounding modes, with all 233 range lines showing 0 misrounded and 0 flag errors. `mpfr_ulp` strict ran 10,000 cases × 82 functions with 0 misrounded, flag or direction errors. `mpfr_diff` ran 100000 cases, 4,803,846 checks, 0 failures. The longest job was `mpfr_ulp` at 33 min 50 s (about 34 min).
- Compile-time cost of the Bessel exact-leading-term fix, from the CI `compile-time` job of run 38676713857 (92396d0 vs the ad368d3 baseline, same runner, CPU s, mean of 5): headers 800,551 vs 797,147 bytes (+3,404). `special` TU: g++-14 6.397 vs 6.247 s (+2.4 %), clang++-19 5.569 vs 5.433 s (+2.5 %). `include_only` and `igamma` changed by less than 1 %: g++-14 1.293 vs 1.283 and 3.700 vs 3.674; clang++-19 1.819 vs 1.803 and 3.345 vs 3.359.
- Fix: `cyl_bessel_j` and `cyl_bessel_i` for integer ν in [1, 60], and `sph_bessel` (n ≥ 1), misrounded in the directed modes for tiny |x| (ilogb(x) < −120) when the leading series term x^n/(2^n·n!) or x^n/(2n+1)!! is exactly a binary256 value or a midpoint. The double-double error (~2^-470 relative) could fall on the wrong side, while the true correction (−x²/(4(n+1)) for J and j_n, + for I) is far smaller. The leading term is now computed exactly in integer arithmetic. When it is exact, the sign of the correction decides the rounding. When it is not, the old path is used: the series for ilogb(x) < −300 and the general evaluation otherwise. Fuzz run 2 found this: seed 4231, 60 cases, range `j_negtiny`, ν = 1, x = −0x12034c7600acf127375e26516800b3db8d1a974c0da8e4f2c819bdaa67a1p-544, misrounded in RZ, RU and RD (3 FAIL lines). The same range now gives 0. The same run also hit `i_negtiny` (shard 1): ν = 1, x = −0x1f233670154885dd2484d28a4a786e867215466c1bf042628e6fb27d05f7p-635, misrounded in RZ, RU and RD (3 FAIL lines). With the fix, seed 4231, 60 cases, it gives 0 misrounded and 0 flag errors (3 before). The bug also hit positive x, ilogb(x) in [−300, −120), and j_n (e.g. j_2(15·2^-300)). Unit tests cover exact and midpoint leading terms for J, I and j_n in all modes. New sf_check ranges `j_tinyex`, `i_tinyex`, `sph_tinyex` (ν or n in [1, 8], x = ±m·2^-k, m odd with at most 6 bits, k in [125, 600]). Before the fix (pre-fix `sf_eval`, `modes` = all 5 rounding modes, 40 cases, seed 1): without `SF_RANGE_SEED`, `j_tinyex` had 18 misrounded, `i_tinyex` 12 and `sph_tinyex` 0. With `SF_RANGE_SEED=1`, `j_tinyex` had 24, `i_tinyex` 30 and `sph_tinyex` 6, with 0 flag errors in every run. After the fix, the same `modes` runs (40 cases, seed 1, with and without `SF_RANGE_SEED=1`) give 0 misrounded and 0 flag errors for all three ranges. After the fix (`SF_RANGE_SEED=1`), RNE 40 cases seed 1: `j_tinyex` max 0.47 ulp, `i_tinyex` 0.40, `sph_tinyex` 0.486, 0 not correctly rounded. `modes` 20 cases seed 2: `j_tinyex` (max RNE 0.40, directed 1), `i_tinyex` (0.429 / 1), `sph_tinyex` (0.466 / 1), `j_negtiny`, `i_negtiny`, `sph_jneg`, `j_tinyx`, `sph_j` all have 0 misrounded and 0 flag errors.
- Compile-time measurement: `tools/compile_time.py` compiles representative TUs (`include_only`: `base256.hpp`; `arith`: `float256.hpp` arithmetic; `math`: `math.hpp` exp/log/sin/pow/atan; `special`: Bessel, Airy, polygamma, Hurwitz ζ, Γ(a, x), ζ; `igamma`: `tgamma_upper` only) with g++-14 and clang++-19, `-std=c++23 -O2 -c`. It reports the mean ± sd of N runs of compiler CPU time (user + sys), the wall-clock mean, and the header bytes. Other git revisions can be measured in the same interleaved run. The baseline for 4.23.0 (`ad368d3`) is in `tools/compile_time_baseline.json`. The new CI job `compile-time` measures the baseline revision next to the pushed tree on the same runner (mean of 5). It prints a `::warning::` and never fails when a TU's mean CPU time grows by more than 10 % or the headers by more than 10 %.
- Compile-time baseline, measured on the shared 8-core dev box (another job was running, so sd ≈ 0.05–0.15 s; one 3904fad cell was disturbed). CPU s, mean of 10:

  | | 3904fad (before Airy table) | 575231d (Airy table) | 1638edf (before Γ(a≤0)) | ad368d3 = 4.23.0 |
  |---|---|---|---|---|
  | header bytes | 734,091 | 784,705 | 792,300 | 797,147 |
  | g++-14 include_only | 1.095 | 1.068 | 1.102 | 1.129 |
  | g++-14 arith | (2.104 ± 2.5, disturbed) | 1.446 | 1.292 | 1.292 |
  | g++-14 math | 2.866 | 2.846 | 2.872 | 2.830 |
  | g++-14 special | 5.023 | 5.107 | 5.359 | 5.336 |
  | g++-14 igamma | 2.904 | 2.971 | 2.954 | 3.082 |
  | clang++-19 include_only | 1.210 | 1.222 | 1.265 | 1.244 |
  | clang++-19 arith | 1.239 | 1.246 | 1.259 | 1.266 |
  | clang++-19 math | 2.206 | 2.199 | 2.202 | 2.184 |
  | clang++-19 special | 3.887 | 3.817 | 4.067 | 4.035 |
  | clang++-19 igamma | 2.351 | 2.308 | 2.316 | 2.411 |

  The Airy 2/π table (575231d) added 50,614 header bytes: 49,779 for `detail/two_over_pi_ext.hpp` and 835 in `mp.hpp` and `special_math_ext.hpp` (784,705 bytes in total; the "+49,779" quoted with 575231d was the table file alone). The earlier include-only measurement, g++-14 -O2 wall time, mean of 8, gave 0.956 s → 0.980 s. In this run the include-only and special TUs moved within about one sd. Γ(a, x) for a ≤ 0 (6bed0c1) added 4,847 header bytes. The `igamma` TU went from 2.954 to 3.082 s with g++-14 (+4.3 %) and from 2.316 to 2.411 s with clang++-19 (+4.1 %), about one sd. Include-only changed by +2.4 % (g++) and −1.7 % (clang), which is noise.
- Nightly fuzz workflow `.github/workflows/fuzz.yml` (cron 03:17 UTC and `workflow_dispatch` with inputs `seed`, `cases`, `mpfr_cases`). Without a seed it uses the UTC date (YYYYMMDD), shown in the run summary. It runs every `sf_check` range (230 at the time, 233 now, including the 4.23.0 ones) in `modes` (all 5 modes, exact flags, 2000-bit mpmath reference), split over 16 matrix jobs, plus `mpfr_ulp` in strict mode and `mpfr_diff` against MPFR. Any misrounded result or flag error fails the job and is printed as a `FAIL` line with the function or range, input bits, mode, expected and got value and flags, and the seed. Logs are uploaded as artifacts and the FAIL lines are collected in the run summary.
- sf_check: `--list`; `SF_RANGE_SEED=1` seeds each range from (seed, range name), so its cases do not depend on the filter or shard; `SF_SHARD=i/K`; `SF_TIMING=1`; `SF_MAXFAIL`. `modes` now prints every failing (case, mode) as a `FAIL` line (up to 50 per range) instead of only the first, reports a reference that raises as `REFERR` (skipped, not a failure), reports an evaluator crash as `CRASH` (exit 2), and ends with a `SUMMARY` line. Without the new variables the cases and figures are unchanged.
- mpfr_ulp: `bound_check` = 2 (strict): any misrounded result, flag error or direction error fails, and each is printed as a `FAIL` line with the seed. New function `tgamma_upper_neg` (Γ(a, x), a ∈ (−60, 0], a third of them 0 or negative integers, against `mpfr_gamma_inc` at 3000 bits), also in the ctest list. Local trial, seed 20261014, 400 cases per function, strict, all 82 functions: 0 misrounded, 0 flag errors, 0 direction errors.
- sf_check: the incomplete-gamma references (`gl`, `gu`, `gp`, `gq`, a > 0) now keep the complementary part when it is below the working precision of the whole (Γ(a) or 1). Before, such a value came back as exactly Γ(a) or 1, and the check counted the library's correct result as wrong. The first local trial of the nightly fuzz found this: seed 20261014, range `gq_near`, Q(2652.0, 283.74) = 1 − 2^-5141.4, where the reference said "exactly 1, no flags", while the library correctly returned 1 − ulp in RZ/RD with `inexact`. The library is unchanged. With the fixed reference, `modes` 40 cases per range, seed 1, over `gl`, `gu`, `gp`, `gq`, `gq_mid`, `gp_mid`, `gu_mid`, `gq_tinya`, `gu_tinya`: 0 misrounded, 0 flag errors.
- sf_check: the same fix for the incomplete-beta references (`ib`, `ibc`). mpmath's `betainc` at 2000 bits returned exactly 1 when the complement was below 2^-2000, so the check counted the library's correct 1 − ulp (RZ/RD, `inexact`) as wrong. The complement is now taken as the lower integral I_{1−x}(b, a), and 1 − complement is formed at a precision that keeps it. The second fuzz run found this: seed 4231, 60 cases, `ibeta_skw` (6 misrounded, 9 flag errors) and `ibeta_big` (2 misrounded, 3 flag errors), all reference errors. The library is unchanged. With the fixed reference, `modes`, seed 4231, 60 cases, `SF_RANGE_SEED=1`: `ibeta_skw` and `ibeta_big` have 0 misrounded and 0 flag errors.

## 4.23.0

- `tgamma_upper(a, x)` now takes a ≤ 0 (a = 0, negative integers, negative non-integers) for x > 0, correctly rounded in all 5 modes with exact flags (previously a domain error). For x < 36 and a ≥ −1000 it uses Γ(a) − γ(a, x) with N = round(−a), e = a + N, and the two poles at a = −N cancelled analytically: Γ(a, x) = (−1)^N x^e expm1(L)/(e·N!) − x^a Σ_{n≠N} (−x)^n/(n!(a+n)), L = lnΓ(1+e) − Σ_{k≤N} ln(1 − e/k) − e ln x (e = 0 gives ψ(N+1) − ln x). Otherwise it uses the Legendre continued fraction, which converges fast there for a ≤ 0. Γ(a, 0) = +∞ with divide-by-zero, Γ(a, +∞) = +0 exactly, x < 0 and a = −∞ are invalid; x > 2^1000 underflows, and for a < −2^1000 the result overflows (x < 1) or underflows (x > 1), except Γ(a, 1) = e^{−1}/(1 − a). `tgamma_lower`, `gamma_p` and `gamma_q` still require a > 0. Cost: 0.2–0.8 ms for most arguments, up to ~7 ms for a within 2^-60 of a negative integer near −1000. New sf_check ranges, round-to-nearest 40 cases each, seed 1, and `modes` 20 cases each, seed 2, all with 0 misrounded and 0 flag errors:
  - `gu_neg` (a ∈ (−60, 0), x = 2^[−20, 5]): max 0.491 ulp; `modes` max RNE 0.476, RZ 0.943, RU 0.982, RD 0.943 ulp.
  - `gu_negsx` (a ∈ (−60, 0), x = 2^[−300, −20]): max 0.489 ulp; `modes` max RNE 0.495, RZ 0.892, RU 0.943, RD 0.892 ulp.
  - `gu_negmid` (a ∈ (−1000, 0), x ∈ (0, 36)): max 0.482 ulp; `modes` max RNE 0.455, RZ 0.982, RU 0.87, RD 0.982 ulp.
  - `gu_negbig` (a ∈ (−300, 0), x ∈ (36, 3000)): max 0.493 ulp; `modes` max RNE 0.482, RZ 0.98, RU 0.998, RD 0.98 ulp.
  - `gu_negha` (a ∈ (−30000, −1000), x ∈ (0.75, 1.5)): max 0.485 ulp; `modes` max RNE 0.494, RZ 0.979, RU 0.869, RD 0.979 ulp.
  - `gu_neghx` (a ∈ (−30000, −1000), x ∈ (36, 30000)): max 0.487 ulp; `modes` max RNE 0.46, RZ 0.945, RU 1, RD 0.945 ulp.
  - `gu_zero` (a = 0, x = 2^[−300, 10]): max 0.497 ulp; `modes` max RNE 0.498, RZ 0.922, RU 0.949, RD 0.922 ulp.
  - `gu_negint` (a ∈ {−1, …, −300}, x = 2^[−20, 6]): max 0.48 ulp; `modes` max RNE 0.475, RZ 0.996, RU 0.968, RD 0.996 ulp.
  - `gu_nearint` (a within 2^-[1, 200] of −1 … −40, x = 2^[−20, 5]): max 0.473 ulp; `modes` max RNE 0.493, RZ 0.996, RU 0.982, RD 0.996 ulp.
  - `gu_negtinya` (a = −2^[−300, −1], x = 2^[−20, 6]): max 0.47 ulp; `modes` max RNE 0.491, RZ 0.913, RU 0.952, RD 0.913 ulp.

  Reference: `mp.gammainc` at two working precisions (700+300 and 700+600 bits, or 2000+300/+600 in `modes`) that must agree to 2^-700 (2^-2000); for `gu_neghx`, where `mp.gammainc` does not converge, the Legendre continued fraction in mpmath at +64 and +200 bits, which agreed with `mp.quad` of x^a e^{−x} ∫_0^∞ e^{−xt}(1+t)^{a−1} dt to about 2^-794 or better in 16 spot checks at 800 bits. Unit tests: 16 reference rows (including a next to −3, both sides of a = −1000, and a = −2^1100 at x = 1), special values, overflow/underflow in all modes, and the recurrence Γ(a+1, x) = aΓ(a, x) + x^a e^{−x}.
- Confirmed sf_check figures for Bessel functions at x < 0 (added in `50c9eb8`), round-to-nearest 40 cases each, seed 1, and `modes` 20 cases each, seed 2, all with 0 misrounded and 0 flag errors: `j_negx` (ν ∈ {0…40}, x ∈ (−100, 0)) max 0.476 ulp, `modes` max RNE 0.491, RZ 0.974, RU 0.966, RD 0.974; `j_negbig` (ν ∈ {0…20}, x = −2^[7, 200]) max 0.486, `modes` RNE 0.484, RZ 0.988, RU 0.844, RD 0.988; `j_negtiny` (ν ∈ {0…10}, x = −2^[−400, −8]) max 0.482, `modes` RNE 0.487, RZ 1, RU 0.984, RD 1; `i_negx` (ν ∈ {0…30}, x ∈ (−150, 0)) max 0.491, `modes` RNE 0.476, RZ 0.911, RU 0.911, RD 0.843; `i_negtiny` (ν ∈ {0…10}, x = −2^[−400, −8]) max 0.482, `modes` RNE 0.487, RZ 0.984, RU 1, RD 0.96; `sph_jneg` (n ∈ {0…20}, x ∈ (−50, 0)) max 0.449, `modes` RNE 0.469, RZ 0.953, RU 0.925, RD 0.953; `sph_yneg` (same) max 0.494, `modes` RNE 0.433, RZ 0.905, RU 0.923, RD 0.915 ulp. Non-integer ν for J/I and all x < 0 for Y/K are NaN with `invalid`.
- Confirmed sf_check figures for `polygamma(n, x)` at x < 0 (added in `1638edf`), round-to-nearest 40 cases each, seed 1, and `modes` 20 cases each, seed 2, all with 0 misrounded and 0 flag errors: `pg_negm` (n ∈ {1…10}, x ∈ (−50, 0)) max 0.476 ulp, `modes` max RNE 0.406, RZ 0.999, RU 0.976, RD 0.999; `pg_negs` (x ∈ (−1, 0)) max 0.495, `modes` RNE 0.459, RZ 0.994, RU 0.767, RD 0.994; `pg_negbig` (x = −2^[1, 230]) max 0.499, `modes` RNE 0.482, RZ 0.955, RU 0.955, RD 0.951; `pg_negtiny` (x = −2^[−300, −1]) max 0.49, `modes` RNE 0.46, RZ 0.882, RU 0.924, RD 0.882; `pg_negn` (n ∈ {10…300}, x ∈ (−20, 0)) max 0.498, `modes` RNE 0.414, RZ 0.995, RU 0.995, RD 0.96; `pg_negpole` (x within 2^-[1, 200] of −1 … −40) max 0.491, `modes` RNE 0.3, RZ 1, RU 1, RD 1; `pg_negz` (x within 2^-[1, 236] relative of the zeros of ψ^(n), n ∈ {2, 4, 6}, one per unit interval in (−12, 0)) max 0.497, `modes` RNE 0.488, RZ 0.985, RU 0.932, RD 0.985 ulp. Poles (unit-tested): at x = 0, −1, −2, … odd n gives +∞ with divide-by-zero; even n gives −∞ at +0, +∞ at −0, and NaN with `invalid` at the negative integers.
- sf_check: the Hurwitz reference for s ≤ −89, a ≤ 4 now uses Hurwitz's formula in mpmath at 1200 and 1800 bits (they must agree to 2^-700, else it falls back to `mp.zeta`), because `mp.zeta` at 3000 bits stalls at some points (case 25 of `hz_vneg` seed 1: s ≈ −519.84, a ≈ 0.541). The formula agrees with `mp.zeta` to about 2^-1828 or better in 6 spot checks.
- Re-ran sf_check range `hz_vneg` (s ∈ [−600, −89], a ∈ [0, 4]) against 4.22.0: round-to-nearest 20 cases, seed 1, max 0.492 ulp, 0 misrounded; round-to-nearest 40 cases, seed 1, max 0.492 ulp, 0 misrounded (needs the new reference below); `modes` 20 cases, seed 2, 0 misrounded, 0 flag errors (max RNE 0.5, RNA 0.5, RZ 0.961, RU 0.857, RD 0.961 ulp).
- sf_check: new Airy ranges `ai_deep`, `aip_deep`, `bi_deep`, `bip_deep` (x = −2^[12000, 175000], up to the 2/π-table design limit), with an mpmath reference built from the DLMF 9.7.9–9.7.12 asymptotic series and ζ reduced mod 2π at 1.5·log2|x| + 900 bits (mpmath's `airyai` takes more than 10 minutes per point at −2^20000). The reference agrees with `mp.airyai` to about 2^-897 relative in 16 spot checks at −2^95…−2^3000. Round-to-nearest 40 cases each, seed 1, max 0.5 / 0.481 / 0.493 / 0.467 ulp, 0 misrounded; `modes` 20 cases each, seed 2, 0 misrounded, 0 flag errors (max directed 0.932 / 0.992 / 0.996 / 1 ulp).
- Airy Ai, Ai', Bi, Bi' are now correctly rounded for every negative finite x (previously wrong past about −2^175000, where the 2/π table ran out). New `detail/two_over_pi_ext.hpp` continues the 2/π table to 6176 words = 395264 bits (generated with mpmath, first 4128 words checked against the existing table); only the multiprecision Airy phase reduction uses it. Past the table the functions would return NaN with `invalid`; for finite binary256 x that cannot happen. Cost: headers +49779 bytes (784705 bytes total); compiling a TU that only includes `base256.hpp` (g++-14 -O2, mean of 8) went from 0.956 s to 0.980 s. New sf_check ranges `ai_full`, `aip_full`, `bi_full`, `bip_full` (x = −2^[175000, 262144], up to the most negative finite value): round-to-nearest 40 cases each, seed 1, max 0.497 / 0.498 / 0.47 / 0.491 ulp, 0 misrounded; `modes` 20 cases each, seed 2, 0 misrounded, 0 flag errors (max directed 0.986 / 0.997 / 0.982 / 0.852 ulp). Unit tests: reference rows for all four functions at the most negative finite value, −3·2^200000 and −5·2^175500.

## 4.22.0

- New `detail/mp.hpp`: small BigInt-mantissa multiprecision (exp, log, sqrt, sincos, π, ln 2, Bernoulli numbers, Stirling lgamma) used as a fallback where the dd evaluation cancels.
- `lgamma` is correctly rounded next to its non-trivial negative zeros (was absolute error ~2^-460). sf_check range `lg_negz`: round-to-nearest 100 cases, seed 1, max 0.496 ulp, 0 misrounded; `modes` 40 cases, seed 1, 0 misrounded, 0 flag errors (max RNE 0.485, RZ 0.991, RU 0.991, RD 0.949 ulp).
- `hurwitz_zeta(s, a)`, s < 1: correctly rounded next to zeros in a, via a cancellation estimate and a multiprecision Euler–Maclaurin / Hurwitz-formula fallback; ζ(−2n, ½) = 0 exactly. sf_check range `hz_nz` (15 zeros, s ∈ {−1, −3, −0.5, 0.5, −2.5, −7.25, −30.75, −100.5}): round-to-nearest 30 cases, seed 1, max 0.497 ulp, 0 misrounded; 40 cases, seed 3, max 0.494 ulp, 0 misrounded; `modes` 15 cases, seed 2, 0 misrounded, 0 flag errors (max RNE 0.498, RZ 0.958, RU 0.73, RD 0.958 ulp).
- Airy Ai, Ai', Bi, Bi' for x < −2^93: exact multiprecision reduction of the phase ζ = (2/3)|x|^{3/2} modulo π/2. Correctly rounded down to about −2^175000; below that the 2/π table runs out. sf_check ranges `ai_huge`, `aip_huge`, `bi_huge`, `bip_huge` (x = −2^[94, 3000]): round-to-nearest 20 cases each, seed 1, max 0.453 / 0.484 / 0.463 / 0.5 ulp, 0 misrounded; `modes` 8 cases each, seed 5, 0 misrounded, 0 flag errors (max directed 0.997 / 0.976 / 0.965 / 0.831 ulp). `ai_vhuge` (x = −2^[3000, 12000]): round-to-nearest 4 cases, seed 4, max 0.494 ulp, 0 misrounded.

## 4.21.0

- **Directed modes and flags checked against mpmath.** `tools/sf_check.py EVAL n seed filter modes` (with `SF_MODES=1` in `sf_eval`) evaluates every case in RNE, RNA, RZ, RU and RD. It compares each result with the correctly rounded value of a 2000-bit mpmath reference in that mode and checks the flags (inexact, overflow, underflow with tininess after rounding). Run over every range of polygamma/trigamma, Airy Ai/Bi/Ai'/Bi' (new overflow/underflow ranges near x = 4205), Lambert W0/W−1, erf_inv/erfc_inv, Fresnel, Dawson, dilog/polylog, all Bessel functions, the elliptic integrals and hurwitz_zeta: 0 misroundings and 0 flag errors (40 cases per range).
- **Fixed:** `polylog(2, x)` and `hurwitz_zeta(s, 1)` forwarded to `dilog`/`riemann_zeta` inside their own rounding environment. The result was then computed in round-to-nearest and the inexact flag was dropped. Both now forward before taking the environment, and a unit test covers it.
- Airy overflow/underflow threshold corrected in the docs: x ≈ 4205, not 4960.

## 4.20.0

- **New header `special_math_ext.hpp`** (included by `base256.hpp`), with the rest of C++17 special math and common extras:
  - C++17: `legendre`, `assoc_legendre`, `sph_legendre`, `laguerre`, `assoc_laguerre`, `hermite`, `expint` (Ei).
  - Exponential and trigonometric integrals: `expint(n, x)` (E_n), `sinint`, `cosint`, `sinhint`, `coshint`, `logint`.
  - Gamma family: `digamma`, `trigamma`, `polygamma(n, x)` (x > 0), `lbeta`, `tgamma_ratio`, `tgamma_delta_ratio`.
  - Airy: `airy_ai`, `airy_bi`, `airy_ai_prime`, `airy_bi_prime`.
  - `lambert_w0`, `lambert_wm1`, `erf_inv`, `erfc_inv`, `dawson`, `fresnel_s`, `fresnel_c`.
  - `dilog`, `polylog(n, x)` (integer n), `dirichlet_eta`, `bernoulli_b2n`.
- **Accuracy:** `sf_check` against mpmath (1500–2000 bits) found 0 misroundings, max 0.5 ulp, in every new range. `mpfr_ulp` checks `expint`/`eint`, `logint`/`li`, `digamma`, `lbeta`, `airy_ai`/`ai`, `dilog`/`li2` and the orthogonal polynomials in all five modes: ≤ 1 ulp directed, 0 misroundings, exact flags. `mpfr_ulp` also emulates round-to-nearest-away ties now. Domains and limits are in README "More special functions".
- **Cost** (`base256_bench`, µs/call): polynomials 18–41, Ei/E_n 324/350, Si/Ci/Shi/Chi 300–330, li 168, digamma 157, polygamma 1.2–1.3 ms, lbeta 127, Airy 590–730, Lambert W 91–104, erf_inv 120, Dawson 290, Fresnel 350, dilog 70, polylog(3, ·) 900, eta 497.

## 4.19.0

- **`hurwitz_zeta(s, a)` now takes all real s ≠ 1** (a > 0). For s < 1 it uses analytic continuation: Euler–Maclaurin at N + a, scaled by (N + a)^-s, while that converges. For s ≤ −89 it uses Hurwitz's formula, with a shifted into (0, 1] by its integer part. s = 0 gives 1/2 − a exactly. With a = +∞ the result is 0 for s > 1 and −∞ for s < 1, and s = −∞ is a domain error. For s < 1 the result is correctly rounded except next to zeros of ζ(s, a) in a, where the error is absolute, about 2^-470 of the largest term.
- **Bessel Y/K for x in 30…120** (µs/call at μ = 3.3, 4.18.0 → 4.19.0):

  | x | Y before | Y after | K before | K after |
  |---|---:|---:|---:|---:|
  | 60 | 1377 | 1154 | 2764 | 835 |
  | 100 | 2305 | 1121 | 1908 | 904 |
  | 119 | 2376 | 1122 | 1817 | 1009 |

  Steed's CF1 and CF2 use a division-free Wallis recurrence (CF2 in complex form). For 30 ≤ x < 120, K comes from the Wronskian I_μK_{μ+1} + I_{μ+1}K_μ = 1/x together with the I power series and the K ratio from CF2. Temme now runs to x < 56 (Y) and 30 (K). The target was well under 1 ms, but Y still takes about 1.1 ms and K near 120 about 1.0 ms. I also tried exponentially improved (hyper-)asymptotics; they reach only about 2^-303 at x = 60, so they are not used.
- **`base256_bench`** (range 0.1…50, mostly unaffected): `cyl_neumann` 711 µs, `cyl_bessel_k` 601 µs, `riemann_zeta` 500 µs, `hurwitz_zeta` 1.87 ms. `hurwitz_zeta` timings for single arguments: (2, 0.5) 2.7 ms, (20, 3) 3.8 ms, (1000, 0.5) 29 µs, (1000, 30) 0.52 ms, (0.5, 2) 2.6 ms, (−10, 0.3) 54 µs, (−150, 0.75) 0.20 ms, (−500, 100) 77 µs. For s > 1 with moderate a, most of the time goes into the roughly 56 − a direct terms that Euler–Maclaurin needs before e^{-2π(N+a)} drops below 2^-506.
- **Verification:**
  - `mpfr_ulp` in all modes, plus the special_math, mpfr_diff and fastpath tests, with GCC and clang `-Werror`.
  - `sf_check` against mpmath found 0 misroundings, max 0.5 ulp.
  - New ranges: y_steed, k_wr, hz_01, hz_neg, hz_vneg, hz_negba.
- **Header size and compile time:** `include/` is 646 KB (636 KB before). Including `base256.hpp` costs about the same as 4.18.0, within noise (best of 3, `-std=c++23`): g++ -O0 1.08 → 1.06 s, -O2 1.03 → 1.18 s; clang++ -O0 1.20 → 1.26 s, -O2 1.19 → 1.23 s.

## 4.18.0

- **Faster special functions** (ns/call, `base256_bench`, 4.17.0 → 4.18.0):
  - `cyl_neumann` 1.29 → 0.67 ms, `cyl_bessel_k` 1.34 → 0.63 ms
  - `cyl_bessel_j` 351 → 255 µs, `cyl_bessel_i` 326 → 238 µs, `jn` 281 → 203 µs, `sph_bessel` 295 → 228 µs
  - `ibeta` 1.09 → 0.55 ms, `beta` 157 → 108 µs
  - `gamma_p`/`gamma_q` 691/702 → 333/344 µs (`gamma_p(0.5, 3)`: 11.8 → 0.27 ms)
  - `riemann_zeta` 835 → 471 µs
  - `comp_ellint_3` 1.81 ms → 39 µs, `ellint_3` 1.75 → 0.59 ms, `ellint_2` 612 → 421 µs, `ellint_1` 276 → 204 µs
- **How:**
  - **Bessel Y/K.** The Temme gamma factors (and Γ(z) for 0.5 ≤ z < 40) come from an MPFR Taylor table of 1/Γ(1+μ) instead of gamma calls. The Temme step needs one dd division instead of four. Temme's series now runs to x < 100 (Y) and 60 (K), where it beats Steed's CF2 at this precision; the measured margin is more than 2^-300 relative. Series and Temme terms below 2^-245 of the sum finish in float256.
  - **ibeta and incomplete gamma.** The continued fractions use a division-free Wallis recurrence on the denominator-cleared fraction instead of modified Lentz. The γ series also covers a + 1 ≤ x < 36, where the Γ fraction converged slowly.
  - **Elliptic.** Complete Π uses Bulirsch's quadratically convergent `cel`. Carlson R_F/R_D/R_J duplicate only to 2^-40 and finish with a degree-12 Taylor tail (DLMF 19.19.7) instead of degree 5 at 2^-90. R_C's series truncates by |t|.
  - **Zeta.** A planned Euler–Maclaurin split (N = 24…64 instead of 128), with tabled B_2k/(2k)! and ln p.
- **New `hurwitz_zeta(s, a)`** for s > 1 and a > 0. It picks between a direct sum and Euler–Maclaurin by estimated cost, and scales by a^-s so that huge and tiny a stay in range.
- **Verification:**
  - Full `mpfr_ulp` in all modes (CI), and the special_math, mpfr_diff and fastpath tests pass, with GCC and clang `-Werror`.
  - `tools/sf_check.py` against mpmath found 0 misroundings over all ranges, max 0.5 ulp.
  - New sf_check ranges: y_temme/k_temme, gp/gq/gu_mid, ell_Pi1, ell_Pibig, hz*.
  - The hz* ranges use mpmath at 3000 bits, because its zeta(s, a) loses hundreds of bits at large a.
- **Fixed during development:** a truncation bound in 0f0598d used ilogb(μ) as an upper bound on log2|μ|, which could lose about 60 bits for |μ| near 1/2. It was fixed in e60960f; the checks had not caught it.
- **New tables:** `tools/gen_tables.c` also emits 1/Γ(1+x) Taylor coefficients (112), B_2k/(2k)! (140) and ln p for p < 64. `include/` grows from 580 KB to 636 KB.

## 4.17.0

- **New `special_math.hpp`** (included by `base256.hpp`), using the C++17 special-math names. Every function evaluates in double-float256, rounds once in the current mode and raises IEEE flags.
  - **Bessel:** `cyl_bessel_j`, `cyl_neumann`, `cyl_bessel_i`, `cyl_bessel_k` (real ν ≥ 0, x ≥ 0), `sph_bessel`, `sph_neumann`, and `jn`/`yn` (any integer n; `jn` takes any real x). Methods: power series, Temme series, Steed's CF1/CF2, and Hankel/large-x expansions with stable recurrences.
  - **`beta`, plus `ibeta`/`ibetac`:** the regularized incomplete beta, by continued fraction.
  - **`riemann_zeta`:** all real x. Euler–Maclaurin for x ≥ 0 and the functional equation for x < 0. ζ(1) is a pole, ζ(0) = −1/2 and ζ(−2n) = +0 exactly.
  - **Incomplete gamma:** `tgamma_lower`, `tgamma_upper`, `gamma_p`, `gamma_q`.
  - **Elliptic:** `comp_ellint_1/2/3`, `ellint_1/2/3`. AGM for K and E, Carlson R_F, R_D and R_J for the rest, in the C++17 sign convention for Π.
- **Accuracy.** No result was misrounded in any check.
  - Against MPFR, 1500 cases × 5 modes each, with exact flags: `jn`, `yn`, `cyl_bessel_j`/`cyl_neumann` at integer orders, `beta`, `riemann_zeta` (including around 1 and down to −2^14) and all four incomplete-gamma functions.
  - Against mpmath (`tools/sf_check.py`, 200 cases per range over 63 ranges, round-to-nearest): every function, including non-integer orders, ν ≥ x, ν near integers, huge or tiny x, and the elliptic integrals. All results were within 0.5 ulp.
  - Domain limits and costs are listed in README "Special functions".
- **New tests:**
  - ctest `special_math`: mpmath reference values in all modes, special arguments, poles, domain errors, overflow/underflow and parity.
  - `mpfr_ulp` adds `jn yn cyl_j cyl_y beta zeta zeta1 zetab1 zetan tgamma_lower tgamma_upper gamma_p gamma_q`; these run in CI.
  - New tools: `sf_eval` (built by CMake), `tools/sf_check.py` and `tools/sf_refgen.py`.
- **`tools/gen_tables.c`** now also emits Euler's γ and the z³ Taylor coefficient of 1/Γ(1+z) (for Temme's series).
- **README:** the `fma` benchmark note now says plainly that the 104–121 ns spread is noise. It also lists the Stirling table as 140 coefficients.
- **Header size and compile time:** `include/` grows from 528 KB to 580 KB (`special_math.hpp` is 51 KB). Compiling a TU that includes `base256.hpp` takes the same time as before: g++ -O0 1.25 s, -O2 2.70 s; clang++ -O0 1.22 s, -O2 2.06 s.

## 4.16.0

- **`tgamma` 1.3×, `lgamma` 1.6× faster** (71.6 → 54.5 µs, 32.3 → 20.7 µs).
  - `tools/gen_tables.c` generates 140 Stirling coefficients B_2k/(2k(2k−1)) with MPFR, via ζ(2k). They go into `detail/tables.hpp` and replace the 60-entry hex-string table that was parsed at startup.
  - The Stirling shift drops from 128 to 64, which halves the rising product. At w = 64 the series needs about 101 terms. Fewer than 26 of them are above 2^-270 of the result and use dd Horner steps. The rest are summed in float256, whose rounding error stays below 2^-500 of the result. The term counts come from log2(w).
  - The shift can be set with `BASE256_STIRLING_SHIFT`. Measured: 128 → 61.7 µs, 96 → 59.4, 80 → 58.5, 64 → 56.6, 56 → 65.5 (tgamma, same binary). The internal ln Γ error matches the shift-128 version (worst 2^-464 relative, set by cancellation near z ≈ 2).
- **`fma` has no regression.** Five repeated runs give 107–110 ns, the same as 4.13.0 (107 ns). Callgrind gives 1406 instructions per call vs 1427 in 4.13.0. The 117/121 ns figures were noise.
- **Header size and compile time:** `include/` grows from 515 KB to 528 KB. Compile time is unchanged within noise: g++ -O0 1.28 → 1.29 s, -O2 2.78 → 2.70 s; clang++ -O0 1.22 → 1.17 s, -O2 2.10 → 1.99 s.
- **Accuracy and flags are unchanged.** Full `mpfr_ulp` (3,000 cases × 5 modes, seeds 61 and 72): 0 misrounded, 0 flag errors, 0 direction errors. `mpfr_diff` 300,000 cases (14.4M checks): 0 failures. The fastpath hash on 3M cases is identical.

## 4.15.0

- **MPFR-generated tables** in the new `detail/tables.hpp`, generated by `tools/gen_tables.c` (`cc -O2 tools/gen_tables.c -lmpfr -lgmp && ./a.out > include/base256/detail/tables.hpp`).
- **`log` 2.4×** (27.1 → 11.4 µs), with `log2`/`log10`/`log1p`/`pow`/`asinh`/`acosh`/`atanh` faster too (`pow` 43.7 → 25.1 µs, `asinh` 27.3 → 10.8 µs, `atanh` 53.1 → 18.2 µs). The old approach was a Newton step on exp. Now there is a two-stage table reduction (1 + j/128, then 1 + k/2^15, using short reciprocals so the products are exact), followed by a 2·atanh(t) Horner series with |t| < 2^-16.
- **`erf` 2.3×, `erfc` 2.9×** (57.2 → 25.0 µs, 113.5 → 39.1 µs). For 1/32 ≤ |x| < 16 they use a Taylor expansion of erfc around c_j = j/16. The tables hold erfc(c_j) and 2/√π·e^(−c_j²), both from MPFR. The coefficients come from the Hermite recurrence and are built lazily per center under `std::call_once`. The term count is chosen from |x − c_j|. Linear-range checks `erf_lin`/`erfc_lin` were added to `mpfr_ulp`.
- **`atan`/`asin`/`acos`/`atan2` 2.2–2.7×** (about 44 → 16–18 µs). The argument is reduced by an atan(j/64) table, then an alternating Horner series runs on |t| < 1/128.
- **`sin`/`cos`/`tan` about 13% faster.** A sin/cos(j/64) table and the angle-addition formulas replace 8 double-angle steps for reduced |r| ≤ 1.
- **`fma` regression in 4.14.0 was noise.** Callgrind shows 1406 instructions per call vs 1427 in 4.13.0.
- **Stirling with more coefficients was not adopted.** At shift w the terms needed are about 102 (w = 64), 80 (80), 70 (96) and 59 (128). The total dd operation count stays roughly flat (about 268 / 240 / 246), so lowering the shift gives no real gain.
- **Header size and compile time:** `include/` grows from 319 KB to 515 KB (the tables are 192 KB). Compiling a TU that includes `base256.hpp` takes g++ -O0 1.18 → 1.28 s, -O2 2.59 → 2.58 s; clang++ -O0 1.07 → 1.15 s, -O2 1.99 → 2.03 s.
- **Accuracy and flags are unchanged.** Full `mpfr_ulp` (all checks × 3,000 cases × 5 modes, seeds 41 and 52): 0 misrounded, 0 flag errors, 0 direction errors. `mpfr_diff` 300,000 cases (14.4M checks): 0 failures. The fastpath hash on 3M cases is identical.

## 4.14.0

- **Fast paths for `+`, `−`, `×`**: when both operands are normal and the result is normal before and after rounding, they skip unpacking and special-case handling. They are about 1.45× faster (add 62 → 42 ns, mul 55 → 38 ns). Every other op benefits from `pack()` checking the common normal case first. The results are bit-identical to the general path in all five modes, flags included. `BASE256_NO_FAST_PATHS` turns the fast paths off.
  - New ctest `fastpath`: `tests/fastpath_diff.cpp` is built with and without the fast paths. Each build hashes the result bits and flags of `+ − × ÷ fma sqrt` in all modes over random encodings, biased toward exponent boundaries, cancellation and carries, and the two hashes must match. It was checked locally on 3M cases and across GCC and clang builds on 2M cases. A deliberately off-by-one boundary is detected.
  - `mpfr_diff`: 300,000 cases in all 5 modes (14.4M checks), 0 failures.
- **`erf` 2.5×, `erfc` 2.3×, `tgamma` 1.6×, `lgamma` 1.2× faster.**
  - The erf series uses Horner on a precomputed 2^n/(2n+1)!! table and now covers a < 6. erfc is computed as 1 − erf there, which loses fewer than 52 of about 474 bits.
  - The Laplace continued fraction for erfc (a ≥ 6) uses the forward convergent recurrence: one division in total instead of one per level. Its depth is 6400/a² + 50, measured to keep the truncation error below 2^-300 (the old fixed depths gave 2^-260).
  - `tgamma` computes exp(Stirling) / rising product, or the reflection formula, as m·2^k rather than exp(lgamma). This saves a log and an exp's worth of work.
  - The rising product pairs its factors, halving the dd multiplications.
  - The Stirling correction uses Horner in 1/w², with the term count chosen from the exponents.
- Accuracy and flags are unchanged. The full `mpfr_ulp` suite (53 checks × 3,000 cases × 5 modes, two seeds, about 1.6M evaluations) gives 0 misrounded results, 0 flag errors and 0 direction errors.

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
