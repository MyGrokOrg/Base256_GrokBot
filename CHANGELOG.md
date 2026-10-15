# Changelog

## Unreleased

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
