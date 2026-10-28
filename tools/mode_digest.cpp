// Digest of (result bits, flags) of the global-mode (fenv::set_round) path, in
// all five rounding modes: arithmetic, ldexp, rint/nearbyint, conversions,
// decimal parse/format, every elementary and special function and the complex
// operations. It uses only the pre-4.26 API, so it can be built against an older
// revision's headers: tools/mode_digest.sh REV builds it twice (REV and the
// working tree) and requires identical output, i.e. the default path unchanged.
// usage: mode_digest [arith_cases=200000] [fn_args=24]
#include <base256/base256.hpp>
#include <base256/complex.hpp>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <map>
#include <random>
#include <string>
#include "../tests/fn_lists.hpp"

using namespace base256;

static std::mt19937_64 rng(4260);
static const round_mode kModes[5] = {round_mode::to_nearest_even, round_mode::to_nearest_away, round_mode::toward_zero,
                                     round_mode::upward, round_mode::downward};

// Random encodings biased toward hard cases (as tests/fastpath_diff.cpp).
static float256 pick(const float256& other) {
    std::uint64_t w[4] = {rng(), rng(), rng(), rng()};
    const unsigned r = rng() % 16;
    int be;
    switch (r) {
        case 0: be = 1 + static_cast<int>(rng() % 300); break;
        case 1: be = 524286 - static_cast<int>(rng() % 300); break;
        case 2: be = 0; break;
        case 3: be = 262143 + static_cast<int>(rng() % 40) - 20; break;
        default: be = 1 + static_cast<int>(rng() % 524286); break;
    }
    if (r >= 12) {
        std::uint64_t o[4];
        std::memcpy(o, &other, sizeof o);
        std::memcpy(w, o, sizeof o);
        // Every draw is its own statement: the order of evaluation of operands is
        // unspecified, and GCC and Clang (and MSVC) differ, which would make the
        // arguments, and so the digests, depend on the compiler.
        { const std::uint64_t v = rng(); w[0] ^= v >> (rng() % 64); }
        if (rng() & 1) { const std::uint64_t v = rng(); w[1] ^= v >> (rng() % 64); }
        if (r == 15) w[3] ^= (1ULL << 63);
        be = static_cast<int>((o[3] >> 44) & 0x7FFFF) + static_cast<int>(rng() % 9) - 4;
        if (be < 0) be = 0;
        if (be > 524286) be = 524286;
    }
    if (rng() % 8 == 0) { w[0] = 0; w[1] &= ~0ULL << (rng() % 64); }
    w[3] = (w[3] & ((1ULL << 44) - 1)) | (static_cast<std::uint64_t>(be) << 44) | (rng() & 1 ? (1ULL << 63) : 0);
    float256 x;
    std::memcpy(&x, w, sizeof w);
    return x;
}
// Moderate arguments for the functions: (−10, 10), (0, 1) or ±2^[−30, 30].
static float256 farg() {
    const double hi = static_cast<double>(rng() >> 11), lo = static_cast<double>(rng() >> 11);
    const float256 u = ldexp(float256(hi), -53) + ldexp(float256(lo), -106);  // [0, 1), 106 bits
    switch (rng() % 3) {
        case 0: return u * float256(20) - float256(10);
        case 1: return u;
        default: {
            const int e = static_cast<int>(rng() % 61) - 30;
            return ldexp(u + float256(1), e) * float256(rng() & 1 ? 1 : -1);
        }
    }
}

static std::map<std::string, std::uint64_t> H;
static void mix(const std::string& k, const void* p, std::size_t n, unsigned fl) {
    auto it = H.try_emplace(k, 1469598103934665603ULL).first;
    std::uint64_t h = it->second;
    const auto* b = static_cast<const unsigned char*>(p);
    for (std::size_t i = 0; i < n; ++i) h = (h ^ b[i]) * 1099511628211ULL;
    it->second = (h ^ fl) * 1099511628211ULL;
}
static void mixf(const std::string& k, float256 y, unsigned fl) { mix(k, &y, sizeof y, fl); }
static void mixc(const std::string& k, complex z, unsigned fl) {
    const float256 p[2] = {z.real(), z.imag()};
    mix(k, p, sizeof p, fl);
}
// long double is the x87 80-bit format on x86 Linux/macOS (10 value bytes),
// binary64 on MSVC and Apple arm64, binary128 on some other targets. Only the
// 80-bit digest is comparable with the Linux golden digests; for binary64 the
// result must equal to_double (bits and flags), so it is hashed under its own
// key exactly as to_double is, and tests/digest_golden.cmake requires the two
// digests to agree.
static void mix_ld(long double l, unsigned fl) {
    constexpr int dig = std::numeric_limits<long double>::digits;
    if constexpr (dig == 64) mix("to_long_double", &l, 10, fl);
    else if constexpr (dig == 53) mix("to_long_double_b64", &l, sizeof(double), fl);
    else mix("to_long_double_p" + std::to_string(dig), &l, sizeof l, fl);
}
template <class F> static void each_mode(F&& f) {
    for (auto m : kModes) {
        fenv::set_round(m);
        f();
    }
    fenv::set_round(round_mode::to_nearest_even);
}
#define RUN(key, expr, MIX) do { fenv::clear_flags(); auto r_ = (expr); MIX(key, r_, fenv::get_flags()); } while (0)

int main(int argc, char** argv) {
    const long n = argc > 1 ? std::atol(argv[1]) : 200000;
    const int k = argc > 2 ? std::atoi(argv[2]) : 24;
    for (long i = 0; i < n; ++i) {
        const float256 a = pick(float256::one()), b = pick(a), c = pick(a);
        const int e = static_cast<int>(rng() % 600000) - 300000;
        each_mode([&] {
            RUN("add", a + b, mixf); RUN("sub", a - b, mixf); RUN("mul", a * b, mixf); RUN("div", a / b, mixf);
            RUN("fma", fma(a, b, c), mixf); RUN("sqrt", sqrt(abs(a)), mixf);
            RUN("ldexp", ldexp(a, e % 600), mixf); RUN("ldexp_far", ldexp(a, e), mixf);
            RUN("rint", rint(a), mixf); RUN("nearbyint", nearbyint(a), mixf);
            fenv::clear_flags(); { const double d = a.to_double(); mix("to_double", &d, sizeof d, fenv::get_flags()); }
            fenv::clear_flags(); { const float f = a.to_float(); mix("to_float", &f, sizeof f, fenv::get_flags()); }
            fenv::clear_flags(); { const long double l = a.to_long_double(); mix_ld(l, fenv::get_flags()); }
        });
    }
    for (long i = 0; i < n / 20; ++i) {
        const float256 a = pick(float256::one());
        const int prec = static_cast<int>(rng() % 90) + 1;
        each_mode([&] {
            fenv::clear_flags(); { const std::string s = a.to_string(prec); mix("to_string", s.data(), s.size(), fenv::get_flags()); }
            fenv::clear_flags(); { const std::string s = a.format(std::chars_format::fixed, prec); mix("format_f", s.data(), s.size(), fenv::get_flags()); }
            fenv::clear_flags(); { const std::string s = a.format(std::chars_format::scientific, prec); mix("format_e", s.data(), s.size(), fenv::get_flags()); }
        });
        const std::string dec = a.to_string(80);  // 80 digits: parse must round
        each_mode([&] { RUN("parse", float256::parse(dec).value(), mixf); });
    }
    for (int i = 0; i < k; ++i) {
        const float256 x = farg(), y = farg(), z = farg();
        const unsigned nn = static_cast<unsigned>(rng() % 12), mm = static_cast<unsigned>(rng() % 4);
        const int ii = static_cast<int>(rng() % 21) - 10;
#define X1(f) RUN(#f, f(x), mixf);
#define X2(f) RUN(#f "/2", f(x, y), mixf);
#define X3(f) RUN(#f, f(x, abs(y) / (abs(y) + float256(1)), abs(z) / (abs(z) + float256(1))), mixf);
#define XN(f) RUN(#f "/n", f(nn, x), mixf);
#define XI(f) RUN(#f, f(ii, x), mixf);
#define XNM(f) RUN(#f, f(nn, mm <= nn ? mm : nn, x), mixf);
        each_mode([&] {
            B256_FN_X(X1) B256_FN_XY(X2) B256_FN_XYZ(X3) B256_FN_NX(XN) B256_FN_IX(XI) B256_FN_NMX(XNM)
            RUN("bernoulli_b2n", bernoulli_b2n(nn * 7), mixf);
            fenv::clear_flags(); { float256 s, co; sincos(x, &s, &co); const float256 p[2] = {s, co}; mix("sincos", p, sizeof p, fenv::get_flags()); }
        });
        const complex zc(x, y), wc(z, x * float256(0.25));
#define CZ(f) RUN("cx_" #f, f(zc), mixc);
#define CR(f) RUN("cx_" #f, f(zc), mixf);
        each_mode([&] {
            B256_CX_Z(CZ) B256_CX_R(CR)
            RUN("cx_mul", zc * wc, mixc); RUN("cx_div", zc / wc, mixc); RUN("cx_rdiv", x / wc, mixc);
            RUN("cx_add", zc + wc, mixc); RUN("cx_sub", zc - wc, mixc); RUN("cx_smul", zc * y, mixc);
            RUN("cx_sdiv", zc / y, mixc);
            RUN("cx_pow", pow(zc, wc), mixc); RUN("cx_powr", pow(zc, y), mixc); RUN("cx_rpow", pow(x, wc), mixc);
            RUN("cx_polar", polar(abs(x), y), mixc);
        });
    }
    for (const auto& [key, h] : H) std::printf("%-16s %016llx\n", key.c_str(), static_cast<unsigned long long>(h));
    return 0;
}
