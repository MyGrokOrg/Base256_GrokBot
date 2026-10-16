// Per-call rounding: every f(args..., m) must give the same bits and raise the
// same flags as f(args...) evaluated with fenv::set_round(m), while the thread's
// own mode (set to a different "decoy" mode here) is left as it was. Also checks
// the constexpr overloads in constant expressions.
#include <base256/base256.hpp>
#include <base256/complex.hpp>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <random>
#include <stdexcept>
#include <string>
#include "fn_lists.hpp"

using namespace base256;

static std::mt19937_64 rng(4261);
static const round_mode kModes[5] = {round_mode::to_nearest_even, round_mode::to_nearest_away, round_mode::toward_zero,
                                     round_mode::upward, round_mode::downward};
static long passed = 0, failed = 0;

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
        w[0] ^= rng() >> (rng() % 64);
        if (rng() & 1) w[1] ^= rng() >> (rng() % 64);
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
static float256 farg() {
    const float256 u = ldexp(float256(static_cast<double>(rng() >> 11)), -53) +
                       ldexp(float256(static_cast<double>(rng() >> 11)), -106);
    switch (rng() % 3) {
        case 0: return u * float256(20) - float256(10);
        case 1: return u;
        default: return ldexp(u + float256(1), static_cast<int>(rng() % 61) - 30) * float256(rng() & 1 ? 1 : -1);
    }
}

template <class T> static bool same_bits(const T& a, const T& b) { return std::memcmp(&a, &b, sizeof(T)) == 0; }
static bool same_bits(const long double& a, const long double& b) { return std::memcmp(&a, &b, 10) == 0; }
static bool same_bits(const complex& a, const complex& b) {
    return same_bits(a.real(), b.real()) && same_bits(a.imag(), b.imag());
}
static bool same_bits(const std::string& a, const std::string& b) { return a == b; }
template <class T> static bool same_bits(const std::complex<T>& a, const std::complex<T>& b) {
    const T x[2] = {a.real(), a.imag()}, y[2] = {b.real(), b.imag()};
    return std::memcmp(x, y, sizeof x) == 0;
}
static bool same_bits(const std::optional<float256>& a, const std::optional<float256>& b) {
    return a.has_value() == b.has_value() && (!a || same_bits(*a, *b));
}

// glob(): the operation under fenv::set_round(m); per(m): the per-call overload.
template <class G, class P>
static void check(const char* name, round_mode m, G glob, P per) {
    fenv::set_round(m);
    fenv::clear_flags();
    const auto want = glob();
    const unsigned wf = fenv::get_flags();
    const round_mode decoy = kModes[(static_cast<int>(m) + 2) % 5];
    fenv::set_round(decoy);
    fenv::clear_flags();
    const auto got = per(m);
    const unsigned gf = fenv::get_flags();
    const bool mode_kept = fenv::get_round() == decoy;
    fenv::set_round(round_mode::to_nearest_even);
    if (same_bits(want, got) && wf == gf && mode_kept) { ++passed; return; }
    if (++failed <= 30)
        std::printf("FAIL %s mode=%d: bits %s, flags want %#x got %#x, mode %s\n", name, static_cast<int>(m),
                    same_bits(want, got) ? "equal" : "DIFFER", wf, gf, mode_kept ? "kept" : "CHANGED");
}
#define CHECK(name, glob_expr, per_expr)                                                           \
    for (round_mode m_ : kModes)                                                                   \
        check(name, m_, [&] { return glob_expr; }, [&](round_mode m) { return per_expr; })

// ---- constant evaluation: the per-call overloads honour m at compile time,
// where the global environment is fixed to round-to-nearest-even.
namespace ce {
constexpr float256 one = float256(1), three = float256(3), two = float256(2);
constexpr float256 third_up = div(one, three, round_mode::upward);
constexpr float256 third_dn = div(one, three, round_mode::downward);
constexpr float256 third_rz = div(one, three, round_mode::toward_zero);
static_assert(third_up == nextup(third_dn));
static_assert(third_rz == third_dn);
static_assert(div(one, three, round_mode::to_nearest_even) == one / three);
static_assert(div(-one, three, round_mode::toward_zero) == -third_dn);
static_assert(add(one, float256::denorm_min(), round_mode::upward) == nextup(one));
static_assert(add(one, float256::denorm_min(), round_mode::to_nearest_even) == one);
static_assert(sub(one, float256::denorm_min(), round_mode::downward) == nextdown(one));
static_assert(mul(third_up, three, round_mode::upward) > one);
static_assert(sqrt(two, round_mode::upward) == nextup(sqrt(two, round_mode::downward)));
static_assert(fma(third_up, three, -one, round_mode::downward) > float256(0));
static_assert(sub(one, one, round_mode::downward).signbit());   // exact zero: −0 downward
static_assert(!sub(one, one, round_mode::upward).signbit());
static_assert(mul(float256::max(), two, round_mode::toward_zero) == float256::max());
static_assert(mul(float256::max(), two, round_mode::upward).isinf());
static_assert(ldexp(float256::denorm_min() * float256(3), -1, round_mode::upward) == float256::denorm_min() * two);
static_assert(rint(float256(2.5), round_mode::to_nearest_away) == three);
static_assert(nearbyint(float256(2.5), round_mode::to_nearest_even) == two);
static_assert(third_up.to_double(round_mode::upward) > 1.0 / 3.0);
static_assert(third_dn.to_double(round_mode::downward) < third_dn.to_double(round_mode::upward));
static_assert(third_up.to_float(round_mode::toward_zero) < 1.0f / 3.0f);
constexpr complex z(one, three);
static_assert(div(z, three, round_mode::upward).imag() == one);
static_assert(div(z, three, round_mode::upward).real() == third_up);
static_assert(add(z, complex(float256::denorm_min(), float256()), round_mode::upward).real() == nextup(one));
static_assert(sub(one, z, round_mode::downward).real().signbit());
}  // namespace ce

int main() {
    const long n = std::getenv("PERCALL_CASES") ? std::atol(std::getenv("PERCALL_CASES")) : 20000;
    // Arithmetic, ldexp, rint, conversions: biased random encodings.
    for (long i = 0; i < n; ++i) {
        const float256 a = pick(float256::one()), b = pick(a), c = pick(a);
        const int e = static_cast<int>(rng() % 600000) - 300000;
        CHECK("add", a + b, add(a, b, m));
        CHECK("sub", a - b, sub(a, b, m));
        CHECK("mul", a * b, mul(a, b, m));
        CHECK("div", a / b, div(a, b, m));
        CHECK("fma", fma(a, b, c), fma(a, b, c, m));
        CHECK("sqrt", sqrt(abs(a)), sqrt(abs(a), m));
        CHECK("ldexp", ldexp(a, e % 600), ldexp(a, e % 600, m));
        CHECK("ldexp_far", ldexp(a, e), ldexp(a, e, m));
        CHECK("scalbn", scalbn(a, e % 600), scalbn(a, e % 600, m));
        CHECK("rint", rint(a), rint(a, m));
        CHECK("nearbyint", nearbyint(a), nearbyint(a, m));
        CHECK("to_double", a.to_double(), a.to_double(m));
        CHECK("to_float", a.to_float(), a.to_float(m));
        CHECK("to_long_double", a.to_long_double(), a.to_long_double(m));
        const complex z(a, b), w(c, a);
        CHECK("cx_add", z + w, add(z, w, m));
        CHECK("cx_sub", z - w, sub(z, w, m));
        CHECK("cx_addr", z + c, add(z, c, m));
        CHECK("cx_radd", c + z, add(c, z, m));
        CHECK("cx_subr", z - c, sub(z, c, m));
        CHECK("cx_rsub", c - z, sub(c, z, m));
        CHECK("cx_smul", z * c, mul(z, c, m));
        CHECK("cx_muls", c * z, mul(c, z, m));
        CHECK("cx_sdiv", z / c, div(z, c, m));
        CHECK("cx_to_std_d", static_cast<std::complex<double>>(z), z.to_std<double>(m));
        CHECK("cx_to_std_f", static_cast<std::complex<float>>(z), z.to_std<float>(m));
    }
    // Decimal conversions.
    for (long i = 0; i < n / 20; ++i) {
        const float256 a = pick(float256::one());
        const int prec = static_cast<int>(rng() % 90) + 1;
        const std::string dec = a.to_string(80);
        CHECK("to_string", a.to_string(prec), a.to_string(prec, m));
        CHECK("format_f", a.format(std::chars_format::fixed, prec), a.format(std::chars_format::fixed, prec, m));
        CHECK("format_e", a.format(std::chars_format::scientific, prec), a.format(std::chars_format::scientific, prec, m));
        CHECK("parse", float256::parse(dec), float256::parse(dec, m));
        auto tc = [&](round_mode* mm) {
            char buf[200];
            auto r = mm ? to_chars(buf, buf + sizeof buf, a, std::chars_format::scientific, prec, *mm)
                        : to_chars(buf, buf + sizeof buf, a, std::chars_format::scientific, prec);
            return std::string(buf, r.ptr);
        };
        CHECK("to_chars", tc(nullptr), tc(&m));
        auto fc = [&](round_mode* mm) {
            float256 v;
            if (mm) (void)from_chars(dec.data(), dec.data() + dec.size(), v, std::chars_format::general, *mm);
            else (void)from_chars(dec.data(), dec.data() + dec.size(), v, std::chars_format::general);
            return v;
        };
        CHECK("from_chars", fc(nullptr), fc(&m));
    }
    // Functions: moderate arguments.
    const int k = std::getenv("PERCALL_FN_ARGS") ? std::atoi(std::getenv("PERCALL_FN_ARGS")) : 6;
    for (int i = 0; i < k; ++i) {
        const float256 x = farg(), y = farg(), zz = farg();
        const float256 p = abs(y) / (abs(y) + float256(1)), q = abs(zz) / (abs(zz) + float256(1));
        const unsigned nn = static_cast<unsigned>(rng() % 12), mm = static_cast<unsigned>(rng() % 4);
        const unsigned mk = mm <= nn ? mm : nn;
        const int ii = static_cast<int>(rng() % 21) - 10;
#define X1(f) CHECK(#f, f(x), f(x, m));
#define X2(f) CHECK(#f, f(x, y), f(x, y, m));
#define X3(f) CHECK(#f, f(x, p, q), f(x, p, q, m));
#define XN(f) CHECK(#f, f(nn, x), f(nn, x, m));
#define XI(f) CHECK(#f, f(ii, x), f(ii, x, m));
#define XNM(f) CHECK(#f, f(nn, mk, x), f(nn, mk, x, m));
        B256_FN_X(X1) B256_FN_XY(X2) B256_FN_XYZ(X3) B256_FN_NX(XN) B256_FN_IX(XI) B256_FN_NMX(XNM)
        CHECK("bernoulli_b2n", bernoulli_b2n(nn * 7), bernoulli_b2n(nn * 7, m));
        auto sc = [&](round_mode* mo) {
            float256 s, co;
            if (mo) sincos(x, &s, &co, *mo); else sincos(x, &s, &co);
            return complex(s, co);
        };
        CHECK("sincos", sc(nullptr), sc(&m));
        const complex z(x, y), w(zz, x * float256(0.25));
#define CZ(f) CHECK("cx_" #f, f(z), f(z, m));
        B256_CX_Z(CZ) B256_CX_R(CZ)
        CHECK("cx_mul", z * w, mul(z, w, m));
        CHECK("cx_div", z / w, div(z, w, m));
        CHECK("cx_rdiv", x / w, div(x, w, m));
        CHECK("cx_pow", pow(z, w), pow(z, w, m));
        CHECK("cx_powr", pow(z, y), pow(z, y, m));
        CHECK("cx_rpow", pow(x, w), pow(x, w, m));
        CHECK("cx_polar", polar(abs(x), y), polar(abs(x), y, m));
    }
    // The caller's mode is restored when the computation throws.
    fenv::set_round(round_mode::downward);
    try {
        detail::with_round(round_mode::upward, []() -> int { throw std::runtime_error("x"); });
    } catch (const std::runtime_error&) {
    }
    if (fenv::get_round() == round_mode::downward) ++passed; else { ++failed; std::printf("FAIL with_round restore on throw\n"); }
    fenv::set_round(round_mode::to_nearest_even);
    std::printf("test_percall: %ld passed, %ld failed\n", passed, failed);
    return failed ? 1 : 0;
}
