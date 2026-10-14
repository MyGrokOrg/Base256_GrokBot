// Fast-path check: prints a digest of (result bits, flags) for millions of
// random +, -, *, /, fma, sqrt operations in all five rounding modes. Built twice,
// with and without BASE256_NO_FAST_PATHS; the two outputs must be identical.
#include <base256/base256.hpp>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <random>

using base256::float256;
namespace fenv = base256::fenv;

static std::mt19937_64 rng(2026);

// Random encodings biased toward interesting cases: equal/near exponents,
// cancellation, carries, subnormal/overflow boundaries.
static float256 pick(const float256& other) {
    std::uint64_t w[4] = {rng(), rng(), rng(), rng()};
    const unsigned r = rng() % 16;
    int be;
    switch (r) {
        case 0: be = 1 + static_cast<int>(rng() % 300); break;           // near min normal
        case 1: be = 524286 - static_cast<int>(rng() % 300); break;      // near overflow
        case 2: be = 0; break;                                           // subnormal / zero
        case 3: be = 262143 + static_cast<int>(rng() % 40) - 20; break;  // near 1
        default: be = 1 + static_cast<int>(rng() % 524286); break;
    }
    if (r >= 12) {  // relative to the other operand: cancellation / carries / alignment
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
    if (rng() % 8 == 0) { w[0] = 0; w[1] &= ~0ULL << (rng() % 64); }  // short significands
    w[3] = (w[3] & ((1ULL << 44) - 1)) | (static_cast<std::uint64_t>(be) << 44) | (rng() & 1 ? (1ULL << 63) : 0);
    float256 x;
    std::memcpy(&x, w, sizeof w);
    return x;
}

int main(int argc, char** argv) {
    const long n = argc > 1 ? std::atol(argv[1]) : 1000000;
    const base256::round_mode modes[5] = {base256::round_mode::to_nearest_even, base256::round_mode::to_nearest_away,
                                          base256::round_mode::toward_zero, base256::round_mode::upward,
                                          base256::round_mode::downward};
    std::uint64_t h[6] = {1469598103934665603ULL, 1469598103934665603ULL, 1469598103934665603ULL,
                          1469598103934665603ULL, 1469598103934665603ULL, 1469598103934665603ULL};
    auto mix = [](std::uint64_t& hh, const float256& y, unsigned fl) {
        std::uint64_t w[4];
        std::memcpy(w, &y, sizeof w);
        for (std::uint64_t v : w) hh = (hh ^ v) * 1099511628211ULL;
        hh = (hh ^ fl) * 1099511628211ULL;
    };
    for (long i = 0; i < n; ++i) {
        const float256 a = pick(float256::one());
        const float256 b = pick(a), c = pick(a);
        for (auto m : modes) {
            fenv::set_round(m);
            float256 y[6];
            unsigned f[6];
            fenv::clear_flags(); y[0] = a + b; f[0] = fenv::get_flags();
            fenv::clear_flags(); y[1] = a - b; f[1] = fenv::get_flags();
            fenv::clear_flags(); y[2] = a * b; f[2] = fenv::get_flags();
            fenv::clear_flags(); y[3] = a / b; f[3] = fenv::get_flags();
            fenv::clear_flags(); y[4] = base256::fma(a, b, c); f[4] = fenv::get_flags();
            fenv::clear_flags(); y[5] = base256::sqrt(base256::abs(a)); f[5] = fenv::get_flags();
            for (int k = 0; k < 6; ++k) mix(h[k], y[k], f[k]);
        }
    }
    fenv::set_round(base256::round_mode::to_nearest_even);
    fenv::clear_flags();
    std::printf("add %016llx\nsub %016llx\nmul %016llx\ndiv %016llx\nfma %016llx\nsqrt %016llx\n",
                (unsigned long long)h[0], (unsigned long long)h[1], (unsigned long long)h[2],
                (unsigned long long)h[3], (unsigned long long)h[4], (unsigned long long)h[5]);
    return 0;
}
