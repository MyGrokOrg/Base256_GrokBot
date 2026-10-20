// 512/256-bit division and the 512-bit integer square root against the bitwise
// reference loops, with random and edge-case operands. Built twice by CMake:
// with the 128-bit integer type where the compiler has one (uint) and with
// -DBASE256_NO_INT128 (uint_noint128), the 64-bit-word path that MSVC and
// clang-cl use.
#include <base256/detail/uint.hpp>
#include <cstdio>
#include <random>

using namespace base256::detail;

static std::mt19937_64 rng(4270);
static std::uint64_t word() {
    switch (rng() % 8) {
        case 0: return 0;
        case 1: return ~0ULL;
        case 2: return 1ULL << (rng() % 64);
        case 3: return ~0ULL >> (rng() % 64);
        case 4: return ~0ULL << (rng() % 64);
        default: return rng();
    }
}

int main() {
    long failed = 0, passed = 0;
    for (long it = 0; it < 200000; ++it) {
        u512 n{};
        u256 d{};
        const int nw = 1 + static_cast<int>(rng() % 8), dw = 1 + static_cast<int>(rng() % 4);
        for (int i = 0; i < nw; ++i) n.d[i] = word();
        for (int i = 0; i < dw; ++i) d.d[i] = word();
        if (rng() % 4 == 0) { for (int i = 0; i < 4; ++i) n.d[4 + i] = d.d[i]; }  // top half == d (qhat clamp)
        if (rng() % 8 == 0) { for (int i = 0; i < 4; ++i) n.d[4 + i] = d.d[i] - (i == 0); }
        if (d.is_zero()) d.d[0] = 1;
        const Div512 a = divmod_512_256(n, d), b = divmod_512_256_bitwise(n, d);
        bool ok = true;
        for (int i = 0; i < 8; ++i) ok = ok && a.quot.d[i] == b.quot.d[i];
        for (int i = 0; i < 4; ++i) ok = ok && a.rem.d[i] == b.rem.d[i];
        if (ok) ++passed;
        else if (++failed <= 10) std::printf("FAIL divmod case %ld\n", it);
        u512 s = n;
        if (rng() & 1) s.d[7] |= 1ULL << 62;
        const u256 r1 = isqrt_512(s), r2 = isqrt_512_bitwise(s);
        bool ok2 = true;
        for (int i = 0; i < 4; ++i) ok2 = ok2 && r1.d[i] == r2.d[i];
        if (ok2) ++passed;
        else if (++failed <= 10) std::printf("FAIL isqrt case %ld\n", it);
    }
    for (std::uint64_t t : {~0ULL, (~0ULL) - 1, 1ULL << 62, (1ULL << 63) + 12345}) {  // isqrt edges
        u512 s{};
        s.d[7] = t;
        s.d[6] = rng();
        bool ok = true;
        const u256 r1 = isqrt_512(s), r2 = isqrt_512_bitwise(s);
        for (int i = 0; i < 4; ++i) ok = ok && r1.d[i] == r2.d[i];
        ok ? ++passed : ++failed;
    }
#if defined(BASE256_HAS_INT128)
    std::printf("test_uint (128-bit integer path): %ld passed, %ld failed\n", passed, failed);
#else
    std::printf("test_uint (64-bit word path): %ld passed, %ld failed\n", passed, failed);
#endif
    return failed ? 1 : 0;
}
