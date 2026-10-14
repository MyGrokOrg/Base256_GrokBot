#pragma once

/// Extra elementary functions layered on math.hpp.
/// Gamma, joint sine/cosine, and IEEE remquo live here.
#include "math.hpp"

namespace base256 {

// nearbyint lives in float256.hpp (current rounding mode, no inexact).

[[nodiscard]] inline float256 fdim(float256 x, float256 y) {
    if (x.isnan() || y.isnan()) {
        if (x.issignaling() || y.issignaling()) detail::raise_flag(fe_invalid);
        return x.isnan() ? x.quieted() : y.quieted();
    }
    if (x > y) return x - y;
    return float256::zero(false);
}

// logb lives in float256.hpp (constexpr). A second definition here made
// <base256/base256.hpp> ill-formed.


namespace detail_spec {
// erf(a) for 0 < a < 4 in dd: 2/√π e^{-a²} Σ 2^n a^{2n+1} / (1·3·…·(2n+1)).
// All terms are positive, so the sum has no cancellation.
[[nodiscard]] inline detail_dd::dd erf_series_dd(float256 a) {
    using detail_dd::dd;
    const dd x2 = detail_dd::two_prod(a, a);
    const dd two_x2 = detail_dd::ldexp(x2, 1);
    dd term{a}, sum{a};
    const float256 rel = ldexp(float256::one(), -490);
    for (int n = 1; n <= 2000; ++n) {
        term = term * two_x2 / dd(float256::from_int(2 * n + 1));
        sum = sum + term;
        if (term.hi < abs(sum.hi) * rel) break;
    }
    const auto& c = detail_dd::consts();
    const dd two_over_sqrtpi = dd(float256::two()) / dd(c.sqrtpi[0], c.sqrtpi[1]);
    return sum * detail_dd::exp(-x2) * two_over_sqrtpi;
}

// erfc(a) for a >= 4 as m·2^k: e^{-a²} / (√π · K(a)), with the Laplace
// continued fraction K = a + (1/2)/(a + 1/(a + (3/2)/(a + …))) evaluated
// backward from a fixed depth (chosen so the truncation error is < 2^-260).
[[nodiscard]] inline detail_dd::Scaled erfc_cf_dd(float256 a) {
    using detail_dd::dd;
    const int depth = a < float256(6) ? 400 : a < float256(8) ? 250 : a < float256(10) ? 200
                    : a < float256(20) ? 120 : 80;
    dd f{a};
    for (int n = depth; n >= 1; --n) f = dd(a) + dd(float256::from_int(n) * float256::half()) / f;
    const auto& c = detail_dd::consts();
    detail_dd::Scaled e = detail_dd::exp_scaled(-detail_dd::two_prod(a, a));
    e.m = e.m / (f * dd(c.sqrtpi[0], c.sqrtpi[1]));
    return e;
}
} // namespace detail_spec

[[nodiscard]] inline float256 erf(float256 x) {
    detail_dd::FnEnv env;
    using detail_dd::dd;
    if (x.isnan()) return env.nan(x);
    if (x.isinf()) return env.exact(copysign(float256::one(), x));
    if (x.iszero()) return env.exact(x);
    const float256 a = abs(x);
    const int s = x.signbit() ? -1 : 1;
    dd y;
    int tail = 2;
    if (detail_math::is_tiny(a)) {
        // 2a/√π (1 - a²/3): the omitted a³ term points toward zero.
        const auto& c = detail_dd::consts();
        y = dd(a + a) / dd(c.sqrtpi[0], c.sqrtpi[1]);
        tail = -s;
    } else if (a < float256(4)) {
        y = detail_spec::erf_series_dd(a);
    } else if (a > float256(16)) {
        y = dd(float256::one());  // erfc(16) < 2^-370
        tail = -s;
    } else {
        const detail_dd::Scaled e = detail_spec::erfc_cf_dd(a);
        y = dd(float256::one()) - detail_dd::ldexp(e.m, e.k);
    }
    if (s < 0) y = -y;
    return env.round(y.hi, y.lo, 0, tail);
}

[[nodiscard]] inline float256 erfc(float256 x) {
    detail_dd::FnEnv env;
    using detail_dd::dd;
    if (x.isnan()) return env.nan(x);
    if (x.isinf()) return env.exact(x.signbit() ? float256::two() : float256::zero());
    if (x.iszero()) return env.exact(float256::one());
    const float256 a = abs(x);
    const auto& c = detail_dd::consts();
    if (x.signbit()) {
        // 1 + erf(a)
        if (a > float256(16)) return env.round(float256::two(), float256::zero(), 0, -1);
        if (detail_math::is_tiny(a))
            return (dd(float256::one()) + dd(a + a) / dd(c.sqrtpi[0], c.sqrtpi[1])).out(env, 0, -1);
        if (a < float256(4)) return (dd(float256::one()) + detail_spec::erf_series_dd(a)).out(env);
        const detail_dd::Scaled e = detail_spec::erfc_cf_dd(a);
        return (dd(float256::two()) - detail_dd::ldexp(e.m, e.k)).out(env);
    }
    if (detail_math::is_tiny(a))
        return (dd(float256::one()) - dd(a + a) / dd(c.sqrtpi[0], c.sqrtpi[1])).out(env, 0, +1);
    if (a < float256(4)) return (dd(float256::one()) - detail_spec::erf_series_dd(a)).out(env);
    if (a > float256(520)) return env.underflow(false);  // erfc(520) < 2^-390000
    const detail_dd::Scaled e = detail_spec::erfc_cf_dd(a);
    return e.m.out(env, e.k);
}

[[nodiscard]] inline float256 asinh(float256 x) {
    detail_dd::FnEnv env;
    using detail_dd::dd;
    if (x.isnan()) return env.nan(x);
    if (x.isinf() || x.iszero()) return env.exact(x);
    const float256 a = abs(x);
    if (detail_math::is_tiny(a)) return detail_math::tiny_odd(env, x, -1, 6, +1);  // x - x^3/6 (+3x^5/40)
    dd y;
    if (a > ldexp(float256::one(), 300)) {
        y = detail_dd::log(dd(a)) + detail_dd::ln2();  // log(2a); 1/(4a^2) is negligible
    } else {
        const dd r = detail_dd::sqrt(detail_dd::two_prod(a, a) + dd(float256::one()));
        y = detail_dd::log(dd(a) + r);
    }
    if (x.signbit()) y = -y;
    return y.out(env);
}

[[nodiscard]] inline float256 acosh(float256 x) {
    detail_dd::FnEnv env;
    using detail_dd::dd;
    if (x.isnan()) return env.nan(x);
    if (x < float256::one()) return env.invalid();
    if (x == float256::one()) return env.exact(float256::zero());
    if (x.isinf()) return env.exact(x);
    if (x > ldexp(float256::one(), 300)) return (detail_dd::log(dd(x)) + detail_dd::ln2()).out(env);
    const dd xm1 = detail_dd::two_sum(x, -float256::one());
    const dd xp1 = detail_dd::two_sum(x, float256::one());
    return detail_dd::log(dd(x) + detail_dd::sqrt(xm1 * xp1)).out(env);
}

[[nodiscard]] inline float256 atanh(float256 x) {
    detail_dd::FnEnv env;
    using detail_dd::dd;
    if (x.isnan()) return env.nan(x);
    if (x.iszero()) return env.exact(x);
    const float256 a = abs(x);
    if (a > float256::one()) return env.invalid();
    if (a == float256::one()) return env.pole(x.signbit());
    if (detail_math::is_tiny(a)) return detail_math::tiny_odd(env, x, 1, 3, +1);  // x + x^3/3 (+x^5/5)
    // (log(1 + a) - log(1 - a)) / 2, with 1 ± a exact as dd.
    const dd lp = detail_dd::log(detail_dd::two_sum(float256::one(), a));
    const dd lm = detail_dd::log(detail_dd::two_sum(float256::one(), -a));
    dd y = detail_dd::ldexp(lp - lm, -1);
    if (x.signbit()) y = -y;
    return y.out(env);
}


namespace detail_gamma {

// Stirling series in dd for ln Γ(w), w >= 128:
//   (w-1/2) ln w - w + ln(2π)/2 + Σ_{k=1}^{60} B_{2k} / (2k(2k-1) w^{2k-1})
// The first omitted term is below 2^-508 at w = 128. Coefficients generated
// with MPFR as hi/lo float256 pairs.
inline constexpr int kStirlingShift = 128;
inline constexpr int kStirlingTerms = 60;

struct StirlingTable {
    detail_dd::dd c[kStirlingTerms];
    StirlingTable() {
        static const char* const t[kStirlingTerms][2] = {
        {"0x1.55555555555555555555555555555555555555555555555555555555555p-4", "0x5.55555555555555555555555555555555555555555555555555555555554p-244"},
        {"-0xb.60b60b60b60b60b60b60b60b60b60b60b60b60b60b60b60b60b60b60b6p-12", "-0xb.60b60b60b60b60b60b60b60b60b60b60b60b60b60b60b60b60b60b60b6p-252"},
        {"0x3.4034034034034034034034034034034034034034034034034034034034p-12", "0x3.4034034034034034034034034034034034034034034034034034034034p-252"},
        {"-0x2.7027027027027027027027027027027027027027027027027027027027p-12", "-0x2.7027027027027027027027027027027027027027027027027027027027p-252"},
        {"0x3.72a3c5631fe46ae1d4e700dca8f158c7f91ab87539c0372a3c5631fe46ap-12", "0xe.1d4e700dca8f158c7f91ab87539c0372a3c5631fe46ae1d4e700dca8f18p-252"},
        {"-0x7.daac36664f1f207daac36664f1f207daac36664f1f207daac36664f1f2p-12", "-0x7.daac36664f1f207daac36664f1f207daac36664f1f207daac36664f1f2p-252"},
        {"0x1.a41a41a41a41a41a41a41a41a41a41a41a41a41a41a41a41a41a41a41a4p-8", "0x1.a41a41a41a41a41a41a41a41a41a41a41a41a41a41a41a41a41a41a41a4p-248"},
        {"-0x7.90a1b2c3d4e5f708192a3b4c5d6e7f90a1b2c3d4e5f708192a3b4c5d6e8p-8", "0x6.f5e4d3c2b1a08f7e6d5c4b3a291806f5e4d3c2b1a08f7e6d5c4b3a2918p-252"},
        {"0x2.dfd2c703c0cfff430edfd2c703c0cfff430edfd2c703c0cfff430edfd2cp-4", "0x7.03c0cfff430edfd2c703c0cfff430edfd2c703c0cfff430edfd2c703c0cp-244"},
        {"-0x1.6476701181f39edbdb9ce625987d4c0e9162983e954e04476701181f39fp+0", "0x2.4246319da6782b3f16e9d67c16ab1fbb898fee7e0c6124246319da6782cp-240"},
        {"0xd.672219167002d3a7a9c886459c00b4e9ea72219167002d3a7a9c88645ap+0", "-0x3.ff4b16158dde6e98ffd2c58563779ba63ff4b16158dde6e98ffd2c58564p-236"},
        {"-0x9.cd9292e6660d55b3f712eb9e07ca39db44a9292e6660d55b3f712eb9e08p+4", "0x3.5c624bb56d6d1999f2aa4c08ed1461f835c624bb56d6d1999f2aa4c08eep-236"},
        {"0x8.911a740da740da740da740da740da740da740da740da740da740da740d8p+8", "0x2.740da740da740da740da740da740da740da740da740da740da740da740ep-228"},
        {"-0x8.d0cc570e255bf59ff6eec24b48ff1b394d92e2fd250842c8d541cdbcbc8p+12", "0x1.74a2142efe472ff95fb55e66930a0c563a1c07c27258304d128c801387ap-224"},
        {"0xa.8d1044d3708d1c219ee4fdc4469ccaedcb00698234d582b96093d5ee56p+16", "0x2.c0846967b29ee369d34f6a57da7b74683dca822e34f5118cb084d3708d2p-220"},
        {"-0xe.8844d8a169abbc406169abbc406169abbc406169abbc406169abbc40618p+20", "0x1.65443bf9e965443bf9e965443bf9e965443bf9e965443bf9e965443bf9fp-216"},
        {"0x1.6d29a0f6433b79890cede62433b79890cede62433b79890cede62433b7ap+28", "-0x7.6f31219dbcc48676f31219dbcc48676f31219dbcc48676f31219dbcc488p-212"},
        {"-0x2.88a233b3c8cddaba9809357125d858d36a2301e4fc1e87e3cc19adfa948p+32", "-0x7.17376ec4c938f344120fe6f9bd78c3ebf39d089684b6284ee4a3280146p-208"},
        {"0x5.0dde6f27500939a85c40939a85c40939a85c40939a85c40939a85c40938p+36", "0x1.a85c40939a85c40939a85c40939a85c40939a85c40939a85c40939a85c4p-200"},
        {"-0xb.4005bde03d4642a243581714af642a243581714af642a243581714af64p+40", "-0x2.a243581714af642a243581714af642a243581714af642a243581714af64p-196"},
        {"0x1.bc8cd6f8f1f755c78753cdb5d5c8db31ddd9593d0ae7da3d55812a5a766p+48", "-0x4.6145057b62d9e0dd77d27807d376ee5f9190f6b559740056a717ec9c2ecp-192"},
        {"-0x4.bbebb143bb94de5a0284fa7ec423225b901be2f0a8b30da65eab49b3d48p+52", "0xe.2d1bc40435e8c5d410a2a009633521fa5fd7b05813bdcdda46fe41d0f58p-188"},
        {"0xe.2e1337f5af0bed90b6b0a352d4f335c83da6597d3226a6f46ba523a04cp+56", "0x2.cb444e986e10dd6c40121349260c05f8a82a488b1d92fbaed2877bfc49cp-180"},
        {"-0x2.e78250162b62405ad3e4bfe61b37ba8827636206b515c5587d429f4c3d2p+64", "-0xc.ea3248261c7e9df0eb8d99c56ceb792b8ee0a8d5afa3634d987e8982b5p-176"},
        {"0xa.5f7eef9e71ac7c80326ab4cc8bf3f7c478f4715b08640e90b3d95ed5188p+68", "0xa.0d9d86ba7f9aad3322fcfdf11e3d1c56c21903a42cf657b546228367618p-172"},
        {"-0x2.83be0395e550213369924971b2192469a537bbc232e73af41c5cc3cf144p+76", "-0xf.e2666cdd91e59ec7076e79befa8d1117883c904971b2192469a537bbc2p-164"},
        {"0xa.8ebfe48da17dd999790760b0ce0256ec758797bf482690066ba7710d48p+80", "0x2.31e832502fffebcc75507f8adddf7f1db5ead95a82d0fe1676439a588c6p-156"},
        {"-0x2.fd60a90cd55bec5c93256ad5b2fb0bbcee86d7bc4a2c51e09026dcdf1e6p+88", "0x4.d04f44311792843b5d3ae1f6fd92320e1a4d04f44311792843b5d3ae1f8p-152"},
        {"0xe.940b3723e6c7d0e7770e67104316dcb45c0eb78e86d69d5edfbfdaa6338p+92", "-0x1.788d7c82d447c216474d29b0942f66218b6e5a2e075bc7436b4eaf6fdffp-144"},
        {"-0x4.c4d0a8344b39331fc661da7a28c1af08ba58c64deb053ea5ae39f40e6bp+100", "0x1.88f22f0fcf32aaca751a436dc3701920cf78d310485cb06c6cb43e08dc9p-136"},
        {"0x1.ab9d9181cb599fdab1f7db2b872005987b1a9448be405987b1a9448be4p+108", "0x5.987b1a9448be405987b1a9448be405987b1a9448be405987b1a9448be4p-132"},
        {"-0xa.00536caa55d946e12e5b9c2738f0d8d4ccac2a22017f7756d4ccac2a22p+112", "-0x1.7f7756d4ccac2a22017f7756d4ccac2a22017f7756d4ccac2a22017f775p-128"},
        {"0x3.fd729201c63411e014e64015789cb3b304d797668630918d0039d618a58p+120", "-0x3.9c98cfb28689979cf6e72ffc629e75a839c98cfb28689979cf6e72ffc62p-120"},
        {"-0x1.b1959bd1fa84a486a39cca5b435920df12e1b004135bf4cac31ed2676b4p+128", "0x2.9c66e8d1e4ffbeca40b353ce12d9894c29c66e8d1e4ffbeca40b353ce12p-112"},
        {"0xc.375bbf4bd004bf46d260a614ffdd36d7d64c29cc515283c3609112f8d1p+132", "-0x3.16f72cd63c31fb4da741e41ef10e764dad56896b61cba5715463ba9c67p-104"},
        {"-0x5.d699b10dfd929649cec7478b59bcc9b2cb12a3bba682049979815d954f8p+140", "0x4.c27d340be732a46131d4846a2189738b23e61ca142781277262eb73d14cp-100"},
        {"0x2.f3fd84838e482d17d43f5875cbfa41aaec3d9804ad012b404ad012b404ap+148", "0xd.012b404ad012b404ad012b404ad012b404ad012b404ad012b404ad012b8p-92"},
        {"-0x1.94155f61c1a55b39f1b8f3b904639efd47fc83ca46164ce9ed57275dfbp+156", "0x1.97c790f3f086b67fe2df75a56ed1ce64945dc0bd53834cafac42723b9eap-84"},
        {"0xe.3e6a8d1ceca9d6a978134fc5fbd6df4767c95d95faac26fe3271a69e848p+160", "-0x2.9e370a3326ebe6232d0bb8076151985c7aa0a6e2686d501482b115bb208p-76"},
        {"-0x8.76f6ffd5a2b432dedba5116dcfeff5687f8b38341b4d93f959d1faa42dp+168", "-0x3.2ef24c494e529790564499e2187031287728a40a45a1490bbc7ff589db2p-68"},
        {"0x5.4b17eb5a4753454551d491d2c2e4dea47b4e533415f512bcd7b9d5a1b18p+176", "-0x1.e6e4cca868489a3e044a4727d1306c3e65c53b8aab12a548f1f979d671ap-60"},
        {"-0x3.7a921a734644fd4de88265ae6aa402cb7e2eeebb6a06cc16ec734baf93ap+184", "-0xa.e39173b6270d83343403611d54eb311766146b58e2cae4baa6bb19d9908p-60"},
        {"0x2.665d4d78b694ad10d0ca2772b0653133498ee569155c918069a9a9a9a9ap+192", "0x9.a9a9a9a9a9a9a9a9a9a9a9a9a9a9a9a9a9a9a9a9a9a9a9a9a9a9a9a9a98p-48"},
        {"-0x1.bc5e774c8df222e95d0264a98f8dbeb539d15d51914e7895c7da8d75c71p+200", "0x3.f1920262d2ecb5146106ee0efee204a9baaf31cde68fe5ccdc22169a51p-40"},
        {"0x1.509fec9a4b705bfe0a40efe4bb2df93fb80541cfd46b1d495bc6db0472p+208", "0x2.4e99535d73852516ef22414bbd57b7b8eccab5e9ce510b0658ce4a46b42p-32"},
        {"-0x1.0acbb1e7c7b94abd2964a0bdc566096b9945d6d26fca858fdf1ad3eec94p+216", "-0x3.ba5c54118274bfceb61c6738d5983e47fa492cbe5ef6468a8956e1af778p-24"},
        {"0xd.d021705199e3afaa551f8fc5e912aaba838c252d3c18fb73f758a18b99p+220", "0x3.ba43e64ee90f993ba43e64ee90f993ba43e64ee90f993ba43e64ee90f9ap-16"},
        {"-0xb.f2b7bc7a7979543c792961104e11d65dbe02d2ff8e925d1105f4f9425p+228", "0x2.8ccd5df4acb261afb993aa008769a9688c5c2666f78e464bfb49532d43ap-12"},
        {"0xa.c8296959290aae445070c1e8417442d7813e4ea540f7fac096e5201debp+236", "-0x1.7fc7b2958da2bfc7b2958da2bfc7b2958da2bfc7b2958da2bfc7b2958dap+0"},
        {"-0xa.242ac78f72ab7455d81e2c0defbd028f327a9138da75eac36a513488ef8p+244", "0x2.de13f59895f3f604e2c97efcefd87fff081441f64900a247ba4db15f9dcp+8"},
        {"0x9.ef20a47f622168fea2cac051ad6e842d9faddc72a871b9cc4dc5303f648p+252", "0x1.06ced85d8623fafe0e85a2afa835cb63fef79c9853797452fb75c05db89p+16"},
        {"-0xa.20540300be3e995ba474e3fcce947ac4f36642ffc9c09cd7f1b1b0081bp+260", "0xd.47c322b22aafde6f8766025b753404e3bc1bc63ce1dbe7407e41024f74p+20"},
        {"0xa.bba57a497305882ca5053aa252000c60ff6fa390bed14118ec2adb6bf5p+268", "-0x3.add2b04f5a63d14dd444252801278c7d90d1dfa8bebbc57490b4b529d3cp+32"},
        {"-0xb.d1e97d63da606f88c9f1fdadb10e1b980d77e7071871d5331260722f518p+276", "-0x3.1249580d47a16fbabed251ffa8bc7f3e48b3a80d69ba480f979ad85781p+40"},
        {"0xd.83e7996aea5faa47f59f90214fc5e2bbde3b066fd60084c921738b1be5p+284", "-0xf.c5800aeffd059f39c5833601f5f09f8d615372d657531f71d5ae9858cep+44"},
        {"-0x1.008b8827e6038e132cf0ece68d643720e2e34e51936d43861a1b36361a3p+296", "-0x4.bbbae529cfe23d5c9ae8441a32ecfd65dcbc07bf9c44e4b586fc3789d18p+56"},
        {"0x1.3b93a2b75bbadae5cc0d4a0ba39bb86d4bf8f0afc5da1ca1a8ef2190074p+304", "-0x4.1f8517e5cf229fec83dfec700d3bbad010af7c3dfec700d3bbad010af7cp+64"},
        {"-0x1.923e37172ee0cbb8e0b8989ab428815df2f5eabeb24d875bf6277aeee0cp+312", "-0x3.ab20d0d9d8628f4ea473299f1d1857501b3f0d663764a1b79018d4e11aap+72"},
        {"0x2.12f06458024339ec28dede411ba46f76b0ed8493da01c91d094320407fcp+320", "-0x7.167a8ab1dd7cf72a102d3a2a4d4f7fa2a4d4f7fa2a4d4f7fa2a4d4f7fa4p+80"},
        {"-0x2.d54ab928061f53f9491e8c85ceae0920d682d771106bf63851f84ca0af6p+328", "0xf.b266a9673858176dca5db98b2dd80b8e1b3defd54910db4d440168e4b38p+88"},
        };
        for (int k = 0; k < kStirlingTerms; ++k)
            c[k] = detail_dd::dd(detail_dd::hexc(t[k][0]), detail_dd::hexc(t[k][1]));
    }
};
[[nodiscard]] inline const StirlingTable& stirling_table() {
    static const StirlingTable tab = [] { const fenv::scoped_default_env g; return StirlingTable{}; }();
    return tab;
}

[[nodiscard]] inline detail_dd::dd stirling_dd(const detail_dd::dd& w) {
    using detail_dd::dd;
    const auto& cs = detail_dd::consts();
    const dd lw = detail_dd::log(w);
    dd sum = (w - dd(float256::half())) * lw - w + dd(cs.half_ln_2pi[0], cs.half_ln_2pi[1]);
    const dd inv = dd(float256::one()) / w;
    const dd inv2 = inv * inv;
    dd p = inv;
    const float256 tiny = ldexp(abs(sum.hi), -500);
    const StirlingTable& tab = stirling_table();
    for (int k = 0; k < kStirlingTerms; ++k) {
        const dd term = tab.c[k] * p;
        sum = sum + term;
        if (abs(term.hi) < tiny) break;
        p = p * inv2;
    }
    return sum;
}

// ln Γ(z) for a positive dd z: shift to z + n >= 128 and divide out the
// rising product z (z+1) … (z+n-1), itself kept as dd.
[[nodiscard]] inline detail_dd::dd lgamma_pos_dd(const detail_dd::dd& z) {
    using detail_dd::dd;
    if (z.hi >= float256(kStirlingShift)) return stirling_dd(z);
    const int n = kStirlingShift - static_cast<int>(floor(z.hi).to_double());
    dd prod = z;
    for (int i = 1; i < n; ++i) prod = prod * (z + dd(float256::from_int(i)));
    return stirling_dd(z + dd(float256::from_int(n))) - detail_dd::log(prod);
}

// sin(π g) for g in (0, 1/2], as dd. Positive.
[[nodiscard]] inline detail_dd::dd sin_pi_dd(float256 g) {
    detail_dd::dd s, c;
    if (g <= float256(0.25)) {
        detail_dd::sincos_reduced(detail_dd::pi() * g, s, c);
        return s;
    }
    detail_dd::sincos_reduced(detail_dd::pi() * (float256::half() - g), s, c);
    return c;
}

// ln|Γ(x)| for finite x that is not a pole, as dd.
[[nodiscard]] inline detail_dd::dd lgamma_dd(float256 x) {
    using detail_dd::dd;
    if (!x.signbit()) return lgamma_pos_dd(dd(x));
    // ln|Γ(x)| = ln π − ln|sin(πx)| − ln Γ(1−x); 1 − x is exact as dd.
    const float256 f = x - round(x);
    const dd sp = sin_pi_dd(abs(f));
    return detail_dd::log(detail_dd::pi()) - detail_dd::log(sp) -
           lgamma_pos_dd(detail_dd::two_sum(float256::one(), -x));
}

// ln Γ(1 + e) and ln Γ(2 + e) for |e| < 2^-100 (three Taylor terms; the
// next one is below 2^-300 relative).
[[nodiscard]] inline detail_dd::dd lgamma_near_1_2(float256 e, bool at_two) {
    using detail_dd::dd;
    using detail_dd::hexc;
    const dd a1 = at_two
        ? dd(hexc("0x6.c3b981c824f385b2e41c07efead34a95e3133c509a33fe6f3fc20cb8f64p-4"),
             hexc("0x1.004271b4a605fc560f112f9b63349defa82eefa9516ecdeca5f71bc4b99p-240"))
        : -dd(hexc("0x9.3c467e37db0c7a4d1be3f810152cb56a1cecc3af65cc0190c03df347098p-4"),
              hexc("0x2.ffbd8e4b59fa03a9f0eed0649ccb621057d11056ae9132135a08e43b468p-240"));
    const float256 a2 = at_two ? hexc("0x5.28d3312983e991873d8912200bace5eb5c88daa81162d8b33582c031a6p-4")
                               : hexc("0xd.28d3312983e991873d8912200bace5eb5c88daa81162d8b33582c031a6p-4");
    const float256 a3 = at_two ? hexc("0x1.13e001a5576068125d0741df7857fb3b2ae7f3b26dbc01f45c3428e5b9dp-4")
                               : hexc("0x6.693556faacb5bd67b25c9734cdad5090803d4907c3115749b1897e3b0f4p-4");
    const float256 e2 = e * e;
    return a1 * e + dd(a2 * e2) - dd(a3 * e2 * e);
}

} // namespace detail_gamma

// Natural log of |Γ(x)|. Poles (±0 and non-positive integers) return +∞ with divide-by-zero.
[[nodiscard]] inline float256 lgamma(float256 x) {
    detail_dd::FnEnv env;
    if (x.isnan()) return env.nan(x);
    if (x.isinf()) return env.exact(float256::inf(false));
    if (x.iszero()) return env.pole(false);
    if (x.signbit() && trunc(x) == x) return env.pole(false); // non-positive integer
    if (x == float256::one() || x == float256::two()) return env.exact(float256::zero());
    const float256 tiny = ldexp(float256::one(), -100);
    // Taylor at 1 and 2; the omitted ε^4 term is positive.
    if (abs(x - float256::one()) < tiny) return detail_gamma::lgamma_near_1_2(x - float256::one(), false).out(env, 0, +1);
    if (abs(x - float256::two()) < tiny) return detail_gamma::lgamma_near_1_2(x - float256::two(), true).out(env, 0, +1);
    if (x > ldexp(float256::one(), 262100)) return env.overflow(false);  // x ln x > max
    return detail_gamma::lgamma_dd(x).out(env);
}

// Γ(x). Negative integers and −∞ are invalid; tgamma(±0) is ±∞ (divide-by-zero).
[[nodiscard]] inline float256 tgamma(float256 x) {
    detail_dd::FnEnv env;
    if (x.isnan()) return env.nan(x);
    if (x.iszero()) return env.pole(x.signbit());
    if (x.isinf()) return x.signbit() ? env.invalid() : env.exact(x);
    if (x.signbit() && trunc(x) == x) return env.invalid();
    if (x > float256(20000)) return env.overflow(false);  // Γ(20000) > 2^262144
    // Γ(n) = (n-1)! is exact when representable: dd product, exact iff lo == 0.
    if (!x.signbit() && trunc(x) == x && x <= float256(300)) {
        detail_dd::dd p{float256::one()};
        const int n = static_cast<int>(x.to_double());
        for (int i = 2; i < n; ++i) p = p * float256::from_int(i);
        return env.round(p.hi, p.lo, 0, p.lo.iszero() ? 0 : 2);
    }
    bool neg = false;
    if (x.signbit()) {
        // sign of Γ(x) = sign of sin(πx) = (-1)^k sign(x - k), k = round(x)
        const float256 k = round(x);
        const bool odd = trunc(k * float256::half()) != k * float256::half();
        neg = odd != (x - k).signbit();
        if (x < float256(-20000)) return env.underflow(neg);  // |Γ(x)| < 2^-262378
    }
    const detail_dd::dd L = detail_gamma::lgamma_dd(x);
    if (L.hi > float256(181800)) return env.overflow(neg);
    if (L.hi < float256(-182000)) return env.underflow(neg);
    detail_dd::Scaled e = detail_dd::exp_scaled(L);
    if (neg) e.m = -e.m;
    return e.m.out(env, e.k);
}

// sin and cos of the same argument, one reduction.
inline void sincos(float256 x, float256* sin_out, float256* cos_out) {
    if (sin_out == nullptr && cos_out == nullptr) return;
    if (x.isnan() || x.isinf() || x.iszero() || detail_math::is_tiny(abs(x))) {
        // Same results and flags as the separate calls.
        if (sin_out) *sin_out = sin(x);
        if (cos_out) *cos_out = cos(x);
        return;
    }
    detail_dd::FnEnv env;
    int quad = 0;
    detail_dd::dd sd, cd;
    detail_math::trig_dd(abs(x), sd, cd, quad);
    detail_dd::dd sv = (quad & 1) ? cd : sd;
    detail_dd::dd cv = (quad & 1) ? -sd : cd;
    if (quad >= 2) { sv = -sv; cv = -cv; }
    if (x.signbit()) sv = -sv;
    const float256 so = sv.out(env);  // restores the caller's environment
    if (sin_out) *sin_out = so;
    if (cos_out) *cos_out = float256::round_dd(cv.hi, cv.lo, 0, env.mode(), 2);
}

// IEEE remainder, plus the low 31 bits of the integer quotient (sign of x/y).
/// IEEE remainder plus the sign and low 31 bits of the rounded quotient.
[[nodiscard]] constexpr float256 remquo(float256 x, float256 y, int* quo) noexcept {
    return float256::remainder_impl(x, y, true, quo);
}

} // namespace base256
