// Special functions: reference values (mpmath, see tools/sf_refgen.py),
// special arguments, domain errors, IEEE flags and rounding modes.
#include <base256/base256.hpp>
#include <functional>
#include <iostream>
#include <string>

using base256::float256;
using base256::round_mode;
namespace fenv = base256::fenv;

static int g_fails = 0, g_pass = 0;
#define CHECK(cond)                                                              \
    do {                                                                         \
        if (cond) ++g_pass;                                                      \
        else { ++g_fails; std::cerr << "FAIL " << __LINE__ << "  " #cond "\n"; } \
    } while (0)

constexpr unsigned kInv = base256::fe_invalid, kDiv = base256::fe_divbyzero, kInx = base256::fe_inexact;
const round_mode kModes[5] = {round_mode::to_nearest_even, round_mode::to_nearest_away, round_mode::toward_zero,
                              round_mode::upward, round_mode::downward};

// Literals are parsed in round-to-nearest (parsing follows the current mode).
static float256 P(const std::string& s) {
    fenv::scoped_round g(round_mode::to_nearest_even);
    return float256::parse(s).value();
}
static bool same(float256 a, float256 b) {
    if (a.isnan() || b.isnan()) return a.isnan() && b.isnan();
    return a == b && a.signbit() == b.signbit();
}

static float256 eval(const std::string& op, const std::string& a, const std::string& b, const std::string& c) {
    using namespace base256;
    if (op == "j") return cyl_bessel_j(P(a), P(b));
    if (op == "y") return cyl_neumann(P(a), P(b));
    if (op == "i") return cyl_bessel_i(P(a), P(b));
    if (op == "k") return cyl_bessel_k(P(a), P(b));
    if (op == "sj") return sph_bessel(static_cast<unsigned>(std::stoul(a)), P(b));
    if (op == "sy") return sph_neumann(static_cast<unsigned>(std::stoul(a)), P(b));
    if (op == "jn") return jn(std::stoi(a), P(b));
    if (op == "yn") return yn(std::stoi(a), P(b));
    if (op == "beta") return beta(P(a), P(b));
    if (op == "zeta") return riemann_zeta(P(a));
    if (op == "hz") return hurwitz_zeta(P(a), P(b));
    if (op == "K") return comp_ellint_1(P(a));
    if (op == "E") return comp_ellint_2(P(a));
    if (op == "Pi") return comp_ellint_3(P(a), P(b));
    if (op == "F1") return ellint_1(P(a), P(b));
    if (op == "E2") return ellint_2(P(a), P(b));
    if (op == "P3") return ellint_3(P(a), P(b), P(c));
    if (op == "gl") return tgamma_lower(P(a), P(b));
    if (op == "gu") return tgamma_upper(P(a), P(b));
    if (op == "gp") return gamma_p(P(a), P(b));
    if (op == "gq") return gamma_q(P(a), P(b));
    if (op == "ibeta") return ibeta(P(a), P(b), P(c));
    if (op == "ibetac") return ibetac(P(a), P(b), P(c));
    if (op == "leg") return legendre(static_cast<unsigned>(std::stoul(a)), P(b));
    if (op == "aleg") return assoc_legendre(static_cast<unsigned>(std::stoul(a)), static_cast<unsigned>(std::stoul(b)), P(c));
    if (op == "sleg") return sph_legendre(static_cast<unsigned>(std::stoul(a)), static_cast<unsigned>(std::stoul(b)), P(c));
    if (op == "lag") return laguerre(static_cast<unsigned>(std::stoul(a)), P(b));
    if (op == "alag") return assoc_laguerre(static_cast<unsigned>(std::stoul(a)), static_cast<unsigned>(std::stoul(b)), P(c));
    if (op == "ei") return expint(P(a));
    if (op == "en") return expint(static_cast<unsigned>(std::stoul(a)), P(b));
    if (op == "si") return sinint(P(a));
    if (op == "ci") return cosint(P(a));
    if (op == "shi") return sinhint(P(a));
    if (op == "chi") return coshint(P(a));
    if (op == "li") return logint(P(a));
    if (op == "psi") return digamma(P(a));
    if (op == "pg") return polygamma(static_cast<unsigned>(std::stoul(a)), P(b));
    if (op == "lbeta") return lbeta(P(a), P(b));
    if (op == "gr") return tgamma_ratio(P(a), P(b));
    if (op == "gdr") return tgamma_delta_ratio(P(a), P(b));
    if (op == "ai") return airy_ai(P(a));
    if (op == "aip") return airy_ai_prime(P(a));
    if (op == "bi") return airy_bi(P(a));
    if (op == "bip") return airy_bi_prime(P(a));
    if (op == "her") return hermite(static_cast<unsigned>(std::stoul(a)), P(b));
    std::cerr << "unknown op " << op << "\n";
    return float256::qnan();
}

struct Ref { const char *op, *a, *b, *c, *want; };
const Ref kRefs[] = {
    {"K", "0.5", "", "", "0x1af8d55d323f78f5001ea4446dce3674ae697dc8d7c5d52cd80426770a3ap-236"},
    {"K", "0.999999", "", "", "0x1fca382344497a73eb5a294515cd163d97df3686e5b47c56fb978da2ffc7p-234"},
    {"E", "0.75", "", "", "0x151876358866b77be0765f50faf4fcd5b31769f46949c398ddf4842c5e44p-236"},
    {"E", "-0.25", "", "", "0x18bc3dacffca36c89b781536fbdf77aace187519bed3b25edc3531d65354p-236"},
    {"Pi", "0.5", "0.25", "", "0x1f4e4cdf1a3696939aa4ecedde9c6454a3083bd54a02f6b735ce542d6c4ap-236"},
    {"Pi", "0.8", "-3", "", "0x1d66708b85fada362005d5c570a00ba05cea9904c2209c45364e91d6b8c7p-237"},
    {"F1", "0.7", "1.2", "", "0x15645e7131092dee3de3a88619977a97ff3c6a865b1ba4e6525722836a31p-236"},
    {"F1", "0.3", "-25.5", "", "-0x1a18c5829f57fed46b39aefe3f652bbee1f1e79cd0e703de3f3c66602f68p-232"},
    {"E2", "0.9", "0.8", "", "0x1780d76b1e3fe956fddf954a2edf0d7f2d97cde911872771f70566e38e4ap-237"},
    {"E2", "1", "10", "", "0x1a2d13debe69ea4b39fec8af80938800253abe147a67f723f5f7ff4293bap-234"},
    {"P3", "0.6", "0.5", "1.1", "0x170dbbce416d0b4e495ad5892983c8e33d4022babe5e9ca5bbc4a6849021p-236"},
    {"P3", "0.4", "-2", "7.75", "0x1290fb3a278aad792f2215a1d9fbd1687e5bdc0887ad648ed03b69c2a3c0p-234"},
    {"gl", "2.5", "1.25", "", "0x1303eafab70f6a49de5ee0e4e006453935df1cf57293775c7a17f3f942f0p-238"},
    {"gu", "2.5", "1.25", "", "0x1083ffae99efb56b3d948b9f68631d54b13311c1f6e5e7e706f2c38f9e3dp-236"},
    {"gp", "0.5", "2", "", "0x1e8b4307d3627a06631767ca639834d36677f268934ad0a9cee41ba91302p-237"},
    {"gq", "0.5", "2", "", "0x174bcf82c9d85f99ce898359c67cb2c99880d976cb52f56311be456ecfe1p-241"},
    {"gq", "1e-20", "0.75", "", "0x101277cc5b0db0d33d28b5222e9b44cf38f6c487bb08647fe1fd0ccfbadcp-304"},
    {"gu", "1e-20", "3", "", "0x1ab91c1085c70053d2d9a2a96940d462f11a279a81aff3a6cc5adeedda4ap-243"},
    {"gp", "300", "250", "", "0x130b6ea8596daac26acf36c027111a4c0e244eac4680d68e3a79d4be2bebp-246"},
    {"gq", "300", "350", "", "0x17b8e00a4441008a1596f3bf6ae9ae24c900ab8231b361aa3792efbc23a6p-245"},
    {"gl", "40", "1e-3", "", "0x1082a1f50a644df2710567db8f103c38ac8bb24aa1d60ef8a08972dc95f7p-640"},
    {"gq", "3", "600", "", "0x1cbfa6e535e8ea450cf2cea46fb57375e624e40f4e180faed273e6e07950p-1085"},
    {"hz", "-3.5", "0.25", "", "0x1066bd23470e8b6a49b0162660f5727dd241d6e9eac70c1cd72cfcdf38ccp-244"},
    {"hz", "0.5", "2", "", "-0x13aece583f39ff71bf41a72115e5abcd114c72921c6311823a213afbdd41p-235"},
    {"hz", "-150", "0.75", "", "0x165f1b587cc51c6028f34199c413c12d530aee7c6b4b5b6d6470b30e1916p237"},
    {"hz", "-200", "50", "", "-0x1f3e4645f20d7fd7d7034e0d99f1a96c38fb05a2c2d0f2e598bf8ee98125p886"},
    {"hz", "-30.25", "7", "", "-0x12640fe4f56f40ce745ed2456c5e0aa014c17f110581ae7f88eb4cd0d5d4p-158"},
    {"hz", "2", "0.5", "", "0x13bd3cc9be45de5a4adc4d9b30118358e10acd47fc1a14450cd044204a79p-234"},
    {"hz", "3.5", "7", "", "0x1e2990b2b49fd69410e6ef47c4bb75901c4f69b1ccf60197ddf67e7ca5a3p-245"},
    {"hz", "1.0001", "3", "", "0x138789e393be41f4190c79dd336dab164b8ae96e45d591432f2554333d85p-223"},
    {"hz", "40", "0.25", "", "0x1000000000000000000000008b61313bdf929d4b549ee9e5bf2a9a6f0672p-156"},
    {"hz", "2", "1e6", "", "0x10c6f82d72bca9b5a8491f31ff1797a8b780cd7cf0cc31671bacf70c2020p-256"},
    {"hz", "1500", "1.5", "", "0x1786ed29b7f963c8ff8f5e96fa005ad0e55456795f00cfa6a909cef00e55p-1114"},
    {"zeta", "2", "", "", "0x1a51a6625307d3230e7b1224401759cbd6b911b55022c5b1666b0580634cp-236"},
    {"zeta", "3", "", "", "0x133ba004f0062138371715c59e6907f1b180b7db17493405dd149c7ab12dp-236"},
    {"zeta", "0.5", "", "", "-0x175d9cb07e73fee37e834e422bcb579a2298e52438c62304744275f7ba82p-236"},
    {"zeta", "-7.5", "", "", "0x1ac7ac41b77b834642801c7330bec3ba10c3bf96bed43e46ad4714af8a9dp-245"},
    {"zeta", "1.0001", "", "", "0x138849e27106b9a5f068fccd023e52200dd513ab9ec545d11ce35a15e6efp-223"},
    {"zeta", "-1", "", "", "-0x155555555555555555555555555555555555555555555555555555555555p-240"},
    {"zeta", "1e-30", "", "", "-0x10000000000000000000000002546cda64479c22ce59de10db702f81474ep-237"},
    {"zeta", "-201.25", "", "", "-0x17b04db2778f8414b485cfd202b78c4fdfbc423df7b8738a2ac32af5c960p+483"},
    {"zeta", "250", "", "", "0x100000000000000000000000000000000000000000000000000000000000p-236"},
    {"j", "0", "1", "", "0x187c7fdbd7b8ef839e3530c22b34cedc3cbaf7ccbfde50538b5581f88245p-237"},
    {"j", "0.5", "3", "", "0x10a460536e4369cdeea10a105e02c4d339d175bb19b68ab3bad4748baf03p-240"},
    {"j", "2.25", "17.5", "", "0x11e90ccdd217718c7f58e6dc5dcebfdb9d0294d3d9976efae66fade762bfp-239"},
    {"j", "10", "1000", "", "-0x191bef1b8729ba9c86f779422721b421d95e3a6403194a55c53953cf2e35p-242"},
    {"j", "150", "100", "", "0x139ede31fcbf6c53b8b620481a3dce2a3666524dd1ceef28d9ae84b0a7d5p-288"},
    {"y", "0", "1", "", "0x16980226f358df29efcc5d79265550eb6a378dcf2886fc1c0a09de8b0620p-240"},
    {"y", "1.75", "0.125", "", "-0x12d27997639523a381cbea362edb935e6ec16ed7340a3198ef7752f0c014p-231"},
    {"y", "3", "60", "", "-0x18464d396d8dd927d80a1f109da42be740af51fb4587a56bd01a4e9ad035p-240"},
    {"y", "0.5", "500", "", "0x1025bd880e0ed349a5ef870416c3903d80906fd332b6b8c5831a98dd6c53p-241"},
    {"i", "0", "1", "", "0x1441ce4b386c2caea109def999685cf4119e467a3fc15a41f6af801536f7p-236"},
    {"i", "2.5", "30", "", "0x1476aea133e68109e4e300d13ba2df7a602516abdda2c78cf1d829447fa8p-197"},
    {"i", "1", "1000", "", "0x14e7fd269bb9b3bdf4586f634914be05b7b5fadaaaa789ac46abcf599d90p+1200"},
    {"k", "0", "1", "", "0x1af2107c43e119d373d510b59feb07626add9584e1b15087a4be61723466p-238"},
    {"k", "0.25", "0.0625", "", "0x1a49a6f8c89a47c067e1da8295e0cb26a06b5d9e81132a11590d599a6b44p-235"},
    {"k", "7", "70", "", "0x1b4d1167a72bfd23711b3d2e3313efbbae0d85314dc2a2c12a9c459bb130p-340"},
    {"k", "1", "1000", "", "0x1913fdf0aa7e0c1fa817ba3ae25f882fe74225a321d00661627c54c55fcfp-1684"},
    {"sj", "3", "2", "", "0x1f16f77f368fbe312b6268f8d9581a93c7f689ff039ef7ffa9ac70ae35b5p-241"},
    {"sy", "2", "1.5", "", "-0x15880a088e972991acc4f3441c8ddb3839cb1bc4a72d1fa60afa386c7fa4p-236"},
    {"jn", "-3", "-2.5", "", "0x1bb98fc5e82abb20d377d0fadde40255df35beecc62ddbed9ed35ab495d0p-239"},
    {"yn", "-3", "2.5", "", "0x18319b4c1093dfe62d671ffdb32be6e7c68c6efecb13bda0c2cd721f8a40p-237"},
    {"beta", "0.5", "0.5", "", "0x1921fb54442d18469898cc51701b839a252049c1114cf98e804177d4c762p-235"},
    {"beta", "3.25", "1e-5", "", "0x1869e67e354ec82c32f8c9d4e6353bad1d5bf1e7dd44298bf37a6edafb34p-220"},
    {"beta", "2000.5", "3000.25", "", "0x1879218301e67217c8570d69842cf4e7b940ac91d918a1dfac085628536bp-5096"},
    {"ibeta", "2.5", "3.5", "0.25", "0x1ac9d438759618f7f2aa3fcae85c861403ee6aefeafcd97b6bcfc6daa460p-239"},
    {"ibeta", "0.001", "50", "0.5", "0x1fffffffffffffffeb58a617f0efef9affeca5e38a288b507323fe50d2ebp-237"},
    {"ibeta", "400", "600", "0.375", "0x1adc195d15bd1fcd36a25c2eb71fa70b106eece9dee55e563613a3bbf209p-241"},
    {"ibetac", "2.5", "3.5", "0.875", "0x12d53224bc5bec7b2c13201d7fcf46f98bb7e15a41318e24e01cc7d51bfdp-244"},
    {"leg", "5", "0.3", "", "0x161aceee0f3cb3e5753a3ec02f2f9873ffac1d29dc725c3dee78183f91e6p-238"},
    {"leg", "200", "-0.71", "", "0x17ed3263df78568352c9736fd56634c21c10eb5e96dba9c5d6c44f8b0585p-242"},
    {"leg", "7", "12.3", "", "0x10d61410641d096f148ad17fb9f482f43ba69f3b9a0a01e1818fb7988824p-206"},
    {"aleg", "6", "3", "0.4", "-0x1e123c04d94c21468d13af52d2272e806361124f88572c95228f5bfbca68p-231"},
    {"aleg", "150", "40", "-0.999", "0x15b8bc409cab7f680deaf38a9a85ac58caf96d3cef50089b8476134829c8p-38"},
    {"sleg", "5", "2", "0.7", "0x19fe0f71d07447e79d612370e5388e843172085755535c988d64ed122991p-238"},
    {"sleg", "120", "60", "2.9", "0x1296b96873f71e5f2be1537f9505e6420d486d6d57d908bcd32b09af853dp-283"},
    {"sleg", "10", "0", "1e-30", "0x14aefbf05f1304afec32ceaba4f6f69ee5452a344a0c1911d9c5d512e6d6p-236"},
    {"lag", "6", "2.5", "", "0x152f1c71c71c71c71c71c71c71c71c71c71c71c71c71c71c71c71c71c71cp-237"},
    {"lag", "400", "1000.5", "", "0x1080ea90776696306e3aa1a56ddbe1acbb7adf2c6e37ca56c500f8f36c08p+478"},
    {"alag", "10", "3", "7.25", "0x1d390c949e9199749cfa2aa85ae0b3bb96bf1c4cca7d02d5ddb8e13e6eedp-233"},
    {"her", "9", "1.1", "", "-0x1930840feb40bc6eecd8ecfcd437b660e0d1a2b0aabe6b465b3f520ba1b1p-223"},
    {"her", "300", "-17.3", "", "-0x11189bb4fd03c1512c112159ed8c053bee4cddbbb0d36f1a98674a5b49a2p+1148"},
    {"ei", "0.372", "", "", "-0x1033ce1ed625f0fd9b6f871ae626add65326b89cebe23fc2337a57d6040ep-245"},
    {"ei", "-3.5", "", "", "-0x1c8cb8abc9972f364d969c13aa13dd4a855d56ac15d749ae88ee11629ec7p-244"},
    {"ei", "250.25", "", "", "0x10d49cd8d00217342e4cd5e25700f1b0251b0b878c2da983dc3845e97345p+117"},
    {"ei", "5000", "", "", "0x123962c79035c18f5e1b6260299879a2098c92c28dd6b13549dd7422cac1p+6965"},
    {"ei", "1e-40", "", "", "-0x16e1ad10a9caeeb1ab6921e6ec30a4e38dbbdbcf36722cc6b33b673dd204p-230"},
    {"en", "1", "45.5", "", "0x1c3961401a4edffb2ac7914bca032333c73f740cd7682cd169f71f0b83bdp-308"},
    {"en", "3", "0.125", "", "0x19858949a0b7dcb149efeca22dd62c963d134028c9d0c2caf204ee59544cp-238"},
    {"en", "40", "12", "", "0x101846da169cf1ceaeb1b6a4ee18a58a0c3ec2f3736104b16680087275f1p-259"},
    {"en", "1000", "2", "", "0x11b88d8cfd04cd96e9f4daf406c3ac6659c2aa739c12ea6e6fc4ecf424d7p-249"},
    {"en", "0", "3", "", "0x10fe75dddf2074863019576d0b00b355a2b63da267aecbb2c5f42dfa153ep-242"},
    {"si", "0.75", "", "", "0x1743358da65566584b58f0fc0df85674d313af6c4bc413e1d348a23ee0a6p-237"},
    {"si", "-150.5", "", "", "-0x19080210c31fe701de07f3a8ad3739a7f5c3baac7a95506dcad20cc6d3b6p-236"},
    {"si", "1e5", "", "", "0x192205cee772ab5235640938bd2250326c906c14039946bc0f0c88c82446p-236"},
    {"ci", "3.5", "", "", "-0x1073273242139b63163f3480a1f019b8f7d7da1e965eecc07439e8838de6p-241"},
    {"ci", "90", "", "", "0x14739ae2697bc7508284fae9e50864b5dedfd20d31d881ffc2fa46b1c0abp-243"},
    {"ci", "1e40", "", "", "-0x13d94cbb44d2f741510a989b4bb67fc1f29427e3488245b0ea6e97b2e1a0p-370"},
    {"shi", "-2.5", "", "", "-0x1c650c9538d7a9affda9e313739da2fab2c2493f62ddacae0b615dca076bp-235"},
    {"shi", "700", "", "", "0x15aa94a53af48e6bc6673697e39e2f5b64539d7bbeaa16556182b566d50ap+763"},
    {"chi", "0.125", "", "", "-0x17f91b55581f79c79f63c8ae4170070befc0a50b45301cdd5d7deee8b426p-236"},
    {"chi", "65.5", "", "", "0x16674fd7556ba30674426b28a443f7eb13975ae093e9bd9a7f94c1da7cc0p-149"},
    {"li", "0.5", "", "", "-0x183c257878eee36d8a46ce3ca4763cbf041d1dfb8e1263e5fc32e9697f3dp-238"},
    {"li", "1.0000001", "", "", "-0x1f14ee3708da5aa5faf4425058f911cccd9b5fdecc432e91fc3e7f2eaedfp-233"},
    {"li", "1e30", "", "", "0x17bca23e192a6b0c865bf71a9def67463c314db4998ce6f00162b6645ca9p-143"},
    {"li", "2", "", "", "0x10b8fda7e918078599e853393f0af3deefad54b9ce533d61088126a05a6bp-236"},
    {"psi", "1.5", "", "", "0x12aed059bd608a0734f4b29ded8e8e95b8955af048e3a717bfc3e428f01ep-241"},
    {"psi", "-3.25", "", "", "0x11dd776c1aa21f7468bec30e8e07957f2d18c64eee70d83c6f9d14b64aeep-234"},
    {"psi", "1e-30", "", "", "-0x193e5939a08ce9dbd4800000012788cfc6fb618f49a37c7f01bfeb9dd9dfp-137"},
    {"psi", "1e10", "", "", "0x17069e2aa2736123e0ea2691300091bf7aed131ab0c1e751ec9b1727e1aap-232"},
    {"pg", "1", "0.5", "", "0x13bd3cc9be45de5a4adc4d9b30118358e10acd47fc1a14450cd044204a79p-234"},
    {"pg", "3", "17.25", "", "0x1bd76f17737bc701a71d9b2d23e23ca4afbc606f0097b69e12f2f9f07895p-248"},
    {"pg", "60", "2.5", "", "-0x168dd83d16e3fec7b7ea22157620ee8e257824978508cda6671c3bd26d13p-45"},
    {"lbeta", "2.5", "3.5", "", "-0x1a6a289c2a217a99525b751804b0e7795f58bfd3314796186efa84e938e3p-235"},
    {"lbeta", "1e6", "3e5", "", "-0x156e7d4810d583b00e57792e5da0ab8be26928e005780d85d47f16c9ef14p-217"},
    {"gr", "7.5", "3.25", "", "0x16f05017023fda00515b1ef457e4121e756720c926c72b9428b44d679102p-227"},
    {"gdr", "100.5", "1e-20", "", "0x1fffffffffffffffe4d0db767be128760f587c065d58eb2ff99707602844p-237"},
    {"ai", "0", "", "", "0x16b8c7962715b85ea5b5eec13993db6a5eac401ac9be7ad90848cce50883p-238"},
    {"ai", "-7.5", "", "", "0x1497f92ca01e60b5223eca905d80e5e373ed7fa7c50c9b1b85b55c6d47f3p-238"},
    {"ai", "14", "", "", "0x1c97d1d62adb7c5d3bcc897107b6ce5b277ea6ebe70e00a0caaf0b420376p-290"},
    {"ai", "200", "", "", "0x1daddbd7c2cbf200815f32de0314780f5f7c7946fc74d22b803f947dd342p-2961"},
    {"ai", "-500.25", "", "", "0x1dc628dc955a40a4814e4b70c8b8411731938197a52e23f54c8e61e00edfp-240"},
    {"aip", "2.5", "", "", "-0x1ae182ce295c4087fb88a9ad46d6e1c1f1a10786ca842cae032b1be7bee1p-242"},
    {"aip", "-31", "", "", "0x13bc9894e5ad575a7eaa85cce3d5245bfbdd9946339fce1a0ebdf8f56f54p-236"},
    {"bi", "-2", "", "", "-0x1a632a64d50c9d88029fcd9d7fd3c4d8c34402deb4eabb22a3ce4c84f45fp-238"},
    {"bi", "50", "", "", "0x1c0e2147632ee9a8c05376a3ad09e3ba2e730be017253b3c259177d32228p+101"},
    {"bip", "0", "", "", "0x1cb0c1a680c8a08b085fbf122266921eb46ceaee6d2e3db957f854cd2e96p-238"},
    {"bip", "-1000", "", "", "0x1c51c7ce3eaffd99757b6e366519025cb962731180150b18864d6b2dc157p-236"},
    {"ibetac", "30", "2", "0.0625", "0x1ffffffffffffffffffffffffffffc5c0000000000000000000000000000p-237"},
};

// Value and flags in every rounding mode.
static void expect(int line, const std::function<float256()>& f, float256 value, unsigned want) {
    for (round_mode m : kModes) {
        fenv::scoped_round g(m);
        fenv::clear_flags();
        const float256 y = f();
        const unsigned got = fenv::test_flags();
        fenv::clear_flags();
        if (same(y, value) && got == want) { ++g_pass; continue; }
        ++g_fails;
        std::cerr << "FAIL " << line << " mode " << int(m) << ": got " << y.to_string(10) << " flags " << got
                  << ", want " << value.to_string(10) << " flags " << want << "\n";
    }
}
#define EXPECT(expr, value, flags) expect(__LINE__, [&] { return (expr); }, (value), (flags))
// Expected value depending on the rounding mode.
static void expect_m(int line, const std::function<float256()>& f, const std::function<float256(round_mode)>& value,
                     unsigned want) {
    for (round_mode m : kModes) {
        fenv::scoped_round g(m);
        fenv::clear_flags();
        const float256 y = f();
        const unsigned got = fenv::test_flags();
        fenv::clear_flags();
        const float256 v = value(m);
        if (same(y, v) && got == want) { ++g_pass; continue; }
        ++g_fails;
        std::cerr << "FAIL " << line << " mode " << int(m) << ": got " << y.to_string(10) << " flags " << got
                  << ", want " << v.to_string(10) << " flags " << want << "\n";
    }
}
#define EXPECT_M(expr, valuefn, flags) expect_m(__LINE__, [&] { return (expr); }, (valuefn), (flags))
static bool down(round_mode m) { return m == round_mode::toward_zero || m == round_mode::downward; }

int main() {
    using namespace base256;
    // Reference values: round-to-nearest result within 1/2 ulp (bitwise equal),
    // the directed results bracket it and only inexact is raised.
    for (const Ref& r : kRefs) {
        const float256 want = P(r.want);
        float256 got[5];
        unsigned fl[5];
        for (int m = 0; m < 5; ++m) {
            fenv::scoped_round g(kModes[m]);
            fenv::clear_flags();
            got[m] = eval(r.op, r.a, r.b, r.c);
            fl[m] = fenv::test_flags();
            fenv::clear_flags();
        }
        const bool ok = same(got[0], want) && got[3] >= want && got[4] <= want &&
                        (got[3] == want || got[3] == nextafter(want, float256::inf(false))) &&
                        (got[4] == want || got[4] == nextafter(want, float256::inf(true)));
        bool flags_ok = true;
        for (unsigned f : fl) flags_ok = flags_ok && f == kInx;
        if (ok && flags_ok) { ++g_pass; continue; }
        ++g_fails;
        std::cerr << "FAIL ref " << r.op << "(" << r.a << ", " << r.b << (r.c[0] ? ", " : "") << r.c << ") = "
                  << got[0].to_hex_string() << " want " << want.to_hex_string() << " flags " << fl[0] << "\n";
    }

    const float256 one = float256::one(), zero = float256::zero(), inf = float256::inf(false);
    const float256 nan = float256::qnan();
    // Bessel: special arguments.
    EXPECT(cyl_bessel_j(zero, zero), one, 0);
    EXPECT(cyl_bessel_j(float256(2.5), zero), zero, 0);
    EXPECT(cyl_bessel_i(zero, zero), one, 0);
    EXPECT(cyl_bessel_i(float256(3), zero), zero, 0);
    EXPECT(cyl_neumann(float256(1), zero), -inf, kDiv);
    EXPECT(cyl_bessel_k(float256(0.5), zero), inf, kDiv);
    EXPECT(cyl_bessel_j(float256(1), inf), zero, 0);
    EXPECT(cyl_neumann(float256(1), inf), zero, 0);
    EXPECT(cyl_bessel_i(float256(1), inf), inf, 0);
    EXPECT(cyl_bessel_k(float256(1), inf), zero, 0);
    EXPECT(cyl_bessel_j(float256(-1), one), nan, kInv);
    EXPECT(cyl_bessel_j(one, float256(-1)), nan, kInv);
    EXPECT(cyl_bessel_k(one, float256(-1)), nan, kInv);
    EXPECT(cyl_bessel_j(inf, one), nan, kInv);
    EXPECT(cyl_bessel_j(nan, one), nan, 0);
    EXPECT(cyl_bessel_j(one, float256::snan()), nan, kInv);
    EXPECT(sph_bessel(0, zero), one, 0);
    EXPECT(sph_bessel(4, zero), zero, 0);
    EXPECT(sph_neumann(2, zero), -inf, kDiv);
    EXPECT(yn(3, zero), -inf, kDiv);
    EXPECT(yn(-3, zero), inf, kDiv);
    EXPECT(yn(2, float256(-1)), nan, kInv);
    EXPECT(jn(-2, zero), zero, 0);
    // Parity: J_n(-x) = (-1)^n J_n(x), J_{-n} = (-1)^n J_n in every mode.
    for (round_mode m : kModes) {
        fenv::scoped_round g(m);
        const float256 x = float256(7.25);
        CHECK(same(jn(3, -x), -[&] { fenv::scoped_round h(m == round_mode::upward ? round_mode::downward
                                                          : m == round_mode::downward ? round_mode::upward : m);
                                       return jn(3, x); }()));
        CHECK(same(jn(4, -x), jn(4, x)));
        CHECK(same(jn(-4, x), jn(4, x)));
        CHECK(same(yn(-2, x), yn(2, x)));
    }
    // Tiny argument: J_0(x) = 1 - x^2/4 is just below 1.
    EXPECT_M(cyl_bessel_j(zero, ldexp(one, -400)),
             [](round_mode m) { return down(m) ? nextafter(float256::one(), float256::zero()) : float256::one(); }, kInx);
    EXPECT_M(cyl_bessel_i(zero, ldexp(one, -400)),
             [](round_mode m) { return m == round_mode::upward ? nextafter(float256::one(), float256::inf(false)) : float256::one(); }, kInx);
    // Overflow / underflow.
    EXPECT_M(cyl_bessel_i(zero, float256(200000)),
             [](round_mode m) { return down(m) ? float256::max() : float256::inf(false); }, base256::fe_overflow | kInx);
    EXPECT_M(cyl_bessel_k(zero, float256(200000)),
             [](round_mode m) { return m == round_mode::upward ? float256::denorm_min() : float256::zero(); },
             base256::fe_underflow | kInx);

    // Beta and incomplete beta: special arguments.
    EXPECT(beta(zero, one), inf, kDiv);
    EXPECT(beta(float256(-1), one), nan, kInv);
    EXPECT(beta(inf, one), zero, 0);
    EXPECT(beta(one, float256(4)), float256(0.25), 0);
    EXPECT(beta(one, one), one, 0);
    EXPECT(beta(nan, one), nan, 0);
    EXPECT(ibeta(float256(2), float256(3), zero), zero, 0);
    EXPECT(ibeta(float256(2), float256(3), one), one, 0);
    EXPECT(ibetac(float256(2), float256(3), zero), one, 0);
    EXPECT(ibeta(zero, float256(3), float256(0.5)), nan, kInv);
    EXPECT(ibeta(one, float256(3), float256(1.5)), nan, kInv);
    EXPECT(ibeta(one, float256(3), float256(-0.5)), nan, kInv);
    // I_x(a, b) + I_{1−x}(b, a) = 1
    for (double xv : {0.1, 0.4, 0.7, 0.95}) {
        const float256 x(xv), a(2.75), b(4.5);
        CHECK(abs(ibeta(a, b, x) - ibetac(b, a, float256::one() - x)) <= ldexp(one, -236));
    }

    // Riemann zeta: special arguments.
    EXPECT(riemann_zeta(one), inf, kDiv);
    EXPECT(riemann_zeta(zero), -float256::half(), 0);
    EXPECT(riemann_zeta(-zero), -float256::half(), 0);
    EXPECT(riemann_zeta(float256(-2)), zero, 0);
    EXPECT(riemann_zeta(float256(-1000)), zero, 0);
    EXPECT(riemann_zeta(inf), one, 0);
    EXPECT(riemann_zeta(-inf), nan, kInv);
    EXPECT(hurwitz_zeta(one, float256(2)), inf, kDiv);
    EXPECT(hurwitz_zeta(zero, float256(0.25)), float256(0.25), 0);
    EXPECT(hurwitz_zeta(-inf, float256(2)), nan, kInv);
    EXPECT(hurwitz_zeta(float256(0.5), inf), -inf, 0);
    EXPECT(hurwitz_zeta(float256(2), zero), nan, kInv);
    EXPECT(hurwitz_zeta(float256(2), inf), zero, 0);
    EXPECT(hurwitz_zeta(inf, float256(0.5)), inf, 0);
    EXPECT(hurwitz_zeta(inf, one), one, 0);
    EXPECT(hurwitz_zeta(inf, float256(2)), zero, 0);
    EXPECT(riemann_zeta(nan), nan, 0);
    EXPECT_M(riemann_zeta(float256(700)),
             [](round_mode m) { return m == round_mode::upward ? nextafter(float256::one(), float256::inf(false)) : float256::one(); }, kInx);
    EXPECT_M(riemann_zeta(float256(-30001.5)),
             [](round_mode m) { return (m == round_mode::toward_zero || m == round_mode::upward) ? -float256::max() : -float256::inf(false); },
             base256::fe_overflow | kInx);

    // Incomplete gamma: special arguments.
    EXPECT(gamma_p(float256(2), zero), zero, 0);
    EXPECT(gamma_q(float256(2), zero), one, 0);
    EXPECT(gamma_p(float256(2), inf), one, 0);
    EXPECT(gamma_q(float256(2), inf), zero, 0);
    EXPECT(tgamma_lower(float256(2), zero), zero, 0);
    EXPECT(tgamma_upper(float256(4), zero), float256(6), 0);
    EXPECT(tgamma_lower(float256(4), inf), float256(6), 0);
    EXPECT(tgamma_upper(float256(4), inf), zero, 0);
    EXPECT(gamma_p(zero, one), nan, kInv);
    EXPECT(gamma_p(float256(-1), one), nan, kInv);
    EXPECT(gamma_q(one, float256(-1)), nan, kInv);
    EXPECT(gamma_q(nan, one), nan, 0);
    for (double xv : {0.25, 1.0, 3.5, 12.0}) {  // P + Q = 1, γ + Γ = Γ(a)
        const float256 a(2.75), x(xv);
        CHECK(abs(gamma_p(a, x) + gamma_q(a, x) - one) <= ldexp(one, -235));
        CHECK(abs((tgamma_lower(a, x) + tgamma_upper(a, x)) / tgamma(a) - one) <= ldexp(one, -234));
    }

    // Elliptic integrals: special arguments.
    EXPECT(comp_ellint_1(one), inf, kDiv);
    EXPECT(comp_ellint_1(float256(1.5)), nan, kInv);
    EXPECT(comp_ellint_2(-one), one, 0);
    EXPECT(comp_ellint_2(float256(-1.5)), nan, kInv);
    EXPECT(comp_ellint_3(float256(0.5), one), inf, kDiv);
    EXPECT(comp_ellint_3(float256(0.5), float256(2)), nan, kInv);
    EXPECT(ellint_1(zero, float256(3.25)), float256(3.25), 0);
    EXPECT(ellint_2(zero, float256(-3.25)), float256(-3.25), 0);
    EXPECT(ellint_1(float256(0.5), zero), zero, 0);
    EXPECT(ellint_1(float256(0.5), -zero), -zero, 0);
    EXPECT(ellint_1(float256(0.5), inf), inf, 0);
    EXPECT(ellint_1(float256(2), one), nan, kInv);
    EXPECT(ellint_1(one, float256(2)), inf, kDiv);
    EXPECT(ellint_3(float256(0.5), float256(4), one), nan, kInv);  // 1 − ν sin²φ < 0
    EXPECT(ellint_1(nan, one), nan, 0);
    EXPECT_M(ellint_1(float256(0.5), ldexp(one, -300)),
             [](round_mode m) { return m == round_mode::upward ? nextafter(ldexp(float256::one(), -300), float256::inf(false))
                                                                : ldexp(float256::one(), -300); }, kInx);
    EXPECT_M(ellint_2(float256(0.5), ldexp(one, -300)),
             [](round_mode m) { return down(m) ? nextafter(ldexp(float256::one(), -300), float256::zero())
                                               : ldexp(float256::one(), -300); }, kInx);
    // Odd in φ, F(φ + π) = F(φ) + 2K.
    CHECK(same(ellint_1(float256(0.6), float256(-1.3)), -ellint_1(float256(0.6), float256(1.3))));
    CHECK(abs(ellint_2(float256(0.6), float256(0.25) + float256::pi()) - ellint_2(float256(0.6), float256(0.25)) -
              float256::two() * comp_ellint_2(float256(0.6))) <= ldexp(one, -232));

    // Orthogonal polynomials: exact values, signs of zero, infinities, domains.
    EXPECT(legendre(0, nan), nan, 0);
    EXPECT(legendre(0, float256(5)), one, 0);
    EXPECT(legendre(9, one), one, 0);
    EXPECT(legendre(9, -one), -one, 0);
    EXPECT(legendre(2, zero), float256(-0.5), 0);
    EXPECT(legendre(3, -zero), -zero, 0);
    EXPECT(legendre(3, -inf), -inf, 0);
    EXPECT(legendre(4, -inf), inf, 0);
    EXPECT(hermite(3, float256(0.5)), float256(-5), 0);  // 8x³ − 12x
    EXPECT(hermite(5, -zero), -zero, 0);
    EXPECT(hermite(4, -inf), inf, 0);
    EXPECT(laguerre(2, float256(2)), float256(-1), 0);    // (x² − 4x + 2)/2
    EXPECT(laguerre(3, inf), -inf, 0);
    EXPECT(laguerre(3, -inf), inf, 0);
    EXPECT(assoc_laguerre(0, 7, float256(3)), one, 0);
    EXPECT(assoc_laguerre(1, 2, float256(3)), zero, 0);  // 1 + m − x
    EXPECT(assoc_legendre(2, 3, float256(0.5)), zero, 0);
    EXPECT(assoc_legendre(3, 1, one), zero, 0);
    EXPECT(assoc_legendre(2, 1, float256(1.5)), nan, kInv);
    EXPECT(assoc_legendre(2, 2, float256(0.5)), float256(2.25), 0);  // 3(1 − x²)
    EXPECT(sph_legendre(3, 1, inf), nan, kInv);
    EXPECT(sph_legendre(3, 4, one), zero, 0);
    EXPECT(sph_legendre(3, 1, zero), -zero, 0);
    EXPECT_M(hermite(20000, float256(1e100)),
             [](round_mode m) { return down(m) ? float256::max() : float256::inf(false); }, base256::fe_overflow | kInx);
    // P_l(−x) = (−1)^l P_l(x), H_n likewise.
    for (unsigned l : {4u, 7u, 30u}) {
        CHECK(same(legendre(l, float256(-0.37)), (l & 1) ? -legendre(l, float256(0.37)) : legendre(l, float256(0.37))));
        CHECK(same(hermite(l, float256(-2.5)), (l & 1) ? -hermite(l, float256(2.5)) : hermite(l, float256(2.5))));
    }

    // Exponential integrals.
    EXPECT(expint(zero), -inf, kDiv);
    EXPECT(expint(-zero), -inf, kDiv);
    EXPECT(expint(inf), inf, 0);
    EXPECT(expint(-inf), -zero, 0);
    EXPECT(expint(nan), nan, 0);
    EXPECT(expint(1u, zero), inf, kDiv);
    EXPECT(expint(0u, zero), inf, kDiv);
    EXPECT(expint(3u, zero), float256(0.5), 0);
    EXPECT(expint(5u, inf), zero, 0);
    EXPECT(expint(2u, -one), nan, kInv);
    EXPECT(expint(1u, -inf), -inf, 0);
    EXPECT_M(expint(float256(190000)),
             [](round_mode m) { return down(m) ? float256::max() : float256::inf(false); }, base256::fe_overflow | kInx);
    EXPECT_M(expint(float256(-185000)),
             [](round_mode m) { return (m == round_mode::downward) ? -float256::denorm_min() : -float256::zero(); },
             base256::fe_underflow | kInx);
    CHECK(same(expint(1u, float256(-2.5)), -expint(float256(2.5))));

    // Sine/cosine integrals and li.
    EXPECT(sinint(-zero), -zero, 0);
    EXPECT(sinint(nan), nan, 0);
    {
        fenv::clear_flags();
        CHECK(same(sinint(-inf), -(float256::pi() * float256::half())) && fenv::test_flags() == kInx);
        fenv::clear_flags();
    }
    EXPECT(cosint(zero), -inf, kDiv);
    EXPECT(cosint(-one), nan, kInv);
    EXPECT(cosint(inf), zero, 0);
    EXPECT(sinhint(-inf), -inf, 0);
    EXPECT(sinhint(zero), zero, 0);
    EXPECT(coshint(zero), -inf, kDiv);
    EXPECT(coshint(-one), nan, kInv);
    EXPECT(coshint(inf), inf, 0);
    EXPECT_M(sinhint(float256(-190000)),
             [](round_mode m) { return (m == round_mode::toward_zero || m == round_mode::upward) ? -float256::max() : -float256::inf(false); },
             base256::fe_overflow | kInx);
    EXPECT(logint(zero), zero, 0);
    EXPECT(logint(one), -inf, kDiv);
    EXPECT(logint(-one), nan, kInv);
    EXPECT(logint(inf), inf, 0);
    CHECK(same(sinhint(float256(-3.5)), -sinhint(float256(3.5))));
    CHECK(same(sinint(float256(-200.25)), -sinint(float256(200.25))));

    // Digamma / polygamma / lbeta / gamma ratios.
    EXPECT(digamma(zero), -inf, kDiv);
    EXPECT(digamma(-zero), inf, kDiv);
    EXPECT(digamma(float256(-3)), nan, kInv);
    EXPECT(digamma(inf), inf, 0);
    EXPECT(digamma(-inf), nan, kInv);
    EXPECT(digamma(nan), nan, 0);
    EXPECT(polygamma(1, zero), inf, kDiv);
    EXPECT(polygamma(2, zero), -inf, kDiv);
    EXPECT(polygamma(2, -one), nan, kInv);
    EXPECT(polygamma(1, inf), zero, 0);
    EXPECT(polygamma(2, inf), -zero, 0);
    CHECK(same(trigamma(float256(3.5)), polygamma(1, float256(3.5))));
    EXPECT(lbeta(one, one), zero, 0);
    EXPECT(lbeta(zero, one), inf, kDiv);
    EXPECT(lbeta(-one, one), nan, kInv);
    EXPECT(lbeta(inf, one), -inf, 0);
    EXPECT(tgamma_ratio(float256(2.5), float256(2.5)), one, 0);
    EXPECT(tgamma_ratio(one, zero), zero, 0);
    EXPECT(tgamma_ratio(zero, one), inf, kDiv);
    EXPECT(tgamma_delta_ratio(float256(3), zero), one, 0);
    EXPECT(tgamma_delta_ratio(float256(3), float256(-4)), nan, kInv);

    // Airy.
    EXPECT(airy_ai(inf), zero, 0);
    EXPECT(airy_ai(-inf), zero, 0);
    EXPECT(airy_ai_prime(inf), -zero, 0);
    EXPECT(airy_ai_prime(-inf), nan, kInv);
    EXPECT(airy_bi(inf), inf, 0);
    EXPECT(airy_bi_prime(-inf), nan, kInv);
    EXPECT(airy_ai(nan), nan, 0);
    EXPECT_M(airy_bi(float256(5000)),
             [](round_mode m) { return down(m) ? float256::max() : float256::inf(false); }, base256::fe_overflow | kInx);
    EXPECT_M(airy_ai(float256(5000)),
             [](round_mode m) { return m == round_mode::upward ? float256::denorm_min() : float256::zero(); },
             base256::fe_underflow | kInx);

    std::cout << g_pass << " passed, " << g_fails << " failed\n";
    return g_fails ? 1 : 0;
}
