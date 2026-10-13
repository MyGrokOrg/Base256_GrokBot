// Compile/run test: the library must coexist with <math.h>/<cmath>
// classification functions included before and after it.
#include <math.h>
#include <base256/base256.hpp>
#include <cmath>

int main() {
    using base256::float256;
    int bad = 0;
    // C/C++ classification on built-in types still works, both spellings.
    bad += isnan(1.0) ? 1 : 0;
    bad += std::isnan(NAN) ? 0 : 1;
    bad += signbit(-1.0) ? 0 : 1;
    bad += std::isinf(HUGE_VAL) ? 0 : 1;
    bad += fpclassify(0.0) == FP_ZERO ? 0 : 1;
    // float256 overloads via ADL and qualified names.
    const float256 n = float256::qnan();
    bad += isnan(n) ? 0 : 1;
    bad += base256::isnan(n) ? 0 : 1;
    bad += n.isnan() ? 0 : 1;
    bad += signbit(-float256::one()) ? 0 : 1;
    bad += isfinite(float256::one()) ? 0 : 1;
    bad += fpclassify(float256::zero()) == FP_ZERO ? 0 : 1;
    return bad;
}
