// Latency of dependent chains of the default-path operators (s = s + x[i] * y[i], ...).
// bench/bench.cpp measures independent calls (throughput); this catches a lost
// inlining of the operators, as in the first 4.26.0 draft (CHANGELOG 4.26.0).
// build: g++-14 -O3 -DNDEBUG -std=c++23 -Iinclude bench/chain.cpp -o chain
// usage: chain muladd|add|mul|div|fma|sqrt   (prints ns per iteration and the
// leading hex digits of the result, which must agree between builds)
#include <base256/float256.hpp>
#include <chrono>
#include <cstdio>
#include <cstring>
using namespace base256;
int main(int argc, char** argv) {
  const char* op = argc > 1 ? argv[1] : "muladd";
  float256 x[64], y[64];
  for (int i = 0; i < 64; ++i) { x[i] = float256(1.0 + i * 0.013) / float256(3.0); y[i] = float256(2.0 - i * 0.007) / float256(7.0); }
  float256 s = float256(0.5);
  const long R = std::strcmp(op, "sqrt") == 0 || std::strcmp(op, "div") == 0 ? 8000 : 40000;
  auto t0 = std::chrono::steady_clock::now();
  for (long r = 0; r < R; ++r) for (int i = 0; i < 64; ++i) {
    if (!std::strcmp(op, "muladd")) s = s + x[i] * y[i];
    else if (!std::strcmp(op, "add")) s = s + x[i];
    else if (!std::strcmp(op, "mul")) s = s * (x[i] + float256(0.75)) * y[i] * float256(2.6);
    else if (!std::strcmp(op, "div")) s = (s + x[i]) / y[i] / float256(5.0);
    else if (!std::strcmp(op, "fma")) s = fma(s, y[i], x[i]);
    else if (!std::strcmp(op, "sqrt")) s = sqrt(s + x[i]);
  }
  double el = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  std::printf("%s %.2f ns/iter %s\n", op, el / (R * 64.0) * 1e9, s.to_hex_string().substr(0, 12).c_str());
}
