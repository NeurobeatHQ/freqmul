// FreqMul on ESP-IDF: a 10 Hz tone through each method, each at its own internal rate, interpolated to 48 kHz.
// Prints the measured output frequency (PASS/FAIL) and the CPU share of real time on this core.
#include <math.h>
#include <stdio.h>

#include "FreqMul.h"
#include "esp_timer.h"

namespace {

constexpr float kInRate = 250.f, kOutRate = 48000.f;
constexpr int kOutBlock = 192;  // 48000 / 250

fm::Upsampler raw;
fm::FreqShifter shifter;
fm::GrainScaler grain;
fm::SpectralBank bank;
fm::Interpolator up;

enum Method { kRaw, kShift, kGrain, kBank };

float runMethod(Method m, float rate, float tone, float seconds, float* cpuPct) {
  const int factor = (int)(kOutRate / rate), n = kOutBlock / factor;
  up.begin(factor);
  fm::ZeroCrossHz zc;
  zc.begin(kOutRate);
  float in[64], out[kOutBlock];
  int64_t busy = 0;
  const int total = (int)(seconds * kInRate), settle = (int)(4 * kInRate);
  for (int i = 0; i < total; ++i) {
    float x = 0.25f * sinf(fm::kTwoPi * tone * i / kInRate);
    int64_t t0 = esp_timer_get_time();
    switch (m) {
      case kRaw: raw.push(x); raw.render(in, n); break;
      case kShift: shifter.push(x); shifter.render(in, n); break;
      case kGrain: grain.push(x); grain.render(in, n); break;
      case kBank: bank.push(x); bank.update(); bank.render(in, n); break;
    }
    up.process(in, n, out);
    int64_t dt = esp_timer_get_time() - t0;
    if (i >= settle) {
      busy += dt;
      for (int k = 0; k < kOutBlock; ++k) zc.add(out[k]);
    }
  }
  *cpuPct = 100.f * (float)busy / (1e6f * (float)(total - settle) / kInRate);
  return zc.hz();
}

void report(const char* name, float hz, float want, float cpu) {
  printf("%s%.2f Hz (want %.2f) %s   cpu %.1f %%\n", name, hz, want, fabsf(hz - want) <= want * 0.006f ? "PASS" : "FAIL", cpu);
}

}  // namespace

extern "C" void app_main(void) {
  printf("FreqMul on ESP-IDF: 10 Hz tone -> method at its internal rate -> interpolator -> 48 kHz\n");
  float cpu;
  raw.begin(kInRate, 4000.f);
  report("raw              4 kHz x12: ", runMethod(kRaw, 4000.f, 10.f, 8.f, &cpu), 10.f, cpu);
  shifter.begin(kInRate, 4000.f);
  shifter.setShift(200.f);
  report("A shift +200 Hz  4 kHz x12: ", runMethod(kShift, 4000.f, 10.f, 8.f, &cpu), 210.f, cpu);
  grain.begin(kInRate, 16000.f);
  grain.setRatio(32.f);
  report("B grain x32     16 kHz x3:  ", runMethod(kGrain, 16000.f, 10.f, 8.f, &cpu), 320.f, cpu);
  bank.begin(kInRate, 8000.f);
  bank.setRatio(32.f);
  report("C bank x32       8 kHz x6:  ", runMethod(kBank, 8000.f, 10.f, 8.f, &cpu), 320.f, cpu);
  printf("done\n");
}
