// BasicRender: a 10 Hz tone through each FreqMul method, each at its own internal rate, interpolated to 48 kHz.
// Prints the measured output frequency (PASS/FAIL) and the CPU share of real time. No DAC needed: it runs on any
// supported board (UNO R4, ESP32-S3, ...) and doubles as a benchmark. Serial at 115200.
#include <FreqMul.h>

namespace {

constexpr float kInRate = 250.f;     // input sample rate
constexpr float kOutRate = 48000.f;  // DAC rate
constexpr int kOutBlock = 192;       // output samples per input sample (48000 / 250)

fm::Upsampler raw;
fm::FreqShifter shifter;
fm::GrainScaler grain;
fm::SpectralBank bank;
fm::Interpolator up;

enum Method { kRaw, kShift, kGrain, kBank };

// Feeds `seconds` of a tone, renders each method at `rate`, interpolates to 48 kHz. Returns the output frequency
// (measured after 4 s of settling) and the CPU share of real time spent on push + render + interpolate.
float runMethod(Method m, float rate, float tone, float seconds, float* cpuPct) {
  const int factor = (int)(kOutRate / rate);
  const int n = kOutBlock / factor;  // internal samples per input sample
  up.begin(factor);
  fm::ZeroCrossHz zc;
  zc.begin(kOutRate);
  float in[64], out[kOutBlock];
  uint32_t busy = 0;
  const int total = (int)(seconds * kInRate), settle = (int)(4 * kInRate);
  for (int i = 0; i < total; ++i) {
    float x = 0.25f * sinf(fm::kTwoPi * tone * i / kInRate);
    uint32_t t0 = micros();
    switch (m) {
      case kRaw: raw.push(x); raw.render(in, n); break;
      case kShift: shifter.push(x); shifter.render(in, n); break;
      case kGrain: grain.push(x); grain.render(in, n); break;
      case kBank: bank.push(x); bank.update(); bank.render(in, n); break;
    }
    up.process(in, n, out);
    uint32_t dt = micros() - t0;
    if (i >= settle) {
      busy += dt;
      for (int k = 0; k < kOutBlock; ++k) zc.add(out[k]);
    }
  }
  *cpuPct = 100.f * busy / (1e6f * (total - settle) / kInRate);
  return zc.hz();
}

void report(const char* name, float hz, float want, float cpu) {
  bool ok = fabsf(hz - want) <= want * 0.006f;
  Serial.print(name);
  Serial.print(hz, 2);
  Serial.print(F(" Hz (want "));
  Serial.print(want, 2);
  Serial.print(ok ? F(") PASS   cpu ") : F(") FAIL   cpu "));
  Serial.print(cpu, 1);
  Serial.println(F(" %"));
}

}  // namespace

void setup() {
  Serial.begin(115200);
  delay(2000);
  Serial.println(F("FreqMul BasicRender: 10 Hz tone -> method at its internal rate -> interpolator -> 48 kHz"));
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

  bank.begin(kInRate, 8000.f);
  bank.setRatio(32.f);
  bank.setQuantize(fm::SpectralBank::kMajor, 110.f);
  report("C bank x32 A maj 8 kHz x6:  ", runMethod(kBank, 8000.f, 10.f, 8.f, &cpu), 329.63f, cpu);
  Serial.println(F("done"));
}

void loop() {}
