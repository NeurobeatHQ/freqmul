// Desktop test for FreqMul (not compiled by the Arduino IDE: it only builds the sketch root and src/).
//   c++ -O2 -std=c++17 -I../../src host_test.cpp ../../src/FreqMul.cpp -o host_test && ./host_test
// Checks every method's output frequency against the expected mapping, two-tone ratio/spacing preservation
// and clock-drift tolerance, then writes demo WAVs (TestSignal -> each method) into out/ next to the binary.
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>

#include <algorithm>
#include <functional>
#include <string>
#include <vector>

#include "FreqMul.h"

static const float IN_RATE = 250.f, AUDIO = 16000.f;
static const int BLOCK = 64;  // = AUDIO / IN_RATE: one input sample per audio block
static int failures = 0;

using Push = std::function<void(float)>;
using Render = std::function<void(float*, int)>;
using Update = std::function<void()>;

// feed `seconds` of input (src(t)) and collect the audio output
static std::vector<float> run(Push push, Render render, Update update, std::function<float(float)> src, float seconds,
                              float inRate = IN_RATE, int block = BLOCK) {
  // integer clock: input sample i is due before audio sample i * AUDIO / inRate (no floating-point ties)
  std::vector<float> out;
  float buf[256];
  const int64_t rateIn = (int64_t)inRate, rateOut = (int64_t)AUDIO;
  int64_t i = 0;
  const int blocks = (int)(seconds * AUDIO / block);
  for (int b = 0; b < blocks; ++b) {
    int64_t a = (int64_t)b * block;
    while (i * rateOut <= a * rateIn) { push(src((float)((double)i / inRate))); ++i; }
    update();
    render(buf, block);
    out.insert(out.end(), buf, buf + block);
  }
  return out;
}

// strongest spectral peaks (Hann window, parabolic interpolation) of out[from:]
static std::vector<std::pair<float, float>> peaks(const std::vector<float>& x, size_t from, int npk) {
  int n = 1;
  while (n * 2 <= (int)(x.size() - from) && n < (1 << 18)) n *= 2;
  std::vector<float> re(n), im(n, 0.f);
  for (int i = 0; i < n; ++i) re[i] = x[from + i] * (0.5f - 0.5f * cosf(fm::kTwoPi * i / n));
  fm::fft(re.data(), im.data(), n);
  std::vector<float> mag(n / 2);
  for (int i = 0; i < n / 2; ++i) mag[i] = sqrtf(re[i] * re[i] + im[i] * im[i]);
  std::vector<std::pair<float, float>> res;
  std::vector<char> used(n / 2, 0);
  for (int k = 0; k < npk; ++k) {
    int best = 2;
    for (int i = 2; i < n / 2 - 1; ++i)
      if (!used[i] && mag[i] > mag[best]) best = i;
    float a = mag[best - 1], b = mag[best], c = mag[best + 1];
    float off = 0.5f * (a - c) / (a - 2 * b + c);
    res.push_back({(best + off) * AUDIO / n, b});
    int guard = (int)(1.5 * n / AUDIO) + 2;  // mask +-1.5 Hz around the peak
    for (int i = best - guard; i <= best + guard; ++i)
      if (i >= 0 && i < n / 2) used[i] = 1;
  }
  return res;
}

static void check(const char* name, float got, float want, float tolPct) {
  bool ok = fabsf(got - want) <= want * tolPct / 100.f;
  if (!ok) ++failures;
  printf("  %-44s got %8.2f Hz  want %8.2f Hz  %s\n", name, got, want, ok ? "PASS" : "FAIL");
}

static float zeroCross(const std::vector<float>& x, size_t from) {
  fm::ZeroCrossHz z;
  z.begin(AUDIO);
  for (size_t i = from; i < x.size(); ++i) z.add(x[i]);
  return z.hz();
}

static void writeWav(const std::string& path, const std::vector<float>& x, float gain) {
  FILE* f = fopen(path.c_str(), "wb");
  if (!f) return;
  uint32_t n = (uint32_t)x.size(), rate = (uint32_t)AUDIO, bytes = n * 2;
  auto u32 = [&](uint32_t v) { fwrite(&v, 4, 1, f); };
  auto u16 = [&](uint16_t v) { fwrite(&v, 2, 1, f); };
  fwrite("RIFF", 1, 4, f); u32(36 + bytes); fwrite("WAVEfmt ", 1, 8, f);
  u32(16); u16(1); u16(1); u32(rate); u32(rate * 2); u16(2); u16(16);
  fwrite("data", 1, 4, f); u32(bytes);
  for (float v : x) { int16_t s = (int16_t)(fm::softClip(v * gain) * 32000.f); fwrite(&s, 2, 1, f); }
  fclose(f);
}

int main(int, char** argv) {
  const float kTone = 10.f;
  auto tone = [&](float t) { return 0.25f * sinf(fm::kTwoPi * kTone * t); };
  auto twoTone = [&](float t) { return 0.2f * sinf(fm::kTwoPi * 6.f * t) + 0.2f * sinf(fm::kTwoPi * 10.f * t); };
  const size_t settle = (size_t)(4 * AUDIO);  // ignore the first 4 s (filters, FFT history, lag control)

  printf("single tone %.1f Hz (zero-crossing frequency of the output):\n", kTone);
  {
    fm::Upsampler u; u.begin(IN_RATE, AUDIO);
    auto o = run([&](float x) { u.push(x); }, [&](float* b, int n) { u.render(b, n); }, [] {}, tone, 8);
    check("pass-through", zeroCross(o, settle), 10.f, 1.f);
  }
  {
    fm::FreqShifter s; s.begin(IN_RATE, AUDIO); s.setShift(200.f);
    auto o = run([&](float x) { s.push(x); }, [&](float* b, int n) { s.render(b, n); }, [] {}, tone, 20);
    check("A shift +200 Hz", zeroCross(o, settle), 210.f, 0.5f);
    auto p = peaks(o, settle, 4);
    float main = p[0].second, worst = 0, at = 0;
    for (size_t i = 1; i < p.size(); ++i) if (p[i].second > worst) { worst = p[i].second; at = p[i].first; }
    float db = 20.f * log10f(worst / main + 1e-12f);
    bool ok = db < -35.f;
    if (!ok) ++failures;
    printf("  %-44s %6.1f dB at %.1f Hz (mirror would be 190 Hz)  %s\n", "A shift: strongest spur vs 210 Hz", db, at, ok ? "PASS" : "FAIL");
  }
  for (float k : {8.f, 32.f}) {
    fm::GrainScaler g; g.begin(IN_RATE, AUDIO); g.setRatio(k);
    auto o = run([&](float x) { g.push(x); }, [&](float* b, int n) { g.render(b, n); }, [] {}, tone, 8);
    check(("B grain x" + std::to_string((int)k)).c_str(), zeroCross(o, settle), kTone * k, 1.f);
  }
  {
    // 10.3 Hz: a 1 s grain holds 10.3 cycles, so plain splices jump 0.3 cycle each time and pull the pitch
    auto tone103 = [](float t) { return 0.25f * sinf(fm::kTwoPi * 10.3f * t); };
    for (bool sync : {false, true}) {
      fm::GrainScaler g; g.begin(IN_RATE, AUDIO); g.setRatio(32); g.setSync(sync);
      auto o = run([&](float x) { g.push(x); }, [&](float* b, int n) { g.render(b, n); }, [] {}, tone103, 10);
      float hz = zeroCross(o, settle);
      if (sync) check("B grain x32, 10.3 Hz, sync (SOLA) on", hz, 329.6f, 0.5f);
      else printf("  %-44s got %8.2f Hz  (plain splices; for reference, not checked)\n", "B grain x32, 10.3 Hz, sync off", hz);
    }
  }
  for (float k : {8.f, 32.f}) {
    fm::SpectralBank c; c.begin(IN_RATE, AUDIO); c.setRatio(k);
    auto o = run([&](float x) { c.push(x); }, [&](float* b, int n) { c.render(b, n); }, [&] { c.update(); }, tone, 8);
    check(("C bank x" + std::to_string((int)k)).c_str(), zeroCross(o, settle), kTone * k, 1.f);
  }
  {
    fm::SpectralBank c; c.begin(IN_RATE, AUDIO); c.setRatio(32); c.setQuantize(fm::SpectralBank::kMajor, 110.f);
    auto o = run([&](float x) { c.push(x); }, [&](float* b, int n) { c.render(b, n); }, [&] { c.update(); }, tone, 8);
    check("C bank x32, A major (320 -> E4 329.63)", zeroCross(o, settle), 329.63f, 0.5f);
  }

  printf("two tones 6 + 10 Hz (two strongest FFT peaks):\n");
  {
    fm::FreqShifter s; s.begin(IN_RATE, AUDIO); s.setShift(200.f);
    auto o = run([&](float x) { s.push(x); }, [&](float* b, int n) { s.render(b, n); }, [] {}, twoTone, 20);
    auto p = peaks(o, settle, 2);
    float lo = fminf(p[0].first, p[1].first), hi = fmaxf(p[0].first, p[1].first);
    check("A shift: lower", lo, 206.f, 0.5f);
    check("A shift: upper (spacing kept: 4 Hz)", hi, 210.f, 0.5f);
  }
  {
    fm::GrainScaler g; g.begin(IN_RATE, AUDIO); g.setRatio(32);
    auto o = run([&](float x) { g.push(x); }, [&](float* b, int n) { g.render(b, n); }, [] {}, twoTone, 20);
    auto p = peaks(o, settle, 2);
    float lo = fminf(p[0].first, p[1].first), hi = fmaxf(p[0].first, p[1].first);
    check("B grain x32: lower", lo, 192.f, 1.f);
    check("B grain x32: upper (ratio kept: 5/3)", hi, 320.f, 1.f);
  }
  {
    fm::SpectralBank c; c.begin(IN_RATE, AUDIO); c.setRatio(32);
    auto o = run([&](float x) { c.push(x); }, [&](float* b, int n) { c.render(b, n); }, [&] { c.update(); }, twoTone, 20);
    auto p = peaks(o, settle, 2);
    float lo = fminf(p[0].first, p[1].first), hi = fmaxf(p[0].first, p[1].first);
    check("C bank x32: lower", lo, 192.f, 1.f);
    check("C bank x32: upper", hi, 320.f, 1.f);
  }

  printf("B at shorter windows (share of output energy within +-2 Hz of the targets, x32, 2 heads, sync):\n");
  for (float win : {1.f, 0.5f, 0.25f}) {
    for (int mix = 0; mix < 3; ++mix) {
      std::vector<float> tones = mix == 0 ? std::vector<float>{10.f} : mix == 1 ? std::vector<float>{6.f} : std::vector<float>{6.f, 10.f};
      fm::GrainScaler g; g.begin(IN_RATE, AUDIO); g.setRatio(32); g.setWindow(win);
      auto o = run([&](float x) { g.push(x); }, [&](float* b, int n) { g.render(b, n); }, [] {},
                   [&](float t) { float x = 0; for (float f : tones) x += 0.2f * sinf(fm::kTwoPi * f * t); return x; }, 14);
      int n = 1 << 17;
      std::vector<float> re(n), im(n, 0.f);
      for (int i = 0; i < n; ++i) re[i] = o[settle + i] * (0.5f - 0.5f * cosf(fm::kTwoPi * i / n));
      fm::fft(re.data(), im.data(), n);
      double tot = 0, on = 0;
      for (int i = 1; i < n / 2; ++i) {
        double e = (double)re[i] * re[i] + (double)im[i] * im[i];
        tot += e;
        for (float f : tones) if (fabsf(i * AUDIO / n - 32.f * f) <= 2.f) on += e;
      }
      float pct = (float)(100.0 * on / tot);
      bool ok = pct >= 95.f;
      if (!ok) ++failures;
      char nm[64];
      snprintf(nm, sizeof nm, "B x32 window %.2f s, %s", win, mix == 0 ? "10 Hz" : mix == 1 ? "6 Hz" : "6+10 Hz mix");
      printf("  %-44s %5.1f %% on target  %s\n", nm, pct, ok ? "PASS" : "FAIL");
    }
  }

  printf("C peak picking: one oscillator per component (strongest component >30 Hz from k*f):\n");
  for (auto kt : std::vector<std::pair<float, float>>{{2, 30}, {3, 30}, {3, 31}, {8, 10}, {32, 10}, {32, 23}, {64, 10}}) {
    float k = kt.first, tn = kt.second;
    fm::SpectralBank c; c.begin(IN_RATE, AUDIO); c.setRatio(k);
    auto o = run([&](float x) { c.push(x); }, [&](float* b, int n) { c.render(b, n); }, [&] { c.update(); },
                 [&](float t) { return 0.25f * sinf(fm::kTwoPi * tn * t); }, 14);
    std::string nm = "C bank x" + std::to_string((int)k) + ", " + std::to_string((int)tn) + " Hz";
    check(nm.c_str(), zeroCross(o, settle), k * tn, 0.5f);
    auto p = peaks(o, settle, 6);
    float main = p[0].second, worst = 0, at = 0;
    for (auto& q : p) if (fabsf(q.first - k * tn) > 30.f && q.second > worst) { worst = q.second; at = q.first; }
    float db = worst > 0 ? 20.f * log10f(worst / main) : -200.f;
    bool ok = db < -40.f;
    if (!ok) ++failures;
    printf("  %-44s %6.1f dB at %.1f Hz  %s\n", "   strongest other component", db, at, ok ? "PASS" : "FAIL");
  }

  printf("Interpolator (x12, 16 taps/phase, pass 0.4):\n");
  {
    const float fin = 4000.f, fout = 48000.f;
    auto interp = [&](std::function<float(int)> x, int n, int blk) {
      fm::Interpolator ip; ip.begin(12);
      std::vector<float> in(n), out(n * 12);
      for (int i = 0; i < n; ++i) in[i] = x(i);
      for (int i = 0; i < n; i += blk) ip.process(in.data() + i, std::min(blk, n - i), out.data() + i * 12);
      return out;
    };
    auto dc = interp([](int) { return 0.7f; }, 400, 400);
    float lo = 1e9, hi = -1e9;
    for (size_t i = 200; i < dc.size(); ++i) { lo = fminf(lo, dc[i]); hi = fmaxf(hi, dc[i]); }
    bool ok = fabsf(lo - 0.7f) < 1e-5f && fabsf(hi - 0.7f) < 1e-5f;
    if (!ok) ++failures;
    printf("  %-44s %.6f .. %.6f  %s\n", "DC 0.7 in -> out (unity per phase)", lo, hi, ok ? "PASS" : "FAIL");
    for (float ft : {600.f, 1200.f, 1600.f}) {  // passband up to 0.4 * 4 kHz
      auto o = interp([&](int i) { return 0.5f * sinf(fm::kTwoPi * ft * i / fin); }, 16000, 400);
      auto pk = [&](const std::vector<float>& x) {
        int n = 1 << 17; std::vector<float> re(n), im(n, 0.f);
        for (int i = 0; i < n; ++i) re[i] = x[2000 + i] * (0.5f - 0.5f * cosf(fm::kTwoPi * i / n));
        fm::fft(re.data(), im.data(), n);
        std::vector<float> m(n / 2); for (int i = 0; i < n / 2; ++i) m[i] = hypotf(re[i], im[i]);
        return m;
      };
      auto m = pk(o);
      int n = 1 << 17;
      auto at = [&](float f) { int b = (int)(f * n / fout + 0.5f); float v = 0; for (int i = b - 3; i <= b + 3; ++i) v = fmaxf(v, m[i]); return v; };
      float main = at(ft), img = 0, imgF = 0;
      for (int r = 1; r <= 11; ++r) for (float f : {r * fin - ft, r * fin + ft}) if (at(f) > img) { img = at(f); imgF = f; }
      double sc = 0, ss = 0; int cnt = 0;  // exact amplitude: projection onto the tone (FFT bins scallop)
      for (size_t q = 2000; q < o.size(); ++q) { double ph = 2 * M_PI * ft * (double)q / fout; sc += o[q] * cos(ph); ss += o[q] * sin(ph); ++cnt; }
      float gain = 20.f * log10f((float)(2.0 * sqrt(sc * sc + ss * ss) / cnt) / 0.5f);
      float idb = 20.f * log10f(img / main);
      bool okg = fabsf(gain) < 0.5f, oki = idb < -50.f;
      if (!okg || !oki) ++failures;
      printf("  tone %4.0f Hz: passband gain %+.2f dB %s, strongest image %.0f Hz at %.1f dB %s\n", ft, gain, okg ? "PASS" : "FAIL", imgF, idb, oki ? "PASS" : "FAIL");
    }
    auto a = interp([&](int i) { return sinf(0.3f * i) + 0.2f * sinf(1.1f * i); }, 1000, 1000);
    auto b = interp([&](int i) { return sinf(0.3f * i) + 0.2f * sinf(1.1f * i); }, 1000, 7);
    double d = 0; for (size_t i = 0; i < a.size(); ++i) d = fmax(d, fabs(a[i] - b[i]));
    ok = d == 0.0;
    if (!ok) ++failures;
    printf("  %-44s max diff %.3g  %s\n", "streaming: 7-sample vs 1000-sample blocks", d, ok ? "PASS" : "FAIL");
  }
  printf("C at 8 kHz internal rate, interpolated x6 to 48 kHz:\n");
  {
    fm::SpectralBank c; c.begin(IN_RATE, 8000.f); c.setRatio(32);
    fm::Interpolator ip; ip.begin(6);
    std::vector<float> o; float b8[32], b48[192];
    int64_t i = 0;
    for (int blk = 0; blk < 12 * 250; ++blk) {   // per input sample: 32 internal samples -> 192 output samples
      c.push(0.25f * sinf(fm::kTwoPi * 10.f * (float)i / IN_RATE)); ++i;
      c.update(); c.render(b8, 32); ip.process(b8, 32, b48);
      if (blk >= 4 * 250) o.insert(o.end(), b48, b48 + 192);
    }
    fm::ZeroCrossHz z; z.begin(48000.f); for (float v : o) z.add(v);
    check("C x32 via 8 kHz -> 48 kHz", z.hz(), 320.f, 0.3f);
  }

  printf("sketch chain: each method at its internal rate -> Interpolator -> 48 kHz (lag 12):\n");
  {
    struct M { const char* name; float rate; int factor; float tone, want; };
    const M ms[] = {{"raw 4 kHz x12", 4000, 12, 10, 10}, {"A shift +200, 4 kHz x12", 4000, 12, 10, 210},
                    {"B grain x32, 16 kHz x3", 16000, 3, 10, 320}, {"C bank x32, 8 kHz x6", 8000, 6, 10, 320}};
    for (int m = 0; m < 4; ++m) {
      fm::Upsampler u; u.begin(IN_RATE, ms[m].rate, 12);
      fm::FreqShifter a; a.begin(IN_RATE, ms[m].rate, 1, 40, 12); a.setShift(200);
      fm::GrainScaler g; g.begin(IN_RATE, ms[m].rate, 12); g.setRatio(32);
      fm::SpectralBank c; c.begin(IN_RATE, ms[m].rate); c.setRatio(32);
      fm::Interpolator ip; ip.begin(ms[m].factor);
      const int n = 192 / ms[m].factor;
      float in[64], out[192];
      std::vector<float> o;
      for (int i = 0; i < 14 * 250; ++i) {
        float x = 0.25f * sinf(fm::kTwoPi * ms[m].tone * i / IN_RATE);
        if (m == 0) { u.push(x); u.render(in, n); }
        if (m == 1) { a.push(x); a.render(in, n); }
        if (m == 2) { g.push(x); g.render(in, n); }
        if (m == 3) { c.push(x); c.update(); c.render(in, n); }
        ip.process(in, n, out);
        if (i >= 4 * 250) o.insert(o.end(), out, out + 192);
      }
      fm::ZeroCrossHz z; z.begin(48000.f); for (float v : o) z.add(v);
      check(ms[m].name, z.hz(), ms[m].want, 0.5f);
    }
  }

  printf("small render blocks (the per-call work must not change the result):\n");
  for (int blk : {4, 7, 16}) {
    std::string b = " (" + std::to_string(blk) + "-sample blocks)";
    {
      fm::FreqShifter x; x.begin(IN_RATE, AUDIO); x.setShift(200.f);
      auto o = run([&](float v) { x.push(v); }, [&](float* q, int n) { x.render(q, n); }, [] {}, tone, 8, IN_RATE, blk);
      check(("A shift +200 Hz" + b).c_str(), zeroCross(o, settle), 210.f, 0.5f);
    }
    {
      fm::GrainScaler x; x.begin(IN_RATE, AUDIO); x.setRatio(32);
      auto o = run([&](float v) { x.push(v); }, [&](float* q, int n) { x.render(q, n); }, [] {}, tone, 8, IN_RATE, blk);
      check(("B grain x32" + b).c_str(), zeroCross(o, settle), 320.f, 1.f);
    }
    {
      fm::SpectralBank x; x.begin(IN_RATE, AUDIO); x.setRatio(32);
      auto o = run([&](float v) { x.push(v); }, [&](float* q, int n) { x.render(q, n); }, [&] { x.update(); }, tone, 8, IN_RATE, blk);
      check(("C bank x32" + b).c_str(), zeroCross(o, settle), 320.f, 1.f);
    }
  }
  {  // C is deterministic: any block size must give the same samples as 64-sample blocks
    auto bankRun = [&](int blk) {
      fm::SpectralBank x; x.begin(IN_RATE, AUDIO); x.setRatio(32);
      return run([&](float v) { x.push(v); }, [&](float* q, int n) { x.render(q, n); }, [&] { x.update(); }, twoTone, 6, IN_RATE, blk);
    };
    auto a = bankRun(64), b = bankRun(4);
    double diff = 0, ref = 0;
    for (size_t i = 0; i < a.size() && i < b.size(); ++i) { diff += (a[i] - b[i]) * (a[i] - b[i]); ref += a[i] * a[i]; }
    float db = 10.f * log10f((float)(diff / ref) + 1e-30f);
    bool ok = db < -100.f;
    if (!ok) ++failures;
    printf("  %-44s %7.1f dB  %s\n", "C bank: 4-sample vs 64-sample blocks, difference", db, ok ? "PASS" : "FAIL");
  }

  printf("clock drift: input converter 0.4 %% fast (251 Hz vs 250 nominal), output must follow the true input:\n");
  {
    fm::GrainScaler g; g.begin(IN_RATE, AUDIO); g.setRatio(32);
    auto o = run([&](float x) { g.push(x); }, [&](float* b, int n) { g.render(b, n); }, [] {}, tone, 60, 251.f);
    // a real 10 Hz tone arrives as 251 samples/s; a cursor stuck at the nominal 250/s would play it at
    // 10 * 250/251 Hz and keep falling behind. The bridge follows the real sample clock: 32 * 10 = 320 Hz.
    check("B grain x32 with 251 Hz input clock", zeroCross(o, (size_t)(40 * AUDIO)), 320.f, 0.3f);
    printf("  bridge lag %.2f samples (target 4), holds %u, step %.5f (nominal %.5f)\n", g.bridge.lagNow(),
           g.bridge.holds, g.bridge.stepNow(), IN_RATE / AUDIO);
    bool ok = fabsf(g.bridge.lagNow() - 4.f) < 2.f;
    if (!ok) ++failures;
    printf("  lag held near target: %s\n", ok ? "PASS" : "FAIL");
  }

  // demo WAVs: TestSignal through the front end and each method, written next to the binary
  const char* slash = strrchr(argv[0], '/');
  const std::string outDir = (slash ? std::string(argv[0], slash - argv[0] + 1) : std::string()) + "out/";
  mkdir(outDir.c_str(), 0755);
  printf("demo WAVs (24 s of TestSignal: a tone swells 4 s out of every 8) -> %s\n", outDir.c_str());
  auto demo = [&](const char* name, Push push, Render render, Update update, float gain) {
    fm::TestSignal sig; sig.begin(IN_RATE, 7);
    fm::FrontEnd fe; fe.begin(IN_RATE, 50.f); fe.setAgc(true, 0.25f, 8.f);
    auto o = run([&](float) { push(fe.process(sig.next())); }, render, update, [](float) { return 0.f; }, 24);
    writeWav(outDir + name, o, gain);
    printf("  %s%s\n", outDir.c_str(), name);
  };
  {
    fm::FreqShifter s; s.begin(IN_RATE, AUDIO); s.setShift(200.f);
    demo("A_shift_plus200Hz.wav", [&](float x) { s.push(x); }, [&](float* b, int n) { s.render(b, n); }, [] {}, 2.f);
  }
  {
    fm::GrainScaler g; g.begin(IN_RATE, AUDIO); g.setRatio(32); g.setHeads(3);
    demo("B_grain_x32.wav", [&](float x) { g.push(x); }, [&](float* b, int n) { g.render(b, n); }, [] {}, 2.f);
  }
  {
    fm::SpectralBank c; c.begin(IN_RATE, AUDIO); c.setRatio(32);
    demo("C_bank_x32.wav", [&](float x) { c.push(x); }, [&](float* b, int n) { c.render(b, n); }, [&] { c.update(); }, 1.5f);
  }
  {
    fm::SpectralBank c; c.begin(IN_RATE, AUDIO); c.setRatio(32); c.setQuantize(fm::SpectralBank::kPentatonicMinor, 110.f);
    demo("C_bank_x32_A_minor_pentatonic.wav", [&](float x) { c.push(x); }, [&](float* b, int n) { c.render(b, n); },
         [&] { c.update(); }, 1.5f);
  }

  printf(failures ? "\n%d FAILED\n" : "\nall passed\n", failures);
  return failures ? 1 : 0;
}
