// FreqMul: real-time moving of low-frequency signals (e.g. 1-40 Hz) into the audio range.
//
// Three methods, all fed at the input sample rate (push) and read at the audio rate (render):
//   FreqShifter   A. single-sideband frequency shift (allpass Hilbert pair): f -> f + shift (Hz spacing kept)
//   GrainScaler   B. delay-line pitch shifter with phase-synced splices: f -> k * f (ratios kept)
//   SpectralBank  C. STFT analysis + oscillator bank: f -> k * f per partial, optionally snapped to a scale
// plus FrontEnd (drift high-pass, mains notches, low-pass, slow AGC), Interpolator (integer-factor polyphase
// upsampler: render a method at a low internal rate, deliver at the DAC rate), Upsampler (plain band-limited
// pass-through for comparison) and TestSignal (a synthetic low-frequency test signal).
//
// Plain C++ (no Arduino.h, no third-party code), so it builds for the RA4M1 and on a desktop for tests.
// push() and render() must run in the same context (e.g. both from loop()); SpectralBank::update() runs
// the FFT and belongs in loop(), never in an interrupt.
#pragma once
#include <math.h>
#include <stdint.h>

#if defined(__GNUC__)
#define FM_INLINE inline __attribute__((always_inline))  // per-sample helpers: inline even at -Os
#else
#define FM_INLINE inline
#endif

namespace fm {

constexpr float kPi = 3.14159265358979f;
constexpr float kTwoPi = 6.28318530717959f;

// ---- shared helpers --------------------------------------------------------------------------------------

// sin(2*pi*phase/2^32) from a 1024-entry table (no interpolation: ~60 dB spurious-free, plenty here)
extern const float kSineTable[1024];
FM_INLINE float sine(uint32_t phase) { return kSineTable[phase >> 22]; }
FM_INLINE float cosine(uint32_t phase) { return kSineTable[(phase + 0x40000000u) >> 22]; }
// phase increment per sample for an oscillator at hz, sampled at rate (negative hz runs backwards)
inline uint32_t phaseInc(float hz, float rate) { return (uint32_t)(int64_t)(hz / rate * 4294967296.0); }

// smooth limiter: unity gain at 0, reaches +-1 with zero slope at +-1.5
inline float softClip(float x) {
  if (x > 1.5f) return 1.f;
  if (x < -1.5f) return -1.f;
  return x - (4.f / 27.f) * x * x * x;
}

// Catmull-Rom interpolation between y0 and y1 (t in [0,1))
FM_INLINE float cubic(float ym1, float y0, float y1, float y2, float t) {
  float c1 = 0.5f * (y1 - ym1);
  float c2 = ym1 - 2.5f * y0 + 2.f * y1 - 0.5f * y2;
  float c3 = 0.5f * (y2 - ym1) + 1.5f * (y0 - y1);
  return ((c3 * t + c2) * t + c1) * t + y0;
}

void fft(float* re, float* im, int n);  // in-place radix-2 complex FFT, n a power of two

struct Biquad {  // RBJ cookbook, transposed direct form II
  float b0 = 1, b1 = 0, b2 = 0, a1 = 0, a2 = 0, z1 = 0, z2 = 0;
  void lowpass(float fs, float f, float q);
  void highpass(float fs, float f, float q);
  void notch(float fs, float f, float q);
  void reset() { z1 = z2 = 0; }
  inline float process(float x) {
    float y = b0 * x + z1;
    z1 = b1 * x - a1 * y + z2;
    z2 = b2 * x - a2 * y;
    return y;
  }
 private:
  void set(float nb0, float nb1, float nb2, float a0, float na1, float na2);
};

// Follows the pushed input samples with a fractional read cursor that advances at the audio rate. The step
// is trimmed slowly (PI control on a ~1 s average of the lag, <= 2 %) to hold a fixed lag behind the newest
// sample, so the input converter and the audio timer may run from different clocks and input may arrive in
// bursts (e.g. DMA blocks): set lag >= burst length + 4 samples.
class RateBridge {
 public:
  void begin(float inRate, float outRate, float lag);  // lag in input samples
  void setLag(float lag) { lag_ = lag; }
  void reset();
  inline void pushed() { ++written; }                   // after storing input sample number `written`
  FM_INLINE void tick() {                               // advance one output sample
    ++ticks_;
    float f = frac + step;
    uint32_t i = idx;
    while (f >= 1.f) { f -= 1.f; ++i; }
    if (i + 2 < written) { idx = i; frac = f; } else { ++holds; }   // never read past the newest sample
  }
  void endBlock();          // once per render() block: lag control and resync
  float lagNow() const { return (float)(written - idx) - frac; }
  float stepNow() const { return step; }
  uint32_t idx = 0, written = 0, holds = 0;             // cursor = idx + frac (input-sample index)
  float frac = 0;
 private:
  float nominal_ = 0, step = 0, lag_ = 4, integ_ = 0, lagAvg_ = 0, tau_ = 16000;
  uint32_t ticks_ = 0;
  bool avgReady_ = false;
};

// ---- front end -------------------------------------------------------------------------------------------

class FrontEnd {
 public:
  // fs: input rate. hpHz removes drift, the notches remove mains and its 2nd harmonic, lpHz band-limits.
  void begin(float fs, float mainsHz = 50.f, float hpHz = 0.7f, float lpHz = 40.f);
  void setGain(float g) { gain_ = g; }                                 // input gain before the filters
  void setAgc(bool on, float targetRms = 0.25f, float seconds = 8.f);  // slow level normalisation
  float process(float x);                                              // one input sample in, one out
  float agcGain() const { return agcGain_; }
 private:
  Biquad hp_, n1_, n2_, lp1_, lp2_;
  bool useN2_ = false, agc_ = false;
  float fs_ = 250, gain_ = 1, agcGain_ = 1, ms_ = 0, agcCoef_ = 0, target_ = 0.25f, agcSeconds_ = 8.f;
  uint32_t agcN_ = 0;
};

// ---- A: single-sideband frequency shift ----------------------------------------------------------------------
// The band-limited input goes through a pair of allpass chains whose outputs are 90 deg apart (Olli Niemitalo's
// 8th-order IIR Hilbert pair: +-0.7 deg from 0.0008 to 0.5 x the sample rate, i.e. >= 44 dB image rejection
// over 0.2-120 Hz at 250 Hz). The analytic pair is interpolated to the audio rate and multiplied by a complex
// oscillator: every component f moves to f + shift (single sideband, Hz spacing kept).
class FreqShifter {
 public:
  // lag: how far (input samples) the read cursor trails the newest pushed sample; raise it when audio is
  // rendered ahead of real time into a FIFO (lag >= 4 + FIFO length / (audio rate / input rate)).
  void begin(float inRate, float audioRate, float bandLo = 1.f, float bandHi = 40.f, float lag = 4.f);
  void setShift(float hz);           // output = input + hz
  float shift() const { return shift_; }
  void push(float x);                // one input sample (at inRate)
  void render(float* out, int n);    // n audio samples (at audioRate)
  void reset();
  RateBridge bridge;
 private:
  struct Allpass {  // H(z) = (a^2 - z^-2) / (1 - a^2 z^-2)
    float a2 = 0, x1 = 0, x2 = 0, y1 = 0, y2 = 0;
    inline float process(float x) {
      float y = a2 * (x + y2) - x2;
      x2 = x1; x1 = x; y2 = y1; y1 = y;
      return y;
    }
  };
  static const int kRing = 64;
  float I_[kRing] = {}, Q_[kRing] = {};
  Biquad hp_, lp1_, lp2_;
  Allpass chainI_[4], chainQ_[4];
  float qDelay_ = 0;
  float inRate_ = 250, audioRate_ = 16000, shift_ = 200.f;
  uint32_t phOut_ = 0, incOut_ = 0;
};

// ---- B: delay-line (granular) pitch shifter ----------------------------------------------------------------
// Read heads sweep the last `window` seconds of input k times faster than real time and cross-fade (Hann
// windows on evenly spaced phases): every component f moves to k * f. Latency ~ window. With sync on (SOLA),
// each head's jump back is moved by up to one period of minSyncHz to where it best matches the head that is
// currently loudest, so splices land in phase: tones keep their exact pitch and bursts don't smear.
class GrainScaler {
 public:
  void begin(float inRate, float audioRate, float lag = 4.f);
  void setRatio(float k);             // 0.25 .. 64 (k = 32: five octaves up)
  void setWindow(float seconds);      // grain length (clamped to the ring buffer)
  void setHeads(int n);               // 2 .. 4 (more heads smooth the grain flutter)
  void setSync(bool on, float minSyncHz = 4.f);  // phase-matched splices, search range = 1 / minSyncHz
  float ratio() const { return k_; }
  void push(float x);                 // one input sample (at inRate)
  void render(float* out, int n);     // n audio samples (at audioRate)
  void reset();
  RateBridge bridge;
 private:
  static const int kRing = 512;       // 2 s at a 250 Hz input rate
  static const int kMaxHeads = 4;
  float buf_[kRing] = {};
  float delay_[kMaxHeads] = {}, ph_[kMaxHeads] = {};
  int heads_ = 2, search_ = 62, corr_ = 12;
  bool sync_ = true;
  float inRate_ = 250, audioRate_ = 16000, k_ = 32, win_ = 250, dmin_ = 3, winReq_ = 1.f;
  void spreadHeads();
  float splice(int h, float d, uint32_t i0);
};

// ---- C: spectral oscillator bank ---------------------------------------------------------------------------
// Every kHop input samples a kN-point Hann-windowed FFT of the newest samples is taken. Spectral peaks in
// [lo, hi] are picked (local maxima above -45 dB of the frame's strongest; sidelobe duplicates merged), each
// with its frequency from the phase advance between frames (phase vocoder) and its amplitude from parabolic
// interpolation. Peaks are matched to the oscillators of the previous frame (nearest frequency), so each
// component is one phase-continuous sine at k * f (or the nearest note of a scale); amplitudes ramp over one
// hop, unmatched oscillators fade out.
class SpectralBank {
 public:
  enum Scale { kOff = 0, kChromatic, kMajor, kMinor, kPentatonicMajor, kPentatonicMinor, kScaleCount };
  void begin(float inRate, float audioRate, float bandLo = 1.f, float bandHi = 40.f);
  void setRatio(float k) { k_ = k; }
  void setQuantize(Scale s, float rootHz = 110.f) { scale_ = s; root_ = rootHz; }
  void setTilt(float t) { tilt_ = t; }          // amplitude *= (f / 10 Hz)^t, flattens a 1/f spectrum
  void setCompress(float g) { gamma_ = g; }     // amplitude = magnitude^g
  void setRelease(float s) { release_ = s; }    // loudness memory: bursts stand out against this many seconds
  float ratio() const { return k_; }
  void push(float x);                           // one input sample (at inRate)
  bool update();                                // FFT when a hop is complete; call from loop()
  void render(float* out, int n);               // n audio samples (at audioRate)
  void reset();
  int partials() const;                         // oscillators currently sounding
  float quantize(float hz) const;               // nearest note of the current scale
  static const char* scaleName(int s);
 private:
  static const int kN = 512, kHop = 64, kMaxBins = 96, kMaxPartials = 40;
  float ring_[kN] = {};
  float re_[kN], im_[kN];
  float prevPhase_[kMaxBins + 2] = {};
  struct Partial { uint32_t phase, inc; float amp, target, dAmp, fIn; bool live; };
  Partial p_[kMaxPartials] = {};
  uint32_t written_ = 0, analysed_ = 0;
  int binLo_ = 3, binHi_ = 81, rampLeft_ = 0, rampLen_ = 4096;
  float inRate_ = 250, audioRate_ = 16000, k_ = 32, tilt_ = 0.5f, gamma_ = 0.7f, level_ = 0, release_ = 10.f;
  Scale scale_ = kOff;
  float root_ = 110.f;
};

// ---- interpolator ------------------------------------------------------------------------------------------
// Integer-factor polyphase FIR upsampler (Kaiser-windowed sinc, each phase normalised to unity DC gain, so CV
// passes exactly). Keeps content below pass * inRate and removes its images above (1 - pass) * inRate; with
// 16 taps per phase and pass = 0.4 the images are ~55 dB down. Cost: `taps` multiply-adds per output sample;
// delay: taps / 2 input samples. Streaming: any block size.
class Interpolator {
 public:
  void begin(int factor, int taps = 16, float pass = 0.4f);  // factor 1..16, taps 4..24 per phase
  int factor() const { return L_; }
  void process(const float* in, int n, float* out);           // n input samples -> n * factor() outputs
  void reset();
 private:
  static const int kMaxFactor = 16, kMaxTaps = 24;
  float h_[kMaxFactor * kMaxTaps] = {};  // phase-major: h_[p * T_ + t] multiplies the input t samples back
  float hist_[2 * kMaxTaps] = {};        // newest-first window at hist_[pos_ .. pos_ + T_ - 1]
  int L_ = 1, T_ = 16, pos_ = 0;
};

// ---- pass-through ------------------------------------------------------------------------------------------
class Upsampler {  // the input itself at the audio rate (cubic interpolation), for scopes and A/B tests
 public:
  void begin(float inRate, float audioRate, float lag = 4.f);
  void push(float x);
  void render(float* out, int n);
  RateBridge bridge;
 private:
  static const int kRing = 64;
  float buf_[kRing] = {};
};

// ---- test signal -------------------------------------------------------------------------------------------
// Synthetic low-frequency test signal (arbitrary units): pink background, a burst tone that swells for 4 s
// every 8 s, two steady tones (6 and 20 Hz), mains hum, slow drift and a transient pulse every ~7 s.
class TestSignal {
 public:
  void begin(float fs, uint32_t seed = 1);
  // burst tone: frequency, amplitude while swelled (high) and between swells (low)
  void setBurst(float hz, float high = 30.f, float low = 5.f) { burstHz_ = hz; burstHigh_ = high; burstLow_ = low; }
  void setMains(float hz, float level = 20.f) { mainsHz_ = hz; mainsLevel_ = level; }
  float next();                       // next sample (at fs)
 private:
  float fs_ = 250, burstHz_ = 10, burstHigh_ = 30, burstLow_ = 5, mainsHz_ = 50, mainsLevel_ = 20;
  float b0_ = 0, b1_ = 0, b2_ = 0;
  uint32_t n_ = 0, pulseAt_ = 0, rng_ = 1;
  uint32_t phBurst_ = 0, phLow_ = 0, phHigh_ = 0, phMains_ = 0, phDrift_ = 0, phWobble_ = 0;
  float white();
};

// ---- measurement -------------------------------------------------------------------------------------------
// Frequency of a tone from interpolated positive-going zero crossings (streaming, no buffer).
class ZeroCrossHz {
 public:
  void begin(float rate, float hysteresis = 0.01f) { rate_ = rate; hyst_ = hysteresis; reset(); }
  void reset() { n_ = 0; first_ = last_ = -1; count_ = 0; prev_ = 0; armed_ = false; }
  void add(float x);
  float hz() const { return count_ > 1 ? (float)((count_ - 1) * rate_ / (last_ - first_)) : 0.f; }
 private:
  float rate_ = 16000, prev_ = 0, hyst_ = 0.01f;
  bool armed_ = false;
  double first_ = -1, last_ = -1;
  int64_t n_ = 0;
  int count_ = 0;
};

}  // namespace fm
