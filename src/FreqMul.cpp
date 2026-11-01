// The Arduino core builds everything at -Os, which leaves the per-sample helpers as calls; the DSP here is
// worth the extra flash at -O2 (the sketch itself keeps -Os).
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC optimize("O2")
#endif

#include "FreqMul.h"

#include <string.h>

#include "SineTable.h"

namespace fm {

// ---- shared helpers --------------------------------------------------------------------------------------

static inline void initTable() {}

static FM_INLINE float hann01(float ph) {  // 0.5 - 0.5 cos(2 pi ph), ph in [0, 1)
  int i = (int)(ph * 1024.f);
  if (i < 0) i = 0;
  if (i > 1023) i = 1023;
  return 0.5f - 0.5f * kSineTable[(i + 256) & 1023];
}

void fft(float* re, float* im, int n) {
  for (int i = 1, j = 0; i < n; ++i) {  // bit reversal
    int bit = n >> 1;
    for (; j & bit; bit >>= 1) j ^= bit;
    j ^= bit;
    if (i < j) {
      float t = re[i]; re[i] = re[j]; re[j] = t;
      t = im[i]; im[i] = im[j]; im[j] = t;
    }
  }
  for (int len = 2; len <= n; len <<= 1) {
    float ang = -kTwoPi / len, wr = cosf(ang), wi = sinf(ang);
    int half = len >> 1;
    for (int i = 0; i < n; i += len) {
      float cr = 1.f, ci = 0.f;
      for (int k = 0; k < half; ++k) {
        int a = i + k, b = a + half;
        float tr = re[b] * cr - im[b] * ci, ti = re[b] * ci + im[b] * cr;
        re[b] = re[a] - tr; im[b] = im[a] - ti;
        re[a] += tr; im[a] += ti;
        float ncr = cr * wr - ci * wi;
        ci = cr * wi + ci * wr;
        cr = ncr;
      }
    }
  }
}

void Biquad::set(float nb0, float nb1, float nb2, float a0, float na1, float na2) {
  b0 = nb0 / a0; b1 = nb1 / a0; b2 = nb2 / a0; a1 = na1 / a0; a2 = na2 / a0;
  reset();
}
void Biquad::lowpass(float fs, float f, float q) {
  float w = kTwoPi * f / fs, c = cosf(w), al = sinf(w) / (2.f * q);
  set((1 - c) / 2, 1 - c, (1 - c) / 2, 1 + al, -2 * c, 1 - al);
}
void Biquad::highpass(float fs, float f, float q) {
  float w = kTwoPi * f / fs, c = cosf(w), al = sinf(w) / (2.f * q);
  set((1 + c) / 2, -(1 + c), (1 + c) / 2, 1 + al, -2 * c, 1 - al);
}
void Biquad::notch(float fs, float f, float q) {
  float w = kTwoPi * f / fs, c = cosf(w), al = sinf(w) / (2.f * q);
  set(1, -2 * c, 1, 1 + al, -2 * c, 1 - al);
}

void RateBridge::begin(float inRate, float outRate, float lag) {
  nominal_ = inRate / outRate;
  tau_ = outRate;  // lag averaging time constant: 1 s of output samples
  lag_ = lag;
  reset();
}

void RateBridge::reset() {
  idx = written = holds = ticks_ = 0;
  frac = integ_ = lagAvg_ = 0;
  avgReady_ = false;
  step = nominal_;
}

void RateBridge::endBlock() {
  if (ticks_ < 32) return;  // run the control at most every 32 output samples: small blocks stay cheap
  float lag = lagNow();
  if (lag > 2.f * lag_ + 48.f) {  // far behind (start-up, stalled renderer): jump to the target
    uint32_t back = (uint32_t)ceilf(lag_);
    if (written > back) { idx = written - back; frac = (float)back - lag_; }
    lag = lagNow();
    avgReady_ = false;
  }
  // average the lag over ~1 s first: bursty input (DMA blocks) makes it a sawtooth that would otherwise
  // become vibrato
  if (!avgReady_) { lagAvg_ = lag; avgReady_ = true; }
  lagAvg_ += (lag - lagAvg_) * ((float)ticks_ / ((float)ticks_ + tau_));
  // PI control on the averaged lag: P ~4 s time constant at 250 Hz, I removes the steady offset a constant
  // clock mismatch would leave
  float err = lagAvg_ - lag_;
  integ_ += err * 1e-8f * (float)ticks_;
  ticks_ = 0;
  if (integ_ > 0.01f) integ_ = 0.01f;
  if (integ_ < -0.01f) integ_ = -0.01f;
  float trim = err * 0.001f + integ_;
  if (trim > 0.02f) trim = 0.02f;
  if (trim < -0.02f) trim = -0.02f;
  step = nominal_ * (1.f + trim);
}

// ---- front end -------------------------------------------------------------------------------------------

void FrontEnd::begin(float fs, float mainsHz, float hpHz, float lpHz) {
  fs_ = fs;
  hp_.highpass(fs, hpHz, 0.7071f);
  n1_.notch(fs, mainsHz, 8.f);
  useN2_ = 2.f * mainsHz < 0.45f * fs;
  if (useN2_) n2_.notch(fs, 2.f * mainsHz, 8.f);
  lp1_.lowpass(fs, lpHz, 0.5412f);  // 4th-order Butterworth
  lp2_.lowpass(fs, lpHz, 1.3066f);
  ms_ = 0;
  agcGain_ = 1;
  agcN_ = 0;
  setAgc(agc_, target_, agcSeconds_);
}

void FrontEnd::setAgc(bool on, float targetRms, float seconds) {
  agc_ = on;
  target_ = targetRms;
  agcSeconds_ = seconds;
  agcCoef_ = 1.f / (seconds * fs_);
}

float FrontEnd::process(float x) {
  float y = hp_.process(x * gain_);
  y = n1_.process(y);
  if (useN2_) y = n2_.process(y);
  y = lp2_.process(lp1_.process(y));
  if (agc_) {
    float c = 1.f / (float)(++agcN_);  // running mean until the time constant takes over
    if (c < agcCoef_) c = agcCoef_;
    ms_ += c * (y * y - ms_);
    float r = sqrtf(ms_);
    agcGain_ = r > 1e-12f ? target_ / r : 1.f;
    if (agcGain_ > 1e7f) agcGain_ = 1e7f;
    y *= agcGain_;
  }
  return y;
}

// ---- A: frequency shifter --------------------------------------------------------------------------------

// Olli Niemitalo's IIR Hilbert pair (yehar.com): chain I and chain Q (+ one sample) differ by 90 deg
static const float kHilbertI[4] = {0.4021921162426f, 0.8561710882420f, 0.9722909545651f, 0.9952884791278f};
static const float kHilbertQ[4] = {0.6923878f, 0.9360654322959f, 0.9882295226860f, 0.9987488452737f};

void FreqShifter::begin(float inRate, float audioRate, float bandLo, float bandHi, float lag) {
  initTable();
  inRate_ = inRate;
  audioRate_ = audioRate;
  hp_.highpass(inRate, bandLo, 0.7071f);
  lp1_.lowpass(inRate, bandHi, 0.5412f);
  lp2_.lowpass(inRate, bandHi, 1.3066f);
  for (int i = 0; i < 4; ++i) {
    chainI_[i].a2 = kHilbertI[i] * kHilbertI[i];
    chainQ_[i].a2 = kHilbertQ[i] * kHilbertQ[i];
  }
  bridge.begin(inRate, audioRate, lag);
  setShift(shift_);
  reset();
}

void FreqShifter::setShift(float hz) {
  shift_ = hz;
  incOut_ = phaseInc(hz, audioRate_);
}

void FreqShifter::reset() {
  memset(I_, 0, sizeof(I_));
  memset(Q_, 0, sizeof(Q_));
  hp_.reset(); lp1_.reset(); lp2_.reset();
  for (int i = 0; i < 4; ++i) {
    float a = chainI_[i].a2, b = chainQ_[i].a2;
    chainI_[i] = Allpass(); chainI_[i].a2 = a;
    chainQ_[i] = Allpass(); chainQ_[i].a2 = b;
  }
  qDelay_ = 0;
  phOut_ = 0;
  bridge.reset();
}

void FreqShifter::push(float x) {
  x = lp2_.process(lp1_.process(hp_.process(x)));
  float i = x, q = x;
  for (int k = 0; k < 4; ++k) { i = chainI_[k].process(i); q = chainQ_[k].process(q); }
  float qd = qDelay_;  // Q lags I by 90 deg: analytic signal I + jQ
  qDelay_ = q;
  uint32_t w = bridge.written & (kRing - 1);
  I_[w] = i;
  Q_[w] = qd;
  bridge.pushed();
}

void FreqShifter::render(float* out, int n) {
  const uint32_t m = kRing - 1;
  for (int k = 0; k < n; ++k) {
    uint32_t i0 = bridge.idx;
    float t = bridge.frac;
    float ii = cubic(I_[(i0 - 1) & m], I_[i0 & m], I_[(i0 + 1) & m], I_[(i0 + 2) & m], t);
    float qq = cubic(Q_[(i0 - 1) & m], Q_[i0 & m], Q_[(i0 + 1) & m], Q_[(i0 + 2) & m], t);
    out[k] = ii * cosine(phOut_) - qq * sine(phOut_);  // Re{(I + jQ) e^(j shift t)}: upper sideband
    phOut_ += incOut_;
    bridge.tick();
  }
  bridge.endBlock();
}

// ---- B: grain scaler -------------------------------------------------------------------------------------

void GrainScaler::begin(float inRate, float audioRate, float lag) {
  initTable();
  inRate_ = inRate;
  audioRate_ = audioRate;
  bridge.begin(inRate, audioRate, lag);
  setSync(sync_, 4.f);
  reset();
}

void GrainScaler::setRatio(float k) {
  if (k < 0.25f) k = 0.25f;
  if (k > 64.f) k = 64.f;
  bool flip = (k < 1.f) != (k_ < 1.f);
  k_ = k;
  if (flip) spreadHeads();
}

void GrainScaler::setSync(bool on, float minSyncHz) {
  sync_ = on;
  search_ = (int)(inRate_ / (minSyncHz > 0.5f ? minSyncHz : 0.5f));
  corr_ = (int)(0.05f * inRate_);  // ~50 ms comparison window
  if (corr_ < 4) corr_ = 4;
  setWindow(winReq_);
}

void GrainScaler::setWindow(float seconds) {
  winReq_ = seconds;
  float w = seconds * inRate_;
  // ring must hold: interpolation taps, bridge lag, window, splice wander (2 x search), correlation window
  float maxW = (float)kRing - dmin_ - 64.f - (sync_ ? 2.f * search_ + 2.f * corr_ : 0.f);
  if (w > maxW) w = maxW;
  if (w < 16.f) w = 16.f;
  win_ = w;
  spreadHeads();
}

void GrainScaler::setHeads(int n) {
  if (n < 2) n = 2;
  if (n > kMaxHeads) n = kMaxHeads;
  heads_ = n;
  spreadHeads();
}

void GrainScaler::spreadHeads() {
  for (int h = 0; h < heads_; ++h) {
    ph_[h] = (float)h / heads_;
    delay_[h] = k_ >= 1.f ? dmin_ + win_ * (1.f - ph_[h]) : dmin_ + win_ * ph_[h];
  }
}

void GrainScaler::reset() {
  memset(buf_, 0, sizeof(buf_));
  spreadHeads();
  bridge.reset();
}

void GrainScaler::push(float x) {
  buf_[bridge.written & (kRing - 1)] = x;
  bridge.pushed();
}

// New delay for head h, whose window just restarted. Plain: jump back one window. Sync: pick the jump within
// +-search samples whose recent history best matches (normalised correlation) the loudest other head's.
float GrainScaler::splice(int h, float d, uint32_t i0) {
  float base = d + win_;
  if (!sync_ || heads_ < 2) return base;
  int dom = h == 0 ? 1 : 0;
  for (int j = 0; j < heads_; ++j)
    if (j != h && hann01(ph_[j]) > hann01(ph_[dom])) dom = j;
  const uint32_t m = kRing - 1;
  // candidates are whole samples away from the loudest head's exact delay, so both read points share the same
  // fraction and the integer-index comparison is exact (rounding them separately biased every splice)
  const float dDom = delay_[dom];
  const uint32_t rDom = i0 - (uint32_t)(int)dDom;
  // the head travels one window per cycle, so its end-of-cycle delay is (candidate - win): keep it in
  // [dmin, dmin + 2S] and the candidate within +-S of the plain jump
  float cLo = base - search_, cHi = base + search_;
  if (cLo < dmin_ + win_) cLo = dmin_ + win_;
  if (cHi > dmin_ + win_ + 2.f * search_) cHi = dmin_ + win_ + 2.f * search_;
  int lo = (int)ceilf(cLo - dDom), hi = (int)floorf(cHi - dDom);  // offsets from the loudest head
  if (lo > hi) return base;
  float best = -1e30f, sPrev = 0, sBest = 0, sNext = 0;
  int bestK = 0;
  float prev = -1e30f;
  for (int k = lo; k <= hi; ++k) {
    uint32_t r = rDom - (uint32_t)k;  // candidate delay dDom + k
    float dot = 0, en = 1e-20f;
    for (int i = -2 * corr_; i <= 0; ++i) {
      float a = buf_[(rDom + i) & m], b = buf_[(r + i) & m];
      dot += a * b;
      en += b * b;
    }
    float sc = dot / sqrtf(en);
    if (sc > best) { best = sc; bestK = k; sPrev = prev; sBest = sc; sNext = -1e30f; }
    else if (k == bestK + 1) sNext = sc;
    prev = sc;
  }
  float frac = 0.f;  // parabolic refinement of the correlation peak
  if (sPrev > -1e29f && sNext > -1e29f) {
    float den = sPrev - 2.f * sBest + sNext;
    if (den < 0.f) frac = 0.5f * (sPrev - sNext) / den;
  }
  return dDom + (float)bestK + frac;
}

void GrainScaler::render(float* out, int n) {
  const uint32_t m = kRing - 1;
  const float gain = 2.f / heads_;  // evenly spaced Hann windows sum to heads/2
  const float dd = (k_ - 1.f) * bridge.stepNow();  // heads move k times faster than the cursor
  const float dph = (dd < 0.f ? -dd : dd) / win_;   // (the step only changes in endBlock)
  for (int k = 0; k < n; ++k) {
    uint32_t i0 = bridge.idx;
    float t = bridge.frac;
    float s = 0.f;
    for (int h = 0; h < heads_; ++h) {  // 1) read every head at its current position
      float d = delay_[h];
      int di = (int)d;
      float rf = t - (d - di);
      uint32_t ri = i0 - (uint32_t)di;
      if (rf < 0.f) { rf += 1.f; --ri; }
      s += hann01(ph_[h]) * cubic(buf_[(ri - 1) & m], buf_[ri & m], buf_[(ri + 1) & m], buf_[(ri + 2) & m], rf);
    }
    for (int h = 0; h < heads_; ++h) { delay_[h] -= dd; ph_[h] += dph; }  // 2) advance all heads together
    for (int h = 0; h < heads_; ++h) {  // 3) splices see every head at the same instant
      if (ph_[h] < 1.f) continue;
      ph_[h] -= 1.f;
      delay_[h] = k_ >= 1.f ? splice(h, delay_[h], i0) : delay_[h] - win_;
    }
    out[k] = s * gain;
    bridge.tick();
  }
  bridge.endBlock();
}

// ---- C: spectral oscillator bank ---------------------------------------------------------------------------

static const uint16_t kScaleMask[] = {0x000, 0xFFF, 0xAB5, 0x5AD, 0x295, 0x4A9};

const char* SpectralBank::scaleName(int s) {
  static const char* names[] = {"off", "chromatic", "major", "minor", "pentatonic major", "pentatonic minor"};
  return (s >= 0 && s < kScaleCount) ? names[s] : "?";
}

void SpectralBank::begin(float inRate, float audioRate, float bandLo, float bandHi) {
  initTable();
  inRate_ = inRate;
  audioRate_ = audioRate;
  float binHz = inRate / kN;
  binLo_ = (int)ceilf(bandLo / binHz);
  if (binLo_ < 2) binLo_ = 2;                       // peaks need a neighbour on each side
  binHi_ = (int)floorf(bandHi / binHz);
  if (binHi_ > kN / 2 - 2) binHi_ = kN / 2 - 2;
  if (binHi_ - binLo_ + 3 > kMaxBins + 2) binHi_ = binLo_ + kMaxBins - 1;
  rampLen_ = (int)(kHop * audioRate / inRate);
  reset();
}

void SpectralBank::reset() {
  memset(ring_, 0, sizeof(ring_));
  memset(prevPhase_, 0, sizeof(prevPhase_));
  memset(p_, 0, sizeof(p_));
  written_ = analysed_ = 0;
  rampLeft_ = 0;
  level_ = 0;
}

int SpectralBank::partials() const {
  int n = 0;
  for (int j = 0; j < kMaxPartials; ++j) n += (p_[j].live || p_[j].amp > 0.f) ? 1 : 0;
  return n;
}

void SpectralBank::push(float x) {
  ring_[written_ & (kN - 1)] = x;
  ++written_;
}

float SpectralBank::quantize(float hz) const {
  if (scale_ == kOff || hz <= 0.f) return hz;
  uint16_t mask = kScaleMask[scale_];
  float semis = 12.f * log2f(hz / root_);
  int base = (int)floorf(semis);
  float best = (float)base, err = 1e9f;
  for (int s = base - 6; s <= base + 7; ++s) {
    int deg = ((s % 12) + 12) % 12;
    if (!((mask >> deg) & 1)) continue;
    float e = fabsf(semis - (float)s);
    if (e < err) { err = e; best = (float)s; }
  }
  return root_ * powf(2.f, best / 12.f);
}

bool SpectralBank::update() {
  if (written_ - analysed_ < (uint32_t)kHop) return false;
  uint32_t hop = written_ - analysed_;  // normally kHop; larger if update() was called late
  analysed_ = written_;
  for (int i = 0; i < kN; ++i) {
    float w = 0.5f - 0.5f * kSineTable[(i * (1024 / kN) + 256) & 1023];  // Hann from the table (exact for kN | 1024)
    re_[i] = w * ring_[(written_ - kN + i) & (kN - 1)];
    im_[i] = 0.f;
  }
  fft(re_, im_, kN);

  // magnitude and phase-vocoder frequency for bins binLo-1 .. binHi+1 (stored from index 0)
  const float binHz = inRate_ / kN;
  const bool phaseOk = hop <= (uint32_t)(kN / 4);  // beyond this the phase advance is ambiguous
  const int nb = binHi_ - binLo_ + 3;
  float mag[kMaxBins + 2], fin[kMaxBins + 2], maxMag = 0.f;
  for (int j = 0; j < nb; ++j) {
    int b = binLo_ - 1 + j;
    float r = re_[b], q = im_[b];
    mag[j] = sqrtf(r * r + q * q) * (4.f / kN);  // sine amplitude (Hann coherent gain 0.5)
    float ph = atan2f(q, r);
    float dev = ph - prevPhase_[j] - kTwoPi * (float)b * (float)hop / kN;
    dev -= kTwoPi * floorf((dev + kPi) / kTwoPi);  // wrap to [-pi, pi)
    prevPhase_[j] = ph;
    fin[j] = ((float)b + (phaseOk ? dev * kN / (kTwoPi * (float)hop) : 0.f)) * binHz;
    if (mag[j] > maxMag) maxMag = mag[j];
  }

  // peaks: local maxima above -45 dB of the strongest, amplitude by parabolic interpolation of log magnitude
  struct Peak { float f, mag; };
  Peak pk[kMaxBins / 2 + 1];
  int np = 0;
  const float floor = maxMag * 0.0056f;
  for (int j = 1; j < nb - 1 && np < kMaxBins / 2; ++j) {
    if (!(mag[j] > mag[j - 1] && mag[j] >= mag[j + 1] && mag[j] > floor)) continue;
    float a = logf(mag[j - 1] + 1e-20f), b = logf(mag[j]), g = logf(mag[j + 1] + 1e-20f);
    float den = a - 2.f * b + g, d = den < 0.f ? 0.5f * (a - g) / den : 0.f;
    pk[np].f = fin[j];
    pk[np].mag = expf(b - 0.25f * (a - g) * d);
    ++np;
  }
  // strongest first; drop peaks within merge distance of a stronger one (sidelobes report the same frequency)
  for (int i = 1; i < np; ++i)
    for (int j = i; j > 0 && pk[j].mag > pk[j - 1].mag; --j) { Peak t = pk[j]; pk[j] = pk[j - 1]; pk[j - 1] = t; }
  int nk = 0;
  for (int i = 0; i < np && nk < kMaxPartials; ++i) {
    bool dup = pk[i].f <= 0.05f;
    for (int j = 0; j < nk && !dup; ++j) dup = fabsf(pk[i].f - pk[j].f) < 1.5f * binHz;
    if (!dup) pk[nk++] = pk[i];
  }

  // amplitudes, output frequencies and loudness normalisation
  float amp[kMaxPartials], fout[kMaxPartials], energy = 0.f;
  for (int i = 0; i < nk; ++i) {
    float a = powf(pk[i].mag, gamma_) * powf(pk[i].f / 10.f, tilt_);
    float fo = quantize(k_ * pk[i].f);
    if (fo < 20.f || fo > 0.45f * audioRate_) a = 0.f;
    amp[i] = a;
    fout[i] = fo;
    energy += a * a;
  }
  float e = sqrtf(energy);
  float rel = expf(-(float)hop / (release_ * inRate_));  // instant attack, release_ seconds
  level_ = e > level_ ? e : level_ * rel;
  float norm = level_ > 1e-20f ? 0.35f / level_ : 0.f;  // ~0.25 rms at the recent peak

  // tracking: each peak continues the nearest sounding oscillator, else takes a silent one
  bool used[kMaxPartials] = {};
  for (int i = 0; i < nk; ++i) {
    float tol = 1.5f * binHz > 0.05f * pk[i].f ? 1.5f * binHz : 0.05f * pk[i].f;
    int best = -1;
    for (int o = 0; o < kMaxPartials; ++o) {
      if (used[o] || !(p_[o].live || p_[o].amp > 0.f)) continue;
      float d = fabsf(p_[o].fIn - pk[i].f);
      if (d < tol) { tol = d; best = o; }
    }
    if (best < 0)
      for (int o = 0; o < kMaxPartials && best < 0; ++o)
        if (!used[o] && !p_[o].live && p_[o].amp == 0.f) best = o;
    if (best < 0) continue;  // pool full
    Partial& p = p_[best];
    used[best] = true;
    p.live = true;
    p.fIn = pk[i].f;
    p.inc = phaseInc(fout[i], audioRate_);
    p.target = amp[i] * norm;
  }
  for (int o = 0; o < kMaxPartials; ++o) {
    Partial& p = p_[o];
    if (!used[o]) { p.live = false; p.target = 0.f; }  // fade out
    p.dAmp = (p.target - p.amp) / (float)rampLen_;
  }
  rampLeft_ = rampLen_;
  return true;
}

// Partials outer (two per pass), samples inner: the partials' state lives in registers for the whole block and
// each out[k] is loaded and stored once per pair; the inner loop is two table lookups and multiply-adds.
void SpectralBank::render(float* out, int n) {
  for (int k = 0; k < n; ++k) out[k] = 0.f;
  const int ramp = rampLeft_ < n ? rampLeft_ : n;  // samples of this block still inside the amplitude ramp
  const bool rampEnds = rampLeft_ > 0 && rampLeft_ <= n;
  int j = 0;
  for (; j + 1 < kMaxPartials; j += 2) {
    Partial& p = p_[j];
    Partial& q = p_[j + 1];
    uint32_t ph1 = p.phase, ph2 = q.phase;
    const uint32_t inc1 = p.inc, inc2 = q.inc;
    float a1 = p.amp, a2 = q.amp;
    const float da1 = p.dAmp, da2 = q.dAmp;
    if (a1 == 0.f && a2 == 0.f && ((da1 == 0.f && da2 == 0.f) || ramp == 0)) {  // silent pair: phases only
      p.phase = ph1 + inc1 * (uint32_t)n;
      q.phase = ph2 + inc2 * (uint32_t)n;
      continue;
    }
    int k = 0;
    for (; k < ramp; ++k) {
      out[k] += a1 * sine(ph1) + a2 * sine(ph2);
      ph1 += inc1; ph2 += inc2; a1 += da1; a2 += da2;
    }
    if (rampEnds) { a1 = p.target; a2 = q.target; p.dAmp = q.dAmp = 0.f; }
    for (; k < n; ++k) {
      out[k] += a1 * sine(ph1) + a2 * sine(ph2);
      ph1 += inc1; ph2 += inc2;
    }
    p.phase = ph1; q.phase = ph2;
    p.amp = a1; q.amp = a2;
  }
  if (j < kMaxPartials) {  // odd count: last partial alone (locals: out[] could alias the struct)
    Partial& p = p_[j];
    uint32_t ph = p.phase;
    const uint32_t inc = p.inc;
    float a = p.amp;
    const float da = p.dAmp;
    int k = 0;
    for (; k < ramp; ++k) { out[k] += a * sine(ph); ph += inc; a += da; }
    if (rampEnds) { a = p.target; p.dAmp = 0.f; }
    for (; k < n; ++k) { out[k] += a * sine(ph); ph += inc; }
    p.phase = ph;
    p.amp = a;
  }
  rampLeft_ -= ramp;
}

// ---- interpolator ------------------------------------------------------------------------------------------

static float besselI0(float x) {  // modified Bessel function of the first kind, order 0 (series)
  float sum = 1.f, term = 1.f, q = 0.25f * x * x;
  for (int k = 1; k < 50 && term > 1e-9f * sum; ++k) { term *= q / ((float)k * (float)k); sum += term; }
  return sum;
}

void Interpolator::begin(int factor, int taps, float pass) {
  initTable();
  if (factor < 1) factor = 1;
  if (factor > kMaxFactor) factor = kMaxFactor;
  if (taps < 4) taps = 4;
  if (taps > kMaxTaps) taps = kMaxTaps;
  if (pass > 0.49f) pass = 0.49f;
  L_ = factor;
  T_ = taps;
  const int n = L_ * T_;
  const float centre = 0.5f * (float)(n - 1);
  const float fc = 0.5f / (float)L_;                                  // half the input rate, in output cycles
  const float transition = (1.f - 2.f * pass) / (float)L_;
  float atten = 2.285f * kTwoPi * transition * (float)(n - 1) + 8.f;  // Kaiser's estimate
  float beta = atten > 50.f ? 0.1102f * (atten - 8.7f) : atten > 21.f ? 0.5842f * powf(atten - 21.f, 0.4f) + 0.07886f * (atten - 21.f) : 0.f;
  float i0b = besselI0(beta);
  for (int m = 0; m < n; ++m) {
    float x = (float)m - centre;
    float s = x == 0.f ? 2.f * fc : sinf(kTwoPi * fc * x) / (kPi * x);
    float r = 2.f * x / (float)(n - 1);
    float w = besselI0(beta * sqrtf(fmaxf(0.f, 1.f - r * r))) / i0b;
    int p = m % L_, t = m / L_;
    h_[p * T_ + t] = s * w;
  }
  for (int p = 0; p < L_; ++p) {  // unity DC gain per phase: no ripple on slow CV
    float sum = 0.f;
    for (int t = 0; t < T_; ++t) sum += h_[p * T_ + t];
    if (sum != 0.f)
      for (int t = 0; t < T_; ++t) h_[p * T_ + t] /= sum;
  }
  reset();
}

void Interpolator::reset() {
  memset(hist_, 0, sizeof(hist_));
  pos_ = 0;
}

void Interpolator::process(const float* in, int n, float* out) {
  for (int i = 0; i < n; ++i) {
    pos_ = pos_ == 0 ? T_ - 1 : pos_ - 1;
    hist_[pos_] = hist_[pos_ + T_] = in[i];  // mirrored: hist_[pos_ .. pos_ + T_ - 1] is contiguous, newest first
    const float* x = hist_ + pos_;
    int p = 0;
    for (; p + 1 < L_; p += 2) {  // two phases per pass: each history value is loaded once for both
      const float* h0 = h_ + p * T_;
      const float* h1 = h0 + T_;
      float s0 = 0.f, s1 = 0.f;
      for (int t = 0; t < T_; ++t) {
        float v = x[t];
        s0 += h0[t] * v;
        s1 += h1[t] * v;
      }
      out[p] = s0;
      out[p + 1] = s1;
    }
    if (p < L_) {
      const float* h0 = h_ + p * T_;
      float s0 = 0.f;
      for (int t = 0; t < T_; ++t) s0 += h0[t] * x[t];
      out[p] = s0;
    }
    out += L_;
  }
}

// ---- pass-through ------------------------------------------------------------------------------------------

void Upsampler::begin(float inRate, float audioRate, float lag) {
  initTable();
  bridge.begin(inRate, audioRate, lag);
  memset(buf_, 0, sizeof(buf_));
}

void Upsampler::push(float x) {
  buf_[bridge.written & (kRing - 1)] = x;
  bridge.pushed();
}

void Upsampler::render(float* out, int n) {
  const uint32_t m = kRing - 1;
  for (int k = 0; k < n; ++k) {
    uint32_t i0 = bridge.idx;
    out[k] = cubic(buf_[(i0 - 1) & m], buf_[i0 & m], buf_[(i0 + 1) & m], buf_[(i0 + 2) & m], bridge.frac);
    bridge.tick();
  }
  bridge.endBlock();
}

// ---- test signal -------------------------------------------------------------------------------------------

void TestSignal::begin(float fs, uint32_t seed) {
  initTable();
  fs_ = fs;
  rng_ = seed ? seed : 1;
  n_ = 0;
  b0_ = b1_ = b2_ = 0;
  pulseAt_ = (uint32_t)(3.f * fs);
  phBurst_ = phLow_ = phHigh_ = phMains_ = phDrift_ = phWobble_ = 0;
}

float TestSignal::white() {
  rng_ ^= rng_ << 13;
  rng_ ^= rng_ >> 17;
  rng_ ^= rng_ << 5;
  return (float)(rng_ >> 8) * (2.f / 16777216.f) - 1.f;
}

float TestSignal::next() {
  float w = white();
  b0_ = 0.99765f * b0_ + w * 0.0990460f;  // economy pink noise (Kellet)
  b1_ = 0.96300f * b1_ + w * 0.2965164f;
  b2_ = 0.57000f * b2_ + w * 1.0526913f;
  float pink = b0_ + b1_ + b2_ + w * 0.1848f;

  // burst envelope: 8 s cycle, swelled for 4 s with 0.5 s raised-cosine edges
  uint32_t cycle = (uint32_t)(8.f * fs_);
  float c = (float)(n_ % cycle) / fs_, env;
  if (c < 0.5f) env = 0.5f - 0.5f * cosf(kPi * c / 0.5f);
  else if (c < 4.f) env = 1.f;
  else if (c < 4.5f) env = 0.5f + 0.5f * cosf(kPi * (c - 4.f) / 0.5f);
  else env = 0.f;

  float wobble = 0.03f * burstHz_ * sine(phWobble_);  // burst frequency drifts +-3 %
  float burst = (burstLow_ + (burstHigh_ - burstLow_) * env) * sine(phBurst_);
  phBurst_ += phaseInc(burstHz_ + wobble, fs_);
  phWobble_ += phaseInc(0.17f, fs_);
  float low = 5.f * sine(phLow_);     phLow_ += phaseInc(6.f, fs_);
  float high = 3.f * sine(phHigh_);   phHigh_ += phaseInc(20.f, fs_);
  float mains = mainsLevel_ * sine(phMains_);  phMains_ += phaseInc(mainsHz_, fs_);
  float drift = 40.f * sine(phDrift_);  phDrift_ += phaseInc(0.05f, fs_);

  float pulse = 0.f;  // 0.15 s Gaussian transient every 6-8 s
  float dt = ((float)n_ - (float)pulseAt_) / fs_;
  if (dt > -0.5f && dt < 0.5f) pulse = 100.f * expf(-(dt * dt) / (2.f * 0.075f * 0.075f));
  if (dt >= 0.5f) pulseAt_ = n_ + (uint32_t)((6.5f + white()) * fs_);
  ++n_;
  return 8.f * pink + burst + low + high + mains + drift + pulse;
}

// ---- measurement -------------------------------------------------------------------------------------------

void ZeroCrossHz::add(float x) {
  if (x < -hyst_) armed_ = true;
  if (armed_ && prev_ < 0.f && x >= 0.f) {
    double t = (double)(n_ - 1) + (double)(-prev_ / (x - prev_));
    if (first_ < 0) first_ = t;
    last_ = t;
    ++count_;
    armed_ = false;
  }
  prev_ = x;
  ++n_;
}

}  // namespace fm
