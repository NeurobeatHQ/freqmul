# FreqMul

Real-time conversion of low-frequency signals (about 1–40 Hz) into audio. Plain C++ (C++11 or later), no
dependencies, no `Arduino.h`. Runs on the RA4M1 (Arduino UNO R4), the ESP32-S3 and desktops. Everything is in
namespace `fm`.

## 0. Installation

The folder is an Arduino library, an ESP-IDF component and a plain CMake library at the same time.

```
FreqMul/
  src/FreqMul.h, FreqMul.cpp, SineTable.h   the library (SineTable.h: a const table that stays in flash)
  library.properties, keywords.txt          Arduino
  CMakeLists.txt, idf_component.yml         ESP-IDF component; plain CMake library elsewhere
  examples/BasicRender/                     Arduino example: every method, frequency check + CPU load
  extras/idf_example/                       the same as an ESP-IDF app
  extras/host_test/                         desktop test suite
```

- **Arduino IDE**: copy or symlink the `FreqMul` folder into your sketchbook's `libraries/` folder, or zip it and
  use *Sketch → Include Library → Add .ZIP Library*. Then `#include <FreqMul.h>`.
- **arduino-cli**: `arduino-cli compile --library /path/to/FreqMul ...`, or install it with
  `arduino-cli lib install --zip-path FreqMul.zip`.
- **PlatformIO**: `lib_deps = symlink:///path/to/FreqMul` in `platformio.ini`.
- **ESP-IDF**: add the folder to the project's component search path, then `REQUIRES FreqMul` in your component:
  ```cmake
  # top-level CMakeLists.txt, before include(project.cmake); APPEND keeps paths PlatformIO adds
  list(APPEND EXTRA_COMPONENT_DIRS "/path/to/FreqMul")
  ```
  The example is `extras/idf_example`: `idf.py set-target esp32s3 && idf.py build flash monitor`.
- **Plain CMake**: `add_subdirectory(FreqMul)` and `target_link_libraries(app PRIVATE FreqMul)`.

## 1. Quick start

```cpp
#include "FreqMul.h"

fm::FrontEnd     fe;     // drift high-pass, mains notches, band-limit, slow AGC
fm::GrainScaler  grain;  // method B: every frequency x k
fm::Interpolator up;     // internal rate -> DAC rate

void setup() {
  fe.begin(250, 50);              // input rate 250 Hz, 50 Hz mains
  fe.setAgc(true);
  grain.begin(250, 16000, 12);    // input rate, internal rate, lag (input samples, see section 5)
  grain.setRatio(32);             // 10 Hz -> 320 Hz (five octaves)
  up.begin(3);                    // 16 kHz x 3 = 48 kHz
}

void onInputSample(float x) { grain.push(fe.process(x)); }   // at the input rate

void fillDacBlock(float* dac48k /* 192 */) {                  // at the output rate, any block size
  float buf[64];
  grain.render(buf, 64);          // 64 samples at 16 kHz
  up.process(buf, 64, dac48k);    // -> 192 samples at 48 kHz
}
```

The pattern is the same for every method: `push()` at the input rate, `render(out, n)` at the method's internal
rate, `Interpolator::process()` up to the DAC rate. Method C also needs `update()` called from the main loop.

## 2. The three methods

| | Class | Mapping | Sounds like | Latency |
|---|---|---|---|---|
| A | `FreqShifter` | f → f + shift (single sideband) | metallic; spacing in Hz kept | ~20 ms + lag |
| B | `GrainScaler` | f → k·f (delay-line pitch shifter) | the input's texture, sped up; ratios kept | ~window/2 + ~0.1 s (0.62 s at the default 1 s window) |
| C | `SpectralBank` | f → k·f per spectral peak, optionally snapped to a scale | organ or pad; can stay in key | ~2 s |

```cpp
fm::FreqShifter a;  a.begin(250, 4000, 1, 40, 12);  a.setShift(200);     // 10 Hz -> 210 Hz
fm::GrainScaler b;  b.begin(250, 16000, 12);        b.setRatio(32);      // setWindow, setHeads, setSync
fm::SpectralBank c; c.begin(250, 8000);             c.setRatio(32);      // setQuantize, setTilt, setCompress
c.setQuantize(fm::SpectralBank::kPentatonicMinor, 110);                  // partials snap to A minor pentatonic
// in loop(): c.update();   // runs the FFT every 64 input samples; never call it from an interrupt
```

- **A**: an allpass Hilbert pair (Niemitalo's coefficients, ±0.7° from 0.0008 to 0.5 × the rate) followed by
  single-sideband mixing. The mirror image is ≥ 44 dB down (measured −45 dB).
- **B**: read heads sweep the last `window` seconds k times faster and crossfade. With sync on (SOLA) each splice
  is moved to where it matches the louder head best, so pitch stays exact: a 10.3 Hz tone at ×32 measures
  329.5 Hz against 329.6 expected, and plain splices give 320.2. `setWindow()` sets B's latency, its CPU and the
  lowest frequency it handles: see *Tuning B's latency* in section 5.
- **C**: a 512-point STFT every 64 samples. Local spectral maxima above −45 dB are picked, duplicates are merged,
  frequencies come from the phase vocoder and amplitudes from parabolic interpolation. Peaks are tracked
  frame to frame, so each component is one phase-continuous sine. Measured stray components ≤ −58 dB at
  ratios ×2 to ×64. It uses about 23 oscillators on noisy input (pool of 40).

`Upsampler` passes the input through unchanged at audio rate, for scopes and A/B comparisons. `TestSignal`
generates a synthetic test input (pink noise, a tone that swells 4 s out of every 8, tones, hum, drift, pulses).
`ZeroCrossHz` measures a tone's frequency.

## 3. Ratios

`setRatio(k)` on B and C takes any real number. Powers of two (2, 4, …, 64) are octave shifts that keep each
note's pitch class; 3, 6 and 12 work as well. Measured accuracy:

| k | 2 | 3 | 4 | 6 | 8 | 12 | 16 | 32 | 64 |
|---|---|---|---|---|---|---|---|---|---|
| B error | 0.03 % | 0.03 % | 0.03 % | 0.02 % | 0.02 % | 0.02 % | 0.02 % | 0.02 % | 0.41 % |
| C | exact | exact | exact | exact | exact | exact | exact | exact | exact |

Limits:
- The output band is k × [bandLo, bandHi].
- C mutes partials that would land below 20 Hz or above 0.45 × its rate.
- B clamps k to 0.25–64, and its splice cost grows with k − 1.
- A is additive (Hz), so ratios don't apply to it.

## 4. Sample rates: internal rate plus Interpolator

Render each method only as fast as its output needs, then interpolate to the DAC rate:

| Method | Highest output | Internal rate | Factor to 48 kHz |
|---|---|---|---|
| A | 40 Hz + shift (≤ 1.5 kHz) | 4 kHz | 12 |
| B | 40·k, plus input images at k·(inRate − 40) | 16 kHz | 3 |
| C | 40·k (pure sines) | 4 kHz up to ×32, 8 kHz up to ×90 | 12 / 6 |

```cpp
fm::Interpolator up;  up.begin(/*factor*/ 6, /*taps*/ 16, /*pass*/ 0.4f);
up.process(in, n, out);   // n inputs -> n * 6 outputs; streaming, any block size
```

- **Interpolator**: a Kaiser-windowed sinc, polyphase, factors 1–16 and 4–24 taps per phase. Each phase has unity
  DC gain, so CV passes exactly. Measured: passband ±0.01 dB up to 0.4 × the input rate, images −55 to −65 dB,
  output identical for any block size. Delay is taps / 2 input samples.
- **Why the DAC rate should be ~48 kHz**: a DAC holds each value until the next, which leaves images at
  fs − f, and without an analog filter they reach the output. At fs = 16 kHz the image of 1 kHz sits at
  15 kHz, −23.5 dB (300 Hz: −34 dB). A wavefolder after the DAC acts like a digital waveshaper running at the DAC
  rate, so harmonics above fs/2 fold back. 48 kHz moves that limit from the 6th harmonic of 1.28 kHz to about
  the 19th. Mutable Instruments' Braids drives a similar string DAC at 96 kHz.
- **Analog alternative**: a single capacitor across the feedback resistor of a non-inverting ×8 stage can't
  attenuate below about 1/gain (−18 dB). An analog reconstruction filter needs to be 2nd order or higher.

## 5. Clocks, latency and bursty input

The input and output clocks may differ, and input may arrive in DMA bursts. Each method follows the input with a
`RateBridge` cursor that trails the newest input sample by `lag` input samples. A PI controller acts on a ~1 s
average of the lag and trims the read speed by at most 2 %.

```
lag >= 4 + burst length + (output FIFO length / (output rate / input rate))     // all in input samples
```

| Situation (input 250 Hz) | lag | Measured |
|---|---|---|
| sample-by-sample input | 4 | lag held within 0.5 sample |
| input clock 0.4 % fast (251 Hz) | 4 | pitch follows the true input (319.5 Hz vs 320) |
| 25-sample bursts every 100 ms, ±2 ms jitter | 17–30 | no stalls after start-up; B within −9/+6 cents |
| bursts with lag 4 (too small) | 4 | cursor stalls at every burst: −30 cent dips |

Latency is lag/inRate plus the method's own delay (A ~20 ms, B ~window/2 + ~0.1 s, C ~2 s) plus output buffering.

### Tuning B's latency

Measured at ×32 with 2 heads, sync on and lag 4. The tone columns give the share of output energy within ±2 Hz of
target. CPU is estimated as splices/s × 62–93 candidate positions per splice × ~480 cycles; on the ESP32-S3,
divide by about 5.

| `setWindow` | Latency | 10 Hz | 6 Hz | 6+10 Hz mix | Splices/s | Splice CPU, R4 |
|---|---|---|---|---|---|---|
| 1 s (default) | 0.62 s | 100 % | 100 % | 100 % | 62 | ~4 % |
| 0.5 s | 0.36 s | 100 % | 100 % | 100 % | 124 | ~12 % |
| 0.25 s | 0.25 s | 100 % | 100 % | 99 % | 248 | ~23 % |
| 0.125 s | 0.15 s | 100 % | 100 % | 3 % | 496 | ~45 % |

- **`setWindow(seconds)` trades latency against CPU and low-frequency coverage.** Each halving of the window
  takes about window/4 off the latency and doubles the splices. Splice cost also grows with k − 1.
- **The window must hold about 1.5 periods of the lowest frequency**: 0.25 s for 6 Hz and up, 0.5 s for 3 Hz and
  up. At 0.125 s a pure 6 Hz tone still passes, but the 6+10 Hz mix falls apart.
- **Total latency** = DMA block length (input waits in the ADC buffer) + lag/inRate (16 ms at lag 4) + ~window/2
  \+ ~0.1 s of splice wander (0.07–0.11 s measured) + output buffering.
- **Don't shrink the splice search to save CPU.** The search spans ±1/minSyncHz seconds (`setSync(true, 4)`:
  ±0.25 s). Capping it at half a window (minSyncHz 8 at a 0.25 s window) dropped the mix to 37 % on target.

## 6. Performance

Cycles are estimated from the compiled Cortex-M4 code. The `test` command in the `freq-dup-test` sketch measures
the real numbers on a UNO R4.

| Part | Cost | UNO R4, 48 MHz | ESP32-S3, 240 MHz, one core |
|---|---|---|---|
| A at 4 kHz | ~100 cycles/sample | ~1 % | < 0.5 % |
| B at 16 kHz, 2 heads | ~150 cycles/sample + splices | ~5 % + splices: ~4 / 12 / 23 % at a 1 / 0.5 / 0.25 s window (∝ (k − 1) / window) | ~2–6 % |
| C at 8 kHz | ~10.5 cycles per oscillator-sample, ~23 oscillators | ~4 % | ~1 % |
| C analysis | FFT + peaks every 256 ms | ~2–3 ms burst, ~1 % average | ~0.5 ms |
| Interpolator to 48 kHz | ~6.5 cycles × taps per output sample | ~10 % at 16 taps | ~2 % |

What matters for performance:
- **Build flags**: `FreqMul.cpp` compiles itself at `-O2` with a pragma, and the per-sample helpers are force-inlined.
  At the Arduino default `-Os` they became function calls, which cost 2–3× more.
- **Block size**: A and B pay about 10 % more at 4-sample blocks than at 64. C sets up its 40-slot oscillator
  pool once per block, which adds about 15 % at 16-sample blocks and about 60 % at 4. Use 16 or more samples per
  block. Output is identical for any block size.
- **Bursts versus average**: B's splice search (~0.6–0.9 ms) and C's analysis (~2–3 ms on the R4) are short bursts.
  The output buffer must cover them, or the bursts must be moved out of the audio interrupt.
- **Levers**: C's cost depends on its oscillator count and internal rate; B's on heads, window and k (each splice
  tries 62–93 positions at ~480 cycles each; a cheaper search with running sums and no sqrt/divide could cut that
  2–4×);
  the interpolator's on taps (12 taps: images −42 dB, 25 % cheaper) and DAC rate.
- **On the ESP32-S3**: plain C runs at roughly Cortex-M4 cycle counts, with 5× the clock and two cores (ESP-DSP
  benchmarks: a 512-point FFT takes 80k cycles in ANSI C and 44.6k in assembly). There is no built-in DAC (use I2S
  or SPI). Put render code in IRAM (`IRAM_ATTR`) and the table in RAM to avoid cache stalls, and pin audio to the
  core that isn't running WiFi.

## 7. Memory (RAM per object)

| FrontEnd | FreqShifter | GrainScaler | SpectralBank | Interpolator | Upsampler |
|---|---|---|---|---|---|
| 176 B | 828 B | 2.2 KB (2 s ring at 250 Hz) | 7.7 KB | 1.7 KB | 304 B |

The sine table (4 KB) is `const` and stays in flash.

## 8. Rules and known limits

- **Contexts**: call `push()` and `render()` from the same context, or make the input counter update atomic.
  `SpectralBank::update()` belongs in the main loop. If `render()` runs in an interrupt while `update()` runs in
  the loop, they share oscillator parameters: double-buffer them first (not done yet).
- **B's grain sidebands**: crossfades leave sidebands around each tone. Measured at a 1 s window with 2 heads:
  −75 dB (×32, 10 Hz) to −49 dB (×64, 23 Hz).
- **B with 4 heads at a 1 s window is broken**: a 6 Hz tone keeps only ~57 % of its energy on target (cause not
  traced yet). Use 2 heads (the default), or 4 heads only at windows ≤ 0.5 s, where they measure clean.
- **C's latency**: C needs a ~2 s analysis window; its time resolution is 256 ms.
- **A's range**: keep the shift ≤ 0.4 × the internal rate (1.5 kHz at 4 kHz).

## 9. Testing

```bash
cd extras/host_test && c++ -O2 -std=c++17 -I../../src host_test.cpp ../../src/FreqMul.cpp -o host_test && ./host_test
```

The desktop suite checks:
- the frequency mapping of every method for single and paired tones, at many ratios;
- B at 1, 0.5 and 0.25 s windows with 10 Hz, 6 Hz and a 6+10 Hz mix (≥ 95 % of energy on target);
- A's image rejection and C's stray components;
- the interpolator's DC, passband, images and streaming;
- the full internal-rate → 48 kHz chain;
- 4-, 7- and 16-sample blocks (C is bit-identical across block sizes);
- a 0.4 % clock mismatch.

It also writes demo WAVs to `extras/host_test/out/` (next to the binary, from any working directory). On hardware, `examples/BasicRender` (Arduino) and `extras/idf_example`
(ESP-IDF) re-check the frequencies of every method and print the CPU load on that chip.
