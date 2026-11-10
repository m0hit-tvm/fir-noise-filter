# FIR noise filter: DSP interview problem

A complete, self-contained Python solution to an entry-level DSP interview problem, written to be read aloud in an interview with the key formulas in comments.

## Setup and run

```bash
pip install -r requirements.txt
python fir_noise_filter.py
```

The script prints the results and saves a plot to `plots/fir_filter.png`.

## Problem: FIR noise filtering of a sensor signal

### Statement

A sensor is sampled at `fs = 1000 Hz`. The wanted signal is a 50 Hz sine of
amplitude 1.0. The measurement also contains a 120 Hz interference tone of
amplitude 0.5 and white Gaussian noise with standard deviation 0.3.

1. Generate 1 s of the measured signal.
2. Compute its amplitude spectrum and identify the tones present.
3. Design a linear-phase FIR low-pass filter with cutoff 80 Hz using the
   windowed-sinc method, by hand (no `firwin`), then verify against
   `scipy.signal.firwin`.
4. Apply the filter by direct convolution and compensate the group delay.
5. Report the attenuation at 50 Hz and 120 Hz and the SNR before and after.

### How the program solves it

* The spectrum uses `rfft` on a Hann-windowed signal, scaled by `2/sum(w)`
  so a sine of amplitude A appears as a peak of height A.
* The FIR taps are `h[n] = 2 fc/fs * sinc(2 fc/fs (n - M/2))` times a
  symmetric Hamming window, normalised to unity DC gain. With 101 taps the
  manual taps match `firwin` to machine precision.
* The filter is applied with `np.convolve(..., "full")` and the first
  `(N-1)/2 = 50` samples are dropped, which removes the 50 ms group delay.
* Typical printed results: about 0 dB at 50 Hz, about -60 dB at 120 Hz, SNR
  rising from roughly 3.6 dB to 15.7 dB.

### Questions an interviewer may ask (and short answers)

1. **Why does a 1000 Hz sample rate work here?** Nyquist: the highest
   frequency present (120 Hz) is well below fs/2 = 500 Hz, so no aliasing.
   An anti-aliasing analog filter is still needed before the ADC for
   real-world broadband noise.
2. **Why apply a Hann window before the FFT?** The record length is not an
   exact multiple of every period, so the rectangular window leaks energy
   into neighbouring bins. Hann trades a slightly wider main lobe for much
   lower side lobes, so small tones are not hidden by leakage.
3. **FIR vs IIR: why FIR here?** FIR is always stable, can be exactly
   linear-phase (symmetric taps), and is easy to implement. IIR achieves the
   same stop-band with far fewer coefficients but has non-linear phase and
   possible numerical stability issues at high order.
4. **What determines the transition width of a windowed-sinc filter?**
   Approximately `window_constant * fs / N`. For Hamming that is about
   `3.3 fs/N`, so 101 taps at 1000 Hz gives about 33 Hz. More taps means a
   sharper filter but more delay and more multiply-accumulates per sample.
5. **What is group delay and how did you remove it?** A symmetric N-tap FIR
   delays every frequency by `(N-1)/2` samples. Trimming that many samples
   from the output re-aligns it. In offline processing `filtfilt` gives zero
   phase by filtering forwards and backwards (but squares the magnitude
   response and cannot be used in real time).
6. **Why normalise the taps to sum to 1?** So the pass-band (DC) gain is
   exactly 0 dB. Truncating and windowing the ideal sinc slightly changes the
   DC gain.
7. **How would you implement this on a microcontroller?** A circular buffer
   of the last N inputs and a multiply-accumulate loop, in fixed point if
   there is no FPU. Exploit tap symmetry to halve the multiplies.
8. **How is SNR defined here?** Signal power over error power in dB, where
   the error is the difference from the clean 50 Hz reference. Edge samples
   affected by the filter start-up transient are excluded.
