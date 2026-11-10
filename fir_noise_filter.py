#!/usr/bin/env python3
"""
FIR noise filter: cleaning a noisy sensor signal with a linear-phase FIR low-pass.

Scenario
--------
A vibration sensor is sampled at fs = 1000 Hz. The signal of interest is a
50 Hz tone (amplitude 1.0). It is corrupted by:
  * a 120 Hz interference tone (amplitude 0.5), and
  * additive white Gaussian noise (sigma = 0.3).

Tasks
-----
1. Generate the signal and plot it in the time domain.
2. Compute the FFT and identify the frequency components (with a Hann window).
3. Design a linear-phase FIR low-pass filter (windowed-sinc, cutoff 80 Hz)
   BY HAND and compare its taps with scipy.signal.firwin.
4. Apply the filter by direct convolution, then remove the group delay.
5. Compare SNR before and after filtering.
6. Plot everything to plots/fir_filter.png.

Run:  python fir_noise_filter.py
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import signal

OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "plots")
os.makedirs(OUT_DIR, exist_ok=True)


# --------------------------------------------------------------------------
# 1. Signal generation
# --------------------------------------------------------------------------
def generate_signal(fs=1000.0, duration=1.0, seed=0):
    """Return time vector, clean signal, and noisy measured signal."""
    rng = np.random.default_rng(seed)
    t = np.arange(0, duration, 1.0 / fs)
    clean = 1.0 * np.sin(2 * np.pi * 50 * t)
    interference = 0.5 * np.sin(2 * np.pi * 120 * t)
    noise = 0.3 * rng.standard_normal(t.size)
    measured = clean + interference + noise
    return t, clean, measured


# --------------------------------------------------------------------------
# 2. Spectrum
# --------------------------------------------------------------------------
def amplitude_spectrum(x, fs, window="hann"):
    """
    Single-sided amplitude spectrum, scaled so a sine of amplitude A shows
    a peak of ~A. A window is applied to reduce spectral leakage; the window
    gain is compensated by dividing by sum(w).
    """
    n = x.size
    w = signal.get_window(window, n) if window else np.ones(n)
    X = np.fft.rfft(x * w)
    freqs = np.fft.rfftfreq(n, d=1.0 / fs)
    amp = 2.0 * np.abs(X) / np.sum(w)
    amp[0] /= 2.0            # DC bin is not doubled
    return freqs, amp


def find_peaks_in_spectrum(freqs, amp, threshold=0.1):
    idx, _ = signal.find_peaks(amp, height=threshold)
    return [(float(freqs[i]), float(amp[i])) for i in idx]


# --------------------------------------------------------------------------
# 3. FIR design (by hand) -- windowed sinc
# --------------------------------------------------------------------------
def design_lowpass_fir(num_taps, cutoff_hz, fs, window="hamming"):
    """
    Windowed-sinc FIR low-pass design.

    Ideal low-pass impulse response:  h[n] = 2*fc/fs * sinc(2*fc/fs * (n - M/2))
    truncated to num_taps samples and multiplied by a window to control the
    ripple of the Gibbs phenomenon. num_taps must be odd for a Type-I linear
    phase filter with an integer group delay of (num_taps-1)/2 samples.
    """
    if num_taps % 2 == 0:
        raise ValueError("Use an odd number of taps for an integer group delay.")
    fc_norm = cutoff_hz / fs                 # cycles / sample
    m = (num_taps - 1) / 2.0
    n = np.arange(num_taps)
    h = 2 * fc_norm * np.sinc(2 * fc_norm * (n - m))   # np.sinc = sin(pi x)/(pi x)
    h *= signal.get_window(window, num_taps, fftbins=False)  # symmetric window
    h /= np.sum(h)                           # unity DC gain
    return h


# --------------------------------------------------------------------------
# 4. Filtering
# --------------------------------------------------------------------------
def apply_fir(x, h):
    """
    Direct convolution y[n] = sum_k h[k] x[n-k], 'full' mode, then trim the
    group delay so the output is aligned with the input.
    """
    y_full = np.convolve(x, h, mode="full")
    delay = (len(h) - 1) // 2
    return y_full[delay: delay + len(x)]


# --------------------------------------------------------------------------
# 5. Metrics
# --------------------------------------------------------------------------
def snr_db(reference, estimate):
    """SNR of 'estimate' relative to 'reference' in dB."""
    err = estimate - reference
    return 10 * np.log10(np.sum(reference ** 2) / np.sum(err ** 2))


def main():
    fs = 1000.0
    t, clean, measured = generate_signal(fs)

    # ---- spectrum of the measured signal
    freqs, amp = amplitude_spectrum(measured, fs)
    print("Detected spectral peaks in the measured signal (Hz, amplitude):")
    for f, a in find_peaks_in_spectrum(freqs, amp, threshold=0.2):
        print(f"   {f:7.1f} Hz   A = {a:.3f}")

    # ---- FIR design, by hand vs scipy
    num_taps, cutoff = 101, 80.0
    h_manual = design_lowpass_fir(num_taps, cutoff, fs, window="hamming")
    h_scipy = signal.firwin(num_taps, cutoff, fs=fs, window="hamming")
    max_tap_diff = np.max(np.abs(h_manual - h_scipy))
    print(f"\nFIR design: {num_taps} taps, cutoff {cutoff} Hz, Hamming window")
    print(f"   max |h_manual - h_firwin| = {max_tap_diff:.2e}  (should be ~1e-16)")
    print(f"   group delay = {(num_taps - 1) // 2} samples = "
          f"{(num_taps - 1) / 2 / fs * 1e3:.1f} ms")

    # ---- frequency response, attenuation at the tones
    w, H = signal.freqz(h_manual, worN=4096, fs=fs)
    H_db = 20 * np.log10(np.abs(H) + 1e-12)
    for f0 in (50, 120):
        print(f"   |H({f0} Hz)| = {H_db[np.argmin(np.abs(w - f0))]:6.2f} dB")

    # ---- filtering
    filtered = apply_fir(measured, h_manual)
    # zero-phase alternative (forward-backward) for comparison
    filtered_zp = signal.filtfilt(h_manual, [1.0], measured)

    # ---- SNR, ignoring filter edge transients
    edge = num_taps
    sl = slice(edge, -edge)
    print("\nSNR relative to the clean 50 Hz tone:")
    print(f"   before filtering : {snr_db(clean[sl], measured[sl]):6.2f} dB")
    print(f"   after  FIR       : {snr_db(clean[sl], filtered[sl]):6.2f} dB")
    print(f"   after  filtfilt  : {snr_db(clean[sl], filtered_zp[sl]):6.2f} dB")

    # ---- plots
    fig, ax = plt.subplots(3, 1, figsize=(10, 10))
    ax[0].plot(t, measured, lw=0.8, label="measured (noisy)")
    ax[0].plot(t, clean, lw=1.2, label="clean 50 Hz")
    ax[0].plot(t, filtered, lw=1.2, label="FIR filtered")
    ax[0].set_xlim(0, 0.2)
    ax[0].set_xlabel("time [s]"); ax[0].set_ylabel("amplitude")
    ax[0].set_title("Time domain"); ax[0].legend(loc="upper right")

    f2, amp_f = amplitude_spectrum(filtered, fs)
    ax[1].plot(freqs, amp, label="measured")
    ax[1].plot(f2, amp_f, label="filtered")
    ax[1].set_xlim(0, 250)
    ax[1].set_xlabel("frequency [Hz]"); ax[1].set_ylabel("amplitude")
    ax[1].set_title("Amplitude spectrum (Hann window)"); ax[1].legend()

    ax[2].plot(w, H_db)
    ax[2].axvline(cutoff, color="k", ls="--", lw=0.8, label="cutoff 80 Hz")
    ax[2].set_ylim(-100, 5); ax[2].set_xlim(0, fs / 2)
    ax[2].set_xlabel("frequency [Hz]"); ax[2].set_ylabel("|H| [dB]")
    ax[2].set_title(f"FIR frequency response ({num_taps} taps, Hamming)")
    ax[2].legend()
    fig.tight_layout()
    path = os.path.join(OUT_DIR, "fir_filter.png")
    fig.savefig(path, dpi=130)
    print(f"\nPlot saved to {path}")


if __name__ == "__main__":
    main()
