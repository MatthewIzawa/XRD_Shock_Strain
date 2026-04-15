"""
pdf.py — Pair distribution function (PDF) analysis from powder XRD data.

Computes the pair distribution function g(r) from experimental
2-theta/intensity data via Fourier transform of the structure factor
S(Q), and estimates crystallite size by fitting envelope functions
to the PDF maxima.

References
----------
Gesing, T. M. & Robben, L. (2024). Crystallite size determination
    from pair distribution function analysis using the envelope
    approach. J. Appl. Cryst. 57, 1466-1476.

Adapted from crystallite_size_calculator
(https://github.com/bafgreat/crystallite_size_calculator) by Dinga Wonanke.
Modified by M.R.M. Izawa.
"""

import numpy as np
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
from scipy.fft import fft


class PDFProcessor:
    """
    Compute the pair distribution function (PDF) from powder XRD data
    via Fourier transform of the structure factor S(Q).

    Parameters
    ----------
    two_theta : np.ndarray
        2-theta diffraction angles in degrees.
    intensity : np.ndarray
        Measured intensities.
    wavelength : float
        X-ray wavelength in angstroms.
    """

    def __init__(self, two_theta, intensity, wavelength):
        self.two_theta = two_theta
        self.intensity = intensity
        self.wavelength = wavelength

    def compute_q(self):
        """Convert 2-theta to scattering vector magnitude Q (inverse angstroms)."""
        theta = np.radians(self.two_theta / 2)
        return (4 * np.pi / self.wavelength) * np.sin(theta)

    def interpolate_data(self, q):
        """Interpolate to uniform Q spacing for FFT."""
        q_uniform = np.linspace(q.min(), q.max(), len(q))
        interp_func = interp1d(q, self.intensity, kind='linear',
                               fill_value="extrapolate")
        return q_uniform, interp_func(q_uniform)

    def compute_structure_factor(self, intensity_uniform):
        """Normalised structure factor S(Q)."""
        background = np.min(intensity_uniform)
        s_q = ((intensity_uniform - background)
               / np.max(intensity_uniform - background))
        return s_q

    def fourier_transform(self, q_uniform, s_q):
        """FFT of Q*S(Q) to obtain g(r)."""
        delta_q = q_uniform[1] - q_uniform[0]
        g_r = fft(s_q * q_uniform) * delta_q
        r = np.fft.fftfreq(len(q_uniform), delta_q)
        return r, g_r.real

    def compute_pdf(self):
        """
        Full pipeline: 2-theta/intensity -> Q -> S(Q) -> g(r).

        Returns
        -------
        r : np.ndarray
            Radial distance array (angstroms).
        g_r : np.ndarray
            Pair distribution function values.
        """
        q = self.compute_q()
        q_uniform, intensity_uniform = self.interpolate_data(q)
        s_q = self.compute_structure_factor(intensity_uniform)
        r, g_r = self.fourier_transform(q_uniform, s_q)
        return r, g_r


# ============================================================
# Envelope function approach for crystallite size from PDF
# ============================================================

def _remove_nan(g_r, radii):
    """Remove NaN entries from paired arrays."""
    if len(g_r) != len(radii):
        raise ValueError("g_r and radii must have the same length.")
    valid = ~np.isnan(g_r)
    return g_r[valid], radii[valid]


def _extract_maxima(g_r, radii, num_intervals=5):
    """Extract one maximum per interval from g(r)."""
    g_r, radii = _remove_nan(g_r, radii)
    interval_length = len(radii) // num_intervals
    max_r = []
    max_g = []
    for i in range(num_intervals):
        start = i * interval_length
        end = (i + 1) * interval_length
        r_int = radii[start:end]
        g_int = g_r[start:end]
        if len(g_int) == 0:
            continue
        idx = np.argmax(g_int)
        max_r.append(r_int[idx])
        max_g.append(g_int[idx])
    return np.array(max_r), np.array(max_g)


def _adaptive_maxima_extraction(g_r, r, min_interval=3,
                                max_interval=50, tolerance=3):
    """Adaptively determine the optimal interval size for maxima extraction."""
    best_interval = min_interval
    max_prom = 0
    no_improvement = 0

    for interval in range(min_interval, max_interval + 1):
        r_max, g_max = _extract_maxima(g_r, r, num_intervals=interval)
        if len(g_max) < 3:
            continue
        peaks, props = find_peaks(g_max, prominence=0.1)
        if len(props["prominences"]) > 0:
            avg_prom = np.mean(props["prominences"])
            if avg_prom > max_prom:
                max_prom = avg_prom
                best_interval = interval
                no_improvement = 0
            else:
                no_improvement += 1
        if no_improvement >= tolerance:
            break

    g_norm = g_r / np.max(np.abs(g_r)) if np.max(np.abs(g_r)) > 0 else g_r
    num_intervals = max(1, len(g_norm) // best_interval)
    return _extract_maxima(g_r, r, num_intervals=num_intervals)


def spherical_envelope(r, d_crys):
    """Spherical crystallite envelope function."""
    return np.heaviside(d_crys - r, 0.5) * (
        1 - (3/2) * (r / d_crys) + (1/2) * (r / d_crys)**3
    )


def cubic_envelope(r, d_crys):
    """Cubic crystallite envelope function."""
    return np.heaviside(d_crys - r, 0.5) * (
        1 - 2 * (r / d_crys) + (r / d_crys)**2
    )


def plate_envelope(r, d_crys):
    """Plate-like crystallite envelope function."""
    return np.heaviside(d_crys - r, 0.5) * (1 - r / d_crys)


def crystallite_size_from_pdf(r, g_r, shape='spherical'):
    """
    Estimate crystallite size by fitting an envelope function to the
    maxima of the pair distribution function g(r).

    Parameters
    ----------
    r : np.ndarray
        Radial distances from PDF.
    g_r : np.ndarray
        PDF values.
    shape : str
        Crystallite shape model: 'spherical', 'cubic', or 'plate'.

    Returns
    -------
    d_crys : float
        Estimated crystallite size in the same units as r.
    """
    envelope_funcs = {
        'spherical': spherical_envelope,
        'cubic': cubic_envelope,
        'plate': plate_envelope,
    }
    if shape not in envelope_funcs:
        raise ValueError(f"Shape must be one of {list(envelope_funcs.keys())}")

    g_norm = g_r / np.max(np.abs(g_r)) if np.max(np.abs(g_r)) > 0 else g_r
    r_max, g_max = _adaptive_maxima_extraction(g_norm, r)

    if len(r_max) < 2:
        return np.nan

    try:
        popt, _ = curve_fit(envelope_funcs[shape], r_max, g_max,
                            bounds=([0], [np.inf]))
        return popt[0]
    except (RuntimeError, ValueError):
        return np.nan
