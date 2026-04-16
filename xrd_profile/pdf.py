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
from scipy.signal import find_peaks, savgol_filter
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


# ============================================================
# Enhanced PDF: Chebyshev background, sine transform, Lorch
# ============================================================

def chebyshev_background(Q, intensity, order=15, n_iter=20, clip_sigma=1.5):
    """
    Iterative Chebyshev polynomial background estimation in Q-space.

    Fits a Chebyshev polynomial to the intensity data, iteratively
    down-weighting points that lie above the current fit (Bragg peaks)
    so that the polynomial converges to the lower envelope (background
    from capillary, Compton scattering, thermal diffuse scattering, etc.).

    Parameters
    ----------
    Q : np.ndarray
        Scattering vector (inverse angstroms).
    intensity : np.ndarray
        Measured intensity on a uniform Q grid.
    order : int
        Chebyshev polynomial order. Typical values: 10-20 for
        synchrotron data, 6-15 for laboratory data.
    n_iter : int
        Number of iterations for peak clipping.
    clip_sigma : float
        Threshold in sigma units above the fit. Points above
        clip_sigma * noise_std are heavily down-weighted.

    Returns
    -------
    background : np.ndarray
        Estimated background array.
    coeffs : np.ndarray
        Chebyshev polynomial coefficients.
    """
    Q_min, Q_max_val = Q.min(), Q.max()
    Q_norm = 2 * (Q - Q_min) / (Q_max_val - Q_min) - 1

    weights = np.ones_like(intensity)

    for _ in range(n_iter):
        coeffs = np.polynomial.chebyshev.chebfit(
            Q_norm, intensity, order, w=weights)
        bg = np.polynomial.chebyshev.chebval(Q_norm, coeffs)
        residuals = intensity - bg

        # Noise estimate from negative residuals only (below background)
        neg_residuals = residuals[residuals < 0]
        std_neg = np.std(neg_residuals) if len(neg_residuals) > 0 else 0
        if std_neg <= 0:
            std_neg = np.std(residuals) / 2

        # Tiered down-weighting
        weights = np.where(residuals > clip_sigma * std_neg, 0.01, 1.0)
        moderate = (residuals > 0) & (residuals <= clip_sigma * std_neg)
        weights[moderate] = 0.5

    background = np.polynomial.chebyshev.chebval(Q_norm, coeffs)
    return background, coeffs


def compute_pdf_sine(two_theta, intensity, wavelength, cheby_order=15,
                     lorch=True, r_max=30.0, n_r=2000):
    """
    Compute the pair distribution function via sine Fourier transform
    with Chebyshev background subtraction and optional Lorch modification.

    This is the recommended PDF method for quantitative peak profile
    analysis, replacing the FFT-based PDFProcessor for studies where
    background control and termination-ripple suppression matter.

    Parameters
    ----------
    two_theta : np.ndarray
        2-theta array (degrees).
    intensity : np.ndarray
        Intensity array.
    wavelength : float
        X-ray wavelength (angstroms).
    cheby_order : int
        Chebyshev polynomial order for background fitting.
    lorch : bool
        Apply Lorch modification function M(Q) = sinc(Q/Q_max) to
        suppress termination ripples. Recommended for synchrotron
        data with large Q_max.
    r_max : float
        Maximum radial distance (angstroms) for the output PDF.
    n_r : int
        Number of points in the r-grid.

    Returns
    -------
    r : np.ndarray
        Radial distance array (angstroms).
    G_r : np.ndarray
        Pair distribution function values.
    Q_max : float
        Maximum Q reached (inverse angstroms). Determines real-space
        resolution: dr = pi / Q_max.
    """
    # Convert to Q
    theta_rad = np.radians(two_theta / 2)
    Q = (4 * np.pi * np.sin(theta_rad)) / wavelength

    # Interpolate to uniform Q grid
    Q_uniform = np.linspace(Q.min(), Q.max(), len(Q))
    I_uniform = np.interp(Q_uniform, Q, intensity)

    # Chebyshev background subtraction
    bg, _ = chebyshev_background(Q_uniform, I_uniform, order=cheby_order)
    S = I_uniform - bg
    mx = np.max(np.abs(S))
    S = S / mx if mx > 0 else S

    # Build F(Q) = Q * S(Q) [* M(Q)]
    if lorch:
        M_Q = np.sinc(Q_uniform / Q_uniform[-1])
        F_Q = Q_uniform * S * M_Q
    else:
        F_Q = Q_uniform * S

    # Sine Fourier transform: G(r) = (2/pi) * integral{ F(Q) * sin(Q*r) dQ }
    r = np.linspace(0.3, r_max, n_r)
    G_r = np.zeros_like(r)
    for i, ri in enumerate(r):
        G_r[i] = (2 / np.pi) * np.trapezoid(F_Q * np.sin(Q_uniform * ri),
                                          Q_uniform)

    return r, G_r, float(Q_uniform[-1])


def measure_pdf_peaks(r, G_r, min_r=1.0, max_r=20.0, n_max=20):
    """
    Detect and characterise peaks in the pair distribution function.

    Parameters
    ----------
    r : np.ndarray
        Radial distance array (angstroms).
    G_r : np.ndarray
        PDF values.
    min_r, max_r : float
        Radial range to search (angstroms).
    n_max : int
        Maximum number of peaks to return.

    Returns
    -------
    list of dict
        Each dict has keys 'r' (peak position in angstroms),
        'height' (G(r) value), and 'fwhm' (angstroms, or NaN if
        not measurable).
    """
    mask = (r >= min_r) & (r <= max_r)
    rs, gs = r[mask], G_r[mask]

    if len(gs) < 10:
        return []

    # Smooth for peak detection
    win = min(11, len(gs) // 3 * 2 + 1)
    if win >= 5 and len(gs) > 15:
        gs_s = savgol_filter(gs, win, 3)
    else:
        gs_s = gs

    peaks_idx, _ = find_peaks(
        gs_s, prominence=0.005 * np.max(np.abs(gs_s)), distance=5)

    results = []
    for idx in peaks_idx[:n_max]:
        height = gs_s[idx]
        half_max = height / 2 if height > 0 else 0

        # FWHM by half-maximum interpolation
        left = np.where(gs_s[:idx] < half_max)[0]
        right = np.where(gs_s[idx:] < half_max)[0]
        if len(left) > 0 and len(right) > 0:
            fwhm = rs[right[0] + idx] - rs[left[-1]]
        else:
            fwhm = np.nan

        results.append({'r': float(rs[idx]),
                        'height': float(height),
                        'fwhm': float(fwhm)})

    return results


def fit_first_pdf_peak(r, G_r, r_min=1.2, r_max=6.0):
    """
    Fit a Gaussian to the first coordination-shell peak in the PDF.

    Searches for the first positive peak in the range [r_min, r_max],
    fits a Gaussian profile, and returns the peak position and FWHM.

    Parameters
    ----------
    r : np.ndarray
        Radial distance array (angstroms).
    G_r : np.ndarray
        PDF values.
    r_min, r_max : float
        Search range (angstroms). Default 1.2-6.0 covers first-shell
        distances for silicates and oxides.

    Returns
    -------
    r_peak : float
        Fitted peak centre (angstroms), or NaN if fitting fails.
    fwhm : float
        Fitted FWHM (angstroms), or NaN if fitting fails.
    """
    mask = (r >= r_min) & (r <= r_max)
    rs, gs = r[mask], G_r[mask]

    if len(gs) < 10:
        return np.nan, np.nan

    # Smooth and find peaks
    win = min(11, len(gs) // 3 * 2 + 1)
    if win >= 5 and len(gs) > 15:
        gs_s = savgol_filter(gs, win, 3)
    else:
        gs_s = gs

    pks, _ = find_peaks(gs_s, prominence=0.02 * np.max(np.abs(gs)),
                         distance=3)

    for idx in pks:
        amp = gs_s[idx]
        if amp <= 0:
            continue
        r0 = rs[idx]

        # Fit Gaussian in a +/- 0.5 A window around the peak
        fit_mask = (rs >= r0 - 0.5) & (rs <= r0 + 0.5)
        rf, gf = rs[fit_mask], gs[fit_mask]
        if len(rf) < 5:
            continue

        try:
            def gauss(x, a, mu, sig):
                return a * np.exp(-(x - mu)**2 / (2 * sig**2))

            popt, _ = curve_fit(
                gauss, rf, gf, p0=[amp, r0, 0.15],
                bounds=([0, r0 - 0.3, 0.01], [amp * 3, r0 + 0.3, 1.0]))
            return float(popt[1]), float(2.355 * popt[2])
        except (RuntimeError, ValueError):
            continue

    return np.nan, np.nan
