"""
Noise characterisation and zero-point offset estimation for XRD patterns.
"""

import numpy as np
from scipy.optimize import minimize_scalar

from .conversions import d_to_two_theta


def estimate_noise(two_theta, intensity, bg_window=151, quiet_window=100):
    """
    Estimate the noise floor of an XRD pattern.

    Uses a large-window median filter to model the background, then
    measures the standard deviation of residuals in rolling windows.
    The minimum rolling sigma (the quietest region) defines the noise
    floor, avoiding contamination from actual diffraction peaks.

    Parameters
    ----------
    two_theta : np.ndarray
        2-theta array (degrees).
    intensity : np.ndarray
        Intensity array.
    bg_window : int
        Median filter window for background estimation (points).
    quiet_window : int
        Rolling window size for sigma estimation (points).

    Returns
    -------
    dict with keys:
        background : np.ndarray
            Estimated background curve.
        sigma : float
            Noise floor (counts).
        residuals : np.ndarray
            intensity - background.
    """
    from scipy.ndimage import median_filter

    background = median_filter(intensity, size=bg_window)
    residuals = intensity - background

    rolling_std = np.array([
        np.std(residuals[i:i + quiet_window])
        for i in range(0, len(residuals) - quiet_window, quiet_window // 2)
    ])
    sigma = np.min(rolling_std) if len(rolling_std) > 0 else np.std(residuals)

    return {
        'background': background,
        'sigma': sigma,
        'residuals': residuals,
    }


def estimate_zero_offset(two_theta, intensity, ref_d, wavelength,
                         search_range=0.3, n_strong=10):
    """
    Estimate 2-theta zero-point offset by maximising the summed intensity
    at the strongest reference peak positions over a range of trial offsets.

    Parameters
    ----------
    two_theta : np.ndarray
        Observed 2-theta array (degrees).
    intensity : np.ndarray
        Observed intensity array.
    ref_d : np.ndarray
        Reference d-spacings (sorted by decreasing intensity).
    wavelength : float
        X-ray wavelength (angstroms).
    search_range : float
        Maximum offset to search (degrees, symmetric about zero).
    n_strong : int
        Number of strongest reference peaks to use.

    Returns
    -------
    offset : float
        Estimated zero-point offset (degrees). Apply as:
        two_theta_corrected = two_theta + offset.
    """
    ref_tt = d_to_two_theta(ref_d[:n_strong], wavelength)
    ref_tt = ref_tt[~np.isnan(ref_tt)]
    ref_tt = ref_tt[(ref_tt > two_theta.min() + 1) &
                    (ref_tt < two_theta.max() - 1)]

    if len(ref_tt) == 0:
        return 0.0

    def neg_score(offset):
        corrected = two_theta + offset
        score = 0.0
        for rtt in ref_tt:
            idx = np.argmin(np.abs(corrected - rtt))
            lo = max(0, idx - 25)
            hi = min(len(intensity), idx + 25)
            score -= np.max(intensity[lo:hi])
        return score

    result = minimize_scalar(neg_score,
                             bounds=(-search_range, search_range),
                             method='bounded')
    return result.x
