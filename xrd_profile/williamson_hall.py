"""
williamson_hall.py — Williamson-Hall peak profile analysis.

Provides conventional (beta*cos vs 4*sin) and reciprocal-space
(DeltaK vs K) Williamson-Hall analysis for separating crystallite
size and microstrain contributions to XRD peak broadening.

References
----------
Williamson, G. K. & Hall, W. H. (1953). X-ray line broadening from
    filed aluminium and wolfram. Acta Metallurgica 1, 22-31.
Das Bakshi, S., Sinha, D., & Chowdhury, S. G. (2018). Anisotropic
    broadening of XRD peaks of alpha-Fe: Williamson-Hall and
    Warren-Averbach analysis using full width at half maximum (FWHM)
    and integral breadth (IB). Materials Characterization 142, 144-153.

Adapted from crystallite_size_calculator
(https://github.com/bafgreat/crystallite_size_calculator) by Dinga Wonanke.
Modified by M.R.M. Izawa.
"""

import numpy as np
from scipy.stats import linregress

from .conversions import two_theta_to_d, fwhm_to_deltaK, two_theta_to_K
from .noise import estimate_noise, estimate_zero_offset
from .peak_detection import estimate_fwhm_simple, estimate_fwhm_voigt, find_peaks_guided


def williamson_hall(two_theta, intensities, wavelength,
                    use_voigt=False, height_threshold=0.05):
    """
    Williamson-Hall analysis: decompose peak broadening into crystallite
    size and microstrain contributions.

    The W-H equation:
        beta * cos(theta) = (K * lambda / D) + 4 * epsilon * sin(theta)

    Parameters
    ----------
    two_theta : np.ndarray
        Full 2-theta array (degrees).
    intensities : np.ndarray
        Full intensity array.
    wavelength : float
        X-ray wavelength in angstroms (1.5406 for Cu K-alpha).
    use_voigt : bool
        If True, use Voigt profile fitting for FWHM; otherwise use
        simple half-maximum interpolation.
    height_threshold : float
        Minimum peak height as fraction of maximum intensity.

    Returns
    -------
    dict with keys:
        strain : float
            Microstrain (epsilon).
        crystallite_size : float
            Crystallite size D in angstroms.
        n_peaks : int
            Number of peaks used in the fit.
        peak_positions : np.ndarray
            2-theta positions of fitted peaks.
        fwhm : np.ndarray
            FWHM values in degrees.
        r_squared : float
            R-squared of the linear fit.
        slope : float
            Slope of the W-H plot.
        intercept : float
            Y-intercept of the W-H plot.
        x : np.ndarray
            4*sin(theta) values for plotting.
        y : np.ndarray
            beta*cos(theta) values for plotting.
    """
    if use_voigt:
        fwhm_data, peak_positions, _, _ = estimate_fwhm_voigt(
            two_theta, intensities, height_threshold
        )
    else:
        fwhm_data, peak_positions = estimate_fwhm_simple(
            two_theta, intensities, height_threshold
        )

    if len(fwhm_data) < 2:
        return {
            'strain': np.nan, 'crystallite_size': np.nan,
            'n_peaks': len(fwhm_data), 'peak_positions': peak_positions,
            'fwhm': fwhm_data, 'r_squared': np.nan,
            'slope': np.nan, 'intercept': np.nan,
            'x': np.array([]), 'y': np.array([])
        }

    fwhm_radians = np.radians(fwhm_data)
    theta_radians = np.radians(peak_positions / 2)
    cos_theta = np.cos(theta_radians)
    sin_theta = np.sin(theta_radians)

    y = fwhm_radians * cos_theta       # beta * cos(theta)
    x = 4 * sin_theta                   # 4 * sin(theta)

    slope, intercept, r_value, _, _ = linregress(x, y)

    strain = slope  # slope = epsilon (some formulations use slope/4)
    Kfactor = 0.9
    if intercept != 0:
        crystallite_size = (Kfactor * wavelength) / intercept  # in angstroms
    else:
        crystallite_size = np.inf

    # d-spacing for each detected peak
    d_spacings = two_theta_to_d(peak_positions, wavelength)

    return {
        'strain': strain,
        'crystallite_size': crystallite_size,
        'n_peaks': len(fwhm_data),
        'peak_positions': peak_positions,
        'd_spacings': d_spacings,
        'fwhm': fwhm_data,
        'r_squared': r_value**2,
        'slope': slope,
        'intercept': intercept,
        'x': x,
        'y': y,
    }


def williamson_hall_reciprocal(two_theta, intensities, wavelength,
                               use_voigt=False, height_threshold=0.05):
    """
    Williamson-Hall analysis in reciprocal space (DeltaK vs K).

    This wavelength-independent formulation plots DeltaK against K
    where K = 2*sin(theta)/lambda and DeltaK = beta*cos(theta)/lambda.
    The relationship is:

        DeltaK = 1/D + epsilon * K

    or equivalently:

        (DeltaK)^2 = (1/D)^2 + (epsilon * K)^2   (quadratic form)

    The reciprocal-space formulation enables direct comparison of data
    collected at different wavelengths (e.g., lab Cu K-alpha vs synchrotron).

    Parameters
    ----------
    two_theta : np.ndarray
        Full 2-theta array (degrees).
    intensities : np.ndarray
        Full intensity array.
    wavelength : float
        X-ray wavelength in angstroms.
    use_voigt : bool
        If True, use Voigt profile fitting for FWHM.
    height_threshold : float
        Minimum peak height as fraction of maximum intensity.

    Returns
    -------
    dict with keys:
        strain : float
            Microstrain (epsilon).
        crystallite_size : float
            Crystallite size D in angstroms.
        n_peaks : int
            Number of peaks used.
        peak_positions : np.ndarray
            2-theta positions (degrees).
        d_spacings : np.ndarray
            d-spacing for each peak (angstroms).
        fwhm : np.ndarray
            FWHM values in degrees 2-theta.
        r_squared : float
            R-squared of the linear fit.
        slope : float
            Slope of DeltaK vs K (= epsilon).
        intercept : float
            Y-intercept (= 1/D).
        K : np.ndarray
            Reciprocal variable for plotting.
        deltaK : np.ndarray
            Peak broadening in reciprocal space for plotting.
    """
    if use_voigt:
        fwhm_data, peak_positions, _, _ = estimate_fwhm_voigt(
            two_theta, intensities, height_threshold
        )
    else:
        fwhm_data, peak_positions = estimate_fwhm_simple(
            two_theta, intensities, height_threshold
        )

    nan_result = {
        'strain': np.nan, 'crystallite_size': np.nan,
        'n_peaks': len(fwhm_data), 'peak_positions': peak_positions,
        'd_spacings': two_theta_to_d(peak_positions, wavelength) if len(peak_positions) > 0 else np.array([]),
        'fwhm': fwhm_data, 'r_squared': np.nan,
        'slope': np.nan, 'intercept': np.nan,
        'K': np.array([]), 'deltaK': np.array([])
    }

    if len(fwhm_data) < 2:
        return nan_result

    K_vals = two_theta_to_K(peak_positions, wavelength)
    deltaK_vals = fwhm_to_deltaK(fwhm_data, peak_positions, wavelength)
    d_spacings = two_theta_to_d(peak_positions, wavelength)

    slope, intercept, r_value, _, _ = linregress(K_vals, deltaK_vals)

    strain = slope
    if intercept > 0:
        crystallite_size = 1.0 / intercept  # angstroms
    elif intercept < 0:
        crystallite_size = 1.0 / intercept  # negative = physically invalid
    else:
        crystallite_size = np.inf

    return {
        'strain': strain,
        'crystallite_size': crystallite_size,
        'n_peaks': len(fwhm_data),
        'peak_positions': peak_positions,
        'd_spacings': d_spacings,
        'fwhm': fwhm_data,
        'r_squared': r_value**2,
        'slope': slope,
        'intercept': intercept,
        'K': K_vals,
        'deltaK': deltaK_vals,
    }


def guided_williamson_hall(two_theta, intensity, ref_d, wavelength,
                           tolerance_d=0.03, n_sigma=3.0,
                           min_fwhm_steps=3, correct_offset=True):
    """
    Reference-guided Williamson-Hall analysis in reciprocal space.

    Combines zero-point offset correction, reference-guided peak detection
    with noise cutoff, and DeltaK-vs-K linear regression.

    Parameters
    ----------
    two_theta : np.ndarray
        2-theta array (degrees).
    intensity : np.ndarray
        Intensity array.
    ref_d : np.ndarray
        Reference mineral d-spacings (angstroms), sorted by
        decreasing intensity.
    wavelength : float
        X-ray wavelength (angstroms).
    tolerance_d : float
        d-spacing matching tolerance (angstroms).
    n_sigma : float
        Noise cutoff in units of sigma.
    min_fwhm_steps : int
        Minimum FWHM in step-size units.
    correct_offset : bool
        If True, estimate and apply zero-point offset correction.

    Returns
    -------
    dict with keys:
        strain : float — microstrain (slope of DeltaK vs K)
        crystallite_size : float — size D in angstroms (1/intercept)
        r_squared : float — R-squared of linear fit
        n_peaks : int — number of peaks used
        zero_offset : float — applied zero-point correction (degrees)
        peaks : dict — full output from find_peaks_guided
        K : np.ndarray — reciprocal variable for each peak
        deltaK : np.ndarray — broadening in reciprocal space
        slope : float
        intercept : float
    """
    # Zero-point offset
    if correct_offset:
        offset = estimate_zero_offset(two_theta, intensity, ref_d, wavelength)
        tt_corr = two_theta + offset
    else:
        offset = 0.0
        tt_corr = two_theta

    # Guided peak detection
    peaks = find_peaks_guided(tt_corr, intensity, ref_d, wavelength,
                              tolerance_d=tolerance_d, n_sigma=n_sigma,
                              min_fwhm_steps=min_fwhm_steps)

    result = {
        'strain': np.nan, 'crystallite_size': np.nan,
        'r_squared': np.nan, 'n_peaks': len(peaks['d_obs']),
        'zero_offset': offset, 'peaks': peaks,
        'K': np.array([]), 'deltaK': np.array([]),
        'slope': np.nan, 'intercept': np.nan,
    }

    if len(peaks['d_obs']) < 3:
        return result

    K_vals = 1.0 / peaks['d_obs']
    deltaK_vals = fwhm_to_deltaK(peaks['fwhm'], peaks['two_theta_obs'],
                                  wavelength)

    slope, intercept, r_value, _, _ = linregress(K_vals, deltaK_vals)

    result['strain'] = slope
    result['crystallite_size'] = (1.0 / intercept) if intercept > 0 else np.nan
    result['r_squared'] = r_value ** 2
    result['slope'] = slope
    result['intercept'] = intercept
    result['K'] = K_vals
    result['deltaK'] = deltaK_vals

    return result
