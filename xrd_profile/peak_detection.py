"""
Peak detection and FWHM estimation for powder XRD data.

Includes simple half-maximum interpolation, Voigt profile fitting,
and reference-guided peak detection with noise cutoff.

The compute_optimal_prominence function is adapted from
crystallite_size_calculator by Dinga Wonanke
(https://github.com/bafgreat/crystallite_size_calculator).
"""

import numpy as np
from scipy.signal import find_peaks
from scipy.special import wofz
from scipy.optimize import curve_fit

from .conversions import two_theta_to_d, fwhm_to_deltaK
from .noise import estimate_noise


def compute_optimal_prominence(intensity, min_prominence=0.01,
                               max_prominence=1.0, tolerance=3, step=0.01):
    """
    Automatically determine the optimal prominence threshold for peak
    detection by iteratively testing values and selecting the one that
    maximises average peak prominence.

    Adapted from crystallite_size_calculator (Wonanke).

    Parameters
    ----------
    intensity : np.ndarray
        Diffraction intensity array.
    min_prominence, max_prominence : float
        Search bounds.
    tolerance : int
        Consecutive iterations without improvement before stopping.
    step : float
        Increment between tested prominence values.

    Returns
    -------
    optimal_prominence : float
    peaks : np.ndarray
        Indices of detected peaks at the optimal prominence.
    """
    max_prominence_avg = 0
    optimal_prominence = min_prominence
    no_improvement_count = 0

    if np.max(intensity) == 0:
        raise ValueError("Maximum intensity is zero; unable to normalise.")

    norm_intensity = intensity / np.max(intensity)

    for current_prominence in np.arange(min_prominence, max_prominence, step):
        peaks, properties = find_peaks(norm_intensity,
                                       prominence=current_prominence)
        prominences = properties["prominences"]
        if len(prominences) > 0:
            current_avg = np.mean(prominences)
            if current_avg > max_prominence_avg:
                max_prominence_avg = current_avg
                optimal_prominence = current_prominence
                no_improvement_count = 0
            else:
                no_improvement_count += 1
        if no_improvement_count >= tolerance:
            break

    peaks, _ = find_peaks(norm_intensity, prominence=optimal_prominence)
    return optimal_prominence, peaks


def estimate_fwhm_simple(two_theta, intensities, height_threshold=0.05):
    """
    Estimate FWHM by direct half-maximum interpolation (no profile fitting).

    Parameters
    ----------
    two_theta : np.ndarray
        2-theta values in degrees.
    intensities : np.ndarray
        Corresponding intensity values.
    height_threshold : float
        Minimum peak height as fraction of maximum intensity.

    Returns
    -------
    fwhm_data : np.ndarray
        FWHM values in degrees 2-theta.
    peak_positions : np.ndarray
        2-theta positions of detected peaks.
    """
    normalised = intensities / np.max(intensities)
    prominence, _ = compute_optimal_prominence(intensities)
    peak_indices, _ = find_peaks(normalised, height=height_threshold,
                                 prominence=prominence)

    fwhms = []
    peak_positions = []

    for peak_index in peak_indices:
        half_max = intensities[peak_index] / 2.0
        left_indices = np.where(intensities[:peak_index] < half_max)[0]
        right_indices = np.where(intensities[peak_index:] < half_max)[0]

        if len(left_indices) == 0 or len(right_indices) == 0:
            continue

        left_idx = left_indices[-1]
        right_idx = right_indices[0] + peak_index
        fwhm = two_theta[right_idx] - two_theta[left_idx]
        fwhms.append(fwhm)
        peak_positions.append(two_theta[peak_index])

    return np.array(fwhms), np.array(peak_positions)


def voigt_profile(x, amplitude, center, sigma, gamma):
    """
    Voigt profile: convolution of Gaussian and Lorentzian components.

    Parameters
    ----------
    x : np.ndarray
        Position values (e.g. 2-theta).
    amplitude : float
        Peak height.
    center : float
        Peak centre position.
    sigma : float
        Gaussian width parameter.
    gamma : float
        Lorentzian width parameter.

    Returns
    -------
    np.ndarray
        Voigt profile evaluated at each x.
    """
    z = ((x - center) + 1j * gamma) / (sigma * np.sqrt(2))
    return amplitude * np.real(wofz(z)) / (sigma * np.sqrt(2 * np.pi))


def estimate_fwhm_voigt(two_theta, intensities, height_threshold=0.01):
    """
    Estimate FWHM by fitting each detected peak with a Voigt profile.

    Parameters
    ----------
    two_theta : np.ndarray
        2-theta values in degrees.
    intensities : np.ndarray
        Corresponding intensity values.
    height_threshold : float
        Minimum peak height as fraction of maximum intensity.

    Returns
    -------
    fwhm_data : np.ndarray
        FWHM values in degrees 2-theta.
    peak_positions : np.ndarray
        2-theta positions of detected peaks.
    sigma_values : np.ndarray
        Fitted Gaussian width for each peak.
    gamma_values : np.ndarray
        Fitted Lorentzian width for each peak.
    """
    prominence, peaks_idx = compute_optimal_prominence(intensities)
    peaks_idx, properties = find_peaks(
        intensities, prominence=prominence,
        height=height_threshold * np.max(intensities)
    )

    fwhm_data = []
    peak_positions = []
    sigma_values = []
    gamma_values = []

    for i, peak in enumerate(peaks_idx):
        left_base = properties["left_bases"][i]
        right_base = properties["right_bases"][i]
        peak_region_x = two_theta[left_base:right_base]
        peak_region_y = intensities[left_base:right_base]

        if len(peak_region_x) < 4:
            continue

        amplitude = peak_region_y.max()
        center = two_theta[peak]

        half_max = amplitude / 2
        close = np.where(
            np.isclose(peak_region_y, half_max, atol=0.1 * half_max)
        )[0]

        if len(close) >= 2:
            estimated_fwhm = peak_region_x[close[-1]] - peak_region_x[close[0]]
        else:
            estimated_fwhm = (peak_region_x[-1] - peak_region_x[0]) / 2

        sigma = estimated_fwhm / (2 * np.sqrt(2 * np.log(2)))
        gamma = estimated_fwhm / 2

        try:
            popt, _ = curve_fit(
                voigt_profile, peak_region_x, peak_region_y,
                p0=[amplitude, center, sigma, gamma],
                method='trf', max_nfev=5000,
                bounds=([0, center - 2, 0, 0],
                        [np.inf, center + 2, np.inf, np.inf])
            )
            fitted_sigma = popt[2]
            fitted_gamma = popt[3]
            fwhm = (0.5346 * (2 * fitted_gamma)
                    + np.sqrt(0.2166 * (2 * fitted_gamma)**2
                              + (2 * fitted_sigma)**2))
        except RuntimeError:
            fwhm = estimated_fwhm
            fitted_sigma = sigma
            fitted_gamma = gamma

        fwhm_data.append(fwhm)
        peak_positions.append(center)
        sigma_values.append(fitted_sigma)
        gamma_values.append(fitted_gamma)

    return (np.array(fwhm_data), np.array(peak_positions),
            np.array(sigma_values), np.array(gamma_values))


def find_peaks_guided(two_theta, intensity, ref_d, wavelength,
                      tolerance_d=0.03, n_sigma=3.0, min_fwhm_steps=3,
                      bg_window=151):
    """
    Reference-guided peak detection with noise cutoff.

    Searches for peaks in the observed pattern near expected reference
    d-spacings. A peak is accepted only if it exceeds the local
    background by n_sigma * noise_floor AND its FWHM exceeds
    min_fwhm_steps * step_size.

    Parameters
    ----------
    two_theta : np.ndarray
        2-theta array (degrees), should be zero-offset corrected.
    intensity : np.ndarray
        Intensity array.
    ref_d : np.ndarray
        Reference d-spacings to search for (angstroms).
    wavelength : float
        X-ray wavelength (angstroms).
    tolerance_d : float
        Maximum allowed d-spacing mismatch (angstroms).
    n_sigma : float
        Minimum peak height above background in noise sigma units.
    min_fwhm_steps : int
        Minimum FWHM in units of the 2-theta step size.
    bg_window : int
        Median filter window for background estimation (points).

    Returns
    -------
    dict with arrays (one entry per accepted peak):
        d_obs, two_theta_obs, fwhm, intensity_obs, d_ref, d_shift, snr
    Plus scalar diagnostics:
        noise_sigma, n_rejected_noise, n_rejected_fwhm, n_rejected_dup
    """
    step_size = np.median(np.diff(two_theta))
    min_fwhm = min_fwhm_steps * step_size

    noise_info = estimate_noise(two_theta, intensity, bg_window=bg_window)
    background = noise_info['background']
    sigma = noise_info['sigma']
    threshold = n_sigma * sigma

    d_pattern = two_theta_to_d(two_theta, wavelength)

    results = {
        'd_obs': [], 'two_theta_obs': [], 'fwhm': [],
        'intensity_obs': [], 'd_ref': [], 'd_shift': [], 'snr': [],
    }
    n_rejected_noise = 0
    n_rejected_fwhm = 0
    n_rejected_dup = 0

    for d_ref_val in ref_d:
        d_lo = d_ref_val - tolerance_d
        d_hi = d_ref_val + tolerance_d
        mask = (d_pattern >= d_lo) & (d_pattern <= d_hi)
        if not np.any(mask):
            continue

        indices = np.where(mask)[0]
        local_i = intensity[indices]
        local_bg = background[indices]

        if len(local_i) < 3:
            continue

        peak_idx_local = np.argmax(local_i)
        peak_intensity = local_i[peak_idx_local]
        peak_bg = local_bg[peak_idx_local]
        global_peak_idx = indices[peak_idx_local]
        peak_tt = two_theta[global_peak_idx]

        peak_height = peak_intensity - peak_bg
        snr = peak_height / sigma if sigma > 0 else 0
        if peak_height < threshold:
            n_rejected_noise += 1
            continue

        half_max = (peak_intensity + peak_bg) / 2.0
        expand = int(2.5 / step_size)
        lo = max(0, global_peak_idx - expand)
        hi = min(len(intensity), global_peak_idx + expand)
        region_tt = two_theta[lo:hi]
        region_i = intensity[lo:hi]

        local_peak_pos = global_peak_idx - lo
        left_below = np.where(region_i[:local_peak_pos] < half_max)[0]
        right_below = np.where(region_i[local_peak_pos:] < half_max)[0]

        if len(left_below) == 0 or len(right_below) == 0:
            n_rejected_fwhm += 1
            continue

        li = left_below[-1]
        ri = right_below[0] + local_peak_pos

        if li + 1 < len(region_i) and region_i[li + 1] != region_i[li]:
            frac_l = ((half_max - region_i[li])
                      / (region_i[li + 1] - region_i[li]))
            tt_left = region_tt[li] + frac_l * (region_tt[li + 1] - region_tt[li])
        else:
            tt_left = region_tt[li]

        if ri - 1 >= 0 and region_i[ri - 1] != region_i[ri]:
            frac_r = ((half_max - region_i[ri])
                      / (region_i[ri - 1] - region_i[ri]))
            tt_right = region_tt[ri] - frac_r * (region_tt[ri] - region_tt[ri - 1])
        else:
            tt_right = region_tt[ri]

        fwhm = tt_right - tt_left

        if fwhm <= min_fwhm or fwhm > 5.0:
            n_rejected_fwhm += 1
            continue

        d_obs_val = two_theta_to_d(peak_tt, wavelength)
        if len(results['two_theta_obs']) > 0:
            if np.min(np.abs(np.array(results['two_theta_obs']) - peak_tt)) < 5 * step_size:
                n_rejected_dup += 1
                continue

        results['d_obs'].append(float(d_obs_val))
        results['two_theta_obs'].append(float(peak_tt))
        results['fwhm'].append(float(fwhm))
        results['intensity_obs'].append(float(peak_intensity))
        results['d_ref'].append(float(d_ref_val))
        results['d_shift'].append(float(d_obs_val - d_ref_val))
        results['snr'].append(float(snr))

    for key in ['d_obs', 'two_theta_obs', 'fwhm', 'intensity_obs',
                'd_ref', 'd_shift', 'snr']:
        results[key] = np.array(results[key])

    results['noise_sigma'] = sigma
    results['n_rejected_noise'] = n_rejected_noise
    results['n_rejected_fwhm'] = n_rejected_fwhm
    results['n_rejected_dup'] = n_rejected_dup

    return results
