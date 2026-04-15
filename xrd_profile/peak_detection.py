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


# ============================================================
# Cross-phase overlap detection and peak quality scoring
# ============================================================

def _measure_peak_asymmetry(two_theta, intensity, peak_tt, fwhm,
                            expand_fwhm=2.5):
    """
    Measure the asymmetry of a single peak as the ratio of left/right
    integrated area above local background.

    Returns
    -------
    asymmetry : float
        0.0 = perfectly symmetric, 1.0 = completely one-sided.
    """
    step = np.median(np.diff(two_theta))
    half_width = int(expand_fwhm * fwhm / step)
    centre_idx = np.argmin(np.abs(two_theta - peak_tt))
    lo = max(0, centre_idx - half_width)
    hi = min(len(intensity), centre_idx + half_width + 1)

    region = intensity[lo:hi].astype(float)
    if len(region) < 5:
        return 0.5

    # Local background from edges
    bg = np.mean([np.mean(region[:3]), np.mean(region[-3:])])
    corrected = np.maximum(region - bg, 0)

    local_peak = centre_idx - lo
    if local_peak <= 0 or local_peak >= len(corrected) - 1:
        return 0.5

    left_area = np.sum(corrected[:local_peak])
    right_area = np.sum(corrected[local_peak + 1:])
    total = left_area + right_area
    if total <= 0:
        return 0.5

    return abs(left_area - right_area) / total


def check_cross_phase_overlap(peak_results, other_phase_d, wavelength,
                              overlap_tol_deg=0.15):
    """
    Flag peaks that overlap with peaks from other mineral phases.

    Parameters
    ----------
    peak_results : dict
        Output from find_peaks_guided.
    other_phase_d : list of np.ndarray
        d-spacing arrays for each interfering phase
        (e.g., [pyroxene_d, olivine_d]).
    wavelength : float
        X-ray wavelength (angstroms).
    overlap_tol_deg : float
        Minimum angular separation (degrees 2-theta) required.
        Peaks closer than max(overlap_tol_deg, 1.5 * own_fwhm)
        to any other-phase peak are flagged.

    Returns
    -------
    dict
        Copy of peak_results with added arrays:
        - 'cross_phase_overlap' : bool array
        - 'nearest_other_phase_deg' : float array
        - 'n_rejected_overlap' : int
    """
    from .conversions import d_to_two_theta

    # Convert all other-phase d-spacings to 2-theta
    other_tt_all = []
    for d_arr in other_phase_d:
        tt_arr = d_to_two_theta(np.asarray(d_arr), wavelength)
        tt_arr = tt_arr[~np.isnan(tt_arr)]
        other_tt_all.append(tt_arr)
    other_tt_combined = np.concatenate(other_tt_all) if other_tt_all else np.array([])

    n_peaks = len(peak_results['two_theta_obs'])
    overlap = np.zeros(n_peaks, dtype=bool)
    nearest_deg = np.full(n_peaks, np.inf)

    if len(other_tt_combined) > 0:
        for i in range(n_peaks):
            pk_tt = peak_results['two_theta_obs'][i]
            pk_fwhm = peak_results['fwhm'][i]
            dists = np.abs(other_tt_combined - pk_tt)
            min_dist = np.min(dists)
            nearest_deg[i] = min_dist

            effective_tol = max(overlap_tol_deg, 1.5 * pk_fwhm)
            if min_dist < effective_tol:
                overlap[i] = True

    result = dict(peak_results)
    result['cross_phase_overlap'] = overlap
    result['nearest_other_phase_deg'] = nearest_deg
    result['n_rejected_overlap'] = int(np.sum(overlap))

    return result


def score_peak_quality(peak_results, two_theta, intensity, wavelength,
                       tolerance_d=0.03,
                       weights=None):
    """
    Compute a per-peak composite quality score.

    Combines signal-to-noise ratio, isolation from neighbours,
    peak symmetry, and d-spacing match accuracy into a single
    score in [0, 1] for each detected peak.

    Parameters
    ----------
    peak_results : dict
        Output from find_peaks_guided (optionally with
        cross_phase_overlap from check_cross_phase_overlap).
    two_theta : np.ndarray
        Full 2-theta array.
    intensity : np.ndarray
        Full intensity array.
    wavelength : float
        X-ray wavelength (angstroms).
    tolerance_d : float
        d-spacing tolerance used in peak detection (angstroms).
    weights : dict or None
        Relative weights for sub-scores. Defaults:
        {'snr': 0.35, 'isolation': 0.25, 'symmetry': 0.20,
         'd_match': 0.20}

    Returns
    -------
    dict
        Copy of peak_results with added arrays:
        - 'quality_score' : float array in [0, 1]
        - 'snr_score', 'isolation_score', 'symmetry_score',
          'd_match_score' : float arrays in [0, 1]
    """
    if weights is None:
        weights = {'snr': 0.35, 'isolation': 0.25,
                   'symmetry': 0.20, 'd_match': 0.20}

    n = len(peak_results['two_theta_obs'])
    snr_scores = np.zeros(n)
    isolation_scores = np.zeros(n)
    symmetry_scores = np.zeros(n)
    d_match_scores = np.zeros(n)

    snr_ref = 20.0  # reference "good" SNR

    for i in range(n):
        # SNR score
        snr_scores[i] = min(1.0, peak_results['snr'][i] / snr_ref)

        # Isolation score
        pk_tt = peak_results['two_theta_obs'][i]
        pk_fwhm = peak_results['fwhm'][i]

        # Distance to nearest same-phase peak
        other_tt = np.delete(peak_results['two_theta_obs'], i)
        if len(other_tt) > 0:
            min_same = np.min(np.abs(other_tt - pk_tt))
        else:
            min_same = np.inf

        # Distance to nearest other-phase peak (if available)
        if 'nearest_other_phase_deg' in peak_results:
            min_other = peak_results['nearest_other_phase_deg'][i]
            min_neighbor = min(min_same, min_other)
        else:
            min_neighbor = min_same

        isolation_scores[i] = min(1.0, min_neighbor / (3 * pk_fwhm)) \
            if pk_fwhm > 0 else 0.0

        # Symmetry score
        asym = _measure_peak_asymmetry(two_theta, intensity, pk_tt, pk_fwhm)
        symmetry_scores[i] = 1.0 - asym

        # d-spacing match score
        d_shift = abs(peak_results['d_shift'][i])
        d_match_scores[i] = max(0.0, 1.0 - d_shift / tolerance_d) \
            if tolerance_d > 0 else 1.0

    # Composite
    composite = (weights['snr'] * snr_scores
                 + weights['isolation'] * isolation_scores
                 + weights['symmetry'] * symmetry_scores
                 + weights['d_match'] * d_match_scores)
    # Normalise to [0, 1]
    w_total = sum(weights.values())
    if w_total > 0:
        composite = composite / w_total

    # Zero out overlapping peaks
    if 'cross_phase_overlap' in peak_results:
        composite[peak_results['cross_phase_overlap']] = 0.0

    result = dict(peak_results)
    result['quality_score'] = composite
    result['snr_score'] = snr_scores
    result['isolation_score'] = isolation_scores
    result['symmetry_score'] = symmetry_scores
    result['d_match_score'] = d_match_scores

    return result
