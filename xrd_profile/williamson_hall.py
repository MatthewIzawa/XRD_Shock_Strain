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
from .peak_detection import (estimate_fwhm_simple, estimate_fwhm_voigt,
                             find_peaks_guided, check_cross_phase_overlap,
                             score_peak_quality)


def _apply_caglioti_correction(fwhm_deg, two_theta_deg, inst_profile):
    """Caglioti-subtract instrumental FWHM from observed FWHM.

    Uses Gaussian-quadrature combination:
        beta_corr^2 = max(beta_obs^2 - beta_inst^2, eps)

    Peaks where beta_obs <= beta_inst are flagged: their corrected FWHM
    is replaced with NaN (caller is expected to filter).

    Parameters
    ----------
    fwhm_deg : np.ndarray
        Observed sample FWHMs (degrees).
    two_theta_deg : np.ndarray
        Peak 2-theta positions (degrees).
    inst_profile : InstrumentalProfile
        Instrumental Caglioti profile.

    Returns
    -------
    fwhm_corr : np.ndarray
        Corrected FWHMs. NaN where the peak was over-corrected.
    """
    fwhm_obs = np.asarray(fwhm_deg, dtype=float)
    tt = np.asarray(two_theta_deg, dtype=float)
    fwhm_inst = np.array([inst_profile.fwhm_at(t) for t in tt])
    diff_sq = fwhm_obs**2 - fwhm_inst**2
    # The np.maximum(diff_sq, 0) inside the where guards against
    # floating-point jitter where diff_sq is ~0+epsilon: sqrt would
    # otherwise occasionally produce NaN on a value mathematically
    # identical to zero.
    fwhm_corr = np.where(diff_sq > 0, np.sqrt(np.maximum(diff_sq, 0)),
                          np.nan)
    return fwhm_corr


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
                           min_fwhm_steps=3, correct_offset=True,
                           other_phase_d=None, other_phase_names=None,
                           overlap_tol_deg=0.15,
                           min_quality=0.3,
                           quality_weights=None,
                           weighted_regression=True,
                           min_peaks_reliable=5,
                           min_r2_reliable=0.3,
                           min_r2_marginal=0.05,
                           sample_flags=None,
                           export_path=None,
                           inst_profile=None):
    """
    Reference-guided Williamson-Hall analysis in reciprocal space
    with cross-phase overlap rejection, peak quality scoring,
    weighted regression, and reliability classification.

    Parameters
    ----------
    two_theta, intensity, ref_d, wavelength, tolerance_d, n_sigma,
    min_fwhm_steps, correct_offset :
        Same as before (see find_peaks_guided).
    other_phase_d : list of np.ndarray or None
        d-spacings for each interfering phase. Peaks overlapping
        with these are excluded.
    other_phase_names : list of str or None
        Names for the interfering phases (for reporting).
    overlap_tol_deg : float
        Minimum angular separation from other-phase peaks (degrees).
    min_quality : float
        Minimum composite quality score [0, 1] for inclusion.
    quality_weights : dict or None
        Weights for quality sub-scores.
    weighted_regression : bool
        Weight the regression by peak quality scores.
    min_peaks_reliable : int
        Minimum peaks for 'reliable' classification.
    min_r2_reliable : float
        Minimum R-squared for 'reliable'.
    min_r2_marginal : float
        Minimum R-squared for 'marginal'.
    sample_flags : dict or None
        Keys: 'multiple_plagioclase_populations' (bool),
        'maskelynite_present' (bool), 'maskelynite_fraction' (float).
    export_path : str or None
        Path for CSV validation export.

    Returns
    -------
    dict with all original keys plus:
        n_peaks_total, n_peaks_used, n_rejected_overlap,
        n_rejected_quality, reliability, reliability_reasons,
        warnings, peak_quality, regression_weights, residuals.
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

    n_total = len(peaks['d_obs'])

    # Cross-phase overlap rejection
    n_rejected_overlap = 0
    if other_phase_d is not None and len(other_phase_d) > 0:
        peaks = check_cross_phase_overlap(peaks, other_phase_d, wavelength,
                                          overlap_tol_deg=overlap_tol_deg)
        n_rejected_overlap = peaks['n_rejected_overlap']

    # Peak quality scoring
    peaks = score_peak_quality(peaks, tt_corr, intensity, wavelength,
                               tolerance_d=tolerance_d,
                               weights=quality_weights)

    # Build inclusion mask
    quality = peaks['quality_score']
    include = quality >= min_quality
    if 'cross_phase_overlap' in peaks:
        include = include & ~peaks['cross_phase_overlap']

    n_rejected_quality = int(np.sum((quality < min_quality)
                                    & (~peaks.get('cross_phase_overlap',
                                                  np.zeros(n_total, dtype=bool)))))
    n_used = int(np.sum(include))

    # Prepare result template
    warnings_list = []
    reliability_reasons = []

    result = {
        'strain': np.nan, 'crystallite_size': np.nan,
        'r_squared': np.nan,
        'n_peaks': n_used,
        'n_peaks_total': n_total,
        'n_peaks_used': n_used,
        'n_rejected_overlap': n_rejected_overlap,
        'n_rejected_quality': n_rejected_quality,
        'zero_offset': offset,
        'peaks': peaks,
        'K': np.array([]), 'deltaK': np.array([]),
        'slope': np.nan, 'intercept': np.nan,
        'reliability': 'unreliable',
        'reliability_reasons': [],
        'warnings': [],
        'peak_quality': quality,
        'regression_weights': np.array([]),
        'residuals': np.array([]),
    }

    # Sample suitability warnings
    if sample_flags:
        if sample_flags.get('multiple_plagioclase_populations', False):
            warnings_list.append(
                "Multiple plagioclase populations detected. W-H assumes "
                "a single population; strain/size values represent a "
                "convolution of distinct histories and may not be "
                "physically meaningful.")
        if sample_flags.get('maskelynite_present', False):
            warnings_list.append(
                "Maskelynite (amorphous plagioclase) is present. Broad "
                "amorphous scattering may contaminate FWHM measurements "
                "of remaining crystalline peaks.")
            frac = sample_flags.get('maskelynite_fraction', 0)
            if frac > 0.5:
                warnings_list.append(
                    f"Maskelynite fraction ({frac:.0%}) exceeds 50%. "
                    "W-H of remaining crystalline peaks is likely not "
                    "representative of the original material.")

    result['warnings'] = warnings_list

    if n_used < 3:
        reliability_reasons.append(
            f"Only {n_used} peaks after filtering (need >= 3 for fit)")
        result['reliability_reasons'] = reliability_reasons
        _export_wh_csv(export_path, peaks, include, result, wavelength)
        return result

    # NEW v0.4.0: Caglioti instrumental correction. Subtract the
    # instrumental FWHM from each peak's measured FWHM in quadrature.
    # Peaks where beta_obs <= beta_inst are flagged in warnings and
    # removed from the analysis.
    if inst_profile is not None:
        fwhm_corr = _apply_caglioti_correction(
            peaks['fwhm'], peaks['two_theta_obs'], inst_profile)
        keep = ~np.isnan(fwhm_corr)
        n_dropped = int(np.sum(~keep))
        if n_dropped > 0:
            warnings_list.append(
                f'{n_dropped} peak(s) over-corrected by instrumental '
                f'subtraction (beta_obs <= beta_inst); excluded from fit.')
            result['warnings'] = warnings_list
        # Filter every parallel array on `peaks` plus the local
        # `quality` and `include` arrays, keeping them aligned.
        peaks = {k: (v[keep] if isinstance(v, np.ndarray) and v.shape
                     and v.shape[0] == len(keep) else v)
                 for k, v in peaks.items()}
        peaks['fwhm'] = fwhm_corr[keep]
        quality = quality[keep]
        include = include[keep]
        n_used = int(np.sum(include))
        result['peaks'] = peaks
        result['peak_quality'] = quality
        result['n_peaks_used'] = n_used
        result['n_peaks'] = n_used          # keep alias in sync
        if n_used < 3:
            reliability_reasons.append(
                f'Only {n_used} peaks after instrumental correction '
                f'(need >= 3 for fit)')
            result['reliability_reasons'] = reliability_reasons
            _export_wh_csv(export_path, peaks, include, result, wavelength)
            return result

    # Extract included peaks
    K_all = 1.0 / peaks['d_obs']
    dK_all = fwhm_to_deltaK(peaks['fwhm'], peaks['two_theta_obs'],
                             wavelength)

    K_vals = K_all[include]
    deltaK_vals = dK_all[include]
    w = quality[include]

    # Regression
    if weighted_regression and np.sum(w) > 0:
        coeffs = np.polyfit(K_vals, deltaK_vals, 1, w=w)
        slope, intercept = coeffs[0], coeffs[1]

        y_pred = slope * K_vals + intercept
        w_mean = np.average(deltaK_vals, weights=w)
        ss_res = np.sum(w * (deltaK_vals - y_pred) ** 2)
        ss_tot = np.sum(w * (deltaK_vals - w_mean) ** 2)
        r_squared = 1 - ss_res / ss_tot if ss_tot > 0 else 0.0
    else:
        slope, intercept, r_value, _, _ = linregress(K_vals, deltaK_vals)
        r_squared = r_value ** 2
        y_pred = slope * K_vals + intercept

    residuals = deltaK_vals - y_pred

    strain = slope
    crystallite_size = (1.0 / intercept) if intercept > 0 else np.nan

    # Reliability classification
    is_reliable = True

    if r_squared < min_r2_reliable:
        is_reliable = False
        reliability_reasons.append(
            f"R-squared = {r_squared:.3f} < {min_r2_reliable} threshold")

    if n_used < min_peaks_reliable:
        is_reliable = False
        reliability_reasons.append(
            f"Only {n_used} peaks used (need {min_peaks_reliable} "
            f"for reliable)")

    if slope < 0:
        is_reliable = False
        reliability_reasons.append(
            f"Negative slope ({slope:.5f}) implies negative strain")

    if intercept <= 0:
        is_reliable = False
        reliability_reasons.append(
            f"Non-positive intercept ({intercept:.5f}) implies "
            "unphysical crystallite size")

    if is_reliable:
        reliability = 'reliable'
    elif r_squared >= min_r2_marginal and n_used >= 3:
        reliability = 'marginal'
    else:
        reliability = 'unreliable'

    if warnings_list:
        if reliability == 'reliable':
            reliability = 'marginal'
            reliability_reasons.append(
                "Downgraded due to sample suitability warnings")

    result.update({
        'strain': strain,
        'crystallite_size': crystallite_size,
        'r_squared': r_squared,
        'slope': slope,
        'intercept': intercept,
        'K': K_vals,
        'deltaK': deltaK_vals,
        'reliability': reliability,
        'reliability_reasons': reliability_reasons,
        'regression_weights': w,
        'residuals': residuals,
    })

    _export_wh_csv(export_path, peaks, include, result, wavelength)

    return result


def _export_wh_csv(export_path, peaks, include_mask, result, wavelength):
    """Write per-peak data to CSV for manual validation."""
    if export_path is None:
        return

    n = len(peaks['d_obs'])
    if n == 0:
        return

    K_all = 1.0 / peaks['d_obs']
    dK_all = fwhm_to_deltaK(peaks['fwhm'], peaks['two_theta_obs'],
                             wavelength)

    with open(export_path, 'w') as f:
        f.write(f"# Williamson-Hall peak validation export\n")
        f.write(f"# Wavelength: {wavelength} angstroms\n")
        f.write(f"# Zero offset: {result['zero_offset']:.4f} deg\n")
        f.write(f"# Reliability: {result['reliability']}\n")
        for reason in result.get('reliability_reasons', []):
            f.write(f"# Reason: {reason}\n")
        for warning in result.get('warnings', []):
            f.write(f"# WARNING: {warning}\n")
        f.write(f"# Peaks total: {result['n_peaks_total']}, "
                f"used: {result['n_peaks_used']}, "
                f"rejected overlap: {result['n_rejected_overlap']}, "
                f"rejected quality: {result['n_rejected_quality']}\n")
        f.write("two_theta_obs,d_obs,d_ref,d_shift,fwhm,snr,"
                "quality_score,snr_score,isolation_score,"
                "symmetry_score,d_match_score,"
                "cross_phase_overlap,included_in_fit,K,deltaK\n")

        overlap = peaks.get('cross_phase_overlap',
                           np.zeros(n, dtype=bool))
        for i in range(n):
            qs = peaks.get('quality_score', np.zeros(n))
            ss = peaks.get('snr_score', np.zeros(n))
            is_ = peaks.get('isolation_score', np.zeros(n))
            sy = peaks.get('symmetry_score', np.zeros(n))
            dm = peaks.get('d_match_score', np.zeros(n))

            f.write(f"{peaks['two_theta_obs'][i]:.4f},"
                    f"{peaks['d_obs'][i]:.4f},"
                    f"{peaks['d_ref'][i]:.4f},"
                    f"{peaks['d_shift'][i]:.5f},"
                    f"{peaks['fwhm'][i]:.5f},"
                    f"{peaks['snr'][i]:.2f},"
                    f"{qs[i]:.4f},"
                    f"{ss[i]:.4f},"
                    f"{is_[i]:.4f},"
                    f"{sy[i]:.4f},"
                    f"{dm[i]:.4f},"
                    f"{int(overlap[i])},"
                    f"{int(include_mask[i])},"
                    f"{K_all[i]:.5f},"
                    f"{dK_all[i]:.6f}\n")
