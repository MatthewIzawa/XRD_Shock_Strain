"""
warren_averbach.py — Warren-Averbach Fourier analysis of XRD peak profiles.

Separates crystallite size and microstrain contributions to peak
broadening via Fourier decomposition of harmonic peak families.

Includes adaptive extraction windows, Tukey tapering, uniform s-grid
interpolation, overlap detection, and Fourier coefficient quality
filtering.

References
----------
Warren, B. E. & Averbach, B. L. (1950). The effect of cold-work
    distortion on X-ray patterns. J. Appl. Phys. 21, 595-599.

Adapted from crystallite_size_calculator
(https://github.com/bafgreat/crystallite_size_calculator) by Dinga Wonanke.
Modified by M.R.M. Izawa.
"""

import numpy as np
from math import gcd
from scipy.stats import linregress

from .conversions import two_theta_to_d, d_to_two_theta, two_theta_to_K
from .noise import estimate_noise, estimate_zero_offset
from .peak_detection import estimate_fwhm_simple


def extract_peak_profile(two_theta, intensity, peak_tt, wavelength,
                         fwhm_est=None, width_fwhm=6.0, bg_window=151,
                         ref_peaks_tt=None, min_isolation=None):
    """
    Extract a single peak profile with background subtraction, adaptive
    window, Tukey tapering, and overlap flagging.

    Parameters
    ----------
    two_theta : np.ndarray
        Full 2-theta array (degrees).
    intensity : np.ndarray
        Full intensity array.
    peak_tt : float
        2-theta position of the peak centre (degrees).
    wavelength : float
        X-ray wavelength (angstroms).
    fwhm_est : float or None
        Estimated FWHM (degrees). If None, measured from the data.
        Used to set adaptive extraction window.
    width_fwhm : float
        Extraction window = +/- width_fwhm * FWHM/2 around peak.
        Default 6.0 (3 FWHM each side).
    bg_window : int
        Median filter window for global background estimation.
    ref_peaks_tt : list/array or None
        2-theta positions of all reference peaks. Used to flag
        overlap: if a reference peak falls within the extraction
        window, overlap_fraction is reported.
    min_isolation : float or None
        Minimum distance (degrees) to nearest reference peak for
        the profile to be considered clean. If None, defaults to
        2 * FWHM.

    Returns
    -------
    dict with keys:
        s : np.ndarray — reciprocal variable on uniform grid
        s0 : float — s at peak centre
        profile : np.ndarray — background-subtracted, tapered, normalised
        two_theta_region : np.ndarray
        raw : np.ndarray
        background : np.ndarray
        overlap_flag : bool — True if a strong neighbour contaminates
        nearest_neighbor_deg : float — distance to nearest ref peak
        fwhm_used : float — FWHM used for window sizing
    """
    from scipy.ndimage import median_filter

    step = np.median(np.diff(two_theta))
    centre_idx = np.argmin(np.abs(two_theta - peak_tt))

    # Adaptive FWHM estimation if not provided
    if fwhm_est is None:
        peak_val = intensity[centre_idx]
        # Use the global background
        bg_full = median_filter(intensity, size=bg_window)
        bg_at_peak = bg_full[centre_idx]
        half_max = (peak_val + bg_at_peak) / 2.0
        # Search left and right for half-max crossing
        left_idx = centre_idx
        while left_idx > 0 and intensity[left_idx] > half_max:
            left_idx -= 1
        right_idx = centre_idx
        while right_idx < len(intensity) - 1 and intensity[right_idx] > half_max:
            right_idx += 1
        fwhm_est = two_theta[right_idx] - two_theta[left_idx]
        fwhm_est = max(fwhm_est, 3 * step)  # floor at 3 steps

    # Extraction window: adaptive, clamped to reasonable range
    half_width_deg = max(width_fwhm * fwhm_est / 2, 0.3)
    half_width_deg = min(half_width_deg, 2.0)  # cap at 2 degrees
    half_width = int(half_width_deg / step)

    lo = max(0, centre_idx - half_width)
    hi = min(len(two_theta), centre_idx + half_width + 1)

    tt_region = two_theta[lo:hi]
    i_region = intensity[lo:hi].astype(float)

    # Background: global median-filter background, extracted for region
    bg_full = median_filter(intensity, size=bg_window)
    bg_region = bg_full[lo:hi].astype(float)

    profile = i_region - bg_region
    profile = np.maximum(profile, 0)

    # Tukey (tapered cosine) window to suppress truncation artifacts
    n_pts = len(profile)
    if n_pts > 4:
        taper_frac = 0.3  # taper the outer 30% on each side
        taper_n = int(n_pts * taper_frac / 2)
        window = np.ones(n_pts)
        if taper_n > 0:
            taper = 0.5 * (1 - np.cos(np.pi * np.arange(taper_n) / taper_n))
            window[:taper_n] = taper
            window[-taper_n:] = taper[::-1]
        profile = profile * window

    # Normalise so area = 1
    area = np.trapz(profile, tt_region)
    if area > 0:
        profile_norm = profile / area
    else:
        profile_norm = profile

    # Convert to reciprocal variable s = 2*sin(theta)/lambda
    theta_rad = np.radians(tt_region / 2)
    s_raw = (2 * np.sin(theta_rad)) / wavelength
    s0 = (2 * np.sin(np.radians(peak_tt / 2))) / wavelength

    # Interpolate to uniform s-grid for clean Fourier transform
    s_uniform = np.linspace(s_raw[0], s_raw[-1], len(s_raw))
    profile_uniform = np.interp(s_uniform, s_raw, profile_norm)

    # Overlap detection
    nearest_deg = np.inf
    overlap_flag = False
    if min_isolation is None:
        min_isolation = 2.0 * fwhm_est

    if ref_peaks_tt is not None:
        for rtt in ref_peaks_tt:
            dist = abs(rtt - peak_tt)
            if dist < 0.05:  # same peak
                continue
            if dist < nearest_deg:
                nearest_deg = dist
            if dist < min_isolation:
                overlap_flag = True

    return {
        's': s_uniform, 's0': s0,
        'profile': profile_uniform,
        'two_theta_region': tt_region,
        'raw': i_region,
        'background': bg_region,
        'overlap_flag': overlap_flag,
        'nearest_neighbor_deg': nearest_deg,
        'fwhm_used': fwhm_est,
    }


def fourier_coefficients(s, profile, s0, n_coeffs=30):
    """
    Compute the real Fourier cosine coefficients A(L) of a peak profile
    on a uniform s-grid.

    A(L_n) = integral{ P(s) * cos(2*pi*L_n*(s - s0)) ds }

    normalised so A(L_0) = 1.

    The column length spacing is L_n = n / (2 * Delta_s), where
    Delta_s is the s-range of the profile.

    Parameters
    ----------
    s : np.ndarray
        Reciprocal variable on uniform grid.
    profile : np.ndarray
        Normalised peak profile (area ~ 1).
    s0 : float
        s at peak centre.
    n_coeffs : int
        Number of Fourier coefficients to compute.

    Returns
    -------
    L : np.ndarray
        Column lengths (angstroms).
    A_L : np.ndarray
        Real Fourier cosine coefficients, normalised so A(0) = 1.
    converged : bool
        True if A(L) decays monotonically without strong oscillation
        over the first n_coeffs/2 values.
    """
    ds = s - s0
    delta_s = s[-1] - s[0]

    if delta_s <= 0 or np.sum(profile) <= 0:
        return np.zeros(n_coeffs), np.zeros(n_coeffs), False

    # Use trapezoidal integration on uniform grid
    A_L = np.zeros(n_coeffs)
    L = np.zeros(n_coeffs)

    for n in range(n_coeffs):
        L[n] = n / (2 * delta_s)
        cos_term = np.cos(2 * np.pi * L[n] * ds)
        A_L[n] = np.trapz(profile * cos_term, s)

    # Normalise so A(0) = 1
    if A_L[0] != 0:
        A_L = A_L / A_L[0]

    # Quality check: A(L) should decay and not oscillate strongly
    half = n_coeffs // 2
    sign_changes = np.sum(np.abs(np.diff(np.sign(A_L[1:half]))) > 0)
    early_negative = np.any(A_L[1:max(3, half // 3)] < -0.1)
    converged = (sign_changes <= half // 3) and not early_negative

    return L, A_L, converged


def warren_averbach(two_theta, intensities, initial_crystallite_size,
                    wavelength, height_threshold=0.05):
    """
    Warren-Averbach method: Fourier analysis to separate size and strain
    broadening contributions.

    Parameters
    ----------
    two_theta : np.ndarray
        Full 2-theta array (degrees).
    intensities : np.ndarray
        Full intensity array.
    initial_crystallite_size : float
        Initial estimate of crystallite size in angstroms.
    wavelength : float
        X-ray wavelength in angstroms.
    height_threshold : float
        Minimum peak height as fraction of maximum intensity.

    Returns
    -------
    dict with keys:
        strain : float
        crystallite_size : float (angstroms)
        n_peaks : int
        peak_positions : np.ndarray
        fwhm : np.ndarray
    """
    fwhm_data, peak_positions = estimate_fwhm_simple(
        two_theta, intensities, height_threshold
    )

    if len(fwhm_data) < 3:
        return {
            'strain': np.nan, 'crystallite_size': np.nan,
            'n_peaks': len(fwhm_data), 'peak_positions': peak_positions,
            'fwhm': fwhm_data,
        }

    fwhm_corrected = np.sqrt(fwhm_data)
    fwhm_radians = np.radians(fwhm_corrected)
    theta = np.radians(peak_positions / 2)
    d_spacing = wavelength / (2 * np.sin(theta))
    l_values = d_spacing
    beta_total = fwhm_radians

    a_l = (np.exp(-l_values / initial_crystallite_size)
           * np.exp(-2 * np.pi**2 * (beta_total**2) * l_values**2))

    log_a_l = np.log(a_l)
    coefficients = np.polyfit(l_values, log_a_l, 2)

    d_crys = -1 / coefficients[1] if coefficients[1] != 0 else np.inf
    strain_sq = -coefficients[0] / (2 * np.pi**2)
    strain = np.sqrt(max(strain_sq, 0))

    return {
        'strain': strain,
        'crystallite_size': d_crys,
        'n_peaks': len(fwhm_data),
        'peak_positions': peak_positions,
        'fwhm': fwhm_data,
    }


def guided_warren_averbach(two_theta, intensity, ref_peaks, wavelength,
                           tolerance_d=0.03, n_sigma=3.0, min_fwhm_steps=3,
                           correct_offset=True, n_coeffs=20,
                           width_fwhm=6.0, min_ref_intensity=1.0,
                           require_clean=False):
    """
    Reference-guided Warren-Averbach analysis using harmonic peak families.

    Improvements over naive implementation:
    - Adaptive extraction window sized to each peak's FWHM
    - Tukey tapering to suppress truncation artifacts
    - Uniform s-grid interpolation for clean Fourier transform
    - Overlap detection: peaks with a strong neighbour within 2*FWHM
      are flagged (and optionally rejected via require_clean)
    - Fourier coefficient quality filtering: families where A(L)
      oscillates or goes negative early are rejected
    - Global median-filter background rather than linear edge fit

    The classical W-A equation for harmonic reflections of order n:

        ln A(L, s_n) = ln A_S(L) - 2*pi^2 * L^2 * n^2 * <e^2(L)> / d^2

    Parameters
    ----------
    two_theta, intensity : np.ndarray
        Observed diffraction pattern.
    ref_peaks : list of dict
        Reference peak list with 'd', 'two_theta', 'intensity',
        'h', 'k', 'l' keys.
    wavelength : float
        X-ray wavelength (angstroms).
    tolerance_d : float
        d-spacing matching tolerance (angstroms).
    n_sigma : float
        Noise cutoff in sigma units.
    min_fwhm_steps : int
        Minimum FWHM in step-size units.
    correct_offset : bool
        Estimate and apply zero-point offset.
    n_coeffs : int
        Number of Fourier coefficients per peak.
    width_fwhm : float
        Extraction window = +/- width_fwhm * FWHM / 2.
    min_ref_intensity : float
        Minimum reference intensity (%) to include a reflection.
    require_clean : bool
        If True, reject peaks flagged as overlapping.

    Returns
    -------
    dict with keys:
        families : list of dict — one per analysed family
        families_rejected : list of dict — families rejected by quality
        mean_crystallite_size, median_crystallite_size, mean_rms_strain : float
        zero_offset : float
        n_families, n_families_rejected : int
    """
    # Zero-point offset
    ref_d_sorted = np.array(sorted(
        [p['d'] for p in ref_peaks if p['intensity'] >= 3.0],
        reverse=True))
    if correct_offset:
        offset = estimate_zero_offset(two_theta, intensity,
                                       ref_d_sorted, wavelength)
        tt_corr = two_theta + offset
    else:
        offset = 0.0
        tt_corr = two_theta

    # Noise estimation
    noise_info = estimate_noise(tt_corr, intensity)
    sigma = noise_info['sigma']
    step_size = np.median(np.diff(tt_corr))

    # All reference peak 2-theta positions (for overlap detection)
    all_ref_tt = []
    for p in ref_peaks:
        if p['intensity'] >= min_ref_intensity:
            rtt = d_to_two_theta(p['d'], wavelength)
            if not np.isnan(rtt):
                all_ref_tt.append(rtt + offset)
    all_ref_tt = np.array(all_ref_tt)

    # Build harmonic families from reference
    families_ref = {}
    for p in ref_peaks:
        if p['intensity'] < min_ref_intensity:
            continue
        h, k, l = p['h'], p['k'], p['l']
        g = gcd(gcd(abs(h), abs(k)), abs(l)) if (h or k or l) else 1
        base = (h // g, k // g, l // g)
        order = g
        if base not in families_ref:
            families_ref[base] = []
        families_ref[base].append({'order': order, 'peak': p})

    multi_families = {k: sorted(v, key=lambda x: x['order'])
                      for k, v in families_ref.items() if len(v) >= 2}

    results_families = []
    rejected_families = []

    for base_hkl, members in multi_families.items():
        family_profiles = []
        any_overlap = False

        for member in members:
            p = member['peak']
            d_ref = p['d']
            order = member['order']

            tt_expected = d_to_two_theta(d_ref, wavelength) + offset
            if np.isnan(tt_expected):
                continue
            if tt_expected < tt_corr.min() + 1 or tt_expected > tt_corr.max() - 1:
                continue

            # Search for peak near expected position
            d_pattern = two_theta_to_d(tt_corr, wavelength)
            mask = np.abs(d_pattern - d_ref) <= tolerance_d
            if not np.any(mask):
                continue

            indices = np.where(mask)[0]
            local_i = intensity[indices]
            peak_idx = indices[np.argmax(local_i)]
            peak_val = intensity[peak_idx]
            bg_val = noise_info['background'][peak_idx]

            if (peak_val - bg_val) < n_sigma * sigma:
                continue

            # FWHM measurement
            half_max = (peak_val + bg_val) / 2
            left_below = np.where(intensity[:peak_idx] < half_max)[0]
            right_below = np.where(intensity[peak_idx:] < half_max)[0]
            if len(left_below) == 0 or len(right_below) == 0:
                continue

            # Interpolated FWHM
            li = left_below[-1]
            ri = right_below[0] + peak_idx
            if li + 1 < len(intensity) and intensity[li + 1] != intensity[li]:
                frac_l = (half_max - intensity[li]) / (intensity[li + 1] - intensity[li])
                tt_left = tt_corr[li] + frac_l * step_size
            else:
                tt_left = tt_corr[li]
            if ri > 0 and intensity[ri - 1] != intensity[ri]:
                frac_r = (half_max - intensity[ri]) / (intensity[ri - 1] - intensity[ri])
                tt_right = tt_corr[ri] - frac_r * step_size
            else:
                tt_right = tt_corr[ri]
            fwhm = tt_right - tt_left

            if fwhm < min_fwhm_steps * step_size or fwhm > 5.0:
                continue

            # Extract profile with adaptive window and overlap detection
            peak_tt_obs = tt_corr[peak_idx]
            prof = extract_peak_profile(
                tt_corr, intensity, peak_tt_obs, wavelength,
                fwhm_est=fwhm, width_fwhm=width_fwhm,
                ref_peaks_tt=all_ref_tt
            )

            if prof['overlap_flag']:
                any_overlap = True
                if require_clean:
                    continue

            # Fourier coefficients with quality check
            L, A_L, converged = fourier_coefficients(
                prof['s'], prof['profile'], prof['s0'],
                n_coeffs=n_coeffs
            )

            if not converged:
                continue

            d_obs = two_theta_to_d(peak_tt_obs, wavelength)
            family_profiles.append({
                'order': order,
                'd_obs': float(d_obs),
                'd_ref': float(d_ref),
                'L': L,
                'A_L': A_L,
                'peak_tt': float(peak_tt_obs),
                'fwhm': float(fwhm),
                'overlap': prof['overlap_flag'],
                'nearest_neighbor': prof['nearest_neighbor_deg'],
            })

        if len(family_profiles) < 2:
            continue

        # --- Warren-Averbach separation ---
        orders = np.array([fp['order'] for fp in family_profiles])
        n_sq = orders.astype(float) ** 2
        d_base = family_profiles[0]['d_obs']
        L_vals = family_profiles[0]['L']

        A_size = np.ones(n_coeffs)
        mean_sq_strain = np.zeros(n_coeffs)

        valid_L_count = 0
        for li in range(1, n_coeffs):
            ln_A_vals = []
            for fp in family_profiles:
                if li < len(fp['A_L']) and fp['A_L'][li] > 0.01:
                    ln_A_vals.append(np.log(fp['A_L'][li]))
                else:
                    ln_A_vals.append(np.nan)

            ln_A_vals = np.array(ln_A_vals)
            valid = ~np.isnan(ln_A_vals)

            if np.sum(valid) >= 2:
                sl, interc, _, _, _ = linregress(n_sq[valid], ln_A_vals[valid])
                A_size[li] = np.exp(interc)
                strain_sq = -sl * d_base**2 / (2 * np.pi**2)
                mean_sq_strain[li] = max(strain_sq, 0)
                valid_L_count += 1
            else:
                A_size[li] = np.nan
                mean_sq_strain[li] = np.nan

        if valid_L_count < 3:
            rejected_families.append({
                'base_hkl': base_hkl,
                'reason': 'insufficient valid L points',
                'valid_L_count': valid_L_count,
            })
            continue

        # Quality check: A_S(L) should decrease monotonically
        valid_as = ~np.isnan(A_size[1:8]) & (A_size[1:8] > 0)
        if np.sum(valid_as) < 3:
            rejected_families.append({
                'base_hkl': base_hkl,
                'reason': 'A_S(L) has too few valid points',
            })
            continue

        As_valid = A_size[1:8][valid_as]
        L_valid = L_vals[1:8][valid_as]

        # Check for strong non-monotonicity (allow small noise)
        increases = np.sum(np.diff(As_valid) > 0.15)
        if increases > len(As_valid) // 2:
            rejected_families.append({
                'base_hkl': base_hkl,
                'reason': 'A_S(L) strongly non-monotonic',
            })
            continue

        # Crystallite size from initial slope of A_S(L)
        sl_as, int_as, r_as, _, _ = linregress(L_valid, As_valid)
        if sl_as < 0 and int_as > 0:
            cryst_size = -int_as / sl_as
        else:
            cryst_size = np.nan

        # RMS strain at representative L values
        valid_strain = ~np.isnan(mean_sq_strain[1:8]) & (mean_sq_strain[1:8] > 0)
        if np.any(valid_strain):
            representative_strain = np.sqrt(
                np.nanmean(mean_sq_strain[1:8][valid_strain])
            )
        else:
            representative_strain = np.nan

        results_families.append({
            'base_hkl': base_hkl,
            'orders': [fp['order'] for fp in family_profiles],
            'd_spacings': [fp['d_obs'] for fp in family_profiles],
            'fwhm_values': [fp['fwhm'] for fp in family_profiles],
            'A_size': A_size,
            'L': L_vals,
            'mean_sq_strain': mean_sq_strain,
            'crystallite_size': cryst_size,
            'rms_strain': representative_strain,
            'A_size_r2': r_as**2,
            'has_overlap': any_overlap,
        })

    # Aggregate — exclude obvious outliers
    sizes = [f['crystallite_size'] for f in results_families
             if not np.isnan(f['crystallite_size'])
             and 0 < f['crystallite_size'] < 10000]
    strains = [f['rms_strain'] for f in results_families
               if not np.isnan(f['rms_strain'])
               and 0 < f['rms_strain'] < 1.0]

    return {
        'families': results_families,
        'families_rejected': rejected_families,
        'mean_crystallite_size': np.mean(sizes) if sizes else np.nan,
        'median_crystallite_size': np.median(sizes) if sizes else np.nan,
        'mean_rms_strain': np.mean(strains) if strains else np.nan,
        'zero_offset': offset,
        'n_families': len(results_families),
        'n_families_rejected': len(rejected_families),
    }
