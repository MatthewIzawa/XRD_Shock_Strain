"""
scherrer.py — Scherrer equation analysis for crystallite size estimation.

Provides standard per-peak Scherrer analysis and the modified Scherrer
method using log-linear regression for an average crystallite size.

References
----------
Scherrer, P. (1918). Bestimmung der Grosse und der inneren Struktur
    von Kolloidteilchen mittels Rontgenstrahlen. Nachrichten von der
    Gesellschaft der Wissenschaften zu Gottingen, Mathematisch-
    Physikalische Klasse, 98-100.
Langford, J. I. & Wilson, A. J. C. (1978). Scherrer after sixty years:
    a survey and some new results in the determination of crystallite
    size. Journal of Applied Crystallography, 11, 102-113.

Adapted from crystallite_size_calculator
(https://github.com/bafgreat/crystallite_size_calculator) by Dinga Wonanke.
Modified by M.R.M. Izawa.
"""

import numpy as np
from scipy.stats import linregress


SCHERRER_K_FOR_SHAPE = {
    'spherical':   0.94,
    'cubic':       0.83,
    'cylindrical': 1.84,
    'platey':      1.0,
}


def _resolve_K(K, shape):
    """Resolve effective K from K and shape kwargs.

    both None  -> 0.9 (v0.2.0 default)
    shape only -> SCHERRER_K_FOR_SHAPE[shape]
    K only     -> K
    both       -> K wins, shape silently ignored
    """
    if K is None and shape is None:
        return 0.9
    if K is None and shape is not None:
        try:
            return SCHERRER_K_FOR_SHAPE[shape]
        except KeyError:
            raise ValueError(
                f'unknown shape {shape!r}; choose from '
                f'{list(SCHERRER_K_FOR_SHAPE)}')
    return K


def scherrer(fwhm_deg, two_theta_positions, wavelength,
             K=None, shape=None):
    """
    Scherrer equation: D = K * lambda / (beta * cos(theta))

    Parameters
    ----------
    fwhm_deg : array-like
        FWHM values in degrees.
    two_theta_positions : array-like
        2-theta peak positions in degrees.
    wavelength : float
        X-ray wavelength in angstroms.
    K : float or None
        Scherrer constant. None = use v0.2.0 default 0.9 unless shape=
        is provided, in which case K is looked up from
        SCHERRER_K_FOR_SHAPE.
    shape : {'spherical', 'cubic', 'cylindrical', 'platey'} or None
        Crystallite shape. Used to look up K when K is None. If both
        are provided, K wins.

    Returns
    -------
    crystallite_sizes : np.ndarray
        Crystallite sizes in angstroms for each peak.
    """
    K_eff = _resolve_K(K, shape)
    fwhm_rad = np.radians(np.asarray(fwhm_deg))
    theta_rad = np.radians(np.asarray(two_theta_positions) / 2)
    return (K_eff * wavelength) / (fwhm_rad * np.cos(theta_rad))


def modified_scherrer(fwhm_deg, two_theta_positions, wavelength,
                       K=None, shape=None):
    """
    Modified Scherrer equation using log-linear regression.

    K, shape resolution: same as scherrer().
    """
    K_eff = _resolve_K(K, shape)
    fwhm_rad = np.radians(np.asarray(fwhm_deg))
    ln_beta = np.log(fwhm_rad)
    ln_1_cos = np.log(1 / np.cos(np.radians(np.asarray(two_theta_positions) / 2)))
    slope, intercept, _, _, _ = linregress(ln_1_cos, ln_beta)
    return (K_eff * wavelength) / np.exp(intercept)


def _filter_to_phase_peaks(fwhm_deg, positions_deg, ref_d, wavelength,
                            tolerance_d=0.03):
    """Keep only peaks whose d-spacing matches a reference d within
    `tolerance_d`. Returns filtered (fwhm, positions) arrays —
    (fwhm, positions) order matches `estimate_fwhm_simple`'s output."""
    from .conversions import two_theta_to_d
    fwhm = np.asarray(fwhm_deg)
    positions = np.asarray(positions_deg)
    if len(positions) == 0:
        return fwhm, positions
    d_obs = two_theta_to_d(positions, wavelength)
    ref_d = np.asarray(ref_d)
    keep_mask = np.zeros(len(d_obs), dtype=bool)
    for i, d in enumerate(d_obs):
        if np.any(np.abs(ref_d - d) < tolerance_d):
            keep_mask[i] = True
    return fwhm[keep_mask], positions[keep_mask]


def _caglioti_subtract_fwhm(fwhm_deg, positions_deg, inst_profile):
    """Caglioti-subtract instrumental FWHM from each measured FWHM.
    Drops peaks where beta_obs <= beta_inst (over-correction).

    Returns (fwhm_corr, positions_kept) — the (fwhm, positions) order
    matches `estimate_fwhm_simple` and the call sites in profile.py.
    """
    fwhm_obs = np.asarray(fwhm_deg, dtype=float)
    pos = np.asarray(positions_deg, dtype=float)
    fwhm_inst = np.array([inst_profile.fwhm_at(t) for t in pos])
    diff_sq = fwhm_obs**2 - fwhm_inst**2
    keep = diff_sq > 0
    fwhm_corr = np.sqrt(np.maximum(diff_sq[keep], 0))
    return fwhm_corr, pos[keep]
