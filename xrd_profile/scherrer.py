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

Adapted from crystallite_size_calculator
(https://github.com/bafgreat/crystallite_size_calculator) by Dinga Wonanke.
Modified by M.R.M. Izawa.
"""

import numpy as np
from scipy.stats import linregress


def scherrer(fwhm_deg, two_theta_positions, wavelength, K=0.9):
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
    K : float
        Shape factor (default 0.9).

    Returns
    -------
    crystallite_sizes : np.ndarray
        Crystallite sizes in angstroms for each peak.
    """
    fwhm_rad = np.radians(np.asarray(fwhm_deg))
    theta_rad = np.radians(np.asarray(two_theta_positions) / 2)
    return (K * wavelength) / (fwhm_rad * np.cos(theta_rad))


def modified_scherrer(fwhm_deg, two_theta_positions, wavelength, K=0.9):
    """
    Modified Scherrer equation using log-linear regression:
        ln(beta) = ln(K*lambda/D) + ln(1/cos(theta))

    Returns a single average crystallite size.

    Parameters
    ----------
    fwhm_deg : array-like
        FWHM values in degrees.
    two_theta_positions : array-like
        2-theta peak positions in degrees.
    wavelength : float
        X-ray wavelength in angstroms.
    K : float
        Shape factor (default 0.9).

    Returns
    -------
    crystallite_size : float
        Average crystallite size in angstroms.
    """
    fwhm_rad = np.radians(np.asarray(fwhm_deg))
    ln_beta = np.log(fwhm_rad)
    ln_1_cos = np.log(1 / np.cos(np.radians(np.asarray(two_theta_positions) / 2)))
    slope, intercept, _, _, _ = linregress(ln_1_cos, ln_beta)
    return (K * wavelength) / np.exp(intercept)
