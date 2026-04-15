"""
Bragg's law coordinate conversions for powder diffraction data.

Provides wavelength-independent conversions between 2-theta, d-spacing,
scattering vector Q, and reciprocal variable K = 1/d.
"""

import numpy as np


def two_theta_to_d(two_theta, wavelength):
    """
    Convert 2-theta (degrees) to d-spacing (angstroms) via Bragg's law.

    d = lambda / (2 * sin(theta))

    Parameters
    ----------
    two_theta : float or np.ndarray
        Diffraction angle(s) in degrees.
    wavelength : float
        X-ray wavelength in angstroms.

    Returns
    -------
    d : float or np.ndarray
        d-spacing(s) in angstroms.
    """
    theta_rad = np.radians(np.asarray(two_theta) / 2)
    sin_theta = np.sin(theta_rad)
    with np.errstate(divide='ignore', invalid='ignore'):
        d = wavelength / (2 * sin_theta)
    return d


def d_to_two_theta(d, wavelength):
    """
    Convert d-spacing (angstroms) to 2-theta (degrees) via Bragg's law.

    Parameters
    ----------
    d : float or np.ndarray
        d-spacing(s) in angstroms.
    wavelength : float
        X-ray wavelength in angstroms.

    Returns
    -------
    two_theta : float or np.ndarray
        Diffraction angle(s) in degrees.
    """
    d = np.asarray(d)
    with np.errstate(invalid='ignore'):
        two_theta = 2 * np.degrees(np.arcsin(wavelength / (2 * d)))
    return two_theta


def two_theta_to_Q(two_theta, wavelength):
    """
    Convert 2-theta (degrees) to scattering vector magnitude Q.

    Q = 4*pi*sin(theta) / lambda  [inverse angstroms]

    Parameters
    ----------
    two_theta : float or np.ndarray
        Diffraction angle(s) in degrees.
    wavelength : float
        X-ray wavelength in angstroms.

    Returns
    -------
    Q : float or np.ndarray
        Scattering vector magnitude(s) in inverse angstroms.
    """
    theta_rad = np.radians(np.asarray(two_theta) / 2)
    return (4 * np.pi * np.sin(theta_rad)) / wavelength


def two_theta_to_K(two_theta, wavelength):
    """
    Convert 2-theta (degrees) to reciprocal variable K = 1/d.

    K = 2*sin(theta)/lambda  [inverse angstroms]

    This is the natural variable for Williamson-Hall analysis in
    reciprocal space (Das Bakshi et al., 2018, Materials
    Characterization 142, 144-153).

    Parameters
    ----------
    two_theta : float or np.ndarray
        Diffraction angle(s) in degrees.
    wavelength : float
        X-ray wavelength in angstroms.

    Returns
    -------
    K : float or np.ndarray
        Reciprocal variable(s) in inverse angstroms.
    """
    theta_rad = np.radians(np.asarray(two_theta) / 2)
    return (2 * np.sin(theta_rad)) / wavelength


def fwhm_to_deltaK(fwhm_deg, two_theta_deg, wavelength):
    """
    Convert FWHM in 2-theta (degrees) to DeltaK in reciprocal space.

    DeltaK = beta * cos(theta) / lambda

    where beta = FWHM in radians.

    Parameters
    ----------
    fwhm_deg : float or np.ndarray
        FWHM in degrees 2-theta.
    two_theta_deg : float or np.ndarray
        Peak position(s) in degrees 2-theta.
    wavelength : float
        X-ray wavelength in angstroms.

    Returns
    -------
    deltaK : float or np.ndarray
        Peak broadening in reciprocal space (inverse angstroms).
    """
    beta_rad = np.radians(np.asarray(fwhm_deg))
    theta_rad = np.radians(np.asarray(two_theta_deg) / 2)
    return (beta_rad * np.cos(theta_rad)) / wavelength
