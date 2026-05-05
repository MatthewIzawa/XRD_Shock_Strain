"""
instrumental.py — Instrumental broadening characterisation and
deconvolution.

Two classes:
    InstrumentalStandard: holds a structural Phase plus a measured
        diffraction pattern of a known standard (LaB6, Si). Supports
        both Caglioti FWHM correction (for W-H, Scherrer) and Stokes
        Fourier deconvolution (for W-A).
    InstrumentalProfile: lightweight Caglioti carrier (U, V, W). No
        measured pattern. Supports W-H, Scherrer; W-A raises a clear
        ValueError.

Caglioti polynomial:
    FWHM(2theta)^2 = U * tan^2(theta) + V * tan(theta) + W

References
----------
Caglioti, G., Paoletti, A., Ricci, F. P. (1958). Choice of collimators
    for a crystal spectrometer for neutron diffraction. Nuclear
    Instruments 3, 223-228.
"""
import json
import warnings
from pathlib import Path

import numpy as np


class InstrumentalProfile:
    """Caglioti-polynomial carrier for instrumental broadening.

    Parameters
    ----------
    U, V, W : float
        Caglioti coefficients (deg^2).
    wavelength : float
        X-ray wavelength (angstroms).
    name : str
        Optional human-readable label.
    """

    def __init__(self, U: float, V: float, W: float,
                 wavelength: float, name: str = ''):
        self.U = float(U)
        self.V = float(V)
        self.W = float(W)
        self.wavelength = float(wavelength)
        self.name = str(name)

    def fwhm_at(self, two_theta_deg) -> float:
        """Caglioti FWHM (degrees) at the given 2-theta (degrees).

        Parameters
        ----------
        two_theta_deg : float or array-like
            2-theta value(s) in degrees. Scalar or numpy-array input is
            accepted; the return shape matches the input.

        Returns
        -------
        float or np.ndarray
            FWHM in degrees, clamped to >= 0. If the Caglioti polynomial
            evaluates negative at any input (an indication that the
            coefficients are outside their valid angular range), a
            UserWarning is emitted before clamping.
        """
        theta = np.deg2rad(np.asarray(two_theta_deg) / 2.0)
        fwhm_sq = (self.U * np.tan(theta)**2
                   + self.V * np.tan(theta)
                   + self.W)
        if np.any(fwhm_sq < 0.0):
            warnings.warn(
                'Caglioti polynomial negative at one or more 2theta '
                'values; FWHM clamped to 0. Check coefficient validity '
                'for this angular range.',
                UserWarning, stacklevel=2)
        return np.sqrt(np.maximum(fwhm_sq, 0.0))
