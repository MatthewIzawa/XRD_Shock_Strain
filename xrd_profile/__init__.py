"""
xrd_profile: XRD peak profile analysis toolkit.

A Python package for quantitative analysis of powder X-ray diffraction
peak profiles, providing Williamson-Hall, Warren-Averbach, Scherrer,
and pair distribution function methods with reference-guided peak
detection and wavelength-independent reciprocal-space formulations.

Adapted in part from crystallite_size_calculator by Dinga Wonanke
(https://github.com/bafgreat/crystallite_size_calculator), with
substantial additions for reference-guided analysis, reciprocal-space
methods, and improved Fourier coefficient extraction.

See README.md for full attribution and citation information.
"""

__version__ = '0.3.0'

from .conversions import (two_theta_to_d, d_to_two_theta, two_theta_to_Q,
                          two_theta_to_K, fwhm_to_deltaK)
from .noise import estimate_noise, estimate_zero_offset
from .peak_detection import (estimate_fwhm_simple, estimate_fwhm_voigt,
                             voigt_profile, find_peaks_guided,
                             check_cross_phase_overlap, score_peak_quality)
from .scherrer import scherrer, modified_scherrer, SCHERRER_K_FOR_SHAPE
from .pdf import (chebyshev_background, compute_pdf_sine,
                  measure_pdf_peaks, fit_first_pdf_peak)
from .phases import Phase, build_reference_peaks
from .instrumental import InstrumentalProfile, InstrumentalStandard
from .profile import XRDProfile

__all__ = [
    'XRDProfile',
    'Phase', 'build_reference_peaks',
    'InstrumentalProfile', 'InstrumentalStandard',
    'two_theta_to_d', 'd_to_two_theta', 'two_theta_to_Q',
    'two_theta_to_K', 'fwhm_to_deltaK',
    'estimate_noise', 'estimate_zero_offset',
    'estimate_fwhm_simple', 'estimate_fwhm_voigt',
    'voigt_profile', 'find_peaks_guided',
    'scherrer', 'modified_scherrer', 'SCHERRER_K_FOR_SHAPE',
    'check_cross_phase_overlap', 'score_peak_quality',
    'chebyshev_background', 'compute_pdf_sine',
    'measure_pdf_peaks', 'fit_first_pdf_peak',
]
