"""
XRDProfile: unified interface for XRD peak profile analysis.

Provides a single class that wraps all analysis methods (Williamson-Hall,
Warren-Averbach, Scherrer, PDF) with coordinate conversion and plotting
utilities.
"""

import numpy as np
import matplotlib.pyplot as plt

from .conversions import (two_theta_to_d, two_theta_to_Q, two_theta_to_K,
                          fwhm_to_deltaK)
from .peak_detection import estimate_fwhm_simple
from .scherrer import scherrer, modified_scherrer


class XRDProfile:
    """
    Unified interface for XRD peak profile analysis.

    Parameters
    ----------
    two_theta : np.ndarray
        2-theta array in degrees.
    intensity : np.ndarray
        Intensity array.
    wavelength : float
        X-ray wavelength in angstroms (default: 1.5406 for Cu K-alpha).
    sample_name : str
        Sample identifier.
    """

    def __init__(self, two_theta, intensity, wavelength=1.5406,
                 sample_name=''):
        self.two_theta = np.asarray(two_theta, dtype=float)
        self.intensity = np.asarray(intensity, dtype=float)
        self.wavelength = wavelength
        self.sample_name = sample_name

    # --- Coordinate conversions ---

    def get_d_spacing(self):
        """Return d-spacing array corresponding to the 2-theta array."""
        return two_theta_to_d(self.two_theta, self.wavelength)

    def get_Q(self):
        """Return scattering vector Q array."""
        return two_theta_to_Q(self.two_theta, self.wavelength)

    def get_K(self):
        """Return reciprocal variable K = 1/d array."""
        return two_theta_to_K(self.two_theta, self.wavelength)

    # --- Peak profile analyses ---

    def williamson_hall(self, use_voigt=False, height_threshold=0.05):
        """Run conventional Williamson-Hall analysis. Returns result dict."""
        from .williamson_hall import williamson_hall
        return williamson_hall(self.two_theta, self.intensity,
                              self.wavelength, use_voigt, height_threshold)

    def williamson_hall_reciprocal(self, use_voigt=False,
                                   height_threshold=0.05):
        """Run reciprocal-space Williamson-Hall (DeltaK vs K)."""
        from .williamson_hall import williamson_hall_reciprocal
        return williamson_hall_reciprocal(self.two_theta, self.intensity,
                                         self.wavelength, use_voigt,
                                         height_threshold)

    def guided_williamson_hall(self, ref_d, tolerance_d=0.03,
                               n_sigma=3.0, min_fwhm_steps=3,
                               correct_offset=True):
        """
        Reference-guided Williamson-Hall in reciprocal space.

        Parameters
        ----------
        ref_d : np.ndarray
            Reference d-spacings sorted by decreasing intensity.
        tolerance_d : float
            d-spacing matching tolerance (angstroms).
        n_sigma : float
            Noise cutoff (multiples of sigma).
        min_fwhm_steps : int
            Minimum FWHM in step-size units.
        correct_offset : bool
            Estimate and apply zero-point offset.

        Returns
        -------
        dict with strain, crystallite_size, r_squared, peaks, etc.
        """
        from .williamson_hall import guided_williamson_hall
        return guided_williamson_hall(
            self.two_theta, self.intensity, ref_d, self.wavelength,
            tolerance_d=tolerance_d, n_sigma=n_sigma,
            min_fwhm_steps=min_fwhm_steps, correct_offset=correct_offset
        )

    def guided_warren_averbach(self, ref_peaks, tolerance_d=0.03,
                                n_sigma=3.0, min_fwhm_steps=3,
                                correct_offset=True, n_coeffs=20,
                                width_fwhm=6.0, require_clean=False):
        """
        Reference-guided Warren-Averbach analysis using harmonic families.

        Parameters
        ----------
        ref_peaks : list of dict
            Reference peaks with 'd', 'two_theta', 'intensity',
            'h', 'k', 'l' keys.
        tolerance_d, n_sigma, min_fwhm_steps, correct_offset, n_coeffs,
        width_fwhm, require_clean : see guided_warren_averbach function.

        Returns
        -------
        dict with families, crystallite sizes, strains, etc.
        """
        from .warren_averbach import guided_warren_averbach
        return guided_warren_averbach(
            self.two_theta, self.intensity, ref_peaks, self.wavelength,
            tolerance_d=tolerance_d, n_sigma=n_sigma,
            min_fwhm_steps=min_fwhm_steps, correct_offset=correct_offset,
            n_coeffs=n_coeffs, width_fwhm=width_fwhm,
            require_clean=require_clean
        )

    def scherrer(self, K=0.9, height_threshold=0.05):
        """Run Scherrer analysis on all detected peaks."""
        fwhm, positions = estimate_fwhm_simple(
            self.two_theta, self.intensity, height_threshold
        )
        if len(fwhm) == 0:
            return {'sizes': np.array([]), 'peak_positions': np.array([]),
                    'd_spacings': np.array([]), 'fwhm': np.array([]),
                    'mean_size': np.nan, 'median_size': np.nan}

        sizes = scherrer(fwhm, positions, self.wavelength, K)
        d_sp = two_theta_to_d(positions, self.wavelength)
        return {
            'sizes': sizes, 'peak_positions': positions,
            'd_spacings': d_sp, 'fwhm': fwhm,
            'mean_size': np.mean(sizes), 'median_size': np.median(sizes),
        }

    def modified_scherrer(self, K=0.9, height_threshold=0.05):
        """Run modified Scherrer equation. Returns average size (angstroms)."""
        fwhm, positions = estimate_fwhm_simple(
            self.two_theta, self.intensity, height_threshold
        )
        if len(fwhm) < 2:
            return np.nan
        return modified_scherrer(fwhm, positions, self.wavelength, K)

    def warren_averbach(self, initial_size=None, height_threshold=0.05):
        """Run simplified (unguided) Warren-Averbach analysis."""
        from .warren_averbach import warren_averbach
        if initial_size is None:
            sch = self.scherrer(height_threshold=height_threshold)
            initial_size = sch['mean_size']
            if np.isnan(initial_size):
                initial_size = 500.0
        return warren_averbach(self.two_theta, self.intensity,
                               initial_size, self.wavelength,
                               height_threshold)

    def compute_pdf(self):
        """
        Compute pair distribution function from the diffraction pattern.

        Returns
        -------
        r : np.ndarray - radial distances (angstroms)
        g_r : np.ndarray - PDF values
        """
        from .pdf import PDFProcessor
        processor = PDFProcessor(self.two_theta, self.intensity,
                                 self.wavelength)
        return processor.compute_pdf()

    def crystallite_size_from_pdf(self, shape='spherical'):
        """Estimate crystallite size from the PDF envelope function."""
        from .pdf import crystallite_size_from_pdf
        r, g_r = self.compute_pdf()
        mask = r > 0
        return crystallite_size_from_pdf(r[mask], g_r[mask], shape)

    def full_analysis(self, use_voigt=False, height_threshold=0.05):
        """Run all unguided analyses and return comprehensive results."""
        from .williamson_hall import (williamson_hall,
                                      williamson_hall_reciprocal)
        wh = williamson_hall(self.two_theta, self.intensity,
                            self.wavelength, use_voigt, height_threshold)
        wh_r = williamson_hall_reciprocal(self.two_theta, self.intensity,
                                          self.wavelength, use_voigt,
                                          height_threshold)
        sch = self.scherrer(height_threshold=height_threshold)
        mod_sch = self.modified_scherrer(height_threshold=height_threshold)
        wa = self.warren_averbach(height_threshold=height_threshold)

        try:
            pdf_size = self.crystallite_size_from_pdf()
        except Exception:
            pdf_size = np.nan

        return {
            'sample_name': self.sample_name,
            'williamson_hall': wh,
            'williamson_hall_reciprocal': wh_r,
            'scherrer': sch,
            'modified_scherrer': mod_sch,
            'warren_averbach': wa,
            'pdf_crystallite_size': pdf_size,
        }

    # --- Plotting ---

    def plot_pattern(self, ax=None, offset=0, x_axis='d_spacing', **kwargs):
        """
        Plot the diffraction pattern.

        Parameters
        ----------
        ax : matplotlib Axes or None
        offset : float - vertical offset for stacked plots
        x_axis : str - 'two_theta', 'd_spacing', 'Q', or 'K'
        """
        if ax is None:
            fig, ax = plt.subplots(figsize=(10, 4))

        if x_axis == 'd_spacing':
            x = self.get_d_spacing()
            xlabel = r'd-spacing ($\AA$)'
        elif x_axis == 'Q':
            x = self.get_Q()
            xlabel = r'Q ($\AA^{-1}$)'
        elif x_axis == 'K':
            x = self.get_K()
            xlabel = r'K = 1/d ($\AA^{-1}$)'
        else:
            x = self.two_theta
            xlabel = r'$2\theta$ (degrees)'

        ax.plot(x, self.intensity + offset,
                label=self.sample_name, **kwargs)
        ax.set_xlabel(xlabel)
        ax.set_ylabel('Intensity')
        return ax

    def plot_williamson_hall(self, ax=None, use_voigt=False, **kwargs):
        """Generate conventional Williamson-Hall plot."""
        wh = self.williamson_hall(use_voigt=use_voigt)
        if len(wh['x']) < 2:
            return None

        if ax is None:
            fig, ax = plt.subplots(figsize=(6, 4))

        ax.scatter(wh['x'], wh['y'], marker='o', s=40,
                   label=self.sample_name, **kwargs)
        x_fit = np.linspace(wh['x'].min(), wh['x'].max(), 100)
        y_fit = wh['slope'] * x_fit + wh['intercept']
        ax.plot(x_fit, y_fit, '--', alpha=0.7)
        ax.set_xlabel(r'$4\sin\theta$')
        ax.set_ylabel(r'$\beta\cos\theta$')
        ax.set_title(f'Williamson-Hall: {self.sample_name}')
        ax.legend()
        return ax

    def plot_williamson_hall_reciprocal(self, ax=None, use_voigt=False,
                                        **kwargs):
        """Generate reciprocal-space Williamson-Hall plot."""
        wh = self.williamson_hall_reciprocal(use_voigt=use_voigt)
        if len(wh['K']) < 2:
            return None

        if ax is None:
            fig, ax = plt.subplots(figsize=(6, 4))

        ax.scatter(wh['K'], wh['deltaK'], marker='o', s=40,
                   label=self.sample_name, **kwargs)
        x_fit = np.linspace(wh['K'].min(), wh['K'].max(), 100)
        y_fit = wh['slope'] * x_fit + wh['intercept']
        ax.plot(x_fit, y_fit, '--', alpha=0.7)
        ax.set_xlabel(r'K = 2sin$\theta$/$\lambda$ ($\AA^{-1}$)')
        ax.set_ylabel(r'$\Delta$K ($\AA^{-1}$)')
        ax.set_title(f'W-H reciprocal: {self.sample_name}')
        ax.legend()
        return ax

    def plot_pdf(self, ax=None, r_max=30, **kwargs):
        """Plot the pair distribution function."""
        r, g_r = self.compute_pdf()
        mask = (r > 0) & (r < r_max)

        if ax is None:
            fig, ax = plt.subplots(figsize=(8, 4))
        ax.plot(r[mask], g_r[mask], label=self.sample_name, **kwargs)
        ax.set_xlabel(r'r ($\AA$)')
        ax.set_ylabel('g(r)')
        ax.set_title(f'PDF: {self.sample_name}')
        ax.legend()
        return ax
