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

    def guided_williamson_hall(self, ref_d=None, tolerance_d=0.03,
                               n_sigma=3.0, min_fwhm_steps=3,
                               correct_offset=True,
                               other_phase_d=None,
                               other_phase_names=None,
                               overlap_tol_deg=0.15,
                               min_quality=0.3,
                               quality_weights=None,
                               weighted_regression=True,
                               sample_flags=None,
                               export_path=None,
                               *,
                               phase=None,
                               instrumental=None,
                               **kwargs):
        """
        Reference-guided Williamson-Hall with cross-phase overlap
        rejection, quality scoring, weighted regression, and
        reliability classification.

        Parameters
        ----------
        ref_d : np.ndarray or None
            Reference d-spacings sorted by decreasing intensity.
            Mutually exclusive with phase=.
        phase : Phase or None
            New in v0.3.0. If provided, ref_d is computed from
            phase.get_ref_d(wavelength, two_theta_range=data_range).
            Mutually exclusive with ref_d=.
        instrumental : InstrumentalStandard, InstrumentalProfile, or None
            New in v0.4.0. If provided, the Caglioti FWHM is subtracted
            in quadrature from each observed peak FWHM before the W-H
            regression. InstrumentalStandard objects are converted via
            caglioti_fit(); InstrumentalProfile objects are used directly.
            Peaks where the subtraction over-corrects (beta_obs <= beta_inst)
            are excluded and flagged in result['warnings'].
        other_phase_d : list of np.ndarray or None
            d-spacings for interfering phases (excluded from fit).
        sample_flags : dict or None
            e.g. {'maskelynite_present': True}. Triggers warnings.
        export_path : str or None
            CSV path for manual validation export.

        See guided_williamson_hall() for all parameters.

        Returns
        -------
        dict with strain, crystallite_size, r_squared, reliability,
        warnings, peak_quality, etc.
        """
        if phase is not None and ref_d is not None:
            raise ValueError('pass either ref_d or phase, not both')
        if phase is not None:
            tt_range = (float(self.two_theta.min()),
                        float(self.two_theta.max()))
            ref_d = phase.get_ref_d(self.wavelength,
                                     two_theta_range=tt_range)
        from .instrumental import (InstrumentalStandard,
                                    InstrumentalProfile)
        if instrumental is None:
            inst_profile = None
        elif isinstance(instrumental, InstrumentalStandard):
            inst_profile = instrumental.caglioti_fit()
        elif isinstance(instrumental, InstrumentalProfile):
            inst_profile = instrumental
        else:
            raise TypeError(
                f'instrumental= must be InstrumentalStandard, '
                f'InstrumentalProfile, or None; got '
                f'{type(instrumental).__name__}')
        if ref_d is None:
            raise ValueError('must pass either ref_d or phase')

        from .williamson_hall import guided_williamson_hall
        return guided_williamson_hall(
            self.two_theta, self.intensity, ref_d, self.wavelength,
            tolerance_d=tolerance_d, n_sigma=n_sigma,
            min_fwhm_steps=min_fwhm_steps, correct_offset=correct_offset,
            other_phase_d=other_phase_d,
            other_phase_names=other_phase_names,
            overlap_tol_deg=overlap_tol_deg,
            min_quality=min_quality,
            quality_weights=quality_weights,
            weighted_regression=weighted_regression,
            sample_flags=sample_flags,
            export_path=export_path,
            inst_profile=inst_profile,
            **kwargs
        )

    def guided_warren_averbach(self, ref_peaks=None, tolerance_d=0.03,
                                n_sigma=3.0, min_fwhm_steps=3,
                                correct_offset=True, n_coeffs=20,
                                width_fwhm=6.0, require_clean=False,
                                *,
                                phase=None,
                                instrumental=None):
        """
        Reference-guided Warren-Averbach analysis using harmonic families.

        Parameters
        ----------
        ref_peaks : list of dict or None
            Reference peaks with 'd', 'two_theta', 'intensity',
            'h', 'k', 'l' keys. Mutually exclusive with phase=.
        phase : Phase or None
            New in v0.3.0. If provided, ref_peaks is computed from
            phase.get_ref_peaks(wavelength, two_theta_range=data_range).
            Mutually exclusive with ref_peaks=.
        instrumental : reserved for Phase 2 / v1.0
            Pass None (default). Any other value raises NotImplementedError.
        tolerance_d, n_sigma, min_fwhm_steps, correct_offset, n_coeffs,
        width_fwhm, require_clean : see guided_warren_averbach function.

        Returns
        -------
        dict with families, crystallite sizes, strains, etc.
        """
        if phase is not None and ref_peaks is not None:
            raise ValueError('pass either ref_peaks or phase, not both')
        if phase is not None:
            tt_range = (float(self.two_theta.min()),
                        float(self.two_theta.max()))
            ref_peaks = phase.get_ref_peaks(self.wavelength,
                                             two_theta_range=tt_range)
        if instrumental is not None:
            raise NotImplementedError(
                'instrumental= is reserved for Phase 2 / v1.0; '
                'see xrd_profile roadmap')
        if ref_peaks is None:
            raise ValueError('must pass either ref_peaks or phase')

        from .warren_averbach import guided_warren_averbach
        return guided_warren_averbach(
            self.two_theta, self.intensity, ref_peaks, self.wavelength,
            tolerance_d=tolerance_d, n_sigma=n_sigma,
            min_fwhm_steps=min_fwhm_steps, correct_offset=correct_offset,
            n_coeffs=n_coeffs, width_fwhm=width_fwhm,
            require_clean=require_clean
        )

    def scherrer(self, K=None, shape=None, height_threshold=0.05):
        """Run Scherrer analysis on all detected peaks.

        K and shape: see xrd_profile.scherrer.scherrer for resolution rules.
        """
        fwhm, positions = estimate_fwhm_simple(
            self.two_theta, self.intensity, height_threshold
        )
        if len(fwhm) == 0:
            return {'sizes': np.array([]), 'peak_positions': np.array([]),
                    'd_spacings': np.array([]), 'fwhm': np.array([]),
                    'mean_size': np.nan, 'median_size': np.nan}

        sizes = scherrer(fwhm, positions, self.wavelength, K=K, shape=shape)
        d_sp = two_theta_to_d(positions, self.wavelength)
        return {
            'sizes': sizes, 'peak_positions': positions,
            'd_spacings': d_sp, 'fwhm': fwhm,
            'mean_size': np.mean(sizes), 'median_size': np.median(sizes),
        }

    def modified_scherrer(self, K=None, shape=None, height_threshold=0.05):
        """Run modified Scherrer equation. Returns average size (angstroms).

        K and shape: see xrd_profile.scherrer.modified_scherrer.
        """
        fwhm, positions = estimate_fwhm_simple(
            self.two_theta, self.intensity, height_threshold
        )
        if len(fwhm) < 2:
            return np.nan
        return modified_scherrer(fwhm, positions, self.wavelength,
                                  K=K, shape=shape)

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
        Compute pair distribution function from the diffraction pattern
        using the basic FFT method.

        For quantitative work, use compute_pdf_sine() instead.

        Returns
        -------
        r : np.ndarray - radial distances (angstroms)
        g_r : np.ndarray - PDF values
        """
        from .pdf import PDFProcessor
        processor = PDFProcessor(self.two_theta, self.intensity,
                                 self.wavelength)
        return processor.compute_pdf()

    def compute_pdf_sine(self, cheby_order=15, lorch=True,
                         r_max=30.0, n_r=2000):
        """
        Compute PDF via sine Fourier transform with Chebyshev
        background subtraction and optional Lorch modification.

        Parameters
        ----------
        cheby_order : int
            Chebyshev polynomial order for background.
        lorch : bool
            Apply Lorch modification to suppress termination ripples.
        r_max : float
            Maximum radial distance (angstroms).
        n_r : int
            Number of points in r-grid.

        Returns
        -------
        r : np.ndarray - radial distances (angstroms)
        G_r : np.ndarray - PDF values
        Q_max : float - maximum Q reached (inverse angstroms)
        """
        from .pdf import compute_pdf_sine
        return compute_pdf_sine(self.two_theta, self.intensity,
                                self.wavelength, cheby_order=cheby_order,
                                lorch=lorch, r_max=r_max, n_r=n_r)

    def measure_pdf_peaks(self, min_r=1.0, max_r=20.0, n_max=20,
                          **pdf_kwargs):
        """
        Detect and characterise peaks in the PDF.

        Parameters
        ----------
        min_r, max_r : float
            Radial search range (angstroms).
        n_max : int
            Maximum peaks to return.
        **pdf_kwargs
            Passed to compute_pdf_sine (cheby_order, lorch, etc.).

        Returns
        -------
        list of dict with 'r', 'height', 'fwhm' per peak.
        """
        from .pdf import measure_pdf_peaks
        r, G_r, _ = self.compute_pdf_sine(**pdf_kwargs)
        return measure_pdf_peaks(r, G_r, min_r=min_r, max_r=max_r,
                                 n_max=n_max)

    def fit_first_pdf_peak(self, r_min=1.2, r_max=6.0, **pdf_kwargs):
        """
        Fit a Gaussian to the first coordination-shell PDF peak.

        Parameters
        ----------
        r_min, r_max : float
            Search range (angstroms).
        **pdf_kwargs
            Passed to compute_pdf_sine (cheby_order, lorch, etc.).

        Returns
        -------
        r_peak : float - peak centre (angstroms)
        fwhm : float - FWHM (angstroms)
        """
        from .pdf import fit_first_pdf_peak
        r, G_r, _ = self.compute_pdf_sine(**pdf_kwargs)
        return fit_first_pdf_peak(r, G_r, r_min=r_min, r_max=r_max)

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

    def run_all(self,
                 methods=None,
                 phases=None,
                 wh=None,
                 wa=None,
                 pdf=None,
                 scherrer=None,
                 instrumental=None):
        """
        Run a configurable bundle of analyses.

        Parameters
        ----------
        methods : list of {'wh', 'wa', 'pdf', 'scherrer'} or None
            Which analyses to run. None = all four.
        phases : list of Phase, single Phase, or None
            For guided W-H and W-A. None = unguided forms run.
        wh, wa, pdf, scherrer : dict or None
            Per-method kwargs.
            e.g. wh={'n_sigma': 3.0, 'tolerance_d': 0.02}.
        instrumental : reserved for Phase 2; raises NotImplementedError.

        Returns
        -------
        dict. With phases:
          {'wh': {phase.name: result, ...},
           'wa': {phase.name: result, ...},
           'pdf': result,
           'scherrer': result}
        Without phases:
          {'wh': result, 'wa': result, 'pdf': result, 'scherrer': result}
        """
        if instrumental is not None:
            raise NotImplementedError(
                'instrumental= is reserved for Phase 2 / v1.0; '
                'see xrd_profile roadmap')

        if methods is None:
            methods = ['wh', 'wa', 'pdf', 'scherrer']
        wh_kwargs = wh or {}
        wa_kwargs = wa or {}
        pdf_kwargs = pdf or {}
        scherrer_kwargs = scherrer or {}

        # Normalise phases to a list (or None)
        if phases is not None and not isinstance(phases, (list, tuple)):
            phases = [phases]

        results = {}

        if 'wh' in methods:
            if phases:
                results['wh'] = {
                    p.name: self.guided_williamson_hall(phase=p, **wh_kwargs)
                    for p in phases
                }
            else:
                from .williamson_hall import williamson_hall
                results['wh'] = williamson_hall(
                    self.two_theta, self.intensity, self.wavelength,
                    **wh_kwargs)

        if 'wa' in methods:
            if phases:
                results['wa'] = {
                    p.name: self.guided_warren_averbach(phase=p, **wa_kwargs)
                    for p in phases
                }
            else:
                results['wa'] = self.warren_averbach(**wa_kwargs)

        if 'pdf' in methods:
            r, G_r, Q_max = self.compute_pdf_sine(**pdf_kwargs)
            results['pdf'] = {'r': r, 'G_r': G_r, 'Q_max': Q_max}

        if 'scherrer' in methods:
            results['scherrer'] = self.scherrer(**scherrer_kwargs)

        return results

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

    def plot_pdf(self, ax=None, r_max=30, method='sine', **kwargs):
        """
        Plot the pair distribution function.

        Parameters
        ----------
        ax : matplotlib Axes or None
        r_max : float
            Maximum r to plot (angstroms).
        method : str
            'sine' (default, Chebyshev + Lorch) or 'fft' (basic).
        """
        if method == 'sine':
            r, g_r, _ = self.compute_pdf_sine(r_max=r_max)
        else:
            r, g_r = self.compute_pdf()
        mask = (r > 0) & (r < r_max)

        if ax is None:
            fig, ax = plt.subplots(figsize=(8, 4))
        ax.plot(r[mask], g_r[mask], label=self.sample_name, **kwargs)
        ax.set_xlabel(r'r ($\AA$)')
        ax.set_ylabel('G(r)')
        ax.set_title(f'PDF: {self.sample_name}')
        ax.axhline(0, color='grey', linewidth=0.3, linestyle=':')
        ax.legend()
        return ax
