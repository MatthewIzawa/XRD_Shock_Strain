"""Unit tests for enhanced PDF analysis functions."""
import numpy as np
import pytest
from xrd_profile.pdf import (chebyshev_background, compute_pdf_sine,
                              measure_pdf_peaks, fit_first_pdf_peak)


CU_KA = 1.5406


def _synthetic_pattern(wavelength=CU_KA, tt_min=10, tt_max=80, n_pts=3500):
    """Generate a synthetic diffraction pattern with known peaks on a
    polynomial background for testing."""
    two_theta = np.linspace(tt_min, tt_max, n_pts)
    # Smooth polynomial background
    bg = 1000 * np.exp(-0.03 * two_theta) + 200
    # Add Gaussian peaks at known 2-theta positions
    peak_positions = [22.0, 28.0, 36.0, 42.0, 50.0, 58.0, 65.0]
    peak_heights = [800, 1200, 600, 400, 900, 300, 500]
    peak_widths = [0.15, 0.12, 0.18, 0.20, 0.14, 0.22, 0.16]
    intensity = bg.copy()
    for pos, h, w in zip(peak_positions, peak_heights, peak_widths):
        intensity += h * np.exp(-(two_theta - pos)**2 / (2 * w**2))
    # Add small noise
    rng = np.random.default_rng(42)
    intensity += rng.normal(0, 5, len(intensity))
    return two_theta, intensity


class TestChebyshevBackground:
    """Tests for iterative Chebyshev background estimation."""

    def test_returns_correct_shapes(self):
        """Background and coeffs should have expected shapes."""
        Q = np.linspace(1.0, 10.0, 1000)
        I = 100 * np.exp(-0.2 * Q) + np.random.default_rng(0).normal(0, 1, 1000)
        bg, coeffs = chebyshev_background(Q, I, order=10)
        assert bg.shape == Q.shape
        assert len(coeffs) == 11  # order + 1

    def test_background_below_peaks(self):
        """Background should converge below the Bragg peaks."""
        Q = np.linspace(1.0, 10.0, 2000)
        bg_true = 50 * np.exp(-0.1 * Q) + 20
        # Add sharp peaks
        peaks = bg_true.copy()
        for pos in [3.0, 5.0, 7.0]:
            peaks += 200 * np.exp(-(Q - pos)**2 / (2 * 0.02**2))
        bg_est, _ = chebyshev_background(Q, peaks, order=10)
        # The estimated background should be closer to the true background
        # than to the peak intensities, on average
        residual_from_true = np.mean(np.abs(bg_est - bg_true))
        residual_from_peaks = np.mean(np.abs(bg_est - peaks))
        assert residual_from_true < residual_from_peaks

    def test_different_orders(self):
        """Higher order should fit more complex backgrounds."""
        Q = np.linspace(1.0, 10.0, 1000)
        I = np.sin(Q) * 50 + 100  # wavy background, no peaks
        bg_low, _ = chebyshev_background(Q, I, order=3)
        bg_high, _ = chebyshev_background(Q, I, order=15)
        # Higher order should fit better
        err_low = np.std(I - bg_low)
        err_high = np.std(I - bg_high)
        assert err_high <= err_low


class TestComputePdfSine:
    """Tests for the sine Fourier transform PDF."""

    def test_returns_correct_types(self):
        """Should return (r, G_r, Q_max) with correct shapes."""
        tt, intensity = _synthetic_pattern()
        r, G_r, Q_max = compute_pdf_sine(tt, intensity, CU_KA, n_r=500)
        assert isinstance(r, np.ndarray)
        assert isinstance(G_r, np.ndarray)
        assert len(r) == 500
        assert len(G_r) == 500
        assert Q_max > 0

    def test_q_max_scales_with_wavelength(self):
        """Shorter wavelength should give larger Q_max."""
        tt, intensity = _synthetic_pattern()
        _, _, Q_max_cu = compute_pdf_sine(tt, intensity, CU_KA, n_r=100)
        # Synchrotron wavelength (shorter)
        _, _, Q_max_syn = compute_pdf_sine(tt, intensity, 0.8265, n_r=100)
        assert Q_max_syn > Q_max_cu

    def test_lorch_suppresses_oscillation(self):
        """Lorch modification should reduce high-frequency oscillation."""
        tt, intensity = _synthetic_pattern()
        _, G_lorch, _ = compute_pdf_sine(tt, intensity, CU_KA,
                                          lorch=True, n_r=500)
        _, G_nolorch, _ = compute_pdf_sine(tt, intensity, CU_KA,
                                            lorch=False, n_r=500)
        # Lorch version should have smaller variance at large r
        # (less termination ripple)
        var_lorch = np.var(G_lorch[-100:])
        var_nolorch = np.var(G_nolorch[-100:])
        assert var_lorch <= var_nolorch


class TestMeasurePdfPeaks:
    """Tests for PDF peak detection."""

    def test_finds_peaks_in_synthetic_pdf(self):
        """Should detect peaks in a synthetic G(r)."""
        r = np.linspace(0.5, 20.0, 2000)
        # Synthetic PDF with known peaks
        G_r = (0.8 * np.exp(-(r - 3.3)**2 / (2 * 0.15**2))
               + 0.5 * np.exp(-(r - 5.2)**2 / (2 * 0.2**2))
               + 0.3 * np.exp(-(r - 7.5)**2 / (2 * 0.25**2)))
        peaks = measure_pdf_peaks(r, G_r, min_r=1.0, max_r=15.0)
        assert len(peaks) >= 3
        # First peak should be near 3.3 A
        assert abs(peaks[0]['r'] - 3.3) < 0.2
        assert peaks[0]['height'] > 0

    def test_returns_empty_for_flat(self):
        """No peaks in a flat G(r)."""
        r = np.linspace(1.0, 20.0, 500)
        G_r = np.zeros_like(r)
        peaks = measure_pdf_peaks(r, G_r)
        assert len(peaks) == 0

    def test_fwhm_is_reasonable(self):
        """Measured FWHM should be close to the known width."""
        r = np.linspace(0.5, 20.0, 4000)
        sigma = 0.15
        G_r = np.exp(-(r - 5.0)**2 / (2 * sigma**2))
        peaks = measure_pdf_peaks(r, G_r, min_r=3.0, max_r=8.0)
        assert len(peaks) >= 1
        expected_fwhm = 2.355 * sigma
        assert abs(peaks[0]['fwhm'] - expected_fwhm) < 0.1


class TestFitFirstPdfPeak:
    """Tests for Gaussian fitting to the first coordination shell."""

    def test_recovers_known_peak(self):
        """Should recover position and width of a known Gaussian peak."""
        r = np.linspace(0.5, 10.0, 2000)
        r0, sigma = 3.3, 0.12
        G_r = 0.8 * np.exp(-(r - r0)**2 / (2 * sigma**2))
        r_peak, fwhm = fit_first_pdf_peak(r, G_r, r_min=1.5, r_max=6.0)
        assert abs(r_peak - r0) < 0.05
        assert abs(fwhm - 2.355 * sigma) < 0.05

    def test_returns_nan_for_no_peak(self):
        """Should return NaN when no positive peak exists."""
        r = np.linspace(1.0, 10.0, 500)
        G_r = -0.1 * np.ones_like(r)
        r_peak, fwhm = fit_first_pdf_peak(r, G_r)
        assert np.isnan(r_peak)
        assert np.isnan(fwhm)

    def test_selects_first_peak(self):
        """When multiple peaks exist, should fit the first one."""
        r = np.linspace(0.5, 10.0, 2000)
        G_r = (0.8 * np.exp(-(r - 2.5)**2 / (2 * 0.1**2))
               + 0.6 * np.exp(-(r - 5.0)**2 / (2 * 0.15**2)))
        r_peak, fwhm = fit_first_pdf_peak(r, G_r, r_min=1.5, r_max=6.0)
        # Should find the first peak near 2.5, not the second at 5.0
        assert abs(r_peak - 2.5) < 0.1
