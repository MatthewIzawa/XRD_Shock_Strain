"""Unit tests for Bragg's law conversions."""
import numpy as np
import pytest
from xrd_profile.conversions import (two_theta_to_d, d_to_two_theta,
                                     two_theta_to_Q, two_theta_to_K,
                                     fwhm_to_deltaK)


CU_KA = 1.5406


def test_two_theta_to_d_known():
    """Cu Ka, 2theta=28.44 deg -> d ~ 3.135 A (anorthite 220)."""
    d = two_theta_to_d(28.44, CU_KA)
    assert abs(d - 3.135) < 0.01


def test_d_to_two_theta_roundtrip():
    """d -> 2theta -> d should be identity."""
    d_orig = np.array([3.18, 2.50, 1.75])
    tt = d_to_two_theta(d_orig, CU_KA)
    d_back = two_theta_to_d(tt, CU_KA)
    np.testing.assert_allclose(d_back, d_orig, atol=1e-10)


def test_two_theta_to_Q():
    """Q = 4*pi*sin(theta)/lambda."""
    Q = two_theta_to_Q(90.0, CU_KA)
    expected = 4 * np.pi * np.sin(np.radians(45)) / CU_KA
    assert abs(Q - expected) < 1e-10


def test_K_equals_inverse_d():
    """K = 1/d."""
    tt = np.array([20.0, 30.0, 45.0, 60.0])
    K = two_theta_to_K(tt, CU_KA)
    d = two_theta_to_d(tt, CU_KA)
    np.testing.assert_allclose(K, 1.0 / d, atol=1e-10)


def test_fwhm_to_deltaK():
    """DeltaK = beta*cos(theta)/lambda."""
    fwhm = 0.2  # degrees
    tt = 30.0
    dK = fwhm_to_deltaK(fwhm, tt, CU_KA)
    expected = np.radians(0.2) * np.cos(np.radians(15)) / CU_KA
    assert abs(dK - expected) < 1e-10


def test_array_inputs():
    """All functions should handle numpy arrays."""
    tt = np.linspace(10, 80, 100)
    d = two_theta_to_d(tt, CU_KA)
    Q = two_theta_to_Q(tt, CU_KA)
    K = two_theta_to_K(tt, CU_KA)
    assert len(d) == 100
    assert len(Q) == 100
    assert len(K) == 100
    assert np.all(d > 0)
    assert np.all(Q > 0)
    assert np.all(K > 0)


def test_wavelength_independence():
    """Same d-spacing from different wavelengths at different 2-theta."""
    d_target = 3.18  # angstroms
    tt_cu = d_to_two_theta(d_target, 1.5406)
    tt_co = d_to_two_theta(d_target, 1.7890)
    d_cu = two_theta_to_d(tt_cu, 1.5406)
    d_co = two_theta_to_d(tt_co, 1.7890)
    assert abs(d_cu - d_co) < 1e-10
