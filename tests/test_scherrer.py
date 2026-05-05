"""Unit tests for Scherrer K and shape factor extension (v0.3.0)."""
import numpy as np
import pytest

from xrd_profile import scherrer, modified_scherrer


CU_KA = 1.5406
FWHM = np.array([0.30, 0.40, 0.55, 0.60])
TT = np.array([20.0, 30.0, 45.0, 60.0])


class TestSchererDefaults:
    def test_no_kwargs_preserves_v020_behavior(self):
        """scherrer(fwhm, tt, wl) without kwargs uses K=0.9."""
        sizes = scherrer(FWHM, TT, CU_KA)
        sizes_explicit_09 = scherrer(FWHM, TT, CU_KA, K=0.9)
        np.testing.assert_allclose(sizes, sizes_explicit_09)

    def test_default_K_is_0_9(self):
        # Hand-computed reference: D = 0.9 * lambda / (beta_rad * cos(theta_rad))
        beta_rad = np.radians(FWHM)
        theta_rad = np.radians(TT / 2)
        expected = 0.9 * CU_KA / (beta_rad * np.cos(theta_rad))
        actual = scherrer(FWHM, TT, CU_KA)
        np.testing.assert_allclose(actual, expected)


class TestShapeKwarg:
    def test_lookup_table_values(self):
        from xrd_profile import SCHERRER_K_FOR_SHAPE
        assert SCHERRER_K_FOR_SHAPE['spherical'] == 0.94
        assert SCHERRER_K_FOR_SHAPE['cubic'] == 0.83
        assert SCHERRER_K_FOR_SHAPE['cylindrical'] == 1.84
        assert SCHERRER_K_FOR_SHAPE['platey'] == 1.0

    def test_shape_spherical_uses_K_094(self):
        sizes_shape = scherrer(FWHM, TT, CU_KA, shape='spherical')
        sizes_explicit = scherrer(FWHM, TT, CU_KA, K=0.94)
        np.testing.assert_allclose(sizes_shape, sizes_explicit)

    def test_shape_cylindrical_uses_K_184(self):
        sizes_shape = scherrer(FWHM, TT, CU_KA, shape='cylindrical')
        sizes_explicit = scherrer(FWHM, TT, CU_KA, K=1.84)
        np.testing.assert_allclose(sizes_shape, sizes_explicit)


class TestKWinsOverShape:
    def test_explicit_K_overrides_shape(self):
        sizes_K = scherrer(FWHM, TT, CU_KA, K=1.0, shape='spherical')
        sizes_K_only = scherrer(FWHM, TT, CU_KA, K=1.0)
        np.testing.assert_allclose(sizes_K, sizes_K_only)


class TestModifiedScherrer:
    def test_default_K_preserves_v020(self):
        size_default = modified_scherrer(FWHM, TT, CU_KA)
        size_explicit = modified_scherrer(FWHM, TT, CU_KA, K=0.9)
        assert np.isclose(size_default, size_explicit)

    def test_shape_works(self):
        size = modified_scherrer(FWHM, TT, CU_KA, shape='cubic')
        size_explicit = modified_scherrer(FWHM, TT, CU_KA, K=0.83)
        assert np.isclose(size, size_explicit)
