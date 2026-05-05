"""
Numerical regression test against frozen v0.2.0 outputs.

Asserts that the v0.2.0 array-based public API produces identical
output (within tight tolerances) to what was generated against the
v0.2.0 codebase. If this test fails after a code change, that change
has perturbed v0.2.0 default behavior -- the strict-additive policy of
Phase 1 forbids that. Either fix the code or regenerate the golden
file with explicit reasoning (see scripts/regenerate_goldens.py).
"""
import json
from pathlib import Path

import numpy as np
import pytest

from xrd_profile import (XRDProfile, scherrer, modified_scherrer,
                         compute_pdf_sine, estimate_fwhm_simple)
from xrd_profile import Phase

FIXTURE_DIR = Path(__file__).parent / 'fixtures'
LAMBDA_I11 = 0.826517

# Same constants as scripts/regenerate_goldens.py -- kept in sync.
ANORTHITE_REF_D = [
    3.20, 3.18, 3.65, 4.04, 6.41, 5.69, 3.74, 3.21, 4.04, 2.94,
]
ANORTHITE_REF_PEAKS = [
    {'d': 3.20, 'two_theta': 14.84, 'intensity': 100.0, 'h': 0, 'k': 4, 'l': 0},
    {'d': 3.18, 'two_theta': 14.93, 'intensity':  85.0, 'h': 2, 'k': 0, 'l': -2},
    {'d': 6.41, 'two_theta':  7.39, 'intensity':  60.0, 'h': 0, 'k': 2, 'l': 0},
    {'d': 4.04, 'two_theta': 11.74, 'intensity':  50.0, 'h': 0, 'k': 0, 'l': 2},
    {'d': 3.65, 'two_theta': 13.00, 'intensity':  45.0, 'h': 1, 'k': 3, 'l': 0},
]
ANORTHITE_CIF = (Path(__file__).parent.parent / 'examples'
                 / 'cifs' / 'Anorthite.cif')


@pytest.fixture(scope='module')
def golden():
    return json.loads((FIXTURE_DIR / 'golden_v0.2.0_results.json').read_text())


@pytest.fixture(scope='module')
def pattern():
    data = np.loadtxt(FIXTURE_DIR / 'tirhert_subset.xy')
    return data[:, 0], data[:, 1]


def _assert_close_scalar(name, actual, expected, rtol=1e-6, atol=1e-10):
    if expected is None:
        assert actual is None or np.isnan(actual), f'{name}: expected None, got {actual}'
        return
    if np.isnan(expected):
        assert np.isnan(actual), f'{name}: expected NaN, got {actual}'
        return
    assert np.isclose(actual, expected, rtol=rtol, atol=atol), (
        f'{name}: expected {expected}, got {actual}')


def _assert_close_array(name, actual, expected, rtol=1e-6, atol=1e-10):
    actual = np.asarray(actual)
    expected = np.asarray(expected)
    assert actual.shape == expected.shape, (
        f'{name}: shape mismatch {actual.shape} vs {expected.shape}')
    np.testing.assert_allclose(actual, expected, rtol=rtol, atol=atol,
                                err_msg=f'{name}: numerical drift')


class TestGuidedWilliamsonHall:
    def test_crystallite_size_matches_v020(self, pattern, golden):
        tt, I = pattern
        profile = XRDProfile(tt, I, wavelength=LAMBDA_I11)
        result = profile.guided_williamson_hall(
            np.array(ANORTHITE_REF_D), n_sigma=3.0, tolerance_d=0.03)
        _assert_close_scalar(
            'crystallite_size',
            result['crystallite_size'],
            golden['guided_williamson_hall']['crystallite_size'])

    def test_strain_matches_v020(self, pattern, golden):
        tt, I = pattern
        profile = XRDProfile(tt, I, wavelength=LAMBDA_I11)
        result = profile.guided_williamson_hall(
            np.array(ANORTHITE_REF_D), n_sigma=3.0, tolerance_d=0.03)
        _assert_close_scalar(
            'strain', result['strain'],
            golden['guided_williamson_hall']['strain'])


class TestGuidedWarrenAverbach:
    def test_median_crystallite_size_matches_v020(self, pattern, golden):
        tt, I = pattern
        profile = XRDProfile(tt, I, wavelength=LAMBDA_I11)
        result = profile.guided_warren_averbach(
            ANORTHITE_REF_PEAKS, n_sigma=3.0, tolerance_d=0.03)
        _assert_close_scalar(
            'median_crystallite_size',
            result['median_crystallite_size'],
            golden['guided_warren_averbach']['median_crystallite_size'])


class TestComputePdfSine:
    def test_pdf_arrays_match_v020(self, pattern, golden):
        tt, I = pattern
        r, G_r, Q_max = compute_pdf_sine(tt, I, LAMBDA_I11,
                                          cheby_order=15, lorch=True)
        _assert_close_array('r', r, golden['compute_pdf_sine']['r'])
        _assert_close_array('G_r', G_r, golden['compute_pdf_sine']['G_r'])
        _assert_close_scalar('Q_max', Q_max, golden['compute_pdf_sine']['Q_max'])


class TestScherrer:
    def test_scherrer_default_K_matches_v020(self, pattern, golden):
        tt, I = pattern
        fwhm, positions = estimate_fwhm_simple(tt, I, height_threshold=0.05)
        if golden['scherrer_default'] is None:
            pytest.skip('No detectable peaks in fixture for default Scherrer')
        sizes = scherrer(fwhm, positions, LAMBDA_I11)
        _assert_close_array('scherrer_sizes', sizes,
                            golden['scherrer_default']['sizes'])

    def test_modified_scherrer_default_K_matches_v020(self, pattern, golden):
        tt, I = pattern
        fwhm, positions = estimate_fwhm_simple(tt, I, height_threshold=0.05)
        if golden['modified_scherrer_default'] is None:
            pytest.skip('Insufficient peaks for modified Scherrer')
        size = modified_scherrer(fwhm, positions, LAMBDA_I11)
        _assert_close_scalar('modified_scherrer', size,
                              golden['modified_scherrer_default'])


# --- v0.3.0 tier (added in v0.4.0) ---
# Asserts: every key in golden_v0.3.0_results.json is reproducible at
# the current tag with byte-equivalent numerical value (key-subset
# value-equality). New top-level keys added at v0.4+ tags are allowed
# in the live result and ignored here.


@pytest.fixture(scope='module')
def golden_v030():
    return json.loads(
        (FIXTURE_DIR / 'golden_v0.3.0_results.json').read_text())


@pytest.fixture(scope='module')
def anorthite_phase():
    return Phase.from_cif(str(ANORTHITE_CIF), name='anorthite')


class TestV030GuidedViaPhase:
    def test_wh_crystallite_size_matches_v030(
            self, pattern, golden_v030, anorthite_phase):
        tt, I = pattern
        profile = XRDProfile(tt, I, wavelength=LAMBDA_I11)
        result = profile.guided_williamson_hall(
            phase=anorthite_phase, n_sigma=3.0, tolerance_d=0.03)
        _assert_close_scalar(
            'wh_crystallite_size',
            result['crystallite_size'],
            golden_v030['guided_wh_via_phase']['crystallite_size'])

    def test_wa_crystallite_size_matches_v030(
            self, pattern, golden_v030, anorthite_phase):
        tt, I = pattern
        profile = XRDProfile(tt, I, wavelength=LAMBDA_I11)
        result = profile.guided_warren_averbach(
            phase=anorthite_phase, n_sigma=3.0, tolerance_d=0.03)
        _assert_close_scalar(
            'wa_mean_crystallite_size',
            result['mean_crystallite_size'],
            golden_v030['guided_wa_via_phase']['mean_crystallite_size'])


class TestV030ScherrerShapeTable:
    def test_scherrer_spherical_matches_v030(self, pattern, golden_v030):
        tt, I = pattern
        profile = XRDProfile(tt, I, wavelength=LAMBDA_I11)
        result = profile.scherrer(shape='spherical')
        _assert_close_scalar(
            'scherrer_spherical_mean',
            result['mean_size'],
            golden_v030['scherrer_spherical']['mean_size'])

    def test_scherrer_cubic_matches_v030(self, pattern, golden_v030):
        tt, I = pattern
        profile = XRDProfile(tt, I, wavelength=LAMBDA_I11)
        result = profile.scherrer(shape='cubic')
        _assert_close_scalar(
            'scherrer_cubic_mean',
            result['mean_size'],
            golden_v030['scherrer_cubic']['mean_size'])

    def test_scherrer_default_matches_v030(self, pattern, golden_v030):
        tt, I = pattern
        profile = XRDProfile(tt, I, wavelength=LAMBDA_I11)
        result = profile.scherrer()
        _assert_close_scalar(
            'scherrer_default_mean',
            result['mean_size'],
            golden_v030['scherrer_default']['mean_size'])


class TestV030RunAll:
    def test_run_all_no_phases_wh_matches_v030(
            self, pattern, golden_v030):
        tt, I = pattern
        profile = XRDProfile(tt, I, wavelength=LAMBDA_I11)
        result = profile.run_all(methods=['wh', 'scherrer'])
        # Unguided W-H returns a single dict (not phase-keyed) when
        # phases is None — assert its crystallite_size matches the
        # golden's snapshot.
        _assert_close_scalar(
            'run_all_wh_crystallite_size',
            result['wh']['crystallite_size'],
            golden_v030['run_all_no_phases']['wh']['crystallite_size'])

    def test_run_all_no_phases_scherrer_matches_v030(
            self, pattern, golden_v030):
        tt, I = pattern
        profile = XRDProfile(tt, I, wavelength=LAMBDA_I11)
        result = profile.run_all(methods=['wh', 'scherrer'])
        _assert_close_scalar(
            'run_all_scherrer_mean',
            result['scherrer']['mean_size'],
            golden_v030['run_all_no_phases']['scherrer']['mean_size'])
