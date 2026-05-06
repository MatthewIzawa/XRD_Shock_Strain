"""Tests for xrd_profile.instrumental — InstrumentalStandard,
InstrumentalProfile, Caglioti fitting, Stokes deconvolution."""
import json
from pathlib import Path

import numpy as np
import pytest

from xrd_profile import InstrumentalProfile

FIXTURE_DIR = Path(__file__).parent / 'fixtures'
SYNTH_LAB6 = FIXTURE_DIR / 'synthetic_lab6.xy'

# Documented synthesis params for the LaB6 fixture.
SYNTH_U, SYNTH_V, SYNTH_W = 5.0e-3, -1.0e-3, 5.0e-3
LAMBDA_CU = 1.5406


class TestInstrumentalProfileBasics:
    def test_construct_with_uvw(self):
        prof = InstrumentalProfile(U=SYNTH_U, V=SYNTH_V, W=SYNTH_W,
                                    wavelength=LAMBDA_CU,
                                    name='test_profile')
        assert prof.U == SYNTH_U
        assert prof.V == SYNTH_V
        assert prof.W == SYNTH_W
        assert prof.wavelength == LAMBDA_CU
        assert prof.name == 'test_profile'

    def test_fwhm_at_recovers_synthesis(self):
        """fwhm_at() with the exact synthesis coefficients should
        return the Caglioti formula's value at any 2theta."""
        prof = InstrumentalProfile(U=SYNTH_U, V=SYNTH_V, W=SYNTH_W,
                                    wavelength=LAMBDA_CU)
        # At theta=22.5deg (2theta=45deg), tan(theta)=tan(22.5deg)
        tt = 45.0
        theta = np.deg2rad(tt / 2.0)
        expected = np.sqrt(SYNTH_U * np.tan(theta)**2
                           + SYNTH_V * np.tan(theta) + SYNTH_W)
        assert np.isclose(prof.fwhm_at(tt), expected, rtol=1e-12)

    def test_fwhm_at_is_positive_over_typical_range(self):
        prof = InstrumentalProfile(U=SYNTH_U, V=SYNTH_V, W=SYNTH_W,
                                    wavelength=LAMBDA_CU)
        for tt in np.linspace(20.0, 100.0, 50):
            assert prof.fwhm_at(tt) > 0.0

    def test_default_name_is_empty_string(self):
        prof = InstrumentalProfile(U=0.005, V=0.0, W=0.005,
                                    wavelength=LAMBDA_CU)
        assert prof.name == ''

    def test_fwhm_at_clamps_negative_polynomial_with_warning(self):
        """When Caglioti coefficients produce a negative polynomial
        (out of valid range), fwhm_at clamps to 0 AND warns."""
        # U=0, V=-1, W=0 gives fwhm_sq = -tan(theta), negative for
        # all theta > 0 (i.e., all 2theta > 0).
        prof = InstrumentalProfile(U=0.0, V=-1.0, W=0.0,
                                    wavelength=LAMBDA_CU)
        with pytest.warns(UserWarning, match='Caglioti polynomial negative'):
            result = prof.fwhm_at(45.0)
        assert result == 0.0


class TestInstrumentalProfileJsonIO:
    def test_to_json_from_json_round_trip(self, tmp_path):
        prof = InstrumentalProfile(U=SYNTH_U, V=SYNTH_V, W=SYNTH_W,
                                    wavelength=LAMBDA_CU,
                                    name='test_round_trip')
        path = tmp_path / 'profile.json'
        prof.to_json(path)
        loaded = InstrumentalProfile.from_json(path)
        assert loaded.U == prof.U
        assert loaded.V == prof.V
        assert loaded.W == prof.W
        assert loaded.wavelength == prof.wavelength
        assert loaded.name == prof.name

    def test_json_file_contains_documented_fields(self, tmp_path):
        prof = InstrumentalProfile(U=0.005, V=-0.001, W=0.005,
                                    wavelength=LAMBDA_CU,
                                    name='lab_bruker_cu_ka')
        path = tmp_path / 'p.json'
        prof.to_json(path)
        contents = json.loads(path.read_text())
        assert set(contents) == {'U', 'V', 'W', 'wavelength',
                                 'name', 'schema_version'}
        assert contents['schema_version'] == '1'


    def test_from_json_wrong_schema_version_raises(self, tmp_path):
        """from_json must reject a JSON file with an unrecognised
        schema_version."""
        bad = tmp_path / 'bad.json'
        bad.write_text(
            '{"schema_version": "99", "U": 0.005, "V": -0.001, '
            '"W": 0.005, "wavelength": 1.5406, "name": "wrong_schema"}')
        with pytest.raises(ValueError, match='schema_version'):
            InstrumentalProfile.from_json(bad)


class TestInstrumentalProfileRegistry:
    def test_from_registry_unknown_name_raises_keyerror(self):
        with pytest.raises(KeyError, match='nonexistent_profile'):
            InstrumentalProfile.from_registry('nonexistent_profile')


from xrd_profile.instrumental import _caglioti_fit


class TestCagliotiFit:
    def test_recovers_synthesis_coefficients_within_5pct(self):
        """Fitting Caglioti to the synthetic LaB6 fixture should
        recover the documented synthesis U, V, W within 5%."""
        data = np.loadtxt(SYNTH_LAB6)
        tt, intensity = data[:, 0], data[:, 1]

        # Reference peak positions (LaB6 cubic, a=4.156825 angstroms,
        # Cu K-alpha) — first 8 reflections in 2-theta order.
        ref_tt = np.array([
            21.358, 30.385, 37.443, 43.527, 48.999, 54.087,
            63.198, 67.494,
        ])

        U, V, W, info = _caglioti_fit(tt, intensity, ref_tt)
        assert abs(U - SYNTH_U) / SYNTH_U < 0.05, \
            f'U: expected {SYNTH_U}, got {U}'
        assert abs(W - SYNTH_W) / SYNTH_W < 0.05, \
            f'W: expected {SYNTH_W}, got {W}'
        # V is small (-1e-3) and noise-prone; absolute tolerance 5e-4
        # corresponds to 50% relative tolerance — V is noise-dominated
        # at this angular coverage, so a tighter tolerance would be
        # flaky on the synthetic fixture.
        assert abs(V - SYNTH_V) < 5.0e-4, \
            f'V: expected {SYNTH_V}, got {V}'
        assert info['n_peaks'] == len(ref_tt)

    def test_fit_info_contains_documented_keys(self):
        data = np.loadtxt(SYNTH_LAB6)
        tt, intensity = data[:, 0], data[:, 1]
        ref_tt = np.array([21.358, 30.385, 37.443, 43.527])
        _, _, _, info = _caglioti_fit(tt, intensity, ref_tt)
        assert set(info) >= {'n_peaks', 'measured_fwhms',
                             'measured_positions', 'cov'}

    def test_raises_when_too_few_peaks_resolvable(self):
        """`_caglioti_fit` must raise ValueError when fewer than 4
        peaks resolve in the search windows."""
        data = np.loadtxt(SYNTH_LAB6)
        tt, intensity = data[:, 0], data[:, 1]
        # Pick reference positions far from any actual LaB6 peak.
        ref_tt = np.array([22.5, 25.5, 28.5])  # 3 phantom positions
        with pytest.raises(ValueError, match='at least 4'):
            _caglioti_fit(tt, intensity, ref_tt)
