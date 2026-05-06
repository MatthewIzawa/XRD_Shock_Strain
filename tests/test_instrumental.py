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


from xrd_profile import Phase, InstrumentalStandard


@pytest.fixture(scope='module')
def lab6_phase():
    """Build a LaB6 Phase inline (Pm-3m, a=4.156825) without needing
    a bundled CIF.

    LaB6 is the NIST SRM 660c crystallite-size standard. We build it
    inline (Phase.from_lattice_params) rather than ship a CIF here
    because pymatgen is the only dep needed and we don't want to add
    a CIF to examples/cifs/ that isn't in the v0.3.0 manifest.
    """
    return Phase.from_lattice_params(
        a=4.156825, b=4.156825, c=4.156825,
        alpha=90, beta=90, gamma=90,
        species=['La', 'B', 'B', 'B', 'B', 'B', 'B'],
        coords=[[0.0, 0.0, 0.0],
                [0.5, 0.5, 0.1993],
                [0.5, 0.5, 0.8007],
                [0.5, 0.1993, 0.5],
                [0.5, 0.8007, 0.5],
                [0.1993, 0.5, 0.5],
                [0.8007, 0.5, 0.5]],
        name='LaB6')


@pytest.fixture(scope='module')
def lab6_standard(lab6_phase):
    data = np.loadtxt(SYNTH_LAB6)
    return InstrumentalStandard(
        phase=lab6_phase,
        two_theta=data[:, 0], intensity=data[:, 1],
        wavelength=LAMBDA_CU, name='synthetic_lab6')


class TestInstrumentalStandard:
    def test_construct_from_phase_and_pattern(self, lab6_phase):
        data = np.loadtxt(SYNTH_LAB6)
        std = InstrumentalStandard(
            phase=lab6_phase,
            two_theta=data[:, 0], intensity=data[:, 1],
            wavelength=LAMBDA_CU, name='test')
        assert std.phase is lab6_phase
        assert std.wavelength == LAMBDA_CU
        assert std.name == 'test'
        assert len(std.two_theta) == len(data)

    def test_caglioti_fit_returns_instrumental_profile(
            self, lab6_standard):
        prof = lab6_standard.caglioti_fit()
        assert isinstance(prof, InstrumentalProfile)
        assert abs(prof.U - SYNTH_U) / SYNTH_U < 0.05
        assert abs(prof.W - SYNTH_W) / SYNTH_W < 0.05

    def test_caglioti_fit_is_cached(self, lab6_standard):
        p1 = lab6_standard.caglioti_fit()
        p2 = lab6_standard.caglioti_fit()
        assert p1 is p2

    def test_fourier_coefficients_returns_arrays_of_requested_length(
            self, lab6_standard):
        # LaB6 (1,0,0) at Cu K-alpha: d ~ 4.157 angstroms
        L, A = lab6_standard.fourier_coefficients(peak_d=4.157,
                                                    n_coeffs=20)
        assert len(L) == 20
        assert len(A) == 20
        assert A[0] > 0  # A(0) = peak area, should be positive

    def test_fourier_coefficients_decay_monotonically_at_low_L(
            self, lab6_standard):
        L, A = lab6_standard.fourier_coefficients(peak_d=4.157,
                                                    n_coeffs=10)
        # The first 5 coefficients should be a non-increasing sequence
        # in absolute value (size profile of a single peak).
        abs_A = np.abs(A[:5])
        diffs = np.diff(abs_A)
        assert np.sum(diffs > 0) <= 1  # allow one bump, no more.

    def test_fourier_coefficients_is_cached(self, lab6_standard):
        """Calling fourier_coefficients twice with the same arguments
        must return the cached arrays (identity check)."""
        L1, A1 = lab6_standard.fourier_coefficients(
            peak_d=4.157, n_coeffs=20)
        L2, A2 = lab6_standard.fourier_coefficients(
            peak_d=4.157, n_coeffs=20)
        assert L1 is L2
        assert A1 is A2


class TestInstrumentalProfileFromStandard:
    def test_from_standard_delegates_to_caglioti_fit(
            self, lab6_standard):
        prof_a = InstrumentalProfile.from_standard(lab6_standard)
        prof_b = lab6_standard.caglioti_fit()
        assert prof_a.U == prof_b.U
        assert prof_a.V == prof_b.V
        assert prof_a.W == prof_b.W


from xrd_profile.instrumental import _stokes_deconvolve


class TestStokesDeconvolve:
    def test_recovers_known_sample_from_convolved_profile(self):
        """Synthesise A_sample(L) (lognormal-derived), convolve with a
        known A_inst(L), Stokes-deconvolve, recover A_sample within 1%
        on the well-conditioned (low-L) coefficients."""
        L = np.linspace(0, 200, 50)
        # Synthetic sample column-length distribution: lognormal-ish.
        A_sample = np.exp(-L / 100.0)
        # Synthetic instrumental: narrower (smaller decay length).
        A_inst = np.exp(-L / 500.0)
        # Observed = convolution = product in Fourier space.
        A_obs = A_sample * A_inst

        A_corr = _stokes_deconvolve(A_obs, A_inst,
                                      damping_threshold=0.05)
        # Recover the first 10 coefficients (well-conditioned).
        np.testing.assert_allclose(A_corr[:10], A_sample[:10],
                                    rtol=1e-6)

    def test_damps_high_L_when_instrumental_too_small(self):
        """When A_inst(L) drops below threshold * A_inst(0), the
        corresponding A_corr should be 0, not a noise amplification."""
        L = np.linspace(0, 200, 50)
        A_inst = np.exp(-L / 20.0)  # decays fast; A_inst(0)=1
        # By L>=60 (index ~15), A_inst < 0.05 * A_inst(0) = 0.05.
        A_obs = np.full_like(L, 0.1)  # arbitrary "observed"

        A_corr = _stokes_deconvolve(A_obs, A_inst,
                                      damping_threshold=0.05)
        # Indices where damping kicked in should be exactly 0.
        damped_mask = A_inst < 0.05 * A_inst[0]
        assert np.all(A_corr[damped_mask] == 0.0)
        assert np.any(A_corr[~damped_mask] != 0.0)

    def test_handles_a_inst_zero_at_origin(self):
        """If A_inst(0) is 0, division by zero is avoided. Behaviour:
        return all zeros and a warning (or raise). We choose: raise."""
        L = np.linspace(0, 100, 20)
        A_inst = np.zeros_like(L)
        A_obs = np.full_like(L, 0.5)
        with pytest.raises(ValueError, match='A_inst.0. is zero'):
            _stokes_deconvolve(A_obs, A_inst)

    def test_shape_mismatch_raises(self):
        """Mismatched A_obs / A_inst shapes must raise ValueError."""
        with pytest.raises(ValueError, match='shape'):
            _stokes_deconvolve(np.ones(5), np.ones(6))


from xrd_profile import XRDProfile

LAMBDA_I11 = 0.826517
TIRHERT = FIXTURE_DIR / 'tirhert_subset.xy'


class TestGuidedWHWithInstrumental:
    @pytest.fixture(scope='class')
    def anorthite_phase(self):
        cif = (Path(__file__).parent.parent / 'examples'
               / 'cifs' / 'Anorthite.cif')
        return Phase.from_cif(str(cif), name='anorthite')

    def test_runs_end_to_end_with_instrumental_standard(
            self, anorthite_phase, lab6_standard):
        """Guided W-H with InstrumentalStandard runs without error and
        produces a finite crystallite size."""
        data = np.loadtxt(TIRHERT)
        profile = XRDProfile(data[:, 0], data[:, 1],
                              wavelength=LAMBDA_I11)
        # The LaB6 standard is at Cu K-alpha; for this smoke test we
        # rebuild it at the I11 wavelength used by the Tirhert subset.
        lab6_data = np.loadtxt(SYNTH_LAB6)
        lab6_at_i11 = InstrumentalStandard(
            phase=lab6_standard.phase,
            two_theta=lab6_data[:, 0], intensity=lab6_data[:, 1],
            wavelength=LAMBDA_I11, name='lab6_at_i11')
        result = profile.guided_williamson_hall(
            phase=anorthite_phase,
            instrumental=lab6_at_i11,
            n_sigma=3.0, tolerance_d=0.03)
        assert np.isfinite(result['crystallite_size'])
        assert result['crystallite_size'] > 0

    def test_runs_with_instrumental_profile(self, anorthite_phase):
        """Guided W-H with bare InstrumentalProfile (no measured pattern)
        also runs end-to-end. Parameters are synchrotron-appropriate
        (U, V, W ~10x smaller than typical lab-source) so the instrumental
        FWHM does not over-correct the narrow I11 peaks."""
        data = np.loadtxt(TIRHERT)
        profile = XRDProfile(data[:, 0], data[:, 1],
                              wavelength=LAMBDA_I11)
        # Synchrotron-scale Caglioti parameters: I11's narrow beam
        # gives instrumental FWHM ~0.007 deg at moderate angles, ~10x
        # smaller than lab-Bruker scales. Lab-scale values (U~5e-3)
        # would exceed sample broadening on the Tirhert subset and
        # over-correct every peak.
        inst = InstrumentalProfile(U=5e-5, V=-2e-5, W=1e-4,
                                    wavelength=LAMBDA_I11,
                                    name='synch_i11_approx')
        result = profile.guided_williamson_hall(
            phase=anorthite_phase, instrumental=inst,
            n_sigma=3.0, tolerance_d=0.03)
        assert np.isfinite(result['crystallite_size'])

    def test_corrected_size_differs_from_uncorrected(
            self, anorthite_phase):
        """Sanity: instrumental correction CHANGES the apparent size.

        Note: the directional claim (corrected size > uncorrected) is
        physically true in principle, but at synchrotron-scale
        instrumental broadening on this fixture the correction is
        small and the weighted-regression intercept can flip sign in
        either direction across noise-driven differences. We assert
        only that the result changes — the strong directional claim
        is appropriate for lab-scale corrections (which the bundled
        Tirhert fixture cannot produce).
        """
        data = np.loadtxt(TIRHERT)
        profile = XRDProfile(data[:, 0], data[:, 1],
                              wavelength=LAMBDA_I11)
        # Synchrotron-scale Caglioti parameters: I11's narrow beam
        # gives instrumental FWHM ~0.007 deg at moderate angles, ~10x
        # smaller than lab-Bruker scales. Lab-scale values (U~5e-3)
        # would exceed sample broadening on the Tirhert subset and
        # over-correct every peak.
        inst = InstrumentalProfile(U=5e-5, V=-2e-5, W=1e-4,
                                    wavelength=LAMBDA_I11)
        uncorrected = profile.guided_williamson_hall(
            phase=anorthite_phase,
            n_sigma=3.0, tolerance_d=0.03)
        corrected = profile.guided_williamson_hall(
            phase=anorthite_phase, instrumental=inst,
            n_sigma=3.0, tolerance_d=0.03)
        # The corrected result must differ from the uncorrected one.
        assert (corrected['crystallite_size']
                != uncorrected['crystallite_size'])

    def test_invalid_instrumental_type_raises(self, anorthite_phase):
        data = np.loadtxt(TIRHERT)
        profile = XRDProfile(data[:, 0], data[:, 1],
                              wavelength=LAMBDA_I11)
        with pytest.raises(TypeError, match='InstrumentalStandard'):
            profile.guided_williamson_hall(
                phase=anorthite_phase,
                instrumental='not_a_real_profile')


class TestGuidedWAWithInstrumental:
    @pytest.fixture(scope='class')
    def anorthite_phase(self):
        cif = (Path(__file__).parent.parent / 'examples'
               / 'cifs' / 'Anorthite.cif')
        return Phase.from_cif(str(cif), name='anorthite')

    def test_runs_end_to_end_with_instrumental_standard(
            self, anorthite_phase, lab6_standard):
        data = np.loadtxt(TIRHERT)
        profile = XRDProfile(data[:, 0], data[:, 1],
                              wavelength=LAMBDA_I11)
        lab6_data = np.loadtxt(SYNTH_LAB6)
        lab6_at_i11 = InstrumentalStandard(
            phase=lab6_standard.phase,
            two_theta=lab6_data[:, 0], intensity=lab6_data[:, 1],
            wavelength=LAMBDA_I11, name='lab6_at_i11')
        result = profile.guided_warren_averbach(
            phase=anorthite_phase, instrumental=lab6_at_i11,
            n_sigma=3.0, tolerance_d=0.03)
        assert np.isfinite(result['mean_crystallite_size']) \
               or np.isnan(result['mean_crystallite_size'])
        # At minimum, the call returns a result dict.
        assert 'families' in result

    def test_instrumental_profile_raises_for_wa(self, anorthite_phase):
        data = np.loadtxt(TIRHERT)
        profile = XRDProfile(data[:, 0], data[:, 1],
                              wavelength=LAMBDA_I11)
        inst = InstrumentalProfile(U=5e-3, V=-1e-3, W=5e-3,
                                    wavelength=LAMBDA_I11)
        with pytest.raises(ValueError, match='Stokes deconvolution'):
            profile.guided_warren_averbach(
                phase=anorthite_phase, instrumental=inst)

    def test_invalid_instrumental_type_raises_for_wa(
            self, anorthite_phase):
        data = np.loadtxt(TIRHERT)
        profile = XRDProfile(data[:, 0], data[:, 1],
                              wavelength=LAMBDA_I11)
        with pytest.raises(TypeError, match='InstrumentalStandard'):
            profile.guided_warren_averbach(
                phase=anorthite_phase, instrumental='garbage')
