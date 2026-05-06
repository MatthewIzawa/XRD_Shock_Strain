"""Unit tests for the Phase class and reference-peak helpers."""
import sys
from unittest.mock import patch

import numpy as np
import pytest


def test_phase_from_cif_without_pymatgen_raises_clear_error():
    """When pymatgen is not importable, Phase.from_cif raises an
    ImportError whose message names the [cif] extra install command."""
    with patch.dict(sys.modules, {'pymatgen': None,
                                    'pymatgen.core.structure': None,
                                    'pymatgen.analysis.diffraction.xrd': None}):
        if 'xrd_profile.phases' in sys.modules:
            del sys.modules['xrd_profile.phases']
        from xrd_profile.phases import Phase
        with pytest.raises(ImportError) as exc_info:
            Phase.from_cif('nonexistent.cif')
        assert 'pip install xrd_profile[cif]' in str(exc_info.value)


# Skip subsequent tests if pymatgen isn't available in the test env.
pytest.importorskip('pymatgen')


# Quartz (low) lattice + atoms — small, well-known structure for tests.
QUARTZ_A = 4.9133
QUARTZ_C = 5.4053
QUARTZ_SPECIES = ['Si', 'Si', 'Si', 'O', 'O', 'O', 'O', 'O', 'O']
QUARTZ_COORDS = [
    [0.4697, 0.0000, 0.0000],
    [0.0000, 0.4697, 0.6667],
    [0.5303, 0.5303, 0.3333],
    [0.4135, 0.2669, 0.1191],
    [0.7331, 0.1466, 0.4524],
    [0.8534, 0.5865, 0.7857],
    [0.2669, 0.4135, 0.5476],
    [0.5865, 0.8534, 0.2143],
    [0.1466, 0.7331, 0.8809],
]
CU_KA = 1.5406


@pytest.fixture
def quartz_phase():
    from xrd_profile.phases import Phase
    return Phase.from_lattice_params(
        QUARTZ_A, QUARTZ_A, QUARTZ_C, 90, 90, 120,
        species=QUARTZ_SPECIES, coords=QUARTZ_COORDS, name='Quartz')


class TestPhaseFromLatticeParams:
    def test_returns_phase_with_name(self, quartz_phase):
        assert quartz_phase.name == 'Quartz'

    def test_structure_is_pymatgen(self, quartz_phase):
        from pymatgen.core.structure import Structure
        assert isinstance(quartz_phase.structure, Structure)

    def test_repr_includes_formula(self, quartz_phase):
        r = repr(quartz_phase)
        assert 'Quartz' in r
        assert 'SiO2' in r or 'O2Si' in r


class TestGetRefPeaks:
    def test_returns_list_of_dicts_with_required_keys(self, quartz_phase):
        peaks = quartz_phase.get_ref_peaks(CU_KA)
        assert len(peaks) > 0
        required = {'d', 'two_theta', 'intensity', 'h', 'k', 'l'}
        for p in peaks:
            assert required.issubset(p.keys())

    def test_min_intensity_filters_weak_peaks(self, quartz_phase):
        all_peaks = quartz_phase.get_ref_peaks(CU_KA, min_intensity=0.0)
        strong = quartz_phase.get_ref_peaks(CU_KA, min_intensity=50.0)
        assert len(strong) <= len(all_peaks)
        for p in strong:
            assert p['intensity'] >= 50.0


class TestGetRefD:
    def test_returns_numpy_array_of_floats(self, quartz_phase):
        ref_d = quartz_phase.get_ref_d(CU_KA)
        assert isinstance(ref_d, np.ndarray)
        assert ref_d.dtype == float

    def test_sorted_by_intensity_descending_by_default(self, quartz_phase):
        ref_d = quartz_phase.get_ref_d(CU_KA, sorted_by_intensity=True)
        peaks = quartz_phase.get_ref_peaks(CU_KA)
        peaks_by_intensity = sorted(peaks, key=lambda p: -p['intensity'])
        expected_d = [p['d'] for p in peaks_by_intensity]
        np.testing.assert_allclose(ref_d, expected_d)


class TestBuildReferencePeaks:
    def test_equivalent_to_phase_method(self, quartz_phase):
        from xrd_profile.phases import build_reference_peaks
        via_func = build_reference_peaks(quartz_phase.structure, CU_KA)
        via_method = quartz_phase.get_ref_peaks(CU_KA)
        assert len(via_func) == len(via_method)
        for a, b in zip(via_func, via_method):
            for key in ('d', 'two_theta', 'intensity', 'h', 'k', 'l'):
                assert a[key] == b[key]


class TestXRDProfileGuidedWilliamsonHallPhaseKwarg:
    """Tests for the phase= kwarg added in v0.3.0."""

    def _make_profile_and_phase(self):
        from xrd_profile import XRDProfile
        # Use the bundled regression-test fixture
        from pathlib import Path
        fix = Path(__file__).parent / 'fixtures' / 'tirhert_subset.xy'
        data = np.loadtxt(fix)
        profile = XRDProfile(data[:, 0], data[:, 1],
                              wavelength=0.826517, sample_name='Tirhert')
        # Anorthite via from_lattice_params (same numbers as the
        # synchrotron_low_shock.py example).
        from xrd_profile.phases import Phase
        an = Phase.from_lattice_params(
            8.1809, 12.881, 7.1101, 93.465, 116.11, 90.369,
            species=['Ca','Al','Al','Si','Si','O','O','O','O','O','O','O','O'],
            coords=[
                [0.269,0.988,0.086],[0.507,0.314,0.621],[0.992,0.815,0.118],
                [0.505,0.320,0.110],[0.006,0.816,0.613],[0.491,0.625,0.487],
                [0.024,0.124,0.995],[0.073,0.488,0.635],[0.576,0.990,0.143],
                [0.298,0.356,0.612],[0.817,0.855,0.142],[0.517,0.179,0.610],
                [0.000,0.680,0.104],
            ],
            name='Anorthite',
        )
        return profile, an

    def test_phase_kwarg_produces_same_result_as_manual_ref_d(self):
        profile, an = self._make_profile_and_phase()
        tt_range = (float(profile.two_theta.min()),
                    float(profile.two_theta.max()))
        manual_ref_d = an.get_ref_d(profile.wavelength,
                                     two_theta_range=tt_range)
        manual = profile.guided_williamson_hall(
            manual_ref_d, n_sigma=3.0, tolerance_d=0.03)
        via_phase = profile.guided_williamson_hall(
            phase=an, n_sigma=3.0, tolerance_d=0.03)
        # Both should yield the same crystallite_size to high tolerance
        assert np.isclose(manual['crystallite_size'],
                          via_phase['crystallite_size'], rtol=1e-10)

    def test_passing_both_phase_and_ref_d_raises_value_error(self):
        profile, an = self._make_profile_and_phase()
        with pytest.raises(ValueError, match='either ref_d or phase'):
            profile.guided_williamson_hall(
                ref_d=np.array([3.2, 4.0]), phase=an)

    def test_invalid_instrumental_type_raises_type_error(self):
        """Now that instrumental= is implemented for guided_williamson_hall,
        passing an unsupported type must raise TypeError, not
        NotImplementedError."""
        profile, an = self._make_profile_and_phase()
        with pytest.raises(TypeError, match='InstrumentalStandard'):
            profile.guided_williamson_hall(phase=an, instrumental='anything')


class TestXRDProfileGuidedWarrenAverbachPhaseKwarg:
    def _make_profile_and_phase(self):
        # Reuse the exact same setup as the W-H integration test
        return TestXRDProfileGuidedWilliamsonHallPhaseKwarg(
            )._make_profile_and_phase()

    def test_phase_kwarg_produces_same_result_as_manual_ref_peaks(self):
        profile, an = self._make_profile_and_phase()
        tt_range = (float(profile.two_theta.min()),
                    float(profile.two_theta.max()))
        manual_peaks = an.get_ref_peaks(profile.wavelength,
                                         two_theta_range=tt_range)
        manual = profile.guided_warren_averbach(
            manual_peaks, n_sigma=3.0, tolerance_d=0.03)
        via_phase = profile.guided_warren_averbach(
            phase=an, n_sigma=3.0, tolerance_d=0.03)
        assert manual['n_families'] == via_phase['n_families']

    def test_passing_both_phase_and_ref_peaks_raises_value_error(self):
        profile, an = self._make_profile_and_phase()
        with pytest.raises(ValueError, match='either ref_peaks or phase'):
            profile.guided_warren_averbach(
                ref_peaks=[{'d': 3.0, 'h': 1, 'k': 0, 'l': 0,
                            'intensity': 100, 'two_theta': 30}],
                phase=an)

    def test_instrumental_kwarg_raises_not_implemented(self):
        profile, an = self._make_profile_and_phase()
        with pytest.raises(NotImplementedError, match='Phase 2'):
            profile.guided_warren_averbach(phase=an, instrumental='x')


class TestExampleCIFsLoad:
    """Smoke test that all bundled example CIFs are parseable."""

    @pytest.mark.parametrize('cif_name', [
        'Forsterite.cif', 'Anorthite.cif', 'Pigeonite.cif',
        'Quartz.cif', 'Hematite.cif',
    ])
    def test_example_cif_loads(self, cif_name):
        from pathlib import Path
        from xrd_profile.phases import Phase
        cif_path = (Path(__file__).parent.parent / 'examples' /
                    'cifs' / cif_name)
        if not cif_path.exists():
            pytest.skip(f'{cif_name} not present in examples/cifs/')
        phase = Phase.from_cif(cif_path)
        assert phase.name == cif_name.replace('.cif', '')
        peaks = phase.get_ref_peaks(1.5406)
        assert len(peaks) > 0
