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
