"""Unit tests for XRDProfile.run_all dispatcher (v0.3.0 + v0.4.0)."""
from pathlib import Path

import numpy as np
import pytest


pytest.importorskip('pymatgen')


FIXTURE = Path(__file__).parent / 'fixtures' / 'tirhert_subset.xy'
LAMBDA_I11 = 0.826517


def _profile_and_anorthite():
    from xrd_profile import XRDProfile
    from xrd_profile.phases import Phase
    data = np.loadtxt(FIXTURE)
    profile = XRDProfile(data[:, 0], data[:, 1],
                          wavelength=LAMBDA_I11, sample_name='Tirhert')
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


class TestRunAllDispatch:
    def test_methods_subset_only_runs_subset(self):
        profile, an = _profile_and_anorthite()
        result = profile.run_all(methods=['pdf'], phases=None)
        assert set(result.keys()) == {'pdf'}

    def test_methods_none_runs_all_four(self):
        profile, an = _profile_and_anorthite()
        result = profile.run_all(methods=None, phases=[an],
                                  wh={'tolerance_d': 0.03},
                                  wa={'tolerance_d': 0.03})
        assert set(result.keys()) == {'wh', 'wa', 'pdf', 'scherrer'}


class TestRunAllPerPhaseKeys:
    def test_wh_result_keyed_by_phase_name_when_phases_given(self):
        profile, an = _profile_and_anorthite()
        result = profile.run_all(methods=['wh'], phases=[an],
                                  wh={'tolerance_d': 0.03})
        assert 'wh' in result
        assert 'Anorthite' in result['wh']

    def test_wh_result_is_flat_dict_when_no_phases(self):
        profile, an = _profile_and_anorthite()
        result = profile.run_all(methods=['wh'], phases=None)
        # unguided W-H returns a dict (no phase wrapping)
        assert 'wh' in result
        # The unguided result is the raw dict from williamson_hall(),
        # not a phase-keyed sub-dict.
        assert 'crystallite_size' in result['wh'] or 'x' in result['wh']

    def test_single_phase_accepted_as_list_or_scalar(self):
        profile, an = _profile_and_anorthite()
        a = profile.run_all(methods=['wh'], phases=an,
                             wh={'tolerance_d': 0.03})
        b = profile.run_all(methods=['wh'], phases=[an],
                             wh={'tolerance_d': 0.03})
        # Both yield the same nested structure
        assert 'Anorthite' in a['wh']
        assert 'Anorthite' in b['wh']


class TestRunAllInstrumental:
    @pytest.fixture(scope='class')
    def anorthite_phase(self):
        from xrd_profile import Phase
        cif = (Path(__file__).parent.parent / 'examples'
               / 'cifs' / 'Anorthite.cif')
        return Phase.from_cif(str(cif), name='anorthite')

    def test_run_all_with_instrumental_profile_runs(
            self, anorthite_phase):
        from xrd_profile import (XRDProfile, InstrumentalProfile)
        FIX = Path(__file__).parent / 'fixtures' / 'tirhert_subset.xy'
        data = np.loadtxt(FIX)
        profile = XRDProfile(data[:, 0], data[:, 1],
                              wavelength=0.826517)
        # Synchrotron-scale instrumental (lab-scale would over-correct).
        inst = InstrumentalProfile(U=5e-5, V=-2e-5, W=1e-4,
                                    wavelength=0.826517)
        # methods=['wh', 'scherrer'] only - W-A would reject Profile.
        result = profile.run_all(
            methods=['wh', 'scherrer'],
            phases=anorthite_phase,
            instrumental=inst)
        assert 'wh' in result
        assert 'scherrer' in result
        # Per-phase dict shape: scherrer keyed by phase name when
        # phases= is given (mirrors W-H per-phase shape).
        assert 'anorthite' in result['scherrer']

    def test_run_all_with_instrumental_profile_rejects_wa(
            self, anorthite_phase):
        from xrd_profile import (XRDProfile, InstrumentalProfile)
        FIX = Path(__file__).parent / 'fixtures' / 'tirhert_subset.xy'
        data = np.loadtxt(FIX)
        profile = XRDProfile(data[:, 0], data[:, 1],
                              wavelength=0.826517)
        inst = InstrumentalProfile(U=5e-5, V=-2e-5, W=1e-4,
                                    wavelength=0.826517)
        with pytest.raises(ValueError, match='Stokes'):
            profile.run_all(methods=['wa'], phases=anorthite_phase,
                             instrumental=inst)
