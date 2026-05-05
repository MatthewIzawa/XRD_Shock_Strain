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
