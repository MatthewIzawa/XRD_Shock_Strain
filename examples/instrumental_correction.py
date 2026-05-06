"""
Demonstrate v0.4.0 instrumental broadening correction.

Loads the synthetic LaB6 fixture as the standard, the bundled Tirhert
subset as the sample, and runs guided Warren-Averbach with and without
Stokes Fourier deconvolution. Prints the apparent crystallite size in
each case.

Run from the repository root:
    python examples/instrumental_correction.py
"""
from pathlib import Path
import sys

import numpy as np

REPO_ROOT = Path(__file__).parent.parent
# Ensure the in-tree xrd_profile is preferred over any site-packages
# install (so v0.4.0 features are picked up when running directly
# from the source tree).
sys.path.insert(0, str(REPO_ROOT))

from xrd_profile import (XRDProfile, Phase,             # noqa: E402
                         InstrumentalStandard)


SAMPLE_PATTERN = REPO_ROOT / 'tests' / 'fixtures' / 'tirhert_subset.xy'
LAB6_PATTERN = REPO_ROOT / 'tests' / 'fixtures' / 'synthetic_lab6.xy'
ANORTHITE_CIF = REPO_ROOT / 'examples' / 'cifs' / 'Anorthite.cif'

LAMBDA_I11 = 0.826517


def main():
    sample_data = np.loadtxt(SAMPLE_PATTERN)
    lab6_data = np.loadtxt(LAB6_PATTERN)

    # Sample profile.
    profile = XRDProfile(sample_data[:, 0], sample_data[:, 1],
                          wavelength=LAMBDA_I11,
                          sample_name='Tirhert_subset')
    anorthite = Phase.from_cif(str(ANORTHITE_CIF), name='anorthite')

    # Instrumental standard. (Synthetic LaB6, used here as if measured
    # at the I11 wavelength so the demo is self-contained.)
    lab6_phase = Phase.from_lattice_params(
        a=4.156825, b=4.156825, c=4.156825,
        alpha=90, beta=90, gamma=90,
        species=['La', 'B', 'B', 'B', 'B', 'B', 'B'],
        coords=[[0.0, 0.0, 0.0],
                [0.5, 0.5, 0.1993], [0.5, 0.5, 0.8007],
                [0.5, 0.1993, 0.5], [0.5, 0.8007, 0.5],
                [0.1993, 0.5, 0.5], [0.8007, 0.5, 0.5]],
        name='LaB6')
    lab6 = InstrumentalStandard(
        phase=lab6_phase,
        two_theta=lab6_data[:, 0], intensity=lab6_data[:, 1],
        wavelength=LAMBDA_I11, name='synthetic_lab6')

    # Without correction.
    wa_uncorrected = profile.guided_warren_averbach(
        phase=anorthite, n_sigma=3.0, tolerance_d=0.03)

    # With Stokes Fourier deconvolution.
    wa_corrected = profile.guided_warren_averbach(
        phase=anorthite, instrumental=lab6,
        n_sigma=3.0, tolerance_d=0.03)

    print('=== Anorthite W-A on Tirhert_subset ===')
    print(f'Uncorrected mean crystallite size: '
          f'{wa_uncorrected["mean_crystallite_size"]:.1f} angstroms')
    print(f'  (instrumental contribution included in the broadening)')
    print(f'Corrected mean crystallite size:   '
          f'{wa_corrected["mean_crystallite_size"]:.1f} angstroms')
    print(f'  (Stokes-deconvolved against synthetic LaB6 standard)')
    print()
    print(f'Per-family size-distribution fits (lognormal D_median):')
    for fam in wa_corrected['families']:
        sd = fam['size_distribution']
        if sd is None:
            print(f'  hkl {fam["base_hkl"]}: insufficient L points')
        else:
            ln = sd['lognormal']
            print(f'  hkl {fam["base_hkl"]}: '
                  f'D_median={ln["D_median"]:.1f} A, '
                  f'sigma={ln["sigma"]:.2f}, '
                  f'R2={ln["fit_r2"]:.3f}')


if __name__ == '__main__':
    main()
