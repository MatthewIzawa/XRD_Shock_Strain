"""
Regenerate tests/fixtures/golden_vX.Y.Z_results.json from the bundled
Tirhert subset using only the public API available at that tag.

Run when the numerical behavior at a given tag intentionally changes
(rare). Each regeneration must be accompanied by explicit reasoning
in the commit message.

Usage:
    python scripts/regenerate_goldens.py --tier v0.2.0
    python scripts/regenerate_goldens.py --tier v0.3.0
"""
import argparse
import json
from pathlib import Path
import numpy as np

from xrd_profile import (XRDProfile, two_theta_to_d,
                         scherrer, modified_scherrer,
                         compute_pdf_sine, fit_first_pdf_peak,
                         estimate_fwhm_simple, Phase)

LAMBDA_I11 = 0.826517
FIXTURE_DIR = Path(__file__).parent.parent / 'tests' / 'fixtures'
PATTERN_FILE = FIXTURE_DIR / 'tirhert_subset.xy'
ANORTHITE_CIF = (Path(__file__).parent.parent / 'examples'
                 / 'cifs' / 'Anorthite.cif')

# Fixed reference d-spacing list for guided W-H, derived from anorthite
# at I11 wavelength. Values frozen here so regeneration is deterministic.
ANORTHITE_REF_D = [
    3.20, 3.18, 3.65, 4.04, 6.41, 5.69, 3.74, 3.21, 4.04, 2.94,
]

# Fixed reference peak list for guided W-A: (d, two_theta, intensity, h, k, l).
ANORTHITE_REF_PEAKS = [
    {'d': 3.20, 'two_theta': 14.84, 'intensity': 100.0,
     'h': 0, 'k': 4, 'l': 0},
    {'d': 3.18, 'two_theta': 14.93, 'intensity':  85.0,
     'h': 2, 'k': 0, 'l': -2},
    {'d': 6.41, 'two_theta':  7.39, 'intensity':  60.0,
     'h': 0, 'k': 2, 'l': 0},
    {'d': 4.04, 'two_theta': 11.74, 'intensity':  50.0,
     'h': 0, 'k': 0, 'l': 2},
    {'d': 3.65, 'two_theta': 13.00, 'intensity':  45.0,
     'h': 1, 'k': 3, 'l': 0},
]


def to_serializable(obj):
    """Recursively convert numpy/tuple types to JSON-serialisable form.

    Dict keys that are tuples (e.g., (h, k, l) Miller indices that
    later v0.5+ goldens will carry) are stringified so json.dumps does
    not raise TypeError.
    """
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, (np.floating, np.integer)):
        return float(obj)
    if isinstance(obj, (bool, np.bool_)):
        return bool(obj)
    if isinstance(obj, dict):
        return {str(k) if isinstance(k, tuple) else k: to_serializable(v)
                for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [to_serializable(x) for x in obj]
    return obj


def regenerate_v020():
    """v0.2.0 tier: exercises the original array-based public API.
    Frozen reference for the v0.2.0 release tag."""
    data = np.loadtxt(PATTERN_FILE)
    tt, I = data[:, 0], data[:, 1]
    profile = XRDProfile(tt, I, wavelength=LAMBDA_I11,
                         sample_name='Tirhert_subset')

    results = {
        'metadata': {
            'fixture': str(PATTERN_FILE.name),
            'wavelength': LAMBDA_I11,
            'n_points': int(len(tt)),
            'tt_range': [float(tt.min()), float(tt.max())],
            'tier': 'v0.2.0',
        },
        'guided_williamson_hall': profile.guided_williamson_hall(
            np.array(ANORTHITE_REF_D), n_sigma=3.0, tolerance_d=0.03),
        'guided_warren_averbach': profile.guided_warren_averbach(
            ANORTHITE_REF_PEAKS, n_sigma=3.0, tolerance_d=0.03),
        'compute_pdf_sine': None,
        'scherrer_default': None,
        'modified_scherrer_default': None,
    }
    r, G_r, Q_max = compute_pdf_sine(tt, I, LAMBDA_I11,
                                      cheby_order=15, lorch=True)
    results['compute_pdf_sine'] = {
        'r': r.tolist(), 'G_r': G_r.tolist(), 'Q_max': float(Q_max),
    }
    fwhm, positions = estimate_fwhm_simple(tt, I, height_threshold=0.05)
    if len(fwhm) > 0:
        sizes = scherrer(fwhm, positions, LAMBDA_I11)
        results['scherrer_default'] = {
            'sizes': sizes.tolist(),
            'mean_size': float(np.mean(sizes)),
            'n_peaks': int(len(sizes)),
        }
        if len(fwhm) >= 2:
            mod_size = modified_scherrer(fwhm, positions, LAMBDA_I11)
            results['modified_scherrer_default'] = float(mod_size)

    out = FIXTURE_DIR / 'golden_v0.2.0_results.json'
    out.write_text(json.dumps(to_serializable(results), indent=2))
    print(f'Wrote {out}')


def regenerate_v030():
    """v0.3.0 tier: exercises Phase API, run_all, Scherrer shape table.
    All v0.4+ kwargs left at their defaults (instrumental=None, etc.)."""
    data = np.loadtxt(PATTERN_FILE)
    tt, I = data[:, 0], data[:, 1]
    profile = XRDProfile(tt, I, wavelength=LAMBDA_I11,
                         sample_name='Tirhert_subset')
    anorthite = Phase.from_cif(str(ANORTHITE_CIF), name='anorthite')

    results = {
        'metadata': {
            'fixture': str(PATTERN_FILE.name),
            'wavelength': LAMBDA_I11,
            'n_points': int(len(tt)),
            'tt_range': [float(tt.min()), float(tt.max())],
            'tier': 'v0.3.0',
        },
        'guided_wh_via_phase': profile.guided_williamson_hall(
            phase=anorthite, n_sigma=3.0, tolerance_d=0.03),
        'guided_wa_via_phase': profile.guided_warren_averbach(
            phase=anorthite, n_sigma=3.0, tolerance_d=0.03),
        'scherrer_spherical': profile.scherrer(shape='spherical'),
        'scherrer_cubic': profile.scherrer(shape='cubic'),
        'scherrer_default': profile.scherrer(),
        'modified_scherrer_default': profile.modified_scherrer(),
        'run_all_no_phases': profile.run_all(
            methods=['wh', 'scherrer']),
    }

    out = FIXTURE_DIR / 'golden_v0.3.0_results.json'
    out.write_text(json.dumps(to_serializable(results), indent=2))
    print(f'Wrote {out}')


def main():
    parser = argparse.ArgumentParser(
        description='Regenerate a versioned golden fixture for '
                    'tests/test_backward_compat.py. Requires explicit '
                    'reasoning in the commit message because the '
                    'fixture is treated as a frozen reference for the '
                    'corresponding tag.')
    parser.add_argument(
        '--tier', choices=['v0.2.0', 'v0.3.0'], default='v0.2.0',
        help='which fixture to overwrite (default: v0.2.0)')
    args = parser.parse_args()
    if args.tier == 'v0.2.0':
        regenerate_v020()
    elif args.tier == 'v0.3.0':
        regenerate_v030()


if __name__ == '__main__':
    main()
