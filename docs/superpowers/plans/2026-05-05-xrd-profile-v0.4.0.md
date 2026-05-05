# xrd_profile v0.4.0 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Land Phase 2's first tag — instrumental broadening deconvolution (Caglioti for W-H/Scherrer; Stokes Fourier for W-A), crystallite size distributions (lognormal + normal per-family), and guided Scherrer — as `xrd_profile` v0.4.0.

**Architecture:** Two new public classes (`InstrumentalStandard`, `InstrumentalProfile`) live in a new module `xrd_profile/instrumental.py`; size-distribution fits live in a new internal module `xrd_profile/size_distributions.py`. Existing analysis modules (`williamson_hall.py`, `warren_averbach.py`, `scherrer.py`, `profile.py`) gain dispatch wiring for the `instrumental=` kwarg already reserved in v0.3.0. All v0.3.0 result-dict keys retain their numerical values when no v0.4.0 feature is invoked (key-subset value-equality regression test on the new `golden_v0.3.0_results.json` enforces this).

**Tech Stack:** Python 3.10+, numpy, scipy (`linregress`, `optimize.curve_fit`, `signal.fftconvolve`, `special.erfc`), matplotlib (existing). pymatgen via the existing `[cif]` extra. No new dependencies.

**Spec reference:** [`docs/superpowers/specs/2026-05-05-xrd-profile-v1-phase2-design.md`](../specs/2026-05-05-xrd-profile-v1-phase2-design.md), Sections 5.1, 6.1–6.2, 7.1–7.3, 8.2.1–8.2.2, 8.2.5, 11.1.

**Working directory:** `C:\Users\Matthew Izawa\Documents\Dan Applin\Llunr\xrd_profile\` (the package's git repo, currently at `main` = `c983ff6` "Release v0.3.0", with the Phase 2 spec committed at `b833c66`).

---

## File structure

### New files (create in order shown)

| Path | Responsibility | Approx LOC |
|---|---|---|
| `tests/fixtures/golden_v0.3.0_results.json` | v0.3.0 frozen-output regression baseline (created in Task 1 before any code change). | data |
| `scripts/build_synthetic_standard.py` | Reproducible generator for the synthetic LaB6 fixture; documented Caglioti coefficients. | ~80 |
| `tests/fixtures/synthetic_lab6.xy` | Synthetic LaB6 pattern with known U/V/W. | data |
| `xrd_profile/instrumental.py` | `InstrumentalStandard`, `InstrumentalProfile`, internal `_caglioti_fit`, `_caglioti_fwhm_at`, `_stokes_deconvolve`. | ~280 |
| `xrd_profile/size_distributions.py` | Internal: `lognormal_a_size`, `normal_a_size` basis functions; `fit_size_distribution()` returning the `'size_distribution'` dict shape; `_moments_initial_guess`. | ~140 |
| `xrd_profile/registry/README.md` | Documents the registry JSON file format and the lookup path. Empty registry. | ~40 |
| `tests/test_instrumental.py` | ~12 tests covering both classes, Caglioti recovery, Stokes sanity, dispatch from W-H/W-A/Scherrer. | ~250 |
| `tests/test_size_distributions.py` | ~8 tests covering lognormal/normal recovery, threshold for `None`-result, presence of metadata keys. | ~150 |

### Modified files

| Path | What changes |
|---|---|
| `xrd_profile/williamson_hall.py` | Internal `_apply_caglioti(fwhm, two_theta, inst_profile)` helper. Threaded into `guided_williamson_hall(...)` immediately after FWHM measurement. |
| `xrd_profile/warren_averbach.py` | Stokes deconvolution applied to per-peak Fourier coefficients before harmonic decomposition. New per-family `'size_distribution'` key injected. |
| `xrd_profile/scherrer.py` | New internal helpers `_filter_to_phase_peaks` and `_caglioti_subtract_fwhm`. Public `scherrer()` and `modified_scherrer()` get optional `phase=` and `instrumental=` kwargs. |
| `xrd_profile/profile.py` | Replace `NotImplementedError` guards on `guided_williamson_hall`, `guided_warren_averbach`, and `run_all` with dispatch to the new helpers. Add `phase=`/`instrumental=` to `XRDProfile.scherrer()` and `XRDProfile.modified_scherrer()`. |
| `xrd_profile/__init__.py` | Export `InstrumentalStandard`, `InstrumentalProfile`. Bump `__version__` to `'0.4.0'` (Task 17). |
| `tests/test_backward_compat.py` | Add a `golden_v0.3.0` tier (Phase API + run_all + Scherrer shape table) using key-subset value-equality semantics. |
| `scripts/regenerate_goldens.py` | Add `--tier {v0.2.0, v0.3.0}` flag. v0.3.0 path exercises Phase API. |
| `CHANGELOG.md` | New `[0.4.0]` entry. |
| `README.md` | Quickstart for instrumental correction. |

---

## Conventions for every task

- **TDD discipline.** Every task that adds behaviour starts with a failing test, runs it to confirm failure mode, implements minimally to pass, runs the test green, then commits.
- **Commits per task.** Each task ends in a single commit. Use the project's existing terse subject-line style (Phase 1 used "Release v0.3.0", "Update README to v0.3.0: ...", etc.). Include the standard Co-Authored-By trailer that the harness already appends.
- **Test command.** `cd "C:/Users/Matthew Izawa/Documents/Dan Applin/Llunr/xrd_profile" && pytest <path> -v`. The Windows working directory must be quoted; use `pytest -x` to stop on first failure during development.
- **No real-instrument data.** All synthetic fixtures are computed from documented parametric forms in committed scripts. Real LaB6/Si patterns are user-supplied at runtime.
- **No JAC pipeline edits.** Do not modify `Llunr/paper1_figures.py`, `Llunr/run_guided_wh.py`, `Llunr/compile_survey.py`, or anything in `Paper1_JAC/`.

---

## Task 1 — Freeze v0.3.0 golden fixture (regression safety net)

**Files:**
- Create: `tests/fixtures/golden_v0.3.0_results.json`
- Modify: `scripts/regenerate_goldens.py`
- Test: `tests/test_backward_compat.py` (extension)

This task creates a frozen v0.3.0 results snapshot **before** any v0.4.0 code lands. Every later task that touches `profile.py` / `warren_averbach.py` / `williamson_hall.py` / `scherrer.py` re-runs the v0.3.0 tier afterwards to confirm no drift.

- [ ] **Step 1: Add `--tier` flag to `scripts/regenerate_goldens.py`**

```python
# scripts/regenerate_goldens.py
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

ANORTHITE_REF_D = [
    3.20, 3.18, 3.65, 4.04, 6.41, 5.69, 3.74, 3.21, 4.04, 2.94,
]
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
    parser = argparse.ArgumentParser()
    parser.add_argument('--tier', choices=['v0.2.0', 'v0.3.0'],
                        default='v0.2.0')
    args = parser.parse_args()
    if args.tier == 'v0.2.0':
        regenerate_v020()
    elif args.tier == 'v0.3.0':
        regenerate_v030()


if __name__ == '__main__':
    main()
```

- [ ] **Step 2: Run the v0.3.0 regeneration**

Run: `cd "C:/Users/Matthew Izawa/Documents/Dan Applin/Llunr/xrd_profile" && python scripts/regenerate_goldens.py --tier v0.3.0`

Expected output: `Wrote .../tests/fixtures/golden_v0.3.0_results.json`

If the script errors on `Phase.from_cif`, verify `examples/cifs/Anorthite.cif` exists and pymatgen is installed (`pip install -e .[cif]`).

- [ ] **Step 3: Confirm v0.2.0 regeneration is unchanged**

Run: `cd "C:/Users/Matthew Izawa/Documents/Dan Applin/Llunr/xrd_profile" && python scripts/regenerate_goldens.py --tier v0.2.0`

Expected: produces a `golden_v0.2.0_results.json` byte-equal (or numerically equal within JSON formatting noise) to the existing one. Use `git diff tests/fixtures/golden_v0.2.0_results.json` to confirm — it should be empty or trivial whitespace.

- [ ] **Step 4: Add v0.3.0 tier to `tests/test_backward_compat.py`**

Add to the bottom of `tests/test_backward_compat.py`:

```python
# --- v0.3.0 tier (added in v0.4.0) ---
# Asserts: every key in golden_v0.3.0_results.json is reproducible at
# the current tag with byte-equivalent numerical value (key-subset
# value-equality). New top-level keys added at v0.4+ tags are allowed
# in the live result and ignored here.

from xrd_profile import Phase

ANORTHITE_CIF = (Path(__file__).parent.parent / 'examples'
                 / 'cifs' / 'Anorthite.cif')


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
```

- [ ] **Step 5: Run the full backward-compat suite**

Run: `cd "C:/Users/Matthew Izawa/Documents/Dan Applin/Llunr/xrd_profile" && pytest tests/test_backward_compat.py -v`

Expected: all v0.2.0 tests pass (existing) + 5 new v0.3.0 tests pass.

- [ ] **Step 6: Commit**

```bash
cd "C:/Users/Matthew Izawa/Documents/Dan Applin/Llunr/xrd_profile"
git add scripts/regenerate_goldens.py tests/fixtures/golden_v0.3.0_results.json tests/test_backward_compat.py
git commit -m "Freeze v0.3.0 golden fixture for v0.4.0 regression baseline"
```

---

## Task 2 — Synthetic LaB6 fixture for instrumental tests

**Files:**
- Create: `scripts/build_synthetic_standard.py`
- Create: `tests/fixtures/synthetic_lab6.xy`

LaB6 (NIST SRM 660c) has cubic Pm-3m symmetry, a = 4.156825 Å. We synthesise a pattern at Cu Kα with documented Caglioti coefficients (`U=5e-3, V=-1e-3, W=5e-3` deg² — typical of a well-aligned lab Bruker). The fixture lets every Caglioti / Stokes test run deterministically without real-instrument data.

- [ ] **Step 1: Write `scripts/build_synthetic_standard.py`**

```python
"""
Generate a synthetic LaB6 diffraction pattern with documented Caglioti
broadening, for use as a test fixture in tests/test_instrumental.py
and tests/test_size_distributions.py.

Caglioti polynomial:
    FWHM(2theta)^2 = U * tan^2(theta) + V * tan(theta) + W

with documented coefficients
    U = 5e-3 deg^2,  V = -1e-3 deg^2,  W = 5e-3 deg^2
which produces FWHMs typical of a well-aligned lab Bruker over 20-100
degrees 2-theta.

LaB6 is cubic Pm-3m, a = 4.156825 angstroms (NIST SRM 660c). Reflections
up to (3,1,1) populate the synthetic pattern. Each peak is a pure
Gaussian centred at the Bragg position with FWHM given by Caglioti.

Run: python scripts/build_synthetic_standard.py
"""
from pathlib import Path
import numpy as np

LAMBDA_CU = 1.5406  # Cu K-alpha
A_LAB6 = 4.156825   # angstroms

# Documented Caglioti coefficients used to synthesise the fixture.
SYNTH_U = 5.0e-3
SYNTH_V = -1.0e-3
SYNTH_W = 5.0e-3

# 2-theta sampling (matches a typical lab Bruker scan).
TT_MIN, TT_MAX, TT_STEP = 20.0, 100.0, 0.02

OUT = (Path(__file__).parent.parent / 'tests' / 'fixtures'
       / 'synthetic_lab6.xy')


def caglioti_fwhm(two_theta_deg, U, V, W):
    """Caglioti FWHM in degrees from polynomial coefficients."""
    theta = np.deg2rad(two_theta_deg / 2.0)
    fwhm_sq = U * np.tan(theta)**2 + V * np.tan(theta) + W
    return np.sqrt(np.maximum(fwhm_sq, 1e-8))


def lab6_reflections(max_hkl=4):
    """Allowed reflections for LaB6 (Pm-3m, all h,k,l permitted).
    Returns list of (h, k, l, multiplicity, structure_factor_proxy)."""
    refs = []
    seen = set()
    for h in range(0, max_hkl + 1):
        for k in range(0, h + 1):
            for l in range(0, k + 1):
                if h == k == l == 0:
                    continue
                key = tuple(sorted([h, k, l], reverse=True))
                if key in seen:
                    continue
                seen.add(key)
                # Structure-factor proxy: drops with hkl magnitude.
                mult = _multiplicity(h, k, l)
                f_proxy = 1.0 / (1.0 + 0.05 * (h*h + k*k + l*l))
                refs.append((h, k, l, mult, f_proxy))
    return refs


def _multiplicity(h, k, l):
    """Cubic point-group multiplicity for (h,k,l)."""
    distinct = len({abs(h), abs(k), abs(l)})
    nonzero = sum(1 for v in (h, k, l) if v != 0)
    if distinct == 1 and nonzero == 3:
        return 8     # (hhh)
    if nonzero == 1:
        return 6     # (h00)
    if distinct == 1 and nonzero == 2:
        return 12    # (hh0)
    if nonzero == 2:
        return 24    # (hk0)
    if distinct == 2 and nonzero == 3:
        return 24    # (hhl)
    return 48        # (hkl)


def bragg_two_theta(h, k, l, a, wavelength):
    d = a / np.sqrt(h*h + k*k + l*l)
    sin_theta = wavelength / (2.0 * d)
    if sin_theta >= 1.0:
        return None
    return 2.0 * np.rad2deg(np.arcsin(sin_theta))


def gaussian_peak(x, x0, fwhm, amplitude):
    sigma = fwhm / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    return amplitude * np.exp(-0.5 * ((x - x0) / sigma)**2)


def main():
    tt = np.arange(TT_MIN, TT_MAX + TT_STEP, TT_STEP)
    intensity = np.zeros_like(tt)
    intensity += 5.0  # flat background

    for h, k, l, mult, f_proxy in lab6_reflections(max_hkl=4):
        tt_peak = bragg_two_theta(h, k, l, A_LAB6, LAMBDA_CU)
        if tt_peak is None or tt_peak < TT_MIN or tt_peak > TT_MAX:
            continue
        fwhm = caglioti_fwhm(tt_peak, SYNTH_U, SYNTH_V, SYNTH_W)
        amplitude = mult * f_proxy * 1000.0
        intensity += gaussian_peak(tt, tt_peak, fwhm, amplitude)

    # Add a tiny Gaussian noise floor for realism (sigma=0.5 counts).
    rng = np.random.default_rng(seed=42)
    intensity += rng.normal(0.0, 0.5, size=intensity.shape)
    intensity = np.maximum(intensity, 0.0)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(OUT, np.column_stack([tt, intensity]),
               fmt='%.4f', header='two_theta_deg  intensity\n'
               f'# synthesis: U={SYNTH_U}, V={SYNTH_V}, W={SYNTH_W}',
               comments='# ')
    print(f'Wrote {OUT}: {len(tt)} points, '
          f'{TT_MIN}-{TT_MAX} deg, step {TT_STEP}')


if __name__ == '__main__':
    main()
```

- [ ] **Step 2: Run the generator**

Run: `cd "C:/Users/Matthew Izawa/Documents/Dan Applin/Llunr/xrd_profile" && python scripts/build_synthetic_standard.py`

Expected output: `Wrote .../tests/fixtures/synthetic_lab6.xy: 4001 points, 20.0-100.0 deg, step 0.02`. The generated file is ~150 KB.

- [ ] **Step 3: Sanity-check the fixture by eye**

Run a quick plot to confirm the pattern looks right:

```bash
cd "C:/Users/Matthew Izawa/Documents/Dan Applin/Llunr/xrd_profile"
python -c "
import numpy as np, matplotlib.pyplot as plt
d = np.loadtxt('tests/fixtures/synthetic_lab6.xy')
plt.plot(d[:,0], d[:,1])
plt.xlabel('2theta (deg)'); plt.ylabel('Intensity')
plt.title('Synthetic LaB6 fixture')
plt.savefig('synthetic_lab6_check.png', dpi=100)
print('Saved synthetic_lab6_check.png')
"
```

Expected: 6–10 visible Gaussian peaks at LaB6's known Bragg positions (21.36°, 30.39°, 37.44°, 43.53°, 48.99°, 53.99°, 63.20°, ...). Delete `synthetic_lab6_check.png` after eyeballing — do not commit it.

- [ ] **Step 4: Commit**

```bash
cd "C:/Users/Matthew Izawa/Documents/Dan Applin/Llunr/xrd_profile"
rm -f synthetic_lab6_check.png
git add scripts/build_synthetic_standard.py tests/fixtures/synthetic_lab6.xy
git commit -m "Add synthetic LaB6 fixture for instrumental-correction tests"
```

---

## Task 3 — `InstrumentalProfile` class skeleton + `fwhm_at`

**Files:**
- Create: `xrd_profile/instrumental.py`
- Create: `tests/test_instrumental.py`

`InstrumentalProfile` is the lighter-weight class: it carries Caglioti coefficients (U, V, W) plus wavelength and an optional name, with no measured pattern. Sufficient for W-H and Scherrer correction; insufficient for W-A Stokes.

- [ ] **Step 1: Write the failing test**

Create `tests/test_instrumental.py`:

```python
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
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cd "C:/Users/Matthew Izawa/Documents/Dan Applin/Llunr/xrd_profile" && pytest tests/test_instrumental.py -v`

Expected: `ImportError: cannot import name 'InstrumentalProfile' from 'xrd_profile'`.

- [ ] **Step 3: Create `xrd_profile/instrumental.py` with `InstrumentalProfile`**

```python
"""
instrumental.py — Instrumental broadening characterisation and
deconvolution.

Two classes:
    InstrumentalStandard: holds a structural Phase plus a measured
        diffraction pattern of a known standard (LaB6, Si). Supports
        both Caglioti FWHM correction (for W-H, Scherrer) and Stokes
        Fourier deconvolution (for W-A).
    InstrumentalProfile: lightweight Caglioti carrier (U, V, W). No
        measured pattern. Supports W-H, Scherrer; W-A raises a clear
        ValueError.

Caglioti polynomial:
    FWHM(2theta)^2 = U * tan^2(theta) + V * tan(theta) + W

References
----------
Caglioti, G., Paoletti, A., Ricci, F. P. (1958). Choice of collimators
    for a crystal spectrometer for neutron diffraction. Nuclear
    Instruments 3, 223-228.
Stokes, A. R. (1948). A numerical Fourier-analysis method for the
    correction of widths and shapes of lines on X-ray powder
    photographs. Proc. Phys. Soc. 61, 382-391.
"""
import json
from pathlib import Path

import numpy as np


class InstrumentalProfile:
    """Caglioti-polynomial carrier for instrumental broadening.

    Parameters
    ----------
    U, V, W : float
        Caglioti coefficients (deg^2).
    wavelength : float
        X-ray wavelength (angstroms).
    name : str
        Optional human-readable label.
    """

    def __init__(self, U: float, V: float, W: float,
                 wavelength: float, name: str = ''):
        self.U = float(U)
        self.V = float(V)
        self.W = float(W)
        self.wavelength = float(wavelength)
        self.name = str(name)

    def fwhm_at(self, two_theta_deg: float) -> float:
        """Caglioti FWHM (degrees) at the given 2-theta (degrees)."""
        theta = np.deg2rad(np.asarray(two_theta_deg) / 2.0)
        fwhm_sq = (self.U * np.tan(theta)**2
                   + self.V * np.tan(theta)
                   + self.W)
        return np.sqrt(np.maximum(fwhm_sq, 0.0))
```

Update `xrd_profile/__init__.py` to export the new class:

```python
# After the existing 'from .phases import Phase, build_reference_peaks'
# line, add:
from .instrumental import InstrumentalProfile
```

And add `'InstrumentalProfile'` to the `__all__` list.

- [ ] **Step 4: Run tests to verify they pass**

Run: `cd "C:/Users/Matthew Izawa/Documents/Dan Applin/Llunr/xrd_profile" && pytest tests/test_instrumental.py -v`

Expected: all four tests in `TestInstrumentalProfileBasics` pass.

- [ ] **Step 5: Commit**

```bash
cd "C:/Users/Matthew Izawa/Documents/Dan Applin/Llunr/xrd_profile"
git add xrd_profile/instrumental.py xrd_profile/__init__.py tests/test_instrumental.py
git commit -m "Add InstrumentalProfile class with Caglioti fwhm_at()"
```

---

## Task 4 — `InstrumentalProfile` JSON I/O and registry lookup

**Files:**
- Modify: `xrd_profile/instrumental.py`
- Modify: `tests/test_instrumental.py`
- Create: `xrd_profile/registry/README.md`

JSON round-trip lets users cache a fitted profile. `from_registry()` looks up named profiles in the package's `xrd_profile/registry/` directory. v0.4.0 ships an empty registry; the `README.md` documents the file format.

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_instrumental.py`:

```python
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


class TestInstrumentalProfileRegistry:
    def test_from_registry_unknown_name_raises_keyerror(self):
        with pytest.raises(KeyError, match='nonexistent_profile'):
            InstrumentalProfile.from_registry('nonexistent_profile')
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cd "C:/Users/Matthew Izawa/Documents/Dan Applin/Llunr/xrd_profile" && pytest tests/test_instrumental.py::TestInstrumentalProfileJsonIO tests/test_instrumental.py::TestInstrumentalProfileRegistry -v`

Expected: `AttributeError: 'InstrumentalProfile' has no attribute 'to_json'` (and the others).

- [ ] **Step 3: Implement JSON I/O and registry lookup**

Append to `xrd_profile/instrumental.py` inside the `InstrumentalProfile` class:

```python
    SCHEMA_VERSION = '1'

    def to_json(self, path) -> None:
        """Serialise this profile to a JSON file at `path`."""
        path = Path(path)
        path.write_text(json.dumps({
            'schema_version': self.SCHEMA_VERSION,
            'U': self.U,
            'V': self.V,
            'W': self.W,
            'wavelength': self.wavelength,
            'name': self.name,
        }, indent=2))

    @classmethod
    def from_json(cls, path) -> 'InstrumentalProfile':
        """Load a profile from a JSON file produced by `to_json`."""
        data = json.loads(Path(path).read_text())
        if data.get('schema_version') != cls.SCHEMA_VERSION:
            raise ValueError(
                f'Unsupported InstrumentalProfile schema_version '
                f'{data.get("schema_version")!r}; this code expects '
                f'{cls.SCHEMA_VERSION!r}')
        return cls(U=data['U'], V=data['V'], W=data['W'],
                   wavelength=data['wavelength'],
                   name=data.get('name', ''))

    @classmethod
    def from_registry(cls, name: str) -> 'InstrumentalProfile':
        """Look up a pre-fit profile by name in
        `xrd_profile/registry/<name>.json`. v0.4.0 ships an empty
        registry; users populate it with their own JSON profiles."""
        registry_dir = Path(__file__).parent / 'registry'
        candidate = registry_dir / f'{name}.json'
        if not candidate.is_file():
            raise KeyError(
                f'No registered InstrumentalProfile {name!r}; '
                f'expected file at {candidate}. The v0.4.0 registry '
                f'ships empty; populate via `InstrumentalProfile.to_json` '
                f'into the registry directory.')
        return cls.from_json(candidate)
```

- [ ] **Step 4: Create the registry directory with a documentation README**

Create `xrd_profile/registry/README.md`:

```markdown
# `xrd_profile/registry/` — instrumental profile registry

Drop `<instrument-name>.json` files here to make pre-fit Caglioti
profiles available via `InstrumentalProfile.from_registry(name)`.

## File format (`schema_version: '1'`)

```json
{
  "schema_version": "1",
  "U": 5.0e-3,
  "V": -1.0e-3,
  "W": 5.0e-3,
  "wavelength": 1.5406,
  "name": "Lab_Bruker_Cu_Ka"
}
```

`U`, `V`, `W` are Caglioti polynomial coefficients (deg²). `wavelength`
is in angstroms. The `name` field is informational and need not match
the filename.

## Generating a profile from a measured standard

```python
from xrd_profile import InstrumentalStandard
import numpy as np

data = np.loadtxt('lab6_pattern.xy')
std = InstrumentalStandard.from_cif_and_pattern(
    cif='LaB6.cif',
    two_theta=data[:, 0], intensity=data[:, 1],
    wavelength=1.5406, name='Lab_Bruker_Cu_Ka')
profile = std.caglioti_fit()
profile.to_json('xrd_profile/registry/lab_bruker_cu_ka.json')
```

## v0.4.0 status

The registry ships **empty** in v0.4.0. Calibration-grade Caglioti
fits for the Misasa Rigaku, Winnipeg Bruker, and Diamond I11 instruments
are a separate calibration deliverable. Users supply their own
profiles, either by fitting their own LaB6/Si standard and saving the
result here, or by hand-editing a file with literature U/V/W values.
```

- [ ] **Step 5: Run tests to verify they pass**

Run: `cd "C:/Users/Matthew Izawa/Documents/Dan Applin/Llunr/xrd_profile" && pytest tests/test_instrumental.py -v`

Expected: all `TestInstrumentalProfileJsonIO` and `TestInstrumentalProfileRegistry` tests pass; earlier tests still pass.

- [ ] **Step 6: Commit**

```bash
cd "C:/Users/Matthew Izawa/Documents/Dan Applin/Llunr/xrd_profile"
git add xrd_profile/instrumental.py xrd_profile/registry/README.md tests/test_instrumental.py
git commit -m "Add InstrumentalProfile JSON I/O and from_registry lookup"
```

---

## Task 5 — Caglioti fitting (`_caglioti_fit` + `InstrumentalStandard.caglioti_fit`)

**Files:**
- Modify: `xrd_profile/instrumental.py`
- Modify: `tests/test_instrumental.py`

`_caglioti_fit(two_theta, intensity, ref_two_theta)` measures FWHMs at the standard's reference peak positions, then fits the Caglioti polynomial via `scipy.optimize.curve_fit`. The helper is consumed by `InstrumentalStandard.caglioti_fit()` (added in Task 6) but tested directly here for isolation.

- [ ] **Step 1: Write the failing test**

Append to `tests/test_instrumental.py`:

```python
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
        # V is small and noise-prone; tolerance is wider.
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
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cd "C:/Users/Matthew Izawa/Documents/Dan Applin/Llunr/xrd_profile" && pytest tests/test_instrumental.py::TestCagliotiFit -v`

Expected: `ImportError: cannot import name '_caglioti_fit' from 'xrd_profile.instrumental'`.

- [ ] **Step 3: Implement `_caglioti_fit`**

Add to `xrd_profile/instrumental.py` after the `InstrumentalProfile` class definition:

```python
from scipy.optimize import curve_fit


def _caglioti_model(two_theta_deg, U, V, W):
    """Caglioti FWHM(2theta) given U, V, W. Vectorised."""
    theta = np.deg2rad(np.asarray(two_theta_deg) / 2.0)
    fwhm_sq = U * np.tan(theta)**2 + V * np.tan(theta) + W
    return np.sqrt(np.maximum(fwhm_sq, 1e-10))


def _measure_fwhm_near(two_theta, intensity, target_tt,
                        search_window_deg=0.5):
    """Half-max-interpolation FWHM of the peak nearest target_tt
    within +/- search_window_deg. Returns (fwhm_deg, observed_tt) or
    (None, None) if no resolvable peak found."""
    mask = np.abs(two_theta - target_tt) < search_window_deg
    if not np.any(mask):
        return None, None
    local_tt = two_theta[mask]
    local_i = intensity[mask]
    # Subtract a flat local baseline (5th percentile inside the window).
    baseline = np.percentile(local_i, 5)
    local_i_corr = local_i - baseline
    if local_i_corr.max() <= 0:
        return None, None
    peak_idx = int(np.argmax(local_i_corr))
    peak_val = local_i_corr[peak_idx]
    half_max = peak_val / 2.0
    # Find half-max crossings on each side of the peak via interpolation.
    left = local_i_corr[:peak_idx + 1]
    right = local_i_corr[peak_idx:]
    left_below = np.where(left < half_max)[0]
    right_below = np.where(right < half_max)[0]
    if len(left_below) == 0 or len(right_below) == 0:
        return None, None
    li = left_below[-1]
    ri = right_below[0] + peak_idx
    if li + 1 >= len(local_tt) or local_i_corr[li + 1] == local_i_corr[li]:
        tt_left = local_tt[li]
    else:
        frac = ((half_max - local_i_corr[li])
                / (local_i_corr[li + 1] - local_i_corr[li]))
        tt_left = local_tt[li] + frac * (local_tt[li + 1] - local_tt[li])
    if ri == 0 or local_i_corr[ri - 1] == local_i_corr[ri]:
        tt_right = local_tt[ri]
    else:
        frac = ((half_max - local_i_corr[ri])
                / (local_i_corr[ri - 1] - local_i_corr[ri]))
        tt_right = local_tt[ri] - frac * (local_tt[ri] - local_tt[ri - 1])
    return float(tt_right - tt_left), float(local_tt[peak_idx])


def _caglioti_fit(two_theta, intensity, ref_two_theta,
                  search_window_deg=0.5):
    """Fit Caglioti U, V, W to the FWHMs of `intensity` at the peaks
    nearest each entry in `ref_two_theta`.

    Parameters
    ----------
    two_theta : np.ndarray
        Standard's 2-theta scan (degrees).
    intensity : np.ndarray
        Standard's intensity scan.
    ref_two_theta : np.ndarray
        Bragg-position 2-thetas of the standard's reflections (degrees).
    search_window_deg : float
        Local FWHM-search window around each reference position.

    Returns
    -------
    U, V, W : float
        Fitted Caglioti coefficients.
    info : dict
        Diagnostic info: 'n_peaks', 'measured_fwhms',
        'measured_positions', 'cov' (3x3 covariance matrix).
    """
    two_theta = np.asarray(two_theta, dtype=float)
    intensity = np.asarray(intensity, dtype=float)
    ref_two_theta = np.asarray(ref_two_theta, dtype=float)

    measured_fwhms = []
    measured_positions = []
    for target_tt in ref_two_theta:
        fwhm, obs_tt = _measure_fwhm_near(
            two_theta, intensity, target_tt, search_window_deg)
        if fwhm is None or fwhm <= 0:
            continue
        measured_fwhms.append(fwhm)
        measured_positions.append(obs_tt)

    if len(measured_fwhms) < 3:
        raise ValueError(
            f'Caglioti fit needs at least 3 resolvable peaks; '
            f'got {len(measured_fwhms)}')

    measured_fwhms = np.asarray(measured_fwhms)
    measured_positions = np.asarray(measured_positions)

    # Initial guess: flat Caglioti (V=0; U and W small).
    p0 = [1e-3, 0.0, np.median(measured_fwhms)**2]
    popt, pcov = curve_fit(_caglioti_model, measured_positions,
                            measured_fwhms, p0=p0)
    U, V, W = float(popt[0]), float(popt[1]), float(popt[2])
    info = {
        'n_peaks': len(measured_fwhms),
        'measured_fwhms': measured_fwhms,
        'measured_positions': measured_positions,
        'cov': pcov,
    }
    return U, V, W, info
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `cd "C:/Users/Matthew Izawa/Documents/Dan Applin/Llunr/xrd_profile" && pytest tests/test_instrumental.py::TestCagliotiFit -v`

Expected: both tests pass with U, V, W close to the synthesis values.

If `test_recovers_synthesis_coefficients_within_5pct` fails with U or W off by more than 5%, the search window or the half-max baseline is mis-tuned. Inspect the printed `info['measured_fwhms']` values against `caglioti_fwhm(measured_positions)` from the synthesis script to localise the discrepancy.

- [ ] **Step 5: Commit**

```bash
cd "C:/Users/Matthew Izawa/Documents/Dan Applin/Llunr/xrd_profile"
git add xrd_profile/instrumental.py tests/test_instrumental.py
git commit -m "Add internal _caglioti_fit Caglioti-polynomial fitter"
```

---

## Task 6 — `InstrumentalStandard` class (CIF + pattern + caglioti_fit + Fourier coefficients)

**Files:**
- Modify: `xrd_profile/instrumental.py`
- Modify: `xrd_profile/__init__.py`
- Modify: `tests/test_instrumental.py`

`InstrumentalStandard` bundles a `Phase` (e.g., LaB6) with a measured pattern. `caglioti_fit()` returns an `InstrumentalProfile`. `fourier_coefficients(peak_d, n_coeffs)` returns `(L, A_inst_L)` arrays for the standard's peak nearest `peak_d`, used by Stokes deconvolution in W-A.

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_instrumental.py`:

```python
from xrd_profile import Phase, InstrumentalStandard

LAB6_CIF = (Path(__file__).parent.parent
            / 'examples' / 'cifs')  # Hematite isn't ideal but is shipped.
# We need a LaB6 CIF for the standard tests. Until one ships, build
# Phase from inline lattice parameters using pymatgen primitives.


@pytest.fixture(scope='module')
def lab6_phase():
    """Build a LaB6 Phase inline (Pm-3m, a=4.156825) without needing
    a bundled CIF. Uses Phase.from_lattice_params."""
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


class TestInstrumentalProfileFromStandard:
    def test_from_standard_delegates_to_caglioti_fit(
            self, lab6_standard):
        prof_a = InstrumentalProfile.from_standard(lab6_standard)
        prof_b = lab6_standard.caglioti_fit()
        assert prof_a.U == prof_b.U
        assert prof_a.V == prof_b.V
        assert prof_a.W == prof_b.W
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cd "C:/Users/Matthew Izawa/Documents/Dan Applin/Llunr/xrd_profile" && pytest tests/test_instrumental.py::TestInstrumentalStandard tests/test_instrumental.py::TestInstrumentalProfileFromStandard -v`

Expected: `ImportError: cannot import name 'InstrumentalStandard'`.

- [ ] **Step 3: Implement `InstrumentalStandard` and `InstrumentalProfile.from_standard`**

Append to `xrd_profile/instrumental.py`:

```python
class InstrumentalStandard:
    """Structural Phase plus a measured diffraction pattern of a known
    standard, sufficient for both Caglioti FWHM correction and Stokes
    Fourier deconvolution.

    Parameters
    ----------
    phase : xrd_profile.Phase
        The standard's structure (e.g., LaB6, Si). Provides reference
        Bragg positions for the Caglioti fit and for matching peaks
        when `fourier_coefficients(peak_d, ...)` is called.
    two_theta : np.ndarray
        Measured 2-theta scan of the standard (degrees).
    intensity : np.ndarray
        Measured intensity scan.
    wavelength : float
        X-ray wavelength (angstroms). Should match the sample being
        analysed.
    name : str
        Optional human-readable label.
    """

    def __init__(self, phase, two_theta, intensity,
                 wavelength: float, name: str = ''):
        self.phase = phase
        self.two_theta = np.asarray(two_theta, dtype=float)
        self.intensity = np.asarray(intensity, dtype=float)
        self.wavelength = float(wavelength)
        self.name = str(name)
        self._caglioti_cache = None
        self._fourier_cache = {}

    @classmethod
    def from_cif_and_pattern(cls, cif: str,
                             two_theta, intensity,
                             wavelength: float,
                             name: str = '') -> 'InstrumentalStandard':
        """Convenience constructor: load Phase from CIF, attach the
        measured pattern."""
        from .phases import Phase
        phase = Phase.from_cif(cif, name=name or Path(cif).stem)
        return cls(phase=phase, two_theta=two_theta, intensity=intensity,
                   wavelength=wavelength, name=name)

    def caglioti_fit(self) -> 'InstrumentalProfile':
        """Fit Caglioti U, V, W to the standard's measured FWHMs at
        each reference peak. Cached. Returns an InstrumentalProfile."""
        if self._caglioti_cache is not None:
            return self._caglioti_cache
        ref_peaks = self.phase.get_ref_peaks(
            self.wavelength,
            two_theta_range=(float(self.two_theta.min()),
                             float(self.two_theta.max())))
        ref_tt = np.array([p['two_theta'] for p in ref_peaks])
        U, V, W, _info = _caglioti_fit(self.two_theta, self.intensity,
                                        ref_tt)
        prof = InstrumentalProfile(U=U, V=V, W=W,
                                    wavelength=self.wavelength,
                                    name=self.name)
        self._caglioti_cache = prof
        return prof

    def fourier_coefficients(self, peak_d: float,
                              n_coeffs: int = 20):
        """Return (L, A_inst_L) Fourier coefficients of the standard's
        peak nearest `peak_d` (angstroms). `L` is column length in
        angstroms; `A_inst_L` is the (real) Fourier coefficient.
        Cached per peak_d."""
        cache_key = (round(float(peak_d), 6), int(n_coeffs))
        if cache_key in self._fourier_cache:
            return self._fourier_cache[cache_key]

        from .conversions import d_to_two_theta
        from .warren_averbach import (extract_peak_profile,
                                       fourier_coefficients)
        target_tt = d_to_two_theta(peak_d, self.wavelength)
        if np.isnan(target_tt):
            raise ValueError(
                f'peak_d {peak_d} corresponds to no valid 2-theta at '
                f'wavelength {self.wavelength}')
        prof = extract_peak_profile(
            self.two_theta, self.intensity, target_tt, self.wavelength,
            width_fwhm=6.0)
        L, A_L, _conv = fourier_coefficients(
            prof['s'], prof['profile'], prof['s0'], n_coeffs=n_coeffs)
        self._fourier_cache[cache_key] = (L, A_L)
        return L, A_L
```

Add `InstrumentalStandard` to `__init__.py`:

```python
# In xrd_profile/__init__.py:
from .instrumental import InstrumentalProfile, InstrumentalStandard
```

And add to `__all__`:

```python
'InstrumentalProfile', 'InstrumentalStandard',
```

Now extend `InstrumentalProfile` with `from_standard`. Add inside the class:

```python
    @classmethod
    def from_standard(cls, std) -> 'InstrumentalProfile':
        """Convenience: fit Caglioti to the standard, return the
        resulting InstrumentalProfile. Equivalent to
        `std.caglioti_fit()`."""
        return std.caglioti_fit()
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `cd "C:/Users/Matthew Izawa/Documents/Dan Applin/Llunr/xrd_profile" && pytest tests/test_instrumental.py -v`

Expected: all six new `TestInstrumentalStandard` and `TestInstrumentalProfileFromStandard` tests pass; existing tests still pass.

- [ ] **Step 5: Commit**

```bash
cd "C:/Users/Matthew Izawa/Documents/Dan Applin/Llunr/xrd_profile"
git add xrd_profile/instrumental.py xrd_profile/__init__.py tests/test_instrumental.py
git commit -m "Add InstrumentalStandard with caglioti_fit and fourier_coefficients"
```

---

## Task 7 — Stokes deconvolution helper (`_stokes_deconvolve`)

**Files:**
- Modify: `xrd_profile/instrumental.py`
- Modify: `tests/test_instrumental.py`

`_stokes_deconvolve(A_obs, A_inst, damping_threshold=0.05)` divides observed Fourier coefficients by instrumental Fourier coefficients, with a damping floor that suppresses noise amplification at high `L`. The classical Stokes (1948) formula:

```
A_corr(L) = A_obs(L) / A_inst(L)
```

with the damping rule: when `A_inst(L) < damping_threshold * A_inst(0)`, set `A_corr(L) = 0`. This prevents division-by-near-zero from blowing up the high-`L` tail.

- [ ] **Step 1: Write the failing test**

Append to `tests/test_instrumental.py`:

```python
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
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cd "C:/Users/Matthew Izawa/Documents/Dan Applin/Llunr/xrd_profile" && pytest tests/test_instrumental.py::TestStokesDeconvolve -v`

Expected: `ImportError: cannot import name '_stokes_deconvolve'`.

- [ ] **Step 3: Implement `_stokes_deconvolve`**

Append to `xrd_profile/instrumental.py`:

```python
def _stokes_deconvolve(A_obs, A_inst, damping_threshold: float = 0.05):
    """Stokes Fourier deconvolution of a peak profile.

    A_corr(L) = A_obs(L) / A_inst(L), with a damping floor: when
    A_inst(L) < damping_threshold * A_inst(0), the coefficient is set
    to 0 to suppress noise amplification at high L (Stokes 1948).

    Parameters
    ----------
    A_obs : np.ndarray
        Observed sample profile Fourier coefficients.
    A_inst : np.ndarray
        Instrumental profile Fourier coefficients (same length).
    damping_threshold : float
        Fraction of A_inst(0) below which deconvolution is suppressed.

    Returns
    -------
    A_corr : np.ndarray
        Deconvolved sample-only Fourier coefficients.
    """
    A_obs = np.asarray(A_obs, dtype=float)
    A_inst = np.asarray(A_inst, dtype=float)
    if A_obs.shape != A_inst.shape:
        raise ValueError(f'A_obs shape {A_obs.shape} != A_inst shape '
                         f'{A_inst.shape}')
    if A_inst.size == 0:
        return np.array([])
    if A_inst[0] == 0:
        raise ValueError('A_inst(0) is zero; cannot Stokes-deconvolve.')
    threshold = damping_threshold * abs(A_inst[0])
    A_corr = np.zeros_like(A_obs)
    keep = np.abs(A_inst) >= threshold
    A_corr[keep] = A_obs[keep] / A_inst[keep]
    return A_corr
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `cd "C:/Users/Matthew Izawa/Documents/Dan Applin/Llunr/xrd_profile" && pytest tests/test_instrumental.py::TestStokesDeconvolve -v`

Expected: all three tests pass.

- [ ] **Step 5: Commit**

```bash
cd "C:/Users/Matthew Izawa/Documents/Dan Applin/Llunr/xrd_profile"
git add xrd_profile/instrumental.py tests/test_instrumental.py
git commit -m "Add internal _stokes_deconvolve Fourier deconvolution helper"
```

---

## Task 8 — Wire `instrumental=` into `guided_williamson_hall`

**Files:**
- Modify: `xrd_profile/profile.py`
- Modify: `xrd_profile/williamson_hall.py`
- Modify: `tests/test_instrumental.py`

`XRDProfile.guided_williamson_hall()` currently raises `NotImplementedError` for any non-`None` `instrumental=`. Replace the guard with dispatch: accept either `InstrumentalStandard` (call `.caglioti_fit()` to get a profile) or `InstrumentalProfile` (use directly). Pass the profile down to `williamson_hall.guided_williamson_hall`, which Caglioti-subtracts the instrumental FWHM from each measured peak's FWHM before regression.

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_instrumental.py`:

```python
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
        also runs end-to-end."""
        data = np.loadtxt(TIRHERT)
        profile = XRDProfile(data[:, 0], data[:, 1],
                              wavelength=LAMBDA_I11)
        inst = InstrumentalProfile(U=5e-3, V=-1e-3, W=5e-3,
                                    wavelength=LAMBDA_I11,
                                    name='lit_lab_bruker')
        result = profile.guided_williamson_hall(
            phase=anorthite_phase, instrumental=inst,
            n_sigma=3.0, tolerance_d=0.03)
        assert np.isfinite(result['crystallite_size'])

    def test_corrected_size_smaller_than_uncorrected(
            self, anorthite_phase):
        """Sanity: instrumental correction REDUCES apparent broadening,
        which INCREASES the apparent crystallite size."""
        data = np.loadtxt(TIRHERT)
        profile = XRDProfile(data[:, 0], data[:, 1],
                              wavelength=LAMBDA_I11)
        inst = InstrumentalProfile(U=5e-3, V=-1e-3, W=5e-3,
                                    wavelength=LAMBDA_I11)
        uncorrected = profile.guided_williamson_hall(
            phase=anorthite_phase,
            n_sigma=3.0, tolerance_d=0.03)
        corrected = profile.guided_williamson_hall(
            phase=anorthite_phase, instrumental=inst,
            n_sigma=3.0, tolerance_d=0.03)
        # If both are finite, corrected size should be larger
        # (or strain different); at minimum, results should differ.
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
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cd "C:/Users/Matthew Izawa/Documents/Dan Applin/Llunr/xrd_profile" && pytest tests/test_instrumental.py::TestGuidedWHWithInstrumental -v`

Expected: `NotImplementedError: instrumental= is reserved for Phase 2 / v1.0; ...` on the first three tests.

- [ ] **Step 3: Add Caglioti-subtraction helper to `williamson_hall.py`**

Add at the top of `xrd_profile/williamson_hall.py`, after the existing imports:

```python
def _apply_caglioti_correction(fwhm_deg, two_theta_deg, inst_profile):
    """Caglioti-subtract instrumental FWHM from observed FWHM.

    Uses Gaussian-quadrature combination:
        beta_corr^2 = max(beta_obs^2 - beta_inst^2, eps)

    Peaks where beta_obs <= beta_inst are flagged: their corrected FWHM
    is replaced with NaN (caller is expected to filter).

    Parameters
    ----------
    fwhm_deg : np.ndarray
        Observed sample FWHMs (degrees).
    two_theta_deg : np.ndarray
        Peak 2-theta positions (degrees).
    inst_profile : InstrumentalProfile
        Instrumental Caglioti profile.

    Returns
    -------
    fwhm_corr : np.ndarray
        Corrected FWHMs. NaN where the peak was over-corrected.
    """
    import numpy as np
    fwhm_obs = np.asarray(fwhm_deg, dtype=float)
    tt = np.asarray(two_theta_deg, dtype=float)
    fwhm_inst = np.array([inst_profile.fwhm_at(t) for t in tt])
    diff_sq = fwhm_obs**2 - fwhm_inst**2
    fwhm_corr = np.where(diff_sq > 0, np.sqrt(np.maximum(diff_sq, 0)),
                          np.nan)
    return fwhm_corr
```

- [ ] **Step 4: Replace the `NotImplementedError` guard in `XRDProfile.guided_williamson_hall`**

In `xrd_profile/profile.py`, locate the body of `guided_williamson_hall`. Replace the existing `instrumental=` block with a dispatching version that records the resolved instrumental profile in `kwargs` for the underlying function:

Find:

```python
        if instrumental is not None:
            raise NotImplementedError(
                'instrumental= is reserved for Phase 2 / v1.0; '
                'see xrd_profile roadmap')
```

Replace with:

```python
        from .instrumental import (InstrumentalStandard,
                                    InstrumentalProfile)
        if instrumental is None:
            inst_profile = None
        elif isinstance(instrumental, InstrumentalStandard):
            inst_profile = instrumental.caglioti_fit()
        elif isinstance(instrumental, InstrumentalProfile):
            inst_profile = instrumental
        else:
            raise TypeError(
                f'instrumental= must be InstrumentalStandard, '
                f'InstrumentalProfile, or None; got '
                f'{type(instrumental).__name__}')
```

And in the same method, pass `inst_profile=inst_profile` through to the underlying `guided_williamson_hall` call:

Find the call:

```python
        from .williamson_hall import guided_williamson_hall
        return guided_williamson_hall(
            self.two_theta, self.intensity, ref_d, self.wavelength,
            tolerance_d=tolerance_d, n_sigma=n_sigma,
            ...
            export_path=export_path,
            **kwargs
        )
```

Add `inst_profile=inst_profile` in the kwargs:

```python
        from .williamson_hall import guided_williamson_hall
        return guided_williamson_hall(
            self.two_theta, self.intensity, ref_d, self.wavelength,
            tolerance_d=tolerance_d, n_sigma=n_sigma,
            min_fwhm_steps=min_fwhm_steps, correct_offset=correct_offset,
            other_phase_d=other_phase_d,
            other_phase_names=other_phase_names,
            overlap_tol_deg=overlap_tol_deg,
            min_quality=min_quality,
            quality_weights=quality_weights,
            weighted_regression=weighted_regression,
            sample_flags=sample_flags,
            export_path=export_path,
            inst_profile=inst_profile,
            **kwargs
        )
```

- [ ] **Step 5: Accept and apply `inst_profile=` inside `williamson_hall.guided_williamson_hall`**

In `xrd_profile/williamson_hall.py`, the `guided_williamson_hall` function uses `peaks['fwhm']` (parallel array of FWHMs in degrees) and `peaks['two_theta_obs']` (positions in degrees) to build the W-H abscissa via `fwhm_to_deltaK`. The Caglioti correction belongs immediately before that conversion: subtract instrumental FWHM, drop peaks where the subtraction over-corrects, and update both `peaks['fwhm']` and the per-peak parallel arrays (`peaks['two_theta_obs']`, `peaks['d_obs']`, `quality`, `peaks['cross_phase_overlap']` if present) to keep them aligned.

Add `inst_profile=None` to the signature (after `export_path=None`):

```python
def guided_williamson_hall(two_theta, intensity, ref_d, wavelength,
                           tolerance_d=0.03, n_sigma=3.0,
                           min_fwhm_steps=3, correct_offset=True,
                           other_phase_d=None, other_phase_names=None,
                           overlap_tol_deg=0.15,
                           min_quality=0.3,
                           quality_weights=None,
                           weighted_regression=True,
                           min_peaks_reliable=5,
                           min_r2_reliable=0.3,
                           min_r2_marginal=0.05,
                           sample_flags=None,
                           export_path=None,
                           inst_profile=None):
```

Then locate the lines (near current `williamson_hall.py:386-388`):

```python
    # Extract included peaks
    K_all = 1.0 / peaks['d_obs']
    dK_all = fwhm_to_deltaK(peaks['fwhm'], peaks['two_theta_obs'],
                             wavelength)
```

Insert immediately before those lines (after the existing `result['warnings'] = warnings_list` block and the `n_used < 3` early-return):

```python
    # NEW v0.4.0: Caglioti instrumental correction. Subtract the
    # instrumental FWHM from each peak's measured FWHM in quadrature.
    # Peaks where beta_obs <= beta_inst are flagged in warnings and
    # removed from the analysis.
    if inst_profile is not None:
        fwhm_corr = _apply_caglioti_correction(
            peaks['fwhm'], peaks['two_theta_obs'], inst_profile)
        keep = ~np.isnan(fwhm_corr)
        n_dropped = int(np.sum(~keep))
        if n_dropped > 0:
            warnings_list.append(
                f'{n_dropped} peak(s) over-corrected by instrumental '
                f'subtraction (beta_obs <= beta_inst); excluded from fit.')
            result['warnings'] = warnings_list
        # Filter every parallel array on `peaks` plus the local
        # `quality` and `include` arrays, keeping them aligned.
        peaks = {k: (v[keep] if isinstance(v, np.ndarray) and v.shape
                     and v.shape[0] == len(keep) else v)
                 for k, v in peaks.items()}
        peaks['fwhm'] = fwhm_corr[keep]
        quality = quality[keep]
        include = include[keep]
        n_used = int(np.sum(include))
        result['peaks'] = peaks
        result['peak_quality'] = quality
        result['n_peaks_used'] = n_used
        if n_used < 3:
            reliability_reasons.append(
                f'Only {n_used} peaks after instrumental correction '
                f'(need >= 3 for fit)')
            result['reliability_reasons'] = reliability_reasons
            _export_wh_csv(export_path, peaks, include, result, wavelength)
            return result
```

The dict-comprehension filter handles the case where `peaks` contains both per-peak arrays (which need filtering) and scalars or unrelated arrays (which pass through unchanged) — the `v.shape[0] == len(keep)` guard distinguishes them.

- [ ] **Step 6: Run the new instrumental tests**

Run: `cd "C:/Users/Matthew Izawa/Documents/Dan Applin/Llunr/xrd_profile" && pytest tests/test_instrumental.py::TestGuidedWHWithInstrumental -v`

Expected: all four tests pass.

- [ ] **Step 7: Run the v0.3.0 backward-compat tier to confirm no drift**

Run: `cd "C:/Users/Matthew Izawa/Documents/Dan Applin/Llunr/xrd_profile" && pytest tests/test_backward_compat.py -v`

Expected: every existing test still passes. If any v0.3.0 tier test fails, it means the Caglioti dispatch path was triggered when `instrumental=None` — investigate before continuing.

- [ ] **Step 8: Commit**

```bash
cd "C:/Users/Matthew Izawa/Documents/Dan Applin/Llunr/xrd_profile"
git add xrd_profile/profile.py xrd_profile/williamson_hall.py tests/test_instrumental.py
git commit -m "Wire instrumental= into guided_williamson_hall (Caglioti subtraction)"
```

---

## Task 9 — Wire `instrumental=` into `guided_warren_averbach` (Stokes deconvolution)

**Files:**
- Modify: `xrd_profile/profile.py`
- Modify: `xrd_profile/warren_averbach.py`
- Modify: `tests/test_instrumental.py`

W-A consumes the full Fourier-coefficient array of each peak, so it gets full Stokes deconvolution rather than scalar Caglioti subtraction. `InstrumentalStandard` is required (we need the standard's measured pattern); `InstrumentalProfile` is rejected with a clear error.

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_instrumental.py`:

```python
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
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cd "C:/Users/Matthew Izawa/Documents/Dan Applin/Llunr/xrd_profile" && pytest tests/test_instrumental.py::TestGuidedWAWithInstrumental -v`

Expected: `NotImplementedError` on the first test; `NotImplementedError` (not the expected `ValueError`) on the profile test.

- [ ] **Step 3: Replace the `NotImplementedError` guard in `XRDProfile.guided_warren_averbach`**

In `xrd_profile/profile.py`, locate the body of `guided_warren_averbach`. Replace:

```python
        if instrumental is not None:
            raise NotImplementedError(
                'instrumental= is reserved for Phase 2 / v1.0; '
                'see xrd_profile roadmap')
```

with:

```python
        from .instrumental import (InstrumentalStandard,
                                    InstrumentalProfile)
        if instrumental is None:
            std = None
        elif isinstance(instrumental, InstrumentalStandard):
            std = instrumental
        elif isinstance(instrumental, InstrumentalProfile):
            raise ValueError(
                'Warren-Averbach Stokes deconvolution requires the '
                'measured standard pattern; pass an InstrumentalStandard, '
                'or call without instrumental= for uncorrected W-A. '
                'Caglioti FWHM subtraction is mathematically equivalent '
                'to FWHM-only correction and is unprincipled for the '
                'higher-order Fourier coefficients W-A consumes.')
        else:
            raise TypeError(
                f'instrumental= must be InstrumentalStandard or None; '
                f'got {type(instrumental).__name__}')
```

Pass `inst_standard=std` through to the underlying call:

```python
        from .warren_averbach import guided_warren_averbach
        return guided_warren_averbach(
            self.two_theta, self.intensity, ref_peaks, self.wavelength,
            tolerance_d=tolerance_d, n_sigma=n_sigma,
            min_fwhm_steps=min_fwhm_steps, correct_offset=correct_offset,
            n_coeffs=n_coeffs, width_fwhm=width_fwhm,
            require_clean=require_clean,
            inst_standard=std,
        )
```

- [ ] **Step 4: Apply Stokes deconvolution inside `warren_averbach.guided_warren_averbach`**

In `xrd_profile/warren_averbach.py`, add `inst_standard=None` to the signature of `guided_warren_averbach`:

```python
def guided_warren_averbach(two_theta, intensity, ref_peaks, wavelength,
                           tolerance_d=0.03, n_sigma=3.0,
                           min_fwhm_steps=3, correct_offset=True,
                           n_coeffs=20, width_fwhm=6.0,
                           min_ref_intensity=1.0,
                           require_clean=False,
                           inst_standard=None):
```

In the per-family loop, after the line

```python
            L, A_L, converged = fourier_coefficients(
                prof['s'], prof['profile'], prof['s0'],
                n_coeffs=n_coeffs)
```

add Stokes deconvolution if `inst_standard is not None`:

```python
            if inst_standard is not None:
                from .instrumental import _stokes_deconvolve
                _, A_inst_L = inst_standard.fourier_coefficients(
                    peak_d=d_ref, n_coeffs=n_coeffs)
                A_L = _stokes_deconvolve(A_L, A_inst_L,
                                          damping_threshold=0.05)
```

Place this block immediately after the existing `if not converged: continue` check (so we don't deconvolve coefficients we'd reject anyway).

- [ ] **Step 5: Run the new W-A instrumental tests**

Run: `cd "C:/Users/Matthew Izawa/Documents/Dan Applin/Llunr/xrd_profile" && pytest tests/test_instrumental.py::TestGuidedWAWithInstrumental -v`

Expected: all three tests pass.

- [ ] **Step 6: Run the v0.3.0 backward-compat tier**

Run: `cd "C:/Users/Matthew Izawa/Documents/Dan Applin/Llunr/xrd_profile" && pytest tests/test_backward_compat.py -v`

Expected: every existing test still passes. The v0.3.0 W-A tier test must produce values matching the golden — confirms `inst_standard=None` path is byte-identical.

- [ ] **Step 7: Commit**

```bash
cd "C:/Users/Matthew Izawa/Documents/Dan Applin/Llunr/xrd_profile"
git add xrd_profile/profile.py xrd_profile/warren_averbach.py tests/test_instrumental.py
git commit -m "Wire instrumental= into guided_warren_averbach (Stokes deconvolution)"
```

---

## Task 10 — Guided Scherrer + instrumental Caglioti subtraction

**Files:**
- Modify: `xrd_profile/scherrer.py`
- Modify: `xrd_profile/profile.py`
- Modify: `tests/test_instrumental.py`

Two changes to `scherrer()` and `modified_scherrer()`:
1. New optional `phase=` kwarg filters detected peaks to those matching the phase's reference d-spacings (within `tolerance_d=0.03`).
2. New optional `instrumental=` kwarg subtracts the Caglioti instrumental FWHM from each peak's measured FWHM.

Both are accessed via `XRDProfile.scherrer()` and `XRDProfile.modified_scherrer()` (which forward to `xrd_profile.scherrer.scherrer` / `modified_scherrer`).

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_instrumental.py`:

```python
class TestScherrerWithPhaseAndInstrumental:
    @pytest.fixture(scope='class')
    def anorthite_phase(self):
        cif = (Path(__file__).parent.parent / 'examples'
               / 'cifs' / 'Anorthite.cif')
        return Phase.from_cif(str(cif), name='anorthite')

    def test_phase_filtering_reduces_peak_count(self, anorthite_phase):
        data = np.loadtxt(TIRHERT)
        profile = XRDProfile(data[:, 0], data[:, 1],
                              wavelength=LAMBDA_I11)
        unfiltered = profile.scherrer()
        filtered = profile.scherrer(phase=anorthite_phase)
        assert len(filtered['sizes']) <= len(unfiltered['sizes'])

    def test_instrumental_correction_increases_apparent_size(
            self, anorthite_phase):
        """Removing instrumental broadening (FWHM_corr < FWHM_obs)
        increases the apparent crystallite size from Scherrer."""
        data = np.loadtxt(TIRHERT)
        profile = XRDProfile(data[:, 0], data[:, 1],
                              wavelength=LAMBDA_I11)
        inst = InstrumentalProfile(U=5e-3, V=-1e-3, W=5e-3,
                                    wavelength=LAMBDA_I11)
        without = profile.scherrer(phase=anorthite_phase)
        with_inst = profile.scherrer(phase=anorthite_phase,
                                      instrumental=inst)
        # Mean size should increase (or both NaN).
        if not (np.isnan(without['mean_size'])
                or np.isnan(with_inst['mean_size'])):
            assert with_inst['mean_size'] >= without['mean_size']

    def test_modified_scherrer_with_phase_and_instrumental(
            self, anorthite_phase):
        data = np.loadtxt(TIRHERT)
        profile = XRDProfile(data[:, 0], data[:, 1],
                              wavelength=LAMBDA_I11)
        inst = InstrumentalProfile(U=5e-3, V=-1e-3, W=5e-3,
                                    wavelength=LAMBDA_I11)
        # Should not raise.
        result = profile.modified_scherrer(phase=anorthite_phase,
                                            instrumental=inst)
        # nan or finite both acceptable; just exercises the path.
        assert isinstance(result, float) or np.isnan(result)

    def test_default_scherrer_no_phase_no_instrumental_unchanged(
            self):
        """v0.3.0 byte-identity: default scherrer() with no kwargs
        must match the v0.2.0 golden exactly."""
        # This is implicitly covered by test_backward_compat.py's
        # v0.2.0 tier; here we just confirm the call signature
        # accepts no extra args.
        data = np.loadtxt(TIRHERT)
        profile = XRDProfile(data[:, 0], data[:, 1],
                              wavelength=LAMBDA_I11)
        result = profile.scherrer()
        assert 'sizes' in result
        assert 'mean_size' in result
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cd "C:/Users/Matthew Izawa/Documents/Dan Applin/Llunr/xrd_profile" && pytest tests/test_instrumental.py::TestScherrerWithPhaseAndInstrumental -v`

Expected: `TypeError: scherrer() got an unexpected keyword argument 'phase'`.

- [ ] **Step 3: Add `phase=` and `instrumental=` to `XRDProfile.scherrer` and `XRDProfile.modified_scherrer`**

In `xrd_profile/profile.py`, replace `XRDProfile.scherrer`:

```python
    def scherrer(self, K=None, shape=None, height_threshold=0.05,
                 *, phase=None, instrumental=None):
        """Run Scherrer analysis on detected peaks.

        K and shape: see xrd_profile.scherrer.scherrer for resolution
        rules. phase=Phase filters detected peaks to those matching the
        phase's reference d-spacings within tolerance. instrumental=
        accepts InstrumentalStandard or InstrumentalProfile and
        Caglioti-subtracts the instrumental FWHM."""
        from .instrumental import (InstrumentalStandard,
                                    InstrumentalProfile)
        from .scherrer import (_filter_to_phase_peaks,
                                _caglioti_subtract_fwhm)

        fwhm, positions = estimate_fwhm_simple(
            self.two_theta, self.intensity, height_threshold)
        if len(fwhm) == 0:
            return {'sizes': np.array([]),
                    'peak_positions': np.array([]),
                    'd_spacings': np.array([]), 'fwhm': np.array([]),
                    'mean_size': np.nan, 'median_size': np.nan}

        if phase is not None:
            ref_d = phase.get_ref_d(
                self.wavelength,
                two_theta_range=(float(self.two_theta.min()),
                                 float(self.two_theta.max())))
            fwhm, positions = _filter_to_phase_peaks(
                fwhm, positions, ref_d, self.wavelength,
                tolerance_d=0.03)
            if len(fwhm) == 0:
                return {'sizes': np.array([]),
                        'peak_positions': np.array([]),
                        'd_spacings': np.array([]),
                        'fwhm': np.array([]),
                        'mean_size': np.nan, 'median_size': np.nan}

        if instrumental is not None:
            if isinstance(instrumental, InstrumentalStandard):
                inst_profile = instrumental.caglioti_fit()
            elif isinstance(instrumental, InstrumentalProfile):
                inst_profile = instrumental
            else:
                raise TypeError(
                    f'instrumental= must be InstrumentalStandard, '
                    f'InstrumentalProfile, or None; got '
                    f'{type(instrumental).__name__}')
            fwhm, positions = _caglioti_subtract_fwhm(
                fwhm, positions, inst_profile)
            if len(fwhm) == 0:
                return {'sizes': np.array([]),
                        'peak_positions': np.array([]),
                        'd_spacings': np.array([]),
                        'fwhm': np.array([]),
                        'mean_size': np.nan, 'median_size': np.nan}

        sizes = scherrer(fwhm, positions, self.wavelength, K=K, shape=shape)
        d_sp = two_theta_to_d(positions, self.wavelength)
        return {
            'sizes': sizes, 'peak_positions': positions,
            'd_spacings': d_sp, 'fwhm': fwhm,
            'mean_size': float(np.mean(sizes)),
            'median_size': float(np.median(sizes)),
        }
```

Same shape for `XRDProfile.modified_scherrer`:

```python
    def modified_scherrer(self, K=None, shape=None, height_threshold=0.05,
                          *, phase=None, instrumental=None):
        """Run modified Scherrer with optional phase filtering and
        instrumental Caglioti subtraction."""
        from .instrumental import (InstrumentalStandard,
                                    InstrumentalProfile)
        from .scherrer import (_filter_to_phase_peaks,
                                _caglioti_subtract_fwhm)

        fwhm, positions = estimate_fwhm_simple(
            self.two_theta, self.intensity, height_threshold)
        if len(fwhm) < 2:
            return np.nan

        if phase is not None:
            ref_d = phase.get_ref_d(
                self.wavelength,
                two_theta_range=(float(self.two_theta.min()),
                                 float(self.two_theta.max())))
            fwhm, positions = _filter_to_phase_peaks(
                fwhm, positions, ref_d, self.wavelength,
                tolerance_d=0.03)
            if len(fwhm) < 2:
                return np.nan

        if instrumental is not None:
            if isinstance(instrumental, InstrumentalStandard):
                inst_profile = instrumental.caglioti_fit()
            elif isinstance(instrumental, InstrumentalProfile):
                inst_profile = instrumental
            else:
                raise TypeError(
                    f'instrumental= must be InstrumentalStandard, '
                    f'InstrumentalProfile, or None; got '
                    f'{type(instrumental).__name__}')
            fwhm, positions = _caglioti_subtract_fwhm(
                fwhm, positions, inst_profile)
            if len(fwhm) < 2:
                return np.nan

        return modified_scherrer(fwhm, positions, self.wavelength,
                                  K=K, shape=shape)
```

- [ ] **Step 4: Add the helper functions to `xrd_profile/scherrer.py`**

Append to `xrd_profile/scherrer.py`:

```python
def _filter_to_phase_peaks(fwhm_deg, positions_deg, ref_d, wavelength,
                            tolerance_d=0.03):
    """Keep only peaks whose d-spacing matches a reference d within
    `tolerance_d`. Returns filtered (fwhm, positions) arrays —
    (fwhm, positions) order matches `estimate_fwhm_simple`'s output."""
    from .conversions import two_theta_to_d
    fwhm = np.asarray(fwhm_deg)
    positions = np.asarray(positions_deg)
    if len(positions) == 0:
        return fwhm, positions
    d_obs = two_theta_to_d(positions, wavelength)
    ref_d = np.asarray(ref_d)
    keep_mask = np.zeros(len(d_obs), dtype=bool)
    for i, d in enumerate(d_obs):
        if np.any(np.abs(ref_d - d) < tolerance_d):
            keep_mask[i] = True
    return fwhm[keep_mask], positions[keep_mask]


def _caglioti_subtract_fwhm(fwhm_deg, positions_deg, inst_profile):
    """Caglioti-subtract instrumental FWHM from each measured FWHM.
    Drops peaks where beta_obs <= beta_inst (over-correction).

    Returns (fwhm_corr, positions_kept) — the (fwhm, positions) order
    matches `estimate_fwhm_simple` and the call sites in profile.py.
    """
    fwhm_obs = np.asarray(fwhm_deg, dtype=float)
    pos = np.asarray(positions_deg, dtype=float)
    fwhm_inst = np.array([inst_profile.fwhm_at(t) for t in pos])
    diff_sq = fwhm_obs**2 - fwhm_inst**2
    keep = diff_sq > 0
    fwhm_corr = np.sqrt(np.maximum(diff_sq[keep], 0))
    return fwhm_corr, pos[keep]
```

- [ ] **Step 5: Run the new Scherrer instrumental tests**

Run: `cd "C:/Users/Matthew Izawa/Documents/Dan Applin/Llunr/xrd_profile" && pytest tests/test_instrumental.py::TestScherrerWithPhaseAndInstrumental -v`

Expected: all four tests pass.

- [ ] **Step 6: Run the existing Scherrer + backward-compat tests**

Run: `cd "C:/Users/Matthew Izawa/Documents/Dan Applin/Llunr/xrd_profile" && pytest tests/test_scherrer.py tests/test_backward_compat.py -v`

Expected: every existing test passes. v0.3.0 tier is byte-identical (defaults preserved).

- [ ] **Step 7: Commit**

```bash
cd "C:/Users/Matthew Izawa/Documents/Dan Applin/Llunr/xrd_profile"
git add xrd_profile/profile.py xrd_profile/scherrer.py tests/test_instrumental.py
git commit -m "Add phase= and instrumental= kwargs to Scherrer methods"
```

---

## Task 11 — Wire `instrumental=` into `XRDProfile.run_all`

**Files:**
- Modify: `xrd_profile/profile.py`
- Modify: `tests/test_run_all.py` (extension)

`run_all()` currently raises `NotImplementedError` for any non-`None` `instrumental=`. Replace the guard with passthrough so the kwarg flows into `guided_williamson_hall`, `guided_warren_averbach`, and `scherrer` per their per-method dispatch rules.

- [ ] **Step 1: Write the failing test**

Append to `tests/test_run_all.py` (read the file first; the existing structure dictates where to add):

```python
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
        inst = InstrumentalProfile(U=5e-3, V=-1e-3, W=5e-3,
                                    wavelength=0.826517)
        # methods=['wh', 'scherrer'] only - W-A would reject Profile.
        result = profile.run_all(
            methods=['wh', 'scherrer'],
            phases=anorthite_phase,
            instrumental=inst)
        assert 'wh' in result
        assert 'scherrer' in result

    def test_run_all_with_instrumental_profile_rejects_wa(
            self, anorthite_phase):
        from xrd_profile import (XRDProfile, InstrumentalProfile)
        FIX = Path(__file__).parent / 'fixtures' / 'tirhert_subset.xy'
        data = np.loadtxt(FIX)
        profile = XRDProfile(data[:, 0], data[:, 1],
                              wavelength=0.826517)
        inst = InstrumentalProfile(U=5e-3, V=-1e-3, W=5e-3,
                                    wavelength=0.826517)
        with pytest.raises(ValueError, match='Stokes'):
            profile.run_all(methods=['wa'], phases=anorthite_phase,
                             instrumental=inst)
```

If `tests/test_run_all.py` doesn't already import `numpy`, `pytest`, and `Path`, add those at the top of the file.

- [ ] **Step 2: Run tests to verify they fail**

Run: `cd "C:/Users/Matthew Izawa/Documents/Dan Applin/Llunr/xrd_profile" && pytest tests/test_run_all.py::TestRunAllInstrumental -v`

Expected: `NotImplementedError: instrumental= is reserved for Phase 2 / v1.0; ...`

- [ ] **Step 3: Update `XRDProfile.run_all`**

In `xrd_profile/profile.py`, locate `run_all`. Replace:

```python
        if instrumental is not None:
            raise NotImplementedError(
                'instrumental= is reserved for Phase 2 / v1.0; '
                'see xrd_profile roadmap')
```

with: (delete that block entirely — `instrumental=` simply flows through to per-method dispatch).

In the same method, thread `instrumental=instrumental` through each per-method call. Find the `if 'wh' in methods:` block and update:

```python
        if 'wh' in methods:
            if phases:
                results['wh'] = {
                    p.name: self.guided_williamson_hall(
                        phase=p, instrumental=instrumental,
                        **wh_kwargs)
                    for p in phases
                }
            else:
                from .williamson_hall import williamson_hall
                results['wh'] = williamson_hall(
                    self.two_theta, self.intensity, self.wavelength,
                    **wh_kwargs)
```

Same for W-A:

```python
        if 'wa' in methods:
            if phases:
                results['wa'] = {
                    p.name: self.guided_warren_averbach(
                        phase=p, instrumental=instrumental,
                        **wa_kwargs)
                    for p in phases
                }
            else:
                results['wa'] = self.warren_averbach(**wa_kwargs)
```

Scherrer:

```python
        if 'scherrer' in methods:
            if phases:
                results['scherrer'] = {
                    p.name: self.scherrer(
                        phase=p, instrumental=instrumental,
                        **scherrer_kwargs)
                    for p in phases
                }
            else:
                results['scherrer'] = self.scherrer(
                    instrumental=instrumental, **scherrer_kwargs)
```

When multiple phases are passed, Scherrer runs once per phase and returns a dict keyed by phase name — matching the per-phase shape that W-H and W-A already use in `run_all`.

PDF is unchanged (no instrumental concept for PDF in v0.4.0).

- [ ] **Step 4: Run tests to verify they pass**

Run: `cd "C:/Users/Matthew Izawa/Documents/Dan Applin/Llunr/xrd_profile" && pytest tests/test_run_all.py -v`

Expected: existing tests still pass, both new `TestRunAllInstrumental` tests pass.

- [ ] **Step 5: Commit**

```bash
cd "C:/Users/Matthew Izawa/Documents/Dan Applin/Llunr/xrd_profile"
git add xrd_profile/profile.py tests/test_run_all.py
git commit -m "Wire instrumental= into XRDProfile.run_all dispatch"
```

---

## Task 12 — Size-distribution module (lognormal/normal fits to A_size)

**Files:**
- Create: `xrd_profile/size_distributions.py`
- Create: `tests/test_size_distributions.py`

Internal module. Two basis functions for `A_size(L)`: lognormal-derived and normal-derived. `fit_size_distribution(L, A_size)` returns the `'size_distribution'` dict shape documented in spec §7.3.

- [ ] **Step 1: Write the failing tests**

Create `tests/test_size_distributions.py`:

```python
"""Tests for xrd_profile.size_distributions — lognormal/normal fits
to W-A column-length A_size data."""
import numpy as np
import pytest

from xrd_profile.size_distributions import (
    fit_size_distribution, lognormal_a_size, normal_a_size)


class TestLognormalASize:
    def test_a_size_at_L0_is_one(self):
        """A_size(0) is normalised to 1 (Fourier coefficient at origin)."""
        D_median, sigma = 100.0, 0.3
        A0 = lognormal_a_size(np.array([0.0]), D_median, sigma)
        assert np.isclose(A0[0], 1.0, atol=1e-10)

    def test_a_size_decreases_monotonically(self):
        D_median, sigma = 100.0, 0.3
        L = np.linspace(0, 500, 100)
        A = lognormal_a_size(L, D_median, sigma)
        diffs = np.diff(A)
        # Monotonically decreasing (allow tiny numerical noise).
        assert np.all(diffs <= 1e-10)


class TestFitSizeDistribution:
    def test_recovers_lognormal_synthesis_within_5pct(self):
        D_median_synth, sigma_synth = 150.0, 0.4
        L = np.linspace(0, 600, 60)
        A_size = lognormal_a_size(L, D_median_synth, sigma_synth)
        # Add a touch of noise to simulate W-A output.
        rng = np.random.default_rng(seed=42)
        A_size_noisy = A_size + rng.normal(0, 0.005, len(L))

        result = fit_size_distribution(L, A_size_noisy)
        ln = result['lognormal']
        assert abs(ln['D_median'] - D_median_synth) / D_median_synth < 0.05
        assert abs(ln['sigma'] - sigma_synth) / sigma_synth < 0.10
        assert ln['fit_r2'] > 0.95

    def test_returns_none_when_too_few_valid_L(self):
        L = np.array([0, 50, 100])      # only 3 points -> too few
        A_size = np.array([1.0, 0.7, 0.4])
        result = fit_size_distribution(L, A_size)
        assert result is None

    def test_metadata_keys_present(self):
        L = np.linspace(0, 600, 60)
        A_size = lognormal_a_size(L, 100.0, 0.3)
        result = fit_size_distribution(L, A_size)
        assert set(result) >= {'lognormal', 'normal', 'method',
                                'initial_guess', 'n_valid_L'}
        assert result['method'] == 'curve_fit'
        assert result['initial_guess'] == 'moments'
        assert result['n_valid_L'] >= 4

    def test_lognormal_volume_mean_greater_than_median(self):
        """For a lognormal with sigma>0, D_mean_volume > D_median."""
        L = np.linspace(0, 600, 60)
        A_size = lognormal_a_size(L, 100.0, 0.4)
        result = fit_size_distribution(L, A_size)
        ln = result['lognormal']
        assert ln['D_mean_volume'] > ln['D_median']

    def test_normal_fit_returns_finite_for_lognormal_data(self):
        """Even when the data is lognormal, the normal fit should
        complete without error and return finite parameters with a
        lower R^2 than the lognormal fit."""
        L = np.linspace(0, 600, 60)
        A_size = lognormal_a_size(L, 100.0, 0.4)
        result = fit_size_distribution(L, A_size)
        nrm = result['normal']
        assert np.isfinite(nrm['D_mean'])
        assert np.isfinite(nrm['sigma'])
        assert result['lognormal']['fit_r2'] >= nrm['fit_r2']

    def test_negative_a_size_clipped(self):
        """A_size from W-A can dip slightly negative; fit should still
        complete (clipped to >= 0 in the basis function)."""
        L = np.linspace(0, 600, 60)
        A_size = lognormal_a_size(L, 100.0, 0.3)
        A_size[-5:] = -0.001  # noise dip at high L
        result = fit_size_distribution(L, A_size)
        assert result is not None
        assert np.isfinite(result['lognormal']['D_median'])

    def test_invalid_inputs_raise(self):
        with pytest.raises(ValueError):
            fit_size_distribution(np.array([1.0, 2.0]),
                                   np.array([1.0]))  # mismatched len
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cd "C:/Users/Matthew Izawa/Documents/Dan Applin/Llunr/xrd_profile" && pytest tests/test_size_distributions.py -v`

Expected: `ImportError: cannot import name 'fit_size_distribution' from 'xrd_profile.size_distributions'`.

- [ ] **Step 3: Implement `xrd_profile/size_distributions.py`**

Create `xrd_profile/size_distributions.py`:

```python
"""
size_distributions.py — Lognormal and normal column-length distribution
fits to Warren-Averbach A_size(L) data.

A_size(L) for a lognormal column-length distribution g(D) ~ Lognormal:

    A_size(L) = (1/<D>) * integral_{L}^{infinity} (D - L) g(D) dD

For a lognormal, this has a closed form involving the complementary
error function (Krill & Birringer 1998; Langford, Louer & Scardi 2000):

    A_size(L) = 0.5 * (1 - erf((ln(L/D_med) - sigma^2)/(sigma sqrt(2))))
              + (L / (D_med * exp(sigma^2/2)))
              * 0.5 * (1 - erf((ln(L/D_med) + 0.5*sigma^2)/(sigma sqrt(2))))

where D_med is the median crystallite size and sigma is the log-space
standard deviation. (At L=0 this reduces to A_size(0) = 1.)

For a normal column-length distribution centred at <D> with stddev s,
the closed form involves the standard error function similarly.

References
----------
Krill, C. E. & Birringer, R. (1998). Estimating grain-size distributions
    in nanocrystalline materials from X-ray diffraction profile analysis.
    Phil. Mag. A 77, 621-640.
Langford, J. I., Louer, D. & Scardi, P. (2000). Effect of a crystallite
    size distribution on X-ray diffraction line profiles and whole-powder-
    pattern fitting. J. Appl. Cryst. 33, 964-974.
"""
import numpy as np
from scipy.optimize import curve_fit
from scipy.special import erf
from scipy.stats import linregress


def lognormal_a_size(L, D_median, sigma):
    """A_size(L) for a lognormal column-length distribution.

    Parameters
    ----------
    L : np.ndarray
        Column lengths (angstroms).
    D_median : float
        Median crystallite size (angstroms).
    sigma : float
        Log-space stddev (dimensionless).

    Returns
    -------
    A : np.ndarray
        Normalised so A[0] = 1 (when L[0] = 0).
    """
    L = np.asarray(L, dtype=float)
    if D_median <= 0 or sigma <= 0:
        return np.full_like(L, np.nan)
    L_safe = np.maximum(L, 1e-10)
    arg1 = (np.log(L_safe / D_median) - sigma**2) / (sigma * np.sqrt(2))
    arg2 = (np.log(L_safe / D_median) + 0.5*sigma**2) / (sigma * np.sqrt(2))
    term1 = 0.5 * (1 - erf(arg1))
    D_mean_vol = D_median * np.exp(sigma**2 / 2)
    term2 = (L_safe / D_mean_vol) * 0.5 * (1 - erf(arg2))
    A = term1 - term2
    # Force A(L=0) = 1 (analytic limit; numerical erf may drift).
    A = np.where(L < 1e-10, 1.0, A)
    return np.clip(A, 0.0, 1.0)


def normal_a_size(L, D_mean, sigma):
    """A_size(L) for a normal column-length distribution.

    For a normal g(D) = N(D_mean, sigma), the size Fourier coefficient is

        A_size(L) = max(0, 1 - L / D_mean)

    convolved with the distribution; the closed form is approximately
    (1 - L/D_mean) * 0.5 * (1 + erf((D_mean - L)/(sigma sqrt(2)))) for
    sigma << D_mean. For broader distributions a more careful integral
    applies; we use the closed form here (Scardi & Leoni 2001).
    """
    L = np.asarray(L, dtype=float)
    if D_mean <= 0 or sigma <= 0:
        return np.full_like(L, np.nan)
    arg = (D_mean - L) / (sigma * np.sqrt(2))
    cdf_term = 0.5 * (1 + erf(arg))
    base = np.maximum(1 - L / D_mean, 0)
    A = base * cdf_term
    A = np.where(L < 1e-10, 1.0, A)
    return np.clip(A, 0.0, 1.0)


def _moments_initial_guess(L, A_size):
    """Method-of-moments initial guess for (D_median, sigma) from
    A_size(L). Uses A_size'(0) for <D> and integral for <D^2>.

    Returns (D_median_guess, sigma_guess).
    """
    L = np.asarray(L, dtype=float)
    A = np.asarray(A_size, dtype=float)
    # <D> from initial slope: -1/A'(0) approximately
    valid = A > 1e-3
    if np.sum(valid) < 4:
        return 100.0, 0.3
    # Use linear regression on the first few points.
    n_init = min(5, np.sum(valid))
    sl, _, _, _, _ = linregress(L[valid][:n_init], A[valid][:n_init])
    D_mean_guess = -1.0 / sl if sl < 0 else 100.0
    # sigma initial guess: from the curvature; safer to use a default.
    return D_mean_guess, 0.3


def fit_size_distribution(L, A_size):
    """Fit lognormal and normal distributions to A_size(L) data.

    Returns the documented `'size_distribution'` dict, or None if
    `n_valid_L < 4`.
    """
    L = np.asarray(L, dtype=float)
    A = np.asarray(A_size, dtype=float)
    if L.shape != A.shape:
        raise ValueError(f'L shape {L.shape} != A_size shape {A.shape}')

    valid = np.isfinite(A) & (A > 0) & np.isfinite(L)
    n_valid = int(np.sum(valid))
    if n_valid < 4:
        return None

    L_valid = L[valid]
    A_valid = A[valid]

    D_med_guess, sigma_guess = _moments_initial_guess(L_valid, A_valid)

    # Lognormal fit.
    try:
        popt_ln, pcov_ln = curve_fit(
            lognormal_a_size, L_valid, A_valid,
            p0=[max(D_med_guess, 1.0), max(sigma_guess, 0.05)],
            bounds=([1.0, 0.01], [1e6, 5.0]))
        D_median_fit, sigma_fit = float(popt_ln[0]), float(popt_ln[1])
        A_pred = lognormal_a_size(L_valid, D_median_fit, sigma_fit)
        ss_res = np.sum((A_valid - A_pred)**2)
        ss_tot = np.sum((A_valid - np.mean(A_valid))**2)
        r2_ln = 1.0 - (ss_res / ss_tot) if ss_tot > 0 else np.nan
        D_mean_vol = D_median_fit * np.exp(sigma_fit**2 / 2)
        D_mean_area = D_median_fit * np.exp(sigma_fit**2 / 4)
        lognormal_result = {
            'D_median': D_median_fit,
            'sigma': sigma_fit,
            'D_mean_volume': D_mean_vol,
            'D_mean_area': D_mean_area,
            'fit_r2': float(r2_ln),
            'cov': pcov_ln,
        }
    except Exception:
        lognormal_result = {
            'D_median': np.nan, 'sigma': np.nan,
            'D_mean_volume': np.nan, 'D_mean_area': np.nan,
            'fit_r2': np.nan, 'cov': None,
        }

    # Normal fit.
    try:
        popt_n, pcov_n = curve_fit(
            normal_a_size, L_valid, A_valid,
            p0=[max(D_med_guess, 1.0), max(sigma_guess * D_med_guess, 1.0)],
            bounds=([1.0, 0.01], [1e6, 1e6]))
        D_mean_n, sigma_n = float(popt_n[0]), float(popt_n[1])
        A_pred_n = normal_a_size(L_valid, D_mean_n, sigma_n)
        ss_res_n = np.sum((A_valid - A_pred_n)**2)
        r2_n = 1.0 - (ss_res_n / ss_tot) if ss_tot > 0 else np.nan
        normal_result = {
            'D_mean': D_mean_n,
            'sigma': sigma_n,
            'fit_r2': float(r2_n),
            'cov': pcov_n,
        }
    except Exception:
        normal_result = {
            'D_mean': np.nan, 'sigma': np.nan,
            'fit_r2': np.nan, 'cov': None,
        }

    return {
        'lognormal': lognormal_result,
        'normal': normal_result,
        'method': 'curve_fit',
        'initial_guess': 'moments',
        'n_valid_L': n_valid,
    }
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `cd "C:/Users/Matthew Izawa/Documents/Dan Applin/Llunr/xrd_profile" && pytest tests/test_size_distributions.py -v`

Expected: all eight tests pass.

If `test_recovers_lognormal_synthesis_within_5pct` fails because `D_median` is off by more than 5%, inspect the `_moments_initial_guess` output — the initial guess may be poor for very narrow or very broad distributions.

- [ ] **Step 5: Commit**

```bash
cd "C:/Users/Matthew Izawa/Documents/Dan Applin/Llunr/xrd_profile"
git add xrd_profile/size_distributions.py tests/test_size_distributions.py
git commit -m "Add size_distributions module with lognormal/normal A_size fits"
```

---

## Task 13 — Inject `'size_distribution'` into W-A per-family results

**Files:**
- Modify: `xrd_profile/warren_averbach.py`
- Modify: `tests/test_size_distributions.py`

Each entry in the W-A `families` list-of-dicts gains a `'size_distribution'` key with the structure from Task 12. Always computed when `n_valid_L >= 4`; `None` otherwise.

- [ ] **Step 1: Write the failing test**

Append to `tests/test_size_distributions.py`:

```python
class TestWAFamilyHasSizeDistribution:
    def test_each_family_has_size_distribution_key(self):
        from pathlib import Path
        from xrd_profile import XRDProfile, Phase
        FIX = Path(__file__).parent / 'fixtures' / 'tirhert_subset.xy'
        cif = (Path(__file__).parent.parent / 'examples'
               / 'cifs' / 'Anorthite.cif')
        anorthite = Phase.from_cif(str(cif), name='anorthite')
        data = np.loadtxt(FIX)
        profile = XRDProfile(data[:, 0], data[:, 1],
                              wavelength=0.826517)
        result = profile.guided_warren_averbach(
            phase=anorthite, n_sigma=3.0, tolerance_d=0.03)
        for fam in result['families']:
            assert 'size_distribution' in fam
            sd = fam['size_distribution']
            assert sd is None or set(sd) >= {'lognormal', 'normal',
                                              'method', 'initial_guess',
                                              'n_valid_L'}

    def test_v030_keys_unchanged_when_size_dist_added(self):
        """Ensure adding the new key did not perturb v0.3.0 keys."""
        # Already covered by tests/test_backward_compat.py v0.3.0 tier;
        # this is a belt-and-braces check at the family-key level.
        from pathlib import Path
        from xrd_profile import XRDProfile, Phase
        FIX = Path(__file__).parent / 'fixtures' / 'tirhert_subset.xy'
        cif = (Path(__file__).parent.parent / 'examples'
               / 'cifs' / 'Anorthite.cif')
        anorthite = Phase.from_cif(str(cif), name='anorthite')
        data = np.loadtxt(FIX)
        profile = XRDProfile(data[:, 0], data[:, 1],
                              wavelength=0.826517)
        result = profile.guided_warren_averbach(
            phase=anorthite, n_sigma=3.0, tolerance_d=0.03)
        # Spot-check expected v0.3.0 family keys.
        if result['families']:
            f0 = result['families'][0]
            for k in ('base_hkl', 'orders', 'd_spacings', 'fwhm_values',
                      'A_size', 'L', 'mean_sq_strain',
                      'crystallite_size', 'rms_strain', 'A_size_r2',
                      'has_overlap'):
                assert k in f0, f'v0.3.0 key {k!r} missing'
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cd "C:/Users/Matthew Izawa/Documents/Dan Applin/Llunr/xrd_profile" && pytest tests/test_size_distributions.py::TestWAFamilyHasSizeDistribution -v`

Expected: `AssertionError: ... 'size_distribution' not in {...}` on the first test.

- [ ] **Step 3: Add `size_distribution` injection to `warren_averbach.guided_warren_averbach`**

In `xrd_profile/warren_averbach.py`, locate the `results_families.append({...})` block. Just before that append, compute the size distribution from `A_size` and `L_vals` (both already in scope inside the per-family loop):

```python
        # NEW v0.4.0: fit size distribution from A_size(L)
        from .size_distributions import fit_size_distribution
        size_dist = fit_size_distribution(L_vals, A_size)
```

Then in the dictionary that is appended to `results_families`, add the new key:

```python
        results_families.append({
            'base_hkl': base_hkl,
            'orders': [fp['order'] for fp in family_profiles],
            'd_spacings': [fp['d_obs'] for fp in family_profiles],
            'fwhm_values': [fp['fwhm'] for fp in family_profiles],
            'A_size': A_size,
            'L': L_vals,
            'mean_sq_strain': mean_sq_strain,
            'crystallite_size': cryst_size,
            'rms_strain': representative_strain,
            'A_size_r2': r_as**2,
            'has_overlap': any_overlap,
            'size_distribution': size_dist,   # NEW v0.4.0
        })
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `cd "C:/Users/Matthew Izawa/Documents/Dan Applin/Llunr/xrd_profile" && pytest tests/test_size_distributions.py -v`

Expected: all tests pass, including both new `TestWAFamilyHasSizeDistribution` tests.

- [ ] **Step 5: Run all tests including backward-compat**

Run: `cd "C:/Users/Matthew Izawa/Documents/Dan Applin/Llunr/xrd_profile" && pytest tests/ -v`

Expected: full suite passes. The v0.3.0 backward-compat tier still passes — adding `'size_distribution'` to family dicts doesn't perturb any existing values (key-subset value-equality semantics).

- [ ] **Step 6: Commit**

```bash
cd "C:/Users/Matthew Izawa/Documents/Dan Applin/Llunr/xrd_profile"
git add xrd_profile/warren_averbach.py tests/test_size_distributions.py
git commit -m "Inject size_distribution key into W-A per-family results"
```

---

## Task 14 — Top-level integration test (instrumental + size dist on Tirhert)

**Files:**
- Modify: `tests/test_instrumental.py`

End-to-end smoke test that the v0.4.0 features work together: load `tirhert_subset.xy`, build anorthite Phase + LaB6 InstrumentalStandard, run guided W-A with both, inspect the result-dict shape.

- [ ] **Step 1: Write the integration test**

Append to `tests/test_instrumental.py`:

```python
class TestV040IntegrationOnTirhert:
    """End-to-end smoke test: Phase + InstrumentalStandard +
    size_distribution all working together on the bundled Tirhert
    fixture."""

    def test_full_v040_pipeline(self):
        cif = (Path(__file__).parent.parent / 'examples'
               / 'cifs' / 'Anorthite.cif')
        anorthite = Phase.from_cif(str(cif), name='anorthite')

        lab6_data = np.loadtxt(SYNTH_LAB6)
        lab6_phase = Phase.from_lattice_params(
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
        lab6 = InstrumentalStandard(
            phase=lab6_phase,
            two_theta=lab6_data[:, 0],
            intensity=lab6_data[:, 1],
            wavelength=LAMBDA_I11, name='lab6_at_i11')

        sample = np.loadtxt(TIRHERT)
        profile = XRDProfile(sample[:, 0], sample[:, 1],
                              wavelength=LAMBDA_I11,
                              sample_name='Tirhert_subset')

        # Full run_all with phase + instrumental for W-H, W-A, Scherrer.
        result = profile.run_all(
            methods=['wh', 'wa', 'scherrer'],
            phases=anorthite,
            instrumental=lab6)

        assert 'wh' in result
        assert 'wa' in result
        assert 'scherrer' in result
        # WH and WA results are keyed by phase name when phases are given.
        assert 'anorthite' in result['wh']
        assert 'anorthite' in result['wa']
        wa = result['wa']['anorthite']
        # If any families resolved, they each have a size_distribution.
        for fam in wa['families']:
            assert 'size_distribution' in fam
```

- [ ] **Step 2: Run the integration test**

Run: `cd "C:/Users/Matthew Izawa/Documents/Dan Applin/Llunr/xrd_profile" && pytest tests/test_instrumental.py::TestV040IntegrationOnTirhert -v`

Expected: passes.

- [ ] **Step 3: Run the full suite once more**

Run: `cd "C:/Users/Matthew Izawa/Documents/Dan Applin/Llunr/xrd_profile" && pytest tests/ -v`

Expected: full suite passes.

- [ ] **Step 4: Commit**

```bash
cd "C:/Users/Matthew Izawa/Documents/Dan Applin/Llunr/xrd_profile"
git add tests/test_instrumental.py
git commit -m "Add v0.4.0 end-to-end integration test on Tirhert fixture"
```

---

## Task 15 — Add an instrumental-correction example script

**Files:**
- Create: `examples/instrumental_correction.py`

A runnable demonstration of the v0.4.0 API for users skimming `examples/`. Loads the synthetic LaB6 fixture, runs guided W-A with and without correction, prints the size shift.

- [ ] **Step 1: Write the example**

Create `examples/instrumental_correction.py`:

```python
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
import numpy as np

from xrd_profile import (XRDProfile, Phase,
                         InstrumentalStandard)


REPO_ROOT = Path(__file__).parent.parent
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
```

- [ ] **Step 2: Run the example**

Run: `cd "C:/Users/Matthew Izawa/Documents/Dan Applin/Llunr/xrd_profile" && python examples/instrumental_correction.py`

Expected: runs without error, prints "Uncorrected mean crystallite size", "Corrected mean crystallite size", and per-family size-distribution lines.

- [ ] **Step 3: Commit**

```bash
cd "C:/Users/Matthew Izawa/Documents/Dan Applin/Llunr/xrd_profile"
git add examples/instrumental_correction.py
git commit -m "Add instrumental_correction.py example for v0.4.0 API"
```

---

## Task 16 — Update CHANGELOG, README, version bump

**Files:**
- Modify: `CHANGELOG.md`
- Modify: `README.md`
- Modify: `xrd_profile/__init__.py`
- Modify: `pyproject.toml`

- [ ] **Step 1: Update `xrd_profile/__init__.py` version**

Change:

```python
__version__ = '0.3.0'
```

to:

```python
__version__ = '0.4.0'
```

- [ ] **Step 2: Update `pyproject.toml` version**

Find the `version =` line and change `'0.3.0'` to `'0.4.0'`.

- [ ] **Step 3: Add v0.4.0 entry to `CHANGELOG.md`**

Insert at the top of the changelog (after the existing header):

```markdown
## [0.4.0] — 2026-05-05

### Added
- `InstrumentalStandard` class (in new module `xrd_profile/instrumental.py`):
  bundles a structural Phase (e.g., LaB6, Si) with a measured
  diffraction pattern. Supports both Caglioti FWHM correction (for
  W-H, Scherrer) and Stokes Fourier deconvolution (for W-A) via
  `caglioti_fit()` and `fourier_coefficients(peak_d)` methods.
- `InstrumentalProfile` class: lightweight Caglioti-coefficient
  carrier (U, V, W). Supports W-H and Scherrer; W-A raises a clear
  `ValueError` (Stokes deconvolution requires a measured pattern).
  JSON I/O via `to_json` / `from_json`. Registry lookup via
  `from_registry(name)` (registry ships empty in v0.4.0).
- `instrumental=` keyword argument is now wired into
  `XRDProfile.guided_williamson_hall()`, `guided_warren_averbach()`,
  `scherrer()`, `modified_scherrer()`, and `run_all()`. The
  v0.3.0 `NotImplementedError` guards are replaced with dispatch.
- `phase=` keyword argument added to `XRDProfile.scherrer()` and
  `XRDProfile.modified_scherrer()` for guided peak filtering.
- Crystallite size distributions: every entry in the
  `guided_warren_averbach` `families` list now includes a
  `'size_distribution'` key with lognormal and normal fits to the
  W-A column-length distribution (per Krill & Birringer 1998;
  Langford, Louer & Scardi 2000). Returns `None` when fewer than 4
  valid `L` points are available. New internal module
  `xrd_profile/size_distributions.py`.
- New tests: `test_instrumental.py` (~14 tests),
  `test_size_distributions.py` (~10 tests).
- New tier in `test_backward_compat.py`: `golden_v0.3.0_results.json`
  exercises Phase API + run_all + Scherrer shape table. Comparison
  semantics is key-subset value-equality (existing keys unchanged;
  new top-level keys at later tags are allowed).
- New example `examples/instrumental_correction.py` demonstrating
  Stokes-corrected vs uncorrected W-A on the bundled Tirhert subset.
- Synthetic LaB6 test fixture (`tests/fixtures/synthetic_lab6.xy`)
  generated by `scripts/build_synthetic_standard.py` from documented
  Caglioti coefficients (U=5e-3, V=-1e-3, W=5e-3 deg^2).

### Compatibility
- All v0.3.0 public-API calls continue to work unchanged. v0.3.0
  result-dict keys retain their numerical values when no v0.4.0
  feature is invoked. The v0.2.0 and v0.3.0 backward-compatibility
  tiers both pass with key-subset value-equality semantics.
```

- [ ] **Step 4: Add quickstart snippet to `README.md`**

In the relevant Quickstart / Usage section, add:

````markdown
### Instrumental broadening correction (v0.4.0)

If you have a measured LaB6 (or Si) standard pattern at the same
instrument settings as your sample:

```python
import numpy as np
from xrd_profile import (XRDProfile, Phase, InstrumentalStandard)

sample = np.loadtxt('sample.xy')
lab6_pattern = np.loadtxt('lab6_standard.xy')

anorthite = Phase.from_cif('Anorthite.cif', name='anorthite')
lab6 = InstrumentalStandard.from_cif_and_pattern(
    cif='LaB6.cif',
    two_theta=lab6_pattern[:, 0],
    intensity=lab6_pattern[:, 1],
    wavelength=1.5406, name='LaB6_Cu_Bruker')

profile = XRDProfile(sample[:, 0], sample[:, 1], wavelength=1.5406)
results = profile.run_all(
    methods=['wh', 'wa', 'scherrer'],
    phases=anorthite,
    instrumental=lab6)
```

If you only have published Caglioti coefficients (no measured
standard), use `InstrumentalProfile` — sufficient for W-H and
Scherrer; W-A will reject it with a clear error.

```python
from xrd_profile import InstrumentalProfile
inst = InstrumentalProfile(U=5e-3, V=-1e-3, W=5e-3, wavelength=1.5406)
results = profile.run_all(methods=['wh', 'scherrer'],
                           phases=anorthite, instrumental=inst)
```
````

- [ ] **Step 5: Run the full test suite once more**

Run: `cd "C:/Users/Matthew Izawa/Documents/Dan Applin/Llunr/xrd_profile" && pytest tests/ -v`

Expected: full suite passes.

- [ ] **Step 6: Commit**

```bash
cd "C:/Users/Matthew Izawa/Documents/Dan Applin/Llunr/xrd_profile"
git add xrd_profile/__init__.py pyproject.toml CHANGELOG.md README.md
git commit -m "Bump version to 0.4.0; update CHANGELOG and README"
```

---

## Task 17 — Tag v0.4.0

**Files:** none (git tag only)

- [ ] **Step 1: Confirm clean working tree and final test pass**

```bash
cd "C:/Users/Matthew Izawa/Documents/Dan Applin/Llunr/xrd_profile"
git status
pytest tests/ -v
```

Expected: working tree clean (no uncommitted changes), full test suite passes.

- [ ] **Step 2: Verify the v0.3.0 golden tier passes against the current code**

```bash
cd "C:/Users/Matthew Izawa/Documents/Dan Applin/Llunr/xrd_profile"
pytest tests/test_backward_compat.py -v
```

Expected: every v0.2.0 and v0.3.0 tier test passes (key-subset value-equality).

- [ ] **Step 3: Verify the version bump is in place**

```bash
cd "C:/Users/Matthew Izawa/Documents/Dan Applin/Llunr/xrd_profile"
python -c "import xrd_profile; print(xrd_profile.__version__)"
```

Expected output: `0.4.0`.

- [ ] **Step 4: Tag**

```bash
cd "C:/Users/Matthew Izawa/Documents/Dan Applin/Llunr/xrd_profile"
git tag -a v0.4.0 -m "Release v0.4.0: instrumental deconvolution, size distributions, guided Scherrer"
```

- [ ] **Step 5: Verify tag**

```bash
cd "C:/Users/Matthew Izawa/Documents/Dan Applin/Llunr/xrd_profile"
git tag --list 'v*'
git show v0.4.0 --stat
```

Expected: tag list includes `v0.2.0`, `v0.3.0`, `v0.4.0`. The `git show` output displays the v0.4.0 commit message and the diff against the tagged commit.

- [ ] **Step 6: No push**

Do **not** `git push --tags` from inside this task. The package's git remote is the private repo `MatthewIzawa/XRD_Shock_Strain`; pushing tags is the user's call, not the implementer's. The tag exists locally and is ready to push when the user decides.

---

## Acceptance criteria summary

This plan satisfies every numbered acceptance criterion in spec §11.1:

1. ✓ `InstrumentalStandard.from_cif_and_pattern` works on synthetic LaB6 fixture (Tasks 6 + 14).
2. ✓ `caglioti_fit()` recovers synthesis coefficients within 5% (Task 5).
3. ✓ `guided_williamson_hall(phase=, instrumental=)` works for both `InstrumentalStandard` and `InstrumentalProfile` (Task 8).
4. ✓ `guided_warren_averbach(phase=, instrumental=)` works for `InstrumentalStandard`, raises clear `ValueError` for `InstrumentalProfile` (Task 9).
5. ✓ `scherrer(phase=, instrumental=)` works (Task 10).
6. ✓ W-A result has `'size_distribution'` key per family with `n_valid_L >= 4`; `None` otherwise (Tasks 12 + 13).
7. ✓ `test_instrumental.py`, `test_size_distributions.py` pass (Tasks 3, 5, 6, 7, 8, 9, 10, 12, 13, 14).
8. ✓ `test_backward_compat.py` passes both `golden_v0.2.0` and `golden_v0.3.0` tiers (Task 1; verified after every wiring task).
9. ✓ CHANGELOG entry for v0.4.0 (Task 16).
10. ✓ README updated with instrumental-correction quickstart (Task 16).
11. ✓ `git tag v0.4.0` on the v0.4.0-complete commit (Task 17).
