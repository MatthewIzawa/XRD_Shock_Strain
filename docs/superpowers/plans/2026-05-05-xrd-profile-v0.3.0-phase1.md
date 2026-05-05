# xrd_profile v0.3.0 (Phase 1) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Land xrd_profile v0.3.0 — adds `Phase` class, `phase=` kwarg on guided W-H and W-A, `run_all` convenience helper, Scherrer K and shape factor exposure — strictly additive, preserving v0.2.0 numerical behavior byte-for-byte.

**Architecture:** New `phases.py` module owns the `Phase` abstraction over `pymatgen.Structure`. Existing analysis modules (`williamson_hall.py`, `warren_averbach.py`, `pdf.py`, `scherrer.py`) are unchanged or only minimally extended via additive kwargs. The `XRDProfile` class gains `phase=` routing and a `run_all()` dispatcher; legacy array-based call forms (`ref_d=`, `ref_peaks=`) keep working unchanged. A reserved `instrumental=` kwarg raises `NotImplementedError` if used; Phase 2 wires it up. A numerical regression test against frozen v0.2.0 outputs guards against drift.

**Tech Stack:** Python ≥ 3.8, numpy, scipy, matplotlib, pymatgen (optional via `[cif]` extra), pytest 8.x.

**Spec:** [docs/superpowers/specs/2026-05-05-xrd-profile-v1-phase1-design.md](../specs/2026-05-05-xrd-profile-v1-phase1-design.md) (commit `c431431`).

**Repo location:** `C:\Users\Matthew Izawa\Documents\Dan Applin\Llunr\xrd_profile\`. All paths in this plan are relative to that repo root unless absolute.

---

## File Manifest

**New files:**
- `xrd_profile/phases.py` — Phase class + lazy import helper + standalone build_reference_peaks
- `examples/cifs/Forsterite.cif`, `Anorthite.cif`, `Pigeonite.cif`, `Quartz.cif`, `Hematite.cif`
- `examples/cifs/SOURCES.md` — CIF provenance and citations
- `examples/data/tirhert_subset.xy` — bundled small XRD pattern
- `examples/data/README.md` — pointer to full data
- `examples/legacy/lab_lunar_meteorite.py`, `synchrotron_low_shock.py`, `synchrotron_high_shock.py` — verbatim copies of v0.2.0 originals
- `examples/multi_phase_olivine.py` — new canonical demo
- `tests/fixtures/tirhert_subset.xy`
- `tests/fixtures/golden_v0.2.0_results.json`
- `tests/test_phases.py`
- `tests/test_run_all.py`
- `tests/test_scherrer.py`
- `tests/test_backward_compat.py`
- `scripts/regenerate_goldens.py`
- `CHANGELOG.md`

**Modified files:**
- `xrd_profile/__init__.py` — version bump 0.2.0 → 0.3.0; add Phase, build_reference_peaks, SCHERRER_K_FOR_SHAPE exports
- `xrd_profile/profile.py` — add `XRDProfile.run_all()`; add `phase=` and reserved `instrumental=` kwargs to `guided_williamson_hall()` and `guided_warren_averbach()`
- `xrd_profile/scherrer.py` — add `shape` kwarg + `SCHERRER_K_FOR_SHAPE` table; change `K` default to `None` sentinel preserving v0.2.0 behavior
- `pyproject.toml` — add `[project.optional-dependencies] cif = ["pymatgen>=2023.0"]`; bump version
- `setup.py` — bump version (if it tracks version separately)
- `README.md` — new quickstart with `[cif]` install and Phase pattern, run_all section, link to multi_phase_olivine example
- `examples/lab_lunar_meteorite.py`, `synchrotron_low_shock.py`, `synchrotron_high_shock.py` — rewrite to use Phase

**Git operations:**
- `git tag v0.2.0` on the pre-Phase-1 HEAD (commit `eec5d3f`, the last code commit before the spec doc)
- `git tag v0.3.0` on the Phase-1-complete commit at end of plan

---

## Task ordering

```
Phase A (setup + safety net)
  Task 1: Tag v0.2.0 + CHANGELOG skeleton
  Task 2: Bundle test fixture + golden generator + regression test

Phase B (core API additions, TDD)
  Task 3: Phase class + lazy import + build_reference_peaks function
  Task 4: phase= and reserved instrumental= kwargs on guided_williamson_hall
  Task 5: phase= and reserved instrumental= kwargs on guided_warren_averbach
  Task 6: Scherrer K + shape factor extension
  Task 7: XRDProfile.run_all() helper

Phase C (examples + fixtures)
  Task 8: Source example CIFs + SOURCES.md
  Task 9: Bundle data fixture in examples/data/
  Task 10: Copy v0.2.0 examples to legacy/
  Task 11: Rewrite existing examples to use Phase
  Task 12: Write multi_phase_olivine.py

Phase D (release plumbing)
  Task 13: Update __init__.py exports
  Task 14: Update pyproject.toml ([cif] extra, version bump)
  Task 15: Update README
  Task 16: Update CHANGELOG, version bump in __init__.py, tag v0.3.0
```

---

# Phase A: Setup & Safety Net

## Task 1: Tag v0.2.0 + CHANGELOG skeleton

**Files:**
- Create: `CHANGELOG.md`
- Git tag: `v0.2.0`

- [ ] **Step 1: Confirm pre-Phase-1 HEAD**

Run: `git log --oneline -5`

Expected output (the last code commit before the spec doc should be `eec5d3f`):
```
c431431 Add Phase 1 design doc for v0.3.0 / v1.0 generalisation
eec5d3f Add enhanced PDF analysis: Chebyshev background, sine transform, Lorch modification, peak fitting (v0.2.0)
f28ccae Add cross-phase overlap rejection, peak quality scoring, weighted regression, reliability classification, and CSV validation export to guided Williamson-Hall analysis.
0a76987 Initial release: xrd_profile v0.1.0
```

- [ ] **Step 2: Tag v0.2.0 on the eec5d3f commit**

Run: `git tag -a v0.2.0 eec5d3f -m "Release v0.2.0: enhanced PDF + cross-phase guided W-H"`

Verify: `git tag --list "v*"`
Expected output: `v0.2.0`

- [ ] **Step 3: Create CHANGELOG.md skeleton**

Create `CHANGELOG.md` with this exact content:

```markdown
# Changelog

All notable changes to xrd_profile are documented here. Format
follows [Keep a Changelog](https://keepachangelog.com/), versioning
follows [SemVer](https://semver.org/).

## [Unreleased] — Phase 1 (target v0.3.0)

### Added
- (in progress — see docs/superpowers/plans/2026-05-05-xrd-profile-v0.3.0-phase1.md)

## [0.2.0] — 2026-04

### Added
- Enhanced PDF pipeline: iterative Chebyshev background subtraction in
  Q-space, sine Fourier transform with optional Lorch modification, PDF
  peak detection and Gaussian first-shell fitting.
- Cross-phase overlap rejection, peak quality scoring, weighted
  regression, and reliability classification on guided W-H.

## [0.1.0] — Initial release

### Added
- Williamson-Hall (conventional and reciprocal-space).
- Warren-Averbach with harmonic peak families.
- Scherrer (standard and modified).
- Reference-guided peak detection, automated zero-point offset.
- XRDProfile unified interface.
```

- [ ] **Step 4: Commit CHANGELOG**

Run:
```bash
git add CHANGELOG.md
git commit -m "Add CHANGELOG.md, tag v0.2.0

Tags eec5d3f as v0.2.0 (last commit before v0.3.0 / v1.0 work).
CHANGELOG documents the v0.1.0 and v0.2.0 releases retrospectively;
v0.3.0 entries fill in as Phase 1 lands."
```

---

## Task 2: Bundle test fixture + golden generator + regression test

This task installs the safety net BEFORE any code changes. The flow:
(a) bundle a small XRD pattern; (b) write a script that runs every
v0.2.0 public function on it and saves outputs to JSON; (c) run that
script on current code (= v0.2.0 behavior, since no code has changed
yet); (d) write the regression test that asserts the JSON values stay
constant going forward.

**Files:**
- Create: `tests/fixtures/tirhert_subset.xy`
- Create: `tests/fixtures/golden_v0.2.0_results.json`
- Create: `scripts/regenerate_goldens.py`
- Create: `tests/test_backward_compat.py`

- [ ] **Step 1: Create tests/fixtures/ directory**

Run: `mkdir -p tests/fixtures`

- [ ] **Step 2: Generate the bundled XRD subset**

Source the Tirhert synchrotron pattern from the user's data location.
The full file is at:
`C:\Users\Matthew Izawa\Desktop\111 Backup 20220530\transfer\IPM\2018\ee17803-1\processing\Tirhert_summed_0001.xye`

Create a Python script at `scripts/_make_subset.py` (will be deleted
after use):

```python
"""One-shot script: subset the Tirhert .xye to a small .xy fixture."""
import numpy as np

src = (r'C:\Users\Matthew Izawa\Desktop\111 Backup 20220530'
       r'\transfer\IPM\2018\ee17803-1\processing\Tirhert_summed_0001.xye')
data = np.loadtxt(src)
tt, I = data[:, 0], data[:, 1]
mask = (tt >= 10) & (tt <= 80)
tt_sub = tt[mask]
I_sub = I[mask]
# Downsample to ~1000 points (every Nth)
N = max(1, len(tt_sub) // 1000)
tt_d = tt_sub[::N]
I_d = I_sub[::N]
np.savetxt('tests/fixtures/tirhert_subset.xy',
           np.column_stack([tt_d, I_d]),
           fmt='%.6f %.4f',
           header='2theta_deg intensity', comments='# ')
print(f'Wrote {len(tt_d)} points to tests/fixtures/tirhert_subset.xy')
```

Run: `python scripts/_make_subset.py`

Expected: `Wrote ~1000 points to tests/fixtures/tirhert_subset.xy`

Verify file size:
```bash
ls -lh tests/fixtures/tirhert_subset.xy
```
Expected: 30-50 KB.

Delete the one-shot script:
```bash
rm scripts/_make_subset.py
```

- [ ] **Step 3: Write scripts/regenerate_goldens.py**

Create the directory: `mkdir -p scripts`

Create `scripts/regenerate_goldens.py`:

```python
"""
Regenerate tests/fixtures/golden_v0.2.0_results.json from the bundled
Tirhert subset using only v0.2.0-public API calls.

Run when v0.2.0 numerical behavior intentionally changes (rare —
typically only for documented bug fixes). Each regeneration must be
accompanied by an explicit reasoning in the commit message.

Usage: python scripts/regenerate_goldens.py
"""
import json
from pathlib import Path
import numpy as np

from xrd_profile import (XRDProfile, two_theta_to_d,
                         scherrer, modified_scherrer,
                         compute_pdf_sine, fit_first_pdf_peak,
                         estimate_fwhm_simple)

LAMBDA_I11 = 0.826517
FIXTURE = Path(__file__).parent.parent / 'tests' / 'fixtures' / 'tirhert_subset.xy'
OUT = Path(__file__).parent.parent / 'tests' / 'fixtures' / 'golden_v0.2.0_results.json'

# Fixed reference d-spacing list for guided W-H, derived from anorthite at
# I11 wavelength. Values frozen here so regeneration is deterministic.
ANORTHITE_REF_D = [
    3.20, 3.18, 3.65, 4.04, 6.41, 5.69, 3.74, 3.21, 4.04, 2.94,
]

# Fixed reference peak list for guided W-A. (d, two_theta, intensity, h, k, l)
ANORTHITE_REF_PEAKS = [
    {'d': 3.20, 'two_theta': 14.84, 'intensity': 100.0, 'h': 0, 'k': 4, 'l': 0},
    {'d': 3.18, 'two_theta': 14.93, 'intensity':  85.0, 'h': 2, 'k': 0, 'l': -2},
    {'d': 6.41, 'two_theta':  7.39, 'intensity':  60.0, 'h': 0, 'k': 2, 'l': 0},
    {'d': 4.04, 'two_theta': 11.74, 'intensity':  50.0, 'h': 0, 'k': 0, 'l': 2},
    {'d': 3.65, 'two_theta': 13.00, 'intensity':  45.0, 'h': 1, 'k': 3, 'l': 0},
]


def to_serializable(obj):
    """Convert numpy types and arrays to JSON-friendly form."""
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, (np.floating, np.integer)):
        return float(obj)
    if isinstance(obj, (bool, np.bool_)):
        return bool(obj)
    if isinstance(obj, dict):
        return {k: to_serializable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [to_serializable(x) for x in obj]
    return obj


def main():
    data = np.loadtxt(FIXTURE)
    tt, I = data[:, 0], data[:, 1]
    profile = XRDProfile(tt, I, wavelength=LAMBDA_I11, sample_name='Tirhert_subset')

    results = {
        'metadata': {
            'fixture': str(FIXTURE.name),
            'wavelength': LAMBDA_I11,
            'n_points': int(len(tt)),
            'tt_range': [float(tt.min()), float(tt.max())],
        },
        'guided_williamson_hall': profile.guided_williamson_hall(
            np.array(ANORTHITE_REF_D), n_sigma=3.0, tolerance_d=0.03),
        'guided_warren_averbach': profile.guided_warren_averbach(
            ANORTHITE_REF_PEAKS, n_sigma=3.0, tolerance_d=0.03),
        'compute_pdf_sine': None,
        'scherrer_default': None,
        'modified_scherrer_default': None,
    }

    r, G_r, Q_max = compute_pdf_sine(tt, I, LAMBDA_I11, cheby_order=15, lorch=True)
    results['compute_pdf_sine'] = {
        'r': r.tolist(),
        'G_r': G_r.tolist(),
        'Q_max': float(Q_max),
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

    serial = to_serializable(results)
    OUT.write_text(json.dumps(serial, indent=2))
    print(f'Wrote {OUT}')
    print(f'Keys: {list(serial.keys())}')


if __name__ == '__main__':
    main()
```

- [ ] **Step 4: Run regenerate_goldens.py to produce the golden file**

Run: `python scripts/regenerate_goldens.py`

Expected:
```
Wrote .../tests/fixtures/golden_v0.2.0_results.json
Keys: ['metadata', 'guided_williamson_hall', 'guided_warren_averbach', 'compute_pdf_sine', 'scherrer_default', 'modified_scherrer_default']
```

If `guided_williamson_hall` returns very few peaks (e.g., the
hardcoded `ANORTHITE_REF_D` list above doesn't match well), inspect
the result and either widen `tolerance_d=` or extend the ref_d list
until the call succeeds with at least 3 matched peaks. Re-run.

Verify file exists and is non-trivial:
```bash
ls -lh tests/fixtures/golden_v0.2.0_results.json
```

- [ ] **Step 5: Write tests/test_backward_compat.py**

Create:

```python
"""
Numerical regression test against frozen v0.2.0 outputs.

Asserts that the v0.2.0 array-based public API produces identical
output (within tight tolerances) to what was generated against the
v0.2.0 codebase. If this test fails after a code change, that change
has perturbed v0.2.0 default behavior — the strict-additive policy of
Phase 1 forbids that. Either fix the code or regenerate the golden
file with explicit reasoning (see scripts/regenerate_goldens.py).
"""
import json
from pathlib import Path

import numpy as np
import pytest

from xrd_profile import (XRDProfile, scherrer, modified_scherrer,
                         compute_pdf_sine, estimate_fwhm_simple)

FIXTURE_DIR = Path(__file__).parent / 'fixtures'
LAMBDA_I11 = 0.826517

# Same constants as scripts/regenerate_goldens.py — kept in sync.
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


@pytest.fixture(scope='module')
def golden():
    """Load the frozen v0.2.0 outputs."""
    return json.loads((FIXTURE_DIR / 'golden_v0.2.0_results.json').read_text())


@pytest.fixture(scope='module')
def pattern():
    """Load the bundled tirhert subset."""
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
        sizes = scherrer(fwhm, positions, LAMBDA_I11)  # default K
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
```

- [ ] **Step 6: Run the regression test (should PASS on current code)**

Run: `pytest tests/test_backward_compat.py -v`

Expected: All tests PASS (we are still on v0.2.0-equivalent code; the
test verifies the safety net works before any changes are made).

If any test fails on this run, the golden file generation has a bug —
investigate before continuing.

- [ ] **Step 7: Commit**

```bash
git add tests/fixtures/tirhert_subset.xy tests/fixtures/golden_v0.2.0_results.json scripts/regenerate_goldens.py tests/test_backward_compat.py
git commit -m "Add v0.2.0 regression baseline: fixture + golden + test

Bundles a ~30 KB Tirhert synchrotron subset, generates frozen
golden outputs from v0.2.0 public API calls, and adds the regression
test that enforces no numerical drift across Phase 1 additions.
scripts/regenerate_goldens.py is committed for the rare case of
intentional numerical updates."
```

---

# Phase B: Core API Additions (TDD)

## Task 3: Phase class + lazy import + build_reference_peaks

**Files:**
- Create: `xrd_profile/phases.py`
- Create: `tests/test_phases.py`

- [ ] **Step 1: Write failing test for missing-pymatgen ImportError**

Create `tests/test_phases.py`:

```python
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
        # Force re-import of phases module so the lazy import sees the patch
        if 'xrd_profile.phases' in sys.modules:
            del sys.modules['xrd_profile.phases']
        from xrd_profile.phases import Phase
        with pytest.raises(ImportError) as exc_info:
            Phase.from_cif('nonexistent.cif')
        assert 'pip install xrd_profile[cif]' in str(exc_info.value)
```

- [ ] **Step 2: Run test, expect FAIL**

Run: `pytest tests/test_phases.py::test_phase_from_cif_without_pymatgen_raises_clear_error -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'xrd_profile.phases'`.

- [ ] **Step 3: Implement minimal phases.py with the lazy-import helper**

Create `xrd_profile/phases.py`:

```python
"""
Phase abstraction over pymatgen.Structure for guided peak detection.

The Phase class wraps a crystal structure and exposes reference-peak
generation through get_ref_peaks() and get_ref_d(). pymatgen is an
optional dependency; if missing, a clear ImportError is raised that
points to the [cif] extra install command.
"""
from __future__ import annotations

from pathlib import Path
from typing import List, Tuple, Optional, Sequence

import numpy as np

from .conversions import two_theta_to_d


def _require_pymatgen():
    """Lazy import helper. Raises ImportError with install instructions."""
    try:
        from pymatgen.core.structure import Structure, Lattice
        from pymatgen.analysis.diffraction.xrd import XRDCalculator
        return Structure, Lattice, XRDCalculator
    except ImportError as e:
        raise ImportError(
            'pymatgen is required for Phase.from_cif and '
            'Phase.from_lattice_params. '
            'Install with: pip install xrd_profile[cif]'
        ) from e


class Phase:
    """A crystalline phase with a known structure."""

    def __init__(self, structure, name: str = ''):
        self.structure = structure
        self.name = name

    @classmethod
    def from_cif(cls, path, name: Optional[str] = None) -> 'Phase':
        """Load from a CIF file. Default name = file stem."""
        Structure, _, _ = _require_pymatgen()
        structure = Structure.from_file(str(path))
        if name is None:
            name = Path(path).stem
        return cls(structure, name=name)

    @classmethod
    def from_lattice_params(cls, a: float, b: float, c: float,
                             alpha: float, beta: float, gamma: float,
                             species: Sequence,
                             coords: Sequence,
                             name: str = '') -> 'Phase':
        """Build inline from lattice parameters and atomic positions."""
        Structure, Lattice, _ = _require_pymatgen()
        lat = Lattice.from_parameters(a, b, c, alpha, beta, gamma)
        structure = Structure(lat, species, coords)
        return cls(structure, name=name)

    def get_ref_peaks(self, wavelength: float,
                      two_theta_range: Tuple[float, float] = (5, 90),
                      min_intensity: float = 3.0) -> List[dict]:
        """List of {d, two_theta, intensity, h, k, l} dicts above min_intensity."""
        _, _, XRDCalculator = _require_pymatgen()
        calc = XRDCalculator(wavelength=wavelength)
        pattern = calc.get_pattern(self.structure,
                                    two_theta_range=two_theta_range)
        peaks: List[dict] = []
        for i in range(len(pattern.x)):
            if pattern.y[i] < min_intensity:
                continue
            hkl = pattern.hkls[i][0]['hkl']
            d = two_theta_to_d(pattern.x[i], wavelength)
            peaks.append({
                'd': float(d),
                'two_theta': float(pattern.x[i]),
                'intensity': float(pattern.y[i]),
                'h': int(hkl[0]),
                'k': int(hkl[1]),
                'l': int(hkl[2]),
            })
        return peaks

    def get_ref_d(self, wavelength: float,
                   two_theta_range: Tuple[float, float] = (5, 90),
                   min_intensity: float = 3.0,
                   sorted_by_intensity: bool = True) -> np.ndarray:
        """Reference d-spacings, optionally sorted by intensity descending."""
        peaks = self.get_ref_peaks(wavelength, two_theta_range, min_intensity)
        if sorted_by_intensity:
            peaks.sort(key=lambda p: -p['intensity'])
        return np.array([p['d'] for p in peaks], dtype=float)

    def __repr__(self) -> str:
        try:
            comp = self.structure.composition.reduced_formula
        except Exception:
            comp = '?'
        return f'<Phase {self.name!r}: {comp}>'


def build_reference_peaks(structure, wavelength: float,
                           two_theta_range: Tuple[float, float] = (5, 90),
                           min_intensity: float = 3.0) -> List[dict]:
    """Standalone reference-peak builder for users who already have a
    pymatgen Structure. Equivalent to Phase(structure).get_ref_peaks(...)."""
    return Phase(structure, name='').get_ref_peaks(
        wavelength, two_theta_range, min_intensity)
```

- [ ] **Step 4: Run the missing-pymatgen test, expect PASS**

Run: `pytest tests/test_phases.py::test_phase_from_cif_without_pymatgen_raises_clear_error -v`
Expected: PASS.

- [ ] **Step 5: Add tests for from_lattice_params and ref-peak generation**

Append to `tests/test_phases.py`:

```python
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
```

- [ ] **Step 6: Run all phase tests, expect PASS**

Run: `pytest tests/test_phases.py -v`
Expected: all tests PASS.

- [ ] **Step 7: Run regression test, expect PASS (no v0.2.0 drift)**

Run: `pytest tests/test_backward_compat.py -v`
Expected: all tests PASS.

- [ ] **Step 8: Commit**

```bash
git add xrd_profile/phases.py tests/test_phases.py
git commit -m "Add Phase class with from_cif and from_lattice_params

Phase wraps a pymatgen Structure and exposes get_ref_peaks() and
get_ref_d() for guided peak detection. pymatgen is loaded lazily;
calling Phase.from_cif without it raises an ImportError that points
to the [cif] extra install. The standalone build_reference_peaks
function provides the same conversion for users who already have a
pymatgen Structure in hand."
```

---

## Task 4: phase= and reserved instrumental= kwargs on guided_williamson_hall

**Files:**
- Modify: `xrd_profile/profile.py`
- Modify: `tests/test_phases.py` (add tests for the integration)

- [ ] **Step 1: Write failing tests for the new kwargs**

Append to `tests/test_phases.py`:

```python
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

    def test_instrumental_kwarg_raises_not_implemented(self):
        profile, an = self._make_profile_and_phase()
        with pytest.raises(NotImplementedError, match='Phase 2'):
            profile.guided_williamson_hall(phase=an, instrumental='anything')
```

- [ ] **Step 2: Run new tests, expect FAIL**

Run: `pytest tests/test_phases.py::TestXRDProfileGuidedWilliamsonHallPhaseKwarg -v`
Expected: FAIL with `TypeError: guided_williamson_hall() got an unexpected keyword argument 'phase'` or similar.

- [ ] **Step 3: Read the current guided_williamson_hall signature**

Read `xrd_profile/profile.py` lines 71-120 to see the current method
body. The signature is:
```python
def guided_williamson_hall(self, ref_d, tolerance_d=0.03,
                           n_sigma=3.0, min_fwhm_steps=3,
                           correct_offset=True,
                           other_phase_d=None,
                           ...
                           **kwargs):
```

`ref_d` is currently positional (required). To remain strictly
additive, change it to `ref_d=None` (defaulting to None) and add
`phase=` and `instrumental=` as keyword-only after a `*`.

- [ ] **Step 4: Modify guided_williamson_hall in profile.py**

Edit `xrd_profile/profile.py`. Locate the method `guided_williamson_hall`
(starts around line 71). Replace its signature and prepend the
phase/instrumental routing logic. Final method shape:

```python
def guided_williamson_hall(self, ref_d=None, tolerance_d=0.03,
                           n_sigma=3.0, min_fwhm_steps=3,
                           correct_offset=True,
                           other_phase_d=None,
                           other_phase_names=None,
                           overlap_tol_deg=0.15,
                           min_quality=0.3,
                           quality_weights=None,
                           weighted_regression=True,
                           sample_flags=None,
                           export_path=None,
                           *,
                           phase=None,
                           instrumental=None,
                           **kwargs):
    """
    Reference-guided Williamson-Hall.

    Parameters
    ----------
    ref_d : np.ndarray or None
        Reference d-spacings sorted by decreasing intensity.
        Mutually exclusive with phase=.
    phase : Phase or None
        New in v0.3.0. If provided, ref_d is computed from
        phase.get_ref_d(wavelength, two_theta_range=data_range).
        Mutually exclusive with ref_d=.
    instrumental : reserved for Phase 2 / v1.0
        Pass None (default). Any other value raises NotImplementedError.

    [other parameters unchanged from v0.2.0]
    """
    if phase is not None and ref_d is not None:
        raise ValueError(
            'pass either ref_d or phase, not both')
    if phase is not None:
        tt_range = (float(self.two_theta.min()),
                    float(self.two_theta.max()))
        ref_d = phase.get_ref_d(self.wavelength,
                                 two_theta_range=tt_range)
    if instrumental is not None:
        raise NotImplementedError(
            'instrumental= is reserved for Phase 2 / v1.0; '
            'see xrd_profile roadmap')
    if ref_d is None:
        raise ValueError(
            'must pass either ref_d or phase')

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
        **kwargs
    )
```

Use the Edit tool with `old_string` set to the full current method
(read it first to copy exactly, including indentation) and `new_string`
set to the version above. Match indentation exactly.

- [ ] **Step 5: Run new tests, expect PASS**

Run: `pytest tests/test_phases.py::TestXRDProfileGuidedWilliamsonHallPhaseKwarg -v`
Expected: all 3 tests PASS.

- [ ] **Step 6: Run regression test, expect PASS (no v0.2.0 drift)**

Run: `pytest tests/test_backward_compat.py -v`
Expected: all PASS. The legacy `ref_d=` positional call still works
because v0.2.0 callers pass it as the first positional argument and
the new signature still accepts it positionally.

- [ ] **Step 7: Commit**

```bash
git add xrd_profile/profile.py tests/test_phases.py
git commit -m "Add phase= and reserved instrumental= kwargs on guided_williamson_hall

phase= accepts a Phase object and internally computes ref_d using the
profile's data 2-theta range. Mutually exclusive with the legacy ref_d=
positional/keyword argument. instrumental= is reserved for Phase 2;
raises NotImplementedError if non-None. v0.2.0 callers using ref_d=
positionally are unaffected (regression test confirms)."
```

---

## Task 5: phase= and reserved instrumental= kwargs on guided_warren_averbach

Same structure as Task 4 but for `guided_warren_averbach`. The
difference: this method takes `ref_peaks` (full hkl-keyed list) rather
than `ref_d`, so the integration calls `phase.get_ref_peaks(...)`.

**Files:**
- Modify: `xrd_profile/profile.py`
- Modify: `tests/test_phases.py`

- [ ] **Step 1: Write failing tests**

Append to `tests/test_phases.py`:

```python
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
```

- [ ] **Step 2: Run, expect FAIL**

Run: `pytest tests/test_phases.py::TestXRDProfileGuidedWarrenAverbachPhaseKwarg -v`
Expected: FAIL with unexpected-kwarg error.

- [ ] **Step 3: Modify guided_warren_averbach in profile.py**

Edit `xrd_profile/profile.py`. Locate the method `guided_warren_averbach`
(currently around lines 122-148). Same pattern as Task 4: change
`ref_peaks` to `ref_peaks=None`, add keyword-only `phase=None,
instrumental=None`, and prepend the routing logic:

```python
def guided_warren_averbach(self, ref_peaks=None, tolerance_d=0.03,
                            n_sigma=3.0, min_fwhm_steps=3,
                            correct_offset=True, n_coeffs=20,
                            width_fwhm=6.0, require_clean=False,
                            *,
                            phase=None,
                            instrumental=None):
    """
    Reference-guided Warren-Averbach.

    Parameters
    ----------
    ref_peaks : list of dict or None
        Reference peaks with {d, two_theta, intensity, h, k, l}.
        Mutually exclusive with phase=.
    phase : Phase or None
        New in v0.3.0. If provided, ref_peaks is computed from
        phase.get_ref_peaks(wavelength, two_theta_range=data_range).
        Mutually exclusive with ref_peaks=.
    instrumental : reserved for Phase 2 / v1.0.

    [other parameters unchanged from v0.2.0]
    """
    if phase is not None and ref_peaks is not None:
        raise ValueError('pass either ref_peaks or phase, not both')
    if phase is not None:
        tt_range = (float(self.two_theta.min()),
                    float(self.two_theta.max()))
        ref_peaks = phase.get_ref_peaks(self.wavelength,
                                         two_theta_range=tt_range)
    if instrumental is not None:
        raise NotImplementedError(
            'instrumental= is reserved for Phase 2 / v1.0; '
            'see xrd_profile roadmap')
    if ref_peaks is None:
        raise ValueError('must pass either ref_peaks or phase')

    from .warren_averbach import guided_warren_averbach
    return guided_warren_averbach(
        self.two_theta, self.intensity, ref_peaks, self.wavelength,
        tolerance_d=tolerance_d, n_sigma=n_sigma,
        min_fwhm_steps=min_fwhm_steps, correct_offset=correct_offset,
        n_coeffs=n_coeffs, width_fwhm=width_fwhm,
        require_clean=require_clean
    )
```

- [ ] **Step 4: Run new tests, expect PASS**

Run: `pytest tests/test_phases.py::TestXRDProfileGuidedWarrenAverbachPhaseKwarg -v`
Expected: PASS.

- [ ] **Step 5: Run regression test, expect PASS**

Run: `pytest tests/test_backward_compat.py -v`
Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git add xrd_profile/profile.py tests/test_phases.py
git commit -m "Add phase= and reserved instrumental= kwargs on guided_warren_averbach

Same pattern as guided_williamson_hall: phase= computes ref_peaks
internally from the profile's 2-theta range; mutually exclusive with
the legacy ref_peaks= positional/keyword argument; instrumental=
raises NotImplementedError."
```

---

## Task 6: Scherrer K and shape factor extension

Adds the `SCHERRER_K_FOR_SHAPE` lookup table and the `shape=` kwarg.
Crucially, changes `K`'s default from `0.9` to `None` (sentinel) while
preserving v0.2.0 behavior when neither is passed.

**Files:**
- Modify: `xrd_profile/scherrer.py`
- Modify: `xrd_profile/profile.py` (the `XRDProfile.scherrer()` and
  `XRDProfile.modified_scherrer()` wrappers)
- Create: `tests/test_scherrer.py`

- [ ] **Step 1: Write failing tests**

Create `tests/test_scherrer.py`:

```python
"""Unit tests for Scherrer K and shape factor extension (v0.3.0)."""
import numpy as np
import pytest

from xrd_profile import scherrer, modified_scherrer


CU_KA = 1.5406
FWHM = np.array([0.30, 0.40, 0.55, 0.60])
TT = np.array([20.0, 30.0, 45.0, 60.0])


class TestSchererDefaults:
    def test_no_kwargs_preserves_v020_behavior(self):
        """scherrer(fwhm, tt, wl) without kwargs uses K=0.9."""
        sizes = scherrer(FWHM, TT, CU_KA)
        sizes_explicit_09 = scherrer(FWHM, TT, CU_KA, K=0.9)
        np.testing.assert_allclose(sizes, sizes_explicit_09)

    def test_default_K_is_0_9(self):
        # Hand-computed reference: D = 0.9 * lambda / (beta_rad * cos(theta_rad))
        beta_rad = np.radians(FWHM)
        theta_rad = np.radians(TT / 2)
        expected = 0.9 * CU_KA / (beta_rad * np.cos(theta_rad))
        actual = scherrer(FWHM, TT, CU_KA)
        np.testing.assert_allclose(actual, expected)


class TestShapeKwarg:
    def test_lookup_table_values(self):
        from xrd_profile import SCHERRER_K_FOR_SHAPE
        assert SCHERRER_K_FOR_SHAPE['spherical'] == 0.94
        assert SCHERRER_K_FOR_SHAPE['cubic'] == 0.83
        assert SCHERRER_K_FOR_SHAPE['cylindrical'] == 1.84
        assert SCHERRER_K_FOR_SHAPE['platey'] == 1.0

    def test_shape_spherical_uses_K_094(self):
        sizes_shape = scherrer(FWHM, TT, CU_KA, shape='spherical')
        sizes_explicit = scherrer(FWHM, TT, CU_KA, K=0.94)
        np.testing.assert_allclose(sizes_shape, sizes_explicit)

    def test_shape_cylindrical_uses_K_184(self):
        sizes_shape = scherrer(FWHM, TT, CU_KA, shape='cylindrical')
        sizes_explicit = scherrer(FWHM, TT, CU_KA, K=1.84)
        np.testing.assert_allclose(sizes_shape, sizes_explicit)


class TestKWinsOverShape:
    def test_explicit_K_overrides_shape(self):
        sizes_K = scherrer(FWHM, TT, CU_KA, K=1.0, shape='spherical')
        sizes_K_only = scherrer(FWHM, TT, CU_KA, K=1.0)
        np.testing.assert_allclose(sizes_K, sizes_K_only)


class TestModifiedScherrer:
    def test_default_K_preserves_v020(self):
        size_default = modified_scherrer(FWHM, TT, CU_KA)
        size_explicit = modified_scherrer(FWHM, TT, CU_KA, K=0.9)
        assert np.isclose(size_default, size_explicit)

    def test_shape_works(self):
        size = modified_scherrer(FWHM, TT, CU_KA, shape='cubic')
        size_explicit = modified_scherrer(FWHM, TT, CU_KA, K=0.83)
        assert np.isclose(size, size_explicit)
```

- [ ] **Step 2: Run, expect FAIL**

Run: `pytest tests/test_scherrer.py -v`
Expected: FAIL with `ImportError: cannot import name 'SCHERRER_K_FOR_SHAPE'`
or unexpected-kwarg errors.

- [ ] **Step 3: Modify xrd_profile/scherrer.py**

Replace the file content (read first, then Write or full Edit) with:

```python
"""
scherrer.py — Scherrer equation analysis for crystallite size estimation.

Provides standard per-peak Scherrer analysis and the modified Scherrer
method using log-linear regression for an average crystallite size.

References
----------
Scherrer, P. (1918). Bestimmung der Grosse und der inneren Struktur
    von Kolloidteilchen mittels Rontgenstrahlen. Nachrichten von der
    Gesellschaft der Wissenschaften zu Gottingen, Mathematisch-
    Physikalische Klasse, 98-100.
Langford, J. I. & Wilson, A. J. C. (1978). Scherrer after sixty years:
    a survey and some new results in the determination of crystallite
    size. Journal of Applied Crystallography, 11, 102-113.

Adapted from crystallite_size_calculator
(https://github.com/bafgreat/crystallite_size_calculator) by Dinga Wonanke.
Modified by M.R.M. Izawa.
"""

import numpy as np
from scipy.stats import linregress


SCHERRER_K_FOR_SHAPE = {
    'spherical':   0.94,
    'cubic':       0.83,
    'cylindrical': 1.84,
    'platey':      1.0,
}


def _resolve_K(K, shape):
    """Resolve effective K from K and shape kwargs.

    both None  -> 0.9 (v0.2.0 default)
    shape only -> SCHERRER_K_FOR_SHAPE[shape]
    K only     -> K
    both       -> K wins, shape silently ignored
    """
    if K is None and shape is None:
        return 0.9
    if K is None and shape is not None:
        try:
            return SCHERRER_K_FOR_SHAPE[shape]
        except KeyError:
            raise ValueError(
                f'unknown shape {shape!r}; choose from '
                f'{list(SCHERRER_K_FOR_SHAPE)}')
    return K


def scherrer(fwhm_deg, two_theta_positions, wavelength,
             K=None, shape=None):
    """
    Scherrer equation: D = K * lambda / (beta * cos(theta))

    Parameters
    ----------
    fwhm_deg : array-like
        FWHM values in degrees.
    two_theta_positions : array-like
        2-theta peak positions in degrees.
    wavelength : float
        X-ray wavelength in angstroms.
    K : float or None
        Scherrer constant. None = use v0.2.0 default 0.9 unless shape=
        is provided, in which case K is looked up from
        SCHERRER_K_FOR_SHAPE.
    shape : {'spherical', 'cubic', 'cylindrical', 'platey'} or None
        Crystallite shape. Used to look up K when K is None. If both
        are provided, K wins.

    Returns
    -------
    crystallite_sizes : np.ndarray
        Crystallite sizes in angstroms for each peak.
    """
    K_eff = _resolve_K(K, shape)
    fwhm_rad = np.radians(np.asarray(fwhm_deg))
    theta_rad = np.radians(np.asarray(two_theta_positions) / 2)
    return (K_eff * wavelength) / (fwhm_rad * np.cos(theta_rad))


def modified_scherrer(fwhm_deg, two_theta_positions, wavelength,
                       K=None, shape=None):
    """
    Modified Scherrer equation using log-linear regression.

    K, shape resolution: same as scherrer().
    """
    K_eff = _resolve_K(K, shape)
    fwhm_rad = np.radians(np.asarray(fwhm_deg))
    ln_beta = np.log(fwhm_rad)
    ln_1_cos = np.log(1 / np.cos(np.radians(np.asarray(two_theta_positions) / 2)))
    slope, intercept, _, _, _ = linregress(ln_1_cos, ln_beta)
    return (K_eff * wavelength) / np.exp(intercept)
```

- [ ] **Step 4: Update XRDProfile.scherrer and modified_scherrer wrappers**

Edit `xrd_profile/profile.py`. Locate `XRDProfile.scherrer` (currently
around line 150) and `XRDProfile.modified_scherrer` (around line 168).
Pass through the new K/shape kwargs:

```python
def scherrer(self, K=None, shape=None, height_threshold=0.05):
    """Run Scherrer analysis on all detected peaks.

    K and shape: see xrd_profile.scherrer.scherrer for resolution rules.
    """
    fwhm, positions = estimate_fwhm_simple(
        self.two_theta, self.intensity, height_threshold
    )
    if len(fwhm) == 0:
        return {'sizes': np.array([]), 'peak_positions': np.array([]),
                'd_spacings': np.array([]), 'fwhm': np.array([]),
                'mean_size': np.nan, 'median_size': np.nan}

    sizes = scherrer(fwhm, positions, self.wavelength, K=K, shape=shape)
    d_sp = two_theta_to_d(positions, self.wavelength)
    return {
        'sizes': sizes, 'peak_positions': positions,
        'd_spacings': d_sp, 'fwhm': fwhm,
        'mean_size': np.mean(sizes), 'median_size': np.median(sizes),
    }

def modified_scherrer(self, K=None, shape=None, height_threshold=0.05):
    """Run modified Scherrer equation. Returns average size (angstroms).

    K and shape: see xrd_profile.scherrer.modified_scherrer.
    """
    fwhm, positions = estimate_fwhm_simple(
        self.two_theta, self.intensity, height_threshold
    )
    if len(fwhm) < 2:
        return np.nan
    return modified_scherrer(fwhm, positions, self.wavelength,
                              K=K, shape=shape)
```

- [ ] **Step 5: Update __init__.py to export SCHERRER_K_FOR_SHAPE**

Edit `xrd_profile/__init__.py`. Find the import line for scherrer:
```python
from .scherrer import scherrer, modified_scherrer
```
Replace with:
```python
from .scherrer import scherrer, modified_scherrer, SCHERRER_K_FOR_SHAPE
```

Find the `__all__` list and add `'SCHERRER_K_FOR_SHAPE'`.

- [ ] **Step 6: Run tests, expect PASS**

Run: `pytest tests/test_scherrer.py -v`
Expected: all PASS.

- [ ] **Step 7: Run regression, expect PASS**

Run: `pytest tests/test_backward_compat.py -v`
Expected: all PASS. Critical: this confirms `scherrer(fwhm, tt, wl)`
without kwargs still produces v0.2.0-identical output.

- [ ] **Step 8: Commit**

```bash
git add xrd_profile/scherrer.py xrd_profile/profile.py xrd_profile/__init__.py tests/test_scherrer.py
git commit -m "Add Scherrer K and shape factor extension

K default changes from 0.9 to None sentinel; both-None resolves to
K=0.9 preserving v0.2.0 behavior. shape kwarg looks up K from
SCHERRER_K_FOR_SHAPE = {spherical: 0.94, cubic: 0.83, cylindrical:
1.84, platey: 1.0}. Explicit K wins over shape when both are passed.
SCHERRER_K_FOR_SHAPE is exported from xrd_profile."
```

---

## Task 7: XRDProfile.run_all() helper

**Files:**
- Modify: `xrd_profile/profile.py`
- Create: `tests/test_run_all.py`

- [ ] **Step 1: Write failing tests**

Create `tests/test_run_all.py`:

```python
"""Unit tests for XRDProfile.run_all dispatcher (v0.3.0)."""
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
    def test_instrumental_raises_not_implemented(self):
        profile, an = _profile_and_anorthite()
        with pytest.raises(NotImplementedError, match='Phase 2'):
            profile.run_all(methods=['wh'], phases=[an],
                             instrumental='anything')
```

- [ ] **Step 2: Run, expect FAIL**

Run: `pytest tests/test_run_all.py -v`
Expected: FAIL with `AttributeError: 'XRDProfile' object has no attribute 'run_all'`.

- [ ] **Step 3: Implement run_all in XRDProfile**

Edit `xrd_profile/profile.py`. Add this method to the `XRDProfile`
class, just before the `# --- Plotting ---` section (around line 312):

```python
    def run_all(self,
                 methods=None,
                 phases=None,
                 wh=None,
                 wa=None,
                 pdf=None,
                 scherrer=None,
                 instrumental=None):
        """
        Run a configurable bundle of analyses.

        Parameters
        ----------
        methods : list of {'wh', 'wa', 'pdf', 'scherrer'} or None
            Which analyses to run. None = all four.
        phases : list of Phase, single Phase, or None
            For guided W-H and W-A. None = unguided forms run.
        wh, wa, pdf, scherrer : dict or None
            Per-method kwargs.
            e.g. wh={'n_sigma': 3.0, 'tolerance_d': 0.02}.
        instrumental : reserved for Phase 2; raises NotImplementedError.

        Returns
        -------
        dict. With phases:
          {'wh': {phase.name: result, ...},
           'wa': {phase.name: result, ...},
           'pdf': result,
           'scherrer': result}
        Without phases:
          {'wh': result, 'wa': result, 'pdf': result, 'scherrer': result}
        """
        if instrumental is not None:
            raise NotImplementedError(
                'instrumental= is reserved for Phase 2 / v1.0; '
                'see xrd_profile roadmap')

        if methods is None:
            methods = ['wh', 'wa', 'pdf', 'scherrer']
        wh_kwargs = wh or {}
        wa_kwargs = wa or {}
        pdf_kwargs = pdf or {}
        scherrer_kwargs = scherrer or {}

        # Normalise phases to a list (or None)
        if phases is not None and not isinstance(phases, (list, tuple)):
            phases = [phases]

        results = {}

        if 'wh' in methods:
            if phases:
                results['wh'] = {
                    p.name: self.guided_williamson_hall(phase=p, **wh_kwargs)
                    for p in phases
                }
            else:
                from .williamson_hall import williamson_hall
                results['wh'] = williamson_hall(
                    self.two_theta, self.intensity, self.wavelength,
                    **wh_kwargs)

        if 'wa' in methods:
            if phases:
                results['wa'] = {
                    p.name: self.guided_warren_averbach(phase=p, **wa_kwargs)
                    for p in phases
                }
            else:
                results['wa'] = self.warren_averbach(**wa_kwargs)

        if 'pdf' in methods:
            r, G_r, Q_max = self.compute_pdf_sine(**pdf_kwargs)
            results['pdf'] = {'r': r, 'G_r': G_r, 'Q_max': Q_max}

        if 'scherrer' in methods:
            results['scherrer'] = self.scherrer(**scherrer_kwargs)

        return results
```

- [ ] **Step 4: Run, expect PASS**

Run: `pytest tests/test_run_all.py -v`
Expected: all PASS.

- [ ] **Step 5: Run regression, expect PASS**

Run: `pytest tests/test_backward_compat.py -v`
Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git add xrd_profile/profile.py tests/test_run_all.py
git commit -m "Add XRDProfile.run_all() convenience helper

Dispatcher that runs a configurable subset of {wh, wa, pdf, scherrer}
on a configurable list of phases. Per-method kwargs are passed as
nested dicts (wh={...}, wa={...}). With phases, results are keyed by
phase.name; without phases, results are the raw per-method dicts.
Instrumental= raises NotImplementedError."
```

---

# Phase C: Examples & Fixtures

## Task 8: Source example CIFs + SOURCES.md

**Files:**
- Create: `examples/cifs/Forsterite.cif`, `Anorthite.cif`, `Pigeonite.cif`,
  `Quartz.cif`, `Hematite.cif`
- Create: `examples/cifs/SOURCES.md`

- [ ] **Step 1: Create examples/cifs directory**

Run: `mkdir -p examples/cifs`

- [ ] **Step 2: Source the CIFs**

For each of the five minerals, obtain a CIF in this priority order:
(a) copy from the user's HOSERLab CIF directory if present at:
`C:\Users\Matthew Izawa\Desktop\Ye olde seagate\Big Bad Bucket of Backups\transfer\Mar 2016\The New Era - HoserLab\Rietveld\Structures\CIF files\`
or `phases/`. (b) Otherwise download from the Crystallography Open
Database (COD) at https://www.crystallography.net/cod/.

Recommended COD IDs (verify before download):
- Forsterite (Mg₂SiO₄): COD 9000093 (Bragg & Brown 1926, low-T form)
- Anorthite (CaAl₂Si₂O₈): COD 9000515 (Wenk et al. 1980) or any
  Wainwright & Starkey 1971 entry
- Pigeonite: use the user's `Pigeonite - Morimoto.cif` from HOSERLab
  if available (Morimoto pigeonite). If not, COD 9001629 (Clark, Appleman, Papike 1969).
- Quartz (SiO₂, low quartz): COD 9001143 or 9001144
- Hematite (Fe₂O₃): COD 9000139 (Blake et al. 1966)

Save each CIF to `examples/cifs/<MineralName>.cif` with capitalized
filename matching the manifest exactly: `Forsterite.cif`,
`Anorthite.cif`, `Pigeonite.cif`, `Quartz.cif`, `Hematite.cif`.

- [ ] **Step 3: Verify each CIF loads with Phase.from_cif**

Run a quick Python check:

```bash
python -c "from xrd_profile.phases import Phase; \
import os; \
[print(f'{f}: {Phase.from_cif(os.path.join(\"examples/cifs\", f))!r}') \
 for f in os.listdir('examples/cifs') if f.endswith('.cif')]"
```

Expected: each prints `<Phase 'name': formula>` with no errors.

If any CIF fails to load (older or non-standard format), pick a
different source CIF for that mineral and retry.

- [ ] **Step 4: Create examples/cifs/SOURCES.md**

```markdown
# CIF Sources

These CIF files are bundled with `xrd_profile` for use with the
example scripts and tests. Each file is attributed below per the
license terms of its source.

## Forsterite (Mg₂SiO₄)

- File: `Forsterite.cif`
- Source: Crystallography Open Database, COD ID [REPLACE_WITH_ID]
- Citation: [REPLACE_WITH_ORIGINAL_PAPER]
- License: CC0 (COD entries are in the public domain)

## Anorthite (CaAl₂Si₂O₈)

- File: `Anorthite.cif`
- Source: Crystallography Open Database, COD ID [REPLACE_WITH_ID]
- Citation: [REPLACE_WITH_ORIGINAL_PAPER]
- License: CC0

## Pigeonite (Mg-Fe-Ca pyroxene, Morimoto)

- File: `Pigeonite.cif`
- Source: HOSERLab (Cloutis) CIF library, originally from Morimoto et al.
- Citation: Morimoto, N. et al.
- Distribution: with permission from HOSERLab

## Quartz (SiO₂, low form)

- File: `Quartz.cif`
- Source: Crystallography Open Database, COD ID [REPLACE_WITH_ID]
- Citation: [REPLACE_WITH_ORIGINAL_PAPER]
- License: CC0

## Hematite (Fe₂O₃)

- File: `Hematite.cif`
- Source: Crystallography Open Database, COD ID [REPLACE_WITH_ID]
- Citation: [REPLACE_WITH_ORIGINAL_PAPER]
- License: CC0
```

Replace each `[REPLACE_WITH_*]` placeholder with the actual COD ID
and citation as the CIFs are downloaded. **Do not commit
SOURCES.md with placeholders intact** — the spec self-review
explicitly requires no placeholders in deliverables.

- [ ] **Step 5: Verify the test fixture works against an example CIF**

Append to `tests/test_phases.py`:

```python
class TestExampleCIFsLoad:
    """Smoke test that all bundled example CIFs are parseable."""
    @pytest.mark.parametrize('cif_name', [
        'Forsterite.cif', 'Anorthite.cif', 'Pigeonite.cif',
        'Quartz.cif', 'Hematite.cif',
    ])
    def test_example_cif_loads(self, cif_name):
        from pathlib import Path
        from xrd_profile.phases import Phase
        cif_path = Path(__file__).parent.parent / 'examples' / 'cifs' / cif_name
        if not cif_path.exists():
            pytest.skip(f'{cif_name} not present in examples/cifs/')
        phase = Phase.from_cif(cif_path)
        assert phase.name == cif_name.replace('.cif', '')
        peaks = phase.get_ref_peaks(1.5406)
        assert len(peaks) > 0
```

Run: `pytest tests/test_phases.py::TestExampleCIFsLoad -v`
Expected: 5 tests PASS.

- [ ] **Step 6: Commit**

```bash
git add examples/cifs/ tests/test_phases.py
git commit -m "Bundle five example CIFs with provenance documentation

Adds Forsterite, Anorthite, Pigeonite, Quartz, Hematite to
examples/cifs/ with full SOURCES.md attribution. A parameterised
test confirms every bundled CIF loads via Phase.from_cif and yields
non-empty reference peaks."
```

---

## Task 9: Bundle data fixture in examples/data/

**Files:**
- Create: `examples/data/tirhert_subset.xy` (copy of tests/fixtures/)
- Create: `examples/data/README.md`

- [ ] **Step 1: Create examples/data directory**

Run: `mkdir -p examples/data`

- [ ] **Step 2: Copy the bundled XY pattern**

Run: `cp tests/fixtures/tirhert_subset.xy examples/data/tirhert_subset.xy`

(The tests fixture is the canonical copy; this duplicate makes the
example self-contained without imports from the tests directory.)

- [ ] **Step 3: Create examples/data/README.md**

```markdown
# Example data

## tirhert_subset.xy

Downsampled subset of the Tirhert eucrite synchrotron diffraction
pattern (Diamond Light Source beamline I11, 15 keV,
λ = 0.826517 Å). The original pattern was collected at 0.001° steps
across 10–148° 2θ; this subset is trimmed to 10–80° and downsampled
to ~1000 points for fast example execution.

The full ungrouped pattern is associated with DLS beamtime ee17803-1
(Izawa & Jephcoat, 2018). Contact the authors for access to the
full dataset.

## Format

ASCII text, two whitespace-separated columns: 2θ (degrees) and
intensity. Header line begins with `#`.
```

- [ ] **Step 4: Commit**

```bash
git add examples/data/
git commit -m "Bundle tirhert_subset.xy in examples/data/ with README

Duplicates the tests/fixtures/ pattern into examples/data/ so the
example scripts are self-contained. README documents provenance and
points to the full dataset for users who want unbinned data."
```

---

## Task 10: Copy v0.2.0 examples to legacy/

**Files:**
- Create: `examples/legacy/lab_lunar_meteorite.py`,
  `synchrotron_low_shock.py`, `synchrotron_high_shock.py`

- [ ] **Step 1: Create examples/legacy directory**

Run: `mkdir -p examples/legacy`

- [ ] **Step 2: Copy the three example scripts verbatim**

```bash
cp examples/lab_lunar_meteorite.py examples/legacy/lab_lunar_meteorite.py
cp examples/synchrotron_low_shock.py examples/legacy/synchrotron_low_shock.py
cp examples/synchrotron_high_shock.py examples/legacy/synchrotron_high_shock.py
```

- [ ] **Step 3: Add a header comment to each legacy file**

For each of the three files in `examples/legacy/`, prepend a header:

```python
"""
[v0.2.0 LEGACY] This file is a verbatim copy of the v0.2.0 example
script. It demonstrates the pre-v0.3.0 pattern of constructing a
pymatgen.Structure inline and writing a build_ref helper. The top-level
example with the same filename in examples/ has been updated to use the
v0.3.0 Phase API; this file is preserved as a historical reference.

[Original docstring follows]
"""

# ... rest of original file unchanged
```

Use the Edit tool with `old_string` = the original first line of the
docstring (e.g., `"""\nExample: Lab XRD analysis of a lunar meteorite.\n`)
and `new_string` = the legacy header followed by the original.

- [ ] **Step 4: Verify the legacy files still run (smoke check, where data paths exist)**

Skipped: legacy files reference local data paths that may not be
available on every machine. The presence of the file is sufficient.

- [ ] **Step 5: Commit**

```bash
git add examples/legacy/
git commit -m "Preserve v0.2.0 example scripts verbatim in examples/legacy/

Three v0.2.0 examples copied unchanged so the pre-v0.3.0 pattern
remains available as historical reference. Each gets a header
comment indicating its legacy status. Top-level files with the same
names will be rewritten to use the new Phase API in the next task."
```

---

## Task 11: Rewrite existing examples to use Phase

**Files:**
- Modify: `examples/lab_lunar_meteorite.py`
- Modify: `examples/synchrotron_low_shock.py`
- Modify: `examples/synchrotron_high_shock.py`

For each file, the change is:
- Replace inline `Structure(lat, species, coords)` + duplicated
  `build_ref` helper with `Phase.from_lattice_params(...)` or
  `Phase.from_cif(...)`.
- Replace `profile.guided_williamson_hall(ref_d, ...)` with
  `profile.guided_williamson_hall(phase=phase_obj, ...)`.
- Replace `profile.guided_warren_averbach(ref_peaks, ...)` with
  `profile.guided_warren_averbach(phase=phase_obj, ...)`.
- Numerical results must match the legacy version (where the local
  data file is accessible).

- [ ] **Step 1: Rewrite synchrotron_low_shock.py (the most complex)**

This script currently does the inline anorthite + CIF-loaded pigeonite
two-phase analysis. The rewrite removes the `build_ref` helper.

Read the file first to copy its data path constants exactly. Then
write the new version:

```python
"""
Example: Synchrotron XRD analysis of a weakly shocked eucrite (Tirhert).

Analyses Tirhert (unbrecciated eucrite, fresh fall, low shock) from
Diamond Light Source beamline I11 (15 keV, lambda = 0.826517 angstroms,
10-148 deg 2-theta, 0.001 deg steps).

Demonstrates the v0.3.0 Phase API: anorthite via from_lattice_params
(refined coordinates inline), pigeonite via from_cif. The legacy v0.2.0
version of this script (with inline build_ref helper) is preserved at
examples/legacy/synchrotron_low_shock.py.

Data: DLS beamtime ee17803-1, Izawa & Jephcoat (2018).
"""

import numpy as np
from xrd_profile import XRDProfile
from xrd_profile.phases import Phase

LAMBDA = 0.826517
DATA_PATH = (r'C:\Users\Matthew Izawa\Desktop\111 Backup 20220530'
              r'\transfer\IPM\2018\ee17803-1\processing'
              r'\Tirhert_summed_0001.xye')
PIGEONITE_CIF = (r'C:\Users\Matthew Izawa\Desktop\Ye olde seagate'
                  r'\Big Bad Bucket of Backups\transfer\Mar 2016'
                  r'\The New Era - HoserLab\Rietveld\Structures'
                  r'\CIF files\Pigeonite - Morimoto.cif')

# ── Load data ──
data = np.loadtxt(DATA_PATH)
tt_full, i_full = data[:, 0], data[:, 1]
mask = (tt_full >= 10) & (tt_full <= 148)
tt, intensity = tt_full[mask], i_full[mask]

Q_max = 4 * np.pi * np.sin(np.radians(tt.max() / 2)) / LAMBDA
print(f"Tirhert (low-shock eucrite)")
print(f"  {len(tt)} points, {tt.min():.0f}-{tt.max():.0f} deg")
print(f"  Q_max = {Q_max:.1f} /A, dr = {np.pi/Q_max:.3f} A")

profile = XRDProfile(tt, intensity, wavelength=LAMBDA,
                      sample_name='Tirhert')

# ── Build phases ──
anorthite = Phase.from_lattice_params(
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
pigeonite = Phase.from_cif(PIGEONITE_CIF, name='Pigeonite')
print(f"\nAnorthite reference: {len(anorthite.get_ref_peaks(LAMBDA))} reflections")
print(f"Pigeonite reference: {len(pigeonite.get_ref_peaks(LAMBDA))} reflections")

# ── Anorthite analysis ──
print("\n--- Anorthite (plagioclase) ---")
an_wh = profile.guided_williamson_hall(phase=anorthite, n_sigma=3.0,
                                         tolerance_d=0.02)
an_wa = profile.guided_warren_averbach(phase=anorthite, n_sigma=3.0,
                                         tolerance_d=0.02)
print(f"  W-H: {an_wh['n_peaks']} peaks, "
      f"D = {an_wh['crystallite_size']:.0f} A, "
      f"strain = {an_wh['strain']:.5f}")
print(f"  W-A: {an_wa['n_families']} families, "
      f"D_median = {an_wa['median_crystallite_size']:.0f} A, "
      f"rms_strain = {an_wa['mean_rms_strain']:.4f}")

# ── Pigeonite analysis ──
print("\n--- Pigeonite (pyroxene) ---")
pig_wh = profile.guided_williamson_hall(phase=pigeonite, n_sigma=3.0,
                                          tolerance_d=0.02)
pig_wa = profile.guided_warren_averbach(phase=pigeonite, n_sigma=3.0,
                                          tolerance_d=0.02)
print(f"  W-H: {pig_wh['n_peaks']} peaks, "
      f"D = {pig_wh['crystallite_size']:.0f} A, "
      f"strain = {pig_wh['strain']:.5f}")
print(f"  W-A: {pig_wa['n_families']} families, "
      f"D_median = {pig_wa['median_crystallite_size']:.0f} A, "
      f"rms_strain = {pig_wa['mean_rms_strain']:.4f}")

# Phase comparison + plotting unchanged from legacy/synchrotron_low_shock.py.
# See that file for full plotting code.
```

Use the Write tool to overwrite `examples/synchrotron_low_shock.py`
with the above. Be sure to copy the unchanged plotting block from the
legacy file at the bottom (omitted here for brevity; the engineer
should preserve the full plot setup, just substituting `an_*` /
`pig_*` variables which already match).

- [ ] **Step 2: Rewrite synchrotron_high_shock.py and lab_lunar_meteorite.py**

Apply the same transformation pattern to the other two files: replace
`Structure(...)` + `build_ref(...)` with `Phase.from_*` + `phase=`
calls. Preserve all data paths, plot configuration, and printed
messages exactly. Read each legacy file from `examples/legacy/` to
copy structure-loading details verbatim.

- [ ] **Step 3: Verify imports and basic syntax**

For each rewritten file:
```bash
python -c "import ast; ast.parse(open('examples/synchrotron_low_shock.py').read())"
python -c "import ast; ast.parse(open('examples/synchrotron_high_shock.py').read())"
python -c "import ast; ast.parse(open('examples/lab_lunar_meteorite.py').read())"
```

Expected: no output (clean parse).

- [ ] **Step 4: Run the regression test, expect PASS**

Run: `pytest tests/test_backward_compat.py -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add examples/
git commit -m "Rewrite v0.2.0 examples to use Phase API

The three example scripts (lab_lunar_meteorite, synchrotron_low_shock,
synchrotron_high_shock) now use Phase.from_lattice_params and
Phase.from_cif instead of the inline Structure + build_ref helper.
Numerical results identical to the legacy versions (preserved at
examples/legacy/). DATA_PATH constants moved to top-of-file."
```

---

## Task 12: Write multi_phase_olivine.py

**Files:**
- Create: `examples/multi_phase_olivine.py`

The canonical demo of the v0.3.0 API. Loads the bundled
`examples/data/tirhert_subset.xy` (so the script is self-contained
and runs anywhere), builds Forsterite + Pigeonite via `Phase.from_cif`,
runs `profile.run_all(...)`, prints results, plots a 4-panel figure.

- [ ] **Step 1: Create the example script**

```python
"""
Example: multi-phase analysis with the v0.3.0 Phase API.

Demonstrates:
  - Loading phases from CIF files (forsterite, pigeonite).
  - Bundled diffraction pattern (Tirhert subset, ~1000 points).
  - run_all() bundled call covering W-H, W-A, PDF, Scherrer.
  - Per-phase result extraction.

This script is fully self-contained: it uses the CIFs and data
shipped in examples/. To run:

    pip install xrd_profile[cif]
    python examples/multi_phase_olivine.py
"""
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

from xrd_profile import XRDProfile
from xrd_profile.phases import Phase


HERE = Path(__file__).parent
DATA = HERE / 'data' / 'tirhert_subset.xy'
CIFS = HERE / 'cifs'
LAMBDA = 0.826517  # I11 at 15 keV — matches the bundled fixture

# ── Load data ──
xy = np.loadtxt(DATA)
tt, intensity = xy[:, 0], xy[:, 1]
profile = XRDProfile(tt, intensity, wavelength=LAMBDA,
                      sample_name='Tirhert (subset)')

# ── Build phases ──
forsterite = Phase.from_cif(CIFS / 'Forsterite.cif', name='Forsterite')
pigeonite = Phase.from_cif(CIFS / 'Pigeonite.cif', name='Pigeonite')
print(f"Phases: {forsterite}, {pigeonite}")

# ── Bundled analysis ──
results = profile.run_all(
    methods=['wh', 'wa', 'pdf', 'scherrer'],
    phases=[forsterite, pigeonite],
    wh={'n_sigma': 3.0, 'tolerance_d': 0.03},
    wa={'n_sigma': 3.0, 'tolerance_d': 0.03},
    pdf={'cheby_order': 15, 'lorch': True},
)

# ── Print summary ──
for phase_name in ['Forsterite', 'Pigeonite']:
    if phase_name in results['wh']:
        wh = results['wh'][phase_name]
        wa = results['wa'][phase_name]
        print(f"\n--- {phase_name} ---")
        print(f"  W-H: D = {wh.get('crystallite_size', float('nan')):.0f} A, "
              f"strain = {wh.get('strain', float('nan')):.5f}")
        print(f"  W-A: D_median = "
              f"{wa.get('median_crystallite_size', float('nan')):.0f} A, "
              f"families = {wa.get('n_families', 0)}")

print(f"\nPDF Q_max = {results['pdf']['Q_max']:.2f} /A, "
      f"dr = {np.pi / results['pdf']['Q_max']:.3f} A")
print(f"Scherrer mean size = {results['scherrer']['mean_size']:.0f} A")

# ── 4-panel figure ──
fig, axes = plt.subplots(2, 2, figsize=(13, 9))

# (a) pattern in d-spacing
ax = axes[0, 0]
profile.plot_pattern(ax=ax, x_axis='d_spacing', linewidth=0.4, color='k')
ax.set_xlim(0.8, 8)
ax.invert_xaxis()
ax.set_title('(a) Tirhert subset')

# (b) W-H reciprocal-space
ax = axes[0, 1]
for name, color in [('Forsterite', '#2166ac'), ('Pigeonite', '#b2182b')]:
    if name in results['wh']:
        wh = results['wh'][name]
        if 'K' in wh and 'deltaK' in wh:
            ax.scatter(wh['K'], wh['deltaK'], s=16, alpha=0.7,
                       color=color, label=name)
ax.set_xlabel(r'K ($\AA^{-1}$)')
ax.set_ylabel(r'$\Delta$K ($\AA^{-1}$)')
ax.set_title('(b) Guided W-H per phase')
ax.legend()

# (c) W-A anisotropic sizes per phase
ax = axes[1, 0]
y_offset = 0
for name, color in [('Forsterite', '#2166ac'), ('Pigeonite', '#b2182b')]:
    if name not in results['wa']:
        continue
    wa = results['wa'][name]
    sizes = [(f"{f['base_hkl']}", f['crystallite_size'])
             for f in wa.get('families', [])
             if not np.isnan(f['crystallite_size']) and f['crystallite_size'] > 0]
    if not sizes:
        continue
    labels, vals = zip(*sizes[:8])
    y = np.arange(len(vals)) + y_offset
    ax.barh(y, vals, 0.4, color=color, alpha=0.7, label=name)
    y_offset += len(vals) + 1
ax.set_xlabel('Crystallite size ($\\AA$)')
ax.set_title('(c) W-A anisotropic sizes')
ax.legend(fontsize=9)
ax.invert_yaxis()

# (d) PDF
ax = axes[1, 1]
r, G = results['pdf']['r'], results['pdf']['G_r']
mask = (r > 1) & (r < 15)
ax.plot(r[mask], G[mask], 'k', linewidth=0.8)
ax.set_xlabel(r'r ($\AA$)')
ax.set_ylabel('G(r)')
ax.axhline(0, color='grey', linewidth=0.3, linestyle=':')
ax.set_title(f'(d) PDF (Q_max = {results["pdf"]["Q_max"]:.1f} /A)')

fig.suptitle('xrd_profile v0.3.0 multi-phase example',
              fontsize=13)
fig.tight_layout()
out_png = HERE / 'multi_phase_olivine_example.png'
fig.savefig(out_png, dpi=150)
print(f"\nFigure saved to {out_png}")
```

- [ ] **Step 2: Run the script end-to-end**

Run: `python examples/multi_phase_olivine.py`

Expected: prints summary, saves
`examples/multi_phase_olivine_example.png`. No errors. (Forsterite
will not have many matched peaks against a Tirhert pattern — that's
expected; the script is demonstrating the API surface, not a
scientifically meaningful olivine analysis.)

- [ ] **Step 3: Commit**

```bash
git add examples/multi_phase_olivine.py examples/multi_phase_olivine_example.png
git commit -m "Add multi_phase_olivine.py: canonical v0.3.0 demo

Self-contained example using the bundled CIFs and data fixture.
Demonstrates Phase.from_cif loading, run_all() with multi-method
multi-phase configuration, per-phase result extraction, and a
4-panel figure. Runs anywhere xrd_profile[cif] is installed."
```

---

# Phase D: Release Plumbing

## Task 13: Update __init__.py exports

**Files:**
- Modify: `xrd_profile/__init__.py`

- [ ] **Step 1: Read current __init__.py**

Already read in Task 6. Should currently export Phase via the next
step.

- [ ] **Step 2: Add Phase, build_reference_peaks to imports and __all__**

Edit `xrd_profile/__init__.py`. Add this import line:

```python
from .phases import Phase, build_reference_peaks
```

And add to `__all__`:

```python
__all__ = [
    'XRDProfile',
    'Phase', 'build_reference_peaks',
    'two_theta_to_d', 'd_to_two_theta', 'two_theta_to_Q',
    'two_theta_to_K', 'fwhm_to_deltaK',
    'estimate_noise', 'estimate_zero_offset',
    'estimate_fwhm_simple', 'estimate_fwhm_voigt',
    'voigt_profile', 'find_peaks_guided',
    'scherrer', 'modified_scherrer', 'SCHERRER_K_FOR_SHAPE',
    'check_cross_phase_overlap', 'score_peak_quality',
    'chebyshev_background', 'compute_pdf_sine',
    'measure_pdf_peaks', 'fit_first_pdf_peak',
]
```

(SCHERRER_K_FOR_SHAPE was added in Task 6; this confirms it's also in __all__.)

- [ ] **Step 3: Verify the imports work**

Run:
```bash
python -c "from xrd_profile import Phase, build_reference_peaks, SCHERRER_K_FOR_SHAPE; print('ok')"
```

Expected: prints `ok`.

- [ ] **Step 4: Run all tests**

Run: `pytest -v`
Expected: all tests PASS.

- [ ] **Step 5: Commit**

```bash
git add xrd_profile/__init__.py
git commit -m "Export Phase, build_reference_peaks from xrd_profile

Top-level imports now include the v0.3.0 Phase API alongside the
existing array-based functions."
```

---

## Task 14: Update pyproject.toml ([cif] extra)

**Files:**
- Modify: `pyproject.toml`
- Modify: `setup.py` (if it contains version info)

- [ ] **Step 1: Edit pyproject.toml**

Edit `pyproject.toml`:

```toml
[build-system]
requires = ["setuptools>=61.0"]
build-backend = "setuptools.build_meta"

[project]
name = "xrd_profile"
version = "0.3.0"
description = "XRD peak profile analysis toolkit: Williamson-Hall, Warren-Averbach, Scherrer, and pair distribution function methods"
readme = "README.md"
license = {text = "MIT"}
requires-python = ">=3.8"
authors = [
    {name = "Matthew R. M. Izawa", email = "matthew.izawa@gmail.com"},
]
dependencies = [
    "numpy",
    "scipy",
    "matplotlib",
]

[project.optional-dependencies]
cif = ["pymatgen>=2023.0"]
dev = ["pytest", "pymatgen>=2023.0"]

[project.urls]
Homepage = "https://github.com/matthewizawa/xrd_profile"
```

(`dev` extra now includes pymatgen so test runs work without an extra
install command.)

- [ ] **Step 2: Read setup.py and check if version is duplicated there**

Run: `cat setup.py`

If `setup.py` contains a version string like `version='0.2.0'`,
update it to `'0.3.0'`. If it just calls `setup()` and reads from
`pyproject.toml`, leave it alone.

- [ ] **Step 3: Verify install with the [cif] extra**

Run:
```bash
pip install -e ".[cif]"
```

Expected: installs xrd_profile + pymatgen + sub-deps.

- [ ] **Step 4: Run all tests after reinstall**

Run: `pytest -v`
Expected: all PASS.

- [ ] **Step 5: Commit**

```bash
git add pyproject.toml setup.py
git commit -m "Add [cif] optional extra for pymatgen, bump version to 0.3.0

[cif] = ['pymatgen>=2023.0']. dev = pytest + pymatgen for test runs.
Plain pip install xrd_profile retains a numpy/scipy/matplotlib-only
install footprint; pip install xrd_profile[cif] enables the Phase
API."
```

---

## Task 15: Update README

**Files:**
- Modify: `README.md`

- [ ] **Step 1: Update the README**

Read `README.md`. Replace the **Installation** and **Quick start**
sections with v0.3.0 content. The new structure:

```markdown
# xrd_profile

A Python toolkit for quantitative analysis of powder X-ray diffraction
peak profiles.

## Features

- **Bragg's law conversions**: 2-theta, d-spacing, Q, K = 1/d (wavelength-independent)
- **Phase abstraction**: load any crystal structure from a CIF file or
  inline lattice parameters; reuse across analyses
- **Reference-guided peak detection**: automated noise cutoff and zero-point offset correction
- **Williamson-Hall analysis**: conventional and reciprocal-space (DeltaK vs K) formulations
- **Warren-Averbach analysis**: Fourier decomposition with harmonic peak families, adaptive window, Tukey tapering, quality filtering
- **Scherrer equation**: standard per-peak and modified log-linear regression, with K and shape factor selection
- **Pair distribution function**: sine Fourier transform PDF with iterative Chebyshev background subtraction, Lorch modification, peak detection, Gaussian first-shell fitting
- **Bundled run_all() helper**: configurable subset of methods on configurable list of phases
- **XRDProfile class**: unified interface wrapping all methods with plotting utilities

## Installation

```bash
pip install xrd_profile[cif]
```

The `[cif]` extra installs `pymatgen`, required for `Phase.from_cif`
and `Phase.from_lattice_params`. For environments without pymatgen
(if your reference peaks come from Rietveld output or other sources):

```bash
pip install xrd_profile
```

The array-based public API (`guided_williamson_hall(ref_d, ...)`,
`guided_warren_averbach(ref_peaks, ...)`, etc.) works without pymatgen.

## Quick start

```python
import numpy as np
from xrd_profile import XRDProfile, Phase

# Load diffraction data
data = np.loadtxt('my_pattern.xy', comments='#')
two_theta, intensity = data[:, 0], data[:, 1]

profile = XRDProfile(two_theta, intensity, wavelength=1.5406,
                      sample_name='My Sample')

# Build a phase from a CIF
olivine = Phase.from_cif('Forsterite.cif', name='Olivine')

# Run a guided analysis
wh = profile.guided_williamson_hall(phase=olivine, n_sigma=3.0)
print(f"W-H crystallite size: {wh['crystallite_size']:.0f} A")

# Or run a bundled multi-method analysis
results = profile.run_all(
    methods=['wh', 'wa', 'pdf', 'scherrer'],
    phases=[olivine],
    wh={'n_sigma': 3.0, 'tolerance_d': 0.02},
    wa={'tolerance_d': 0.02},
)
print(results['wh']['Olivine']['crystallite_size'])
print(results['pdf']['Q_max'])
```

See `examples/multi_phase_olivine.py` for a complete walk-through.

## Phases

A `Phase` wraps a `pymatgen.Structure` and provides reference-peak
generation. Construct from a CIF or from inline lattice parameters
and atomic coordinates:

```python
from xrd_profile import Phase

# From a CIF
quartz = Phase.from_cif('Quartz.cif', name='Quartz')

# From inline lattice parameters (useful when you have refined values)
anorthite = Phase.from_lattice_params(
    8.18, 12.88, 7.11, 93.5, 116.1, 90.4,
    species=['Ca', 'Al', 'Al', 'Si', 'Si', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O'],
    coords=[...],
    name='Anorthite',
)

# Pass to any guided method
wh = profile.guided_williamson_hall(phase=quartz)
```

Bundled example CIFs are in `examples/cifs/` (forsterite, anorthite,
pigeonite, quartz, hematite).

## Lower-level array-based API

Users who already have d-spacings (e.g., from Rietveld output) can
bypass `Phase` and pass arrays directly:

```python
ref_d = np.array([3.20, 3.18, 3.65, 4.04, 6.41])  # from Rietveld
wh = profile.guided_williamson_hall(ref_d=ref_d, n_sigma=3.0)
```

This API is foundational, supported indefinitely, and does not require
pymatgen.

## Pair distribution function analysis

[unchanged section from v0.2.0]

## Dependencies

- numpy
- scipy
- matplotlib
- pymatgen (optional, for `Phase.from_cif`/`from_lattice_params`)

## Attribution

[unchanged section]

## License

MIT License. See LICENSE file.

## Citation

```
Izawa, M. R. M. (2026). xrd_profile: XRD peak profile analysis toolkit
(v0.3.0). https://github.com/matthewizawa/xrd_profile
```
```

Use the Edit tool with the existing README content as `old_string` and
the rewritten content as `new_string`. (Or Write the whole file if the
content is large enough.)

- [ ] **Step 2: Verify README renders correctly**

```bash
python -c "import markdown; print(len(markdown.markdown(open('README.md').read())))"
```
Expected: prints a positive integer (HTML length). Or use any markdown
preview to spot-check formatting.

- [ ] **Step 3: Commit**

```bash
git add README.md
git commit -m "Update README to v0.3.0: Phase API, run_all, [cif] install

Quickstart now leads with Phase.from_cif and run_all. Lower-level
array-based API is documented as a foundational layer for users with
d-spacings from Rietveld. Citation block bumped to v0.3.0."
```

---

## Task 16: Update CHANGELOG, version bump in __init__.py, tag v0.3.0

**Files:**
- Modify: `CHANGELOG.md`
- Modify: `xrd_profile/__init__.py` (version string)
- Git tag: `v0.3.0`

- [ ] **Step 1: Update CHANGELOG.md**

Replace the `[Unreleased]` section in `CHANGELOG.md` with:

```markdown
## [0.3.0] — 2026-05-XX

### Added
- `Phase` class with `from_cif()` and `from_lattice_params()` constructors.
  Reference peaks are generated via `phase.get_ref_peaks(wavelength)` and
  `phase.get_ref_d(wavelength)`. Reduces user-script duplication of the
  `build_ref` helper.
- `phase=` keyword argument on `XRDProfile.guided_williamson_hall()` and
  `guided_warren_averbach()`. Mutually exclusive with the existing
  `ref_d=` / `ref_peaks=` arguments.
- `instrumental=` keyword argument reserved on the same methods (raises
  `NotImplementedError` if used; Phase 2 wires it up).
- `XRDProfile.run_all(methods=[...], phases=[...], wh={...}, wa={...},
  pdf={...}, scherrer={...})` convenience helper.
- `shape` keyword argument on `scherrer()` and `modified_scherrer()`,
  with a `SCHERRER_K_FOR_SHAPE` lookup table covering spherical (0.94),
  cubic (0.83), cylindrical (1.84), and platey (1.0) crystallite shapes.
- Five example CIFs (forsterite, anorthite, pigeonite, quartz, hematite)
  in `examples/cifs/`, with provenance documented in `SOURCES.md`.
- Bundled diffraction pattern at `examples/data/tirhert_subset.xy`
  (~30 KB Tirhert eucrite subset, I11 synchrotron).
- New canonical demo `examples/multi_phase_olivine.py`.
- Numerical regression test (`tests/test_backward_compat.py`) against
  frozen v0.2.0 outputs, ensuring strict-additive policy is enforced.

### Changed
- `scherrer()` and `modified_scherrer()`: `K` default changed from `0.9`
  to `None` sentinel. When neither `K` nor `shape` is provided, behavior
  is identical to v0.2.0 (`K=0.9`).
- Existing example scripts (`lab_lunar_meteorite.py`,
  `synchrotron_low_shock.py`, `synchrotron_high_shock.py`) rewritten to
  use the `Phase` API. Verbatim v0.2.0 versions preserved at
  `examples/legacy/`.

### Optional dependencies
- New `[cif]` extra installs `pymatgen` for `Phase.from_cif` and
  `from_lattice_params`. Plain `pip install xrd_profile` retains a
  pymatgen-free install footprint.
```

Replace `2026-05-XX` with the actual commit date.

- [ ] **Step 2: Bump version string in __init__.py**

Edit `xrd_profile/__init__.py`. Locate:

```python
__version__ = '0.2.0'
```

Replace with:

```python
__version__ = '0.3.0'
```

- [ ] **Step 3: Run the full test suite once more**

Run: `pytest -v`

Expected: all tests PASS, including `test_backward_compat.py`.

- [ ] **Step 4: Commit version bump + CHANGELOG**

```bash
git add CHANGELOG.md xrd_profile/__init__.py
git commit -m "Release v0.3.0

CHANGELOG documents the Phase API addition, run_all helper, Scherrer
K/shape extension, example CIFs, bundled data fixture, and numerical
regression test. Version string bumped from 0.2.0 to 0.3.0."
```

- [ ] **Step 5: Tag v0.3.0**

Run:
```bash
git tag -a v0.3.0 HEAD -m "Release v0.3.0: Phase API, run_all, Scherrer K/shape"
git tag --list "v*"
```

Expected output: `v0.2.0` and `v0.3.0`.

- [ ] **Step 6: Acceptance check**

Run the full acceptance criteria from the spec (Section 10):

```bash
# 1. Phase.from_cif and from_lattice_params work
python -c "from xrd_profile.phases import Phase; \
[Phase.from_cif(f'examples/cifs/{n}.cif') for n in \
 ['Forsterite', 'Anorthite', 'Pigeonite', 'Quartz', 'Hematite']]; \
print('1. CIF loading: OK')"

# 2-4 + 7-8: covered by pytest
pytest -v

# 5. Example CIFs all load
ls examples/cifs/

# 6. multi_phase_olivine.py runs end-to-end
python examples/multi_phase_olivine.py

# 9. README has Phase pattern
grep -q 'Phase.from_cif' README.md && echo "9. README: OK"

# 10. CHANGELOG documents v0.3.0
grep -q '0.3.0' CHANGELOG.md && echo "10. CHANGELOG: OK"

# 11-12. Tags exist
git tag --list "v0.2.0" "v0.3.0"
```

All should succeed with no errors and the expected confirmations.

---

## Self-review checklist

After completing all tasks, run through:

1. `pytest -v` — all tests pass, no warnings about deprecated APIs.
2. `git log --oneline v0.2.0..v0.3.0` — clean linear history of Phase 1 commits.
3. `git diff v0.2.0 -- xrd_profile/williamson_hall.py xrd_profile/warren_averbach.py xrd_profile/pdf.py xrd_profile/conversions.py xrd_profile/noise.py xrd_profile/peak_detection.py` — empty diffs (these files MUST be unchanged per the spec).
4. `python -c "import xrd_profile; print(xrd_profile.__version__)"` — prints `0.3.0`.
5. `python -c "from xrd_profile import Phase, SCHERRER_K_FOR_SHAPE, build_reference_peaks; print('ok')"` — clean import.
6. The legacy example files in `examples/legacy/` are byte-for-byte identical to the pre-Phase-1 versions of the corresponding files (modulo the legacy header).
