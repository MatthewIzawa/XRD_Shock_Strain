# v0.3.0 vs v0.4.0 Survey Comparison Harness — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a parameterised single-shot harness that re-runs the published 29-sample JAC survey through `xrd_profile` v0.4.0 — once with `instrumental=None` for a strict-additive regression check (Column B vs the saved CSV's Column A) and once with a per-instrument synthetic-LaB6 `InstrumentalStandard` for sensitivity (Column C) — producing a merged comparison CSV, a per-family size-distribution dump, a `run_log.md`, and 5 publication-track figures.

**Architecture:** Four Python files under `Llunr/comparison_v040/` (`config.py`, `synthetic_standards.py`, `run_comparison.py`, `plots.py`) plus `README.md` and one pytest test. `run_comparison.py` reads `Paper1_JAC/survey_results_29samples.csv` (Column A) and the source XY patterns, runs `XRDProfile.guided_warren_averbach` twice per sample (Column B with `instrumental=None`, Column C with `instrumental=<synthetic LaB6 InstrumentalStandard>` per instrument), writes the merged CSV / size-distribution dump / `run_log.md`, then calls `plots.py` to render the 5 figures.

**Tech Stack:** Python 3.13.9 (Anaconda base env at `C:\Users\Matthew Izawa\anaconda3\python.exe`); numpy, scipy, matplotlib, pymatgen, pytest (already in the base env per `xrd_profile/CLAUDE.md`); `xrd_profile` v0.4.0 importable from `Llunr/xrd_profile/`.

**Repository note:** `Llunr/` is **not** a git repo; only `xrd_profile/` is. The spec at `xrd_profile/docs/superpowers/specs/2026-05-06-v030-v040-survey-comparison-design.md` and this plan are in the package's git tree, but the harness implementation files at `Llunr/comparison_v040/*.py` are project-local and outside git. "Commit" steps in this plan apply only to the spec and plan in `xrd_profile/`; implementation files are saved in place. A single end-of-plan commit (Task 9) lands spec + plan in the package repo.

**Test command (harness-internal pytest):**
```
"/c/Users/Matthew Izawa/anaconda3/python.exe" -m pytest "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/comparison_v040/tests/" -v
```

**Run command (full harness end-to-end):**
```
"/c/Users/Matthew Izawa/anaconda3/python.exe" "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/comparison_v040/run_comparison.py"
```

---

## File Structure

| File | Role |
|---|---|
| `Llunr/comparison_v040/README.md` | How to invoke, expected outputs, prerequisites |
| `Llunr/comparison_v040/config.py` | Pure data: `SAMPLES` (29 entries), `INSTRUMENTAL_CAGLIOTI`, paths, `PHASE_REFS_SPEC` |
| `Llunr/comparison_v040/synthetic_standards.py` | `build_synthetic_lab6(wavelength, U, V, W, output_path) → (xy_path, InstrumentalStandard)` |
| `Llunr/comparison_v040/run_comparison.py` | Orchestrator: startup, per-sample loop, output, run_log.md, dispatches plots.py |
| `Llunr/comparison_v040/plots.py` | The 5 figures — one function per figure |
| `Llunr/comparison_v040/tests/test_synthetic_standards.py` | One pytest test for FWHM consistency |
| `Llunr/comparison_v040/output/` | All artifacts: CSVs, figures, synthetic standards, run_log.md |

---

## Reference: spec sections used per task

| Task | Spec sections |
|---|---|
| 1. Setup directory + README | §4.1 |
| 2. synthetic_standards.py | §6.3 |
| 3. config.py | §5.1, §5.2, §5.3, §5.4 |
| 4. run_comparison.py — startup phase | §5.5, §6.1 |
| 5. run_comparison.py — per-sample loop | §6.2, §6.4, §6.5 |
| 6. run_comparison.py — full run + run_log.md | §7.1, §7.2, §8.2, §9.1, §9.3 |
| 7. plots.py — figs 01 and 02 (scatter) | §7.3 |
| 8. plots.py — figs 03, 04, 05 | §7.3 |
| 9. Commit spec + plan to xrd_profile git | — |

---

## Task 1: Set up directory structure and README

**Files:**
- Create: `Llunr/comparison_v040/README.md`
- Create directories: `Llunr/comparison_v040/output/`, `Llunr/comparison_v040/output/synthetic_standards/`, `Llunr/comparison_v040/tests/`

- [ ] **Step 1: Create directories**

```bash
mkdir -p "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/comparison_v040/output/synthetic_standards"
mkdir -p "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/comparison_v040/tests"
```

- [ ] **Step 2: Verify directory tree**

```bash
ls "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/comparison_v040/"
ls "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/comparison_v040/output/"
```

Expected: `output/`, `tests/` visible; `synthetic_standards/` inside `output/`.

- [ ] **Step 3: Write README.md**

File: `Llunr/comparison_v040/README.md`

````markdown
# v0.3.0 vs v0.4.0 survey comparison harness

Validation and sensitivity comparison of `xrd_profile` v0.4.0 against
the saved 29-sample survey table at
`Paper1_JAC/survey_results_29samples.csv`.

See:
- Spec: `xrd_profile/docs/superpowers/specs/2026-05-06-v030-v040-survey-comparison-design.md`
- Plan: `xrd_profile/docs/superpowers/plans/2026-05-06-v030-v040-survey-comparison.md`

## Run

```
"/c/Users/Matthew Izawa/anaconda3/python.exe" run_comparison.py
```

(from this directory; or pass the absolute path to `run_comparison.py`).

Expected runtime: 10–20 minutes for the full 29-sample loop.

## Outputs

All artifacts land in `output/`:

- `comparison.csv` — 29-row merged A | B | C table (Sections 7.1 of spec)
- `size_distributions.csv` — per-family lognormal/normal fits
- `synthetic_standards/lab6_<key>.xy` — generated synthetic LaB6 patterns
- `fig01_D_median_AvB.{png,svg}` through `fig05_size_distribution_panels.{png,svg}`
- `run_log.md` — exact Caglioti values, A-vs-B pass/fail summary, top-5 deltas

## Prerequisites

- Anaconda Python 3.13.9 with numpy, scipy, matplotlib, pymatgen, pytest
- `Llunr/xrd_profile/` v0.4.0 importable
- `Paper1_JAC/survey_results_29samples.csv` exists
- Source XY data accessible at the locations encoded in `config.py`
- Co Kα data path filled into `config.py` before first run (see Task 3)

## Tests

```
"/c/Users/Matthew Izawa/anaconda3/python.exe" -m pytest tests/ -v
```
````

- [ ] **Step 4: Save and verify**

```bash
ls "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/comparison_v040/README.md"
```

Expected: file exists.

---

## Task 2: synthetic_standards.py — parameterised LaB6 builder

**Files:**
- Create: `Llunr/comparison_v040/synthetic_standards.py`
- Create: `Llunr/comparison_v040/tests/test_synthetic_standards.py`

The package's existing helper at `xrd_profile/scripts/build_synthetic_standard.py` is hard-coded to Cu Kα 1.5406 Å and lab-scale Caglioti `(U=5e-3, V=-1e-3, W=5e-3)`. This task generalises that math to arbitrary wavelength + Caglioti and returns both the written XY path AND a constructed `InstrumentalStandard` ready to pass to W-A.

The W-A code path requires `InstrumentalStandard` (not bare `InstrumentalProfile`) — verified at `xrd_profile/xrd_profile/profile.py:208–214`. The constructor signature is `InstrumentalStandard(phase, two_theta, intensity, wavelength, name='')` (verified at `xrd_profile/xrd_profile/instrumental.py:334–340`).

- [ ] **Step 1: Write the failing pytest test**

File: `Llunr/comparison_v040/tests/test_synthetic_standards.py`

```python
"""Tests for build_synthetic_lab6."""
import sys
from pathlib import Path

import numpy as np
import pytest

# Make the comparison_v040 package importable.
HERE = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(HERE))

# Make xrd_profile importable.
PKG = Path(r'C:\Users\Matthew Izawa\Documents\Dan Applin\Llunr\xrd_profile')
sys.path.insert(0, str(PKG))

from synthetic_standards import build_synthetic_lab6


def test_build_synthetic_lab6_returns_instrumental_standard(tmp_path):
    """Generated LaB6 wraps as an InstrumentalStandard."""
    from xrd_profile.instrumental import InstrumentalStandard
    out_xy = tmp_path / 'lab6_test.xy'
    xy_path, std = build_synthetic_lab6(
        wavelength=1.5406, U=5e-3, V=-1e-3, W=5e-3,
        output_path=str(out_xy))
    assert Path(xy_path).resolve() == out_xy.resolve()
    assert out_xy.exists()
    assert isinstance(std, InstrumentalStandard)
    assert std.wavelength == pytest.approx(1.5406)


def test_build_synthetic_lab6_fwhm_matches_caglioti(tmp_path):
    """Strongest peak below 35° 2θ has FWHM matching Caglioti to within 5%."""
    U, V, W = 5e-3, -1e-3, 5e-3
    out_xy = tmp_path / 'lab6_test.xy'
    build_synthetic_lab6(wavelength=1.5406, U=U, V=V, W=W,
                         output_path=str(out_xy))
    data = np.loadtxt(out_xy)
    tt, intensity = data[:, 0], data[:, 1]
    # Find strongest peak below 35° 2θ.
    mask = tt < 35
    masked_intensity = np.where(mask, intensity, -np.inf)
    peak_idx = int(np.argmax(masked_intensity))
    peak_tt = tt[peak_idx]
    # Expected Caglioti FWHM at that 2θ.
    theta = np.deg2rad(peak_tt / 2)
    expected_fwhm_sq = U * np.tan(theta)**2 + V * np.tan(theta) + W
    expected_fwhm = np.sqrt(max(expected_fwhm_sq, 1e-8))
    # Measure FWHM by half-max width on the contiguous segment around the peak.
    half_max = intensity[peak_idx] / 2.0
    above = intensity > half_max
    left = peak_idx
    while left > 0 and above[left - 1]:
        left -= 1
    right = peak_idx
    while right < len(tt) - 1 and above[right + 1]:
        right += 1
    measured_fwhm = tt[right] - tt[left]
    rel_err = abs(measured_fwhm - expected_fwhm) / expected_fwhm
    assert rel_err < 0.05, (
        f'measured FWHM {measured_fwhm:.4f}° vs expected '
        f'{expected_fwhm:.4f}° (rel err {rel_err:.3f})')


def test_build_synthetic_lab6_synchrotron_wavelength(tmp_path):
    """At 0.8265 Å (Diamond I11), reflections are accessible and peaks
    are an order of magnitude sharper than lab-scale Caglioti."""
    out_xy = tmp_path / 'lab6_synch.xy'
    xy_path, std = build_synthetic_lab6(
        wavelength=0.826517, U=1e-4, V=1e-5, W=1e-4,
        output_path=str(out_xy))
    data = np.loadtxt(out_xy)
    tt, intensity = data[:, 0], data[:, 1]
    # Should have peaks below ~35° 2θ even at this short wavelength.
    above_threshold = (tt < 35) & (intensity > 0.1 * intensity.max())
    assert above_threshold.sum() > 5, (
        f'too few peaks below 35° at synchrotron wavelength: '
        f'{above_threshold.sum()}')
```

- [ ] **Step 2: Run test to verify it fails**

```bash
"/c/Users/Matthew Izawa/anaconda3/python.exe" -m pytest "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/comparison_v040/tests/test_synthetic_standards.py" -v
```

Expected: ImportError on `from synthetic_standards import build_synthetic_lab6` — the file doesn't exist yet.

- [ ] **Step 3: Implement synthetic_standards.py**

File: `Llunr/comparison_v040/synthetic_standards.py`

```python
"""
Parameterised synthetic LaB6 pattern builder.

Generalises xrd_profile/scripts/build_synthetic_standard.py (which is
hard-coded to Cu Kα + lab-scale Caglioti) to arbitrary wavelength and
Caglioti coefficients. Used by run_comparison.py to build per-instrument
InstrumentalStandard objects from literature Caglioti values, since
W-A Stokes deconvolution requires a measured (or synthetic) standard
pattern, not just U/V/W (see xrd_profile/profile.py:208-214).

LaB6 is cubic Pm-3m, a = 4.156825 Å (NIST SRM 660c). Reflections up
to (3,1,1) populate the synthetic pattern. Each peak is a pure
Gaussian centred at the Bragg 2θ with FWHM given by the Caglioti
polynomial.
"""
from pathlib import Path
import sys

import numpy as np

# Make xrd_profile importable.
_PKG = Path(r'C:\Users\Matthew Izawa\Documents\Dan Applin\Llunr\xrd_profile')
if str(_PKG) not in sys.path:
    sys.path.insert(0, str(_PKG))

A_LAB6 = 4.156825  # NIST SRM 660c lattice parameter, Å


def caglioti_fwhm(two_theta_deg, U, V, W):
    """Caglioti FWHM in degrees from polynomial coefficients."""
    theta = np.deg2rad(two_theta_deg / 2.0)
    fwhm_sq = U * np.tan(theta)**2 + V * np.tan(theta) + W
    return np.sqrt(np.maximum(fwhm_sq, 1e-8))


def _multiplicity(h, k, l):
    """Cubic point-group multiplicity for (h,k,l)."""
    distinct = len({abs(h), abs(k), abs(l)})
    nonzero = sum(1 for v in (h, k, l) if v != 0)
    if distinct == 1 and nonzero == 3:
        return 8
    if nonzero == 1:
        return 6
    if distinct == 1 and nonzero == 2:
        return 12
    if nonzero == 2:
        return 24
    if distinct == 2 and nonzero == 3:
        return 24
    return 48


def _lab6_reflections(max_hkl=4):
    """Allowed reflections for LaB6 (Pm-3m, all h,k,l permitted)."""
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
                mult = _multiplicity(h, k, l)
                f_proxy = 1.0 / (1.0 + 0.05 * (h*h + k*k + l*l))
                refs.append((h, k, l, mult, f_proxy))
    return refs


def _build_lab6_phase():
    """Build a minimal LaB6 Phase suitable for InstrumentalStandard."""
    from xrd_profile import Phase
    return Phase.from_lattice_params(
        a=A_LAB6, b=A_LAB6, c=A_LAB6,
        alpha=90.0, beta=90.0, gamma=90.0,
        species=['La', 'B', 'B', 'B', 'B', 'B', 'B'],
        coords=[
            [0.0, 0.0, 0.0],
            [0.5, 0.7, 0.5], [0.5, 0.3, 0.5],
            [0.7, 0.5, 0.5], [0.3, 0.5, 0.5],
            [0.5, 0.5, 0.7], [0.5, 0.5, 0.3],
        ],
        name='LaB6_synthetic',
    )


def build_synthetic_lab6(wavelength, U, V, W, output_path,
                          tt_min=None, tt_max=None, tt_step=0.02):
    """
    Generate a synthetic LaB6 pattern with Caglioti-Gaussian peaks
    and wrap it as an InstrumentalStandard.

    Parameters
    ----------
    wavelength : float
        X-ray wavelength in angstroms.
    U, V, W : float
        Caglioti polynomial coefficients in deg², deg², deg².
    output_path : str or Path
        Where to write the synthetic XY pattern.
    tt_min, tt_max : float, optional
        2θ range in degrees. Defaults: tt_min = 5°, tt_max = the highest
        accessible LaB6 reflection 2θ + 5°, capped at 148°.
    tt_step : float
        2θ sampling step in degrees.

    Returns
    -------
    output_path : str
        Absolute path to the written XY file.
    standard : InstrumentalStandard
        Wrapped via the InstrumentalStandard constructor; ready to pass
        to XRDProfile.guided_warren_averbach as instrumental=.
    """
    refs = _lab6_reflections(max_hkl=4)
    peak_centres = []
    for h, k, l, mult, f_proxy in refs:
        d = A_LAB6 / np.sqrt(h*h + k*k + l*l)
        sin_theta = wavelength / (2.0 * d)
        if sin_theta >= 1.0:
            continue
        two_theta = 2.0 * np.degrees(np.arcsin(sin_theta))
        peak_centres.append((two_theta, mult * f_proxy * f_proxy))

    if not peak_centres:
        raise ValueError(
            f'No LaB6 reflections accessible at wavelength {wavelength} Å '
            f'(all sin(θ) >= 1).')

    if tt_min is None:
        tt_min = 5.0
    if tt_max is None:
        max_peak = max(c[0] for c in peak_centres)
        tt_max = min(148.0, max_peak + 5.0)

    tt = np.arange(tt_min, tt_max + tt_step, tt_step)
    intensity = np.zeros_like(tt)

    for centre, amp in peak_centres:
        if centre < tt_min or centre > tt_max:
            continue
        fwhm = caglioti_fwhm(centre, U, V, W)
        sigma = fwhm / (2.0 * np.sqrt(2.0 * np.log(2.0)))
        intensity += amp * np.exp(-((tt - centre) ** 2) / (2.0 * sigma ** 2))

    intensity += 0.001 * intensity.max()  # small flat baseline

    output_path = str(Path(output_path).resolve())
    np.savetxt(output_path, np.column_stack([tt, intensity]),
               fmt='%.6f %.6f')

    from xrd_profile.instrumental import InstrumentalStandard
    phase = _build_lab6_phase()
    std = InstrumentalStandard(
        phase=phase, two_theta=tt, intensity=intensity,
        wavelength=wavelength, name='synthetic_lab6')
    return output_path, std
```

- [ ] **Step 4: Run test to verify it passes**

```bash
"/c/Users/Matthew Izawa/anaconda3/python.exe" -m pytest "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/comparison_v040/tests/test_synthetic_standards.py" -v
```

Expected: `3 passed`.

- [ ] **Step 5: Visual sanity check**

```bash
"/c/Users/Matthew Izawa/anaconda3/python.exe" -c "
import sys
sys.path.insert(0, r'C:\Users\Matthew Izawa\Documents\Dan Applin\Llunr\comparison_v040')
sys.path.insert(0, r'C:\Users\Matthew Izawa\Documents\Dan Applin\Llunr\xrd_profile')
from synthetic_standards import build_synthetic_lab6
import tempfile, os, numpy as np
with tempfile.TemporaryDirectory() as td:
    xy, std = build_synthetic_lab6(0.826517, 1e-4, 1e-5, 1e-4, os.path.join(td, 'i11.xy'))
    data = np.loadtxt(xy)
    n_peaks = (data[:,1] > 0.1*data[:,1].max()).sum()
    print(f'I11 synthetic: 2θ {data[0,0]:.1f} to {data[-1,0]:.1f}°, n_pts {len(data)}, n_peaks (intensity > 10% max) {n_peaks}')
    print(f'wavelength on standard: {std.wavelength}')
"
```

Expected: 2θ range covering ~5° to ~70° at synchrotron wavelength, several thousand sample points, several distinct peaks above the 10% threshold; `wavelength: 0.826517`.

- [ ] **Step 6: Save in place (no commit — outside git)**

The two files (`synthetic_standards.py` and `tests/test_synthetic_standards.py`) live under `Llunr/comparison_v040/` which is outside any git repo. Verify they exist:

```bash
ls "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/comparison_v040/synthetic_standards.py"
ls "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/comparison_v040/tests/test_synthetic_standards.py"
```

Expected: both files visible.

---

## Task 3: config.py — sample list, Caglioti, paths

**Files:**
- Create: `Llunr/comparison_v040/config.py`

The 25 non-Co samples are copied verbatim from `Llunr/Paper1_JAC/compile_survey.py:118–218`. The 4 Co Kα samples are new entries; the user supplies the data path and filenames before first run.

`INSTRUMENTAL_CAGLIOTI` initial values are TYPICAL placeholders matching the package's existing synthetic-standard fixture (Cu Kα), with Co Kα equal to Cu Kα at TYPICAL magnitude and I11 set to 0.1× lab Caglioti. Each entry has a `'flag'` field set to `'TYPICAL'` until literature lookup pins exact values per Spec §8.

- [ ] **Step 1: Write config.py with all 29 samples and PHASE_REFS_SPEC**

File: `Llunr/comparison_v040/config.py`

```python
"""
Pure-data configuration for the v0.3.0 vs v0.4.0 survey comparison.

Replace the Caglioti values in INSTRUMENTAL_CAGLIOTI with literature
values (per Spec §8) before publication-track interpretation. Initial
values are TYPICAL fallbacks (lab Bruker for Cu/Co; 0.1× lab for I11)
and are flagged as such in run_log.md.

When measured calibration patterns become available, replace the
relevant entry's {'U','V','W','flag','citation'} dict with
{'measured_xy_path': ..., 'cif_path': ..., 'flag': 'MEASURED',
'citation': ...} and update the synthetic-standard build branch in
run_comparison.py accordingly.
"""
import os

# ── Wavelengths ──
CU_KA = 1.54056
CO_KA = 1.78897
SYNCH_I11 = 0.826517

# ── Source data directories ──
LUNAR_XY_DIR = r'C:\Users\Matthew Izawa\Documents\Conferences\LPSC 2025\Lunar\XY files'
SYNCH_DIR = (r'C:\Users\Matthew Izawa\Desktop\111 Backup 20220530'
             r'\transfer\IPM\2018\ee17803-1\processing')
# USER INPUT: fill in the Winnipeg Co Kα directory before first run.
CO_DIR = r'<USER INPUT NEEDED — Winnipeg Co Kα data directory>'

# ── CIF directory used by build_ref (mirror compile_survey.py) ──
CIF_DIR = (r'C:\Users\Matthew Izawa\Desktop\Ye olde seagate'
           r'\Big Bad Bucket of Backups\transfer\Mar 2016'
           r'\The New Era - HoserLab\Rietveld\Structures\CIF files')

# ── Saved survey table (Column A) ──
SURVEY_CSV_29 = (r'C:\Users\Matthew Izawa\Documents\Dan Applin\Llunr'
                 r'\Paper1_JAC\survey_results_29samples.csv')

# ── Output directory ──
OUTPUT_DIR = (r'C:\Users\Matthew Izawa\Documents\Dan Applin\Llunr'
              r'\comparison_v040\output')

# ── Per-instrument Caglioti (TYPICAL placeholders; replace per Spec §8) ──
INSTRUMENTAL_CAGLIOTI = {
    'cu_misasa': {
        'U': 5.0e-3, 'V': -1.0e-3, 'W': 5.0e-3,
        'wavelength': CU_KA,
        'flag': 'TYPICAL',
        'citation': 'TYPICAL Bragg-Brentano Cu Kα; magnitude matches package synthetic LaB6 fixture (xrd_profile/scripts/build_synthetic_standard.py).',
    },
    'co_winnipeg': {
        'U': 5.0e-3, 'V': -1.0e-3, 'W': 5.0e-3,
        'wavelength': CO_KA,
        'flag': 'TYPICAL',
        'citation': 'TYPICAL Bragg-Brentano Co Kα; magnitude matches Cu Kα BB fallback.',
    },
    'i11': {
        'U': 5.0e-4, 'V': -1.0e-4, 'W': 5.0e-4,
        'wavelength': SYNCH_I11,
        'flag': 'TYPICAL',
        'citation': 'TYPICAL synchrotron I11; 0.1× lab BB Caglioti pending Thompson et al. 2009 lookup.',
    },
}

# ── Phase reference build specification ──
# Used by run_comparison.py to populate PHASE_REFS at startup.
# 'name'    : phase identifier referenced by SAMPLES entries
# 'source'  : either 'lattice' (build inline) or 'cif' (load from CIF_DIR)
# 'wavelengths' : list of wavelengths to build references at
# For the inline anorthite, lattice + species + coords mirror compile_survey.py:49-56.
PHASE_REFS_SPEC = [
    {
        'name': 'anorthite',
        'source': 'lattice',
        'lattice': (8.1809, 12.881, 7.1101, 93.465, 116.11, 90.369),
        'species': ['Ca', 'Al', 'Al', 'Si', 'Si',
                    'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O'],
        'coords': [
            [0.269, 0.988, 0.086],
            [0.507, 0.314, 0.621],
            [0.992, 0.815, 0.118],
            [0.505, 0.320, 0.110],
            [0.006, 0.816, 0.613],
            [0.491, 0.625, 0.487],
            [0.024, 0.124, 0.995],
            [0.073, 0.488, 0.635],
            [0.576, 0.990, 0.143],
            [0.298, 0.356, 0.612],
            [0.817, 0.855, 0.142],
            [0.517, 0.179, 0.610],
            [0.000, 0.680, 0.104],
        ],
        'wavelengths': [CU_KA, CO_KA, SYNCH_I11],
    },
    {'name': 'enstatite',  'source': 'cif',
     'cif_filename': 'Enstatite - Hugh-Jones.cif',
     'wavelengths': [SYNCH_I11]},
    {'name': 'forsterite', 'source': 'cif',
     'cif_filename': 'Fo100 - Smyth.cif',
     'wavelengths': [SYNCH_I11]},
    {'name': 'pigeonite',  'source': 'cif',
     'cif_filename': 'Pigeonite - Morimoto.cif',
     'wavelengths': [SYNCH_I11]},
]

# ── Sample list (29 entries; align 1:1 with rows of SURVEY_CSV_29) ──
SAMPLES = [
    # ───── Lunar (7) — Cu Kα Misasa ─────
    {'name': 'Gadamis 004', 'instrument_key': 'cu_misasa',
     'wavelength': CU_KA, 'phase_ref': 'anorthite',
     'file': os.path.join(LUNAR_XY_DIR, 'gadamis004.xy'),
     'tt_range': (10, 60)},
    {'name': 'NWA 11182a', 'instrument_key': 'cu_misasa',
     'wavelength': CU_KA, 'phase_ref': 'anorthite',
     'file': os.path.join(LUNAR_XY_DIR, 'nwa11182a.xy'),
     'tt_range': (10, 60)},
    {'name': 'NWA 11182b', 'instrument_key': 'cu_misasa',
     'wavelength': CU_KA, 'phase_ref': 'anorthite',
     'file': os.path.join(LUNAR_XY_DIR, 'nwa11182b.xy'),
     'tt_range': (10, 60)},
    {'name': 'NWA 11788', 'instrument_key': 'cu_misasa',
     'wavelength': CU_KA, 'phase_ref': 'anorthite',
     'file': os.path.join(LUNAR_XY_DIR, 'nwa11788.xy'),
     'tt_range': (10, 60)},
    {'name': 'NWA 13788', 'instrument_key': 'cu_misasa',
     'wavelength': CU_KA, 'phase_ref': 'anorthite',
     'file': os.path.join(LUNAR_XY_DIR, 'nwa13788.xy'),
     'tt_range': (10, 60)},
    {'name': 'NWA 12593', 'instrument_key': 'cu_misasa',
     'wavelength': CU_KA, 'phase_ref': 'anorthite',
     'file': os.path.join(LUNAR_XY_DIR, 'nwa12594.xy'),  # mislabelled file
     'tt_range': (10, 60)},
    {'name': 'NWA 10401', 'instrument_key': 'cu_misasa',
     'wavelength': CU_KA, 'phase_ref': 'anorthite',
     'file': os.path.join(LUNAR_XY_DIR, 'nwa10401.xy'),
     'tt_range': (10, 60)},

    # ───── HED synchrotron (12) — I11 ─────
    {'name': 'Tirhert', 'instrument_key': 'i11',
     'wavelength': SYNCH_I11, 'phase_ref': 'anorthite',
     'file': os.path.join(SYNCH_DIR, 'Tirhert_summed_0001.xye'),
     'tt_range': (10, 148)},
    {'name': 'Bereba', 'instrument_key': 'i11',
     'wavelength': SYNCH_I11, 'phase_ref': 'anorthite',
     'file': os.path.join(SYNCH_DIR, 'Bereba_summed_0001.xye'),
     'tt_range': (10, 148)},
    {'name': 'NWA 6477', 'instrument_key': 'i11',
     'wavelength': SYNCH_I11, 'phase_ref': 'anorthite',
     'file': os.path.join(SYNCH_DIR, 'NWA_6477_summed_0001.xye'),
     'tt_range': (10, 148)},
    {'name': 'NWA 1836', 'instrument_key': 'i11',
     'wavelength': SYNCH_I11, 'phase_ref': 'anorthite',
     'file': os.path.join(SYNCH_DIR, 'NWA_1836_summed_0001.xye'),
     'tt_range': (10, 148)},
    {'name': 'NWA 6711', 'instrument_key': 'i11',
     'wavelength': SYNCH_I11, 'phase_ref': 'anorthite',
     'file': os.path.join(SYNCH_DIR, 'NWA_6711_summed_0001.xye'),
     'tt_range': (10, 148)},
    {'name': 'Millbillillie', 'instrument_key': 'i11',
     'wavelength': SYNCH_I11, 'phase_ref': 'anorthite',
     'file': os.path.join(SYNCH_DIR, 'Millbillillie_summed_0001.xye'),
     'tt_range': (10, 148)},
    {'name': 'JaH 626', 'instrument_key': 'i11',
     'wavelength': SYNCH_I11, 'phase_ref': 'anorthite',
     'file': os.path.join(SYNCH_DIR, 'JaH_626_summed_0001.xye'),
     'tt_range': (10, 148)},
    {'name': 'NWA 1943', 'instrument_key': 'i11',
     'wavelength': SYNCH_I11, 'phase_ref': 'anorthite',
     'file': os.path.join(SYNCH_DIR, 'NWA_1943_summed_0001.xye'),
     'tt_range': (10, 148)},
    {'name': 'NWA 1942', 'instrument_key': 'i11',
     'wavelength': SYNCH_I11, 'phase_ref': 'anorthite',
     'file': os.path.join(SYNCH_DIR, 'NWA_1942_summed_0001.xye'),
     'tt_range': (10, 148)},
    {'name': 'Dhofar 485', 'instrument_key': 'i11',
     'wavelength': SYNCH_I11, 'phase_ref': 'anorthite',
     'file': os.path.join(SYNCH_DIR, 'Dho_485_summed_0001.xye'),
     'tt_range': (10, 148)},
    {'name': 'NWA 2968', 'instrument_key': 'i11',
     'wavelength': SYNCH_I11, 'phase_ref': 'enstatite',
     'file': os.path.join(SYNCH_DIR, 'NWA_2968_summed_0001.xye'),
     'tt_range': (10, 148)},

    # ───── Non-HED synchrotron (6) — I11 ─────
    {'name': 'NWA 869', 'instrument_key': 'i11',
     'wavelength': SYNCH_I11, 'phase_ref': 'enstatite',
     'file': os.path.join(SYNCH_DIR, 'NWA_869_summed_0001.xye'),
     'tt_range': (10, 148)},
    {'name': 'NWA 801', 'instrument_key': 'i11',
     'wavelength': SYNCH_I11, 'phase_ref': 'forsterite',
     'file': os.path.join(SYNCH_DIR, 'NWA_801_summed_0001.xye'),
     'tt_range': (10, 148)},
    {'name': 'NWA 4872', 'instrument_key': 'i11',
     'wavelength': SYNCH_I11, 'phase_ref': 'forsterite',
     'file': os.path.join(SYNCH_DIR, 'NWA_4872_summed_0001.xye'),
     'tt_range': (10, 148)},
    {'name': 'NWA 7042', 'instrument_key': 'i11',
     'wavelength': SYNCH_I11, 'phase_ref': 'pigeonite',
     'file': os.path.join(SYNCH_DIR, 'NWA_7042_summed_0001.xye'),
     'tt_range': (10, 148)},
    {'name': 'Tagish Lake', 'instrument_key': 'i11',
     'wavelength': SYNCH_I11, 'phase_ref': 'forsterite',
     'file': os.path.join(SYNCH_DIR, 'Tagish_Lake_summed_0001.xye'),
     'tt_range': (10, 148)},
    {'name': 'NWA 753', 'instrument_key': 'i11',
     'wavelength': SYNCH_I11, 'phase_ref': 'forsterite',
     'file': os.path.join(SYNCH_DIR, 'NWA_753_summed_0001.xye'),
     'tt_range': (10, 148)},

    # ───── Co Kα Winnipeg (4) — USER INPUT NEEDED on filenames ─────
    {'name': 'NWA 5751', 'instrument_key': 'co_winnipeg',
     'wavelength': CO_KA, 'phase_ref': 'anorthite',
     'file': os.path.join(CO_DIR, '<USER INPUT NEEDED filename>'),
     'tt_range': (15, 90)},
    {'name': 'NWA 6013', 'instrument_key': 'co_winnipeg',
     'wavelength': CO_KA, 'phase_ref': 'anorthite',
     'file': os.path.join(CO_DIR, '<USER INPUT NEEDED filename>'),
     'tt_range': (15, 90)},
    {'name': 'Talampaya', 'instrument_key': 'co_winnipeg',
     'wavelength': CO_KA, 'phase_ref': 'anorthite',
     'file': os.path.join(CO_DIR, '<USER INPUT NEEDED filename>'),
     'tt_range': (15, 90)},
    {'name': 'Tatahouine', 'instrument_key': 'co_winnipeg',
     'wavelength': CO_KA, 'phase_ref': 'anorthite',
     'file': os.path.join(CO_DIR, '<USER INPUT NEEDED filename>'),
     'tt_range': (15, 90)},
]
```

- [ ] **Step 2: Sanity-check config imports cleanly**

```bash
"/c/Users/Matthew Izawa/anaconda3/python.exe" -c "
import sys
sys.path.insert(0, r'C:\Users\Matthew Izawa\Documents\Dan Applin\Llunr\comparison_v040')
import config
print(f'SAMPLES: {len(config.SAMPLES)} entries')
print(f'INSTRUMENTAL_CAGLIOTI keys: {list(config.INSTRUMENTAL_CAGLIOTI.keys())}')
print(f'PHASE_REFS_SPEC: {len(config.PHASE_REFS_SPEC)} entries')
print(f'OUTPUT_DIR: {config.OUTPUT_DIR}')
"
```

Expected output:
```
SAMPLES: 29 entries
INSTRUMENTAL_CAGLIOTI keys: ['cu_misasa', 'co_winnipeg', 'i11']
PHASE_REFS_SPEC: 4 entries
OUTPUT_DIR: C:\Users\Matthew Izawa\Documents\Dan Applin\Llunr\comparison_v040\output
```

- [ ] **Step 3: Cross-check 29-entry alignment with the saved CSV**

```bash
"/c/Users/Matthew Izawa/anaconda3/python.exe" -c "
import csv, sys
sys.path.insert(0, r'C:\Users\Matthew Izawa\Documents\Dan Applin\Llunr\comparison_v040')
import config
with open(config.SURVEY_CSV_29, encoding='utf-8') as f:
    r = csv.DictReader(f)
    csv_names = {row['Sample'] for row in r}
sample_names = {s['name'] for s in config.SAMPLES}
missing_in_config = csv_names - sample_names
missing_in_csv = sample_names - csv_names
print(f'csv_names: {len(csv_names)}; sample_names: {len(sample_names)}')
print(f'missing in config: {missing_in_config}')
print(f'missing in csv: {missing_in_csv}')
assert not missing_in_config, 'Some saved-CSV rows have no SAMPLES entry'
assert not missing_in_csv,    'Some SAMPLES entries are not in the saved CSV'
print('OK: 29 samples align 1:1.')
"
```

Expected: `OK: 29 samples align 1:1.`

If the cross-check fails on Co Kα names, double-check the spelling in the saved CSV vs `SAMPLES`.

- [ ] **Step 4: Save in place**

```bash
ls "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/comparison_v040/config.py"
```

Expected: file exists.

---

## Task 4: run_comparison.py — startup phase (config validation, PHASE_REFS, INSTRUMENTAL_STANDARDS)

**Files:**
- Create: `Llunr/comparison_v040/run_comparison.py` (this task adds startup; later tasks extend it)

This task gets the script to the point where it has all the prerequisites loaded and validated, but does not yet run any per-sample W-A. The end-of-task smoke test runs `python run_comparison.py` and expects it to print "Startup OK" and exit cleanly (or fail with a clear error message if Co Kα paths are still placeholders).

- [ ] **Step 1: Implement run_comparison.py with startup phase only**

File: `Llunr/comparison_v040/run_comparison.py`

```python
"""
v0.3.0 vs v0.4.0 survey comparison harness.

See xrd_profile/docs/superpowers/specs/2026-05-06-v030-v040-survey-comparison-design.md
for the design.

Run:
    "/c/Users/Matthew Izawa/anaconda3/python.exe" run_comparison.py
"""
import csv
import os
import sys
from pathlib import Path
from datetime import datetime

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

# Make xrd_profile importable.
PKG = Path(r'C:\Users\Matthew Izawa\Documents\Dan Applin\Llunr\xrd_profile')
sys.path.insert(0, str(PKG))

import config
from synthetic_standards import build_synthetic_lab6


def build_phase_refs():
    """Build PHASE_REFS dict from config.PHASE_REFS_SPEC.

    Returns dict keyed by (phase_name, wavelength) → list of {d, two_theta,
    intensity, h, k, l} reference peak dicts."""
    from xrd_profile import Phase
    from pymatgen.core.structure import Structure
    refs = {}
    for spec in config.PHASE_REFS_SPEC:
        if spec['source'] == 'lattice':
            phase = Phase.from_lattice_params(
                a=spec['lattice'][0], b=spec['lattice'][1], c=spec['lattice'][2],
                alpha=spec['lattice'][3], beta=spec['lattice'][4],
                gamma=spec['lattice'][5],
                species=spec['species'], coords=spec['coords'],
                name=spec['name'])
        elif spec['source'] == 'cif':
            cif_path = os.path.join(config.CIF_DIR, spec['cif_filename'])
            if not os.path.exists(cif_path):
                raise FileNotFoundError(
                    f'CIF for phase "{spec["name"]}" not found at {cif_path}')
            phase = Phase.from_cif(cif_path, name=spec['name'])
        else:
            raise ValueError(f'Unknown source for phase {spec["name"]}: '
                             f'{spec["source"]}')
        for w in spec['wavelengths']:
            ref_peaks = phase.get_ref_peaks(
                wavelength=w, two_theta_range=(5, 148), min_intensity=3.0)
            refs[(spec['name'], w)] = ref_peaks
    return refs


def build_instrumental_standards():
    """Build per-instrument synthetic LaB6 standards.

    Returns dict keyed by instrument_key → InstrumentalStandard."""
    standards = {}
    synth_dir = Path(config.OUTPUT_DIR) / 'synthetic_standards'
    synth_dir.mkdir(parents=True, exist_ok=True)
    for key, cag in config.INSTRUMENTAL_CAGLIOTI.items():
        out_xy = synth_dir / f'lab6_{key}.xy'
        _, std = build_synthetic_lab6(
            wavelength=cag['wavelength'],
            U=cag['U'], V=cag['V'], W=cag['W'],
            output_path=str(out_xy))
        standards[key] = std
    return standards


def load_saved_column_a():
    """Load Column A from the saved 29-sample CSV into a name-keyed dict.

    Returns dict: sample_name → row_dict (string-valued entries from csv.DictReader).
    """
    rows = {}
    with open(config.SURVEY_CSV_29, encoding='utf-8') as f:
        for row in csv.DictReader(f):
            rows[row['Sample']] = row
    return rows


def validate_startup(saved_a, phase_refs, instrumental_standards):
    """Spec §5.5 startup assertions."""
    samples = config.SAMPLES
    assert len(samples) == 29, (
        f'Expected 29 SAMPLES; got {len(samples)}')
    sample_names = [s['name'] for s in samples]
    assert len(set(sample_names)) == 29, (
        f'SAMPLES has duplicate names: '
        f'{sorted(n for n in sample_names if sample_names.count(n) > 1)}')
    csv_names = set(saved_a.keys())
    missing_in_csv = set(sample_names) - csv_names
    if missing_in_csv:
        raise AssertionError(
            f'SAMPLES names absent from saved CSV: {sorted(missing_in_csv)}')
    missing_in_samples = csv_names - set(sample_names)
    if missing_in_samples:
        raise AssertionError(
            f'Saved-CSV sample names absent from SAMPLES: '
            f'{sorted(missing_in_samples)}')
    for s in samples:
        if not os.path.exists(s['file']):
            raise FileNotFoundError(
                f"Sample {s['name']!r}: XY file not found at {s['file']}")
        if s['instrument_key'] not in config.INSTRUMENTAL_CAGLIOTI:
            raise KeyError(
                f"Sample {s['name']!r}: instrument_key "
                f"{s['instrument_key']!r} has no Caglioti entry")
        if s['instrument_key'] not in instrumental_standards:
            raise KeyError(
                f"Sample {s['name']!r}: instrument_key "
                f"{s['instrument_key']!r} has no synthetic standard")
        if (s['phase_ref'], s['wavelength']) not in phase_refs:
            raise KeyError(
                f"Sample {s['name']!r}: ({s['phase_ref']!r}, "
                f"{s['wavelength']}) has no PHASE_REFS entry")
        if not phase_refs[(s['phase_ref'], s['wavelength'])]:
            raise ValueError(
                f"Sample {s['name']!r}: "
                f"({s['phase_ref']!r}, {s['wavelength']}) "
                f"resolves to an empty reference peak list")


def main():
    print('=' * 60)
    print('xrd_profile v0.3.0 vs v0.4.0 survey comparison harness')
    print('=' * 60)
    print(f'Run start:    {datetime.now().isoformat()}')
    import xrd_profile
    print(f'xrd_profile:  v{xrd_profile.__version__}')
    print(f'Output dir:   {config.OUTPUT_DIR}')
    Path(config.OUTPUT_DIR).mkdir(parents=True, exist_ok=True)

    print('\n[1/4] Loading saved 29-sample CSV (Column A)...')
    saved_a = load_saved_column_a()
    print(f'      {len(saved_a)} rows.')

    print('\n[2/4] Building PHASE_REFS...')
    phase_refs = build_phase_refs()
    print(f'      {len(phase_refs)} (phase, wavelength) entries.')

    print('\n[3/4] Building synthetic LaB6 InstrumentalStandards...')
    instrumental_standards = build_instrumental_standards()
    for k, std in instrumental_standards.items():
        flag = config.INSTRUMENTAL_CAGLIOTI[k]['flag']
        print(f'      {k:<14} λ={std.wavelength:.4f} Å  flag={flag}')

    print('\n[4/4] Validating startup assertions...')
    validate_startup(saved_a, phase_refs, instrumental_standards)
    print('      OK — all 29 SAMPLES aligned with saved CSV; all XY files '
          'present; all instruments and phase refs resolved.')

    print('\nStartup OK. (Per-sample loop and output writing are added in '
          'subsequent tasks.)')


if __name__ == '__main__':
    main()
```

- [ ] **Step 2: Run the startup phase as a smoke test**

```bash
"/c/Users/Matthew Izawa/anaconda3/python.exe" "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/comparison_v040/run_comparison.py"
```

If Co Kα paths are still `<USER INPUT NEEDED ...>` placeholders, expect a clear `FileNotFoundError` naming `NWA 5751` (the first Co Kα entry):

```
FileNotFoundError: Sample 'NWA 5751': XY file not found at C:\...\<USER INPUT NEEDED — Winnipeg Co Kα data directory>\<USER INPUT NEEDED filename>
```

That is the expected first-run failure. Either:
- (a) supply the Co Kα path + filenames in `config.py` and re-run, OR
- (b) (acceptable interim) temporarily comment out the four Co Kα entries in `SAMPLES` and adjust the `len(samples) == 29` assertion to `>= 25`, to validate the rest of the harness; restore before the final Task 6 run.

If the four Co Kα paths are filled in correctly (or the entries are temporarily disabled), expect:

```
Startup OK. (Per-sample loop and output writing are added in subsequent tasks.)
```

If a *different* error appears (PHASE_REFS build fails, synthetic standard import fails, etc.), the implementer fixes that before continuing.

- [ ] **Step 3: Save in place**

```bash
ls "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/comparison_v040/run_comparison.py"
```

Expected: file exists.

---

## Task 5: run_comparison.py — per-sample loop (Columns B and C, A-vs-B tolerance)

**Files:**
- Modify: `Llunr/comparison_v040/run_comparison.py`

Add the per-sample loop, A-vs-B tolerance check, and output-CSV writing logic. Defer plot rendering and run_log.md to Task 6.

- [ ] **Step 1: Add the tolerance helper, sample-result extraction, CSV writers, and per-sample loop**

Open `Llunr/comparison_v040/run_comparison.py` and add the following helpers ABOVE `def main()`:

```python
# ──────────────────────────────────────────────────────────
# Per-sample analysis
# ──────────────────────────────────────────────────────────

# Tolerance for A vs B equality, per Spec §6.4.
TOL_FLOAT_D     = 1.0e-6   # |Δ D_median| in Å for Cu Kα + I11
TOL_FLOAT_E     = 1.0e-9   # |Δ rms_strain| for Cu Kα + I11
TOL_ROUNDED_D   = 0.5      # |Δ D_median| in Å for Co Kα (saved CSV is integer-rounded)
TOL_ROUNDED_E   = 5.0e-5   # |Δ rms_strain| for Co Kα (saved CSV is 4-decimal)


def _is_blank(s):
    return s is None or s == '' or str(s).strip() == ''


def _to_float(s):
    if _is_blank(s):
        return float('nan')
    return float(s)


def check_a_vs_b(saved_row, wa_b, instrument_key):
    """Compare saved Column A row against v0.4.0 instrumental=None Column B.

    Returns dict with status ('PASS', 'FAIL', 'N/A') and per-field
    abs-deltas."""
    a_n = _to_float(saved_row['WA_families'])
    a_d = _to_float(saved_row['WA_D_median_A'])
    a_e = _to_float(saved_row['WA_rms_strain'])
    b_n = wa_b.get('n_families', float('nan')) if wa_b else float('nan')
    b_d = wa_b.get('median_crystallite_size', float('nan')) if wa_b else float('nan')
    b_e = wa_b.get('mean_rms_strain', float('nan')) if wa_b else float('nan')

    if np.isnan(a_d) and np.isnan(a_e) and a_n == 0:
        return {'status': 'N/A', 'd_delta': float('nan'),
                'e_delta': float('nan'), 'n_delta': 0}

    if instrument_key == 'co_winnipeg':
        tol_d, tol_e = TOL_ROUNDED_D, TOL_ROUNDED_E
    else:
        tol_d, tol_e = TOL_FLOAT_D, TOL_FLOAT_E

    d_delta = abs(a_d - b_d) if not (np.isnan(a_d) or np.isnan(b_d)) else float('nan')
    e_delta = abs(a_e - b_e) if not (np.isnan(a_e) or np.isnan(b_e)) else float('nan')
    n_delta = abs(int(a_n) - int(b_n)) if not (np.isnan(a_n) or np.isnan(b_n)) else 0

    if (np.isnan(d_delta) or np.isnan(e_delta) or
            d_delta > tol_d or e_delta > tol_e or n_delta != 0):
        status = 'FAIL'
    else:
        status = 'PASS'
    return {'status': status, 'd_delta': d_delta,
            'e_delta': e_delta, 'n_delta': n_delta}


def load_xy(file_path, tt_range):
    """Load a 2-column XY (or .xye) file and slice to tt_range."""
    data = np.loadtxt(file_path, comments='#')
    tt = data[:, 0]
    intensity = data[:, 1]
    mask = (tt >= tt_range[0]) & (tt <= tt_range[1])
    return tt[mask], intensity[mask]


def run_one_sample(s, ref_peaks, std):
    """Run W-A on one sample, twice: with and without instrumental.

    Returns (wa_B, wa_C). Either may be None if a load or fit error occurs."""
    from xrd_profile import XRDProfile
    try:
        tt, intensity = load_xy(s['file'], s['tt_range'])
    except Exception as exc:
        print(f'    [load error] {s["name"]}: {exc}')
        return None, None
    profile = XRDProfile(tt, intensity, s['wavelength'],
                          sample_name=s['name'])
    try:
        wa_B = profile.guided_warren_averbach(
            ref_peaks, n_sigma=3.0, tolerance_d=0.02)
    except Exception as exc:
        print(f'    [W-A B error] {s["name"]}: {exc}')
        wa_B = None
    try:
        wa_C = profile.guided_warren_averbach(
            ref_peaks, n_sigma=3.0, tolerance_d=0.02,
            instrumental=std)
    except Exception as exc:
        print(f'    [W-A C error] {s["name"]}: {exc}')
        wa_C = None
    return wa_B, wa_C


def comparison_row(s, saved_row, wa_B, wa_C, ab_check):
    """Build one row dict for comparison.csv."""
    def field(d, k, default=float('nan')):
        return d.get(k, default) if d else default

    return {
        'Sample': s['name'],
        'Group': saved_row['Group'],
        'Type': saved_row['Type'],
        'Shock': saved_row['Shock'],
        'Source': saved_row['Source'],
        'Q_max': _to_float(saved_row['Q_max']),
        'A_WA_families': _to_float(saved_row['WA_families']),
        'A_WA_D_median_A': _to_float(saved_row['WA_D_median_A']),
        'A_WA_rms_strain': _to_float(saved_row['WA_rms_strain']),
        'B_WA_families': field(wa_B, 'n_families'),
        'B_WA_D_median_A': field(wa_B, 'median_crystallite_size'),
        'B_WA_rms_strain': field(wa_B, 'mean_rms_strain'),
        'C_WA_families': field(wa_C, 'n_families'),
        'C_WA_D_median_A': field(wa_C, 'median_crystallite_size'),
        'C_WA_rms_strain': field(wa_C, 'mean_rms_strain'),
        'A_vs_B_status': ab_check['status'],
        'delta_BC_D_median_A': (field(wa_C, 'median_crystallite_size')
                                  - field(wa_B, 'median_crystallite_size')),
        'delta_BC_D_median_pct': (
            100.0 * (field(wa_C, 'median_crystallite_size')
                      - field(wa_B, 'median_crystallite_size'))
            / field(wa_B, 'median_crystallite_size')
            if field(wa_B, 'median_crystallite_size') and
               not np.isnan(field(wa_B, 'median_crystallite_size')) and
               field(wa_B, 'median_crystallite_size') != 0
            else float('nan')),
        'delta_BC_rms_strain': (field(wa_C, 'mean_rms_strain')
                                 - field(wa_B, 'mean_rms_strain')),
        'A_PDF_pk1_r_A': _to_float(saved_row.get('PDF_pk1_r_A', '')),
        'A_PDF_pk1_FWHM_A': _to_float(saved_row.get('PDF_pk1_FWHM_A', '')),
    }


def size_distribution_rows(s, wa_C):
    """Build rows for size_distributions.csv from wa_C['families']."""
    if not wa_C:
        return []
    rows = []
    for fam_idx, fam in enumerate(wa_C.get('families', [])):
        sd = fam.get('size_distribution', {})
        ln = sd.get('lognormal', {})
        nm = sd.get('normal', {})
        # Representative reflection: lowest-order family member.
        peaks = fam.get('peaks', [])
        if peaks:
            rep = min(peaks, key=lambda p: (
                p.get('h', 0)**2 + p.get('k', 0)**2 + p.get('l', 0)**2))
            h, k, l = rep.get('h', ''), rep.get('k', ''), rep.get('l', '')
        else:
            h = k = l = ''
        rows.append({
            'Sample': s['name'],
            'instrument_key': s['instrument_key'],
            'family_index': fam_idx,
            'h': h, 'k': k, 'l': l,
            'LN_D_median_A': ln.get('D_median', float('nan')),
            'LN_sigma':       ln.get('sigma', float('nan')),
            'LN_R2':          ln.get('R2', float('nan')),
            'N_D_mean_A':     nm.get('D_mean', float('nan')),
            'N_sigma':        nm.get('sigma', float('nan')),
            'N_R2':           nm.get('R2', float('nan')),
            'n_valid_L':      sd.get('n_valid_L', 0),
        })
    return rows


COMPARISON_FIELDS = [
    'Sample', 'Group', 'Type', 'Shock', 'Source', 'Q_max',
    'A_WA_families', 'A_WA_D_median_A', 'A_WA_rms_strain',
    'B_WA_families', 'B_WA_D_median_A', 'B_WA_rms_strain',
    'C_WA_families', 'C_WA_D_median_A', 'C_WA_rms_strain',
    'A_vs_B_status',
    'delta_BC_D_median_A', 'delta_BC_D_median_pct',
    'delta_BC_rms_strain',
    'A_PDF_pk1_r_A', 'A_PDF_pk1_FWHM_A',
]

SIZE_DIST_FIELDS = [
    'Sample', 'instrument_key', 'family_index', 'h', 'k', 'l',
    'LN_D_median_A', 'LN_sigma', 'LN_R2',
    'N_D_mean_A',    'N_sigma',  'N_R2',
    'n_valid_L',
]


def write_csv(path, fieldnames, rows):
    with open(path, 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow(r)
```

- [ ] **Step 2: Update main() to call the per-sample loop**

Replace the existing `main()` function in `run_comparison.py` with:

```python
def main():
    print('=' * 60)
    print('xrd_profile v0.3.0 vs v0.4.0 survey comparison harness')
    print('=' * 60)
    print(f'Run start:    {datetime.now().isoformat()}')
    import xrd_profile
    print(f'xrd_profile:  v{xrd_profile.__version__}')
    print(f'Output dir:   {config.OUTPUT_DIR}')
    Path(config.OUTPUT_DIR).mkdir(parents=True, exist_ok=True)

    print('\n[1/5] Loading saved 29-sample CSV (Column A)...')
    saved_a = load_saved_column_a()
    print(f'      {len(saved_a)} rows.')

    print('\n[2/5] Building PHASE_REFS...')
    phase_refs = build_phase_refs()
    print(f'      {len(phase_refs)} (phase, wavelength) entries.')

    print('\n[3/5] Building synthetic LaB6 InstrumentalStandards...')
    instrumental_standards = build_instrumental_standards()
    for k, std in instrumental_standards.items():
        flag = config.INSTRUMENTAL_CAGLIOTI[k]['flag']
        print(f'      {k:<14} λ={std.wavelength:.4f} Å  flag={flag}')

    print('\n[4/5] Validating startup assertions...')
    validate_startup(saved_a, phase_refs, instrumental_standards)
    print('      OK.')

    print('\n[5/5] Running per-sample loop (29 samples × 2 W-A)...')
    comparison_rows = []
    sd_rows = []
    for i, s in enumerate(config.SAMPLES, start=1):
        print(f'  [{i:>2}/29] {s["name"]:<16} '
              f'instrument={s["instrument_key"]:<12}', end=' ')
        ref_peaks = phase_refs[(s['phase_ref'], s['wavelength'])]
        std = instrumental_standards[s['instrument_key']]
        wa_B, wa_C = run_one_sample(s, ref_peaks, std)
        ab = check_a_vs_b(saved_a[s['name']], wa_B, s['instrument_key'])
        comparison_rows.append(
            comparison_row(s, saved_a[s['name']], wa_B, wa_C, ab))
        sd_rows.extend(size_distribution_rows(s, wa_C))
        n_b = wa_B.get('n_families', '?') if wa_B else 'err'
        n_c = wa_C.get('n_families', '?') if wa_C else 'err'
        print(f'B fams={n_b}  C fams={n_c}  A-vs-B={ab["status"]}')

    out_dir = Path(config.OUTPUT_DIR)
    write_csv(out_dir / 'comparison.csv',
              COMPARISON_FIELDS, comparison_rows)
    write_csv(out_dir / 'size_distributions.csv',
              SIZE_DIST_FIELDS, sd_rows)
    print(f'\nWrote {out_dir / "comparison.csv"}')
    print(f'Wrote {out_dir / "size_distributions.csv"}')
    print('\nLoop complete. (run_log.md and figures are added in '
          'subsequent tasks.)')


if __name__ == '__main__':
    main()
```

- [ ] **Step 3: Single-sample smoke test before full run**

For a fast iteration, temporarily edit `config.py` `SAMPLES` to keep only `Tirhert` (synchrotron, 37 families per the saved CSV — exercises the most W-A code paths). Run:

```bash
"/c/Users/Matthew Izawa/anaconda3/python.exe" "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/comparison_v040/run_comparison.py"
```

Expected output (line for Tirhert):
```
  [ 1/29] Tirhert          instrument=i11           B fams=37  C fams=<int 1-37>  A-vs-B=PASS
```

`A-vs-B=PASS` confirms strict-additive holds for that one sample.
`C fams` should be in the range 1–37; if it's 0, the synthetic LaB6 for I11 is incompatible with the sample's data range.

If the smoke test fails, iterate before restoring the full SAMPLES list.

- [ ] **Step 4: Restore full SAMPLES list**

Restore the 29-entry `SAMPLES` list in `config.py`. Verify by re-running the cross-check:

```bash
"/c/Users/Matthew Izawa/anaconda3/python.exe" -c "
import sys
sys.path.insert(0, r'C:\Users\Matthew Izawa\Documents\Dan Applin\Llunr\comparison_v040')
import config
print(f'SAMPLES: {len(config.SAMPLES)} entries')
"
```

Expected: `SAMPLES: 29 entries`.

---

## Task 6: run_comparison.py — full survey run + run_log.md

**Files:**
- Modify: `Llunr/comparison_v040/run_comparison.py` (add `write_run_log()` and call it from `main()`)

This task adds the human-readable `run_log.md` per Spec §8.2 and §9.3, and runs the harness end-to-end on all 29 samples.

- [ ] **Step 1: Add write_run_log()**

In `run_comparison.py`, ABOVE `def main()` (alongside the other helpers), add:

```python
import subprocess


def _git_sha_of_xrd_profile():
    try:
        out = subprocess.check_output(
            ['git', 'rev-parse', '--short', 'HEAD'],
            cwd=str(PKG), stderr=subprocess.DEVNULL)
        return out.decode().strip()
    except Exception:
        return '(git SHA unavailable)'


def write_run_log(start_time, comparison_rows):
    """Emit run_log.md per Spec §8.2."""
    import xrd_profile
    out_path = Path(config.OUTPUT_DIR) / 'run_log.md'

    # A-vs-B summary by subset.
    by_subset = {'cu_misasa': [], 'co_winnipeg': [], 'i11': []}
    for r in comparison_rows:
        # Determine instrument_key from the SAMPLES list.
        for s in config.SAMPLES:
            if s['name'] == r['Sample']:
                by_subset[s['instrument_key']].append(r)
                break

    def summarise(rows):
        n = len(rows)
        passed = sum(1 for r in rows if r['A_vs_B_status'] == 'PASS')
        failed = sum(1 for r in rows if r['A_vs_B_status'] == 'FAIL')
        na     = sum(1 for r in rows if r['A_vs_B_status'] == 'N/A')
        return n, passed, failed, na

    cu_n, cu_p, cu_f, cu_na = summarise(by_subset['cu_misasa'])
    co_n, co_p, co_f, co_na = summarise(by_subset['co_winnipeg'])
    i11_n, i11_p, i11_f, i11_na = summarise(by_subset['i11'])

    # Top-5 |delta_BC_D_median_pct| across all rows.
    def abs_pct(r):
        v = r['delta_BC_D_median_pct']
        return abs(v) if v == v and v != float('inf') else -1.0
    top5 = sorted(comparison_rows, key=abs_pct, reverse=True)[:5]

    lines = []
    lines.append('# v0.3.0 vs v0.4.0 survey comparison — run log\n')
    lines.append(f'- **Run timestamp:** {start_time.isoformat()}')
    lines.append(f'- **xrd_profile version:** {xrd_profile.__version__}')
    lines.append(f'- **xrd_profile git SHA:** {_git_sha_of_xrd_profile()}')
    lines.append(f'- **Sample count:** {len(comparison_rows)}\n')

    lines.append('## Caglioti per instrument\n')
    lines.append('| Instrument | λ (Å) | U (deg²) | V (deg²) | W (deg²) | Flag | Citation |')
    lines.append('|---|---|---|---|---|---|---|')
    for k, c in config.INSTRUMENTAL_CAGLIOTI.items():
        lines.append(
            f'| `{k}` | {c["wavelength"]:.4f} | {c["U"]:.2e} | {c["V"]:.2e} | '
            f'{c["W"]:.2e} | {c["flag"]} | {c["citation"]} |')
    lines.append('')

    lines.append('## A vs B summary\n')
    lines.append(f'- Cu Kα Misasa: {cu_p}/{cu_n} PASS, {cu_f} FAIL, '
                 f'{cu_na} N/A')
    lines.append(f'- Diamond I11:  {i11_p}/{i11_n} PASS, {i11_f} FAIL, '
                 f'{i11_na} N/A')
    lines.append(f'- Co Kα Winnipeg: {co_p}/{co_n} PASS, {co_f} FAIL, '
                 f'{co_na} N/A')
    lines.append('')
    lines.append(f'Tolerances applied (per Spec §6.4): '
                 f'Cu+I11 |Δ D_median|<{TOL_FLOAT_D:.0e} Å, '
                 f'|Δ strain|<{TOL_FLOAT_E:.0e}; '
                 f'Co Kα |Δ D_median|<{TOL_ROUNDED_D} Å, '
                 f'|Δ strain|<{TOL_ROUNDED_E:.0e}.\n')

    lines.append('## Top-5 |B → C| fractional shifts in D_median\n')
    lines.append('| Sample | Instrument | B D_median (Å) | C D_median (Å) | Δ_BC % |')
    lines.append('|---|---|---|---|---|')
    for r in top5:
        ik = next((s['instrument_key'] for s in config.SAMPLES
                   if s['name'] == r['Sample']), '?')
        lines.append(
            f'| {r["Sample"]} | {ik} | '
            f'{r["B_WA_D_median_A"]:.1f} | {r["C_WA_D_median_A"]:.1f} | '
            f'{r["delta_BC_D_median_pct"]:.2f}% |')
    lines.append('')

    lines.append('## Caveat\n')
    lines.append(
        '> Column C is conditional on the Caglioti values listed above. '
        'TYPICAL values in the table are placeholders; LITERATURE values '
        'tighten this from "magnitude estimate" to "literature-grade '
        'estimate"; MEASURED LaB6 calibration on each instrument would '
        'tighten further. A FAIL in the A vs B summary on the Cu+I11 '
        'subset is a strict-additive contract violation and a real bug '
        'in the package — must be fixed before column C is used for '
        'science interpretation.\n')

    out_path.write_text('\n'.join(lines), encoding='utf-8')
    print(f'Wrote {out_path}')
```

- [ ] **Step 2: Wire write_run_log into main()**

Modify the end of `main()` to capture the start time and call `write_run_log()`:

Replace the `print('\nLoop complete...'` line at the end of `main()` with:

```python
    write_run_log(start_time=start_time,
                   comparison_rows=comparison_rows)
    print('\nLoop and run_log.md complete. (Figures are added in '
          'subsequent tasks.)')
```

…and at the **top** of `main()`, immediately after the heading prints, capture the start time:

```python
    start_time = datetime.now()
```

(The existing `print(f'Run start: ...')` line uses `datetime.now()` directly — replace that line to use `start_time` instead so both prints and run_log show the same instant.)

The relevant snippet at the top of main() should now read:

```python
def main():
    start_time = datetime.now()
    print('=' * 60)
    print('xrd_profile v0.3.0 vs v0.4.0 survey comparison harness')
    print('=' * 60)
    print(f'Run start:    {start_time.isoformat()}')
    import xrd_profile
    print(f'xrd_profile:  v{xrd_profile.__version__}')
    print(f'Output dir:   {config.OUTPUT_DIR}')
    Path(config.OUTPUT_DIR).mkdir(parents=True, exist_ok=True)
    # ... rest unchanged ...
```

- [ ] **Step 3: Run the full survey end-to-end**

```bash
"/c/Users/Matthew Izawa/anaconda3/python.exe" "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/comparison_v040/run_comparison.py"
```

Expected runtime: 10–20 minutes. Expected final output:

```
[5/5] Running per-sample loop (29 samples × 2 W-A)...
  [ 1/29] Gadamis 004      instrument=cu_misasa     B fams=16  C fams=<int>  A-vs-B=PASS
  [ 2/29] NWA 11182a       instrument=cu_misasa     B fams=11  C fams=<int>  A-vs-B=PASS
  ...
  [29/29] Tatahouine       instrument=co_winnipeg   B fams=9   C fams=<int>  A-vs-B=PASS

Wrote .../output/comparison.csv
Wrote .../output/size_distributions.csv
Wrote .../output/run_log.md

Loop and run_log.md complete. (Figures are added in subsequent tasks.)
```

- [ ] **Step 4: Verify pass criteria (Spec §9.3)**

Open `output/run_log.md` and confirm:

```
- Cu Kα Misasa:   7/7 PASS, 0 FAIL, 0 N/A
- Diamond I11:    17/18 PASS, 0 FAIL, 1 N/A    ← NWA 7042 is N/A (saved CSV has families=0)
- Co Kα Winnipeg: 4/4 PASS, 0 FAIL, 0 N/A
```

If any FAIL appears on the Cu+I11 subset (8/25 expected to be PASS plus 17/18 = 24/25 plus 1 N/A — actually 24 PASS and 1 N/A across Cu+I11), this is a strict-additive contract violation in the package. Stop and investigate before continuing to figures.

If a FAIL appears on the Co Kα subset, this is informational — the saved CSV likely came from a different historical processing pipeline. Note the deltas in the run log and continue.

- [ ] **Step 5: Inspect the top-5 delta table in run_log.md**

The largest |B→C| fractional shifts should be on lab Cu Kα and Co Kα samples (where instrumental contribution is biggest). If the largest shifts are on synchrotron I11 samples, that's a red flag — the synthetic LaB6 for I11 is too wide; check the I11 Caglioti in `config.py`.

Expected qualitative pattern (per Spec §8.1):
- Lab Cu+Co lines in the top-5: −20% to −50% shifts.
- I11 lines (if any in top-5): −1% to −5% shifts.

If the pattern is reversed, edit the I11 Caglioti to make it sharper (smaller U/W) and re-run. The TYPICAL flag stays in run_log.md.

- [ ] **Step 6: Save in place**

```bash
ls -la "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/comparison_v040/output/"
```

Expected: `comparison.csv`, `size_distributions.csv`, `run_log.md`, `synthetic_standards/` all present.

---

## Task 7: plots.py — figs 01 and 02 (D_median scatter plots)

**Files:**
- Create: `Llunr/comparison_v040/plots.py`

The two scatter plots are the simplest figures and share a lot of code (same axes, same colour-by-instrument scheme, just different x/y data). Implement them together.

- [ ] **Step 1: Create plots.py with figs 01 and 02**

File: `Llunr/comparison_v040/plots.py`

```python
"""
Figure rendering for the v0.3.0 vs v0.4.0 comparison harness.

Reads output/comparison.csv and output/size_distributions.csv and
writes the five publication-track figures specified in Spec §7.3.

Run:
    python plots.py            # render all five figures
    python plots.py 01 02      # render only figs 01 and 02
"""
import csv
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import config

OUTPUT_DIR = Path(config.OUTPUT_DIR)
INSTRUMENT_COLOURS = {
    'cu_misasa':   '#1f77b4',  # blue
    'co_winnipeg': '#d62728',  # red
    'i11':         '#2ca02c',  # green
}
INSTRUMENT_LABELS = {
    'cu_misasa':   'Cu Kα (Misasa)',
    'co_winnipeg': 'Co Kα (Winnipeg)',
    'i11':         'Diamond I11',
}


def _load_comparison():
    rows = []
    with open(OUTPUT_DIR / 'comparison.csv', encoding='utf-8') as f:
        for r in csv.DictReader(f):
            rows.append(r)
    return rows


def _instrument_for(name):
    """Look up instrument_key from config.SAMPLES by sample name."""
    for s in config.SAMPLES:
        if s['name'] == name:
            return s['instrument_key']
    raise KeyError(name)


def _to_float_safe(s):
    if s is None or s == '':
        return float('nan')
    try:
        return float(s)
    except ValueError:
        return float('nan')


def _save(fig, base_name):
    for ext in ('png', 'svg'):
        path = OUTPUT_DIR / f'{base_name}.{ext}'
        fig.savefig(path, dpi=300)
        print(f'Wrote {path}')
    plt.close(fig)


def fig01_d_median_a_vs_b():
    """A vs B regression scatter: should hug y = x."""
    rows = _load_comparison()
    fig, ax = plt.subplots(figsize=(6, 6), constrained_layout=True)
    for ik in INSTRUMENT_COLOURS:
        a_vals, b_vals = [], []
        for r in rows:
            if _instrument_for(r['Sample']) != ik:
                continue
            a = _to_float_safe(r['A_WA_D_median_A'])
            b = _to_float_safe(r['B_WA_D_median_A'])
            if np.isnan(a) or np.isnan(b):
                continue
            a_vals.append(a)
            b_vals.append(b)
        ax.scatter(a_vals, b_vals, c=INSTRUMENT_COLOURS[ik],
                   label=INSTRUMENT_LABELS[ik], s=50, alpha=0.8,
                   edgecolors='k', linewidth=0.5)
    # y = x diagonal.
    a_all = [_to_float_safe(r['A_WA_D_median_A']) for r in rows]
    a_all = [v for v in a_all if not np.isnan(v)]
    if a_all:
        lo, hi = 0, max(a_all) * 1.1
        ax.plot([lo, hi], [lo, hi], 'k--', linewidth=0.8, alpha=0.5,
                label='y = x')
    ax.set_xlabel('A: saved CSV W-A D$_{median}$ (Å)')
    ax.set_ylabel('B: v0.4.0 instrumental=None W-A D$_{median}$ (Å)')
    ax.set_title('Fig 01 — Regression check: A vs B')
    ax.legend(loc='lower right')
    _save(fig, 'fig01_D_median_AvB')


def fig02_d_median_b_vs_c():
    """B vs C sensitivity scatter: distance from y=x is the correction effect."""
    rows = _load_comparison()
    fig, ax = plt.subplots(figsize=(6, 6), constrained_layout=True)
    for ik in INSTRUMENT_COLOURS:
        b_vals, c_vals = [], []
        for r in rows:
            if _instrument_for(r['Sample']) != ik:
                continue
            b = _to_float_safe(r['B_WA_D_median_A'])
            c = _to_float_safe(r['C_WA_D_median_A'])
            if np.isnan(b) or np.isnan(c):
                continue
            b_vals.append(b)
            c_vals.append(c)
        ax.scatter(b_vals, c_vals, c=INSTRUMENT_COLOURS[ik],
                   label=INSTRUMENT_LABELS[ik], s=50, alpha=0.8,
                   edgecolors='k', linewidth=0.5)
    b_all = [_to_float_safe(r['B_WA_D_median_A']) for r in rows]
    b_all = [v for v in b_all if not np.isnan(v)]
    if b_all:
        lo, hi = 0, max(b_all) * 1.1
        ax.plot([lo, hi], [lo, hi], 'k--', linewidth=0.8, alpha=0.5,
                label='y = x (no correction)')
    ax.set_xlabel('B: v0.4.0 instrumental=None D$_{median}$ (Å)')
    ax.set_ylabel('C: v0.4.0 with synthetic LaB6 D$_{median}$ (Å)')
    ax.set_title('Fig 02 — Instrumental correction effect: B vs C')
    ax.legend(loc='lower right')
    _save(fig, 'fig02_D_median_BvC')


def main(args):
    selected = set(args) if args else {'01', '02', '03', '04', '05'}
    if '01' in selected:
        fig01_d_median_a_vs_b()
    if '02' in selected:
        fig02_d_median_b_vs_c()
    # Figures 03, 04, 05 are added in Task 8.


if __name__ == '__main__':
    main(sys.argv[1:])
```

- [ ] **Step 2: Render figs 01 and 02**

```bash
"/c/Users/Matthew Izawa/anaconda3/python.exe" "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/comparison_v040/plots.py" 01 02
```

Expected output:
```
Wrote .../output/fig01_D_median_AvB.png
Wrote .../output/fig01_D_median_AvB.svg
Wrote .../output/fig02_D_median_BvC.png
Wrote .../output/fig02_D_median_BvC.svg
```

- [ ] **Step 3: Visual sanity check**

Open `output/fig01_D_median_AvB.png` in an image viewer.

Expected: 28 points (29 minus NWA 7042 N/A) all sitting essentially exactly on the y=x diagonal — this is the strict-additive contract visualised. Any visible off-diagonal point means a real regression in the package.

Open `output/fig02_D_median_BvC.png`.

Expected: lab Cu (blue) and Co (red) points sit *below* y=x (corrected D < uncorrected D — instrumental subtraction makes peaks intrinsically narrower → smaller D). Synchrotron I11 (green) points should sit on or just below y=x (small correction). If green points sit far below y=x, the I11 synthetic standard is too wide; tune the I11 Caglioti in `config.py` and re-run Task 6.

- [ ] **Step 4: Save in place**

```bash
ls "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/comparison_v040/plots.py"
ls -la "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/comparison_v040/output/"*.png
```

Expected: `plots.py` + 4 PNG files (fig01 and fig02 in two formats each).

---

## Task 8: plots.py — figs 03 (boxplot), 04 (FWHM panels), 05 (size-distribution panels)

**Files:**
- Modify: `Llunr/comparison_v040/plots.py`

These three figures are more complex than the scatters. Implement one at a time, with a visual sanity check after each.

- [ ] **Step 1: Add fig 03 — Δ_BC % vs shock stage, boxplot per instrument**

Add to `plots.py` (above `def main(args)`):

```python
SHOCK_ORDER = ['S1', 'S1-S2', 'S2', 'S2-S3', 'S3', 'S3-S4',
                'S4', 'S4-S5', 'S5', 'S5-S6', 'S6']


def fig03_delta_bc_by_shock():
    """Boxplot of Δ_BC (%) on D_median, grouped by shock stage, one box per instrument."""
    rows = _load_comparison()
    # Collect Δ_BC % per (shock, instrument).
    by_key = {}
    for r in rows:
        ik = _instrument_for(r['Sample'])
        d = _to_float_safe(r['delta_BC_D_median_pct'])
        if np.isnan(d):
            continue
        by_key.setdefault((r['Shock'], ik), []).append(d)

    # Order shocks present in the data by SHOCK_ORDER.
    shocks_in_data = sorted({k[0] for k in by_key},
                             key=lambda s: (SHOCK_ORDER.index(s)
                                              if s in SHOCK_ORDER
                                              else len(SHOCK_ORDER)))

    fig, ax = plt.subplots(figsize=(10, 5), constrained_layout=True)
    width = 0.25
    instruments = list(INSTRUMENT_COLOURS.keys())
    for j, ik in enumerate(instruments):
        positions = []
        data = []
        for i, sh in enumerate(shocks_in_data):
            vals = by_key.get((sh, ik), [])
            if vals:
                positions.append(i + (j - 1) * width)
                data.append(vals)
        if data:
            bp = ax.boxplot(data, positions=positions, widths=width * 0.8,
                            patch_artist=True, manage_ticks=False)
            for box in bp['boxes']:
                box.set_facecolor(INSTRUMENT_COLOURS[ik])
                box.set_alpha(0.6)
            for med in bp['medians']:
                med.set_color('k')

    ax.set_xticks(range(len(shocks_in_data)))
    ax.set_xticklabels(shocks_in_data)
    ax.set_xlabel('Shock stage')
    ax.set_ylabel('Δ$_{BC}$ % on D$_{median}$  (= 100·(C − B)/B)')
    ax.set_title('Fig 03 — Instrumental correction effect by shock stage and instrument')
    ax.axhline(0, color='k', linewidth=0.5, alpha=0.5)
    ax.grid(True, axis='y', alpha=0.3)
    handles = [plt.Rectangle((0, 0), 1, 1, color=INSTRUMENT_COLOURS[ik],
                              alpha=0.6, label=INSTRUMENT_LABELS[ik])
               for ik in instruments]
    ax.legend(handles=handles, loc='upper right')
    _save(fig, 'fig03_delta_BC_by_shock')
```

Update `main()`:

```python
def main(args):
    selected = set(args) if args else {'01', '02', '03', '04', '05'}
    if '01' in selected:
        fig01_d_median_a_vs_b()
    if '02' in selected:
        fig02_d_median_b_vs_c()
    if '03' in selected:
        fig03_delta_bc_by_shock()
    # Figures 04, 05 are added below.
```

Render and inspect:

```bash
"/c/Users/Matthew Izawa/anaconda3/python.exe" "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/comparison_v040/plots.py" 03
```

Expected: Δ_BC % is negative (column C ≤ column B by physics — correction reduces D); blue/red boxes ~−20% to −50%; green boxes ~−1% to −5%. If the boxes are positive, sign convention is wrong (review the Δ definition in `comparison_row()`).

- [ ] **Step 2: Add fig 04 — instrumental vs sample FWHM, three-panel**

For each instrument, pick a sample whose B-vs-C delta is roughly median for that instrument (representative, not extreme). The implementer chooses; suggested defaults: `'NWA 11182a'` (Cu), `'Tirhert'` (I11), `'Talampaya'` (Co).

Add to `plots.py`:

```python
REPRESENTATIVE_SAMPLES = {
    'cu_misasa':   'NWA 11182a',   # Cu Kα Misasa
    'i11':         'Tirhert',      # synchrotron
    'co_winnipeg': 'Talampaya',    # Co Kα Winnipeg
}


def fig04_instrumental_contribution():
    """Three-panel: Caglioti FWHM(2θ) vs measured-peak FWHMs for one
    representative sample per instrument."""
    from synthetic_standards import caglioti_fwhm
    from xrd_profile import XRDProfile

    fig, axes = plt.subplots(1, 3, figsize=(15, 4),
                              constrained_layout=True, sharey=True)
    for ax, (ik, sample_name) in zip(axes, REPRESENTATIVE_SAMPLES.items()):
        cag = config.INSTRUMENTAL_CAGLIOTI[ik]
        # Caglioti curve over relevant 2θ range.
        s = next(s for s in config.SAMPLES if s['name'] == sample_name)
        tt_range = s['tt_range']
        tt_curve = np.linspace(tt_range[0] + 1, tt_range[1] - 1, 200)
        fwhm_curve = caglioti_fwhm(tt_curve, cag['U'], cag['V'], cag['W'])
        ax.plot(tt_curve, fwhm_curve, 'k-', linewidth=1.5,
                label=f'Caglioti (synthetic LaB6, {cag["flag"]})')

        # Sample's measured peak FWHMs from a quick guided W-H pass.
        tt, intensity = np.loadtxt(s['file'], comments='#').T[:2]
        mask = (tt >= tt_range[0]) & (tt <= tt_range[1])
        profile = XRDProfile(tt[mask], intensity[mask], s['wavelength'],
                              sample_name=s['name'])
        # Use guided_williamson_hall to get peak list with FWHMs.
        # build_phase_refs result is not in scope here; rebuild ad-hoc.
        from xrd_profile import Phase
        if s['phase_ref'] == 'anorthite':
            phase = Phase.from_lattice_params(
                a=8.1809, b=12.881, c=7.1101,
                alpha=93.465, beta=116.11, gamma=90.369,
                species=['Ca', 'Al', 'Al', 'Si', 'Si',
                          'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O'],
                coords=[[0.269,0.988,0.086],[0.507,0.314,0.621],
                        [0.992,0.815,0.118],[0.505,0.320,0.110],
                        [0.006,0.816,0.613],[0.491,0.625,0.487],
                        [0.024,0.124,0.995],[0.073,0.488,0.635],
                        [0.576,0.990,0.143],[0.298,0.356,0.612],
                        [0.817,0.855,0.142],[0.517,0.179,0.610],
                        [0.000,0.680,0.104]],
                name='anorthite')
            ref_peaks = phase.get_ref_peaks(s['wavelength'],
                                              two_theta_range=tt_range,
                                              min_intensity=3.0)
        else:
            ref_peaks = []
        if ref_peaks:
            wh = profile.guided_williamson_hall(
                ref_peaks, n_sigma=3.0, tolerance_d=0.02)
            peaks = wh.get('peaks', {})
            tt_obs = peaks.get('two_theta_obs', [])
            fwhm_obs = peaks.get('fwhm', [])
            if len(tt_obs) > 0:
                ax.scatter(tt_obs, fwhm_obs,
                           c=INSTRUMENT_COLOURS[ik], s=40,
                           edgecolors='k', linewidth=0.5,
                           label=f'{sample_name} measured FWHM')
        ax.set_xlabel('2θ (°)')
        ax.set_title(f'{INSTRUMENT_LABELS[ik]}')
        ax.legend(loc='upper left', fontsize=8)
        ax.grid(True, alpha=0.3)
    axes[0].set_ylabel('FWHM (°)')
    fig.suptitle('Fig 04 — Instrumental vs measured peak FWHM '
                 'per instrument')
    _save(fig, 'fig04_instrumental_contribution')
```

(Note: the `phase_ref` lookup in fig 04 is anorthite-only because all three representative samples use anorthite. If the implementer chooses a non-anorthite representative, mirror the lookup pattern from `run_comparison.build_phase_refs`.)

Update `main()`:

```python
    if '04' in selected:
        fig04_instrumental_contribution()
```

Render and inspect:

```bash
"/c/Users/Matthew Izawa/anaconda3/python.exe" "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/comparison_v040/plots.py" 04
```

Expected: in each panel, the Caglioti curve is monotonically increasing with 2θ (W dominates near 0°; U·tan² dominates at high 2θ); measured-FWHM points sit ABOVE the Caglioti curve (sample broadening adds in quadrature on top of instrumental). Margin between Caglioti and measured points is the sample contribution. Lab panels show ~30–50% gap; I11 panel shows much wider gap (most observed broadening is sample-driven).

- [ ] **Step 3: Add fig 05 — size-distribution panels**

This figure picks one sample (default: Tirhert) and shows 4 representative families' size distributions.

Add to `plots.py`:

```python
SIZE_DIST_SAMPLE = 'Tirhert'   # representative for fig 05


def fig05_size_distribution_panels():
    """4 panels: A_size(L) data + lognormal + normal fits for 4 families
    of one representative sample."""
    # Re-run W-A on this one sample to get the A_size(L) arrays
    # (size_distributions.csv has only the fitted params, not the underlying L/A_size).
    # Build the sample, ref_peaks, and InstrumentalStandard.
    sample = next(s for s in config.SAMPLES if s['name'] == SIZE_DIST_SAMPLE)
    # Reuse run_comparison's helpers.
    from run_comparison import (build_phase_refs,
                                  build_instrumental_standards, load_xy)
    phase_refs = build_phase_refs()
    standards = build_instrumental_standards()
    ref_peaks = phase_refs[(sample['phase_ref'], sample['wavelength'])]
    std = standards[sample['instrument_key']]
    tt, intensity = load_xy(sample['file'], sample['tt_range'])

    from xrd_profile import XRDProfile
    profile = XRDProfile(tt, intensity, sample['wavelength'],
                          sample_name=sample['name'])
    wa = profile.guided_warren_averbach(
        ref_peaks, n_sigma=3.0, tolerance_d=0.02, instrumental=std)

    families = wa.get('families', [])
    # Pick 4 families with the most valid L points (most informative fits).
    families_sorted = sorted(
        families,
        key=lambda fm: -fm.get('size_distribution', {}).get('n_valid_L', 0))
    chosen = families_sorted[:4]
    if len(chosen) < 4:
        print(f'WARNING: only {len(chosen)} families available for fig 05')

    fig, axes = plt.subplots(2, 2, figsize=(10, 8), constrained_layout=True)
    from synthetic_standards import A_LAB6   # reuse for any constants if needed
    for ax, fam in zip(axes.flat, chosen):
        sd = fam.get('size_distribution', {})
        L = np.asarray(fam.get('L', []))
        A = np.asarray(fam.get('A_size', []))
        ax.scatter(L, A, c='k', s=20, label='A_size(L) data', zorder=3)
        # Lognormal fit.
        ln = sd.get('lognormal', {})
        if ln and not np.isnan(ln.get('D_median', float('nan'))):
            from xrd_profile.size_distributions import lognormal_a_size
            L_smooth = np.linspace(L.min(), L.max(), 200)
            ax.plot(L_smooth,
                    lognormal_a_size(L_smooth, ln['D_median'], ln['sigma']),
                    'b-', linewidth=1.5,
                    label=f'lognormal R²={ln["R2"]:.2f}')
        nm = sd.get('normal', {})
        if nm and not np.isnan(nm.get('D_mean', float('nan'))):
            from xrd_profile.size_distributions import normal_a_size
            L_smooth = np.linspace(L.min(), L.max(), 200)
            ax.plot(L_smooth,
                    normal_a_size(L_smooth, nm['D_mean'], nm['sigma']),
                    'r--', linewidth=1.5,
                    label=f'normal R²={nm["R2"]:.2f}')
        peaks = fam.get('peaks', [])
        if peaks:
            rep = min(peaks, key=lambda p: (
                p.get('h', 0)**2 + p.get('k', 0)**2 + p.get('l', 0)**2))
            ax.set_title(f'family ({rep["h"]} {rep["k"]} {rep["l"]})')
        ax.set_xlabel('Column length L (Å)')
        ax.set_ylabel('A$_{size}$(L)')
        ax.legend(loc='upper right', fontsize=8)
        ax.grid(True, alpha=0.3)
    fig.suptitle(f'Fig 05 — Size-distribution fits, {SIZE_DIST_SAMPLE} '
                 f'(top-4 families by n_valid_L)')
    _save(fig, 'fig05_size_distribution_panels')
```

Update `main()`:

```python
    if '05' in selected:
        fig05_size_distribution_panels()
```

Render and inspect:

```bash
"/c/Users/Matthew Izawa/anaconda3/python.exe" "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/comparison_v040/plots.py" 05
```

Expected: 4 panels each showing scattered black points (A_size(L) data) decaying from 1.0 at L=0 toward 0 at high L; blue solid line (lognormal fit) and red dashed line (normal fit) overlaid. R² values typically > 0.95 for well-behaved families. If a panel is empty, the family's W-A fit lacked enough valid L points (n_valid_L < ~5).

- [ ] **Step 4: Render all 5 figures end-to-end**

```bash
"/c/Users/Matthew Izawa/anaconda3/python.exe" "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/comparison_v040/plots.py"
```

Expected: 5 figures × 2 formats = 10 files written to `output/`.

- [ ] **Step 5: Wire plots.py into run_comparison.py main()**

At the end of `run_comparison.py`'s `main()` (after `write_run_log`), add:

```python
    print('\nRendering figures...')
    import plots
    plots.main([])
    print('All figures rendered.')
```

- [ ] **Step 6: Final end-to-end run from a clean output dir**

Move/rename the existing output dir to confirm a fresh run produces everything:

```bash
mv "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/comparison_v040/output" \
   "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/comparison_v040/output_prev"
"/c/Users/Matthew Izawa/anaconda3/python.exe" "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/comparison_v040/run_comparison.py"
ls "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/comparison_v040/output/"
```

Expected: full 13+ files present after one `python run_comparison.py` invocation
(comparison.csv, size_distributions.csv, run_log.md, synthetic_standards/<3 xy>,
fig01-05 × {png, svg}).

- [ ] **Step 7: Save in place**

```bash
ls "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/comparison_v040/plots.py"
ls "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/comparison_v040/output/"*.svg
```

Expected: `plots.py` exists; 5 .svg figures present.

---

## Task 9: Commit spec + plan to xrd_profile git repo

**Files:**
- Modify (in git): `xrd_profile/docs/superpowers/specs/2026-05-06-v030-v040-survey-comparison-design.md` (already written by brainstorming skill)
- Modify (in git): `xrd_profile/docs/superpowers/plans/2026-05-06-v030-v040-survey-comparison.md` (already written by writing-plans skill)

This task lands the design doc and this implementation plan in the package's git history. Implementation files at `Llunr/comparison_v040/*.py` remain outside git per the repository note in the header.

- [ ] **Step 1: Inspect git status in xrd_profile**

```bash
cd "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/xrd_profile"
git status
```

Expected: untracked files
```
Untracked files:
  docs/superpowers/specs/2026-05-06-v030-v040-survey-comparison-design.md
  docs/superpowers/plans/2026-05-06-v030-v040-survey-comparison.md
```

- [ ] **Step 2: Stage the spec and plan**

```bash
cd "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/xrd_profile"
git add docs/superpowers/specs/2026-05-06-v030-v040-survey-comparison-design.md \
        docs/superpowers/plans/2026-05-06-v030-v040-survey-comparison.md
git status
```

Expected: both files now under "Changes to be committed".

- [ ] **Step 3: Commit**

```bash
cd "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/xrd_profile"
git commit -m "$(cat <<'EOF'
Add v0.3.0 vs v0.4.0 survey comparison design and plan

Validation harness for the 29-sample JAC survey under v0.4.0:
column B (instrumental=None) checks the strict-additive contract on
real data; column C (synthetic LaB6 InstrumentalStandard built from
literature Caglioti per instrument) quantifies the new methods'
effect. Implementation lives at Llunr/comparison_v040/ (outside git,
per repo conventions). Spec + plan landed here as the package-level
record of design intent.
EOF
)"
```

Expected: commit succeeds; new commit visible in `git log -1`.

- [ ] **Step 4: Verify commit**

```bash
cd "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/xrd_profile"
git log -1 --stat
```

Expected: shows the two new files as added.

- [ ] **Step 5: Do NOT push without explicit user approval**

Per `xrd_profile/CLAUDE.md`: "Never push to remote without explicit user approval." Stop here. The user decides whether to push the spec + plan to `origin/main`.

---

## Self-review notes

**Spec coverage** — every section of the spec maps to one or more tasks:

| Spec section | Task |
|---|---|
| §1 Background, §2 Goals, §3 Non-goals | Header of plan |
| §4.1 Layout, §4.2 Components, §4.3 Read-only contracts | Task 1 (dirs), Tasks 2-5 (files) |
| §5.1 Wavelengths/paths | Task 3 |
| §5.2 INSTRUMENTAL_CAGLIOTI | Task 3 |
| §5.3 SAMPLES (29 entries) | Task 3 |
| §5.4 PHASE_REFS | Task 3 (declaration), Task 4 (build) |
| §5.5 Startup assertions | Task 4 |
| §6.1 Startup phase | Task 4 |
| §6.2 Per-sample loop | Task 5 |
| §6.3 W-A requires InstrumentalStandard (not InstrumentalProfile) | Task 2 (synthetic standard) |
| §6.4 A vs B tolerance | Task 5 (`check_a_vs_b`) |
| §6.5 NaN handling | Task 5 (NaN-check returns 'N/A') |
| §7.1 comparison.csv schema | Task 5 |
| §7.2 size_distributions.csv schema | Task 5 |
| §7.3 Five figures | Tasks 7 (figs 01, 02), 8 (figs 03, 04, 05) |
| §8 Caglioti per instrument | Task 3 (TYPICAL placeholders), Task 6 (run_log records flag) |
| §8.1 Anticipated qualitative result | Task 6 (sanity check after first run) |
| §8.2 run_log.md per-run capture | Task 6 |
| §9.1 Failure modes | Tasks 4, 5 |
| §9.2 Reproducibility | Task 6 (run_log captures version + SHA + Caglioti) |
| §9.3 Testing pass criteria | Task 6 step 4 |
| §9.4 Performance | Header, Task 6 (10–20 min) |
| §10 Open items | Header, Task 3 (USER INPUT NEEDED placeholders) |

**Placeholder scan** — all `<USER INPUT NEEDED>` sentinels are sample-data paths the user fills in once before first run; the harness fails loudly with a clear message if they remain. No "TBD" / "implement later" / vague-error-handling phrases in the plan body. Every code step shows complete code.

**Type consistency** — function signatures match across tasks:
- `build_synthetic_lab6(wavelength, U, V, W, output_path) → (xy_path, InstrumentalStandard)` (Task 2, used in Task 4)
- `check_a_vs_b(saved_row, wa_b, instrument_key) → dict` (Task 5)
- `run_one_sample(s, ref_peaks, std) → (wa_B, wa_C)` (Task 5)
- `comparison_row(s, saved_row, wa_B, wa_C, ab_check) → dict` (Task 5)
- `write_run_log(start_time, comparison_rows)` (Task 6)
- `INSTRUMENT_COLOURS` and `INSTRUMENT_LABELS` consistent across figs 01-04 (Tasks 7, 8)

---

## Execution outcome (2026-05-07, γ-scope)

The plan was executed via subagent-driven-development; the actual outcome diverged from the original task list because of two findings discovered during Tasks 4–6. Recorded here for the on-disk record. See spec §12 for the prose narrative.

### Task-by-task status

| # | Task | Status | Notes |
|---|---|---|---|
| 1 | Setup directory + README | ✅ shipped as planned | |
| 2 | `synthetic_standards.py` + pytest | ✅ shipped with two reviewer fixes | LaB6 boron coords corrected to NIST SRM 660c (0.1993/0.8007 not 0.2/0.8); zero-intensity guard added; `tests/conftest.py` factored out the `sys.path` boilerplate; `tt_step` default `0.005` not `0.02` (FWHM half-max walk required finer sampling). |
| 3 | `config.py` | ✅ shipped with sample-count correction | The saved CSV `survey_results_29samples.csv` actually has **28** sample rows, not 29 (filename misleading). 7 Cu + 17 I11 + 4 Co Kα. The 4 Co Kα entries' `file` paths and source dir are still `<USER INPUT NEEDED>` placeholders. |
| 4–6 | `run_comparison.py` startup + per-sample loop + run_log | ✅ shipped with min_intensity correction + γ-scope run_log narrative | One-line bug found: `Phase.get_ref_peaks(min_intensity=3.0)` filters anorthite to ~96 reflections; `compile_survey.py` passes the unfiltered XRDCalculator output (~17 000 reflections) directly. Changed harness to `min_intensity=0.0`. Result: column B reproduces saved CSV row-for-row at floating-point precision. Also: graceful Co Kα skip-with-warning when paths are placeholders; UTF-8 stdout wrapper; sys.path order fixed so local `plots.py` shadows `xrd_profile/plots.py`. |
| 7 | `plots.py` figs 01 + 02 | ⚠️ partial — fig 01 only | Fig 01 (A vs B regression scatter) ships. Figs 02 (B vs C scatter) deferred — column C is empty across the board. |
| 8 | `plots.py` figs 03, 04, 05 | ⏸️ deferred to v0.5.0 | All three figures depend on column C being populated. |
| 9 | Commit spec + plan to xrd_profile git | ✅ shipped (this commit) | Spec §12 (execution result) and this section added. |

### Why figs 02–05 are deferred — the column C blocker

The Stokes deconvolution path in `xrd_profile/instrumental.py:375–412` extracts from the standard pattern at the *sample's* 2θ, not the standard's nearest peak. With our synchrotron-scale synthetic LaB6 (narrow Gaussians on near-zero baseline), every off-peak extraction yields zero area → `A_inst[0] = 0` → `cannot Stokes-deconvolve`. This affected 24/24 samples in the run.

Spec §12.3 enumerates the v0.5.0 fix options. Preferred: synthesise a Gaussian at the sample's 2θ from the fitted Caglioti U/V/W. Once that lands, `Llunr/comparison_v040/run_comparison.py` re-runs without code changes and column C populates.

### What's on disk in the harness output dir

`Llunr/comparison_v040/output/` after the 2026-05-07 run:
- `comparison.csv` (24 rows; columns A and B populated, column C all `nan`)
- `size_distributions.csv` (header only — column C produces no families)
- `synthetic_standards/lab6_{cu_misasa,co_winnipeg,i11}.xy` (built but unused once column C aborts)
- `fig01_D_median_AvB.{png,svg}` (the validation visual — all 23 points on y=x)
- `run_log.md` (the human-readable summary, including the v0.5.0 ticket text)
