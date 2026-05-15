# JAC N=32 Expansion Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Integrate the 24 Apr 2026 Diamond Light Source I11 data (4 new samples, 2 lab→synchrotron upgrades, 1 blank capillary) into the JAC manuscript and produce a v2 .docx bundle for a single co-author re-circulation round, ahead of JAC submission.

**Architecture:** Five sequential phases with explicit author-review gates. Phase 1 validates capillary subtraction on 3 test eucrites. Phase 2 produces `survey_results_32samples.csv` and triggers a sanity-check gate. Phase 3 regenerates all data-bearing figures into `Paper1_JAC/paper1_figures/`. Phase 4 updates manuscript prose, tables, cover letter, and submission documents. Phase 5 packages v2 .docx with track changes and a cover note.

**Tech Stack:** Python 3.13.9 (Anaconda base env at `C:\Users\Matthew Izawa\anaconda3\python.exe`); numpy, scipy, matplotlib, pymatgen, pytest, python-docx (already in the base env per `xrd_profile/CLAUDE.md`); `xrd_profile` v0.3.x importable from `Llunr/xrd_profile/` (read-only API consumer — no package changes); pandoc for the track-changes .docx diff.

**Repository note:** `Llunr/` is **not** a git repo; only `xrd_profile/` is. The spec at `xrd_profile/docs/superpowers/specs/2026-05-15-jac-n33-expansion-design.md` and this plan are in the package's git tree, but all analysis code at `Llunr/n33_expansion/*.py` and all manuscript work at `Paper1_JAC/` are project-local and outside git. "Commit" steps in this plan apply only to the spec and plan in `xrd_profile/`. The spec was already committed at `3ea4b17` on `main`; this plan is committed at the end of Task 1.

**TDD scope:** Tasks 2, 4, 8, 12, 13 ship genuine new Python code and follow TDD (write failing test → run red → minimal impl → run green → commit-equivalent). All other tasks are procedural — they read or edit existing scripts, regenerate figures, or modify Markdown prose; verification is by inspection or by running the script and inspecting output. The plan flags each task with `[TDD]` or `[PROCEDURAL]` or `[AUTHOR TASK]`.

**Test command (n33_expansion-internal pytest):**
```
"/c/Users/Matthew Izawa/anaconda3/python.exe" -m pytest "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/n33_expansion/tests/" -v
```

**Anaconda Python path** (used throughout):
```
"/c/Users/Matthew Izawa/anaconda3/python.exe"
```

---

## File Structure

| Path | Role |
|---|---|
| `Llunr/n33_expansion/README.md` | How to run end-to-end, expected outputs, prerequisites |
| `Llunr/n33_expansion/scripts/reduce_dat_to_xy.py` | Convert DLS `.dat` (header + tab-separated 2θ-counts-error) to xrd_profile-readable `.xy` |
| `Llunr/n33_expansion/scripts/capillary_subtract.py` | High-Q-tail scaling and capillary subtraction; importable as a module |
| `Llunr/n33_expansion/scripts/phase1_validation.py` | Phase 1 driver: process capillary, subtract from 3 test eucrites, generate diagnostic.png |
| `Llunr/n33_expansion/scripts/process_all_samples.py` | Phase 2 driver: run xrd_profile on 6 new samples and re-fit PDF for 17 existing |
| `Llunr/n33_expansion/scripts/build_survey_csv.py` | Aggregate per-sample results into `survey_results_32samples.csv` |
| `Llunr/n33_expansion/scripts/build_cross_instrument.py` | Build `cross_instrument_pairs.csv` for Table S2 |
| `Llunr/n33_expansion/tests/test_reduce_dat.py` | Unit tests for DLS .dat parser |
| `Llunr/n33_expansion/tests/test_capillary_subtract.py` | Unit tests for scaling and subtraction |
| `Llunr/n33_expansion/tests/test_build_survey_csv.py` | Unit tests for CSV aggregation row construction |
| `Llunr/n33_expansion/data/i11_2026Apr_reduced/` | Reduced `.xy` from new `.dat` (7 files) |
| `Llunr/n33_expansion/data/i11_corrected/` | Capillary-subtracted `.xy` (23 files: 6 new + 17 existing) |
| `Llunr/n33_expansion/phase1_validation/diagnostic.png` | Phase 1 GO/PAUSE evidence figure |
| `Llunr/n33_expansion/phase1_validation/decision_note.md` | Phase 1 decision record |
| `Llunr/n33_expansion/results/per_sample/*.json` | Per-sample xrd_profile output for 6 new analyses |
| `Llunr/n33_expansion/results/per_sample_pdf_refits/*.json` | PDF-only refit output for 17 existing I11 samples |
| `Llunr/n33_expansion/results/cross_instrument_pairs.csv` | Lab Co + I11 W-A values for NWA 5751 and Talampaya |
| `Paper1_JAC/survey_results_32samples.csv` | Final canonical per-sample CSV |
| `Paper1_JAC/paper1_figures/*.png` | Regenerated figures |
| `Paper1_JAC/manuscript_paper1_JAC.md` | Updated manuscript |
| `Paper1_JAC/table1_sample_inventory.md` | Updated Table 1 source |
| `Paper1_JAC/Table_S2_cross_instrument.md` | New Table S2 source |
| `Paper1_JAC/submission/cover_letter.md` | Updated cover letter |
| `Paper1_JAC/submission/suggested_reviewers.md` | Light updates |
| `Paper1_JAC/submission/build_auxiliary_docx.py` | Inline-content updates for cover letter |
| `Paper1_JAC/manuscript_paper1_JAC_v2_clean.docx` | Clean v2 .docx |
| `Paper1_JAC/manuscript_paper1_JAC_v2_trackchanges.docx` | Track-changes vs PJAM-MI v1 |
| `Paper1_JAC/submission/cover_letter_v2.docx` | v2 cover letter |
| `Paper1_JAC/submission/v2_change_summary.md` | One-page cover note for co-authors |

---

## Reference: spec sections used per task

| Task | Spec sections |
|---|---|
| 1 | (Plan-setup; commits the spec + plan) |
| 2 | §4.1, §5.2 |
| 3 | §4.1, §5.2 |
| 4 | §5.3, §3.3 (decisions: high-Q-tail scaling) |
| 5 | §5.4–5.6 |
| 6 | §5.7 (decision gate) |
| 7 | §6.1 |
| 8 | §6.2 |
| 9 | §6.3 |
| 10 | §6.4 |
| 11 | §6.5, §10 items 2–5 |
| 12 | §6.6 |
| 13 | §6.7 |
| 14 | §6.8 (decision gate) |
| 15–24 | §7 |
| 25–40 | §8 |
| 41–46 | §9 |

---

## Phase 1 — Methods Validation

### Task 1: Confirm beamtime wavelength, set up `n33_expansion/` directory, commit this plan

**Files:**
- Create: `Llunr/n33_expansion/README.md`
- Create directories: `Llunr/n33_expansion/{scripts,tests,data,phase1_validation,results}/`
- Create: `Llunr/n33_expansion/wavelength.txt` (one-line record of the confirmed Apr 2026 I11 MAC λ)
- Commit: `xrd_profile/docs/superpowers/plans/2026-05-15-jac-n33-expansion.md`

[PROCEDURAL + AUTHOR TASK]

- [ ] **Step 0 (AUTHOR TASK): Confirm Apr 2026 I11 MAC wavelength**

The I11 MAC wavelength was recorded by the Diamond beamline team at the time
of analysis: λ = 0.824883 Å, identical for the 2018 ee17803-1 and Apr 2026
beamtimes (author-confirmed 2026-05-15). Write this value into
`Llunr/n33_expansion/wavelength.txt`. **Provenance note:** the
`2026-05-11-reference-phase-sensitivity-design.md` spec and the existing
Rietveld project files (`HED_XRD_Shock/Rietveld/Tirhert/iter12/`,
`phase_manifest.yaml`, etc.) cite λ = 0.826517 Å, which is an incorrect value
whose source of error is pending retrace. The manuscript Table 1 footnote
(currently `λ = 0.8265 Å` at `Paper1_JAC/manuscript_paper1_JAC.md` line 237)
and the §2.1.1 example code (line 61) also show the incorrect 0.8265; Phase 4
of this plan must fix both. The Rietveld retrace and any consequences for
existing refinements are out of scope for this plan.

```bash
mkdir -p "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/n33_expansion"
echo "0.824883  # I11 MAC wavelength recorded by Diamond beamline team at analysis. Identical for the 2018 ee17803-1 and Apr 2026 (Izawa_achondrites_20260424) beamtimes; author-confirmed 2026-05-15. NB: HED_XRD_Shock/Rietveld/ project files and 2026-05-11-reference-phase-sensitivity-design.md cite 0.826517 — incorrect value, source of error pending retrace." > \
    "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/n33_expansion/wavelength.txt"
cat "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/n33_expansion/wavelength.txt"
```

- [ ] **Step 1: Create directory tree**

```bash
mkdir -p "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/n33_expansion/scripts"
mkdir -p "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/n33_expansion/tests"
mkdir -p "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/n33_expansion/data/i11_2026Apr_reduced"
mkdir -p "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/n33_expansion/data/i11_corrected"
mkdir -p "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/n33_expansion/phase1_validation"
mkdir -p "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/n33_expansion/results/per_sample"
mkdir -p "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/n33_expansion/results/per_sample_pdf_refits"
```

- [ ] **Step 2: Verify tree**

```bash
ls -la "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/n33_expansion/"
```

Expected: `scripts/`, `tests/`, `data/`, `phase1_validation/`, `results/` all visible.

- [ ] **Step 3: Write `README.md`**

File: `Llunr/n33_expansion/README.md`

```markdown
# n33_expansion — JAC N=33 manuscript revision analysis

Standalone analysis directory for integrating the 24 Apr 2026 DLS I11
data into the JAC manuscript. Lives outside git (Llunr/ is not
git-tracked).

## Goal

Produce `Paper1_JAC/survey_results_32samples.csv` and supporting
figures for a v2 manuscript with capillary subtraction, 4 new
samples, and 2 lab→synchrotron upgrades.

## Entry points

- Phase 1 (validation): `python scripts/phase1_validation.py`
- Phase 2 (full processing): `python scripts/process_all_samples.py`
- Phase 2 (CSV build): `python scripts/build_survey_csv.py`

## Prerequisites

- xrd_profile v0.3.x importable
- numpy, scipy, matplotlib, pymatgen, pytest, python-docx
- Anaconda Python at `C:\Users\Matthew Izawa\anaconda3\python.exe`

## Source data

Raw DLS data: `Llunr/HED_XRD_Shock/New_DLS_20260513/`
- 1437211–1437215, 1437218: 6 sample scans
- 1437220: blank capillary

Existing I11 data: per `Paper1_JAC/survey_results_28samples.csv`.

## Outputs

See `xrd_profile/docs/superpowers/specs/2026-05-15-jac-n33-expansion-design.md`
§13 for the full file inventory.
```

- [ ] **Step 4: Verify the plan is already committed**

The plan file was committed by the plan author at the end of `writing-plans`.
Confirm:

```bash
cd "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/xrd_profile"
git log --oneline -1 -- docs/superpowers/plans/2026-05-15-jac-n33-expansion.md
```

Expected: a single recent commit line with subject
"Add JAC N=32 expansion plan: …". If absent, run

```bash
git add docs/superpowers/plans/2026-05-15-jac-n33-expansion.md
git commit -m "Add JAC N=32 expansion plan: 5-phase validate-first integration"
```

---

### Task 2: Implement `.dat` → `.xy` reducer

**Files:**
- Create: `Llunr/n33_expansion/scripts/reduce_dat_to_xy.py`
- Create: `Llunr/n33_expansion/tests/test_reduce_dat.py`

[TDD]

The DLS MAC `.dat` file has a header block delimited by `&DLS` ... `&END`,
followed by a tab-separated table `tth\tcounts\terror` with ~150 020 rows.
The reducer parses the header for metadata (RunNumber, ScanTime, Date,
Wavelength), reads the numeric table, and emits a 2-column `.xy` file
(2θ, counts) matching the existing I11 file convention. The wavelength
defaults to 0.824883 Å (from Task 1 Step 0; override at call site if
`wavelength.txt` records a different value) but is overridable per call.

- [ ] **Step 1: Write the failing test**

File: `Llunr/n33_expansion/tests/test_reduce_dat.py`

```python
"""Tests for DLS .dat → .xy reducer."""
from pathlib import Path

import numpy as np
import pytest

from n33_expansion.scripts.reduce_dat_to_xy import (
    parse_dls_header,
    read_dat,
    write_xy,
)

NEW_DLS_DIR = Path(
    "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/HED_XRD_Shock/New_DLS_20260513"
).expanduser()


def test_parse_header_extracts_run_number():
    header_bytes = (
        b"&DLS\n"
        b"CarouselNo=Not Set\n"
        b"RunNumber=1437211\n"
        b"ScanTime=1800.0\n"
        b"Date=Apr 24, 2026\n"
        b"Wavelength=Not Set\n"
        b"&END\n"
    )
    meta = parse_dls_header(header_bytes.decode("latin1"))
    assert meta["RunNumber"] == "1437211"
    assert meta["ScanTime"] == "1800.0"


def test_read_dat_returns_two_arrays_same_length():
    dat = NEW_DLS_DIR / "1437211-mac-summed.dat"
    tth, counts = read_dat(dat)
    assert isinstance(tth, np.ndarray)
    assert isinstance(counts, np.ndarray)
    assert tth.shape == counts.shape
    assert tth.size > 100_000
    assert np.all(np.diff(tth) > 0)  # monotonic increasing


def test_write_xy_round_trip(tmp_path):
    tth = np.array([0.001, 0.002, 0.003])
    counts = np.array([11.2, 9.4, 9.5])
    out = tmp_path / "test.xy"
    write_xy(out, tth, counts)
    assert out.exists()
    loaded = np.loadtxt(out)
    assert loaded.shape == (3, 2)
    np.testing.assert_allclose(loaded[:, 0], tth)
    np.testing.assert_allclose(loaded[:, 1], counts)
```

- [ ] **Step 2: Run test to verify it fails**

```bash
"/c/Users/Matthew Izawa/anaconda3/python.exe" -m pytest "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/n33_expansion/tests/test_reduce_dat.py" -v
```

Expected: FAIL with `ModuleNotFoundError: No module named 'n33_expansion'` or `cannot import name`.

- [ ] **Step 3: Write minimal implementation**

File: `Llunr/n33_expansion/scripts/reduce_dat_to_xy.py`

```python
"""Reduce DLS I11 MAC .dat files to xrd_profile-readable .xy.

The .dat file format:
  &DLS
  key=value (header lines)
  ...
  &END
  tth\tcounts\terror
  <numeric rows, tab-separated>

This reducer parses the header into a dict, reads the numeric table,
and writes a 2-column .xy (2θ in degrees, counts).
"""
from __future__ import annotations

from pathlib import Path

import numpy as np


HEADER_START = "&DLS"
HEADER_END = "&END"


def parse_dls_header(text: str) -> dict[str, str]:
    """Parse the &DLS ... &END block into a metadata dict."""
    meta: dict[str, str] = {}
    in_header = False
    for line in text.splitlines():
        if line.strip().startswith(HEADER_START):
            in_header = True
            continue
        if line.strip().startswith(HEADER_END):
            break
        if in_header and "=" in line:
            key, _, value = line.partition("=")
            meta[key.strip()] = value.strip()
    return meta


def read_dat(path: Path | str) -> tuple[np.ndarray, np.ndarray]:
    """Read a DLS MAC .dat and return (2θ, counts) arrays."""
    path = Path(path)
    raw = path.read_bytes().decode("latin1")
    # Header parsing (kept for traceability; not returned here)
    _ = parse_dls_header(raw)
    # Locate the start of the numeric block
    end_marker = f"\n{HEADER_END}\n"
    end_idx = raw.find(end_marker)
    if end_idx < 0:
        raise ValueError(f"{path}: &END marker not found")
    numeric_block = raw[end_idx + len(end_marker):]
    # Skip the column-header line "tth\tcounts\terror"
    first_nl = numeric_block.find("\n")
    table_text = numeric_block[first_nl + 1:]
    data = np.loadtxt(table_text.splitlines())
    # Columns: tth, counts, error
    return data[:, 0], data[:, 1]


def write_xy(path: Path | str, tth: np.ndarray, counts: np.ndarray) -> None:
    """Write a 2-column .xy file (2θ in degrees, counts)."""
    path = Path(path)
    arr = np.column_stack([tth, counts])
    np.savetxt(path, arr, fmt="%.6f\t%.6f")


def reduce_one(dat_path: Path | str, out_dir: Path | str) -> Path:
    """Convert one .dat to .xy in `out_dir`. Returns output path."""
    dat_path = Path(dat_path)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    tth, counts = read_dat(dat_path)
    out_path = out_dir / (dat_path.stem + ".xy")
    write_xy(out_path, tth, counts)
    return out_path


if __name__ == "__main__":
    import sys
    for dat in sys.argv[1:]:
        out = reduce_one(dat, "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/n33_expansion/data/i11_2026Apr_reduced")
        print(f"{dat} → {out}")
```

- [ ] **Step 4: Add `__init__.py` files for importability**

```bash
touch "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/n33_expansion/__init__.py"
touch "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/n33_expansion/scripts/__init__.py"
touch "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/n33_expansion/tests/__init__.py"
```

- [ ] **Step 5: Run test to verify it passes**

```bash
cd "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr"
"/c/Users/Matthew Izawa/anaconda3/python.exe" -m pytest n33_expansion/tests/test_reduce_dat.py -v
```

Expected: 3 passed.

---

### Task 3: Reduce the 7 new `.dat` files and sanity-check capillary

**Files:**
- Run: `Llunr/n33_expansion/scripts/reduce_dat_to_xy.py` (as script)
- Verify: 7 `.xy` files in `Llunr/n33_expansion/data/i11_2026Apr_reduced/`
- Verify: capillary `.xy` shows broad amorphous hump only, no Bragg peaks

[PROCEDURAL]

- [ ] **Step 1: Run reducer over all 7 .dat files**

```bash
cd "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr"
"/c/Users/Matthew Izawa/anaconda3/python.exe" n33_expansion/scripts/reduce_dat_to_xy.py \
  HED_XRD_Shock/New_DLS_20260513/1437211-mac-summed.dat \
  HED_XRD_Shock/New_DLS_20260513/1437212-mac-summed.dat \
  HED_XRD_Shock/New_DLS_20260513/1437213-mac-summed.dat \
  HED_XRD_Shock/New_DLS_20260513/1437214-mac-summed.dat \
  HED_XRD_Shock/New_DLS_20260513/1437215-mac-summed.dat \
  HED_XRD_Shock/New_DLS_20260513/1437218-mac-summed.dat \
  HED_XRD_Shock/New_DLS_20260513/1437220-mac-summed.dat
```

Expected output: 7 lines, each `<dat> → <xy>`.

- [ ] **Step 2: Verify output files**

```bash
ls -la "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/n33_expansion/data/i11_2026Apr_reduced/"
```

Expected: 7 `.xy` files, each ~5 MB (150 020 rows of `2θ\tcounts`).

- [ ] **Step 3: Plot capillary alone for sanity check**

File: `Llunr/n33_expansion/scripts/plot_capillary.py`

```python
"""One-off: plot capillary .xy and save as PNG."""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


CAP = Path(
    "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/n33_expansion/data/i11_2026Apr_reduced/1437220-mac-summed.xy"
)
OUT = Path(
    "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/n33_expansion/phase1_validation/capillary_alone.png"
)


def main() -> None:
    data = np.loadtxt(CAP)
    fig, ax = plt.subplots(figsize=(8, 4), constrained_layout=True)
    ax.plot(data[:, 0], data[:, 1], lw=0.6, color="k")
    ax.set_xlabel("2θ (°)")
    ax.set_ylabel("Counts")
    ax.set_title("Blank capillary 1437220 — sanity check (broad hump only)")
    fig.savefig(OUT, dpi=150)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
```

Run:

```bash
"/c/Users/Matthew Izawa/anaconda3/python.exe" "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/n33_expansion/scripts/plot_capillary.py"
```

- [ ] **Step 4: Inspect `phase1_validation/capillary_alone.png`**

Expected: a single broad hump around 2θ ≈ 10–25° (amorphous borosilicate) with no narrow Bragg peaks. If sharp peaks are visible, the file is not what we think; halt and investigate.

---

### Task 4: Implement capillary subtraction with high-Q-tail scaling

**Files:**
- Create: `Llunr/n33_expansion/scripts/capillary_subtract.py`
- Create: `Llunr/n33_expansion/tests/test_capillary_subtract.py`

[TDD]

The subtraction module:
1. Converts 2θ to Q via `Q = (4π/λ) sin(θ)` (λ from caller).
2. Interpolates capillary onto the sample's Q-grid (they share a 2θ grid by construction, but the function should not assume that).
3. Fits a scalar `s` minimising `Σ_(Q ∈ [Q_lo, Q_hi]) (I_sample(Q) − s · I_cap(Q))²` via closed-form least-squares: `s = Σ I_s I_c / Σ I_c²` over the window.
4. Returns `I_corrected = I_sample − s · I_cap` and the scale factor `s`.

Default window: `Q_lo = 12.0`, `Q_hi = 14.0` Å⁻¹.

- [ ] **Step 1: Write the failing test**

File: `Llunr/n33_expansion/tests/test_capillary_subtract.py`

```python
"""Tests for high-Q-tail capillary subtraction."""
from __future__ import annotations

import numpy as np
import pytest

from n33_expansion.scripts.capillary_subtract import (
    tth_to_q,
    fit_scale_factor,
    subtract_capillary,
)


WAVELENGTH = 0.824883  # Å, I11 MAC


def test_tth_to_q_round():
    # At 2θ = 90° and λ = 0.824883 Å, Q = (4π/λ) sin(45°) = (4π/0.824883)(√2/2)
    expected = (4 * np.pi / WAVELENGTH) * np.sin(np.deg2rad(45.0))
    assert tth_to_q(np.array([90.0]), WAVELENGTH)[0] == pytest.approx(expected, rel=1e-9)


def test_fit_scale_factor_recovers_known_scale():
    rng = np.random.default_rng(0)
    Q = np.linspace(0.5, 14.5, 5000)
    I_cap = 100.0 * np.exp(-((Q - 1.6) / 1.0) ** 2) + 5.0 + rng.normal(0, 0.01, Q.size)
    s_true = 0.73
    I_sample = s_true * I_cap.copy()
    I_sample[(Q > 3) & (Q < 8)] += 200.0 * np.exp(-((Q - 5) / 0.3) ** 2)  # add Bragg peaks far from the high-Q window
    s_hat = fit_scale_factor(Q, I_sample, Q, I_cap, q_lo=12.0, q_hi=14.0)
    assert s_hat == pytest.approx(s_true, rel=5e-3)


def test_subtract_capillary_returns_corrected_and_scale():
    Q = np.linspace(0.5, 14.5, 5000)
    I_cap = np.ones_like(Q)
    I_sample = 0.5 * np.ones_like(Q) + 1.0  # constant offset on the high-Q tail
    corrected, s = subtract_capillary(
        Q, I_sample, Q, I_cap, q_lo=12.0, q_hi=14.0
    )
    assert s == pytest.approx(1.5, rel=1e-9)
    assert corrected.shape == I_sample.shape
```

- [ ] **Step 2: Run test to verify it fails**

```bash
cd "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr"
"/c/Users/Matthew Izawa/anaconda3/python.exe" -m pytest n33_expansion/tests/test_capillary_subtract.py -v
```

Expected: FAIL with `ImportError` or `ModuleNotFoundError`.

- [ ] **Step 3: Write minimal implementation**

File: `Llunr/n33_expansion/scripts/capillary_subtract.py`

```python
"""Capillary background subtraction via high-Q-tail least-squares scaling."""
from __future__ import annotations

import numpy as np


def tth_to_q(tth_deg: np.ndarray, wavelength_A: float) -> np.ndarray:
    """Convert 2θ in degrees to Q in inverse Å. Q = (4π/λ) sin(θ)."""
    theta = np.deg2rad(tth_deg) / 2.0
    return (4.0 * np.pi / wavelength_A) * np.sin(theta)


def fit_scale_factor(
    q_sample: np.ndarray,
    i_sample: np.ndarray,
    q_cap: np.ndarray,
    i_cap: np.ndarray,
    *,
    q_lo: float = 12.0,
    q_hi: float = 14.0,
) -> float:
    """Find s minimising Σ_(Q∈[q_lo, q_hi]) (I_sample − s · I_cap)² .

    Capillary intensities are interpolated onto the sample Q-grid first.
    Closed-form solution: s = Σ I_s I_c / Σ I_c² over the window.
    """
    i_cap_on_sample = np.interp(q_sample, q_cap, i_cap)
    mask = (q_sample >= q_lo) & (q_sample <= q_hi)
    if mask.sum() < 10:
        raise ValueError(f"Window [{q_lo}, {q_hi}] contains too few points ({mask.sum()})")
    num = float(np.sum(i_sample[mask] * i_cap_on_sample[mask]))
    den = float(np.sum(i_cap_on_sample[mask] ** 2))
    if den == 0.0:
        raise ValueError("Capillary intensity is zero in the fit window")
    return num / den


def subtract_capillary(
    q_sample: np.ndarray,
    i_sample: np.ndarray,
    q_cap: np.ndarray,
    i_cap: np.ndarray,
    *,
    q_lo: float = 12.0,
    q_hi: float = 14.0,
) -> tuple[np.ndarray, float]:
    """Return (I_corrected, scale_factor)."""
    s = fit_scale_factor(q_sample, i_sample, q_cap, i_cap, q_lo=q_lo, q_hi=q_hi)
    i_cap_on_sample = np.interp(q_sample, q_cap, i_cap)
    return i_sample - s * i_cap_on_sample, s
```

- [ ] **Step 4: Run test to verify it passes**

```bash
"/c/Users/Matthew Izawa/anaconda3/python.exe" -m pytest n33_expansion/tests/test_capillary_subtract.py -v
```

Expected: 3 passed.

---

### Task 5: Generate Phase 1 validation diagnostic figure

**Files:**
- Create: `Llunr/n33_expansion/scripts/phase1_validation.py`
- Output: `Llunr/n33_expansion/phase1_validation/diagnostic.png`

[PROCEDURAL]

This is a one-shot driver: load capillary + 3 test eucrites (Tirhert S1, NWA 6477 S3, JaH 626 S6), apply subtraction, run xrd_profile PDF on each, plot a 4-panel diagnostic figure.

**Author input needed before running:** the I11 `.xy` paths for the 3 test eucrites. Use `Paper1_JAC/survey_results_28samples.csv` or inspect the existing I11 data directory.

- [ ] **Step 1: Identify existing I11 .xy paths for Tirhert, NWA 6477, JaH 626**

```bash
"/c/Users/Matthew Izawa/anaconda3/python.exe" -c "
import csv
with open('/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/Paper1_JAC/survey_results_28samples.csv') as f:
    r = csv.DictReader(f)
    for row in r:
        if row.get('sample_id') in ('Tirhert', 'NWA 6477', 'JaH 626'):
            print(row.get('sample_id'), '→', row.get('source_path', row))
"
```

Expected: 3 lines, one per sample, showing the source `.xy` path (or whatever field name the CSV uses for the input pattern).

If `source_path` is not a column, inspect the CSV header:

```bash
head -1 "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/Paper1_JAC/survey_results_28samples.csv"
```

and adapt the field name. The plan author has not seen the CSV columns directly; the implementing agent should record the actual paths in `Llunr/n33_expansion/phase1_validation/decision_note.md`.

- [ ] **Step 2: Write the driver**

File: `Llunr/n33_expansion/scripts/phase1_validation.py`

```python
"""Phase 1 driver: validate capillary subtraction on 3 test eucrites.

Outputs:
  Llunr/n33_expansion/phase1_validation/diagnostic.png  (4 panels)
  Llunr/n33_expansion/phase1_validation/test_fwhm.csv   (before/after FWHM table)
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from xrd_profile import XRDProfile

from n33_expansion.scripts.capillary_subtract import (
    tth_to_q,
    subtract_capillary,
)
from n33_expansion.scripts.reduce_dat_to_xy import write_xy


WAVELENGTH = 0.824883  # Å — confirm via Llunr/n33_expansion/wavelength.txt (Task 1 Step 0)
N33 = Path("/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/n33_expansion")
CAP_XY = N33 / "data/i11_2026Apr_reduced/1437220-mac-summed.xy"
OUT_DIR = N33 / "phase1_validation"

# === FILL IN from Step 1 ===
TEST_SAMPLES = {
    "Tirhert": {"path": "[AUTHOR INPUT NEEDED]", "shock": "S1"},
    "NWA 6477": {"path": "[AUTHOR INPUT NEEDED]", "shock": "S3"},
    "JaH 626": {"path": "[AUTHOR INPUT NEEDED]", "shock": "S6"},
}


def load_xy(path: Path) -> tuple[np.ndarray, np.ndarray]:
    data = np.loadtxt(path)
    return data[:, 0], data[:, 1]


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    tth_cap, I_cap = load_xy(CAP_XY)
    q_cap = tth_to_q(tth_cap, WAVELENGTH)

    fig, axes = plt.subplots(2, 2, figsize=(12, 9), constrained_layout=True)
    rows = []

    # Panel (a): capillary I(Q)
    ax = axes[0, 0]
    ax.plot(q_cap, I_cap, lw=0.6, color="k")
    ax.set_xlabel("Q (Å⁻¹)")
    ax.set_ylabel("Counts")
    ax.set_title("(a) Blank capillary")

    # Panel (b): sample I(Q) before vs after, low-Q emphasis
    ax = axes[0, 1]
    for name, info in TEST_SAMPLES.items():
        tth_s, I_s = load_xy(Path(info["path"]))
        q_s = tth_to_q(tth_s, WAVELENGTH)
        I_corr, s = subtract_capillary(q_s, I_s, q_cap, I_cap)
        ax.plot(q_s, I_s, lw=0.4, alpha=0.5, label=f"{name} raw")
        ax.plot(q_s, I_corr, lw=0.6, label=f"{name} −cap (s={s:.2f})")
    ax.set_xlim(0.5, 5.0)
    ax.set_xlabel("Q (Å⁻¹)")
    ax.set_ylabel("Counts")
    ax.set_title("(b) Sample I(Q): raw vs capillary-subtracted (low-Q)")
    ax.legend(fontsize=7)

    # Panels (c) and (d) below: PDF before/after and FWHM trend
    pdf_axis = axes[1, 0]
    fwhm_axis = axes[1, 1]
    pdf_axis.set_title("(c) PDF G(r) before (dashed) vs after (solid)")
    pdf_axis.set_xlabel("r (Å)")
    pdf_axis.set_ylabel("G(r)")
    fwhm_axis.set_title("(d) First-peak FWHM vs shock — old (○) vs new (●)")
    fwhm_axis.set_xlabel("Shock stage (literature)")
    fwhm_axis.set_ylabel("FWHM (Å)")

    for name, info in TEST_SAMPLES.items():
        tth_s, I_s = load_xy(Path(info["path"]))
        q_s = tth_to_q(tth_s, WAVELENGTH)
        I_corr, s = subtract_capillary(q_s, I_s, q_cap, I_cap)

        # Write corrected .xy
        corr_path = N33 / "data/i11_corrected" / Path(info["path"]).name.replace(".xy", "_cap.xy")
        corr_path.parent.mkdir(parents=True, exist_ok=True)
        write_xy(corr_path, tth_s, I_corr)

        # Run xrd_profile PDF on raw and corrected
        prof_raw = XRDProfile.from_file(str(Path(info["path"])), wavelength=WAVELENGTH)
        pdf_raw = prof_raw.pdf(q_max=14.6)
        prof_corr = XRDProfile.from_file(str(corr_path), wavelength=WAVELENGTH)
        pdf_corr = prof_corr.pdf(q_max=14.6)

        pdf_axis.plot(pdf_raw["r"], pdf_raw["G"], "--", alpha=0.5, label=f"{name} raw")
        pdf_axis.plot(pdf_corr["r"], pdf_corr["G"], "-", label=f"{name} −cap")

        fwhm_old = pdf_raw["first_peak_fwhm"]
        fwhm_new = pdf_corr["first_peak_fwhm"]
        rows.append({
            "sample": name, "shock": info["shock"],
            "fwhm_raw": fwhm_old, "fwhm_corr": fwhm_new, "scale": s,
        })

        fwhm_axis.plot([info["shock"]], [fwhm_old], "o", markerfacecolor="none", color="k")
        fwhm_axis.plot([info["shock"]], [fwhm_new], "o", color="C0")

    pdf_axis.set_xlim(0, 10)
    pdf_axis.legend(fontsize=7)

    fig.savefig(OUT_DIR / "diagnostic.png", dpi=180)
    print(f"wrote {OUT_DIR / 'diagnostic.png'}")

    import csv
    with open(OUT_DIR / "test_fwhm.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["sample", "shock", "fwhm_raw", "fwhm_corr", "scale"])
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print(f"wrote {OUT_DIR / 'test_fwhm.csv'}")


if __name__ == "__main__":
    main()
```

- [ ] **Step 3: Confirm wavelength and fill in TEST_SAMPLES**

[AUTHOR TASK] — confirm `WAVELENGTH = 0.824883` Å from Apr 2026 I11 MAC log book or proposal; fill in the 3 `.xy` paths from Step 1.

- [ ] **Step 4: Run the driver**

```bash
cd "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr"
"/c/Users/Matthew Izawa/anaconda3/python.exe" n33_expansion/scripts/phase1_validation.py
```

Expected: prints `wrote …/diagnostic.png` and `wrote …/test_fwhm.csv`.

- [ ] **Step 5: Open `diagnostic.png` for visual inspection**

The diagnostic figure has 4 panels; inspect each. Note any anomaly in `decision_note.md` (next task).

---

### Task 6: Phase 1 decision gate — AUTHOR REVIEW

**Files:**
- Create: `Llunr/n33_expansion/phase1_validation/decision_note.md`

[AUTHOR TASK]

- [ ] **Step 1: Author reviews `diagnostic.png` and `test_fwhm.csv`**

Decision criteria (from spec §5.7):
- GO: monotonic broadening trend preserved or strengthened; FWHM rank-order with shock stage intact (Tirhert FWHM < NWA 6477 FWHM < JaH 626 FWHM).
- PAUSE: trend disrupted.

- [ ] **Step 2: Write decision note**

File: `Llunr/n33_expansion/phase1_validation/decision_note.md`

Template:

```markdown
# Phase 1 validation decision

Date: YYYY-MM-DD
Author: M. R. M. Izawa

## Wavelength

Apr 2026 I11 MAC λ = ___ Å (source: ___).

## Capillary scale factors from the 3 test samples

| Sample | Shock | FWHM raw (Å) | FWHM corrected (Å) | Scale s |
|---|---|---|---|---|
| Tirhert | S1 | … | … | … |
| NWA 6477 | S3 | … | … | … |
| JaH 626 | S6 | … | … | … |

## Trend assessment

[GO/PAUSE]. Reasoning:

…

## Outstanding follow-ups

- [ ] …
```

- [ ] **Step 3: Halt if PAUSE; proceed to Phase 2 if GO**

If PAUSE: the implementing agent halts and asks the user how to proceed (spec §11 abort conditions). If GO: continue to Task 7.

---

## Phase 2 — Data Processing and Survey Expansion

### Task 7: Verify pre-expansion inventory against the canonical CSV

**Files:**
- Read: `Paper1_JAC/survey_results_28samples.csv`
- Output: a printed inventory table (no file written)

[PROCEDURAL]

The spec assumes 7 Lab Cu + 17 I11 + 4 Lab Co = 28 (NOTE: original spec said "29" — that was a prose miscount; the CSV has 28 rows; see decision_note.md), with the I11 breakdown of 11 HED + 6 non-HED. Verify this against the CSV before generating the 32-sample CSV.

- [ ] **Step 1: Read CSV column headers**

```bash
head -1 "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/Paper1_JAC/survey_results_28samples.csv"
```

Record the exact column names that identify `sample`, `instrument`, and `group` (e.g., HED / Lunar / OC / CC / Achondrite / Martian).

- [ ] **Step 2: Run inventory check**

```bash
"/c/Users/Matthew Izawa/anaconda3/python.exe" -c "
import csv
from collections import Counter
with open('/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/Paper1_JAC/survey_results_28samples.csv') as f:
    rows = list(csv.DictReader(f))
print('Total rows:', len(rows))
print('Instruments:', Counter(r.get('instrument') for r in rows))
print('Groups:', Counter(r.get('group') for r in rows))
for r in rows:
    if r.get('instrument', '').startswith(('Cu', 'Co')):
        print(r.get('sample_id'), r.get('instrument'), r.get('group'))
"
```

Replace `instrument` / `group` / `sample_id` with the actual column names from Step 1.

Expected: total 28; instruments break into 7 Lab Cu + 17 I11 + 4 Lab Co; the 4 Lab Co samples are NWA 5751, NWA 6013, Talampaya, Tatahouine.

- [ ] **Step 3: Resolve any discrepancy**

If totals don't match, halt and update the spec inventory before proceeding. The 32-sample CSV inherits whatever the 28-sample CSV says.

---

### Task 8: Apply capillary subtraction to all 23 I11 samples

**Files:**
- Create: `Llunr/n33_expansion/scripts/apply_capillary_all.py`
- Output: 24 `.xy` files in `Llunr/n33_expansion/data/i11_corrected/`

[TDD] (for the path-list logic — small but worth a test)

- [ ] **Step 1: Write the failing test**

File: append to `Llunr/n33_expansion/tests/test_capillary_subtract.py`

```python
def test_apply_capillary_writes_one_xy_per_input(tmp_path):
    """Smoke test that the driver iterates over all inputs."""
    from n33_expansion.scripts.apply_capillary_all import process_sample

    # Make a fake sample
    sample_in = tmp_path / "fake_sample.xy"
    cap_in = tmp_path / "fake_cap.xy"
    out_dir = tmp_path / "out"
    tth = np.linspace(0.1, 60.0, 200)
    np.savetxt(sample_in, np.column_stack([tth, tth * 0 + 100.0]))
    np.savetxt(cap_in, np.column_stack([tth, tth * 0 + 50.0]))

    out_path = process_sample(
        sample_in, cap_in, out_dir, wavelength=0.824883
    )
    assert out_path.exists()
    assert out_path.parent == out_dir
```

- [ ] **Step 2: Run to verify it fails**

```bash
"/c/Users/Matthew Izawa/anaconda3/python.exe" -m pytest n33_expansion/tests/test_capillary_subtract.py::test_apply_capillary_writes_one_xy_per_input -v
```

Expected: FAIL.

- [ ] **Step 3: Write minimal implementation**

File: `Llunr/n33_expansion/scripts/apply_capillary_all.py`

```python
"""Apply capillary subtraction to all I11 samples (new + existing)."""
from __future__ import annotations

from pathlib import Path

import numpy as np

from n33_expansion.scripts.capillary_subtract import (
    tth_to_q,
    subtract_capillary,
)
from n33_expansion.scripts.reduce_dat_to_xy import write_xy


def process_sample(
    sample_xy: Path,
    cap_xy: Path,
    out_dir: Path,
    *,
    wavelength: float,
    q_lo: float = 12.0,
    q_hi: float = 14.0,
) -> Path:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    s = np.loadtxt(sample_xy)
    c = np.loadtxt(cap_xy)
    q_s = tth_to_q(s[:, 0], wavelength)
    q_c = tth_to_q(c[:, 0], wavelength)
    I_corr, scale = subtract_capillary(
        q_s, s[:, 1], q_c, c[:, 1], q_lo=q_lo, q_hi=q_hi
    )
    out_path = out_dir / (Path(sample_xy).stem + "_cap.xy")
    write_xy(out_path, s[:, 0], I_corr)
    return out_path


def main() -> None:
    """Iterate over the 23 I11 samples (6 new + 17 existing).

    The 18 existing paths come from `survey_results_28samples.csv`;
    the 6 new from `i11_2026Apr_reduced/`. See the inventory in
    `decision_note.md`.
    """
    import csv

    N33 = Path("/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/n33_expansion")
    PAPER = Path("/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/Paper1_JAC")
    CAP = N33 / "data/i11_2026Apr_reduced/1437220-mac-summed.xy"
    OUT = N33 / "data/i11_corrected"
    WAVELENGTH = 0.824883

    # 6 new (skip the capillary itself)
    new_dir = N33 / "data/i11_2026Apr_reduced"
    new_samples = [
        p for p in sorted(new_dir.glob("*.xy"))
        if "1437220" not in p.name
    ]

    # 18 existing — read from CSV
    existing = []
    with open(PAPER / "survey_results_28samples.csv") as f:
        for row in csv.DictReader(f):
            if row.get("instrument", "").strip() == "I11":  # field name placeholder; adjust to actual column name from Task 7 Step 1
                existing.append(Path(row["source_path"]))  # field name placeholder; adjust to actual column name from Task 7 Step 1

    print(f"Processing {len(new_samples)} new + {len(existing)} existing = "
          f"{len(new_samples) + len(existing)} samples")
    for s in list(new_samples) + existing:
        out = process_sample(s, CAP, OUT, wavelength=WAVELENGTH)
        print(f"  {s.name} → {out.name}")


if __name__ == "__main__":
    main()
```

- [ ] **Step 4: Run test to verify it passes**

```bash
"/c/Users/Matthew Izawa/anaconda3/python.exe" -m pytest n33_expansion/tests/test_capillary_subtract.py -v
```

Expected: 4 passed (3 existing + 1 new).

- [ ] **Step 5: Adapt the CSV column names in `apply_capillary_all.main()`**

Replace `row.get("instrument")` and `row["source_path"]` with the actual field names from Task 7 Step 1. The plan author did not have direct CSV access; the implementing agent fills these in.

- [ ] **Step 6: Run the full pipeline**

```bash
cd "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr"
"/c/Users/Matthew Izawa/anaconda3/python.exe" n33_expansion/scripts/apply_capillary_all.py
```

Expected: prints `Processing 6 new + 17 existing = 23 samples`, then 23 lines `<in> → <out>`.

- [ ] **Step 7: Verify output count**

```bash
ls "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/n33_expansion/data/i11_corrected/" | wc -l
```

Expected: 24.

---

### Task 9: Run `xrd_profile` on the 6 new analyses

**Files:**
- Create: `Llunr/n33_expansion/scripts/process_new_samples.py`
- Output: 6 JSON files in `Llunr/n33_expansion/results/per_sample/`

[PROCEDURAL]

For each of the 6 new I11 samples (4 new + 2 upgrades), run the standard `XRDProfile.run_all` with the appropriate reference phase (anorthite for plagioclase-bearing eucrites/howardites, enstatite/pigeonite for pyroxene-dominated, forsterite for olivine-rich samples).

- [ ] **Step 1: Map each new sample to its primary reference phase**

| Scan | Sample | Type | Primary reference phase | Notes |
|---|---|---|---|---|
| 1437211 | NWA 7465 | Monomict eucrite | anorthite | also pigeonite secondary |
| 1437212 | NWA 5478 | Shocked polymict eucrite | anorthite | also pigeonite secondary |
| 1437213 | NWA 5751 | Howardite | anorthite | pigeonite, augite secondary |
| 1437214 | NWA 7831 | Diogenite | enstatite | minor olivine |
| 1437215 | NWA 6693 | Ungrouped achondrite | forsterite | Fe-rich olivine dominant; pyroxene secondary |
| 1437218 | Talampaya | Cumulate eucrite | anorthite | also pigeonite secondary |

- [ ] **Step 2: Write the driver**

File: `Llunr/n33_expansion/scripts/process_new_samples.py`

```python
"""Run xrd_profile.run_all on the 6 new I11 sample patterns."""
from __future__ import annotations

import json
from pathlib import Path

from xrd_profile import Phase, XRDProfile


N33 = Path("/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/n33_expansion")
PAPER = Path("/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/Paper1_JAC")
CIF_DIR = Path("/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/xrd_profile/examples/cifs")
WAVELENGTH = 0.824883

PRIMARY_REFS = {
    "1437211-mac-summed_cap.xy": ("NWA 7465", "Anorthite.cif"),
    "1437212-mac-summed_cap.xy": ("NWA 5478", "Anorthite.cif"),
    "1437213-mac-summed_cap.xy": ("NWA 5751", "Anorthite.cif"),
    "1437214-mac-summed_cap.xy": ("NWA 7831", "Enstatite.cif"),
    "1437215-mac-summed_cap.xy": ("NWA 6693", "Forsterite.cif"),
    "1437218-mac-summed_cap.xy": ("Talampaya", "Anorthite.cif"),
}


def main() -> None:
    out = N33 / "results/per_sample"
    out.mkdir(parents=True, exist_ok=True)
    corrected = N33 / "data/i11_corrected"

    for fname, (sample, cif_name) in PRIMARY_REFS.items():
        xy = corrected / fname
        cif = CIF_DIR / cif_name
        if not xy.exists():
            raise FileNotFoundError(xy)
        if not cif.exists():
            raise FileNotFoundError(cif)
        print(f"Running {sample} (ref={cif_name}, file={xy.name})")
        phase = Phase.from_cif(str(cif))
        prof = XRDProfile.from_file(str(xy), wavelength=WAVELENGTH)
        result = prof.run_all(phase=phase, q_max=14.6)
        # Serialise (xrd_profile result objects expose .to_dict() per package convention)
        with open(out / f"{sample.replace(' ', '_')}.json", "w") as f:
            json.dump(result.to_dict(), f, indent=2, default=str)
        print(f"  → wrote {out / sample.replace(' ', '_')}.json")


if __name__ == "__main__":
    main()
```

- [ ] **Step 3: Verify CIF availability**

```bash
ls "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/xrd_profile/examples/cifs/" | grep -E "Anorthite|Enstatite|Forsterite"
```

Expected: at least `Anorthite.cif`, `Enstatite.cif`, `Forsterite.cif` present. If absent, locate the project's CIF library and update `CIF_DIR`.

- [ ] **Step 4: Run the driver**

```bash
cd "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr"
"/c/Users/Matthew Izawa/anaconda3/python.exe" n33_expansion/scripts/process_new_samples.py
```

Expected: 6 lines `Running …` + 6 lines `→ wrote …`. Total 6 JSON files in `results/per_sample/`.

- [ ] **Step 5: Verify per-sample JSON content**

```bash
"/c/Users/Matthew Izawa/anaconda3/python.exe" -c "
import json
from pathlib import Path
p = Path('/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/n33_expansion/results/per_sample')
for f in sorted(p.glob('*.json')):
    d = json.load(open(f))
    keys = ', '.join(list(d.keys())[:6])
    print(f.stem, '→ keys:', keys, '...')
"
```

Expected: 6 lines, each showing top-level result keys (W-A, PDF, W-H, etc.).

---

### Task 10: Re-run PDF first-peak fit on the 17 existing I11 samples

**Files:**
- Create: `Llunr/n33_expansion/scripts/refit_existing_pdfs.py`
- Output: 17 JSON files in `Llunr/n33_expansion/results/per_sample_pdf_refits/`

[PROCEDURAL]

Same logic as Task 9, but PDF-only (skip W-A and W-H) and using capillary-subtracted `.xy` of existing samples.

- [ ] **Step 1: Write the driver**

File: `Llunr/n33_expansion/scripts/refit_existing_pdfs.py`

```python
"""Re-run PDF first-peak fit on the 17 existing I11 samples with
capillary-subtracted intensities."""
from __future__ import annotations

import csv
import json
from pathlib import Path

from xrd_profile import XRDProfile


N33 = Path("/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/n33_expansion")
PAPER = Path("/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/Paper1_JAC")
WAVELENGTH = 0.824883


def main() -> None:
    out = N33 / "results/per_sample_pdf_refits"
    out.mkdir(parents=True, exist_ok=True)
    corrected = N33 / "data/i11_corrected"

    # Read existing I11 sample list from the 28-sample CSV
    with open(PAPER / "survey_results_28samples.csv") as f:
        rows = [r for r in csv.DictReader(f) if r.get("instrument", "").strip() == "I11"]

    print(f"Refitting PDFs on {len(rows)} existing I11 samples")
    for row in rows:
        sample = row["sample_id"]  # field name placeholder; adjust to actual column name from Task 7 Step 1
        # The corrected file shares the existing sample's stem with `_cap` suffix.
        # Map sample → file via the CSV's source_path field.
        src = Path(row["source_path"])
        corr = corrected / (src.stem + "_cap.xy")
        if not corr.exists():
            raise FileNotFoundError(corr)
        prof = XRDProfile.from_file(str(corr), wavelength=WAVELENGTH)
        pdf_result = prof.pdf(q_max=14.6)
        with open(out / f"{sample.replace(' ', '_')}.json", "w") as f:
            json.dump(pdf_result, f, indent=2, default=str)
        print(f"  {sample} → first_peak_fwhm = {pdf_result.get('first_peak_fwhm')}")


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Adapt CSV field names**

As in Task 8 Step 5, replace `sample_id` / `source_path` / `instrument` with the actual CSV field names.

- [ ] **Step 3: Run the driver**

```bash
cd "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr"
"/c/Users/Matthew Izawa/anaconda3/python.exe" n33_expansion/scripts/refit_existing_pdfs.py
```

Expected: 18 lines, each `<sample> → first_peak_fwhm = <value>`.

- [ ] **Step 4: Verify count**

```bash
ls "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/n33_expansion/results/per_sample_pdf_refits/" | wc -l
```

Expected: 18.

---

### Task 11: Sample metadata research — AUTHOR TASK

**Files:**
- Update: `Llunr/n33_expansion/results/new_samples_metadata.md` (create)

[AUTHOR TASK]

For each of the 4 new samples, look up Meteoritical Bulletin entry plus any classification or shock-petrography reference.

- [ ] **Step 1: Author fills in the metadata table**

File: `Llunr/n33_expansion/results/new_samples_metadata.md`

```markdown
# Metadata for 4 new samples

| Sample | Group | Type | Major minerals | Lit. shock stage | Primary reference | DOI/URL |
|---|---|---|---|---|---|---|
| NWA 7465 | HED | Monomict eucrite | Anorthite + pigeonite + augite | … | … | … |
| NWA 5478 | HED | Shocked polymict eucrite | Anorthite + pigeonite + augite | … (S5? S6?) | … | … |
| NWA 7831 | HED | Diogenite | Orthopyroxene (+ minor olivine?) | … | … | … |
| NWA 6693 | Ungrouped achondrite (8th group) | CC-associated impact melt, paired with NWA 6704 | Fe-rich olivine + Fe-rich pyroxene + awaruite (Ni₃Fe, trace) | low (S1-S2?) | Warren et al. … | … |
```

- [ ] **Step 2: Note any `[AUTHOR INPUT NEEDED]` items still outstanding**

`[AUTHOR INPUT NEEDED]` placeholders are acceptable in this file and propagate into Table 1 / Table 2 for v2 co-author circulation; they must be resolved before final JAC submission.

---

### Task 12: Build `survey_results_32samples.csv` via aggregation script

**Files:**
- Create: `Llunr/n33_expansion/scripts/build_survey_csv.py`
- Create: `Llunr/n33_expansion/tests/test_build_survey_csv.py`
- Output: `Paper1_JAC/survey_results_32samples.csv`

[TDD]

The aggregator:
1. Loads the 28-sample CSV as a baseline.
2. Drops the 2 Lab Co rows for NWA 5751 and Talampaya.
3. Updates the PDF FWHM column for 17 existing I11 rows from `per_sample_pdf_refits/`.
4. Adds 4 new rows (NWA 7465, NWA 5478, NWA 7831, NWA 6693) from `per_sample/`.
5. Adds 2 new rows (NWA 5751, Talampaya) from `per_sample/` (their upgrade entries).
6. Writes 32 rows to `survey_results_32samples.csv`.

- [ ] **Step 1: Write the failing test**

File: `Llunr/n33_expansion/tests/test_build_survey_csv.py`

```python
"""Test the aggregation logic for the 32-sample CSV build."""
from __future__ import annotations

import csv

import pytest

from n33_expansion.scripts.build_survey_csv import (
    drop_lab_co_upgrades,
    update_pdf_fwhm,
    add_new_rows,
)


def test_drop_lab_co_upgrades_removes_two_rows():
    rows = [
        {"sample_id": "NWA 5751", "instrument": "Lab Co"},
        {"sample_id": "Talampaya", "instrument": "Lab Co"},
        {"sample_id": "Tirhert", "instrument": "I11"},
        {"sample_id": "Tatahouine", "instrument": "Lab Co"},
    ]
    result = drop_lab_co_upgrades(rows)
    ids = [r["sample_id"] for r in result]
    assert "NWA 5751" not in ids
    assert "Talampaya" not in ids
    assert "Tirhert" in ids
    assert "Tatahouine" in ids


def test_update_pdf_fwhm_replaces_for_named_samples():
    rows = [
        {"sample_id": "Tirhert", "instrument": "I11", "pdf_fwhm": "0.468"},
        {"sample_id": "Bereba", "instrument": "I11", "pdf_fwhm": "0.543"},
    ]
    refits = {"Tirhert": 0.520, "Bereba": 0.560}
    result = update_pdf_fwhm(rows, refits)
    assert float(result[0]["pdf_fwhm"]) == pytest.approx(0.520)
    assert float(result[1]["pdf_fwhm"]) == pytest.approx(0.560)


def test_add_new_rows_appends_correct_count():
    rows = [{"sample_id": "Tirhert"}]
    new = [
        {"sample_id": "NWA 7465"},
        {"sample_id": "NWA 5478"},
        {"sample_id": "NWA 7831"},
        {"sample_id": "NWA 6693"},
        {"sample_id": "NWA 5751"},
        {"sample_id": "Talampaya"},
    ]
    result = add_new_rows(rows, new)
    assert len(result) == 7
    assert result[-1]["sample_id"] == "Talampaya"
```

- [ ] **Step 2: Run test to verify it fails**

```bash
"/c/Users/Matthew Izawa/anaconda3/python.exe" -m pytest n33_expansion/tests/test_build_survey_csv.py -v
```

Expected: FAIL with `ImportError`.

- [ ] **Step 3: Write minimal implementation**

File: `Llunr/n33_expansion/scripts/build_survey_csv.py`

```python
"""Build survey_results_32samples.csv from the 28-sample baseline + new analyses."""
from __future__ import annotations

import csv
import json
from pathlib import Path


N33 = Path("/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/n33_expansion")
PAPER = Path("/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/Paper1_JAC")

UPGRADE_SAMPLES = {"NWA 5751", "Talampaya"}
NEW_SAMPLES = ["NWA 7465", "NWA 5478", "NWA 7831", "NWA 6693", "NWA 5751", "Talampaya"]


def drop_lab_co_upgrades(rows: list[dict]) -> list[dict]:
    """Remove the Lab Co rows for NWA 5751 and Talampaya."""
    return [
        r for r in rows
        if not (r["sample_id"] in UPGRADE_SAMPLES and r["instrument"] == "Lab Co")
    ]


def update_pdf_fwhm(rows: list[dict], refits: dict[str, float]) -> list[dict]:
    """Update the pdf_fwhm column for existing I11 samples that have a refit."""
    out = []
    for r in rows:
        if r["sample_id"] in refits and r["instrument"] == "I11":
            r = dict(r)
            r["pdf_fwhm"] = f"{refits[r['sample_id']]:.3f}"
        out.append(r)
    return out


def add_new_rows(rows: list[dict], new_rows: list[dict]) -> list[dict]:
    return rows + new_rows


def load_refits() -> dict[str, float]:
    refit_dir = N33 / "results/per_sample_pdf_refits"
    refits = {}
    for f in refit_dir.glob("*.json"):
        sample = f.stem.replace("_", " ")
        d = json.load(open(f))
        refits[sample] = d["first_peak_fwhm"]
    return refits


def build_new_row(sample: str, json_path: Path, template: dict) -> dict:
    d = json.load(open(json_path))
    row = dict(template)  # inherit column structure
    row["sample_id"] = sample
    row["instrument"] = "I11"
    row["q_max"] = "14.6"
    # Field names are subject to actual CSV columns — adjust per Task 7
    row["wa_median_column_length"] = d.get("warren_averbach", {}).get("median_column_length")
    row["pdf_fwhm"] = d.get("pdf", {}).get("first_peak_fwhm")
    return row


def main() -> None:
    src = PAPER / "survey_results_28samples.csv"
    dst = PAPER / "survey_results_32samples.csv"

    with open(src) as f:
        reader = csv.DictReader(f)
        fields = reader.fieldnames
        rows = list(reader)
    assert len(rows) == 28, f"expected 28 baseline rows, got {len(rows)}"

    template = dict.fromkeys(fields, "")  # for new rows

    rows = drop_lab_co_upgrades(rows)
    assert len(rows) == 27, f"after upgrade-drop expected 27, got {len(rows)}"

    refits = load_refits()
    rows = update_pdf_fwhm(rows, refits)

    new_rows = []
    per_sample = N33 / "results/per_sample"
    for s in NEW_SAMPLES:
        new_rows.append(build_new_row(s, per_sample / f"{s.replace(' ', '_')}.json", template))
    rows = add_new_rows(rows, new_rows)
    assert len(rows) == 32, f"final expected 32, got {len(rows)}"

    with open(dst, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print(f"wrote {dst} with {len(rows)} rows")


if __name__ == "__main__":
    main()
```

- [ ] **Step 4: Adapt field names per Task 7**

Replace `instrument`, `sample_id`, `pdf_fwhm`, `wa_median_column_length`, `q_max` with the actual CSV column names.

- [ ] **Step 5: Run test to verify it passes**

```bash
"/c/Users/Matthew Izawa/anaconda3/python.exe" -m pytest n33_expansion/tests/test_build_survey_csv.py -v
```

Expected: 3 passed.

- [ ] **Step 6: Run the aggregator**

```bash
cd "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr"
"/c/Users/Matthew Izawa/anaconda3/python.exe" n33_expansion/scripts/build_survey_csv.py
```

Expected: `wrote …/survey_results_32samples.csv with 32 rows`. Sanity-check via:

```bash
wc -l "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/Paper1_JAC/survey_results_32samples.csv"
```

Expected: 33 (32 rows + 1 header).

---

### Task 13: Build `cross_instrument_pairs.csv` for Table S2

**Files:**
- Create: `Llunr/n33_expansion/scripts/build_cross_instrument.py`
- Output: `Llunr/n33_expansion/results/cross_instrument_pairs.csv`

[PROCEDURAL]

For NWA 5751 and Talampaya: pair the (pre-existing) Lab Co W-A column length against the (new) I11 W-A column length.

- [ ] **Step 1: Write the script**

File: `Llunr/n33_expansion/scripts/build_cross_instrument.py`

```python
"""Build cross_instrument_pairs.csv: Lab Co vs I11 W-A for NWA 5751 + Talampaya."""
from __future__ import annotations

import csv
import json
from pathlib import Path


N33 = Path("/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/n33_expansion")
PAPER = Path("/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/Paper1_JAC")
SAMPLES = ["NWA 5751", "Talampaya"]


def main() -> None:
    rows = []
    # Lab Co values from the 28-sample CSV
    with open(PAPER / "survey_results_28samples.csv") as f:
        for r in csv.DictReader(f):
            if r["sample_id"] in SAMPLES and r["instrument"] == "Lab Co":
                rows.append({
                    "sample": r["sample_id"],
                    "instrument": "Lab Co",
                    "q_max": r["q_max"],
                    "wa_median_column_length": r["wa_median_column_length"],
                })
    # I11 values from the new per_sample/ JSON
    for s in SAMPLES:
        d = json.load(open(N33 / "results/per_sample" / f"{s.replace(' ', '_')}.json"))
        rows.append({
            "sample": s,
            "instrument": "I11",
            "q_max": "14.6",
            "wa_median_column_length": d["warren_averbach"]["median_column_length"],
        })

    out = N33 / "results/cross_instrument_pairs.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["sample", "instrument", "q_max", "wa_median_column_length"])
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print(f"wrote {out} ({len(rows)} rows)")


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Run**

```bash
"/c/Users/Matthew Izawa/anaconda3/python.exe" n33_expansion/scripts/build_cross_instrument.py
```

Expected: `wrote …/cross_instrument_pairs.csv (4 rows)`.

---

### Task 14: Phase 2 sanity-check gate — AUTHOR REVIEW

**Files:**
- Create: `Llunr/n33_expansion/results/phase2_sanity_check.md`

[AUTHOR TASK]

Inspect the 6 new sample analyses against predictions from spec §6.8.

- [ ] **Step 1: Print the predictions table**

```bash
"/c/Users/Matthew Izawa/anaconda3/python.exe" -c "
import json
from pathlib import Path
p = Path('/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/n33_expansion/results/per_sample')
samples = ['NWA_7465', 'NWA_5478', 'NWA_5751', 'NWA_7831', 'NWA_6693', 'Talampaya']
for s in samples:
    f = p / f'{s}.json'
    if not f.exists():
        print(s, 'MISSING')
        continue
    d = json.load(open(f))
    wa = d.get('warren_averbach', {})
    pdf = d.get('pdf', {})
    print(f'{s:12s}  WA median={wa.get(\"median_column_length\")}  Plag/Pyx={wa.get(\"plag_pyx_ratio\")}  PDF FWHM={pdf.get(\"first_peak_fwhm\")}')
"
```

Expected: 6 lines of analysis output.

- [ ] **Step 2: Author writes the sanity-check note**

File: `Llunr/n33_expansion/results/phase2_sanity_check.md`

```markdown
# Phase 2 sanity check

Date: YYYY-MM-DD
Author: M. R. M. Izawa

## Predictions vs observations

| Sample | Prediction | Observed | Verdict |
|---|---|---|---|
| NWA 5478 (shocked polymict eucrite) | Plag/Pyx < 1 (high-shock side of reversal) | … | … |
| Talampaya (cumulate eucrite S2) | I11 W-A > Lab Co 310 Å; Plag/Pyx > 1 | … | … |
| NWA 7831 (diogenite) | W-A similar to Tatahouine, NWA 6013, NWA 2968 | … | … |
| NWA 6693 (ungrouped, oxidised, low shock) | W-A in upper range; PDF FWHM at low-shock end | … | … |

## Verdict

[GO / PAUSE]. Reasoning: …
```

- [ ] **Step 3: Halt if PAUSE; proceed to Phase 3 if GO**

If any prediction fails, halt and investigate (per spec §11). If all pass, proceed.

---

## Phase 3 — Figure Regeneration

### Task 15: Adapt figure scripts to read the 32-sample CSV

**Files:**
- Inspect (read-only): `Paper1_JAC/paper1_figures.py`, `Paper1_JAC/regen_fig2.py`, `Paper1_JAC/fig7_pdf_examples.py`, `Paper1_JAC/figS1_column_length_schematic.py`
- Modify: each script's CSV path reference

[PROCEDURAL]

- [ ] **Step 1: Identify which scripts reference `survey_results_28samples.csv`**

```bash
cd "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/Paper1_JAC"
grep -l "survey_results_29samples" *.py
```

Expected: a list of figure-generating scripts.

- [ ] **Step 2: For each script, replace the CSV path**

For each script in the grep output, use `Edit` to replace:

```
survey_results_28samples.csv
```

with:

```
survey_results_32samples.csv
```

Also search for hardcoded `29` counts in titles/labels/legends and replace with `33` (and `seven`/`7` with `eight`/`8` if relevant).

- [ ] **Step 3: Verify**

```bash
grep -c "survey_results_29samples" "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/Paper1_JAC/"*.py
```

Expected: 0 occurrences after the change.

---

### Task 16: Regenerate Fig 1

**Files:**
- Run: relevant script(s) producing `Paper1_JAC/paper1_figures/fig1_*.png`
- Manual decision: whether to add an 8th panel for NWA 6693 (if Fig 1 is one-panel-per-group)

[PROCEDURAL]

- [ ] **Step 1: Inspect Fig 1's current panel structure**

```bash
"/c/Users/Matthew Izawa/anaconda3/python.exe" -c "
import matplotlib.image as mpimg
img = mpimg.imread('/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/Paper1_JAC/paper1_figures/fig1_families_vs_column_length.png')
print('shape:', img.shape)
"
```

Open the PNG visually. Note whether panels are per-group or per-shock-stage.

- [ ] **Step 2: Regenerate**

```bash
cd "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/Paper1_JAC"
"/c/Users/Matthew Izawa/anaconda3/python.exe" paper1_figures.py  # or whichever script makes Fig 1
```

- [ ] **Step 3: Inspect output**

Open `paper1_figures/fig1_*.png`. Verify it includes the new samples and the count in any caption text annotation.

- [ ] **Step 4: If 8th panel for NWA 6693 is needed, AUTHOR DECISION**

Note in `Llunr/n33_expansion/results/phase3_figure_decisions.md`.

---

### Task 17: Regenerate Fig 2 (anpig crossover)

**Files:**
- Run: `Paper1_JAC/regen_fig2.py`
- Output: `Paper1_JAC/paper1_figures/fig2_anpig_crossover.png`

[PROCEDURAL]

- [ ] **Step 1: Run the regeneration script**

```bash
cd "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/Paper1_JAC"
"/c/Users/Matthew Izawa/anaconda3/python.exe" regen_fig2.py
```

- [ ] **Step 2: Inspect output**

Open `paper1_figures/fig2_anpig_crossover.png`. Verify:
- 3 new eucrite points are visible: NWA 7465 (mid shock), NWA 5478 (high shock), Talampaya (S2 low shock)
- Plag/Pyx reversal is still visible — and now anchored by 3 high-shock eucrites instead of 2
- Symbol/colour scheme is consistent with the prior version
- Caption text (if rendered) says "33" not "29"

- [ ] **Step 3: Spot-check the new high-shock anchor**

Print the actual Plag/Pyx ratio for NWA 5478:

```bash
"/c/Users/Matthew Izawa/anaconda3/python.exe" -c "
import json
d = json.load(open('/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/n33_expansion/results/per_sample/NWA_5478.json'))
print('NWA 5478 Plag/Pyx:', d['warren_averbach'].get('plag_pyx_ratio'))
"
```

Expected: < 1 (high-shock side of reversal).

---

### Task 18: Regenerate Fig 3 (PDF FWHM vs W-A column length, synchrotron)

**Files:**
- Run: relevant script in `Paper1_JAC/` (likely `paper1_figures.py`)
- Output: `Paper1_JAC/paper1_figures/fig3_pdf_fwhm_vs_column_length.png`

[PROCEDURAL]

- [ ] **Step 1: Run**

```bash
cd "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/Paper1_JAC"
"/c/Users/Matthew Izawa/anaconda3/python.exe" paper1_figures.py
```

- [ ] **Step 2: Inspect**

Open `paper1_figures/fig3_pdf_fwhm_vs_column_length.png`. Verify:
- 24 total synchrotron points (18 existing + 6 new)
- Symbol for NWA 6693 reflects the 8th parent-body group (new legend entry)
- FWHM values overall consistent with the capillary-corrected Phase 1 outcome

---

### Task 19: Regenerate Fig 4 (W-H reliability classification)

**Files:**
- Run: relevant script
- Output: `Paper1_JAC/paper1_figures/fig4_wh_reliability.png`

[PROCEDURAL]

- [ ] **Step 1: Run**

```bash
cd "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/Paper1_JAC"
"/c/Users/Matthew Izawa/anaconda3/python.exe" paper1_figures.py
```

- [ ] **Step 2: Inspect**

Open `paper1_figures/fig4_wh_reliability.png`. Verify the 6 new samples appear and the qualitative reliability classification (likely-OK / likely-unreliable / inconclusive) still tells the same story.

---

### Task 20: Regenerate Fig 5 (families vs Q_max + cross-instrument bonus)

**Files:**
- Modify: relevant script for Fig 5 (cross-instrument bonus addition)
- Output: `Paper1_JAC/paper1_figures/fig5_cross_instrument.png` (or `fig5_families_vs_qmax.png`)

[PROCEDURAL]

Spec-approved bonus: plot NWA 5751 and Talampaya on both panel (b) Lab Co and panel (c) I11, with connecting lines between the same-specimen pair.

- [ ] **Step 1: Identify the script and panel logic**

```bash
grep -n "Fig.*5\|families_vs_qmax\|cross_instrument" "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/Paper1_JAC/"*.py
```

- [ ] **Step 2: Modify the script to add cross-instrument lines**

Open the script. Find the panel-c rendering. Add:

```python
# Cross-instrument pairs for NWA 5751 and Talampaya
CROSS_INSTRUMENT = ["NWA 5751", "Talampaya"]
for sample in CROSS_INSTRUMENT:
    # Look up both Lab Co (panel b) and I11 (panel c) family counts
    lab_co = lookup(sample, instrument="Lab Co")
    i11 = lookup(sample, instrument="I11")
    # Plot both with a connecting line
    ax_b.plot([lab_co["q_max"]], [lab_co["families"]], "o", markerfacecolor="C2", label=sample)
    ax_c.plot([i11["q_max"]], [i11["families"]], "o", markerfacecolor="C2")
    # Connecting line across the panels (requires a parent figure-level annotation)
```

Note: implementing the cross-panel line requires `fig.add_artist(ConnectionPatch(…))`. Alternative: render both points on a SINGLE panel showing families vs Q_max for all samples on one axis, with cross-instrument pairs highlighted. Implementing agent picks the cleaner approach.

- [ ] **Step 3: Run the script**

```bash
cd "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/Paper1_JAC"
"/c/Users/Matthew Izawa/anaconda3/python.exe" paper1_figures.py
```

- [ ] **Step 4: Inspect**

Open `fig5_*.png`. Verify:
- Lab Co panel shows 2 samples (Tatahouine, NWA 6013)
- I11 panel shows 23 samples
- Cross-instrument pairs for NWA 5751 + Talampaya are visible and visually obvious
- Caption updated for new counts

---

### Task 21: Regenerate Fig 6 (W-A column length vs first PDF peak position)

**Files:**
- Run: relevant script
- Output: `Paper1_JAC/paper1_figures/fig6_column_length_vs_pdf_position.png`

[PROCEDURAL]

- [ ] **Step 1: Run**

```bash
cd "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/Paper1_JAC"
"/c/Users/Matthew Izawa/anaconda3/python.exe" paper1_figures.py
```

- [ ] **Step 2: Inspect**

Open `fig6_*.png`. Verify:
- Tagish Lake serpentine outlier still separates cleanly
- NWA 6711 anomalous first-shell distance still separates
- 6 new points visible

---

### Task 22: Regenerate Fig 7 (PDF examples / G(r) curves)

**Files:**
- Modify: `Paper1_JAC/fig7_pdf_examples.py` (point at corrected `.xy` files)
- Output: `Paper1_JAC/paper1_figures/fig7_pdf_examples.png`

[PROCEDURAL]

- [ ] **Step 1: Update the script to read capillary-subtracted `.xy`**

In `fig7_pdf_examples.py`, replace any `.xy` path pointing at raw I11 data with the corresponding `_cap.xy` in `Llunr/n33_expansion/data/i11_corrected/`.

- [ ] **Step 2: Optionally add NWA 5478 and/or NWA 6693 panels**

Author decision. Default: keep the existing 3-4 representative samples. Note any change in `phase3_figure_decisions.md`.

- [ ] **Step 3: Run**

```bash
"/c/Users/Matthew Izawa/anaconda3/python.exe" "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/Paper1_JAC/fig7_pdf_examples.py"
```

- [ ] **Step 4: Inspect**

Open `fig7_pdf_examples.png`. Verify G(r) curves are smoother / less wavy at low r than the v1 version (capillary contribution removed).

---

### Task 23: Generate Fig S2 (capillary before/after demo) — OPTIONAL

**Files:**
- Create: `Paper1_JAC/figS2_capillary_demo.py`
- Output: `Paper1_JAC/paper1_figures/figS2_capillary_demo.png`

[PROCEDURAL]

Per the user's brainstorm decision: "Fig 7 I like but might end up as supplement." So the capillary before/after demonstration goes here as Fig S2.

- [ ] **Step 1: Write the script**

File: `Paper1_JAC/figS2_capillary_demo.py`

```python
"""Fig S2: before/after capillary subtraction demonstration."""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from xrd_profile import XRDProfile


PAPER = Path("/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/Paper1_JAC")
N33 = Path("/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/n33_expansion")
WAVELENGTH = 0.824883

# Pick one representative low-r-sensitive sample
DEMO_SAMPLE = "Tirhert"  # or whatever the author prefers


def main() -> None:
    raw_xy = "[AUTHOR INPUT NEEDED: source path for Tirhert raw I11 .xy]"
    cap_subtracted_xy = N33 / "data/i11_corrected" / "[fill in matching _cap.xy filename]"
    cap_alone_xy = N33 / "data/i11_2026Apr_reduced/1437220-mac-summed.xy"

    prof_raw = XRDProfile.from_file(raw_xy, wavelength=WAVELENGTH)
    prof_corr = XRDProfile.from_file(str(cap_subtracted_xy), wavelength=WAVELENGTH)
    pdf_raw = prof_raw.pdf(q_max=14.6)
    pdf_corr = prof_corr.pdf(q_max=14.6)

    fig, axes = plt.subplots(2, 1, figsize=(8, 8), constrained_layout=True)
    ax = axes[0]
    cap = np.loadtxt(cap_alone_xy)
    ax.plot(cap[:, 0], cap[:, 1], lw=0.6, color="k")
    ax.set_xlabel("2θ (°)")
    ax.set_ylabel("Counts")
    ax.set_title("(a) Blank capillary alone")

    ax = axes[1]
    ax.plot(pdf_raw["r"], pdf_raw["G"], "--", color="0.5", label=f"{DEMO_SAMPLE} raw")
    ax.plot(pdf_corr["r"], pdf_corr["G"], "-", color="k", label=f"{DEMO_SAMPLE} capillary-subtracted")
    ax.set_xlim(0, 10)
    ax.set_xlabel("r (Å)")
    ax.set_ylabel("G(r)")
    ax.set_title(f"(b) PDF before vs after capillary subtraction ({DEMO_SAMPLE})")
    ax.legend(fontsize=9)

    out = PAPER / "paper1_figures/figS2_capillary_demo.png"
    fig.savefig(out, dpi=180)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Author fills in the [AUTHOR INPUT NEEDED] path**

- [ ] **Step 3: Run**

```bash
"/c/Users/Matthew Izawa/anaconda3/python.exe" "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/Paper1_JAC/figS2_capillary_demo.py"
```

- [ ] **Step 4: Inspect**

Open `figS2_capillary_demo.png`. Verify the capillary hump is visible in (a) and the subtraction-induced change in (b) is clearly demonstrable.

---

### Task 24: Phase 3 sanity-check — AUTHOR REVIEW

**Files:**
- Create: `Llunr/n33_expansion/results/phase3_figure_decisions.md`

[AUTHOR TASK]

- [ ] **Step 1: Author reviews all regenerated figures**

Side-by-side compare each new `paper1_figures/figN_*.png` against the corresponding `submission/figures/figN_*.png` (the v1 Inkscape-edited version).

- [ ] **Step 2: Sanity-check the three core claims (spec §7 gate):**

- Fig 2: Plag/Pyx reversal anchored by ≥3 high-shock eucrites — visible?
- Fig 3 / Fig 6: PDF FWHM monotonic-with-shock trend preserved or strengthened?
- Fig 5 cross-instrument bonus: NWA 5751 + Talampaya pairs visually obvious?

- [ ] **Step 3: Write the decision note**

File: `Llunr/n33_expansion/results/phase3_figure_decisions.md`

```markdown
# Phase 3 figure decisions

Date: YYYY-MM-DD

## Per-figure status

| Figure | Approve | Inkscape rework needed | Notes |
|---|---|---|---|
| Fig 1 | … | … | … |
| Fig 2 | … | … | … |
| Fig 3 | … | … | … |
| Fig 4 | … | … | … |
| Fig 5 | … | … | … |
| Fig 6 | … | … | … |
| Fig 7 | … | … | … |
| Fig S2 | … | … | … |

## Sanity-check core claims

- Plag/Pyx reversal anchored by 3 high-shock eucrites: [yes/no]
- PDF FWHM monotonic with shock: [yes/no]
- Cross-instrument bonus visible: [yes/no]

## Punch list for Inkscape rework

…

## Verdict

[GO to Phase 4 / PAUSE]
```

---

## Phase 4 — Text and Metadata Updates

### Task 25: Update Abstract

**Files:**
- Modify: `Paper1_JAC/manuscript_paper1_JAC.md`

[PROCEDURAL]

- [ ] **Step 1: Locate the Abstract section**

```bash
grep -n "^## Abstract\|^# Abstract\|28-meteorite" "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/Paper1_JAC/manuscript_paper1_JAC.md" | head -5
```

- [ ] **Step 2: Apply edits**

Replace:
- `28-meteorite suite spanning seven parent-body groups` → `32-meteorite suite spanning eight parent-body groups`
- Any hardcoded `29` → `33`
- Any `seven` → `eight` where it refers to group count (do not blindly replace other uses)

Do NOT add capillary subtraction language to the abstract (per spec §3.3 brainstorm decision: this is a software paper).

- [ ] **Step 3: Read back the updated Abstract to verify**

```bash
"/c/Users/Matthew Izawa/anaconda3/python.exe" -c "
import re
with open('/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/Paper1_JAC/manuscript_paper1_JAC.md') as f:
    text = f.read()
m = re.search(r'## Abstract\n+(.+?)\n## ', text, re.DOTALL)
print(m.group(1) if m else 'not found')
"
```

Expected: abstract reads "32-meteorite suite spanning eight parent-body groups" with the asymmetric-amorphisation language preserved.

---

### Task 26: Update §2.5 PDF methods (capillary paragraph)

**Files:**
- Modify: `Paper1_JAC/manuscript_paper1_JAC.md`

[PROCEDURAL]

- [ ] **Step 1: Locate §2.5**

```bash
grep -n "^### 2\.5\|^## 2\.5\|pair distribution function" "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/Paper1_JAC/manuscript_paper1_JAC.md" | head -5
```

- [ ] **Step 2: Apply the edit**

Find the existing PDF methods paragraph and immediately after it, insert:

```markdown
For synchrotron measurements, a blank borosilicate capillary scan acquired
under identical beam conditions is subtracted from each sample intensity
prior to PDF computation. A per-sample scale factor is determined by
least-squares matching in the high-Q tail (Q = 12-14 Å⁻¹), a region
dominated by amorphous and incoherent scattering rather than crystalline
Bragg peaks. The corrected intensity is then passed to the iterative
Chebyshev background fit and Lorch-modified Fourier transform described
above.
```

(Edit the surrounding sentence to reference the capillary step before the Chebyshev background, not after, so the order in prose matches the order in code.)

- [ ] **Step 3: Verify**

```bash
grep -A 5 "blank borosilicate capillary" "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/Paper1_JAC/manuscript_paper1_JAC.md"
```

---

### Task 27: Update §2.6 Samples and data collection

**Files:**
- Modify: `Paper1_JAC/manuscript_paper1_JAC.md`

[PROCEDURAL]

- [ ] **Step 1: Locate §2.6**

```bash
grep -n "^### 2\.6\|Samples and data" "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/Paper1_JAC/manuscript_paper1_JAC.md" | head -3
```

- [ ] **Step 2: Apply the edits**

- Update the sample count from 29 → 33
- Add a sentence listing the 4 new samples and 2 upgrades, e.g.:

```markdown
Four additional samples (NWA 7465 — monomict eucrite; NWA 5478 — shocked
polymict eucrite; NWA 7831 — diogenite; NWA 6693 — CC-associated impact
melt, paired with NWA 6704 and treated here as a distinct ungrouped
achondrite parent body) were measured at I11 in April 2026 alongside
synchrotron re-measurements of NWA 5751 and Talampaya (originally
measured at Lab Co). A blank borosilicate capillary background was
acquired under identical conditions during the same beamtime.
```

- [ ] **Step 3: Verify**

```bash
grep -B 1 -A 8 "NWA 7465\|NWA 5478\|NWA 7831\|NWA 6693" "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/Paper1_JAC/manuscript_paper1_JAC.md" | head -30
```

---

### Task 28: Update §3.3 Plag/Pyx reversal narrative

**Files:**
- Modify: `Paper1_JAC/manuscript_paper1_JAC.md`

[PROCEDURAL]

- [ ] **Step 1: Locate §3.3**

```bash
grep -n "^### 3\.3\|Plag.*Pyx\|crossover" "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/Paper1_JAC/manuscript_paper1_JAC.md" | head -5
```

- [ ] **Step 2: Apply edits**

Update the narrative to mention:
- NWA 5478 as a third high-shock anchor alongside Millbillillie (S5) and JaH 626 (S6)
- Talampaya as a low-shock anchor (S2 cumulate eucrite) joining Tirhert (S1) and Bereba (S2-S3)
- NWA 7465 as a mid-ladder eucrite

Reuse the §4.2 asymmetric-amorphisation framing without expanding it.

Example sentence to add:

```markdown
NWA 5478 (shocked polymict eucrite, with shock-melt textures in the
Meteoritical Bulletin entry) provides a third independent high-shock
data point with Plag/Pyx < 1, alongside Millbillillie and JaH 626. The
upgrade of Talampaya to synchrotron resolution places a second
low-shock cumulate eucrite on the high-Plag/Pyx side of the crossover.
```

- [ ] **Step 3: Verify the new Plag/Pyx values are consistent**

```bash
"/c/Users/Matthew Izawa/anaconda3/python.exe" -c "
import json
from pathlib import Path
p = Path('/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/n33_expansion/results/per_sample')
for s in ['NWA_5478', 'NWA_7465', 'Talampaya']:
    d = json.load(open(p / f'{s}.json'))
    print(s, 'Plag/Pyx:', d['warren_averbach'].get('plag_pyx_ratio'))
"
```

Use the actual numbers in the prose if appropriate (table-style or in-text).

---

### Task 29: Update §4.2 asymmetric-amorphisation rescue

**Files:**
- Modify: `Paper1_JAC/manuscript_paper1_JAC.md`

[PROCEDURAL]

- [ ] **Step 1: Locate §4.2**

```bash
grep -n "^### 4\.2\|asymmetric.*amorph\|maskelynite" "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/Paper1_JAC/manuscript_paper1_JAC.md" | head -5
```

- [ ] **Step 2: Apply edits**

Strengthen the §4.2 narrative with NWA 5478 as evidence:

```markdown
The reversal pattern is now sampled by three independent high-shock
eucrites (Millbillillie, JaH 626, NWA 5478), each consistent with the
asymmetric-amorphisation interpretation: plagioclase column lengths
collapse as the structure approaches its ~25-30 GPa amorphisation
threshold (Stöffler 2018) while pyroxene retains crystalline coherence
well below its higher ~50 GPa threshold (Cao 2025).
```

- [ ] **Step 3: Verify**

```bash
grep -A 5 "three independent high-shock" "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/Paper1_JAC/manuscript_paper1_JAC.md"
```

---

### Task 30: Retire §4.3 capillary "would enable" sentence

**Files:**
- Modify: `Paper1_JAC/manuscript_paper1_JAC.md`

[PROCEDURAL]

- [ ] **Step 1: Locate the sentence**

```bash
grep -n "would enable\|blank capillary measurement" "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/Paper1_JAC/manuscript_paper1_JAC.md"
```

Expected: one match at approximately line 163.

- [ ] **Step 2: Replace with a present-tense description**

The old sentence:

```
A blank capillary measurement at the same beam conditions would enable
subtraction of this contribution.
```

is replaced with:

```
A blank borosilicate capillary measurement acquired in the April 2026
beamtime is subtracted from each synchrotron pattern prior to PDF
computation (§2.5), removing the overlapping amorphous capillary
contribution from the first coordination shell and yielding the
capillary-corrected G(r) profiles reported in §3.x.
```

(Adjust §3.x to the actual section number once Phase 4 prose is consolidated.)

- [ ] **Step 3: Verify**

```bash
grep -c "would enable" "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/Paper1_JAC/manuscript_paper1_JAC.md"
```

Expected: 0.

---

### Task 31: Update §3 Results (other sections) with new sample counts

**Files:**
- Modify: `Paper1_JAC/manuscript_paper1_JAC.md`

[PROCEDURAL]

- [ ] **Step 1: Search for any hardcoded "29" or "seven groups"**

```bash
grep -n "\b29\b\|seven parent-body\|seven groups" "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/Paper1_JAC/manuscript_paper1_JAC.md"
```

- [ ] **Step 2: Replace each occurrence with 33 / eight as appropriate**

Use `Edit` per occurrence; do not blindly replace because "29" might appear in unrelated contexts (e.g., year 2029 — unlikely but possible).

- [ ] **Step 3: Verify zero hits**

```bash
grep -c "28-meteorite\|29 meteorites\|survey of 29\|seven parent-body" "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/Paper1_JAC/manuscript_paper1_JAC.md"
```

Expected: 0.

---

### Task 32: Update §5 Conclusions

**Files:**
- Modify: `Paper1_JAC/manuscript_paper1_JAC.md`

[PROCEDURAL]

- [ ] **Step 1: Locate §5**

```bash
grep -n "^## 5\.\|^## Conclusion" "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/Paper1_JAC/manuscript_paper1_JAC.md"
```

- [ ] **Step 2: Apply edits**

- Update any sample counts (29 → 33, seven → eight)
- Strengthen the asymmetric-amorphisation conclusion (three high-shock anchors now)
- Add a conclusion or modify an existing one to mention capillary subtraction's role, e.g.:

```markdown
8. A blank borosilicate capillary measurement, subtracted with a high-Q-tail
   scaling, removes the overlapping amorphous contribution from synchrotron
   PDFs and improves the resolution of the first coordination shell, enabling
   the PDF FWHM trend with shock stage to be quantified across the full
   synchrotron suite.
```

- [ ] **Step 3: Verify**

Inspect §5 in the file.

---

### Task 33: Update Table 1 inventory

**Files:**
- Modify: `Paper1_JAC/table1_sample_inventory.md`

[PROCEDURAL]

- [ ] **Step 1: Read current Table 1**

```bash
cat "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/Paper1_JAC/table1_sample_inventory.md"
```

- [ ] **Step 2: Apply edits**

- Caption: "29 meteorites" → "33 meteorites"
- Remove the existing Lab Co rows for NWA 5751 and Talampaya
- Add new I11 rows for NWA 5751 and Talampaya with their updated W-A median column length, Q_max = 14.6, and PDF FWHM
- Add 4 new I11 rows: NWA 7465, NWA 5478, NWA 7831, NWA 6693 with metadata from Task 11 (use `[AUTHOR INPUT NEEDED]` where lit data is pending)
- Update PDF FWHM column for the 17 existing I11 rows using values from `n33_expansion/results/per_sample_pdf_refits/*.json`

- [ ] **Step 3: Verify row count**

```bash
grep -c "^\|" "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/Paper1_JAC/table1_sample_inventory.md"
```

Expected: 34 (33 data rows + 1 header row; or different counting depending on Markdown table style).

---

### Task 34: Update Table 2 (Plag/Pyx ratios for eucrites)

**Files:**
- Locate and modify the Table 2 source file in `Paper1_JAC/`

[PROCEDURAL]

Table 2 was inserted inline in the manuscript at the §3.3 rewrite during the McCausland response work. Find it:

- [ ] **Step 1: Locate Table 2 markdown**

```bash
grep -n "Table 2\|table2" "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/Paper1_JAC/manuscript_paper1_JAC.md"
```

- [ ] **Step 2: Apply edits**

Add 3 new rows: NWA 7465, NWA 5478, Talampaya, with their Plag/Pyx ratios from `n33_expansion/results/per_sample/`. Order rows by shock stage. NWA 5478 enters as the third high-shock eucrite.

- [ ] **Step 3: Verify**

Inspect the updated Table 2 and confirm 3 new rows present.

---

### Task 35: Create Table S2 (cross-instrument comparison)

**Files:**
- Create: `Paper1_JAC/Table_S2_cross_instrument.md`

[PROCEDURAL]

- [ ] **Step 1: Write the table**

File: `Paper1_JAC/Table_S2_cross_instrument.md`

```markdown
# Table S2. Cross-instrument comparison for NWA 5751 and Talampaya

Two cumulate eucrites originally measured at Lab Co K-α (Q_max = 5.4 Å⁻¹)
were re-measured at I11 (Q_max = 14.6 Å⁻¹) during the April 2026 beamtime.
The Lab Co rows are reported here to demonstrate the apparent-microstructure
effect of restricted Q-range on the same physical specimen (cf. §4.3).

| Sample | Instrument | Q_max (Å⁻¹) | W-A median column length (Å) | PDF first-peak FWHM (Å) |
|---|---|---|---|---|
| NWA 5751 | Lab Co | 5.4 | … | — |
| NWA 5751 | I11    | 14.6 | … | … |
| Talampaya | Lab Co | 5.4 | 310 | — |
| Talampaya | I11    | 14.6 | … | … |

Values pulled from `Llunr/n33_expansion/results/cross_instrument_pairs.csv`.
```

Use actual values from `cross_instrument_pairs.csv` and `per_sample/*.json`.

- [ ] **Step 2: Reference Table S2 in §4.3**

In §4.3 of `manuscript_paper1_JAC.md`, add a sentence:

```markdown
The lab-Q_max effect is illustrated quantitatively in Table S2 for NWA 5751
and Talampaya, both measured at Lab Co and I11.
```

---

### Task 36: Update cover letter

**Files:**
- Modify: `Paper1_JAC/submission/cover_letter.md`

[PROCEDURAL]

- [ ] **Step 1: Apply edits**

Replace in `cover_letter.md`:
- `28-meteorite suite` → `32-meteorite suite`
- `seven parent-body groups` → `eight parent-body groups`

Do NOT mention capillary subtraction (consistent with abstract decision).

- [ ] **Step 2: Verify**

```bash
grep "32-meteorite\|eight parent-body" "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/Paper1_JAC/submission/cover_letter.md"
```

---

### Task 37: Update suggested reviewers (light)

**Files:**
- Modify: `Paper1_JAC/submission/suggested_reviewers.md`

[PROCEDURAL]

The suggested reviewers list is unchanged from v1. Only update if any introductory text references the survey size.

- [ ] **Step 1: Check for sample-count references**

```bash
grep "29\|seven" "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/Paper1_JAC/submission/suggested_reviewers.md"
```

- [ ] **Step 2: If hits, update to 33/eight; if not, no action**

---

### Task 38: Update acknowledgments and data availability statement

**Files:**
- Modify: `Paper1_JAC/manuscript_paper1_JAC.md` (acknowledgments section)
- Modify: `Paper1_JAC/submission/build_auxiliary_docx.py` (data availability builder)

[PROCEDURAL]

- [ ] **Step 1: Locate acknowledgments**

```bash
grep -n "^## Acknowledg\|## Data availability\|Diamond Light Source\|DLS proposal" "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/Paper1_JAC/manuscript_paper1_JAC.md"
```

- [ ] **Step 2: Add Apr 2026 Diamond proposal number**

Add a sentence in acknowledgments:

```markdown
Additional I11 beamtime in April 2026 (proposal [AUTHOR INPUT NEEDED:
Apr 2026 proposal number]) is gratefully acknowledged.
```

- [ ] **Step 3: Update data availability**

Edit `build_auxiliary_docx.py:build_data_availability()` (or the inline content) to add the new raw `.dat` files (1437211–1437220 from 24 Apr 2026) and the capillary scan to the archived dataset description.

- [ ] **Step 4: Verify**

Spot-check the relevant sections.

---

### Task 39: Add references for new samples + capillary citation (if any)

**Files:**
- Modify: `Paper1_JAC/manuscript_paper1_JAC.md` (references section)

[PROCEDURAL]

- [ ] **Step 1: Locate references section**

```bash
grep -n "^## References\|^# References" "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/Paper1_JAC/manuscript_paper1_JAC.md"
```

- [ ] **Step 2: Add Meteoritical Bulletin / classification refs for the 4 new samples**

Insert refs from Task 11 metadata table. Use `[AUTHOR INPUT NEEDED]` placeholders where verifications are still pending.

- [ ] **Step 3: Add a citation for capillary subtraction methodology (optional)**

Candidates: Billinge group / Wright et al. methodology papers. Author decision on whether this is warranted. If yes, add to references and cite in §2.5.

---

### Task 40: Update `build_auxiliary_docx.py` inline content

**Files:**
- Modify: `Paper1_JAC/submission/build_auxiliary_docx.py`

[PROCEDURAL]

The previous cover-letter rewrite (Greek α, no corporate-speak, no date) was not yet propagated into the inline Python content. Re-sync now.

- [ ] **Step 1: Read the current cover-letter Python content**

```bash
grep -n "def build_cover_letter\|28-meteorite\|seven parent-body" "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/Paper1_JAC/submission/build_auxiliary_docx.py"
```

- [ ] **Step 2: Update the inline strings**

- `28-meteorite suite` → `32-meteorite suite`
- `seven parent-body groups` → `eight parent-body groups`
- Verify Greek α is preserved (Cu Kα, Co Kα).

- [ ] **Step 3: Verify**

```bash
grep "32-meteorite\|eight parent-body" "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/Paper1_JAC/submission/build_auxiliary_docx.py"
```

---

## Phase 5 — Co-author Re-circulation

### Task 41: Build v2 clean `.docx`

**Files:**
- Output: `Paper1_JAC/manuscript_paper1_JAC_v2_clean.docx`

[PROCEDURAL]

- [ ] **Step 1: Run the existing docx builder**

```bash
cd "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/Paper1_JAC"
"/c/Users/Matthew Izawa/anaconda3/python.exe" build_docx.py
```

(Or whichever script is the canonical builder. Check filenames first.)

- [ ] **Step 2: Verify output**

```bash
ls -la "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/Paper1_JAC/manuscript_paper1_JAC.docx"
```

- [ ] **Step 3: Rename to v2 clean**

```bash
mv "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/Paper1_JAC/manuscript_paper1_JAC.docx" \
   "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/Paper1_JAC/manuscript_paper1_JAC_v2_clean.docx"
```

- [ ] **Step 4: Open in Word/LibreOffice and spot-check**

Verify:
- Headings are black (not blue)
- No em-dashes in prose
- Captions/footnotes at 10pt
- Greek α renders correctly
- Sample count says 33 throughout

---

### Task 42: Generate v2 track-changes `.docx` vs PJAM-MI v1

**Files:**
- Output: `Paper1_JAC/manuscript_paper1_JAC_v2_trackchanges.docx`

[PROCEDURAL]

- [ ] **Step 1: Confirm the v1 baseline file exists**

```bash
ls "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/Paper1_JAC/manuscript_paper1_JAC_submission_EAC may 9 2026-PJAM-MI.docx"
```

- [ ] **Step 2: Generate the diff using Word's Compare feature OR pandoc**

**Option A (Word):** Open Word → Review → Compare → select v1 baseline and v2 clean as Original/Revised → save as `manuscript_paper1_JAC_v2_trackchanges.docx`.

**Option B (pandoc):** Run:

```bash
pandoc \
  "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/Paper1_JAC/manuscript_paper1_JAC_submission_EAC may 9 2026-PJAM-MI.docx" \
  -o /tmp/v1.md
pandoc \
  "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/Paper1_JAC/manuscript_paper1_JAC_v2_clean.docx" \
  -o /tmp/v2.md
# Manual review: diff /tmp/v1.md /tmp/v2.md
```

[AUTHOR TASK] — pandoc-based track changes are not natively round-trip-perfect. The Word Compare flow is recommended.

- [ ] **Step 3: Spot-check that the most consequential edits show as tracked**

Open the result in Word. Verify:
- §4.3 capillary "would enable" sentence shows as deleted, new sentence as inserted
- Abstract `28-meteorite suite` → `32-meteorite suite` shown as edit
- New rows in Table 1 / Table 2 shown as inserts

---

### Task 43: Build `cover_letter_v2.docx`

**Files:**
- Output: `Paper1_JAC/submission/cover_letter_v2.docx`

[PROCEDURAL]

- [ ] **Step 1: Run the auxiliary builder**

```bash
cd "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/Paper1_JAC/submission"
"/c/Users/Matthew Izawa/anaconda3/python.exe" build_auxiliary_docx.py
```

- [ ] **Step 2: Rename and verify**

```bash
mv "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/Paper1_JAC/submission/cover_letter.docx" \
   "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/Paper1_JAC/submission/cover_letter_v2.docx"
ls -la "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/Paper1_JAC/submission/cover_letter_v2.docx"
```

Open in Word/LibreOffice. Verify:
- Greek α renders (Cu Kα, Co Kα)
- No em-dashes; no coloured headings
- Sample count says 33; group count says 8
- No date stamp
- 5 named reviewers (Ungár, Billinge, Langenhorst, Tsuchiyama, Toby)

---

### Task 44: Generate `figure_thumbnails.pdf` for the bundle

**Files:**
- Create: `Llunr/n33_expansion/scripts/build_figure_thumbnails.py`
- Output: `Paper1_JAC/submission/figure_thumbnails.pdf`

[PROCEDURAL]

- [ ] **Step 1: Write the script**

File: `Llunr/n33_expansion/scripts/build_figure_thumbnails.py`

```python
"""Assemble all v2 figures into a single PDF for skim review by co-authors."""
from __future__ import annotations

from pathlib import Path

import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


PAPER = Path("/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/Paper1_JAC")
OUT = PAPER / "submission/figure_thumbnails.pdf"


def main() -> None:
    figs = sorted((PAPER / "paper1_figures").glob("fig*.png"))
    with PdfPages(OUT) as pdf:
        # Pack 4 figures per page
        for i in range(0, len(figs), 4):
            fig, axes = plt.subplots(2, 2, figsize=(11, 8.5), constrained_layout=True)
            for ax, fpath in zip(axes.flat, figs[i:i + 4]):
                ax.imshow(mpimg.imread(fpath))
                ax.set_title(fpath.stem, fontsize=8)
                ax.axis("off")
            for ax in axes.flat[len(figs[i:i + 4]):]:
                ax.axis("off")
            pdf.savefig(fig)
            plt.close(fig)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Run**

```bash
"/c/Users/Matthew Izawa/anaconda3/python.exe" "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/n33_expansion/scripts/build_figure_thumbnails.py"
```

Expected: `wrote …/figure_thumbnails.pdf`.

---

### Task 45: Write `v2_change_summary.md` cover note

**Files:**
- Create: `Paper1_JAC/submission/v2_change_summary.md`

[PROCEDURAL]

- [ ] **Step 1: Write the cover note**

File: `Paper1_JAC/submission/v2_change_summary.md`

```markdown
# JAC paper v2 — summary of changes since the v1 (PJAM-MI) round

## What's new in v2

- **4 new I11 samples** added from the April 2026 Diamond Light Source
  beamtime: NWA 7465 (monomict eucrite), NWA 5478 (shocked polymict
  eucrite), NWA 7831 (diogenite), NWA 6693 (CC-associated impact melt
  paired with NWA 6704, treated as a distinct ungrouped achondrite
  parent body — bringing the survey total to 8 groups).
- **2 samples upgraded** from Lab Co K-α to I11 synchrotron: NWA 5751
  (howardite) and Talampaya (cumulate eucrite). Their Lab Co
  measurements are retained in a new supplementary Table S2 as a
  cross-instrument comparison.
- **Blank capillary subtraction** is now integrated into the
  synchrotron PDF pipeline (§2.5). The validation behind this change
  is summarised in `Llunr/n33_expansion/phase1_validation/decision_note.md`.
- **All data-bearing figures regenerated** (Figs 1–7). One new
  supplementary figure: Fig S2 shows the capillary before/after
  demonstration. Inkscape rework punch list is in
  `Llunr/n33_expansion/results/phase3_figure_decisions.md`.
- The §4.3 "would enable a blank capillary measurement" sentence
  is replaced by a present-tense description of the subtraction
  in place.

## What strengthened (not changed)

- The Plag/Pyx reversal claim is now anchored by **three independent
  high-shock eucrites** (Millbillillie, JaH 626, NWA 5478) instead
  of two, and a second low-shock cumulate eucrite (Talampaya) on
  the high-Plag/Pyx side.
- The PDF FWHM monotonic-with-shock trend is now capillary-corrected
  across the entire synchrotron suite.
- The cross-instrument apparent-microstructure caveat (§4.3) is now
  backed by same-specimen measurements (Table S2).

## What's unchanged

- Authorship and order (11 authors, as in PJAM-MI).
- Suggested reviewers (Ungár, Billinge, Langenhorst, Tsuchiyama, Toby).
- The asymmetric-amorphisation interpretation in §4.2.
- The `xrd_profile` methodology and package version.

## Outstanding `[AUTHOR INPUT NEEDED]` items

- Apr 2026 DLS proposal number (acknowledgments + data availability)
- Lit shock stage + primary references for the 4 new samples
- Reviewer email / affiliation confirmations
- Final decision on capillary-subtraction citation (Billinge group?
  Wright et al.?)

## Deadline for co-author comments

[author to fill in]
```

- [ ] **Step 2: Verify**

Inspect the cover note. Add or remove sections as needed for the specific co-author audience.

---

### Task 46: Final bundle assembly and email — AUTHOR TASK

**Files:**
- Bundle (assemble by hand): `manuscript_paper1_JAC_v2_clean.docx`, `..._v2_trackchanges.docx`, `cover_letter_v2.docx`, `suggested_reviewers.docx`, `figure_thumbnails.pdf`, `v2_change_summary.md`

[AUTHOR TASK]

- [ ] **Step 1: Confirm all bundle files exist**

```bash
ls -la "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/Paper1_JAC/manuscript_paper1_JAC_v2_clean.docx"
ls -la "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/Paper1_JAC/manuscript_paper1_JAC_v2_trackchanges.docx"
ls -la "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/Paper1_JAC/submission/cover_letter_v2.docx"
ls -la "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/Paper1_JAC/submission/suggested_reviewers.docx"
ls -la "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/Paper1_JAC/submission/figure_thumbnails.pdf"
ls -la "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/Paper1_JAC/submission/v2_change_summary.md"
```

- [ ] **Step 2: Pin v1 as the approved snapshot**

```bash
cp "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/Paper1_JAC/manuscript_paper1_JAC_submission_EAC may 9 2026-PJAM-MI.docx" \
   "/c/Users/Matthew Izawa/Documents/Dan Applin/Llunr/Paper1_JAC/manuscript_paper1_JAC_v1_approved.docx"
```

- [ ] **Step 3: Email co-authors** (author task — no script)

Use the cover-note text as the email body. Attach the 4 .docx files + 1 PDF + 1 .md.

Set a deadline for comments (typically 1-2 weeks).

- [ ] **Step 4: Track co-author feedback**

Use a single document (`Llunr/n33_expansion/results/coauthor_feedback_v2.md`) to consolidate comments as they arrive. Address each in a final pass before JAC submission.

---

## Plan complete

Tasks 1–46 cover the full v1 → v2 manuscript revision. Phase 1 and Phase 2 each have explicit author-review gates that can abort the plan (Tasks 6 and 14). Phase 3 has a per-figure decision step (Task 24). Phase 4 is mechanical prose updates. Phase 5 produces the co-author bundle and ends in an author task (email).

**Post-plan workflow:** After v2 co-author feedback (1-2 weeks) and integration, invoke `superpowers:finishing-a-development-branch` to complete the work and submit to JAC.
