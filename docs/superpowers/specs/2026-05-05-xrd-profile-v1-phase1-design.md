---
title: xrd_profile v1.0 — Phase 1 design
date: 2026-05-05
status: approved
target_version: 0.3.0 (Phase 1), 1.0.0 (Phase 2)
author: Matthew R. M. Izawa
---

# xrd_profile v1.0 — Phase 1 design

## 1. Background and motivation

`xrd_profile` is a Python toolkit for quantitative analysis of powder X-ray
diffraction peak profiles, currently at v0.2.0. The package is the
methodological subject of a Journal of Applied Crystallography (JAC)
manuscript (draft complete, with co-authors at the time of writing) and the
analysis backbone for a follow-on Meteoritics & Planetary Science (MAPS)
study of HED shock systematics.

In its current form the package handles plagioclase and pyroxene multi-phase
analyses through ad-hoc, per-script reference-peak builders: each user
script constructs a `pymatgen.Structure` (sometimes from inline lattice
parameters and atomic coordinates, sometimes from a CIF file path) and
duplicates a `build_ref` helper to convert pymatgen output into the dict
format the analysis functions expect. There is no library-owned `Phase`
abstraction; the duplication is a maintenance burden, and adding a new
mineral (olivine, hematite, quartz) is a copy-paste-and-edit operation
rather than a single-line call.

The MAPS follow-on will need olivine support (olivine-bearing diogenites,
mesosiderites). More broadly, a v1.0 push positions the package as a
general-purpose tool that any group with powder XRD data can pick up,
load arbitrary phases from CIF files, and run a configurable bundle of
analyses (Williamson-Hall, Warren-Averbach, pair distribution function,
Scherrer) with sensible toggling.

## 2. Goals

- Make `xrd_profile` capable of representing **any crystalline phase** via
  a uniform `Phase` abstraction loaded from CIF or built inline from
  lattice parameters and atomic coordinates.
- Provide a **single-call convenience helper** (`run_all`) for running
  configurable subsets of methods (W-H, W-A, PDF, Scherrer) on
  configurable lists of phases, returning a structured nested-dict result.
- Expose **Scherrer constant K and shape factor** as first-class kwargs,
  with a small lookup table for common crystallite shapes.
- Lay the **foundations for v1.0** so that Phase 2 (instrumental
  deconvolution, size distributions, anisotropic reporting, modified
  Williamson-Hall plug-in interface, YAML/CLI frontend, Sphinx
  documentation, PyPI release) lands as additive diffs without retrofitting
  Phase 1 architecture.

## 3. Non-goals

- Phase 1 does not implement instrumental broadening deconvolution. The
  `instrumental=` kwarg is reserved on guided-method signatures and
  raises `NotImplementedError` if a non-None value is passed; Phase 2
  wires it up.
- Phase 1 does not implement crystallite size distributions
  (lognormal/normal fits to W-A column-length distributions). Deferred
  to Phase 2.
- Phase 1 does not implement anisotropic / direction-resolved size
  reporting. Deferred to Phase 2.
- Phase 1 does not implement the modified Williamson-Hall plug-in
  interface (Ungár-Borbély dislocation contrast factors). Deferred to
  Phase 2 as droppable scope (silicate contrast factors are a research
  problem in their own right; if Phase 2 scope tightens, this drops).
- Phase 1 does not introduce a YAML/TOML config or CLI. The Python API
  is the only interface in v0.3.0.
- Phase 1 does not introduce Sphinx/mkdocs documentation. Inline
  docstrings are sufficient until v1.0.
- Phase 1 does not change any existing public API signatures, default
  values, or numerical results. All additions are strictly additive
  optional kwargs.
- Phase 1 does not touch the JAC manuscript draft, the figures in
  `Paper1_JAC/figures/`, or the analysis pipeline scripts in `Llunr/`
  (e.g., `paper1_figures.py`, `run_guided_wh.py`, `compile_survey.py`).

## 4. Constraints

- **Strictly additive policy.** No public API removed or renamed. New
  functionality lands as new optional kwargs. Default behavior identical
  to v0.2.0 byte-for-byte on the array-based API. A numerical regression
  test (Section 8.4) enforces this.
- **`pymatgen` stays optional.** Required only for `Phase.from_cif()`
  and `Phase.from_lattice_params()`. Array-based functions retain no
  pymatgen dependency. Installed via `pip install xrd_profile[cif]`.
- **Phase 2 must land without Phase 1 rewrites.** The kwarg hook pattern
  adopted in Phase 1 (`phase=`, plus a placeholder-but-rejecting
  `instrumental=` accepted on guided methods) is the same pattern Phase 2
  fills in. No retrofitting.
- **JAC paper reproducibility preserved.** A `git tag v0.2.0` is created
  on the pre-Phase-1 HEAD before any Phase 1 work begins. The paper's
  data-availability link can pin to `v0.2.0` (frozen) or `v0.3.0` /
  `main` (evolving) at the author's discretion at submission time.

## 5. Sequencing: Phase 1 vs Phase 2

### Phase 1 (this design, v0.3.0)

1. `Phase` class with `from_cif()` and `from_lattice_params()` constructors.
2. `phase=` kwarg added to `XRDProfile.guided_williamson_hall()` and
   `guided_warren_averbach()`, alongside existing `ref_d=` / `ref_peaks=`.
3. `instrumental=` kwarg reserved on the same methods (no-op,
   `NotImplementedError` if used).
4. `XRDProfile.run_all(methods=[...], phases=[...], wh={...}, wa={...},
   pdf={...}, scherrer={...})` convenience helper.
5. Scherrer K and shape factor exposed via `K` and `shape` kwargs on
   `scherrer()`, `modified_scherrer()`, `XRDProfile.scherrer()`,
   `XRDProfile.modified_scherrer()`. Sentinel default `K=None` preserves
   v0.2.0 behavior.
6. Five example CIFs in `examples/cifs/` (forsterite, anorthite,
   pigeonite, quartz, hematite).
7. One small bundled fixture pattern (`examples/data/tirhert_subset.xy`,
   ~30-50 KB).
8. New example script `examples/multi_phase_olivine.py`.
9. Existing example scripts (`lab_lunar_meteorite.py`,
   `synchrotron_low_shock.py`, `synchrotron_high_shock.py`) updated to
   use `Phase`. Originals preserved verbatim in `examples/legacy/`.
10. Test additions: `test_phases.py`, `test_run_all.py`, `test_scherrer.py`,
    `test_backward_compat.py`. Fixtures in `tests/fixtures/`.
11. README updated with new quickstart pointing to Phase API.

### Phase 2 (deferred to a separate design, v1.0.0)

1. Instrumental broadening deconvolution (LaB6/Si standard input,
   wired into the reserved `instrumental=` kwarg).
2. Crystallite size distributions (lognormal/normal fits to W-A
   column-length distributions, exposed as new keys in W-A result dict).
3. Anisotropic / direction-resolved size reporting (per-axis sizes from
   W-A families, exposed as new keys).
4. Modified Williamson-Hall plug-in interface (user-supplied contrast
   factor function as `contrast_factor=` kwarg). Droppable.
5. YAML/TOML config + CLI frontend (`xrd_profile run experiment.yaml`).
6. Sphinx or mkdocs API reference, tutorial notebooks for representative
   mineral systems, contribution guide.
7. PyPI release of v1.0.0.

## 6. Architecture

### 6.1 New file

`xrd_profile/phases.py` — contains the `Phase` class, the standalone
`build_reference_peaks(structure, wavelength, ...)` function for users
who already have a `pymatgen.Structure`, and a lazy-import helper that
gives a clear `ImportError` if pymatgen is not installed.

### 6.2 Modified files

- `xrd_profile/profile.py` — adds `XRDProfile.run_all()` method;
  adds optional `phase=` kwarg to `guided_williamson_hall()` and
  `guided_warren_averbach()`; reserves `instrumental=` kwarg with
  `NotImplementedError` guard.
- `xrd_profile/scherrer.py` — adds `shape` kwarg with
  `SCHERRER_K_FOR_SHAPE` lookup table; changes `K` default to a
  sentinel (`None`) that resolves to `0.9` when both `K` and `shape`
  are `None` (preserves v0.2.0 behavior). Same treatment for
  `modified_scherrer()`.
- `xrd_profile/__init__.py` — exports `Phase`, `build_reference_peaks`,
  `SCHERRER_K_FOR_SHAPE`.

### 6.3 Unchanged files

`conversions.py`, `noise.py`, `peak_detection.py`, `pdf.py`,
`warren_averbach.py`, `williamson_hall.py` — no internal changes.
Phase routing happens at the `XRDProfile`-method level; the analysis
modules continue to receive `ref_d` / `ref_peaks` arrays exactly as today.

### 6.4 Dependencies

`pyproject.toml` adds an optional extra:

```toml
[project.optional-dependencies]
cif = ["pymatgen>=2023.0"]
```

- `pip install xrd_profile` → no pymatgen, array-based API only.
- `pip install xrd_profile[cif]` → pymatgen installed, `Phase` works.

The README quickstart shows the `[cif]` install in the recommended
path so new users do not hit the import error.

### 6.5 Directory layout after Phase 1

```
xrd_profile/                                 # repo root
├── CITATION.cff
├── LICENSE
├── README.md
├── pyproject.toml
├── setup.py
├── docs/
│   └── superpowers/
│       └── specs/
│           └── 2026-05-05-xrd-profile-v1-phase1-design.md
├── examples/
│   ├── cifs/
│   │   ├── Forsterite.cif
│   │   ├── Anorthite.cif
│   │   ├── Pigeonite.cif
│   │   ├── Quartz.cif
│   │   ├── Hematite.cif
│   │   └── SOURCES.md
│   ├── data/
│   │   ├── tirhert_subset.xy
│   │   └── README.md
│   ├── legacy/
│   │   ├── lab_lunar_meteorite.py
│   │   ├── synchrotron_low_shock.py
│   │   └── synchrotron_high_shock.py
│   ├── lab_lunar_meteorite.py              # updated to use Phase
│   ├── synchrotron_low_shock.py            # updated to use Phase
│   ├── synchrotron_high_shock.py           # updated to use Phase
│   ├── multi_phase_olivine.py              # NEW
│   └── *.png                               # existing figures unchanged
├── scripts/
│   └── regenerate_goldens.py               # NEW
├── tests/
│   ├── fixtures/
│   │   ├── tirhert_subset.xy
│   │   ├── forsterite_for_test.cif
│   │   └── golden_v0.2.0_results.json
│   ├── test_conversions.py                 # existing
│   ├── test_pdf.py                         # existing
│   ├── test_phases.py                      # NEW
│   ├── test_run_all.py                     # NEW
│   ├── test_scherrer.py                    # NEW
│   └── test_backward_compat.py             # NEW
└── xrd_profile/
    ├── __init__.py                         # updated exports
    ├── conversions.py                      # unchanged
    ├── noise.py                            # unchanged
    ├── pdf.py                              # unchanged
    ├── peak_detection.py                   # unchanged
    ├── phases.py                           # NEW
    ├── profile.py                          # updated
    ├── scherrer.py                         # updated
    ├── warren_averbach.py                  # unchanged
    └── williamson_hall.py                  # unchanged
```

### 6.6 Versioning

`__version__ = '0.3.0'` at end of Phase 1 (semver MINOR bump =
additive new feature, no breaks). `v1.0.0` at end of Phase 2.

## 7. API design

### 7.1 `Phase` class

```python
# xrd_profile/phases.py

class Phase:
    """A crystalline phase with a known structure, used to generate
    reference peak positions for guided peak detection."""

    def __init__(self, structure, name: str = ''):
        self.structure = structure  # pymatgen Structure
        self.name = name

    @classmethod
    def from_cif(cls, path, name=None) -> 'Phase':
        """Load from a CIF file. Default name = file stem."""

    @classmethod
    def from_lattice_params(cls, a, b, c, alpha, beta, gamma,
                            species, coords, name='') -> 'Phase':
        """Build inline from lattice parameters + atomic positions.
        Preserves the existing inline-anorthite pattern in user scripts."""

    def get_ref_peaks(self, wavelength,
                      two_theta_range=(5, 90),
                      min_intensity=3.0) -> list[dict]:
        """List of dicts with keys d, two_theta, intensity, h, k, l.
        Suitable for ref_peaks= argument to guided_warren_averbach()."""

    def get_ref_d(self, wavelength,
                  two_theta_range=(5, 90),
                  min_intensity=3.0,
                  sorted_by_intensity=True) -> np.ndarray:
        """d-spacing array sorted by intensity. Suitable for ref_d=."""


def build_reference_peaks(structure, wavelength,
                          two_theta_range=(5, 90),
                          min_intensity=3.0) -> list[dict]:
    """Standalone function for users who already have a pymatgen
    Structure in hand. Equivalent to Phase(structure).get_ref_peaks(...)."""
```

A lazy-import helper in `phases.py` raises a clear `ImportError` if
pymatgen is missing:

> `pymatgen is required for Phase.from_cif and from_lattice_params.
> Install with: pip install xrd_profile[cif]`

### 7.2 Integration with `XRDProfile`

```python
# In xrd_profile/profile.py:

def guided_williamson_hall(self, ref_d=None, *, phase=None,
                           instrumental=None, **kwargs):
    if phase is not None and ref_d is not None:
        raise ValueError("pass either ref_d or phase, not both")
    if phase is not None:
        tt_range = (float(self.two_theta.min()),
                    float(self.two_theta.max()))
        ref_d = phase.get_ref_d(self.wavelength, two_theta_range=tt_range)
    if instrumental is not None:
        raise NotImplementedError(
            "instrumental= is reserved for Phase 2 / v1.0; "
            "see xrd_profile roadmap")
    # ... existing v0.2.0 logic with ref_d
```

Same pattern for `guided_warren_averbach(ref_peaks=None, *, phase=None,
instrumental=None, ...)`.

API contracts:

- `phase=` and `ref_d=` (or `ref_peaks=`) are mutually exclusive; raises
  `ValueError` if both provided.
- When `phase=` is used, `two_theta_range` defaults to the profile's
  actual data range (e.g., `(10, 148)` for I11 synchrotron, `(5, 80)`
  for Co Kα Bruker). Users wanting a different range compute `ref_d`
  manually with `phase.get_ref_d(wavelength, two_theta_range=(...))`
  and pass that.
- `min_intensity=3.0` default in `get_ref_d` / `get_ref_peaks` matches
  existing user practice.

### 7.3 `XRDProfile.run_all()`

```python
def run_all(self,
            methods: list[str] | None = None,
            phases: list[Phase] | Phase | None = None,
            wh: dict | None = None,
            wa: dict | None = None,
            pdf: dict | None = None,
            scherrer: dict | None = None,
            instrumental=None) -> dict:
    """
    Run a configurable bundle of analyses.

    Parameters
    ----------
    methods : list of {'wh', 'wa', 'pdf', 'scherrer'} or None
        Which analyses to run. None = all four.
    phases : list of Phase, single Phase, or None
        For guided W-H and W-A. None = unguided forms run.
        A single Phase is wrapped in a list internally.
    wh, wa, pdf, scherrer : dict or None
        Per-method kwargs. e.g. wh={'n_sigma': 3.0, 'tolerance_d': 0.02}.
    instrumental : reserved for Phase 2; raises NotImplementedError.

    Returns
    -------
    dict. Examples:
      {'wh': {phase.name: result_dict, ...},   # if phases provided
       'wa': {phase.name: result_dict, ...},
       'pdf': pdf_result_dict,
       'scherrer': scherrer_result_dict}

      {'wh': result_dict,                       # if phases=None
       'wa': result_dict,
       ...}
    """
```

Per-method results retain their existing v0.2.0 dict schemas — `run_all`
is purely a dispatcher, not a result-shape redesign. The unified
results-object question (typed dataclass vs nested dict) is deferred to
Phase 2.

### 7.4 Scherrer K and shape factor

```python
# xrd_profile/scherrer.py

SCHERRER_K_FOR_SHAPE = {
    'spherical':   0.94,   # Langford & Wilson (1978)
    'cubic':       0.83,
    'cylindrical': 1.84,
    'platey':      1.0,
}

def scherrer(fwhm_deg, two_theta_deg, wavelength,
             K: float | None = None,
             shape: str | None = None):
    """
    Scherrer crystallite size.

    K, shape resolution:
      both None  → K = 0.9 (v0.2.0 default, identical behavior)
      shape only → K = SCHERRER_K_FOR_SHAPE[shape]
      K only     → K = K (shape ignored)
      both       → K wins, shape silently ignored

    The two-None case preserves v0.2.0 numerical results exactly.
    """
```

Same change to `modified_scherrer()`, `XRDProfile.scherrer()`,
`XRDProfile.modified_scherrer()`. The `SCHERRER_K_FOR_SHAPE` table is
exported from `__init__.py` so users can introspect what values get
applied.

## 8. Examples and tests

### 8.1 Example CIFs

`examples/cifs/` ships five CIFs:

| File | Mineral | Source | Rationale |
|---|---|---|---|
| `Forsterite.cif` | Mg₂SiO₄ olivine | Crystallography Open Database | Phase 2 driver (HED follow-on) |
| `Anorthite.cif` | CaAl₂Si₂O₈ plagioclase | Crystallography Open Database | Matches existing Paper 1 analyses |
| `Pigeonite.cif` | Mg-Fe-Ca pyroxene (Morimoto) | HOSERLab CIF library | Matches existing Paper 1 analyses |
| `Quartz.cif` | SiO₂ | Crystallography Open Database | Universal calibration phase |
| `Hematite.cif` | Fe₂O₃ | Crystallography Open Database | Mars / OSIRIS-REx relevance |

`examples/cifs/SOURCES.md` lists the COD ID and citation for each
external CIF, plus license attribution where applicable.

### 8.2 Bundled data fixture

`examples/data/tirhert_subset.xy` is a ~30-50 KB downsampled real XRD
pattern (Tirhert subset trimmed to ~1000 points across 10-80°), with
`examples/data/README.md` pointing to the full data location for users
who want the ungrouped pattern. This bundled file lets the example
scripts and tests run without external dependencies.

### 8.3 Example scripts

- **`examples/multi_phase_olivine.py` (new):** demonstrates the Phase
  API on a simulated olivine-bearing diogenite scenario. Loads the
  bundled fixture pattern, builds Forsterite + Pigeonite via
  `Phase.from_cif()`, runs
  `profile.run_all(methods=['wh', 'wa', 'pdf'], phases=[...], ...)`,
  prints results, produces a 4-panel figure (pattern, W-H comparison,
  W-A anisotropic sizes per phase, PDF). Canonical "what does the new
  API look like in practice" reference.
- **Updated `lab_lunar_meteorite.py`, `synchrotron_low_shock.py`,
  `synchrotron_high_shock.py`:** rewritten to use
  `Phase.from_lattice_params()` (for anorthite, where refined lattice
  parameters are inline) and `Phase.from_cif()` (for pigeonite, where
  it loads from a CIF). The duplicated `build_ref` helper is removed.
  Numerical results identical to v0.2.0 (regression-tested). Hardcoded
  data paths parameterized as a top-of-file `DATA_PATH = ...` for
  clarity but stay local to the author's machine.
- **`examples/legacy/lab_lunar_meteorite.py`,
  `synchrotron_low_shock.py`, `synchrotron_high_shock.py`:** v0.2.0
  originals copied verbatim. Frozen as historical artifacts.

### 8.4 Tests

```
tests/
├── fixtures/
│   ├── tirhert_subset.xy
│   ├── forsterite_for_test.cif
│   └── golden_v0.2.0_results.json
├── test_conversions.py             # existing, unchanged
├── test_pdf.py                     # existing, unchanged
├── test_phases.py                  # NEW (~6 tests)
├── test_run_all.py                 # NEW (~5 tests)
├── test_scherrer.py                # NEW (~4 tests)
└── test_backward_compat.py         # NEW (~5 regression tests)
```

**`test_phases.py`:** from_cif loads expected formula; from_lattice_params
produces equivalent ref_d to from_cif on the same structure;
get_ref_d returns sorted by intensity when sorted_by_intensity=True;
min_intensity filtering works as documented; required keys
(d, two_theta, intensity, h, k, l) are present in get_ref_peaks output;
missing-pymatgen monkeypatch produces a clear ImportError with install
instructions in the message.

**`test_run_all.py`:** methods=None runs all four methods; methods subset
runs only the subset; per-phase result keys appear when phases are
given; unguided result shape (single dict, not nested) when phases=None;
instrumental=anything raises NotImplementedError with the documented
message text.

**`test_scherrer.py`:** default `scherrer(fwhm, tt, wl)` produces
v0.2.0-identical output (compares to a frozen array); shape='spherical'
uses K=0.94; explicit K wins over shape; SCHERRER_K_FOR_SHAPE has the
documented values.

**`test_backward_compat.py` (most important):** loads
`tirhert_subset.xy`, runs the v0.2.0 array-based API
(`guided_williamson_hall(ref_d, ...)`,
`guided_warren_averbach(ref_peaks, ...)`,
`compute_pdf_sine(...)`, `scherrer(fwhm, tt, wl)`), asserts each result
matches the frozen `golden_v0.2.0_results.json` within tight numerical
tolerances (1e-10 for direct outputs, 1e-6 for derived quantities).
A `scripts/regenerate_goldens.py` is committed for the rare case of
intentional numerical updates (e.g., bug fixes), with explicit
documentation that regeneration requires explicit reasoning in the
commit message.

### 8.5 README updates

- Quickstart shows `pip install xrd_profile[cif]` and the Phase-based
  pattern as primary, with a "lower-level array-based API" subsection
  for users with d-spacings from Rietveld.
- New "Phases" section documents `Phase.from_cif`, `from_lattice_params`,
  and the example CIFs.
- New "run_all" section shows the multi-phase, multi-method bundled call.
- Link to `examples/multi_phase_olivine.py` as the canonical demo.

## 9. Phase 2 forward-compatibility

What Phase 1's design buys Phase 2:

- **Instrumental deconvolution (item c):** `instrumental=` kwarg already
  in `guided_williamson_hall` and `guided_warren_averbach` signatures.
  Phase 2 replaces the `NotImplementedError` with the actual deconvolution
  logic. No signature change. `Phase` (or a similar `InstrumentalStandard`)
  represents LaB6 / Si standards via the same `from_cif()` pathway.
- **Size distributions (b):** lognormal/normal fits to W-A column-length
  distributions added as new keys in the existing W-A result dict
  (e.g., `'column_length_distribution': {'mean': ..., 'sigma': ...,
  'shape': 'lognormal'}`). Existing keys remain untouched. Strictly
  additive.
- **Anisotropic reporting (d):** new keys in W-A result for per-axis
  sizes (e.g., `'sizes_by_axis': {'a': ..., 'b': ..., 'c': ...}`). Same
  pattern.
- **Modified W-H plug-in (e):** new optional kwarg on
  `guided_williamson_hall(contrast_factor=callable_or_none, ...)`. Same
  kwarg pattern as `phase=` and `instrumental=`. Droppable without
  breaking anything.
- **YAML/CLI:** the CLI is a thin frontend that builds `Phase` objects,
  constructs `XRDProfile`, and calls `run_all`. No new internal API
  needed — Phase 1's `run_all` is already the right ingestion point.
- **Sphinx/mkdocs:** every public symbol added in Phase 1 has a
  docstring suitable for autodoc. No retrofit required.

The strict-additive policy plus the kwarg-reservation pattern means
Phase 2 lands as a series of "add code paths" diffs rather than
"modify code paths" diffs. Smaller review surface, less regression risk.

## 10. Phase 1 acceptance criteria

Phase 1 is done when:

1. `Phase.from_cif()` and `Phase.from_lattice_params()` both work
   end-to-end on the five bundled CIFs and on inline parameters.
2. `XRDProfile.guided_williamson_hall(phase=phase, ...)` produces
   identical numerical results to
   `guided_williamson_hall(ref_d=phase.get_ref_d(wavelength), ...)`.
3. `XRDProfile.run_all(...)` returns the documented nested-dict
   structure for all method/phase combinations.
4. `scherrer(fwhm, tt, wl)` with no kwargs produces v0.2.0-identical
   output (regression test passes).
5. All five example CIFs are present, licensed, and parse via
   `Phase.from_cif()`.
6. `multi_phase_olivine.py` runs end-to-end on the bundled fixture.
7. The three updated existing examples run end-to-end with identical
   numbers vs. the legacy versions.
8. New tests pass; existing tests still pass; `test_backward_compat.py`
   confirms no numerical drift.
9. README quickstart shows the Phase pattern; legacy patterns
   documented separately.
10. `CHANGELOG.md` documents the v0.2.0 → v0.3.0 changes.
11. `git tag v0.2.0` exists pointing at the pre-Phase-1 HEAD.
12. `git tag v0.3.0` created on the Phase-1-complete commit.

## 11. Decisions deferred to Phase 2 design

- **Unified results object (typed dataclass).** Phase 1 keeps per-method
  dicts. Phase 2 may introduce a hierarchical dataclass
  (`results.williamson_hall.crystallite_size`) once the field set has
  stabilized through real Phase 2 use.
- **`pymatgen` as hard vs. optional dependency.** Phase 1 keeps it
  optional. If real-world Phase 2 usage shows that the overwhelming
  majority of installs use `[cif]`, the v1.0 PyPI release may promote
  pymatgen to a hard dependency at that point.
- **`InstrumentalStandard` as a separate class vs. reuse of `Phase`.**
  Phase 2 design will determine whether the LaB6 / Si standard gets
  its own class or rides on `Phase`.
- **Guided Scherrer (per-phase peak sets).** Phase 1 leaves Scherrer
  unguided. Phase 2 may add a `phase=` kwarg to `XRDProfile.scherrer()`.

## 12. Open questions

None at this time. All clarifications from the 2026-05-05 brainstorming
session have been incorporated.
