---
title: xrd_profile v1.0 — Phase 2 design
date: 2026-05-05
status: approved
target_versions: 0.4.0, 0.5.0, 1.0.0
predecessor: 2026-05-05-xrd-profile-v1-phase1-design.md
author: Matthew R. M. Izawa
---

# xrd_profile v1.0 — Phase 2 design

## 1. Background and motivation

`xrd_profile` is a Python toolkit for quantitative analysis of powder
X-ray diffraction peak profiles. Phase 1 (v0.3.0, 2026-05-05) shipped a
`Phase` abstraction with `from_cif()` / `from_lattice_params()`
constructors, an `XRDProfile.run_all()` convenience helper, exposed
Scherrer K and shape factor as first-class kwargs, and reserved (but
did not implement) the `instrumental=` kwarg on guided W-H, guided W-A,
and `run_all`. The Phase 1 spec deliberately laid hooks so Phase 2 work
lands as additive diffs rather than retrofits.

Phase 2 fills those hooks in and brings the package to v1.0:
instrumental broadening deconvolution, crystallite size distributions,
direction-resolved size reporting, a modified-Williamson-Hall plug-in
interface, a YAML config + CLI frontend, an mkdocs documentation site,
and the infrastructure for a PyPI release. The PyPI public upload
itself is sequenced behind the JAC manuscript decision; the v1.0.0 tag
ships the build infrastructure and a TestPyPI dry-run, with the
canonical repository remaining private until the broader project is
ready for public consumption.

The work is structured as three tagged releases — v0.4.0 (instrumental
deconvolution + size distributions + guided Scherrer), v0.5.0
(anisotropic reporting + modified-W-H plug-in), and v1.0.0 (CLI + docs
+ release infrastructure) — rather than as a single monolithic v1.0
push. The stepwise sequencing matches the precedent set in Phase 1
(`git tag v0.2.0` frozen for JAC reproducibility before any Phase 1
work began): each Phase 2 tag is a citation-stable artifact so a
reviewer asking "which version had instrumental correction" can be
answered with a single tag.

## 2. Goals

- Wire up the **instrumental broadening deconvolution** Phase 1
  reserved, using a Caglioti polynomial for W-H/Scherrer (which consume
  scalar FWHM) and Stokes Fourier deconvolution for W-A (which already
  lives in Fourier space). Each method gets the deconvolution that
  matches its mathematical basis.
- Add **crystallite size distributions** (lognormal and normal fits to
  the W-A column-length distribution) as additive keys on each W-A
  family entry.
- Add **direction-resolved size reporting** as new top-level keys
  (`sizes_by_hkl`, `size_distributions_by_hkl`) on the W-A result, with
  optional crystal-axis projection (`sizes_by_axis`) when the lattice
  symmetry is at least orthorhombic.
- Add the **modified-Williamson-Hall plug-in interface** (a
  `contrast_factor=` callable kwarg) so users supplying their own
  contrast-factor function can run Ungár-Borbély-style analyses without
  the package committing to a contrast-factor library for
  poorly-characterised systems (silicates).
- Provide a **YAML/CLI frontend** (`xrd_profile run experiment.yaml`)
  as a thin layer over `XRDProfile.run_all()`, gated behind a new
  `[cli]` optional extra so the core library install stays minimal.
- Ship an **mkdocs-material documentation site** (theory pages,
  per-mineral-system tutorial notebooks, auto-generated API reference,
  contribution guide) under the package's existing private repo, with
  a manually-triggered GitHub Actions workflow ready to flip to
  automatic deployment when the repo goes public.
- Land the **PyPI release infrastructure** (build configuration, twine
  wiring, complete project metadata, TestPyPI dry-run) so that the
  public upload to PyPI is a single command when the JAC manuscript and
  repo-visibility decisions are made.

## 3. Non-goals

- Phase 2 does not introduce a typed-dataclass results object. Phase 1
  Section 11 deferred this to "Phase 2 once the field set has
  stabilized through real Phase 2 use." The field set stays mutating
  through v0.4 and v0.5; locking a dataclass schema mid-Phase-2 is
  premature. Deferred to v2.x.
- Phase 2 does not promote `pymatgen` from optional to hard dependency.
  Stays under the `[cif]` extra. v2.x decision once PyPI usage data
  exists.
- Phase 2 does not ship a contrast-factor library. The plug-in
  interface accepts a user-supplied callable; one worked example
  (`examples/modified_wh_cubic.py` with the Ungár-Tichy cubic form) is
  provided as a template. Materials-science contrast factors for
  silicates remain a research problem outside the package's scope.
- Phase 2 does not ship a built-in instrumental-profile registry with
  pre-fit values for the Misasa Rigaku, Winnipeg Bruker, or Diamond
  I11. The infrastructure (`InstrumentalProfile.from_registry()`, JSON
  format, lookup path) lands in v0.4.0; the registry itself ships
  empty and is populated as a separate calibration deliverable.
- Phase 2 does not add new analysis methods (PDF variants beyond what
  ships in v0.3.0; Rietveld integration; pole-figure inversion; etc.).
- Phase 2 does not publish the package on PyPI. The public upload
  waits for the JAC manuscript decision and the
  repo-visibility decision (Section 10). v1.0.0 ships the
  infrastructure and the git tag; the repo remains private.
- Phase 2 does not deploy a public documentation site. The
  GitHub Actions workflow is scaffolded with `on: workflow_dispatch`
  only (manual trigger). When the repo goes public, the trigger flips
  to `on: push` and Pages activates.
- Phase 2 does not touch the JAC manuscript draft, the figures in
  `Paper1_JAC/figures/`, or the analysis pipeline scripts in `Llunr/`
  (`paper1_figures.py`, `run_guided_wh.py`, `compile_survey.py`).
- Phase 2 does not remove the `Hematite.cif` example shipped in v0.3.0.
  Removing a public artifact would break the strict-additivity
  contract. The CIF stays; the v1.0 documentation positions the
  package's validated scope as silicates only (Section 9.4).

## 4. Constraints

- **Tiered strict-additivity policy.** Each Phase 2 tag (v0.4.0,
  v0.5.0, v1.0.0) ships a new golden-results fixture. All earlier-tier
  goldens (v0.2.0, v0.3.0, ...) continue to pass at every later
  release under key-subset value-equality semantics: every key in the
  earlier-tier golden retains its numerical value to within
  1e-10 / 1e-6 tolerances, while the live result may carry additional
  keys introduced at later tiers (Section 8.2.5). The numerical
  regression test (`tests/test_backward_compat.py`) runs all tiers
  each release; the newest tier locks the API surface added at that
  tag.
- **Defaults preserve numerical results.** `instrumental=None`,
  `axis_projection=False`, `contrast_factor=None`, and `phase=None`
  (on Scherrer) all reduce to v0.3.0 numerics on existing keys. A user
  upgrading from v0.3.0 → v0.5.0 without invoking any v0.4+ feature
  sees the same values for every result-dict key that existed in
  v0.3.0.
- **Optional extras kept compositional.** `[cif]` (v0.3.0, pymatgen),
  `[cli]` (v1.0.0, pyyaml + pydantic). Pure-numpy install stays
  minimal: `pip install xrd_profile` works without either extra and
  retains the array-based public API.
- **No JAC pipeline / manuscript edits.** `Llunr/paper1_figures.py`,
  `Llunr/run_guided_wh.py`, `Llunr/compile_survey.py`, and
  `Paper1_JAC/` are untouched. Phase 2 is package-internal.
- **No real-instrument standard data committed.** The synthetic LaB6
  fixture (`tests/fixtures/synthetic_lab6.xy`) is generated by a
  committed script with documented Caglioti coefficients. Real LaB6 /
  Si patterns are user-supplied at runtime.
- **Phase 2 must not retrofit Phase 1 architecture.** All v0.4.0
  changes either (a) wire up a kwarg already reserved in v0.3.0, or
  (b) add new keys to existing result dicts, or (c) add new public
  classes/modules. No v0.3.0 public signature is modified beyond
  removing the `NotImplementedError` guard on `instrumental=`.

## 5. Sequencing: three tagged releases

### 5.1 v0.4.0 — instrumental deconvolution + size distributions + guided Scherrer

1. New module `xrd_profile/instrumental.py` containing
   `InstrumentalStandard` and `InstrumentalProfile` classes (Section 7.1).
2. Caglioti polynomial fitting / evaluation (`fit_caglioti`,
   `caglioti_fwhm_at`) implemented inside `instrumental.py` for use by
   W-H and Scherrer.
3. Stokes Fourier deconvolution (`stokes_deconvolve`) implemented
   inside `instrumental.py` for use by W-A. Operates on per-peak
   Fourier coefficients.
4. `XRDProfile.guided_williamson_hall(..., instrumental=)` removes the
   `NotImplementedError` guard, accepts `InstrumentalStandard` or
   `InstrumentalProfile`, and applies the Caglioti FWHM correction
   (Section 7.2.1).
5. `XRDProfile.guided_warren_averbach(..., instrumental=)` accepts
   `InstrumentalStandard` only and applies Stokes Fourier deconvolution
   to each peak family before harmonic decomposition. Passing an
   `InstrumentalProfile` raises a clear `ValueError` (Section 7.2.2).
6. `XRDProfile.scherrer(..., phase=, instrumental=)` adds both kwargs.
   `phase=` filters detected peaks to those matching the phase's
   reference d-spacings within tolerance (Section 7.2.3).
   `instrumental=` accepts the same types as W-H.
7. W-A result dict gains a per-family `'size_distribution'` key
   (Section 7.3) with lognormal and normal fits to the column-length
   distribution. Always computed when `n_valid_L >= 4`; `None`
   otherwise.
8. `xrd_profile/__init__.py` exports `InstrumentalStandard`,
   `InstrumentalProfile`.
9. Synthetic LaB6 fixture at `tests/fixtures/synthetic_lab6.xy`,
   produced by `scripts/build_synthetic_standard.py`.
10. Tests: `test_instrumental.py`, `test_size_distributions.py`,
    extension of `test_backward_compat.py` to add a `golden_v0.3.0`
    tier.
11. README updated with a Phase-2 quickstart showing instrumental
    correction. CHANGELOG entry.
12. `git tag v0.4.0` on the v0.4.0-complete commit.

### 5.2 v0.5.0 — anisotropic reporting + modified-W-H plug-in

1. W-A result dict gains top-level keys `'sizes_by_hkl'` and
   `'size_distributions_by_hkl'` aggregating the per-family data
   (Section 7.4). Always computed when families exist.
2. `XRDProfile.guided_warren_averbach(..., axis_projection=False)`
   adds the `axis_projection` kwarg. Default `False` → no
   `'sizes_by_axis'` key emitted (the always-on `'sizes_by_hkl'` and
   `'size_distributions_by_hkl'` from item 1 are still added; v0.4.0
   keys retain their values). `True` → adds `'sizes_by_axis'` key when
   the phase lattice has orthorhombic-or-higher symmetry; raises a
   clear `ValueError` for monoclinic / triclinic with the documented
   reason.
3. `XRDProfile.guided_williamson_hall(..., contrast_factor=None)` adds
   the `contrast_factor` kwarg. Default `None` → W-H result
   value-identical to v0.4.0 on every v0.4.0 key. Callable →
   modified-W-H ordinate `(sinθ/λ)·√C̄_hkl` (Section 7.5).
4. New example `examples/modified_wh_cubic.py` with the Ungár-Tichy
   cubic form and a user-set `q` parameter, demonstrating the callable
   contract.
5. Tests: `test_anisotropic.py`, `test_modified_wh.py`, extension of
   `test_backward_compat.py` to add a `golden_v0.4.0` tier.
6. README and CHANGELOG updates.
7. `git tag v0.5.0` on the v0.5.0-complete commit.

### 5.3 v1.0.0 — YAML CLI + docs site + PyPI infrastructure

1. New module `xrd_profile/cli.py` (entry point `main()`).
2. New module `xrd_profile/config_schema.py` (Pydantic v2 models for
   the YAML config — Section 7.6).
3. `pyproject.toml`: add `[cli]` extra (`pyyaml>=6.0`,
   `pydantic>=2.0`); add `[project.scripts]` entry for the
   `xrd_profile` console script; complete `[project.urls]` metadata
   (Section 10).
4. mkdocs-material site under `docs/`: index, quickstart, four theory
   concept pages, three tutorial notebooks (silicate-only:
   olivine-diogenite, pyroxene-eucrite, quartz-calibration), API
   reference auto-generated via mkdocstrings, scope page, migration
   guide, contribution guide.
5. `mkdocs.yml` at repo root; `docs/` source tree; `site/`
   gitignored.
6. `.github/workflows/docs.yml` with `on: workflow_dispatch` only
   (manual trigger). Flipping to `on: push` is a one-line edit when
   the repo goes public.
7. `CONTRIBUTING.md` and standard `.github/` issue / PR templates.
8. TestPyPI dry-run via `scripts/release.sh` (or `Makefile` target):
   builds the sdist + wheel, twine-checks them, uploads to TestPyPI,
   verifies install. Public-PyPI upload is a separate command, not
   part of the v1.0.0 tag.
9. Tests: extension of `test_backward_compat.py` to add a
   `golden_v0.5.0` tier; `test_cli.py` covering the YAML
   load + dispatch path.
10. `git tag v1.0.0` on the v1.0.0-complete commit.

### 5.4 Sequencing rationale

- **v0.4.0 first.** Instrumental correction is the largest scientific
  payload; it benefits most from a tag of its own (cite-able as "the
  version with instrumental correction"). Size distributions ride
  along because they require zero new module surface (W-A already
  produces A_size; size-distribution fits are pure additive keys) and
  because the Caglioti work pulls us into `scherrer.py` anyway, where
  guided Scherrer is a five-line addition.
- **v0.5.0 second.** Anisotropic / per-(hkl) reporting is small but
  builds on the v0.4.0 size-distribution keys (each (hkl) family now
  carries a distribution; the v0.5.0 aggregation lifts that into
  top-level dicts). Modified-W-H plug-in is independent of
  v0.4.0 work and could land in either tier; it lands in v0.5.0
  because the v0.4.0 surface is already substantial.
- **v1.0.0 last.** CLI is purposefully a thin layer over `run_all`,
  so it has to land *after* the analysis API stabilises. Docs work is
  large but mostly orthogonal to code changes — could be done in
  parallel with v0.4.0/v0.5.0 implementation, but the tutorial
  notebooks need v0.4.0 + v0.5.0 features to demonstrate, so the
  natural order is code-first, docs-last.

## 6. Architecture

### 6.1 New files (per release)

**v0.4.0:**
- `xrd_profile/instrumental.py` — `InstrumentalStandard`,
  `InstrumentalProfile`, internal helpers `_caglioti_fit`,
  `_caglioti_fwhm_at`, `_stokes_deconvolve`,
  `_lognormal_fit_to_a_size`, `_normal_fit_to_a_size`. Targeting
  ~280 LOC.
- `tests/test_instrumental.py` — Caglioti fit recovers synthesis
  parameters; Stokes inverse-FT is stable; both classes round-trip
  correctly.
- `tests/test_size_distributions.py` — lognormal / normal recovery on
  synthetic A_size data; threshold for `None`-result; coverage of the
  edge cases.
- `tests/fixtures/synthetic_lab6.xy` — generated, ~50 KB.
- `tests/fixtures/golden_v0.3.0_results.json` — frozen v0.3.0 outputs
  on the same `tirhert_subset.xy` fixture.
- `scripts/build_synthetic_standard.py` — reproducibly generates the
  synthetic LaB6 fixture from documented Caglioti coefficients.

**v0.5.0:**
- `examples/modified_wh_cubic.py` — Ungár-Tichy cubic form, q
  parameter, ~80 LOC.
- `tests/test_anisotropic.py` — `sizes_by_hkl` aggregation;
  `axis_projection` with cubic, orthorhombic, hexagonal, and
  monoclinic phases; correct `ValueError` on low-symmetry.
- `tests/test_modified_wh.py` — identity callable reproduces
  `contrast_factor=None` results; non-trivial callable shifts
  abscissa as documented.
- `tests/fixtures/golden_v0.4.0_results.json`.

**v1.0.0:**
- `xrd_profile/cli.py` — Click-or-argparse command dispatcher
  (argparse, no new dependency), ~120 LOC.
- `xrd_profile/config_schema.py` — Pydantic v2 models, ~200 LOC.
- `xrd_profile/registry/` — package-data directory for
  `InstrumentalProfile.from_registry()` lookups; ships empty (a
  `README.md` documenting the JSON file format and the `XDG_DATA_HOME`
  fallback path for user-supplied profiles).
- `docs/` tree:
    - `index.md`, `quickstart.md`, `migration.md`, `scope.md`,
      `changelog.md` (include of CHANGELOG.md).
    - `concepts/williamson-hall.md`, `concepts/warren-averbach.md`,
      `concepts/instrumental-correction.md`, `concepts/pdf.md`.
    - `tutorials/olivine-diogenite.ipynb`,
      `tutorials/pyroxene-eucrite.ipynb`,
      `tutorials/quartz-calibration.ipynb`.
    - `api/` (auto-generated by mkdocstrings; populated at build time,
      no committed files beyond a `nav` stub in `mkdocs.yml`).
- `mkdocs.yml` at repo root.
- `.github/workflows/docs.yml` (workflow_dispatch only).
- `.github/ISSUE_TEMPLATE/`, `.github/PULL_REQUEST_TEMPLATE.md`.
- `CONTRIBUTING.md`.
- `tests/test_cli.py`.
- `tests/fixtures/golden_v0.5.0_results.json`.
- `scripts/release.sh` (or `Makefile` target).

### 6.2 Modified files

**v0.4.0:**
- `xrd_profile/profile.py` — three methods updated:
    - `guided_williamson_hall`: replace `NotImplementedError` guard
      with dispatch to `_apply_caglioti_correction`.
    - `guided_warren_averbach`: replace `NotImplementedError` guard
      with dispatch to `_apply_stokes_deconvolution` (only for
      `InstrumentalStandard`; `InstrumentalProfile` raises clear
      `ValueError`).
    - `scherrer` and `modified_scherrer`: add `phase=` and
      `instrumental=` kwargs (Section 7.2.3).
    - `run_all`: replace `NotImplementedError` guard with passthrough
      to per-method dispatch.
- `xrd_profile/williamson_hall.py` — internal `_apply_caglioti` helper
  invoked when an `InstrumentalProfile` is in scope; modifies the
  fitted FWHM array via `β_corr² = β_obs² − β_inst²` (Gaussian-quadrature
  combination, since W-H assumes Gaussian-broadening dominance).
- `xrd_profile/warren_averbach.py` — Stokes deconvolution invoked
  inside the per-family loop, between Fourier-coefficient extraction
  and harmonic separation. Per-family `size_distribution` key added
  to the result entries.
- `xrd_profile/scherrer.py` — accept `phase=` and `instrumental=`,
  filter peaks by phase, apply Caglioti subtraction.
- `xrd_profile/__init__.py` — add `InstrumentalStandard`,
  `InstrumentalProfile` to imports and `__all__`.

**v0.5.0:**
- `xrd_profile/warren_averbach.py` — top-level aggregation of
  `sizes_by_hkl`, `size_distributions_by_hkl` after the per-family
  loop; conditional `sizes_by_axis` when `axis_projection=True`.
- `xrd_profile/williamson_hall.py` — `contrast_factor` kwarg threaded
  through to the abscissa computation.

**v1.0.0:**
- `pyproject.toml` — `[cli]` extra; `[project.scripts]` entry;
  `[project.urls]` complete; classifiers.
- `README.md` — link to the docs site (when public); CLI quickstart;
  full extras documentation.
- `xrd_profile/__init__.py` — bump `__version__ = '1.0.0'`.

### 6.3 Unchanged files (across all of Phase 2)

- `xrd_profile/conversions.py`, `xrd_profile/noise.py`,
  `xrd_profile/peak_detection.py`, `xrd_profile/pdf.py`,
  `xrd_profile/phases.py`. Phase 2 does not modify any v0.3.0
  internal logic; all changes are at the `XRDProfile`-method
  boundary, in new modules, or as additive result-dict keys.
- `Llunr/paper1_figures.py`, `Llunr/run_guided_wh.py`,
  `Llunr/compile_survey.py`, `Paper1_JAC/`. JAC pipeline untouched.

### 6.4 Dependencies

`pyproject.toml` after v1.0.0:

```toml
[project]
dependencies = [
    "numpy>=1.24",
    "scipy>=1.11",
    "matplotlib>=3.7",
]

[project.optional-dependencies]
cif = ["pymatgen>=2023.0"]
cli = ["pyyaml>=6.0", "pydantic>=2.0"]
docs = ["mkdocs-material>=9.5", "mkdocstrings[python]>=0.24",
        "mkdocs-jupyter>=0.24"]
dev = ["pytest>=7.0", "pytest-cov>=4.0", "twine>=5.0", "build>=1.0"]

[project.scripts]
xrd_profile = "xrd_profile.cli:main"
```

`pip install xrd_profile` → core (numpy, scipy, matplotlib) only.
`pip install xrd_profile[cif]` → adds pymatgen for `Phase`.
`pip install xrd_profile[cli]` → adds the CLI.
`pip install xrd_profile[cli,cif]` → both.
`pip install xrd_profile[docs]` → for documentation builders.
`pip install -e .[dev,cif,cli,docs]` → contributor install.

### 6.5 Versioning

`__version__` progresses `'0.3.0'` → `'0.4.0'` → `'0.5.0'` → `'1.0.0'`.
Each version bump lands in the same commit as the corresponding
`CHANGELOG.md` entry and the corresponding `git tag`.

## 7. API design

### 7.1 Instrumental classes

```python
# xrd_profile/instrumental.py

class InstrumentalStandard:
    """Full instrumental characterisation: a Phase (LaB6 / Si CIF)
    plus a measured pattern of that standard at the same instrument
    settings as the sample. Supports both Caglioti FWHM correction
    (for W-H, Scherrer) and Stokes Fourier deconvolution (for W-A)."""

    def __init__(self, phase: 'Phase',
                 two_theta: np.ndarray, intensity: np.ndarray,
                 wavelength: float, name: str = ''):
        self.phase = phase
        self.two_theta = np.asarray(two_theta, dtype=float)
        self.intensity = np.asarray(intensity, dtype=float)
        self.wavelength = wavelength
        self.name = name
        self._caglioti_cache = None
        self._fourier_cache = None

    @classmethod
    def from_cif_and_pattern(cls, cif: str,
                             two_theta, intensity,
                             wavelength: float,
                             name: str = '') -> 'InstrumentalStandard':
        """Convenience constructor: load Phase from CIF, attach the
        measured pattern."""

    def caglioti_fit(self) -> 'InstrumentalProfile':
        """Fit U·tan²θ + V·tanθ + W to the standard's measured FWHMs
        at each reference peak. Cached. Returns an InstrumentalProfile
        wrapping the fit (for downstream W-H / Scherrer use)."""

    def fourier_coefficients(self, peak_d: float, n_coeffs: int = 20):
        """Return (L, A_inst_L) for the standard's peak nearest peak_d.
        Cached per-peak. Used by Stokes deconvolution in W-A."""

    def to_json(self, path) -> None: ...
    @classmethod
    def from_json(cls, path) -> 'InstrumentalStandard': ...


class InstrumentalProfile:
    """Lightweight: Caglioti coefficients (U, V, W) plus wavelength and
    optional name. Sufficient for W-H and Scherrer correction. W-A
    raises a clear ValueError if given an InstrumentalProfile, since
    Stokes Fourier deconvolution requires the measured pattern."""

    def __init__(self, U: float, V: float, W: float,
                 wavelength: float, name: str = ''):
        self.U = U
        self.V = V
        self.W = W
        self.wavelength = wavelength
        self.name = name

    @classmethod
    def from_standard(cls, std: InstrumentalStandard) -> 'InstrumentalProfile':
        """Convenience: fit Caglioti to the standard, return the profile."""
        return std.caglioti_fit()

    @classmethod
    def from_registry(cls, name: str) -> 'InstrumentalProfile':
        """Look up a pre-fit profile by name. v0.4.0 ships an empty
        registry; raises KeyError until users populate it via
        `to_json` to the registry path."""

    def fwhm_at(self, two_theta_deg: float) -> float:
        """Caglioti FWHM at a given 2θ. Used by W-H / Scherrer
        broadening correction."""
        theta = np.deg2rad(two_theta_deg / 2.0)
        return np.sqrt(self.U * np.tan(theta)**2
                       + self.V * np.tan(theta)
                       + self.W)

    def to_json(self, path) -> None: ...
    @classmethod
    def from_json(cls, path) -> 'InstrumentalProfile': ...
```

### 7.2 `XRDProfile` integration

#### 7.2.1 W-H

```python
def guided_williamson_hall(self, ref_d=None, ..., *,
                            phase=None, instrumental=None,
                            contrast_factor=None,    # NEW v0.5.0
                            **kwargs):
    if phase is not None and ref_d is not None:
        raise ValueError("pass either ref_d or phase, not both")
    if phase is not None:
        ref_d = phase.get_ref_d(self.wavelength,
                                 two_theta_range=(self.two_theta.min(),
                                                  self.two_theta.max()))
    # NEW v0.4.0: instrumental correction
    if instrumental is None:
        inst_profile = None
    elif isinstance(instrumental, InstrumentalStandard):
        inst_profile = instrumental.caglioti_fit()
    elif isinstance(instrumental, InstrumentalProfile):
        inst_profile = instrumental
    else:
        raise TypeError(
            f"instrumental= must be InstrumentalStandard, "
            f"InstrumentalProfile, or None; got {type(instrumental)}")
    # ... existing logic, with FWHM correction applied if inst_profile
    # is not None, and abscissa modified if contrast_factor is not None.
```

The Caglioti correction is applied as
β²_corr = max(β²_obs − β²_inst, ε)
with ε = (machine epsilon) to guard against under-corrected peaks.
Peaks where β²_obs ≤ β²_inst are flagged in the result dict's
`'warnings'` list and excluded from regression.

#### 7.2.2 W-A

```python
def guided_warren_averbach(self, ref_peaks=None, ..., *,
                            phase=None, instrumental=None,
                            axis_projection=False):  # NEW v0.5.0
    # ... phase / ref_peaks dispatch as v0.3.0 ...

    # NEW v0.4.0: Stokes deconvolution requires measured pattern
    if instrumental is None:
        std = None
    elif isinstance(instrumental, InstrumentalStandard):
        std = instrumental
    elif isinstance(instrumental, InstrumentalProfile):
        raise ValueError(
            "Warren-Averbach Stokes deconvolution requires the measured "
            "standard pattern; pass an InstrumentalStandard, or call "
            "without instrumental= for uncorrected W-A. (Caglioti "
            "FWHM subtraction is mathematically equivalent to FWHM-only "
            "correction and is unprincipled for the higher-order "
            "Fourier coefficients W-A consumes.)")
    else:
        raise TypeError(...)
    # ... per-family loop, with Stokes deconvolution applied to
    # A_obs(L) → A_corr(L) when std is not None ...
```

The Stokes step is, for each family member peak:
1. Get sample profile Fourier coefficients A_obs(L) (already computed
   in v0.3.0).
2. Get standard profile Fourier coefficients A_inst(L) at the matching
   peak via `std.fourier_coefficients(d_ref, n_coeffs)`.
3. A_corr(L) = A_obs(L) / A_inst(L), with a damping floor:
   if A_inst(L) < threshold * A_inst(0), set A_corr(L) = 0 to suppress
   noise amplification at high L.
4. Continue with the v0.3.0 harmonic decomposition using A_corr.

#### 7.2.3 Scherrer (v0.4.0)

```python
def scherrer(self, K=None, shape=None, height_threshold=0.05,
             *, phase=None, instrumental=None):  # NEW v0.4.0
    fwhm, positions = estimate_fwhm_simple(
        self.two_theta, self.intensity, height_threshold)

    # NEW v0.4.0: phase filtering
    if phase is not None:
        ref_d = phase.get_ref_d(self.wavelength, ...)
        positions, fwhm = _filter_to_phase_peaks(positions, fwhm,
                                                  ref_d, self.wavelength,
                                                  tolerance_d=0.03)

    # NEW v0.4.0: instrumental correction
    if instrumental is not None:
        inst_profile = (instrumental.caglioti_fit()
                        if isinstance(instrumental, InstrumentalStandard)
                        else instrumental)
        fwhm_inst = np.array([inst_profile.fwhm_at(tt) for tt in positions])
        fwhm = np.sqrt(np.maximum(fwhm**2 - fwhm_inst**2, 0))

    sizes = scherrer(fwhm, positions, self.wavelength, K=K, shape=shape)
    # ... existing return shape ...
```

`modified_scherrer` gets the same kwargs, same dispatch.

### 7.3 Crystallite size distributions (v0.4.0)

Per-family in the W-A result `families` list, each entry gains:

```python
{
  # ... existing v0.3.0 keys: base_hkl, orders, d_spacings,
  #     fwhm_values, A_size, L, mean_sq_strain, crystallite_size,
  #     rms_strain, A_size_r2, has_overlap ...

  'size_distribution': {
      'lognormal': {
          'D_median': float,           # angstroms
          'sigma': float,              # log-space stddev
          'D_mean_volume': float,      # volume-weighted mean
          'D_mean_area': float,        # area-weighted mean
          'fit_r2': float,
          'cov': np.ndarray,           # 2x2 parameter covariance
      },
      'normal': {
          'D_mean': float,
          'sigma': float,
          'fit_r2': float,
          'cov': np.ndarray,
      },
      'method': 'curve_fit',
      'initial_guess': 'moments',
      'n_valid_L': int,
  },
}
```

`'size_distribution'` is `None` when `n_valid_L < 4` (insufficient
points for a 2-parameter nonlinear fit with residual DOF > 1).

Mathematical basis: A_size(L) for a lognormal column-length
distribution has a closed-form expression involving the complementary
error function (Krill & Birringer 1998; Langford, Louër, Scardi
2000). Method-of-moments estimators from A_size'(0) and ∫A_size dL
provide initial values; `scipy.optimize.curve_fit` refines them.
Normal fit follows the same machinery with a different basis function.

### 7.4 Direction-resolved size reporting (v0.5.0)

Top-level keys added to the W-A result dict:

```python
{
  'families': [...],                                   # unchanged from v0.4.0

  # NEW v0.5.0: aggregation across families
  'sizes_by_hkl': {(2,0,0): 350.0, (0,2,0): 380.0, (0,0,2): 410.0},
  'size_distributions_by_hkl': {(2,0,0): {'lognormal': {...},
                                           'normal': {...}}, ...},

  # NEW v0.5.0: optional crystal-axis projection
  'sizes_by_axis': {                                   # only if axis_projection=True
      'a': 348.0,
      'b': 380.0,
      'c': 411.0,
      'method': 'direct',                              # only value emitted in v0.5.0
      'lattice_system': 'orthorhombic',
  },

  # ... v0.4.0 top-level keys unchanged ...
}
```

`sizes_by_hkl` and `size_distributions_by_hkl` are always present in
v0.5+ W-A results — empty dicts (`{}`) when no families were resolved,
populated otherwise. `sizes_by_axis` is present only when
`axis_projection=True` was requested.

`axis_projection=True` is honoured only for cubic / hexagonal /
tetragonal / orthorhombic lattices, where each family direction maps
unambiguously to a single crystal axis (cubic: all → isotropic, so
`a == b == c`; tetragonal / hexagonal: in-plane axes equivalent by
symmetry, so `a == b ≠ c`; orthorhombic: three distinct values).
The dict always reports three named axes for uniformity; equality is
by lattice symmetry, not by happy coincidence. For lower-symmetry
systems (monoclinic, triclinic) the projection is an under-determined
linear inversion (`'lsq_inversion'`) and is out of scope for v0.5.0;
the kwarg raises `ValueError` with the documented reason. The
`'method'` field is always `'direct'` in v0.5.0; `'lsq_inversion'` is
reserved for a future tag and not emitted here. Lattice-system
detection uses `Phase.structure.lattice` properties (available via
pymatgen).

### 7.5 Modified Williamson-Hall plug-in (v0.5.0)

`guided_williamson_hall` accepts `contrast_factor=None | callable`.
When non-`None`, the W-H abscissa changes from `4·sinθ` (or its
reciprocal-space variant) to `4·sinθ·√C̄_hkl`, where C̄_hkl is the
contrast factor returned by `contrast_factor(h, k, l, lattice)`.

Callable contract:

```python
def contrast_factor(h: int, k: int, l: int,
                    lattice: 'pymatgen.core.lattice.Lattice') -> float:
    """Return the dislocation contrast factor for reflection (h,k,l).
    Caller supplies a function appropriate to their crystal system,
    elastic constants, dislocation type, and slip system."""
```

`examples/modified_wh_cubic.py` provides the Ungár-Tichy (1999)
analytical cubic form as a worked template:

```python
def ungar_tichy_cubic(q: float):
    """Returns a contrast_factor callable for cubic crystals
    parameterised by a single elastic-anisotropy parameter q.
    Reference: Ungár & Tichy (1999), Phys. Stat. Sol. A 171, 425-434."""
    def C_hkl(h, k, l, lattice):
        H_sq = (h*h*k*k + k*k*l*l + l*l*h*h) / (h*h + k*k + l*l)**2
        Ch00 = 0.30  # representative; user-tunable
        return Ch00 * (1 - q * H_sq)
    return C_hkl
```

The package itself does not commit to specific Ch00 / q values for any
crystal system. The callable is the user's responsibility; the example
shows the shape.

### 7.6 YAML config + CLI (v1.0.0)

#### 7.6.1 Schema (Pydantic v2 models in `config_schema.py`)

```python
class SampleConfig(BaseModel):
    pattern: str                           # path to .xy / .csv / .dat
    wavelength: float
    name: str = ''

class PhaseConfig(BaseModel):
    name: str
    cif: str | None = None                 # path to .cif
    lattice_params: dict | None = None     # {a, b, c, alpha, beta, gamma,
                                           #  species, coords}
    # exactly one of cif / lattice_params

class InstrumentalConfig(BaseModel):
    type: Literal['standard', 'profile', 'none'] = 'none'
    cif: str | None = None
    pattern: str | None = None
    U: float | None = None
    V: float | None = None
    W: float | None = None
    # 'standard' requires cif + pattern
    # 'profile' requires U + V + W
    # 'none' uses no fields

class OutputConfig(BaseModel):
    results: str | None = None             # JSON output path
    figures: str | None = None             # PNG output dir

class ExperimentConfig(BaseModel):
    sample: SampleConfig
    phases: list[PhaseConfig] = []
    instrumental: InstrumentalConfig = InstrumentalConfig()
    methods: list[Literal['wh', 'wa', 'pdf', 'scherrer']] = \
        ['wh', 'wa', 'pdf', 'scherrer']
    wh: dict = Field(default_factory=dict)         # method kwargs
    wa: dict = Field(default_factory=dict)
    pdf: dict = Field(default_factory=dict)
    scherrer: dict = Field(default_factory=dict)
    output: OutputConfig = OutputConfig()
```

The four method-kwarg dicts (`wh`, `wa`, `pdf`, `scherrer`) live
directly on `ExperimentConfig` to match the YAML structure shown in
Section 1.6.1 above (`wh: {n_sigma: 3.0, ...}` at root level rather
than nested under `method_kwargs:`). The unknown-keys behaviour for
each dict is `Extra.allow` so users can pass any kwarg that
`run_all`'s per-method dispatch accepts — Pydantic does not validate
the inner kwarg names against the method signatures (validation
happens at the `run_all` call site instead, with Python's standard
TypeError if a kwarg is unknown).

#### 7.6.2 CLI (`xrd_profile/cli.py`)

```bash
xrd_profile run experiment.yaml
xrd_profile init [path/to/scaffold.yaml]
xrd_profile --version
xrd_profile --help
```

`run` flow: load YAML → validate via `ExperimentConfig` → resolve
relative paths → build `Phase` objects via `Phase.from_cif` or
`Phase.from_lattice_params` → build `InstrumentalStandard` /
`InstrumentalProfile` / `None` → load pattern via `numpy.loadtxt`
(extension dispatch) → construct `XRDProfile` → call
`profile.run_all(methods=cfg.methods, phases=phases,
instrumental=instrumental, wh=cfg.wh, wa=cfg.wa, pdf=cfg.pdf,
scherrer=cfg.scherrer)` → serialise
results to JSON (custom encoder converts numpy arrays to lists, dicts
keyed by tuples are converted to lists of {key, value} pairs) →
generate figures if `output.figures` is set.

`init` writes a starter YAML to the given path (default
`./experiment.yaml`) with comments documenting each field.

argparse handles the verb dispatch; no Click dependency needed
(keeping the optional-extras footprint to PyYAML + Pydantic only).

## 8. Examples and tests

### 8.1 Example scripts (cumulative through v1.0)

- `examples/multi_phase_olivine.py` — v0.3.0, unchanged.
- `examples/lab_lunar_meteorite.py`, `synchrotron_low_shock.py`,
  `synchrotron_high_shock.py` — v0.3.0, unchanged.
- `examples/legacy/*.py` — v0.3.0 originals, unchanged.
- `examples/instrumental_correction.py` — NEW v0.4.0. Shows
  `InstrumentalStandard.from_cif_and_pattern(...)` with the synthetic
  LaB6 fixture, runs guided W-A with and without correction, prints
  the size shift.
- `examples/modified_wh_cubic.py` — NEW v0.5.0 (Section 7.5).
- `examples/cli_experiment.yaml` — NEW v1.0.0. Example YAML config
  exercising the full pipeline.

### 8.2 Tests (cumulative through v1.0)

```
tests/
├── fixtures/
│   ├── tirhert_subset.xy                       # v0.3.0
│   ├── forsterite_for_test.cif                 # v0.3.0
│   ├── synthetic_lab6.xy                       # NEW v0.4.0
│   ├── golden_v0.2.0_results.json              # v0.3.0
│   ├── golden_v0.3.0_results.json              # NEW v0.4.0
│   ├── golden_v0.4.0_results.json              # NEW v0.5.0
│   └── golden_v0.5.0_results.json              # NEW v1.0.0
├── test_conversions.py                         # v0.3.0
├── test_pdf.py                                 # v0.3.0
├── test_phases.py                              # v0.3.0
├── test_run_all.py                             # v0.3.0
├── test_scherrer.py                            # v0.3.0
├── test_backward_compat.py                     # v0.3.0; tiered each release
├── test_instrumental.py                        # NEW v0.4.0
├── test_size_distributions.py                  # NEW v0.4.0
├── test_anisotropic.py                         # NEW v0.5.0
├── test_modified_wh.py                         # NEW v0.5.0
└── test_cli.py                                 # NEW v1.0.0
```

#### 8.2.1 `test_instrumental.py` (~12 tests)

- `InstrumentalStandard.from_cif_and_pattern` constructs correctly
  given the synthetic LaB6 fixture.
- Caglioti fit recovers the synthesis parameters (`U`, `V`, `W`)
  within 5%.
- `caglioti_fwhm_at(2θ)` is monotonic in 2θ above the V-term-induced
  minimum, as the polynomial dictates.
- `fourier_coefficients(d, n_coeffs)` returns arrays of the requested
  length.
- `InstrumentalProfile.from_standard(std)` round-trips: the returned
  profile's `fwhm_at` matches the standard's measured FWHMs.
- `InstrumentalProfile.from_registry('nonexistent')` raises clear
  `KeyError` (registry empty in v0.4.0).
- `to_json` / `from_json` round-trip both classes.
- `guided_williamson_hall(phase=phase, instrumental=std)` runs end to
  end on `tirhert_subset.xy` and produces a finite crystallite size.
- `guided_warren_averbach(phase=phase, instrumental=profile)` raises
  `ValueError` with the documented "Stokes requires measured pattern"
  message.
- `guided_warren_averbach(phase=phase, instrumental=std)` runs end to
  end and produces a finite crystallite size.
- Stokes deconvolution reduces the apparent crystallite size compared
  to uncorrected W-A on the same data (qualitative sanity).
- `scherrer(phase=phase, instrumental=std)` filters peaks correctly
  and produces the expected count.

#### 8.2.2 `test_size_distributions.py` (~8 tests)

- Synthetic A_size from a known lognormal `(D_median, σ)` is fit by
  the package within 5% of the synthesis parameters.
- Synthetic A_size from a known normal distribution is fit within 5%
  of synthesis parameters.
- Family with `n_valid_L < 4` returns `'size_distribution': None`.
- `'method'`, `'initial_guess'`, `'n_valid_L'` keys are present.
- Lognormal `D_mean_volume > D_median > D_mean_area`-ish ordering
  for σ > 0 (sanity).
- Normal fit on a clearly-lognormal distribution gives lower R²
  than the lognormal fit (sanity).
- Distribution fit does not modify the existing `crystallite_size`
  scalar.
- Frozen v0.3.0 W-A output (no distribution key) compared against
  v0.4.0 W-A output (distribution key present): all v0.3.0 keys
  retain their numerical values (key-subset value-equality).

#### 8.2.3 `test_anisotropic.py` (~6 tests)

- `sizes_by_hkl` has one entry per family in `families` and the
  values match the per-family `crystallite_size`.
- `size_distributions_by_hkl` mirrors the per-family
  `size_distribution`.
- `axis_projection=True` on a cubic phase produces
  `{'a': x, 'b': x, 'c': x}` with all three equal (within float
  tolerance).
- `axis_projection=True` on an orthorhombic phase produces three
  distinct values matching the (100), (010), (001) family sizes.
- `axis_projection=True` on a hexagonal phase produces `a`, `b`, `c`
  with `a == b ≠ c`.
- `axis_projection=True` on a monoclinic phase raises `ValueError`
  with the documented message.

#### 8.2.4 `test_modified_wh.py` (~5 tests)

- `contrast_factor=None` produces W-H output value-identical to v0.4.0
  on every key.
- A constant callable returning `1.0` everywhere produces W-H output
  value-identical to v0.4.0 on every key (mathematical equivalence).
- A non-trivial callable (e.g., `lambda h,k,l,lat: 0.3 * (h+k+l)`)
  shifts the W-H slope as expected.
- The Ungár-Tichy callable from `examples/modified_wh_cubic.py`
  imports cleanly and runs end to end.
- Invalid callable signature raises clear `TypeError` at first
  invocation (not a deferred numpy error).

#### 8.2.5 `test_backward_compat.py` (tiered)

- `golden_v0.2.0`: as in v0.3.0. Exercises the array-based public API.
- `golden_v0.3.0`: NEW v0.4.0. Exercises `Phase`, `run_all`, Scherrer
  shape table; all v0.4+ kwargs left at their defaults
  (`instrumental=None`, etc.).
- `golden_v0.4.0`: NEW v0.5.0. Exercises v0.4.0 features
  (instrumental correction, size distributions); all v0.5+ kwargs left
  at their defaults (`axis_projection=False`, `contrast_factor=None`).
- `golden_v0.5.0`: NEW v1.0.0. Exercises v0.5.0 features.

**Comparison semantics — key-subset value-equality, not literal byte-identity.**
Each tier asserts that, for every key K present in the golden file,
the live result contains K with a numerically equal value (1e-10
tolerance for direct outputs, 1e-6 for derived quantities) on
`tirhert_subset.xy`. The live result *may* carry additional keys
introduced at later tiers — for example, a v0.5.0-built result run
against the v0.3.0 golden will pass the comparison even though the
v0.5.0 result has the new top-level `'sizes_by_hkl'` and
`'size_distributions_by_hkl'` keys that the v0.3.0 golden never knew
about. The comparison loops over the golden's keys, never the live
result's.

This semantics is what "strict additivity" actually means: existing
documented behaviour does not change; new behaviour is allowed to add
new fields. Adding a later tier never modifies an earlier tier's
golden file — that is the additivity contract.

`scripts/regenerate_goldens.py` accepts `--tier v0.X.0` to regenerate a
single tier; regeneration requires explicit reasoning in the commit
message, same gate as v0.3.0. A regeneration that *adds keys* to the
golden file is the normal case for a new tier; a regeneration that
*changes existing key values* requires a load-bearing justification in
the commit message and counts as a numerical break (semver-major
caveat — Phase 2 commits to never doing this without an explicit user
decision).

#### 8.2.6 `test_cli.py` (~10 tests)

- `xrd_profile init` writes a valid YAML scaffold.
- The scaffold passes `ExperimentConfig` validation.
- `xrd_profile run examples/cli_experiment.yaml` returns exit 0
  on the bundled fixture.
- Missing-file YAML produces a clear validation error.
- Mutually-exclusive fields (e.g., `cif` + `lattice_params` on the
  same phase) produce a clear validation error.
- Output JSON loads cleanly and contains the expected top-level keys.
- Output figures are produced at the configured paths.
- `xrd_profile --version` prints the package version.
- `xrd_profile --help` prints the help text.
- Unknown verb produces a clear error.

### 8.3 Bundled fixture for the CLI tutorial

`examples/cli_experiment.yaml` references `examples/data/tirhert_subset.xy`
(already shipped in v0.3.0) and `examples/cifs/Forsterite.cif` /
`Pigeonite.cif` (also v0.3.0). No new bundled data needed for the CLI
example.

## 9. Documentation site (v1.0.0)

### 9.1 Tooling

mkdocs-material with the following plugin set:

- `mkdocstrings[python]` — auto-API reference from docstrings.
- `mkdocs-jupyter` — render the three tutorial notebooks inline.
- `mkdocs-include-markdown-plugin` — include `CHANGELOG.md` directly
  into `docs/changelog.md`.

Site source under `docs/`. Build output `site/` (gitignored). Live
preview via `mkdocs serve`.

### 9.2 Site structure

```
docs/
├── index.md                          # landing: install, 5-line example
├── quickstart.md                     # README quickstart, expanded
├── concepts/
│   ├── williamson-hall.md
│   ├── warren-averbach.md
│   ├── instrumental-correction.md
│   ├── pdf.md
│   └── scope.md                      # silicate-only positioning
├── tutorials/
│   ├── olivine-diogenite.ipynb
│   ├── pyroxene-eucrite.ipynb
│   └── quartz-calibration.ipynb
├── api/                              # mkdocstrings-generated
├── migration.md                      # v0.2 → v0.3 → v0.4 → v0.5 → v1.0
└── changelog.md                      # include of CHANGELOG.md
```

### 9.3 Hosting

GitHub Actions workflow at `.github/workflows/docs.yml` with
`on: workflow_dispatch` only — no automatic triggers. The workflow
builds the site and deploys to the `gh-pages` branch when manually
invoked. Pages serving stays disabled at the repo level until the
repo goes public.

When the repo goes public:
1. Flip `on: workflow_dispatch` → `on: push: branches: [main]`.
2. Enable Pages serving from `gh-pages` in repo settings.
3. Add the resulting URL to `pyproject.toml` `[project.urls]`
   `Documentation`.

Steps 1-3 are a single PR when the time comes. v1.0.0 ships in the
"private repo, infrastructure ready" state.

### 9.4 Scope page (`docs/concepts/scope.md`)

Explicit positioning paragraph: v1.0 is validated for
mesodesmic-bonded silicate phases (Si-O frameworks: olivine,
pyroxene, plagioclase, quartz, feldspathoids). Ionic systems
(oxides like hematite, magnetite, ilmenite) and metallic systems
work mechanically through the same `Phase` API but are not
included in the package's published validation; results require
independent verification before publication. Plans to extend
validation to ionic and metallic systems are documented as future
work.

The same statement appears (abridged) in `examples/cifs/SOURCES.md`
next to `Hematite.cif`, in the README, and as a docstring note on
`Phase.from_cif`.

## 10. PyPI release path (v1.0.0 infrastructure; public upload deferred)

### 10.1 What lands in v1.0.0

- `pyproject.toml` complete:
    - `[project]` with `name`, `version`, `description`, `readme`,
      `requires-python`, `license`, `authors`, `keywords`,
      `classifiers`.
    - `[project.urls]`: `Homepage = "https://matthewizawa.github.io/bio/"`
      (author bio); `Documentation` left unset until public docs go
      live; `Repository` left unset until repo public.
    - `[project.optional-dependencies]` (Section 6.4).
    - `[project.scripts]` for the CLI.
    - `[build-system]` with `setuptools` (current) or `hatchling`
      (decision made during implementation; either works).
- `scripts/release.sh`:
    1. `python -m build` → produces `dist/*.tar.gz` and `dist/*.whl`.
    2. `twine check dist/*` → metadata + README rendering.
    3. `twine upload --repository testpypi dist/*` (interactive
       authentication via `~/.pypirc` or `TWINE_*` env vars).
    4. `pip install --index-url https://test.pypi.org/simple/
       --extra-index-url https://pypi.org/simple/ xrd_profile` in a
       fresh virtualenv → smoke-test the install.
    5. Print the URL of the TestPyPI listing.
- The script does **not** upload to public PyPI. That is a separate,
  manual command (`twine upload dist/*`) when the user decides.

### 10.2 What does not land in v1.0.0

- Public PyPI listing. The `xrd_profile` name on PyPI is reserved as
  part of the TestPyPI dry-run (TestPyPI uploads do not reserve PyPI
  names; if name reservation matters, register the name on PyPI
  separately as a 0.0.0 placeholder — out of scope for v1.0.0 unless
  the user requests).
- Public documentation site. See Section 9.3.
- Public repo visibility. See Section 11.

### 10.3 Sequence to public release (post-v1.0.0)

1. JAC manuscript decision arrives.
2. Repo visibility decision: keep private, or split package to its
   own public repo, or make the existing repo public.
3. Flip the docs workflow to automatic deploy (Section 9.3).
4. Run `twine upload dist/*` against public PyPI.
5. Tag a `v1.0.0-public` release on GitHub (or rely on the v1.0.0
   tag created at v1.0.0 ship time).

This sequence is documented in `docs/migration.md` as a forward-looking
note; the user makes each decision when ready.

## 11. Phase 2 acceptance criteria (per tag)

### 11.1 v0.4.0

1. `InstrumentalStandard.from_cif_and_pattern` works end-to-end on
   the synthetic LaB6 fixture.
2. `caglioti_fit()` recovers the synthesis Caglioti coefficients
   within 5%.
3. `guided_williamson_hall(phase=, instrumental=)` works for both
   `InstrumentalStandard` and `InstrumentalProfile` arguments.
4. `guided_warren_averbach(phase=, instrumental=)` works for
   `InstrumentalStandard`; raises clear `ValueError` for
   `InstrumentalProfile`.
5. `scherrer(phase=, instrumental=)` works.
6. W-A result has `'size_distribution'` key on each family with
   `n_valid_L >= 4`; `None` otherwise.
7. `test_instrumental.py`, `test_size_distributions.py` pass.
8. `test_backward_compat.py` passes both `golden_v0.2.0` and
   `golden_v0.3.0` tiers.
9. CHANGELOG entry for v0.4.0.
10. README updated with instrumental-correction quickstart.
11. `git tag v0.4.0` on the v0.4.0-complete commit.

### 11.2 v0.5.0

1. `sizes_by_hkl` and `size_distributions_by_hkl` always present on
   v0.5+ W-A results: empty dicts when `len(families) == 0`,
   populated otherwise.
2. `axis_projection=True` produces correct values for cubic,
   hexagonal, tetragonal, and orthorhombic test phases.
3. `axis_projection=True` raises `ValueError` for monoclinic,
   triclinic.
4. `contrast_factor=None` produces W-H output value-identical to
   v0.4.0 on every v0.4.0 key.
5. The Ungár-Tichy cubic callable from
   `examples/modified_wh_cubic.py` runs end-to-end.
6. `test_anisotropic.py`, `test_modified_wh.py` pass.
7. `test_backward_compat.py` passes all of `golden_v0.2.0`,
   `golden_v0.3.0`, `golden_v0.4.0`.
8. CHANGELOG entry for v0.5.0.
9. README updated.
10. `git tag v0.5.0` on the v0.5.0-complete commit.

### 11.3 v1.0.0

1. `xrd_profile run examples/cli_experiment.yaml` exits 0 on the
   bundled fixture.
2. `xrd_profile init` writes a valid YAML scaffold that passes
   `ExperimentConfig` validation.
3. `xrd_profile --version` prints `1.0.0`.
4. `pip install -e .[cli]` succeeds and the `xrd_profile` console
   script is on PATH.
5. `mkdocs build` succeeds with no warnings.
6. `mkdocs serve` renders the three tutorial notebooks correctly.
7. `python -m build` produces `dist/*.tar.gz` and `dist/*.whl`
   without error.
8. `twine check dist/*` passes.
9. `scripts/release.sh` completes the TestPyPI dry-run successfully
   (manual run; not in CI).
10. `test_cli.py` passes.
11. `test_backward_compat.py` passes all four golden tiers.
12. `CONTRIBUTING.md` and `.github/` templates present.
13. CHANGELOG entry for v1.0.0.
14. `git tag v1.0.0` on the v1.0.0-complete commit.

## 12. Decisions deferred to v2.x

- **Typed-dataclass results object.** Field set is still mutating
  through Phase 2; lock after real-world v1.0 usage informs which
  fields stabilise.
- **pymatgen as hard dependency.** Promotion decision once PyPI usage
  data shows the `[cif]`-vs-bare-install split.
- **Built-in instrumental-profile registry.** v0.4.0 ships
  infrastructure; populating registry entries for the Misasa Rigaku,
  Winnipeg Bruker, and Diamond I11 is a separate calibration
  deliverable on its own timeline.
- **Contrast-factor library.** The plug-in interface is the v0.5.0
  commitment. A library of cubic / hexagonal / orthorhombic factors
  could ship in v2.x; silicate factors require new research.
- **Public PyPI upload and public repo visibility.** Coordinated with
  the JAC manuscript decision (Section 10.3).
- **Validation extension to ionic and metallic systems.** Future work
  documented in `docs/concepts/scope.md`.

## 13. Open questions

None at this time. All clarifications from the 2026-05-05 Phase 2
brainstorming session have been incorporated.
