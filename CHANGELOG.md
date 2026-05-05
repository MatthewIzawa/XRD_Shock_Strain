# Changelog

All notable changes to xrd_profile are documented here. Format
follows [Keep a Changelog](https://keepachangelog.com/), versioning
follows [SemVer](https://semver.org/).

## [0.3.0] — 2026-05-05

### Added
- `Phase` class with `from_cif()` and `from_lattice_params()`
  constructors. Reference peaks are generated via
  `phase.get_ref_peaks(wavelength)` and `phase.get_ref_d(wavelength)`.
  Removes the per-script `build_ref` helper duplication.
- `phase=` keyword argument on `XRDProfile.guided_williamson_hall()`
  and `guided_warren_averbach()`. Mutually exclusive with the
  legacy `ref_d=` / `ref_peaks=` arguments.
- `instrumental=` keyword argument reserved on the same methods and
  on `XRDProfile.run_all()` (raises `NotImplementedError` if used;
  Phase 2 wires it up).
- `XRDProfile.run_all(methods=[...], phases=[...], wh={...},
  wa={...}, pdf={...}, scherrer={...})` convenience helper.
- `shape` keyword argument on `scherrer()` and `modified_scherrer()`,
  with a `SCHERRER_K_FOR_SHAPE` lookup table covering spherical
  (0.94), cubic (0.83), cylindrical (1.84), and platey (1.0)
  crystallite shapes.
- Five example CIFs (forsterite, anorthite, pigeonite, quartz,
  hematite) in `examples/cifs/`, with provenance documented in
  `SOURCES.md`.
- Bundled diffraction pattern at `examples/data/tirhert_subset.xy`
  (~460 KB Tirhert eucrite subset, I11 synchrotron).
- New canonical demo `examples/multi_phase_olivine.py`.
- Numerical regression test (`tests/test_backward_compat.py`)
  against frozen v0.2.0 outputs, enforcing the strict-additive
  policy across Phase 1 changes.
- `[cif]` optional dependency extra: `pip install xrd_profile[cif]`
  installs `pymatgen` for the Phase API.

### Changed
- `scherrer()` and `modified_scherrer()`: `K` default changed from
  `0.9` to `None` sentinel. When neither `K` nor `shape` is provided,
  behaviour is identical to v0.2.0 (`K=0.9`).
- Existing example scripts (`lab_lunar_meteorite.py`,
  `synchrotron_low_shock.py`, `synchrotron_high_shock.py`)
  rewritten to use the `Phase` API. Verbatim v0.2.0 versions
  preserved at `examples/legacy/`.

### Compatibility
- All v0.2.0 array-based public API calls continue to work
  unchanged. Numerical results are byte-for-byte identical when no
  v0.3.0 features are invoked (regression test enforces this).

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
