# xrd_profile

A Python toolkit for quantitative analysis of powder X-ray diffraction peak profiles.

## Features

- **Bragg's law conversions**: 2-theta, d-spacing, Q, K = 1/d (wavelength-independent)
- **`Phase` abstraction**: load any crystal structure from a CIF file or inline lattice parameters; reuse across analyses (new in v0.3.0)
- **Reference-guided peak detection**: automated noise cutoff and zero-point offset correction
- **Williamson-Hall analysis**: conventional and reciprocal-space (ΔK vs K) formulations
- **Warren-Averbach analysis**: Fourier decomposition with harmonic peak families, adaptive window, Tukey tapering, and quality filtering
- **Scherrer equation**: standard per-peak and modified log-linear regression, with `K` and crystallite shape factor selection (new in v0.3.0)
- **Pair distribution function**: sine Fourier transform PDF with iterative Chebyshev background subtraction, Lorch modification, peak detection, and Gaussian first-shell fitting
- **Bundled `run_all()` helper**: configurable subset of methods on configurable list of phases (new in v0.3.0)
- **`XRDProfile` class**: unified interface wrapping all methods with plotting utilities

## Installation

Recommended (includes `pymatgen` for CIF loading and `Phase.from_lattice_params`):

```bash
pip install -e ".[cif]"
```

Minimum (no `Phase`; the array-based public API works without `pymatgen`, useful when reference d-spacings come from Rietveld output or another source):

```bash
pip install -e .
```

## Quick start

```python
import numpy as np
from xrd_profile import XRDProfile, Phase

# Load diffraction data (2-theta in degrees, intensity in counts)
data = np.loadtxt('my_pattern.xy', comments='#')
two_theta, intensity = data[:, 0], data[:, 1]

profile = XRDProfile(two_theta, intensity, wavelength=1.5406,
                     sample_name='My Sample')

# Build a phase from a CIF
olivine = Phase.from_cif('Forsterite.cif', name='Olivine')

# Run a single guided analysis
wh = profile.guided_williamson_hall(phase=olivine, n_sigma=3.0)
print(f"W-H crystallite size: {wh['crystallite_size']:.0f} A")

# Or run a bundled multi-method, multi-phase analysis
results = profile.run_all(
    methods=['wh', 'wa', 'pdf', 'scherrer'],
    phases=[olivine],
    wh={'n_sigma': 3.0, 'tolerance_d': 0.02},
    wa={'tolerance_d': 0.02},
)
print(results['wh']['Olivine']['crystallite_size'])
print(results['pdf']['Q_max'])
```

See `examples/multi_phase_olivine.py` for a complete walk-through using the bundled CIFs and data fixture.

## Phases

A `Phase` wraps a `pymatgen.Structure` and provides reference-peak generation. Construct from a CIF or from inline lattice parameters and atomic coordinates:

```python
from xrd_profile import Phase

# From a CIF
quartz = Phase.from_cif('Quartz.cif', name='Quartz')

# From inline lattice parameters (useful when refined values are at hand)
anorthite = Phase.from_lattice_params(
    8.18, 12.88, 7.11, 93.5, 116.1, 90.4,
    species=['Ca','Al','Al','Si','Si','O','O','O','O','O','O','O','O'],
    coords=[...],
    name='Anorthite',
)

# Pass to any guided method
wh = profile.guided_williamson_hall(phase=quartz)
wa = profile.guided_warren_averbach(phase=anorthite)
```

Bundled example CIFs are in `examples/cifs/` (forsterite, anorthite, pigeonite, quartz, hematite) with full provenance in `examples/cifs/SOURCES.md`.

## Lower-level array-based API

Users who already have d-spacings (e.g., from Rietveld output) can bypass `Phase` and pass arrays directly:

```python
ref_d = np.array([3.20, 3.18, 3.65, 4.04, 6.41])  # from Rietveld
wh = profile.guided_williamson_hall(ref_d=ref_d, n_sigma=3.0)
```

This API is foundational, supported indefinitely, and does not require `pymatgen`.

## Pair distribution function analysis

The enhanced PDF pipeline uses iterative Chebyshev background subtraction in Q-space and a sine Fourier transform with optional Lorch modification to suppress termination ripples:

```python
from xrd_profile import XRDProfile

profile = XRDProfile(two_theta, intensity, wavelength=0.826517,
                     sample_name='Millbillillie')

# Compute PDF with Chebyshev background and Lorch modification
r, G_r, Q_max = profile.compute_pdf_sine(cheby_order=20, lorch=True)

# Measure all PDF peaks
peaks = profile.measure_pdf_peaks(min_r=1.0, max_r=15.0)
for p in peaks[:5]:
    print(f"  r = {p['r']:.3f} A, FWHM = {p['fwhm']:.3f} A")

# Gaussian fit to the first coordination shell
r_peak, fwhm = profile.fit_first_pdf_peak()
print(f"First shell: r = {r_peak:.3f} A, FWHM = {fwhm:.4f} A")
```

The functions are also available standalone for use outside the `XRDProfile` class:

```python
from xrd_profile import (chebyshev_background, compute_pdf_sine,
                          measure_pdf_peaks, fit_first_pdf_peak)

r, G_r, Q_max = compute_pdf_sine(two_theta, intensity, wavelength,
                                  cheby_order=15, lorch=True)
peaks = measure_pdf_peaks(r, G_r)
r_pk, fwhm = fit_first_pdf_peak(r, G_r)
```

## Scherrer K and shape factor (v0.3.0)

The Scherrer constant `K` and a crystallite shape factor are exposed as kwargs:

```python
from xrd_profile import scherrer, SCHERRER_K_FOR_SHAPE

# Default (K=0.9) — backward-compatible with v0.2.0
sizes = scherrer(fwhm, two_theta, wavelength=1.5406)

# Use a shape preset
sizes = scherrer(fwhm, two_theta, 1.5406, shape='spherical')   # K=0.94
sizes = scherrer(fwhm, two_theta, 1.5406, shape='cylindrical') # K=1.84

# Or set K explicitly (wins over shape if both are passed)
sizes = scherrer(fwhm, two_theta, 1.5406, K=1.0)
```

`SCHERRER_K_FOR_SHAPE` is exported for introspection: `{'spherical': 0.94, 'cubic': 0.83, 'cylindrical': 1.84, 'platey': 1.0}`.

## Dependencies

- numpy
- scipy
- matplotlib
- pymatgen (optional, installed via the `[cif]` extra; required for `Phase.from_cif` and `Phase.from_lattice_params`)

## Attribution

This package incorporates and adapts code from the following sources:

### crystallite_size_calculator

The Williamson-Hall, Warren-Averbach, Scherrer, pair distribution function, and envelope function implementations are adapted from [crystallite_size_calculator](https://github.com/bafgreat/crystallite_size_calculator) by Dinga Wonanke. The original package required CIF file inputs via ase and pymatgen; this adaptation works directly with experimental 2-theta/intensity arrays and removes those dependencies from the core library.

### Methods references

- **Williamson-Hall method**: Williamson, G. K. & Hall, W. H. (1953). X-ray line broadening from filed aluminium and wolfram. *Acta Metallurgica*, 1, 22-31.
- **Modified Williamson-Hall**: Das Bakshi, S., Sinha, D. & Ghosh Chowdhury, S. (2018). Anisotropic broadening of XRD peaks of α'-Fe: Williamson-Hall and Warren-Averbach analysis. *Materials Characterization*, 142, 144-153.
- **Warren-Averbach method**: Warren, B. E. & Averbach, B. L. (1950). The effect of cold-work distortion on X-ray patterns. *Journal of Applied Physics*, 21, 595-599.
- **Scherrer equation**: Scherrer, P. (1918). Bestimmung der Größe und der inneren Struktur von Kolloidteilchen mittels Röntgenstrahlen. *Nachrichten von der Gesellschaft der Wissenschaften zu Göttingen*, 26, 98-100.
- **Scherrer K shape factors**: Langford, J. I. & Wilson, A. J. C. (1978). Scherrer after sixty years: a survey and some new results in the determination of crystallite size. *Journal of Applied Crystallography*, 11, 102-113.
- **Envelope function approach**: Gesing, T. M. & Robben, L. (2024). Accurate determination of crystallite sizes and crystallite size distributions. *Journal of Applied Crystallography*, 57, 1466-1476.

### Original contributions

The following features are original to this package:

- Reciprocal-space Williamson-Hall formulation (ΔK vs K) for wavelength-independent cross-source comparison
- Reference-guided peak detection using computed mineral reference positions with noise cutoff
- Automated zero-point offset estimation
- Improved Warren-Averbach Fourier coefficient extraction with adaptive windowing, Tukey tapering, uniform s-grid interpolation, overlap detection, and quality filtering
- Multi-phase analysis capability (separate guides for different mineral phases in the same pattern)
- Enhanced PDF pipeline: iterative Chebyshev polynomial background subtraction in Q-space, sine Fourier transform with optional Lorch modification, PDF peak detection and Gaussian first-shell fitting
- `Phase` abstraction with `from_cif` / `from_lattice_params` constructors and library-owned reference-peak generation (v0.3.0)
- `XRDProfile.run_all()` configurable multi-method, multi-phase dispatcher (v0.3.0)

## License

MIT License. See LICENSE file.

## Citation

If you use this package in published research, please cite:

> Izawa, M. R. M. (2026). xrd_profile: XRD peak profile analysis toolkit (v0.3.0). https://github.com/matthewizawa/xrd_profile

and acknowledge the crystallite_size_calculator package by Wonanke (see Attribution above).
