# xrd_profile

A Python toolkit for quantitative analysis of powder X-ray diffraction peak profiles.

## Features

- **Bragg's law conversions**: 2-theta, d-spacing, Q, K = 1/d (wavelength-independent)
- **Reference-guided peak detection**: use computed mineral reference positions with automated noise cutoff and zero-point offset correction
- **Williamson-Hall analysis**: conventional and reciprocal-space (DeltaK vs K) formulations
- **Warren-Averbach analysis**: Fourier decomposition with harmonic peak families, adaptive window, Tukey tapering, and quality filtering
- **Scherrer equation**: standard per-peak and modified log-linear regression
- **Pair distribution function**: sine Fourier transform PDF with iterative Chebyshev background subtraction, Lorch modification, peak detection, and Gaussian first-shell fitting; plus FFT-based PDF with envelope function crystallite size estimation
- **XRDProfile class**: unified interface wrapping all methods with plotting utilities

## Installation

```bash
cd xrd_profile
pip install -e .
```

## Quick start

```python
import numpy as np
from xrd_profile import XRDProfile

# Load diffraction data (2-theta in degrees, intensity in counts)
data = np.loadtxt('my_pattern.xy', comments='#')
two_theta, intensity = data[:, 0], data[:, 1]

# Create profile object (specify your wavelength)
profile = XRDProfile(two_theta, intensity, wavelength=1.5406,
                     sample_name='My Sample')

# Run unguided analysis
results = profile.full_analysis()
print(f"W-H crystallite size: {results['williamson_hall']['crystallite_size']:.0f} A")
print(f"Scherrer mean size: {results['scherrer']['mean_size']:.0f} A")

# Plot in d-spacing
profile.plot_pattern(x_axis='d_spacing')
```

### Reference-guided analysis

For multi-phase samples, use computed reference peak positions from crystal structure data to guide the analysis:

```python
from pymatgen.core.structure import Structure
from pymatgen.analysis.diffraction.xrd import XRDCalculator
from xrd_profile import XRDProfile, two_theta_to_d

# Build reference peak list from CIF
structure = Structure.from_file('anorthite.cif')
xrd_calc = XRDCalculator(wavelength='CuKa')
pattern = xrd_calc.get_pattern(structure, two_theta_range=(5, 90))

ref_peaks = []
for i in range(len(pattern.x)):
    hkl = pattern.hkls[i][0]['hkl']
    d = two_theta_to_d(pattern.x[i], 1.5406)
    ref_peaks.append({
        'd': float(d), 'two_theta': float(pattern.x[i]),
        'intensity': float(pattern.y[i]),
        'h': hkl[0], 'k': hkl[1], 'l': hkl[2],
    })

# d-spacings sorted by intensity for W-H
ref_d = sorted([p['d'] for p in ref_peaks if p['intensity'] >= 3.0],
               reverse=True)

# Guided Williamson-Hall
wh = profile.guided_williamson_hall(ref_d, n_sigma=3.0)

# Guided Warren-Averbach (needs full hkl for harmonic families)
wa = profile.guided_warren_averbach(ref_peaks, n_sigma=3.0)
```

### Pair distribution function analysis

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

## Dependencies

- numpy
- scipy
- matplotlib
- pymatgen (optional, for computing reference peak positions from CIF files)

## Attribution

This package incorporates and adapts code from the following sources:

### crystallite_size_calculator

The Williamson-Hall, Warren-Averbach, Scherrer, pair distribution function, and envelope function implementations are adapted from [crystallite_size_calculator](https://github.com/bafgreat/crystallite_size_calculator) by Dinga Wonanke. The original package required CIF file inputs via ase and pymatgen; this adaptation works directly with experimental 2-theta/intensity arrays and removes those dependencies from the core library.

### Methods references

- **Williamson-Hall method**: Williamson, G. K. & Hall, W. H. (1953). X-ray line broadening from filed aluminium and wolfram. *Acta Metallurgica*, 1, 22-31.
- **Modified Williamson-Hall**: Das Bakshi, S., Sinha, D. & Ghosh Chowdhury, S. (2018). Anisotropic broadening of XRD peaks of alpha'-Fe: Williamson-Hall and Warren-Averbach analysis. *Materials Characterization*, 142, 144-153.
- **Warren-Averbach method**: Warren, B. E. & Averbach, B. L. (1950). The effect of cold-work distortion on X-ray patterns. *Journal of Applied Physics*, 21, 595-599.
- **Scherrer equation**: Scherrer, P. (1918). Bestimmung der Grosse und der inneren Struktur von Kolloidteilchen mittels Rontgenstrahlen. *Nachrichten von der Gesellschaft der Wissenschaften zu Gottingen*, 26, 98-100.
- **Envelope function approach**: Gesing, T. M. & Robben, L. (2024). Accurate determination of crystallite sizes and crystallite size distributions. *Journal of Applied Crystallography*, 57, 1466-1476.

### Original contributions

The following features are original to this package:

- Reciprocal-space Williamson-Hall formulation (DeltaK vs K) for wavelength-independent cross-source comparison
- Reference-guided peak detection using computed mineral reference positions with noise cutoff
- Automated zero-point offset estimation
- Improved Warren-Averbach Fourier coefficient extraction with adaptive windowing, Tukey tapering, uniform s-grid interpolation, overlap detection, and quality filtering
- Multi-phase analysis capability (separate guides for different mineral phases in the same pattern)
- Enhanced PDF pipeline: iterative Chebyshev polynomial background subtraction in Q-space, sine Fourier transform with optional Lorch modification, PDF peak detection and Gaussian first-shell fitting

## License

MIT License. See LICENSE file.

## Citation

If you use this package in published research, please cite:

> Izawa, M. R. M. (2026). xrd_profile: XRD peak profile analysis toolkit (v0.2.0). https://github.com/matthewizawa/xrd_profile

and acknowledge the crystallite_size_calculator package by Wonanke (see Attribution above).
