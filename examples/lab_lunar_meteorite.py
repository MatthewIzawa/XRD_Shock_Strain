"""
Example: Laboratory XRD analysis of a lunar meteorite.

Analyses NWA 11182a (weakly shocked lunar meteorite) from a Rigaku
SmartLab with Cu K-alpha radiation (10-60 deg 2-theta, 0.02 deg steps).
Demonstrates reference-guided Williamson-Hall and Warren-Averbach
analysis using the v0.3.0 Phase API.

The legacy v0.2.0 version (with inline Structure + build_ref helper)
is preserved at examples/legacy/lab_lunar_meteorite.py.

Data: University of Winnipeg HOSERLab / Okayama University IPM.
"""

import numpy as np
from xrd_profile import XRDProfile
from xrd_profile.phases import Phase

# ── Instrument parameters ──
CU_KA = 1.54056  # Cu K-alpha wavelength (angstroms)

DATA_PATH = (r'C:\Users\Matthew Izawa\Documents\Conferences\LPSC 2025'
              r'\Lunar\XY files\nwa11182a.xy')

# ── Load diffraction data ──
# Format: space-separated 2-theta (deg), intensity (counts)
# Adapt this path to your own data
data = np.loadtxt(DATA_PATH, comments='#')
two_theta, intensity = data[:, 0], data[:, 1]

print(f"Loaded: {len(two_theta)} points, "
      f"{two_theta.min():.1f}-{two_theta.max():.1f} deg 2-theta")

# ── Create profile object ──
profile = XRDProfile(two_theta, intensity, wavelength=CU_KA,
                     sample_name='NWA 11182a')

# ── Step 1: Unguided analysis (quick overview) ──
print("\n--- Unguided analysis ---")
results = profile.full_analysis()
wh = results['williamson_hall']
print(f"Peaks detected: {wh['n_peaks']}")
print(f"W-H crystallite size: {wh['crystallite_size']:.0f} angstroms")
print(f"Scherrer mean size: {results['scherrer']['mean_size']:.0f} angstroms")

# ── Step 2: Build anorthite phase ──
# Anorthite (CaAl2Si2O8), P-1, RRUFF R060193 cell parameters
anorthite = Phase.from_lattice_params(
    8.1809, 12.881, 7.1101, 93.465, 116.11, 90.369,
    species=['Ca', 'Al', 'Al', 'Si', 'Si',
             'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O'],
    coords=[
        [0.269, 0.988, 0.086],  # Ca
        [0.507, 0.314, 0.621],  # Al
        [0.992, 0.815, 0.118],  # Al
        [0.505, 0.320, 0.110],  # Si
        [0.006, 0.816, 0.613],  # Si
        [0.491, 0.625, 0.487],  # O
        [0.024, 0.124, 0.995],  # O
        [0.073, 0.488, 0.635],  # O
        [0.576, 0.990, 0.143],  # O
        [0.298, 0.356, 0.612],  # O
        [0.817, 0.855, 0.142],  # O
        [0.517, 0.179, 0.610],  # O
        [0.000, 0.680, 0.104],  # O
    ],
    name='Anorthite',
)

ref_peaks = anorthite.get_ref_peaks(CU_KA, two_theta_range=(5, 90))
print(f"\nAnorthite reference: {len(ref_peaks)} reflections")

# ── Step 3: Guided Williamson-Hall ──
print("\n--- Guided Williamson-Hall (reciprocal space) ---")
gwh = profile.guided_williamson_hall(phase=anorthite, n_sigma=3.0,
                                      tolerance_d=0.03)

print(f"Zero-point offset: {gwh['zero_offset']:.4f} deg")
print(f"Anorthite peaks found: {gwh['n_peaks']}")
print(f"Crystallite size: {gwh['crystallite_size']:.0f} angstroms")
print(f"Microstrain: {gwh['strain']:.5f}")
print(f"R-squared: {gwh['r_squared']:.4f}")

# Peak details in d-spacing
pk = gwh['peaks']
print(f"\n  {'d_obs':>8} {'d_ref':>8} {'shift':>8} {'FWHM':>8} {'SNR':>6}")
for i in range(min(10, len(pk['d_obs']))):
    print(f"  {pk['d_obs'][i]:8.3f} {pk['d_ref'][i]:8.3f} "
          f"{pk['d_shift'][i]:8.4f} {pk['fwhm'][i]:8.4f} "
          f"{pk['snr'][i]:6.1f}")
if len(pk['d_obs']) > 10:
    print(f"  ... and {len(pk['d_obs']) - 10} more peaks")

# ── Step 4: Guided Warren-Averbach ──
print("\n--- Guided Warren-Averbach ---")
gwa = profile.guided_warren_averbach(phase=anorthite, n_sigma=3.0,
                                      tolerance_d=0.03)

print(f"Harmonic families analysed: {gwa['n_families']}")
print(f"Families rejected: {gwa['n_families_rejected']}")
print(f"Median crystallite size: {gwa['median_crystallite_size']:.0f} angstroms")
print(f"Mean RMS strain: {gwa['mean_rms_strain']:.4f}")

# Per-family details
print(f"\n  {'hkl':>14} {'orders':<12} {'D (A)':>8} {'rms_e':>8}")
for fam in gwa['families']:
    hkl = fam['base_hkl']
    hkl_str = f"({hkl[0]},{hkl[1]},{hkl[2]})"
    orders = ','.join(str(o) for o in fam['orders'])
    D = f"{fam['crystallite_size']:.0f}" if not np.isnan(fam['crystallite_size']) else 'N/A'
    e = f"{fam['rms_strain']:.4f}" if not np.isnan(fam['rms_strain']) else 'N/A'
    print(f"  {hkl_str:>14} n=[{orders}]{'':<6} {D:>8} {e:>8}")

# ── Step 5: Pair distribution function ──
print("\n--- Pair distribution function ---")
r, g_r = profile.compute_pdf()
mask = r > 0
print(f"PDF computed: {len(r[mask])} points, "
      f"r = {r[mask].min():.1f} to {r[mask].max():.1f} angstroms")

# ── Plots (optional) ──
try:
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(15, 4))

    profile.plot_pattern(ax=axes[0], x_axis='d_spacing')
    axes[0].set_xlim(1.5, 8)
    axes[0].invert_xaxis()

    profile.plot_williamson_hall_reciprocal(ax=axes[1])

    profile.plot_pdf(ax=axes[2], r_max=20)

    fig.tight_layout()
    fig.savefig('lab_lunar_example.png', dpi=150)
    print("\nFigure saved to lab_lunar_example.png")
except Exception as e:
    print(f"\nPlotting skipped: {e}")
