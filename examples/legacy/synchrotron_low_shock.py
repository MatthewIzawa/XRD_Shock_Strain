"""
[v0.2.0 LEGACY] This file is a verbatim copy of the v0.2.0 example
script. It demonstrates the pre-v0.3.0 pattern of constructing a
pymatgen.Structure inline and writing a build_ref helper. The top-level
example with the same filename in examples/ has been updated to use
the v0.3.0 Phase API; this file is preserved as a historical reference.

(Original docstring follows.)
"""

"""
Example: Synchrotron XRD analysis of a weakly shocked eucrite (Tirhert).

Analyses Tirhert (unbrecciated eucrite, fresh fall, low shock) from
Diamond Light Source beamline I11 (15 keV, lambda = 0.826517 angstroms,
10-148 deg 2-theta, 0.001 deg steps).

Demonstrates multi-phase guided analysis using both anorthite and
pigeonite reference structures, plus pair distribution function
analysis at high Q-range.

Data: DLS beamtime ee17803-1, Izawa & Jephcoat (2018).
"""

import numpy as np
from pymatgen.core.structure import Structure, Lattice
from pymatgen.analysis.diffraction.xrd import XRDCalculator
from xrd_profile import XRDProfile, two_theta_to_d

# ── Instrument parameters ──
LAMBDA = 0.826517  # I11 at 15 keV

# ── Load data ──
data = np.loadtxt(
    r'C:\Users\Matthew Izawa\Desktop\111 Backup 20220530\transfer\IPM'
    r'\2018\ee17803-1\processing\Tirhert_summed_0001.xye')
tt_full, i_full, err_full = data[:, 0], data[:, 1], data[:, 2]

# Use the refined angular range
mask = (tt_full >= 10) & (tt_full <= 148)
tt, intensity = tt_full[mask], i_full[mask]

Q_max = 4 * np.pi * np.sin(np.radians(tt.max() / 2)) / LAMBDA
print(f"Tirhert (low-shock eucrite)")
print(f"  {len(tt)} points, {tt.min():.0f}-{tt.max():.0f} deg")
print(f"  Q_max = {Q_max:.1f} /A, dr = {np.pi/Q_max:.3f} A")

profile = XRDProfile(tt, intensity, wavelength=LAMBDA,
                     sample_name='Tirhert')

# ── Build anorthite reference ──
an_lat = Lattice.from_parameters(8.1809, 12.881, 7.1101,
                                  93.465, 116.11, 90.369)
an_species = ['Ca','Al','Al','Si','Si','O','O','O','O','O','O','O','O']
an_coords = [
    [0.269,0.988,0.086],[0.507,0.314,0.621],[0.992,0.815,0.118],
    [0.505,0.320,0.110],[0.006,0.816,0.613],[0.491,0.625,0.487],
    [0.024,0.124,0.995],[0.073,0.488,0.635],[0.576,0.990,0.143],
    [0.298,0.356,0.612],[0.817,0.855,0.142],[0.517,0.179,0.610],
    [0.000,0.680,0.104],
]

# ── Build pigeonite reference ──
# Pigeonite CIF from Morimoto, loaded via pymatgen
cif_dir = (r'C:\Users\Matthew Izawa\Desktop\Ye olde seagate'
           r'\Big Bad Bucket of Backups\transfer\Mar 2016'
           r'\The New Era - HoserLab\Rietveld\Structures\CIF files')

xrd_calc = XRDCalculator(wavelength=LAMBDA)

def build_ref(structure, two_theta_range=(5, 148), min_intensity=3.0):
    """Build reference peak list from a pymatgen Structure."""
    pattern = xrd_calc.get_pattern(structure, two_theta_range=two_theta_range)
    ref_peaks = []
    for i in range(len(pattern.x)):
        hkl = pattern.hkls[i][0]['hkl']
        d = two_theta_to_d(pattern.x[i], LAMBDA)
        ref_peaks.append({
            'd': float(d), 'two_theta': float(pattern.x[i]),
            'intensity': float(pattern.y[i]),
            'h': hkl[0], 'k': hkl[1], 'l': hkl[2],
        })
    ref_d = np.array(sorted(
        [p['d'] for p in ref_peaks if p['intensity'] >= min_intensity],
        reverse=True))
    return ref_peaks, ref_d

an_struct = Structure(an_lat, an_species, an_coords)
an_peaks, an_d = build_ref(an_struct)
print(f"\nAnorthite reference: {len(an_peaks)} reflections, "
      f"{len(an_d)} peaks (I >= 3%)")

pig_struct = Structure.from_file(f'{cif_dir}\\Pigeonite - Morimoto.cif')
pig_peaks, pig_d = build_ref(pig_struct)
print(f"Pigeonite reference: {len(pig_peaks)} reflections, "
      f"{len(pig_d)} peaks (I >= 3%)")

# ── Anorthite analysis ──
print("\n--- Anorthite (plagioclase) ---")
an_wh = profile.guided_williamson_hall(an_d, n_sigma=3.0, tolerance_d=0.02)
an_wa = profile.guided_warren_averbach(an_peaks, n_sigma=3.0, tolerance_d=0.02)

print(f"  W-H: {an_wh['n_peaks']} peaks, D = {an_wh['crystallite_size']:.0f} A, "
      f"strain = {an_wh['strain']:.5f}")
print(f"  W-A: {an_wa['n_families']} families, "
      f"D_median = {an_wa['median_crystallite_size']:.0f} A, "
      f"rms_strain = {an_wa['mean_rms_strain']:.4f}")

# ── Pigeonite analysis ──
print("\n--- Pigeonite (pyroxene) ---")
pig_wh = profile.guided_williamson_hall(pig_d, n_sigma=3.0, tolerance_d=0.02)
pig_wa = profile.guided_warren_averbach(pig_peaks, n_sigma=3.0, tolerance_d=0.02)

print(f"  W-H: {pig_wh['n_peaks']} peaks, D = {pig_wh['crystallite_size']:.0f} A, "
      f"strain = {pig_wh['strain']:.5f}")
print(f"  W-A: {pig_wa['n_families']} families, "
      f"D_median = {pig_wa['median_crystallite_size']:.0f} A, "
      f"rms_strain = {pig_wa['mean_rms_strain']:.4f}")

# ── Phase comparison ──
an_s = an_wh['crystallite_size']
pig_s = pig_wh['crystallite_size']
if an_s > 0 and pig_s > 0:
    ratio = an_s / pig_s
    print(f"\n--- Phase comparison ---")
    print(f"  Anorthite / Pigeonite size ratio: {ratio:.2f}")
    print(f"  (Ratio near 1.0 expected for unshocked material)")

# ── PDF ──
print("\n--- Pair distribution function ---")
r, g_r = profile.compute_pdf()
print(f"  Q_max = {Q_max:.1f} /A gives real-space resolution "
      f"dr = {np.pi/Q_max:.3f} A")
print(f"  Compare: lab Cu Ka gives dr ~ 0.77 A, "
      f"lab Co Ka gives dr ~ 0.70 A")

# ── Plots ──
try:
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # Pattern in d-spacing
    ax = axes[0, 0]
    profile.plot_pattern(ax=ax, x_axis='d_spacing', linewidth=0.3, color='k')
    ax.set_xlim(0.6, 10)
    ax.invert_xaxis()
    ax.set_title('Tirhert: diffraction pattern')

    # W-H comparison
    ax = axes[0, 1]
    ax.scatter(an_wh['K'], an_wh['deltaK'], s=15, alpha=0.6,
               color='#2166ac', label='Anorthite')
    ax.scatter(pig_wh['K'], pig_wh['deltaK'], s=15, alpha=0.6,
               color='#b2182b', label='Pigeonite')
    ax.set_xlabel(r'K ($\AA^{-1}$)')
    ax.set_ylabel(r'$\Delta$K ($\AA^{-1}$)')
    ax.set_title('Guided W-H: two-phase comparison')
    ax.legend()

    # PDF
    ax = axes[1, 0]
    mask_pdf = (r > 1) & (r < 20)
    ax.plot(r[mask_pdf], g_r[mask_pdf], 'k', linewidth=0.7)
    ax.set_xlabel(r'r ($\AA$)')
    ax.set_ylabel('G(r)')
    ax.set_title(f'PDF (Q_max = {Q_max:.1f} /A)')

    # W-A anisotropic sizes
    ax = axes[1, 1]
    for phase, wa, color, label in [
        ('an', an_wa, '#2166ac', 'Anorthite'),
        ('pig', pig_wa, '#b2182b', 'Pigeonite')
    ]:
        sizes = [(f"({f['base_hkl'][0]},{f['base_hkl'][1]},{f['base_hkl'][2]})",
                  f['crystallite_size'])
                 for f in wa['families']
                 if not np.isnan(f['crystallite_size']) and f['crystallite_size'] > 0]
        if sizes:
            labels, vals = zip(*sizes[:8])
            y = np.arange(len(vals))
            ax.barh(y, vals, 0.4, color=color, alpha=0.7, label=label)
            ax.set_yticks(y)
            ax.set_yticklabels(labels, fontsize=7)
    ax.set_xlabel('Crystallite size ($\\AA$)')
    ax.set_title('W-A anisotropic sizes')
    ax.legend(fontsize=8)
    ax.invert_yaxis()

    fig.suptitle('Tirhert (low-shock eucrite): synchrotron multi-phase analysis',
                 fontsize=13)
    fig.tight_layout()
    fig.savefig('synchrotron_low_shock_example.png', dpi=150)
    print("\nFigure saved to synchrotron_low_shock_example.png")
except Exception as e:
    print(f"\nPlotting skipped: {e}")
