"""
Example: Synchrotron XRD analysis of a weakly shocked eucrite (Tirhert).

Analyses Tirhert (unbrecciated eucrite, fresh fall, low shock) from
Diamond Light Source beamline I11 (15 keV, lambda = 0.826517 angstroms,
10-148 deg 2-theta, 0.001 deg steps).

Demonstrates the v0.3.0 Phase API: anorthite via from_lattice_params
(refined coordinates inline), pigeonite via from_cif. The legacy
v0.2.0 version (with inline build_ref helper) is preserved at
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
