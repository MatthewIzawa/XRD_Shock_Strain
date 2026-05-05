"""
Example: Synchrotron XRD analysis of a heavily shocked eucrite (JaH 626).

Analyses JaH 626 (shocked polymict eucrite with impact melt glass)
from Diamond Light Source beamline I11. Compares plagioclase and
pyroxene shock response, and demonstrates the PDF approach for
characterising structural disorder.

JaH 626 is one of the most heavily shocked eucrites known, containing
extensive impact melt glass. The plagioclase is largely converted to
maskelynite (amorphous). This sample represents the high-shock
endmember for comparison with Tirhert (low shock).

Demonstrates the v0.3.0 Phase API: anorthite via from_lattice_params
(refined coordinates inline), pigeonite via from_cif. The legacy
v0.2.0 version (with inline build_ref helper) is preserved at
examples/legacy/synchrotron_high_shock.py.

Data: DLS beamtime ee17803-1, Izawa & Jephcoat (2018).
Mineralogy: Cloutis, Izawa et al. (2013) Icarus 223, 850-877.
  40 wt% plagioclase, 59 wt% pigeonite, 1 wt% ilmenite.
"""

import numpy as np
from xrd_profile import XRDProfile
from xrd_profile.phases import Phase

LAMBDA = 0.826517
DATA_PATH = (r'C:\Users\Matthew Izawa\Desktop\111 Backup 20220530'
              r'\transfer\IPM\2018\ee17803-1\processing'
              r'\JaH_626_summed_0001.xye')
PIGEONITE_CIF = (r'C:\Users\Matthew Izawa\Desktop\Ye olde seagate'
                  r'\Big Bad Bucket of Backups\transfer\Mar 2016'
                  r'\The New Era - HoserLab\Rietveld\Structures'
                  r'\CIF files\Pigeonite - Morimoto.cif')

# ── Load data ──
data = np.loadtxt(DATA_PATH)
mask = (data[:, 0] >= 10) & (data[:, 0] <= 148)
tt, intensity = data[mask, 0], data[mask, 1]

Q_max = 4 * np.pi * np.sin(np.radians(tt.max() / 2)) / LAMBDA
print(f"JaH 626 (very high shock eucrite)")
print(f"  {len(tt)} points, Q_max = {Q_max:.1f} /A")

profile = XRDProfile(tt, intensity, wavelength=LAMBDA,
                     sample_name='JaH 626')

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

print(f"\nReferences: anorthite ({len(anorthite.get_ref_peaks(LAMBDA))} reflections), "
      f"pigeonite ({len(pigeonite.get_ref_peaks(LAMBDA))} reflections)")

# ── Multi-phase analysis ──
print("\n--- Anorthite (plagioclase / maskelynite) ---")
an_wh = profile.guided_williamson_hall(phase=anorthite, n_sigma=3.0,
                                        tolerance_d=0.02)
an_wa = profile.guided_warren_averbach(phase=anorthite, n_sigma=3.0,
                                        tolerance_d=0.02)

print(f"  W-H: {an_wh['n_peaks']} peaks, D = {an_wh['crystallite_size']:.0f} A, "
      f"strain = {an_wh['strain']:.5f}")
print(f"  W-A: {an_wa['n_families']} families, "
      f"D_median = {an_wa['median_crystallite_size']:.0f} A, "
      f"rms_strain = {an_wa['mean_rms_strain']:.4f}")
print(f"  Note: fewer peaks and families than Tirhert reflects")
print(f"  loss of crystalline anorthite to maskelynite (amorphous)")

print("\n--- Pigeonite (pyroxene) ---")
pig_wh = profile.guided_williamson_hall(phase=pigeonite, n_sigma=3.0,
                                          tolerance_d=0.02)
pig_wa = profile.guided_warren_averbach(phase=pigeonite, n_sigma=3.0,
                                          tolerance_d=0.02)

print(f"  W-H: {pig_wh['n_peaks']} peaks, D = {pig_wh['crystallite_size']:.0f} A, "
      f"strain = {pig_wh['strain']:.5f}")
print(f"  W-A: {pig_wa['n_families']} families, "
      f"D_median = {pig_wa['median_crystallite_size']:.0f} A, "
      f"rms_strain = {pig_wa['mean_rms_strain']:.4f}")

# ── Differential response ──
an_s = an_wh['crystallite_size']
pig_s = pig_wh['crystallite_size']
print(f"\n--- Differential shock response ---")
if an_s > 0 and pig_s > 0:
    ratio = an_s / pig_s
    print(f"  Anorthite D = {an_s:.0f} A")
    print(f"  Pigeonite D = {pig_s:.0f} A")
    print(f"  An/Pig ratio = {ratio:.2f}")
    if ratio < 1:
        print(f"  Plagioclase has SMALLER crystallites than pyroxene,")
        print(f"  consistent with onset of maskelynite formation.")
    else:
        print(f"  Plagioclase still has larger crystallites than pyroxene.")

print(f"\n  Compare with Tirhert (low shock): An/Pig ~ 0.96 (equal)")
print(f"  The reversal of the An/Pig ratio marks the transition from")
print(f"  pyroxene-dominated strain accommodation to plagioclase")
print(f"  amorphisation (maskelynite formation).")

# ── PDF ──
print(f"\n--- Pair distribution function ---")
r, g_r = profile.compute_pdf()
print(f"  PDF at Q_max = {Q_max:.1f} /A (dr = {np.pi/Q_max:.3f} A)")
print(f"  The PDF includes contributions from both crystalline pyroxene")
print(f"  and amorphous maskelynite. The capillary glass also")
print(f"  contributes an amorphous signal that cannot be separated")
print(f"  without a blank capillary measurement.")

# ── Plots ──
try:
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # Pattern
    ax = axes[0, 0]
    profile.plot_pattern(ax=ax, x_axis='d_spacing', linewidth=0.3, color='k')
    ax.set_xlim(0.6, 10)
    ax.invert_xaxis()
    ax.set_title('JaH 626: diffraction pattern')

    # Two-phase W-H
    ax = axes[0, 1]
    if len(an_wh['K']) > 0:
        ax.scatter(an_wh['K'], an_wh['deltaK'], s=15, alpha=0.6,
                   color='#2166ac', label=f"Anorthite ({an_wh['n_peaks']} pk)")
    if len(pig_wh['K']) > 0:
        ax.scatter(pig_wh['K'], pig_wh['deltaK'], s=15, alpha=0.6,
                   color='#b2182b', label=f"Pigeonite ({pig_wh['n_peaks']} pk)")
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

    # W-A sizes by phase
    ax = axes[1, 1]
    for wa, color, label in [
        (an_wa, '#2166ac', 'Anorthite'),
        (pig_wa, '#b2182b', 'Pigeonite')
    ]:
        sizes = [(f"({f['base_hkl'][0]},{f['base_hkl'][1]},{f['base_hkl'][2]})",
                  f['crystallite_size'])
                 for f in wa['families']
                 if not np.isnan(f['crystallite_size']) and f['crystallite_size'] > 0]
        if sizes:
            labels, vals = zip(*sizes[:6])
            y = np.arange(len(vals))
            ax.barh(y, vals, 0.4, color=color, alpha=0.7, label=label)
            ax.set_yticks(y)
            ax.set_yticklabels(labels, fontsize=7)
    ax.set_xlabel('Crystallite size ($\\AA$)')
    ax.set_title('W-A anisotropic sizes')
    ax.legend(fontsize=8)
    ax.invert_yaxis()

    fig.suptitle('JaH 626 (very high shock eucrite): synchrotron multi-phase analysis',
                 fontsize=13)
    fig.tight_layout()
    fig.savefig('synchrotron_high_shock_example.png', dpi=150)
    print("\nFigure saved to synchrotron_high_shock_example.png")
except Exception as e:
    print(f"\nPlotting skipped: {e}")
