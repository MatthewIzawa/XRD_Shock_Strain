"""
Example: multi-phase analysis with the v0.3.0 Phase API.

Demonstrates:
  - Loading phases from CIF files (forsterite, pigeonite).
  - Bundled diffraction pattern (Tirhert subset, ~23k points).
  - run_all() bundled call covering W-H, W-A, PDF, Scherrer.
  - Per-phase result extraction.

This script is fully self-contained: it uses the CIFs and data
shipped in examples/. To run:

    pip install xrd_profile[cif]
    python examples/multi_phase_olivine.py
"""
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

from xrd_profile import XRDProfile
from xrd_profile.phases import Phase


HERE = Path(__file__).parent
DATA = HERE / 'data' / 'tirhert_subset.xy'
CIFS = HERE / 'cifs'
LAMBDA = 0.826517  # I11 at 15 keV — matches the bundled fixture

# ── Load data ──
xy = np.loadtxt(DATA)
tt, intensity = xy[:, 0], xy[:, 1]
profile = XRDProfile(tt, intensity, wavelength=LAMBDA,
                      sample_name='Tirhert (subset)')

# ── Build phases ──
forsterite = Phase.from_cif(CIFS / 'Forsterite.cif', name='Forsterite')
pigeonite = Phase.from_cif(CIFS / 'Pigeonite.cif', name='Pigeonite')
print(f"Phases: {forsterite}, {pigeonite}")

# ── Bundled analysis ──
results = profile.run_all(
    methods=['wh', 'wa', 'pdf', 'scherrer'],
    phases=[forsterite, pigeonite],
    wh={'n_sigma': 3.0, 'tolerance_d': 0.03},
    wa={'n_sigma': 3.0, 'tolerance_d': 0.03},
    pdf={'cheby_order': 15, 'lorch': True},
)

# ── Print summary ──
for phase_name in ['Forsterite', 'Pigeonite']:
    if phase_name in results['wh']:
        wh = results['wh'][phase_name]
        wa = results['wa'][phase_name]
        print(f"\n--- {phase_name} ---")
        print(f"  W-H: D = {wh.get('crystallite_size', float('nan')):.0f} A, "
              f"strain = {wh.get('strain', float('nan')):.5f}")
        print(f"  W-A: D_median = "
              f"{wa.get('median_crystallite_size', float('nan')):.0f} A, "
              f"families = {wa.get('n_families', 0)}")

print(f"\nPDF Q_max = {results['pdf']['Q_max']:.2f} /A, "
      f"dr = {np.pi / results['pdf']['Q_max']:.3f} A")
print(f"Scherrer mean size = {results['scherrer']['mean_size']:.0f} A")

# ── 4-panel figure ──
fig, axes = plt.subplots(2, 2, figsize=(13, 9))

# (a) pattern in d-spacing
ax = axes[0, 0]
profile.plot_pattern(ax=ax, x_axis='d_spacing', linewidth=0.4, color='k')
ax.set_xlim(0.8, 8)
ax.invert_xaxis()
ax.set_title('(a) Tirhert subset')

# (b) W-H reciprocal-space
ax = axes[0, 1]
for name, color in [('Forsterite', '#2166ac'), ('Pigeonite', '#b2182b')]:
    if name in results['wh']:
        wh = results['wh'][name]
        if 'K' in wh and 'deltaK' in wh:
            ax.scatter(wh['K'], wh['deltaK'], s=16, alpha=0.7,
                       color=color, label=name)
ax.set_xlabel(r'K ($\AA^{-1}$)')
ax.set_ylabel(r'$\Delta$K ($\AA^{-1}$)')
ax.set_title('(b) Guided W-H per phase')
ax.legend()

# (c) W-A anisotropic sizes per phase
ax = axes[1, 0]
y_offset = 0
for name, color in [('Forsterite', '#2166ac'), ('Pigeonite', '#b2182b')]:
    if name not in results['wa']:
        continue
    wa = results['wa'][name]
    sizes = [(f"{f['base_hkl']}", f['crystallite_size'])
             for f in wa.get('families', [])
             if not np.isnan(f['crystallite_size']) and f['crystallite_size'] > 0]
    if not sizes:
        continue
    labels, vals = zip(*sizes[:8])
    y = np.arange(len(vals)) + y_offset
    ax.barh(y, vals, 0.4, color=color, alpha=0.7, label=name)
    y_offset += len(vals) + 1
ax.set_xlabel('Crystallite size ($\\AA$)')
ax.set_title('(c) W-A anisotropic sizes')
ax.legend(fontsize=9)
ax.invert_yaxis()

# (d) PDF
ax = axes[1, 1]
r, G = results['pdf']['r'], results['pdf']['G_r']
mask = (r > 1) & (r < 15)
ax.plot(r[mask], G[mask], 'k', linewidth=0.8)
ax.set_xlabel(r'r ($\AA$)')
ax.set_ylabel('G(r)')
ax.axhline(0, color='grey', linewidth=0.3, linestyle=':')
ax.set_title(f'(d) PDF (Q_max = {results["pdf"]["Q_max"]:.1f} /A)')

fig.suptitle('xrd_profile v0.3.0 multi-phase example',
              fontsize=13)
fig.tight_layout()
out_png = HERE / 'multi_phase_olivine_example.png'
fig.savefig(out_png, dpi=150)
print(f"\nFigure saved to {out_png}")
