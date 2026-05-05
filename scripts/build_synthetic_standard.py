"""
Generate a synthetic LaB6 diffraction pattern with documented Caglioti
broadening, for use as a test fixture in tests/test_instrumental.py
and tests/test_size_distributions.py.

Caglioti polynomial:
    FWHM(2theta)^2 = U * tan^2(theta) + V * tan(theta) + W

with documented coefficients
    U = 5e-3 deg^2,  V = -1e-3 deg^2,  W = 5e-3 deg^2
which produces FWHMs typical of a well-aligned lab Bruker over 20-100
degrees 2-theta.

LaB6 is cubic Pm-3m, a = 4.156825 angstroms (NIST SRM 660c). Reflections
up to (3,1,1) populate the synthetic pattern. Each peak is a pure
Gaussian centred at the Bragg position with FWHM given by Caglioti.

Run: python scripts/build_synthetic_standard.py
"""
from pathlib import Path
import numpy as np

LAMBDA_CU = 1.5406  # Cu K-alpha
A_LAB6 = 4.156825   # angstroms

# Documented Caglioti coefficients used to synthesise the fixture.
SYNTH_U = 5.0e-3
SYNTH_V = -1.0e-3
SYNTH_W = 5.0e-3

# 2-theta sampling (matches a typical lab Bruker scan).
TT_MIN, TT_MAX, TT_STEP = 20.0, 100.0, 0.02

OUT = (Path(__file__).parent.parent / 'tests' / 'fixtures'
       / 'synthetic_lab6.xy')


def caglioti_fwhm(two_theta_deg, U, V, W):
    """Caglioti FWHM in degrees from polynomial coefficients."""
    theta = np.deg2rad(two_theta_deg / 2.0)
    fwhm_sq = U * np.tan(theta)**2 + V * np.tan(theta) + W
    return np.sqrt(np.maximum(fwhm_sq, 1e-8))


def lab6_reflections(max_hkl=4):
    """Allowed reflections for LaB6 (Pm-3m, all h,k,l permitted).
    Returns list of (h, k, l, multiplicity, structure_factor_proxy)."""
    refs = []
    seen = set()
    for h in range(0, max_hkl + 1):
        for k in range(0, h + 1):
            for l in range(0, k + 1):
                if h == k == l == 0:
                    continue
                key = tuple(sorted([h, k, l], reverse=True))
                if key in seen:
                    continue
                seen.add(key)
                # Structure-factor proxy: drops with hkl magnitude.
                mult = _multiplicity(h, k, l)
                f_proxy = 1.0 / (1.0 + 0.05 * (h*h + k*k + l*l))
                refs.append((h, k, l, mult, f_proxy))
    return refs


def _multiplicity(h, k, l):
    """Cubic point-group multiplicity for (h,k,l)."""
    distinct = len({abs(h), abs(k), abs(l)})
    nonzero = sum(1 for v in (h, k, l) if v != 0)
    if distinct == 1 and nonzero == 3:
        return 8     # (hhh)
    if nonzero == 1:
        return 6     # (h00)
    if distinct == 1 and nonzero == 2:
        return 12    # (hh0)
    if nonzero == 2:
        return 24    # (hk0)
    if distinct == 2 and nonzero == 3:
        return 24    # (hhl)
    return 48        # (hkl)


def bragg_two_theta(h, k, l, a, wavelength):
    d = a / np.sqrt(h*h + k*k + l*l)
    sin_theta = wavelength / (2.0 * d)
    if sin_theta >= 1.0:
        return None
    return 2.0 * np.rad2deg(np.arcsin(sin_theta))


def gaussian_peak(x, x0, fwhm, amplitude):
    sigma = fwhm / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    return amplitude * np.exp(-0.5 * ((x - x0) / sigma)**2)


def main():
    tt = np.arange(TT_MIN, TT_MAX + TT_STEP, TT_STEP)
    intensity = np.zeros_like(tt)
    intensity += 5.0  # flat background

    for h, k, l, mult, f_proxy in lab6_reflections(max_hkl=4):
        tt_peak = bragg_two_theta(h, k, l, A_LAB6, LAMBDA_CU)
        if tt_peak is None or tt_peak < TT_MIN or tt_peak > TT_MAX:
            continue
        fwhm = caglioti_fwhm(tt_peak, SYNTH_U, SYNTH_V, SYNTH_W)
        amplitude = mult * f_proxy * 1000.0
        intensity += gaussian_peak(tt, tt_peak, fwhm, amplitude)

    # Add a tiny Gaussian noise floor for realism (sigma=0.5 counts).
    rng = np.random.default_rng(seed=42)
    intensity += rng.normal(0.0, 0.5, size=intensity.shape)
    intensity = np.maximum(intensity, 0.0)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(OUT, np.column_stack([tt, intensity]),
               fmt='%.4f', header='two_theta_deg  intensity\n'
               f'# synthesis: U={SYNTH_U}, V={SYNTH_V}, W={SYNTH_W}',
               comments='# ')
    print(f'Wrote {OUT}: {len(tt)} points, '
          f'{TT_MIN}-{TT_MAX} deg, step {TT_STEP}')


if __name__ == '__main__':
    main()
