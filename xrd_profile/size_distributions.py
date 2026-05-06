"""
size_distributions.py — Lognormal and normal column-length distribution
fits to Warren-Averbach A_size(L) data.

A_size(L) for a lognormal column-length distribution g(D) ~ Lognormal:

    A_size(L) = (1/<D>) * integral_{L}^{infinity} (D - L) g(D) dD

For a lognormal, this has a closed form involving the complementary
error function (Krill & Birringer 1998; Langford, Louer & Scardi 2000):

    A_size(L) = 0.5 * (1 - erf((ln(L/D_med) - sigma^2)/(sigma sqrt(2))))
              + (L / (D_med * exp(sigma^2/2)))
              * 0.5 * (1 - erf((ln(L/D_med) + 0.5*sigma^2)/(sigma sqrt(2))))

where D_med is the median crystallite size and sigma is the log-space
standard deviation. (At L=0 this reduces to A_size(0) = 1.)

For a normal column-length distribution centred at <D> with stddev s,
the closed form involves the standard error function similarly.

References
----------
Krill, C. E. & Birringer, R. (1998). Estimating grain-size distributions
    in nanocrystalline materials from X-ray diffraction profile analysis.
    Phil. Mag. A 77, 621-640.
Langford, J. I., Louer, D. & Scardi, P. (2000). Effect of a crystallite
    size distribution on X-ray diffraction line profiles and whole-powder-
    pattern fitting. J. Appl. Cryst. 33, 964-974.
"""
import numpy as np
from scipy.optimize import curve_fit
from scipy.special import erf
from scipy.stats import linregress


def lognormal_a_size(L, D_median, sigma):
    """A_size(L) for a lognormal column-length distribution.

    Parameters
    ----------
    L : np.ndarray
        Column lengths (angstroms).
    D_median : float
        Median crystallite size (angstroms).
    sigma : float
        Log-space stddev (dimensionless).

    Returns
    -------
    A : np.ndarray
        Normalised so A[0] = 1 (when L[0] = 0).
    """
    L = np.asarray(L, dtype=float)
    if D_median <= 0 or sigma <= 0:
        return np.full_like(L, np.nan)
    L_safe = np.maximum(L, 1e-10)
    arg1 = (np.log(L_safe / D_median) - sigma**2) / (sigma * np.sqrt(2))
    arg2 = (np.log(L_safe / D_median) + 0.5*sigma**2) / (sigma * np.sqrt(2))
    term1 = 0.5 * (1 - erf(arg1))
    D_mean_vol = D_median * np.exp(sigma**2 / 2)
    term2 = (L_safe / D_mean_vol) * 0.5 * (1 - erf(arg2))
    A = term1 - term2
    # Force A(L=0) = 1 (analytic limit; numerical erf may drift).
    A = np.where(L < 1e-10, 1.0, A)
    return np.clip(A, 0.0, 1.0)


def normal_a_size(L, D_mean, sigma):
    """A_size(L) for a normal column-length distribution.

    For a normal g(D) = N(D_mean, sigma), the size Fourier coefficient is

        A_size(L) = max(0, 1 - L / D_mean)

    convolved with the distribution; the closed form is approximately
    (1 - L/D_mean) * 0.5 * (1 + erf((D_mean - L)/(sigma sqrt(2)))) for
    sigma << D_mean. For broader distributions a more careful integral
    applies; we use the closed form here (Scardi & Leoni 2001).
    """
    L = np.asarray(L, dtype=float)
    if D_mean <= 0 or sigma <= 0:
        return np.full_like(L, np.nan)
    arg = (D_mean - L) / (sigma * np.sqrt(2))
    cdf_term = 0.5 * (1 + erf(arg))
    base = np.maximum(1 - L / D_mean, 0)
    A = base * cdf_term
    A = np.where(L < 1e-10, 1.0, A)
    return np.clip(A, 0.0, 1.0)


def _moments_initial_guess(L, A_size):
    """Initial guess for (D_median, sigma) from A_size(L).

    Estimates a characteristic column length from the initial slope
    A_size'(0) ~= -1/<D> via linear regression over the first few
    points (Warren 1969; Bertaut 1949), and uses this as the initial
    guess for D_median. The sigma guess is a fixed default of 0.3
    (mid-range for crystallite-size lognormals); curve_fit refines
    both. Falls back to (100.0, 0.3) when there are too few valid
    points or the initial slope is non-negative.

    Returns (D_median_guess, sigma_guess).
    """
    L = np.asarray(L, dtype=float)
    A = np.asarray(A_size, dtype=float)
    valid = A > 1e-3
    if np.sum(valid) < 4:
        return 100.0, 0.3
    # Linear regression on the first few points -> initial slope.
    n_init = min(5, np.sum(valid))
    sl, _, _, _, _ = linregress(L[valid][:n_init], A[valid][:n_init])
    D_median_guess = -1.0 / sl if sl < 0 else 100.0
    return D_median_guess, 0.3


def fit_size_distribution(L, A_size):
    """Fit lognormal and normal distributions to A_size(L) data.

    Returns the documented `'size_distribution'` dict, or None if
    `n_valid_L < 4`.
    """
    L = np.asarray(L, dtype=float)
    A = np.asarray(A_size, dtype=float)
    if L.shape != A.shape:
        raise ValueError(f'L shape {L.shape} != A_size shape {A.shape}')

    valid = np.isfinite(A) & (A > 0) & np.isfinite(L)
    n_valid = int(np.sum(valid))
    if n_valid < 4:
        return None

    L_valid = L[valid]
    A_valid = A[valid]

    D_med_guess, sigma_guess = _moments_initial_guess(L_valid, A_valid)

    # Lognormal fit.
    try:
        popt_ln, pcov_ln = curve_fit(
            lognormal_a_size, L_valid, A_valid,
            p0=[max(D_med_guess, 1.0), max(sigma_guess, 0.05)],
            bounds=([1.0, 0.01], [1e6, 5.0]))
        D_median_fit, sigma_fit = float(popt_ln[0]), float(popt_ln[1])
        A_pred = lognormal_a_size(L_valid, D_median_fit, sigma_fit)
        ss_res = np.sum((A_valid - A_pred)**2)
        ss_tot = np.sum((A_valid - np.mean(A_valid))**2)
        r2_ln = 1.0 - (ss_res / ss_tot) if ss_tot > 0 else np.nan
        D_mean_vol = D_median_fit * np.exp(sigma_fit**2 / 2)
        D_mean_area = D_median_fit * np.exp(sigma_fit**2 / 4)
        lognormal_result = {
            'D_median': D_median_fit,
            'sigma': sigma_fit,
            'D_mean_volume': D_mean_vol,
            'D_mean_area': D_mean_area,
            'fit_r2': float(r2_ln),
            'cov': pcov_ln,
        }
    except Exception:
        lognormal_result = {
            'D_median': np.nan, 'sigma': np.nan,
            'D_mean_volume': np.nan, 'D_mean_area': np.nan,
            'fit_r2': np.nan, 'cov': None,
        }

    # Normal fit.
    try:
        popt_n, pcov_n = curve_fit(
            normal_a_size, L_valid, A_valid,
            p0=[max(D_med_guess, 1.0), max(sigma_guess * D_med_guess, 1.0)],
            bounds=([1.0, 0.01], [1e6, 1e6]))
        D_mean_n, sigma_n = float(popt_n[0]), float(popt_n[1])
        A_pred_n = normal_a_size(L_valid, D_mean_n, sigma_n)
        ss_res_n = np.sum((A_valid - A_pred_n)**2)
        r2_n = 1.0 - (ss_res_n / ss_tot) if ss_tot > 0 else np.nan
        normal_result = {
            'D_mean': D_mean_n,
            'sigma': sigma_n,
            'fit_r2': float(r2_n),
            'cov': pcov_n,
        }
    except Exception:
        normal_result = {
            'D_mean': np.nan, 'sigma': np.nan,
            'fit_r2': np.nan, 'cov': None,
        }

    return {
        'lognormal': lognormal_result,
        'normal': normal_result,
        'method': 'curve_fit',
        'initial_guess': 'moments',
        'n_valid_L': n_valid,
    }
