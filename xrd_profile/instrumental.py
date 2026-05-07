"""
instrumental.py — Instrumental broadening characterisation and
deconvolution.

Two classes:
    InstrumentalStandard: holds a structural Phase plus a measured
        diffraction pattern of a known standard (LaB6, Si). Supports
        both Caglioti FWHM correction (for W-H, Scherrer) and Stokes
        Fourier deconvolution (for W-A).
    InstrumentalProfile: lightweight Caglioti carrier (U, V, W). No
        measured pattern. Supports W-H, Scherrer; W-A raises a clear
        ValueError.

Caglioti polynomial:
    FWHM(2theta)^2 = U * tan^2(theta) + V * tan(theta) + W

References
----------
Caglioti, G., Paoletti, A., Ricci, F. P. (1958). Choice of collimators
    for a crystal spectrometer for neutron diffraction. Nuclear
    Instruments 3, 223-228.
"""
import json
import warnings
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.optimize import curve_fit


class InstrumentalProfile:
    """Caglioti-polynomial carrier for instrumental broadening.

    Parameters
    ----------
    U, V, W : float
        Caglioti coefficients (deg^2).
    wavelength : float
        X-ray wavelength (angstroms).
    name : str
        Optional human-readable label.
    """

    def __init__(self, U: float, V: float, W: float,
                 wavelength: float, name: str = ''):
        self.U = float(U)
        self.V = float(V)
        self.W = float(W)
        self.wavelength = float(wavelength)
        self.name = str(name)

    def fwhm_at(self, two_theta_deg) -> float:
        """Caglioti FWHM (degrees) at the given 2-theta (degrees).

        Parameters
        ----------
        two_theta_deg : float or array-like
            2-theta value(s) in degrees. Scalar or numpy-array input is
            accepted; the return shape matches the input.

        Returns
        -------
        float or np.ndarray
            FWHM in degrees, clamped to >= 0. If the Caglioti polynomial
            evaluates negative at any input (an indication that the
            coefficients are outside their valid angular range), a
            UserWarning is emitted before clamping.
        """
        theta = np.deg2rad(np.asarray(two_theta_deg) / 2.0)
        fwhm_sq = (self.U * np.tan(theta)**2
                   + self.V * np.tan(theta)
                   + self.W)
        if np.any(fwhm_sq < 0.0):
            warnings.warn(
                'Caglioti polynomial negative at one or more 2theta '
                'values; FWHM clamped to 0. Check coefficient validity '
                'for this angular range.',
                UserWarning, stacklevel=2)
        return np.sqrt(np.maximum(fwhm_sq, 0.0))

    SCHEMA_VERSION = '1'

    def to_json(self, path) -> None:
        """Serialise this profile to a JSON file at `path`."""
        path = Path(path)
        path.write_text(json.dumps({
            'schema_version': self.SCHEMA_VERSION,
            'U': self.U,
            'V': self.V,
            'W': self.W,
            'wavelength': self.wavelength,
            'name': self.name,
        }, indent=2))

    @classmethod
    def from_json(cls, path) -> 'InstrumentalProfile':
        """Load a profile from a JSON file produced by `to_json`."""
        data = json.loads(Path(path).read_text())
        if data.get('schema_version') != cls.SCHEMA_VERSION:
            raise ValueError(
                f'Unsupported InstrumentalProfile schema_version '
                f'{data.get("schema_version")!r}; this code expects '
                f'{cls.SCHEMA_VERSION!r}')
        return cls(U=data['U'], V=data['V'], W=data['W'],
                   wavelength=data['wavelength'],
                   name=data.get('name', ''))

    @classmethod
    def from_standard(cls, std) -> 'InstrumentalProfile':
        """Convenience: fit Caglioti to the standard, return the
        resulting InstrumentalProfile. Equivalent to
        `std.caglioti_fit()`."""
        return std.caglioti_fit()

    @classmethod
    def from_registry(cls, name: str) -> 'InstrumentalProfile':
        """Look up a pre-fit profile by name in
        `xrd_profile/registry/<name>.json`. v0.4.0 ships an empty
        registry; users populate it with their own JSON profiles."""
        registry_dir = Path(__file__).parent / 'registry'
        candidate = registry_dir / f'{name}.json'
        if not candidate.is_file():
            raise KeyError(
                f'No registered InstrumentalProfile {name!r}; '
                f'expected file at {candidate}. The v0.4.0 registry '
                f'ships empty; populate via `InstrumentalProfile.to_json` '
                f'into the registry directory.')
        return cls.from_json(candidate)


def _stokes_deconvolve(A_obs, A_inst, damping_threshold: float = 0.05):
    """Stokes Fourier deconvolution of a peak profile.

    A_corr(L) = A_obs(L) / A_inst(L), with a damping floor: when
    A_inst(L) < damping_threshold * A_inst(0), the coefficient is set
    to 0 to suppress noise amplification at high L (Stokes 1948).

    Parameters
    ----------
    A_obs : np.ndarray
        Observed sample profile Fourier coefficients.
    A_inst : np.ndarray
        Instrumental profile Fourier coefficients (same length).
    damping_threshold : float
        Fraction of A_inst(0) below which deconvolution is suppressed.

    Returns
    -------
    A_corr : np.ndarray
        Deconvolved sample-only Fourier coefficients.
    """
    A_obs = np.asarray(A_obs, dtype=float)
    A_inst = np.asarray(A_inst, dtype=float)
    if A_obs.shape != A_inst.shape:
        raise ValueError(f'A_obs shape {A_obs.shape} != A_inst shape '
                         f'{A_inst.shape}')
    if A_inst.size == 0:
        return np.array([])
    if A_inst[0] == 0:
        raise ValueError('A_inst(0) is zero; cannot Stokes-deconvolve.')
    # 5% floor is the conventional Stokes (1948) damping heuristic;
    # see Warren & Averbach (1952) J. Appl. Phys. 23, 497 for usage.
    threshold = damping_threshold * abs(A_inst[0])
    A_corr = np.zeros_like(A_obs)
    keep = np.abs(A_inst) >= threshold
    A_corr[keep] = A_obs[keep] / A_inst[keep]
    return A_corr


def _caglioti_model(two_theta_deg, U, V, W):
    """Caglioti FWHM(2theta) given U, V, W. Vectorised."""
    theta = np.deg2rad(np.asarray(two_theta_deg) / 2.0)
    fwhm_sq = U * np.tan(theta)**2 + V * np.tan(theta) + W
    # Internal-only floor at 1e-10 (vs 0.0 in fwhm_at): curve_fit needs
    # a strictly positive sqrt argument to avoid a zero-derivative wall
    # during optimisation. The public fwhm_at uses 0.0 + a UserWarning.
    return np.sqrt(np.maximum(fwhm_sq, 1e-10))


def _measure_fwhm_near(two_theta, intensity, target_tt,
                        search_window_deg=0.5, n_dense=5000):
    """Half-max-interpolation FWHM of the peak nearest target_tt
    within +/- search_window_deg. Returns (fwhm_deg, observed_tt) or
    (None, None) if no resolvable peak found.

    A cubic-spline is fitted to the local window before finding the
    half-max crossings. This is equivalent to half-max interpolation
    but with sub-pixel precision — necessary when the scan step is
    only a few times smaller than the FWHM.
    """
    mask = np.abs(two_theta - target_tt) < search_window_deg
    if not np.any(mask):
        return None, None
    local_tt = two_theta[mask]
    local_i = intensity[mask]
    if len(local_tt) < 4:
        return None, None
    # Subtract a flat local baseline (5th percentile inside the window).
    baseline = np.percentile(local_i, 5)
    local_i_corr = local_i - baseline
    if local_i_corr.max() <= 0:
        return None, None
    # Fit a cubic spline for sub-pixel half-max crossing resolution.
    try:
        cs = CubicSpline(local_tt, local_i_corr)
    except ValueError:
        return None, None
    dense_tt = np.linspace(local_tt[0], local_tt[-1], n_dense)
    dense_i = np.maximum(cs(dense_tt), 0.0)
    peak_idx = int(np.argmax(dense_i))
    peak_val = dense_i[peak_idx]
    if peak_val <= 0:
        return None, None
    half_max = peak_val / 2.0
    obs_tt = float(dense_tt[peak_idx])
    # Find half-max crossings on each side of the peak.
    left = dense_i[:peak_idx + 1]
    right = dense_i[peak_idx:]
    left_below = np.where(left < half_max)[0]
    right_below = np.where(right < half_max)[0]
    if len(left_below) == 0 or len(right_below) == 0:
        return None, None
    li = left_below[-1]
    ri = right_below[0] + peak_idx
    if li + 1 >= len(dense_tt) or dense_i[li + 1] == dense_i[li]:
        tt_left = dense_tt[li]
    else:
        frac = ((half_max - dense_i[li])
                / (dense_i[li + 1] - dense_i[li]))
        tt_left = dense_tt[li] + frac * (dense_tt[li + 1] - dense_tt[li])
    if ri == 0 or dense_i[ri - 1] == dense_i[ri]:
        tt_right = dense_tt[ri]
    else:
        frac = ((half_max - dense_i[ri])
                / (dense_i[ri - 1] - dense_i[ri]))
        tt_right = dense_tt[ri] - frac * (dense_tt[ri] - dense_tt[ri - 1])
    return float(tt_right - tt_left), obs_tt


def _caglioti_fit(two_theta, intensity, ref_two_theta,
                  search_window_deg=0.5):
    """Fit Caglioti U, V, W to the FWHMs of `intensity` at the peaks
    nearest each entry in `ref_two_theta`.

    Parameters
    ----------
    two_theta : np.ndarray
        Standard's 2-theta scan (degrees).
    intensity : np.ndarray
        Standard's intensity scan.
    ref_two_theta : np.ndarray
        Bragg-position 2-thetas of the standard's reflections (degrees).
    search_window_deg : float
        Local FWHM-search window around each reference position.

    Returns
    -------
    U, V, W : float
        Fitted Caglioti coefficients.
    info : dict
        Diagnostic info: 'n_peaks', 'measured_fwhms',
        'measured_positions', 'cov' (3x3 covariance matrix).
    """
    two_theta = np.asarray(two_theta, dtype=float)
    intensity = np.asarray(intensity, dtype=float)
    ref_two_theta = np.asarray(ref_two_theta, dtype=float)

    measured_fwhms = []
    measured_positions = []
    for target_tt in ref_two_theta:
        fwhm, obs_tt = _measure_fwhm_near(
            two_theta, intensity, target_tt, search_window_deg)
        if fwhm is None or fwhm <= 0:
            continue
        measured_fwhms.append(fwhm)
        measured_positions.append(obs_tt)

    if len(measured_fwhms) < 4:
        raise ValueError(
            f'Caglioti fit needs at least 4 resolvable peaks for a '
            f'meaningful 3-parameter fit; got {len(measured_fwhms)}')

    measured_fwhms = np.asarray(measured_fwhms)
    measured_positions = np.asarray(measured_positions)

    # Initial guess: U=1e-3 is appropriate for typical lab Bruker
    # broadening; for synchrotron data with much narrower peaks,
    # consider a data-driven starting value (e.g. var(fwhms)).
    # V=0 is a sensible neutral; W is initialised from the median FWHM^2.
    p0 = [1e-3, 0.0, np.median(measured_fwhms)**2]
    popt, pcov = curve_fit(
        _caglioti_model, measured_positions, measured_fwhms,
        p0=p0,
        bounds=([0.0, -np.inf, 0.0], [np.inf, np.inf, np.inf]),
    )
    U, V, W = float(popt[0]), float(popt[1]), float(popt[2])
    info = {
        'n_peaks': len(measured_fwhms),
        'measured_fwhms': measured_fwhms,
        'measured_positions': measured_positions,
        'cov': pcov,
    }
    return U, V, W, info


class InstrumentalStandard:
    """Structural Phase plus a measured diffraction pattern of a known
    standard, sufficient for both Caglioti FWHM correction and Stokes
    Fourier deconvolution.

    This object is treated as immutable; mutating `two_theta`,
    `intensity`, or `phase` after construction yields undefined behaviour
    because the cached Caglioti fit and Fourier coefficients will not be
    invalidated.

    Parameters
    ----------
    phase : xrd_profile.Phase
        The standard's structure (e.g., LaB6, Si). Provides reference
        Bragg positions for the Caglioti fit and for matching peaks
        when `fourier_coefficients(peak_d, ...)` is called.
    two_theta : np.ndarray
        Measured 2-theta scan of the standard (degrees).
    intensity : np.ndarray
        Measured intensity scan.
    wavelength : float
        X-ray wavelength (angstroms). Should match the sample being
        analysed.
    name : str
        Optional human-readable label.
    """

    def __init__(self, phase, two_theta, intensity,
                 wavelength: float, name: str = ''):
        self.phase = phase
        self.two_theta = np.asarray(two_theta, dtype=float)
        self.intensity = np.asarray(intensity, dtype=float)
        self.wavelength = float(wavelength)
        self.name = str(name)
        self._caglioti_cache = None
        self._fourier_cache = {}

    @classmethod
    def from_cif_and_pattern(cls, cif: str,
                             two_theta, intensity,
                             wavelength: float,
                             name: str = '') -> 'InstrumentalStandard':
        """Convenience constructor: load Phase from CIF, attach the
        measured pattern."""
        from .phases import Phase
        resolved_name = name or Path(cif).stem
        phase = Phase.from_cif(cif, name=resolved_name)
        return cls(phase=phase, two_theta=two_theta, intensity=intensity,
                   wavelength=wavelength, name=resolved_name)

    def caglioti_fit(self) -> 'InstrumentalProfile':
        """Fit Caglioti U, V, W to the standard's measured FWHMs at
        each reference peak. Cached. Returns an InstrumentalProfile."""
        if self._caglioti_cache is not None:
            return self._caglioti_cache
        ref_peaks = self.phase.get_ref_peaks(
            self.wavelength,
            two_theta_range=(float(self.two_theta.min()),
                             float(self.two_theta.max())))
        ref_tt = np.array([p['two_theta'] for p in ref_peaks])
        U, V, W, _info = _caglioti_fit(self.two_theta, self.intensity,
                                        ref_tt)
        prof = InstrumentalProfile(U=U, V=V, W=W,
                                    wavelength=self.wavelength,
                                    name=self.name)
        self._caglioti_cache = prof
        return prof

    def fourier_coefficients(self, peak_d: float,
                              n_coeffs: int = 20,
                              width_fwhm: float = 6.0):
        """Return (L, A_inst_L) Fourier coefficients of the standard's
        peak shape at the 2-theta corresponding to `peak_d`. `L` is
        column length in angstroms; `A_inst_L` is the (real) Fourier
        cosine coefficient, normalised so A(L=0) = 1. Cached per
        (peak_d, n_coeffs, width_fwhm).

        Two paths:

        1. **Extraction** (preferred): pull a peak window from the
           standard's measured pattern at the target 2-theta and take
           its Fourier coefficients. Captures real-instrument peak-shape
           detail when the standard has measurable signal at target_tt.

        2. **Caglioti synthesis** (fallback): if extraction returns a
           degenerate profile (zero area; A(0) == 0), synthesise a
           Gaussian peak at target_tt with FWHM given by the standard's
           Caglioti fit, then take Fourier coefficients of the synthetic
           peak. Required when target_tt falls on a near-zero baseline
           of a sparse standard pattern — the typical mode for synthetic
           LaB6 patterns generated from literature U/V/W (where peaks
           are discrete Gaussians and most of the 2-theta range is
           empty).

        New in v0.4.1: the fallback path. Prior v0.4.0 always took the
        extraction path and raised `A_inst(0) is zero` when the standard
        was sparse at target_tt.

        Parameters
        ----------
        peak_d : float
            d-spacing (angstroms) of the target reflection.
        n_coeffs : int, default 20
            Number of Fourier coefficients to return.
        width_fwhm : float, default 6.0
            Window width in units of FWHM. Used by the extraction path
            (passed to warren_averbach.extract_peak_profile) and the
            synthesis path (s-grid spans width_fwhm * fwhm_s around s0).
        """
        cache_key = (round(float(peak_d), 6), int(n_coeffs),
                     float(width_fwhm))
        if cache_key in self._fourier_cache:
            return self._fourier_cache[cache_key]

        from .conversions import d_to_two_theta
        target_tt = d_to_two_theta(peak_d, self.wavelength)
        if np.isnan(target_tt):
            raise ValueError(
                f'peak_d {peak_d} corresponds to no valid 2-theta at '
                f'wavelength {self.wavelength}')

        L, A_L, converged = self._fourier_via_extraction(
            target_tt, n_coeffs, width_fwhm)
        # Fall back to synthesis when extraction yields a degenerate
        # result. Two failure modes are seen in practice:
        # (a) zero-area profile (target_tt fell on a region with no
        #     measurable signal): A_L[0] == 0 from the W-A
        #     fourier_coefficients zero-area early return.
        # (b) noise-driven profile (median-filter background subtraction
        #     leaves only baseline noise): A_L[0] = 1 after
        #     normalisation but A_L[1:] oscillates because the
        #     "peak" was actually noise. fourier_coefficients reports
        #     this via converged=False.
        if A_L.size == 0 or A_L[0] == 0 or not converged:
            L, A_L = self._fourier_via_caglioti(
                target_tt, n_coeffs, width_fwhm)
        self._fourier_cache[cache_key] = (L, A_L)
        return L, A_L

    def _fourier_via_extraction(self, target_tt: float,
                                 n_coeffs: int, width_fwhm: float):
        """Path 1: extract a peak window from the standard's measured
        pattern at target_tt, take its Fourier coefficients. Returns
        (L, A_L, converged); A_L[0] == 0 or converged == False signals
        a degenerate extraction (caller falls back to synthesis)."""
        from .warren_averbach import (extract_peak_profile,
                                       fourier_coefficients)
        prof = extract_peak_profile(
            self.two_theta, self.intensity, target_tt, self.wavelength,
            width_fwhm=width_fwhm)
        L, A_L, converged = fourier_coefficients(
            prof['s'], prof['profile'], prof['s0'], n_coeffs=n_coeffs)
        return L, A_L, converged

    def _fourier_via_caglioti(self, target_tt: float,
                                n_coeffs: int, width_fwhm: float):
        """Path 2 (v0.4.1 fallback): synthesise a Gaussian at target_tt
        with FWHM from the standard's Caglioti fit, then take Fourier
        coefficients of the synthetic peak.

        Conversion of FWHM from 2-theta (degrees) to s = 2 sin(theta)/lambda
        (inverse angstroms): around the peak centre,
            ds/d(2-theta)_radians = cos(theta) / wavelength,
        so
            FWHM_s = FWHM_2theta * cos(theta) * (pi/180) / wavelength.
        """
        from .warren_averbach import fourier_coefficients
        cag = self.caglioti_fit()
        fwhm_2theta = float(cag.fwhm_at(target_tt))  # degrees
        theta_rad = np.radians(target_tt / 2.0)
        s0 = (2.0 * np.sin(theta_rad)) / self.wavelength
        fwhm_s = (fwhm_2theta * np.cos(theta_rad)
                  * np.pi / (180.0 * self.wavelength))
        if fwhm_s <= 0:
            # Caglioti fit returned non-positive FWHM at this 2-theta —
            # cannot synthesise a meaningful Gaussian.
            return np.zeros(n_coeffs), np.zeros(n_coeffs)
        sigma_s = fwhm_s / (2.0 * np.sqrt(2.0 * np.log(2.0)))
        half_window = width_fwhm * fwhm_s / 2.0
        n_pts = max(401, n_coeffs * 8)
        s = np.linspace(s0 - half_window, s0 + half_window, n_pts)
        profile = np.exp(-((s - s0) ** 2) / (2.0 * sigma_s ** 2))
        area = np.trapezoid(profile, s)
        if area > 0:
            profile = profile / area
        L, A_L, _ = fourier_coefficients(
            s, profile, s0, n_coeffs=n_coeffs)
        return L, A_L
