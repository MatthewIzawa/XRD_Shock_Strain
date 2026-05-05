"""
Phase abstraction over pymatgen.Structure for guided peak detection.

The Phase class wraps a crystal structure and exposes reference-peak
generation through get_ref_peaks() and get_ref_d(). pymatgen is an
optional dependency; if missing, a clear ImportError is raised that
points to the [cif] extra install command.
"""
from __future__ import annotations

from pathlib import Path
from typing import List, Tuple, Optional, Sequence

import numpy as np

from .conversions import two_theta_to_d


def _require_pymatgen():
    """Lazy import helper. Raises ImportError with install instructions."""
    try:
        from pymatgen.core.structure import Structure, Lattice
        from pymatgen.analysis.diffraction.xrd import XRDCalculator
        return Structure, Lattice, XRDCalculator
    except ImportError as e:
        raise ImportError(
            'pymatgen is required for Phase.from_cif and '
            'Phase.from_lattice_params. '
            'Install with: pip install xrd_profile[cif]'
        ) from e


class Phase:
    """A crystalline phase with a known structure."""

    def __init__(self, structure, name: str = ''):
        self.structure = structure
        self.name = name

    @classmethod
    def from_cif(cls, path, name: Optional[str] = None) -> 'Phase':
        """Load from a CIF file. Default name = file stem."""
        Structure, _, _ = _require_pymatgen()
        structure = Structure.from_file(str(path))
        if name is None:
            name = Path(path).stem
        return cls(structure, name=name)

    @classmethod
    def from_lattice_params(cls, a: float, b: float, c: float,
                             alpha: float, beta: float, gamma: float,
                             species: Sequence,
                             coords: Sequence,
                             name: str = '') -> 'Phase':
        """Build inline from lattice parameters and atomic positions."""
        Structure, Lattice, _ = _require_pymatgen()
        lat = Lattice.from_parameters(a, b, c, alpha, beta, gamma)
        structure = Structure(lat, species, coords)
        return cls(structure, name=name)

    def get_ref_peaks(self, wavelength: float,
                      two_theta_range: Tuple[float, float] = (5, 90),
                      min_intensity: float = 3.0) -> List[dict]:
        """List of {d, two_theta, intensity, h, k, l} dicts with intensity >= min_intensity."""
        _, _, XRDCalculator = _require_pymatgen()
        calc = XRDCalculator(wavelength=wavelength)
        pattern = calc.get_pattern(self.structure,
                                    two_theta_range=two_theta_range)
        peaks: List[dict] = []
        for i in range(len(pattern.x)):
            if pattern.y[i] < min_intensity:
                continue
            hkl = pattern.hkls[i][0]['hkl']
            d = two_theta_to_d(pattern.x[i], wavelength)
            peaks.append({
                'd': float(d),
                'two_theta': float(pattern.x[i]),
                'intensity': float(pattern.y[i]),
                'h': int(hkl[0]),
                'k': int(hkl[1]),
                'l': int(hkl[2]),
            })
        return peaks

    def get_ref_d(self, wavelength: float,
                   two_theta_range: Tuple[float, float] = (5, 90),
                   min_intensity: float = 3.0,
                   sorted_by_intensity: bool = True) -> np.ndarray:
        """Reference d-spacings, optionally sorted by intensity descending."""
        peaks = self.get_ref_peaks(wavelength, two_theta_range, min_intensity)
        if sorted_by_intensity:
            peaks.sort(key=lambda p: -p['intensity'])
        return np.array([p['d'] for p in peaks], dtype=float)

    def __repr__(self) -> str:
        try:
            comp = self.structure.composition.reduced_formula
        except Exception:
            comp = '?'
        return f'<Phase {self.name!r}: {comp}>'


def build_reference_peaks(structure, wavelength: float,
                           two_theta_range: Tuple[float, float] = (5, 90),
                           min_intensity: float = 3.0) -> List[dict]:
    """Standalone reference-peak builder for users who already have a
    pymatgen Structure. Equivalent to Phase(structure).get_ref_peaks(...)."""
    return Phase(structure, name='').get_ref_peaks(
        wavelength, two_theta_range, min_intensity)
