# `xrd_profile/registry/` — instrumental profile registry

Drop `<instrument-name>.json` files here to make pre-fit Caglioti
profiles available via `InstrumentalProfile.from_registry(name)`.

## File format (`schema_version: '1'`)

```json
{
  "schema_version": "1",
  "U": 5.0e-3,
  "V": -1.0e-3,
  "W": 5.0e-3,
  "wavelength": 1.5406,
  "name": "Lab_Bruker_Cu_Ka"
}
```

`U`, `V`, `W` are Caglioti polynomial coefficients (deg²). `wavelength`
is in angstroms. The `name` field is informational and need not match
the filename.

## Generating a profile from a measured standard

```python
from xrd_profile import InstrumentalStandard
import numpy as np

data = np.loadtxt('lab6_pattern.xy')
std = InstrumentalStandard.from_cif_and_pattern(
    cif='LaB6.cif',
    two_theta=data[:, 0], intensity=data[:, 1],
    wavelength=1.5406, name='Lab_Bruker_Cu_Ka')
profile = std.caglioti_fit()
profile.to_json('xrd_profile/registry/lab_bruker_cu_ka.json')
```

## v0.4.0 status

The registry ships **empty** in v0.4.0. Calibration-grade Caglioti
fits for the Misasa Rigaku, Winnipeg Bruker, and Diamond I11 instruments
are a separate calibration deliverable. Users supply their own
profiles, either by fitting their own LaB6/Si standard and saving the
result here, or by hand-editing a file with literature U/V/W values.
