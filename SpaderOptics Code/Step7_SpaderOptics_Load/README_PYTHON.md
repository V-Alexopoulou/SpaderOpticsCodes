# SpaderOptics database: access with Python

SpaderOptics contains optical, solar-interaction, thermal-emission, equilibrium-temperature and
slab-temperature data for metallic and non-metallic materials, organized by material and by
literature source.

For access with MATLAB, see `README_MATLAB.md`.

The data are stored in the file `SpaderOptics.mat` (MATLAB `-v7.3` format). A `-v7.3` file is an
**HDF5 file**, so it can be read with open-source tools. This README explains how to read it with
Python (`h5py`, `numpy`) using the script `SpaderOptics_load.py` that is included in the deposit, and how
to find every type of data.

| File | Content |
|---|---|
| `SpaderOptics.mat` | the database (about 1.4 GB on disk, about 2.8 GB in memory when fully loaded) |
| `SpaderOptics_load.py` | Python script that loads the whole database |
| `README_PYTHON.md` | this file |

---

## 1. Requirements

- Python 3
- `numpy` and `h5py`

Installation (Linux, Debian/Ubuntu):

```bash
sudo apt install python3-numpy python3-h5py
```

or with pip: `python3 -m pip install numpy h5py`.

No MATLAB license is needed.

## 2. Loading the database

Place `SpaderOptics.mat` in the current folder (or pass its path as an argument), then run

```bash
python3 -i SpaderOptics_load.py [path/to/SpaderOptics.mat]
```

The option `-i` keeps Python open afterwards. The script creates two variables:

- `db[material][author][field]`: the data of every literature source
- `spread[material][field]`: the between-source statistics of every material (see section 5)

Loading takes a few minutes and about 3 GB of memory. If the memory is not sufficient, load a single
material instead (section 6). Copying the file from an external drive to the internal disk first
makes the loading faster.

```python
print(list(db.keys()))              # all materials
print(list(db['Al'].keys()))        # all sources (authors) of aluminum
src = db['Al']['Rakic']             # one source
print(list(src.keys()))             # all fields of this source
```

## 3. Structure of the database

```
SpaderOptics
├─ Material 1
│  ├─ dataset            -> one record per literature source (author)
│  └─ source_spread      -> statistics of the emissivity across the sources
├─ Material 2
└─ ...
```

Sources whose name ends with `-Tdepend` (for example `Ferrera-Tdepend`) combine the datasets that
the authors measured at several temperatures into one temperature-dependent entry. These entries
have an additional temperature dimension in some arrays (see section 4).

## 4. Data of one source

Notation: `λ` wavelength, `θ` incident/emission angle, `d` thickness, `T` temperature,
`C` solar concentration multiplier. Arrays have the same shape as in MATLAB. Dimensions are listed
in order.

### 4.1 Bibliographic information

| Field | Content |
|---|---|
| `author` | name of the source |
| `reference` | full reference of the optical data |
| `source_url` | DOI or URL |
| `upper_temperature_limit_K` | highest temperature of the data (melting or decomposition point) |

```python
print(src['author'], src['reference'], src['source_url'])
```

### 4.2 Original data from the literature

| Field | Content | Dimensions |
|---|---|---|
| `original_wavelength_um` | wavelength grid of the source [µm] | λ_original |
| `original_n_complex` | complex refractive index n + ik | λ_original |
| `original_temperature_K` | temperature of the source (number, or text if not available) | |
| `original_thickness_um` | thickness of the source [µm] (empty text if not available) | |

```python
wl = src['original_wavelength_um'].flatten()
n  = src['original_n_complex'].flatten()
n_real, k = n.real, n.imag
```

### 4.3 Resampled optical constants

| Field | Content | Dimensions |
|---|---|---|
| `resampled_wavelength_um` | solar wavelength grid used in the calculations [µm] | λ |
| `resampled_n_complex` | refractive index on this grid | λ |

### 4.4 Calculation grids

| Field | Content |
|---|---|
| `calculation_incident_angles_deg` | angles [deg] |
| `calculation_thickness_um` | thicknesses [µm] |
| `calculation_temperature_K` | temperatures [K] |
| `calculation_solar_concentration_multipliers` | solar concentration multipliers |

### 4.5 Optical properties

| Field | Dimensions (not T-dependent) | Dimensions (`-Tdepend`) |
|---|---|---|
| `calculated_R_directional`, `calculated_T_directional`, `calculated_A_directional` | λ × θ × d | λ × θ × d × T |
| `calculated_R_hemispherical`, `calculated_T_hemispherical`, `calculated_A_hemispherical` | λ × d | λ × d × T |
| `calculated_E_directional` | T × θ × d | T × θ × d |
| `calculated_E_hemispherical` | T × d | T × d |

(R reflectivity, T transmissivity, A absorptivity, E emissivity.)

### 4.6 Fluxes (units W m⁻²)

| Field | Dimensions (not T-dependent) | Dimensions (`-Tdepend`) |
|---|---|---|
| `calculated_directional_Reflected_flux_W_m2`, `..._Transmitted_...`, `..._Absorbed_...` | θ × d | T × θ × d |
| `calculated_hemispherical_Reflected_flux_W_m2`, `..._Transmitted_...`, `..._Absorbed_...` | d | T × d |
| `calculated_hemispherical_Emitted_flux_W_m2` | T × d | T × d |

### 4.7 Equilibrium and slab temperatures, conductive flux

| Field | Dimensions |
|---|---|
| `calculated_lumped_equilibrium_temperature_dirA_hemiE_K` | θ × d × C |
| `calculated_lumped_equilibrium_temperature_hemiA_hemiE_K` | 1 × d × C |
| `calculated_Conductive_flux_dirA_hemiE_W_m2`, `calculated_Conductive_flux_hemiA_hemiE_W_m2` | θ × d × C, and 1 × d × C |
| `calculated_slab_equilibrium_temperature_{front,rear,mid}_dirA_hemiE_K` | θ × d × C |
| `calculated_slab_equilibrium_temperature_{front,rear,mid}_hemiA_hemiE_K` | 1 × d × C |

`dirA_hemiE` uses the directional absorptivity and the hemispherical emissivity; `hemiA_hemiE`
uses the hemispherical absorptivity and the hemispherical emissivity.

### 4.8 Thermal conductivity

| Field | Content |
|---|---|
| `thermal_conductivity_temperature_grid_K` | temperatures [K] |
| `thermal_conductivity_W_mK` | thermal conductivity [W m⁻¹ K⁻¹] |
| `thermal_conductivity_reference` | reference of these data |

### 4.9 Coverage metrics

| Field | Content |
|---|---|
| `metrics_solar_coverage_ratio` | fraction of the AM0 solar spectrum covered by the wavelength range of the source (a number) |
| `metrics_thermal_coverage_ratio` | thermal coverage ratio for every temperature (vector over T) |

### 4.10 Examples

```python
import numpy as np

src = db['Al']['Rakic']
T = src['calculation_temperature_K'].flatten()
d = src['calculation_thickness_um'].flatten()
a = src['calculation_incident_angles_deg'].flatten()
wl = src['resampled_wavelength_um'].flatten()

# Hemispherical emissivity at 500 K and 10000 um (nearest grid points)
iT = np.abs(T - 500).argmin()
id_ = np.abs(d - 10000).argmin()
E = src['calculated_E_hemispherical']             # T x d
print(T[iT], d[id_], E[iT, id_])

# Directional reflectivity vs wavelength at 5 deg incidence and 100 um (not T-dependent source)
R = src['calculated_R_directional']               # lambda x theta x d
ia = np.abs(a - 5).argmin()
id_ = np.abs(d - 100).argmin()
R_vs_wl = R[:, ia, id_]

# Temperature-dependent source: the temperature is the last dimension
srcT = db['Ag']['Ferrera-Tdepend']
RT = srcT['calculated_R_directional']             # lambda x theta x d x T
TT = srcT['calculation_temperature_K'].flatten()
R_vs_wl_at_T = RT[:, 0, 0, np.abs(TT - 500).argmin()]

# Lumped-body equilibrium temperature (hemispherical A and E) vs concentration, for the first thickness
C = src['calculation_solar_concentration_multipliers'].flatten()
Teq = src['calculated_lumped_equilibrium_temperature_hemiA_hemiE_K']   # 1 x d x C
Teq_vs_C = Teq[0, 0, :]
```

## 5. Spread between sources (`source_spread`)

Every material has a structure `source_spread`, which quantifies the variability of the calculated
hemispherical emissivity between the sources of that material. All statistics are calculated
separately for every temperature and thickness. The temperature-dependent (`-Tdepend`) entries are
not included. `spread[material]` is `None` if a material has fewer than two sources, or if its
sources are not on the same temperature and thickness grids.

| Field | Content | Dimensions |
|---|---|---|
| `description` | text description | |
| `authors` | list of the sources included | |
| `temperature_K` | temperature grid [K] | T |
| `thickness_um` | thickness grid [µm] | d |
| `min`, `max` | lowest and highest emissivity of any source | T × d |
| `mean`, `median` | mean and median emissivity across the sources | T × d |
| `std` | standard deviation across the sources (normalized by N − 1) | T × d |
| `mad` | median absolute deviation from the median | T × d |
| `range_abs` | `max − min` | T × d |
| `range_rel_pct` | `100 × (max − min) / median` [%] | T × d |

```python
sp = spread['Al']
print(sp['authors'])

T = sp['temperature_K'].flatten()
d = sp['thickness_um'].flatten()
iT = np.abs(T - 500).argmin()
id_ = np.abs(d - 10000).argmin()
print(sp['min'][iT, id_], sp['median'][iT, id_], sp['max'][iT, id_])
print(sp['range_rel_pct'][iT, id_])               # relative range at 500 K and 10000 um [%]
```

The relative range becomes very large when the median emissivity is small, so read it together with
the absolute range.

## 6. Loading a single material (low memory)

If the memory is limited, read only one material. Silicon is the largest (about 1 GB); most
materials need far less.

```python
import h5py
import numpy as np

path = 'SpaderOptics.mat'

def to_str(f, ref):
    return ''.join(chr(int(c)) for c in np.array(f[ref]).flatten())

def get(f, ref):
    arr = np.array(f[ref])
    if arr.dtype == np.uint16:
        return ''.join(chr(int(c)) for c in arr.flatten())
    return arr.T

def load_material(path, name):
    with h5py.File(path, 'r') as f:
        root = f['SpaderOptics']
        names = [to_str(f, r) for r in root['material'][:].flatten()]
        i = names.index(name)
        ds = f[root['dataset'][:].flatten()[i]]
        out = {}
        for j in range(ds['author'].size):
            rec = {fld: get(f, ds[fld][:].flatten()[j]) for fld in ds.keys()}
            out[rec['author']] = rec
        return out

al = load_material(path, 'Al')
```

The materials in the database are listed by `list(db.keys())` after running `SpaderOptics_load.py`.

## 7. Other tools

- **MATLAB:** see `README_MATLAB.md`. In short, `load('SpaderOptics.mat')` loads everything, and
  `matfile('SpaderOptics.mat')` reads parts of the file. For example,
  `m = matfile('SpaderOptics.mat'); material = m.SpaderOptics(1,i);` reads only material number `i`.
- **Other languages:** because the file is HDF5, it can be read with any HDF5 library.

## 8. Notes

- Units are given in the field names (`_um`, `_K`, `_W_m2`, `_W_mK`, `_deg`).
- A source is identified by its `author` (for example `Rakic`, `Ferrera-298K`, `Ferrera-Tdepend`).
- The data describe materials in the conditions of the original sources, so a difference between
  sources does not mean that one of them is wrong.
