# SpaderOptics database: access with MATLAB

This file explains how to read `SpaderOptics.mat` in MATLAB. For access with open-source tools
(Python), see `README_PYTHON.md`. The field names, units and array dimensions are the same in both.

No toolbox is needed, only base MATLAB (the file is in `-v7.3` format, MATLAB R2006b or later).

## 1. Loading

```matlab
load('SpaderOptics.mat')        % creates the variable SpaderOptics
```

The complete database needs about 2.8 GB of memory (the file is about 1.4 GB on disk). To read
only part of the file, see section 6.

## 2. Structure

`SpaderOptics` is a 1 x N structure array with one element per material:

```
SpaderOptics(i).material        name of the material (text)
SpaderOptics(i).dataset(j)      the j-th literature source of the material (structure array)
SpaderOptics(i).source_spread   statistics of the emissivity across the sources (structure)
```

```matlab
{SpaderOptics.material}                          % list of all materials

i  = find(strcmp({SpaderOptics.material}, 'Al')); % index of a material
ds = SpaderOptics(i).dataset;                     % all sources of this material
{ds.author}                                       % list of the sources (authors)

j   = find(strcmp({ds.author}, 'Rakic'));         % index of a source
src = ds(j);                                      % all data of this source
fieldnames(src)                                   % list of all fields
```

Sources whose name ends with `-Tdepend` (for example `Ferrera-Tdepend`) combine the datasets that
the authors measured at several temperatures. In these entries some arrays have an additional
temperature dimension.

Notation below: `λ` wavelength, `θ` angle, `d` thickness, `T` temperature, `C` solar
concentration multiplier.

## 3. Data of one source

The examples use `src` as defined above.

### 3.1 Bibliographic information

```matlab
src.author                      % name of the source
src.reference                   % full reference
src.source_url                  % DOI or URL
src.upper_temperature_limit_K   % highest temperature of the data [K]
```

### 3.2 Original data and resampled optical constants

```matlab
% original data from the literature
wl_orig = src.original_wavelength_um;       % [um]
n_orig  = src.original_n_complex;           % complex refractive index n + ik
n = real(n_orig);   k = imag(n_orig);

src.original_temperature_K    % number, or text if the source gives no numerical temperature
src.original_thickness_um     % number, or '' if the source gives no thickness

% resampled to the solar wavelength grid used in all calculations
wl = src.resampled_wavelength_um;           % [um]
nk = src.resampled_n_complex;
```

For `-Tdepend` sources, `resampled_n_complex` has one column per source temperature.

### 3.3 Calculation grids

```matlab
a = src.calculation_incident_angles_deg;               % theta [deg]
d = src.calculation_thickness_um;                      % thickness [um]
T = src.calculation_temperature_K;                     % temperature [K]
C = src.calculation_solar_concentration_multipliers;   % concentration multipliers
```

### 3.4 Optical properties

| Field | Dimensions (not T-dependent) | Dimensions (`-Tdepend`) |
|---|---|---|
| `calculated_R_directional`, `calculated_T_directional`, `calculated_A_directional` | λ × θ × d | λ × θ × d × T |
| `calculated_R_hemispherical`, `calculated_T_hemispherical`, `calculated_A_hemispherical` | λ × d | λ × d × T |
| `calculated_E_directional` | T × θ × d | T × θ × d |
| `calculated_E_hemispherical` | T × d | T × d |

(R reflectivity, T transmissivity, A absorptivity, E emissivity.)

```matlab
% hemispherical emissivity at about 500 K and 10000 um
[~, iT] = min(abs(T - 500));
[~, id] = min(abs(d - 10000));
E = src.calculated_E_hemispherical;        % T x d
E(iT, id)

% hemispherical emissivity vs temperature for all thicknesses
plot(T, E); xlabel('Temperature (K)'); ylabel('Hemispherical emissivity')

% directional reflectivity vs wavelength at 5 deg and 100 um (not T-dependent source)
[~, ia] = min(abs(a - 5));
[~, id] = min(abs(d - 100));
R = src.calculated_R_directional;          % lambda x theta x d
plot(wl, R(:, ia, id))

% T-dependent source: temperature is the last dimension
srcT = SpaderOptics(find(strcmp({SpaderOptics.material}, 'Ag'))).dataset;
srcT = srcT(strcmp({srcT.author}, 'Ferrera-Tdepend'));
RT   = srcT.calculated_R_directional;      % lambda x theta x d x T
[~, iT] = min(abs(srcT.calculation_temperature_K - 500));
plot(srcT.resampled_wavelength_um, RT(:, 1, 1, iT))
```

### 3.5 Fluxes (W m⁻²)

| Field | Dimensions (not T-dependent) | Dimensions (`-Tdepend`) |
|---|---|---|
| `calculated_directional_Reflected_flux_W_m2`, `..._Transmitted_...`, `..._Absorbed_...` | θ × d | T × θ × d |
| `calculated_hemispherical_Reflected_flux_W_m2`, `..._Transmitted_...`, `..._Absorbed_...` | d | T × d |
| `calculated_hemispherical_Emitted_flux_W_m2` | T × d | T × d |

```matlab
qabs = src.calculated_directional_Absorbed_flux_W_m2;   % theta x d
plot(a, qabs(:, id))
```

### 3.6 Equilibrium and slab temperatures, conductive flux

| Field | Dimensions |
|---|---|
| `calculated_lumped_equilibrium_temperature_dirA_hemiE_K` | θ × d × C |
| `calculated_lumped_equilibrium_temperature_hemiA_hemiE_K` | 1 × d × C |
| `calculated_Conductive_flux_dirA_hemiE_W_m2`, `calculated_Conductive_flux_hemiA_hemiE_W_m2` | θ × d × C, and 1 × d × C |
| `calculated_slab_equilibrium_temperature_{front,rear,mid}_dirA_hemiE_K` | θ × d × C |
| `calculated_slab_equilibrium_temperature_{front,rear,mid}_hemiA_hemiE_K` | 1 × d × C |

`dirA_hemiE` uses the directional absorptivity and the hemispherical emissivity; `hemiA_hemiE`
uses the hemispherical absorptivity and the hemispherical emissivity.

```matlab
% lumped-body equilibrium temperature vs concentration (hemispherical A and E), first thickness
Teq = src.calculated_lumped_equilibrium_temperature_hemiA_hemiE_K;   % 1 x d x C
plot(C, squeeze(Teq(1, 1, :)))
xlabel('Solar concentration multiplier'); ylabel('Temperature (K)')

% slab front temperature, directional A, normal incidence, first thickness
Tfront = src.calculated_slab_equilibrium_temperature_front_dirA_hemiE_K;   % theta x d x C
plot(C, squeeze(Tfront(1, 1, :)))
```

### 3.7 Thermal conductivity and coverage metrics

```matlab
plot(src.thermal_conductivity_temperature_grid_K, src.thermal_conductivity_W_mK)
src.thermal_conductivity_reference

src.metrics_solar_coverage_ratio     % fraction of the AM0 spectrum covered (a number)
src.metrics_thermal_coverage_ratio   % thermal coverage ratio for every temperature (vector over T)
```

## 4. Spread between sources (`source_spread`)

Every material has a structure `source_spread` with statistics of the calculated hemispherical
emissivity across its sources. The statistics are calculated separately for every temperature and
thickness. The temperature-dependent (`-Tdepend`) entries are not included. The structure is
empty if the material has fewer than two sources, or if its sources are not on the same
temperature and thickness grids.

| Field | Content | Dimensions |
|---|---|---|
| `description` | text description | |
| `authors` | sources included | |
| `temperature_K` | temperature grid [K] | T |
| `thickness_um` | thickness grid [µm] | d |
| `min`, `max` | lowest and highest emissivity of any source | T × d |
| `mean`, `median` | mean and median across the sources | T × d |
| `std` | standard deviation across the sources (normalized by N − 1) | T × d |
| `mad` | median absolute deviation from the median | T × d |
| `range_abs` | `max − min` | T × d |
| `range_rel_pct` | `100 × (max − min) / median` [%] | T × d |

```matlab
sp = SpaderOptics(i).source_spread;

if isfield(sp, 'min')               % check that a spread exists for this material
    sp.authors                      % sources included

    [~, iT] = min(abs(sp.temperature_K - 500));
    [~, id] = min(abs(sp.thickness_um - 10000));
    sp.min(iT, id), sp.median(iT, id), sp.max(iT, id)
    sp.range_rel_pct(iT, id)        % relative range at 500 K and 10000 um [%]

    % envelope (min-max) and median vs temperature for one thickness
    Ts = sp.temperature_K(:);
    fill([Ts; flipud(Ts)], [sp.min(:, id); flipud(sp.max(:, id))], ...
         [0.2 0.4 0.8], 'FaceAlpha', 0.25, 'EdgeColor', 'none'); hold on
    plot(Ts, sp.median(:, id), '--')
    xlabel('Temperature (K)'); ylabel('Hemispherical emissivity')
end
```

The relative range becomes very large when the median emissivity is small, so read it together
with the absolute range.

## 5. Interactive viewer

The folder `Step6_SpaderOptics_App_Development` (if included in the deposit) contains the MATLAB
App `SpaderOpticsApp`, which plots the quantities of the database for a chosen material and
source, and shows the spread between sources.

## 6. Reading only a part of the file

The `matfile` function reads parts of a `-v7.3` file without loading the rest:

```matlab
m = matfile('SpaderOptics.mat');
size(m, 'SpaderOptics')           % size of the structure array
mat = m.SpaderOptics(1, i);       % loads only material number i (all its sources)
```

The structure array has one row, so the index is `(1, i)`. The element contains all sources of
that material. Silicon is the largest (about 1 GB); most materials need far less.

## 7. Notes

- Units are given in the field names (`_um`, `_K`, `_W_m2`, `_W_mK`, `_deg`).
- A source is identified by its `author` (for example `Rakic`, `Ferrera-298K`, `Ferrera-Tdepend`).
- The data describe materials in the conditions of the original sources, so a difference between
  sources does not mean that one of them is wrong.
