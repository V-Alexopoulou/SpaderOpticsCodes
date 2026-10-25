# SpaderOptics: code

Code used to build the **SpaderOptics database** (optical, solar-interaction, thermal-emission and
equilibrium-temperature data of materials for space applications) and the **SpaderOptics App**.

- Database (Zenodo): "10.5281/zenodo.23643910"
- This code (archived on Zenodo): <software DOI>
- Paper: <reference / DOI>

> The database file `SpaderOptics.mat` (about 1.4 GB) is **not** stored in this repository.
> Download it from Zenodo.

## Two ways to use this repository

**A. Use the database only (no calculation).** Download `SpaderOptics.mat` from Zenodo.
- With the App: copy it into `Step6_SpaderOptics_App_Development` and run `SpaderOpticsApp.mlapp`.
- With scripts: follow `Step7_SpaderOptics_Load/README_MATLAB.md` or `README_PYTHON.md`.

**B. Rebuild the database (Steps 1-5).** See below.

## Contents

| Folder | Content | External data to add |
|---|---|---|
| `Step1_Inorganic_Materials` | import of the refractive indices of inorganic materials | RefractiveIndex.INFO folders |
| `Step2_Organic_Materials` | import of the refractive indices of the other materials (Ge, Si, SiC, SiO2, TiC, ZrO2, C, Kapton, PMMA) | RefractiveIndex.INFO folders |
| `Step3_Solar_Flux` | solar-matter interaction (`Step3_Solar_Matter_Interaction.m`, `LoadSolar.m`) | AM0 spectrum `solar_spectrum.xlsx` |
| `Step4_Emissivity` | emissivity and equilibrium temperatures | none |
| `Step5_SpaderOptics_Database_Development` | assembles the final database `SpaderOptics.mat` | none |
| `Step6_SpaderOptics_App_Development` | the SpaderOptics App (`SpaderOpticsApp.mlapp`) | the database |
| `Step7_SpaderOptics_Load` | Python loader and MATLAB/Python readmes | the database |

Each folder has its own `README.txt` with the exact instructions.

## Requirements

- MATLAB R2025b (developed and tested with this version)
- Python 3 with `numpy` and `h5py` (only for loading the database with Python)

## Rebuilding the database (Steps 1-5)

1. Add the external data (next section) to the folders of Steps 1, 2 and 3.
2. Open MATLAB and run, **in this order**: Step 1, Step 2, Step 3, Step 4, Step 5. Run each script by itself.
3. **Do not clear the workspace between the steps.** Each script uses the arrays created by the previous ones.
4. Step 5 creates `SpaderOptics.mat` in its own folder.

## External data (not included)

The following data were created by other authors and are therefore not redistributed here.
Please download them from the original sources and cite them.

| Data | Source | Used in |
|---|---|---|
| Refractive indices (YAML files) | RefractiveIndex.INFO database, https://github.com/polyanskiy/refractiveindex.info-database | Steps 1, 2 |
| AM0 solar spectrum | ASTM E490 (<edition>), data file from NREL, https://www.nrel.gov/grid/solar-resource/spectra-astm-e490 | Step 3 |

The deposited database was created with the data retrieved on **4 May 2026**. The original
repositories may have changed since then, so a rebuild with newer data may give slightly different results.
To reproduce the deposited database exactly, use the refractiveindex.info-database commit
`<e5203c860cf87b9b0daf95831b8d436e5e63c04d>` (last commit on 3 May 2026).

## Licence

The code is released under the GNU General Public License v3.0 (see `LICENSE`).
The database is released separately on Zenodo under CC BY 4.0.

## How to cite

Please cite the paper, the database and this code (see `CITATION.cff`).

## Funding

This work was funded by the European Union under grant agreement number 101223931. Its contents are
the sole responsibility of the SPADER consortium and do not necessarily reflect the views of the European Commission.
