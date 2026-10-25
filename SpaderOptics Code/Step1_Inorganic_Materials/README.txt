STEP 1 - Inorganic materials: import of the refractive indices

Code:  Step1_Import_RefractiveIndices_for_InorganicMaterials.m

EXTERNAL DATA (not included; created by other authors)
Download the folders below from the RefractiveIndex.INFO database:
  https://github.com/polyanskiy/refractiveindex.info-database
  path: database/data/main/

  Ag, Al, Al2O3, Au, Cu, Cu2O, Fe, Fe2O3, Fe3O4, Ni, Ti, TiO2

(Use the folder names exactly as the code expects them; check the list at the top of the .m file.)

HOW TO ARRANGE THEM
In the repository each material folder contains two subfolders, "nk" and "n2".
Only the data in "nk" are used. Delete the intermediate "nk" folder and copy its
contents (one YAML file per author) directly into the material folder:

   repository:   Ag/nk/<author>.yml
   this folder:  Ag/<author>.yml

Do the same for every material above. The folders are placed inside this folder
(Step1_Inorganic_Materials), next to the .m file.

The database in the Zenodo deposit was built from the RefractiveIndex.INFO data
retrieved on 4 May 2026 (last commit on 3 May 2026). The upstream repository may have changed since then.
For an exact reproduction, use the commit of that date:
  <e5203c860cf87b9b0daf95831b8d436e5e63c04d>

RUN
Open MATLAB in this folder and run Step1_Import_RefractiveIndices_for_InorganicMaterials.m.
Do not clear the workspace afterwards: Step 5 uses the arrays created here.

PLEASE CITE, additionally, the RefractiveIndex.INFO database and the original sources of the data.
