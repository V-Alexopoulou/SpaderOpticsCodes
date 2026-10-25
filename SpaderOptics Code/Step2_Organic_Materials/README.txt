STEP 2 - Other materials (semiconductors, ceramics, polymers, carbon): import of the refractive indices

Code:  Step2_Import_RefractiveIndices_for_OrganicMaterials.m

EXTERNAL DATA (not included; created by other authors)
Download the folders below from the RefractiveIndex.INFO database:
  https://github.com/polyanskiy/refractiveindex.info-database

  from database/data/main/ :
     Ge, Si, SiC, SiO2, TiC, ZrO2
     C   (carbon)  -> keep ONLY the amorphous-carbon sources: Arakawa, Hagemann, Larruquert
                      (the C folder also contains other forms, e.g. diamond, which are not used)

  from database/data/other/commercial polymers/ :
     Kapton

  from database/data/organic/ :
     (C5H8O2)n-poly(methyl methacrylate)   -> used as PMMA
     <rename the folder to the name expected by the code, if different>

HOW TO ARRANGE THEM
Only the data in the "nk" subfolders are used. Delete the intermediate "nk" folder and copy
its contents (one YAML file per author) directly into the material folder:

   repository:   Si/nk/<author>.yml
   this folder:  Si/<author>.yml

The folders are placed inside this folder (Step2_Organic_Materials), next to the .m file.

The database in the Zenodo deposit was built from the RefractiveIndex.INFO data
retrieved on 4 May 2026 (last commit on 3 May 2026). For an exact reproduction use the commit of that date:
  <e5203c860cf87b9b0daf95831b8d436e5e63c04d>

RUN
Open MATLAB in this folder and run Step2_Import_RefractiveIndices_for_OrganicMaterials.m.
Do not clear the workspace afterwards.

PLEASE CITE, additionally, the RefractiveIndex.INFO database and the original sources of the data.
