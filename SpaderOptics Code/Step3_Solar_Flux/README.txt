STEP 3 - Solar-matter interaction

Code:  Step3_Solar_Matter_Interaction.m   (main script)
       LoadSolar.m                        (function called by the main script)

EXTERNAL DATA (not included; ASTM standard, not redistributed)
The extraterrestrial (AM0) reference solar spectrum, ASTM E490 (<edition, ASTM E490-00a(2019)>).
A data file is available from NREL:
  https://www.nrel.gov/grid/solar-resource/spectra-astm-e490
Save it in this folder as an Excel file named

  solar_spectrum.xlsx

Retrieved on: 4 May 2026.

RUN
Run Step3_Solar_Matter_Interaction.m only. It loads the spectrum automatically by calling
LoadSolar.m. Do not clear the workspace afterwards.

PLEASE CITE, additionally, the ASTM E490 standard and the source from which you downloaded the data file.
