%% Step4_Emissivity_EquilibriumTemperature.m
%  =========================================================================
%  STEP 4 — Emissivity, lumped equilibrium temperature, and slab
%           (two-surface, conduction-coupled) temperature computation
%
%  Uses the optical data from the `solar_flux` struct (produced by
%  compute_solar_flux.m) to compute, for every material / literature
%  source:
%
%  1. SPECTRAL EMISSIVITY (Kirchhoff's law)
%        eps_lambda(lambda) = alpha_lambda(lambda)
%
%  2. TOTAL (DIRECTIONAL AND HEMISPHERICAL) EMISSIVITY vs TEMPERATURE
%     Planck-weighted ratio of real-body to blackbody spectral radiance,
%     evaluated over the material-specific temperature grid described
%     below.
%
%  3. LUMPED EQUILIBRIUM TEMPERATURE
%     Single-temperature energy balance  Q_abs = Q_emit,hemispherical,
%     solved by a 1-D sign-change search over the temperature grid.
%
%  4. SLAB (TWO-SURFACE) TEMPERATURE
%     Steady-state 1-D conduction between a front surface (absorbing
%     solar flux) and a rear surface, both radiating hemispherically,
%     for a free-standing plane-parallel slab of given thickness.
%
%  PHYSICAL SCOPE — WHAT IS AND IS NOT COMPUTED
%  ----------------------------------------------------------------------
%  Directional *emission* is not used anywhere in this script to build a
%  standalone emitted flux. eps_directional(theta) is a well-defined
%  angle-resolved material property (Kirchhoff's law applied at a given
%  angle), but eps_directional(theta)*sigma*T^4 is NOT a physically
%  complete flux: sigma*T^4 is a whole-hemisphere-integrated quantity,
%  while eps_directional is a per-steradian (radiance-type) ratio.
%  Turning a directional emissivity into an actual emitted flux requires
%  specifying a receiving geometry (a target surface, distance, and
%  solid angle) that is undefined for an isolated fragment radiating
%  into the omnidirectional deep-space background. Accordingly, only
%  HEMISPHERICAL emission is used in every energy balance below; the
%  only "directional" quantity that ever appears in an equilibrium
%  calculation is the *absorbed* solar flux (which is well-posed,
%  because there genuinely is one real incident beam direction), and it
%  is explicitly projected by cos(theta_i) onto the real, tilted slab
%  surface before being compared against a per-real-surface-area
%  emitted flux.
%
%  TEMPERATURE GRID
%  ----------------------------------------------------------------------
%  Each material's temperature grid runs from the lowest temperature
%  available in its literature thermal-conductivity data up to
%  min(T_max, melting/degradation point). Properties are therefore never
%  evaluated, extended, or frozen above the melting/degradation point:
%  if the true equilibrium would lie above it, the corresponding
%  root-finding search simply returns NaN, which should be read as
%  "this material melts/degrades under this condition" (check
%  melting_point_K for the relevant threshold).
%
%  PREREQUISITES
%  ----------------------------------------------------------------------
%  `solar_flux`, `solar_flux_valid`, `solar_flux_best` must exist in the
%  workspace (or on disk as .mat files), produced by compute_solar_flux.m.
%  `wl_sun`, `I_sun` (the AM0 grid/spectrum) must also be in the workspace.
%  =========================================================================

clc;

%% ==== USER SETTINGS ======================================================

T_min  = 5;      % Absolute floor of any temperature grid [K]
T_max  = 5000;   % Absolute ceiling of any temperature grid [K]
N_temp = 50;     % Number of temperature points per material-specific grid

% Solar concentration multipliers (suns), as used throughout the paper.
multipliers = [1, 2:2:20, 30:10:100, 200:100:500, 1000, 2000];

% ==========================================================================

% --- Physical constants ---
sigma_SB = 5.670374419e-8;      % Stefan-Boltzmann constant [W/m^2/K^4]

%% ==== CHECK SOLAR FLUX DATASETS ==========================================

if ~exist('solar_flux', 'var') || isempty(solar_flux)
    if exist('solar_flux.mat', 'file')
        S = load('solar_flux.mat', 'solar_flux');
        solar_flux = S.solar_flux;
    else
        error(['`solar_flux` not found. Run compute_solar_flux.m first.']);
    end
end

%% ==== MELTING / DEGRADATION POINTS ========================================
%
% Values in Kelvin. Materials without a defined melting point are
% assigned NaN; for those, the temperature grid runs all the way to
% T_max (see the "material-specific temperature range" block below),
% and this absence of a known validity limit should be treated as a
% known limitation for that material.

melting_points = struct();

melting_points.C          = 3925;
melting_points.Ge         = 1210.6;
melting_points.Kapton     = 873.2;
melting_points.PMMA       = 332.2;
melting_points.Si         = 1685.2;
melting_points.SiC        = 3100;
melting_points.SiO2       = 1986.0;
melting_points.TiC        = 3523;
melting_points.ZrO2       = 2973.0;
melting_points.Ag         = 1234.0;
melting_points.Al         = 933.2;
melting_points.Al2O3      = 2323;
melting_points.Au         = 1336.2;
melting_points.Cu         = 1356;
melting_points.Cu2O       = 423;
melting_points.Fe         = 1810;
melting_points.Fe2O3      = 950;
melting_points.Fe3O4      = 858;
melting_points.Ni         = 1726;
melting_points.Ti         = 1953;
melting_points.TiO2       = 2133;

%% ==== THERMAL CONDUCTIVITY DATA ===========================================
%
% Literature temperature-conductivity pairs (T [K], k [W/m/K]), retained
% exactly as tabulated. Used only by the Fourier (two-surface) model;
% interpolation to the model's temperature grid is performed later and
% is NOT extrapolated beyond the tabulated range (see
% evaluateThermalConductivity).

thermal_conductivity = struct();

thermal_conductivity.C.T = [10:10:100 150:50:250 273.2 300:50:400 500:100:1000 1200:200:3800 melting_points.C];
thermal_conductivity.C.k = [0.49 2.5 6 11 17 24 32 41 49 58 94 120 131 131 129 124 118 106 95 85 77 70 64 55 49 45 42 40 38 36 35 33 31 28 25 20 11.3 11.3];
thermal_conductivity.C.ref = "Touloukian, Y.S.; Powell, R.W.; Ho, C.Y.; Klemens, P.G. Thermophysical Properties of Matter - The TPCR Data Series -- Vol 2. Thermal Conductivity - Nonmetallic Solids. IFI/Plenum, New York - Washington (1970)";

thermal_conductivity.Ge.T = [5:16 18 20:5:50 60:10:100 150:50:250 273.2 300 350 400:100:1100 melting_points.Ge];
thermal_conductivity.Ge.k = [1160 1390 1550 1660 1730 1770 1790 1800 1790 1770 1730 1690 1590 1490 1270 1080 920 788 695 619 487 393 325 270 232 132 96.8 74.9 66.7 59.9 49.5 43.2 33.8 27.3 22.7 19.8 18.2 17.4 17.0 17.4];
thermal_conductivity.Ge.ref = "Touloukian, Y.S.; Powell, R.W.; Ho, C.Y.; Klemens, P.G. Thermophysical Properties of Matter - The TPCR Data Series -- Vol 1. Thermal Conductivity - Metallic Elements and Alloys. IFI/Plenum, New York - Washington (1970)";

thermal_conductivity.Kapton.T = [5 773.2 melting_points.Kapton];
thermal_conductivity.Kapton.k = [0.25 0.25 0.65];
thermal_conductivity.Kapton.ref = "S. Venkatachalam, M. Depriester, A. Hadj Sahraoui, B. Capoen, M.R. Ammar, D. Hourlier. Thermal conductivity of Kapton-derived carbon. Carbon 114, 134 -- 140 (2017). https://doi.org/10.1016/j.carbon.2016.11.072";

thermal_conductivity.PMMA.T = [5 6 9 15 19 21 24 183.3 184.4 244.8 245.6 298.2 300.3 332.1 melting_points.PMMA];
thermal_conductivity.PMMA.k = [0.06 0.06 0.06 0.065 0.07 0.077 0.09 0.105 0.106 0.134 0.135 0.192 0.154 0.158 0.157];
thermal_conductivity.PMMA.ref = "Touloukian, Y.S.; Powell, R.W.; Ho, C.Y.; Klemens, P.G. Thermophysical Properties of Matter - The TPCR Data Series -- Vol 2. Thermal Conductivity - Nonmetallic Solids. IFI/Plenum, New York - Washington (1970)";

thermal_conductivity.Si.T = [5:16 18 20:5:50 60:10:100 150:50:250 273.2 300 350 400:100:1600 melting_points.Si];
thermal_conductivity.Si.k = [424 686 991 1340 1720 2110 2480 2870 3250 3600 3930 4220 4670 4940 5140 4810 4130 3530 3060 2680 2110 1680 1340 1080 884 409 264 190 168 146 119 98.9 76.2 61.9 50.8 42.2 35.9 31.2 27.9 25.7 24.4 23.5 22.7 22.1 22];
thermal_conductivity.Si.ref = "Touloukian, Y.S.; Powell, R.W.; Ho, C.Y.; Klemens, P.G. Thermophysical Properties of Matter - The TPCR Data Series -- Vol 1. Thermal Conductivity - Metallic Elements and Alloys. IFI/Plenum, New York - Washington (1970)";

thermal_conductivity.SiC.T = [5 7.85 8.22 9.82 11.2 13.8 16.2 19 20.7 23.2 28.3 41.7 48.5 100 102 104 115 145 152 165 173 185 270 873.2 1073.2 1273.2 1473.2 1673.2 melting_points.SiC];
thermal_conductivity.SiC.k = [2.75 7.07 9.08 10.3 14.7 21 30.1 45.5 51 62.5 93 115 155 89.5 105 89.7 98.7 98.5 93.7 89.7 85.7 83.6 65.5 18.4 15.9 13.8 12.1 10.9 0.71];
thermal_conductivity.SiC.ref = "Touloukian, Y.S.; Powell, R.W.; Ho, C.Y.; Klemens, P.G. Thermophysical Properties of Matter - The TPCR Data Series -- Vol 2. Thermal Conductivity - Nonmetallic Solids. IFI/Plenum, New York - Washington (1970)";

thermal_conductivity.SiO2.T = [5 7:1:13 15:5:50 60:10:100 150:50:250 273 300:50:500 600:100:800 melting_points.SiO2];
thermal_conductivity.SiO2.k = [400 900 1210 1500 1650 1680 1630 1520 1250 720 460 318 233 179 143 118 83 66 54 45 39 23.1 16.4 12.7 11.6 10.4 8.8 7.6 6.7 6 5 4.47 4.2 4.2];
thermal_conductivity.SiO2.ref = "Touloukian, Y.S.; Powell, R.W.; Ho, C.Y.; Klemens, P.G. Thermophysical Properties of Matter - The TPCR Data Series -- Vol 2. Thermal Conductivity - Nonmetallic Solids. IFI/Plenum, New York - Washington (1970)";

thermal_conductivity.TiC.T = [373.2 471.2 513.2 558.2 653.2 758.2 793.2 843.2 872.2 873.2 933.2 943.2 963.2 978.2 988.2 1018.2 1028.2 1048.2 1073.2 1174.2 1238.2 1263.2 1388.2 1448.2 1518 1602 1813 1913 1994 2074 2125 2254 2303 2333 melting_points.TiC];
thermal_conductivity.TiC.k = [24.1 20.9 19.8 17.8 15.6 34.9 34.3 37 37.4 38.5 38.7 39.1 38.9 39.1 40 41.8 39.1 41.4 41.8 39.3 43.5 38.3 39.7 42.7 41.4 43.9 47.7 45.2 49 48.1 47.7 47.7 47.3 48.1 48.1];
thermal_conductivity.TiC.ref = "Touloukian, Y.S.; Powell, R.W.; Ho, C.Y.; Klemens, P.G. Thermophysical Properties of Matter - The TPCR Data Series -- Vol 2. Thermal Conductivity - Nonmetallic Solids. IFI/Plenum, New York - Washington (1970)";

thermal_conductivity.ZrO2.T = [369.2 747.2 649.2 701.2 749.2 818.2 883.2 949.2 1001.2 1061.2 1125.2 1189.2 1243.2 1283.2 1307.2 1325.2 1333.2 1353.2 1389.2 1403.2 1425.2 1433.2 1447.2 1487.2 1503.2 1553.2 1569.2 melting_points.ZrO2];
thermal_conductivity.ZrO2.k = [1.63 1.7 1.72 1.76 1.77 1.8 1.82 1.84 1.86 1.88 1.91 1.92 1.94 1.97 1.97 1.97 1.97 1.99 2 2.01 2.01 2.01 2.02 2.03 2.05 2.05 2.06 2.06];
thermal_conductivity.ZrO2.ref = "Touloukian, Y.S.; Powell, R.W.; Ho, C.Y.; Klemens, P.G. Thermophysical Properties of Matter - The TPCR Data Series -- Vol 2. Thermal Conductivity - Nonmetallic Solids. IFI/Plenum, New York - Washington (1970)";

thermal_conductivity.Ag.T = [5:16 18 20:5:50 60:10:100 150:50:250 273.2 300 350 400:100:1200 melting_points.Ag];
thermal_conductivity.Ag.k = [17200 18700 19300 19000 18100 16800 15400 13900 12400 10900 9600 8500 6600 5100 2950 1930 1370 1050 840 700 550 497 471 460 450 432 430 428 428 427 424 420 413 405 397 389 382 374 366 358 355];
thermal_conductivity.Ag.ref = "Touloukian, Y.S.; Powell, R.W.; Ho, C.Y.; Klemens, P.G. Thermophysical Properties of Matter - The TPCR Data Series -- Vol 1. Thermal Conductivity - Metallic Elements and Alloys. IFI/Plenum, New York - Washington (1970)";

thermal_conductivity.Al.T = [5:16 18 20:5:50 60:10:100 150:50:250 273.2 300 350 400:100:900 melting_points.Al];
thermal_conductivity.Al.k = [18800 21300 22900 23800 24000 23600 22700 2500 20200 19900 17600 16300 13800 11700 7740 3180 3490 2380 1680 1230 734 532 414 344 302 248 237 235 236 237 240 240 237 232 226 220 213 211];
thermal_conductivity.Al.ref = "Touloukian, Y.S.; Powell, R.W.; Ho, C.Y.; Klemens, P.G. Thermophysical Properties of Matter - The TPCR Data Series -- Vol 1. Thermal Conductivity - Metallic Elements and Alloys. IFI/Plenum, New York - Washington (1970)";

thermal_conductivity.Al2O3.T = [100:50:250 273.2 300:50:400 500:100:2100 melting_points.Al2O3];
thermal_conductivity.Al2O3.k = [133 77 55 43.4 39.7 36 30.7 26.4 20.2 15.8 12.6 10.4 8.9 7.85 7.1 6.55 6.13 5.83 5.66 5.56 5.54 5.59 5.74 6 6.44 6.44];
thermal_conductivity.Al2O3.ref = "Touloukian, Y.S.; Powell, R.W.; Ho, C.Y.; Klemens, P.G. Thermophysical Properties of Matter - The TPCR Data Series -- Vol 2. Thermal Conductivity - Nonmetallic Solids. IFI/Plenum, New York - Washington (1970)";

thermal_conductivity.Au.T = [5:16 18 20:5:50 60:10:100 150:50:250 273.2 300 350 400:100:1300 melting_points.Au];
thermal_conductivity.Au.k = [2070 2370 2600 2750 2820 2820 2770 2670 2550 2410 2260 2090 1770 1500 1020 760 610 520 460 420 380 358 352 348 345 335 327 320 318 315 313 312 309 304 298 292 285 278 271 262 251 247];
thermal_conductivity.Au.ref = "Touloukian, Y.S.; Powell, R.W.; Ho, C.Y.; Klemens, P.G. Thermophysical Properties of Matter - The TPCR Data Series -- Vol 1. Thermal Conductivity - Metallic Elements and Alloys. IFI/Plenum, New York - Washington (1970)";

thermal_conductivity.Cu.T = [5:16 18 20:5:50 60:10:100 150:50:250 273.2 300 350 400:100:1300 melting_points.Cu];
thermal_conductivity.Cu.k = [13800 15900 17700 18900 19500 19600 19300 18300 17600 16600 15600 14500 12400 10500 6800 4300 2900 2030 1530 1220 850 670 570 514 483 428 413 404 401 398 394 392 388 383 377 371 364 357 350 342 334 330];
thermal_conductivity.Cu.ref = "Touloukian, Y.S.; Powell, R.W.; Ho, C.Y.; Klemens, P.G. Thermophysical Properties of Matter - The TPCR Data Series -- Vol 1. Thermal Conductivity - Metallic Elements and Alloys. IFI/Plenum, New York - Washington (1970)";

thermal_conductivity.Cu2O.T = [125.2 137 139.3 152.4 155.4 174 177.2 196.9 198.2 219.1 234.7 235 237 239.7 239.8 244.3 296.6 321.6 357 361.8 melting_points.Cu2O];
thermal_conductivity.Cu2O.k = [3.46 5.36 4.98 6.42 6.23 6.62 6.27 6.25 6.57 6.24 6.62 6.45 5.78 6.44 5.95 5.94 5.44 4.97 4.65 4.48 4.48];
thermal_conductivity.Cu2O.ref = "Touloukian, Y.S.; Powell, R.W.; Ho, C.Y.; Klemens, P.G. Thermophysical Properties of Matter - The TPCR Data Series -- Vol 2. Thermal Conductivity - Nonmetallic Solids. IFI/Plenum, New York - Washington (1970)";

thermal_conductivity.Fe.T = [5:16 18 20:5:50 60:10:100 150:50:250 273.2 300 350 400:100:1000 1059 1100 1183 1200:100:1600 1673 1700 1800 melting_points.Fe];
thermal_conductivity.Fe.k = [373 442 513 580 645 705 762 813 834 897 930 956 988 997 936 814 681 555 450 372 265 204 168 146 132 104 94 86.5 83.5 80.3 74.4 69.4 61.3 54.7 48.7 43.3 38 32.6 29.6 29.7 29.9 28.2 29.9 30.9 31.8 32.1 33.4 33.6 34.5 34.6];
thermal_conductivity.Fe.ref = "Touloukian, Y.S.; Powell, R.W.; Ho, C.Y.; Klemens, P.G. Thermophysical Properties of Matter - The TPCR Data Series -- Vol 1. Thermal Conductivity - Metallic Elements and Alloys. IFI/Plenum, New York - Washington (1970)";

thermal_conductivity.Fe2O3.T = [5 8.49 11.71 15.04 21.54 29.56 65.23 81.02 90.26 198.6 261.63 300.51 317.1 335.7 353.9 385.6 453.2 melting_points.Fe2O3];
thermal_conductivity.Fe2O3.k = [13.4 25 37.6 53.4 72.5 80.8 43.6 28.6 24.5 7.7 7.4 7.1 4.44 4.35 4.35 4.31 4.14 4.14];
thermal_conductivity.Fe2O3.ref = "Touloukian, Y.S.; Powell, R.W.; Ho, C.Y.; Klemens, P.G. Thermophysical Properties of Matter - The TPCR Data Series -- Vol 2. Thermal Conductivity - Nonmetallic Solids. IFI/Plenum, New York - Washington (1970)";

thermal_conductivity.Fe3O4.T = [5 8.49 11.71 15.04 21.54 29.56 65.23 81.02 90.26 198.6 261.63 300.51 317.1 335.7 353.9 385.6 453.2 melting_points.Fe3O4];
thermal_conductivity.Fe3O4.k = [13.4 25 37.6 53.4 72.5 80.8 43.6 28.6 24.5 7.7 7.4 7.1 4.44 4.35 4.35 4.31 4.14 4.14];
thermal_conductivity.Fe3O4.ref = "Touloukian, Y.S.; Powell, R.W.; Ho, C.Y.; Klemens, P.G. Thermophysical Properties of Matter - The TPCR Data Series -- Vol 2. Thermal Conductivity - Nonmetallic Solids. IFI/Plenum, New York - Washington (1970)";

thermal_conductivity.Ni.T = [5:16 18 20:5:50 60:10:100 150:50:250 273.2 300 350 400:100:600 631 700:100:1500 melting_points.Ni];
thermal_conductivity.Ni.k = [316 377 436 494 549 600 648 691 730 764 792 815 845 856 813 695 562 463 391 330 263 221 193 172 158 121 106 97 94 90.5 85.0 80.1 72.1 65.5 63.8 65.3 67.4 69.6 71.8 73.9 76.1 78.2 80.4 82.5 82.5];
thermal_conductivity.Ni.ref = "Touloukian, Y.S.; Powell, R.W.; Ho, C.Y.; Klemens, P.G. Thermophysical Properties of Matter - The TPCR Data Series -- Vol 1. Thermal Conductivity - Metallic Elements and Alloys. IFI/Plenum, New York - Washington (1970)";

thermal_conductivity.Ti.T = [5:16 18 20:5:50 60:10:100 150:50:250 273.2 300 350 400:100:1900 melting_points.Ti];
thermal_conductivity.Ti.k = [7.19 8.63 10.1 11.5 12.9 14.4 15.8 17.2 18.6 20 21.4 22.7 25.4 27.9 33.7 38.2 41.1 42.2 41.6 40.1 37.7 35.6 33.9 32.4 31.2 27 24.5 22.9 22.4 21.9 21 20.4 17.7 19.4 19.4 19.7 20.2 20.7 21.3 22 22.8 23.6 24.5 25.3 26.2 27.1 28 28];
thermal_conductivity.Ti.ref = "Touloukian, Y.S.; Powell, R.W.; Ho, C.Y.; Klemens, P.G. Thermophysical Properties of Matter - The TPCR Data Series -- Vol 1. Thermal Conductivity - Metallic Elements and Alloys. IFI/Plenum, New York - Washington (1970)";

thermal_conductivity.TiO2.T = [5 10:1:15 20:5:30 40:10:100 150:50:250 273.2 300:50:400 500:100:1400 melting_points.TiO2];
thermal_conductivity.TiO2.k = [313 1790 2000 2060 2020 1910 1770 1000 540 285 117 66 45 35 30 26.4 23.5 16.8 13.7 11.8 11.1 10.4 9.4 8.5 5.88 5.02 4.39 3.94 3.65 3.46 3.35 3.28 3.23 3.21 3.21];
thermal_conductivity.TiO2.ref = "Touloukian, Y.S.; Powell, R.W.; Ho, C.Y.; Klemens, P.G. Thermophysical Properties of Matter - The TPCR Data Series -- Vol 2. Thermal Conductivity - Nonmetallic Solids. IFI/Plenum, New York - Washington (1970)";

%% ==== MAIN LOOP ==========================================================
%
% Compute emissivity/equilibrium temperatures, then Fourier temperatures,
% for every source in solar_flux (no coverage-based filtering or "best
% source per material" selection is applied).

fprintf('\n============================================================\n');
fprintf('Computing emissivity and equilibrium temperatures\n');
fprintf('============================================================\n');

emissivity = computeEmissivityDataset( ...
    solar_flux, wl_sun, I_sun, multipliers, sigma_SB, ...
    melting_points, thermal_conductivity, T_min, T_max, N_temp);

fourier = computeFourierDataset( ...
    emissivity, thermal_conductivity, multipliers, sigma_SB);

%% ==== SUBSTRATE / ADHESIVE-SYSTEM EXCLUSIONS ==============================
%
% A few literature sources report optical constants for a measured STACK
% (a thin film on a substrate, and/or with an adhesive layer) rather than
% for the pure bulk material. The thermal-conductivity data used by the
% Fourier (two-surface, conduction-coupled) model above is the pure
% material's, so the resulting conductive flux and Fourier front/rear/
% mid-plane temperatures are not physically meaningful for these specific
% sources and are overwritten with NaN below.
%
% This does NOT affect the lumped equilibrium temperatures
% (T_eq_dir_hemi, T_eq_hemi_hemi, stored in `emissivity`): those are a
% pure radiative energy balance and do not use thermal-conductivity data
% at all, so they remain valid for these sources.
%
% Exclusions are specified by [material_index, source_index] — i.e. the
% position of that source within fourier(material_index).sources (the
% `fourier` and `emissivity` struct arrays share identical material and
% source ordering, so the same indices apply to both). Run this script
% once with print_index_map = true to print every material/source name
% next to its index, then fill in the list below using those indices.
 
print_index_map = false;
 
if print_index_map
    printMaterialSourceIndexMap(fourier);
end
 
% Each row: [material_index, source_index]   % material / author (comment only, for reference)
substrate_adhesive_exclusions = [ ...
      18, 3;   % Ag / Ciesielski-Ge
      18, 4;   % Ag / Ciesielski-Ni
      18, 5;   % Ag / Ciesielski
      21, 2;   % Au / Ciesielski-Ge
      21, 3;   % Au / Ciesielski
];
 
fourier = applySubstrateAdhesiveExclusions(fourier, substrate_adhesive_exclusions);



%% #########################################################################
%%  LOCAL FUNCTIONS
%% #########################################################################

function emissivity_results = computeEmissivityDataset( ...
    solar_flux_dataset, wl_sun, I_sun, multipliers, sigma_SB, ...
    melting_points, thermal_conductivity, T_min, T_max, N_temp)
%COMPUTEEMISSIVITYDATASET
%
% For every material/source in solar_flux_dataset, computes:
%
%   1. Directional total emissivity eps_directional(T,theta,d) (Eq. 41
%      of the Methods) — a material property, Kirchhoff's law applied
%      at each angle, retained in the database but NOT converted into a
%      standalone emitted flux (see file header).
%
%   2. Hemispherical total emissivity eps_hemispherical(T,d), via the
%      standard Lambertian-weighted angular integral:
%         eps_hemi(d) = 2*integral[ eps_lambda(lambda,theta,d)
%                                    * cos(theta)*sin(theta) dtheta ]
%
%   3. TWO lumped-body equilibrium-temperature combinations (both using
%      hemispherical emission, per the file header):
%         I.  directional absorption (cos(theta_i)-projected onto the
%             real surface) -- hemispherical emission
%         II. hemispherical absorption -- hemispherical emission
%
% Each material's temperature grid spans [T_lower, min(T_max, Tm)],
% where T_lower is the lowest temperature in that material's literature
% thermal-conductivity data and Tm is its melting/degradation point (if
% known). Properties are therefore never evaluated above Tm.

wl_sun = wl_sun(:);
I_sun  = I_sun(:);

emissivity_results = struct([]);

for m = 1:numel(solar_flux_dataset)

    mat = solar_flux_dataset(m);
    emissivity_results(m).material = mat.material;
    sources_out = struct([]);

    %% ====================================================================
    % MELTING POINT AND MATERIAL-SPECIFIC TEMPERATURE GRID
    % ====================================================================

    Tm = getMeltingPoint(mat.material, melting_points);

    thermal_data = getThermalConductivity(mat.material, thermal_conductivity);

    if isempty(thermal_data.T)
        warning('%s: thermal conductivity data unavailable. Source skipped.', ...
            mat.material);
        continue
    end

    T_lower = thermal_data.T(1);

    if isfinite(Tm)
        T_upper = min(T_max, Tm);
    else
        % No known melting/degradation point: fall back to the global
        % ceiling. This is a known limitation for this material (no
        % upper validity bound can be enforced).
        T_upper = T_max;
    end

    if T_lower >= T_upper
        warning('%s: invalid temperature range [%g, %g] K. Source skipped.', ...
            mat.material, T_lower, T_upper);
        continue
    end

    T_range = linspace(T_lower, T_upper, N_temp);

    % Make sure the melting point itself is exactly represented, when
    % it falls within range.
    if isfinite(Tm) && Tm <= T_max
        T_range(end) = Tm;
    end

    %% ====================================================================
    % SOURCE LOOP
    % ====================================================================

    for s = 1:numel(mat.sources)

        src = mat.sources(s);
        fprintf('\nEmissivity: %s / %s\n', mat.material, src.name);

        isTdepend = isfield(src,'temperature_K') && numel(src.temperature_K) > 1;

        %% ----------------------------------------------------------------
        % Wavelength grid and matching solar-spectrum samples
        %% ----------------------------------------------------------------

        wl_um = src.wavelength_um(:);
        wl_m  = wl_um * 1e-6;
        Nw    = numel(wl_um);

        [~, idx_sun_common] = ismember(wl_um, wl_sun);
        if any(idx_sun_common == 0)
            error('Some wavelengths in src.wavelength_um are not present in wl_sun.');
        end
        I_sun_common = I_sun(idx_sun_common);

        %% ----------------------------------------------------------------
        % Directional spectral absorptivity, normalized to [Nw x Na x Nt]
        %% ----------------------------------------------------------------

        A_dir = src.A_lambda;

        if ndims(A_dir) == 2
            % Fresnel (opaque, no thickness dependence): [Nw x Na]
            Na = size(A_dir,2);
            Nt = 1;
            A_dir = reshape(A_dir,Nw,Na,1);
        else
            % Slab: [Nw x Na x Nt]
            Na = size(A_dir,2);
            Nt = size(A_dir,3);
        end

        %% ----------------------------------------------------------------
        % Angle grid and Lambertian hemispherical weighting
        %% ----------------------------------------------------------------

        theta = deg2rad(src.angles_deg(:));
        angular_weight = 2*cos(theta).*sin(theta);

        %% ----------------------------------------------------------------
        % TEMPERATURE-DEPENDENT OPTICAL PROPERTIES
        %
        % For T-dependent sources: author temperatures -> T_range, via
        % linear interpolation inside the literature range and linear
        % extrapolation outside it (both directions), followed by
        % clipping to the physical bounds 0<=A,R,T<=1 and enforcement of
        % A+R+T=1 (A treated as primary; R and T share the remainder in
        % their locally interpolated/extrapolated proportion).
        %
        % Because T_range never extends above the melting/degradation
        % point (see above), this extrapolation is never evaluated past
        % that point; it may still be evaluated below the lowest
        % reported literature temperature, which remains a solid-state
        % (same-phase) approximation.
        %% ----------------------------------------------------------------

        T_author = [];

        if isTdepend

            T_author = src.temperature_K(:);

            A_author     = src.A_lambda;
            R_author     = src.R_lambda;
            T_author_opt = src.T_lambda;

            A_dir_T = zeros(Nw,Na,Nt,numel(T_range));
            R_dir_T = zeros(Nw,Na,Nt,numel(T_range));
            T_dir_T = zeros(Nw,Na,Nt,numel(T_range));

            for iw = 1:Nw
                for ia = 1:Na
                    for it = 1:Nt

                        A_values = squeeze(A_author(iw,ia,it,:));
                        R_values = squeeze(R_author(iw,ia,it,:));
                        T_values = squeeze(T_author_opt(iw,ia,it,:));

                        A_values = A_values(:);
                        R_values = R_values(:);
                        T_values = T_values(:);

                        % --- Absorptivity ---
                        A_extrapolated = interp1(T_author, A_values, T_range, 'linear', 'extrap');
                        A_extrapolated = max(A_extrapolated,0);
                        A_extrapolated = min(A_extrapolated,1);

                        A_dir_T(iw,ia,it,:) = reshape(A_extrapolated,1,1,1,[]);

                        % --- Reflectance / transmittance, with A+R+T=1 enforced ---
                        R_raw = interp1(T_author, R_values, T_range, 'linear', 'extrap');
                        T_raw = interp1(T_author, T_values, T_range, 'linear', 'extrap');

                        R_raw = max(min(R_raw,1),0);
                        T_raw = max(min(T_raw,1),0);

                        remaining = 1 - A_extrapolated;

                        RT_sum = R_raw + T_raw;
                        R_fraction = zeros(size(R_raw));
                        T_fraction = zeros(size(T_raw));
                        valid_RT = RT_sum > 0;

                        R_fraction(valid_RT) = R_raw(valid_RT) ./ RT_sum(valid_RT);
                        T_fraction(valid_RT) = T_raw(valid_RT) ./ RT_sum(valid_RT);
                        R_fraction(~valid_RT) = 0.5;
                        T_fraction(~valid_RT) = 0.5;

                        R_extrapolated = min(max(remaining .* R_fraction,0),1);
                        T_extrapolated = min(max(remaining .* T_fraction,0),1);

                        R_dir_T(iw,ia,it,:) = reshape(R_extrapolated,1,1,1,[]);
                        T_dir_T(iw,ia,it,:) = reshape(T_extrapolated,1,1,1,[]);

                    end
                end
            end

        else
            A_dir_T = [];
            R_dir_T = [];
            T_dir_T = [];
        end

        %% ----------------------------------------------------------------
        % Planck spectral radiance B_lambda(lambda,T) [W/m^2/sr/m],
        % evaluated with wavelength in metres (SI), for every T in T_range
        %% ----------------------------------------------------------------

        h_planck = 6.62607015e-34;
        c_light  = 2.99792458e8;
        kB       = 1.380649e-23;

        coeff1 = 2*h_planck*c_light^2;
        coeff2 = h_planck*c_light/kB;

        x = coeff2 ./ (wl_m*T_range);
        x = min(x,500);   % avoid overflow in exp()

        B = coeff1 ./ (wl_m.^5 .* (exp(x)-1));

        %% ----------------------------------------------------------------
        % Thermal spectral coverage ratio: fraction of the full blackbody
        % spectrum captured within this source's available wavelength
        % range, at each temperature. Diagnostic only.
        %% ----------------------------------------------------------------

        thermal_coverage_ratio = zeros(numel(T_range),1);
        for iT = 1:numel(T_range)
            captured = trapz(wl_m, B(:,iT));
            total    = sigma_SB/pi * T_range(iT)^4;
            thermal_coverage_ratio(iT) = captured / total;
        end

        intB = trapz(wl_m,B,1);   % Planck denominator, [1 x Ntemp]

        %% ----------------------------------------------------------------
        % DIRECTIONAL TOTAL EMISSIVITY (Eq. 41) — stored as a material
        % property only; not used to build an emitted flux (see header).
        %% ----------------------------------------------------------------

        eps_dir = zeros(numel(T_range),Na,Nt);

        for it = 1:Nt
            for ia = 1:Na

                if isTdepend
                    alpha_dir_T = squeeze(A_dir_T(:,ia,it,:));
                else
                    alpha_dir_T = A_dir(:,ia,it);
                end

                if isTdepend
                    intEB = zeros(numel(T_range),1);
                    for iT = 1:numel(T_range)
                        intEB(iT) = trapz(wl_m, alpha_dir_T(:,iT).*B(:,iT));
                    end
                else
                    intEB = trapz(wl_m, alpha_dir_T .* B, 1);
                    intEB = intEB(:);
                end

                eps_dir(:,ia,it) = intEB ./ intB(:);

            end
        end

        %% ----------------------------------------------------------------
        % Temperature-dependent solar-flux arrays (only populated for
        % T-dependent sources, where absorptivity itself varies with T).
        %% ----------------------------------------------------------------

        q_absorbed_directional_T    = zeros(numel(T_range), Na, Nt);
        q_reflected_directional_T   = zeros(numel(T_range), Na, Nt);
        q_transmitted_directional_T = zeros(numel(T_range), Na, Nt);

        q_absorbed_hemispherical_T    = zeros(numel(T_range), Nt);
        q_reflected_hemispherical_T   = zeros(numel(T_range), Nt);
        q_transmitted_hemispherical_T = zeros(numel(T_range), Nt);

        if isTdepend

            W = reshape(angular_weight,1,[],1,1);

            A_lambda_hemi_T = reshape(trapz(theta, A_dir_T .* W, 2), Nw, Nt, numel(T_range));
            R_lambda_hemi_T = reshape(trapz(theta, R_dir_T .* W, 2), Nw, Nt, numel(T_range));
            T_lambda_hemi_T = reshape(trapz(theta, T_dir_T .* W, 2), Nw, Nt, numel(T_range));

            for it = 1:Nt
                for iT = 1:numel(T_range)
                    q_absorbed_hemispherical_T(iT,it)    = trapz(wl_um, I_sun_common .* A_lambda_hemi_T(:,it,iT));
                    q_reflected_hemispherical_T(iT,it)   = trapz(wl_um, I_sun_common .* R_lambda_hemi_T(:,it,iT));
                    q_transmitted_hemispherical_T(iT,it) = trapz(wl_um, I_sun_common .* T_lambda_hemi_T(:,it,iT));
                end
            end

            for it = 1:Nt
                for ia = 1:Na
                    for iT = 1:numel(T_range)
                        q_absorbed_directional_T(iT,ia,it)    = trapz(wl_um, I_sun_common .* A_dir_T(:,ia,it,iT));
                        q_reflected_directional_T(iT,ia,it)   = trapz(wl_um, I_sun_common .* R_dir_T(:,ia,it,iT));
                        q_transmitted_directional_T(iT,ia,it) = trapz(wl_um, I_sun_common .* T_dir_T(:,ia,it,iT));
                    end
                end
            end

            eps_lambda_hemi = [];  % not used in the T-dependent branch

        else

            W = reshape(angular_weight,1,[],1);
            eps_lambda_hemi = squeeze(trapz(theta, A_dir .* W, 2));
            eps_lambda_hemi = reshape(eps_lambda_hemi, Nw, Nt);

            A_lambda_hemi_T = [];

        end

        %% ----------------------------------------------------------------
        % PLANCK-WEIGHTED TOTAL HEMISPHERICAL EMISSIVITY
        %% ----------------------------------------------------------------

        eps_hemi = zeros(numel(T_range),Nt);

        for it = 1:Nt
            if isTdepend
                alpha_hemi_T = squeeze(A_lambda_hemi_T(:,it,:));
                for iT = 1:numel(T_range)
                    intEB_hemi = trapz(wl_m, alpha_hemi_T(:,iT).*B(:,iT));
                    eps_hemi(iT,it) = intEB_hemi / intB(iT);
                end
            else
                alpha_hemi_vec = eps_lambda_hemi(:,it);
                intEB_hemi = trapz(wl_m, alpha_hemi_vec .* B, 1);
                eps_hemi(:,it) = intEB_hemi(:) ./ intB(:);
            end
        end

        %% ----------------------------------------------------------------
        % Sanity checks
        %% ----------------------------------------------------------------

        if any(~isfinite(eps_dir(:)) | eps_dir(:) < 0 | eps_dir(:) > 1)
            warning('%s/%s: directional emissivity contains values outside [0,1].', mat.material, src.name);
        end

        if any(~isfinite(eps_hemi(:)) | eps_hemi(:) < 0 | eps_hemi(:) > 1)
            warning('%s/%s: hemispherical emissivity contains values outside [0,1].', mat.material, src.name);
        end

        %% ----------------------------------------------------------------
        % Hemispherical emitted flux [W/m^2]. A standalone "directional
        % emitted flux" is intentionally not computed (see file header).
        %% ----------------------------------------------------------------

        q_emitted_hemispherical = zeros(numel(T_range), Nt);
        for it = 1:Nt
            q_emitted_hemispherical(:,it) = eps_hemi(:,it) .* sigma_SB .* T_range(:).^4;
        end

        %% ----------------------------------------------------------------
        % Solar absorbed flux (angle/thickness-resolved, as computed by
        % compute_solar_flux.m)
        %% ----------------------------------------------------------------

        q_abs_dir  = src.absorbed_flux;              % [Na x Nt]
        q_abs_hemi = reshape(src.absorbed_flux_hemi,1,[]);  % [1 x Nt]

        %% ----------------------------------------------------------------
        % EQUILIBRIUM TEMPERATURE
        %
        % Two combinations only, both balanced against hemispherical
        % emission (see file header):
        %   I.  directional absorption, projected by cos(theta_i) onto
        %       the real surface, vs hemispherical emission
        %   II. hemispherical absorption vs hemispherical emission
        %% ----------------------------------------------------------------

        Nmult = numel(multipliers);

        T_eq_dir_hemi  = NaN(Na,Nt,Nmult);
        T_eq_hemi_hemi = NaN(1,Nt,Nmult);

        for it = 1:Nt

            q_em_hemi = eps_hemi(:,it) .* sigma_SB .* T_range(:).^4;

            if isTdepend
                qa_hemi = q_absorbed_hemispherical_T(:,it);
            else
                qa_hemi = repmat(q_abs_hemi(min(it,numel(q_abs_hemi))), numel(T_range), 1);
            end

            for ia = 1:Na

                if isTdepend
                    qa_dir = q_absorbed_directional_T(:,ia,it);
                else
                    qa_dir = repmat(q_abs_dir(ia, min(it,size(q_abs_dir,2))), numel(T_range), 1);
                end

                % Project the directional (per-beam-area) absorbed flux
                % onto the real, tilted surface so it is on the same
                % per-unit-real-surface-area basis as q_em_hemi.
                qa_dir_proj = qa_dir .* cos(theta(ia));

                for im = 1:Nmult

                    T_eq_dir_hemi(ia,it,im) = findEquilibriumTemperature( ...
                        T_range, q_em_hemi, multipliers(im)*qa_dir_proj);

                end
            end

            for im = 1:Nmult
                T_eq_hemi_hemi(1,it,im) = findEquilibriumTemperature( ...
                    T_range, q_em_hemi, multipliers(im)*qa_hemi);
            end

        end

        %% ----------------------------------------------------------------
        % STORE RESULTS
        %% ----------------------------------------------------------------

        sOut.name            = src.name;
        sOut.model           = src.model;
        sOut.wavelength_um   = src.wavelength_um;
        sOut.angles_deg      = src.angles_deg;
        sOut.thicknesses_um  = src.thicknesses_um;
        sOut.melting_point_K = Tm;

        sOut.eps_directional   = eps_dir;
        sOut.eps_hemispherical = eps_hemi;

        sOut.q_absorbed_directional   = q_abs_dir;
        sOut.q_absorbed_hemispherical = q_abs_hemi;

        if isTdepend
            sOut.T_range              = T_range;
            sOut.temperature_K_author = T_author;

            sOut.A_lambda_T = A_dir_T;
            sOut.R_lambda_T = R_dir_T;
            sOut.T_lambda_T = T_dir_T;

            sOut.A_lambda_hemi_T = A_lambda_hemi_T;
            sOut.R_lambda_hemi_T = R_lambda_hemi_T;
            sOut.T_lambda_hemi_T = T_lambda_hemi_T;

            sOut.q_absorbed_directional_T    = q_absorbed_directional_T;
            sOut.q_reflected_directional_T   = q_reflected_directional_T;
            sOut.q_transmitted_directional_T = q_transmitted_directional_T;

            sOut.q_absorbed_hemispherical_T    = q_absorbed_hemispherical_T;
            sOut.q_reflected_hemispherical_T   = q_reflected_hemispherical_T;
            sOut.q_transmitted_hemispherical_T = q_transmitted_hemispherical_T;
        else
            sOut.T_range              = T_range;
            sOut.temperature_K_author = [];

            sOut.A_lambda_T = [];
            sOut.R_lambda_T = [];
            sOut.T_lambda_T = [];

            sOut.A_lambda_hemi_T = [];
            sOut.R_lambda_hemi_T = [];
            sOut.T_lambda_hemi_T = [];

            sOut.q_absorbed_directional_T    = [];
            sOut.q_reflected_directional_T   = [];
            sOut.q_transmitted_directional_T = [];

            sOut.q_absorbed_hemispherical_T    = [];
            sOut.q_reflected_hemispherical_T   = [];
            sOut.q_transmitted_hemispherical_T = [];
        end

        sOut.q_emitted_hemispherical = q_emitted_hemispherical;

        sOut.multipliers = multipliers;

        sOut.T_eq_dir_hemi  = T_eq_dir_hemi;
        sOut.T_eq_hemi_hemi = T_eq_hemi_hemi;

        sOut.thermal_coverage_ratio = thermal_coverage_ratio;

        if isempty(sources_out)
            sources_out = sOut;
        else
            sources_out(end+1) = sOut;
        end

    end

    emissivity_results(m).sources = sources_out;

end

end


function fourier_results = computeFourierDataset( ...
    emissivity_results, thermal_conductivity, multipliers, sigma_SB)
%COMPUTEFOURIERDATASET
%
% For each material / source / thickness / solar concentration in
% emissivity_results, solves the steady-state, two-surface (front/rear),
% 1-D conduction energy balance for a free-standing, plane-parallel slab:
%
%   front surface: q_abs(Tf) = eps_front(Tf)*sigma*Tf^4 + q_cond
%   rear surface:  q_cond     = eps_rear(Tb)*sigma*Tb^4
%   conduction:    q_cond = (1/d) * integral_{Tb}^{Tf} k(T) dT
%
% Front and rear are assumed to share the same (hemispherical)
% emissivity function, since the underlying optical data characterize
% one homogeneous slab, not two distinct faces.
%
% Only TWO absorption/emission combinations are computed, mirroring the
% lumped-body model and for the same reason (see file header of
% compute_emissivity_equilibrium.m — a standalone "directional emitted
% flux" is not a physically complete quantity):
%
%   I.  directional absorption (cos(theta_i)-projected) -- hemispherical
%       emission from both faces
%   II. hemispherical absorption -- hemispherical emission from both
%       faces
%
% NUMERICAL NOTE: the conductive-flux sensitivity k/d spans many orders
% of magnitude across the thickness grid (sub-micron to cm). The solver
% (solveFourierEquilibrium) addresses this by first finding the
% physically-motivated isothermal-limit temperature (Tf=Tb) via a 1-D
% search, then using it as the initial guess for a scaled 2-D solve; for
% thin/high-conductivity cases where conduction resistance is genuinely
% negligible, the isothermal result is returned directly rather than
% relying on a poorly-conditioned 2-D optimization to resolve a
% physically negligible front/rear split.

fourier_results = struct([]);

for m = 1:numel(emissivity_results)

    material_name = emissivity_results(m).material;
    fourier_results(m).material = material_name;

    thermal_data = getThermalConductivity(material_name, thermal_conductivity);

    sources_in  = emissivity_results(m).sources;
    sources_out = struct([]);

    for s = 1:numel(sources_in)

        src = sources_in(s);
        fprintf('\nFourier: %s / %s\n', material_name, src.name);

        T_range = src.T_range(:);
        thicknesses_um = src.thicknesses_um(:);
        Nt = numel(thicknesses_um);
        Na = numel(src.angles_deg);
        Nmult = numel(multipliers);

        theta = deg2rad(src.angles_deg(:));

        eps_hemi = src.eps_hemispherical;

        q_abs_dir  = src.q_absorbed_directional;
        q_abs_hemi = src.q_absorbed_hemispherical(:);

        isTdepend = ~isempty(src.temperature_K_author);

        if isTdepend
            q_abs_dir_T  = src.q_absorbed_directional_T;
            q_abs_hemi_T = src.q_absorbed_hemispherical_T;
        end

        %% ----------------------------------------------------------------
        % OUTPUT ARRAYS
        %% ----------------------------------------------------------------

        T_front_dir_hemi = NaN(Na,Nt,Nmult);
        T_rear_dir_hemi  = NaN(Na,Nt,Nmult);
        T_mid_dir_hemi   = NaN(Na,Nt,Nmult);
        q_cond_dir_hemi  = NaN(Na,Nt,Nmult);

        T_front_hemi_hemi = NaN(1,Nt,Nmult);
        T_rear_hemi_hemi  = NaN(1,Nt,Nmult);
        T_mid_hemi_hemi   = NaN(1,Nt,Nmult);
        q_cond_hemi_hemi  = NaN(1,Nt,Nmult);

        %% ----------------------------------------------------------------
        % THICKNESS LOOP
        %% ----------------------------------------------------------------

        for it = 1:Nt

            d_um = thicknesses_um(it);

            epsH = eps_hemi(:,it);
            eps_hemi_fun = @(T) interp1(T_range, epsH, T, 'linear', 'extrap');

            if isTdepend
                qa_hemi_all = q_abs_hemi_T(:,it);
                qa_hemi_fun = @(T) interp1(T_range, qa_hemi_all, T, 'linear', 'extrap');
            else
                qa_hemi_fun = @(T) repmat(q_abs_hemi(it), size(T));
            end

            %% ------------------------------------------------------------
            % I. DIRECTIONAL ABSORPTION (projected) -- HEMISPHERICAL EMISSION
            %% ------------------------------------------------------------

            for ia = 1:Na

                if isTdepend
                    qa_dir_all = q_abs_dir_T(:,ia,it);
                    qa_dir_proj = qa_dir_all .* cos(theta(ia));
                    qa_dir_fun = @(T) interp1(T_range, qa_dir_proj, T, 'linear', 'extrap');
                else
                    qa_dir_proj = q_abs_dir(ia,it) .* cos(theta(ia));
                    qa_dir_fun = @(T) repmat(qa_dir_proj, size(T));
                end

                for im = 1:Nmult

                    C = multipliers(im);
                    q_abs_C_fun = @(T) C .* qa_dir_fun(T);

                    [Tf,Tb,Tc,qc] = solveFourierEquilibrium( ...
                        d_um, q_abs_C_fun, eps_hemi_fun, eps_hemi_fun, ...
                        sigma_SB, thermal_data, T_range);

                    T_front_dir_hemi(ia,it,im) = Tf;
                    T_rear_dir_hemi(ia,it,im)  = Tb;
                    T_mid_dir_hemi(ia,it,im)   = Tc;
                    q_cond_dir_hemi(ia,it,im)  = qc;

                end
            end

            %% ------------------------------------------------------------
            % II. HEMISPHERICAL ABSORPTION -- HEMISPHERICAL EMISSION
            %% ------------------------------------------------------------

            for im = 1:Nmult

                C = multipliers(im);
                q_abs_C_fun = @(T) C .* qa_hemi_fun(T);

                [Tf,Tb,Tc,qc] = solveFourierEquilibrium( ...
                    d_um, q_abs_C_fun, eps_hemi_fun, eps_hemi_fun, ...
                    sigma_SB, thermal_data, T_range);

                T_front_hemi_hemi(1,it,im) = Tf;
                T_rear_hemi_hemi(1,it,im)  = Tb;
                T_mid_hemi_hemi(1,it,im)   = Tc;
                q_cond_hemi_hemi(1,it,im)  = qc;

            end

        end

        %% ----------------------------------------------------------------
        % STORE RESULTS
        %% ----------------------------------------------------------------

        sOut.name            = src.name;
        sOut.model           = src.model;
        sOut.wavelength_um   = src.wavelength_um;
        sOut.angles_deg      = src.angles_deg;
        sOut.thicknesses_um  = src.thicknesses_um;
        sOut.T_range         = src.T_range;
        sOut.melting_point_K = src.melting_point_K;

        sOut.T_thermal_conductivity   = thermal_data.T;
        sOut.k_thermal_conductivity   = thermal_data.k;
        sOut.ref_thermal_conductivity = thermal_data.ref;

        sOut.T_front_dir_hemi = T_front_dir_hemi;
        sOut.T_rear_dir_hemi  = T_rear_dir_hemi;
        sOut.T_mid_dir_hemi   = T_mid_dir_hemi;
        sOut.q_cond_dir_hemi  = q_cond_dir_hemi;

        sOut.T_front_hemi_hemi = T_front_hemi_hemi;
        sOut.T_rear_hemi_hemi  = T_rear_hemi_hemi;
        sOut.T_mid_hemi_hemi   = T_mid_hemi_hemi;
        sOut.q_cond_hemi_hemi  = q_cond_hemi_hemi;

        sOut.multipliers = multipliers;

        if isempty(sources_out)
            sources_out = sOut;
        else
            sources_out(end+1) = sOut;
        end

    end

    fourier_results(m).sources = sources_out;

end

end


function [Tf, Tb, Tc, q_cond] = solveFourierEquilibrium( ...
    d_um, q_abs_fun, eps_front_fun, eps_rear_fun, sigma_SB, thermal_data, T_range)
%SOLVEFOURIEREQUILIBRIUM
%
% See the header of computeFourierDataset for the governing equations.
%
% Solves the 2-D front/rear system on NON-DIMENSIONALIZED variables
% (temperature scaled by the isothermal-limit temperature T0, residuals
% scaled by sigma*T0^4), seeded from T0 itself, because the raw
% conductive sensitivity k(T)/d varies by many orders of magnitude
% across the thickness grid and would otherwise leave an unscaled
% solver started from an arbitrary fixed guess poorly conditioned,
% especially for thin films. If the 2-D solve fails to converge, or
% resolves no meaningful front/rear split, the isothermal-limit result
% is returned directly — which is also the physically correct answer
% whenever conduction resistance across the slab is negligible.

Tf = NaN; Tb = NaN; Tc = NaN; q_cond = NaN;

if isempty(thermal_data.T) || isempty(thermal_data.k)
    return
end

T_low  = max(T_range(1),   thermal_data.T(1));
T_high = min(T_range(end), thermal_data.T(end));

if ~isfinite(T_low) || ~isfinite(T_high) || T_low >= T_high
    return
end

d = d_um * 1e-6;
if d <= 0 || ~isfinite(d)
    return
end

%% ------------------------------------------------------------
% Isothermal-limit initial guess (also the fallback result)
%% ------------------------------------------------------------

T0 = estimateIsothermalTemperature( ...
    T_range, q_abs_fun, eps_front_fun, eps_rear_fun, sigma_SB, T_low, T_high);

if ~isfinite(T0)
    return
end

%% ------------------------------------------------------------
% Scaled 2-D solve
%% ------------------------------------------------------------

Tscale = max(T0, 1);
Fscale = sigma_SB * Tscale^4;

y0   = [T0, T0] ./ Tscale;
lb_T = [T_low, T_low];
ub_T = [T_high, T_high];

options = optimoptions( ...
    'fsolve', ...
    'Display', 'off', ...
    'FunctionTolerance', 1e-12, ...
    'StepTolerance', 1e-12, ...
    'MaxIterations', 500, ...
    'MaxFunctionEvaluations', 2000);

scaledResiduals = @(y) fourierResiduals( ...
    min(max(y .* Tscale, lb_T), ub_T), ...
    d, q_abs_fun, eps_front_fun, eps_rear_fun, sigma_SB, thermal_data) ./ Fscale;

Tf_candidate = NaN;
Tb_candidate = NaN;

try
    [y_sol, ~, exitflag] = fsolve(scaledResiduals, y0, options);
    if exitflag > 0
        x_sol = min(max(y_sol .* Tscale, lb_T), ub_T);
        Tf_candidate = x_sol(1);
        Tb_candidate = x_sol(2);
    end
catch
    % fall through to the isothermal fallback below
end

if ~isfinite(Tf_candidate) || ~isfinite(Tb_candidate) || Tf_candidate < Tb_candidate
    % 2-D solve failed, or did not resolve a physical (front >= rear)
    % split: fall back to the isothermal-limit result.
    Tf_candidate = T0;
    Tb_candidate = T0;
end

%% ------------------------------------------------------------
% Conductive flux and mid-plane temperature
%% ------------------------------------------------------------

if Tf_candidate == Tb_candidate

    q_cond_candidate = 0;
    Tc_candidate = Tf_candidate;

else

    k_integral = integral( ...
        @(T) evaluateThermalConductivity(T, thermal_data), ...
        Tb_candidate, Tf_candidate, 'ArrayValued', true);

    q_cond_candidate = k_integral / d;

    target_integral = 0.5 * k_integral;

    try
        Tc_candidate = fzero( ...
            @(T) integral( ...
                    @(TT) evaluateThermalConductivity(TT, thermal_data), ...
                    Tb_candidate, T, 'ArrayValued', true) ...
                 - target_integral, ...
            [Tb_candidate, Tf_candidate]);
    catch
        Tc_candidate = 0.5*(Tf_candidate + Tb_candidate);
    end

end

Tf = Tf_candidate;
Tb = Tb_candidate;
Tc = Tc_candidate;
q_cond = q_cond_candidate;

end


function T0 = estimateIsothermalTemperature( ...
    T_range, q_abs_fun, eps_front_fun, eps_rear_fun, sigma_SB, T_low, T_high)
%ESTIMATEISOTHERMALTEMPERATURE
%
% 1-D sign-change search (reusing findEquilibriumTemperature) for the
% isothermal (zero conductive-gradient) balance:
%
%   q_abs(T) = [eps_front(T) + eps_rear(T)] * sigma * T^4
%
% This is the physically correct limit whenever conduction resistance
% across the slab is negligible, and is used both as the initial guess
% for, and — when that does not converge to a physically meaningful
% split — the reported result of, the full 2-D front/rear solve.

T_scan = T_range(T_range >= T_low & T_range <= T_high);

if numel(T_scan) < 2
    T0 = NaN;
    return
end

T_scan = T_scan(:);
q_em  = (eps_front_fun(T_scan) + eps_rear_fun(T_scan)) .* sigma_SB .* T_scan.^4;
q_abs = q_abs_fun(T_scan);

T0 = findEquilibriumTemperature(T_scan, q_em, q_abs(:));

end


function F = fourierResiduals( ...
    x, d, q_abs_fun, eps_front_fun, eps_rear_fun, sigma_SB, thermal_data)
%FOURIERRESIDUALS  Physical-units residuals of the front/rear energy
%balance, evaluated at x = [Tf, Tb] (both in Kelvin).

Tf = x(1);
Tb = x(2);

try
    k_integral = integral( ...
        @(T) evaluateThermalConductivity(T, thermal_data), ...
        Tb, Tf, 'ArrayValued', true);
catch
    F = [1e20; 1e20];
    return
end

if ~isfinite(k_integral)
    F = [1e20; 1e20];
    return
end

q_cond = k_integral / d;

q_abs = q_abs_fun(Tf);
eps_front = eps_front_fun(Tf);
eps_rear  = eps_rear_fun(Tb);

if ~isfinite(q_abs) || ~isfinite(eps_front) || ~isfinite(eps_rear)
    F = [1e20; 1e20];
    return
end

F_rear  = q_cond - eps_rear  .* sigma_SB .* Tb.^4;
F_front = q_abs  - eps_front .* sigma_SB .* Tf.^4 - q_cond;

F = [F_front; F_rear];

end


function T_eq = findEquilibriumTemperature(T_range, q_emitted, target)
%FINDEQUILIBRIUMTEMPERATURE
%
% Finds the root of q_emitted(T) - target(T) = 0 by scanning T_range
% for a sign change and linearly interpolating between the two
% bracketing points. Returns NaN if no sign change is found in range
% (i.e. no equilibrium exists within [T_range(1), T_range(end)]).

T_range   = T_range(:);
q_emitted = q_emitted(:);

if isscalar(target)
    target_flux = repmat(target, numel(T_range), 1);
else
    target_flux = target(:);
end

diff_q = q_emitted - target_flux;

idx = find(diff_q == 0, 1, 'first');
if ~isempty(idx)
    T_eq = T_range(idx);
    return
end

idx = find(diff_q(1:end-1).*diff_q(2:end) < 0, 1, 'first');

if ~isempty(idx)
    T_eq = interp1(diff_q(idx:idx+1), T_range(idx:idx+1), 0);
else
    T_eq = NaN;
end

end


function normalized_name = normalizeMaterialName(material_name)
%NORMALIZEMATERIALNAME
%
% Strips common descriptive suffixes (film/thickness, crystal phase,
% polytype) from a material or struct-field name and reduces it to a
% lowercase, alphanumeric-only string, for matching database material
% names against melting_points / thermal_conductivity struct fields.

base_name = material_name;

base_name = regexprep(base_name, '\s+(Thick|Thin)\s+Film$', '', 'ignorecase');
base_name = regexprep(base_name, '\s+(Amorphous|amorphous|Crystal|Fused|Nanoparticles)$', '', 'ignorecase');
base_name = regexprep(base_name, '\s+(3C|4H|6H)$', '', 'ignorecase');

normalized_name = lower(regexprep(base_name, '[^a-zA-Z0-9]', ''));

end


function Tm = getMeltingPoint(material_name, melting_points)
%GETMELTINGPOINT  Looks up the melting/degradation point [K] for
%material_name by normalized-name matching against melting_points'
%field names. Returns NaN if no match is found.

Tm = NaN;
normalized_name = normalizeMaterialName(material_name);

fields = fieldnames(melting_points);
for i = 1:numel(fields)
    if strcmp(normalized_name, normalizeMaterialName(fields{i}))
        Tm = melting_points.(fields{i});
        return
    end
end

end


function thermal_data = getThermalConductivity(material_name, thermal_conductivity)
%GETTHERMALCONDUCTIVITY  Looks up and validates the literature
%thermal-conductivity data (T, k, ref) for material_name by
%normalized-name matching against thermal_conductivity's field names.
%Returns empty fields if no match is found or the data fail validation.

thermal_data = struct('T', [], 'k', [], 'ref', []);

normalized_name = normalizeMaterialName(material_name);

matched_field = '';
fields = fieldnames(thermal_conductivity);
for i = 1:numel(fields)
    if strcmp(normalized_name, normalizeMaterialName(fields{i}))
        matched_field = fields{i};
        break
    end
end

if isempty(matched_field)
    return
end

data = thermal_conductivity.(matched_field);

if ~isfield(data,'T') || ~isfield(data,'k') || ~isfield(data,'ref')
    return
end

T = data.T(:);
k = data.k(:);
ref = data.ref(:);

if isempty(T) || isempty(k) || isempty(ref)
    return
end

if numel(T) ~= numel(k)
    error('Thermal conductivity data for %s have inconsistent dimensions.', material_name);
end

valid = isfinite(T) & isfinite(k);
T = T(valid);
k = k(valid);

if any(T <= 0)
    error('Thermal conductivity data for %s contain non-positive temperatures.', material_name);
end

if any(k <= 0)
    error('Thermal conductivity data for %s contain non-positive values.', material_name);
end

[T, sort_idx] = sort(T);
k = k(sort_idx);

if numel(unique(T)) ~= numel(T)
    error('Thermal conductivity data for %s contain duplicate temperatures.', material_name);
end

thermal_data.T = T;
thermal_data.k = k;
thermal_data.ref = ref;

end


function k = evaluateThermalConductivity(T, thermal_data)
%EVALUATETHERMALCONDUCTIVITY  Interpolates k(T) from literature data.
%No extrapolation: temperatures outside the literature range return NaN.

k = NaN(size(T));

if isempty(thermal_data.T) || isempty(thermal_data.k)
    return
end

T_data = thermal_data.T;
k_data = thermal_data.k;

valid = (T >= T_data(1)) & (T <= T_data(end));

k(valid) = interp1(T_data, k_data, T(valid), 'linear');

end

function printMaterialSourceIndexMap(fourier_results)
%PRINTMATERIALSOURCEINDEXMAP
%
% Prints every [material_index, source_index] pair in fourier_results
% next to its material and source name, so that entries for
% substrate_adhesive_exclusions can be looked up by name and then
% recorded by index.
 
fprintf('\n============================================================\n');
fprintf('Material / source index map (for substrate_adhesive_exclusions)\n');
fprintf('============================================================\n');
 
for m_idx = 1:numel(fourier_results)
 
    mat_name = fourier_results(m_idx).material;
    srcs = fourier_results(m_idx).sources;
 
    for s_idx = 1:numel(srcs)
 
        fprintf('  [%2d, %2d]   %-10s   %s\n', ...
            m_idx, s_idx, mat_name, srcs(s_idx).name);
 
    end
 
end
 
fprintf('============================================================\n\n');
 
end
 
 
function fourier_results = applySubstrateAdhesiveExclusions(fourier_results, exclusions)
%APPLYSUBSTRATEADHESIVEEXCLUSIONS
%
% Overwrites the Fourier-model outputs (conductive flux and front/rear/
% mid-plane temperatures, for both the dir_hemi and hemi_hemi
% combinations) with NaN for the [material_index, source_index] pairs
% listed in `exclusions`, because the thermal-conductivity data used by
% the Fourier model refers to the pure material, not to the substrate-
% or adhesive-containing stack actually measured by these sources.
%
% `exclusions` is an [N x 2] array of [material_index, source_index]
% rows, referring to positions within fourier_results(material_index).
% sources. Out-of-range indices are skipped with a warning rather than
% causing an error, so a typo in the list does not halt the run.
 
fourier_fields = { ...
    'T_front_dir_hemi', 'T_rear_dir_hemi', 'T_mid_dir_hemi', 'q_cond_dir_hemi', ...
    'T_front_hemi_hemi','T_rear_hemi_hemi','T_mid_hemi_hemi','q_cond_hemi_hemi'};
 
if isempty(exclusions)
    return
end
 
for r = 1:size(exclusions,1)
 
    m_idx = exclusions(r,1);
    s_idx = exclusions(r,2);
 
    if m_idx < 1 || m_idx > numel(fourier_results)
        warning('Substrate/adhesive exclusion row %d: material index %d out of range.', ...
            r, m_idx);
        continue
    end
 
    if s_idx < 1 || s_idx > numel(fourier_results(m_idx).sources)
        warning('Substrate/adhesive exclusion row %d: source index %d out of range for material %s.', ...
            r, s_idx, fourier_results(m_idx).material);
        continue
    end
 
    src = fourier_results(m_idx).sources(s_idx);
 
    fprintf('Excluding Fourier results (substrate/adhesive system): %s / %s [material %d, source %d]\n', ...
        fourier_results(m_idx).material, src.name, m_idx, s_idx);
 
    for f = 1:numel(fourier_fields)
        fname = fourier_fields{f};
        if isfield(src, fname)
            src.(fname) = NaN(size(src.(fname)));
        end
    end
 
    fourier_results(m_idx).sources(s_idx) = src;
 
end
 
end
