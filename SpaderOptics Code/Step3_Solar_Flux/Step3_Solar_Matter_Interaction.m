%% Step3_Solar_Matter_Interaction.m
%  =========================================================================
%  Computes — for every material/source imported by the Step 1/2 import
%  scripts — the optical response (reflectivity, transmissivity,
%  absorptivity) and the TOTAL solar flux (W/m^2) absorbed, reflected, and
%  transmitted, at each incidence angle and slab thickness, given the AM0
%  E-490 solar spectrum.
%
%  PREREQUISITE
%  ------------
%  The struct array `materials` (and/or `Inorganic_materials`, which is
%  merged into it below) produced by the import scripts must already be
%  present in the base workspace, or loadable from materials.mat.
%
%  OPTICAL MODEL
%  --------------
%  Every material/source is treated with the same plane-parallel-slab,
%  incoherent-multiple-internal-reflection model (see calcSlabOptics
%  below), for every angle in angles_deg and every thickness in
%  thicknesses_um. There is no separate opaque/Fresnel-only special case:
%  the slab formulas reduce to the single-interface Fresnel result on
%  their own in the opaque (large absorption, large thickness) limit, so
%  applying them uniformly is both simpler and correct.
%
%  SOLAR SPECTRUM
%  --------------
%  The AM0 E-490 (ASTM E490) standard spectrum is loaded via LoadSolar()
%  (a separate utility function/script that must be on the MATLAB path;
%  not included in this file). Only wavelength (um) and spectral
%  irradiance (W/m^2/um) are needed.
%
%  METHOD
%  ------
%  For each material, source, incidence angle theta and slab thickness d:
%
%     Spectral absorbed irradiance  P_abs(lambda) = I_sun(lambda).*A(lambda,theta,d)   [W/m^2/um]
%     Spectral reflected irradiance P_ref(lambda) = I_sun(lambda).*R(lambda,theta,d)   [W/m^2/um]
%     Spectral transmitted          P_tr (lambda) = I_sun(lambda).*T(lambda,theta,d)   [W/m^2/um]
%
%  The TOTAL fluxes are the wavelength integrals (trapezoidal rule):
%
%     F_abs(theta,d) = integral P_abs(lambda) dlambda    [W/m^2]
%     F_ref(theta,d) = integral P_ref(lambda) dlambda    [W/m^2]
%     F_tr (theta,d) = integral P_tr (lambda) dlambda    [W/m^2]
%
%  Integration is performed over the OVERLAP of the material's wavelength
%  grid and the solar spectrum range — outside that overlap there is no
%  optical data, so the contribution is zero. The incident solar flux
%  used for normalisation (F_sun_total) is integrated over the same
%  overlap range so that F_abs + F_ref + F_tr = F_sun_total holds to
%  numerical accuracy (checked explicitly below).
%
%  Hemispherical quantities (R, T, A, and the corresponding fluxes) are
%  obtained from the directional ones via the standard Lambertian
%  angular-averaging integral:
%
%     X_hemispherical = 2 * integral[ X_directional(theta) * cos(theta) * sin(theta) dtheta ]
%
%  OUTPUT
%  ------
%  struct array  solar_flux(m)  with fields:
%     .material
%     .sources(s) :
%        .name, .source_url, .temperature, .temperature_K, .film_thickness_um
%        .model                'slab' (always — see OPTICAL MODEL above)
%        .angles_deg           [1 x Na]
%        .thicknesses_um       [1 x Nt]
%        .wavelength_range_um  [lambda_min lambda_max]  overlap actually integrated
%        .wavelength_um        [Nw x 1]  wavelength grid actually used
%        .F_sun_total          scalar — total incident solar flux (W/m^2) over the overlap
%        .coverage_ratio       scalar — fraction of the full AM0 spectrum covered by this overlap
%        .n_complex_resampled  [Nw x 1] or [Nw x NT]  n+ik resampled onto wavelength_um
%        .R_lambda, .T_lambda, .A_lambda       directional spectral optical properties
%        .R_hemi_lambda, .T_hemi_lambda, .A_hemi_lambda   hemispherical spectral optical properties
%        .absorbed_flux, .reflected_flux, .transmitted_flux             directional fluxes [W/m^2]
%        .absorbed_flux_hemi, .reflected_flux_hemi, .transmitted_flux_hemi   hemispherical fluxes [W/m^2]
%
%  No coverage-based filtering or "best source per material" selection is
%  performed here: every source that has at least some wavelength overlap
%  with the AM0 spectrum is kept, and its coverage_ratio is stored so a
%  user can filter downstream if desired.
%  =========================================================================

clc;

%% ==== USER INPUTS =======================================================
angleStep = 5;                % Angular resolution [deg]

%% ==== ANGLE GRID =========================================================
angles_deg = [0:angleStep:(90-angleStep), 89.9];
angles_rad = deg2rad(angles_deg);

%% ==== MERGE MATERIAL DATABASES ===========================================
%
% Step 1 (Inorganic_materials) and Step 2 (materials) are merged here into
% a single `materials` struct array before computing optical/solar
% quantities.

if exist('Inorganic_materials','var') && ~isempty(Inorganic_materials)

    if ~exist('materials','var') || isempty(materials)
        materials = Inorganic_materials;
    else
        end_mats = length(materials);
        for i = 1:length(Inorganic_materials)
            materials(end_mats+i) = Inorganic_materials(i);
        end
    end

end

%% ==== LOAD MATERIAL DATABASE =============================================
if ~exist('materials','var') || isempty(materials)

    if exist('materials.mat','file')
        fprintf('Loading materials.mat ...\n');
        S = load('materials.mat','materials');
        materials = S.materials;
    else
        error(['materials struct not found. Run the Step 1/2 import ', ...
               'scripts first.']);
    end

end

%% ==== LOAD SOLAR SPECTRUM ================================================

[I_sun,wl_sun,wl_sun_min,wl_sun_max,AM0_total] = LoadSolar();

fprintf('\nSolar spectrum loaded:\n');
fprintf('   lambda = %.4f - %.4f um\n', wl_sun_min, wl_sun_max);
fprintf('   Integrated AM0 flux = %.2f W/m^2\n\n', AM0_total);

%% ==== OUTPUT STRUCTURES ==================================================

solar_flux = struct([]);
noOverlap  = struct([]);

iKeep      = 0;
iNoOverlap = 0;

%% ==== COMPUTE FLUXES PER MATERIAL / SOURCE ==============================

for m = 1:numel(materials)

    mat = materials(m);
    sources_all = struct([]);

    for s = 1:numel(mat.sources)

        src = mat.sources(s);

        %% ===============================================================
        %  MATERIAL WAVELENGTH OVERLAP WITH SOLAR SPECTRUM
        %% ===============================================================

        wl_mat = src.wavelength(:);

        % Remove duplicated wavelengths
        [wl_mat, ia_u] = unique(wl_mat,'stable');

        % Detect temperature-dependent source
        isTdepend = isfield(src,'temperature_K') && ...
                    numel(src.temperature_K) > 1 && ...
                    size(src.n_complex,2) == numel(src.temperature_K);

        if isTdepend
            temperature_K = src.temperature_K(:).';
            n_raw = real(src.n_complex(ia_u,:));   % wavelength x temperature
            k_raw = imag(src.n_complex(ia_u,:));
        else
            n_raw = real(src.n_complex(ia_u));
            k_raw = imag(src.n_complex(ia_u));
        end

        % Solar overlap range
        loLim = max(min(wl_mat),min(wl_sun));
        hiLim = min(max(wl_mat),max(wl_sun));

        if hiLim <= loLim

            fprintf('NO OVERLAP : %-15s | %s\n', mat.material, src.name);

            iNoOverlap = iNoOverlap + 1;
            noOverlap(iNoOverlap).material = mat.material;
            noOverlap(iNoOverlap).source   = src.name;

            continue

        end

        %% ===============================================================
        %  SOLAR GRID USED FOR INTEGRATION
        %% ===============================================================

        solarMask = wl_sun>=loLim & wl_sun<=hiLim;
        wl_int = wl_sun(solarMask);
        I_int  = I_sun(solarMask);

        %% ===============================================================
        %  SOLAR COVERAGE (diagnostic metric; not used for filtering)
        %% ===============================================================

        F_sun_covered = trapz(wl_int,I_int);
        coverageRatio = F_sun_covered/AM0_total;
        F_sun_total   = F_sun_covered;

        %% ===============================================================
        %  RESAMPLE n AND k (NOT R,T,A) ONTO THE SOLAR WAVELENGTH GRID
        %% ===============================================================

        n_int = interp1(wl_mat, n_raw, wl_int, 'linear', 'extrap');
        k_int = interp1(wl_mat, k_raw, wl_int, 'linear', 0);

        N_int = complex(n_int,k_int);

        %% ===============================================================
        %  COMPUTE OPTICAL PROPERTIES (plane-parallel slab, all angles
        %  and thicknesses; see OPTICAL MODEL in the file header)
        %% ===============================================================

        thk_vec = src.thicknesses_um;

        if isTdepend

            NT = numel(temperature_K);
            Nt = numel(thk_vec);
            Nw = numel(wl_int);
            Na = numel(angles_rad);

            % Output dimensions: wavelength x angle x thickness x temperature
            R = zeros(Nw,Na,Nt,NT);
            T = zeros(Nw,Na,Nt,NT);
            A = zeros(Nw,Na,Nt,NT);

            for iT = 1:NT
                N_T = N_int(:,iT);
                [R(:,:,:,iT), T(:,:,:,iT), A(:,:,:,iT)] = ...
                    calcSlabOptics(N_T, wl_int, angles_rad, thk_vec);
            end

        else

            [R,T,A] = calcSlabOptics(N_int, wl_int, angles_rad, thk_vec);

        end

        %% ===============================================================
        % SANITY CHECKS
        %% ===============================================================

        if any(R(:)<0 | R(:)>1)
            warning('%s/%s: invalid reflectance corrected.', mat.material, src.name);
            R = max(min(R,1),0);
        end

        if any(T(:)<0 | T(:)>1)
            warning('%s/%s: invalid transmissivity corrected.', mat.material, src.name);
            T = max(min(T,1),0);
        end

        if any(A(:)<0 | A(:)>1)
            warning('%s/%s: invalid absorptivity corrected.', mat.material, src.name);
            A = max(min(A,1),0);
        end

        energy_error = max(abs(R+T+A-1),[],'all');
        if energy_error > 1e-6
            warning('%s/%s: energy conservation error %.3e.', mat.material, src.name, energy_error);
        end

        %% ===============================================================
        %  HEMISPHERICAL OPTICAL PROPERTIES
        %
        %  X_hemispherical(lambda) = 2 * integral[ X(lambda,theta) * cos(theta)*sin(theta) dtheta ]
        %% ===============================================================

        theta = angles_rad(:);
        angular_weight = 2*cos(theta).*sin(theta);

        if isTdepend
            W = reshape(angular_weight,1,[],1,1);
        else
            W = reshape(angular_weight,1,[],1);
        end

        R_hemi_lambda = squeeze(trapz(theta, R .* W, 2));
        T_hemi_lambda = squeeze(trapz(theta, T .* W, 2));
        A_hemi_lambda = squeeze(trapz(theta, A .* W, 2));

        %% ===============================================================
        %  DIRECTIONAL SOLAR FLUX INTEGRATION
        %% ===============================================================

        Na = size(R,2);
        Nt = size(R,3);

        if isTdepend
            NT = size(R,4);
            I_bc = reshape(I_int,[],1,1,1);
        else
            I_bc = reshape(I_int,[],1,1);
        end

        spec_abs = I_bc .* A;
        spec_ref = I_bc .* R;
        spec_tr  = I_bc .* T;

        F_abs = squeeze(trapz(wl_int,spec_abs,1));
        F_ref = squeeze(trapz(wl_int,spec_ref,1));
        F_tr  = squeeze(trapz(wl_int,spec_tr,1));

        if isTdepend
            % angle x thickness x temperature
            F_abs = reshape(F_abs,Na,Nt,NT);
            F_ref = reshape(F_ref,Na,Nt,NT);
            F_tr  = reshape(F_tr ,Na,Nt,NT);
        else
            % angle x thickness
            F_abs = reshape(F_abs,Na,Nt);
            F_ref = reshape(F_ref,Na,Nt);
            F_tr  = reshape(F_tr ,Na,Nt);
        end

        %% ===============================================================
        %  ENERGY CONSERVATION CHECK (directional fluxes)
        %% ===============================================================

        energy_flux_error = max(abs(F_abs + F_ref + F_tr - F_sun_total),[],'all');

        if energy_flux_error > 1e-3
            warning('%s/%s: solar energy balance error %.3e W/m^2', ...
                mat.material, src.name, energy_flux_error);
        end

        %% ===============================================================
        %  HEMISPHERICAL SOLAR FLUX
        %% ===============================================================

        if isTdepend

            % R/T/A_hemi_lambda: wavelength x thickness x temperature
            F_ref_hemi = trapz(wl_int, reshape(I_int,[],1,1) .* R_hemi_lambda, 1);
            F_tr_hemi  = trapz(wl_int, reshape(I_int,[],1,1) .* T_hemi_lambda, 1);
            F_abs_hemi = trapz(wl_int, reshape(I_int,[],1,1) .* A_hemi_lambda, 1);

            F_ref_hemi = reshape(F_ref_hemi,Nt,NT);
            F_tr_hemi  = reshape(F_tr_hemi ,Nt,NT);
            F_abs_hemi = reshape(F_abs_hemi,Nt,NT);

        else

            F_ref_hemi = trapz(wl_int, I_int .* R_hemi_lambda, 1);
            F_tr_hemi  = trapz(wl_int, I_int .* T_hemi_lambda, 1);
            F_abs_hemi = trapz(wl_int, I_int .* A_hemi_lambda, 1);

        end

        %% ===============================================================
        %  STORE RESULTS
        %% ===============================================================

        sOut.name = src.name;

        if isfield(src,'source_url')
            sOut.source_url = src.source_url;
        else
            sOut.source_url = '';
        end

        if isfield(src,'temperature')
            sOut.temperature = src.temperature;
        else
            sOut.temperature = '';
        end

        if isfield(src,'temperature_K')
            sOut.temperature_K = src.temperature_K;
        else
            sOut.temperature_K = NaN;
        end

        if isfield(src,'film_thickness_um')
            sOut.film_thickness_um = src.film_thickness_um;
        else
            sOut.film_thickness_um = NaN;
        end

        sOut.model = 'slab';

        sOut.wavelength_range_um = [loLim hiLim];
        sOut.coverage_ratio      = coverageRatio;
        sOut.wavelength_um       = wl_int;
        sOut.F_sun_total         = F_sun_total;
        sOut.n_complex_resampled = N_int;

        sOut.angles_deg     = angles_deg;
        sOut.thicknesses_um = thk_vec;

        sOut.R_lambda = R;
        sOut.T_lambda = T;
        sOut.A_lambda = A;

        sOut.absorbed_flux    = F_abs;
        sOut.reflected_flux   = F_ref;
        sOut.transmitted_flux = F_tr;

        sOut.R_hemi_lambda = R_hemi_lambda;
        sOut.T_hemi_lambda = T_hemi_lambda;
        sOut.A_hemi_lambda = A_hemi_lambda;

        sOut.absorbed_flux_hemi    = F_abs_hemi;
        sOut.reflected_flux_hemi   = F_ref_hemi;
        sOut.transmitted_flux_hemi = F_tr_hemi;

        if isempty(sources_all)
            sources_all = sOut;
        else
            sources_all(end+1) = sOut;
        end

    end

    if ~isempty(sources_all)
        iKeep = iKeep + 1;
        solar_flux(iKeep).material = mat.material;
        solar_flux(iKeep).sources  = sources_all;
    end

end

%% ==== SUMMARY ===========================================================
fprintf('\n========== Solar flux computation complete ==========\n');
for m = 1:numel(solar_flux)
    nSrc = numel(solar_flux(m).sources);
    fprintf('  %-12s : %d source(s)\n', solar_flux(m).material, nSrc);
    for s = 1:nSrc
        sf = solar_flux(m).sources(s);

        iAng = 1;
        iThk = 1;

        fprintf('      %-25s [%s] lambda=[%.3f-%.3f] um  Fsun=%.2f W/m^2\n', ...
            sf.name, sf.model, sf.wavelength_range_um(1), sf.wavelength_range_um(2), sf.F_sun_total);

        fprintf(['      %25s theta=%.0f deg, d=%.3g um : ', ...
                 'Abs=%.2f  Ref=%.2f  Trans=%.2f W/m^2\n'], ...
            '', sf.angles_deg(iAng), sf.thicknesses_um(iThk), ...
            sf.absorbed_flux(iAng,iThk), sf.reflected_flux(iAng,iThk), sf.transmitted_flux(iAng,iThk));

    end
end

fprintf('\n=============================================\n');
fprintf('Sources excluded due to zero solar overlap:\n');
fprintf('=============================================\n');
for i = 1:numel(noOverlap)
    fprintf('%-15s   %s\n', noOverlap(i).material, noOverlap(i).source);
end


%% ============ LOCAL FUNCTIONS  =======================

function [R_slab,T_slab,A_slab] = calcSlabOptics( ...
        n_complex, wavelength_um, angles_rad, thicknesses_um)
%CALCSLABOPTICS
%
% Incoherent multiple-internal-reflection model for a plane-parallel slab
% surrounded by vacuum. Calculates R(lambda,theta,d), T(lambda,theta,d),
% A(lambda,theta,d).
%
% For each polarisation (s,p), with single-interface Fresnel reflectance
% R_if, single-pass intensity attenuation tau = exp(-2*Im(kz)*d), and the
% round-trip (two-pass) attenuation tau2 = tau^2:
%
%   T_pol     = (1-R_if)^2 * tau   / (1 - R_if^2 * tau2)
%   R_pol_tot = R_if + R_if*(1-R_if)^2 * tau2 / (1 - R_if^2 * tau2)
%
% (tau2, not tau, appears in the denominator and in the internal-
% reflection buildup term of R_pol_tot: each additional internally-
% reflected contribution corresponds to one full round trip through the
% slab, i.e. two single passes.)
%
% Non-polarised: average of s and p.  Absorptance: A = 1 - R - T.

Nw = numel(n_complex);
Na = numel(angles_rad);
Nt = numel(thicknesses_um);

%% Expand dimensions

N      = repmat(n_complex(:),1,Na);
lambda = repmat(wavelength_um(:),1,Na);

sin1 = repmat(sin(angles_rad),Nw,1);
cos1 = repmat(cos(angles_rad),Nw,1);

%% Snell law

cos2 = sqrt(1-(sin1./N).^2);

%% Fresnel reflectance

rs = (cos1-N.*cos2)./(cos1+N.*cos2);
rp = (N.*cos1-cos2)./(N.*cos1+cos2);

Rs = abs(rs).^2;
Rp = abs(rp).^2;

%% Propagation inside material
%
% kz = (2*pi/lambda)*sqrt(N^2-sin^2(theta))

kz = (2*pi./lambda).*sqrt(N.^2-sin1.^2);
alpha = imag(kz);

%% Output arrays

R_slab = zeros(Nw,Na,Nt);
T_slab = zeros(Nw,Na,Nt);

%% Thickness loop

for it = 1:Nt

    d = thicknesses_um(it);

    % Single-pass attenuation, and the round-trip (two-pass) attenuation
    % used in the internal-reflection buildup terms.
    tau  = exp(-2*alpha*d);
    tau2 = tau.^2;

    % avoid numerical singularities
    denomS = max(1-Rs.^2.*tau2,1e-30);
    denomP = max(1-Rp.^2.*tau2,1e-30);

    %% s polarization
    Ts = ((1-Rs).^2 .* tau)./denomS;
    Rs_tot = Rs + Rs.*(1-Rs).^2.*tau2./denomS;

    %% p polarization
    Tp = ((1-Rp).^2 .* tau)./denomP;
    Rp_tot = Rp + Rp.*(1-Rp).^2.*tau2./denomP;

    %% Unpolarised
    R_slab(:,:,it) = 0.5*(Rs_tot+Rp_tot);
    T_slab(:,:,it) = 0.5*(Ts+Tp);

end

%% Absorptance

A_slab = 1-R_slab-T_slab;

%% Final numerical corrections

R_slab(~isfinite(R_slab)) = 0;
T_slab(~isfinite(T_slab)) = 0;
A_slab(~isfinite(A_slab)) = 0;

R_slab = max(min(R_slab,1),0);
T_slab = max(min(T_slab,1),0);
A_slab = max(min(A_slab,1),0);

end