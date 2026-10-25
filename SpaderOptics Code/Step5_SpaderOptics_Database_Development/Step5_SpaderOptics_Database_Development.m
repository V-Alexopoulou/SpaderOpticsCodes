%% Step5_SpaderOptics_Database_Development

%% Remove sources that were excluded from the calculations (per material)
% {material, author} pairs. Only these exact combinations are removed.
exclude_pairs = { ...
    'Al', 'Larruquert'; ...
    'Cu', 'Brimhall' };

for p = 1:size(exclude_pairs,1)

    mat_name    = exclude_pairs{p,1};
    author_name = exclude_pairs{p,2};

    imat = find(strcmpi({materials.material}, mat_name), 1);

    if isempty(imat)
        warning('Material %s not found. Skipping.', mat_name)
        continue
    end

    names  = {materials(imat).sources.name};
    remove = contains(names, author_name, 'IgnoreCase', true);

    if any(remove)
        fprintf('%s: removing %s\n', materials(imat).material, ...
            strjoin(names(remove), ', '));
        materials(imat).sources(remove) = [];
    else
        fprintf('%s: %s not present (already removed?)\n', ...
            materials(imat).material, author_name);
    end

end

% Safety check: every struct must now have the same sources per material
for i = 1:numel(materials)

    n = [numel(materials(i).sources), ...
         numel(solar_flux(i).sources), ...
         numel(emissivity(i).sources), ...
         numel(fourier(i).sources)];

    if any(n ~= n(1))
        error('Source count mismatch for %s: materials=%d, solar_flux=%d, emissivity=%d, fourier=%d', ...
            materials(i).material, n(1), n(2), n(3), n(4));
    end

end

%% Create the database struct
SpaderOptics = struct( ...
    'material', {}, ...
    'dataset', {} );

for i = 1:numel(materials)

    % -------------------------------------------------------------
    % Material name
    % -------------------------------------------------------------
    SpaderOptics(i).material = materials(i).material;

    sources_materials = materials(i).sources;
    sources_solarflux = solar_flux(i).sources;
    sources_emissivity = emissivity(i).sources;
    sources_fourier = fourier(i).sources;
    

    % -------------------------------------------------------------
    % Preallocate publication dataset
    % -------------------------------------------------------------
    dataset = struct( ...
        'author', {}, ...
        'reference', {}, ...
        'source_url', {}, ...
        'upper_temperature_limit_K', {}, ...
        'thermal_conductivity_temperature_grid_K', {}, ...
        'thermal_conductivity_W_mK', {}, ...
        'thermal_conductivity_reference', {}, ...
        'original_temperature_K', {}, ...
        'original_thickness_um', {}, ...
        'original_wavelength_um', {}, ...
        'original_n_complex', {}, ...
        'resampled_wavelength_um', {}, ...
        'resampled_n_complex', {}, ...
        'calculation_incident_angles_deg', {}, ...
        'calculation_thickness_um', {}, ...
        'calculation_temperature_K', {}, ...
        'calculation_solar_concentration_multipliers', {}, ...
        'calculated_R_directional', {}, ...
        'calculated_T_directional', {}, ...
        'calculated_A_directional', {}, ...
        'calculated_E_directional', {}, ...
        'calculated_R_hemispherical', {}, ...
        'calculated_T_hemispherical', {}, ...
        'calculated_A_hemispherical', {}, ...
        'calculated_E_hemispherical', {}, ...
        'calculated_directional_Reflected_flux_W_m2', {}, ...
        'calculated_directional_Transmitted_flux_W_m2', {}, ...
        'calculated_directional_Absorbed_flux_W_m2', {}, ...
        'calculated_hemispherical_Reflected_flux_W_m2', {}, ...
        'calculated_hemispherical_Transmitted_flux_W_m2', {}, ...
        'calculated_hemispherical_Absorbed_flux_W_m2', {}, ...
        'calculated_hemispherical_Emitted_flux_W_m2', {}, ...
        'calculated_lumped_equilibrium_temperature_dirA_hemiE_K', {}, ...
        'calculated_lumped_equilibrium_temperature_hemiA_hemiE_K', {}, ...
        'calculated_Conductive_flux_dirA_hemiE_W_m2', {}, ...
        'calculated_Conductive_flux_hemiA_hemiE_W_m2', {}, ...
        'calculated_slab_equilibrium_temperature_front_dirA_hemiE_K', {}, ...
        'calculated_slab_equilibrium_temperature_front_hemiA_hemiE_K', {}, ...
        'calculated_slab_equilibrium_temperature_rear_dirA_hemiE_K', {}, ...
        'calculated_slab_equilibrium_temperature_rear_hemiA_hemiE_K', {}, ...
        'calculated_slab_equilibrium_temperature_mid_dirA_hemiE_K', {}, ...
        'calculated_slab_equilibrium_temperature_mid_hemiA_hemiE_K', {}, ...
        'metrics_solar_coverage_ratio', {}, ...
        'metrics_thermal_coverage_ratio', {});

    for j = 1:numel(sources_materials)

        src_materials = sources_materials(j);
        

        % ---------------------------------------------------------
        % Author
        % ---------------------------------------------------------
        dataset(j).author = src_materials.name;

        % ---------------------------------------------------------
        % Reference
        % ---------------------------------------------------------
        dataset(j).reference = src_materials.reference;

        % ---------------------------------------------------------
        % DOI / source URL
        % ---------------------------------------------------------
        dataset(j).source_url = src_materials.source_url;


        % ---------------------------------------------------------
        % Original temperature
        %
        % If temperature_K exists, use it.
        % If temperature_K is NaN, use exactly what is stored in
        % the original 'temperature' field.
        % ---------------------------------------------------------
        if isfield(src_materials, 'temperature_K') && ~isempty(src_materials.temperature_K)
        
            T_K = src_materials.temperature_K;
        
            if isnumeric(T_K) && any(isfinite(T_K(:)))
        
                % At least one valid numerical temperature exists
                dataset(j).original_temperature_K = T_K;
        
            else
        
                % No valid numerical temperature -> preserve original text
                dataset(j).original_temperature_K = src_materials.temperature;
        
            end
        
        else
        
            % temperature_K does not exist or is empty
            dataset(j).original_temperature_K = src_materials.temperature;
        
        end

        % ---------------------------------------------------------
        % Original thickness
        %
        % NaN -> empty string
        % ---------------------------------------------------------
        if isfield(src_materials, 'film_thickness_um') && ...
                ~isempty(src_materials.film_thickness_um) && ...
                ~isnan(src_materials.film_thickness_um)

            dataset(j).original_thickness_um = ...
                src_materials.film_thickness_um;

        else

            dataset(j).original_thickness_um = '';

        end

        % ---------------------------------------------------------
        % Original wavelength grid
        % ---------------------------------------------------------
        dataset(j).original_wavelength_um = ...
            src_materials.wavelength;

        % ---------------------------------------------------------
        % Original complex refractive index
        % ---------------------------------------------------------
        dataset(j).original_n_complex = src_materials.n_complex;

    end

    for j = 1:numel(sources_solarflux)

        src_solarflux = sources_solarflux(j);
        
        % ---------------------------------------------------------
        % Resampled Wavelength
        % ---------------------------------------------------------
        
        dataset(j).resampled_wavelength_um = src_solarflux.wavelength_um;

        % ---------------------------------------------------------
        % Resampled n_complex
        % ---------------------------------------------------------
        
        dataset(j).resampled_n_complex = src_solarflux.n_complex_resampled;

        % ---------------------------------------------------------
        % Incident angles for calculations
        % ---------------------------------------------------------
        
        dataset(j).calculation_incident_angles_deg = src_solarflux.angles_deg;

        % ---------------------------------------------------------
        % Foil thicknesses for calculations
        % ---------------------------------------------------------
        
        dataset(j).calculation_thickness_um = src_solarflux.thicknesses_um;

        % ---------------------------------------------------------
        % Directional Reflectivity (for non- T dependent cases)
        % ---------------------------------------------------------
        
        dataset(j).calculated_R_directional = src_solarflux.R_lambda;

        % ---------------------------------------------------------
        % Directional Transmissivity (for non- T dependent cases)
        % ---------------------------------------------------------
        
        dataset(j).calculated_T_directional = src_solarflux.T_lambda;

        % ---------------------------------------------------------
        % Directional Absorptivity (for non- T dependent cases)
        % ---------------------------------------------------------
        
        dataset(j).calculated_A_directional = src_solarflux.A_lambda;

        % ---------------------------------------------------------
        % Hemispherical Reflectivity (for non- T dependent cases)
        % ---------------------------------------------------------
        
        dataset(j).calculated_R_hemispherical = src_solarflux.R_hemi_lambda;

        % ---------------------------------------------------------
        % Hemispherical Transmissivity (for non- T dependent cases)
        % ---------------------------------------------------------
        
        dataset(j).calculated_T_hemispherical = src_solarflux.T_hemi_lambda;

        % ---------------------------------------------------------
        % Hemispherical Absorptivity (for non- T dependent cases)
        % ---------------------------------------------------------
        
        dataset(j).calculated_A_hemispherical = src_solarflux.A_hemi_lambda;

        % ---------------------------------------------------------
        % Directional Reflected Flux (for non- T dependent cases)
        % ---------------------------------------------------------
        
        dataset(j).calculated_directional_Reflected_flux_W_m2 = src_solarflux.reflected_flux;

        % ---------------------------------------------------------
        % Directional Transmitted Flux (for non- T dependent cases)
        % ---------------------------------------------------------
        
        dataset(j).calculated_directional_Transmitted_flux_W_m2 = src_solarflux.transmitted_flux;

        % ---------------------------------------------------------
        % Directional Absorbed Flux (for non- T dependent cases)
        % ---------------------------------------------------------
        
        dataset(j).calculated_directional_Absorbed_flux_W_m2 = src_solarflux.absorbed_flux;

        % ---------------------------------------------------------
        % Hemispherical Reflected Flux (for non- T dependent cases)
        % ---------------------------------------------------------
        
        dataset(j).calculated_hemispherical_Reflected_flux_W_m2 = src_solarflux.reflected_flux_hemi;

        % ---------------------------------------------------------
        % Hemispherical Transmitted Flux (for non- T dependent cases)
        % ---------------------------------------------------------
        
        dataset(j).calculated_hemispherical_Transmitted_flux_W_m2 = src_solarflux.transmitted_flux_hemi;

        % ---------------------------------------------------------
        % Hemispherical Absorbed Flux (for non- T dependent cases)
        % ---------------------------------------------------------
        
        dataset(j).calculated_hemispherical_Absorbed_flux_W_m2 = src_solarflux.absorbed_flux_hemi;
        
        % ---------------------------------------------------------
        % Solar coverage ratio (metrics)
        % ---------------------------------------------------------
        
        dataset(j).metrics_solar_coverage_ratio = src_solarflux.coverage_ratio; 

    end
    
    for j = 1:numel(sources_emissivity)

        src_emissivity = sources_emissivity(j);

        
        % ---------------------------------------------------------
        % Melting or Degradation points
        % ---------------------------------------------------------
        
        dataset(j).upper_temperature_limit_K = src_emissivity.melting_point_K;

        % ---------------------------------------------------------
        % Temperatures for calculations
        % ---------------------------------------------------------
        
        dataset(j).calculation_temperature_K = src_emissivity.T_range;

        % ---------------------------------------------------------
        % Solar concentration multipliers
        % ---------------------------------------------------------
        
        dataset(j).calculation_solar_concentration_multipliers = src_emissivity.multipliers;

        
        % ---------------------------------------------------------
        % Calculated directional emissivity
        % ---------------------------------------------------------
        
        dataset(j).calculated_E_directional = src_emissivity.eps_directional;

        % ---------------------------------------------------------
        % Calculated hemispherical emissivity
        % ---------------------------------------------------------
        
        dataset(j).calculated_E_hemispherical = src_emissivity.eps_hemispherical;

        % ---------------------------------------------------------
        % Calculated hemispherical emitted flux
        % ---------------------------------------------------------
        
        dataset(j).calculated_hemispherical_Emitted_flux_W_m2 = src_emissivity.q_emitted_hemispherical;

        % ---------------------------------------------------------
        % Calculated lumped equilibrium temperature (directional absorption - hemispherical emittance)
        % ---------------------------------------------------------
        
        dataset(j).calculated_lumped_equilibrium_temperature_dirA_hemiE_K = src_emissivity.T_eq_dir_hemi;

        % ---------------------------------------------------------
        % Calculated lumped equilibrium temperature (hemispherical absorption - hemispherical emittance)
        % ---------------------------------------------------------
        
        dataset(j).calculated_lumped_equilibrium_temperature_hemiA_hemiE_K = src_emissivity.T_eq_hemi_hemi;

                % ---------------------------------------------------------
        % Thermal coverage ratio (metrics), capped at 1
        % ---------------------------------------------------------

        tcr = src_emissivity.thermal_coverage_ratio;

        if isnumeric(tcr)

            above = tcr > 1;          % NaN > 1 is false, so NaNs are kept as NaN

            if any(above(:))
                fprintf('%s | %s: thermal coverage ratio above 1 in %d of %d values (max %.4f) -> set to 1\n', ...
                    materials(i).material, src_emissivity.name, ...
                    nnz(above), numel(tcr), max(tcr(above)));

                tcr(above) = 1;
            end

        end

        dataset(j).metrics_thermal_coverage_ratio = tcr;

       
    end

    % -------------------------------------------------------------
    % Replace temperature dependent quantities for T-dependent sources
    % -------------------------------------------------------------
    
    for j = 1:numel(sources_emissivity)
    
        src_emissivity = sources_emissivity(j);
    
        % Check whether this is a T-dependent source
        if endsWith(src_emissivity.name, '-Tdepend')
    
            % Store T-dependent author name
            author_Tdepend = src_emissivity.name;
           
    
            % Store T-dependent directional reflectivity temporarily
            R_lambda_T_temp = src_emissivity.R_lambda_T;
            T_lambda_T_temp = src_emissivity.T_lambda_T;
            A_lambda_T_temp = src_emissivity.A_lambda_T;
            R_lambda_hemi_T_temp = src_emissivity.R_lambda_hemi_T;
            T_lambda_hemi_T_temp = src_emissivity.T_lambda_hemi_T;
            A_lambda_hemi_T_temp = src_emissivity.A_lambda_hemi_T;
            reflected_flux_T_temp = src_emissivity.q_reflected_directional_T;
            transmitted_flux_T_temp = src_emissivity.q_transmitted_directional_T;
            absorbed_flux_T_temp = src_emissivity.q_absorbed_directional_T;
            reflected_flux_hemi_T_temp = src_emissivity.q_reflected_hemispherical_T;
            transmitted_flux_hemi_T_temp = src_emissivity.q_transmitted_hemispherical_T;
            absorbed_flux_hemi_T_temp = src_emissivity.q_absorbed_hemispherical_T;
    
            % -----------------------------------------------------
            % Find the corresponding source in SpaderOptics.dataset
            % -----------------------------------------------------
            idx_dataset = find( ...
                strcmp({dataset.author}, author_Tdepend), ...
                1, 'first');
    
            if ~isempty(idx_dataset)
    
                % Replace the solar_flux quantity with the emissivity quantity
                dataset(idx_dataset).calculated_R_directional = R_lambda_T_temp;
                dataset(idx_dataset).calculated_T_directional = T_lambda_T_temp;
                dataset(idx_dataset).calculated_A_directional = A_lambda_T_temp;
                dataset(idx_dataset).calculated_R_hemispherical = R_lambda_hemi_T_temp;
                dataset(idx_dataset).calculated_T_hemispherical = T_lambda_hemi_T_temp;
                dataset(idx_dataset).calculated_A_hemispherical = A_lambda_hemi_T_temp;
                dataset(idx_dataset).calculated_directional_Reflected_flux_W_m2 = reflected_flux_T_temp;
                dataset(idx_dataset).calculated_directional_Transmitted_flux_W_m2 = transmitted_flux_T_temp;
                dataset(idx_dataset).calculated_directional_Absorbed_flux_W_m2 = absorbed_flux_T_temp;
                dataset(idx_dataset).calculated_hemispherical_Reflected_flux_W_m2 = reflected_flux_hemi_T_temp;
                dataset(idx_dataset).calculated_hemispherical_Transmitted_flux_W_m2 = transmitted_flux_hemi_T_temp;
                dataset(idx_dataset).calculated_hemispherical_Absorbed_flux_W_m2 = absorbed_flux_hemi_T_temp;
            end
    
        end

    end

    for j = 1:numel(sources_fourier)

        src_fourier = sources_fourier(j);

        % ---------------------------------------------------------
        % Thermal conductivity data
        % ---------------------------------------------------------
        
        dataset(j).thermal_conductivity_temperature_grid_K = src_fourier.T_thermal_conductivity; 
        dataset(j).thermal_conductivity_W_mK = src_fourier.k_thermal_conductivity;
        dataset(j).thermal_conductivity_reference = char(src_fourier.ref_thermal_conductivity);

        % ---------------------------------------------------------
        % Conductive heat flux
        % --------------------------------------------------------- 
        dataset(j).calculated_Conductive_flux_dirA_hemiE_W_m2 = src_fourier.q_cond_dir_hemi;
        dataset(j).calculated_Conductive_flux_hemiA_hemiE_W_m2 = src_fourier.q_cond_hemi_hemi;

        % ---------------------------------------------------------
        % Equilibrium temperature at the front surface
        % ---------------------------------------------------------
        dataset(j).calculated_slab_equilibrium_temperature_front_dirA_hemiE_K = src_fourier.T_front_dir_hemi;
        dataset(j).calculated_slab_equilibrium_temperature_front_hemiA_hemiE_K = src_fourier.T_front_hemi_hemi;

        % ---------------------------------------------------------
        % Equilibrium temperature at the rear surface
        % ---------------------------------------------------------
        dataset(j).calculated_slab_equilibrium_temperature_rear_dirA_hemiE_K = src_fourier.T_rear_dir_hemi;
        dataset(j).calculated_slab_equilibrium_temperature_rear_hemiA_hemiE_K = src_fourier.T_rear_hemi_hemi;

        % ---------------------------------------------------------
        % Equilibrium temperature at the middle of the slab
        % ---------------------------------------------------------
        dataset(j).calculated_slab_equilibrium_temperature_mid_dirA_hemiE_K = src_fourier.T_mid_dir_hemi;
        dataset(j).calculated_slab_equilibrium_temperature_mid_hemiA_hemiE_K = src_fourier.T_mid_hemi_hemi;

    end


        % -------------------------------------------------------------
    % Spread of calculated hemispherical emissivity between sources
    % -------------------------------------------------------------
    exclude_from_spread = {'-Tdepend'};     % author-name patterns to leave out

    keep = true(1, numel(dataset));
    for k = 1:numel(exclude_from_spread)
        keep = keep & ~contains({dataset.author}, exclude_from_spread{k});
    end

    Eh  = {dataset(keep).calculated_E_hemispherical};
    Tg  = {dataset(keep).calculation_temperature_K};
    Thk = {dataset(keep).calculation_thickness_um};

    spread = struct();

    same = numel(Eh) > 1 && ...
           all(cellfun(@(x) isequal(size(x), size(Eh{1})), Eh)) && ...
           all(cellfun(@(x) isequal(x, Tg{1}),  Tg)) && ...
           all(cellfun(@(x) isequal(x, Thk{1}), Thk));

    if same

        E3 = cat(3, Eh{:});                        % nT x nThickness x nSources

        spread.description   = ['Spread of calculated hemispherical emissivity ' ...
                                'across sources (rows: temperature, columns: thickness)'];
        spread.authors       = {dataset(keep).author};
        spread.temperature_K = Tg{1}(:);
        spread.thickness_um  = Thk{1}(:).';
        spread.min           = min(E3, [], 3);
        spread.max           = max(E3, [], 3);
        spread.mean          = mean(E3, 3, 'omitnan');
        spread.median        = median(E3, 3, 'omitnan');
        spread.std           = std(E3, 0, 3, 'omitnan');
        spread.mad           = median(abs(E3 - spread.median), 3, 'omitnan');
        spread.range_abs     = spread.max - spread.min;
        spread.range_rel_pct = 100*spread.range_abs ./ spread.median;

    else

        warning('%s: fewer than 2 comparable sources or different grids, spread not computed.', ...
            materials(i).material)

    end

    SpaderOptics(i).source_spread = spread;

    % -------------------------------------------------------------
    % Store nested source dataset
    % -------------------------------------------------------------
    SpaderOptics(i).dataset = dataset;

end

% 1. Save in the current folder (folder Step5_SpaderOptics_Database_Development)
save('SpaderOptics.mat', 'SpaderOptics', '-v7.3');

% 2. Save in the next folder (folder Step6_SpaderOptics_App_Development)(using a relative path)
% '..' moves up one level to the parent folder, then we enter 'Step6_SpaderOptics_App_Development'
pathTo_Step6_SpaderOptics_App_Development = fullfile('..', 'Step6_SpaderOptics_App_Development', 'SpaderOptics.mat');
save(pathTo_Step6_SpaderOptics_App_Development, 'SpaderOptics', '-v7.3');
