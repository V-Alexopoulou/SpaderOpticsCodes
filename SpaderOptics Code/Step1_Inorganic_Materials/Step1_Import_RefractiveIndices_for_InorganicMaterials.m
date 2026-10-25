%% Step1_Import_RefractiveIndices_for_InorganicMaterials.m
%  =========================================================================
%  STEP 1 — Import inorganic-material optical constants
%
%  Parses every .yml file from the RefractiveIndex.INFO database that lives
%  under <rootDir>, grouped by material sub-folders (one folder per
%  material, one .yml file per literature source).
%
%  For each file, this script:
%    1. parses the YAML DATA block(s) into wavelength-dependent n(lambda)
%       and k(lambda) (tabulated or evaluated from a dispersion formula),
%    2. validates and cleans the resulting (lambda, n, k) triples (removes
%       non-finite/negative values, enforces a monotonically increasing
%       wavelength grid, merges duplicate wavelength points by averaging),
%    3. extracts bibliographic/metadata fields (reference, source URL,
%       temperature, film thickness),
%    4. groups same-author sources measured at different temperatures into
%       an additional combined "-Tdepend" entry.
%
%  This script does NOT compute any optical response quantity (Fresnel or
%  slab reflectance/transmittance/absorptance). That calculation is
%  performed once, from the raw (lambda, n_complex) data stored here, by
%  Step3_Solar_Matter_Interaction.m — on its own resampled wavelength and angle grid.
%  Computing it here as well would duplicate that calculation on a
%  different grid and leave results that are never read by any downstream
%  script.
%
%  OUTPUT — struct array  Inorganic_materials(m)  with fields:
%     .material       — folder (material) name  (e.g. 'Al', 'SiO2')
%     .sources(s)     — sub-struct for each literature source:
%
%        --- metadata ---
%        .name              — source filename without extension
%        .reference         — bibliographic reference string from the file
%        .source_url        — URL extracted from the href tag in REFERENCES
%        .temperature       — temperature string extracted from COMMENTS
%                             (e.g. 'room temperature', '300 K', '25 C')
%        .temperature_K     — temperature in Kelvin, parsed from the source
%                             filename (e.g. "Ferrera-298K" -> 298), NaN if
%                             not encoded in the name
%        .film_thickness_um — film thickness in micrometres, extracted from
%                             the CONDITIONS field (NaN if not specified)
%
%        --- optical constants ---
%        .wavelength        — [Nw x 1]   wavelength in micrometres
%        .n_complex         — [Nw x 1]   complex refractive index  n + i*k
%
%        --- grid for the later slab-optics calculation ---
%        .thicknesses_um    — [1 x Nt]  slab thicknesses [um], passed
%                             through unchanged for compute_solar_flux.m
%
%  Supported YML DATA types
%  -------------------------
%     tabulated nk  — three-column table (lambda, n, k)
%     tabulated n   — two-column table   (lambda, n)
%     tabulated k   — two-column table   (lambda, k)
%     formula 1-9   — all RefractiveIndex.INFO dispersion formulas
%
%  A single .yml file may contain MULTIPLE data blocks (e.g. one formula
%  block for n and one tabulated-k block for k). The script merges them by
%  interpolating the shorter dataset onto the wavelength grid of the first.
%  =========================================================================

clear; clc; close all;

%% ==== USER SETTINGS ======================================================
rootDir        = pwd;                                     % root folder with material sub-dirs
defaultNw      = 1000;                                     % wavelength points for formulas
thicknesses_um = [0.001 0.01 0.1 1 10 100 1000 10000];      % slab thicknesses [um]
% ==========================================================================

%% ==== DISCOVER MATERIAL FOLDERS ==========================================
entries = dir(rootDir);
entries = entries([entries.isdir]);
entries = entries(~startsWith({entries.name}, '.'));

Inorganic_materials = struct([]);
mIdx = 0;

for d = 1:numel(entries)
    folderName = entries(d).name;
    folderPath = fullfile(rootDir, folderName);

    ymlFiles = [dir(fullfile(folderPath, '*.yml')); ...
                dir(fullfile(folderPath, '*.yaml'))];
    if isempty(ymlFiles), continue; end

    mIdx = mIdx + 1;
    Inorganic_materials(mIdx).material = folderName;
    sources = struct([]);
    sIdx = 0;

    for f = 1:numel(ymlFiles)
        filePath = fullfile(folderPath, ymlFiles(f).name);
        [~, srcName] = fileparts(ymlFiles(f).name);
        fprintf('Processing: %s / %s\n', folderName, ymlFiles(f).name);

        try
            [wl, n_re, k_im, refStr, comStr, condStr] = ...
                    parseYML(filePath, defaultNw, folderName, srcName);
        catch ME
            warning('  Skipped  %s : %s', ymlFiles(f).name, ME.message);
            continue
        end
        if isempty(wl), continue; end

        n_cplx = complex(n_re, k_im);

        %% ----------------------------------------------------------------
        % Metadata
        %% ----------------------------------------------------------------

        sIdx = sIdx + 1;
        sources(sIdx).name      = srcName;
        sources(sIdx).reference = refStr;

        hrefTok = regexp(refStr, 'href\s*=\s*"([^"]+)"', 'tokens', 'once');
        if ~isempty(hrefTok)
            sources(sIdx).source_url = hrefTok{1};
        else
            sources(sIdx).source_url = '';
        end

        sources(sIdx).temperature       = extractTemperature(comStr);
        sources(sIdx).film_thickness_um = extractFilmThickness(condStr);
        sources(sIdx).temperature_K     = extractTemperatureFromName(srcName);

        %% ----------------------------------------------------------------
        % Optical constants and the thickness grid used later by
        % compute_solar_flux.m
        %% ----------------------------------------------------------------

        sources(sIdx).wavelength     = wl(:);
        sources(sIdx).n_complex      = n_cplx(:);
        sources(sIdx).thicknesses_um = thicknesses_um;
    end

    %% ====================================================================
    % CREATE TEMPERATURE-DEPENDENT SOURCES
    %% ====================================================================
    sources = createTemperatureDependentSources(sources);
    Inorganic_materials(mIdx).sources = sources;
end

%% ==== SUMMARY ============================================================
fprintf('\n========== Step 1 - Import complete: %d Inorganic material(s) ==========\n', numel(Inorganic_materials));
for m = 1:numel(Inorganic_materials)
    nSrc = numel(Inorganic_materials(m).sources);
    fprintf('  %-12s : %d source(s)\n', Inorganic_materials(m).material, nSrc);
    for s = 1:nSrc
        src = Inorganic_materials(m).sources(s);

        tempStr = src.temperature;
        if isempty(tempStr)
            tempStr = '(not specified)';
        end

        if isnan(src.temperature_K)
            tempKStr = '(no T)';
        else
            tempKStr = sprintf('%.2f K', src.temperature_K);
        end

        if isnan(src.film_thickness_um)
            thkStr = '';
        else
            thkStr = sprintf('  d=%.4g um', src.film_thickness_um);
        end

        fprintf('      %-25s  wl=[%.4f-%.4f] um  (%d pts)  T=%s  nameT=%s%s\n', ...
            src.name, min(src.wavelength), max(src.wavelength), ...
            numel(src.wavelength), tempStr, tempKStr, thkStr);
    end
end

clearvars -except Inorganic_materials;

return


%% #########################################################################
%%  LOCAL FUNCTIONS
%% #########################################################################

function [wl, n, k, refStr, comStr, condStr] = parseYML(filePath, defaultNw, materialName, sourceName)
%PARSEYML  Read one RefractiveIndex.INFO .yml file.
%  Returns wavelength (um), real part n, extinction coefficient k,
%  the REFERENCES string, the COMMENTS string, and the CONDITIONS string.
%  Also validates and cleans (lambda, n, k): removes non-finite/negative
%  values, sorts and de-duplicates the wavelength grid.

    raw = fileread(filePath);

    refStr  = extractMultilineField(raw, 'REFERENCES');
    comStr  = extractMultilineField(raw, 'COMMENTS');
    condStr = extractMultilineField(raw, 'CONDITIONS');

    blocks = splitDataBlocks(raw);
    if isempty(blocks)
        error('No DATA blocks found in file.');
    end

    wl = []; n = []; k = [];

    for b = 1:numel(blocks)
        blk = blocks{b};
        dtype = lower(strtrim(extractScalarField(blk, 'type')));
        if isempty(dtype), continue; end

        if contains(dtype, 'tabulated')
            [wl_b, n_b, k_b] = parseTabulated(blk, dtype);
        elseif startsWith(dtype, 'formula')
            [wl_b, n_b, k_b] = parseFormula(blk, dtype, defaultNw);
        else
            continue
        end
        if isempty(wl_b), continue; end

        if isempty(wl)
            wl = wl_b;  n = n_b;  k = k_b;
        else
            if any(n_b ~= 0) && all(n == 0)
                n = interp1(wl_b, n_b, wl, 'linear', 'extrap');
            end
            if any(k_b ~= 0) && all(k == 0)
                k = interp1(wl_b, k_b, wl, 'linear', 0);
            end
        end
    end

    %% ============================================================
    %  SANITY CHECKS & AUTOMATIC CLEANING
    %% ============================================================

    % ---------- 1. Check n and k ---------------------------------
    valid = isfinite(n) & isfinite(k) & (n >= 0) & (k >= 0);

    if any(~valid)

        fprintf('\nProblem: %s / %s\n', materialName, sourceName);

        bad = find(~valid);

        for ii = 1:numel(bad)

            idx = bad(ii);

            if ~isfinite(n(idx))
                fprintf('  NaN/Inf n : ');
            elseif ~isfinite(k(idx))
                fprintf('  NaN/Inf k : ');
            elseif n(idx) < 0
                fprintf('  negative n : ');
            else
                fprintf('  negative k : ');
            end

            fprintf('lambda = %.8g um   n = %.8g   k = %.8g\n', ...
                wl(idx), n(idx), k(idx));

        end

    end

    % Remove invalid rows
    wl = wl(valid);
    n  = n(valid);
    k  = k(valid);

    % ---------- 2. Non-monotonic wavelengths ---------------------

    dwl = diff(wl);

    if any(dwl <= 0)

        fprintf('\nProblem: %s / %s\n', materialName, sourceName);
        fprintf('  Non-monotonic wavelength grid:\n');

        idx = find(dwl <= 0);

        for ii = 1:numel(idx)

            fprintf('    %.8g -> %.8g um\n', ...
                wl(idx(ii)), wl(idx(ii)+1));

        end

    end

    % ---------- 3. Duplicate wavelengths -------------------------

    [u,~,ic] = unique(wl);

    counts = accumarray(ic,1);

    dupGroups = find(counts > 1);

    if ~isempty(dupGroups)

        fprintf('\nProblem: %s / %s\n', materialName, sourceName);
        fprintf('  Duplicate wavelength(s):\n');

        for g = dupGroups.'

            rows = find(ic==g);

            fprintf('    lambda = %.8g um\n', u(g));

            for r = rows.'

                fprintf('       n = %.8g   k = %.8g\n', ...
                    n(r),k(r));

            end

        end

        % Average duplicates
        n = accumarray(ic,n,[],@mean);
        k = accumarray(ic,k,[],@mean);
        wl = u;

    end

    % ---------- 4. Sort wavelengths ------------------------------

    [wl,idx] = sort(wl);

    n = n(idx);
    k = k(idx);

    if isempty(wl)
        error('Could not extract optical data from file.');
    end
end


function tempStr = extractTemperature(comStr)
%EXTRACTTEMPERATURE  Find temperature info in the COMMENTS string.
%  Handles: "room temperature", "293 K", "100 K (-173 C)", "T = 300 K",
%  "293K" (no space), etc.

    tempStr = '';
    if isempty(comStr), return; end

    if contains(lower(comStr), 'room temp')
        tempStr = 'room temperature';
        return
    end

    safe = comStr;
    safe(safe > char(127)) = ' ';
    safe = strrep(safe, '=', ' ');
    safe = strrep(safe, ':', ' ');
    safe = strrep(safe, '(', ' ');
    safe = strrep(safe, ')', ' ');
    safe = strrep(safe, ',', ' ');
    safe = strrep(safe, ';', ' ');

    tokens = strsplit(strtrim(safe));

    for i = 1:numel(tokens)
        tok = tokens{i};

        if numel(tok) == 1 && any(tok == 'K' | tok == 'C' | tok == 'F')
            if i > 1
                num = str2double(tokens{i-1});
                if ~isnan(num)
                    tempStr = [tokens{i-1} ' ' tok];
                    return
                end
            end
        end

        if numel(tok) >= 2
            lastChar = tok(end);
            if any(lastChar == 'K' | lastChar == 'C' | lastChar == 'F')
                numPart = tok(1:end-1);
                num = str2double(numPart);
                if ~isnan(num)
                    tempStr = [numPart ' ' lastChar];
                    return
                end
            end
        end
    end
end


function thickness_um = extractFilmThickness(condStr)
%EXTRACTFILMTHICKNESS  Extract film_thickness from the CONDITIONS block.
%  The value in the YML is in metres (e.g. 20e-9 for 20 nm).
%  Returns the thickness converted to micrometres, or NaN if not found.

    thickness_um = NaN;
    if isempty(condStr), return; end

    lines = splitlines(condStr);
    for i = 1:numel(lines)
        ln = strtrim(lines{i});
        if startsWith(ln, 'film_thickness')
            colonIdx = strfind(ln, ':');
            if isempty(colonIdx), continue; end
            valStr = strtrim(ln(colonIdx(1)+1:end));
            val = str2double(valStr);
            if ~isnan(val)
                thickness_um = val * 1e6;   % metres -> micrometres
            end
            return
        end
    end
end


function blocks = splitDataBlocks(raw)
%SPLITDATABLOCKS  Locate individual "- type: ..." sub-blocks inside DATA.
    blocks = {};
    expr = '(?m)^\s*-\s*type\s*:';
    starts = regexp(raw, expr);
    if isempty(starts), return; end
    for i = 1:numel(starts)
        if i < numel(starts)
            blocks{end+1} = raw(starts(i):starts(i+1)-1); %#ok<AGROW>
        else
            blocks{end+1} = raw(starts(i):end); %#ok<AGROW>
        end
    end
end


function val = extractScalarField(txt, field)
%EXTRACTSCALARFIELD  Get a single-line YAML field value.
    tok = regexp(txt, ['(?m)' field '\s*:\s*"?([^"\r\n]+)"?'], 'tokens', 'once');
    if isempty(tok), val = ''; else, val = strtrim(tok{1}); end
end


function val = extractMultilineField(txt, field)
%EXTRACTMULTILINEFIELD  Get a YAML block-scalar (field: |) or bare-mapping
%  (field:) value.

    val = '';
    lines = splitlines(txt);
    collecting = false;
    parts = {};

    for i = 1:numel(lines)
        ln = lines{i};
        stripped = strtrim(ln);

        if collecting
            if ~isempty(ln) && (ln(1) == ' ' || ln(1) == char(9))
                parts{end+1} = stripped; %#ok<AGROW>
            else
                break
            end
        else
            if startsWith(stripped, field)
                rest = strtrim(stripped(numel(field)+1:end));
                if startsWith(rest, ':')
                    rest2 = strtrim(rest(2:end));
                    if isempty(rest2) || startsWith(rest2, '|')
                        collecting = true;
                    end
                end
            end
        end
    end

    if ~isempty(parts)
        val = strjoin(parts, newline);
    end
end


function [wl, n, k] = parseTabulated(blk, dtype)
%PARSETABULATED  Handle tabulated nk / tabulated n / tabulated k.

    tok = regexp(blk, 'data\s*:\s*\|?\s*\n([\s\S]+)', 'tokens', 'once');
    if isempty(tok), wl=[]; n=[]; k=[]; return; end

    lines = splitlines(strtrim(tok{1}));
    nums = [];
    for i = 1:numel(lines)
        ln = strtrim(lines{i});
        if isempty(ln) || ln(1) == '#', continue; end
        if ~isempty(regexp(ln, '^[a-zA-Z\-]', 'once')), break; end
        row = sscanf(ln, '%f')';
        if ~isempty(row)
            nums(end+1, 1:numel(row)) = row; %#ok<AGROW>
        end
    end

    if isempty(nums), wl=[]; n=[]; k=[]; return; end

    wl = nums(:,1);

    if contains(dtype, 'nk')
        n = nums(:,2);
        k = nums(:,3);
    elseif endsWith(strtrim(dtype), 'n')
        n = nums(:,2);
        k = zeros(size(wl));
    elseif endsWith(strtrim(dtype), 'k')
        n = zeros(size(wl));
        k = nums(:,2);
    else
        n = nums(:,2);
        k = zeros(size(wl));
    end
end


function [wl, n, k] = parseFormula(blk, dtype, defaultNw)
%PARSEFORMULA  Evaluate dispersion formulas 1-9.

    fnum = sscanf(dtype, 'formula %d');
    if isempty(fnum), wl=[]; n=[]; k=[]; return; end

    rangeStr = extractScalarField(blk, 'wavelength_range');
    if isempty(rangeStr)
        rangeStr = extractScalarField(blk, 'range');
    end
    if ~isempty(rangeStr)
        rng = sscanf(rangeStr, '%f %f');
        if numel(rng) >= 2
            wl = linspace(rng(1), rng(2), defaultNw).';
        else
            wl = linspace(0.2, 2.0, defaultNw).';
        end
    else
        wl = linspace(0.2, 2.0, defaultNw).';
    end

    coeffStr = extractScalarField(blk, 'coefficients');
    C = sscanf(coeffStr, '%f')';
    if isempty(C), wl=[]; n=[]; k=[]; return; end

    n = evalDispersion(fnum, C, wl);
    k = zeros(size(wl));
end


function n = evalDispersion(fnum, C, wl)
%EVALDISPERSION  Compute n(lambda) for formula types 1-9.

    lam  = wl(:);
    lam2 = lam.^2;
    Nc   = numel(C);

    switch fnum

        case 1   % Sellmeier (modified form)
            nsq = 1 + C(1);
            for j = 2:2:Nc-1
                nsq = nsq + C(j) .* lam2 ./ (lam2 - C(j+1).^2);
            end
            n = sqrt(abs(nsq));

        case 2   % Sellmeier-2
            nsq = 1 + C(1);
            for j = 2:2:Nc-1
                nsq = nsq + C(j) .* lam2 ./ (lam2 - C(j+1));
            end
            n = sqrt(abs(nsq));

        case 3   % Polynomial
            nsq = C(1);
            for j = 2:2:Nc-1
                nsq = nsq + C(j) .* lam.^C(j+1);
            end
            n = sqrt(abs(nsq));

        case 4   % RefractiveIndex.INFO
            nsq = C(1);
            idx = 2;
            while idx + 3 <= Nc
                nsq = nsq + C(idx) .* lam.^C(idx+1) ./ ...
                            (lam2 - C(idx+2).^C(idx+3));
                idx = idx + 4;
                if idx > 9, break; end
            end
            while idx + 1 <= Nc
                nsq = nsq + C(idx) .* lam.^C(idx+1);
                idx = idx + 2;
            end
            n = sqrt(abs(nsq));

        case 5   % Cauchy
            n = C(1) * ones(size(lam));
            for j = 2:2:Nc-1
                n = n + C(j) .* lam.^C(j+1);
            end

        case 6   % Gases
            n = 1 + C(1);
            for j = 2:2:Nc-1
                n = n + C(j) ./ (C(j+1) - lam.^(-2));
            end

        case 7   % Herzberger
            L = 1 ./ (lam2 - 0.028);
            Cp = [C, zeros(1, max(0, 6-Nc))];
            n = Cp(1) + Cp(2).*L + Cp(3).*L.^2 + ...
                Cp(4).*lam2 + Cp(5).*lam.^4 + Cp(6).*lam.^6;

        case 8   % Retro
            Cp = [C, zeros(1, max(0, 4-Nc))];
            a = Cp(1) + Cp(2).*lam2./(lam2 - Cp(3)) + Cp(4).*lam2;
            n = sqrt(abs( (2.*a + 1) ./ (1 - a) ));

        case 9   % Exotic
            Cp = [C, zeros(1, max(0, 6-Nc))];
            n = Cp(1) + Cp(2)./(lam2 - Cp(3)) + ...
                Cp(4).*(lam - Cp(5)) ./ ((lam - Cp(5)).^2 + Cp(6));

        otherwise
            warning('Formula %d not implemented. Returning n = 1.', fnum);
            n = ones(size(lam));
    end

    n = real(n);
end


function temperature_K = extractTemperatureFromName(sourceName)
%EXTRACTTEMPERATUREFROMNAME  Parse a temperature in Kelvin from a source
%  filename such as "Ferrera-298K" or "Author-20C". Returns NaN if no
%  temperature is encoded in the name.

    temperature_K = NaN;

    sourceName = char(sourceName);

    tok = regexp(sourceName, '(\d+(?:\.\d+)?)\s*K', 'tokens', 'once');
    if ~isempty(tok)
        temperature_K = str2double(tok{1});
        return
    end

    tok = regexp(sourceName, '(\d+(?:\.\d+)?)\s*C', 'tokens', 'once');
    if ~isempty(tok)
        temperature_C = str2double(tok{1});
        temperature_K = temperature_C + 273.15;
        return
    end

end


function sources = createTemperatureDependentSources(sources)
%CREATETEMPERATUREDEPENDENTSOURCES
%
% Finds sources belonging to the same author but measured at different
% temperatures and creates one additional combined source called
% "Author-Tdepend" (e.g. Ferrera-298K, Ferrera-404K, ... -> also adds
% Ferrera-Tdepend), with n(lambda,T) and k(lambda,T) stacked along a new
% temperature dimension. Sources whose name encodes a thickness rather
% than a temperature (e.g. "Ferrera-10nm") are ignored for this grouping.

    if isempty(sources)
        return
    end

    temperature_K = [sources.temperature_K];
    hasTemperature = isfinite(temperature_K);

    if ~any(hasTemperature)
        return
    end

    tempSources = sources(hasTemperature);

    % Determine author name by stripping the temperature suffix, e.g.
    % Ferrera-298K -> Ferrera.
    authorNames = cell(size(tempSources));
    for i = 1:numel(tempSources)
        name = tempSources(i).name;
        name = regexprep(name, '[-_ ]?\d+(?:\.\d+)?\s*[KC]$', '');
        authorNames{i} = name;
    end

    uniqueAuthors = unique(authorNames);

    newSources = sources([]);
    nNew = 0;

    for a = 1:numel(uniqueAuthors)

        author = uniqueAuthors{a};
        idx = strcmp(authorNames, author);
        groupedSources = tempSources(idx);

        % Need at least two temperatures to form a T-dependent entry.
        if numel(groupedSources) < 2
            continue
        end

        [~, order] = sort([groupedSources.temperature_K]);
        groupedSources = groupedSources(order);

        % Wavelength grids must match exactly across temperatures.
        wl1 = groupedSources(1).wavelength(:);
        compatible = true;
        for j = 2:numel(groupedSources)
            wlj = groupedSources(j).wavelength(:);
            if numel(wlj) ~= numel(wl1) || any(abs(wlj - wl1) > 1e-12)
                compatible = false;
                break
            end
        end

        if ~compatible
            warning(['Cannot create %s-Tdepend because the ', ...
                     'temperature datasets do not have identical ', ...
                     'wavelength grids.'], author);
            continue
        end

        nNew = nNew + 1;
        Tsrc = groupedSources(1);
        Tsrc.name = [author '-Tdepend'];
        Tsrc.temperature_K = [groupedSources.temperature_K];
        Tsrc.wavelength = wl1;

        Nw = numel(wl1);
        Nt = numel(groupedSources);
        n_T = zeros(Nw, Nt);
        k_T = zeros(Nw, Nt);

        for j = 1:Nt
            N_complex = groupedSources(j).n_complex(:);
            n_T(:,j) = real(N_complex);
            k_T(:,j) = imag(N_complex);
        end

        Tsrc.n_complex = complex(n_T, k_T);
        Tsrc.temperature = 'temperature dependent';
        Tsrc.film_thickness_um = NaN;

        newSources(nNew) = Tsrc;

    end

    if ~isempty(newSources)
        sources = [sources newSources];
    end

end