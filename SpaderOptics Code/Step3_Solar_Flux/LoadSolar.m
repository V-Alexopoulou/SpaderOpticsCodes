%load as
%% ==== Load Solar Wavelenght/Irradiance ==============================
% [I_sun,wl_sun,wl_sun_min,wl_sun_max,AM0_total] = LoadSolar();
function [I_sun,wl_sun,wl_sun_min,wl_sun_max,AM0_total] = LoadSolar()
    %% ==== USER SETTINGS ======================================================
    solarFile = fullfile(pwd, 'solar_spectrum.xlsx');
    solarSheet = 'NewAM0';
    % ==========================================================================
    
    
    %% ==== LOAD THE SOLAR SPECTRUM ============================================
    fprintf('Reading solar spectrum: %s\n', solarFile);
    raw = readcell(solarFile, 'Sheet', solarSheet);
    % Columns A and B; skip header rows until we hit numeric data
    wl_sun  = [];
    I_sun   = [];
    for r = 1:size(raw, 1)
        a = raw{r, 1};
        b = raw{r, 2};
        if isnumeric(a) && isnumeric(b) && ~isnan(a) && ~isnan(b)
            wl_sun(end+1, 1) = a;    %#ok<SAGROW>   wavelength  (µm)
            I_sun(end+1, 1)  = b;    %#ok<SAGROW>   irradiance  (W/m²/µm)
        end
    end
    
    if isempty(wl_sun)
        error('Could not extract numeric (wavelength, irradiance) pairs from %s.', solarFile);
    end
    %clear solarSheet solarFile;
    %clear r a b raw;
    % Enforce monotonic ascending wavelength
    [wl_sun, sortIdx] = sort(wl_sun);
    I_sun = I_sun(sortIdx);
    
    wl_sun_min = wl_sun(1);
    wl_sun_max = wl_sun(end);
    fprintf('   Solar spectrum: %d points,  λ = [%.4f – %.4f] µm\n', ...
        numel(wl_sun), wl_sun_min, wl_sun_max);
    
    % Total solar constant (informational)
    AM0_total = trapz(wl_sun, I_sun);
    fprintf('   Integrated solar irradiance over full range: %.2f W/m²\n', AM0_total);
   % clear  sortIdx AM0_total

   tiledlayout(2,1);

%% Full spectrum
nexttile
plot(wl_sun,I_sun,'LineWidth',1.5,'Color','k');
xscale log;
yscale log;
xlabel('Wavelength (\mum)');
ylabel({'Solar spectral irradiance', ...
        '(Wm^{-2}\mum^{-1})'});
title('a) Full ASTM E-490 spectrum');
grid on;

%% Zoom where most energy lies
nexttile
plot(wl_sun,I_sun,'LineWidth',1.5,'Color','k');
xlim([0 5]);
xlabel('Wavelength (\mum)');
ylabel({'Solar spectral irradiance', ...
        '(Wm^{-2}\mum^{-1})'});
title('b) 0-5 \mum region');
grid on;

sgtitle('Extraterrestrial solar spectrum (ASTM E-490)');

%% ==== EXPORT FIGURE ==========================================
%exportgraphics(gcf, 'SolarSpectrum.pdf', 'ContentType','vector');

end

  
   
   
   
   
  
