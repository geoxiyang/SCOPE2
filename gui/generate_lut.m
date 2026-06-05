%% generate_lut.m
% Pre-compute SCOPE2 sensitivity lookup table for the Streamlit GUI.
%
% Run ONCE from the SCOPE2 root directory:
%   >> generate_lut
%
% Output:  gui/lut_data/sensitivity_<param>.csv   (Actot, aPAR, Rd per run)
%          gui/lut_data/reflectance_<param>.csv    (reflectance spectrum per run)
%          gui/lut_data/wavelengths.csv            (wavelength axis, nm)
%
% Settings: 7 sweep parameters x 3 layer profiles x 40 values = 840 SCOPE runs.
% Expected runtime: 15-30 minutes.

clear all %#ok<CLALL>
restoredefaultpath

% Navigate to SCOPE2 root regardless of where this script was called from.
% This script lives in <root>/gui/, so the root is one level up.
script_dir = fileparts(mfilename('fullpath'));
root_dir   = fileparts(script_dir);
cd(root_dir);

addpath src/RTMs src/supporting src/fluxes src/IO

%% 1. Constants and spectral bands
constants = define_constants;
spectral  = define_bands;
path_input = 'input/';

%% 2. Read parameter files (mirrors SCOPE.m startup)
fid = fopen('set_parameter_filenames.csv','r');
parameter_file = textscan(fid,'%s','Delimiter',',');
fclose(fid);

% Options
fid = fopen([path_input parameter_file{1}{1}],'r');
Ni  = textscan(fid,'%d%s','Delimiter',',');
fclose(fid);
N   = double(Ni{1});

options.lite               = N(1);
options.calc_fluor         = 0;   % disabled for speed
options.calc_planck        = 0;
options.calc_xanthophyllabs = 0;
options.soilspectrum       = N(5);
options.Fluorescence_model = N(6);
options.apply_T_corr       = N(7);
options.verify             = 0;
options.saveCSV            = 0;
options.mSCOPE             = 1;   % always use mSCOPE
options.simulation         = 0;
options.calc_directional   = 0;
options.calc_vert_profiles = 0;
options.soil_heat_method   = N(14);
options.calc_rss_rbs       = N(15);
options.Cca_function_of_Cab = 0;
options.calc_zo            = 0;

% File names
f_names = {'Simulation_Name','soil_file','optipar_file','atmos_file','Dataset_dir',...
    'meteo_ec_csv','vegetation_retrieved_csv','LIDF_file','verification_dir',...
    'mSCOPE_csv','nly'};
cols = {'t','year','Rin','Rli','p','Ta','ea','u','RH','VPD','tts','tto','psi',...
    'Cab','Cca','Cdm','Cw','Cs','Cant','N','SMC','BSMBrightness','BSMlat','BSMlon',...
    'LAI','hc','LIDFa','LIDFb','z','Ca','Vcmo','m','atmos_names'};
fnc = [f_names, cols];
F   = struct('FileID', fnc);

fid = fopen([path_input parameter_file{1}{2}],'r');
while ~feof(fid)
    line = fgetl(fid);
    if ~isempty(line) && ~(line(1) == '%')
        X = textscan(line,'%s%s','Delimiter',',','Whitespace','\t');
        if ~isempty(X{1}) && ~isempty(X{2})
            k = find(strcmp(fnc, X{1}{1}));
            if ~isempty(k), F(k).FileName = X{2}{1}; end
        end
    end
end
fclose(fid);

% Input data
k = 1;
fid = fopen([path_input parameter_file{1}{3}],'r');
clear X varnames
while ~feof(fid)
    line = fgetl(fid);
    if ischar(line) && ~isempty(line)
        y = textscan(line,'%s','Delimiter',',','TreatAsEmpty',' ');
        varnames(k) = string(y{1}{1}); %#ok<SAGROW>
        X(k).Val    = str2double(y{:}); %#ok<SAGROW>
        k = k+1;
    end
end
fclose(fid);

V = assignvarnames();
for i = 1:length(V)
    j = find(strcmp(varnames, V(i).Name));
    if isempty(j)
        if i == 2, options.Cca_function_of_Cab = 1; end
    else
        kk = find(~isnan(X(j).Val));
        if ~isempty(kk), V(i).Val = X(j).Val(kk); else V(i).Val = -999; end
    end
end

%% 3. Load spectral data (done once)
load([path_input,'fluspect_parameters/', F(3).FileName]);
rsfile = load([path_input,'soil_spectra/', F(2).FileName]);
atmo   = load_atmo(fullfile(path_input,'radiationdata', F(4).FileName), spectral.SCOPEspec);

%% 4. Canopy structure constants
canopy.nlincl  = 13;
canopy.nlazi   = 36;
canopy.litab   = [5:10:75 81:2:89]';
canopy.lazitab = (5:10:355);

%% 5. Default parameter structs
vi  = ones(length(V),1);
xyt_def.t = 0; xyt_def.year = 0;
xyt_def.startDOY = 20060618; xyt_def.endDOY = 20300101;
xyt_def.LAT = 51.55; xyt_def.LON = 5.55; xyt_def.timezn = 1;
soil_dummy = struct;
[soil0, leafbio0, canopy0, meteo0, angles0, ~] = ...
    select_input(V, vi, canopy, options, constants, xyt_def, soil_dummy);

%% 6. Layer profile weights (nly = 10, top-to-bottom)
nly = 10;
idx = 1:nly;
profiles = {'Uniform','Top-heavy','Bottom-heavy'};
weights  = zeros(3, nly);
weights(1,:) = ones(1,nly) / nly;                                    % Uniform
w = exp(-0.3*(idx-1));  weights(2,:) = w / sum(w);                   % Top-heavy
w = exp(-0.3*(nly-idx)); weights(3,:) = w / sum(w);                  % Bottom-heavy

%% 7. Sweep definitions
sweeps.Rin    = linspace(0.1, 1200, 40);
sweeps.Ta     = linspace(-5,  45,   40);
sweeps.LAI    = linspace(0.1, 8.0,  40);
sweeps.Vcmo   = linspace(5,   200,  40);
sweeps.Ca     = linspace(200, 800,  40);
sweeps.ea     = linspace(5,   40,   40);
sweeps.LIDFa  = linspace(-0.5, 0.5, 40);
sweep_names   = fieldnames(sweeps);

%% 8. Wavelength decimation: every 5 nm over optical range (400-2400 nm)
wl_all   = spectral.wlP;         % 400:1:2400  (2001 values)
wl_idx   = 1:5:length(wl_all);   % every 5th → 401 wavelengths
wl_decim = wl_all(wl_idx);

mkdir('gui/lut_data');
writematrix(wl_decim, 'gui/lut_data/wavelengths.csv');
n_wl = length(wl_decim);

%% 9. Main loop
fprintf('\nStarting LUT generation: %d params x %d profiles x 40 points = %d runs\n', ...
    length(sweep_names), length(profiles), length(sweep_names)*length(profiles)*40);

for s = 1:length(sweep_names)
    param  = sweep_names{s};
    values = sweeps.(param);
    n      = length(values);
    total  = n * length(profiles);

    all_profile = repmat("", total, 1);
    all_value   = zeros(total, 1);
    all_Actot   = zeros(total, 1);
    all_aPAR    = zeros(total, 1);
    all_Rd      = zeros(total, 1);
    all_refl    = zeros(total, n_wl);

    row = 0;
    for p = 1:length(profiles)
        fprintf('  %s | %s ...\n', param, profiles{p});
        for i = 1:n
            row = row + 1;

            % Copy defaults and override swept parameter
            leafbio  = leafbio0;
            meteo    = meteo0;
            canopy_r = canopy0;
            soil_r   = soil0;
            angles_r = angles0;

            switch param
                case 'Rin',   meteo.Rin    = values(i);
                case 'Ta',    meteo.Ta     = values(i);
                case 'LAI',   canopy_r.LAI = max(1e-9, values(i));
                case 'Vcmo',  leafbio.Vcmo = values(i);
                case 'Ca',    meteo.Ca     = values(i);
                case 'ea',    meteo.ea     = values(i);
                case 'LIDFa', canopy_r.LIDFa = values(i);
            end

            % Update derived quantities that depend on LAI / canopy geometry
            [canopy_r.zo, canopy_r.d] = zo_and_d(soil_r, canopy_r, constants);

            % Build mly struct for this profile
            mly.nly    = nly;
            mly.totLAI = canopy_r.LAI;
            mly.pLAI   = canopy_r.LAI * weights(p,:);
            mly.pCab   = leafbio.Cab  * ones(1,nly);
            mly.pCca   = leafbio.Cca  * ones(1,nly);
            mly.pCdm   = leafbio.Cdm  * ones(1,nly);
            mly.pCw    = leafbio.Cw   * ones(1,nly);
            mly.pCs    = leafbio.Cs   * ones(1,nly);
            mly.pN     = leafbio.N    * ones(1,nly);

            % RTM sublayer grid: must be large enough that every mSCOPE layer
            % gets >= 1 sublayer even for the most skewed weight profile.
            % For nly=10, Bottom/Top-heavy min weight ≈ 0.018 → need nl >= 55.
            % Use nly*7 = 70 as a floating-point-safe minimum (avoids 2/min_w rounding).
            nl_rtm  = max(ceil(10*canopy_r.LAI/canopy_r.Cv), nly * 7);
            nl_rtm  = max(nl_rtm, 2);

            canopy_r.nlayers = nl_rtm;
            x_grid           = (-1/nl_rtm : -1/nl_rtm : -1)';
            canopy_r.xl      = [0; x_grid];
            canopy_r.lidf    = leafangles(canopy_r.LIDFa, canopy_r.LIDFb);

            % Soil reflectance
            soil_r.refl               = rsfile(:, soil_r.spectrum+1);
            soil_r.refl(spectral.IwlT) = soil_r.rs_thermal;
            soil_r.Tsold              = meteo.Ta * ones(12,2);

            % Leaf optical properties
            leafbio.emis = 1 - leafbio.rho_thermal - leafbio.tau_thermal;
            leafbio.V2Z  = 0;
            leafopt      = fluspect_mSCOPE(mly, spectral, leafbio, optipar, nl_rtm);
            leafopt.refl(:, spectral.IwlT) = leafbio.rho_thermal;
            leafopt.tran(:, spectral.IwlT) = leafbio.tau_thermal;

            % fluspect_mSCOPE uses floor(cumsum(...)*nl_rtm) internally; floating-point
            % rounding on the normalised weights can make the last index land at nl_rtm-1
            % instead of nl_rtm.  Sync canopy.nlayers to the actual output size so RTMo
            % never receives a mismatched nl.
            actual_nl = size(leafopt.refl, 1);
            if actual_nl ~= nl_rtm
                canopy_r.nlayers = actual_nl;
                x_grid_a         = (-1/actual_nl : -1/actual_nl : -1)';
                canopy_r.xl      = [0; x_grid_a];
            end

            % Radiative transfer
            [rad, gap, ~] = RTMo(spectral, atmo, soil_r, leafopt, canopy_r, angles_r, ...
                                  constants, meteo, options);

            % Energy balance
            xyt_r.t = 0; xyt_r.year = 0;
            [~, rad, ~, soil_r, bcu, bch, ~] = ebal(constants, options, rad, gap, ...
                meteo, soil_r, canopy_r, leafbio, 1, xyt_r);

            % Canopy-level aggregation (same as SCOPE.m lines 365-384)
            nl   = canopy_r.nlayers;
            Ps   = gap.Ps(1:nl);
            Ph   = 1 - Ps;
            if options.lite == 0
                integr = 'angles_and_layers';
            else
                integr = 'layers';
            end

            Actot = canopy_r.LAI * (meanleaf(canopy_r, bch.A,  'layers', Ph) + ...
                                     meanleaf(canopy_r, bcu.A,  integr,  Ps));
            aPAR  = canopy_r.LAI * (meanleaf(canopy_r, rad.Pnh,'layers', Ph) + ...
                                     meanleaf(canopy_r, rad.Pnu, integr, Ps));
            Rd    = canopy_r.LAI * (meanleaf(canopy_r, bch.Rd, 'layers', Ph) + ...
                                     meanleaf(canopy_r, bcu.Rd, integr,  Ps));
            irrad = rad.Esun_ + rad.Esky_;
            refl  = zeros(size(rad.Lo_));
            ok    = irrad > 0;
            refl(ok) = pi * rad.Lo_(ok) ./ irrad(ok);

            all_profile(row) = string(profiles{p});
            all_value(row)   = values(i);
            all_Actot(row)   = Actot;
            all_aPAR(row)    = aPAR;
            all_Rd(row)      = Rd;
            all_refl(row,:)  = refl(wl_idx)';
        end
    end

    % Save sensitivity CSV
    T_sens = table(all_profile, all_value, all_Actot, all_aPAR, all_Rd, ...
        'VariableNames', {'profile','value','Actot','aPAR','Rd'});
    writetable(T_sens, sprintf('gui/lut_data/sensitivity_%s.csv', param));

    % Save reflectance CSV
    wl_headers = arrayfun(@(w) sprintf('wl_%d', w), wl_decim', 'UniformOutput', false);
    T_refl = [table(all_profile, all_value, 'VariableNames', {'profile','value'}), ...
               array2table(all_refl, 'VariableNames', wl_headers)];
    writetable(T_refl, sprintf('gui/lut_data/reflectance_%s.csv', param));

    fprintf('    Saved sensitivity_%s.csv and reflectance_%s.csv\n', param, param);
end

fprintf('\nDone. Files written to gui/lut_data/\n');
