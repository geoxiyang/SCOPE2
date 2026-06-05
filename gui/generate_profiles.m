%% generate_profiles.m
% Pre-compute vertical profile data for the Streamlit GUI.
%
% Run ONCE from the SCOPE2 root directory:
%   >> generate_profiles
%
% Output:  gui/lut_data/vertical_profiles.csv
%
% Three leaf-angle scenarios (all with uniform LAI distribution):
%   Uniform          : spherical (LIDFa = 0) at every layer
%   Vertical-on-top  : LIDFa gradient from +0.4 (top) to -0.4 (bottom)
%   Horizontal-on-top: LIDFa gradient from -0.4 (top) to +0.4 (bottom)
%
% Sweeps: 3 scenarios x 5 total-LAI values + 3 scenarios x 5 Rin values = 30 runs.
% Expected runtime: 2-5 minutes.

clear all %#ok<CLALL>
restoredefaultpath

% Navigate to SCOPE2 root regardless of where this script was called from.
script_dir = fileparts(mfilename('fullpath'));
root_dir   = fileparts(script_dir);
cd(root_dir);

addpath src/RTMs src/supporting src/fluxes src/IO

%% 1. Constants and spectral bands
constants = define_constants;
spectral  = define_bands;
path_input = 'input/';

%% 2. Read parameter files (mirrors SCOPE.m / generate_lut.m startup)
fid = fopen('set_parameter_filenames.csv','r');
parameter_file = textscan(fid,'%s','Delimiter',',');
fclose(fid);

fid = fopen([path_input parameter_file{1}{1}],'r');
Ni  = textscan(fid,'%d%s','Delimiter',',');
fclose(fid);
N   = double(Ni{1});

options.lite               = N(1);
options.calc_fluor         = 0;
options.calc_planck        = 0;
options.calc_xanthophyllabs = 0;
options.soilspectrum       = N(5);
options.Fluorescence_model = N(6);
options.apply_T_corr       = N(7);
options.verify             = 0;
options.saveCSV            = 0;
options.mSCOPE             = 1;
options.simulation         = 0;
options.calc_directional   = 0;
options.calc_vert_profiles = 0;
options.soil_heat_method   = N(14);
options.calc_rss_rbs       = N(15);
options.Cca_function_of_Cab = 0;
options.calc_zo            = 0;

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

k2 = 1;
fid = fopen([path_input parameter_file{1}{3}],'r');
clear X varnames
while ~feof(fid)
    line = fgetl(fid);
    if ischar(line) && ~isempty(line)
        y = textscan(line,'%s','Delimiter',',','TreatAsEmpty',' ');
        varnames(k2) = string(y{1}{1}); %#ok<SAGROW>
        X(k2).Val    = str2double(y{:}); %#ok<SAGROW>
        k2 = k2+1;
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

%% 3. Load spectral data
load([path_input,'fluspect_parameters/', F(3).FileName]);
rsfile = load([path_input,'soil_spectra/', F(2).FileName]);
atmo   = load_atmo(fullfile(path_input,'radiationdata', F(4).FileName), spectral.SCOPEspec);

%% 4. Canopy structure
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

%% 6. Leaf-angle scenarios (nly=10 layers, top→bottom)
nly = 10;
idx_layers = 1:nly;
scenario_names = {'Uniform','Vertical-on-top','Horizontal-on-top'};
% Per-layer LIDFa: linspace gives smooth gradient top→bottom
lidf_a_profiles = {
    zeros(1, nly),                    % spherical throughout
    linspace(+0.4, -0.4, nly),        % erectophile→planophile
    linspace(-0.4, +0.4, nly),        % planophile→erectophile
};
% Uniform LAI distribution for all scenarios
weights_uniform = ones(1, nly) / nly;

%% 7. Sweep values
sweep_params = {'LAI', 'Rin'};
sweep_values = {[0.5, 1.5, 3.0, 5.0, 8.0], ...
                [100, 300, 600, 900, 1200]};

n_total = length(scenario_names) * (length(sweep_values{1}) + length(sweep_values{2}));
fprintf('\nStarting profile generation: %d runs total\n', n_total);

%% 8. Pre-allocate output table
out_scenario   = {};
out_sweep_param = {};
out_sweep_value = [];
out_depth_rel  = [];
out_A          = [];
out_aPAR       = [];
out_Rd         = [];
out_Ps         = [];
out_Kn         = [];

%% 9. Main loop
for sc = 1:length(scenario_names)
    lidf_a_vec = lidf_a_profiles{sc};   % [1 x nly] LIDFa per mSCOPE layer

    for sp = 1:length(sweep_params)
        param  = sweep_params{sp};
        values = sweep_values{sp};

        for vi_idx = 1:length(values)
            fprintf('  %s | %s = %.1f ...\n', scenario_names{sc}, param, values(vi_idx));

            % Copy defaults and override swept parameter
            leafbio  = leafbio0;
            meteo    = meteo0;
            canopy_r = canopy0;
            soil_r   = soil0;
            angles_r = angles0;

            switch param
                case 'Rin',  meteo.Rin    = values(vi_idx);
                case 'LAI',  canopy_r.LAI = max(1e-9, values(vi_idx));
            end

            [canopy_r.zo, canopy_r.d] = zo_and_d(soil_r, canopy_r, constants);

            % Build mly struct (uniform LAI distribution)
            mly.nly    = nly;
            mly.totLAI = canopy_r.LAI;
            mly.pLAI   = canopy_r.LAI * weights_uniform;
            mly.pCab   = leafbio.Cab  * ones(1,nly);
            mly.pCca   = leafbio.Cca  * ones(1,nly);
            mly.pCdm   = leafbio.Cdm  * ones(1,nly);
            mly.pCw    = leafbio.Cw   * ones(1,nly);
            mly.pCs    = leafbio.Cs   * ones(1,nly);
            mly.pN     = leafbio.N    * ones(1,nly);

            % RTM sublayer grid (floating-point-safe minimum)
            nl_rtm  = max(ceil(10*canopy_r.LAI/canopy_r.Cv), nly * 7);
            nl_rtm  = max(nl_rtm, 2);

            canopy_r.nlayers = nl_rtm;
            x_grid           = (-1/nl_rtm : -1/nl_rtm : -1)';
            canopy_r.xl      = [0; x_grid];
            canopy_r.lidf    = leafangles(canopy_r.LIDFa, canopy_r.LIDFb);  % temp scalar

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

            % Sync nlayers to actual fluspect output (floating-point guard)
            actual_nl = size(leafopt.refl, 1);
            if actual_nl ~= nl_rtm
                canopy_r.nlayers = actual_nl;
                x_grid_a         = (-1/actual_nl : -1/actual_nl : -1)';
                canopy_r.xl      = [0; x_grid_a];
            end
            nl = canopy_r.nlayers;

            % Build per-layer LIDF matrix [13 x nl] using mSCOPE→RTM layer mapping
            pLAI_norm   = mly.pLAI / sum(mly.pLAI);
            indStar_raw = floor(cumsum(pLAI_norm) * nl);
            indStar_raw(end) = nl;           % force last index to nl
            for ii = 2:nly                   % ensure strictly increasing
                if indStar_raw(ii) <= indStar_raw(ii-1)
                    indStar_raw(ii) = indStar_raw(ii-1) + 1;
                end
            end
            indStar = [1, min(indStar_raw, nl)];

            lidf_mat = zeros(13, nl);
            for jj = 1:nly
                j1 = indStar(jj);
                j2 = indStar(jj+1);
                lidf_j = leafangles(lidf_a_vec(jj), canopy_r.LIDFb);  % [13 x 1]
                lidf_mat(:, j1:j2) = repmat(lidf_j, 1, j2-j1+1);
            end
            canopy_r.lidf = lidf_mat;   % [13 x nl] — activates per-layer path in RTMo

            % Radiative transfer
            [rad, gap, ~] = RTMo(spectral, atmo, soil_r, leafopt, canopy_r, angles_r, ...
                                  constants, meteo, options);

            % Energy balance
            xyt_r.t = 0; xyt_r.year = 0;
            [~, rad, ~, soil_r, bcu, bch, ~] = ebal(constants, options, rad, gap, ...
                meteo, soil_r, canopy_r, leafbio, 1, xyt_r);

            % Per-layer metrics (sunlit+shaded weighted mean per unit leaf area)
            A_layer    = zeros(nl, 1);
            aPAR_layer = zeros(nl, 1);
            Rd_layer   = zeros(nl, 1);
            Kn_layer   = zeros(nl, 1);
            Ps_layer   = gap.Ps(1:nl);

            for j = 1:nl
                Ps_j = Ps_layer(j);
                Ph_j = 1 - Ps_j;
                if options.lite
                    A_layer(j)    = Ph_j*bch.A(j)   + Ps_j*bcu.A(j);
                    Rd_layer(j)   = Ph_j*bch.Rd(j)  + Ps_j*bcu.Rd(j);
                    Kn_layer(j)   = Ph_j*bch.Kn(j)  + Ps_j*bcu.Kn(j);
                    aPAR_layer(j) = Ph_j*rad.Pnh(j)  + Ps_j*rad.Pnu(j);
                else
                    A_layer(j)    = Ph_j*bch.A(j)   + Ps_j*mean(bcu.A(:,:,j),   [1 2]);
                    Rd_layer(j)   = Ph_j*bch.Rd(j)  + Ps_j*mean(bcu.Rd(:,:,j),  [1 2]);
                    Kn_layer(j)   = Ph_j*bch.Kn(j)  + Ps_j*mean(bcu.Kn(:,:,j),  [1 2]);
                    aPAR_layer(j) = Ph_j*rad.Pnh(j)  + Ps_j*mean(rad.Pnu(:,:,j), [1 2]);
                end
            end

            % Relative depth: midpoint of each RTM sublayer [0=top, 1=bottom]
            depth_rel = ((1:nl)' - 0.5) / nl;

            % Accumulate rows
            n_new = nl;
            out_scenario   = [out_scenario;   repmat({scenario_names{sc}}, n_new, 1)]; %#ok<AGROW>
            out_sweep_param = [out_sweep_param; repmat({param}, n_new, 1)]; %#ok<AGROW>
            out_sweep_value = [out_sweep_value; repmat(values(vi_idx), n_new, 1)]; %#ok<AGROW>
            out_depth_rel   = [out_depth_rel;  depth_rel]; %#ok<AGROW>
            out_A           = [out_A;          A_layer];   %#ok<AGROW>
            out_aPAR        = [out_aPAR;       aPAR_layer]; %#ok<AGROW>
            out_Rd          = [out_Rd;         Rd_layer];   %#ok<AGROW>
            out_Ps          = [out_Ps;         Ps_layer];   %#ok<AGROW>
            out_Kn          = [out_Kn;         Kn_layer];   %#ok<AGROW>
        end
    end
end

%% 10. Save CSV
mkdir('gui/lut_data');
T = table(string(out_scenario), string(out_sweep_param), out_sweep_value, ...
          out_depth_rel, out_A, out_aPAR, out_Rd, out_Ps, out_Kn, ...
          'VariableNames', {'scenario','sweep_param','sweep_value','depth_rel',...
                            'A','aPAR','Rd','Ps','Kn'});
writetable(T, 'gui/lut_data/vertical_profiles.csv');
fprintf('\nDone. Saved %d rows to gui/lut_data/vertical_profiles.csv\n', height(T));
