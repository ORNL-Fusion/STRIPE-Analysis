clc; clear; close all;

%% ----------------------------
% CONFIG
%% ----------------------------
tilts = [0 45 85];  % <<< set your three tilts here

charge_index = 1;          % 1 = sum over all Ta charge states
SKIP_FIRST_NSURFACES = 3;  % skip first 3 surfaces

E_MAX_BIN = 100;           % plot up to this energy bin
CLIM_MIN  = 1e8;

scale = 2.331579550458116e+09;

% Storage
spect_inc  = {};  % each cell: 1 x nE (angle-averaged incident)
spect_refl = {};  % each cell: 1 x nE (angle-averaged reflected)
nE_list    = zeros(size(tilts));

%% ----------------------------
% LOOP OVER TILTS
%% ----------------------------
for it = 1:numel(tilts)
    tilt = tilts(it);
    file = sprintf('../tilted_targets/%d_degrees/MPEX_runs/gitrm-surface.nc', tilt);

    % --- Read ---
    surfEDist    = ncread(file,'surfEDist');
    surfReflDist = ncread(file,'surfReflDist');

    % Expected dims: (nAngles, nEnergies, nSurfaces, nSpecies)
    surfEDist    = surfEDist(:,:,:,charge_index);
    surfReflDist = surfReflDist(:,:,:,charge_index);

    % --- Remove first N surfaces ---
    surfEDist    = surfEDist(:,:,SKIP_FIRST_NSURFACES+1:end);
    surfReflDist = surfReflDist(:,:,SKIP_FIRST_NSURFACES+1:end);

    % --- Sum over surfaces ---
    surfEDist    = sum(surfEDist,3);
    surfReflDist = sum(surfReflDist,3);

    % --- Squeeze to 2D (Angle x Energy) ---
    surfEDist    = squeeze(surfEDist);
    surfReflDist = squeeze(surfReflDist);

    % --- Scale + clean ---
    surfEDist    = scale .* surfEDist;
    surfReflDist = scale .* surfReflDist;

    surfEDist(surfEDist<=0)       = NaN;
    surfReflDist(surfReflDist<=0) = NaN;

    % --- Angle-averaged spectra vs energy (mean over angle dimension) ---
    % surf*(nAngles x nEnergies) -> nanmean over dim=1 => (1 x nEnergies)
    inc_E  = nanmean(surfEDist,    1);
    refl_E = nanmean(surfReflDist, 1);

    spect_inc{it}  = inc_E;
    spect_refl{it} = refl_E;
    nE_list(it)    = numel(inc_E);
end

%% ----------------------------
% PLOT: 1D spectra vs energy for all tilts
%% ----------------------------
% Use common energy range across tilts
nE_common = min(nE_list);
E_MAX_BIN = min(E_MAX_BIN, nE_common);
E = 1:nE_common;

figure('Color','w','Position',[200 120 900 700]);

% --- Incident panel ---
subplot(2,1,1); hold on; box on; grid on;
for it = 1:numel(tilts)
    y = spect_inc{it}(1:nE_common);
    semilogy(E(1:E_MAX_BIN), y(1:E_MAX_BIN), 'LineWidth', 2);
end
xlabel('Energy bin');
ylabel('\langle surfEDist \rangle_{angle}');
title('\Sigma Ta_x^{+}: Angle-averaged incident spectrum vs energy');
legend(arrayfun(@(t)sprintf('%d°',t), tilts, 'UniformOutput', false), ...
       'Location','northeast');
ylim([CLIM_MIN, max(cellfun(@(v) max(v(1:E_MAX_BIN),[],'omitnan'), spect_inc))]);

% --- Reflected panel ---
subplot(2,1,2); hold on; box on; grid on;
for it = 1:numel(tilts)
    y = spect_refl{it}(1:nE_common);
    semilogy(E(1:E_MAX_BIN), y(1:E_MAX_BIN), 'LineWidth', 2);
end
xlabel('Energy bin');
ylabel('\langle surfReflDist \rangle_{angle}');
title('\Sigma Ta_x^{+}: Angle-averaged reflected spectrum vs energy');
legend(arrayfun(@(t)sprintf('%d°',t), tilts, 'UniformOutput', false), ...
       'Location','northeast');
ylim([CLIM_MIN, max(cellfun(@(v) max(v(1:E_MAX_BIN),[],'omitnan'), spect_refl))]);