%% ieads_plot_full_resample_energy_pcolor.m
% Single-file script: read energy axis, resample in energy-space, pcolor w/edges.
% Paste into MATLAB and run.

clc; clear; close all;

%% ----------------------------
% CONFIG
%% ----------------------------
tilt = 0;   % degrees
file = sprintf('../tilted_targets/test/%d_degrees/MPEX_runs/gitrm-surface.nc', tilt);

% surface selection
SKIP_FIRST_NSURFACES = 3;
DROP_LAST_NSURFACES  = 3;

% resample / coarsen
TARGET_N_E_BINS = 5000;      % set [] to skip resample (good for quick test)
COARSEN_ANGLE_FACTOR  = 2;
COARSEN_ENERGY_FACTOR = 20;  % energy block-avg after resample

% angle marker
TARGET_ANGLE_BIN = 180;

% scale (set to 1 to inspect raw iEADS if you expect ~0..60)
scale = 1.165799141267602e+16 ;

% charge bins
BIN_COMBINED = 1;
BIN_TA = 3:6;  % Ta1+..Ta4+
if numel(BIN_TA) ~= 4, error('BIN_TA must be exactly 4 bins (3:6).'); end

% ----------------------------
% REAL energy axis (physical) & plotting behavior
% ----------------------------
E_PHYS_MIN = 0;     % eV (fallback)
E_PHYS_MAX = 500;   % eV (fallback)

USE_ENERGY_ZOOM_EV = true;
E_ZOOM_MIN_EV = 0;
E_ZOOM_MAX_EV = 500;

FORCE_LOG_X_ON_ENERGY = true;   % log-x on REAL energy
E_LOG_MIN_POSITIVE = 1e-3;      % eV floor to avoid log(0)

% color limits behavior
USE_AUTO_CLIM_2D = true;
AUTO_CLIM_PCTS = [2 98];        % percentiles for CLIM
CLIM_2D_FIXED = [1e8 1e10];     % fallback if auto disabled

%% ----------------------------
% READ DATA
%% ----------------------------
surfEDist_all    = ncread(file,'surfEDist');     % [A,E,S,species,bins]
surfReflDist_all = ncread(file,'surfReflDist');  % [A,E,S,species]

[nA0, nE0, nSurf0, nSpecies0, nBins0] = size(surfEDist_all);
fprintf('surfEDist dims: A=%d, E=%d, S=%d, species=%d, bins=%d\n', nA0, nE0, nSurf0, nSpecies0, nBins0);
if nBins0 < 6, error('Need at least 6 bins for Ta1+..Ta4+.'); end

% try to read native energy grid (length = nE0 or nE0+1 edges)
[E_native, E_name] = try_read_energy_axis(file, nE0, E_PHYS_MIN, E_PHYS_MAX);
fprintf('Energy axis source: %s (n=%d, min=%.6g, max=%.6g)\n', E_name, numel(E_native), min(E_native), max(E_native));

%% ----------------------------
% PREP reflected species-summed (no charge dimension)
%% ----------------------------
refl_AES = squeeze(sum(surfReflDist_all,4)); % [A,E,S]
[refl_AE, ~, ~] = prep_AES_to_AE_energy(refl_AES, E_native, SKIP_FIRST_NSURFACES, DROP_LAST_NSURFACES, ...
                                        TARGET_N_E_BINS, COARSEN_ANGLE_FACTOR, COARSEN_ENERGY_FACTOR, scale);

%% ----------------------------
% PREP incident AE for combined + Ta bins
%% ----------------------------
bins_to_plot = [BIN_COMBINED, BIN_TA]; % [1,3,4,5,6]
nPlot = numel(bins_to_plot);

inc_AE_list = cell(nPlot,1);
nA_list = zeros(nPlot,1);
nE_list = zeros(nPlot,1);

for k = 1:nPlot
    b = bins_to_plot(k);
    if ~isscalar(b) || b~=floor(b) || b < 1 || b > nBins0
        error('Invalid bin index bins_to_plot(%d) = %s; valid range 1..%d', k, mat2str(b), nBins0);
    end

    tmp = surfEDist_all(:,:,:,:,b); % may be 4-D or squeezed
    expected4 = [nA0, nE0, nSurf0, nSpecies0];
    if ~isequal(size(tmp), expected4)
        if numel(tmp) == prod(expected4)
            tmp = reshape(tmp, expected4);
        else
            tmp = reshape(tmp, [size(tmp,1), size(tmp,2), 1, size(tmp,3)]);
        end
    end

    AES = squeeze(sum(tmp,4)); % [A,E,S]
    if ndims(AES) ~= 3
        error('After sum over species, AES must be 3-D. size(AES)=%s', mat2str(size(AES)));
    end

    [inc_AE, nA_k, nE_k] = prep_AES_to_AE_energy(AES, E_native, SKIP_FIRST_NSURFACES, DROP_LAST_NSURFACES, ...
                                                TARGET_N_E_BINS, COARSEN_ANGLE_FACTOR, COARSEN_ENERGY_FACTOR, scale);
    inc_AE_list{k} = inc_AE;
    nA_list(k) = nA_k;
    nE_list(k) = nE_k;
end

% unify dims across panels
nA_min = min(nA_list);
nE_min = min(nE_list);
for k = 1:nPlot
    inc_AE_list{k} = inc_AE_list{k}(1:nA_min, 1:nE_min);
end
refl_AE = refl_AE(1:nA_min, 1:nE_min);

%% ----------------------------
% BUILD REAL energy axis (E_plot) matching final energy dimension
%% ----------------------------
% After resample, we expect nE_min bins in the prepared arrays.
% Build E_new (centers) consistent with TARGET_N_E_BINS if used; otherwise use E_native mapping.
if ~isempty(TARGET_N_E_BINS)
    % The resample inside prep_AES_to_AE_energy used E_resample vector of length TARGET_N_E_BINS.
    % We must reconstruct the E_plot used there (it matches TARGET_N_E_BINS), then clip to nE_min.
    E_resample = build_resample_energy_axis(E_native, TARGET_N_E_BINS, E_PHYS_MIN, E_PHYS_MAX);
    E_plot_full = E_resample;
else
    % No resample: use E_native (centers). If E_native length > nE_min, take first nE_min
    E_plot_full = E_native(:).';
end

% Now clip to whatever nE_min we ended up with (prep AES may have coarsened)
if numel(E_plot_full) < nE_min
    warning('E_plot_full shorter than nE_min; truncating nE_min to match.');
    nE_min = numel(E_plot_full);
end
E_plot_full = E_plot_full(1:nE_min);
E_plot = E_plot_full;   % plotted centers
nE_plot = numel(E_plot);

angle_idx = min(max(1, TARGET_ANGLE_BIN), nA_min);

%% ----------------------------
% Build edges (for pcolor) consistent with centers
%% ----------------------------
E_edges = zeros(1, nE_plot+1);
if FORCE_LOG_X_ON_ENERGY
    for i = 1:(nE_plot-1)
        E_edges(i+1) = sqrt(E_plot(i) * E_plot(i+1));
    end
    if E_edges(2) <= 0
        E_edges(2) = E_plot(1)*1.0001;
    end
    E_edges(1) = max(E_plot(1)^2 / E_edges(2), E_LOG_MIN_POSITIVE);
    E_edges(nE_plot+1) = E_plot(end)^2 / E_edges(nE_plot);
else
    for i = 1:(nE_plot-1)
        E_edges(i+1) = 0.5*(E_plot(i)+E_plot(i+1));
    end
    dfirst = E_plot(2)-E_plot(1);
    dlast  = E_plot(end)-E_plot(end-1);
    E_edges(1) = E_plot(1) - 0.5*dfirst;
    E_edges(nE_plot+1) = E_plot(end) + 0.5*dlast;
end

% Angle edges
A_edges = 0:nA_min;

% meshgrid for pcolor
[Xe, Ye] = meshgrid(E_edges, A_edges);    % sizes (nA_min+1 x nE_plot+1)

%% ----------------------------
% DEBUG prints
%% ----------------------------
AE_test = inc_AE_list{1}(:, 1:nE_plot);
AE_f = AE_test(isfinite(AE_test));
fprintf('E_plot: min=%.6g eV, med=%.6g eV, max=%.6g eV (n=%d), XScale=%s\n', ...
    min(E_plot), median(E_plot), max(E_plot), numel(E_plot), ternary(FORCE_LOG_X_ON_ENERGY,'log','linear'));
if isempty(AE_f)
    fprintf('Combined AE: NO FINITE VALUES in selected range!\n');
else
    fprintf('Combined AE: finite=%d / %d, min=%.6g, max=%.6g\n', ...
        nnz(isfinite(AE_test)), numel(AE_test), min(AE_f), max(AE_f));
    fprintf('log10(Combined AE): min=%.3g, max=%.3g\n', log10(min(AE_f)), log10(max(AE_f)));
end

%% ----------------------------
% Choose CLIM for 2D panels (linear units)
%% ----------------------------
if USE_AUTO_CLIM_2D
    AE_all = inc_AE_list{1}(:);
    AE_all = AE_all(isfinite(AE_all));
    if isempty(AE_all)
        warning('Auto CLIM: no finite values found. Using a default range.');
        CLIM_2D = [1 10]; % fallback
    else
        p_lo = prctile(AE_all, AUTO_CLIM_PCTS(1));
        p_hi = prctile(AE_all, AUTO_CLIM_PCTS(2));
        p_lo = max(p_lo, 1e-20);
        p_hi = max(p_hi, p_lo*10);
        CLIM_2D = [p_lo, p_hi];
    end
else
    CLIM_2D = CLIM_2D_FIXED;
end
fprintf('Using CLIM_2D (linear units) = [%.3g, %.3g]\n', CLIM_2D(1), CLIM_2D(2));

%% ----------------------------
% PLOTTING (2x3) - inline pcolor panels
%% ----------------------------
figure('Color','w','Position',[80 80 1400 900]);
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
colormap(jet);
sgtitle(sprintf('MPEX target (tilt = %d^\\circ)', tilt), 'FontSize',16,'FontWeight','bold');

% panel 1: Combined
nexttile(1);
AE_plot = inc_AE_list{1}(:, 1:nE_plot);  % [nA_min x nE_plot]
Z = nan(size(Ye)); Z(1:size(AE_plot,1), 1:size(AE_plot,2)) = AE_plot;
pcolor(Xe, Ye, Z); shading flat;
set(gca,'XScale', ternary(FORCE_LOG_X_ON_ENERGY,'log','linear'), 'YDir','normal', 'ColorScale','log');
caxis(CLIM_2D); cb = colorbar; format_colorbar_decades(cb, CLIM_2D);
xlabel('Energy (eV)'); ylabel('Angle bin');
title(sprintf('Combined iEADS (bin=%d)', BIN_COMBINED));
xlim([min(E_plot) max(E_plot)]); ylim([1 nA_min]);
hold on; plot([min(E_plot) max(E_plot)], [angle_idx angle_idx], 'w--', 'LineWidth', 1.5);

% panel 2: Ta1+
nexttile(2);
AE_plot = inc_AE_list{2}(:, 1:nE_plot);
Z = nan(size(Ye)); Z(1:size(AE_plot,1), 1:size(AE_plot,2)) = AE_plot;
pcolor(Xe, Ye, Z); shading flat;
set(gca,'XScale', ternary(FORCE_LOG_X_ON_ENERGY,'log','linear'), 'YDir','normal', 'ColorScale','log');
caxis(CLIM_2D); cb = colorbar; format_colorbar_decades(cb, CLIM_2D);
xlabel('Energy (eV)'); ylabel('Angle bin'); title('Ta1^+');
xlim([min(E_plot) max(E_plot)]); ylim([1 nA_min]);
hold on; plot([min(E_plot) max(E_plot)], [angle_idx angle_idx], 'w--', 'LineWidth', 1.2);

% panel 3: Ta2+
nexttile(3);
AE_plot = inc_AE_list{3}(:, 1:nE_plot);
Z = nan(size(Ye)); Z(1:size(AE_plot,1), 1:size(AE_plot,2)) = AE_plot;
pcolor(Xe, Ye, Z); shading flat;
set(gca,'XScale', ternary(FORCE_LOG_X_ON_ENERGY,'log','linear'), 'YDir','normal', 'ColorScale','log');
caxis(CLIM_2D); cb = colorbar; format_colorbar_decades(cb, CLIM_2D);
xlabel('Energy (eV)'); ylabel('Angle bin'); title('Ta2^+');
xlim([min(E_plot) max(E_plot)]); ylim([1 nA_min]);
hold on; plot([min(E_plot) max(E_plot)], [angle_idx angle_idx], 'w--', 'LineWidth', 1.2);

% panel 4: 1D spectra (angle-averaged)
nexttile(4); hold on; box on;
spec_comb = nanmean(inc_AE_list{1}(:, 1:nE_plot), 1);
plot(E_plot, spec_comb, 'k', 'LineWidth', 2.5);
ta_labels = {'Ta1^+','Ta2^+','Ta3^+','Ta4^+'};
for j = 1:4
    specj = nanmean(inc_AE_list{1+j}(:, 1:nE_plot), 1);
    plot(E_plot, specj, 'LineWidth', 1.6);
end
grid on; xlabel('Energy (eV)'); ylabel('\langle iEADS \rangle_\theta (linear units)');
title('Angle-averaged iEADS'); set(gca,'XScale', ternary(FORCE_LOG_X_ON_ENERGY,'log','linear'));
legend(['Combined', ta_labels], 'Location','northeast'); xlim([min(E_plot) max(E_plot)]);

% panel 5: Ta3+
nexttile(5);
AE_plot = inc_AE_list{4}(:, 1:nE_plot);
Z = nan(size(Ye)); Z(1:size(AE_plot,1), 1:size(AE_plot,2)) = AE_plot;
pcolor(Xe, Ye, Z); shading flat;
set(gca,'XScale', ternary(FORCE_LOG_X_ON_ENERGY,'log','linear'), 'YDir','normal', 'ColorScale','log');
caxis(CLIM_2D); cb = colorbar; format_colorbar_decades(cb, CLIM_2D);
xlabel('Energy (eV)'); ylabel('Angle bin'); title('Ta3^+');
xlim([min(E_plot) max(E_plot)]); ylim([1 nA_min]);
hold on; plot([min(E_plot) max(E_plot)], [angle_idx angle_idx], 'w--', 'LineWidth', 1.2);

% panel 6: Ta4+
nexttile(6);
AE_plot = inc_AE_list{5}(:, 1:nE_plot);
Z = nan(size(Ye)); Z(1:size(AE_plot,1), 1:size(AE_plot,2)) = AE_plot;
pcolor(Xe, Ye, Z); shading flat;
set(gca,'XScale', ternary(FORCE_LOG_X_ON_ENERGY,'log','linear'), 'YDir','normal', 'ColorScale','log');
caxis(CLIM_2D); cb = colorbar; format_colorbar_decades(cb, CLIM_2D);
xlabel('Energy (eV)'); ylabel('Angle bin'); title('Ta4^+');
xlim([min(E_plot) max(E_plot)]); ylim([1 nA_min]);
hold on; plot([min(E_plot) max(E_plot)], [angle_idx angle_idx], 'w--', 'LineWidth', 1.2);

%% ----------------------------
% Reflected spectrum (REAL energy axis + optional log-x)
%% ----------------------------
figure('Color','w','Position',[240 240 1000 460]); hold on; box on;
refl_spec = nanmean(refl_AE(:, 1:nE_plot), 1);
plot(E_plot, refl_spec, 'LineWidth', 2.2);
grid on; xlabel('Energy (eV)'); ylabel('\langle Flux \rangle_{\theta} (linear units)');
title('Reflected spectrum (species-summed; NOT charge-resolved)');
set(gca,'XScale', ternary(FORCE_LOG_X_ON_ENERGY,'log','linear'));
xlim([min(E_plot) max(E_plot)]);

fprintf('Done.\n');

%% ========================================================================
% Helper: prep_AES_to_AE_energy
% - inputs: AESmat [nA x nE_old x nSurf]
% - E_old: native energy centers (length nE_old)
% - optionally resamples to TARGET_N_E_BINS by energy interpolation
% - returns AEmat [nA' x nE_out]
%% ========================================================================
function [AEmat, nA_out, nE_out] = prep_AES_to_AE_energy( ...
    AESmat, E_old, SKIP_FIRST_NSURFACES, DROP_LAST_NSURFACES, TARGET_N_E_BINS, ...
    COARSEN_ANGLE_FACTOR, COARSEN_ENERGY_FACTOR, scale)

    % surfaces select
    nSurf = size(AESmat,3);
    i0 = SKIP_FIRST_NSURFACES + 1;
    i1 = nSurf - DROP_LAST_NSURFACES;
    if i0 > i1
        error('Surface selection invalid: i0=%d, i1=%d. Check SKIP_FIRST_NSURFACES/DROP_LAST_NSURFACES.', i0, i1);
    end
    AESsel = AESmat(:,:,i0:i1);

    % sum surfaces -> [A,E_old]
    AE = squeeze(sum(AESsel,3));
    [nA, nE_old] = size(AE);

    % build target energy centers if resampling requested
    if ~isempty(TARGET_N_E_BINS) && TARGET_N_E_BINS ~= nE_old
        E_new = build_resample_energy_axis(E_old, TARGET_N_E_BINS, min(E_old), max(E_old));
        AEi = nan(nA, numel(E_new));
        for ia = 1:nA
            y = AE(ia,:);
            if all(isnan(y))
                AEi(ia,:) = NaN;
            else
                mask = ~isnan(y);
                % energy-based interpolation (physical)
                yi = interp1(E_old(mask), y(mask), E_new, 'linear', NaN);
                AEi(ia,:) = yi;
            end
        end
        AE = AEi;
    end

    % coarsen (block-average)
    if (COARSEN_ANGLE_FACTOR>1) || (COARSEN_ENERGY_FACTOR>1)
        fa = COARSEN_ANGLE_FACTOR;
        fe = COARSEN_ENERGY_FACTOR;
        nA2 = floor(size(AE,1)/fa)*fa;
        nE2 = floor(size(AE,2)/fe)*fe;
        if nA2 < 1 || nE2 < 1
            warning('Coarsen factors too large for AE. Skipping coarsen.');
        else
            AE = AE(1:nA2, 1:nE2);
            AE = reshape(AE, fa, nA2/fa, fe, nE2/fe);
            AE = squeeze(mean(mean(AE,1,'omitnan'),3,'omitnan'));
            AE = reshape(AE, nA2/fa, nE2/fe);
        end
    end

    % scale & clean
    AE = scale .* AE;
    AE(AE <= 0) = NaN;

    AEmat = AE;
    nA_out = size(AE,1);
    nE_out = size(AE,2);
end

%% ========================================================================
% Helper: build_resample_energy_axis
% - construct E_new centers (length nNew) from E_old centers
% - uses linear spacing across [min(E_old), max(E_old)] unless E_old suggests edges
%% ========================================================================
function E_new = build_resample_energy_axis(E_old, nNew, Emin, Emax)
    % If E_old looks log-spaced, build geometric linspace; otherwise linear.
    ratios = diff(E_old) ./ E_old(1:end-1);
    if all(ratios > 0) && (std(log10(E_old(2:end)./E_old(1:end-1))) < 1e-6)
        % approximately geometric spacing -> make geometric resample
        E_new = logspace(log10(Emin + eps), log10(Emax + eps), nNew);
    else
        % fallback linear in energy (but good enough)
        E_new = linspace(Emin, Emax, nNew);
    end
end

%% ========================================================================
% Helper: try_read_energy_axis
% Tries several common variable names; returns centers (not edges) and var name.
%% ========================================================================
function [E, used_name] = try_read_energy_axis(ncfile, nE, Emin, Emax)
    candidates = { ...
        'energy','Energy','E','egrid','Egrid','EGrid', ...
        'surfEGrid','surfEgrid','surfEnergy','energyGrid', ...
        'eBins','EBins','Ebins','energy_bins','energyBins', ...
        'e_edges','E_edges','energy_edges','energyEdges', ...
        'energy_centers','E_centers' ...
    };

    E = [];
    used_name = '';
    for i = 1:numel(candidates)
        v = candidates{i};
        try
            tmp = ncread(ncfile, v);
            tmp = double(tmp(:));
            if numel(tmp) == nE
                E = tmp(:).';
                used_name = v;
                break;
            elseif numel(tmp) == nE+1
                % edges -> centers
                E = 0.5*(tmp(1:end-1) + tmp(2:end));
                used_name = v;
                break;
            end
        catch
            % ignore
        end
    end

    if isempty(E)
        used_name = 'FALLBACK_linspace';
        E = linspace(Emin, Emax, nE);
    end

    % If not strictly increasing -> fallback
    if any(diff(E) <= 0)
        warning('Energy axis "%s" not strictly increasing. Falling back to uniform linspace.', used_name);
        used_name = 'FALLBACK_linspace_nonmonotonic';
        E = linspace(Emin, Emax, nE);
    end
end

%% ========================================================================
% Helper: colorbar formatting (decades)
%% ========================================================================
function format_colorbar_decades(cb, clim)
    lo = clim(1); hi = clim(2);
    if lo <= 0 || hi <= 0 || ~isfinite(lo) || ~isfinite(hi) || lo >= hi
        return;
    end
    e1 = floor(log10(lo));
    e2 = ceil(log10(hi));
    exps = e1:e2;
    cb.Ticks = 10.^exps;
    cb.TickLabels = arrayfun(@(e) sprintf('10^{%d}', e), exps, 'UniformOutput', false);
    cb.Label.String = 'iEADS (linear units)';
end

%% ========================================================================
% Helper: tiny ternary (truthy string)
%% ========================================================================
function out = ternary(cond, a, b)
    if cond, out = a; else out = b; end
end