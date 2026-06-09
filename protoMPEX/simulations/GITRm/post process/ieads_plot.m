%% revised_script_bins_layout_fixed_with_reflected.m
% Full script with reflected EAD plotted in the same 2x3 layout as incident.
% Uses energy BIN INDICES on x-axis by default (falls back to indices).
clc; clear; close all;

%% ----------------------------
% CONFIG
%% ----------------------------
tilt = 0;   % degrees
file = sprintf('../tilted_targets/test/%d_degrees/MPEX_runs_test/gitrm-surface.nc', tilt);

% surface selection
SKIP_FIRST_NSURFACES = 3; % 6 for 85 degrees case
DROP_LAST_NSURFACES  = 3;

% resample / coarsen
TARGET_N_E_BINS = 5000;  % set [] to skip
COARSEN_ANGLE_FACTOR  = 2;
COARSEN_ENERGY_FACTOR = 10;

% angle marker
TARGET_ANGLE_BIN = 180;

% scale (incident scaling)
scale = 1;

% charge bins (incident)
BIN_COMBINED = 1;
BIN_TA = 3:6;  % Ta1+..Ta4+
if numel(BIN_TA) ~= 4, error('BIN_TA must be exactly 4 bins (3:6).'); end

% color limits
CLIM_2D_INC  = [1e-5 1e1];   % incident default
CLIM_2D_REFL = [1e-2 1e-1];  % reflected default (your suggestion)

% energy zoom in BIN INDEX units (integer indices)
E_ZOOM_MIN = 0;
E_ZOOM_MAX = 100;

% x-axis mode
FORCE_LOG_X_ON_BINS = false;   % <-- keep false unless you want log-x visual

% reflected options
PLOT_REFLECTED = true;
NORMALIZE_REFLECTED = true;    % normalize reflected 2D to its total flux (sum over E,angle)
PLOT_REFLECTED_CHARGE_PANELS = false; % false: mark Ta1..Ta4 as N/A (no charge-resolved refl avail)

%% ----------------------------
% READ DATA
%% ----------------------------
surfEDist_all    = ncread(file,'surfEDist');     % [A,E,S,species,bins]
surfReflDist_all = ncread(file,'surfReflDist');  % [A,E,S,species]

[nA0, nE0, nSurf0, nSpecies0, nBins0] = size(surfEDist_all);
fprintf('surfEDist dims: A=%d, E=%d, S=%d, species=%d, bins=%d\n', nA0, nE0, nSurf0, nSpecies0, nBins0);
if nBins0 < 6, error('Need at least 6 bins for Ta1+..Ta4+.'); end

%% ----------------------------
% PREP reflected species-summed (no charge dimension)
%% ----------------------------
refl_AES = squeeze(sum(surfReflDist_all,4)); % [A,E,S]
[refl_AE, nA_r, nE_r] = prep_AES_to_AE(refl_AES, SKIP_FIRST_NSURFACES, DROP_LAST_NSURFACES, ...
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
    % validate b
    if ~isscalar(b) || b~=floor(b) || b < 1 || b > nBins0
        error('Invalid bin index bins_to_plot(%d) = %s; valid range 1..%d', k, mat2str(b), nBins0);
    end

    % safe extraction (preserve expected dims)
    tmp = surfEDist_all(:,:,:,:,b); % may be 4-D or squeezed
    expected4 = [nA0, nE0, nSurf0, nSpecies0];
    if ~isequal(size(tmp), expected4)
        if numel(tmp) == prod(expected4)
            tmp = reshape(tmp, expected4);
        else
            tmp = reshape(tmp, [size(tmp,1), size(tmp,2), 1, size(tmp,3)]);
        end
    end

    % sum over species -> AES [A x E x S]
    AES = squeeze(sum(tmp,4));
    if ndims(AES) ~= 3
        error('After sum over species, AES must be 3-D. size(AES)=%s', mat2str(size(AES)));
    end

    [inc_AE, nA_k, nE_k] = prep_AES_to_AE(AES, SKIP_FIRST_NSURFACES, DROP_LAST_NSURFACES, ...
                                         TARGET_N_E_BINS, COARSEN_ANGLE_FACTOR, COARSEN_ENERGY_FACTOR, scale);
    inc_AE_list{k} = inc_AE;
    nA_list(k) = nA_k;
    nE_list(k) = nE_k;
end

% unify dims
nA_min = min(nA_list);
nE_min = min(nE_list);
for k = 1:nPlot
    inc_AE_list{k} = inc_AE_list{k}(1:nA_min, 1:nE_min);
end
refl_AE = refl_AE(1:nA_min, 1:nE_min);

% energy bins to plot (integer indices, clipped to available)
E_ZOOM_MIN = max(1, round(E_ZOOM_MIN));
E_ZOOM_MAX = max(E_ZOOM_MIN, round(E_ZOOM_MAX));
E_ZOOM_MAX = min(E_ZOOM_MAX, nE_min);
Ebins = E_ZOOM_MIN:E_ZOOM_MAX;       % integer index vector for columns

% numeric index vector for safe indexing
idxE = Ebins;                        % guaranteed integer indices within 1..nE_min

angle_idx = min(max(1,TARGET_ANGLE_BIN), nA_min);

%% ----------------------------
% INCIDENT PLOTTING (original 2x3 layout)
%% ----------------------------
figure('Color','w','Position',[80 80 1400 900]);
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
colormap(jet);
sgtitle(sprintf('IEADs at MPEX target (tilt = %d^\\circ) - INCIDENT', tilt), 'FontSize',16,'FontWeight','bold');

apply_xaxis = @(ax) set_xaxis_mode(ax, FORCE_LOG_X_ON_BINS, E_ZOOM_MIN, E_ZOOM_MAX);

% (1) Combined 2D
nexttile(1);
AE_comb = inc_AE_list{1};
pcolor(Ebins, 1:nA_min, AE_comb(:, idxE)); shading flat;
set(gca,'YDir','normal','ColorScale','linear');
caxis(CLIM_2D_INC); colorbar;
xlabel('Energy [eV]'); ylabel('Angle [\theta]');
title(sprintf('Combined IEADs (bin=%d)', BIN_COMBINED));
apply_xaxis(gca);
hold on; plot([E_ZOOM_MIN E_ZOOM_MAX],[angle_idx angle_idx],'w--','LineWidth',1.5);

% (2) Ta1+ 2D
nexttile(2);
AE = inc_AE_list{2};
pcolor(Ebins, 1:nA_min, AE(:, idxE)); shading flat;
set(gca,'YDir','normal','ColorScale','linear');
caxis(CLIM_2D_INC); colorbar;
xlabel('Energy [eV]'); ylabel('Angle [\theta]');
title('Ta1^+');
apply_xaxis(gca);
hold on; plot([E_ZOOM_MIN E_ZOOM_MAX],[angle_idx angle_idx],'w--','LineWidth',1.2);

% (3) Ta2+ 2D
nexttile(3);
AE = inc_AE_list{3};
pcolor(Ebins, 1:nA_min, AE(:, idxE)); shading flat;
set(gca,'YDir','normal','ColorScale','linear');
caxis(CLIM_2D_INC); colorbar;
xlabel('Energy [eV]'); ylabel('Angle [\theta]');
title('Ta2^+');
apply_xaxis(gca);
hold on; plot([E_ZOOM_MIN E_ZOOM_MAX],[angle_idx angle_idx],'w--','LineWidth',1.2);

% (4) 1D spectra: Combined + Ta1..Ta4
nexttile(4); hold on; box on;
spec_comb = nanmean(AE_comb(:, idxE),1);
if FORCE_LOG_X_ON_BINS
    semilogx(Ebins, spec_comb, 'k', 'LineWidth', 2.5);
else
    plot(Ebins, spec_comb, 'k', 'LineWidth', 2.5);
end

ta_labels = {'Ta1^+','Ta2^+','Ta3^+','Ta4^+'};
for j = 1:4
    AEj = inc_AE_list{1+j};
    specj = nanmean(AEj(:, idxE),1);
    if FORCE_LOG_X_ON_BINS
        semilogx(Ebins, specj, 'LineWidth', 1.6);
    else
        plot(Ebins, specj, 'LineWidth', 1.6);
    end
end
grid on;
xlabel('Energy [eV]');
ylabel('\langle iEADS \rangle_\theta');
title('Angle-averaged iEADS (incident)');
apply_xaxis(gca);
legend(['Combined',ta_labels],'Location','northeast');

% (5) Ta3+ 2D
nexttile(5);
AE = inc_AE_list{4};
pcolor(Ebins, 1:nA_min, AE(:, idxE)); shading flat;
set(gca,'YDir','normal','ColorScale','linear');
caxis(CLIM_2D_INC); colorbar;
xlabel('Energy [eV]'); ylabel('Angle [\theta]');
title('Ta3^+');
apply_xaxis(gca);
hold on; plot([E_ZOOM_MIN E_ZOOM_MAX],[angle_idx angle_idx],'w--','LineWidth',1.2);

% (6) Ta4+ 2D
nexttile(6);
AE = inc_AE_list{5};
pcolor(Ebins, 1:nA_min, AE(:, idxE)); shading flat;
set(gca,'YDir','normal','ColorScale','linear');
caxis(CLIM_2D_INC); colorbar;
xlabel('Energy [eV]'); ylabel('Angle [\theta]');
title('Ta4^+');
apply_xaxis(gca);
hold on; plot([E_ZOOM_MIN E_ZOOM_MAX],[angle_idx angle_idx],'w--','LineWidth',1.2);

fprintf('Incident plotting done.\n');

%% ----------------------------
% REFLECTED EAD (single 2D panel)
%% ----------------------------
NORMALIZE_REFLECTED = true;     % normalize by total reflected flux (sum over E and angle)
CLIM_2D_REFL = [1e-5 1e-2];     % your suggested range (adjust if needed)

refl2 = refl_AE(:, idxE);       % [angle x energy] subset

if NORMALIZE_REFLECTED
    tot_ref = nansum(refl2(:));
    if tot_ref > 0
        refl2 = refl2 ./ tot_ref;
    end
end

figure('Color','w','Position',[200 200 900 650]);
pcolor(Ebins, 1:nA_min, refl2); shading flat;
set(gca,'YDir','normal','ColorScale','log');
caxis(CLIM_2D_REFL); colorbar;
colormap("jet")

xlabel('Energy [eV]');           % use "Energy bin" since you're using bin indices
ylabel('Angle [\theta]');
title(sprintf('Reflected EAD (species-summed) | tilt=%d^\\circ', tilt));

set_xaxis_mode(gca, FORCE_LOG_X_ON_BINS, E_ZOOM_MIN, E_ZOOM_MAX);
hold on;
plot([E_ZOOM_MIN E_ZOOM_MAX],[angle_idx angle_idx],'w--','LineWidth',1.5);

%% ========================================================================
% Local function: prep_AES_to_AE
% Same as original: AESmat = [nAngles x nEnergies x nSurfaces]
% Output: AEmat  = [nAngles' x nEnergies'] after surface select, resample, coarsen, scale, clean
%% ========================================================================
function [AEmat, nA_out, nE_out] = prep_AES_to_AE( ...
    AESmat, SKIP_FIRST_NSURFACES, DROP_LAST_NSURFACES, TARGET_N_E_BINS, ...
    COARSEN_ANGLE_FACTOR, COARSEN_ENERGY_FACTOR, scale)

    % surfaces select
    nSurf = size(AESmat,3);
    i0 = SKIP_FIRST_NSURFACES + 1;
    i1 = nSurf - DROP_LAST_NSURFACES;
    if i0 > i1
        error('Surface selection invalid: i0=%d, i1=%d. Check SKIP_FIRST_NSURFACES/DROP_LAST_NSURFACES.', i0, i1);
    end
    AESsel = AESmat(:,:,i0:i1);

    % sum surfaces -> [A,E]
    AE = squeeze(sum(AESsel,3));
    [nA, nE_orig] = size(AE);

    % resample energy (index-based) if needed
    if ~isempty(TARGET_N_E_BINS) && TARGET_N_E_BINS ~= nE_orig
        x_old = 1:nE_orig;
        x_new = linspace(1,nE_orig,TARGET_N_E_BINS);
        AEi = nan(nA, numel(x_new));
        for ia = 1:nA
            y = AE(ia,:);
            if all(isnan(y))
                AEi(ia,:) = NaN;
            else
                mask = ~isnan(y);
                y(~mask) = 0; % set NaNs to 0 for interp
                yi = interp1(x_old, y, x_new, 'linear', 0);
                wi = interp1(x_old, double(mask), x_new, 'linear', 0);
                yi(wi < 0.5) = NaN;
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
    AE(AE <= 0.05) = NaN;

    AEmat = AE;
    nA_out = size(AE,1);
    nE_out = size(AE,2);
end

%% ========================================================================
% Local function: set_xaxis_mode
%% ========================================================================
function set_xaxis_mode(ax, force_log, xmin, xmax)
    if force_log
        set(ax, 'XScale', 'log');
        if xmin <= 0, xmin = 1; end
        xlim(ax, [max(1,xmin), max(2,xmax)]);
    else
        set(ax, 'XScale', 'linear');
        xlim(ax, [xmin, xmax]);
    end
end

%% ========================================================================
% Local helper: ternary
%% ========================================================================
function out = ternary(cond, a, b)
    if cond, out = a; else out = b; end
end