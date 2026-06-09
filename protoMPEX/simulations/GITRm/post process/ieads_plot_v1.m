%% ieads_plot_v1.m
% Per-charge (Ta0+..Ta4+) per-nBins 2D IEADS plots (incident surfEDist)
% Robust CLim (percentile) so higher charge states are visible.
% Compatible with older MATLAB versions (no colorbar(t,'eastoutside') syntax).

clc; clear; close all;

%% ----------------------------
% CONFIG (edit as needed)
%% ----------------------------
tilt = 45;
file = sprintf('../tilted_targets/test/%d_degrees/MPEX_runs/gitrm-surface.nc', tilt);

SKIP_FIRST_NSURFACES = 3;
DROP_LAST_NSURFACES  = 3;

E_YMAX = 500;                  % show up to this energy-bin index
COARSEN_ANGLE_FACTOR  = 5;     % set 1 to disable
COARSEN_ENERGY_FACTOR = 20;    % set 1 to disable

scale = 2.331579550458116e+09; % your scale

MAX_CHARGE_TO_SHOW = 4;        % plots Ta^0 .. Ta^4+ => species indices 1..5

% CLim options:
%   'per_charge_percentile' : robust CLim per charge (recommended)
%   'global_percentile'     : one CLim across all charges+bins (absolute comparison)
%   'normalize_per_charge'  : normalize each charge to its own max (shape-only)
CLIM_METHOD = 'per_charge_percentile';

P_LO = 0.5;    % lower percentile for CLim (ignore extreme low outliers)
P_HI = 99.5;   % upper percentile for CLim (ignore extreme hot pixels)
NORM_CLIM_LO = 1e-4;  % if normalize_per_charge, show down to this fraction of max

%% ----------------------------
% READ DATA
%% ----------------------------
surfEDist_all = ncread(file,'surfEDist');   % nA x nE x nSurf x nSpecies x nBins
[nA, nE, nSurfTotal, nSpecies, nBins] = size(surfEDist_all);

fprintf('Loaded surfEDist: nA=%d, nE=%d, nSurf=%d, nSpecies=%d, nBins=%d\n', ...
    nA, nE, nSurfTotal, nSpecies, nBins);

% Surface selection
i0 = SKIP_FIRST_NSURFACES + 1;
i1 = nSurfTotal - DROP_LAST_NSURFACES;
if i0 > i1
    error('Invalid surface selection: i0=%d i1=%d. Check SKIP_FIRST_NSURFACES/DROP_LAST_NSURFACES.', i0, i1);
end

% Slice desired surfaces
surfEDist_all = surfEDist_all(:,:,i0:i1,:,:);

% Sum over surfaces -> angle x energy x species x bins
E_sumsurf = squeeze(sum(surfEDist_all, 3));

% Charges to show: Ta0..Ta4+ => indices 1..5
nShowCharges = min(nSpecies, MAX_CHARGE_TO_SHOW+1);
if nShowCharges < (MAX_CHARGE_TO_SHOW+1)
    warning('File has nSpecies=%d; plotting only species 1..%d (Ta^0..Ta^{%d}+).', ...
        nSpecies, nShowCharges, nShowCharges-1);
end

%% ----------------------------
% BUILD PANELS (scale + clean + coarsen) in a cell array
%% panels{s}{b} = 2D matrix (angle x energy) for species s and bin b
%% ----------------------------
panels = cell(nShowCharges, 1);
for s = 1:nShowCharges
    panels{s} = cell(nBins, 1);
    for b = 1:nBins
        A = squeeze(E_sumsurf(:,:,s,b));                  % nA x nE
        A = local_coarsen2D(A, COARSEN_ANGLE_FACTOR, COARSEN_ENERGY_FACTOR);
        A = scale .* A;
        A(A<=0) = NaN;
        panels{s}{b} = A;
    end
end

%% ----------------------------
% COMPUTE CLim(s)
%% ----------------------------
switch lower(CLIM_METHOD)
    case 'per_charge_percentile'
        CLims = nan(nShowCharges,2);
        for s = 1:nShowCharges
            vals = [];
            for b = 1:nBins
                A = panels{s}{b};
                v = A(~isnan(A) & A>0);
                vals = [vals; v(:)];
            end
            if isempty(vals)
                continue
            end
            lo = prctile(vals, P_LO);
            hi = prctile(vals, P_HI);
            lo = max(lo, min(vals));     % guard
            hi = max(hi, lo*10);         % guard
            clo = 10^floor(log10(max(lo, eps)));
            chi = 10^ceil(log10(max(hi, clo*10)));
            CLims(s,:) = [clo chi];
        end

    case 'global_percentile'
        vals_all = [];
        for s = 1:nShowCharges
            for b = 1:nBins
                A = panels{s}{b};
                v = A(~isnan(A) & A>0);
                vals_all = [vals_all; v(:)];
            end
        end
        if isempty(vals_all)
            CLims = nan(nShowCharges,2);
        else
            lo = prctile(vals_all, P_LO);
            hi = prctile(vals_all, P_HI);
            clo = 10^floor(log10(max(lo, eps)));
            chi = 10^ceil(log10(max(hi, clo*10)));
            CLims = repmat([clo chi], nShowCharges, 1);
        end

    case 'normalize_per_charge'
        CLims = repmat([NORM_CLIM_LO 1], nShowCharges, 1);

    otherwise
        error('Unknown CLIM_METHOD: %s', CLIM_METHOD);
end

%% ----------------------------
% PLOT: one figure per charge, each with 12 tiles (3x4) for nBins
%% ----------------------------
for s = 1:nShowCharges
    figure('Color','w','Units','normalized','Position',[0.03 0.03 0.94 0.90]);
    tiledlayout(3,4,'Padding','compact','TileSpacing','compact');
    colormap(jet)

    % For normalize_per_charge, find per-charge max across bins
    if strcmpi(CLIM_METHOD,'normalize_per_charge')
        vmax = -inf;
        for b = 1:nBins
            vmax = max(vmax, nanmax(panels{s}{b}(:)));
        end
        if ~isfinite(vmax) || vmax<=0
            vmax = NaN;
        end
    else
        vmax = NaN;
    end

    axHandles = gobjects(nBins,1);

    for b = 1:nBins
        nexttile
        A = panels{s}{b};

        if strcmpi(CLIM_METHOD,'normalize_per_charge') && isfinite(vmax)
            Aplot = A ./ vmax;
        else
            Aplot = A;
        end

        ymaxE = min(E_YMAX, size(Aplot,2));
        imagesc(1:ymaxE, 1:size(Aplot,1), Aplot(:,1:ymaxE));
        set(gca,'YDir','normal','ColorScale','log');

        title(sprintf('Bin %d', b), 'FontWeight','bold');
        if b==1
            xlabel('Energy bin'); ylabel('Angle bin');
        else
            set(gca,'XTickLabel',[],'YTickLabel',[]);
        end

        axHandles(b) = gca;
    end

    % Apply CLim across all tiles for this charge
    clim_use = CLims(s,:);
    if all(isfinite(clim_use))
        for k = 1:numel(axHandles)
            set(axHandles(k), 'CLim', clim_use);
        end
    else
        warning('No finite CLim for species %d (Ta^%d+). Leaving autoscale.', s, s-1);
    end

    % Single shared colorbar (compatible syntax)
    cb = colorbar;
    cb.Location = 'eastoutside';
    cb.Label.String = 'Scaled flux (incident, per-bin)';
    cb.Label.FontSize = 10;

    % Optional: set log-spaced ticks if CLim known
    if all(isfinite(clim_use))
        cb.Ticks = logspace(log10(clim_use(1)), log10(clim_use(2)), 6);
    end

    sgtitle(sprintf('Incident IEADS per nBins | Ta^{%d+} (species %d) | tilt=%d° | surfaces %d..%d', ...
        s-1, s, tilt, i0, i1), 'FontSize', 14);
end

fprintf('Done. Plotted %d species (requested up to Ta^%d+).\n', nShowCharges, MAX_CHARGE_TO_SHOW);

%% ----------------------------
% local function: coarsen 2D by block averaging (angle x energy)
%% ----------------------------
function Aout = local_coarsen2D(A, fa, fe)
    if nargin < 2, fa = 1; end
    if nargin < 3, fe = 1; end
    if fa<=1 && fe<=1
        Aout = A;
        return;
    end
    [nA, nE] = size(A);
    nAtrim = floor(nA/fa)*fa;
    nEtrim = floor(nE/fe)*fe;
    if nAtrim < 1 || nEtrim < 1
        warning('Coarsen factors too large; returning original.');
        Aout = A;
        return;
    end
    A = A(1:nAtrim, 1:nEtrim);
    A = reshape(A, fa, nAtrim/fa, fe, nEtrim/fe);
    A = squeeze(mean(mean(A,1,'omitnan'),3,'omitnan'));
    Aout = reshape(A, nAtrim/fa, nEtrim/fe);
end