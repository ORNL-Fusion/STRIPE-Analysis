clc; clear; close all; format long;

%% ==============================
%  1) EFIT grid & limiter geometry
% ==============================
load('extrapolated_data_196154.mat');  % provides X, Y, g (incl. g.lim)
% load('extrapolated_data_200882.mat');  % provides X, Y, g (incl. g.lim)

%% ==============================
%  2) Settings
% ==============================
psiN_mask = 1.0;                   % evaluate inside ψ_N ≤ 1
dt = 1e-9;                         % [s]
nP = 5e6;                          % number of launched particles
% source_strength = 7.176454593335885e+20; % [1/s] for 200882
source_strength = 9.442619262858435e+17; % # 196154 % particles/sec from GITRm


%% ==============================
%  3) Load GITRm spec data
% ==============================
path1 = "../diiid-helicon/DIII-D_helicon_runs_full_196154/";
% path1 = "../diiid-helicon/DIII-D_helicon_runs_full_200882/";
% path1 = "../diiid-helicon/DIII-D_helicon_runs_high_200882/";

file1 = "gitrm-spec.nc";

gridr_2d = ncread(path1 + file1, 'gridr_2d');  % [nR] cell centers
gridz_2d = ncread(path1 + file1, 'gridz_2d');  % [nZ] cell centers
[x, y]   = meshgrid(gridr_2d, gridz_2d);       % [nZ × nR] centers

data_1   = ncread(path1 + file1, 'n_2d');      % [nR × nZ × 1 × nSpec] (counts/bin)
data_trim = squeeze(data_1(:,:,1,1))';         % -> [nZ × nR]
data_trim(data_trim <= 0) = NaN;

%% ---- Symmetric cell-edge construction (prevents outer-cell volume inflation)
% R edges
dr      = diff(gridr_2d);
drL     = dr(1);         drR = dr(end);
r_edges = [gridr_2d(1)-drL/2;  (gridr_2d(1:end-1)+gridr_2d(2:end))/2;  gridr_2d(end)+drR/2];

% Z edges
dz      = diff(gridz_2d);
dzL     = dz(1);         dzR = dz(end);
z_edges = [gridz_2d(1)-dzL/2;  (gridz_2d(1:end-1)+gridz_2d(2:end))/2;  gridz_2d(end)+dzR/2];

% Cell widths
dR = diff(r_edges);      % [nR]
dZ = diff(z_edges);      % [nZ]

% Toroidal volume per cell (ring): V = 2π ΔZ * ( (r_out^2 - r_in^2)/2 )
% Build arrays [nZ×nR]
r_in  = r_edges(1:end-1);
r_out = r_edges(2:end);
ring_area = 0.5*(r_out.^2 - r_in.^2);                  % [nR]
vol_RZ = 2*pi * (dZ(:)) .* ring_area(:)';              % [nZ × nR]

%% ---- Convert raw counts to physical density [m^-3]
density_gitr = source_strength .* data_trim .* dt ./ (nP .* vol_RZ);

%% ==============================
%  4) ψ_N on the SAME (R,Z) grid
% ==============================
[Xg, Yg] = meshgrid(gridr_2d, gridz_2d);  % [nZ x nR]
psiN_flat = calc_psiN(g, Xg(:), Yg(:), []);
psiN_gitr = reshape(psiN_flat, size(Xg)); % [nZ x nR]

% Inclusive limiter mask: inside OR on boundary
[in_lim, on_lim] = inpolygon(Xg, Yg, g.lim(1,:), g.lim(2,:));
mask_limiter_inc = in_lim | on_lim;

% Keep ψN mask but do NOT wipe the last inboard cell near wall
mask_inside = (psiN_gitr <= psiN_mask) & mask_limiter_inc;

% Outside vessel → NaN (not 0) so log plots don't force a fake drop
density_gitr(~mask_limiter_inc) = NaN;

%% ==============================
%  5) Load BACKGROUND multi-species; total carbon + n_e
% ==============================
bkgfile = 'profilesDIIID_196154_multi.nc';
% bkgfile = 'profilesDIIID_200882_multi.nc';
R_bkg   = ncread(bkgfile, 'x');                    % [nR_b]
Z_bkg   = ncread(bkgfile, 'z');                    % [nZ_b]
psiN_b  = ncread(bkgfile, 'psiN');                 % [nR_b x nZ_b]
ne_b    = ncread(bkgfile, 'ne');                   % [nR_b x nZ_b]
ni_all  = ncread(bkgfile, 'ni_all');               % [nR_b x nZ_b x ns]
Z_all   = ncread(bkgfile, 'atomic_number');        % [ns]

idxC = find(Z_all == 6);
if isempty(idxC), error('No carbon species (Z==6) found in background NetCDF.'); end

niC_total_b = sum(ni_all(:,:,idxC), 3);            % [nR_b x nZ_b]

% Transpose to [nZ_b x nR_b] for interp2 with meshgrid(R_bkg,Z_bkg)
psiN_b_T = psiN_b';   ne_b_T = ne_b';   niC_b_T = niC_total_b';

% Interpolate to GITR grid
[Rb2, Zb2] = meshgrid(R_bkg, Z_bkg); % [nZ_b x nR_b]
psiN_b_onG = interp2(Rb2, Zb2, psiN_b_T, Xg, Yg, 'linear', NaN);
ne_onG     = interp2(Rb2, Zb2, ne_b_T,   Xg, Yg, 'linear', NaN);
nC_onG     = interp2(Rb2, Zb2, niC_b_T,  Xg, Yg, 'linear', NaN);

% Apply vessel mask consistently
nC_onG(~mask_limiter_inc) = NaN;  ne_onG(~mask_limiter_inc) = NaN;

% Carbon fraction
fracC_over_ne = nC_onG ./ ne_onG;

%% ==============================
%  6) Midplane extraction and ψ_N-based comparison
% ==============================
[~, mid_idx] = min(abs(gridz_2d));    % z ~ 0
z_mid = gridz_2d(mid_idx);

psiN_mid     = psiN_gitr(mid_idx, :);
nC_b_mid     = nC_onG(mid_idx, :);
ne_mid       = ne_onG(mid_idx, :);
frac_mid     = fracC_over_ne(mid_idx, :);
nC_gitr_mid  = density_gitr(mid_idx, :);

% Clean invalid, sort by ψN
valid = isfinite(psiN_mid) & isfinite(nC_b_mid) & isfinite(nC_gitr_mid);
psiN_mid_f = psiN_mid(valid);  nC_b_mid_f = nC_b_mid(valid);
nC_g_mid_f = nC_gitr_mid(valid); frac_mid_f = frac_mid(valid);

[psiN_mid_s, sidx] = sort(psiN_mid_f);
nC_b_mid_s = nC_b_mid_f(sidx);
nC_g_mid_s = nC_g_mid_f(sidx);
frac_mid_s = frac_mid_f(sidx);

psiN_common = linspace(0, 1.2, 300);
nC_b_vspsi    = interp1(psiN_mid_s, nC_b_mid_s, psiN_common, 'pchip', 'extrap');
nC_gitr_vspsi = interp1(psiN_mid_s, nC_g_mid_s, psiN_common, 'pchip', 'extrap');
frac_vspsi    = interp1(psiN_mid_s, frac_mid_s, psiN_common, 'pchip', 'extrap');

%% ==============================
%  7) Plots
% ==============================
% (a) Map sanity
figure('Name','GITR Carbon Density Map');
h = pcolor(Xg, Yg, density_gitr); set(h,'EdgeColor','none');
xlabel('R [m]'); ylabel('Z [m]'); axis equal tight;
xlim([min(gridr_2d) max(gridr_2d)]); ylim([min(gridz_2d) max(gridz_2d)]);
cb = colorbar; ylabel(cb,'n_C [m^{-3}]');
title('GITR Carbon Density (physical)'); 
set(gca,'ColorScale','log','FontSize',14);
hold on;
contour(Xg, Yg, psiN_gitr, [1 1], 'LineStyle', '--', ...
    'LineColor', 'k', 'LineWidth', 1.2);
plot(g.lim(1,:), g.lim(2,:), 'k', 'LineWidth', 1.0);

%% (b) Midplane densities vs ψ_N  (C and n_e)
%% --- Robust ψ_N-binned midplane profiles (reduces SOL oscillations)
psiN_mid_raw = psiN_gitr(mid_idx, :);
nC_g_raw     = density_gitr(mid_idx, :);
nC_b_raw     = nC_onG(mid_idx, :);
ne_raw       = ne_onG(mid_idx, :);

% Keep finite samples only
ok = isfinite(psiN_mid_raw) & isfinite(nC_g_raw) & isfinite(nC_b_raw) & isfinite(ne_raw);
psiN_mid_raw = psiN_mid_raw(ok);
nC_g_raw     = nC_g_raw(ok);
nC_b_raw     = nC_b_raw(ok);
ne_raw       = ne_raw(ok);

% Define ψ bins focused on edge region (you can widen if needed)
psi_edges = [0.00:0.01:1.20];              % 0.01-wide bins
[~,~,bin] = histcounts(psiN_mid_raw, psi_edges);

% Bin-aggregate using median (robust) — fallback to mean if single sample
nbins = numel(psi_edges)-1;
psi_bin_centers = 0.5*(psi_edges(1:end-1)+psi_edges(2:end));
nC_g_binned = nan(1,nbins); nC_b_binned = nan(1,nbins); ne_b_binned = nan(1,nbins);
for k = 1:nbins
    idx = (bin==k);
    if any(idx)
        nC_g_binned(k) = median(nC_g_raw(idx), 'omitnan');
        nC_b_binned(k) = median(nC_b_raw(idx), 'omitnan');
        ne_b_binned(k) = median(ne_raw(idx),   'omitnan');
    end
end

% Optional gentle smoothing to kill single-bin spikes
nC_g_binned = smoothdata(nC_g_binned, 'movmedian', 3);
nC_b_binned = smoothdata(nC_b_binned, 'movmedian', 3);
ne_b_binned = smoothdata(ne_b_binned, 'movmedian', 3);

%% --- Plot ONLY ψ_N > 0.88
mask_zoom = psi_bin_centers > 0.88;      % your requested zoom range

figure('Name','Midplane (ψ_N>0.88): GITR & Background Densities vs ψ_N');
semilogy(psi_bin_centers(mask_zoom), nC_g_binned(mask_zoom), 'r-', 'LineWidth', 2); hold on;
semilogy(psi_bin_centers(mask_zoom), nC_b_binned(mask_zoom), 'b--', 'LineWidth', 2);
semilogy(psi_bin_centers(mask_zoom), ne_b_binned(mask_zoom), 'k-.', 'LineWidth', 2);
% --- Add separatrix line at ψ_N = 1 ---
xline(1, 'k--', 'LineWidth', 1.2, ...
    'Label', '\psi_N = 1', ...
    'LabelVerticalAlignment', 'bottom', ...
    'LabelHorizontalAlignment', 'left');

xlabel('\psi_N'); ylabel('Density [m^{-3}]');
title(sprintf('Midplane (Z = %.3f m): GITR & Background Densities vs \\psi_N', z_mid));
legend({'GITR Carbon','Background Carbon','Background n_e','\psi_N = 1'}, 'Location','best');
grid on; set(gca,'FontSize',14);
xlim([0.9 1.1]);  % full range or adjust to [0.88 1.2] if zoo



%% ==============================
%  9) Midplane profiles vs R (with n_e, LCFS, and LFS wall)
% ==============================
figure('Name','Midplane Densities vs R (LFS separatrix & wall)');

R_mid       = gridr_2d(:)';                     
nC_gitr_mid = density_gitr(mid_idx, :);         
nC_bkg_mid  = nC_onG(mid_idx, :);               
ne_bkg_mid  = ne_onG(mid_idx, :);               
psiN_mid_full = psiN_gitr(mid_idx, :);          

% Plot
semilogy(R_mid, nC_gitr_mid, 'r-', 'LineWidth', 2); hold on;
semilogy(R_mid, nC_bkg_mid,  'b--', 'LineWidth', 2);
semilogy(R_mid, ne_bkg_mid,  'k-.', 'LineWidth', 2);

% Magnetic axis to define LFS
R_axis = isfield(g,'rmaxis') * g.rmaxis + ~isfield(g,'rmaxis') * mean(R_mid(~isnan(psiN_mid_full)));

% LFS separatrix (ψ_N=1 with R>R_axis)
valid_idx = find(isfinite(psiN_mid_full) & R_mid > R_axis);
if numel(valid_idx) > 2
    psiN_LFS = psiN_mid_full(valid_idx);
    R_LFS    = R_mid(valid_idx);
    [psiN_LFS, sI] = sort(psiN_LFS); R_LFS = R_LFS(sI);
    try
        R_sep = interp1(psiN_LFS, R_LFS, 1, 'linear', 'extrap');
        xline(R_sep, 'k--', 'LineWidth', 1.5, 'Label','\psi_N=1 (LFS)', ...
            'LabelVerticalAlignment','bottom', 'LabelHorizontalAlignment','left');
    catch
        warning('LFS separatrix interpolation failed.');
    end
else
    warning('Insufficient valid ψ_N points to locate LFS separatrix.');
end

% OUTER wall (LFS) from g.lim near midplane
R_limiter = g.lim(1,:);  Z_limiter = g.lim(2,:);
z_tol = 0.05;  % ±5 cm around Z=0
mask_LFS_wall = (abs(Z_limiter) < z_tol) & (R_limiter > R_axis);
if any(mask_LFS_wall)
    R_wall = max(R_limiter(mask_LFS_wall));  % furthest-out point on LFS midplane
    xline(R_wall, 'k-', 'LineWidth', 1.5, 'Label','Wall (LFS)', ...
        'LabelVerticalAlignment','bottom', 'LabelHorizontalAlignment','left');
else
    warning('No limiter points near Z=0 on LFS.');
end

xlabel('R [m]'); ylabel('Density [m^{-3}]');
title(sprintf('Midplane Densities vs R (z = %.3f m, LFS)', z_mid));
legend({'GITR Carbon','Background Carbon','Background n_e','\psi_N=1 (LFS)','Wall (LFS)'}, 'Location','best');
grid on; set(gca,'FontSize',14);
xlim([2.2845 2.37]);           % your preferred zoom
ylim([1e16 1e18]);             % edge-focused y-range

% Print numeric results
if exist('R_sep','var'),  fprintf('LFS separatrix at R = %.4f m\n', R_sep); end
if exist('R_wall','var'), fprintf('LFS outer wall (from g.lim) at R = %.4f m\n', R_wall); end

% --- Wall-clipped outer bin correction (midplane only, stable) ---
r_in  = r_edges(1:end-1);
r_out = r_edges(2:end);
dZmid = mean(diff(z_edges));           % average cell height (avoids indexing error)

% Clip ring area exactly at the wall radius
r_in_eff  = r_in;
r_out_eff = min(r_out, R_wall);
ring_area_eff = 0.5 * (r_out_eff.^2 - r_in_eff.^2);

% Cells entirely outside the wall -> NaN (not zero)
ring_area_eff(r_in >= R_wall) = NaN;

% Toroidal volume for midplane rings
vol_mid_eff = 2*pi * dZmid * ring_area_eff;     % [1 × nR]

% Recompute midplane GITR carbon only up to wall
counts_mid = data_trim(mid_idx, :);
density_clip = source_strength .* counts_mid .* dt ./ (nP .* vol_mid_eff);

% Clean bad points
density_clip(~isfinite(density_clip)) = NaN;

% Smooth a bit to avoid numerical spikes
density_clip = smoothdata(density_clip, 'movmean', 3);

% Overplot corrected curve
% semilogy(R_mid, density_clip, 'm-', 'LineWidth', 2.5);
legend({'GITR Carbon (orig)','Background Carbon','Background n_e', ...
        '\psi_N=1 (LFS)','Wall (LFS)','GITR Carbon (wall-clipped)'}, ...
        'Location','best');