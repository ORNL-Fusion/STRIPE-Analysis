%% SOLPS-ITER (.mat) -> native (R,Z) -> Euclidean (X,Y) extrapolation -> GITR-style NetCDF
% Uses SOLPS Geo/State quads, paints onto native rectangular grid, adds EFIT B,
% decomposes U|| -> (uR,uZ,uT), then extrapolates to XY using:
%   Option-B clamp on native grid
%   interp2 interior
%   Euclidean KNN-IDW + exp-decay fill ONLY where interp2 is NaN
% Neon enforced like electrons with smooth blend across edge.
%
% IMPORTANT PATCH (flows):
%   Extrapolate U|| on XY with controlled decay+cutoff, then project using b-hat(XY).
%   (Optional) cap U|| to M_cap * c_s for selected species (D+ and all Neon states).
%
% Output: profiles_iter_multi_GITRstyle_XY.nc

close all; clear; clc;

%% ------------------- INPUTS -------------------
matFile   = 'solps_iter.mat';          % Geo/State
eqdskFile = 'MOB-348s_eqdsk.txt';      % used inside read_efit_data
% outnc     = 'profiles_iter_multi_GITRstyle_XY.nc';
outputMode = 'multiFluid';   % 'singleFluid' | 'multiFluid'

%% ------------------- User knobs -------------------
% Native rectangular grid used for painting SOLPS quads
nR_native = 600;      % recommend 400–1200
nZ_native = 400;       % recommend 200–800

doSmooth   = true;
sigma_pix  = 1.0;      % 0.8–1.5

% XY extrapolation grid
num_points     = 500;
R_outer_target = 8.6;

% Midplane window used for slope diagnostics
r_min = 8.01295;
r_max = 8.33246;

% Option-B clamp parameters (native grid)
K_edge = 1;
Ldecay = 0.03;     % [m]
n_min  = 1e10;
T_min  = 3;

% Euclidean fill controls (stripe-free Euclidean method)
K_fill = 12;       % set 1 to recover Voronoi stripes; 4–12 recommended
p_idw  = 2;       % IDW power

% ===== PATCH: FLOW extrapolation controls =====
Lflow_decay   = 0.05;   % [m] decay length for U|| outside SOLPS support (0.02–0.08)
flow_cutoff_m = 0.12;   % [m] beyond this distance from support -> force U|| = 0
Upar_cap      = inf;    % optional cap (e.g. 2e4) to prevent wild values

% ===== PATCH: SOUND-SPEED (Mach) CAP for U|| on XY =====
use_cs_cap  = true;     % cap U|| to M_cap*c_s for D+ and all Neon states
M_cap       = 1.0;      % Mach cap
gamma_e     = 1.0;      % polytropic factors (common choice = 1)
gamma_i     = 1.0;
eC          = 1.602176634e-19;   % J/eV
amu         = 1.66053906660e-27; % kg

% psiN masking (auto-disabled if EFIT g invalid)
apply_psi_mask = true;
psiN_mask = 0.86;

% Neon enforcement controls (2D)
dR_anchor_in   = 0.05;
dR_anchor_win  = 0.005;
% dR_transition  = 0.01;
% dR_transition = 2*(x_xy(2)-x_xy(1));
epsv           = 1e-40;

% Species indexing (your ordering)
iDplus     = 2;     % D+
idxNe_rule = [];    % set after ns known

%% ------------------- LOAD SOLPS .mat -------------------
S = load(matFile);
Geo   = S.Geo;
State = S.State;

r = double(Geo.pr);   % [4 x Nc]
z = double(Geo.pz);   % [4 x Nc]
Nc = size(r,2);

ne_cv = double(State.ne(:));        % [Nc x 1]
na_cv = double(State.na);           % [Nc x ns]
ua_cv = double(State.ua);           % [Nc x ns] parallel flow U||

te_cv = double(State.te(:))./1.602e-19;  % eV
ti_cv = double(State.ti(:))./1.602e-19;  % eV

zn = double(State.zn(:));           % [ns x 1]
am = double(State.am(:));           % [ns x 1]  (amu)
ns = size(na_cv,2);

fprintf('Loaded %d cells and %d species from %s.\n', Nc, ns, matFile);

idxNe_rule = 6:min(16, ns);   % your Ne0..Ne10+ rule (adjust if needed)

%% ------------------- DEFINE NATIVE (R,Z) GRID FOR PAINTING -------------------
rmin = min(r(:)); rmax = max(r(:));
zmin = min(z(:)); zmax = max(z(:));

rS = linspace(rmin, rmax, nR_native);
zS = linspace(zmin, zmax, nZ_native);
[Xn, Zn] = meshgrid(rS, zS);          % [nZ x nR]

% Paint helper
paint = @(vals) paintCells(vals, r, z, Xn, Zn);

%% ------------------- PAINT SOLPS FIELDS ON NATIVE GRID -------------------
disp('Painting SOLPS fields onto native rectangular (R,Z) grid ...');

mask_native = paint(ones(Nc,1)) > 0;

neS = paint(ne_cv);
TeS = paint(te_cv);
TiS = paint(ti_cv);

niS = zeros(nZ_native, nR_native, ns);
UaS = zeros(nZ_native, nR_native, ns);
for s = 1:ns
    niS(:,:,s) = paint(na_cv(:,s));
    UaS(:,:,s) = paint(ua_cv(:,s));
end

% enforce mask -> NaN outside SOLPS mesh
neS(~mask_native) = NaN;
TeS(~mask_native) = NaN;
TiS(~mask_native) = NaN;
for s = 1:ns
    tmp = niS(:,:,s); tmp(~mask_native) = NaN; niS(:,:,s) = tmp;
    tmp = UaS(:,:,s); tmp(~mask_native) = NaN; UaS(:,:,s) = tmp;
end

% optional smoothing (mask-aware)
if doSmooth
    neS = masked_gauss_smooth(neS, mask_native, sigma_pix);
    TeS = masked_gauss_smooth(TeS, mask_native, sigma_pix);
    TiS = masked_gauss_smooth(TiS, mask_native, sigma_pix);
    for s = 1:ns
        niS(:,:,s) = masked_gauss_smooth(niS(:,:,s), mask_native, sigma_pix);
        UaS(:,:,s) = masked_gauss_smooth(UaS(:,:,s), mask_native, sigma_pix);
    end
end

%% ------------------- DATA CLEANING -------------------
neS(~isfinite(neS) | neS<=0) = NaN;
TeS(~isfinite(TeS) | TeS<=0) = NaN;
TiS(~isfinite(TiS) | TiS<=0) = NaN;
TeS(TeS < T_min) = T_min;
TiS(TiS < T_min) = T_min;

%% ------------------- READ EFIT / WALL / PSI -------------------
% Expected from read_efit_data:
%   Br, Bt, Bz on EFIT grid
%   r_efit, z_efit
%   g struct (optional; for wall and psiN)
read_efit_data;

% ---- Defensive checks on EFIT outputs ----
has_g_struct = exist('g','var') && isstruct(g);

hasWall = has_g_struct && isfield(g,'lim') && ~isempty(g.lim) && ...
          size(g.lim,1) >= 2 && size(g.lim,2) >= 3 && ...
          any(isfinite(g.lim(:)));

hasPsi  = has_g_struct && all(isfield(g, {'r','z','psirz','ssimag','ssibry'})) && ...
          ~isempty(g.r) && ~isempty(g.z) && ~isempty(g.psirz);

hasB = exist('Br','var') && exist('Bt','var') && exist('Bz','var') && ...
       exist('r_efit','var') && exist('z_efit','var') && ...
       ~isempty(Br) && ~isempty(Bt) && ~isempty(Bz) && ...
       ~isempty(r_efit) && ~isempty(z_efit);

if ~hasPsi
    warning('EFIT g missing/invalid -> disabling psiN mask.');
    apply_psi_mask = false;
else
    if ~isfield(g,'dR') || isempty(g.dR)
        g.dR = g.r(2) - g.r(1);
    end
    if ~isfield(g,'dZ') || isempty(g.dZ)
        g.dZ = g.z(2) - g.z(1);
    end
end

% ---- Sanitize wall polygon once, use everywhere ----
hasWall_poly = false;
r_wall = [];
z_wall = [];

if hasWall
    r_wall = double(g.lim(1,:));
    z_wall = double(g.lim(2,:));

    okw = isfinite(r_wall) & isfinite(z_wall);
    r_wall = r_wall(okw);
    z_wall = z_wall(okw);

    if numel(r_wall) >= 3
        % close polygon if needed
        if r_wall(1) ~= r_wall(end) || z_wall(1) ~= z_wall(end)
            r_wall(end+1) = r_wall(1);
            z_wall(end+1) = z_wall(1);
        end
        hasWall_poly = true;
    else
        r_wall = [];
        z_wall = [];
        hasWall_poly = false;
    end
end

%% ------------------- MAP EFIT B TO NATIVE GRID -------------------
if ~hasB
    warning('EFIT B-field missing/invalid -> using zero B-field everywhere.');
    BrS = zeros(size(Xn));
    BtS = zeros(size(Xn));
    BzS = zeros(size(Xn));
else
    FBr = griddedInterpolant({double(z_efit), double(r_efit)}, double(Br), 'linear', 'nearest');
    FBt = griddedInterpolant({double(z_efit), double(r_efit)}, double(Bt), 'linear', 'nearest');
    FBz = griddedInterpolant({double(z_efit), double(r_efit)}, double(Bz), 'linear', 'nearest');

    BrS = FBr(Zn, Xn);   % [nZ x nR]
    BtS = FBt(Zn, Xn);
    BzS = FBz(Zn, Xn);

    BrS(~isfinite(BrS)) = 0;
    BtS(~isfinite(BtS)) = 0;
    BzS(~isfinite(BzS)) = 0;
end

% ---- Unit vectors on native grid ----
Bmag = sqrt(BrS.^2 + BtS.^2 + BzS.^2);
epsB = 1e-30;
bRhat = BrS ./ max(Bmag, epsB);
bZhat = BzS ./ max(Bmag, epsB);
bThat = BtS ./ max(Bmag, epsB);

%% ------------------- DECOMPOSE U|| -> (uR,uZ,uT) ON NATIVE GRID -------------------
% For diagnostics only
uR_S = zeros(nZ_native, nR_native, ns);
uZ_S = zeros(nZ_native, nR_native, ns);
uT_S = zeros(nZ_native, nR_native, ns);

for s = 1:ns
    Upar = UaS(:,:,s);
    uR_S(:,:,s) = Upar .* bRhat;
    uZ_S(:,:,s) = Upar .* bZhat;
    uT_S(:,:,s) = Upar .* bThat;
end

%% ------------------- QUICK SANITY ON NATIVE GRID -------------------
figure;
imagesc(rS, zS, neS);
set(gca,'YDir','normal');
colorbar;
title('n_e (painted native grid)');
xlabel('R [m]');
ylabel('Z [m]');
hold on;
if hasWall_poly
    plot(r_wall, z_wall, 'r', 'LineWidth', 1.2);
end
hold off;

%% ------------------- NATIVE MESH + MASKS -------------------
plasma_mask = isfinite(neS) & (neS > n_min);

inside_limiter_mask = true(size(plasma_mask));
if hasWall_poly
    inside_limiter_mask = inpolygon(Xn, Zn, r_wall, z_wall);
end

%% ------------------- OMP edge at Z≈0 -------------------
[~, iz0] = min(abs(zS - 0));

edge_idx_mid = find(plasma_mask(iz0,:), 1, 'last');
if isempty(edge_idx_mid)
    edge_idx_mid = find(isfinite(neS(iz0,:)) & neS(iz0,:) > 0, 1, 'last');
end
if isempty(edge_idx_mid)
    error('No valid SOLPS midplane edge point found at Z≈0 on painted grid.');
end

R_edge_omp = rS(edge_idx_mid);
fprintf('OMP edge (midplane): R_edge_omp = %.5f m\n', R_edge_omp);

%% ------------------- XY extrapolation grid -------------------
% Extend slightly beyond the outer wall rather than using a hard-coded target
dR_wall_buffer = 0.02;   % [m], tune 0.01–0.05

% --- XY grid definition (force outer radius to R_outer_target) ---
if hasWall_poly
    R_xy_min = min(r_wall);
    Z_xy_min = min(z_wall);
    Z_xy_max = max(z_wall);
else
    R_xy_min = min(rS);
    Z_xy_min = min(zS);
    Z_xy_max = max(zS);
end

% Force grid to extend to desired radius
R_xy_max = R_outer_target;
[X, Y] = meshgrid(linspace(R_xy_min, R_xy_max, num_points), ...
                  linspace(Z_xy_min, Z_xy_max, num_points));   % [nZout x nRout]

fprintf('XY grid: R=[%.3f, %.3f], Z=[%.3f, %.3f], N=%dx%d\n', ...
    R_xy_min, R_xy_max, Z_xy_min, Z_xy_max, size(X,1), size(X,2));

x_xy = X(1,:);
z_xy = Y(:,1);
[nZout, nRout] = size(X);

dx_xy = x_xy(2) - x_xy(1);
dR_transition = 2 * dx_xy;   % 2 grid-cell blend width
extrap_coords = [X(:), Y(:)];

%% ------------------- SOLPS-derived decay (Option B) + XY extrapolation -------------------
% Strategy:
%   1) Build 1D midplane diagnostics only for checking fits.
%   2) On native (R,Z), Option-B modifies ONLY points beyond the last SOLPS cell.
%   3) On XY, use interp2 wherever SOLPS support exists.
%   4) Only where interp2 is NaN, use Euclidean KNN-IDW + exp-decay fill.

%% =========================
%% (A) Midplane diagnostics + fallback slopes
%% =========================
mpfx = linspace(r_min, r_max, 1000);
mpx_max = max(R_xy_max, R_outer_target);
mpx  = linspace(r_min, mpx_max, 1000);

% define clamp bounds before use
Lmin = 0.005;   % [m]
Lmax = 0.20;    % [m]

densityAtMidplane = interp2(rS, zS, neS, mpx, 0*mpx, 'linear', NaN);
TeAtMidplane      = interp2(rS, zS, TeS, mpx, 0*mpx, 'linear', NaN);
TiAtMidplane      = interp2(rS, zS, TiS, mpx, 0*mpx, 'linear', NaN);

fitDensityAtMidplane = interp2(rS, zS, neS, mpfx, 0*mpfx, 'linear', NaN);
fitTeAtMidplane      = interp2(rS, zS, TeS, mpfx, 0*mpfx, 'linear', NaN);
fitTiAtMidplane      = interp2(rS, zS, TiS, mpfx, 0*mpfx, 'linear', NaN);

if nnz(isfinite(fitDensityAtMidplane) & fitDensityAtMidplane > 0) >= 3
    p_ne = polyfit(mpfx(isfinite(fitDensityAtMidplane) & fitDensityAtMidplane > 0), ...
                   log(fitDensityAtMidplane(isfinite(fitDensityAtMidplane) & fitDensityAtMidplane > 0)), 1);
else
    p_ne = [-10, log(n_min)];
end

if nnz(isfinite(fitTeAtMidplane) & fitTeAtMidplane > 0) >= 3
    p_Te = polyfit(mpfx(isfinite(fitTeAtMidplane) & fitTeAtMidplane > 0), ...
                   log(fitTeAtMidplane(isfinite(fitTeAtMidplane) & fitTeAtMidplane > 0)), 1);
else
    p_Te = [-5, log(max(T_min,1))];
end

if nnz(isfinite(fitTiAtMidplane) & fitTiAtMidplane > 0) >= 3
    p_Ti = polyfit(mpfx(isfinite(fitTiAtMidplane) & fitTiAtMidplane > 0), ...
                   log(fitTiAtMidplane(isfinite(fitTiAtMidplane) & fitTiAtMidplane > 0)), 1);
else
    p_Ti = p_Te;
end

% slopes used only for outer Euclidean fill
s_ne = min(p_ne(1), -1e-4);
s_Te = min(p_Te(1), -1e-4);
s_Ti = min(p_Ti(1), -1e-4);

Lne_mid = max(Lmin, min(Lmax, -1 / s_ne));
LTe_mid = max(Lmin, min(Lmax, -1 / s_Te));
LTi_mid = max(Lmin, min(Lmax, -1 / s_Ti));

extrapolatedne1d = NaN(size(mpx));
extrapolatedTe1d = NaN(size(mpx));
extrapolatedTi1d = NaN(size(mpx));

idx_last_ne = find(isfinite(densityAtMidplane) & densityAtMidplane > 0, 1, 'last');
idx_last_Te = find(isfinite(TeAtMidplane) & TeAtMidplane > 0, 1, 'last');
idx_last_Ti = find(isfinite(TiAtMidplane) & TiAtMidplane > 0, 1, 'last');

if ~isempty(idx_last_ne)
    extrapolatedne1d(1:idx_last_ne) = densityAtMidplane(1:idx_last_ne);
    R0 = mpx(idx_last_ne);
    n0 = max(densityAtMidplane(idx_last_ne), n_min);
    beyond = mpx > R0;
    extrapolatedne1d(beyond) = n0 .* exp(-(mpx(beyond) - R0) / Lne_mid);
end

if ~isempty(idx_last_Te)
    extrapolatedTe1d(1:idx_last_Te) = TeAtMidplane(1:idx_last_Te);
    R0 = mpx(idx_last_Te);
    T0 = max(TeAtMidplane(idx_last_Te), T_min);
    beyond = mpx > R0;
    extrapolatedTe1d(beyond) = max(T0 .* exp(-(mpx(beyond) - R0) / LTe_mid), T_min);
end

if ~isempty(idx_last_Ti)
    extrapolatedTi1d(1:idx_last_Ti) = TiAtMidplane(1:idx_last_Ti);
    R0 = mpx(idx_last_Ti);
    T0 = max(TiAtMidplane(idx_last_Ti), T_min);
    beyond = mpx > R0;
    extrapolatedTi1d(beyond) = max(T0 .* exp(-(mpx(beyond) - R0) / LTi_mid), T_min);
end

figure;
semilogy(mpx, max(extrapolatedne1d,1), 'b-', 'LineWidth', 1.8); hold on;
semilogy(mpx, max(extrapolatedTe1d,1), 'r-', 'LineWidth', 1.8);
semilogy(mpx, max(extrapolatedTi1d,1), 'k-', 'LineWidth', 1.8);
semilogy(mpx, max(densityAtMidplane,1), 'b.', 'MarkerSize', 8);
semilogy(mpx, max(TeAtMidplane,1), 'r.', 'MarkerSize', 8);
semilogy(mpx, max(TiAtMidplane,1), 'k.', 'MarkerSize', 8);
xline(R_edge_omp, 'k--', 'LineWidth', 1.2);
legend('n_e extrap','T_e extrap','T_i extrap', ...
       'n_e midplane','T_e midplane','T_i midplane','R_{edge,OMP}');
xlabel('R [m]');
ylabel('Value');
title('ITER midplane diagnostics: native + fitted extrapolation');
grid on; box on;

%% =========================
%% (B) Option-B clamp on native grid using SOLPS local gradients: L(z)
%%     ONLY beyond last SOLPS-supported cell
%%     PATCH: smooth L(z) before applying clamp
%% =========================
neS_clamp = neS;
TeS_clamp = TeS;
TiS_clamp = TiS;

Ngrad_in = 6;

Lne_z = nan(nZ_native,1);
LTe_z = nan(nZ_native,1);
LTi_z = nan(nZ_native,1);

for jz = 1:nZ_native
    row = plasma_mask(jz,:);
    edge_idx = find(row, 1, 'last');
    if isempty(edge_idx) || edge_idx < 4
        continue;
    end

    ie = max(3, edge_idx - 1);
    i1 = max(2, ie - Ngrad_in);
    i2 = ie;

    Lne_z(jz) = local_decay_from_grad(rS(i1:i2), neS(jz,i1:i2), Lmin, Lmax);
    LTe_z(jz) = local_decay_from_grad(rS(i1:i2), TeS(jz,i1:i2), Lmin, Lmax);
    LTi_z(jz) = local_decay_from_grad(rS(i1:i2), TiS(jz,i1:i2), Lmin, Lmax);

    if ~isfinite(Lne_z(jz)), Lne_z(jz) = Lne_mid; end
    if ~isfinite(LTe_z(jz)), LTe_z(jz) = LTe_mid; end
    if ~isfinite(LTi_z(jz)), LTi_z(jz) = LTi_mid; end
end

% fill missing rows
good = find(isfinite(Lne_z));
if numel(good) >= 2
    Lne_z = interp1(good, Lne_z(good), (1:nZ_native)', 'nearest', 'extrap');
else
    Lne_z(:) = Lne_mid;
end

good = find(isfinite(LTe_z));
if numel(good) >= 2
    LTe_z = interp1(good, LTe_z(good), (1:nZ_native)', 'nearest', 'extrap');
else
    LTe_z(:) = LTe_mid;
end

good = find(isfinite(LTi_z));
if numel(good) >= 2
    LTi_z = interp1(good, LTi_z(good), (1:nZ_native)', 'nearest', 'extrap');
else
    LTi_z(:) = LTi_mid;
end

% smooth L(z) to suppress row-to-row patchiness
Lsmooth_win = 7;   % try 5, 7, or 9
kerL = ones(Lsmooth_win,1) / Lsmooth_win;

Lne_z = conv(Lne_z, kerL, 'same');
LTe_z = conv(LTe_z, kerL, 'same');
LTi_z = conv(LTi_z, kerL, 'same');

Lne_z = max(Lmin, min(Lmax, Lne_z));
LTe_z = max(Lmin, min(Lmax, LTe_z));
LTi_z = max(Lmin, min(Lmax, LTi_z));

% now apply outward clamp row-by-row
for jz = 1:nZ_native
    row = plasma_mask(jz,:);
    edge_idx = find(row, 1, 'last');
    if isempty(edge_idx) || edge_idx < 1
        continue;
    end

    i0 = max(1, edge_idx - K_edge + 1);

    edge_ne = max(neS(jz, i0:edge_idx), [], 'omitnan');
    edge_Te = max(TeS(jz, i0:edge_idx), [], 'omitnan');
    edge_Ti = max(TiS(jz, i0:edge_idx), [], 'omitnan');

    if ~isfinite(edge_ne), edge_ne = neS(jz, edge_idx); end
    if ~isfinite(edge_Te), edge_Te = max(TeS(jz, edge_idx), T_min); end
    if ~isfinite(edge_Ti), edge_Ti = max(TiS(jz, edge_idx), T_min); end

    if ~isfinite(edge_ne), continue; end

    for ir = edge_idx+1:nR_native
        dr = rS(ir) - rS(edge_idx);
        neS_clamp(jz,ir) = max(edge_ne * exp(-dr / Lne_z(jz)), 0);
        TeS_clamp(jz,ir) = max(edge_Te * exp(-dr / LTe_z(jz)), T_min);
        TiS_clamp(jz,ir) = max(edge_Ti * exp(-dr / LTi_z(jz)), T_min);
    end
end

fprintf('Option-B SOLPS-grad L(z) [smoothed]: Lne=[%.3f, %.3f] m, LTe=[%.3f, %.3f] m, LTi=[%.3f, %.3f] m\n', ...
    min(Lne_z), max(Lne_z), min(LTe_z), max(LTe_z), min(LTi_z), max(LTi_z));

%% =========================
%% (C) Mapping to XY for electrons
%%     interp2 inside support, edge-band Euclidean fill only outside
%%     PATCH: use outer-anchor band, not all native valid points
%% =========================
coords_ok = [Xn(:), Zn(:)];

% preserve SOLPS interior first
val_ne_interp = interp2(rS, zS, neS, X, Y, 'linear', NaN);
val_Te_interp = interp2(rS, zS, TeS, X, Y, 'linear', NaN);
val_Ti_interp = interp2(rS, zS, TiS, X, Y, 'linear', NaN);

% ---- build outer-edge anchor mask ----
Nanchor = 3;   % number of last valid native cells per row used as anchors
outer_anchor_mask = false(size(plasma_mask));

for jz = 1:nZ_native
    edge_idx = find(plasma_mask(jz,:), 1, 'last');
    if isempty(edge_idx)
        continue;
    end
    i1 = max(1, edge_idx - Nanchor + 1);
    outer_anchor_mask(jz, i1:edge_idx) = true;
end

% anchor points for scalar outer fill
ok_ne = outer_anchor_mask & isfinite(neS_clamp) & (neS_clamp > 0);
ok_Te = outer_anchor_mask & isfinite(TeS_clamp) & (TeS_clamp > 0);
ok_Ti = outer_anchor_mask & isfinite(TiS_clamp) & (TiS_clamp > 0);

idx_ne = find(ok_ne(:));
idx_Te = find(ok_Te(:));
idx_Ti = find(ok_Ti(:));

coords_ne = coords_ok(idx_ne,:); vals_ne = neS_clamp(idx_ne);
coords_Te = coords_ok(idx_Te,:); vals_Te = TeS_clamp(idx_Te);
coords_Ti = coords_ok(idx_Ti,:); vals_Ti = TiS_clamp(idx_Ti);

num_points_all = numel(X);
pool = gcp('nocreate');
if isempty(pool), pool = parpool; end
num_workers = pool.NumWorkers;

batch_size = ceil(num_points_all / (10 * num_workers));
extrap_coords = [X(:), Y(:)];
nBatches = ceil(num_points_all / batch_size);

% parfor-safe locals
coords_ne_par = coords_ne;
coords_Te_par = coords_Te;
coords_Ti_par = coords_Ti;
vals_ne_par   = vals_ne;
vals_Te_par   = vals_Te;
vals_Ti_par   = vals_Ti;

s_ne_par = s_ne;
s_Te_par = s_Te;
s_Ti_par = s_Ti;
K_fill_par = K_fill;
p_idw_par  = p_idw;
T_min_par  = T_min;

val_ne_cell = cell(nBatches,1);
val_Te_cell = cell(nBatches,1);
val_Ti_cell = cell(nBatches,1);

parfor b = 1:nBatches
    ib = (b-1)*batch_size + 1;
    ie = min(b*batch_size, num_points_all);
    subset = extrap_coords(ib:ie,:);

    % ne
    k_ne = min(K_fill_par, size(coords_ne_par,1));
    [iN, dN] = knnsearch(coords_ne_par, subset, 'K', k_ne);
    if k_ne == 1
        iN = iN(:); dN = dN(:);
        ne_out = vals_ne_par(iN) .* exp(s_ne_par .* dN);
    else
        wN = 1 ./ max(dN, 1e-12).^p_idw_par;
        wN = wN ./ sum(wN, 2);
        vN = vals_ne_par(iN) .* exp(s_ne_par .* dN);
        ne_out = sum(wN .* vN, 2);
    end

    % Te
    k_te = min(K_fill_par, size(coords_Te_par,1));
    [iT, dT] = knnsearch(coords_Te_par, subset, 'K', k_te);
    if k_te == 1
        iT = iT(:); dT = dT(:);
        Te_out = vals_Te_par(iT) .* exp(s_Te_par .* dT);
    else
        wT = 1 ./ max(dT, 1e-12).^p_idw_par;
        wT = wT ./ sum(wT, 2);
        vT = vals_Te_par(iT) .* exp(s_Te_par .* dT);
        Te_out = sum(wT .* vT, 2);
    end
    Te_out = max(Te_out, T_min_par);

    % Ti
    k_ti = min(K_fill_par, size(coords_Ti_par,1));
    [iI, dI] = knnsearch(coords_Ti_par, subset, 'K', k_ti);
    if k_ti == 1
        iI = iI(:); dI = dI(:);
        Ti_out = vals_Ti_par(iI) .* exp(s_Ti_par .* dI);
    else
        wI = 1 ./ max(dI, 1e-12).^p_idw_par;
        wI = wI ./ sum(wI, 2);
        vI = vals_Ti_par(iI) .* exp(s_Ti_par .* dI);
        Ti_out = sum(wI .* vI, 2);
    end
    Ti_out = max(Ti_out, T_min_par);

    val_ne_cell{b} = ne_out;
    val_Te_cell{b} = Te_out;
    val_Ti_cell{b} = Ti_out;
end

val_ne_extrap = reshape(cell2mat(val_ne_cell), size(X));
val_Te_extrap = reshape(cell2mat(val_Te_cell), size(X));
val_Ti_extrap = reshape(cell2mat(val_Ti_cell), size(X));

% keep interp2 values where available, fill only where missing
val_ne = val_ne_interp;
val_Te = val_Te_interp;
val_Ti = val_Ti_interp;

mask_ne_out = isnan(val_ne_interp);
mask_Te_out = isnan(val_Te_interp);
mask_Ti_out = isnan(val_Ti_interp);

val_ne(mask_ne_out) = val_ne_extrap(mask_ne_out);
val_Te(mask_Te_out) = val_Te_extrap(mask_Te_out);
val_Ti(mask_Ti_out) = val_Ti_extrap(mask_Ti_out);
%% =========================

%% ------------------- Map B-fields to XY -------------------
Br_xy = interp2(rS, zS, BrS, X, Y, 'linear', 0);
Bt_xy = interp2(rS, zS, BtS, X, Y, 'linear', 0);
Bz_xy = interp2(rS, zS, BzS, X, Y, 'linear', 0);

Bmag_xy = sqrt(Br_xy.^2 + Bt_xy.^2 + Bz_xy.^2);
bR_xy = Br_xy ./ max(Bmag_xy, 1e-30);
bZ_xy = Bz_xy ./ max(Bmag_xy, 1e-30);
bT_xy = Bt_xy ./ max(Bmag_xy, 1e-30);


%% ------------------- Multi-species mapping (ni + velocities) -------------------
% NEON-AWARE REVISED BLOCK
% Goals:
%  - Non-Ne species: Option-B native clamp + edge-anchor XY mapping (as before)
%  - Neon: map total Ne smoothly, renormalize to follow val_ne outside support,
%          reconstruct charge states from smooth Z-dependent fractions
%  - Apply light 2D smoothing to suppress striping
%  - Freeze last radial columns to avoid artificial drop at outermost grid

val_ni_mspecies = zeros(nZout, nRout, ns);
val_uR_mspecies = zeros(nZout, nRout, ns);
val_uZ_mspecies = zeros(nZout, nRout, ns);
val_uT_mspecies = zeros(nZout, nRout, ns);

src_mask_u = plasma_mask & inside_limiter_mask;   % native support mask for flows

% ---- support distance on XY (for "outside SOLPS" detection) ----
sup_native = isfinite(neS) & (neS > 0);
idx_sup = find(sup_native(:));
coords_sup = [Xn(idx_sup), Zn(idx_sup)];

iSup = knnsearch(coords_sup, extrap_coords, 'K', 1);
dSup = vecnorm(coords_sup(iSup,:) - extrap_coords, 2, 2);
dSup = reshape(dSup, size(X));                     % [nZout x nRout]

Lblend = max(dR_transition, 2*(x_xy(2)-x_xy(1)));   % smooth transition thickness
% corrected outside-support blending: 0 at support, ->1 outside
w_out = 1 - exp(-(dSup ./ max(Lblend,1e-6)).^2);

% ---- sanitized wall polygon for inpolygon (prevents r_wall undefined crash) ----
hasWall_poly = exist('g','var') && isstruct(g) && isfield(g,'lim') && ~isempty(g.lim) ...
               && size(g.lim,1) >= 2 && size(g.lim,2) >= 3;
if hasWall_poly
    r_wall = g.lim(1,:);
    z_wall = g.lim(2,:);
else
    r_wall = [];
    z_wall = [];
end

% ---- Neon helper ----
isNe = false(1,ns);
if ~isempty(idxNe_rule)
    isNe(idxNe_rule) = true;
end

% ---- outer-edge anchor mask (same idea as electrons) ----
Nanchor_species = 3;   % try 3 or 4
outer_anchor_mask_species = false(size(plasma_mask));
for jz = 1:nZ_native
    edge_idx = find(plasma_mask(jz,:), 1, 'last');
    if isempty(edge_idx), continue; end
    i1 = max(1, edge_idx - Nanchor_species + 1);
    outer_anchor_mask_species(jz, i1:edge_idx) = true;
end

% ---- Precompute fNeZ(Z): total-Ne / ne at an anchor band just inside edge ----
fNeZ = zeros(nZout,1);
fNe_max = 1.2;  % cap total Neon fraction of ne (tune 0.01–0.10)

if any(isNe)
    for iz = 1:nZout
        Zq = z_xy(iz);
        [~, jz0] = min(abs(zS - Zq));

        row = plasma_mask(jz0,:);
        edge_idx = find(row, 1, 'last');
        if isempty(edge_idx)
            fNeZ(iz) = 0;
            continue
        end

        Ranchor = rS(edge_idx) - dR_anchor_in;
        RL = Ranchor - dR_anchor_win;
        RU = Ranchor + dR_anchor_win;

        iL = find(rS >= RL, 1, 'first'); if isempty(iL), iL = 1; end
        iU = find(rS <= RU, 1, 'last');  if isempty(iU), iU = nR_native; end
        if iU < iL
            [~,iC] = min(abs(rS - Ranchor));
            iL = iC; iU = iC;
        end

        ne_anchor = median(neS_clamp(jz0,iL:iU), 'omitnan');

        nNe_anchor = 0;
        for j = 1:numel(idxNe_rule)
            kNe = idxNe_rule(j);
            tmp = squeeze(niS(jz0,iL:iU,kNe));
            nNe_anchor = nNe_anchor + median(tmp(:), 'omitnan');
        end

        if ~isfinite(ne_anchor) || ne_anchor <= 0 || ~isfinite(nNe_anchor) || nNe_anchor <= 0
            fNeZ(iz) = 0;
        else
            fNeZ(iz) = nNe_anchor / max(ne_anchor, epsv);
        end
    end

    % cap and smooth vertical behaviour
    fNeZ = min(fNeZ, fNe_max);
    if nZout >= 9
        fNeZ = smoothdata(fNeZ, 'movmean', 11);   % reduces Z striping
    end
    fNeZ = max(fNeZ, 0);
end

% ----------------------------
% MAIN: non-Ne densities + flows
% ----------------------------
for k = 1:ns
    ni_k = niS(:,:,k);

    % cap only for D+ and Neon
    apply_cap_this_species = (k == iDplus) || (any(isNe) && isNe(k));

    % DENSITY: skip individual Neon charge-states here (handled later)
    if ~(any(isNe) && isNe(k))
        % --- outward clamp on native (Option-B style using local gradients) ---
        ni_ex = ni_k;
        Lni_z = nan(nZ_native,1);
        for jz = 1:nZ_native
            row = plasma_mask(jz,:);
            edge_idx = find(row, 1, 'last');
            if isempty(edge_idx) || edge_idx < 4, continue; end
            ie = max(3, edge_idx - 1);
            i1 = max(2, ie - Ngrad_in);
            i2 = ie;
            Lni_z(jz) = local_decay_from_grad(rS(i1:i2), ni_k(jz,i1:i2), Lmin, Lmax);
        end
        good = find(isfinite(Lni_z));
        if numel(good) >= 2
            Lni_z = interp1(good, Lni_z(good), (1:nZ_native)', 'nearest', 'extrap');
        elseif any(isfinite(Lni_z))
            Lni_z(~isfinite(Lni_z)) = median(Lni_z(isfinite(Lni_z)));
        else
            Lni_z(:) = Ldecay;
        end
        Lni_z = conv(Lni_z, ones(7,1)/7, 'same');    % light smoothing
        Lni_z = max(Lmin, min(Lmax, Lni_z));

        for jz = 1:nZ_native
            row = plasma_mask(jz,:);
            edge_idx = find(row, 1, 'last');
            if isempty(edge_idx), continue; end
            i0 = max(1, edge_idx - K_edge + 1);
            e_ni = max(ni_k(jz, i0:edge_idx), [], 'omitnan');
            if ~isfinite(e_ni), e_ni = ni_k(jz,edge_idx); end
            if ~isfinite(e_ni), continue; end
            for ir = edge_idx+1:nR_native
                dr = rS(ir) - rS(edge_idx);
                ni_ex(jz,ir) = max(e_ni * exp(-dr / Lni_z(jz)), 0);
            end
        end

        % --- XY mapping for non-Ne species (interior via interp2, outside via KNN-IDW+decay) ---
        val_ni_interp = interp2(rS, zS, ni_k, X, Y, 'linear', NaN);
        ok_ni = outer_anchor_mask_species & isfinite(ni_ex) & (ni_ex > 0);
        idx_ni = find(ok_ni(:));

        if isempty(idx_ni)
            val_ni_extrap = zeros(size(X));
        else
            coords_ni = coords_ok(idx_ni,:);
            vals_ni   = ni_ex(idx_ni);

            ni_mid = interp2(rS, zS, ni_k, mpfx, 0*mpfx, 'linear', NaN);
            if nnz(isfinite(ni_mid) & ni_mid > 0) >= 3
                p_ni = polyfit(mpfx(isfinite(ni_mid) & ni_mid > 0), ...
                               log(ni_mid(isfinite(ni_mid) & ni_mid > 0)), 1);
            else
                p_ni = p_ne;
            end
            s_use = min(p_ni(1), -1e-4);

            K_use = min(max(4, K_fill), size(coords_ni,1));
            p_use = p_idw;

            ni_cell = cell(nBatches,1);
            parfor b = 1:nBatches
                s = (b-1)*batch_size + 1;
                e = min(b*batch_size, num_points_all);
                subset = extrap_coords(s:e,:);
                [iN, dN] = knnsearch(coords_ni, subset, 'K', K_use);
                if K_use == 1
                    iN = iN(:); dN = dN(:);
                    ni_out = vals_ni(iN) .* exp(s_use .* dN);
                else
                    w = 1 ./ max(dN, 1e-12).^p_use; w = w ./ sum(w, 2);
                    valsK  = vals_ni(iN); decayK = exp(s_use .* dN);
                    ni_out = sum(w .* valsK .* decayK, 2);
                end
                ni_cell{b} = ni_out;
            end
            val_ni_extrap = reshape(cell2mat(ni_cell), size(X));
        end

        val_ni = val_ni_interp;
        mask_ni_out = isnan(val_ni_interp) | (val_ni_interp <= 0);
        val_ni(mask_ni_out) = val_ni_extrap(mask_ni_out);

        % light 2D smoothing (reduces checkerboard/row-striping when increasing resolution)
        if exist('imgaussfilt','file') == 2
            val_ni = imgaussfilt(val_ni, 0.8);
        else
            val_ni = local_imgaussfilt_fallback(val_ni, 0.8);
        end

        val_ni_mspecies(:,:,k) = max(val_ni, 0);
    end

    % ---------------------------
    % =====================================================================
    % (3) FLOWS: preserve interp2 interior, fill ONLY outside support
    %     using edge-anchor Euclidean KNN-IDW, then decay with distance
    % =====================================================================
    if any(src_mask_u(:))
        Upar_native = UaS(:,:,k);

        % --- interior support for flows ---
        supU = src_mask_u & isfinite(Upar_native);

        % interp2 keeps interior smooth
        Upar_interp = interp2(rS, zS, Upar_native, X, Y, 'linear', NaN);

        % ---- build outer-edge anchors for flows (same idea as densities) ----
        outer_anchor_mask_u = false(size(src_mask_u));
        Nanchor_u = 3;   % try 3 or 4

        for jz = 1:nZ_native
            edge_idx = find(src_mask_u(jz,:), 1, 'last');
            if isempty(edge_idx)
                continue;
            end
            i1 = max(1, edge_idx - Nanchor_u + 1);
            outer_anchor_mask_u(jz, i1:edge_idx) = true;
        end

        supU_anchor = outer_anchor_mask_u & isfinite(Upar_native);
        idxU = find(supU_anchor(:));

        if isempty(idxU)
            Upar_xy = zeros(size(X));
        else
            coordsU = [Xn(idxU), Zn(idxU)];
            valsU   = Upar_native(idxU);

            % start from smooth interior
            Upar_xy = Upar_interp;

            % fill only where interp2 failed
            fillU = ~isfinite(Upar_interp);

            if any(fillU(:))
                subset = extrap_coords(fillU(:), :);

                K_use_u = min(max(4, K_fill), size(coordsU,1));
                p_use_u = p_idw;

                [iU, dU] = knnsearch(coordsU, subset, 'K', K_use_u);

                if K_use_u == 1
                    iU = iU(:);
                    dU = dU(:);
                    Ufill = valsU(iU);
                    dmin  = dU;
                else
                    wU = 1 ./ max(dU, 1e-12).^p_use_u;
                    wU = wU ./ sum(wU, 2);
                    valsUK = valsU(iU);
                    Ufill  = sum(wU .* valsUK, 2);
                    dmin   = dU(:,1);
                end

                % decay only outside support
                Lflow_use = max(Lflow_decay, 1e-6);
                Ufill = Ufill .* exp(-dmin ./ Lflow_use);

                if isfinite(flow_cutoff_m) && flow_cutoff_m > 0
                    Ufill(dmin >= flow_cutoff_m) = 0;
                end

                Upar_xy(fillU) = Ufill;
            end
        end

        % optional light smoothing ONLY on filled region to suppress stripes
        fillU_all = ~isfinite(Upar_interp);
        if any(fillU_all(:))
            Utmp = Upar_xy;
            if exist('imgaussfilt','file') == 2
                Usm = imgaussfilt(Utmp, 0.8);   % try 0.6–1.0
            else
                Usm = local_imgaussfilt_fallback(Utmp, 0.8);
            end
            Upar_xy(fillU_all) = Usm(fillU_all);
        end

        % suppress outside limiter polygon safely
        if hasWall_poly
            inside_xy = inpolygon(X, Y, r_wall, z_wall);
            Upar_xy(~inside_xy) = 0;
        end

        Upar_xy(~isfinite(Upar_xy)) = 0;

        % hard cap
        if isfinite(Upar_cap)
            Upar_xy = max(min(Upar_xy, Upar_cap), -Upar_cap);
        end

        % Mach cap for selected species
        if use_cs_cap && apply_cap_this_species
            m_i = max(am(k) * amu, 1e-40);
            cs_xy = sqrt(max((gamma_e * val_Te + gamma_i * val_Ti) * eC ./ m_i, 0));
            Ulim = M_cap .* cs_xy;
            Upar_xy = max(min(Upar_xy, Ulim), -Ulim);
        end

        % project onto XY magnetic-field unit vectors
        val_uR_mspecies(:,:,k) = Upar_xy .* bR_xy;
        val_uZ_mspecies(:,:,k) = Upar_xy .* bZ_xy;
        val_uT_mspecies(:,:,k) = Upar_xy .* bT_xy;
    end
end
%% =============================================================================
%% === REBUILD NEON: keep SOLPS charge-state structure inside support,
%% === but outside support force total-Ne to follow electron fall smoothly
%% =============================================================================
if any(isNe)
    fprintf('\n=== Rebuilding Neon: SOLPS fractions inside, electron-fall outside ===\n');

    % ---------------------------------------------------------------------
    % 1) Native total Neon from SOLPS
    % ---------------------------------------------------------------------
    nNe_native = zeros(size(neS));
    for j = 1:numel(idxNe_rule)
        kNe = idxNe_rule(j);
        tmp = niS(:,:,kNe);
        tmp(~isfinite(tmp)) = 0;
        nNe_native = nNe_native + tmp;
    end
    nNe_native(~isfinite(nNe_native)) = 0;

    % ---------------------------------------------------------------------
    % 2) Map TOTAL Neon from SOLPS to XY (inside structure)
    % ---------------------------------------------------------------------
    val_nNe_in = interp2(rS, zS, nNe_native, X, Y, 'linear', NaN);
    val_nNe_in(~isfinite(val_nNe_in)) = 0;

    % optional very light smoothing of the interior only
    val_nNe_in = local_imgaussfilt_fallback(val_nNe_in, 0.5);
    val_nNe_in = max(val_nNe_in, 0);

    % ---------------------------------------------------------------------
    % 3) Build OUTSIDE target directly from electrons
    %    This guarantees Neon follows electron fall beyond support
    % ---------------------------------------------------------------------
    fNeZ_use = fNeZ;
    good_f = isfinite(fNeZ_use) & (fNeZ_use > 0);

    if any(good_f)
        fref = median(fNeZ_use(good_f), 'omitnan');
    else
        fref = 0.01;
    end

    fNeZ_use(~isfinite(fNeZ_use) | fNeZ_use < 0) = fref;

    if nZout >= 9
        fNeZ_use = smoothdata(fNeZ_use, 'movmean', 21);
        fNeZ_use = smoothdata(fNeZ_use, 'gaussian', 15);
    end
    fNeZ_use = max(fNeZ_use, 0);

    val_nNe_out = fNeZ_use .* val_ne;

    % cap Ne/ne if desired
    fNe_max = 0.05;
    val_nNe_out = min(val_nNe_out, fNe_max .* max(val_ne,0));
    val_nNe_out(~isfinite(val_nNe_out)) = 0;

    % ---------------------------------------------------------------------
    % 4) Smooth blend: inside uses SOLPS Neon, outside uses electron-following Neon
    %    w_out = 0 at support, -> 1 outside support
    % ---------------------------------------------------------------------
    val_nNe_tot = (1 - w_out) .* val_nNe_in + w_out .* val_nNe_out;
    val_nNe_tot = max(val_nNe_tot, 0);

    % light smoothing to suppress any remaining banding
    val_nNe_tot = local_imgaussfilt_fallback(val_nNe_tot, 0.7);
    val_nNe_tot = max(val_nNe_tot, 0);

    % ---------------------------------------------------------------------
    % 5) Build SOLPS-based Neon charge fractions on native grid
    %    Keep spatial charge-state structure from SOLPS
    % ---------------------------------------------------------------------
    fq_native = zeros(nZ_native, nR_native, numel(idxNe_rule));

    for j = 1:numel(idxNe_rule)
        kNe = idxNe_rule(j);
        tmp = niS(:,:,kNe);
        tmp(~isfinite(tmp)) = 0;
        fq_native(:,:,j) = tmp;
    end

    nNe_native_safe = nNe_native;
    nNe_native_safe(~isfinite(nNe_native_safe) | nNe_native_safe <= 0) = epsv;

    for j = 1:numel(idxNe_rule)
        fq_native(:,:,j) = fq_native(:,:,j) ./ nNe_native_safe;
    end

    fq_native(~isfinite(fq_native)) = 0;
    fq_native = max(fq_native, 0);

    for j = 1:numel(idxNe_rule)
        tmp = fq_native(:,:,j);
        tmp = masked_gauss_smooth(tmp, mask_native, 1.0);
        tmp(~isfinite(tmp)) = 0;
        fq_native(:,:,j) = max(tmp, 0);
    end

    fq_sum_native = sum(fq_native, 3);
    good_native = fq_sum_native > 0;
    for j = 1:numel(idxNe_rule)
        tmp = fq_native(:,:,j);
        tmp(good_native) = tmp(good_native) ./ fq_sum_native(good_native);
        fq_native(:,:,j) = tmp;
    end

    % ---------------------------------------------------------------------
    % 6) Map fractions to XY
    %    Inside: preserve SOLPS fraction structure
    %    Outside: smoothly freeze last-edge fractions
    % ---------------------------------------------------------------------
    fq_xy = zeros(nZout, nRout, numel(idxNe_rule));

    for j = 1:numel(idxNe_rule)
        fqj_native = fq_native(:,:,j);

        fqj_interp = interp2(rS, zS, fqj_native, X, Y, 'linear', NaN);

        ok_fq = outer_anchor_mask_species & isfinite(fqj_native) & (fqj_native >= 0);
        idx_fq = find(ok_fq(:));

        if isempty(idx_fq)
            fqj_extrap = zeros(size(X));
        else
            coords_fq = coords_ok(idx_fq,:);
            vals_fq   = fqj_native(idx_fq);

            K_use = min(max(4, K_fill), size(coords_fq,1));
            p_use = p_idw;

            fq_cell = cell(nBatches,1);

            parfor b = 1:nBatches
                s = (b-1)*batch_size + 1;
                e = min(b*batch_size, num_points_all);
                subset = extrap_coords(s:e,:);

                [iN, dN] = knnsearch(coords_fq, subset, 'K', K_use);

                if K_use == 1
                    iN = iN(:);
                    fq_out = vals_fq(iN);
                else
                    w = 1 ./ max(dN, 1e-12).^p_use;
                    w = w ./ sum(w, 2);
                    valsK = vals_fq(iN);
                    fq_out = sum(w .* valsK, 2);
                end

                fq_cell{b} = fq_out;
            end

            fqj_extrap = reshape(cell2mat(fq_cell), size(X));
        end

        fqj = fqj_interp;
        mask_fq_out = isnan(fqj_interp);
        fqj(mask_fq_out) = fqj_extrap(mask_fq_out);

        fqj = local_imgaussfilt_fallback(fqj, 0.8);
        fqj(~isfinite(fqj)) = 0;
        fq_xy(:,:,j) = max(fqj, 0);
    end

    fq_sum_xy = sum(fq_xy, 3);
    good_xy = fq_sum_xy > 0;
    for j = 1:numel(idxNe_rule)
        tmp = fq_xy(:,:,j);
        tmp(good_xy) = tmp(good_xy) ./ fq_sum_xy(good_xy);
        fq_xy(:,:,j) = tmp;
    end

    % ---------------------------------------------------------------------
    % 7) Reconstruct Neon charge states
    % ---------------------------------------------------------------------
    for j = 1:numel(idxNe_rule)
        kNe = idxNe_rule(j);
        val_ni_mspecies(:,:,kNe) = fq_xy(:,:,j) .* val_nNe_tot;
        val_ni_mspecies(:,:,kNe) = max(val_ni_mspecies(:,:,kNe), 0);
    end

    % exact renormalization
    nNe_sum = zeros(size(val_ne));
    for j = 1:numel(idxNe_rule)
        kNe = idxNe_rule(j);
        nNe_sum = nNe_sum + val_ni_mspecies(:,:,kNe);
    end

    renorm = ones(size(val_ne));
    m = nNe_sum > 0;
    renorm(m) = val_nNe_tot(m) ./ max(nNe_sum(m), epsv);

    for j = 1:numel(idxNe_rule)
        kNe = idxNe_rule(j);
        val_ni_mspecies(:,:,kNe) = val_ni_mspecies(:,:,kNe) .* renorm;
    end

    fprintf('✅ Neon rebuilt: no sharp handoff at last SOLPS grid; outside follows electrons.\n');
else
    fprintf('\n(idxNe_rule empty) -> skipping Neon rebuild.\n');
end
%% ------------------- psiN on (X,Y) and mask scalars -------------------
if apply_psi_mask
    try
        psiN = reshape(calc_psiN(g, X(:), Y(:), 0), size(X));
        mask = psiN < psiN_mask;
    catch MEpsi
        warning('calc_psiN failed (%s). Disabling psiN mask.', MEpsi.message);
        psiN = NaN(size(X));
        mask = false(size(X));
        apply_psi_mask = false;
    end
else
    psiN = NaN(size(X));
    mask = false(size(X));
end
fprintf('psiN mask: enabled=%d, threshold=%.2f, masked fraction=%.3f\n', ...
    apply_psi_mask, psiN_mask, nnz(mask)/numel(mask));

val_ne_masked = val_ne; val_ne_masked(mask) = 0;
val_Te_masked = val_Te; val_Te_masked(mask) = 0;
val_Ti_masked = val_Ti; val_Ti_masked(mask) = 0;

val_ni_mspecies_masked = val_ni_mspecies;
for k = 1:ns
    tmp = val_ni_mspecies_masked(:,:,k);
    tmp(mask) = 0;
    val_ni_mspecies_masked(:,:,k) = tmp;
end

% flows: keep as mapped (do not psi-mask)
val_uR_mspecies_masked = val_uR_mspecies;
val_uZ_mspecies_masked = val_uZ_mspecies;
val_uT_mspecies_masked = val_uT_mspecies;

% =========================
% Separatrix location at OMP (Z ≈ 0)
% =========================
[~, izMid] = min(abs(z_xy - 0));

psi_mid = psiN(izMid,:);   % psiN along OMP
Rline   = x_xy(:)';

% find where psiN crosses 1
R_sep_omp = NaN;

if any(isfinite(psi_mid))
    try
        R_sep_omp = interp1(psi_mid, Rline, 1.0, 'linear', 'extrap');
    catch
        % fallback: nearest point
        [~, idx_sep] = min(abs(psi_mid - 1));
        R_sep_omp = Rline(idx_sep);
    end
end

fprintf('OMP separatrix location: R_sep_omp = %.5f m\n', R_sep_omp);

%% =============================================================================
%% --- Load ITER antenna geometry (for OMP radial marker only) ---
geom_mat = 'final_ITER_data.mat';   % <-- set your actual file
z_centroid = [];
r_centroid = [];

if exist(geom_mat,'file')
    try
        load(geom_mat,'centroid');   % centroid is [Ntri x 3]
        if exist('centroid','var') && ~isempty(centroid) && size(centroid,2) >= 3
            z_centroid = centroid(:,3);
            r_centroid = hypot(centroid(:,1), centroid(:,2));
        else
            warning('centroid not found/invalid in %s', geom_mat);
        end
    catch MEg
        warning('Failed loading %s: %s', geom_mat, MEg.message);
    end
else
    warning('geom_mat not found: %s (OMP antenna markers will fall back).', geom_mat);
end

% --- OMP antenna radial marker (not antenna-center Z cuts anymore) ---
if ~isempty(r_centroid) && all(isfinite(r_centroid))
    R_ant_omp_min = min(r_centroid);
    R_ant_omp_max = max(r_centroid);
    R_ant_edge    = R_ant_omp_min;   % keep old naming if used elsewhere
else
    if hasWall
        R_ant_omp_min = max(g.lim(1,:)) - 0.05;
        R_ant_omp_max = max(g.lim(1,:));
        R_ant_edge    = R_ant_omp_min;
    else
        R_ant_omp_min = max(x_xy) - 0.05;
        R_ant_omp_max = max(x_xy);
        R_ant_edge    = R_ant_omp_min;
    end
end

dR_ant = 0.08;  % retained if needed elsewhere


%% =============================================================================
%% === FULL PLOT PACK: native grid + XY grid + flows + species + neon charge states
doPlotAll = true;

if doPlotAll
    fprintf('\n=== PLOT PACK: native grid + XY grid diagnostics (OMP-based 1D cuts) ===\n');

    % ------------------- Common helpers -------------------
    hasJet = true;
    if hasJet
        cmap = jet;
    else
        cmap = parula;
    end

    % Wall polygon for overlay (fully defensive)
    r_wall = []; z_wall = [];
    hasWall_local = false;

    if exist('g','var') && isstruct(g) && isfield(g,'lim') && ~isempty(g.lim) ...
            && size(g.lim,1) >= 2 && size(g.lim,2) >= 3
        hasWall_local = true;
        r_wall = g.lim(1,:);
        z_wall = g.lim(2,:);
    else
        if exist('lim','var') && isnumeric(lim) && size(lim,1) >= 2 && size(lim,2) >= 3
            hasWall_local = true;
            r_wall = lim(1,:);
            z_wall = lim(2,:);
        end
    end

    % Convenience: total Neon on native + XY
    nNe_total_native = zeros(size(neS));
    for j = 1:numel(idxNe_rule)
        nNe_total_native = nNe_total_native + niS(:,:,idxNe_rule(j));
    end

    nNe_total_XY = zeros(size(val_ne_masked));
    for j = 1:numel(idxNe_rule)
        nNe_total_XY = nNe_total_XY + val_ni_mspecies_masked(:,:,idxNe_rule(j));
    end

    % =========================================================================
    % 1) Native grid: painted fields
    % =========================================================================
    figure('Color','w','Position',[80 80 1600 900]);
    tiledlayout(2,3,'TileSpacing','compact','Padding','compact');

    nexttile;
    imagesc(rS, zS, neS); set(gca,'YDir','normal'); axis equal tight; colorbar;
    title('Native (painted) n_e'); xlabel('R [m]'); ylabel('Z [m]'); colormap(gca,cmap);
    hold on; if hasWall_local, plot(r_wall, z_wall, 'k','LineWidth',1.2); end; hold off;

    nexttile;
    imagesc(rS, zS, TeS); set(gca,'YDir','normal'); axis equal tight; colorbar;
    title('Native (painted) T_e [eV]'); xlabel('R [m]'); ylabel('Z [m]'); colormap(gca,cmap);
    hold on; if hasWall_local, plot(r_wall, z_wall, 'k','LineWidth',1.2); end; hold off;

    nexttile;
    imagesc(rS, zS, TiS); set(gca,'YDir','normal'); axis equal tight; colorbar;
    title('Native (painted) T_i [eV]'); xlabel('R [m]'); ylabel('Z [m]'); colormap(gca,cmap);
    hold on; if hasWall_local, plot(r_wall, z_wall, 'k','LineWidth',1.2); end; hold off;

    nexttile;
    imagesc(rS, zS, double(plasma_mask)); set(gca,'YDir','normal'); axis equal tight; colorbar;
    title(sprintf('plasma\\_mask (n_e > %.1e)', n_min)); xlabel('R [m]'); ylabel('Z [m]');
    colormap(gca,gray);
    hold on; if hasWall_local, plot(r_wall, z_wall, 'r','LineWidth',1.2); end; hold off;

    nexttile;
    imagesc(rS, zS, nNe_total_native); set(gca,'YDir','normal'); axis equal tight; colorbar;
    title('Native total Neon density'); xlabel('R [m]'); ylabel('Z [m]'); colormap(gca,cmap);
    set(gca,'ColorScale','log');
    hold on; if hasWall_local, plot(r_wall, z_wall, 'k','LineWidth',1.2); end; hold off;

    nexttile;
    imagesc(rS, zS, sqrt(BrS.^2+BtS.^2+BzS.^2)); set(gca,'YDir','normal'); axis equal tight; colorbar;
    title('|B| on native grid'); xlabel('R [m]'); ylabel('Z [m]'); colormap(gca,cmap);
    hold on; if hasWall_local, plot(r_wall, z_wall, 'k','LineWidth',1.2); end; hold off;

    sgtitle('Native-grid sanity: SOLPS-painted fields + masks + |B|','FontSize',14);

    % =========================================================================
    % 2) Native grid: D+, Ne, and a sample U||
    % =========================================================================
    figure('Color','w','Position',[90 90 1600 500]);
    tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

    nexttile;
    if iDplus <= ns
        imagesc(rS, zS, niS(:,:,iDplus)); set(gca,'YDir','normal'); axis equal tight; colorbar;
        title(sprintf('Native n_{D^+} (species %d)', iDplus)); xlabel('R'); ylabel('Z'); colormap(gca,cmap);
        set(gca,'ColorScale','log');
        hold on; if hasWall_local, plot(r_wall, z_wall, 'k','LineWidth',1.2); end; hold off;
    else
        axis off; text(0.1,0.5,'iDplus exceeds ns','FontSize',13);
    end

    nexttile;
    imagesc(rS, zS, nNe_total_native); set(gca,'YDir','normal'); axis equal tight; colorbar;
    title('Native total Neon'); xlabel('R'); ylabel('Z'); colormap(gca,cmap);
    set(gca,'ColorScale','log');
    hold on; if hasWall_local, plot(r_wall, z_wall, 'k','LineWidth',1.2); end; hold off;

    nexttile;
    kU = min(iDplus, ns);
    imagesc(rS, zS, UaS(:,:,kU)); set(gca,'YDir','normal'); axis equal tight; colorbar;
    title(sprintf('Native U_{||} (species %d)', kU)); xlabel('R'); ylabel('Z'); colormap(gca,cmap);
    hold on; if hasWall_local, plot(r_wall, z_wall, 'k','LineWidth',1.2); end; hold off;

    sgtitle('Native species + flow sanity','FontSize',14);

    % =========================================================================
    % 3) XY grid: scalars + psiN + mask
    % =========================================================================
    figure('Color','w','Position',[70 70 1700 900]);
    tiledlayout(2,3,'TileSpacing','compact','Padding','compact');

    nexttile;
    imagesc(x_xy, z_xy, val_ne_masked); set(gca,'YDir','normal'); axis equal tight; colorbar;
    title('XY n_e (masked)'); xlabel('R'); ylabel('Z'); colormap(gca,cmap);
    set(gca,'ColorScale','log');
    hold on; if hasWall_local, plot(r_wall, z_wall, 'k','LineWidth',1.2); end; hold off;

    nexttile;
    imagesc(x_xy, z_xy, val_Te_masked); set(gca,'YDir','normal'); axis equal tight; colorbar;
    title('XY T_e (masked)'); xlabel('R'); ylabel('Z'); colormap(gca,cmap);
    set(gca,'ColorScale','log');
    hold on; if hasWall_local, plot(r_wall, z_wall, 'k','LineWidth',1.2); end; hold off;

    nexttile;
    imagesc(x_xy, z_xy, val_Ti_masked); set(gca,'YDir','normal'); axis equal tight; colorbar;
    title('XY T_i (masked)'); xlabel('R'); ylabel('Z'); colormap(gca,cmap);
    set(gca,'ColorScale','log');
    hold on; if hasWall_local, plot(r_wall, z_wall, 'k','LineWidth',1.2); end; hold off;

    nexttile;
    imagesc(x_xy, z_xy, sqrt(Br_xy.^2 + Bt_xy.^2 + Bz_xy.^2)); set(gca,'YDir','normal');
    axis equal tight; colorbar; title('XY |B|'); xlabel('R'); ylabel('Z'); colormap(gca,cmap);
    hold on; if hasWall_local, plot(r_wall, z_wall, 'k','LineWidth',1.2); end; hold off;

    nexttile;
    if exist('psiN','var') && any(isfinite(psiN(:)))
        imagesc(x_xy, z_xy, psiN); set(gca,'YDir','normal'); axis equal tight; colorbar;
        title('XY \psi_N'); xlabel('R'); ylabel('Z'); colormap(gca,cmap);
        hold on; if hasWall_local, plot(r_wall, z_wall, 'k','LineWidth',1.2); end; hold off;
    else
        axis off; text(0.1,0.5,'psiN not available','FontSize',13);
    end

    nexttile;
    imagesc(x_xy, z_xy, double(mask)); set(gca,'YDir','normal'); axis equal tight; colorbar;
    title(sprintf('psi-mask (psiN < %.2f)', psiN_mask)); xlabel('R'); ylabel('Z');
    colormap(gca,gray);
    hold on; if hasWall_local, plot(r_wall, z_wall, 'r','LineWidth',1.2); end; hold off;

    sgtitle('XY scalars: n_e, T_e, T_i, |B|, \psi_N, mask','FontSize',14);

    % =========================================================================
    % 4) XY grid: D+, total Ne, Ne/n_e
    % =========================================================================
    ne_ratio_floor = 1e15;
    Ne_over_ne = nNe_total_XY ./ max(val_ne_masked, ne_ratio_floor);
    Ne_over_ne(val_ne_masked < ne_ratio_floor) = NaN;

    figure('Color','w','Position',[90 90 1650 520]);
    tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

    nexttile;
    if iDplus <= ns
        imagesc(x_xy, z_xy, val_ni_mspecies_masked(:,:,iDplus)); set(gca,'YDir','normal');
        axis equal tight; colorbar; title('XY n_{D^+}'); xlabel('R'); ylabel('Z'); colormap(gca,cmap);
        set(gca,'ColorScale','log');
        hold on; if hasWall_local, plot(r_wall, z_wall, 'k','LineWidth',1.2); end; hold off;
    else
        axis off; text(0.1,0.5,'iDplus exceeds ns','FontSize',13);
    end

    nexttile;
    imagesc(x_xy, z_xy, nNe_total_XY); set(gca,'YDir','normal');
    axis equal tight; colorbar; title('XY total Neon'); xlabel('R'); ylabel('Z'); colormap(gca,cmap);
    set(gca,'ColorScale','log');
    hold on; if hasWall_local, plot(r_wall, z_wall, 'k','LineWidth',1.2); end; hold off;

    nexttile;
    imagesc(x_xy, z_xy, Ne_over_ne); set(gca,'YDir','normal');
    axis equal tight; colorbar; title('XY (total Ne) / n_e'); xlabel('R'); ylabel('Z'); colormap(gca,cmap);
    set(gca,'ColorScale','log');
    hold on; if hasWall_local, plot(r_wall, z_wall, 'k','LineWidth',1.2); end; hold off;

    sgtitle('XY species: D+, total Ne, Ne/n_e','FontSize',14);

    % =========================================================================
    % 5) XY grid: flows (uR/uZ/uT) for D+
    % =========================================================================
    kFlow = min(iDplus, ns);
    uRk = val_uR_mspecies_masked(:,:,kFlow);
    uZk = val_uZ_mspecies_masked(:,:,kFlow);
    uTk = val_uT_mspecies_masked(:,:,kFlow);
    Umag = sqrt(uRk.^2 + uZk.^2 + uTk.^2);

    figure('Color','w','Position',[70 70 1700 520]);
    tiledlayout(1,4,'TileSpacing','compact','Padding','compact');

    nexttile;
    imagesc(x_xy, z_xy, uRk); set(gca,'YDir','normal'); axis equal tight; colorbar;
    title(sprintf('XY u_R (species %d)', kFlow)); xlabel('R'); ylabel('Z'); colormap(gca,cmap);

    nexttile;
    imagesc(x_xy, z_xy, uZk); set(gca,'YDir','normal'); axis equal tight; colorbar;
    title(sprintf('XY u_Z (species %d)', kFlow)); xlabel('R'); ylabel('Z'); colormap(gca,cmap);

    nexttile;
    imagesc(x_xy, z_xy, uTk); set(gca,'YDir','normal'); axis equal tight; colorbar;
    title(sprintf('XY u_T (species %d)', kFlow)); xlabel('R'); ylabel('Z'); colormap(gca,cmap);

    nexttile;
    imagesc(x_xy, z_xy, Umag); set(gca,'YDir','normal'); axis equal tight; colorbar;
    title(sprintf('XY |u| (species %d)', kFlow)); xlabel('R'); ylabel('Z'); colormap(gca,cmap);

    sgtitle('XY flows (projected from extrapolated U||)','FontSize',14);

 %% =========================================================================
%% 6) OMP 1D CUTS: full radial range + separatrix + shaded antenna region
%% =========================================================================
[~, izMid] = min(abs(z_xy - 0));
Rline = x_xy(:);

% --- OMP scalar/flow profiles ---
ne_mid = val_ne_masked(izMid,:).';
Te_mid = val_Te_masked(izMid,:).';
Ti_mid = val_Ti_masked(izMid,:).';

uR_mid = squeeze(val_uR_mspecies_masked(izMid,:,kFlow)).';
uZ_mid = squeeze(val_uZ_mspecies_masked(izMid,:,kFlow)).';
uT_mid = squeeze(val_uT_mspecies_masked(izMid,:,kFlow)).';

nD_mid = squeeze(val_ni_mspecies_masked(izMid,:,iDplus)).';

nNe_q_mid = zeros(numel(Rline), numel(idxNe_rule));
for j = 1:numel(idxNe_rule)
    kNe = idxNe_rule(j);
    nNe_q_mid(:,j) = squeeze(val_ni_mspecies_masked(izMid,:,kNe)).';
end
nNe_tot_mid = sum(nNe_q_mid, 2);
% OMP Neon charge-state fractions
f_q_mid = nNe_q_mid ./ max(nNe_tot_mid, epsv);
f_q_mid(~isfinite(f_q_mid)) = 0;

% --- separatrix location at OMP from psiN(R, Z~0) ---
R_sep_omp = NaN;
if exist('psiN','var') && any(isfinite(psiN(:)))
    psi_mid = psiN(izMid,:);
    goodPsi = isfinite(psi_mid) & isfinite(Rline.');
    if nnz(goodPsi) >= 2
        psi_good = psi_mid(goodPsi);
        R_good   = Rline(goodPsi);

        % sort by R for safety
        [R_good, isrt] = sort(R_good);
        psi_good = psi_good(isrt);

        % use nearest crossing to psiN = 1
        sgn = psi_good - 1.0;
        idx_cross = find(sgn(1:end-1).*sgn(2:end) <= 0, 1, 'last');

        if ~isempty(idx_cross)
            R_sep_omp = interp1(psi_good(idx_cross:idx_cross+1), ...
                                R_good(idx_cross:idx_cross+1), ...
                                1.0, 'linear', 'extrap');
        else
            [~, idx_sep] = min(abs(psi_good - 1.0));
            R_sep_omp = R_good(idx_sep);
        end
    end
end

fprintf('OMP separatrix location: R_sep_omp = %.5f m\n', R_sep_omp);

% --- antenna radial region from centroid-based R extent ---
if ~exist('R_ant_omp_min','var') || ~isfinite(R_ant_omp_min)
    R_ant_omp_min = NaN;
end
if ~exist('R_ant_omp_max','var') || ~isfinite(R_ant_omp_max)
    R_ant_omp_max = NaN;
end

has_ant_band = isfinite(R_ant_omp_min) && isfinite(R_ant_omp_max) && (R_ant_omp_max > R_ant_omp_min);

% helper for shaded antenna region
add_ant_patch = @(ax, y1, y2) patch(ax, ...
    [R_ant_omp_min R_ant_omp_max R_ant_omp_max R_ant_omp_min], ...
    [y1 y1 y2 y2], ...
    [1 0 1], 'FaceAlpha', 0.08, 'EdgeColor', 'none', 'HandleVisibility', 'off');

% =========================
% FIGURE 1: densities + temperatures (2 rows)
% =========================
figure('Color','w','Position',[80 80 1650 850]);
tl = tiledlayout(2,3,'TileSpacing','compact','Padding','compact');

% ---- Row 1, Col 1: all densities
ax1 = nexttile;
semilogy(Rline, max(ne_mid,1), 'k-', 'LineWidth', 1.8); hold on; grid on; box on;
semilogy(Rline, max(nD_mid,1), 'b--', 'LineWidth', 1.8);
semilogy(Rline, max(nNe_tot_mid,1), 'r-', 'LineWidth', 1.8);

yl = ylim;
if has_ant_band
    add_ant_patch(ax1, yl(1), yl(2));
end
xline(R_edge_omp,'k:','LineWidth',1.5);
if isfinite(R_sep_omp)
    xline(R_sep_omp,'g--','LineWidth',1.5);
end
if has_ant_band
    xline(R_ant_omp_min,'m--','LineWidth',1.2);
    xline(R_ant_omp_max,'m--','LineWidth',1.2);
end
title('OMP densities across full radial range');
xlabel('R [m]');
ylabel('Density [m^{-3}]');
legend_entries = {'n_e','n_{D^+}','n_{Ne,total}','R_{edge,OMP}'};
if isfinite(R_sep_omp)
    legend_entries{end+1} = 'R_{sep}';
end
if has_ant_band
    legend_entries{end+1} = 'R_{ant,min}';
    legend_entries{end+1} = 'R_{ant,max}';
end
legend(legend_entries,'Location','best');

% ---- Row 1, Col 2: ne only
ax2 = nexttile;
semilogy(Rline, max(ne_mid,1), 'k-', 'LineWidth', 1.8); hold on; grid on; box on;
yl = ylim;
if has_ant_band
    add_ant_patch(ax2, yl(1), yl(2));
end
xline(R_edge_omp,'k:','LineWidth',1.5);
if isfinite(R_sep_omp)
    xline(R_sep_omp,'g--','LineWidth',1.5);
end
if has_ant_band
    xline(R_ant_omp_min,'m--','LineWidth',1.2);
    xline(R_ant_omp_max,'m--','LineWidth',1.2);
end
title('OMP n_e');
xlabel('R [m]');
ylabel('n_e [m^{-3}]');

% ---- Row 1, Col 3: total Neon only
ax3 = nexttile;
semilogy(Rline, max(nNe_tot_mid,1), 'r-', 'LineWidth', 1.8); hold on; grid on; box on;
yl = ylim;
if has_ant_band
    add_ant_patch(ax3, yl(1), yl(2));
end
xline(R_edge_omp,'k:','LineWidth',1.5);
if isfinite(R_sep_omp)
    xline(R_sep_omp,'g--','LineWidth',1.5);
end
if has_ant_band
    xline(R_ant_omp_min,'m--','LineWidth',1.2);
    xline(R_ant_omp_max,'m--','LineWidth',1.2);
end
title('OMP total Neon');
xlabel('R [m]');
ylabel('n_{Ne,total} [m^{-3}]');

% ---- Row 2, Col 1: Te
ax4 = nexttile;
plot(Rline, Te_mid, 'r-', 'LineWidth', 1.8); hold on; grid on; box on;
yl = ylim;
if has_ant_band
    add_ant_patch(ax4, yl(1), yl(2));
end
xline(R_edge_omp,'k:','LineWidth',1.5);
if isfinite(R_sep_omp)
    xline(R_sep_omp,'g--','LineWidth',1.5);
end
if has_ant_band
    xline(R_ant_omp_min,'m--','LineWidth',1.2);
    xline(R_ant_omp_max,'m--','LineWidth',1.2);
end
title('OMP T_e');
xlabel('R [m]');
ylabel('T_e [eV]');

% ---- Row 2, Col 2: Ti
ax5 = nexttile;
plot(Rline, Ti_mid, 'b-', 'LineWidth', 1.8); hold on; grid on; box on;
yl = ylim;
if has_ant_band
    add_ant_patch(ax5, yl(1), yl(2));
end
xline(R_edge_omp,'k:','LineWidth',1.5);
if isfinite(R_sep_omp)
    xline(R_sep_omp,'g--','LineWidth',1.5);
end
if has_ant_band
    xline(R_ant_omp_min,'m--','LineWidth',1.2);
    xline(R_ant_omp_max,'m--','LineWidth',1.2);
end
title('OMP T_i');
xlabel('R [m]');
ylabel('T_i [eV]');

% ---- Row 2, Col 3: Te and Ti together
ax6 = nexttile;
plot(Rline, Te_mid, 'r-', 'LineWidth', 1.6); hold on; grid on; box on;
plot(Rline, Ti_mid, 'b--', 'LineWidth', 1.6);
yl = ylim;
if has_ant_band
    add_ant_patch(ax6, yl(1), yl(2));
end
xline(R_edge_omp,'k:','LineWidth',1.5);
if isfinite(R_sep_omp)
    xline(R_sep_omp,'g--','LineWidth',1.5);
end
if has_ant_band
    xline(R_ant_omp_min,'m--','LineWidth',1.2);
    xline(R_ant_omp_max,'m--','LineWidth',1.2);
end
title('OMP temperatures');
xlabel('R [m]');
ylabel('T [eV]');
legT = {'T_e','T_i','R_{edge,OMP}'};
if isfinite(R_sep_omp)
    legT{end+1} = 'R_{sep}';
end
if has_ant_band
    legT{end+1} = 'R_{ant,min}';
    legT{end+1} = 'R_{ant,max}';
end
legend(legT,'Location','best');

sgtitle(tl,'OMP densities and temperatures across full radial range','FontSize',16);

% =========================
% FIGURE 2: flows (1 row)
% =========================
figure('Color','w','Position',[100 100 1650 380]);
tl2 = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

ax7 = nexttile;
plot(Rline, uR_mid, 'LineWidth', 1.6); hold on; grid on; box on;
yl = ylim;
if has_ant_band
    add_ant_patch(ax7, yl(1), yl(2));
end
xline(R_edge_omp,'k:','LineWidth',1.5);
if isfinite(R_sep_omp)
    xline(R_sep_omp,'g--','LineWidth',1.5);
end
if has_ant_band
    xline(R_ant_omp_min,'m--','LineWidth',1.2);
    xline(R_ant_omp_max,'m--','LineWidth',1.2);
end
title(sprintf('OMP u_R (species %d)', kFlow));
xlabel('R [m]');
ylabel('u_R [m/s]');

ax8 = nexttile;
plot(Rline, uZ_mid, 'LineWidth', 1.6); hold on; grid on; box on;
yl = ylim;
if has_ant_band
    add_ant_patch(ax8, yl(1), yl(2));
end
xline(R_edge_omp,'k:','LineWidth',1.5);
if isfinite(R_sep_omp)
    xline(R_sep_omp,'g--','LineWidth',1.5);
end
if has_ant_band
    xline(R_ant_omp_min,'m--','LineWidth',1.2);
    xline(R_ant_omp_max,'m--','LineWidth',1.2);
end
title(sprintf('OMP u_Z (species %d)', kFlow));
xlabel('R [m]');
ylabel('u_Z [m/s]');

ax9 = nexttile;
plot(Rline, uT_mid, 'LineWidth', 1.6); hold on; grid on; box on;
yl = ylim;
if has_ant_band
    add_ant_patch(ax9, yl(1), yl(2));
end
xline(R_edge_omp,'k:','LineWidth',1.5);
if isfinite(R_sep_omp)
    xline(R_sep_omp,'g--','LineWidth',1.5);
end
if has_ant_band
    xline(R_ant_omp_min,'m--','LineWidth',1.2);
    xline(R_ant_omp_max,'m--','LineWidth',1.2);
end
title(sprintf('OMP u_T (species %d)', kFlow));
xlabel('R [m]');
ylabel('u_T [m/s]');

sgtitle(tl2,'OMP flows across full radial range','FontSize',16);



% =========================
% FIGURE 3: Mean Neon abundance in outer OMP band
% =========================
maskR_outer = (Rline >= (R_edge_omp - 0.08));
fq_mean = mean(f_q_mid(maskR_outer,:), 1, 'omitnan');
fq_mean = fq_mean ./ max(sum(fq_mean), epsv);

figure('Color','w','Position',[140 140 850 450]);
bar(0:(numel(idxNe_rule)-1), 100*fq_mean, 'FaceColor','flat');
grid on; box on;
xlabel('Ne charge state q');
ylabel('Mean abundance [%]');
title(sprintf('Ne charge-state abundance (avg over R >= R_{edge}-0.08, Z≈0)'));
set(gca,'XTick',0:(numel(idxNe_rule)-1));

fprintf('✅ OMP 1D plot section complete.\n');

   %% === NEW FIG: 2x2 OMP summary (row1: n_e & n_{Ne,total}, row2: T_e & T_i) ===
% Assumes Rline, ne_mid, nNe_tot_mid, Te_mid, Ti_mid, R_edge_omp are in workspace.
% Optional: R_sep_omp, R_ant_omp_min, R_ant_omp_max may be present.

% defensive checks / fallbacks
if ~exist('Rline','var') || ~exist('ne_mid','var')
    error('Required variables Rline and ne_mid not found in workspace.');
end
if ~exist('nNe_tot_mid','var')
    nNe_tot_mid = zeros(size(ne_mid));
end
if ~exist('Te_mid','var'), Te_mid = zeros(size(ne_mid)); end
if ~exist('Ti_mid','var'), Ti_mid = zeros(size(ne_mid)); end
if ~exist('R_edge_omp','var'), R_edge_omp = NaN; end
if ~exist('R_sep_omp','var'), R_sep_omp = NaN; end
if ~exist('R_ant_omp_min','var'), R_ant_omp_min = NaN; end
if ~exist('R_ant_omp_max','var'), R_ant_omp_max = NaN; end

has_ant_band = isfinite(R_ant_omp_min) && isfinite(R_ant_omp_max) && (R_ant_omp_max > R_ant_omp_min);

% small helper to draw antenna patch (works per axis)
draw_ant_patch = @(ax, ylo, yhi) patch(ax, ...
    [R_ant_omp_min R_ant_omp_max R_ant_omp_max R_ant_omp_min], ...
    [ylo ylo yhi yhi], ...
    [1 0 1], 'FaceAlpha', 0.08, 'EdgeColor', 'none', 'HandleVisibility','off');

%% === OMP summary (2 panels): densities + temperatures (LOG SCALE) ===
figure('Color','w','Position',[120 120 1000 700]);
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

has_ant_band = isfinite(R_ant_omp_min) && isfinite(R_ant_omp_max);

% =========================
% TOP: densities (LOG)
% =========================
ax1 = nexttile;
semilogy(Rline, max(ne_mid,1),      'k-', 'LineWidth', 1.8); hold on;
semilogy(Rline, max(nNe_tot_mid,1), 'r-', 'LineWidth', 1.8);

grid on; box on;
xlabel('R [m]');
ylabel('Density [m^{-3}]');
title('OMP densities (Z \approx 0)');

xlim([8.035 8.59]);
ylim([1e12 5e20]);

yl = ylim(ax1);

% ---- antenna shading ----
if has_ant_band
    patch([R_ant_omp_min R_ant_omp_max R_ant_omp_max R_ant_omp_min], ...
          [yl(1) yl(1) yl(2) yl(2)], ...
          [1 0 1], 'FaceAlpha',0.08,'EdgeColor','none');

    % label antenna region
    text(mean([R_ant_omp_min R_ant_omp_max]), yl(2)/3, ...
        'Antenna region', 'Color','m', 'HorizontalAlignment','center');
end

% ---- reference lines + labels ----
if isfinite(R_edge_omp)
    xline(R_edge_omp,'k:','LineWidth',1.6);
    text(R_edge_omp, yl(2)/2, ' OMP', 'Color','k','FontWeight','bold');
end

if isfinite(R_sep_omp)
    xline(R_sep_omp,'g--','LineWidth',1.6);
    text(R_sep_omp, yl(2)/4, ' Separatrix', 'Color','g','FontWeight','bold');
end

if has_ant_band
    xline(R_ant_omp_min,'m--','LineWidth',1.2);
    xline(R_ant_omp_max,'m--','LineWidth',1.2);
end

legend('n_e','n_{Ne,total}','Location','best');


% =========================
% BOTTOM: temperatures (LOG)
% =========================
ax2 = nexttile;

Te_plot = max(Te_mid, 1e-3);
Ti_plot = max(Ti_mid, 1e-3);

semilogy(Rline, Te_plot, 'r-',  'LineWidth', 1.6); hold on;
semilogy(Rline, Ti_plot, 'b--', 'LineWidth', 1.6);

grid on; box on;
xlabel('R [m]');
ylabel('T [eV]');
title('OMP temperatures (Z \approx 0)');

xlim([8.035 8.59]);
ylim([1 4e3]);

yl = ylim(ax2);

% antenna shading
if has_ant_band
    patch([R_ant_omp_min R_ant_omp_max R_ant_omp_max R_ant_omp_min], ...
          [yl(1) yl(1) yl(2) yl(2)], ...
          [1 0 1], 'FaceAlpha',0.08,'EdgeColor','none');

    text(mean([R_ant_omp_min R_ant_omp_max]), yl(2)/3, ...
        'Antenna region', 'Color','m','HorizontalAlignment','center');
end

% reference lines + labels
if isfinite(R_edge_omp)
    xline(R_edge_omp,'k:','LineWidth',1.6);
    text(R_edge_omp, yl(2)/2, ' OMP', 'Color','k','FontWeight','bold');
end

if isfinite(R_sep_omp)
    xline(R_sep_omp,'g--','LineWidth',1.6);
    text(R_sep_omp, yl(2)/4, ' Separatrix', 'Color','g','FontWeight','bold');
end

if has_ant_band
    xline(R_ant_omp_min,'m--','LineWidth',1.2);
    xline(R_ant_omp_max,'m--','LineWidth',1.2);
end

legend('T_e','T_i','Location','best');

sgtitle('OMP profiles (log scale) with antenna region, OMP, and separatrix','FontSize',14);


%% === OMP Neon charge-state fractional abundance ===
figure('Color','w','Position',[120 120 1100 520]);
hold on; box on; grid on;

cmap = jet(numel(idxNe_rule));

NeNames = arrayfun(@(q) sprintf('Ne^{%d+}', q), 0:numel(idxNe_rule)-1, 'UniformOutput', false);

for j = 1:numel(idxNe_rule)
    plot(Rline, f_q_mid(:,j), 'LineWidth', 1.6, ...
        'Color', cmap(j,:), 'DisplayName', NeNames{j});
end

xlim([8.02 8.6]);
yl = ylim;

% antenna shading + label
if has_ant_band
    patch([R_ant_omp_min R_ant_omp_max R_ant_omp_max R_ant_omp_min], ...
          [yl(1) yl(1) yl(2) yl(2)], ...
          [1 0 1], 'FaceAlpha',0.08,'EdgeColor','none');

    text(mean([R_ant_omp_min R_ant_omp_max]), 0.8*yl(2), ...
        'Antenna region', 'Color','m','HorizontalAlignment','center');
end

% reference lines + labels
if isfinite(R_edge_omp)
    xline(R_edge_omp,'k:','LineWidth',1.8);
    text(R_edge_omp, 0.9*yl(2), 'OMP', 'Color','k','FontWeight','bold');
end

if isfinite(R_sep_omp)
    xline(R_sep_omp,'g--','LineWidth',1.6);
    text(R_sep_omp, 0.7*yl(2), 'Separatrix', 'Color','g','FontWeight','bold');
end

if has_ant_band
    xline(R_ant_omp_min,'m--','LineWidth',1.2);
    xline(R_ant_omp_max,'m--','LineWidth',1.2);
end

xlabel('R [m]');
ylabel('n_{Ne,q} / \Sigma_q n_{Ne,q}');
title('OMP Neon charge-state fractional abundance');

legend('Location','bestoutside');
legend boxoff;

end
%% =============================================================================
%% ------------------- NetCDF Export (switchable) -------------------
xw = double(x_xy(:));   % R
zw = double(z_xy(:));   % Z
nX = numel(xw);
nZ = numel(zw);

switch outputMode

    case 'singleFluid'
        % ------------------- NetCDF Export (single-fluid only) -------------------
        x  = X(1,:);          % R
        z  = Y(:,1);          % Z
        nX = length(x);
        nZ = length(z);

        outnc = 'profiles_iter_single_GITRstyle_XY.nc';
        ncid = netcdf.create(outnc, ...
            bitor(netcdf.getConstant('NETCDF4'), netcdf.getConstant('CLOBBER')));

        % ------------------- Dimensions -------------------
        dimR = netcdf.defDim(ncid,'nX', nX);
        dimZ = netcdf.defDim(ncid,'nZ', nZ);

        % ------------------- Coordinate variables -------------------
        varidX = netcdf.defVar(ncid,'x','double',dimR);
        varidZ = netcdf.defVar(ncid,'z','double',dimZ);

        % ------------------- 2D plasma profile variables -------------------
        varidNe     = netcdf.defVar(ncid,'ne',     'double',[dimR dimZ]);
        varidNi     = netcdf.defVar(ncid,'ni',     'double',[dimR dimZ]);
        varidTe     = netcdf.defVar(ncid,'te',     'double',[dimR dimZ]);
        varidTi     = netcdf.defVar(ncid,'ti',     'double',[dimR dimZ]);
        varidGradTi = netcdf.defVar(ncid,'gradTi', 'double',[dimR dimZ]);
        varidGradTe = netcdf.defVar(ncid,'gradTe', 'double',[dimR dimZ]);
        varidVr     = netcdf.defVar(ncid,'vr',     'double',[dimR dimZ]);
        varidVt     = netcdf.defVar(ncid,'vt',     'double',[dimR dimZ]);
        varidVz     = netcdf.defVar(ncid,'vz',     'double',[dimR dimZ]);
        varidPsi    = netcdf.defVar(ncid,'psiN',   'double',[dimR dimZ]);

        % ------------------- Optional B-field variables -------------------
        varidBr = netcdf.defVar(ncid,'br','double',[dimR dimZ]);
        varidBt = netcdf.defVar(ncid,'bt','double',[dimR dimZ]);
        varidBz = netcdf.defVar(ncid,'bz','double',[dimR dimZ]);

        netcdf.endDef(ncid);

        % ------------------- Write coordinates -------------------
        netcdf.putVar(ncid, varidX, x);
        netcdf.putVar(ncid, varidZ, z);

        % ------------------- Single-fluid fields -------------------
        % If you want strict single-fluid closure:
        %   ni = ne
        %   Ti = Te
        % Otherwise replace with your preferred single-ion fields.

        netcdf.putVar(ncid, varidNe,     permute(double(val_ne_masked), [2 1]));
        netcdf.putVar(ncid, varidNi,     permute(double(val_ne_masked), [2 1]));
        netcdf.putVar(ncid, varidTe,     permute(double(val_Te_masked), [2 1]));
        netcdf.putVar(ncid, varidTi,     permute(double(val_Ti_masked), [2 1]));

        % gradTi / gradTe:
        % If val_gradTi_masked and val_gradTe_masked exist, use them.
        % Otherwise fall back safely.
        if exist('val_gradTi_masked','var')
            netcdf.putVar(ncid, varidGradTi, permute(double(val_gradTi_masked), [2 1]));
        else
            netcdf.putVar(ncid, varidGradTi, zeros(nX,nZ,'double'));
        end

        if exist('val_gradTe_masked','var')
            netcdf.putVar(ncid, varidGradTe, permute(double(val_gradTe_masked), [2 1]));
        elseif exist('val_gradTi_masked','var')
            netcdf.putVar(ncid, varidGradTe, permute(double(val_gradTi_masked), [2 1]));
        else
            netcdf.putVar(ncid, varidGradTe, zeros(nX,nZ,'double'));
        end

        % flows:
        % Prefer single-fluid velocity fields if they exist.
        if exist('val_vr_masked','var') && exist('val_vt_masked','var') && exist('val_vz_masked','var')
            netcdf.putVar(ncid, varidVr, permute(double(val_vr_masked), [2 1]));
            netcdf.putVar(ncid, varidVt, permute(double(val_vt_masked), [2 1]));
            netcdf.putVar(ncid, varidVz, permute(double(val_vz_masked), [2 1]));
        else
            % fallback to D+ from multi-species if single-fluid vars are absent
            kSingle = iDplus;
            netcdf.putVar(ncid, varidVr, permute(double(val_uR_mspecies_masked(:,:,kSingle)), [2 1]));
            netcdf.putVar(ncid, varidVt, permute(double(val_uT_mspecies_masked(:,:,kSingle)), [2 1]));
            netcdf.putVar(ncid, varidVz, permute(double(val_uZ_mspecies_masked(:,:,kSingle)), [2 1]));
        end

        % psiN
        if exist('psiN','var') && any(isfinite(psiN(:)))
            netcdf.putVar(ncid, varidPsi, permute(double(psiN), [2 1]));
        else
            netcdf.putVar(ncid, varidPsi, zeros(nX,nZ,'double'));
        end

        % B-fields
        netcdf.putVar(ncid, varidBr, permute(double(Br_xy), [2 1]));
        netcdf.putVar(ncid, varidBt, permute(double(Bt_xy), [2 1]));
        netcdf.putVar(ncid, varidBz, permute(double(Bz_xy), [2 1]));

        netcdf.close(ncid);
        fprintf('✅ Wrote %s\n', outnc);


    case 'multiFluid'
        outnc = 'profiles_iter_multi_GITRstyle_XY.nc';

        if exist(outnc,'file'), delete(outnc); end
        ncid = netcdf.create(outnc, bitor(netcdf.getConstant('NETCDF4'), ...
                                          netcdf.getConstant('CLOBBER')));

        % Dims
        dimX = netcdf.defDim(ncid,'nX', nX);
        dimZ = netcdf.defDim(ncid,'nZ', nZ);
        dimS = netcdf.defDim(ncid,'species', ns);

        % Coords / meta
        vid_x   = netcdf.defVar(ncid,'x','double',dimX);
        vid_z   = netcdf.defVar(ncid,'z','double',dimZ);
        vid_Z   = netcdf.defVar(ncid,'atomic_number','double',dimS);
        vid_q   = netcdf.defVar(ncid,'charge_number','double',dimS);
        vid_psi = netcdf.defVar(ncid,'psiN','double',[dimX dimZ]);

        % Single-fluid 2D
        vid_ne  = netcdf.defVar(ncid,'ne','double',[dimX dimZ]);
        vid_te  = netcdf.defVar(ncid,'te','double',[dimX dimZ]);
        vid_ti  = netcdf.defVar(ncid,'ti','double',[dimX dimZ]);

        vid_br  = netcdf.defVar(ncid,'br','double',[dimX dimZ]);
        vid_bt  = netcdf.defVar(ncid,'bt','double',[dimX dimZ]);
        vid_bz  = netcdf.defVar(ncid,'bz','double',[dimX dimZ]);

        % Multi-species 3D
        vid_niA = netcdf.defVar(ncid,'ni_all','double',[dimX dimZ dimS]);
        vid_uRA = netcdf.defVar(ncid,'uR_all','double',[dimX dimZ dimS]);
        vid_uZA = netcdf.defVar(ncid,'uZ_all','double',[dimX dimZ dimS]);
        vid_uTA = netcdf.defVar(ncid,'uT_all','double',[dimX dimZ dimS]);

        netcdf.endDef(ncid);

        % Write coords/meta
        netcdf.putVar(ncid, vid_x, xw);
        netcdf.putVar(ncid, vid_z, zw);
        netcdf.putVar(ncid, vid_Z, double(zn(:)));
        netcdf.putVar(ncid, vid_q, double(zn(:)));   % replace later if explicit charge numbers exist

        if exist('psiN','var') && any(isfinite(psiN(:)))
            netcdf.putVar(ncid, vid_psi, permute(double(psiN), [2 1]));
        else
            netcdf.putVar(ncid, vid_psi, zeros(nX, nZ, 'double'));
        end

        % Write single-fluid fields
        netcdf.putVar(ncid, vid_ne, permute(double(val_ne_masked), [2 1]));
        netcdf.putVar(ncid, vid_te, permute(double(val_Te_masked), [2 1]));
        netcdf.putVar(ncid, vid_ti, permute(double(val_Ti_masked), [2 1]));

        netcdf.putVar(ncid, vid_br, permute(double(Br_xy), [2 1]));
        netcdf.putVar(ncid, vid_bt, permute(double(Bt_xy), [2 1]));
        netcdf.putVar(ncid, vid_bz, permute(double(Bz_xy), [2 1]));

        % Write multi-species fields
        netcdf.putVar(ncid, vid_niA, permute(double(val_ni_mspecies_masked), [2 1 3]));
        netcdf.putVar(ncid, vid_uRA, permute(double(val_uR_mspecies_masked), [2 1 3]));
        netcdf.putVar(ncid, vid_uZA, permute(double(val_uZ_mspecies_masked), [2 1 3]));
        netcdf.putVar(ncid, vid_uTA, permute(double(val_uT_mspecies_masked), [2 1 3]));

        netcdf.close(ncid);
        fprintf('✅ Wrote %s\n', outnc);


    otherwise
        error('Unknown outputMode = %s. Use ''singleFluid'' or ''multiSpecies''.', outputMode);
end

%% =============================================================================
%% === OPTIONAL: Overlay/compare flows from an existing GITR NetCDF on this XY grid ===
doPlotGITRflows = true;
gitr_nc = 'profiles_iter_multi_GITRstyle_XY.nc';   % <-- set this (can be another file)

if doPlotGITRflows && exist(gitr_nc,'file')
    fprintf('\n=== Loading GITR flows for visualization: %s ===\n', gitr_nc);

    infoG = ncinfo(gitr_nc);
    vnames = string({infoG.Variables.Name});
    has = @(nm) any(vnames == string(nm));
    rd  = @(nm) ncread(gitr_nc, nm);

    xG = rd('x');   % [nX]
    zG = rd('z');   % [nZ]
    nXg = numel(xG);
    nZg = numel(zG);

    [XGm, ZGm] = meshgrid(xG, zG);   % [nZ x nX]

    fix2D_ZxX = @(A) local_force_ZxX(A, nXg, nZg);

    % Convention A (single-field)
    vrG = []; vzG = []; vtG = [];
    if has('vr'), vrG = rd('vr'); end
    if has('vz'), vzG = rd('vz'); end
    if has('vt'), vtG = rd('vt'); end

    % Convention B (multi-species)
    uR_G = []; uZ_G = []; uT_G = [];
    if isempty(vrG) && has('uR_all'), uR_G = rd('uR_all'); end
    if isempty(vzG) && has('uZ_all'), uZ_G = rd('uZ_all'); end
    if isempty(vtG) && has('uT_all'), uT_G = rd('uT_all'); end

    if ~isempty(vrG)
        vrG = fix2D_ZxX(vrG);
        vzG = fix2D_ZxX(vzG);
        vtG = fix2D_ZxX(vtG);

        vR_xy = interp2(XGm, ZGm, vrG, X, Y, 'linear', NaN);
        vZ_xy = interp2(XGm, ZGm, vzG, X, Y, 'linear', NaN);
        vT_xy = interp2(XGm, ZGm, vtG, X, Y, 'linear', NaN);
        flowTag = 'GITR (vr/vz/vt)';

    elseif ~isempty(uR_G)
        kPlot = iDplus;   % D+ by default

        uRk = local_force_ZxX_slice(uR_G, nXg, nZg, kPlot);
        uZk = local_force_ZxX_slice(uZ_G, nXg, nZg, kPlot);
        uTk = local_force_ZxX_slice(uT_G, nXg, nZg, kPlot);

        vR_xy = interp2(XGm, ZGm, uRk, X, Y, 'linear', NaN);
        vZ_xy = interp2(XGm, ZGm, uZk, X, Y, 'linear', NaN);
        vT_xy = interp2(XGm, ZGm, uTk, X, Y, 'linear', NaN);
        flowTag = sprintf('GITR (uR/uZ/uT), species %d', kPlot);

    else
        warning('No flow vars found in %s. Expected vr/vz/vt or uR_all/uZ_all/uT_all.', gitr_nc);
        vR_xy = []; vZ_xy = []; vT_xy = []; flowTag = '';
    end

    if ~isempty(vR_xy)
        [~, izMid] = min(abs(z_xy - 0));

        figure('Color','w','Position',[120 120 1500 420]);
        tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

        nexttile; imagesc(x_xy, z_xy, vR_xy); set(gca,'YDir','normal'); axis equal tight; colorbar;
        title(['v_R on XY: ' flowTag]); xlabel('R'); ylabel('Z');

        nexttile; imagesc(x_xy, z_xy, vZ_xy); set(gca,'YDir','normal'); axis equal tight; colorbar;
        title(['v_Z on XY: ' flowTag]); xlabel('R'); ylabel('Z');

        nexttile; imagesc(x_xy, z_xy, vT_xy); set(gca,'YDir','normal'); axis equal tight; colorbar;
        title(['v_T on XY: ' flowTag]); xlabel('R'); ylabel('Z');

        % Midplane comparison (same species)
        kCmp = iDplus;
        uR_my = squeeze(val_uR_mspecies_masked(izMid,:,kCmp));
        uZ_my = squeeze(val_uZ_mspecies_masked(izMid,:,kCmp));
        uT_my = squeeze(val_uT_mspecies_masked(izMid,:,kCmp));

        figure('Color','w','Position',[140 140 1200 450]);
        plot(x_xy, uR_my, 'LineWidth', 1.8); hold on; grid on; box on;
        plot(x_xy, uZ_my, 'LineWidth', 1.8);
        plot(x_xy, uT_my, 'LineWidth', 1.8);
        plot(x_xy, vR_xy(izMid,:), '--', 'LineWidth', 1.6);
        plot(x_xy, vZ_xy(izMid,:), '--', 'LineWidth', 1.6);
        plot(x_xy, vT_xy(izMid,:), '--', 'LineWidth', 1.6);
        xlabel('R [m]'); ylabel('Flow [m/s]');
        title(sprintf('Midplane flows: your XY (solid) vs %s (dashed)', flowTag));
        legend('uR yours','uZ yours','uT yours','uR GITR','uZ GITR','uT GITR','Location','best');
    end

else
    if doPlotGITRflows
        warning('GITR file not found: %s (skipping GITR flow visualization).', gitr_nc);
    end
end

%% =============================================================================
%% ------------------- Local helpers (orientation-safe) -------------------
function V = local_force_ZxX(A, nX, nZ)
% Return V as [nZ x nX] (for meshgrid(x,z)->[nZ x nX])
    sz = size(A);

    if isequal(sz, [nZ nX])
        V = A; return
    elseif isequal(sz, [nX nZ])
        V = A.'; return
    elseif numel(A)==nX*nZ
        V = reshape(A, [nX nZ]).'; return
    else
        error('local_force_ZxX: cannot coerce size %s to [nZ x nX]=[%d %d].', ...
              mat2str(sz), nZ, nX);
    end
end

function Vk = local_force_ZxX_slice(A3, nX, nZ, k)
% Extract species slice and return as [nZ x nX]
    sz = size(A3);
    if numel(sz) ~= 3
        error('local_force_ZxX_slice: expected 3D array, got %s', mat2str(sz));
    end

    % Common NetCDF layouts:
    %   [nX x nZ x ns]  or  [nZ x nX x ns]
    if sz(1)==nX && sz(2)==nZ
        Vk = A3(:,:,k).';
        return
    elseif sz(1)==nZ && sz(2)==nX
        Vk = A3(:,:,k);
        return
    end

    % fallback: try permutes
    perms_all = perms(1:3);
    for i = 1:size(perms_all,1)
        P = perms_all(i,:);
        B = permute(A3, P);
        s2 = size(B);
        if s2(1)==nX && s2(2)==nZ
            Vk = B(:,:,k).';
            return
        elseif s2(1)==nZ && s2(2)==nX
            Vk = B(:,:,k);
            return
        end
    end

    error('local_force_ZxX_slice: cannot coerce size %s with nX=%d nZ=%d.', mat2str(sz), nX, nZ);
end

%% =============================================================================
%% ------------------- Helper functions -------------------
function Q = paintCells(vals, r, z, X, Z)
% Fill each SOLPS quad cell onto (X,Z) grid with constant value.
% r,z: [4 x Nc], vals: [Nc x 1], X,Z: [nZ x nR]
    [nZ, nR] = size(X);
    Q = NaN(nZ, nR);
    Nc = size(r,2);
    for i = 1:Nc
        ri = r(:,i);
        zi = z(:,i);
        if any(isnan(ri)) || any(isnan(zi)), continue; end
        [in,on] = inpolygon(X, Z, ri, zi);
        m = (in | on);
        if any(m(:)), Q(m) = vals(i); end
    end
end

function Aout = masked_gauss_smooth(Ain, mask, sigma_pix)
% Gaussian smoothing using conv2, normalized by convolved mask.
    Ain2 = Ain;
    Ain2(~mask) = 0;

    rad = max(1, ceil(3*sigma_pix));
    x = (-rad:rad);
    g = exp(-(x.^2)/(2*sigma_pix^2));
    g = g / sum(g);
    K = g' * g;

    num = conv2(Ain2, K, 'same');
    den = conv2(double(mask), K, 'same');

    Aout = num ./ max(den, 1e-12);
    Aout(~mask) = NaN;
end
function L = local_decay_from_grad(R, f, Lmin, Lmax)
% L = -1/(d ln(f)/dR) from a short stencil; clamps output to [Lmin,Lmax].
    R = double(R(:));
    f = double(f(:));

    ok = isfinite(R) & isfinite(f) & (f > 0);
    if nnz(ok) < 3
        L = NaN; return;
    end
    R = R(ok); f = f(ok);

    dln_dR = gradient(log(f), R);
    s = dln_dR(end);

    if ~isfinite(s) || s >= -1e-6
        s = -1e-4; % fallback: weak/positive gradient
    end

    L = -1 / s;
    L = max(Lmin, min(Lmax, L));
end

function A = local_imgaussfilt_fallback(Ain, sigma)
    if sigma <= 0
        A = Ain;
        return;
    end

    if exist('imgaussfilt','file') == 2
        A = imgaussfilt(Ain, sigma);
        return;
    end

    rad = max(1, ceil(3*sigma));
    x = -rad:rad;
    g = exp(-(x.^2)/(2*sigma^2));
    g = g / sum(g);
    K = g' * g;
    A = conv2(double(Ain), K, 'same');
end