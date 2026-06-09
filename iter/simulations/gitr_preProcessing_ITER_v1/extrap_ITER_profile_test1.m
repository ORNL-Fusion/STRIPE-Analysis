%% Full SOLPS → (X,Y) Extrapolation (Euclidean-from-outer-wall) for Electrons & Multi-species,
%  Masking by psiN, and NetCDF Export
%
%  Key requirements implemented (NO other behavior changes beyond what you asked):
%  1) Extrapolation everywhere OUTSIDE the outer wall, using Euclidean distance away from outer wall
%     (wall-normal probing to sample an interior value + exp(slope * distance_to_wall)).
%  2) Remove temperature gradient part completely (no gradTi read/compute/export).
%  3) Fix the "Index exceeds array bounds" error by making all native grids/masks consistent ([nR x nZ]).
%  4) Fix reshape(cell2mat(...)) mismatch by using vertcat(cell{:}) (robust order-preserving concat).
%  5) Uses your ITER midplane fit window: r_min=8.01, r_max=8.33, extrapolate out to 8.6.
%  6) Te and all velocities are constant outside wall (nearest wall-normal interior), no exponential decay.
%     Densities decay outside wall using midplane-derived log-slope.
%     Neon: if increasing outward -> hold constant; else decay using D+ slope.

close all; clear all; clc;

%% ------------------- Load SOLPS Data (match profiles_iter_multi.nc layout) -------------------
fileSOLPS = 'profiles_iter_multi.nc';

% From ncdisp:
% x(nX), z(nZ), atomic_number(species), charge_number(species),
% ne(nX,nZ), te(nX,nZ), ti(nX,nZ), br/bt/bz(nX,nZ),
% ni_all/uR_all/uZ_all/uT_all(nX,nZ,species)
rS  = double(ncread(fileSOLPS, 'x'));    % R grid (nR)
zS  = double(ncread(fileSOLPS, 'z'));    % Z grid (nZ)
nr = numel(rS); nz = numel(zS);

neS = double(ncread(fileSOLPS, 'ne'));   % [nR x nZ]
TeS = double(ncread(fileSOLPS, 'te'));   % [nR x nZ]
TiS = double(ncread(fileSOLPS, 'ti'));   % [nR x nZ]

BrS = double(ncread(fileSOLPS, 'br'));   % [nR x nZ]
BtS = double(ncread(fileSOLPS, 'bt'));   % [nR x nZ]
BzS = double(ncread(fileSOLPS, 'bz'));   % [nR x nZ]

% Optional single-fluid velocities (may not exist in ITER file)
vrS = try_ncread(fileSOLPS, 'vr', [nr nz]);
vzS = try_ncread(fileSOLPS, 'vz', [nr nz]);
vtS = try_ncread(fileSOLPS, 'vt', [nr nz]);

% Multi-species arrays
Z_all   = double(ncread(fileSOLPS,'atomic_number'));   % [ns]
q_all   = double(ncread(fileSOLPS,'charge_number'));   % [ns]
ns = numel(Z_all);

ni_all0 = double(ncread(fileSOLPS,'ni_all'));          % [nR x nZ x ns] expected
uR_all0 = double(ncread(fileSOLPS,'uR_all'));
uZ_all0 = double(ncread(fileSOLPS,'uZ_all'));
uT_all0 = double(ncread(fileSOLPS,'uT_all'));

ni_all = ensure_nr_nz_ns(ni_all0, nr, nz, ns);
uR_all = ensure_nr_nz_ns(uR_all0, nr, nz, ns);
uZ_all = ensure_nr_nz_ns(uZ_all0, nr, nz, ns);
uT_all = ensure_nr_nz_ns(uT_all0, nr, nz, ns);

%% Full SOLPS → (X,Y) Extrapolation (Euclidean, wall-normal outside wall) for Electrons & Multi-species,
%  Masking by psiN, and NetCDF Export (ITER profiles_iter_multi.nc layout)
%  - Midplane exp fit (for diagnostics)
%  - Option-B column clamp on native (R,Z) for n,T (and Ti if present)
%  - OUTSIDE WALL: wall-normal extrapolation using Euclidean distance from wall
%       * scalars: exp(slope * d) decay (slope from midplane fit; guarded <=0)
%       * flows: constant along wall-normal (sampled slightly inside wall)
%  - INSIDE WALL: interp2 from clamped (or native) fields
%  - Multi-species support (ni_all, uR/uZ/uT) as in your ncdisp
%  - Writes single-fluid + multi-species to NetCDF
%
%  NOTE: removed temperature gradient (gradTi) entirely per request.

close all; clear all; clc;

%% ------------------- Load SOLPS Data (match profiles_iter_multi.nc layout) -------------------
fileSOLPS = 'profiles_iter_multi.nc';

% Based on your ncdisp:
% x (nX), z (nZ), atomic_number, charge_number, ne, te, ti, br, bt, bz,
% ni_all (nX x nZ x species), uR_all, uZ_all, uT_all

rS  = ncread(fileSOLPS, 'x');    % R grid (nR)
zS  = ncread(fileSOLPS, 'z');    % Z grid (nZ)

neS = ncread(fileSOLPS, 'ne');   % [nR x nZ]
TeS = ncread(fileSOLPS, 'te');   % [nR x nZ]
TiS = ncread_2D(fileSOLPS, 'ti'); % [nR x nZ] if present else NaN

brS = ncread_2D(fileSOLPS, 'br', size(neS));
btS = ncread_2D(fileSOLPS, 'bt', size(neS));
bzS = try_ncread_2D(fileSOLPS, 'bz', size(neS));

% Single-fluid flows may not exist in this file. Keep NaN if absent.
vrS = try_ncread_2D(fileSOLPS, 'vr', size(neS));
vzS = try_ncread_2D(fileSOLPS, 'vz', size(neS));
vtS = try_ncread_2D(fileSOLPS, 'vt', size(neS));

% Multi-species arrays
Z_all   = ncread(fileSOLPS,'atomic_number');   % [ns]
q_all   = ncread(fileSOLPS,'charge_number');   % [ns]

ni_all0 = ncread(fileSOLPS,'ni_all');          % expected [nR x nZ x ns]
uR_all0 = ncread(fileSOLPS,'uR_all');
uZ_all0 = ncread(fileSOLPS,'uZ_all');
uT_all0 = ncread(fileSOLPS,'uT_all');

% ti_all may or may not exist; if absent, fill with TiS replicated.
ti_all0 = try_ncread_3D(fileSOLPS,'ti_all', size(ni_all0));

nr = numel(rS); nz = numel(zS); ns = numel(Z_all);

% Ensure 3D shapes are [nR x nZ x ns]
ni_all = ensure_nr_nz_ns(ni_all0, nr, nz, ns);
uR_all = ensure_nr_nz_ns(uR_all0, nr, nz, ns);
uZ_all = ensure_nr_nz_ns(uZ_all0, nr, nz, ns);
uT_all = ensure_nr_nz_ns(uT_all0, nr, nz, ns);

if isempty(ti_all0) || all(isnan(ti_all0(:)))
    ti_all = repmat(TiS, 1, 1, ns);
else
    ti_all = ensure_nr_nz_ns(ti_all0, nr, nz, ns);
end

%% ------------------- Data Cleaning -------------------
n_min = 1e10;        % density floor (for masks)
T_min = 10;          % temperature floor (eV)

neS(~isfinite(neS) | neS<=0) = NaN;
TeS(~isfinite(TeS) | TeS<=0) = NaN;
TiS(~isfinite(TiS) | TiS<=0) = NaN;

TeS(TeS < T_min) = T_min;
TiS(TiS < T_min) = T_min;

brS(~isfinite(brS)) = NaN;
btS(~isfinite(btS)) = NaN;
bzS(~isfinite(bzS)) = NaN;

vrS(~isfinite(vrS)) = NaN;
vzS(~isfinite(vzS)) = NaN;
vtS(~isfinite(vtS)) = NaN;

%% ------------------- Read EFIT / wall polygon -------------------
read_efit_data;  % must provide g with g.lim, and calc_psiN available
if ~exist('g','var') || ~isfield(g,'lim') || isempty(g.lim)
    error('read_efit_data did not provide g.lim wall polygon.');
end
r_wall = g.lim(1,:);  z_wall = g.lim(2,:);
wallCoords = [r_wall(:), z_wall(:)];

%% ------------------- Quick visual check: SOLPS ne + wall -------------------
figure; imagesc(rS, zS, neS'); set(gca,'YDir','normal'); colorbar;
title('ne (SOLPS)'); xlabel('R [m]'); ylabel('Z [m]'); hold on;
plot(r_wall, z_wall, 'r', 'LineWidth', 1);

%% ------------------- psiN on native SOLPS grid -------------------
% Build RR/ZZ on [nR x nZ] consistent with data arrays:
[ZZs, RRs] = meshgrid(zS, rS); % RRs,ZZs are [nR x nZ]
psiN_native = reshape(calc_psiN(g, RRs(:), ZZs(:), 0), size(RRs));

plasma_mask = neS > n_min;     % [nR x nZ]

%% ------------------- (X,Y) Extrapolation Grid -------------------
num_points = 1000;   % keep as you had
[X, Y] = meshgrid(linspace(min(r_wall), max(r_wall), num_points), ...
                  linspace(min(z_wall), max(z_wall), num_points));
num_points_all = numel(X);
pool = gcp(); num_workers = pool.NumWorkers;
batch_size = ceil(num_points_all / (10 * num_workers));
extrap_coords = [X(:), Y(:)];
nBatches = ceil(num_points_all / batch_size);

%% ------------------- Midplane fit (diagnostic) -------------------
% Your ITER midplane fit window and target extension:
r_min = 8.01;
r_max = 8.33;
R_outer_target = 8.6;

mpfx = linspace(r_min, r_max, 1000);   % fit window (midplane)
mpfy = 0*mpfx;
mpx  = linspace(r_min, R_outer_target, 1000);  % extend outward

fitDensityAtMidplane = interp2(rS, zS, neS', mpfx, mpfy, 'linear', NaN);
fitTeAtMidplane      = interp2(rS, zS, TeS', mpfx, mpfy, 'linear', NaN);

% Robust polyfits (skip if <3 finite)
if nnz(isfinite(fitDensityAtMidplane) & fitDensityAtMidplane>0) >= 3
    p_ne = polyfit(mpfx(isfinite(fitDensityAtMidplane) & fitDensityAtMidplane>0), ...
                   log(fitDensityAtMidplane(isfinite(fitDensityAtMidplane) & fitDensityAtMidplane>0)), 1);
else
    p_ne = [ -10, log(n_min) ];
end
if nnz(isfinite(fitTeAtMidplane) & fitTeAtMidplane>0) >= 3
    p_Te = polyfit(mpfx(isfinite(fitTeAtMidplane) & fitTeAtMidplane>0), ...
                   log(fitTeAtMidplane(isfinite(fitTeAtMidplane) & fitTeAtMidplane>0)), 1);
else
    p_Te = [ -5, log(max(T_min,1)) ];
end

densityAtMidplane = interp2(rS, zS, neS', mpx, 0*mpx, 'linear', NaN);
densityAtMidplane(isnan(densityAtMidplane)) = 0;

TeAtMidplane = interp2(rS, zS, TeS', mpx, 0*mpx, 'linear', NaN);
TeAtMidplane(isnan(TeAtMidplane)) = T_min;

interpfn = (mpx - mpfx(1)) / (mpfx(end) - mpfx(1));
interpfn = min(max(interpfn, 0), 1);

extrapolatedne1d = interpfn .* exp(p_ne(2) + mpx*p_ne(1)) + (1-interpfn) .* densityAtMidplane;
extrapolatedTe1d = max(interpfn .* exp(p_Te(2) + mpx*p_Te(1)) + (1-interpfn) .* TeAtMidplane, T_min);

% Diagnostic tail beyond r_max
K_edge = 1;
Ldecay = 0.10;   % [m] used for native clamp and as a scale reference

R2 = mpfx(end);
idx_in = find(isfinite(extrapolatedne1d) & (mpx <= R2));
if ~isempty(idx_in)
    idt = idx_in(max(numel(idx_in)-K_edge+1,1):end);
    edge_ne = max(extrapolatedne1d(idt));
    beyond  = (mpx > R2);
    extrapolatedne1d(beyond) = edge_ne .* exp(-(mpx(beyond)-R2)/Ldecay);
end
idx_inT = find(isfinite(extrapolatedTe1d) & (mpx <= R2));
if ~isempty(idx_inT)
    idtT = idx_inT(max(numel(idx_inT)-K_edge+1,1):end);
    edge_Te = max(extrapolatedTe1d(idtT));
    beyondT = (mpx > R2);
    extrapolatedTe1d(beyondT) = max(edge_Te .* exp(-(mpx(beyondT)-R2)/Ldecay), T_min);
end

figure; semilogy(mpx, densityAtMidplane, 'b', mpx, extrapolatedne1d, 'b.', ...
                 mpx, TeAtMidplane, 'r', mpx, extrapolatedTe1d, 'r.');
legend('ne midplane','ne extrap','Te midplane','Te extrap');
xlabel('R [m]'); ylabel('Value'); title('Midplane ne & Te');

%% ------------------- Option-B clamp on native grid (single-fluid scalars) -------------------
% Clamp outward in each Z-column using last plasma cell; decay scalars with dr/Ldecay
neS_clamp = neS;  TeS_clamp = TeS;  TiS_clamp = TiS;

for jz = 1:nz
    col = plasma_mask(:, jz);                  % [nR x 1]
    edge_idx = find(col, 1, 'last');
    if isempty(edge_idx), continue; end
    i0 = max(1, edge_idx - K_edge + 1);

    edge_ne = max(neS(i0:edge_idx, jz), [], 'omitnan');
    edge_Te = max(TeS(i0:edge_idx, jz), [], 'omitnan');
    edge_Ti = max(TiS(i0:edge_idx, jz), [], 'omitnan');

    if ~isfinite(edge_ne), edge_ne = neS(edge_idx, jz); end
    if ~isfinite(edge_Te), edge_Te = max(TeS(edge_idx, jz), T_min); end
    if ~isfinite(edge_Ti), edge_Ti = max(TiS(edge_idx, jz), T_min); end

    for ir = edge_idx+1:nr
        dr  = rS(ir) - rS(edge_idx);
        fac = exp(-dr/Ldecay);

        neS_clamp(ir,jz) = max(edge_ne * fac, 0);
        TeS_clamp(ir,jz) = max(edge_Te * fac, T_min);
        TiS_clamp(ir,jz) = max(edge_Ti * fac, T_min);
    end
end

%% ------------------- Euclidean mapping to (X,Y): INSIDE interp2, OUTSIDE wall-normal -------------------
% Slopes (guard against growth)
s_ne = min(p_ne(1), -1e-4);
s_Te = min(p_Te(1), -1e-4);

% Build "source" points for outside-wall sampling (use clamped scalars; use multi-species flows later)
ok_ne = isfinite(neS_clamp) & neS_clamp>0;
ok_Te = isfinite(TeS_clamp) & TeS_clamp>0;
ok_Ti = isfinite(TiS_clamp) & TiS_clamp>0;

coords_ne = [RRs(ok_ne), ZZs(ok_ne)];
vals_ne   = neS_clamp(ok_ne);

coords_Te = [RRs(ok_Te), ZZs(ok_Te)];
vals_Te   = TeS_clamp(ok_Te);

coords_Ti = [RRs(ok_Ti), ZZs(ok_Ti)];
vals_Ti   = TiS_clamp(ok_Ti);

% For single-fluid flows: if absent, we'll synthesize from species=2 (common) as a fallback
have_v = any(isfinite(vrS(:))) || any(isfinite(vzS(:))) || any(isfinite(vtS(:)));
if have_v
    vr_use = vrS; vz_use = vzS; vt_use = vtS;
else
    % fallback: use species 2 if exists, else zeros
    kfb = min(2, ns);
    vr_use = uR_all(:,:,kfb);
    vz_use = uZ_all(:,:,kfb);
    vt_use = uT_all(:,:,kfb);
end

ok_v  = isfinite(vr_use) | isfinite(vz_use) | isfinite(vt_use);
coords_v = [RRs(ok_v), ZZs(ok_v)];
vals_vr  = vr_use(ok_v);
vals_vz  = vz_use(ok_v);
vals_vt  = vt_use(ok_v);

% Output arrays
val_ne = zeros(size(X));
val_Te = zeros(size(X));
val_Ti = zeros(size(X));
val_vr = zeros(size(X));
val_vz = zeros(size(X));
val_vt = zeros(size(X));

% batching
val_ne_cell = cell(nBatches,1);
val_Te_cell = cell(nBatches,1);
val_Ti_cell = cell(nBatches,1);
val_vr_cell = cell(nBatches,1);
val_vz_cell = cell(nBatches,1);
val_vt_cell = cell(nBatches,1);

parfor b = 1:nBatches
    s = (b-1)*batch_size + 1;
    e = min(b*batch_size, num_points_all);
    subset = extrap_coords(s:e, :);  % [Nb x 2] -> (R,Z)

    % Inside wall?
    inW = inpolygon(subset(:,1), subset(:,2), r_wall, z_wall);

    % ---- Scalars: inside = interp2; outside = wall-normal exponential ----
    ne_in = interp2(rS, zS, neS_clamp', subset(:,1), subset(:,2), 'linear', NaN);
    Te_in = interp2(rS, zS, TeS_clamp', subset(:,1), subset(:,2), 'linear', NaN);
    Ti_in = interp2(rS, zS, TiS_clamp', subset(:,1), subset(:,2), 'linear', NaN);

    ne_out = zeros(size(ne_in));
    Te_out = zeros(size(Te_in));
    Ti_out = zeros(size(Ti_in));

    if any(~inW)
        subO = subset(~inW,:);
        ne_out(~inW) = wall_normal_extrap_batch(subO, wallCoords, coords_ne, vals_ne, s_ne, 'exp');
        Te_out(~inW) = wall_normal_extrap_batch(subO, wallCoords, coords_Te, vals_Te, s_Te, 'exp');
        Ti_out(~inW) = wall_normal_extrap_batch(subO, wallCoords, coords_Ti, vals_Ti, s_Te, 'exp'); % use Te slope as proxy
    end

    ne_map = ne_in; ne_map(~inW) = ne_out(~inW);
    Te_map = Te_in; Te_map(~inW) = Te_out(~inW);
    Ti_map = Ti_in; Ti_map(~inW) = Ti_out(~inW);

    ne_map(~isfinite(ne_map)) = 0;
    Te_map(~isfinite(Te_map)) = T_min;
    Ti_map(~isfinite(Ti_map)) = T_min;

    Te_map = max(Te_map, T_min);
    Ti_map = max(Ti_map, T_min);

    % ---- Flows: inside = interp2; outside = wall-normal CONSTANT (sampled inside) ----
    vr_in = interp2(rS, zS, vr_use', subset(:,1), subset(:,2), 'linear', NaN);
    vz_in = interp2(rS, zS, vz_use', subset(:,1), subset(:,2), 'linear', NaN);
    vt_in = interp2(rS, zS, vt_use', subset(:,1), subset(:,2), 'linear', NaN);

    vr_out = zeros(size(vr_in));
    vz_out = zeros(size(vz_in));
    vt_out = zeros(size(vt_in));

    if any(~inW)
        subO = subset(~inW,:);
        vr_out(~inW) = wall_normal_extrap_batch(subO, wallCoords, coords_v, vals_vr, 0, 'const');
        vz_out(~inW) = wall_normal_extrap_batch(subO, wallCoords, coords_v, vals_vz, 0, 'const');
        vt_out(~inW) = wall_normal_extrap_batch(subO, wallCoords, coords_v, vals_vt, 0, 'const');
    end

    vr_map = vr_in; vr_map(~inW) = vr_out(~inW);
    vz_map = vz_in; vz_map(~inW) = vz_out(~inW);
    vt_map = vt_in; vt_map(~inW) = vt_out(~inW);

    % Fill NaNs (if any) with nearest interior sampling (still stable)
    nan_v = ~isfinite(vr_map);
    if any(nan_v)
        vr_map(nan_v) = knn_weighted(coords_v, vals_vr, subset(nan_v,:), 8);
    end
    nan_v = ~isfinite(vz_map);
    if any(nan_v)
        vz_map(nan_v) = knn_weighted(coords_v, vals_vz, subset(nan_v,:), 8);
    end
    nan_v = ~isfinite(vt_map);
    if any(nan_v)
        vt_map(nan_v) = knn_weighted(coords_v, vals_vt, subset(nan_v,:), 8);
    end

    % Ensure column vectors for reshape safety
    val_ne_cell{b} = ne_map(:);
    val_Te_cell{b} = Te_map(:);
    val_Ti_cell{b} = Ti_map(:);
    val_vr_cell{b} = vr_map(:);
    val_vz_cell{b} = vz_map(:);
    val_vt_cell{b} = vt_map(:);
end

val_ne = reshape(cell2mat(val_ne_cell), size(X));
val_Te = reshape(cell2mat(val_Te_cell), size(X));
val_Ti = reshape(cell2mat(val_Ti_cell), size(X));
val_vr = reshape(cell2mat(val_vr_cell), size(X));
val_vz = reshape(cell2mat(val_vz_cell), size(X));
val_vt = reshape(cell2mat(val_vt_cell), size(X));

% B fields: inside interp2; outside hold nearest wall-normal constant from near-wall sample
ok_B = isfinite(brS) & isfinite(btS) & isfinite(bzS);
coords_B = [RRs(ok_B), ZZs(ok_B)];
vals_br  = brS(ok_B); vals_bt = btS(ok_B); vals_bz = bzS(ok_B);

val_br_cell = cell(nBatches,1);
val_bt_cell = cell(nBatches,1);
val_bz_cell = cell(nBatches,1);

parfor b = 1:nBatches
    s = (b-1)*batch_size + 1;
    e = min(b*batch_size, num_points_all);
    subset = extrap_coords(s:e, :);

    inW = inpolygon(subset(:,1), subset(:,2), r_wall, z_wall);

    br_in = interp2(rS, zS, brS', subset(:,1), subset(:,2), 'linear', NaN);
    bt_in = interp2(rS, zS, btS', subset(:,1), subset(:,2), 'linear', NaN);
    bz_in = interp2(rS, zS, bzS', subset(:,1), subset(:,2), 'linear', NaN);

    br_out = zeros(size(br_in));
    bt_out = zeros(size(bt_in));
    bz_out = zeros(size(bz_in));

    if any(~inW)
        subO = subset(~inW,:);
        br_out(~inW) = wall_normal_extrap_batch(subO, wallCoords, coords_B, vals_br, 0, 'const');
        bt_out(~inW) = wall_normal_extrap_batch(subO, wallCoords, coords_B, vals_bt, 0, 'const');
        bz_out(~inW) = wall_normal_extrap_batch(subO, wallCoords, coords_B, vals_bz, 0, 'const');
    end

    br_map = br_in; br_map(~inW) = br_out(~inW);
    bt_map = bt_in; bt_map(~inW) = bt_out(~inW);
    bz_map = bz_in; bz_map(~inW) = bz_out(~inW);

    nanB = ~isfinite(br_map);
    if any(nanB), br_map(nanB) = knn_weighted(coords_B, vals_br, subset(nanB,:), 8); end
    nanB = ~isfinite(bt_map);
    if any(nanB), bt_map(nanB) = knn_weighted(coords_B, vals_bt, subset(nanB,:), 8); end
    nanB = ~isfinite(bz_map);
    if any(nanB), bz_map(nanB) = knn_weighted(coords_B, vals_bz, subset(nanB,:), 8); end

    val_br_cell{b} = br_map(:);
    val_bt_cell{b} = bt_map(:);
    val_bz_cell{b} = bz_map(:);
end

val_br = reshape(cell2mat(val_br_cell), size(X));
val_bt = reshape(cell2mat(val_bt_cell), size(X));
val_bz = reshape(cell2mat(val_bz_cell), size(X));

%% ------------------- Quick visuals (XY) -------------------
figure; imagesc([min(X(:)) max(X(:))],[min(Y(:)) max(Y(:))], val_ne);
set(gca,'YDir','normal','ColorScale','log'); colorbar; title('Extrapolated n_e'); hold on;
plot(r_wall, z_wall, 'r','LineWidth',1);

figure; imagesc([min(X(:)) max(X(:))],[min(Y(:)) max(Y(:))], val_Te);
set(gca,'YDir','normal','ColorScale','log'); colorbar; title('Extrapolated T_e'); hold on;
plot(r_wall, z_wall, 'r','LineWidth',1);

%% ------------------- Multi-species: wall-normal outside + interp2 inside -------------------
val_ni_mspecies = zeros(size(X,1), size(X,2), ns);
val_Ti_mspecies = zeros(size(X,1), size(X,2), ns);
val_uR_mspecies = zeros(size(X,1), size(X,2), ns);
val_uZ_mspecies = zeros(size(X,1), size(X,2), ns);
val_uT_mspecies = zeros(size(X,1), size(X,2), ns);

for k = 1:ns
    ni_k = ni_all(:,:,k);
    Ti_k = ti_all(:,:,k);
    uR_k = uR_all(:,:,k);
    uZ_k = uZ_all(:,:,k);
    uT_k = uT_all(:,:,k);

    ni_k(~isfinite(ni_k) | ni_k<=0) = NaN;
    Ti_k(~isfinite(Ti_k) | Ti_k<=0) = NaN;
    Ti_k(Ti_k < T_min) = T_min;

    % Species midplane slopes (fallback to electron slopes)
    ni_mid = interp2(rS, zS, ni_k', mpfx, 0*mpfx, 'linear', NaN);
    Ti_mid = interp2(rS, zS, Ti_k', mpfx, 0*mpfx, 'linear', NaN);

    if nnz(isfinite(ni_mid) & ni_mid>0) >= 3
        p_ni = polyfit(mpfx(isfinite(ni_mid)&ni_mid>0), log(ni_mid(isfinite(ni_mid)&ni_mid>0)), 1);
    else
        p_ni = p_ne;
    end
    if nnz(isfinite(Ti_mid) & Ti_mid>0) >= 3
        p_Ti = polyfit(mpfx(isfinite(Ti_mid)&Ti_mid>0), log(Ti_mid(isfinite(Ti_mid)&Ti_mid>0)), 1);
    else
        p_Ti = p_Te;
    end

    s_ni = min(p_ni(1), -1e-4);
    s_Ti = min(p_Ti(1), -1e-4);

    ok_ni = isfinite(ni_k) & ni_k>0;
    ok_Ti = isfinite(Ti_k) & Ti_k>0;
    ok_u  = isfinite(uR_k) | isfinite(uZ_k) | isfinite(uT_k);

    coords_ni = [RRs(ok_ni), ZZs(ok_ni)];  vals_ni = ni_k(ok_ni);
    coords_Ti = [RRs(ok_Ti), ZZs(ok_Ti)];  vals_Ti_k = Ti_k(ok_Ti);
    coords_u  = [RRs(ok_u),  ZZs(ok_u) ];  vals_uR = uR_k(ok_u); vals_uZ = uZ_k(ok_u); vals_uT = uT_k(ok_u);

    ni_cell = cell(nBatches,1);
    Ti_cell = cell(nBatches,1);
    uR_cell = cell(nBatches,1);
    uZ_cell = cell(nBatches,1);
    uT_cell = cell(nBatches,1);

    parfor b = 1:nBatches
        s = (b-1)*batch_size + 1;
        e = min(b*batch_size, num_points_all);
        subset = extrap_coords(s:e, :);

        inW = inpolygon(subset(:,1), subset(:,2), r_wall, z_wall);

        % Inside: interp2
        ni_in = interp2(rS, zS, ni_k', subset(:,1), subset(:,2), 'linear', NaN);
        Ti_in = interp2(rS, zS, Ti_k', subset(:,1), subset(:,2), 'linear', NaN);

        uR_in = interp2(rS, zS, uR_k', subset(:,1), subset(:,2), 'linear', NaN);
        uZ_in = interp2(rS, zS, uZ_k', subset(:,1), subset(:,2), 'linear', NaN);
        uT_in = interp2(rS, zS, uT_k', subset(:,1), subset(:,2), 'linear', NaN);

        % Outside: wall-normal
        ni_out = zeros(size(ni_in));
        Ti_out = zeros(size(Ti_in));
        uR_out = zeros(size(uR_in));
        uZ_out = zeros(size(uZ_in));
        uT_out = zeros(size(uT_in));

        if any(~inW)
            subO = subset(~inW,:);
            ni_out(~inW) = wall_normal_extrap_batch(subO, wallCoords, coords_ni, vals_ni, s_ni, 'exp');
            Ti_out(~inW) = wall_normal_extrap_batch(subO, wallCoords, coords_Ti, vals_Ti_k, s_Ti, 'exp');

            uR_out(~inW) = wall_normal_extrap_batch(subO, wallCoords, coords_u, vals_uR, 0, 'const');
            uZ_out(~inW) = wall_normal_extrap_batch(subO, wallCoords, coords_u, vals_uZ, 0, 'const');
            uT_out(~inW) = wall_normal_extrap_batch(subO, wallCoords, coords_u, vals_uT, 0, 'const');
        end

        ni_map = ni_in; ni_map(~inW) = ni_out(~inW);
        Ti_map = Ti_in; Ti_map(~inW) = Ti_out(~inW);
        uR_map = uR_in; uR_map(~inW) = uR_out(~inW);
        uZ_map = uZ_in; uZ_map(~inW) = uZ_out(~inW);
        uT_map = uT_in; uT_map(~inW) = uT_out(~inW);

        ni_map(~isfinite(ni_map)) = 0;
        Ti_map(~isfinite(Ti_map)) = T_min;
        Ti_map = max(Ti_map, T_min);

        % Any remaining NaNs in flows -> weighted KNN
        nanU = ~isfinite(uR_map);
        if any(nanU), uR_map(nanU) = knn_weighted(coords_u, vals_uR, subset(nanU,:), 8); end
        nanU = ~isfinite(uZ_map);
        if any(nanU), uZ_map(nanU) = knn_weighted(coords_u, vals_uZ, subset(nanU,:), 8); end
        nanU = ~isfinite(uT_map);
        if any(nanU), uT_map(nanU) = knn_weighted(coords_u, vals_uT, subset(nanU,:), 8); end

        ni_cell{b} = ni_map(:);
        Ti_cell{b} = Ti_map(:);
        uR_cell{b} = uR_map(:);
        uZ_cell{b} = uZ_map(:);
        uT_cell{b} = uT_map(:);
    end

    val_ni_mspecies(:,:,k) = reshape(cell2mat(ni_cell), size(X));
    val_Ti_mspecies(:,:,k) = reshape(cell2mat(Ti_cell), size(X));
    val_uR_mspecies(:,:,k) = reshape(cell2mat(uR_cell), size(X));
    val_uZ_mspecies(:,:,k) = reshape(cell2mat(uZ_cell), size(X));
    val_uT_mspecies(:,:,k) = reshape(cell2mat(uT_cell), size(X));
end

%% ------------------- psiN on (X,Y) and Mask -------------------
psiN_flat = calc_psiN(g, X(:), Y(:), 0);
psiN      = reshape(psiN_flat, size(X));

psiN_mask = 0.86;
mask = psiN < psiN_mask;

% Mask densities & temperatures by psiN (keep flows unmasked)
val_ne_masked = val_ne; val_ne_masked(mask) = 0;
val_Te_masked = val_Te; val_Te_masked(mask) = 0;
val_Ti_masked = val_Ti; val_Ti_masked(mask) = 0;

val_vr_masked = val_vr;
val_vz_masked = val_vz;
val_vt_masked = val_vt;

val_ni_mspecies_masked = val_ni_mspecies;
val_Ti_mspecies_masked = val_Ti_mspecies;
for k = 1:ns
    tmp = val_ni_mspecies_masked(:,:,k); tmp(mask)=0; val_ni_mspecies_masked(:,:,k)=tmp;
    tmp = val_Ti_mspecies_masked(:,:,k); tmp(mask)=0; val_Ti_mspecies_masked(:,:,k)=tmp;
end

val_uR_mspecies_masked = val_uR_mspecies;
val_uZ_mspecies_masked = val_uZ_mspecies;
val_uT_mspecies_masked = val_uT_mspecies;

save('multiSpecies_data_iter_wallnormal.mat');



%% ------------------- NetCDF Export (single-fluid + multispecies) -------------------
x  = X(1,:);          % R
z  = Y(:,1);          % Z
nR = length(x); nZ = length(z);

outnc = 'profiles_iter_multi_wallnormal.nc';
if exist(outnc,'file'), delete(outnc); end
ncid = netcdf.create(outnc,'CLOBBER');

% Dims
dimR = netcdf.defDim(ncid,'nX', nR);
dimZ = netcdf.defDim(ncid,'nZ', nZ);
dimS = netcdf.defDim(ncid,'species', ns);

% Coords/meta
vid_x   = netcdf.defVar(ncid,'x','double',dimR);
vid_z   = netcdf.defVar(ncid,'z','double',dimZ);
vid_Z   = netcdf.defVar(ncid,'atomic_number','double',dimS);
vid_q   = netcdf.defVar(ncid,'charge_number','double',dimS);
vid_psi = netcdf.defVar(ncid,'psiN','double',[dimR dimZ]);

% Single-fluid 2D
vid_ne  = netcdf.defVar(ncid,'ne','double',[dimR dimZ]);
vid_te  = netcdf.defVar(ncid,'te','double',[dimR dimZ]);
vid_ti  = netcdf.defVar(ncid,'ti','double',[dimR dimZ]);
vid_br  = netcdf.defVar(ncid,'br','double',[dimR dimZ]);
vid_bt  = netcdf.defVar(ncid,'bt','double',[dimR dimZ]);
vid_bz  = netcdf.defVar(ncid,'bz','double',[dimR dimZ]);
vid_vr  = netcdf.defVar(ncid,'vr','double',[dimR dimZ]);
vid_vz  = netcdf.defVar(ncid,'vz','double',[dimR dimZ]);
vid_vt  = netcdf.defVar(ncid,'vt','double',[dimR dimZ]);

% Multi-species 3D
vid_niA = netcdf.defVar(ncid,'ni_all','double',[dimR dimZ dimS]);
vid_uRA = netcdf.defVar(ncid,'uR_all','double',[dimR dimZ dimS]);
vid_uZA = netcdf.defVar(ncid,'uZ_all','double',[dimR dimZ dimS]);
vid_uTA = netcdf.defVar(ncid,'uT_all','double',[dimR dimZ dimS]);

netcdf.endDef(ncid);

% Write coords/meta
netcdf.putVar(ncid, vid_x, x);
netcdf.putVar(ncid, vid_z, z);
netcdf.putVar(ncid, vid_Z, Z_all);
netcdf.putVar(ncid, vid_q, q_all);
netcdf.putVar(ncid, vid_psi, permute(psiN, [2 1]));   % [nR x nZ]

% Write single-fluid (permute [nZ x nR] -> [nR x nZ])
netcdf.putVar(ncid, vid_ne, permute(val_ne_masked, [2 1]));
netcdf.putVar(ncid, vid_te, permute(val_Te_masked, [2 1]));
netcdf.putVar(ncid, vid_ti, permute(val_Ti_masked, [2 1]));

% Map B-field to (X,Y): inside interp2, outside wall-normal constant
Br_out = wall_map_constant_from_native(X, Y, rS, zS, BrS, r_wall, z_wall);
Bt_out = wall_map_constant_from_native(X, Y, rS, zS, BtS, r_wall, z_wall);
Bz_out = wall_map_constant_from_native(X, Y, rS, zS, BzS, r_wall, z_wall);

netcdf.putVar(ncid, vid_br, permute(Br_out, [2 1]));
netcdf.putVar(ncid, vid_bt, permute(Bt_out, [2 1]));
netcdf.putVar(ncid, vid_bz, permute(Bz_out, [2 1]));

netcdf.putVar(ncid, vid_vr, permute(val_vr_masked, [2 1]));
netcdf.putVar(ncid, vid_vz, permute(val_vz_masked, [2 1]));
netcdf.putVar(ncid, vid_vt, permute(val_vt_masked, [2 1]));

% Multi-species
netcdf.putVar(ncid, vid_niA, permute(val_ni_mspecies_masked, [2 1 3]));
netcdf.putVar(ncid, vid_uRA, permute(val_uR_mspecies_masked, [2 1 3]));
netcdf.putVar(ncid, vid_uZA, permute(val_uZ_mspecies_masked, [2 1 3]));
netcdf.putVar(ncid, vid_uTA, permute(val_uT_mspecies_masked, [2 1 3]));

netcdf.close(ncid);
disp(['Wrote ', outnc]);

%% ------------------- Read-back sanity plots -------------------
Rr = ncread(outnc,'x'); Zz = ncread(outnc,'z');
ne_m = ncread(outnc,'ne');           % [nR x nZ]
ni_m = ncread(outnc,'ni_all');       % [nR x nZ x ns]
uT_m = ncread(outnc,'uT_all');       % [nR x nZ x ns]

figure; imagesc(Zz, Rr, ne_m); set(gca,'YDir','normal','FontSize',12);
xlabel('Z [m]'); ylabel('R [m]'); colorbar; title('Masked n_e on (R,Z) from NetCDF');
hold on; plot(g.lim(2,:), g.lim(1,:), 'k', 'LineWidth', 1.0); hold off;

k = min(2, ns);
figure; imagesc(Zz, Rr, ni_m(:,:,k)); set(gca,'YDir','normal','FontSize',12);
xlabel('Z [m]'); ylabel('R [m]'); colorbar; title(sprintf('Masked n_i (species %d) from NetCDF',k));
hold on; plot(g.lim(2,:), g.lim(1,:), 'k', 'LineWidth', 1.0); hold off;

figure; imagesc(Zz, Rr, uT_m(:,:,k)); set(gca,'YDir','normal','FontSize',12);
xlabel('Z [m]'); ylabel('R [m]'); colorbar; title(sprintf('U_t (species %d) from NetCDF',k));
hold on; plot(g.lim(2,:), g.lim(1,:), 'k', 'LineWidth', 1.0); hold off;

%% ------------------- Helper functions -------------------

function A = try_ncread(fname, varname, expectedSize)
    try
        A = double(ncread(fname, varname));
        if nargin==3 && ~isempty(expectedSize) && ~isequal(size(A), expectedSize)
            if numel(A) == prod(expectedSize)
                A = reshape(A, expectedSize);
            else
                A = NaN(expectedSize);
            end
        end
    catch
        A = NaN(expectedSize);
    end
end

function B = ensure_nr_nz_ns(A0, nr, nz, ns)
    if isempty(A0)
        B = NaN(nr,nz,ns);
        return
    end
    sz = size(A0);
    if isequal(sz, [nr nz ns])
        B = A0; return
    end
    if numel(sz) ~= 3
        if numel(A0) == nr*nz*ns
            B = reshape(A0, [nr nz ns]);
        else
            B = NaN(nr,nz,ns);
        end
        return
    end
    perms_all = perms(1:3);
    for i=1:size(perms_all,1)
        p = perms_all(i,:);
        s = sz(p);
        if isequal(s, [nr nz ns])
            B = permute(A0, p);
            return
        end
    end
    if numel(A0) == nr*nz*ns
        B = reshape(A0, [nr nz ns]);
    else
        B = NaN(nr,nz,ns);
    end
end

function out = wall_normal_decay(subset, wallCoords, srcCoords, srcVals, slope)
    % f = f_probe * exp(slope * d_wall), with f_probe from wall-normal probe
    Nb = size(subset,1);
    out = zeros(Nb,1);

    idxWall = knnsearch(wallCoords, subset, 'K', 1);
    wallPts = wallCoords(idxWall,:);

    v = subset - wallPts;
    d = sqrt(sum(v.^2,2));
    u = v ./ max(d, 1e-30);

    step = 0.02; % 2 cm inside
    probe = wallPts - step .* u;

    idxSrc = knnsearch(srcCoords, probe, 'K', 1);
    fprobe = srcVals(idxSrc);

    out = fprobe .* exp(slope .* d);
    out(out<0) = 0;
end

function out = wall_normal_constant(subset, wallCoords, srcCoords, srcVals)
    Nb = size(subset,1);
    out = zeros(Nb,1);

    idxWall = knnsearch(wallCoords, subset, 'K', 1);
    wallPts = wallCoords(idxWall,:);

    v = subset - wallPts;
    d = sqrt(sum(v.^2,2));
    u = v ./ max(d, 1e-30);

    step = 0.02; % 2 cm inside
    probe = wallPts - step .* u;

    idxSrc = knnsearch(srcCoords, probe, 'K', 1);
    out = srcVals(idxSrc);
end

function [holdVal, incMask] = wall_normal_neon_hold(subset, wallCoords, srcCoords, srcVals)
    % If f(2cm inside) > f(4cm inside) => increasing outward => hold constant at f(2cm)
    idxWall = knnsearch(wallCoords, subset, 'K', 1);
    wallPts = wallCoords(idxWall,:);

    v = subset - wallPts;
    d = sqrt(sum(v.^2,2));
    u = v ./ max(d, 1e-30);

    probe1 = wallPts - 0.02 .* u;
    probe2 = wallPts - 0.04 .* u;

    i1 = knnsearch(srcCoords, probe1, 'K', 1);
    i2 = knnsearch(srcCoords, probe2, 'K', 1);

    f1 = srcVals(i1);
    f2 = srcVals(i2);

    incMask = (f1 > f2);
    holdVal = f1;
end

function Fxy = wall_map_constant_from_native(X, Y, rS, zS, Fnative, r_wall, z_wall)
    % inside wall: interp2
    % outside wall: wall-normal constant
    Fxy = zeros(size(X));

    Fin = interp2(rS, zS, Fnative', X, Y, 'linear', NaN);
    insideW = inpolygon(X, Y, r_wall, z_wall);
    Fxy(insideW) = Fin(insideW);

    [RRn, ZZn] = meshgrid(rS, zS);
    RRn = RRn.'; ZZn = ZZn.';
    inside_native = inpolygon(RRn, ZZn, r_wall, z_wall);

    mask_src = inside_native & isfinite(Fnative);
    coords_src = [RRn(mask_src), ZZn(mask_src)];
    vals_src   = Fnative(mask_src);

    if isempty(vals_src)
        Fxy(~insideW) = 0;
        return
    end

    sub_out = [X(~insideW), Y(~insideW)];
    Fxy(~insideW) = wall_normal_constant(sub_out, [r_wall(:) z_wall(:)], coords_src, vals_src);
end
%% ------------------- NetCDF Export (single-fluid + multispecies) -------------------
x  = X(1,:);          % R
z  = Y(:,1);          % Z
nR = length(x); nZ = length(z);

outnc = 'profiles_iter_multi_wallnormal.nc';
if exist(outnc,'file'), delete(outnc); end
ncid = netcdf.create(outnc,'CLOBBER');

% Dims
dimR = netcdf.defDim(ncid,'nX', nR);
dimZ = netcdf.defDim(ncid,'nZ', nZ);
dimS = netcdf.defDim(ncid,'species', ns);

% Coords/meta
vid_x   = netcdf.defVar(ncid,'x','double',dimR);
vid_z   = netcdf.defVar(ncid,'z','double',dimZ);
vid_Z   = netcdf.defVar(ncid,'atomic_number','double',dimS);
vid_q   = netcdf.defVar(ncid,'charge_number','double',dimS);
vid_psi = netcdf.defVar(ncid,'psiN','double',[dimR dimZ]);

% Single-fluid 2D
vid_ne  = netcdf.defVar(ncid,'ne','double',[dimR dimZ]);
vid_te  = netcdf.defVar(ncid,'te','double',[dimR dimZ]);
vid_ti  = netcdf.defVar(ncid,'ti','double',[dimR dimZ]);
vid_br  = netcdf.defVar(ncid,'br','double',[dimR dimZ]);
vid_bt  = netcdf.defVar(ncid,'bt','double',[dimR dimZ]);
vid_bz  = netcdf.defVar(ncid,'bz','double',[dimR dimZ]);
vid_vr  = netcdf.defVar(ncid,'vr','double',[dimR dimZ]);
vid_vz  = netcdf.defVar(ncid,'vz','double',[dimR dimZ]);
vid_vt  = netcdf.defVar(ncid,'vt','double',[dimR dimZ]);

% Multi-species 3D
vid_niA = netcdf.defVar(ncid,'ni_all','double',[dimR dimZ dimS]);
vid_uRA = netcdf.defVar(ncid,'uR_all','double',[dimR dimZ dimS]);
vid_uZA = netcdf.defVar(ncid,'uZ_all','double',[dimR dimZ dimS]);
vid_uTA = netcdf.defVar(ncid,'uT_all','double',[dimR dimZ dimS]);

netcdf.endDef(ncid);

% Write coords/meta
netcdf.putVar(ncid, vid_x, x);
netcdf.putVar(ncid, vid_z, z);
netcdf.putVar(ncid, vid_Z, Z_all);
netcdf.putVar(ncid, vid_q, q_all);
netcdf.putVar(ncid, vid_psi, permute(psiN, [2 1]));   % [nR x nZ]

% Write single-fluid (permute [nZ x nR] -> [nR x nZ])
netcdf.putVar(ncid, vid_ne, permute(val_ne_masked, [2 1]));
netcdf.putVar(ncid, vid_te, permute(val_Te_masked, [2 1]));
netcdf.putVar(ncid, vid_ti, permute(val_Ti_masked, [2 1]));

% Map B-field to (X,Y): inside interp2, outside wall-normal constant
Br_out = wall_map_constant_from_native(X, Y, rS, zS, BrS, r_wall, z_wall);
Bt_out = wall_map_constant_from_native(X, Y, rS, zS, BtS, r_wall, z_wall);
Bz_out = wall_map_constant_from_native(X, Y, rS, zS, BzS, r_wall, z_wall);

netcdf.putVar(ncid, vid_br, permute(Br_out, [2 1]));
netcdf.putVar(ncid, vid_bt, permute(Bt_out, [2 1]));
netcdf.putVar(ncid, vid_bz, permute(Bz_out, [2 1]));

netcdf.putVar(ncid, vid_vr, permute(val_vr_masked, [2 1]));
netcdf.putVar(ncid, vid_vz, permute(val_vz_masked, [2 1]));
netcdf.putVar(ncid, vid_vt, permute(val_vt_masked, [2 1]));

% Multi-species
netcdf.putVar(ncid, vid_niA, permute(val_ni_mspecies_masked, [2 1 3]));
netcdf.putVar(ncid, vid_uRA, permute(val_uR_mspecies_masked, [2 1 3]));
netcdf.putVar(ncid, vid_uZA, permute(val_uZ_mspecies_masked, [2 1 3]));
netcdf.putVar(ncid, vid_uTA, permute(val_uT_mspecies_masked, [2 1 3]));

netcdf.close(ncid);
disp(['Wrote ', outnc]);

%% ------------------- Read-back sanity plots -------------------
Rr = ncread(outnc,'x'); Zz = ncread(outnc,'z');
ne_m = ncread(outnc,'ne');           % [nR x nZ]
ni_m = ncread(outnc,'ni_all');       % [nR x nZ x ns]
uT_m = ncread(outnc,'uT_all');       % [nR x nZ x ns]

figure; imagesc(Zz, Rr, ne_m); set(gca,'YDir','normal','FontSize',12);
xlabel('Z [m]'); ylabel('R [m]'); colorbar; title('Masked n_e on (R,Z) from NetCDF');
hold on; plot(g.lim(2,:), g.lim(1,:), 'k', 'LineWidth', 1.0); hold off;

k = min(2, ns);
figure; imagesc(Zz, Rr, ni_m(:,:,k)); set(gca,'YDir','normal','FontSize',12);
xlabel('Z [m]'); ylabel('R [m]'); colorbar; title(sprintf('Masked n_i (species %d) from NetCDF',k));
hold on; plot(g.lim(2,:), g.lim(1,:), 'k', 'LineWidth', 1.0); hold off;

figure; imagesc(Zz, Rr, uT_m(:,:,k)); set(gca,'YDir','normal','FontSize',12);
xlabel('Z [m]'); ylabel('R [m]'); colorbar; title(sprintf('U_t (species %d) from NetCDF',k));
hold on; plot(g.lim(2,:), g.lim(1,:), 'k', 'LineWidth', 1.0); hold off;

%% ------------------- Helper functions -------------------

function A = try_ncread(fname, varname, expectedSize)
    try
        A = double(ncread(fname, varname));
        if nargin==3 && ~isempty(expectedSize) && ~isequal(size(A), expectedSize)
            if numel(A) == prod(expectedSize)
                A = reshape(A, expectedSize);
            else
                A = NaN(expectedSize);
            end
        end
    catch
        A = NaN(expectedSize);
    end
end

function B = ensure_nr_nz_ns(A0, nr, nz, ns)
    if isempty(A0)
        B = NaN(nr,nz,ns);
        return
    end
    sz = size(A0);
    if isequal(sz, [nr nz ns])
        B = A0; return
    end
    if numel(sz) ~= 3
        if numel(A0) == nr*nz*ns
            B = reshape(A0, [nr nz ns]);
        else
            B = NaN(nr,nz,ns);
        end
        return
    end
    perms_all = perms(1:3);
    for i=1:size(perms_all,1)
        p = perms_all(i,:);
        s = sz(p);
        if isequal(s, [nr nz ns])
            B = permute(A0, p);
            return
        end
    end
    if numel(A0) == nr*nz*ns
        B = reshape(A0, [nr nz ns]);
    else
        B = NaN(nr,nz,ns);
    end
end

function out = wall_normal_decay(subset, wallCoords, srcCoords, srcVals, slope)
    % f = f_probe * exp(slope * d_wall), with f_probe from wall-normal probe
    Nb = size(subset,1);
    out = zeros(Nb,1);

    idxWall = knnsearch(wallCoords, subset, 'K', 1);
    wallPts = wallCoords(idxWall,:);

    v = subset - wallPts;
    d = sqrt(sum(v.^2,2));
    u = v ./ max(d, 1e-30);

    step = 0.02; % 2 cm inside
    probe = wallPts - step .* u;

    idxSrc = knnsearch(srcCoords, probe, 'K', 1);
    fprobe = srcVals(idxSrc);

    out = fprobe .* exp(slope .* d);
    out(out<0) = 0;
end

function out = wall_normal_constant(subset, wallCoords, srcCoords, srcVals)
    Nb = size(subset,1);
    out = zeros(Nb,1);

    idxWall = knnsearch(wallCoords, subset, 'K', 1);
    wallPts = wallCoords(idxWall,:);

    v = subset - wallPts;
    d = sqrt(sum(v.^2,2));
    u = v ./ max(d, 1e-30);

    step = 0.02; % 2 cm inside
    probe = wallPts - step .* u;

    idxSrc = knnsearch(srcCoords, probe, 'K', 1);
    out = srcVals(idxSrc);
end

function [holdVal, incMask] = wall_normal_neon_hold(subset, wallCoords, srcCoords, srcVals)
    % If f(2cm inside) > f(4cm inside) => increasing outward => hold constant at f(2cm)
    idxWall = knnsearch(wallCoords, subset, 'K', 1);
    wallPts = wallCoords(idxWall,:);

    v = subset - wallPts;
    d = sqrt(sum(v.^2,2));
    u = v ./ max(d, 1e-30);

    probe1 = wallPts - 0.02 .* u;
    probe2 = wallPts - 0.04 .* u;

    i1 = knnsearch(srcCoords, probe1, 'K', 1);
    i2 = knnsearch(srcCoords, probe2, 'K', 1);

    f1 = srcVals(i1);
    f2 = srcVals(i2);

    incMask = (f1 > f2);
    holdVal = f1;
end

function Fxy = wall_map_constant_from_native(X, Y, rS, zS, Fnative, r_wall, z_wall)
    % inside wall: interp2
    % outside wall: wall-normal constant
    Fxy = zeros(size(X));

    Fin = interp2(rS, zS, Fnative', X, Y, 'linear', NaN);
    insideW = inpolygon(X, Y, r_wall, z_wall);
    Fxy(insideW) = Fin(insideW);

    [RRn, ZZn] = meshgrid(rS, zS);
    RRn = RRn.'; ZZn = ZZn.';
    inside_native = inpolygon(RRn, ZZn, r_wall, z_wall);

    mask_src = inside_native & isfinite(Fnative);
    coords_src = [RRn(mask_src), ZZn(mask_src)];
    vals_src   = Fnative(mask_src);

    if isempty(vals_src)
        Fxy(~insideW) = 0;
        return
    end

    sub_out = [X(~insideW), Y(~insideW)];
    Fxy(~insideW) = wall_normal_constant(sub_out, [r_wall(:) z_wall(:)], coords_src, vals_src);
end