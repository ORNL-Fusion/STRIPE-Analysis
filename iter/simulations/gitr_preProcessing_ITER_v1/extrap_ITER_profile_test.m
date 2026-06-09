%% Full SOLPS → (X,Y) Extrapolation (Euclidean) for Electrons & Multi-species (ITER file),
%  Masking by psiN, and NetCDF Export (GITR-style) + Cross-check plots
%
%  KEY POINTS (kept consistent with your workflow):
%  - Reads profiles_iter_multi.nc (x,z, ne,te,ti, br,bt,bz, ni_all, uR_all,uZ_all,uT_all)
%  - Midplane exp-fit for diagnostics in window [r_min,r_max], extrapolate to R_outer_target
%  - Native Option-B column clamp (radial tail) for scalars; velocities kept constant at edge
%  - (X,Y) mapping: scalars use exp(slope*d) with Euclidean distance to nearest anchor
%  - Velocities: **NO knn K=1 Voronoi tiles**. Uses scatteredInterpolant (smooth) + NN fallback
%  - Removes gradTi entirely (per request)
%  - Writes outnc with: x,z, atomic_number, charge_number, psiN, ne,te,ti, br,bt,bz,
%                       ni_all, uR_all,uZ_all,uT_all
%  - Adds sanity + read-back + midplane comparisons + Ne fractions / 2D maps checks
%
%  You set: r_min≈8.01, r_max≈8.33, extrapolate to 8.6

close all; clear; clc;

%% ------------------- Load SOLPS Data (match profiles_iter_multi.nc layout) -------------------
fileSOLPS = 'profiles_iter_wideGrid.nc';

% ncdisp shows:
% x (nX=400), z (nZ=800), atomic_number (16), charge_number (16),
% ne, te, ti, br, bt, bz, ni_all (nX x nZ x species),
% uR_all, uZ_all, uT_all

rgrid = ncread(fileSOLPS,'x');   % [nR]
zgrid = ncread(fileSOLPS,'z');   % [nZ]
ne_q  = ncread(fileSOLPS,'ne');  % [nR x nZ]
Te_q  = ncread(fileSOLPS,'te');  % [nR x nZ]
Ti_q  = ncread(fileSOLPS,'ti');  % [nR x nZ]
Br_q  = ncread(fileSOLPS,'br');  % [nR x nZ]
Bt_q  = ncread(fileSOLPS,'bt');  % [nR x nZ]
Bz_q  = ncread(fileSOLPS,'bz');  % [nR x nZ]

Z_all = ncread(fileSOLPS,'atomic_number');   % [ns]
q_all = ncread(fileSOLPS,'charge_number');   % [ns]

ni_all0 = ncread(fileSOLPS,'ni_all');        % expected [nR x nZ x ns]
uR_all0 = ncread(fileSOLPS,'uR_all');
uZ_all0 = ncread(fileSOLPS,'uZ_all');
uT_all0 = ncread(fileSOLPS,'uT_all');

nr = numel(rgrid);
nz = numel(zgrid);
ns = numel(Z_all);

% Ensure 3D arrays are [nr x nz x ns]
ni_all = ensure_nr_nz_ns(ni_all0, nr, nz, ns);
uR_all = ensure_nr_nz_ns(uR_all0, nr, nz, ns);
uZ_all = ensure_nr_nz_ns(uZ_all0, nr, nz, ns);
uT_all = ensure_nr_nz_ns(uT_all0, nr, nz, ns);

% convenience aliases used below (keep naming similar to prior script)
rS  = rgrid;
zS  = zgrid;
neS = ne_q;
TeS = Te_q;
TiS = Ti_q;

%% ------------------- Data Cleaning -------------------
n_min = 1e10;        % density floor (for masks)
T_min = 3;          % temperature floor (eV)
apply_psi_mask = true;
psiN_mask = 0.86;    % mask core region with psi_N < 0.86
adas_ne_dir = '/Users/78k/ORNL Dropbox/Atul Kumar/work/STRIPE-Analysis/iter/simulations/solps_iter/ADAS';
adas_scd_ne = fullfile(adas_ne_dir, 'scd89_ne.dat');
adas_acd_ne = fullfile(adas_ne_dir, 'acd89_ne.dat');

neS(~isfinite(neS) | neS<=0) = NaN;
TeS(~isfinite(TeS) | TeS<=0) = NaN;
TiS(~isfinite(TiS) | TiS<=0) = NaN;
TeS(TeS < T_min) = T_min;
TiS(TiS < T_min) = T_min;

Br_q(~isfinite(Br_q)) = 0;
Bt_q(~isfinite(Bt_q)) = 0;
Bz_q(~isfinite(Bz_q)) = 0;

%% ------------------- Read EFIT / wall polygon -------------------
read_efit_data;  % must create g with g.lim, g.r,g.z,g.psirz,g.ssimag,g.ssibry

if ~exist('g','var') || ~isfield(g,'lim') || isempty(g.lim)
    error('EFIT wall polygon not found. Ensure read_efit_data provides g.lim.');
end

%% ------------------- Quick visual sanity on SOLPS grid -------------------
figure; imagesc(rgrid, zgrid, neS'); set(gca,'YDir','normal'); colorbar;
title('ne (SOLPS)'); xlabel('R [m]'); ylabel('Z [m]'); hold on;
plot(g.lim(1,:), g.lim(2,:), 'r'); hold off;

%% ------------------- psiN on native SOLPS grid (optional diagnostics) -------------------
% Build psiN on native grid (not used for clamp, but used later for XY mask)
[ZZn, RRn] = meshgrid(zS, rS);             % size [nr x nz]
psiN_native = reshape(calc_psiN(g, RRn(:), ZZn(:), 0), size(RRn));  % [nr x nz]

%% ------------------- Masks on native grid -------------------
plasma_mask = isfinite(neS) & (neS > n_min);  % [nr x nz]
inside_limiter_mask = inpolygon(RRn, ZZn, g.lim(1,:), g.lim(2,:));  % [nr x nz]

%% ------------------- (X,Y) Extrapolation Grid -------------------
num_points = 1000;   % keep as in your last script
R_outer_target = 8.6;  % enforce extrapolation domain to this outer radius
R_xy_min = min(g.lim(1,:));
R_xy_max = max(max(g.lim(1,:)), R_outer_target);
Z_xy_min = min(g.lim(2,:));
Z_xy_max = max(g.lim(2,:));
[X, Y] = meshgrid(linspace(R_xy_min, R_xy_max, num_points), ...
                  linspace(Z_xy_min, Z_xy_max, num_points));
fprintf('XY extrapolation grid: R=[%.3f, %.3f], Z=[%.3f, %.3f]\n', ...
        R_xy_min, R_xy_max, Z_xy_min, Z_xy_max);
num_points_all = numel(X);

pool = gcp();
num_workers = pool.NumWorkers;
batch_size = ceil(num_points_all / (10 * num_workers));
extrap_coords = [X(:), Y(:)];
nBatches = ceil(num_points_all / batch_size);

%% ------------------- Midplane fit (diagnostic) -------------------
% ITER midplane window & outward extrapolation target
r_min = 8.01295;
r_max = 8.33246;

mpfx = linspace(r_min, r_max, 1000);   % fit window (midplane)
mpfy = 0*mpfx;
mpx  = linspace(r_min, R_outer_target, 1000);  % extend outward (antenna recess)

fitDensityAtMidplane = interp2(rS, zS, neS', mpfx, mpfy, 'linear', NaN);
fitTeAtMidplane      = interp2(rS, zS, TeS', mpfx, mpfy, 'linear', NaN);

% Robust polyfits (skip if <3 finite)
if nnz(isfinite(fitDensityAtMidplane) & fitDensityAtMidplane>0) >= 3
    p_ne = polyfit(mpfx(isfinite(fitDensityAtMidplane) & fitDensityAtMidplane>0), ...
                   log(fitDensityAtMidplane(isfinite(fitDensityAtMidplane) & fitDensityAtMidplane>0)), 1);
else
    p_ne = [ -10, log(n_min) ];   % safe fallback
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

extrapolatedne1d = interpfn .* exp(p_ne(2) + mpx*p_ne(1)) + (1 - interpfn) .* densityAtMidplane;
extrapolatedTe1d = max(interpfn .* exp(p_Te(2) + mpx*p_Te(1)) + (1 - interpfn) .* TeAtMidplane, T_min);

% Option-B tail on midplane (diagnostic)
K_edge = 1;       % robust edge value from last K points
Ldecay = 0.01;    % [m] clamp decay length

R2 = mpfx(end);   % r_max
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
xlabel('R [m]'); ylabel('Value'); title('Midplane ne & Te (diagnostic)');

%% ------------------- Option-B clamp on native grid (single-fluid scalars; no gradTi) -------------------
neS_clamp = neS;
TeS_clamp = TeS;
TiS_clamp = TiS;

% For single-fluid velocities: not in file; we keep only species velocities.
% (If you want electron/ion bulk, define later from species; not done here.)

for jz = 1:nz
    col = plasma_mask(:, jz);
    edge_idx = find(col, 1, 'last');
    if isempty(edge_idx), continue; end
    i0 = max(1, edge_idx - K_edge + 1);

    edge_ne = max(neS(i0:edge_idx, jz), [], 'omitnan');
    edge_Te = max(TeS(i0:edge_idx, jz), [], 'omitnan');
    edge_Ti = max(TiS(i0:edge_idx, jz), [], 'omitnan');

    if ~isfinite(edge_ne), edge_ne = neS(edge_idx,jz); end
    if ~isfinite(edge_Te), edge_Te = max(TeS(edge_idx,jz), T_min); end
    if ~isfinite(edge_Ti), edge_Ti = max(TiS(edge_idx,jz), T_min); end

    for ir = edge_idx+1:nr
        dr = rS(ir) - rS(edge_idx);
        fac = exp(-dr/Ldecay);
        neS_clamp(ir,jz) = max(edge_ne * fac, 0);
        TeS_clamp(ir,jz) = max(edge_Te * fac, T_min);
        TiS_clamp(ir,jz) = max(edge_Ti * fac, T_min);
    end
end

%% ------------------- Euclidean extrapolation to (X,Y) (single-fluid scalars only) -------------------
coords_ok = [RRn(:), ZZn(:)];  % native (R,Z) anchor coords, size [nr*nz x 2]

ok_ne_mask = isfinite(neS_clamp) & (neS_clamp>0);
ok_Te_mask = isfinite(TeS_clamp) & (TeS_clamp>0);
ok_Ti_mask = isfinite(TiS_clamp) & (TiS_clamp>0);
ok_T_common = ok_Te_mask & ok_Ti_mask;  % enforce same Te/Ti spatial support

idx_ne = find(ok_ne_mask(:)); coords_ne = coords_ok(idx_ne,:); vals_ne = neS_clamp(idx_ne);
idx_Te = find(ok_T_common(:)); coords_Te = coords_ok(idx_Te,:); vals_Te = TeS_clamp(idx_Te);
idx_Ti = find(ok_T_common(:)); coords_Ti = coords_ok(idx_Ti,:); vals_Ti = TiS_clamp(idx_Ti);

% slopes from midplane fit (guard against accidental growth)
s_ne = min(p_ne(1), -1e-4);
s_Te = min(p_Te(1), -1e-4);
% For Ti use Te slope as a reasonable fallback unless you want a separate fit
s_Ti = s_Te;

val_ne_cell = cell(nBatches,1);
val_Te_cell = cell(nBatches,1);
val_Ti_cell = cell(nBatches,1);

parfor b = 1:nBatches
    ibeg = (b-1)*batch_size + 1;
    iend = min(b*batch_size, num_points_all);
    subset = extrap_coords(ibeg:iend, :);  % [Nb x 2]

    % ne: nearest anchor + exp(slope*d)
    iN = knnsearch(coords_ne, subset, 'K', 1);
    dN = vecnorm(coords_ne(iN,:) - subset, 2, 2);
    ne_out = vals_ne(iN) .* exp(s_ne .* dN);

    % Te
    iT = knnsearch(coords_Te, subset, 'K', 1);
    dT = vecnorm(coords_Te(iT,:) - subset, 2, 2);
    Te_out = max(vals_Te(iT) .* exp(s_Te .* dT), T_min);

    % Ti
    iI = knnsearch(coords_Ti, subset, 'K', 1);
    dI = vecnorm(coords_Ti(iI,:) - subset, 2, 2);
    Ti_out = max(vals_Ti(iI) .* exp(s_Ti .* dI), T_min);

    val_ne_cell{b} = ne_out;
    val_Te_cell{b} = Te_out;
    val_Ti_cell{b} = Ti_out;
end

ne_full = cell2mat(val_ne_cell);
Te_full = cell2mat(val_Te_cell);
Ti_full = cell2mat(val_Ti_cell);

if numel(ne_full) ~= num_points_all
    error('Scalar mapping size mismatch: ne_full has %d, expected %d.', numel(ne_full), num_points_all);
end

val_ne = reshape(ne_full, size(X));
val_Te = reshape(Te_full, size(X));
val_Ti = reshape(Ti_full, size(X));

% Keep identical Te/Ti coverage footprint after extrapolation.
maskT = (val_Te > 0) & (val_Ti > 0) & isfinite(val_Te) & isfinite(val_Ti);
val_Te(~maskT) = 0;
val_Ti(~maskT) = 0;

figure; imagesc([min(X(:)) max(X(:))],[min(Y(:)) max(Y(:))], val_ne);
set(gca,'YDir','normal','ColorScale','log'); colorbar; title('Extrapolated n_e (X,Y)');
hold on; plot(g.lim(1,:), g.lim(2,:), 'r'); hold off;

figure; imagesc([min(X(:)) max(X(:))],[min(Y(:)) max(Y(:))], val_Te);
set(gca,'YDir','normal','ColorScale','log'); colorbar; title('Extrapolated T_e (X,Y)');
hold on; plot(g.lim(1,:), g.lim(2,:), 'r'); hold off;

%% ------------------- Multi-species: Option-B clamp + (X,Y) mapping -------------------
% Clamp scalars ni,Ti with radial decay.
% For velocities, do NOT pre-clamp radially (that creates horizontal bands).
% Instead, map from interior source points only using smooth interpolation.
val_ni_mspecies = zeros(size(X,1), size(X,2), ns);
val_uR_mspecies = zeros(size(X,1), size(X,2), ns);
val_uZ_mspecies = zeros(size(X,1), size(X,2), ns);
val_uT_mspecies = zeros(size(X,1), size(X,2), ns);

% wall polygon (for your own reference; not used in smoothing directly)
r_wall = g.lim(1,:); z_wall = g.lim(2,:);
wallCoords = [r_wall(:), z_wall(:)];
idxNe_rule = 6:min(16, ns);  % Ne0..Ne10+ in your current indexing
R_last_SOLPS = max(rS);      % hold Neon constant for R beyond this
iDplus = 2;                  % D+ species index (current ordering)

for k = 1:ns
    ni_k = ni_all(:,:,k);
    uR_k = uR_all(:,:,k);
    uZ_k = uZ_all(:,:,k);
    uT_k = uT_all(:,:,k);

    % Species scalar "Ti_k": not in ITER file; use single-fluid TiS as fallback
    Ti_k = TiS;

    % Species-specific midplane slopes (fallback to electron fits)
    ni_mid = interp2(rS, zS, ni_k', mpfx, 0*mpfx, 'linear', NaN);
    if nnz(isfinite(ni_mid) & ni_mid>0) >= 3
        p_ni = polyfit(mpfx(isfinite(ni_mid)&ni_mid>0), log(ni_mid(isfinite(ni_mid)&ni_mid>0)), 1);
    else
        p_ni = p_ne;
    end
    s_ni = min(p_ni(1), -1e-4);

    use_const_neon = any(k == idxNe_rule);
    % Ti slope fallback
    s_Ti_k = s_Te;

    % Column clamp on native scalars only
    ni_ex = ni_k;
    Ti_ex = Ti_k;

    for jz = 1:nz
        col = plasma_mask(:, jz);
        edge_idx = find(col, 1, 'last');
        if isempty(edge_idx), continue; end
        i0 = max(1, edge_idx - K_edge + 1);

        e_ni = max(ni_k(i0:edge_idx, jz), [], 'omitnan');
        e_Ti = max(Ti_k(i0:edge_idx, jz), [], 'omitnan');

        if ~isfinite(e_ni), e_ni = ni_k(edge_idx,jz); end
        if ~isfinite(e_Ti), e_Ti = max(Ti_k(edge_idx,jz), T_min); end

        for ir = edge_idx+1:nr
            dr = rS(ir) - rS(edge_idx);
            fac = exp(-dr/Ldecay);

            % scalars
            ni_ex(ir,jz) = max(e_ni * fac, 0);
            Ti_ex(ir,jz) = max(e_Ti * fac, T_min);
        end
    end

    % ---- Scalars on (X,Y): nearest anchor + exp(slope*d) ----
    ok_mask_ni = isfinite(ni_ex) & ni_ex>0;
    idx_ni = find(ok_mask_ni(:));
    if isempty(idx_ni)
        warning('Species k=%d has no finite ni points; skipping.', k);
        continue;
    end
    coords_ni = [RRn(idx_ni), ZZn(idx_ni)];
    vals_ni   = ni_ex(idx_ni);
    if use_const_neon
        F_ni_const = scatteredInterpolant(coords_ni(:,1), coords_ni(:,2), vals_ni, 'natural', 'nearest');
    end

    ni_cell = cell(nBatches,1);
    parfor b = 1:nBatches
        ibeg = (b-1)*batch_size + 1;
        iend = min(b*batch_size, num_points_all);
        subset = extrap_coords(ibeg:iend, :);

        if use_const_neon
            % Keep Neon density constant beyond last SOLPS R by clamping query R.
            rq = subset(:,1);
            zq = subset(:,2);
            rq(rq > R_last_SOLPS) = R_last_SOLPS;
            ni_out = F_ni_const(rq, zq);
        else
            iN = knnsearch(coords_ni, subset, 'K', 1);
            dN = vecnorm(coords_ni(iN,:) - subset, 2, 2);
            ni_out = vals_ni(iN) .* exp(s_ni .* dN);
        end

        ni_cell{b} = ni_out;
    end
    ni_full = cell2mat(ni_cell);
    if numel(ni_full) ~= num_points_all
        error('Species %d ni mapping size mismatch: got %d expected %d', k, numel(ni_full), num_points_all);
    end
    val_ni_mspecies(:,:,k) = reshape(ni_full, size(X));

    % ---- Velocities on (X,Y): smooth interpolation from interior anchors only ----
    src_mask_u = plasma_mask & inside_limiter_mask;
    if ~any(src_mask_u(:))
        % fallback: constant zeros
        val_uR_mspecies(:,:,k) = 0;
        val_uZ_mspecies(:,:,k) = 0;
        val_uT_mspecies(:,:,k) = 0;
    else
        uR_full = map_velocity_component(uR_k, src_mask_u, RRn, ZZn, extrap_coords);
        uZ_full = map_velocity_component(uZ_k, src_mask_u, RRn, ZZn, extrap_coords);
        uT_full = map_velocity_component(uT_k, src_mask_u, RRn, ZZn, extrap_coords);

        val_uR_mspecies(:,:,k) = reshape(uR_full, size(X));
        val_uZ_mspecies(:,:,k) = reshape(uZ_full, size(X));
        val_uT_mspecies(:,:,k) = reshape(uT_full, size(X));
    end
end
%% % Enforce constant Neon fraction beyond an anchor taken slightly inside
% the SOLPS edge to avoid edge spikes:
% n_Ne(R>R_anchor,z) = [n_Ne/ne]_anchor(z) * ne(R,z)
if ~isempty(idxNe_rule)
    iR_last_XY = find(X(1,:) <= R_last_SOLPS, 1, 'last');
    if isempty(iR_last_XY)
        iR_last_XY = 1;
    end
    neon_anchor_inward_cells = 12;   % anchor offset inside edge
    neon_anchor_halfwindow   = 4;    % radial averaging half-window around anchor
    iR_anchor = max(1, iR_last_XY - neon_anchor_inward_cells);
    iL = max(1, iR_anchor - neon_anchor_halfwindow);
    iU = min(size(X,2), iR_anchor + neon_anchor_halfwindow);

    outerCols = find(X(1,:) > X(1,iR_anchor));
    if ~isempty(outerCols)
        eps_ne = 1e-40;
        for kNe = idxNe_rule
            if kNe > ns, continue; end
            frac_win = val_ni_mspecies(:, iL:iU, kNe) ./ max(val_ne(:, iL:iU), eps_ne);
            frac_edge_z = median(frac_win, 2, 'omitnan');
            frac_edge_z(~isfinite(frac_edge_z)) = 0;
            frac_edge_z = max(frac_edge_z, 0);
            for j = outerCols
                val_ni_mspecies(:, j, kNe) = frac_edge_z .* val_ne(:, j);
            end
        end
        fprintf('Neon fraction anchor: R=%.4f m (window %.4f..%.4f m)\n', ...
                X(1,iR_anchor), X(1,iL), X(1,iU));
    end
end

% Enforce quasineutrality by adjusting D+:
% ne = nD+ + sum_q(q * nNe_q), using Neon charge states 0..10 at idx 6..16.
if iDplus <= ns && ~isempty(idxNe_rule)
    qNe = (0:numel(idxNe_rule)-1);               % Ne0..Ne10+
    ne_charge = zeros(size(val_ne));             % charged-Ne contribution
    for j = 1:numel(idxNe_rule)
        kNe = idxNe_rule(j);
        if kNe > ns, continue; end
        ne_charge = ne_charge + qNe(j) .* val_ni_mspecies(:,:,kNe);
    end

    % If charged Neon exceeds ne, scale Neon states down locally first.
    % This avoids large residuals after clipping D+ at zero.
    eps_ne = 1e-40;
    overMask = ne_charge > val_ne;
    if any(overMask(:))
        alpha = ones(size(val_ne));
        alpha(overMask) = val_ne(overMask) ./ max(ne_charge(overMask), eps_ne);
        alpha = min(max(alpha,0),1);
        for j = 1:numel(idxNe_rule)
            kNe = idxNe_rule(j);
            if kNe > ns, continue; end
            tmp = val_ni_mspecies(:,:,kNe);
            tmp(overMask) = tmp(overMask) .* alpha(overMask);
            val_ni_mspecies(:,:,kNe) = tmp;
        end

        % Recompute Neon charge after scaling
        ne_charge = zeros(size(val_ne));
        for j = 1:numel(idxNe_rule)
            kNe = idxNe_rule(j);
            if kNe > ns, continue; end
            ne_charge = ne_charge + qNe(j) .* val_ni_mspecies(:,:,kNe);
        end
    end

    nD_new = val_ne - ne_charge;
    nD_new(~isfinite(nD_new)) = 0;
    nD_new = max(nD_new, 0);
    val_ni_mspecies(:,:,iDplus) = nD_new;

    qn_err = val_ne - (val_ni_mspecies(:,:,iDplus) + ne_charge);
    fprintf('Quasi-neutrality after D+ reset: max|err|=%.3e, mean|err|=%.3e\n', ...
            max(abs(qn_err(:)),[],'omitnan'), mean(abs(qn_err(:)),'omitnan'));
end
%%
% ADAS-based coronal comparison at OUTER SOLPS MIDPLANE point (diagnostic only).
if ~isempty(idxNe_rule)
    [~, iz_mid_solps] = min(abs(zS - 0));
    ne_mid_solps = neS(:, iz_mid_solps);
    Te_mid_solps = TeS(:, iz_mid_solps);
    valid_mid = isfinite(ne_mid_solps) & isfinite(Te_mid_solps) & (ne_mid_solps > 0) & (Te_mid_solps > 0);

    if any(valid_mid)
        iR_outer_solps = find(valid_mid, 1, 'last');
        R_outer_solps = rS(iR_outer_solps);
        ne_outer = ne_mid_solps(iR_outer_solps);
        Te_outer = Te_mid_solps(iR_outer_solps);

        nNe_outer = zeros(1, numel(idxNe_rule));
        for j = 1:numel(idxNe_rule)
            kNe = idxNe_rule(j);
            if kNe <= ns
                nNe_outer(j) = ni_all(iR_outer_solps, iz_mid_solps, kNe);
            end
        end
        f_solps_ne = nNe_outer ./ max(sum(nNe_outer), 1e-40);

        try
            % Pure coronal equilibrium from screenshot Eq. (3): SCD/ACD only.
            f_cor_ne = get_ne_coronal_from_adas(Te_outer, ne_outer, adas_scd_ne, adas_acd_ne);
            f_cor_ne = f_cor_ne(1:numel(idxNe_rule));
            f_cor_ne = f_cor_ne ./ max(sum(f_cor_ne), 1e-40);

            devL1_cor = sum(abs(f_solps_ne - f_cor_ne));
            devL2_cor = norm(f_solps_ne - f_cor_ne);
            fprintf('Neon coronal comparison at outer SOLPS midplane: R=%.4f m, Z=%.4f m, Te=%.3f eV, ne=%.3e\n', ...
                    R_outer_solps, zS(iz_mid_solps), Te_outer, ne_outer);
            fprintf('SOLPS vs coronal (SCD/ACD only): L1=%.3e, L2=%.3e\n', devL1_cor, devL2_cor);

            q_plot = 0:(numel(idxNe_rule)-1);
            figure('Color','w','Position',[100 100 900 500]);
            hold on; box on; grid on;
            plot(q_plot, f_solps_ne, 'o-', 'LineWidth', 2, 'DisplayName', 'SOLPS Neon fraction');
            plot(q_plot, f_cor_ne, 's--', 'LineWidth', 2, 'DisplayName', 'Coronal (SCD/ACD only)');
            xlabel('Ne charge state q');
            ylabel('Fractional abundance');
            title(sprintf('Outer-midplane Neon fraction comparison (R=%.3f m, Z=%.3f m, T_e=%.3f eV, n_e=%.3e m^{-3})', ...
                R_outer_solps, zS(iz_mid_solps), Te_outer, ne_outer));
            legend('Location','best');
            text(0.02, 0.98, sprintf('Used for coronal: T_e=%.3f eV, n_e=%.3e m^{-3}', Te_outer, ne_outer), ...
                'Units','normalized', 'HorizontalAlignment','left', 'VerticalAlignment','top', ...
                'BackgroundColor','w', 'Margin',4, 'EdgeColor',[0.7 0.7 0.7]);

            % Additional coronal-only scan at fixed Te and selected ne values.
            ne_scan = [1e16, 1e17, 1e18];
            cscan = lines(numel(ne_scan));
            figure('Color','w','Position',[120 120 900 500]);
            hold on; box on; grid on;
            for ii = 1:numel(ne_scan)
                f_scan = get_ne_coronal_from_adas(Te_outer, ne_scan(ii), adas_scd_ne, adas_acd_ne);
                f_scan = f_scan(1:numel(idxNe_rule));
                f_scan = f_scan ./ max(sum(f_scan), 1e-40);
                plot(q_plot, f_scan, 'o-', 'LineWidth', 2, 'Color', cscan(ii,:), ...
                    'DisplayName', sprintf('Coronal n_e=%.0e m^{-3}', ne_scan(ii)));
            end
            xlabel('Ne charge state q');
            ylabel('Fractional abundance');
            title(sprintf('Coronal Neon equilibrium at fixed T_e=%.3f eV', Te_outer));
            legend('Location','best');
        catch ME_cor
            fprintf('Coronal check skipped: %s\n', ME_cor.message);
        end
    else
        fprintf('Coronal check skipped: no valid SOLPS midplane point found.\n');
    end
else
    fprintf('Coronal check skipped: Neon index set is empty.\n');
end

%% ------------------- psiN on (X,Y) and Mask -------------------
psiN = reshape(calc_psiN(g, X(:), Y(:), 0), size(X));

if apply_psi_mask
    mask = psiN < psiN_mask;
else
    mask = false(size(psiN));
end
fprintf('psiN mask: enabled=%d, threshold=%.2f, masked fraction=%.3f\n', ...
        apply_psi_mask, psiN_mask, nnz(mask)/numel(mask));

% Mask scalars (keep velocities unmasked)
val_ne_masked = val_ne; val_ne_masked(mask) = 0;
val_Te_masked = val_Te; val_Te_masked(mask) = 0;
val_Ti_masked = val_Ti; val_Ti_masked(mask) = 0;

val_ni_mspecies_masked = val_ni_mspecies;
for k = 1:ns
    tmp = val_ni_mspecies_masked(:,:,k); tmp(mask)=0; val_ni_mspecies_masked(:,:,k)=tmp;
end

% velocities (multi-species) unmasked
val_uR_mspecies_masked = val_uR_mspecies;
val_uZ_mspecies_masked = val_uZ_mspecies;
val_uT_mspecies_masked = val_uT_mspecies;

%% ------------------- NetCDF Export (GITR-style) -------------------
x = X(1,:);     % R axis (nRout)
z = Y(:,1);     % Z axis (nZout)
nRout = numel(x);
nZout = numel(z);

outnc = 'profiles_iter_multi_GITRstyle_XY.nc';
if exist(outnc,'file'), delete(outnc); end
ncid = netcdf.create(outnc,'CLOBBER');

% dims
dimR = netcdf.defDim(ncid,'x', nRout);
dimZ = netcdf.defDim(ncid,'z', nZout);
dimS = netcdf.defDim(ncid,'species', ns);

% vars
vid_x   = netcdf.defVar(ncid,'x','double',dimR);
vid_z   = netcdf.defVar(ncid,'z','double',dimZ);
vid_Z   = netcdf.defVar(ncid,'atomic_number','double',dimS);
vid_q   = netcdf.defVar(ncid,'charge_number','double',dimS);
vid_psi = netcdf.defVar(ncid,'psiN','double',[dimR dimZ]);

% single-fluid fields
vid_ne  = netcdf.defVar(ncid,'ne','double',[dimR dimZ]);
vid_te  = netcdf.defVar(ncid,'te','double',[dimR dimZ]);
vid_ti  = netcdf.defVar(ncid,'ti','double',[dimR dimZ]);

% B-fields (map native B to XY by interp2 on native; zero outside)
Br_xy = interp2(rS, zS, Br_q', X, Y, 'linear', 0);
Bt_xy = interp2(rS, zS, Bt_q', X, Y, 'linear', 0);
Bz_xy = interp2(rS, zS, Bz_q', X, Y, 'linear', 0);

vid_br  = netcdf.defVar(ncid,'br','double',[dimR dimZ]);
vid_bt  = netcdf.defVar(ncid,'bt','double',[dimR dimZ]);
vid_bz  = netcdf.defVar(ncid,'bz','double',[dimR dimZ]);

% multispecies
vid_niA = netcdf.defVar(ncid,'ni_all','double',[dimR dimZ dimS]);
vid_uRA = netcdf.defVar(ncid,'uR_all','double',[dimR dimZ dimS]);
vid_uZA = netcdf.defVar(ncid,'uZ_all','double',[dimR dimZ dimS]);
vid_uTA = netcdf.defVar(ncid,'uT_all','double',[dimR dimZ dimS]);

netcdf.endDef(ncid);

% write coords/meta
netcdf.putVar(ncid, vid_x, x);
netcdf.putVar(ncid, vid_z, z);
netcdf.putVar(ncid, vid_Z, Z_all);
netcdf.putVar(ncid, vid_q, q_all);
netcdf.putVar(ncid, vid_psi, permute(psiN, [2 1]));  % [dimR dimZ]

% write single-fluid masked (permute to [nRout x nZout])
netcdf.putVar(ncid, vid_ne, permute(val_ne_masked, [2 1]));
netcdf.putVar(ncid, vid_te, permute(val_Te_masked, [2 1]));
netcdf.putVar(ncid, vid_ti, permute(val_Ti_masked, [2 1]));

% B-fields
netcdf.putVar(ncid, vid_br, permute(Br_xy, [2 1]));
netcdf.putVar(ncid, vid_bt, permute(Bt_xy, [2 1]));
netcdf.putVar(ncid, vid_bz, permute(Bz_xy, [2 1]));

% write multi-species (permute to [nRout x nZout x ns])
netcdf.putVar(ncid, vid_niA, permute(val_ni_mspecies_masked, [2 1 3]));
netcdf.putVar(ncid, vid_uRA, permute(val_uR_mspecies_masked, [2 1 3]));
netcdf.putVar(ncid, vid_uZA, permute(val_uZ_mspecies_masked, [2 1 3]));
netcdf.putVar(ncid, vid_uTA, permute(val_uT_mspecies_masked, [2 1 3]));

netcdf.close(ncid);
fprintf('Wrote %s\n', outnc);

%% =============================================================================
%% === QUICK PLOTS (sanity) on SOLPS native grid ===
figure; imagesc(rgrid, zgrid, sqrt(Br_q.^2+Bt_q.^2+Bz_q.^2)');
set(gca,'YDir','normal'); axis equal tight; colorbar
title('|B| on SOLPS grid'); xlabel('R [m]'); ylabel('Z [m]'); hold on;
plot(g.lim(1,:), g.lim(2,:), 'k','LineWidth',1); hold off;

kplot = min(2,ns);
figure; imagesc(rgrid, zgrid, uT_all(:,:,kplot)');
set(gca,'YDir','normal'); axis equal tight; colorbar
title(sprintf('u_T (species %d) on SOLPS grid',kplot)); xlabel('R [m]'); ylabel('Z [m]');
hold on; plot(g.lim(1,:), g.lim(2,:), 'k','LineWidth',1); hold off;

%% =============================================================================
%% === VISUALIZATION CHECK: Read-back from GITR-style NetCDF ===
fprintf('\n=== Visualizing data from %s ===\n', outnc);

Rcheck = ncread(outnc,'x');
Zcheck = ncread(outnc,'z');

ne_chk  = ncread(outnc,'ne');
te_chk  = ncread(outnc,'te');
ti_chk  = ncread(outnc,'ti');
br_chk  = ncread(outnc,'br');
bt_chk  = ncread(outnc,'bt');
bz_chk  = ncread(outnc,'bz');

ni_chk  = ncread(outnc,'ni_all');
uR_chk  = ncread(outnc,'uR_all');
uZ_chk  = ncread(outnc,'uZ_all');
uT_chk  = ncread(outnc,'uT_all');

Z_all_chk = ncread(outnc,'atomic_number');
q_all_chk = ncread(outnc,'charge_number');

[nRchk, nZchk, ns_chk] = size(ni_chk);
fprintf('GITR file grid: nR=%d, nZ=%d, ns=%d\n', nRchk, nZchk, ns_chk);

Bmag_chk = sqrt(br_chk.^2 + bt_chk.^2 + bz_chk.^2);
Umag_chk = sqrt(uR_chk.^2 + uZ_chk.^2 + uT_chk.^2);

figure('Name','n_e from GITR NetCDF');
imagesc(Rcheck, Zcheck, ne_chk'); set(gca,'YDir','normal');
xlabel('R [m]'); ylabel('Z [m]'); title('Electron density n_e [m^{-3}]');
axis equal tight; colorbar;

figure('Name','|B| from GITR NetCDF');
imagesc(Rcheck, Zcheck, Bmag_chk'); set(gca,'YDir','normal');
xlabel('R [m]'); ylabel('Z [m]'); title('|B| [T]');
axis equal tight; colorbar;

k = min(7,ns_chk);
figure('Name','U_T from GITR NetCDF');
imagesc(Rcheck, Zcheck, uT_chk(:,:,k)'); set(gca,'YDir','normal');
xlabel('R [m]'); ylabel('Z [m]'); title(sprintf('U_T (species %d)',k));
axis equal tight; colorbar; hold on; plot(g.lim(1,:), g.lim(2,:), 'k'); hold off;

figure('Name','U_Z from GITR NetCDF');
imagesc(Rcheck, Zcheck, uZ_chk(:,:,k)'); set(gca,'YDir','normal');
xlabel('R [m]'); ylabel('Z [m]'); title(sprintf('U_Z (species %d)',k));
axis equal tight; colorbar; hold on; plot(g.lim(1,:), g.lim(2,:), 'k'); hold off;

figure('Name','U_R from GITR NetCDF');
imagesc(Rcheck, Zcheck, uR_chk(:,:,k)'); set(gca,'YDir','normal');
xlabel('R [m]'); ylabel('Z [m]'); title(sprintf('U_R (species %d)',k));
axis equal tight; colorbar; hold on; plot(g.lim(1,:), g.lim(2,:), 'k'); hold off;

fprintf('✅ GITR NetCDF visualization complete. Check figures for consistency.\n');

%% =============================================================================
%% === MIDPLANE COMPARISON: SOLPS vs GITR NetCDF (Z≈0) ===
fprintf('\n=== Comparing SOLPS and GITR data at midplane (Z=0) ===\n');

[~, iz_SOLPS] = min(abs(zgrid - 0));
[~, iz_GITR ] = min(abs(Zcheck - 0));

fprintf('Closest Z index in SOLPS = %d (Z=%.3f m)\n', iz_SOLPS, zgrid(iz_SOLPS));
fprintf('Closest Z index in GITR  = %d (Z=%.3f m)\n', iz_GITR, Zcheck(iz_GITR));

ne_SOLPS_mid   = ne_q(:, iz_SOLPS);   % note: native is [R x Z]
Te_SOLPS_mid   = Te_q(:, iz_SOLPS);
Ti_SOLPS_mid   = Ti_q(:, iz_SOLPS);
Bmag_SOLPS_mid = sqrt(Br_q(:,iz_SOLPS).^2 + Bt_q(:,iz_SOLPS).^2 + Bz_q(:,iz_SOLPS).^2);

ne_GITR_mid   = ne_chk(:, iz_GITR);
Te_GITR_mid   = te_chk(:, iz_GITR);
Ti_GITR_mid   = ti_chk(:, iz_GITR);
Bmag_GITR_mid = sqrt(br_chk(:,iz_GITR).^2 + bt_chk(:,iz_GITR).^2 + bz_chk(:,iz_GITR).^2);

ne_SOLPS_mid(ne_SOLPS_mid<=0)=NaN;
ne_GITR_mid(ne_GITR_mid<=0)=NaN;

figure('Color','w','Position',[100 100 1000 700]);
subplot(2,2,1)
plot(rgrid, ne_SOLPS_mid, 'b-', 'LineWidth', 2); hold on;
plot(Rcheck, ne_GITR_mid, 'r--', 'LineWidth', 2);
xlabel('R [m]'); ylabel('n_e [m^{-3}]'); legend('SOLPS','GITR','Location','best');
title('Electron Density at Z=0'); grid on;

subplot(2,2,2)
plot(rgrid, Te_SOLPS_mid, 'b-', 'LineWidth', 2); hold on;
plot(Rcheck, Te_GITR_mid, 'r--', 'LineWidth', 2);
xlabel('R [m]'); ylabel('T_e [eV]'); legend('SOLPS','GITR','Location','best');
title('Electron Temperature at Z=0'); grid on;

subplot(2,2,3)
plot(rgrid, Ti_SOLPS_mid, 'b-', 'LineWidth', 2); hold on;
plot(Rcheck, Ti_GITR_mid, 'r--', 'LineWidth', 2);
xlabel('R [m]'); ylabel('T_i [eV]'); legend('SOLPS','GITR','Location','best');
title('Ion Temperature at Z=0'); grid on;

subplot(2,2,4)
plot(rgrid, Bmag_SOLPS_mid, 'b-', 'LineWidth', 2); hold on;
plot(Rcheck, Bmag_GITR_mid, 'r--', 'LineWidth', 2);
xlabel('R [m]'); ylabel('|B| [T]'); legend('SOLPS','GITR','Location','best');
title('|B| at Z=0'); grid on;

sgtitle('SOLPS vs GITR Midplane Profiles (Z≈0)','FontSize',14);

fprintf('Mean |Δn_e|:  %.3e\n', nanmean(abs(ne_SOLPS_mid - interp1(Rcheck,ne_GITR_mid,rgrid,'linear','extrap'))));
fprintf('Mean |ΔT_e|:  %.3e\n', nanmean(abs(Te_SOLPS_mid - interp1(Rcheck,Te_GITR_mid,rgrid,'linear','extrap'))));
fprintf('Mean |ΔT_i|:  %.3e\n', nanmean(abs(Ti_SOLPS_mid - interp1(Rcheck,Ti_GITR_mid,rgrid,'linear','extrap'))));
fprintf('Mean |Δ|B||:  %.3e\n', nanmean(abs(Bmag_SOLPS_mid - interp1(Rcheck,Bmag_GITR_mid,rgrid,'linear','extrap'))));
fprintf('✅ Midplane comparison complete.\n');

%% =============================================================================
%% === MIDPLANE FRACTIONS: Neon on EXTRAPOLATED XY grid (to outer R) ===
fprintf('\n=== Plotting EXTRAPOLATED midplane n_i/n_e for Ne0–Ne10+ + total (Z≈0) ===\n');

idxNe = 6:16;
NeNames = {'Ne^{0}','Ne^{1+}','Ne^{2+}','Ne^{3+}','Ne^{4+}', ...
           'Ne^{5+}','Ne^{6+}','Ne^{7+}','Ne^{8+}','Ne^{9+}','Ne^{10+}'};

% Use extrapolated XY arrays so profiles extend to outer target radius.
% val_* arrays are [nZout x nRout], with axes z (rows), x (cols).
[~, iz_XY_mid] = min(abs(z - 0));
R_xy = x(:);
ne_mid = val_ne_masked(iz_XY_mid, :).';
ne_mid(ne_mid <= 0) = NaN;

% Right-hand region focus
R_focus_min = 8.0;
R_focus_max = max(R_xy);
maskR = (R_xy >= R_focus_min) & (R_xy <= R_focus_max);
Rzoom = R_xy(maskR);

% Compute total Ne density
nNe_total = zeros(size(ne_mid));
for j = 1:numel(idxNe)
    kNe = idxNe(j);
    if kNe > ns, continue; end
    nNe_total = nNe_total + squeeze(val_ni_mspecies_masked(iz_XY_mid, :, kNe)).';
end

figure('Color','w','Position',[100 100 1100 500]);
hold on; box on; grid on;
cmap = jet(numel(idxNe));

for j = 1:numel(idxNe)
    kNe = idxNe(j);
    if kNe > ns, continue; end
    n_i = squeeze(val_ni_mspecies_masked(iz_XY_mid, :, kNe)).';
    frac = n_i ./ ne_mid;
    plot(Rzoom, frac(maskR), 'LineWidth', 1.6, 'Color', cmap(j,:), 'DisplayName', NeNames{j});
end
frac_total = nNe_total ./ ne_mid;
plot(Rzoom, frac_total(maskR), 'k-', 'LineWidth', 2.8, 'DisplayName', 'Total Ne/n_e');

xlabel('R [m]'); ylabel('n_i / n_e');
title('Ne charge-state fractions at Z≈0 (EXTRAPOLATED XY)');
legend('Location','bestoutside'); legend boxoff; axis tight;

%% =============================================================================
%% === 2-D MAPS (SOLPS native) + overlays, similar to your diagnostics ===
fprintf('\n=== Plotting SOLPS 2-D maps (native): log n_e, log T_e, |B|, total Ne, avg Ne uT ===\n');

log_ne = log10(max(ne_q, 1));
log_Te = log10(max(Te_q, 1));
Bmag_native = sqrt(Br_q.^2 + Bt_q.^2 + Bz_q.^2);

% total Ne on native (R x Z)
nNe_total_2D = zeros(size(ne_q));
for j = 1:numel(idxNe)
    kNe = idxNe(j);
    if kNe > ns, continue; end
    nNe_total_2D = nNe_total_2D + ni_all(:,:,kNe);
end
log_nNe = log10(max(nNe_total_2D, 1));

% average Neon uT on native (use indices 7..16 like your snippet, if available)
Ua_Ne_sum = zeros(size(ne_q));
count = 0;
for kNe = 7:16
    if kNe > ns, continue; end
    Ua_Ne_sum = Ua_Ne_sum + uT_all(:,:,kNe);
    count = count + 1;
end
uT_Ne_avg = Ua_Ne_sum ./ max(count,1);

figure('Color','w','Position',[100 100 1500 800]);
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');

nexttile(1);
imagesc(rgrid, zgrid, log_ne'); set(gca,'YDir','normal'); axis equal tight; colorbar;
title('log_{10}(n_e)'); xlabel('R [m]'); ylabel('Z [m]');
hold on; plot(g.lim(1,:), g.lim(2,:), 'k','LineWidth',1.2); hold off;

nexttile(2);
imagesc(rgrid, zgrid, log_Te'); set(gca,'YDir','normal'); axis equal tight; colorbar;
title('log_{10}(T_e)'); xlabel('R [m]'); ylabel('Z [m]');
hold on; plot(g.lim(1,:), g.lim(2,:), 'k','LineWidth',1.2); hold off;

nexttile(3);
imagesc(rgrid, zgrid, Bmag_native'); set(gca,'YDir','normal'); axis equal tight; colorbar;
title('|B|'); xlabel('R [m]'); ylabel('Z [m]');
hold on; plot(g.lim(1,:), g.lim(2,:), 'k','LineWidth',1.2); hold off;

nexttile(4);
imagesc(rgrid, zgrid, log_nNe'); set(gca,'YDir','normal'); axis equal tight; colorbar;
title('log_{10}(Total Ne)'); xlabel('R [m]'); ylabel('Z [m]');
hold on; plot(g.lim(1,:), g.lim(2,:), 'k','LineWidth',1.2); hold off;

nexttile(5);
imagesc(rgrid, zgrid, Ti_q'); set(gca,'YDir','normal'); axis equal tight; colorbar;
title('T_i [eV]'); xlabel('R [m]'); ylabel('Z [m]');
hold on; plot(g.lim(1,:), g.lim(2,:), 'k','LineWidth',1.2); hold off;

nexttile(6);
imagesc(rgrid, zgrid, uT_Ne_avg'); set(gca,'YDir','normal'); axis equal tight; colorbar;
title('Avg Ne u_T (7..16)'); xlabel('R [m]'); ylabel('Z [m]');
hold on; plot(g.lim(1,:), g.lim(2,:), 'k','LineWidth',1.2); hold off;

sgtitle('SOLPS native diagnostics (log n_e, log T_e, |B|, total Ne, T_i, avg Ne u_T)','FontSize',16);
fprintf('✅ 2-D SOLPS native maps complete.\n');

%% =============================================================================
%% === 2-D MAPS (GITR NetCDF on XY grid) for the "flow looks weird" check ===
fprintf('\n=== Plotting GITR NetCDF (XY) maps to verify smooth velocity extrapolation ===\n');
kcheck = min(7,ns_chk);

figure('Color','w','Position',[100 100 1400 500]);
tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

nexttile(1);
imagesc(Rcheck, Zcheck, uR_chk(:,:,kcheck)'); set(gca,'YDir','normal'); axis equal tight; colorbar;
title(sprintf('u_R (species %d) from NetCDF',kcheck)); xlabel('R'); ylabel('Z');
hold on; plot(g.lim(1,:), g.lim(2,:), 'k','LineWidth',1); hold off;

nexttile(2);
imagesc(Rcheck, Zcheck, uZ_chk(:,:,kcheck)'); set(gca,'YDir','normal'); axis equal tight; colorbar;
title(sprintf('u_Z (species %d) from NetCDF',kcheck)); xlabel('R'); ylabel('Z');
hold on; plot(g.lim(1,:), g.lim(2,:), 'k','LineWidth',1); hold off;

nexttile(3);
imagesc(Rcheck, Zcheck, uT_chk(:,:,kcheck)'); set(gca,'YDir','normal'); axis equal tight; colorbar;
title(sprintf('u_T (species %d) from NetCDF',kcheck)); xlabel('R'); ylabel('Z');
hold on; plot(g.lim(1,:), g.lim(2,:), 'k','LineWidth',1); hold off;

fprintf('✅ If you still see Voronoi/triangular tiles, you are not using scatteredInterpolant for u*. This script does.\n');

%% =============================================================================
%% ------------------- Helper functions -------------------
function B = ensure_nr_nz_ns(A0, nr, nz, ns)
    if isempty(A0)
        B = NaN(nr,nz,ns);
        return;
    end
    sz = size(A0);
    if numel(sz) ~= 3
        error('Expected 3D array; got %s', mat2str(sz));
    end
    if isequal(sz, [nr nz ns])
        B = A0;
        return;
    end

    % try permutations that match [nr nz ns]
    perms_all = perms(1:3);
    for i = 1:size(perms_all,1)
        p = perms_all(i,:);
        if isequal(sz(p), [nr nz ns])
            B = permute(A0, p);
            return;
        end
    end

    % last resort: reshape if element count matches
    if numel(A0) == nr*nz*ns
        B = reshape(A0, [nr nz ns]);
    else
        error('Cannot coerce size %s into [%d %d %d].', mat2str(sz), nr, nz, ns);
    end
end

function v_full = map_velocity_component(v_native, src_mask_u, RRn, ZZn, extrap_coords)
    % Use only finite interior values for this component.
    mask = src_mask_u & isfinite(v_native);
    idx = find(mask(:));
    if isempty(idx)
        v_full = zeros(size(extrap_coords,1),1);
        return
    end

    coords = [RRn(idx), ZZn(idx)];
    vals = v_native(idx);

    % Smooth interpolation in-domain, NN fallback outside hull.
    F = scatteredInterpolant(coords(:,1), coords(:,2), vals, 'natural', 'none');
    v_full = F(extrap_coords(:,1), extrap_coords(:,2));

    nanMask = isnan(v_full);
    if any(nanMask)
        idxNN = knnsearch(coords, extrap_coords(nanMask,:), 'K', 1);
        v_full(nanMask) = vals(idxNN);
    end
end

function f = get_ne_coronal_from_adas(Te_eV, ne_m3, scd_file, acd_file)
    if ~isfile(scd_file)
        error('Missing ADAS SCD file: %s', scd_file);
    end
    if ~isfile(acd_file)
        error('Missing ADAS ACD file: %s', acd_file);
    end

    scd = read_adf11_blocks_local(scd_file);
    acd = read_adf11_blocks_local(acd_file);

    ne_cm3 = ne_m3 * 1e-6;
    logTe = log10(max(Te_eV, 1e-12));
    logNe = log10(max(ne_cm3, 1e-30));

    S = zeros(1,10);      % q=0..9
    A = zeros(1,10);      % q+1=1..10
    for q = 0:9
        z1 = q + 1;
        S(q+1) = interp_adf11_block_local(scd, z1, logTe, logNe);      % cm^3/s
        A(q+1) = interp_adf11_block_local(acd, z1, logTe, logNe);       % cm^3/s
    end

    % coronal fractions for Ne0..Ne10+
    g = zeros(1,11);
    g(1) = 1;
    epsv = 1e-80;
    for q = 0:9
        g(q+2) = g(q+1) * S(q+1) / max(A(q+1), epsv);
    end
    s = sum(g);
    if ~isfinite(s) || s <= 0
        f = zeros(1,11);
        f(1) = 1;
    else
        f = g / s;
    end
end

function val = interp_adf11_block_local(adf, z1, logTe, logNe)
    if z1 > numel(adf.blocks) || isempty(adf.blocks{z1})
        avail = find(~cellfun(@isempty, adf.blocks));
        error('ADF11 block Z1=%d missing in %s (available: %s)', z1, adf.filename, mat2str(avail));
    end
    M = adf.blocks{z1}; % log10(rate)
    v = bilinear_local(adf.logTe, adf.logNe, M, logTe, logNe);
    val = 10.^v;
end

function v = bilinear_local(xg, yg, M, x, y)
    [ix0, ix1, tx] = bracket_local(xg, x);
    [iy0, iy1, ty] = bracket_local(yg, y);
    v00 = M(ix0, iy0); v01 = M(ix0, iy1);
    v10 = M(ix1, iy0); v11 = M(ix1, iy1);
    v0 = v00*(1-ty) + v01*ty;
    v1 = v10*(1-ty) + v11*ty;
    v = v0*(1-tx) + v1*tx;
end

function [i0, i1, t] = bracket_local(g, x)
    if x <= g(1), i0 = 1; i1 = 1; t = 0; return; end
    if x >= g(end), i0 = numel(g); i1 = numel(g); t = 0; return; end
    i1 = find(g >= x, 1, 'first');
    i0 = i1 - 1;
    t = (x - g(i0)) / max(g(i1) - g(i0), 1e-30);
end

function adf = read_adf11_blocks_local(filename)
    txt = readlines(filename);
    txt = cellstr(txt(:));

    nums = sscanf(txt{1}, '%f');
    if numel(nums) < 5
        error('Cannot parse ADF11 header in %s', filename);
    end
    nBlocks = round(nums(1));
    nNe = round(nums(2));
    nTe = round(nums(3));

    i = 2;
    while i <= numel(txt) && (isempty(strtrim(txt{i})) || contains(txt{i}, '-'))
        i = i + 1;
    end

    logNe = [];
    while numel(logNe) < nNe
        if i > numel(txt), error('Unexpected EOF reading logNe in %s', filename); end
        vals = sscanf(txt{i}, '%f').';
        logNe = [logNe vals]; %#ok<AGROW>
        i = i + 1;
    end
    logNe = logNe(1:nNe);

    logTe = [];
    while numel(logTe) < nTe
        if i > numel(txt), error('Unexpected EOF reading logTe in %s', filename); end
        vals = sscanf(txt{i}, '%f').';
        logTe = [logTe vals]; %#ok<AGROW>
        i = i + 1;
    end
    logTe = logTe(1:nTe);

    blocks = cell(1, max(16, nBlocks));
    loaded = 0;
    while i <= numel(txt) && loaded < nBlocks
        ln = txt{i};
        if ~contains(ln, 'Z1=')
            i = i + 1;
            continue;
        end
        z1 = sscanf(extractAfter(ln, 'Z1='), '%f', 1);
        if isempty(z1) || ~isfinite(z1)
            tk = regexp(ln, 'Z1=\s*([0-9]+)', 'tokens', 'once');
            if isempty(tk)
                i = i + 1;
                continue;
            end
            z1 = str2double(tk{1});
        end
        z1 = round(z1);

        i = i + 1;
        M = zeros(nTe, nNe);
        for it = 1:nTe
            row = [];
            while numel(row) < nNe
                if i > numel(txt)
                    error('Unexpected EOF in block Z1=%d (%s)', z1, filename);
                end
                if contains(txt{i}, 'Z1=')
                    error('Malformed block Z1=%d in %s (early next header)', z1, filename);
                end
                vals = sscanf(txt{i}, '%f').';
                if ~isempty(vals)
                    row = [row vals]; %#ok<AGROW>
                end
                i = i + 1;
            end
            M(it,:) = row(1:nNe);
        end
        if z1 > numel(blocks), blocks{z1} = []; end
        blocks{z1} = M;
        loaded = loaded + 1;
    end

    adf.filename = filename;
    adf.logTe = logTe;
    adf.logNe = logNe;
    adf.blocks = blocks;
end
