%% Full SOLPS → (X,Y) Extrapolation (Euclidean) for Electrons & Multi-species,
%  Masking by psiN, and NetCDF Export
%  - Midplane exp fit (for diagnostics)
%  - Option-B column clamp on native (r,z)
%  - Euclidean exp(-d/Ldecay) propagation to (X,Y) for n,T and uR,uZ,uT
%  - Multi-species support (species dim first in file)
%  - Writes single-fluid + multi-species to NetCDF

close all; clear all; clc;

%% ------------------- Load SOLPS Data -------------------
fileSOLPS = 'interpolated_values_200882_multi.nc';

rS  = ncread(fileSOLPS, 'gridr');    % [nr]
zS  = ncread(fileSOLPS, 'gridz');    % [nz]
neS = ncread(fileSOLPS, 'ne');       % [nr x nz]
TeS = ncread(fileSOLPS, 'te');       % [nr x nz]
vrS = ncread(fileSOLPS, 'vr');       % [nr x nz]
vtS = ncread(fileSOLPS, 'vt');       % [nr x nz]
vzS = ncread(fileSOLPS, 'vz');       % [nr x nz]
gradTiS = -ncread(fileSOLPS,'gradTi'); % keep sign convention as in your scripts

% ---- Multi-species raw (species-first in file) ----
Z_all   = ncread(fileSOLPS,'atomic_number');   % [ns]
q_all   = ncread(fileSOLPS,'charge_number');   % [ns]
ni_all0 = ncread(fileSOLPS,'ni_all');          % 3D, species dim present
ti_all0 = ncread(fileSOLPS,'ti_all');
uR_all0 = ncread(fileSOLPS,'uR_all');
uZ_all0 = ncread(fileSOLPS,'uZ_all');
uT_all0 = ncread(fileSOLPS,'uT_all');

nr = numel(rS); nz = numel(zS); ns = numel(Z_all);

% Reorder to [nr x nz x ns] regardless of file ordering
perm_to_nr_nz_ns = @(A) permute(A, order_to_nr_nz_ns(size(A), nr, nz, ns));
ni_all = perm_to_nr_nz_ns(ni_all0);
ti_all = perm_to_nr_nz_ns(ti_all0);
uR_all = perm_to_nr_nz_ns(uR_all0);
uZ_all = perm_to_nr_nz_ns(uZ_all0);
uT_all = perm_to_nr_nz_ns(uT_all0);

%% ------------------- Data Cleaning -------------------
n_min = 1e10;        % density floor (for masks)
T_min = 10;          % temperature floor (eV)
clip0 = @(A) (A.*(A>0));     % helper to zero-out negatives

neS(~isfinite(neS) | neS<=0) = NaN;
TeS(~isfinite(TeS) | TeS<=0) = NaN;
vrS(~isfinite(vrS)) = NaN; vzS(~isfinite(vzS)) = NaN; vtS(~isfinite(vtS)) = NaN;
gradTiS(~isfinite(gradTiS)) = NaN;
TeS(TeS < T_min) = T_min;

% Small visuals (unchanged)
figure; imagesc(rS, zS, neS');  set(gca,'YDir','normal'); colorbar; title('ne (SOLPS)');  hold on; 
read_efit_data; plot(g.lim(1,:), g.lim(2,:), 'r');

%% ------------------- psiN on native SOLPS grid -------------------
[RR_native, ZZ_native] = meshgrid(rS, zS);   % [nz x nr] after meshgrid
R_vec = RR_native(:).'; Z_vec = ZZ_native(:).';
[psiN_vec, ~] = calc_psiN(g, R_vec, Z_vec, 0);
psiN_native = reshape(psiN_vec, size(RR_native));   % [nz x nr]

% Build masks used for column-wise edge detection
[ZZ_nr, RR_nr] = meshgrid(zS, rS);     % RR_nr,ZZ_nr are [nr x nz] (matching array shapes below)
plasma_mask = neS > n_min;             % [nr x nz]
inside_limiter_mask = inpolygon(RR_nr, ZZ_nr, g.lim(1,:), g.lim(2,:));

%% ------------------- (X,Y) Extrapolation Grid -------------------
num_points = 1000;   % keep as in your last script
[X, Y] = meshgrid(linspace(min(g.lim(1,:)), max(g.lim(1,:)), num_points), ...
                  linspace(min(g.lim(2,:)), max(g.lim(2,:)), num_points));
num_points_all = numel(X);
pool = gcp(); num_workers = pool.NumWorkers;
batch_size = ceil(num_points_all / (10 * num_workers));
extrap_coords = [X(:), Y(:)];
nBatches = ceil(num_points_all / batch_size);

%% ------------------- Midplane fit (diagnostic) -------------------
% Shot 196154 midplane window & full R range (as in your scripts)
mpfx = linspace(2.28122, 2.3156, 1000); % DIII-D # 200882
% mpfx = linspace(2.14333, 2.20688, 1000); % DIII-D 196154   
mpfy = 0*mpfx;
% mpx  = linspace(2.14333, 2.5,    1000);
mpx  = linspace(2.28122, 2.5,    1000);

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

densityAtMidplane = interp2(rS, zS, neS', mpx, 0*mpx, 'linear', NaN); densityAtMidplane(isnan(densityAtMidplane)) = 0;
TeAtMidplane      = interp2(rS, zS, TeS', mpx, 0*mpx, 'linear', NaN); TeAtMidplane(isnan(TeAtMidplane)) = T_min;

interpfn = (mpx - mpfx(1)) / (mpfx(end) - mpfx(1)); interpfn = min(max(interpfn,0),1);
extrapolatedne1d = interpfn .* exp(p_ne(2) + mpx*p_ne(1)) + (1-interpfn) .* densityAtMidplane;
extrapolatedTe1d = max(interpfn .* exp(p_Te(2) + mpx*p_Te(1)) + (1-interpfn) .* TeAtMidplane, T_min);

% Option-B tail on midplane (diagnostic)
K_edge = 1;       % use last K points for robust edge value
Ldecay = 0.10;    % [m] Euclidean decay length (used below as well)

% OLD:
% Ldecay_v = 1;    % [m] Euclidean decay length (used below as well)
% NEW: infinite/constant velocity by construction (we'll also remove exp() later)
Ldecay_v = 1e6;    % [m] (kept for reference; no longer used in exp for v)

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

%% ------------------- Option-B clamp on native grid (single-fluid) -------------------
% Column-wise: clamp outward (flat at edge), then exponential tail (dr / Ldecay)
neS_clamp = neS;  TeS_clamp = TeS;
vrS_clamp = vrS;  vzS_clamp = vzS;  vtS_clamp = vtS;
gradTi_clamp = gradTiS;

for jz = 1:nz
    col = plasma_mask(:, jz);
    edge_idx = find(col, 1, 'last');
    if isempty(edge_idx), continue; end
    i0 = max(1, edge_idx - K_edge + 1);

    % Edge values (robust)
    edge_ne  = max(neS(i0:edge_idx, jz), [], 'omitnan');
    edge_Te  = max(TeS(i0:edge_idx, jz), [], 'omitnan');

    % OLD: mean over last K cells (can be ~0 at sheath)
    % e_vr     = mean(vrS(i0:edge_idx, jz), 'omitnan');
    % e_vz     = mean(vzS(i0:edge_idx, jz), 'omitnan');
    % e_vt     = mean(vtS(i0:edge_idx, jz), 'omitnan');

    % NEW: use last finite edge-cell values
    e_vr = vrS(edge_idx, jz);
    e_vz = vzS(edge_idx, jz);
    e_vt = vtS(edge_idx, jz);

    e_gTi    = mean(gradTiS(i0:edge_idx, jz), 'omitnan');

    if ~isfinite(edge_ne), edge_ne = neS(edge_idx,jz); end
    if ~isfinite(edge_Te), edge_Te = max(TeS(edge_idx,jz), T_min); end
    if ~isfinite(e_vr),    e_vr = vrS(edge_idx,jz); end
    if ~isfinite(e_vz),    e_vz = vzS(edge_idx,jz); end
    if ~isfinite(e_vt),    e_vt = vtS(edge_idx,jz); end
    if ~isfinite(e_gTi),   e_gTi = gradTiS(edge_idx,jz); end

    for ir = edge_idx+1:nr
        % OLD: skip outside limiter
        % if ~inside_limiter_mask(ir, jz), continue; end
        dr = rS(ir) - rS(edge_idx);

        % Scalars still decay:
        fac = exp(-dr/Ldecay);

        % OLD (decayed velocities):
        % vrS_clamp(ir,jz)  = e_vr  * exp(-dr/Ldecay_v);
        % vzS_clamp(ir,jz)  = e_vz  * exp(-dr/Ldecay_v);
        % vtS_clamp(ir,jz)  = e_vt  * exp(-dr/Ldecay_v);

        % NEW (constant velocities to the wall):
        vrS_clamp(ir,jz)  = e_vr;
        vzS_clamp(ir,jz)  = e_vz;
        vtS_clamp(ir,jz)  = e_vt;

        % Scalars (unchanged behavior):
        neS_clamp(ir,jz)  = max(edge_ne * fac, 0);
        TeS_clamp(ir,jz)  = max(edge_Te * fac, T_min);
        gradTi_clamp(ir,jz) = e_gTi * fac;
    end
end

%% ------------------- Euclidean extrapolation to (X,Y) (single-fluid) -------------------
% Build ok points from clamped fields
[ZZs, RRs] = meshgrid(zS, rS); % [nr x nz] shaped mapping
coords_ok = [RRs(:), ZZs(:)];  % native (R,Z) anchor coords

ok_ne_mask = isfinite(neS_clamp) & (neS_clamp>0);
ok_Te_mask = isfinite(TeS_clamp) & (TeS_clamp>0);
ok_v_mask  = isfinite(vrS_clamp) | isfinite(vzS_clamp) | isfinite(vtS_clamp);

idx_ne = find(ok_ne_mask(:));   coords_ne = coords_ok(idx_ne,:);   okValues_ne = neS_clamp(idx_ne);
idx_Te = find(ok_Te_mask(:));   coords_Te = coords_ok(idx_Te,:);   okValues_Te = TeS_clamp(idx_Te);
idx_v  = find(ok_v_mask(:));    coords_v  = coords_ok(idx_v,:);    vals_vr = vrS_clamp(idx_v); vals_vz = vzS_clamp(idx_v); vals_vt = vtS_clamp(idx_v);

% Use midplane exponential slopes (1/m) from earlier robust fits
% Guard against accidental growth (in case a fit returns >= 0)
s_ne = min(p_ne(1), -1e-4);    % slope of ln(ne) vs R
s_Te = min(p_Te(1), -1e-4);    % slope of ln(Te) vs R

% Parallel batches
val_ne_cell = cell(nBatches,1);
val_Te_cell = cell(nBatches,1);
val_vr_cell = cell(nBatches,1);
val_vz_cell = cell(nBatches,1);
val_vt_cell = cell(nBatches,1);

parfor b = 1:nBatches
    s = (b-1)*batch_size + 1; 
    e = min(b*batch_size, num_points_all);
    subset = extrap_coords(s:e, :);  % [Nb x 2] target (R,Z)

    % ---- electrons: slope-based exponential with Euclidean distance ----
    % ne
    idxN = knnsearch(coords_ne, subset, 'K', 1);
    dN   = vecnorm(coords_ne(idxN,:) - subset, 2, 2);   % Euclidean distance [m]
    ne_out = okValues_ne(idxN) .* exp(s_ne .* dN);

    % Te
    idxT = knnsearch(coords_Te, subset, 'K', 1);
    dT   = vecnorm(coords_Te(idxT,:) - subset, 2, 2);
    Te_out = max(okValues_Te(idxT) .* exp(s_Te .* dT), T_min);

    % ---- velocities: keep your constant-nearest behavior (no exponential) ----
    idxV  = knnsearch(coords_v, subset, 'K', 1);
    vr_out = vals_vr(idxV);
    vz_out = vals_vz(idxV);
    vt_out = vals_vt(idxV);

    val_ne_cell{b} = ne_out;
    val_Te_cell{b} = Te_out;
    val_vr_cell{b} = vr_out;
    val_vz_cell{b} = vz_out;
    val_vt_cell{b} = vt_out;
end

val_ne = reshape(cell2mat(val_ne_cell), size(X));
val_Te = reshape(cell2mat(val_Te_cell), size(X));
val_vr = reshape(cell2mat(val_vr_cell), size(X));
val_vz = reshape(cell2mat(val_vz_cell), size(X));
val_vt = reshape(cell2mat(val_vt_cell), size(X));

% gradTi: already clamped on native; map linearly
val_gradTi = interp2(rS, zS, gradTi_clamp', X, Y, 'linear', 0);

% Quick visuals (unchanged style)
figure; imagesc([min(X(:)) max(X(:))],[min(Y(:)) max(Y(:))], val_ne);
set(gca,'YDir','normal','ColorScale','log'); colorbar; title('Extrapolated n_e'); hold on; plot(g.lim(1,:), g.lim(2,:), 'r');

figure; imagesc([min(X(:)) max(X(:))],[min(Y(:)) max(Y(:))], val_Te);
set(gca,'YDir','normal','ColorScale','log'); colorbar; title('Extrapolated T_e'); hold on; plot(g.lim(1,:), g.lim(2,:), 'r');

%% ------------------- Multi-species: Option-B clamp + Euclidean to (X,Y) -------------------
val_ni_mspecies = zeros(size(X,1), size(X,2), ns);
val_Ti_mspecies = zeros(size(X,1), size(X,2), ns);
val_uR_mspecies = zeros(size(X,1), size(X,2), ns);
val_uZ_mspecies = zeros(size(X,1), size(X,2), ns);
val_uT_mspecies = zeros(size(X,1), size(X,2), ns);

for k = 1:ns
    ni_k = ni_all(:,:,k); Ti_k = ti_all(:,:,k);
    uR_k = uR_all(:,:,k); uZ_k = uZ_all(:,:,k); uT_k = uT_all(:,:,k);

    % Midplane (diagnostic polyfit only)
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

    % Column-wise Option-B clamp on native
    ni_ex = ni_k; Ti_ex = Ti_k;
    uR_ex = uR_k; uZ_ex = uZ_k; uT_ex = uT_k;
    for jz = 1:nz
        col = plasma_mask(:, jz);
        edge_idx = find(col, 1, 'last');
        if isempty(edge_idx), continue; end
        i0 = max(1, edge_idx - K_edge + 1);

        e_ni = max(ni_k(i0:edge_idx, jz), [], 'omitnan');
        e_Ti = max(Ti_k(i0:edge_idx, jz), [], 'omitnan');

        % OLD (means):
        % e_uR = mean(uR_k(i0:edge_idx, jz), 'omitnan');
        % e_uZ = mean(uZ_k(i0:edge_idx, jz), 'omitnan');
        % e_uT = mean(uT_k(i0:edge_idx, jz), 'omitnan');

        % NEW (last finite edge cell):
        e_uR = uR_k(edge_idx, jz);
        e_uZ = uZ_k(edge_idx, jz);
        e_uT = uT_k(edge_idx, jz);

        if ~isfinite(e_ni), e_ni = ni_k(edge_idx,jz); end
        if ~isfinite(e_Ti), e_Ti = max(Ti_k(edge_idx,jz), T_min); end
        if ~isfinite(e_uR), e_uR = uR_k(edge_idx,jz); end
        if ~isfinite(e_uZ), e_uZ = uZ_k(edge_idx,jz); end
        if ~isfinite(e_uT), e_uT = uT_k(edge_idx,jz); end

        for ir = edge_idx+1:nr
            % OLD: limiter check
            % if ~inside_limiter_mask(ir, jz), continue; end
            dr = rS(ir) - rS(edge_idx);
            fac = exp(-dr/Ldecay);  % for ni, Ti

            % OLD (decay velocities):
            % uR_ex(ir, jz) = e_uR * exp(-dr/Ldecay_v);
            % uZ_ex(ir, jz) = e_uZ * exp(-dr/Ldecay_v);
            % uT_ex(ir, jz) = e_uT * exp(-dr/Ldecay_v);

            % NEW (constant velocities):
            uR_ex(ir, jz) = e_uR;
            uZ_ex(ir, jz) = e_uZ;
            uT_ex(ir, jz) = e_uT;

            % Scalars:
            ni_ex(ir, jz) = max(e_ni * fac, 0);
            Ti_ex(ir, jz) = max(e_Ti * fac, T_min);
        end
    end

       % ===== Euclidean propagation to (X,Y): species-specific slopes =====
    ok_mask_ni = isfinite(ni_ex) & ni_ex>0;
    ok_mask_Ti = isfinite(Ti_ex) & Ti_ex>0;
    ok_mask_u  = isfinite(uR_ex) | isfinite(uZ_ex) | isfinite(uT_ex);

    idx_ni = find(ok_mask_ni(:)); 
    idx_Ti = find(ok_mask_Ti(:)); 
    idx_u  = find(ok_mask_u(:));

    if isempty(idx_ni)
        warning('Species k=%d has no finite ni points; skipping.', k); 
        continue;
    end

    coords_ok = [RRs(idx_ni), ZZs(idx_ni)];          % anchor points for ni
    vals_ni   = ni_ex(idx_ni);

    coords_T  = [RRs(idx_Ti), ZZs(idx_Ti)];          % anchor points for Ti
    vals_Ti   = Ti_ex(idx_Ti);

    coords_u  = [RRs(idx_u),  ZZs(idx_u)];           % anchor points for u
    vals_uR   = uR_ex(idx_u); 
    vals_uZ   = uZ_ex(idx_u); 
    vals_uT   = uT_ex(idx_u);

    % Use species-specific mid-plane slopes (guard against ≥0)
    s_ni = min(p_ni(1), -1e-4);   % slope of ln(ni) vs R for this species
    s_Ti = min(p_Ti(1), -1e-4);   % slope of ln(Ti) vs R for this species

    ni_cell = cell(nBatches,1); 
    Ti_cell = cell(nBatches,1);
    uR_cell = cell(nBatches,1); 
    uZ_cell = cell(nBatches,1); 
    uT_cell = cell(nBatches,1);

    parfor b = 1:nBatches
        s = (b-1)*batch_size + 1; 
        e = min(b*batch_size, num_points_all);
        subset = extrap_coords(s:e, :);    % (R,Z)

        % ---- ni: slope-based exponential in Euclidean distance ----
        iN = knnsearch(coords_ok, subset, 'K', 1);
        dN = vecnorm(coords_ok(iN,:) - subset, 2, 2);
        ni_out = vals_ni(iN) .* exp(s_ni .* dN);     % decays with distance
        %% 
        %% 

        % ---- Ti: slope-based exponential, clamped to T_min ----
        iT = knnsearch(coords_T, subset, 'K', 1);
        dT = vecnorm(coords_T(iT,:) - subset, 2, 2);
        Ti_out = max(vals_Ti(iT) .* exp(s_Ti .* dT), T_min);

        % ---- velocities: constant (nearest) — no exponential ----
        iU = knnsearch(coords_u, subset, 'K', 1);
        uR_out = vals_uR(iU);
        uZ_out = vals_uZ(iU);
        uT_out = vals_uT(iU);

        ni_cell{b} = ni_out;
        Ti_cell{b} = Ti_out;
        uR_cell{b} = uR_out;
        uZ_cell{b} = uZ_out;
        uT_cell{b} = uT_out;
    end

    val_ni_mspecies(:,:,k) = reshape(cell2mat(ni_cell), size(X));
    val_Ti_mspecies(:,:,k) = reshape(cell2mat(Ti_cell), size(X));
    val_uR_mspecies(:,:,k) = reshape(cell2mat(uR_cell), size(X));
    val_uZ_mspecies(:,:,k) = reshape(cell2mat(uZ_cell), size(X));
    val_uT_mspecies(:,:,k) = reshape(cell2mat(uT_cell), size(X));
end

% Optional quick check for carbon (Z=6)
isC = (Z_all == 6);
if any(isC)
    kc = find(isC,1,'first');
    figure; imagesc([min(X(:)) max(X(:))],[min(Y(:)) max(Y(:))], val_uT_mspecies(:,:,kc));
    set(gca,'YDir','normal'); colorbar; title(sprintf('n_i on (X,Y), Carbon (k=%d)',kc)); hold on; plot(g.lim(1,:), g.lim(2,:), 'r');
end

%% ------------------- psiN on (X,Y) and Mask -------------------
psiN_flat = calc_psiN(g, X(:), Y(:), []);
psiN      = reshape(psiN_flat, size(X));

psiN_mask = 0.86;      % your earlier choice (change if needed)
mask = psiN < psiN_mask;

% Single-fluid masked
val_ne_masked      = val_ne;      val_ne_masked(mask)      = 0;
val_Te_masked      = val_Te;      val_Te_masked(mask)      = 0;
val_gradTi_masked  = val_gradTi;  val_gradTi_masked(mask)  = 0;

% OLD: mask velocities too
% val_vr_masked      = val_vr;      val_vr_masked(mask)      = 0;
% val_vz_masked      = val_vz;      val_vz_masked(mask)      = 0;
% val_vt_masked      = val_vt;      val_vt_masked(mask)      = 0;

% NEW: keep velocities unmasked by psiN
val_vr_masked = val_vr;
val_vz_masked = val_vz;
val_vt_masked = val_vt;

% Multi-species masked
val_ni_mspecies_masked = val_ni_mspecies;
val_Ti_mspecies_masked = val_Ti_mspecies;

% OLD: mask multi-species velocities
% val_uR_mspecies_masked = val_uR_mspecies;
% val_uZ_mspecies_masked = val_uZ_mspecies;
% val_uT_mspecies_masked = val_uT_mspecies;
% for k = 1:ns
%     tmp = val_uR_mspecies_masked(:,:,k); tmp(mask)=0; val_uR_mspecies_masked(:,:,k)=tmp;
%     tmp = val_uZ_mspecies_masked(:,:,k); tmp(mask)=0; val_uZ_mspecies_masked(:,:,k)=tmp;
%     tmp = val_uT_mspecies_masked(:,:,k); tmp(mask)=0; val_uT_mspecies_masked(:,:,k)=tmp;
% end

% NEW: do not mask multi-species velocities
val_uR_mspecies_masked = val_uR_mspecies;
val_uZ_mspecies_masked = val_uZ_mspecies;
val_uT_mspecies_masked = val_uT_mspecies;

% Still mask densities and temperatures by psiN
for k = 1:ns
    tmp = val_ni_mspecies_masked(:,:,k); tmp(mask)=0; val_ni_mspecies_masked(:,:,k)=tmp;
    tmp = val_Ti_mspecies_masked(:,:,k); tmp(mask)=0; val_Ti_mspecies_masked(:,:,k)=tmp;
end

save('multiSpecies_data_200882_test.mat')
%% ------------------- NetCDF Export (single-fluid + multispecies) -------------------
x  = X(1,:);          % R
z  = Y(:,1);          % Z
nR = length(x); nZ = length(z);

outnc = 'profilesDIIID_200882_multi_test.nc';
ncid = netcdf.create(outnc,'CLOBBER');

% Dims
dimR = netcdf.defDim(ncid,'x', nR);
dimZ = netcdf.defDim(ncid,'z', nZ);
dimS = netcdf.defDim(ncid,'species', ns);

% Coords / meta
vid_x   = netcdf.defVar(ncid,'x','double',dimR);
vid_z   = netcdf.defVar(ncid,'z','double',dimZ);
vid_Z   = netcdf.defVar(ncid,'atomic_number','double',dimS);
vid_q   = netcdf.defVar(ncid,'charge_number','double',dimS);
vid_psi = netcdf.defVar(ncid,'psiN','double',[dimR dimZ]);

% Single-fluid 2D
vid_ne  = netcdf.defVar(ncid,'ne',     'double',[dimR dimZ]);
vid_te  = netcdf.defVar(ncid,'te',     'double',[dimR dimZ]);
vid_gti = netcdf.defVar(ncid,'gradTi', 'double',[dimR dimZ]);
vid_vr  = netcdf.defVar(ncid,'vr',     'double',[dimR dimZ]);
vid_vz  = netcdf.defVar(ncid,'vz',     'double',[dimR dimZ]);
vid_vt  = netcdf.defVar(ncid,'vt',     'double',[dimR dimZ]);

% Multi-species 3D (R,Z,S)
vid_niA = netcdf.defVar(ncid,'ni_all', 'double',[dimR dimZ dimS]);
vid_tiA = netcdf.defVar(ncid,'ti_all', 'double',[dimR dimZ dimS]);
vid_uRA = netcdf.defVar(ncid,'uR_all', 'double',[dimR dimZ dimS]);
vid_uZA = netcdf.defVar(ncid,'uZ_all', 'double',[dimR dimZ dimS]);
vid_uTA = netcdf.defVar(ncid,'uT_all', 'double',[dimR dimZ dimS]);

netcdf.endDef(ncid);

% Write coords/meta
netcdf.putVar(ncid, vid_x, x);
netcdf.putVar(ncid, vid_z, z);
netcdf.putVar(ncid, vid_Z, Z_all);
netcdf.putVar(ncid, vid_q, q_all);
% note: permute (X,Y)->(R,Z) = (columns,rows) i.e., [nZ x nR] -> [nR x nZ]
netcdf.putVar(ncid, vid_psi, permute(psiN, [2 1]));

% Write single-fluid masked (permute to [nR x nZ])
netcdf.putVar(ncid, vid_ne,  permute(val_ne_masked,     [2 1]));
netcdf.putVar(ncid, vid_te,  permute(val_Te_masked,     [2 1]));
netcdf.putVar(ncid, vid_gti, permute(val_gradTi_masked, [2 1]));
netcdf.putVar(ncid, vid_vr,  permute(val_vr_masked,     [2 1]));
netcdf.putVar(ncid, vid_vz,  permute(val_vz_masked,     [2 1]));
netcdf.putVar(ncid, vid_vt,  permute(val_vt_masked,     [2 1]));

% Write multi-species masked (permute to [nR x nZ x ns])
netcdf.putVar(ncid, vid_niA, permute(val_ni_mspecies_masked, [2 1 3]));
netcdf.putVar(ncid, vid_tiA, permute(val_Ti_mspecies_masked, [2 1 3]));
netcdf.putVar(ncid, vid_uRA, permute(val_uR_mspecies_masked, [2 1 3]));
netcdf.putVar(ncid, vid_uZA, permute(val_uZ_mspecies_masked, [2 1 3]));
netcdf.putVar(ncid, vid_uTA, permute(val_uT_mspecies_masked, [2 1 3]));

netcdf.close(ncid);
disp(['Wrote ', outnc]);

%% ------------------- Read-back sanity plots -------------------
Rr = ncread(outnc,'x'); Zz = ncread(outnc,'z');
ne_m = ncread(outnc,'ne');
te_m = ncread(outnc,'te');
ni_all = ncread(outnc,'ni_all');
uT_all = ncread(outnc,'uT_all');

figure; imagesc(Zz, Rr, ne_m); set(gca,'YDir','normal','FontSize',12);
xlabel('Z [m]'); ylabel('R [m]'); colorbar; title('Masked n_e on (R,Z) from NetCDF');

figure; imagesc(Rr, Zz, ni_all(:,:,6)'); set(gca,'YDir','normal','FontSize',12);
xlabel('Z [m]'); ylabel('R [m]'); colorbar; title('Masked n_e on (R,Z) from NetCDF');
 hold on;
    plot(g.lim(1,:), g.lim(2,:), 'k', 'LineWidth', 1);
    contour(X, Y, psiN, [psiN_mask psiN_mask], 'k--', 'LineWidth', 1.2);

    % 2) Example species slice (e.g., k = 6) — no transpose needed
k = 6;
figure; imagesc(Zz, Rr, uT_all(:,:,k)); axis xy; set(gca,'FontSize',12);
xlabel('Z [m]'); ylabel('R [m]'); colorbar;
title(sprintf('Masked U_t (species %d) on (R,Z) from NetCDF', k));

hold on;
plot(g.lim(2,:), g.lim(1,:), 'k', 'LineWidth', 1.0);
contour(Zz, Rr, psi_m, [psiN_mask psiN_mask], 'k--', 'LineWidth', 1.2);
hold off;

%%
%% ------------------- Read-back sanity plots (axes-consistent) -------------------
Rr = ncread(outnc,'x');                 % R (nR)
Zz = ncread(outnc,'z');                 % Z (nZ)
ne_m = ncread(outnc,'ne');              % [nR x nZ]
te_m = ncread(outnc,'te');              % [nR x nZ]
ni_m = ncread(outnc,'ni_all');          % [nR x nZ x ns]
uR_m = ncread(outnc,'uR_all');          % [nR x nZ x ns]
uZ_m = ncread(outnc,'uZ_all');          % [nR x nZ x ns]
uT_m = ncread(outnc,'uT_all');          % [nR x nZ x ns]

% Build psiN on the SAME grid you plot (Z horizontal, R vertical)
[ZZplot, RRplot] = meshgrid(Zz, Rr);    % ZZplot=X (Z), RRplot=Y (R)
psi_plot = reshape(calc_psiN(g, RRplot(:), ZZplot(:), 0), size(RRplot));

% Example: species k=6, plot U_t
k = 6;
figure;
imagesc(Zz, Rr, uT_m(:,:,k));           % <-- X=Z, Y=R; no transpose needed
set(gca,'YDir','normal','FontSize',12);
xlabel('Z [m]'); ylabel('R [m]');
colorbar; title(sprintf('Masked U_t (species %d) on (R,Z) from NetCDF',k));
hold on;

% IMPORTANT: limiter is stored as [R;Z]; swap to (Z,R) when plotting
plot(g.lim(2,:), g.lim(1,:), 'k', 'LineWidth', 1.0);

% LCFS/psiN mask on the same (Z,R) axes
psiN_mask = 0.86;
contour(Zz, Rr, psi_plot, [psiN_mask psiN_mask], 'k--', 'LineWidth', 1.2);
hold off;

% (Optional) make the coarse grid less blocky
axis tight; axis xy; % already normal
% shading flat; % if you switch to pcolor(RZ) style
% or:
set(gca,'Layer','top');                 % keep overlays visible
%% ------------------- Helper: dimension re-ordering -------------------
function order = order_to_nr_nz_ns(sz, nr, nz, ns)
    if numel(sz) ~= 3
        error('Expected 3D array; got size %s', mat2str(sz));
    end
    idx_r  = find(sz == nr, 1, 'first');
    idx_z  = find(sz == nz, 1, 'first');
    idx_s  = find(sz == ns, 1, 'first');
    if isempty(idx_r) || isempty(idx_z) || isempty(idx_s)
        error('Could not identify (nr,nz,ns) in size %s (nr=%d,nz=%d,ns=%d).', ...
              mat2str(sz), nr, nz, ns);
    end
    order = [idx_r, idx_z, idx_s];
end