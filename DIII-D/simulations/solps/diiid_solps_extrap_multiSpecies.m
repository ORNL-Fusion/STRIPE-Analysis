%% Initialization
close all;
clear all;
clc;

%% Load SOLPS Data
% fileSOLPS = 'interpolated_values_multiSpecies_196154.nc';
fileSOLPS = 'interpolated_values.nc';
rS = ncread(fileSOLPS, 'gridr');    % Radial grid (Nr)
zS = ncread(fileSOLPS, 'gridz');    % Vertical grid (Nz)
neS = ncread(fileSOLPS, 'ne');      % Electron density (Nr x Nz)
TeS = ncread(fileSOLPS, 'te');      % Electron temperature (Nr x Nz)
gradTiS = -ncread(fileSOLPS,'gradTi');
gradTirS =  ncread(fileSOLPS,'gradTir'); %#ok<NASGU>
vrS = ncread(fileSOLPS,'vr');
vtS = ncread(fileSOLPS,'vt');
vzS = ncread(fileSOLPS,'vz');

%% Data Cleaning
n_min = 1e10;  % Minimum density floor
T_min = 10;    % Minimum temperature floor

neS(neS <= 0 | ~isfinite(neS)) = NaN;
TeS(TeS <= 0 | ~isfinite(TeS)) = NaN;
vrS(~isfinite(vrS)) = NaN;
vzS(~isfinite(vzS)) = NaN;
vtS(~isfinite(vtS)) = NaN;
gradTiS(~isfinite(gradTiS)) = NaN;

TeS(TeS < T_min) = T_min;

% Read EFIT data
read_efit_data;

%%
figure; imagesc(rS, zS, neS');  set(gca, 'YDir', 'normal'); colorbar; set(gca, 'ColorScale', 'linear');
title('ne (input)'); hold on; plot(g.lim(1,:), g.lim(2,:), 'r');

figure; imagesc(rS, zS, vzS');  set(gca, 'YDir', 'normal'); colorbar; set(gca, 'ColorScale', 'linear');
title('vz (input)'); hold on; plot(g.lim(1,:), g.lim(2,:), 'r');

figure; imagesc(rS, zS, vtS');  set(gca, 'YDir', 'normal'); colorbar; set(gca, 'ColorScale', 'linear');
title('vt (input)'); hold on; plot(g.lim(1,:), g.lim(2,:), 'r');

figure; imagesc(rS, zS, gradTiS');  set(gca, 'YDir', 'normal'); colorbar; set(gca, 'ColorScale', 'linear');
title('gradTi (input)'); hold on; plot(g.lim(1,:), g.lim(2,:), 'r');

%% 2) Compute psi_N on SOLPS grid (fixed)


%% Grid Preparation
[ZZ, RR] = meshgrid(zS, rS);  % Create a grid matching SOLPS data dimensions

% --- Find Outer SOLPS Boundary for Density ---
plasma_mask = neS > n_min;  % Valid plasma region for density
max_r_plasma = max(RR(plasma_mask));  % maximum r where plasma exists
% (Other boundaries are available if needed)

% --- Determine Edge Threshold for Density ---
outer_edge_mask = plasma_mask & (RR == max_r_plasma);
edge_ne_values = neS(outer_edge_mask);
threshold = min(edge_ne_values);

% --- Limiter Boundary Mask ---
inside_limiter_mask = inpolygon(RR, ZZ, g.lim(1,:), g.lim(2,:));

% --- Extrapolation Mask for Density ---
extrapolation_mask = (neS <= threshold) & inside_limiter_mask;

% --- Exponential Decay Parameters ---
% decay_length = 0.2;  % Adjust decay rate as needed

[RR, ZZ] = meshgrid(rS, zS);        % RR,ZZ are [Nz×Nr]
R_vec = RR(:).';   % 1×(Nz*Nr)
Z_vec = ZZ(:).';
[psiN_vec, ~] = calc_psiN(g, R_vec, Z_vec, 0);
psiN = reshape(psiN_vec, size(RR));   % [Nz×Nr]  %#ok<NASGU>

%% Grid Preparation
rkron = kron(rS, ones(size(zS))')';
zkron = kron(zS, ones(size(rS))');

% Euclidean target grid (keep your resolution)
num_points = 20;
[X, Y] = meshgrid(linspace(min(g.lim(1,:)), max(g.lim(1,:)), num_points), ...
                   linspace(min(g.lim(2,:)), max(g.lim(2,:)), num_points));

% Identify valid data points for knnsearch (density mask)
idx = find(~isnan(neS'));
rcoords = rkron(idx);
zcoords = zkron(idx);
coords  = [rcoords, zcoords];  % N x 2 for knnsearch

okValues_ne = neS'; okValues_ne = okValues_ne(idx);
okValues_Te = TeS'; okValues_Te = okValues_Te(idx);

%% Midplane Density & Temperature Extrapolation (Polynomial Fit)
mpfx = linspace(2.14547, 2.2051, 1000); % DIII-D 196154
mpfy = mpfx * 0;

fitDensityAtMidplane = interp2(rS, zS, neS', mpfx, mpfy, 'linear', NaN);
fitTeAtMidplane      = interp2(rS, zS, TeS', mpfx, mpfy, 'linear', NaN);

% Safe-guard the fits
gn = isfinite(fitDensityAtMidplane) & (fitDensityAtMidplane > 0);
gt = isfinite(fitTeAtMidplane)      & (fitTeAtMidplane      > 0);
if nnz(gn) >= 3
    p_ne = polyfit(mpfx(gn), log(fitDensityAtMidplane(gn)), 1);
else
    p_ne = [0, log(nanmean(fitDensityAtMidplane(gn)))];
end
if nnz(gt) >= 3
    p_Te = polyfit(mpfx(gt), log(fitTeAtMidplane(gt)), 1);
else
    p_Te = [0, log(max(nanmean(fitTeAtMidplane(gt)), T_min))];
end

% Full radial range
mpx = linspace(2.14547, 2.5, 1000);
densityAtMidplane = interp2(rS, zS, neS', mpx, 0*mpx, 'linear', NaN); densityAtMidplane(isnan(densityAtMidplane)) = 0;
TeAtMidplane      = interp2(rS, zS, TeS', mpx, 0*mpx, 'linear', NaN);  TeAtMidplane(isnan(TeAtMidplane)) = T_min;

interpfn = (mpx - 2.14547) / (2.2051 - 2.14547);
interpfn = min(max(interpfn, 0), 1);

extrapolatedne1d = interpfn .* exp(p_ne(2) + mpx * p_ne(1)) + (1 - interpfn) .* densityAtMidplane;
extrapolatedTe1d = max(interpfn .* exp(p_Te(2) + mpx * p_Te(1)) + (1 - interpfn) .* TeAtMidplane, T_min);

% Plot midplane
figure;
semilogy(mpx, densityAtMidplane, 'b');  hold on;
semilogy(mpx, extrapolatedne1d, 'b.');
semilogy(mpx, TeAtMidplane, 'r');
semilogy(mpx, extrapolatedTe1d, 'r.');
legend('Density Midplane', 'Extrapolated Density', 'Temperature Midplane', 'Extrapolated Temperature');
xlabel('R [m]'); ylabel('Value'); title('Midplane Density and Temperature Extrapolation'); hold off;

% Save midplane data
writematrix([mpx; densityAtMidplane], 'nemidplane.csv');
writematrix([mpx; TeAtMidplane], 'Temidplane.csv');

%% Efficient Euclidean Extrapolation on 2D Grid (ne, Te) Using parfor
num_points_all = numel(X);

% Start/attach pool if available
try
    pool = gcp('nocreate'); if isempty(pool), pool = parpool; end
catch
    pool = []; % serial fallback
end
num_workers = max(1, ternary(~isempty(pool), pool.NumWorkers, 1));
batch_size  = ceil(num_points_all / (10 * num_workers));

extrap_coords = [X(:), Y(:)];
nBatches = ceil(num_points_all / batch_size);
val_ne_cell = cell(nBatches, 1);
val_Te_cell = cell(nBatches, 1);

parfor (b = 1:nBatches, num_workers)
    s = (b-1)*batch_size + 1; e = min(b*batch_size, num_points_all);
    subset = extrap_coords(s:e, :);
    idx_nearest = knnsearch(coords, subset, 'K', 1);
    d = vecnorm(coords(idx_nearest, :) - subset, 2, 2);
    val_ne_cell{b} = okValues_ne(idx_nearest) .* exp(p_ne(1) * d);
    val_Te_cell{b} = max(okValues_Te(idx_nearest) .* exp(p_Te(1) * d), T_min);
end

val_ne = reshape(cell2mat(val_ne_cell), size(X));
val_Te = reshape(cell2mat(val_Te_cell), size(X));

%% Visualization of Extrapolated Electron Density
figure; imagesc([min(X(:)) max(X(:))], [min(Y(:)) max(Y(:))], val_ne);
set(gca, 'YDir', 'normal'); colorbar; set(gca, 'ColorScale', 'log');
title('Extrapolated Electron Density (n_e)'); clim([10^1, 10^20]); hold on; plot(g.lim(1,:), g.lim(2,:), 'r');

%% Visualization of Extrapolated Electron Temperature
figure; imagesc([min(X(:)) max(X(:))], [min(Y(:)) max(Y(:))], val_Te);
set(gca, 'YDir', 'normal'); colorbar; set(gca, 'ColorScale', 'log');
title('Extrapolated Electron Temperature (Te)'); clim([1, max(val_Te(:))]); hold on; plot(g.lim(1,:), g.lim(2,:), 'r');

%% ---- Euclidean Extrapolation for velocities and gradTi (REVISED) ----
% Build "ok" sets for velocities/gradTi (finite points)
vr_ok_mask     = isfinite(vrS');
vz_ok_mask     = isfinite(vzS');
vt_ok_mask     = isfinite(vtS');
gradTi_ok_mask = isfinite(gradTiS');

vr_coords     = [rkron(vr_ok_mask'),    zkron(vr_ok_mask')];
vz_coords     = [rkron(vz_ok_mask'),    zkron(vz_ok_mask')];
vt_coords     = [rkron(vt_ok_mask'),    zkron(vt_ok_mask')];
gradTi_coords = [rkron(gradTi_ok_mask'), zkron(gradTi_ok_mask')];

vr_values_ok     = vrS';      vr_values_ok     = vr_values_ok(vr_ok_mask);
vz_values_ok     = vzS';      vz_values_ok     = vz_values_ok(vz_ok_mask);
vt_values_ok     = vtS';      vt_values_ok     = vt_values_ok(vt_ok_mask);
gradTi_values_ok = gradTiS';  gradTi_values_ok = gradTi_values_ok(gradTi_ok_mask);



% Use a Euclidean e-folding length for velocities/gradTi to avoid blow-up
decay_length = 0.2;  % keep your number for consistency
L_vel  = decay_length;
L_grad = decay_length;

val_vr_cell     = cell(nBatches, 1);
val_vz_cell     = cell(nBatches, 1);
val_vt_cell     = cell(nBatches, 1);
val_gradTi_cell = cell(nBatches, 1);

parfor (b = 1:nBatches, num_workers)
    s = (b-1)*batch_size + 1; e = min(b*batch_size, num_points_all);
    subset = extrap_coords(s:e, :);

    % vr
    idx_vr = knnsearch(vr_coords, subset, 'K', 1);
    d_vr   = vecnorm(vr_coords(idx_vr, :) - subset, 2, 2);
    val_vr_cell{b} = vr_values_ok(idx_vr) .* exp(-d_vr / L_vel);

    % vz
    idx_vz = knnsearch(vz_coords, subset, 'K', 1);
    d_vz   = vecnorm(vz_coords(idx_vz, :) - subset, 2, 2);
    val_vz_cell{b} = vz_values_ok(idx_vz) .* exp(-d_vz / L_vel);

    % vt
    idx_vt = knnsearch(vt_coords, subset, 'K', 1);
    d_vt   = vecnorm(vt_coords(idx_vt, :) - subset, 2, 2);
    val_vt_cell{b} = vt_values_ok(idx_vt) .* exp(-d_vt / L_vel);

    % gradTi
    idx_gt = knnsearch(gradTi_coords, subset, 'K', 1);
    d_gt   = vecnorm(gradTi_coords(idx_gt, :) - subset, 2, 2);
    val_gradTi_cell{b} = gradTi_values_ok(idx_gt) .* exp(-d_gt / L_grad);
end

val_vr     = reshape(cell2mat(val_vr_cell),     size(X));
val_vz     = reshape(cell2mat(val_vz_cell),     size(X));
val_vt     = reshape(cell2mat(val_vt_cell),     size(X));
val_gradTi = reshape(cell2mat(val_gradTi_cell), size(X));

%% 5) Visualization of 2D Extrapolated Data (Euclidean grid)
vars2D = {val_ne, val_Te, val_vr, val_vz, val_vt, val_gradTi};
varNames2D = {'Electron Density (neS)', 'Electron Temperature (TeS)', ...
              'Radial Velocity (vrS)', 'Poloidal Velocity (vzS)', ...
              'Toroidal Velocity (vtS)', 'Ion Temperature Gradient (gradTiS)'};

for i = 1:length(vars2D)
    figure;
    imagesc([min(X(:)) max(X(:))], [min(Y(:)) max(Y(:))], vars2D{i});
    set(gca, 'YDir', 'normal'); colorbar; set(gca, 'ColorScale', 'linear');
    title(['Extrapolated ', varNames2D{i}]);
    hold on; plot(g.lim(1,:), g.lim(2,:), 'r');
end


%{
% %% === Scheme-1: Multi-species (e.g., Carbon) extrapolation: SAME method as deuterium ===
% % We will:
% %  1) read species fields
% %  2) for each species, make the SAME midplane fit (exp in R) for n_i and T_i
% %  3) do the SAME Euclidean 2-D extrapolation using knnsearch for n_i, T_i
% %  4) apply the SAME exponential decay beyond the edge for velocities (uR,uZ,uT)
% % Notes:
% %  - uses your existing rS,zS, g, plasma_mask, inside_limiter_mask, decay_length
% %  - keeps deuterium fields untouched
% %  - handles files where the species dimension is first
% 
% % ---------- 0) Read species arrays ----------
% Z_all   = ncread(fileSOLPS,'atomic_number');  % [ns]
% q_all   = ncread(fileSOLPS,'charge_number');  %#ok<NASGU>
% ni_all0 = ncread(fileSOLPS,'ni_all');         % expect [ns x nr x nz] from your file
% ti_all0 = ncread(fileSOLPS,'ti_all');
% uR_all0 = ncread(fileSOLPS,'uR_all');
% uZ_all0 = ncread(fileSOLPS,'uZ_all');
% uT_all0 = ncread(fileSOLPS,'uT_all');
% 
% ns = numel(Z_all);
% nr = numel(rS);
% nz = numel(zS);
% 
% % ---------- 1) Normalize to [nr x nz x ns] ----------
% perm_to_nr_nz_ns = @(A) permute(A, order_to_nr_nz_ns(size(A), nr, nz, ns));
% ni_all = perm_to_nr_nz_ns(ni_all0);
% ti_all = perm_to_nr_nz_ns(ti_all0);
% uR_all = perm_to_nr_nz_ns(uR_all0);
% uZ_all = perm_to_nr_nz_ns(uZ_all0);
% uT_all = perm_to_nr_nz_ns(uT_all0);
% 
% % Precompute coords for knnsearch once (like you did for deuterium)
% [ZZ_s, RR_s] = meshgrid(zS, rS);     % RR_s,ZZ_s are [nr x nz]
% rkron_s = RR_s(:);                   % length = nr*nz
% zkron_s = ZZ_s(:);
% coords_s = [rkron_s, zkron_s];       % N x 2
% 
% % ---------- 2) Loop species with SAME midplane & Euclidean extrap ----------
% % T_min = evalin('base','T_min');      % reuse your thresholds
% % n_min = evalin('base','n_min');
% 
% % Midplane line you used for D
% % mpfx = linspace(2.14333, 2.20675, 1000);  % DIII-D 196154
% mpfx = linspace(2.28122, 2.3156, 1000); % DIII-D # 200882
% 
% 
% mpfy = mpfx*0;
% % mpx  = linspace(2.14333, 2.5, 1000);
% mpx = linspace(2.28122, 2.5, 1000); % DIII-D # 200882
% 
% % For Euclidean 2-D extrap grid (reuse your X,Y)
% % If X,Y not in scope yet, comment next two lines and use your existing X,Y.
% % num_points = 1000; [X, Y] = meshgrid(linspace(min(g.lim(1,:)), max(g.lim(1,:)), num_points), ...
% %                                      linspace(min(g.lim(2,:)), max(g.lim(2,:)), num_points));
% 
% num_points_all = numel(X);
% pool = gcp(); num_workers = pool.NumWorkers;
% batch_size = ceil(num_points_all / (10 * num_workers));
% extrap_coords = [X(:), Y(:)];
% nBatches = ceil(num_points_all / batch_size);
% 
% % Edge masks from your bulk ne (SAME edge used for all species)
% % plasma_mask : [nr x nz], inside_limiter_mask : [nr x nz], decay_length set above
% 
% % Output holders (optional: if you want Euclidean fields on X,Y)
% val_ni_mspecies = zeros(size(X,1), size(X,2), ns);
% val_Ti_mspecies = zeros(size(X,1), size(X,2), ns);
% 
% for k = 1:ns
%     % ----- 2a) Midplane fit for species k (same as you did for ne/Te) -----
%     ni_k = ni_all(:,:,k);          % [nr x nz]
%     Ti_k = ti_all(:,:,k);
% 
%     % Along midplane (Z=0)
%     ni_mid = interp2(rS, zS, ni_k', mpfx, mpfy, 'linear', NaN);
%     Ti_mid = interp2(rS, zS, Ti_k', mpfx, mpfy, 'linear', NaN);
% 
%     % Safe-guard: if too few finite points, skip polyfit to avoid warnings
%     good_n = isfinite(ni_mid) & (ni_mid > 0);
%     good_T = isfinite(Ti_mid) & (Ti_mid > 0);
% 
%     if nnz(good_n) >= 3
%         p_ni = polyfit(mpfx(good_n), log(ni_mid(good_n)), 1);
%     else
%         % fallback: copy deuterium p_ne if available, else flat
%         p_ni = p_ne;  % re-use bulk electron fit slope as proxy
%     end
% 
%     if nnz(good_T) >= 3
%         p_Ti = polyfit(mpfx(good_T), log(Ti_mid(good_T)), 1);
%     else
%         % fallback: copy deuterium p_Te if available, else clamp to T_min
%         p_Ti = p_Te;
%     end
% 
%     % Baseline along full mpx from inside profile
%     ni_base = interp2(rS, zS, ni_k', mpx, 0*mpx, 'linear', NaN); ni_base(~isfinite(ni_base)) = 0;
%     Ti_base = interp2(rS, zS, Ti_k', mpx, 0*mpx, 'linear', NaN); Ti_base(~isfinite(Ti_base)) = T_min;
% 
%     interpfn = (mpx - 2.281) / (2.3156 - 2.28122); interpfn = min(max(interpfn,0),1);
% 
%     ni_1d = interpfn .* exp(p_ni(2) + mpx*p_ni(1)) + (1-interpfn) .* ni_base;
%     Ti_1d = max(interpfn .* exp(p_Ti(2) + mpx*p_Ti(1)) + (1-interpfn) .* Ti_base, T_min);
% 
%     % --- Plot Midplane Density and Temperature for Carbon (Z=6) ---
% isC = (Z_all == 6);
% if any(isC)
%     kc = find(isC,1,'first');
% 
%     % Midplane R array you defined earlier
%     Rmid = mpx;
% 
%     % Get carbon midplane data (already interpolated earlier in your loop)
%     ni_mid_C = interp2(rS, zS, ni_all(:,:,kc)', mpfx, mpfy, 'linear', NaN);
%     Ti_mid_C = interp2(rS, zS, ti_all(:,:,kc)', mpfx, mpfy, 'linear', NaN);
% 
%     % Extrapolated 1D profiles you constructed
%     ni_1d_C = ni_1d;   % from your loop, when k == kc
%     Ti_1d_C = Ti_1d;
% 
%     % Plot
%     figure;
%     semilogy(Rmid, ni_mid_C, 'b'); hold on;
%     semilogy(Rmid, ni_1d_C, 'b.');  % extrapolated density
%     semilogy(Rmid, Ti_mid_C, 'r');
%     semilogy(Rmid, Ti_1d_C, 'r.');
% 
%     legend('Carbon Density Midplane', 'Extrapolated Carbon Density', ...
%            'Carbon Temperature Midplane', 'Extrapolated Carbon Temperature');
%     xlabel('R [m]');
%     ylabel('Value');
%     title('Carbon Midplane Density and Temperature Extrapolation');
%     hold off;
% end
% 
%     % ----- 2b) Euclidean 2D extrap (same as for deuterium) for n_i and T_i -----
%     % Build okValues using in-plasma points (finite & positive)
%     ok_mask = isfinite(ni_k) & ni_k > 0;
%     ok_idx = find(ok_mask(:));
%     if isempty(ok_idx)
%         warning('Species k=%d has no positive finite ni points; skipping Euclidean extrap.', k);
%         continue;
%     end
%     coords_ok = coords_s(ok_idx, :);
%     ok_ni = ni_k(ok_idx);
%     ok_Ti = Ti_k(ok_idx);
% 
%     % parfor batching identical to your D block
%     ni_cell = cell(nBatches,1);
%     Ti_cell = cell(nBatches,1);
%     parfor b = 1:nBatches
%         s = (b-1)*batch_size + 1; e = min(b*batch_size, num_points_all);
%         subset = extrap_coords(s:e, :);
%         idx_nearest = knnsearch(coords_ok, subset, 'K', 1);
%         d = vecnorm(coords_ok(idx_nearest, :) - subset, 2, 2);
% 
%         % Use species-specific slopes p_ni(1), p_Ti(1) exactly as for D
%         ni_cell{b} = ok_ni(idx_nearest) .* exp(-p_ni(1) * d);
%         Ti_cell{b} = max(ok_Ti(idx_nearest) .* exp(-p_Ti(1) * d), T_min);
%     end
% 
%     val_ni = reshape(cell2mat(ni_cell), size(X));
%     val_Ti = reshape(cell2mat(Ti_cell), size(X));
% 
%     % Store
%     val_ni_mspecies(:,:,k) = val_ni;
%     val_Ti_mspecies(:,:,k) = val_Ti;
% 
%     % ----- 2c) Edge-following exponential decay for velocities (same as your D) -----
%     uR_k = uR_all(:,:,k); uZ_k = uZ_all(:,:,k); uT_k = uT_all(:,:,k);
%     for jz = 1:nz
%         col = plasma_mask(:, jz);
%         edge_idx = find(col, 1, 'last');
%         if ~isempty(edge_idx) && edge_idx < nr
%             R_edge   = rS(edge_idx);
%             uR_edge  = uR_k(edge_idx, jz);
%             uZ_edge  = uZ_k(edge_idx, jz);
%             uT_edge  = uT_k(edge_idx, jz);
%             for ir = edge_idx+1:nr
%                 if inside_limiter_mask(ir, jz)
%                     fac = exp(-(rS(ir) - R_edge)/decay_length); % SAME as your velocity rule
%                     uR_k(ir, jz) = uR_edge * fac;
%                     uZ_k(ir, jz) = uZ_edge * fac;
%                     uT_k(ir, jz) = uT_edge * fac;
%                 end
%             end
%         end
%     end
%     uR_all(:,:,k) = uR_k; uZ_all(:,:,k) = uZ_k; uT_all(:,:,k) = uT_k;
% end
% 
% % Optional: quick carbon check (Z=6)
% isC = (Z_all == 6);
% if any(isC)
%     kc = find(isC,1,'first');
%     figure; imagesc(rS, zS, (val_Ti_mspecies(:,:,6))'); set(gca,'YDir','normal'); colorbar;
%     title(sprintf('n_i (post-extrap), Z=6 species k=%d', kc)); hold on; plot(g.lim(1,:), g.lim(2,:), 'r');
%     figure; imagesc(rS, zS, (uR_all(:,:,kc))'); set(gca,'YDir','normal'); colorbar;
%     title(sprintf('u_R (post-extrap), Z=6 species k=%d', kc)); hold on; plot(g.lim(1,:), g.lim(2,:), 'r');
% end
% 
% % ---------- helper: figure out permutation to [nr x nz x ns] ----------
% function order = order_to_nr_nz_ns(sz, nr, nz, ns)
%     if numel(sz) ~= 3
%         error('Expected 3D array for species fields; got size %s', mat2str(sz));
%     end
%     idx_r  = find(sz == nr, 1, 'first');
%     idx_z  = find(sz == nz, 1, 'first');
%     idx_ns = find(sz == ns, 1, 'first');
%     if isempty(idx_r) || isempty(idx_z) || isempty(idx_ns)
%         error('Could not identify (nr,nz,ns) in size %s (nr=%d, nz=%d, ns=%d).', ...
%               mat2str(sz), nr, nz, ns);
%     end
%     order = [idx_r, idx_z, idx_ns];
% end
% 

%}
%% === Scheme 2: Multi-species (e.g., Carbon) extrapolation: SAME method as deuterium ===
% We will:
%  1) read species fields
%  2) for each species, make the SAME midplane fit (exp in R) for n_i and T_i
%  3) (REVISED) 2-D profiles also use Option B: edge clamp + exponential tail on native grid
%  4) (REVISED) velocities (uR,uZ,uT) also use Option B on native grid, then map to (X,Y)
% Notes:
%  - uses your existing rS,zS, g, plasma_mask, inside_limiter_mask, decay_length
%  - keeps deuterium fields untouched
%  - handles files where the species dimension is first

% ---------- 0) Read species arrays ----------
Z_all   = ncread(fileSOLPS,'atomic_number');  % [ns]
q_all   = ncread(fileSOLPS,'charge_number');  %#ok<NASGU>
ni_all0 = ncread(fileSOLPS,'ni_all');         % expect 3D from your file
ti_all0 = ncread(fileSOLPS,'ti_all');
uR_all0 = ncread(fileSOLPS,'uR_all');
uZ_all0 = ncread(fileSOLPS,'uZ_all');
uT_all0 = ncread(fileSOLPS,'uT_all');

ns = numel(Z_all);
nr = numel(rS);
nz = numel(zS);

% ---------- 1) Normalize to [nr x nz x ns] ----------
perm_to_nr_nz_ns = @(A) permute(A, order_to_nr_nz_ns(size(A), nr, nz, ns));
ni_all = perm_to_nr_nz_ns(ni_all0);
ti_all = perm_to_nr_nz_ns(ti_all0);
uR_all = perm_to_nr_nz_ns(uR_all0);
uZ_all = perm_to_nr_nz_ns(uZ_all0);
uT_all = perm_to_nr_nz_ns(uT_all0);

% (coords precompute not needed for Option B, but harmless if present above)

% ---------- 2) Loop species with SAME midplane & 2-D Option B ----------
% (We reuse your global p_ne, p_Te, T_min, n_min, plasma_mask, inside_limiter_mask)
% Use your DIII-D shot midplane lines:
mpfx = linspace(2.14547, 2.2051, 1000); % DIII-D 196154 % for 196154
% mpfx = linspace(2.28122, 2.3156, 1000);     % for 200882
mpfy = mpfx*0;

mpx  = linspace(2.14547, 2.5, 1000);      % for 196154
% mpx = linspace(2.28122, 2.5, 1000);         % for 200882

% Option B controls (keep same style/scale as density/Ti)
K      = 1;          % robust edge from last K inside samples
Ldecay = 0.1;       % [m] decay length for tails (adjust if needed)

% Output holders on (X,Y)
val_ni_mspecies = zeros(size(X,1), size(X,2), ns);
val_Ti_mspecies = zeros(size(X,1), size(X,2), ns);
val_uR_mspecies = zeros(size(X,1), size(X,2), ns);
val_uZ_mspecies = zeros(size(X,1), size(X,2), ns);
val_uT_mspecies = zeros(size(X,1), size(X,2), ns);

for k = 1:ns
    % ----- 2a) Midplane fit for species k (same as you did for ne/Te) -----
    ni_k = ni_all(:,:,k);          % [nr x nz]
    Ti_k = ti_all(:,:,k);

    % Along midplane (Z=0)
    ni_mid = interp2(rS, zS, ni_k', mpfx, mpfy, 'linear', NaN);
    Ti_mid = interp2(rS, zS, Ti_k', mpfx, mpfy, 'linear', NaN);

    % Safe-guard: if too few finite points, skip polyfit to avoid warnings
    good_n = isfinite(ni_mid) & (ni_mid > 0);
    good_T = isfinite(Ti_mid) & (Ti_mid > 0);

    if nnz(good_n) >= 3
        p_ni = polyfit(mpfx(good_n), log(ni_mid(good_n)), 1);
    else
        p_ni = p_ne;   % reuse your bulk slope
    end

    if nnz(good_T) >= 3
        p_Ti = polyfit(mpfx(good_T), log(Ti_mid(good_T)), 1);
    else
        p_Ti = p_Te;
    end

    % Baseline along full mpx from inside profile
    ni_base = interp2(rS, zS, ni_k', mpx, 0*mpx, 'linear', NaN); ni_base(~isfinite(ni_base)) = 0;
    Ti_base = interp2(rS, zS, Ti_k', mpx, 0*mpx, 'linear', NaN); Ti_base(~isfinite(Ti_base)) = T_min;

    % interpfn = (mpx - 2.281) / (2.3156 - 2.28122); interpfn = min(max(interpfn,0),1); % 200882
    interpfn = (mpx - 2.14547) / (2.2051 - 2.14547); interpfn = min(max(interpfn,0),1); % 196154
    ni_1d = interpfn .* exp(p_ni(2) + mpx*p_ni(1)) + (1-interpfn) .* ni_base;
    Ti_1d = max(interpfn .* exp(p_Ti(2) + mpx*p_Ti(1)) + (1-interpfn) .* Ti_base, T_min);

    % Optional: quick carbon check (Z=6)
isC = (Z_all == 6);
if any(isC)
    kc = find(isC,1,'first');
    figure; imagesc(rS, zS, (ni_all(:,:,kc)')); set(gca,'YDir','normal'); colorbar;
    title(sprintf('n_i (post-extrap), Z=6 species k=%d', kc)); hold on; plot(g.lim(1,:), g.lim(2,:), 'r');
    figure; imagesc(rS, zS, (uR_all(:,:,kc)')); set(gca,'YDir','normal'); colorbar;
    title(sprintf('u_R (post-extrap to XY), Z=6 species k=%d', kc)); hold on; plot(g.lim(1,:), g.lim(2,:), 'r');
end

    % ===== Option B (midplane): edge clamp + exponential tail =====
    R2 = mpfx(end);
    % density midplane tail
    idx_in = find(isfinite(ni_1d) & (mpx <= R2));
    if ~isempty(idx_in)
        idx_tail = idx_in(max(numel(idx_in)-K+1,1):end);
        edge_val_ni = max(ni_1d(idx_tail));
        beyond = (mpx > R2);
        ni_1d(beyond) = edge_val_ni .* exp(-(mpx(beyond)-R2)/Ldecay);
    end
    % temperature midplane tail
    idx_inT = find(isfinite(Ti_1d) & (mpx <= R2));
    if ~isempty(idx_inT)
        idx_tailT = idx_inT(max(numel(idx_inT)-K+1,1):end);
        edge_val_Ti = max(Ti_1d(idx_tailT));
        beyondT = (mpx > R2);
        Ti_1d(beyondT) = max(edge_val_Ti .* exp(-(mpx(beyondT)-R2)/Ldecay), T_min);
    end
    % ===== end Option B (midplane) =====

    % (Optional quick carbon plot)
    if any(Z_all == 6) && k == find(Z_all==6,1,'first')
        Rmid = mpx;
        ni_mid_C = interp2(rS, zS, ni_k', mpfx, mpfy, 'linear', NaN);
        Ti_mid_C = interp2(rS, zS, Ti_k', mpfx, mpfy, 'linear', NaN);
        figure; semilogy(Rmid, ni_mid_C, 'b'); hold on;
        semilogy(Rmid, ni_1d, 'b.');
        semilogy(Rmid, Ti_mid_C, 'r');
        semilogy(Rmid, Ti_1d, 'r.');
        legend('Carbon Density Midplane','Extrapolated Carbon Density', ...
               'Carbon Temperature Midplane','Extrapolated Carbon Temperature');
        xlabel('R [m]'); ylabel('Value'); title('Carbon Midplane Extrapolation'); hold off;
    end

    % ----- 2b) 2-D Option B on native (rS,zS) for n_i and T_i -----
    ni_k_ex = ni_k;
    Ti_k_ex = Ti_k;

    for jz = 1:nz
        col = plasma_mask(:, jz);                 % in-plasma mask for this z
        edge_idx = find(col, 1, 'last');          % last valid r index

        if ~isempty(edge_idx)
            % robust edge values from last K points inside the edge
            i0 = max(1, edge_idx - K + 1);
            edge_vals_ni = ni_k(i0:edge_idx, jz);
            edge_vals_Ti = Ti_k(i0:edge_idx, jz);

            edge_val_ni_col = max(edge_vals_ni, [], 'omitnan');
            edge_val_Ti_col = max(edge_vals_Ti, [], 'omitnan');

            if ~isfinite(edge_val_ni_col), edge_val_ni_col = ni_k(edge_idx, jz); end
            if ~isfinite(edge_val_Ti_col), edge_val_Ti_col = max(Ti_k(edge_idx, jz), T_min); end

            % apply exponential tail beyond the edge, only inside limiter
            for ir = edge_idx+1:nr
                if inside_limiter_mask(ir, jz)
                    dr = rS(ir) - rS(edge_idx);
                    ni_k_ex(ir, jz) = max(edge_val_ni_col * exp(-dr / Ldecay), 0);
                    Ti_k_ex(ir, jz) = max(edge_val_Ti_col * exp(-dr / Ldecay), T_min);
                end
            end
        end
    end

    % Map Option-B-extrapolated n_i, T_i to your (X,Y) grid
    val_ni = interp2(rS, zS, ni_k_ex', X, Y, 'linear', 0);
    val_Ti = interp2(rS, zS, Ti_k_ex', X, Y, 'linear', T_min);

    % Store
    val_ni_mspecies(:,:,k) = val_ni;
    val_Ti_mspecies(:,:,k) = val_Ti;

    % ----- 2c) (REVISED) Velocities: Option B on native (rS,zS) then map to (X,Y) -----
    uR_k = uR_all(:,:,k);
    uZ_k = uZ_all(:,:,k);
    uT_k = uT_all(:,:,k);

    uR_k_ex = uR_k;
    uZ_k_ex = uZ_k;
    uT_k_ex = uT_k;

    for jz = 1:nz
        col = plasma_mask(:, jz);
        edge_idx = find(col, 1, 'last');          % last valid r index inside plasma

        if ~isempty(edge_idx)
            % robust edge values (mean over last K points) reduces noise
            i0 = max(1, edge_idx - K + 1);
            uR_edge = mean(uR_k(i0:edge_idx, jz), 'omitnan');
            uZ_edge = mean(uZ_k(i0:edge_idx, jz), 'omitnan');
            uT_edge = mean(uT_k(i0:edge_idx, jz), 'omitnan');

            if ~isfinite(uR_edge), uR_edge = uR_k(edge_idx, jz); end
            if ~isfinite(uZ_edge), uZ_edge = uZ_k(edge_idx, jz); end
            if ~isfinite(uT_edge), uT_edge = uT_k(edge_idx, jz); end

            % exponential tail for each component beyond R_edge, inside limiter
            for ir = edge_idx+1:nr
                if inside_limiter_mask(ir, jz)
                    dr = rS(ir) - rS(edge_idx);
                    fac = exp(-dr / Ldecay);
                    uR_k_ex(ir, jz) = uR_edge * fac;
                    uZ_k_ex(ir, jz) = uZ_edge * fac;
                    uT_k_ex(ir, jz) = uT_edge * fac;
                end
            end
        end
    end

    % Map velocities to (X,Y)
    val_uR = interp2(rS, zS, uR_k_ex', X, Y, 'linear', 0);
    val_uZ = interp2(rS, zS, uZ_k_ex', X, Y, 'linear', 0);
    val_uT = interp2(rS, zS, uT_k_ex', X, Y, 'linear', 0);

    % Store
    val_uR_mspecies(:,:,k) = val_uR;
    val_uZ_mspecies(:,:,k) = val_uZ;
    val_uT_mspecies(:,:,k) = val_uT;

    % keep updated copies if you later want to write them back
    uR_all(:,:,k) = uR_k_ex;
    uZ_all(:,:,k) = uZ_k_ex;
    uT_all(:,:,k) = uT_k_ex;
end
%%
% Optional: quick carbon check (Z=6)
isC = (Z_all == 6);
if any(isC)
    kc = find(isC,1,'first');
    figure; imagesc(rS, zS, (val_ni_mspecies(:,:,kc))); set(gca,'YDir','normal'); colorbar;
    title(sprintf('n_i (post-extrap), Z=6 species k=%d', kc)); hold on; plot(g.lim(1,:), g.lim(2,:), 'r');
    figure; imagesc(rS, zS, (val_uR_mspecies(:,:,kc))); set(gca,'YDir','normal'); colorbar;
    title(sprintf('u_R (post-extrap to XY), Z=6 species k=%d', kc)); hold on; plot(g.lim(1,:), g.lim(2,:), 'r');
end


%% Save the data

% save('extrapolated_data_200882_wallFix.mat'); % DIII-D # 200882

save('extrapolated_data_multiSpecies_196154.mat'); % DIII-D # 196154

%% ---------- helper: figure out permutation to [nr x nz x ns] ----------
function order = order_to_nr_nz_ns(sz, nr, nz, ns)
    if numel(sz) ~= 3
        error('Expected 3D array for species fields; got size %s', mat2str(sz));
    end
    idx_r  = find(sz == nr, 1, 'first');
    idx_z  = find(sz == nz, 1, 'first');
    idx_ns = find(sz == ns, 1, 'first');
    if isempty(idx_r) || isempty(idx_z) || isempty(idx_ns)
        error('Could not identify (nr,nz,ns) in size %s (nr=%d, nz=%d, ns=%d).', ...
              mat2str(sz), nr, nz, ns);
    end
    order = [idx_r, idx_z, idx_ns];
end

%% ---- helpers ----
function y = ternary(cond, a, b)
    if cond, y = a; else, y = b; end
end