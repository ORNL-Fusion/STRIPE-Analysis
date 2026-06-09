% close all;
% clear all;

%% Constants (kept from your snippet)
ME   = 9.10938356e-31;
MI   = 1.6737236e-27;
Q    = 1.60217662e-19;
EPS0 = 8.854187e-12;
amu  = 18;

%% Files
ncfile  = 'profilesDIIID_196154_multi.nc';  % multispecies on (R,z,ns) from our writer
% ncfile1 = 'profilesDIIID_196154.nc';               % single-fluid for baseline plots

%% --- Single-fluid read (unchanged flow) ---
R = ncread(ncfile,'x');     % grid R
z = ncread(ncfile,'z');     % grid z

ne = ncread(ncfile,'ne');
te = ncread(ncfile,'te');
vr = ncread(ncfile,'vr');
vt = ncread(ncfile,'vt');
vz = ncread(ncfile,'vz');
if any(strcmp({ncinfo(ncfile).Variables.Name},'ni')), ni = ncread(ncfile,'ni'); else, ni = []; end
if any(strcmp({ncinfo(ncfile).Variables.Name},'ti')), ti = ncread(ncfile,'ti'); else, ti = []; end

%% ---- NEW: multispecies metadata & fields ----
Z_all = ncread(ncfile,'atomic_number');   % [ns]
q_all = ncread(ncfile,'charge_number');   % [ns]

% 3D arrays as written: [nR x nZ x ns]; but we add a safeguard
ni_all_raw = ncread(ncfile,'ni_all');     % size unknown orientation
ti_all_raw = ncread(ncfile,'ti_all');
uR_all_raw = ncread(ncfile,'uR_all');
uZ_all_raw = ncread(ncfile,'uZ_all');
uT_all_raw = ncread(ncfile,'uT_all');

% ---- Surface geometry (assumed provided upstream) ----
% Variables used: centroid [N x 3], norm_vec [N x 3], X,Y,Z for patch, g.lim for outline
% Make radial coordinate of centroids:
r_centroid = sqrt(centroid(:,1).^2 + centroid(:,2).^2);

%% Ensure ascending axes and reindex fields
Rvec = R(:); zvec = z(:);
[Rvec, iR] = sort(Rvec, 'ascend');
[zvec, iz] = sort(zvec, 'ascend');

% Reindex single-fluid onto sorted axes
ne = ne(iR, iz); te = te(iR, iz);
vr = vr(iR, iz); vt = vt(iR, iz); vz = vz(iR, iz);
if ~isempty(ni), ni = ni(iR, iz); end
if ~isempty(ti), ti = ti(iR, iz); end

% Normalize multispecies dims to [nR x nZ x ns], then sort the first two dims
to_RZNS = @(A) normalize_species_dims(A, numel(R), numel(z), numel(Z_all));  % helper at bottom
ni_all = to_RZNS(ni_all_raw); ni_all = ni_all(iR, iz, :);
ti_all = to_RZNS(ti_all_raw); ti_all = ti_all(iR, iz, :);
uR_all = to_RZNS(uR_all_raw); uR_all = uR_all(iR, iz, :);
uZ_all = to_RZNS(uZ_all_raw); uZ_all = uZ_all(iR, iz, :);
uT_all = to_RZNS(uT_all_raw); uT_all = uT_all(iR, iz, :);

% If single-fluid arrays were stored as (z,R) originally, swap once
if size(ne,1)==numel(zvec) && size(ne,2)==numel(Rvec)
    ne = ne.'; te = te.'; vr = vr.'; vt = vt.'; vz = vz.';
    if ~isempty(ni), ni = ni.'; end
    if ~isempty(ti), ti = ti.'; end
end

%% -----------------------------------------------
%% Select CARBON (Z=6) — per-species (sum over all charge states) by default
Z_target = 6;
idxC_all = find(Z_all == Z_target);
if isempty(idxC_all)
    error('No species with atomic_number == %d found.', Z_target);
end

% --- Choose mode: sum all carbon charge states (per-species), OR pick one q ---
use_single_charge = true;   % <- set true to pick one charge state instead of summing
q_sel = 6;                   % if use_single_charge=true, e.g. C4+

if ~use_single_charge
    % Per-species CARBON: sum over all charge states where Z==6
    niC = sum(ni_all(:,:,idxC_all), 3, 'omitnan');              % carbon density (all charges)
    % Density-weighted T and velocities:
    tiny = 1e-60;
    wC   = ni_all(:,:,idxC_all);                                % [nR x nZ x nC]
    denom = sum(wC, 3, 'omitnan');
    denom(denom<tiny) = tiny;
    tiC = sum(wC .* ti_all(:,:,idxC_all), 3, 'omitnan') ./ denom;
    uR_C = sum(wC .* uR_all(:,:,idxC_all), 3, 'omitnan') ./ denom;
    uZ_C = sum(wC .* uZ_all(:,:,idxC_all), 3, 'omitnan') ./ denom;
    uT_C = sum(wC .* uT_all(:,:,idxC_all), 3, 'omitnan') ./ denom;

    labelC = sprintf('Carbon (Z=%d) total (sum over charges)', Z_target);
else
    % Per-charge: pick a single charge state of carbon
    idx_q = find((Z_all==Z_target) & (q_all==q_sel), 1, 'first');
    if isempty(idx_q)
        error('Requested carbon charge state q=%d not found among Z=%d species.', q_sel, Z_target);
    end
    % niC = ni_all(:,:,idx_q);
    niC = ni_all(:,:,2); % Deuterium
    tiC = ti_all(:,:,2); % TiC=TiD
    % uR_C = uR_all(:,:,idx_q);
    % uZ_C = uZ_all(:,:,idx_q);
    % uT_C = uT_all(:,:,idx_q);

    uR_C = uR_all(:,:,2);
    uZ_C = uZ_all(:,:,2);
    uT_C = uT_all(:,:,2);

    labelC = sprintf('Carbon (Z=%d, q=%d)', Z_target, q_sel);
end

%% -----------------------------------------------
%% Interpolation settings (Euclidean R–z, same as you use elsewhere)
interp_mode = 'linear';   % try 'nearest' to preserve sharp peaks
extrap_mode = 'none';     % or 'extrap' if you want off-grid values

% Build interpolants for single-fluid (kept)
F_ne = griddedInterpolant({Rvec,zvec}, ne, interp_mode, extrap_mode);
F_te = griddedInterpolant({Rvec,zvec}, te, interp_mode, extrap_mode);
F_vr = griddedInterpolant({Rvec,zvec}, vr, interp_mode, extrap_mode);
F_vt = griddedInterpolant({Rvec,zvec}, vt, interp_mode, extrap_mode);
F_vz = griddedInterpolant({Rvec,zvec}, vz, interp_mode, extrap_mode);
if ~isempty(ni), F_ni = griddedInterpolant({Rvec,zvec}, ni, interp_mode, extrap_mode); end
if ~isempty(ti), F_ti = griddedInterpolant({Rvec,zvec}, ti, interp_mode, extrap_mode); end

% Interpolants for CARBON per-species (or per-charge)
F_niC = griddedInterpolant({Rvec,zvec}, niC, interp_mode, extrap_mode);
F_tiC = griddedInterpolant({Rvec,zvec}, tiC, interp_mode, extrap_mode);
F_uR_C = griddedInterpolant({Rvec,zvec}, uR_C, interp_mode, extrap_mode);
F_uZ_C = griddedInterpolant({Rvec,zvec}, uZ_C, interp_mode, extrap_mode);
F_uT_C = griddedInterpolant({Rvec,zvec}, uT_C, interp_mode, extrap_mode);

%% -----------------------------------------------
%% Quick 2-D maps for context
figure;
imagesc(Rvec, zvec, niC.'); set(gca,'YDir','normal');
xlabel('R [m]'); ylabel('Z [m]'); title([labelC ' density n_C']);
colorbar; axis image; hold on; plot(g.lim(1,:), g.lim(2,:), 'k', 'LineWidth', 1);

figure;
imagesc(Rvec, zvec, tiC.'); set(gca,'YDir','normal');
xlabel('R [m]'); ylabel('Z [m]'); title([labelC ' temperature T_C']);
colorbar; axis image; hold on; plot(g.lim(1,:), g.lim(2,:), 'k', 'LineWidth', 1);

figure;
imagesc(Rvec, zvec, uT_C.'); set(gca,'YDir','normal');
xlabel('R [m]'); ylabel('Z [m]'); title([labelC ' temperature T_C']);
colorbar; axis image; hold on; plot(g.lim(1,:), g.lim(2,:), 'k', 'LineWidth', 1);

%% -----------------------------------------------
%% Interpolate onto your 3-D geometry (centroid)
Rq = r_centroid; 
Zq = centroid(:,3);

% Single-fluid samples (still available if you want them)
ne_surf = F_ne(Rq, Zq);
te_surf = F_te(Rq, Zq);
vr_surf = F_vr(Rq, Zq);
vt_surf = F_vt(Rq, Zq);
vz_surf = F_vz(Rq, Zq);
if exist('F_ni','var'), ni_surf = F_ni(Rq, Zq); end
if exist('F_ti','var'), ti_surf = F_ti(Rq, Zq); end

% Carbon per-species (or per-charge) samples on the geometry
niC_surf = F_niC(Rq, Zq);
tiC_surf = F_tiC(Rq, Zq);
uR_C_surf = F_uR_C(Rq, Zq);
uZ_C_surf = F_uZ_C(Rq, Zq);
uT_C_surf = F_uT_C(Rq, Zq);
vpC_surf  = sqrt(uR_C_surf.^2 + uZ_C_surf.^2 + uT_C_surf.^2);

%% -----------------------------------------------
%% 3-D surface plots at geometry for CARBON Z=6
figure;
patch(transpose(X), transpose(Y), transpose(Z), niC_surf, ...
      'FaceAlpha', 1, 'EdgeAlpha', 0.25);
title([labelC ' density on geometry']);
xlabel('x [mm]'); ylabel('y [mm]'); zlabel('z [mm]');
axis equal; view(3); colorbar; colormap(parula);

figure;
patch(transpose(X), transpose(Y), transpose(Z), tiC_surf, ...
      'FaceAlpha', 1, 'EdgeAlpha', 0.25);
title([labelC ' temperature on geometry']);
xlabel('x [mm]'); ylabel('y [mm]'); zlabel('z [mm]');
axis equal; view(3); colorbar; colormap(parula);

figure;
patch(transpose(X), transpose(Y), transpose(Z), vpC_surf, ...
      'FaceAlpha', 1, 'EdgeAlpha', 0.25);
title([labelC ' speed |u| on geometry']);
xlabel('x [mm]'); ylabel('y [mm]'); zlabel('z [mm]');
axis equal; view(3); colorbar; colormap(parula);

%% -----------------------------------------------
%% Diagnostics
fprintf('%s |n_C| range on grid: [%.3e, %.3e]\n', labelC, min(niC(:)), max(niC(:)));
mC = [min(niC_surf), median(niC_surf,'omitnan'), max(niC_surf)];
fprintf('%s on geometry: min=%.3e, median=%.3e, max=%.3e\n', labelC, mC);

% Example probe
test_R = 2.29169; test_Z = 0.502665;
fprintf('n_C(grid) at [R=%.5f,Z=%.5f] ~ %.3e\n', test_R, test_Z, F_niC(test_R,test_Z));
[~,kpt] = min(hypot(Rq - test_R, Zq - test_Z));
fprintf('Nearest centroid d=%.3g m, n_C(surf)=%.3e\n', hypot(Rq(kpt)-test_R, Zq(kpt)-test_Z), niC_surf(kpt));

%% -----------------------------------------------
%% (Unchanged) EFIT B-field → centroids
z_efit = z_efit(:); r_efit = r_efit(:);
[z_efit, izEF] = sort(z_efit, 'ascend');
[r_efit, iREF] = sort(r_efit, 'ascend');
Br = Br(izEF, iREF);  Bt = Bt(izEF, iREF);  Bz = Bz(izEF, iREF);

F_Br = griddedInterpolant({z_efit, r_efit}, Br, 'linear', 'none');
F_Bt = griddedInterpolant({z_efit, r_efit}, Bt, 'linear', 'none');
F_Bz = griddedInterpolant({z_efit, r_efit}, Bz, 'linear', 'none');

br_surf = F_Br(centroid(:,3), r_centroid);
bt_surf = F_Bt(centroid(:,3), r_centroid);
bz_surf = F_Bz(centroid(:,3), r_centroid);

phi_centroid = atan2(centroid(:,2), centroid(:,1));
bx = double(br_surf.*cos(phi_centroid) - bt_surf.*sin(phi_centroid));
by = double(br_surf.*sin(phi_centroid) + bt_surf.*cos(phi_centroid));
bz_comp = double(bz_surf);

b_mag = sqrt(bx.^2 + by.^2 + bz_comp.^2);
ubx   = bx ./ b_mag;  uby = by ./ b_mag;  ubz = bz_comp ./ b_mag;

norm_vec_mag = sqrt(norm_vec(:,1).^2 + norm_vec(:,2).^2 + norm_vec(:,3).^2);
unorm_vec    = norm_vec ./ norm_vec_mag;

theta = acos( unorm_vec(:,1).*ubx + unorm_vec(:,2).*uby + unorm_vec(:,3).*ubz );
ii = find(theta > pi/2); theta(ii) = abs(theta(ii) - pi);

figure; histogram(theta);
xlabel('\theta [rad]'); ylabel('count'); title('\theta between normal and B');

%% ---------------- Helpers ----------------
function A = normalize_species_dims(Ain, nR, nZ, nS)
% Normalize a 3D species array to size [nR x nZ x nS].
% Accepts input as [nR x nZ x nS], [nS x nR x nZ], or [nR x nS x nZ].
    sz = size(Ain);
    if numel(sz)~=3
        error('Species array must be 3D, got size %s', mat2str(sz));
    end
    % If already [nR x nZ x nS]
    if isequal(sz, [nR, nZ, nS])
        A = Ain;
        return;
    end
    % Try permutations
    candidates = {
        [1 2 3], [nR nZ nS];   % identity
        [2 3 1], [nZ nS nR];
        [3 1 2], [nS nR nZ];
        [2 1 3], [nZ nR nS];
        [3 2 1], [nS nZ nR];
        [1 3 2], [nR nS nZ];
    };
    A = [];
    for k = 1:2:size(candidates,1)
        perm = candidates{k};
        tgt  = candidates{k+1};
        if all(sz(perm) == tgt)
            % We want final to be [nR x nZ x nS]
            % figure out which perm takes Ain -> [nR x nZ x nS]
            if all(tgt == [nR nZ nS])
                A = permute(Ain, perm);
                return;
            end
        end
    end
    % Fallback heuristic: find dims by size match
    [~, idxR] = min(abs(sz - nR));
    [~, idxZ] = min(abs(sz - nZ));
    idxS = setdiff(1:3, [idxR idxZ], 'stable');
    A = permute(Ain, [idxR idxZ idxS]);
    if ~isequal(size(A), [nR nZ nS])
        error('Could not coerce species array to [nR x nZ x nS]. Got %s', mat2str(size(A)));
    end
end