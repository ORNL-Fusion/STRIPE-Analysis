%===============================================================
% Surface sampling of SOLPS multi-species fields (R–z) → 3-D geometry
% Fix for "zero flow on surface": build velocity interpolants only from
% valid interior cells and use nearest extrapolation outside.
%===============================================================

% close all;
% clear all;

%% -------------------- Constants (unchanged) --------------------
ME   = 9.10938356e-31;
MI   = 1.6737236e-27;
Q    = 1.60217662e-19;
EPS0 = 8.854187e-12;
amu  = 18;

%% -------------------- Files --------------------
ncfile  = 'profilesDIIID_196154_multi_test.nc';  % multi-species on (R,z,ns)

%% -------------------- Expect geometry in workspace --------------------
% Required variables (must already be in workspace):
%   centroid [N x 3]  - triangle centroids (x,y,z)
%   norm_vec [N x 3]  - triangle outward normals
%   X,Y,Z             - patch vertices for plotting (same topology as centroids)
%   g.lim (2 x M)     - limiter polygon (R;Z) for overlays
%   Br,Bt,Bz,z_efit,r_efit  - EFIT fields on (z,r) for B at centroids
assert(exist('centroid','var')==1 && size(centroid,2)==3, ...
  'centroid [N x 3] must be in workspace.');
assert(exist('norm_vec','var')==1 && size(norm_vec,2)==3, ...
  'norm_vec [N x 3] must be in workspace.');
assert(exist('X','var')==1 && exist('Y','var')==1 && exist('Z','var')==1, ...
  'X,Y,Z patch arrays must be in workspace.');
assert(exist('g','var')==1 && isfield(g,'lim'), 'EFIT limiter g.lim must be present.');

%% -------------------- Read single-fluid & multispecies --------------------
info = ncinfo(ncfile);
varnames = {info.Variables.Name};

R  = ncread(ncfile,'x');      % [nR]
z  = ncread(ncfile,'z');      % [nZ]
ne = ncread(ncfile,'ne');     % [nR x nZ]
te = ncread(ncfile,'te');     % [nR x nZ]
vr = ncread(ncfile,'vr');     % [nR x nZ]
vt = ncread(ncfile,'vt');     % [nR x nZ]
vz = ncread(ncfile,'vz');     % [nR x nZ]

has_ni = any(strcmp(varnames,'ni'));
has_ti = any(strcmp(varnames,'ti'));
if has_ni, ni = ncread(ncfile,'ni'); else, ni = []; end
if has_ti, ti = ncread(ncfile,'ti'); else, ti = []; end

% Multispecies metadata and fields
Z_all = ncread(ncfile,'atomic_number');   % [ns]
q_all = ncread(ncfile,'charge_number');   % [ns]
ni_all_raw = ncread(ncfile,'ni_all');     % unknown orientation
ti_all_raw = ncread(ncfile,'ti_all');
uR_all_raw = ncread(ncfile,'uR_all');
uZ_all_raw = ncread(ncfile,'uZ_all');
uT_all_raw = ncread(ncfile,'uT_all');

% psiN for masking diagnostics / interior detection
psiN_grid_raw = ncread(ncfile,'psiN');    % [nR x nZ] as written by exporter

nR = numel(R); nZ = numel(z); nS = numel(Z_all);

%% -------------------- Enforce (R,z) orientation & sort axes --------------------
Rvec = R(:); zvec = z(:);
[Rvec, iR] = sort(Rvec, 'ascend');
[zvec, iz] = sort(zvec, 'ascend');

ne = ensure_RZ_order(ne, nR, nZ); ne = ne(iR, iz);
te = ensure_RZ_order(te, nR, nZ); te = te(iR, iz);
vr = ensure_RZ_order(vr, nR, nZ); vr = vr(iR, iz);
vt = ensure_RZ_order(vt, nR, nZ); vt = vt(iR, iz);
vz = ensure_RZ_order(vz, nR, nZ); vz = vz(iR, iz);
if ~isempty(ni), ni = ensure_RZ_order(ni, nR, nZ); ni = ni(iR, iz); end
if ~isempty(ti), ti = ensure_RZ_order(ti, nR, nZ); ti = ti(iR, iz); end

% Species arrays → [nR x nZ x nS], then sort first two dims
to_RZNS = @(A) normalize_species_dims(A, nR, nZ, nS);
ni_all = to_RZNS(ni_all_raw); ni_all = ni_all(iR, iz, :);
ti_all = to_RZNS(ti_all_raw); ti_all = ti_all(iR, iz, :);
uR_all = to_RZNS(uR_all_raw); uR_all = uR_all(iR, iz, :);
uZ_all = to_RZNS(uZ_all_raw); uZ_all = uZ_all(iR, iz, :);
uT_all = to_RZNS(uT_all_raw); uT_all = uT_all(iR, iz, :);

psiN_grid = ensure_RZ_order(psiN_grid_raw, nR, nZ); psiN_grid = psiN_grid(iR, iz);

%% -------------------- Select species: Carbon example --------------------
Z_target = 6;                % Carbon
use_single_charge = true;    % true → pick specific charge; false → sum all carbon charges
q_sel = 4;                   % choose a carbon charge state if use_single_charge=true

idxZ = find(Z_all == Z_target);
if isempty(idxZ)
    error('No species with atomic_number == %d found.', Z_target);
end

if ~use_single_charge
    % Sum all carbon charges (density-weighted averages for T and u)
    wC    = ni_all(:,:,idxZ);                    % [nR x nZ x nC]
    tiny  = 1e-60;
    denom = sum(wC, 3, 'omitnan'); denom(denom<tiny) = tiny;

    niC  = sum(wC, 3, 'omitnan');
    tiC  = sum(wC .* ti_all(:,:,idxZ), 3, 'omitnan') ./ denom;
    uR_C = sum(wC .* uR_all(:,:,idxZ), 3, 'omitnan') ./ denom;
    uZ_C = sum(wC .* uZ_all(:,:,idxZ), 3, 'omitnan') ./ denom;
    uT_C = sum(wC .* uT_all(:,:,idxZ), 3, 'omitnan') ./ denom;

    labelC = sprintf('Carbon (Z=%d) total over charges', Z_target);
else
    idx_q = find((Z_all==Z_target) & (q_all==q_sel), 1, 'first');
    if isempty(idx_q)
        error('Requested Z=%d, q=%d not found.', Z_target, q_sel);
    end
    niC  = ni_all(:,:,idx_q);
    tiC  = ti_all(:,:,idx_q);
    uR_C = uR_all(:,:,idx_q);
    uZ_C = uZ_all(:,:,idx_q);
    uT_C = uT_all(:,:,idx_q);

    labelC = sprintf('Carbon (Z=%d, q=%d)', Z_target, q_sel);
end

%% -------------------- 2-D context maps (poloidal R–z) --------------------
figure;
imagesc(zvec, Rvec, niC); set(gca,'YDir','normal');
xlabel('Z [m]'); ylabel('R [m]'); title([labelC ' density n_C']);
colorbar; axis tight; hold on; plot(g.lim(2,:), g.lim(1,:), 'k', 'LineWidth', 1);

figure;
imagesc(zvec, Rvec, tiC); set(gca,'YDir','normal');
xlabel('Z [m]'); ylabel('R [m]'); title([labelC ' temperature T_C (eV)']);
colorbar; axis tight; hold on; plot(g.lim(2,:), g.lim(1,:), 'k', 'LineWidth', 1);

figure;
imagesc(zvec, Rvec, uT_C); set(gca,'YDir','normal');
xlabel('Z [m]'); ylabel('R [m]'); title([labelC ' toroidal velocity u_T (m/s)']);
colorbar; axis tight; hold on; plot(g.lim(2,:), g.lim(1,:), 'k', 'LineWidth', 1);

%% -------------------- Build interpolants (key fix for velocities) --------------------
% Gridded interpolants for scalars (OK to use masked zeros if desired)
F_ne = griddedInterpolant({Rvec,zvec}, ne, 'linear', 'none');
F_te = griddedInterpolant({Rvec,zvec}, te, 'linear', 'none');
if ~isempty(ni), F_ni = griddedInterpolant({Rvec,zvec}, ni, 'linear', 'none'); end
if ~isempty(ti), F_ti = griddedInterpolant({Rvec,zvec}, ti, 'linear', 'none'); end

% -------- This is the important part --------
% Use ONLY interior, nonzero cells to build velocity interpolants.
% -------- Scalars (same approach as flows) --------
[ZZm, RRm] = meshgrid(zvec, Rvec);
psi_cut = 0.86;
is_plasma = (psiN_grid >= psi_cut) & isfinite(psiN_grid);

valid_nC = is_plasma & isfinite(niC) & (niC > 0);
valid_tC = is_plasma & isfinite(tiC) & (tiC > 0);



R_n = RRm(valid_nC);  Z_n = ZZm(valid_nC);  V_n = niC(valid_nC);
R_t = RRm(valid_tC);  Z_t = ZZm(valid_tC);  V_t = tiC(valid_tC);

F_nC = scatteredInterpolant(R_n, Z_n, V_n, 'natural', 'nearest');
F_tC = scatteredInterpolant(R_t, Z_t, V_t, 'natural', 'nearest');

valid_uR = is_plasma & isfinite(uR_C) & (abs(uR_C) > 0);
valid_uZ = is_plasma & isfinite(uZ_C) & (abs(uZ_C) > 0);
valid_uT = is_plasma & isfinite(uT_C) & (abs(uT_C) > 0);

% Flatten valid interior points
R_u = RRm(valid_uR);  Z_u = ZZm(valid_uR);  V_uR = uR_C(valid_uR);
R_w = RRm(valid_uZ);  Z_w = ZZm(valid_uZ);  V_uZ = uZ_C(valid_uZ);
R_t = RRm(valid_uT);  Z_t = ZZm(valid_uT);  V_uT = uT_C(valid_uT);

% Scatter interpolants over interior; 'nearest' extrapolation outside
F_uR_C = scatteredInterpolant(R_u, Z_u, V_uR, 'natural', 'nearest');
F_uZ_C = scatteredInterpolant(R_w, Z_w, V_uZ, 'natural', 'nearest');
F_uT_C = scatteredInterpolant(R_t, Z_t, V_uT, 'natural', 'nearest');

%% -------------------- Interpolate to 3-D geometry surface --------------------
r_centroid = hypot(centroid(:,1), centroid(:,2));
Rq = r_centroid;
Zq = centroid(:,3);



% Scalars (still gridded)
ne_surf = F_ne(Rq, Zq);
te_surf = F_te(Rq, Zq);
if exist('F_ni','var'), ni_surf = F_ni(Rq, Zq); end
if exist('F_ti','var'), ti_surf = F_ti(Rq, Zq); end

% Species (velocity) using scatter interpolants built above
% niC_surf  = scattered_sample_grid(niC,  Rvec, zvec, Rq, Zq);  % helper below
% tiC_surf  = scattered_sample_grid(tiC,  Rvec, zvec, Rq, Zq);

niC_surf = F_nC(Rq, Zq);
tiC_surf = F_tC(Rq, Zq);

uR_C_surf = F_uR_C(Rq, Zq);
uZ_C_surf = F_uZ_C(Rq, Zq);
uT_C_surf = F_uT_C(Rq, Zq);
vpC_surf  = sqrt(uR_C_surf.^2 + uZ_C_surf.^2 + uT_C_surf.^2);

% uR_C_surf = F_uR_C(Rq, Zq);
% uZ_C_surf = F_uZ_C(Rq, Zq);
% uT_C_surf = F_uT_C(Rq, Zq);
% vpC_surf  = sqrt(uR_C_surf.^2 + uZ_C_surf.^2 + uT_C_surf.^2);

%% -------------------- Surface plots --------------------
figure;
patch(transpose(X), transpose(Y), transpose(Z), niC_surf, ...
      'FaceAlpha', 1, 'EdgeAlpha', 0.25);
title([labelC ' density on geometry']);
xlabel('x [mm]'); ylabel('y [mm]'); zlabel('z [mm]');
axis equal; view(3); colorbar; colormap(parula);

figure;
patch(transpose(X), transpose(Y), transpose(Z), tiC_surf, ...
      'FaceAlpha', 1, 'EdgeAlpha', 0.25);
title([labelC ' temperature on geometry (eV)']);
xlabel('x [mm]'); ylabel('y [mm]'); zlabel('z [mm]');
axis equal; view(3); colorbar; colormap(parula);

figure;
patch(transpose(X), transpose(Y), transpose(Z), vpC_surf, ...
      'FaceAlpha', 1, 'EdgeAlpha', 0.25);
title([labelC ' speed |u| on geometry (m/s)']);
xlabel('x [mm]'); ylabel('y [mm]'); zlabel('z [mm]');
axis equal; view(3); colorbar; colormap(parula);

%% -------------------- Diagnostics --------------------
fprintf('%s |n_C| range on grid: [%.3e, %.3e]\n', labelC, min(niC(:)), max(niC(:)));
mC = [min(niC_surf), median(niC_surf,'omitnan'), max(niC_surf)];
fprintf('%s on geometry: min=%.3e, median=%.3e, max=%.3e\n', labelC, mC);

% How many centroids lie outside psi_cut?
F_psi = griddedInterpolant({Rvec,zvec}, psiN_grid, 'linear', 'nearest');
psi_centroid = F_psi(Rq, Zq);
fprintf('Centroids with psiN < %.2f: %d / %d\n', psi_cut, sum(psi_centroid < psi_cut), numel(psi_centroid));

% Example probe
test_R = 2.29169; test_Z = 0.502665;
fprintf('n_C(grid) at [R=%.5f,Z=%.5f] ~ %.3e\n', test_R, test_Z, ...
    scattered_sample_grid(niC, Rvec, zvec, test_R, test_Z));
[~,kpt] = min(hypot(Rq - test_R, Zq - test_Z));
fprintf('Nearest centroid d=%.3g m, n_C(surf)=%.3e\n', hypot(Rq(kpt)-test_R, Zq(kpt)-test_Z), niC_surf(kpt));

%% -------------------- EFIT B-field → centroids (unchanged) --------------------
assert(exist('Br','var')==1 && exist('Bt','var')==1 && exist('Bz','var')==1 ...
    && exist('z_efit','var')==1 && exist('r_efit','var')==1, ...
    'EFIT fields Br,Bt,Bz,z_efit,r_efit must be in workspace.');

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

figure; histogram(57.296*theta);
xlabel('\theta [rad]'); ylabel('count'); title('\theta between normal and B');



%===============================================================
% -------------------- Helpers --------------------
%===============================================================
function A = ensure_RZ_order(Ain, nR, nZ)
% Coerce a 2-D array to [nR x nZ]. If it is [nZ x nR], transpose.

    sz = size(Ain);
    if numel(sz) ~= 2
        error('Array must be 2-D; got %s', mat2str(sz));
    end
    if isequal(sz, [nR, nZ])
        A = Ain;
        return;
    elseif isequal(sz, [nZ, nR])
        A = Ain.';
        return;
    else
        error('Cannot coerce array of size %s to [nR x nZ] or [nZ x nR].', mat2str(sz));
    end
end

function A = normalize_species_dims(Ain, nR, nZ, nS)
% Normalize a 3D species array to size [nR x nZ x nS].
% Accepts input as [nR x nZ x nS], [nS x nR x nZ], [nR x nS x nZ], etc.

    sz = size(Ain);
    if numel(sz)~=3
        error('Species array must be 3D, got size %s', mat2str(sz));
    end

    % Already correct
    if isequal(sz, [nR, nZ, nS])
        A = Ain;
        return;
    end

    permsList = [
        1 2 3;
        1 3 2;
        2 1 3;
        2 3 1;
        3 1 2;
        3 2 1];
    targets = [
        nR nZ nS;
        nR nS nZ;
        nZ nR nS;
        nZ nS nR;
        nS nR nZ;
        nS nZ nR];

    % Try simple single permute cases
    for k = 1:size(permsList,1)
        if all(sz(permsList(k,:)) == targets(k,:))
            A = permute(Ain, permsList(k,:));
            % If this didn't land exactly on [nR nZ nS], try one more perm
            szA = size(A);
            if ~isequal(szA, [nR nZ nS])
                for kk = 1:size(permsList,1)
                    szA = size(A);                % recompute after any change
                    if all(szA(permsList(kk,:)) == [nR nZ nS])
                        A = permute(A, permsList(kk,:));
                        break;
                    end
                end
            end
            if isequal(size(A), [nR nZ nS]), return; end
        end
    end

    % Fallback heuristic: choose dims closest to nR and nZ
    [~, idxR] = min(abs(sz - nR));
    [~, idxZ] = min(abs(sz - nZ));
    idxS = setdiff(1:3, [idxR idxZ], 'stable');
    A = permute(Ain, [idxR idxZ idxS]);

    if ~isequal(size(A), [nR nZ nS])
        error('Could not coerce species array to [nR x nZ x nS]. Got %s', mat2str(size(A)));
    end
end

function vq = scattered_sample_grid(A_RZ, Rvec, zvec, Rq, Zq)
% Convenience: sample a 2-D grid A(R,z) at query points using scatteredInterpolant
% (natural/nearest), which behaves better outside domain than griddedInterpolant
% when the original array may have masked zeros.

    [ZZ, RR] = meshgrid(zvec, Rvec);
    mask = isfinite(A_RZ);
    F = scatteredInterpolant(RR(mask), ZZ(mask), A_RZ(mask), 'natural', 'nearest');
    vq = F(Rq, Zq);
end