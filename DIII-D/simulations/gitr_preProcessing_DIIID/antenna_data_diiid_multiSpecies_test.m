%===============================================================
% SOLPS (R,z) → 3D geometry surface sampling
% - Electrons: ne, Te, vR_e, vZ_e, vT_e are 2D fields → surface via scatteredInterpolant
% - Ions: pick species/charge from 3D cubes → same scatteredInterpolant strategy
% - Interior-only (psiN>=psi_cut) training, nearest extrapolation outside
%===============================================================

% close all; clear all;

%% -------------------- Files --------------------
% ncfile  = 'profilesDIIID_196154_multi_test.nc';  % has x,z, ne,te, vr,vz,vt, psiN, and *_all cubes
ncfile  = 'profilesDIIID_200882.nc';  % has x,z, ne,te, vr,vz,vt, psiN, and *_all cubes
% ncfile  = 'profilesDIIID_196154.nc';  % has x,z, ne,te, vr,vz,vt, psiN, and *_all cubes




%% -------------------- Expect geometry & EFIT in workspace --------------------
assert(exist('centroid','var')==1 && size(centroid,2)==3, 'centroid [N x 3] needed');
assert(exist('norm_vec','var')==1 && size(norm_vec,2)==3, 'norm_vec [N x 3] needed');
assert(exist('X','var')==1 && exist('Y','var')==1 && exist('Z','var')==1, 'X,Y,Z patches needed');
assert(exist('g','var')==1 && isfield(g,'lim'), 'Limiter g.lim needed');

%% -------------------- Read single-fluid & multispecies --------------------
info = ncinfo(ncfile); varnames = {info.Variables.Name};

R  = ncread(ncfile,'x');      % [nR]
z  = ncread(ncfile,'z');      % [nZ]
ne = ncread(ncfile,'ne');     % [nR x nZ]   (2D electrons)
te = ncread(ncfile,'te');     % [nR x nZ]   (2D electrons)
vr = ncread(ncfile,'vr');     % [nR x nZ]   (treat as electron flow fields)
vt = ncread(ncfile,'vt');     % [nR x nZ]
vz = ncread(ncfile,'vz');     % [nR x nZ]

Z_all = ncread(ncfile,'atomic_number');   % [ns]
q_all = ncread(ncfile,'charge_number');   % [ns]
ni_all_raw = ncread(ncfile,'ni_all');     % 3D
ti_all_raw = ncread(ncfile,'ti_all');
uR_all_raw = ncread(ncfile,'uR_all');
uZ_all_raw = ncread(ncfile,'uZ_all');
uT_all_raw = ncread(ncfile,'uT_all');

psiN_grid_raw = ncread(ncfile,'psiN');    % [nR x nZ]

nR = numel(R); nZ = numel(z); nS = numel(Z_all);

%% -------------------- Ensure (R,z) ordering & sort axes --------------------
Rvec = R(:); zvec = z(:);
[Rvec, iR] = sort(Rvec,'ascend');
[zvec, iz] = sort(zvec,'ascend');

ne = ensure_RZ_order(ne,nR,nZ); ne = ne(iR,iz);
te = ensure_RZ_order(te,nR,nZ); te = te(iR,iz);
vr = ensure_RZ_order(vr,nR,nZ); vr = vr(iR,iz);
vt = ensure_RZ_order(vt,nR,nZ); vt = vt(iR,iz);
vz = ensure_RZ_order(vz,nR,nZ); vz = vz(iR,iz);

to_RZNS = @(A) normalize_species_dims(A, nR, nZ, nS);
ni_all = to_RZNS(ni_all_raw); ni_all = ni_all(iR,iz,:);
ti_all = to_RZNS(ti_all_raw); ti_all = ti_all(iR,iz,:);
uR_all = to_RZNS(uR_all_raw); uR_all = uR_all(iR,iz,:);
uZ_all = to_RZNS(uZ_all_raw); uZ_all = uZ_all(iR,iz,:);
uT_all = to_RZNS(uT_all_raw); uT_all = uT_all(iR,iz,:);

psiN_grid = ensure_RZ_order(psiN_grid_raw,nR,nZ); psiN_grid = psiN_grid(iR,iz);

%% -------------------- Species choice (ions) --------------------
Z_target = 6;          % Carbon
use_single_charge = true;
q_sel = 6;

idxZ = find(Z_all==Z_target);
assert(~isempty(idxZ),'No Z=%d in file.',Z_target);

if ~use_single_charge
    wC = ni_all(:,:,idxZ); tiny=1e-60;
    denom = sum(wC,3,'omitnan'); denom(denom<tiny)=tiny;
    niC  = sum(wC,3,'omitnan');
    tiC  = sum(wC.*ti_all(:,:,idxZ),3,'omitnan')./denom;
    uR_C = sum(wC.*uR_all(:,:,idxZ),3,'omitnan')./denom;
    uZ_C = sum(wC.*uZ_all(:,:,idxZ),3,'omitnan')./denom;
    uT_C = sum(wC.*uT_all(:,:,idxZ),3,'omitnan')./denom;
    labelC = sprintf('Carbon (Z=%d) all charges',Z_target);
else
    idx_q = find((Z_all==Z_target)&(q_all==q_sel),1,'first');
    assert(~isempty(idx_q),'No Z=%d, q=%d in file.',Z_target,q_sel);
    niC  = ni_all(:,:,idx_q);
    tiC  = ti_all(:,:,idx_q);
    uR_C = uR_all(:,:,idx_q);
    uZ_C = uZ_all(:,:,idx_q);
    uT_C = uT_all(:,:,idx_q);
    labelC = sprintf('Carbon (Z=%d, q=%d)',Z_target,q_sel);
end
%% -------------------- 2-D context maps (poloidal R–z) --------------------
figure;
imagesc(zvec, Rvec, ne); set(gca,'YDir','normal');
xlabel('Z [m]'); ylabel('R [m]'); title([labelC ' density n_C']);
colorbar; axis tight; hold on; plot(g.lim(2,:), g.lim(1,:), 'k', 'LineWidth', 1);

figure;
imagesc(zvec, Rvec, te); set(gca,'YDir','normal');
xlabel('Z [m]'); ylabel('R [m]'); title([labelC ' temperature T_C (eV)']);
colorbar; axis tight; hold on; plot(g.lim(2,:), g.lim(1,:), 'k', 'LineWidth', 1);

figure;
imagesc(zvec, Rvec, vt); set(gca,'YDir','normal');
xlabel('Z [m]'); ylabel('R [m]'); title([labelC ' toroidal velocity u_T (m/s)']);
colorbar; axis tight; hold on; plot(g.lim(2,:), g.lim(1,:), 'k', 'LineWidth', 1);
%% -------------------- Build interior-only interpolants (same recipe for e- and C) --------------------
[ZZm, RRm] = meshgrid(zvec, Rvec);
psi_cut   = 1;
is_plasma = (psiN_grid >= psi_cut) & isfinite(psiN_grid);

% ---- ION (selected species) : keep your working version ----
valid_nC = is_plasma & isfinite(niC) & (niC > 0);
valid_tC = is_plasma & isfinite(tiC) & (tiC > 0);
valid_uR = is_plasma & isfinite(uR_C) & (abs(uR_C) > 0);
valid_uZ = is_plasma & isfinite(uZ_C) & (abs(uZ_C) > 0);
valid_uT = is_plasma & isfinite(uT_C) & (abs(uT_C) > 0);

F_nC = scatteredInterpolant(RRm(valid_nC), ZZm(valid_nC), niC(valid_nC), 'natural','nearest');
F_tC = scatteredInterpolant(RRm(valid_tC), ZZm(valid_tC), tiC(valid_tC), 'natural','nearest');
F_uR_C = scatteredInterpolant(RRm(valid_uR), ZZm(valid_uR), uR_C(valid_uR), 'natural','nearest');
F_uZ_C = scatteredInterpolant(RRm(valid_uZ), ZZm(valid_uZ), uZ_C(valid_uZ), 'natural','nearest');
F_uT_C = scatteredInterpolant(RRm(valid_uT), ZZm(valid_uT), uT_C(valid_uT), 'natural','nearest');

% ---- ELECTRON / BULK fields (2D): use the SAME interior-only, scatteredInterpolant approach ----
valid_ne = is_plasma & isfinite(ne) & (ne > 0);
valid_te = is_plasma & isfinite(te) & (te > 0);
valid_vr = is_plasma & isfinite(vr);    % allow zeros if physics says so
valid_vz = is_plasma & isfinite(vz);
valid_vt = is_plasma & isfinite(vt);

F_ne = scatteredInterpolant(RRm(valid_ne), ZZm(valid_ne), ne(valid_ne), 'natural','nearest');
F_te = scatteredInterpolant(RRm(valid_te), ZZm(valid_te), te(valid_te), 'natural','nearest');
F_vr = scatteredInterpolant(RRm(valid_vr), ZZm(valid_vr), vr(valid_vr), 'natural','nearest');
F_vz = scatteredInterpolant(RRm(valid_vz), ZZm(valid_vz), vz(valid_vz), 'natural','nearest');
F_vt = scatteredInterpolant(RRm(valid_vt), ZZm(valid_vt), vt(valid_vt), 'natural','nearest');

%% -------------------- Interpolate to 3-D geometry surface --------------------
r_centroid = hypot(centroid(:,1), centroid(:,2));
Rq = r_centroid;   Zq = centroid(:,3);

% ---- Electrons (2D) on surface ----
ne_surf = F_ne(Rq, Zq);
te_surf = F_te(Rq, Zq);
vr_e_surf = F_vr(Rq, Zq);
vz_e_surf = F_vz(Rq, Zq);
vt_e_surf = F_vt(Rq, Zq);
vp_e_surf = sqrt(vr_e_surf.^2 + vz_e_surf.^2 + vt_e_surf.^2);

% ---- Selected ion (species/charge) on surface (unchanged, your working method) ----
niC_surf  = F_nC(Rq, Zq);
tiC_surf  = F_tC(Rq, Zq);
uR_C_surf = F_uR_C(Rq, Zq);
uZ_C_surf = F_uZ_C(Rq, Zq);
uT_C_surf = F_uT_C(Rq, Zq);
vpC_surf  = sqrt(uR_C_surf.^2 + uZ_C_surf.^2 + uT_C_surf.^2);


%% -------------------- Surface plots --------------------
figure;
patch(transpose(X), transpose(Y), transpose(Z), niC_surf, ...
      'FaceAlpha', 1, 'EdgeAlpha', 0.25);
title([labelC ' n_C on geometry']);
xlabel('x [mm]'); ylabel('y [mm]'); zlabel('z [mm]');
axis equal; view(3); colorbar; colormap(parula);

figure;
patch(transpose(X), transpose(Y), transpose(Z), tiC_surf, ...
      'FaceAlpha', 1, 'EdgeAlpha', 0.25);
title([labelC ' T_C on geometry (eV)']);
xlabel('x [mm]'); ylabel('y [mm]'); zlabel('z [mm]');
axis equal; view(3); colorbar; colormap(parula);

figure;
patch(transpose(X), transpose(Y), transpose(Z), vpC_surf, ...
      'FaceAlpha', 1, 'EdgeAlpha', 0.25);
title([labelC ' speed |u_C| on geometry (m/s)']);
xlabel('x [mm]'); ylabel('y [mm]'); zlabel('z [mm]');
axis equal; view(3); colorbar; colormap(parula);

figure;
patch(transpose(X), transpose(Y), transpose(Z), ne_surf, ...
      'FaceAlpha', 1, 'EdgeAlpha', 0.25);
title([labelC ' n_e on geometry']);
xlabel('x [mm]'); ylabel('y [mm]'); zlabel('z [mm]');
axis equal; view(3); colorbar; colormap(parula);

figure;
patch(transpose(X), transpose(Y), transpose(Z), te_surf, ...
      'FaceAlpha', 1, 'EdgeAlpha', 0.25);
title([labelC ' T_e on geometry (eV)']);
xlabel('x [mm]'); ylabel('y [mm]'); zlabel('z [mm]');
axis equal; view(3); colorbar; colormap(parula);

figure;
patch(transpose(X), transpose(Y), transpose(Z), vp_e_surf, ...
      'FaceAlpha', 1, 'EdgeAlpha', 0.25);
title([labelC ' speed |u_e| on geometry (m/s)']);
xlabel('x [mm]'); ylabel('y [mm]'); zlabel('z [mm]');
axis equal; view(3); colorbar; colormap(parula);
%% -------------------- Diagnostics --------------------
psi_cut = 0.86;
F_psi = griddedInterpolant({Rvec,zvec}, psiN_grid, 'linear','nearest');
psi_centroid = F_psi(Rq,Zq);
fprintf('Centroids with psiN < %.2f: %d / %d\n', psi_cut, sum(psi_centroid<psi_cut), numel(psi_centroid));

fprintf('Electron n_e surf:   min=%.3e  med=%.3e  max=%.3e\n', min(ne_surf), median(ne_surf,'omitnan'), max(ne_surf));
fprintf('Electron T_e surf:   min=%.3e  med=%.3e  max=%.3e\n', min(te_surf), median(te_surf,'omitnan'), max(te_surf));
fprintf('Electron |v| surf:   min=%.3e  med=%.3e  max=%.3e\n', min(vp_e_surf), median(vp_e_surf,'omitnan'), max(vp_e_surf));

fprintf('%s n_i surf:         min=%.3e  med=%.3e  max=%.3e\n', labelC, min(niC_surf), median(niC_surf,'omitnan'), max(niC_surf));
fprintf('%s T_i surf:         min=%.3e  med=%.3e  max=%.3e\n', labelC, min(tiC_surf), median(tiC_surf,'omitnan'), max(tiC_surf));
fprintf('%s |u| surf:         min=%.3e  med=%.3e  max=%.3e\n', labelC, min(vpC_surf), median(vpC_surf,'omitnan'), max(vpC_surf));

% Example probe
test_R = mean(Rvec); test_Z = mean(zvec);
[~,kpt] = min(hypot(Rq - test_R, Zq - test_Z));
fprintf('Probe near [R=%.3f,Z=%.3f]: n_e=%.3e, T_e=%.3e, |v_e|=%.3e, %s n_i=%.3e, T_i=%.3e, |u|=%.3e\n', ...
    test_R, test_Z, ne_surf(kpt), te_surf(kpt), vp_e_surf(kpt), ...
    labelC, niC_surf(kpt), tiC_surf(kpt), vpC_surf(kpt));

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
    sz = size(Ain);
    if numel(sz)~=2, error('Need 2-D array'); end
    if isequal(sz,[nR,nZ]), A=Ain; return;
    elseif isequal(sz,[nZ,nR]), A=Ain.'; return;
    else, error('Cannot coerce %s to [nR x nZ] or [nZ x nR]', mat2str(sz));
    end
end

function A = normalize_species_dims(Ain, nR, nZ, nS)
    sz = size(Ain);
    if numel(sz)~=3, error('Species array must be 3D, got %s', mat2str(sz)); end
    if isequal(sz,[nR,nZ,nS]), A=Ain; return; end
    permsList = [1 2 3;1 3 2;2 1 3;2 3 1;3 1 2;3 2 1];
    for k=1:6
        if all(size(Ain,permsList(k,:))==[nR nZ nS])
            A = permute(Ain,permsList(k,:)); return;
        end
    end
    [~,iR]=min(abs(sz-nR)); [~,iZ]=min(abs(sz-nZ)); iS=setdiff(1:3,[iR iZ],'stable');
    A = permute(Ain,[iR iZ iS]);
    if ~isequal(size(A),[nR nZ nS]), error('Could not coerce species array'); end
end