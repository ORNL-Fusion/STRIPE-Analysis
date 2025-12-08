load('multiSpecies_data_200882_test.mat')

%% ------------------- psiN on (X,Y) and Masks -------------------
psiN_flat = calc_psiN(g, X(:), Y(:), []);
psiN      = reshape(psiN_flat, size(X));

psiN_mask = 0.86;                      % your earlier choice (change if needed)
mask_psi  = psiN < psiN_mask;          % psiN-based mask (inside LCFS)

% NEW: mask outside the outer boundary (limiter/vessel)
% NOTE: X==R, Y==Z in your pipeline; g.lim is [R;Z]
mask_outer = ~inpolygon(X, Y, g.lim(1,:), g.lim(2,:));

% NEW: combined mask if you want BOTH effects at once
mask_any = mask_psi | mask_outer;

% If you ONLY want to mask beyond the outer boundary, use mask_outer below.
% If you want both, use mask_any. (I'll use mask_any to satisfy "mask beyond
% the outer boundary" AND preserve your ψN mask.)
mask_use = mask_any;

% ------------------- Single-fluid masked -------------------
val_ne_masked      = val_ne;      
val_Te_masked      = val_Te;      
val_gradTi_masked  = val_gradTi;  

% OLD: mask velocities too by psiN
% val_vr_masked      = val_vr;      val_vr_masked(mask)      = 0;
% val_vz_masked      = val_vz;      val_vz_masked(mask)      = 0;
% val_vt_masked      = val_vt;      val_vt_masked(mask)      = 0;

% OLD: keep velocities unmasked by psiN
% val_vr_masked = val_vr;
% val_vz_masked = val_vz;
% val_vt_masked = val_vt;

% NEW: mask ALL single-fluid profiles beyond the outer boundary (and ψN if desired)
val_ne_masked(mask_use)     = 0;
val_Te_masked(mask_use)     = 0;
val_gradTi_masked(mask_use) = 0;

val_vr_masked = val_vr;  val_vr_masked(mask_use) = 0;
val_vz_masked = val_vz;  val_vz_masked(mask_use) = 0;
val_vt_masked = val_vt;  val_vt_masked(mask_use) = 0;

% ------------------- Multi-species masked -------------------
if ~exist('ns','var') || isempty(ns), ns = size(val_ni_mspecies,3); end

val_ni_mspecies_masked = val_ni_mspecies;
val_Ti_mspecies_masked = val_Ti_mspecies;

% OLD: do not mask multi-species velocities
% val_uR_mspecies_masked = val_uR_mspecies;
% val_uZ_mspecies_masked = val_uZ_mspecies;
% val_uT_mspecies_masked = val_uT_mspecies;

% NEW: mask multi-species velocities beyond the outer boundary (and ψN if desired)
val_uR_mspecies_masked = val_uR_mspecies;
val_uZ_mspecies_masked = val_uZ_mspecies;
val_uT_mspecies_masked = val_uT_mspecies;

for k = 1:ns
    % densities & temps (were already ψN-masked; now also outer-wall masked)
    tmp = val_ni_mspecies_masked(:,:,k); tmp(mask_use)=0; val_ni_mspecies_masked(:,:,k)=tmp;
    tmp = val_Ti_mspecies_masked(:,:,k); tmp(mask_use)=0; val_Ti_mspecies_masked(:,:,k)=tmp;

    % velocities: NEW — also mask outside limiter (and ψN if you keep mask_use)
    tmp = val_uR_mspecies_masked(:,:,k); tmp(mask_use)=0; val_uR_mspecies_masked(:,:,k)=tmp;
    tmp = val_uZ_mspecies_masked(:,:,k); tmp(mask_use)=0; val_uZ_mspecies_masked(:,:,k)=tmp;
    tmp = val_uT_mspecies_masked(:,:,k); tmp(mask_use)=0; val_uT_mspecies_masked(:,:,k)=tmp;
end

save('multiSpecies_data_200882_outerMAsk_test.mat')

%% --------------- NetCDF Export stays unchanged below ---------------
% (your existing NetCDF writing code goes here with *_masked variables)

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
hold on
 % LCFS/psiN mask on the same (Z,R) axes
psiN_mask = 0.95;
contour(Zz, Rr, psi_plot, [psiN_mask psiN_mask], 'k--', 'LineWidth', 1.2);
hold off;

    % 2) Example species slice (e.g., k = 6) — no transpose needed
k = 6;
figure; imagesc(Zz, Rr, uT_all(:,:,k)); axis xy; set(gca,'FontSize',12);
xlabel('Z [m]'); ylabel('R [m]'); colorbar;
title(sprintf('Masked U_t (species %d) on (R,Z) from NetCDF', k));
hold on;
% IMPORTANT: limiter is stored as [R;Z]; swap to (Z,R) when plotting
plot(g.lim(2,:), g.lim(1,:), 'k', 'LineWidth', 1.0);


hold off;

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
imagesc(Zz, Rr, ni_m(:,:,k));           % <-- X=Z, Y=R; no transpose needed
set(gca,'YDir','normal','FontSize',12);
xlabel('Z [m]'); ylabel('R [m]');
colorbar; title(sprintf('Masked U_t (species %d) on (R,Z) from NetCDF',k));
hold on;

% IMPORTANT: limiter is stored as [R;Z]; swap to (Z,R) when plotting
plot(g.lim(2,:), g.lim(1,:), 'k', 'LineWidth', 1.0);

% LCFS/psiN mask on the same (Z,R) axes
psiN_mask = 1;
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