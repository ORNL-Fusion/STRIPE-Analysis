%% Full SOLPS→GITR Processing, Masking, NetCDF Export, and Profile Visualization

% This script:
%   1) Loads extrapolated SOLPS data on the GITR grid
%   2) Cleans NaNs
%   3) Computes normalized poloidal flux (ψ_N)
%   4) Applies a mask (ψ_N < 0.8) to all fields
%   5) Plots ψ_N and masked fields
%   6) Exports masked profiles to a NetCDF file with explicit variable definitions
%   7) Reads back the NetCDF and visualizes 2D maps and 1D slices

close all;
clear all;

%% 1) Load data
load('extrapolated_data_multiSpecies_196154.mat');  
% expects: X, Y, val_ne, val_Te, val_gradTi, val_gradTe, val_vr, val_vz, val_vt, g

%% 2) Replace NaNs with zeros
val_ne(isnan(val_ne))         = 0;
val_Te(isnan(val_Te))         = 0;
val_gradTi(isnan(val_gradTi)) = 0;
val_vr(isnan(val_vr))         = 0;
val_vz(isnan(val_vz))         = 0;
val_vt(isnan(val_vt))         = 0;
% Uncomment if you have magnetic-field arrays:
% val_br(isnan(val_br))     = 0;
% val_bz(isnan(val_bz))     = 0;
% val_bt(isnan(val_bt))     = 0;

%% 3) Compute ψ_N on the 2D grid
psiN_flat = calc_psiN(g, X(:), Y(:), []);
psiN      = reshape(psiN_flat, size(X));

%% 4) Apply mask (ψ_N < 0.8) to all variables
psiN_mask=0.86;
mask = psiN < psiN_mask;



val_ne_masked      = val_ne;      val_ne_masked(mask)      = 0;
val_Te_masked      = val_Te;      val_Te_masked(mask)      = 0;
val_gradTi_masked  = val_gradTi;  val_gradTi_masked(mask)  = 0;
% val_gradTe_masked  = val_gradTe;  val_gradTe_masked(mask)  = 0;
val_vr_masked      = val_vr;      val_vr_masked(mask)      = 0;
val_vz_masked      = val_vz;      val_vz_masked(mask)      = 0;
val_vt_masked      = val_vt;      val_vt_masked(mask)      = 0;
% val_br_masked    = val_br;    val_br_masked(mask)    = 0;
% val_bz_masked    = val_bz;    val_bz_masked(mask)    = 0;
% val_bt_masked    = val_bt;    val_bt_masked(mask)    = 0;

%% 5a) Plot 2D map of ψ_N
figure;
imagesc([min(X(:)) max(X(:))], [min(Y(:)) max(Y(:))], psiN);
set(gca, 'YDir', 'normal');
axis equal tight;
xlabel('R [m]'); ylabel('Z [m]');
colorbar; colormap jet;
title('Normalized Poloidal Flux \psi_N');
hold on;
plot(g.lim(1,:), g.lim(2,:), 'k', 'LineWidth', 1);

%% 5b) Plot masked electron density
figure;
imagesc([min(X(:)) max(X(:))], [min(Y(:)) max(Y(:))], val_ne_masked);
set(gca, 'YDir', 'normal'); axis equal tight;
xlabel('R [m]'); ylabel('Z [m]');
colorbar; colormap parula;
title('Electron Density (n_e) Masked for \psi_N < 0.9');
hold on;
plot(g.lim(1,:), g.lim(2,:), 'k', 'LineWidth', 1);
contour(X, Y, psiN, [psiN_mask psiN_mask], 'k--', 'LineWidth', 1.5);
legend({'Limiter','\psi_N = 0.9'}, 'Location','best');

%% 5c) Plot other masked fields
fields = {...
    val_Te_masked,    'T_e [eV]',                 'Electron Temperature'; ...
    val_gradTi_masked,'\nabla_\parallel T_i',     'Ion Temperature Gradient'; ...
    val_vr_masked,    'v_r [m/s]',                 'Radial Velocity'; ...
    val_vz_masked,    'v_z [m/s]',                 'Poloidal Velocity'; ...
    val_vt_masked,    'v_\phi [m/s]',              'Toroidal Velocity'  ...
};
for k = 1:size(fields,1)
    figure;
    imagesc([min(X(:)) max(X(:))], [min(Y(:)) max(Y(:))], fields{k,1});
    set(gca, 'YDir', 'normal'); axis equal tight;
    xlabel('R [m]'); ylabel('Z [m]');
    colorbar; colormap parula;
    title([fields{k,3}, ' Masked (\psi_N < 0.8)']);
    hold on;
    plot(g.lim(1,:), g.lim(2,:), 'k', 'LineWidth', 1);
    contour(X, Y, psiN, [psiN_mask psiN_mask], 'k--', 'LineWidth', 1.2);
end

%% 6) Export masked profiles to NetCDF
x   = X(1,:);    % radial grid
z   = Y(:,1);    % axial grid
nR  = length(x);
nZ  = length(z);

ncid = netcdf.create('profilesDIIID_196154.nc','NC_WRITE');

% Define dimensions
dimR = netcdf.defDim(ncid,'x', nR);
dimZ = netcdf.defDim(ncid,'z', nZ);

% Define coordinate variables
varidR = netcdf.defVar(ncid,'x','double',dimR);
varidZ = netcdf.defVar(ncid,'z','double',dimZ);

% Define 2D profile variables
varidNe     = netcdf.defVar(ncid,'ne',     'double',[dimR dimZ]);
varidNi     = netcdf.defVar(ncid,'ni',     'double',[dimR dimZ]);
varidTe     = netcdf.defVar(ncid,'te',     'double',[dimR dimZ]);
varidTi     = netcdf.defVar(ncid,'ti',     'double',[dimR dimZ]);
varidGradTi = netcdf.defVar(ncid,'gradTi', 'double',[dimR dimZ]);
varidGradTe = netcdf.defVar(ncid,'gradTe', 'double',[dimR dimZ]);
varidVr     = netcdf.defVar(ncid,'vr',     'double',[dimR dimZ]);
varidVt     = netcdf.defVar(ncid,'vt',     'double',[dimR dimZ]);
varidVz     = netcdf.defVar(ncid,'vz',     'double',[dimR dimZ]);

% (Optional) magnetic field
% varidBr = netcdf.defVar(ncid,'br','double',[dimR dimZ]);
% varidBt = netcdf.defVar(ncid,'bt','double',[dimR dimZ]);
% varidBz = netcdf.defVar(ncid,'bz','double',[dimR dimZ]);

netcdf.endDef(ncid);

% Write data
netcdf.putVar(ncid,varidR,         x);
netcdf.putVar(ncid,varidZ,         z);
netcdf.putVar(ncid,varidNe,        val_ne_masked);
netcdf.putVar(ncid,varidNi,        val_ne_masked);      % n_i = n_e
netcdf.putVar(ncid,varidTe,        val_Te_masked);
netcdf.putVar(ncid,varidTi,        val_Te_masked);      % T_i = T_e
netcdf.putVar(ncid,varidGradTi,    val_gradTi_masked);
netcdf.putVar(ncid,varidGradTe,    val_gradTi_masked);
netcdf.putVar(ncid,varidVr,        val_vr_masked);
netcdf.putVar(ncid,varidVt,        val_vt_masked);
netcdf.putVar(ncid,varidVz,        val_vz_masked);

% (Optional) write B-fields
% netcdf.putVar(ncid,varidBr, val_br_masked);
% netcdf.putVar(ncid,varidBt, val_bt_masked);
% netcdf.putVar(ncid,varidBz, val_bz_masked);

netcdf.close(ncid);

%% === 6b) (NEW) Multispecies masking & NetCDF export ===
% Goal: mask and export ni/ti/uR/uZ/uT for all ion species on your (X,Y) grid.
% Assumptions:
%   - Preferred: 'extrapolated_data_196154.mat' already contains
%       val_ni_mspecies, val_Ti_mspecies, val_uR_mspecies, val_uZ_mspecies, val_uT_mspecies,
%       Z_all, q_all
%   - Fallback: read raw species arrays from 'interpolated_values.nc' and
%       map them to (X,Y) by linear interp (no extra extrapolation here).

fprintf('\n=== Multispecies export: preparing species cubes ===\n');

have_cubes = exist('val_ni_mspecies','var')==1 && ...
             exist('val_Ti_mspecies','var')==1 && ...
             exist('val_uR_mspecies','var')==1 && ...
             exist('val_uZ_mspecies','var')==1 && ...
             exist('val_uT_mspecies','var')==1;

have_meta  = exist('Z_all','var')==1 && exist('q_all','var')==1;

if ~have_cubes || ~have_meta
    % ---- Fallback path: read species from SOLPS netCDF and map to (X,Y) ----
    fileSOLPS = 'interpolated_values.nc';
    fprintf('Species cubes not found in MAT. Reading from %s ...\n', fileSOLPS);

    % metadata
    Z_all = ncread(fileSOLPS,'atomic_number');      % [ns]
    q_all = ncread(fileSOLPS,'charge_number');      % [ns]
    ns    = numel(Z_all);

    % native SOLPS grid for species
    rS_nc = ncread(fileSOLPS,'gridr');              % [nr]
    zS_nc = ncread(fileSOLPS,'gridz');              % [nz]
    nr_nc = numel(rS_nc); nz_nc = numel(zS_nc);

    % species 3D arrays as written by your C++ (dims = [nz nr ns] in file)
    ni_nc = ncread(fileSOLPS,'ni_all');             % [nz nr ns]
    ti_nc = ncread(fileSOLPS,'ti_all');             % [nz nr ns]
    uR_nc = ncread(fileSOLPS,'uR_all');             % [nz nr ns]
    uZ_nc = ncread(fileSOLPS,'uZ_all');             % [nz nr ns]
    uT_nc = ncread(fileSOLPS,'uT_all');             % [nz nr ns]

    % allocate cubes on (X,Y): sizes [nR x nZ x ns] to match your 2D vars layout
    nR = size(X,2);
    nZ = size(Y,1);
    val_ni_mspecies = zeros(nR, nZ, ns);
    val_Ti_mspecies = zeros(nR, nZ, ns);
    val_uR_mspecies = zeros(nR, nZ, ns);
    val_uZ_mspecies = zeros(nR, nZ, ns);
    val_uT_mspecies = zeros(nR, nZ, ns);

    % map each species to (X,Y) by 2-D interp on the native (rS_nc,zS_nc) grid
    % NOTE: arrays for interp2 are [nr x nz] → transpose to [nz x nr].
    for k = 1:ns
        % gather [nr x nz] for each field
        ni_k = squeeze(permute(ni_nc(:,:,k), [2 1]));   % [nr x nz]
        ti_k = squeeze(permute(ti_nc(:,:,k), [2 1]));   % [nr x nz]
        uR_k = squeeze(permute(uR_nc(:,:,k), [2 1]));   % [nr x nz]
        uZ_k = squeeze(permute(uZ_nc(:,:,k), [2 1]));   % [nr x nz]
        uT_k = squeeze(permute(uT_nc(:,:,k), [2 1]));   % [nr x nz]

        % map to (X,Y). Your convention elsewhere: interp2(r, z, field', X, Y, ...)
        val_ni_mspecies(:,:,k) = interp2(rS_nc, zS_nc, ni_k', X, Y, 'linear', 0);
        val_Ti_mspecies(:,:,k) = interp2(rS_nc, zS_nc, ti_k', X, Y, 'linear', 0);
        val_uR_mspecies(:,:,k) = interp2(rS_nc, zS_nc, uR_k', X, Y, 'linear', 0);
        val_uZ_mspecies(:,:,k) = interp2(rS_nc, zS_nc, uZ_k', X, Y, 'linear', 0);
        val_uT_mspecies(:,:,k) = interp2(rS_nc, zS_nc, uT_k', X, Y, 'linear', 0);
    end
else
    ns = size(val_ni_mspecies,3);
    fprintf('Using species cubes already in memory (ns=%d).\n', ns);
end

% ---- Apply the SAME ψ_N mask to all species channels ----
mask3 = repmat(mask, [1 1 ns]);    % (nZ x nR x ns) or (nR x nZ x ns)? — ensure orientation
% Your 2D arrays are [nR x nZ]; mask is [nZ x nR] from reshape(psiN,size(X)).
% We used imagesc([minX maxX],[minY maxY], val_*) with val_* sized like X (nZ x nR) originally,
% but above you stored val_* as [nR x nZ]. So transpose mask to [nR x nZ] first:
mask_T = mask';                     % [nR x nZ]
mask3  = repmat(mask_T, [1 1 ns]);  % [nR x nZ x ns]

val_ni_mspecies_masked = val_ni_mspecies;
val_Ti_mspecies_masked = val_Ti_mspecies;
val_uR_mspecies_masked = val_uR_mspecies;
val_uZ_mspecies_masked = val_uZ_mspecies;
val_uT_mspecies_masked = val_uT_mspecies;

val_ni_mspecies_masked(mask3) = 0;
val_Ti_mspecies_masked(mask3) = 0;
val_uR_mspecies_masked(mask3) = 0;
val_uZ_mspecies_masked(mask3) = 0;
val_uT_mspecies_masked(mask3) = 0;

% ---- Write a multispecies NetCDF (keeps your original single-species file unchanged) ----
outfile_ms = 'profilesDIIID_196154_multispecies.nc';
if exist(outfile_ms,'file'), delete(outfile_ms); end
ncid2 = netcdf.create(outfile_ms, 'NC_WRITE');

% dimensions
dimR = netcdf.defDim(ncid2,'x', size(X,2));  % nR
dimZ = netcdf.defDim(ncid2,'z', size(Y,1));  % nZ
dimS = netcdf.defDim(ncid2,'ns', ns);

% coords
varidR = netcdf.defVar(ncid2,'x','double',dimR);
varidZ = netcdf.defVar(ncid2,'z','double',dimZ);

% species meta
varidZnum = netcdf.defVar(ncid2,'atomic_number','double',dimS);
varidQnum = netcdf.defVar(ncid2,'charge_number','double',dimS);

% multi-species fields (use [x z ns] ordering → arrays sized [nR x nZ x ns])
varid_ni = netcdf.defVar(ncid2,'ni_all','double',[dimR dimZ dimS]);
varid_ti = netcdf.defVar(ncid2,'ti_all','double',[dimR dimZ dimS]);
varid_uR = netcdf.defVar(ncid2,'uR_all','double',[dimR dimZ dimS]);
varid_uZ = netcdf.defVar(ncid2,'uZ_all','double',[dimR dimZ dimS]);
varid_uT = netcdf.defVar(ncid2,'uT_all','double',[dimR dimZ dimS]);

% also store your masked bulk fields for convenience
varidNe = netcdf.defVar(ncid2,'ne','double',[dimR dimZ]);
varidTe = netcdf.defVar(ncid2,'te','double',[dimR dimZ]);
varidVr = netcdf.defVar(ncid2,'vr','double',[dimR dimZ]);
varidVz = netcdf.defVar(ncid2,'vz','double',[dimR dimZ]);
varidVt = netcdf.defVar(ncid2,'vt','double',[dimR dimZ]);
varidGi = netcdf.defVar(ncid2,'gradTi','double',[dimR dimZ]);
% optional: store mask threshold
varidPsi = netcdf.defVar(ncid2,'psiN_threshold','double',[]);

netcdf.endDef(ncid2);

% write coords
netcdf.putVar(ncid2,varidR, X(1,:));    % x = columns of X
netcdf.putVar(ncid2,varidZ, Y(:,1));    % z = rows of Y

% write meta
netcdf.putVar(ncid2,varidZnum, Z_all(:));
netcdf.putVar(ncid2,varidQnum, q_all(:));

% write multi-species masked cubes
netcdf.putVar(ncid2,varid_ni, val_ni_mspecies_masked);
netcdf.putVar(ncid2,varid_ti, val_Ti_mspecies_masked);
netcdf.putVar(ncid2,varid_uR, val_uR_mspecies_masked);
netcdf.putVar(ncid2,varid_uZ, val_uZ_mspecies_masked);
netcdf.putVar(ncid2,varid_uT, val_uT_mspecies_masked);

% write bulk masked fields you already computed
netcdf.putVar(ncid2,varidNe, val_ne_masked);
netcdf.putVar(ncid2,varidTe, val_Te_masked);
netcdf.putVar(ncid2,varidVr, val_vr_masked);
netcdf.putVar(ncid2,varidVz, val_vz_masked);
netcdf.putVar(ncid2,varidVt, val_vt_masked);
netcdf.putVar(ncid2,varidGi, val_gradTi_masked);
netcdf.putVar(ncid2,varidPsi, psiN_mask);

netcdf.close(ncid2);
fprintf('Wrote multispecies file: %s\n', outfile_ms);

%% 7b) (Optional) quick check: plot carbon (Z=6) n_i after masking
isC = find(Z_all==6, 1, 'first');
if ~isempty(isC)
    figure;
    imagesc([min(X(:)) max(X(:))],[min(Y(:)) max(Y(:))], val_ni_mspecies_masked(:,:,isC));
    set(gca,'YDir','normal'); axis equal tight; colorbar;
    title('Masked Carbon n_i on (X,Y)');
    hold on; plot(g.lim(1,:), g.lim(2,:), 'k'); contour(X,Y,psiN,[psiN_mask psiN_mask],'k--');
end

%% 7) Read back & visualize some profile slices
disp('Reading back masked profiles...');
Rr = ncread('profilesDIIID_196154.nc','x');
Zz = ncread('profilesDIIID_196154.nc','z');

ne_m = ncread('profilesDIIID_196154.nc','ne');
te_m = ncread('profilesDIIID_196154.nc','te');

% 2D density
figure; imagesc(Zz, Rr, ne_m);
set(gca,'YDir','normal','FontName','Times','FontSize',24);
xlabel('$z$ [m]','Interpreter','latex');
ylabel('$r$ [m]','Interpreter','latex');
cb = colorbar; ylabel(cb,'$n_e\,[\mathrm{m}^{-3}]$','Interpreter','latex');
title('Masked $n_e$ (\psi_N<0.8)','Interpreter','latex');
pbaspect([2 1 1]);

% 1D slices at inner radius and mid-plane
figure; plot(Zz, ne_m(2,:), 'LineWidth',1.5);
xlabel('$z$ [m]','Interpreter','latex');
ylabel('$n_e$','Interpreter','latex');
title('Axial $n_e$ at innermost radius','Interpreter','latex');

figure; plot(Rr, ne_m(:,round(nZ/2)), 'LineWidth',1.5);
xlabel('$r$ [m]','Interpreter','latex');
ylabel('$n_e$','Interpreter','latex');
title('Radial $n_e$ at mid-plane','Interpreter','latex');