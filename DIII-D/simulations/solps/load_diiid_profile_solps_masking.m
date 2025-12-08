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
load('extrapolated_data_196154.mat');  
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