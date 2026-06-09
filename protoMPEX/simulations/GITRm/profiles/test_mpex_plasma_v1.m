%% test_mpex_plasma_v1.m
% ProtoMPEX / MPEX toy plasma + Ta profiles for impurity transport challenge
% - Geometry: 10 cm x 10 cm x 20 cm box, plasma radius a = 2 cm
% - Uniform axial Bz = 0.5 T, Br = Bt = 0
% - Super-Gaussian profiles in r: exp(-(r/a)^mOrder)
% - Te = Ti, T_Ta = 5 eV, v_||,Ta = 1000 m/s
% - Ta flux + density + temperature + velocities included
%
% Outputs:
%   1) profilesProtoMPEX.nc  : ne, ni, te, ti, gradTi, gradTe, vr, vt, vz,
%                              + Ta density/temp/velocities and Ta fluxes
%   2) BfieldProtoMPEX.nc    : x, z, br, bt, bz

clear; clc; close all;

%% ----------------- Basic geometry & grids -----------------
Lz_box_m  = 0.20;        % [m] 20 cm
R_box_m   = 0.05;        % [m] 5 cm (half-width)
a_m       = 0.02;        % [m] plasma radius (2 cm)
mOrder    = 12;          % super-Gaussian order

nR = 301;
nZ = 401;
x  = linspace(0, R_box_m, nR);   % radial coordinate
z  = linspace(0, Lz_box_m, nZ);  % axial coordinate

[XZ, R] = meshgrid(z, x);        % R(i,j)=r_i, Z(i,j)=z_j

SG1D = exp(-(x./a_m).^mOrder);   % 1D shape in r
SG2D = exp(-(R./a_m).^mOrder);   % 2D shape in (r,z)

%% ----------------- Plasma profiles -----------------
ne_center = 1e19;        % [m^-3]
Te_center = 5.0;         % [eV]  (Te = Ti = 5 eV on axis)

SG2D_norm = SG2D / max(SG2D(:));

ne0 = ne_center * SG2D;        % [m^-3]
ni0 = ne0;                     % [m^-3] quasi-neutral

% Enforce Te = Ti everywhere (same profile)
te0 = Te_center * SG2D_norm;   % [eV]
ti0 = te0;                     % [eV] exactly equal to Te

% Radial gradients
dr = x(2) - x(1);
gradTi0 = zeros(nR, nZ);
gradTe0 = zeros(nR, nZ);
gradTi0(1:end-1,:) = diff(ti0,1,1)/dr;
gradTe0(1:end-1,:) = diff(te0,1,1)/dr;
gradTi0(end,:)     = gradTi0(end-1,:);
gradTe0(end,:)     = gradTe0(end-1,:);

%% ----------------- Flow and B-field -----------------
Bz_val = 0.5;     % [T]
Vr_val = 0.0;     % [m/s]
Vt_val = 0.0;     % [m/s]
Vz_val = 0.0;     % [m/s]

vr0 = Vr_val * zeros(nR,nZ);
vt0 = Vt_val * zeros(nR,nZ);
vz0 = Vz_val * ones(nR,nZ);    % uniform (here 0)

br0 = zeros(nR,nZ);
bt0 = zeros(nR,nZ);
bz0 = Bz_val * ones(nR,nZ);

%% ----------------- Tantalum profiles -----------------
% Charge-state fractions
fTa2 = 0.08; fTa3 = 0.62; fTa4 = 0.30;
S_f  = fTa2 + fTa3 + fTa4;
fTa2 = fTa2 / S_f; fTa3 = fTa3 / S_f; fTa4 = fTa4 / S_f;

% Total Ta flux (given) and charge-state split
GammaTa0_1D   = 1e20 * SG1D;                 % [m^-2 s^-1]
GammaTa_tot0  = repmat(GammaTa0_1D(:),1,nZ); % [nR x nZ]
GammaTa2_0    = fTa2 * GammaTa_tot0;
GammaTa3_0    = fTa3 * GammaTa_tot0;
GammaTa4_0    = fTa4 * GammaTa_tot0;

% Enforce v_||,Ta = 1000 m/s
VTa_par = 1000.0;                             % [m/s] *** fixed ***
nTa_1D  = GammaTa0_1D ./ VTa_par;             % [m^-3]
nTa_2D  = repmat(nTa_1D(:),1,nZ);

% Enforce T_Ta = 5 eV everywhere (constant)
T_Ta_const = 5.0;                              % [eV] *** fixed ***
tTa_2D     = T_Ta_const * ones(nR,nZ);         % flat in r,z

% Tantalum velocity components
vrTa_2D = zeros(nR,nZ);
vtTa_2D = zeros(nR,nZ);
vzTa_2D = VTa_par * ones(nR,nZ);               % uniform 1000 m/s in +z

%% ----------------- 1D visualizations (pre-NC) -----------------
r_cm = x*100; z_cm = z*100;
midZ = round(nZ/2);

figure(10); clf;
subplot(2,3,1);
plot(r_cm, ne0(:,midZ)/1e19,'LineWidth',1.6); grid on;
xlabel('r [cm]'); ylabel('n_e [10^{19} m^{-3}]'); title('n_e(r) midplane');

subplot(2,3,2);
plot(r_cm, te0(:,midZ),'LineWidth',1.6); grid on;
xlabel('r [cm]'); ylabel('T_e [eV]'); title('T_e(r) midplane');

subplot(2,3,3);
plot(r_cm, ti0(:,midZ),'LineWidth',1.6); grid on;
xlabel('r [cm]'); ylabel('T_i [eV]'); title('T_i(r) midplane (Te=Ti)');

subplot(2,3,4);
plot(r_cm, gradTe0(:,midZ),'LineWidth',1.6); grid on;
xlabel('r [cm]'); ylabel('dT_e/dr [eV/m]'); title('gradT_e(r) midplane');

subplot(2,3,5);
plot(r_cm, gradTi0(:,midZ),'LineWidth',1.6); grid on;
xlabel('r [cm]'); ylabel('dT_i/dr [eV/m]'); title('gradT_i(r) midplane');

subplot(2,3,6);
plot(r_cm, vz0(:,midZ),'LineWidth',1.6); grid on;
xlabel('r [cm]'); ylabel('v_z [m/s]'); title('v_z(r) midplane');

sgtitle('Plasma 1D profiles (pre-NetCDF)');

figure(11); clf;
subplot(2,2,1);
plot(r_cm, GammaTa0_1D,'LineWidth',1.6); grid on;
xlabel('r [cm]'); ylabel('\Gamma_{Ta,tot} [m^{-2}s^{-1}]');
title('\Gamma_{Ta,tot}(r)');

subplot(2,2,2);
plot(r_cm, nTa_1D,'LineWidth',1.6); grid on;
xlabel('r [cm]'); ylabel('n_{Ta} [m^{-3}]');
title('n_{Ta}(r)');

subplot(2,2,3);
plot(r_cm, tTa_2D(:,midZ),'LineWidth',1.6); grid on;
xlabel('r [cm]'); ylabel('T_{Ta} [eV]');
title('T_{Ta}(r) midplane (5 eV flat)');

subplot(2,2,4);
plot(r_cm, vzTa_2D(:,midZ),'LineWidth',1.6); grid on;
xlabel('r [cm]'); ylabel('v_{z,Ta} [m/s]');
title('v_{z,Ta}(r) midplane (1000 m/s flat)');

sgtitle('Ta 1D profiles (pre-NetCDF)');

%% ----------------- 2D visualizations (pre-NC) -----------------
figure(20); clf;
subplot(2,3,1);
imagesc(z_cm, r_cm, ne0/1e19); set(gca,'YDir','normal'); colorbar;
xlabel('z [cm]'); ylabel('r [cm]'); title('n_e [10^{19} m^{-3}]');

subplot(2,3,2);
imagesc(z_cm, r_cm, te0); set(gca,'YDir','normal'); colorbar;
xlabel('z [cm]'); ylabel('r [cm]'); title('T_e [eV]');

subplot(2,3,3);
imagesc(z_cm, r_cm, ti0); set(gca,'YDir','normal'); colorbar;
xlabel('z [cm]'); ylabel('r [cm]'); title('T_i [eV] (Te=Ti)');

subplot(2,3,4);
imagesc(z_cm, r_cm, gradTe0); set(gca,'YDir','normal'); colorbar;
xlabel('z [cm]'); ylabel('r [cm]'); title('gradT_e [eV/m]');

subplot(2,3,5);
imagesc(z_cm, r_cm, gradTi0); set(gca,'YDir','normal'); colorbar;
xlabel('z [cm]'); ylabel('r [cm]'); title('gradT_i [eV/m]');

subplot(2,3,6);
imagesc(z_cm, r_cm, vz0); set(gca,'YDir','normal'); colorbar;
xlabel('z [cm]'); ylabel('r [cm]'); title('v_z [m/s]');

sgtitle('Plasma 2D profiles (pre-NetCDF)');

figure(21); clf;
subplot(2,3,1);
imagesc(z_cm, r_cm, GammaTa_tot0); set(gca,'YDir','normal'); colorbar;
xlabel('z [cm]'); ylabel('r [cm]'); title('\Gamma_{Ta,tot}');

subplot(2,3,2);
imagesc(z_cm, r_cm, GammaTa2_0); set(gca,'YDir','normal'); colorbar;
xlabel('z [cm]'); ylabel('r [cm]'); title('\Gamma_{Ta^{2+}}');

subplot(2,3,3);
imagesc(z_cm, r_cm, GammaTa3_0); set(gca,'YDir','normal'); colorbar;
xlabel('z [cm]'); ylabel('r [cm]'); title('\Gamma_{Ta^{3+}}');

subplot(2,3,4);
imagesc(z_cm, r_cm, GammaTa4_0); set(gca,'YDir','normal'); colorbar;
xlabel('z [cm]'); ylabel('r [cm]'); title('\Gamma_{Ta^{4+}}');

subplot(2,3,5);
imagesc(z_cm, r_cm, nTa_2D); set(gca,'YDir','normal'); colorbar;
xlabel('z [cm]'); ylabel('r [cm]'); title('n_{Ta}');

subplot(2,3,6);
imagesc(z_cm, r_cm, tTa_2D); set(gca,'YDir','normal'); colorbar;
xlabel('z [cm]'); ylabel('r [cm]'); title('T_{Ta} [eV] (5 eV flat)');

sgtitle('Ta 2D profiles (pre-NetCDF)');

%% ============================================================
%   WRITE profilesProtoMPEX.nc (plasma vars + Ta)
%% ============================================================
ncfile_profiles = 'profilesProtoMPEX.nc';
if exist(ncfile_profiles,'file'), delete(ncfile_profiles); end

ncid = netcdf.create(ncfile_profiles, ...
    bitor(netcdf.getConstant('NETCDF4'), netcdf.getConstant('CLOBBER')));

% Dimensions
dimR = netcdf.defDim(ncid,'nX', nR);
dimZ = netcdf.defDim(ncid,'nZ', nZ);

% Coordinates
varidR = netcdf.defVar(ncid,'x','double',dimR);
varidZ = netcdf.defVar(ncid,'z','double',dimZ);

% Plasma profiles (DIII-D-style)
varidNe     = netcdf.defVar(ncid,'ne',     'double',[dimR dimZ]);
varidNi     = netcdf.defVar(ncid,'ni',     'double',[dimR dimZ]);
varidTe     = netcdf.defVar(ncid,'te',     'double',[dimR dimZ]);
varidTi     = netcdf.defVar(ncid,'ti',     'double',[dimR dimZ]);
varidGradTi = netcdf.defVar(ncid,'gradTi', 'double',[dimR dimZ]);
varidGradTe = netcdf.defVar(ncid,'gradTe', 'double',[dimR dimZ]);
varidVr     = netcdf.defVar(ncid,'vr',     'double',[dimR dimZ]);
varidVt     = netcdf.defVar(ncid,'vt',     'double',[dimR dimZ]);
varidVz     = netcdf.defVar(ncid,'vz',     'double',[dimR dimZ]);

% Tantalum fields
varidNTa    = netcdf.defVar(ncid,'nTa',    'double',[dimR dimZ]);
varidTTa    = netcdf.defVar(ncid,'tTa',    'double',[dimR dimZ]);
varidVrTa   = netcdf.defVar(ncid,'vrTa',   'double',[dimR dimZ]);
varidVtTa   = netcdf.defVar(ncid,'vtTa',   'double',[dimR dimZ]);
varidVzTa   = netcdf.defVar(ncid,'vzTa',   'double',[dimR dimZ]);

varidGamTaTot = netcdf.defVar(ncid,'GammaTa_tot','double',[dimR dimZ]);
varidGamTa2   = netcdf.defVar(ncid,'GammaTa2',   'double',[dimR dimZ]);
varidGamTa3   = netcdf.defVar(ncid,'GammaTa3',   'double',[dimR dimZ]);
varidGamTa4   = netcdf.defVar(ncid,'GammaTa4',   'double',[dimR dimZ]);

netcdf.endDef(ncid);

% Write coordinates
netcdf.putVar(ncid,varidR,x);
netcdf.putVar(ncid,varidZ,z);

% Write plasma
netcdf.putVar(ncid,varidNe,     ne0);
netcdf.putVar(ncid,varidNi,     ni0);
netcdf.putVar(ncid,varidTe,     te0);
netcdf.putVar(ncid,varidTi,     ti0);
netcdf.putVar(ncid,varidGradTi, gradTi0);
netcdf.putVar(ncid,varidGradTe, gradTe0);
netcdf.putVar(ncid,varidVr,     vr0);
netcdf.putVar(ncid,varidVt,     vt0);
netcdf.putVar(ncid,varidVz,     vz0);

% Write Tantalum
netcdf.putVar(ncid,varidNTa,    nTa_2D);
netcdf.putVar(ncid,varidTTa,    tTa_2D);
netcdf.putVar(ncid,varidVrTa,   vrTa_2D);
netcdf.putVar(ncid,varidVtTa,   vtTa_2D);
netcdf.putVar(ncid,varidVzTa,   vzTa_2D);

netcdf.putVar(ncid,varidGamTaTot, GammaTa_tot0);
netcdf.putVar(ncid,varidGamTa2,   GammaTa2_0);
netcdf.putVar(ncid,varidGamTa3,   GammaTa3_0);
netcdf.putVar(ncid,varidGamTa4,   GammaTa4_0);

netcdf.close(ncid);
fprintf('Wrote plasma + Ta profiles to %s\n', ncfile_profiles);

%% ======================================
%   WRITE BfieldProtoMPEX.nc (B only)
%% ======================================
ncfile_B = 'BfieldProtoMPEX.nc';
if exist(ncfile_B,'file'), delete(ncfile_B); end

ncidB = netcdf.create(ncfile_B, ...
    bitor(netcdf.getConstant('NETCDF4'), netcdf.getConstant('CLOBBER')));

dimRB = netcdf.defDim(ncidB,'nX', nR);
dimZB = netcdf.defDim(ncidB,'nZ', nZ);

varidRB = netcdf.defVar(ncidB,'x','double',dimRB);
varidZB = netcdf.defVar(ncidB,'z','double',dimZB);

varidBr = netcdf.defVar(ncidB,'br','double',[dimRB dimZB]);
varidBt = netcdf.defVar(ncidB,'bt','double',[dimRB dimZB]);
varidBz = netcdf.defVar(ncidB,'bz','double',[dimRB dimZB]);

netcdf.endDef(ncidB);

netcdf.putVar(ncidB,varidRB,x);
netcdf.putVar(ncidB,varidZB,z);
netcdf.putVar(ncidB,varidBr,br0);
netcdf.putVar(ncidB,varidBt,bt0);
netcdf.putVar(ncidB,varidBz,bz0);

netcdf.close(ncidB);
fprintf('Wrote B-field to %s\n', ncfile_B);

%% ----------------- Post-write cross-check -----------------
% Read back from profiles file
x_nc       = ncread(ncfile_profiles,'x');
z_nc       = ncread(ncfile_profiles,'z');
ne_nc      = ncread(ncfile_profiles,'ne');
te_nc      = ncread(ncfile_profiles,'te');
ti_nc      = ncread(ncfile_profiles,'ti');
gradTe_nc  = ncread(ncfile_profiles,'gradTe');
nTa_nc     = ncread(ncfile_profiles,'nTa');
tTa_nc     = ncread(ncfile_profiles,'tTa');
GammaTa_nc = ncread(ncfile_profiles,'GammaTa_tot');

% Read back from Bfile
bz_nc = ncread(ncfile_B,'bz');

x_cm_nc = x_nc*100; z_cm_nc = z_nc*100;
midZ_nc = round(numel(z_nc)/2);

figure(30); clf;
subplot(2,2,1);
plot(x_cm_nc, te0(:,midZ),'k-', x_cm_nc, ti0(:,midZ),'b-.','LineWidth',1.4);
grid on; xlabel('r [cm]'); ylabel('T [eV]');
legend('Te orig','Ti orig'); title('Check Te = Ti (orig)');

subplot(2,2,2);
plot(x_cm_nc, te_nc(:,midZ_nc),'r--', x_cm_nc, ti_nc(:,midZ_nc),'g:','LineWidth',1.4);
grid on; xlabel('r [cm]'); ylabel('T [eV]');
legend('Te nc','Ti nc'); title('Check Te = Ti (from nc)');

subplot(2,2,3);
plot(x_cm_nc, tTa_2D(:,midZ),'k-', x_cm_nc, tTa_nc(:,midZ_nc),'r--','LineWidth',1.4);
grid on; xlabel('r [cm]'); ylabel('T_{Ta} [eV]');
legend('orig','nc'); title('T_{Ta}(r) (should be 5 eV flat)');

subplot(2,2,4);
plot(x_cm_nc, nTa_2D(:,midZ),'k-', x_cm_nc, nTa_nc(:,midZ_nc),'r--','LineWidth',1.4);
grid on; xlabel('r [cm]'); ylabel('n_{Ta} [m^{-3}]');
legend('orig','nc'); title('n_{Ta}(r)');

sgtitle('Te=Ti and Ta checks (orig vs nc)');

figure(31); clf;
subplot(1,2,1);
imagesc(z_cm, x_cm_nc, bz0); set(gca,'YDir','normal'); colorbar;
xlabel('z [cm]'); ylabel('r [cm]'); title('B_z original');

subplot(1,2,2);
imagesc(z_cm_nc, x_cm_nc, bz_nc); set(gca,'YDir','normal'); colorbar;
xlabel('z [cm]'); ylabel('r [cm]'); title('B_z from BfieldProtoMPEX.nc');

sgtitle('B-field cross-check (separate NC)');