clc; clear all; close all;
outnc   = fullfile(pwd, 'profiles_iter_wideGrid.nc'); 
%% === VISUALIZATION CHECK: Read-back from GITR-style NetCDF ===
fprintf('\n=== Visualizing data from %s ===\n', outnc);

% --- Read back main fields ---
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

% --- Derived quantities ---
Bmag_chk = sqrt(br_chk.^2 + bt_chk.^2 + bz_chk.^2);
Umag_chk = sqrt(uR_chk.^2 + uZ_chk.^2 + uT_chk.^2);

% === Basic 2D contour sanity plots ===
figure('Name','n_e from GITR NetCDF');
imagesc(Rcheck, Zcheck, ne_chk'); set(gca,'YDir','normal');
xlabel('R [m]'); ylabel('Z [m]'); title('Electron density n_e [m^{-3}]');
axis equal tight; colorbar; colormap(turbo);

figure('Name','|B| from GITR NetCDF');
imagesc(Rcheck, Zcheck, Bmag_chk'); set(gca,'YDir','normal');
xlabel('R [m]'); ylabel('Z [m]'); title('|B| [T]');
axis equal tight; colorbar; colormap(turbo);

figure('Name','n_e from GITR NetCDF');
imagesc(Rcheck, Zcheck, ne_chk'); set(gca,'YDir','normal');
xlabel('R [m]'); ylabel('Z [m]'); title('Electron density n_e [m^{-3}]');
axis equal tight; colorbar; colormap(turbo);

figure('Name','U_T from GITR NetCDF');
imagesc(Rcheck, Zcheck, uT_chk(:,:,2)'); set(gca,'YDir','normal');
xlabel('R [m]'); ylabel('Z [m]'); title('|B| [T]');
axis equal tight; colorbar; colormap(turbo);

figure('Name','U_Z from GITR NetCDF');
imagesc(Rcheck, Zcheck, uZ_chk(:,:,2)'); set(gca,'YDir','normal');
xlabel('R [m]'); ylabel('Z [m]'); title('|B| [T]');
axis equal tight; colorbar; colormap(turbo);

figure('Name','U_R from GITR NetCDF');
imagesc(Rcheck, Zcheck, uR_chk(:,:,2)'); set(gca,'YDir','normal');
xlabel('R [m]'); ylabel('Z [m]'); title('|B| [T]');
axis equal tight; colorbar; colormap(turbo);



fprintf('✅ GITR NetCDF visualization complete. Check figures for consistency.\n');