%% SOLPS-ITER → Cartesian grid → NetCDF export (multi-species D + Ne)
% Adds EFIT B-field, interpolates to the SOLPS (R,Z) grid, and decomposes
% species parallel flow Ua along (bR, bZ, bT) unit vectors → (uR,uZ,uT).
% -------------------------------------------------------------------------

% clear; clc; close all;

%% === INPUTS ===
matFile = 'solps_iter.mat';                 % SOLPS input
eqdskFile = 'MOB-348s_eqdsk.txt'; 

% eqdskFile= 'g900003.00230_ITER_15MA_eqdsk16VVHR.txt';% for read_EFIT_data(...)
outnc   = fullfile(pwd, 'profiles_iter_wideGrid.nc');     % NetCDF output

if exist(outnc,'file'), delete(outnc); end
fprintf('Writing multiSpecies NetCDF to %s\n', outnc);

%% === LOAD SOLPS DATA ===
S = load(matFile);
Geo   = S.Geo;
State = S.State;

% geometry (cell corners)
r = double(Geo.pr);   % [4 x Nc]
z = double(Geo.pz);   % [4 x Nc]
Nc = size(r,2);

% cell-centered plasma
ne_cv = double(State.ne(:));                % [Nc x 1]
na_cv = double(State.na);                   % [Nc x ns]
ua_cv = double(State.ua);                   % [Nc x ns] (parallel flow U||)
% Convert J→eV if needed (keep your original conversion)
te_cv = 0.25.*double(State.te(:))./1.602e-19;     % [Nc x 1]
ti_cv = 0.25.*double(State.ti(:))./1.602e-19;     % [Nc x 1]

zn = double(State.zn(:));                   % [ns x 1] atomic number
am = double(State.am(:));                   % [ns x 1] atomic mass
ns = size(na_cv,2);

fprintf('Loaded %d cells and %d species.\n', Nc, ns);

%% === DEFINE CARTESIAN GRID (SOLPS target grid) ===
rmin = min(r(:)); rmax = max(r(:));
zmin = min(z(:)); zmax = max(z(:));
nR = 400; nZ = 800;                              % keep your current res
rgrid = linspace(rmin, rmax, nR);
zgrid = linspace(zmin, zmax, nZ);
[X, Z] = meshgrid(rgrid, zgrid);                 % X→R, Z→Z

% mapper: paint polygon cells onto Cartesian grid
fill_by_cells = @(vals) paintCells(vals, r, z, X, Z);

%% === MAP SOLPS FIELDS TO (R,Z) ===
disp('Interpolating SOLPS fields to (R,Z) grid ...');
ne_q = fill_by_cells(ne_cv);
Te_q = fill_by_cells(te_cv);
Ti_q = fill_by_cells(ti_cv);

ni_q = cell(1, ns);
Ua_q = cell(1, ns);
for s = 1:ns
    ni_q{s} = fill_by_cells(na_cv(:,s));    % [nZ x nR]
    Ua_q{s} = fill_by_cells(ua_cv(:,s));    % [nZ x nR]
end

%% === EFIT B-FIELD → interpolate to SOLPS (R,Z) grid ===
fprintf('Reading EFIT and computing B-field on SOLPS grid ...\n');

% This script defines Br, Bt, Bz, r_efit, z_efit
read_efit_data;   % your EFIT routine

% Interpolants on EFIT grid (EFIT arrays are [nZb x nRb])
FBr = griddedInterpolant({z_efit, r_efit}, Br, 'linear', 'nearest');
FBt = griddedInterpolant({z_efit, r_efit}, Bt, 'linear', 'nearest');
FBz = griddedInterpolant({z_efit, r_efit}, Bz, 'linear', 'nearest');

% Evaluate directly on the SOLPS (Z,R) grid (Z:800, R:400)
Br_q = FBr(Z, X);      % size: [800 x 400]
Bt_q = FBt(Z, X);
Bz_q = FBz(Z, X);

% Normalize to unit vectors
Bmag = sqrt(Br_q.^2 + Bt_q.^2 + Bz_q.^2);
epsB = 1e-30;
bRhat = Br_q ./ max(Bmag, epsB);
bZhat = Bz_q ./ max(Bmag, epsB);
bThat = Bt_q ./ max(Bmag, epsB);

% --- Use SOLPS grid sizes here ---
nZ = size(Z,1); nR = size(X,2);
uR_all = zeros(nZ, nR, ns);
uZ_all = zeros(nZ, nR, ns);
uT_all = zeros(nZ, nR, ns);

for s = 1:ns
    Upar = Ua_q{s};                     % [800 x 400] = [nZ x nR]
    uR_all(:,:,s) = Upar .* bRhat;      % OK: 800x400 .* 800x400
    uZ_all(:,:,s) = Upar .* bZhat;
    uT_all(:,:,s) = Upar .* bThat;
end

%% === DECOMPOSE U_parallel INTO (uR,uZ,uT) FOR EACH SPECIES ===
uR_all = zeros(nZ, nR, ns);
uZ_all = zeros(nZ, nR, ns);
uT_all = zeros(nZ, nR, ns);

for s = 1:ns
    Upar = Ua_q{s};          % [nZ x nR]
    uR_all(:,:,s) = Upar .* bRhat;
    uZ_all(:,:,s) = Upar .* bZhat;
    uT_all(:,:,s) = Upar .* bThat;
end

%% === STACK ARRAYS FOR NETCDF ===
ni_all = zeros(nZ, nR, ns);
Ti_all = zeros(nZ, nR, ns);
for s = 1:ns
    ni_all(:,:,s) = ni_q{s};
    Ti_all(:,:,s) = Ti_q;    % same Ti for all species (keep your choice)
end

%% === WRITE MULTISPECIES NETCDF (now with B and uR/uZ/uT) ===
disp('Writing multiSpecies NetCDF...');
ncid = netcdf.create(outnc, bitor(netcdf.getConstant('NETCDF4'), netcdf.getConstant('CLOBBER')));

% dims
dimR = netcdf.defDim(ncid, 'nX', nR);
dimZ = netcdf.defDim(ncid, 'nZ', nZ);
dimS = netcdf.defDim(ncid, 'species', ns);

% coords/meta
vid_x   = netcdf.defVar(ncid,'x','double',dimR);
vid_z   = netcdf.defVar(ncid,'z','double',dimZ);
vid_Z   = netcdf.defVar(ncid,'atomic_number','double',dimS);
vid_q   = netcdf.defVar(ncid,'charge_number','double',dimS);

% single-fluid 2D
vid_ne  = netcdf.defVar(ncid,'ne','double',[dimR dimZ]);
vid_te  = netcdf.defVar(ncid,'te','double',[dimR dimZ]);
vid_ti  = netcdf.defVar(ncid,'ti','double',[dimR dimZ]);

% B-field on grid
vid_br  = netcdf.defVar(ncid,'br','double',[dimR dimZ]);
vid_bt  = netcdf.defVar(ncid,'bt','double',[dimR dimZ]);
vid_bz  = netcdf.defVar(ncid,'bz','double',[dimR dimZ]);

% multispecies 3D
vid_niA = netcdf.defVar(ncid,'ni_all','double',[dimR dimZ dimS]);
vid_uRA = netcdf.defVar(ncid,'uR_all','double',[dimR dimZ dimS]);
vid_uZA = netcdf.defVar(ncid,'uZ_all','double',[dimR dimZ dimS]);
vid_uTA = netcdf.defVar(ncid,'uT_all','double',[dimR dimZ dimS]);

netcdf.endDef(ncid);

% write coords/meta
netcdf.putVar(ncid, vid_x, rgrid);
netcdf.putVar(ncid, vid_z, zgrid);
netcdf.putVar(ncid, vid_Z, zn);
netcdf.putVar(ncid, vid_q, zn);    % if you have explicit charge numbers, swap here

% write 2D (transpose to [nR x nZ])
netcdf.putVar(ncid, vid_ne, ne_q');
netcdf.putVar(ncid, vid_te, Te_q');
netcdf.putVar(ncid, vid_ti, Ti_q');

% write B (transpose to [nR x nZ])
netcdf.putVar(ncid, vid_br, Br_q');
netcdf.putVar(ncid, vid_bt, Bt_q');
netcdf.putVar(ncid, vid_bz, Bz_q');

% write 3D ([nZ x nR x ns] → [nR x nZ x ns])
netcdf.putVar(ncid, vid_niA, permute(ni_all, [2 1 3]));
netcdf.putVar(ncid, vid_uRA, permute(uR_all, [2 1 3]));
netcdf.putVar(ncid, vid_uZA, permute(uZ_all, [2 1 3]));
netcdf.putVar(ncid, vid_uTA, permute(uT_all, [2 1 3]));

netcdf.close(ncid);
fprintf('✅ Wrote %s (multiSpecies + B + velocity components)\n', outnc);

%% === QUICK PLOTS (sanity) ===
figure; imagesc(rgrid, zgrid, sqrt(Br_q.^2+Bt_q.^2+Bz_q.^2));
set(gca,'YDir','normal'); axis equal tight; colorbar
title('|B| on SOLPS grid'); xlabel('R [m]'); ylabel('Z [m]');

k = min(2,ns);
figure; imagesc(rgrid, zgrid, uT_all(:,:,2));
set(gca,'YDir','normal'); axis equal tight; colorbar
title(sprintf('u_T (species %d)',k)); xlabel('R [m]'); ylabel('Z [m]');

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
%% === MIDPLANE COMPARISON: SOLPS vs GITR NetCDF ===
fprintf('\n=== Comparing SOLPS and GITR data at midplane (Z=0) ===\n');

% --- Find nearest Z index to midplane ---
[~, iz_SOLPS] = min(abs(zgrid - 0));
[~, iz_GITR ] = min(abs(Zcheck - 0));

fprintf('Closest Z index in SOLPS = %d (Z=%.3f m)\n', iz_SOLPS, zgrid(iz_SOLPS));
fprintf('Closest Z index in GITR  = %d (Z=%.3f m)\n', iz_GITR, Zcheck(iz_GITR));

% --- Extract SOLPS midplane profiles ---
ne_SOLPS_mid  = ne_q(iz_SOLPS, :);
Te_SOLPS_mid  = Te_q(iz_SOLPS, :);
Ti_SOLPS_mid  = Ti_q(iz_SOLPS, :);
Bmag_SOLPS_mid = sqrt(Br_q(iz_SOLPS,:).^2 + Bt_q(iz_SOLPS,:).^2 + Bz_q(iz_SOLPS,:).^2);

% --- Extract GITR midplane profiles (read back from NetCDF) ---
ne_GITR_mid  = ne_chk(:, iz_GITR);
Te_GITR_mid  = te_chk(:, iz_GITR);
Ti_GITR_mid  = ti_chk(:, iz_GITR);
Bmag_GITR_mid = sqrt(br_chk(:,iz_GITR).^2 + bt_chk(:,iz_GITR).^2 + bz_chk(:,iz_GITR).^2);

% --- Optional normalization for plotting clarity ---
ne_SOLPS_mid(ne_SOLPS_mid<=0)=NaN;
ne_GITR_mid(ne_GITR_mid<=0)=NaN;

% --- Plot midplane overlays ---
figure('Color','w','Position',[100 100 1000 700]);
subplot(2,2,1)
plot(rgrid, ne_SOLPS_mid, 'b-', 'LineWidth', 2); hold on;
plot(Rcheck, ne_GITR_mid, 'r--', 'LineWidth', 2);
xlabel('R [m]'); ylabel('n_e [m^{-3}]');
legend('SOLPS','GITR','Location','best');
title('Electron Density at Z=0');

subplot(2,2,2)
plot(rgrid, Te_SOLPS_mid, 'b-', 'LineWidth', 2); hold on;
plot(Rcheck, Te_GITR_mid, 'r--', 'LineWidth', 2);
xlabel('R [m]'); ylabel('T_e [eV]');
legend('SOLPS','GITR','Location','best');
title('Electron Temperature at Z=0');

subplot(2,2,3)
plot(rgrid, Ti_SOLPS_mid, 'b-', 'LineWidth', 2); hold on;
plot(Rcheck, Ti_GITR_mid, 'r--', 'LineWidth', 2);
xlabel('R [m]'); ylabel('T_i [eV]');
legend('SOLPS','GITR','Location','best');
title('Ion Temperature at Z=0');

subplot(2,2,4)
plot(rgrid, Bmag_SOLPS_mid, 'b-', 'LineWidth', 2); hold on;
plot(Rcheck, Bmag_GITR_mid, 'r--', 'LineWidth', 2);
xlabel('R [m]'); ylabel('|B| [T]');
legend('SOLPS','GITR','Location','best');
title('Magnetic Field Magnitude at Z=0');

sgtitle('SOLPS vs GITR Midplane Profiles','FontSize',14);

% --- Quantitative check (mean absolute deviation) ---
fprintf('Mean |Δn_e|:  %.3e\n', nanmean(abs(ne_SOLPS_mid - interp1(Rcheck,ne_GITR_mid,rgrid,'linear','extrap'))));
fprintf('Mean |ΔT_e|:  %.3e\n', nanmean(abs(Te_SOLPS_mid - interp1(Rcheck,Te_GITR_mid,rgrid,'linear','extrap'))));
fprintf('Mean |ΔT_i|:  %.3e\n', nanmean(abs(Ti_SOLPS_mid - interp1(Rcheck,Ti_GITR_mid,rgrid,'linear','extrap'))));
fprintf('Mean |Δ|B||:  %.3e\n', nanmean(abs(Bmag_SOLPS_mid - interp1(Rcheck,Bmag_GITR_mid,rgrid,'linear','extrap'))));

fprintf('✅ Midplane comparison complete.\n');

%% === MIDPLANE FRACTIONAL DENSITIES (He0–He2+, Ne0–Ne10+ + total, vertical layout) ===
fprintf('\n=== Plotting SOLPS midplane n_i/n_e (He0–2+, Ne0–10+ + total, vertical layout) ===\n');

% --- Midplane setup ---
R_focus_min = 8.0;                    % start of right-side (antenna) region [m]
R_focus_max = max(rgrid);             % outer edge
[~, iz_SOLPS] = min(abs(zgrid - 0));  % Z ≈ 0 slice
fprintf('Using SOLPS midplane index %d (Z = %.3f m)\n', iz_SOLPS, zgrid(iz_SOLPS));

% --- Electron density along midplane ---
ne_mid = ne_q(iz_SOLPS,:);  
ne_mid(ne_mid <= 0) = NaN;            % avoid division by zero

% --- Define species indices and names ---
idxHe = [3 4 5];        % He0, He1+, He2+
HeNames = {'He^{0}','He^{1+}','He^{2+}'};
idxNe = 6:16;           % Ne0, Ne1+ ... Ne10+
NeNames = {'Ne^{0}','Ne^{1+}','Ne^{2+}','Ne^{3+}','Ne^{4+}', ...
           'Ne^{5+}','Ne^{6+}','Ne^{7+}','Ne^{8+}','Ne^{9+}','Ne^{10+}'};

% --- Limit to right-hand region ---
rightMask = rgrid >= R_focus_min & rgrid <= R_focus_max;
Rzoom = rgrid(rightMask);

% --- Figure: two stacked panels (He top, Ne bottom) ---
figure('Color','w','Position',[100 100 1100 900]);
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

%% === Helium fractions ===
nexttile(1);
hold on; box on;
cmap = jet(length(idxHe));
for j = 1:length(idxHe)
    i = idxHe(j);
    if i > numel(ni_q), continue; end
    n_i = ni_q{i}(iz_SOLPS,:);
    frac = n_i ./ ne_mid;            % n_i / n_e
    plot(Rzoom, frac(rightMask), 'LineWidth', 2, 'Color', cmap(j,:), ...
         'DisplayName', HeNames{j});
end
xlabel('R [m]');
ylabel('n_i / n_e');
title('He charge-state fractions (SOLPS midplane)');
legend show; legend boxoff;
colormap(jet); grid on; axis tight;

% === Neon fractions + total ===
nexttile(2);
hold on; box on;
cmap = jet(length(idxNe));
nNe_total = zeros(size(ne_mid));

for j = 1:length(idxNe)
    i = idxNe(j);
    if i > numel(ni_q), continue; end
    n_i = ni_q{i}(iz_SOLPS,:);
    nNe_total = nNe_total + n_i;      % accumulate for total
    frac = n_i ./ ne_mid;
    plot(Rzoom, frac(rightMask), 'LineWidth', 1.5, 'Color', cmap(j,:), ...
         'DisplayName', NeNames{j});
end

% --- Add total Ne / n_e ---
frac_total = nNe_total ./ ne_mid;
plot(Rzoom, frac_total(rightMask), 'k-', 'LineWidth', 2.5, ...
     'DisplayName', 'Total Ne/n_e');

xlabel('R [m]');
ylabel('n_i / n_e');
title('Ne charge-state fractions (SOLPS midplane)');
legend show; legend boxoff;
colormap(jet); grid on; axis tight;

%% === MIDPLANE ELECTRON DENSITY AND TEMPERATURE ===
fprintf('\n=== Plotting SOLPS midplane n_e and T_e ===\n');

% --- Choose Z ≈ 0 midplane slice ---
[~, iz_SOLPS] = min(abs(zgrid - 0));      % index for Z=0
fprintf('Using midplane index %d (Z = %.3f m)\n', iz_SOLPS, zgrid(iz_SOLPS));

% --- Extract midplane data ---
ne_mid = ne_q(iz_SOLPS,:);                 % electron density [m^-3]
Te_mid = Te_q(iz_SOLPS,:);                 % electron temperature [eV]
ne_mid(ne_mid <= 0) = NaN;                 % mask invalid
Te_mid(Te_mid <= 0) = NaN;

% --- Restrict to right-hand region (if desired) ---
R_focus_min = 8.0;                         % adjust for your case
R_focus_max = max(rgrid);
maskR = rgrid >= R_focus_min & rgrid <= R_focus_max;
Rzoom = rgrid(maskR);

% --- Reference radii for vertical lines ---
R_line1 = 8.333;    % outer target
R_line2 = 8.173;   % separatrix or inner reference

% --- Plot n_e and T_e vs R ---
figure('Color','w','Position',[100 100 1100 450]);
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

% ==================== Electron density ====================
nexttile(1);
semilogy(Rzoom, ne_mid(maskR), 'b-', 'LineWidth', 2);
xlabel('R [m]');
ylabel('n_e [m^{-3}]');
title('Midplane Electron Density n_e (Z≈0)');
grid on; box on; axis tight;
set(gca,'FontSize',12);

% Add vertical lines
xline(R_line1,'--k','LineWidth',1.3,'DisplayName','R = 8.33 m');
xline(R_line2,':r','LineWidth',1.3,'DisplayName','R = 8.173 m');
legend({'n_e','R=8.33 m','R=8.173 m'},'Location','best','Box','off','FontSize',10);

% ==================== Electron temperature ====================
nexttile(2);
plot(Rzoom, Te_mid(maskR), 'r-', 'LineWidth', 2);
xlabel('R [m]');
ylabel('T_e [eV]');
title('Midplane Electron Temperature T_e (Z≈0)');
grid on; box on; axis tight;
set(gca,'FontSize',12);

% Add vertical lines
xline(R_line1,'--k','LineWidth',1.3,'DisplayName','R = 8.33 m');
xline(R_line2,':b','LineWidth',1.3,'DisplayName','R = 8.173 m');
legend({'T_e','R=8.33 m','R=8.173 m'},'Location','best','Box','off','FontSize',10);

sgtitle('SOLPS Midplane Profiles (Z≈0): n_e and T_e with Reference Radii','FontSize',14);

%% === NEON DENSITY AND CHARGE-STATE FRACTIONS AT R = 8.33 m (MIDPLANE) ===
fprintf('\n=== Calculating Neon charge-state fractions at R = 8.33 m (midplane) ===\n');

R_target = 8.33;                              % [m] chosen location
[~, iR_target] = min(abs(rgrid - R_target));  % find nearest grid index

% --- Define Neon charge states ---
idxNe = 6:16;
NeNames = {'Ne^{0}','Ne^{1+}','Ne^{2+}','Ne^{3+}','Ne^{4+}', ...
           'Ne^{5+}','Ne^{6+}','Ne^{7+}','Ne^{8+}','Ne^{9+}','Ne^{10+}'};

[~, iz_SOLPS] = min(abs(zgrid - 0));          % Z ≈ 0 slice
fprintf('Using midplane Z index %d (Z = %.3f m)\n', iz_SOLPS, zgrid(iz_SOLPS));
fprintf('Nearest R grid index %d (R = %.3f m)\n', iR_target, rgrid(iR_target));

% --- Extract Neon densities at this (Z≈0, R=8.33) point ---
nNe_atR = zeros(length(idxNe),1);
for j = 1:length(idxNe)
    if j > numel(ni_q), continue; end
    n_i = ni_q{idxNe(j)}(iz_SOLPS,:);     % along R
    nNe_atR(j) = n_i(iR_target);          % value at target R
end

% --- Compute total and fractions -----------------------------------------
totalNe = nansum(nNe_atR);
fracNe  = nNe_atR ./ totalNe;

fprintf('\nTotal Neon density at R = %.3f m, Z = 0: %.3e m^-3\n', rgrid(iR_target), totalNe);
fprintf('Charge-state fractions:\n');
for j = 1:length(idxNe)
    fprintf('%10s : %.3e (%.2f%%)\n', NeNames{j}, fracNe(j), fracNe(j)*100);
end

% --- Optional bar plot ---------------------------------------------------
figure('Color','w','Position',[100 100 700 450]);
bar(1:length(idxNe), fracNe*100,'FaceColor','flat');
set(gca,'XTick',1:length(idxNe),'XTickLabel',NeNames,'XTickLabelRotation',45);
ylabel('Fraction [%]');
title(sprintf('Neon charge-state fractions at R = %.2f m (Z≈0)', rgrid(iR_target)));
grid on; box on;

%% === MIDPLANE ELECTRON DENSITY AND TEMPERATURE (Z = 0.0 m) ===
fprintf('\n=== Plotting SOLPS midplane n_e and T_e (Z = %.3f m) ===\n', 0.0);

Z_target = 0.0;    % target height (m)

% --- Find nearest Z index in SOLPS grid -----------------------------------
[~, iz_SOLPS] = min(abs(zgrid - Z_target));
fprintf('Using SOLPS index %d (Z = %.3f m)\n', iz_SOLPS, zgrid(iz_SOLPS));

% --- Extract SOLPS midplane data ------------------------------------------
ne_mid = ne_q(iz_SOLPS,:);                 % electron density [m^-3]
Te_mid = Te_q(iz_SOLPS,:);                 % electron temperature [eV]
ne_mid(ne_mid <= 0) = NaN;
Te_mid(Te_mid <= 0) = NaN;

% --- Restrict to right-hand region (LFS) ----------------------------------
R_focus_min = 8.0;                         
R_focus_max = max(rgrid);
maskR = rgrid >= R_focus_min & rgrid <= R_focus_max;
Rzoom = rgrid(maskR);

%% === Find LFS separatrix crossing at Z ≈ Z_target using ψ_N = 1 contour ===
psiN = (g.psirz - g.ssimag) / (g.ssibry - g.ssimag);
R_sep_LFS = NaN;  Z_sep_LFS = NaN;

if max(psiN(:)) >= 1
    C = contourc(g.r, g.z, psiN', [1 1]);  % ψ_N = 1 contour
    segs = {}; k = 1;
    while k < size(C,2)
        lev = C(1,k); npt = C(2,k);
        if abs(lev - 1) < 1e-9, segs{end+1} = C(:,k+1:k+npt); end %#ok<AGROW>
        k = k + npt + 1;
    end
    if ~isempty(segs)
        % Pick LFS segment (largest R)
        maxR_all = -inf; seg_best = 1;
        for s = 1:numel(segs)
            if max(segs{s}(1,:)) > maxR_all, seg_best = s; end
        end
        Rlcfs = segs{seg_best}(1,:);  Zlcfs = segs{seg_best}(2,:);
        [Zlcfs_sorted, idxU] = unique(Zlcfs,'sorted');
        Rlcfs_sorted = Rlcfs(idxU);
        if any(Zlcfs_sorted >= Z_target) && any(Zlcfs_sorted <= Z_target)
            R_sep_LFS = interp1(Zlcfs_sorted, Rlcfs_sorted, Z_target,'linear','extrap');
            Z_sep_LFS = Z_target;
        else
            [R_sep_LFS, iMax] = max(Rlcfs_sorted);
            Z_sep_LFS = Zlcfs_sorted(iMax);
        end
        fprintf('LFS separatrix at R = %.3f m, Z = %.3f m\n',R_sep_LFS,Z_sep_LFS);
    end
else
    warning('ψ_N never reaches 1 inside domain (max ψ_N = %.3f).', max(psiN(:)));
end

%% --- Find OUTER (LFS) wall intersection at Z ≈ Z_target ------------------
R_wall_LFS = NaN;
if isfield(g,'lim') && ~isempty(g.lim)
    wall_r = g.lim(1,:); wall_z = g.lim(2,:);
    good = isfinite(wall_r) & isfinite(wall_z);
    wall_r = wall_r(good); wall_z = wall_z(good);
    zTol = 0.05;
    nearMask = abs(wall_z - Z_target) <= zTol;
    if any(nearMask)
        R_wall_LFS = max(wall_r(nearMask));
        Z_wall_LFS = mean(wall_z(nearMask));
    else
        [R_wall_LFS, iMaxRw] = max(wall_r);
        Z_wall_LFS = wall_z(iMaxRw);
    end
    fprintf('Outer wall at R = %.3f m, Z ≈ %.3f m\n', R_wall_LFS, Z_wall_LFS);
end

%% === Exponential falloff of n_e, T_e to the wall =========================
validMask = ~isnan(ne_mid);
R_valid = rgrid(validMask);
ne_valid = ne_mid(validMask);
Te_valid = Te_mid(validMask);

if isempty(R_valid)
    warning('No valid n_e points for extrapolation.');
    Rplot = rgrid; ne_plot = ne_mid; Te_plot = Te_mid;
else
    % --- Last grid point with valid data ---
    R_last = max(R_valid);
    ne_last = ne_valid(end);
    Te_last = Te_valid(end);

    % --- Only extrapolate if wall is beyond last grid point ---
    if R_wall_LFS > R_last
        % Define exponential decay lengths
        lambda_ne = 0.017;   % [m] 17 mm
        lambda_Te = 0.017;   % [m] 17 mm

        % Build exponential extension (30 points)
        R_ext = linspace(R_last, R_wall_LFS, 30);
        ne_ext = ne_last .* exp(-(R_ext - R_last)/lambda_ne);
        Te_ext = Te_last .* exp(-(R_ext - R_last)/lambda_Te);

        % Safely select only values up to R_last (avoid index mismatch)
        insideMask = rgrid < R_last + 1e-6;
        R_core = rgrid(insideMask);
        ne_core = ne_mid(insideMask);
        Te_core = Te_mid(insideMask);

        % Merge core + extrapolated region
        Rplot = [R_core, R_ext];
        ne_plot = [ne_core, ne_ext];
        Te_plot = [Te_core, Te_ext];

        fprintf('Exponential falloff applied: λ=%.1f mm from R=%.3f→%.3f m\n', ...
                lambda_ne*1e3, R_last, R_wall_LFS);
    else
        % Wall inside grid → no extrapolation
        Rplot = rgrid;
        ne_plot = ne_mid;
        Te_plot = Te_mid;
    end
end

%% === Plot n_e and T_e with markers ======================================
figure('Color','w','Position',[100 100 1200 700]);
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

% ---- n_e ----
nexttile(1);
semilogy(Rplot, ne_plot, 'b-', 'LineWidth', 2);
xlabel('R [m]'); ylabel('n_e [m^{-3}]');
title(sprintf('n_e at Z = %.3f m (λ_n = 17 mm exponential falloff)', Z_target));
grid on; box on; axis tight; hold on;
set(gca,'FontSize',12);
if ~isnan(R_sep_LFS), xline(R_sep_LFS,'--k','LineWidth',1.5,'DisplayName','Separatrix'); end
if ~isnan(R_wall_LFS), xline(R_wall_LFS,':r','LineWidth',1.5,'DisplayName','Wall'); end
legend show; legend boxoff;

% ---- T_e ----
nexttile(2);
plot(Rplot, Te_plot, 'r-', 'LineWidth', 2);
xlabel('R [m]'); ylabel('T_e [eV]');
title(sprintf('T_e at Z = %.3f m (λ_T = 17 mm exponential falloff)', Z_target));
grid on; box on; axis tight; hold on;
set(gca,'FontSize',12);
if ~isnan(R_sep_LFS), xline(R_sep_LFS,'--k','LineWidth',1.5,'DisplayName','Separatrix'); end
if ~isnan(R_wall_LFS), xline(R_wall_LFS,':b','LineWidth',1.5,'DisplayName','Wall'); end
legend show; legend boxoff;

% ---- Annotation ----
annotation('textbox',[0.5 0.94 0 0],'String',...
    sprintf('Beyond R ≈ %.2f m, the low-density region is extrapolated exponentially (λ = 17 mm)', ...
    R_wall_LFS),'FontSize',11,'HorizontalAlignment','center','EdgeColor','none');

sgtitle(sprintf('SOLPS Midplane Profiles (Z = %.3f m) — Exponential Falloff to Wall', Z_target),'FontSize',14);

%% --- Global title ---
sgtitle(sprintf('SOLPS Midplane (Z≈0) — Fractional Densities  %.1f ≤ R ≤ %.1f m', ...
                R_focus_min, R_focus_max), 'FontSize',14);

fprintf('✅ SOLPS midplane fractional density plots (He0–2+, Ne0–10+ + total, vertical layout) complete.\n');
%% === 2-D SOLPS MAPS (log-scale n_e, T_e; linear fluxes) & (Total Ne, T_i, Average Neon Ua) with EFIT wall + jet ===
fprintf('\n=== Plotting 2-D SOLPS maps (n_e, T_e, flux_e, total Ne, T_i, avg Neon Ua) using jet colormap ===\n');

% --- EFIT wall outline (from read_efit_data or gfile) ---
if exist('rW','var') && exist('zW','var')
    r_wall = rW; z_wall = zW;
elseif exist('g','var') && isfield(g,'lim')
    r_wall = g.lim(1,:); z_wall = g.lim(2,:);
else
    warning('No EFIT wall data found. Skipping wall overlay.');
    r_wall = []; z_wall = [];
end

% --- Neon species indices (Ne0 ... Ne10+) ---
idxNe = 6:16;
fprintf('Averaging Neon parallel flow Ua over indices %d–%d.\n', idxNe(1), idxNe(end));

% --- Compute total Neon density ---
nNe_total = zeros(size(ne_q));
for j = idxNe
    if j > numel(ni_q), continue; end
    nNe_total = nNe_total + ni_q{j};
end

% --- Compute average Neon parallel flow Ua (indices 7–16) ---
Ua_Ne_sum = zeros(size(Ua_q{1}));
count = 0;
for j = 7:16
    if j > numel(Ua_q) || isempty(Ua_q{j}), continue; end
    Ua_Ne_sum = Ua_Ne_sum + Ua_q{j};
    count = count + 1;
end
Ua_Ne_avg = Ua_Ne_sum ./ max(count,1);
Ua_Ne_avg(isnan(Ua_Ne_avg)) = 0;

% --- Electron flow and flux (use Ua_q{2}) ---
fprintf('Using Ua_q{2} for electron parallel flow.\n');
Ua_e = Ua_q{2};
flux_e = ne_q .* abs(Ua_e);     % electron particle flux (n_e * |U_parallel|)

% --- Apply log10 safely to densities and temperatures ---
log_ne = log10(max(ne_q,0));
log_Te = log10(max(Te_q,0));
log_Ti = log10(max(Ti_q,0));
log_nNe_total = log10(max(nNe_total,0));

% --- Figure layout ---
figure('Color','w','Position',[100 100 1500 900]);
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');

% ======== Row 1: Electrons ========
% log10(n_e)
nexttile(1);
imagesc(rgrid,zgrid,log_ne); set(gca,'YDir','normal');
xlabel('R [m]'); ylabel('Z [m]');
title('log_{10}(Electron Density)  log_{10}(n_e [m^{-3}])');
colormap(gca,jet); colorbar; axis equal tight;
if ~isempty(r_wall), hold on; plot(r_wall,z_wall,'k','LineWidth',1.5); end

% log10(T_e)
nexttile(2);
imagesc(rgrid,zgrid,log_Te); set(gca,'YDir','normal');
xlabel('R [m]'); ylabel('Z [m]');
title('log_{10}(Electron Temperature)  log_{10}(T_e [eV])');
colormap(gca,jet); colorbar; axis equal tight;
if ~isempty(r_wall), hold on; plot(r_wall,z_wall,'k','LineWidth',1.5); end

% Electron flux (n_e × |U_parallel|)
nexttile(3);
imagesc(rgrid,zgrid,flux_e); set(gca,'YDir','normal');
xlabel('R [m]'); ylabel('Z [m]');
title('Electron Flux  n_e × |U_{||,e}| [m^{-2}s^{-1}]');
colormap(gca,jet); colorbar; axis equal tight;
if ~isempty(r_wall), hold on; plot(r_wall,z_wall,'k','LineWidth',1.5); end

% ======== Row 2: Total Neon + Average Flow ========
% log10(total Ne density)
nexttile(4);
imagesc(rgrid,zgrid,log_nNe_total); set(gca,'YDir','normal');
xlabel('R [m]'); ylabel('Z [m]');
title('log_{10}(Total Neon Density)  log_{10}(∑ n_{Ne^q} [m^{-3}])');
colormap(gca,jet); colorbar; axis equal tight;
if ~isempty(r_wall), hold on; plot(r_wall,z_wall,'k','LineWidth',1.5); end

% T_i (linear)
nexttile(5);
imagesc(rgrid,zgrid,log_Ti); set(gca,'YDir','normal');
xlabel('R [m]'); ylabel('Z [m]');
title('Ion Temperature T_i [eV]');
colormap(gca,jet); colorbar; axis equal tight;
if ~isempty(r_wall), hold on; plot(r_wall,z_wall,'k','LineWidth',1.5); end

% Average Neon parallel flow ⟨U_a,Ne⟩
nexttile(6);
imagesc(rgrid,zgrid,Ua_Ne_avg); set(gca,'YDir','normal');
xlabel('R [m]'); ylabel('Z [m]');
title('Average Neon Parallel Flow  ⟨U_{a,Ne}⟩ [m/s]  (indices 7–16)');
colormap(gca,jet); colorbar; axis equal tight;
if ~isempty(r_wall), hold on; plot(r_wall,z_wall,'k','LineWidth',1.5); end

sgtitle('SOLPS 2-D Fields (Electrons & Total Neon) with Average Neon Parallel Flow — log_{10}(n, T_e) shown','FontSize',16);
fprintf('✅ 2-D SOLPS maps (log n_e, log T_e, flux_e, total Ne, T_i, avg Ne Ua) plotted successfully using jet colormap.\n');

%% === Neon charge-state abundance at antenna region (region-averaged) ===
% Define antenna region (edit these bounds to match your antenna footprint)
R_ant_min = 8.20;
R_ant_max = 8.36;
Z_ant_min = -0.15;
Z_ant_max = +0.15;

maskAnt = (X >= R_ant_min) & (X <= R_ant_max) & (Z >= Z_ant_min) & (Z <= Z_ant_max);

% Neon indices in your current mapping (you used Ne0..Ne10+ = 6:16)
idxNe = 6:16;
NeNames = {'Ne^{0}','Ne^{1+}','Ne^{2+}','Ne^{3+}','Ne^{4+}', ...
           'Ne^{5+}','Ne^{6+}','Ne^{7+}','Ne^{8+}','Ne^{9+}','Ne^{10+}'};

% Total Ne density in region (sum over charge states)
nNe_total_map = zeros(size(ne_q));
for j = 1:numel(idxNe)
    nNe_total_map = nNe_total_map + ni_q{idxNe(j)};
end

% Region-averaged abundance (simple average over antenna mask)
% (If you prefer density-weighted, see option below.)
fracNe_region = zeros(numel(idxNe),1);
for j = 1:numel(idxNe)
    nq = ni_q{idxNe(j)};
    frac_map = nq ./ max(nNe_total_map, 1e-60);
    fracNe_region(j) = mean(frac_map(maskAnt), 'omitnan');
end

% Normalize to sum=1 (sometimes mean of ratios won't sum to 1 exactly)
fracNe_region = fracNe_region ./ nansum(fracNe_region);

fprintf('\n=== Neon charge-state abundance (region-averaged over antenna box) ===\n');
fprintf('Antenna box: R=[%.3f, %.3f] m, Z=[%.3f, %.3f] m\n', R_ant_min, R_ant_max, Z_ant_min, Z_ant_max);
for j = 1:numel(idxNe)
    fprintf('%10s : %.4f (%.2f%%)\n', NeNames{j}, fracNe_region(j), 100*fracNe_region(j));
end

figure('Color','w','Position',[100 100 750 450]);
bar(0:10, 100*fracNe_region);
xlabel('Ne charge state q (0..10)'); ylabel('Abundance [%]');
title('Neon charge-state abundance at antenna (region-averaged)');
grid on; box on;
%% === Helper: polygon "paint" mapper ===
function grid_vals = paintCells(cv_vals, r4xN, z4xN, X, Z)
    grid_vals = zeros(size(X));
    N = size(r4xN,2);
    for i = 1:N
        [in,on] = inpolygon(X,Z,r4xN(:,i),z4xN(:,i));
        if any(in(:)) || any(on(:))
            grid_vals(in|on) = cv_vals(i);
        end
    end
end