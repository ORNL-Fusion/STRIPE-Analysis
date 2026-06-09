clc; clear; close all; format long

%% --- Load GITRm per-face data ---
% path1 = '../diiid-helicon/DIII-D_helicon_runs_200882/gitrm-surface2.nc';
path1 = '../diiid-helicon/DIII-D_helicon_runs_196154/gitrm-surface2.nc';
ncid1 = netcdf.open(path1,'NOWRITE');
coords      = netcdf.getVar(ncid1, 2);   % flattened [9*sz x 1]
grossPlasma = netcdf.getVar(ncid1, 3);   % sz x 1
grossDepo   = netcdf.getVar(ncid1, 4);   % sz x 1
grossEroRaw = netcdf.getVar(ncid1, 5);   % sz x 1
netcdf.close(ncid1);

sz = size(grossEroRaw,1);

%% --- Unpack triangle vertices ---
coords_mat = reshape(coords, [9, sz]).';
x1 = coords_mat(:,1); y1 = coords_mat(:,2); z1 = coords_mat(:,3);
x2 = coords_mat(:,4); y2 = coords_mat(:,5); z2 = coords_mat(:,6);
x3 = coords_mat(:,7); y3 = coords_mat(:,8); z3 = coords_mat(:,9);

%% --- Triangle area (3D) ---
P1 = [x1 y1 z1];  P2 = [x2 y2 z2];  P3 = [x3 y3 z3];
v1 = P2 - P1;     v2 = P3 - P1;
cp = cross(v1, v2, 2);
ElementArea = 0.5 * sqrt(sum(cp.^2,2));
ElementArea(ElementArea==0) = eps;

%% --- Erosion per-area fluxes ---
% erosion_rate = 7.176454593335885e+20; % from GITR
erosion_rate = 4.090184441246757e+17;  % #196154 (units depend on your tally definition)
nP = 5e6;
erosionPP = erosion_rate / nP;

grossPlasmaA   = erosionPP .* grossPlasma   ./ ElementArea;
grossDepoA     = erosionPP .* grossDepo     ./ ElementArea;
grossEroRawA   = erosionPP .* grossEroRaw   ./ ElementArea;

grossEroA = grossPlasmaA + grossEroRawA;
netEroA   = grossEroA   - grossDepoA;
netDepoA  = grossDepoA  - grossEroRawA;

%% === NEW SECTION: Area-integrated erosion and deposition ===
grossEro_integrated = sum(grossEroA .* ElementArea);   % total gross erosion [particles/s]
grossDepo_integrated = sum(grossDepoA .* ElementArea); % total gross deposition [particles/s]
netEro_integrated = sum(netEroA .* ElementArea);       % total net erosion [particles/s]
netDepo_integrated = sum(netDepoA .* ElementArea);     % total net redeposition [particles/s]

fprintf('\n--- Integrated Results ---\n');
fprintf('Total Gross Erosion     = %.4e particles/s\n', grossEro_integrated);
fprintf('Total Gross Deposition  = %.4e particles/s\n', grossDepo_integrated);
fprintf('Total Net Erosion       = %.4e particles/s\n', netEro_integrated);
fprintf('Total Net Redeposition  = %.4e particles/s\n', netDepo_integrated);
fprintf('--------------------------\n\n');