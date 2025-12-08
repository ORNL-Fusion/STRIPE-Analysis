%% Initialization
close all;
clear all;
clc;

%% Load SOLPS Data
fileSOLPS = 'interpolated_values_196154.nc';
rS = ncread(fileSOLPS, 'gridr');    % Radial grid (Nr)
zS = ncread(fileSOLPS, 'gridz');    % Vertical grid (Nz)
neS = ncread(fileSOLPS, 'ne');      % Electron density (Nr x Nz)
TeS = ncread(fileSOLPS, 'te');      % Electron temperature (Nr x Nz)
gradTiS = -ncread(fileSOLPS, 'gradTi');
vrS = ncread(fileSOLPS, 'vr');
vzS = ncread(fileSOLPS, 'vz');
vtS = ncread(fileSOLPS, 'vt');

%% Data Cleaning
n_min = 1e10;
T_min = 10;
neS(neS <= 0 | ~isfinite(neS)) = NaN;
TeS(TeS <= 0 | ~isfinite(TeS)) = NaN;
vrS(~isfinite(vrS)) = NaN;
vzS(~isfinite(vzS)) = NaN;
vtS(~isfinite(vtS)) = NaN;
gradTiS(~isfinite(gradTiS)) = NaN;
TeS(TeS < T_min) = T_min;

%% Load Experimental Data at Z = 0.27
load('exp_profile_Z027.mat');  % R_exp, ne_exp

%% Compare SOLPS with Experimental Data
Z_target = 0.27;
[~, z_idx] = min(abs(zS - Z_target));
ne_solps_z027 = neS(:, z_idx);
ne_solps_interp = interp1(rS, ne_solps_z027, R_exp, 'linear', NaN);

figure;
semilogy(R_exp, ne_exp, 'ko-', 'DisplayName', 'Experimental');
hold on;
semilogy(R_exp, ne_solps_interp, 'r.-', 'DisplayName', 'SOLPS @ Z=0.27');
xlabel('R [m]');
ylabel('n_e [m^{-3}]');
legend;
title('Comparison of SOLPS vs. Experimental Density at Z = 0.27 m');
grid on;

%% Fit Experimental Tail
idx_fit = R_exp > 1.0;
p_exp = polyfit(R_exp(idx_fit), log(ne_exp(idx_fit)), 1);
ne_fit_func = @(r) exp(polyval(p_exp, r));

%% Extrapolate on 1D Midplane Grid
mpx = linspace(2.14333, 2.5, 1000);  % DIII-D #196154
mpy = ones(size(mpx)) * Z_target;
densityAtPlane = interp2(rS, zS, neS', mpx, mpy, 'linear', NaN);
densityAtPlane(isnan(densityAtPlane)) = 0;

% Transition from SOLPS to experimental fit
r_sep = 2.20675;  % Approximate separatrix
interpfn = (mpx - rS(1)) / (r_sep - rS(1));
interpfn = min(max(interpfn, 0), 1);
extrapolatedne1d = (1 - interpfn) .* densityAtPlane + interpfn .* ne_fit_func(mpx);

%% Plot Final 1D Profile
figure;
semilogy(mpx, densityAtPlane, 'b-', 'DisplayName', 'SOLPS @ Z=0.27');
hold on;
semilogy(mpx, extrapolatedne1d, 'r--', 'DisplayName', 'Extrapolated (Exp Fit)');
xlabel('R [m]');
ylabel('n_e [m^{-3}]');
title('SOLPS vs. Experimental-Guided Extrapolation');
legend;
grid on;

%% Save for 2D Use
save('ne_extrapolated_z027.mat', 'mpx', 'extrapolatedne1d');