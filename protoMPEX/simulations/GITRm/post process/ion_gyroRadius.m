% Deuterium ion gyro-radius calculation
% Assumes v_perp ~ ion sound speed sqrt(Te/mi)

clear; clc;

%% ---- Inputs ----
Te_eV = 5;        % Electron temperature [eV]
B     = 0.5;      % Magnetic field [T]
Z     = 1;        % Ion charge state (D+)

%% ---- Constants ----
e  = 1.602176634e-19;   % Elementary charge [C]
mp = 1.67262192369e-27; % Proton mass [kg]

%% ---- Deuterium mass ----
mi = 2 * mp;            % Deuterium ion mass [kg]

%% ---- Convert temperature ----
Te_J = Te_eV * e;       % [J]

%% ---- Ion sound speed (perpendicular velocity estimate) ----
cs = sqrt(Te_J / mi);   % [m/s]

%% ---- Ion gyro-radius ----
rho_i = mi * cs / (Z * e * B);   % [m]

%% ---- Two gyro-radii ----
two_rho_i = 2 * rho_i;           % [m]

%% ---- Output ----
fprintf('Deuterium ion gyro-radius rho_i = %.4e m (%.3f mm)\n', rho_i, rho_i*1e3);
fprintf('Two gyro-radii 2*rho_i         = %.4e m (%.3f mm)\n', two_rho_i, two_rho_i*1e3);