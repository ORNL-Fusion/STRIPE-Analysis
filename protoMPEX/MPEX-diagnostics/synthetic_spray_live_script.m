%% MATLAB Script: Bremsstrahlung Emission via Kunze (Particle-Based Gaunt Factor)
% Implements analytic free–free formula from Kunze Chapter 6 (Eq. 6.81)
% Gaunt factor is computed using log-based expression from classical approximation
% Units: [m], [s], [J], [eV], [W]

clear; clc; close all;

%% Physical Constants
e     = 1.602176634e-19;      % Elementary charge [C]
h     = 6.62607015e-34;       % Planck constant [J·s]
hbar  = h/(2*pi);             % Reduced Planck constant [J·s]
c     = 2.99792458e8;         % Speed of light [m/s]
kB    = 1.380649e-23;         % Boltzmann constant [J/K]
alpha = 7.2973525693e-3;      % Fine-structure constant [unitless]
a0    = 5.29177210903e-11;    % Bohr radius [m]
ER_eV = 13.60569193;          % Rydberg energy [eV]
ER_J  = ER_eV * e;            % Rydberg energy [J]
me    = 9.10938356e-31;       % Electron mass [kg]

%% Plasma & Geometry Parameters
Z     = 1;                    % Ion charge (hydrogen)
n0    = 1e19;                 % Electron/ion density [m^-3]
eV0   = 1e4;                  % Electron temperature [eV]
Te_eV = eV0;
Te_J  = Te_eV * e;            % Electron temperature in joules

% Cylinder geometry
r     = 0.06;                 % Radius [m]
L     = 5.0;                  % Length [m]
V_cyl = pi * r^2 * L;         % Volume [m^3]

%% Prefactor from Kunze Eq. (6.82)
prefactor = 4.108e-46;        % [W·m^4]

%% Spectral Emissivity per wavelength ε_λ [W·m^-3·m^-1]
% Photon energy range 0.1–20 keV
E_keV = linspace(0.02, 20, 200);       % [keV]
E_J   = E_keV * 1e3 * e;              % [J]
lambda = (h * c) ./ E_J;             % Wavelength [m]

%% Gaunt Factor (velocity-averaged, classical expression)
gamma = 0.577215664901532;   % Euler–Mascheroni constant

% Estimate using v ~ c for η_nu, η_Z definitions
eta_Z  = Z * e^2 ./ (hbar * c) .* lambda / (2 * pi);     % [unitless]
eta_nu = h ./ (2 * pi * me * c .* lambda);               % [unitless]

Gff = (sqrt(3)/pi) .* log(1 ./ (eta_nu .* max(1, exp(gamma) .* eta_Z)));
Gff(Gff < 1) = 1;  % Enforce physical lower bound

%% Kunze Eq. 6.81: ε_λ = C * n0^2 * Z^2 * sqrt(ER / Te) * 1/λ^2 * exp(-hc/(λ Te)) * g_ff
sqrt_term = sqrt(ER_J / Te_J);  % sqrt(E_R / k_B T)
epsilon_lambda = prefactor .* n0^2 .* Z^2 .* sqrt_term ...
    .* (1 ./ lambda.^2) .* exp(- (h * c) ./ (lambda .* Te_J)) .* Gff;

%% Convert to spectral emissivity per energy ε_E [W·m^-3·J^-1]
% Use dλ/dE = (h c)/E^2
dlambda_dE = (h * c) ./ (E_J.^2);
epsilon_energy = epsilon_lambda .* abs(dlambda_dE);

%% Integrate total emission
power_density         = trapz(E_J, epsilon_energy);               % [W/m^3]
total_power           = power_density * V_cyl;                   % [W]
photon_rate_density   = trapz(E_J, epsilon_energy ./ E_J);       % [photons/m^3/s]
total_photon_rate     = photon_rate_density * V_cyl;             % [photons/s]

%% Display Results
fprintf('=== Kunze Bremsstrahlung (T_e = %.1f eV) ===\n', Te_eV);
fprintf('Cylinder volume       : %.3f m^3\n', V_cyl);
fprintf('Density (n0)          : %.3e m^-3\n', n0);
fprintf('Power density         : %.3e W/m^3\n', power_density);
fprintf('Total power           : %.3e W\n', total_power);
fprintf('Photon rate density   : %.3e photons/m^3/s\n', photon_rate_density);
fprintf('Total photon rate     : %.3e photons/s\n', total_photon_rate);

%% Plot Spectra
figure('Name','Spectral Emissivity vs Wavelength','Position',[100 100 600 400]);
plot(lambda * 1e9, epsilon_lambda, 'LineWidth', 1.5);
xlabel('Wavelength (nm)');
ylabel('\epsilon_\lambda [W·m^{-3}·nm^{-1}]');
title('Bremsstrahlung Emissivity per Wavelength');
grid on;

figure('Name','Spectral Emissivity vs Energy','Position',[100 100 600 400]);
plot(E_keV, epsilon_energy, 'LineWidth', 1.5);
xlabel('Photon Energy (keV)');
ylabel('\epsilon_E [W·m^{-3}·J^{-1}]');
title('Bremsstrahlung Emissivity per Energy');
grid on;