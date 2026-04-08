%% MATLAB Script: Bremsstrahlung Emission via Kunze (Standard + Kappa-Corrected)
% Implements Kunze Eq. (6.81) with classical Gaunt factor
% Adds correction for kappa-distributed nonthermal electrons
% Units: [m], [s], [J], [eV], [W]

clear; clc; close all;

%% Physical Constants
e     = 1.602176634e-19;       % Elementary charge [C]
h     = 6.62607015e-34;        % Planck constant [J·s]
hbar  = h/(2*pi);              % Reduced Planck constant [J·s]
c     = 2.99792458e8;          % Speed of light [m/s]
kB    = 1.380649e-23;          % Boltzmann constant [J/K]
alpha = 7.2973525693e-3;       % Fine-structure constant
a0    = 5.29177210903e-11;     % Bohr radius [m]
ER_eV = 13.60569193;           % Rydberg energy [eV]
ER_J  = ER_eV * e;             % Rydberg energy [J]
me    = 9.10938356e-31;        % Electron mass [kg]
eps0  = 8.8541878128e-12;      % Vacuum permittivity [F/m]

%% Plasma & Geometry Parameters
Z     = 1;                     % Ion charge
n0    = 1e19;                  % Electron/ion density [m^-3]
Te_eV = 1e4;                   % Reference temperature [eV]
Te_J  = Te_eV * e;             % [J]
kappa = 3;                     % Kappa index for nonthermal correction

% Cylinder geometry
r     = 0.06;                  % Radius [m]
L     = 5.0;                   % Length [m]
V_cyl = pi * r^2 * L;          % Volume [m^3]

%% Prefactor from Kunze Eq. (6.82)
prefactor = 4.108e-46;         % [W·m^4]

%% Photon Energy Grid and Wavelength
E_keV  = linspace(0.02, 20, 200);     % [keV]
E_J    = E_keV * 1e3 * e;             % [J]
lambda = (h * c) ./ E_J;              % Wavelength [m]

%% Gaunt Factor (velocity-averaged, classical log formula)
gamma = 0.577215664901532;
eta_Z  = Z * e^2 ./ (hbar * c) .* lambda / (2 * pi);
eta_nu = h ./ (2 * pi * me * c .* lambda);
Gff    = (sqrt(3)/pi) .* log(1 ./ (eta_nu .* max(1, exp(gamma) .* eta_Z)));
Gff(Gff < 1) = 1;

%% Maxwellian-Based Kunze Formula
sqrt_term = sqrt(ER_J / Te_J);  % sqrt(E_R / k_B T)
exp_term  = exp(- (h * c) ./ (lambda .* Te_J));

epsilon_lambda_maxwell = prefactor .* n0^2 .* Z^2 .* sqrt_term ...
    .* (1 ./ lambda.^2) .* exp_term .* Gff;

%% Kappa-Corrected Kunze Formula
exp_term_kappa = (1 + (h * c) ./ (lambda * kappa * Te_J)).^(-(kappa + 1));

epsilon_lambda_kappa = prefactor .* n0^2 .* Z^2 .* sqrt_term ...
    .* (1 ./ lambda.^2) .* exp_term_kappa .* Gff;

%% Convert Both to Spectral Emissivity per Energy ε_E [W·m^-3·J^-1]
dlambda_dE = (h * c) ./ (E_J.^2);
epsilon_energy_maxwell = epsilon_lambda_maxwell .* abs(dlambda_dE);
epsilon_energy_kappa   = epsilon_lambda_kappa   .* abs(dlambda_dE);

%% Integrate Total Emission (Maxwellian)
power_density_maxwell       = trapz(E_J, epsilon_energy_maxwell);
photon_rate_density_maxwell = trapz(E_J, epsilon_energy_maxwell ./ E_J);
total_power_maxwell         = power_density_maxwell * V_cyl;
total_photon_rate_maxwell   = photon_rate_density_maxwell * V_cyl;

%% Integrate Total Emission (Kappa-Corrected)
power_density_kappa       = trapz(E_J, epsilon_energy_kappa);
photon_rate_density_kappa = trapz(E_J, epsilon_energy_kappa ./ E_J);
total_power_kappa         = power_density_kappa * V_cyl;
total_photon_rate_kappa   = photon_rate_density_kappa * V_cyl;

%% Display Results
fprintf('\n=== Kunze Bremsstrahlung: Maxwellian ===\n');
fprintf('T_e                     : %.1f eV\n', Te_eV);
fprintf('Power density           : %.3e W/m^3\n', power_density_maxwell);
fprintf('Total power             : %.3e W\n', total_power_maxwell);
fprintf('Photon rate density     : %.3e photons/m^3/s\n', photon_rate_density_maxwell);
fprintf('Total photon rate       : %.3e photons/s\n', total_photon_rate_maxwell);

fprintf('\n=== Kunze Bremsstrahlung: Kappa-Corrected (κ = %.1f) ===\n', kappa);
fprintf('Power density           : %.3e W/m^3\n', power_density_kappa);
fprintf('Total power             : %.3e W\n', total_power_kappa);
fprintf('Photon rate density     : %.3e photons/m^3/s\n', photon_rate_density_kappa);
fprintf('Total photon rate       : %.3e photons/s\n', total_photon_rate_kappa);

%% Plot Spectral Emissivity Comparison
figure('Name','Spectral Emissivity Comparison','Position',[100 100 600 400]);
plot(E_keV, epsilon_energy_maxwell, 'b-', 'LineWidth', 1.5); hold on;
plot(E_keV, epsilon_energy_kappa,   'm--', 'LineWidth', 1.5);
xlabel('Photon Energy (keV)');
ylabel('\epsilon_E [W·m^{-3}·J^{-1}]');
legend('Maxwellian', sprintf('Kappa-corrected (κ = %.1f)', kappa));
title('Bremsstrahlung Emissivity Spectrum');
grid on;