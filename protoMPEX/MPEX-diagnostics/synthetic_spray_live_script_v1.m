%% Simple Bremsstrahlung Script (Kunze Formula, 1D–2V Particle Data)
clear; clc;

%% Physical Constants
e   = 1.602e-19;      % Elementary charge [C]
h   = 6.626e-34;      % Planck constant [J·s]
hbar= h/(2*pi);       % Reduced Planck constant
c   = 3e8;            % Speed of light [m/s]
me  = 9.11e-31;       % Electron mass [kg]
kB  = 1.38e-23;       % Boltzmann constant
ER  = 13.6 * e;       % Rydberg energy [J]

%% Plasma Setup
Z     = 1;            % Ion charge
N_particles = 1e6;
r_max = 0.06; L = 5.0;
V_cyl = pi * r_max^2 * L;  % [m^3]
ni    = N_particles / V_cyl;

% Generate sample velocities (replace with real data from PICOS++)
v_th  = 5e6;  % thermal scale [m/s]
v_par  = (2*rand(N_particles,1) - 1) * 5*v_th;
v_perp = rand(N_particles,1) * 5*v_th;
V2_particles = v_par.^2 + v_perp.^2;
Vmag = sqrt(V2_particles);

%% Photon Energy Grid
E_keV = linspace(0.1, 20, 200);
E_J   = E_keV * 1e3 * e;
lambda = (h * c) ./ E_J;

%% Gaunt Factor (Classical Log-Based)
gamma = 0.5772;
eta_Z  = Z * e^2 ./ (hbar * Vmag);
eps_omega = zeros(size(E_J));

for j = 1:length(E_J)
    omega = E_J(j);
    eta_nu = hbar * omega ./ (2 * me * Vmag.^2);
    g_ff = (sqrt(3)/pi) .* log(1 ./ (eta_nu .* max(1, exp(gamma) * eta_Z)));
    g_ff(g_ff < 1) = 1;

    % Born cross-section (approximate)
    alpha = 1/137;  % Fine-structure constant

% Born cross-section (approximate)
E_part = 0.5 * me * V2_particles;
log_term = log(2 * E_part / (me * c^2) + 1);
sigma = (16/3) * Z^2 * alpha * (2.82e-15)^2 * ...
        (me * c^2 ./ omega) .* log_term;

    % Particle-based emissivity [J/s per particle]
    emissivity = sigma .* Vmag .* g_ff .* (h * omega);

    % Average and scale by density
    eps_omega(j) = ni * mean(emissivity);
end

%% Plot
figure;
plot(E_keV, eps_omega, 'LineWidth', 1.5);
xlabel('Photon Energy (keV)');
ylabel('\epsilon_\omega [W·m^{-3}·Hz^{-1}]');
title('Bremsstrahlung Emissivity (Particle-Based)');
grid on;