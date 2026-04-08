%% Bremsstrahlung from Power-Law Electrons (Bernhardi 1982 Inspired)
clear; clc;

%% Physical Constants
e     = 1.602e-19;      % Elementary charge [C]
me    = 9.109e-31;      % Electron mass [kg]
c     = 2.998e8;        % Speed of light [m/s]
eps0  = 8.854e-12;      % Vacuum permittivity [F/m]
alpha = 1/137;          % Fine-structure constant
r_e   = e^2 / (4*pi*eps0*me*c^2);  % Classical electron radius [m]

%% Plasma Parameters
nhot = 1e19;            % Hot electron density [m^-3]
ni   = nhot;            % Assume quasi-neutrality
s    = 1.9;             % Spectral index from paper

%% Energy Grid (10 keV to 1 MeV)
E_keV = linspace(1, 1000, 10000);    % Electron kinetic energy [keV]
E_J   = E_keV * e;                   % Energy [J]
dE    = diff(E_J); dE(end+1) = dE(end); % Energy spacing [J]

%% Power-law Distribution Function (not normalized)
f_E = E_J.^(-s);
% Normalize so that ∫f(E)dE = nhot
normF = trapz(E_J, f_E);
f_E   = f_E * (nhot / normF);  % Now in units of [m^-3/J]

%% Photon Grid (10 to 350 keV)
K_eV = linspace(1, 350, 3000);      % Photon energy [keV]
K_J  = K_eV * e;

%% Bremsstrahlung Cross-Section (Elwert corrected)
vE   = @(E) sqrt(2*E/me);
eta  = @(E) alpha * c ./ vE(E);
g_ff = @(E) pi*eta(E) ./ (1 - exp(-2*pi*eta(E)));  % Elwert factor

sigmaB = @(E,K) (16/3)*alpha*r_e^2 .* log(max(2*E./K, 1.01)) ./ K;   % Born
sigmaE = @(E,K) sigmaB(E,K) .* g_ff(E);                             % Elwert

%% Compute Spectral Emissivity ε(K)
epsilon = zeros(size(K_J));  % [W/m^3/keV]

for i = 1:length(K_J)
    K = K_J(i);
    mask = E_J > K;
    E_sample = E_J(mask);
    f_sample = f_E(mask);
    sigma_vals = sigmaE(E_sample, K);
    v_vals     = vE(E_sample);

    epsilon(i) = ni * trapz(E_sample, sigma_vals .* f_sample .* v_vals);
end

epsilon = epsilon / e;  % Convert from per Joule to per keV

%% Plot Spectral Emissivity
figure;
semilogy(K_eV, epsilon, 'LineWidth', 1.5);
xlabel('Photon Energy [keV]');
ylabel('\epsilon(K) [W/m^3/keV]');
title('X-ray Bremsstrahlung from Power-Law Electrons');
grid on;

%% Total Power Estimate
eps_total = trapz(K_J, epsilon * e);  % Integrate over K [W/m^3]
fprintf('Volume Bremsstrahlung Power Density: %.3e W/m^3\n', eps_total);

% Assume a small cylindrical plasma
r_cyl = 0.03; z_len = 1.0;
V = pi*r_cyl^2 * z_len;
P_total = eps_total * V;
fprintf('Total Radiated X-ray Power: %.3e W\n', P_total);
