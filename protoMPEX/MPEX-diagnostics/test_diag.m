%% MATLAB Script: Analytic κ-Corrected Kunze Bremsstrahlung (Eq. 6.81)

clear; clc; close all;

%% MATLAB Script: Analytic κ-Corrected Kunze Bremsstrahlung (Eq. 6.81)

clear; clc; close all;

%% 1) Physical constants
e     = 1.602176634e-19;    % Elementary charge [C]
h     = 6.62607015e-34;     % Planck constant [J·s]
hbar  = h/(2*pi);          % Reduced Planck constant [J·s]
c     = 2.99792458e8;       % Speed of light [m/s]
kB    = 1.380649e-23;       % Boltzmann constant [J/K]
me    = 9.10938356e-31;     % Electron mass [kg]
ER_eV = 13.60569193;        % Rydberg energy [eV]
ER_J  = ER_eV * e;          % Rydberg energy [J]
alpha = 7.2973525693e-3;    % Fine-structure constant


%% 2) Plasma & geometry
Z     = 1;                  % ion charge
n0    = 1e19;               % electron+ion density [m^-3]
Te_eV = 1e4;                % eV
Te_J  = Te_eV * e;          % J
kappa = 3;                  % κ index
r     = 0.06; L = 5.0;      % m
V_cyl = pi*r^2*L;           % m^3

%% 3) Photon grid & wavelength
E_keV = linspace(0.02, 20, 500);   % [keV]
E_J   = E_keV*1e3*e;               % [J]
lambda= (h*c)./E_J;                % [m]

%% 4) Gaunt factor (classical log form)
gamma_E = 0.577215664901532;
eta_Z   = Z*e^2./(hbar*c).*lambda/(2*pi);
eta_nu  = h./(2*pi*me*c.*lambda);
Gff     = (sqrt(3)/pi).*log(1./(eta_nu.*max(1,exp(gamma_E).*eta_Z)));
Gff(Gff<1) = 1;

%% 5) Kunze prefactor (Eq. 6.82)
pref = 4.108e-46;   % [W·m^4]

%% 6) Compute emissivities
eps_lambda_maxwell = zeros(size(lambda));
eps_lambda_kappa   = zeros(size(lambda));
for i = 1:numel(lambda)
    Lmb = lambda(i);
    % Maxwellian factor
    exp_M = exp(-h*c/(Lmb*Te_J));
    % κ‐generalized factor
    exp_k = (1 + h*c/(Lmb*kappa*Te_J))^(-(kappa+1));
    % Kunze Eq.6.81
    common = pref * n0^2 * Z^2 * sqrt(ER_J/Te_J) * (1/Lmb^2) * Gff(i);
    eps_lambda_maxwell(i) = common * exp_M;
    eps_lambda_kappa(i)   = common * exp_k;
end

%% 7) Convert to per‐energy [W·m^-3·J^-1]
dl_dE = (h*c)./E_J.^2;
eps_E_M = eps_lambda_maxwell .* abs(dl_dE);
eps_E_K = eps_lambda_kappa   .* abs(dl_dE);

%% 8) Integrate for power density and total power
P_den_M = trapz(E_J, eps_E_M);
P_den_K = trapz(E_J, eps_E_K);
P_tot_M = P_den_M * V_cyl;
P_tot_K = P_den_K * V_cyl;

%% 9) Print results
idx1  = find(abs(E_keV-1)==min(abs(E_keV-1)),1);
idx5  = find(abs(E_keV-5)==min(abs(E_keV-5)),1);
idx10 = find(abs(E_keV-10)==min(abs(E_keV-10)),1);

fprintf('\n=== Kunze Bremsstrahlung ===\n');
fprintf('Maxwellian: Power density = %.3e W/m^3, Total = %.3e W\n', P_den_M, P_tot_M);
fprintf(' κ-corrected: Power density = %.3e W/m^3, Total = %.3e W\n', P_den_K, P_tot_K);
fprintf('ε_E(1 keV): M = %.3e, K = %.3e [W/m^3/J]\n', eps_E_M(idx1), eps_E_K(idx1));
fprintf('ε_E(5 keV): M = %.3e, K = %.3e [W/m^3/J]\n', eps_E_M(idx5), eps_E_K(idx5));
fprintf('ε_E(10 keV): M=%.3e, K=%.3e [W/m^3/J]\n',  eps_E_M(idx10),eps_E_K(idx10));

%% 10) Plot comparison
figure;
plot(E_keV, eps_E_M, 'b-', 'LineWidth',1.5); hold on;
plot(E_keV, eps_E_K, 'r--','LineWidth',1.5); hold off;
xlabel('Photon Energy (keV)');
ylabel('\epsilon_E [W·m^{-3}·J^{-1}]');
legend('Maxwellian','κ-corrected');
title('Bremsstrahlung Emissivity (Kunze Eq.6.81)');
grid on;