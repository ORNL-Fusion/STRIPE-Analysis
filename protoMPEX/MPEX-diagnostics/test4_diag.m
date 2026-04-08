%% Bremsstrahlung from 1D–2V κ-Distribution with Power-Balance and Kunze Comparison
% Units carefully tracked throughout: [W/m³], [J], [keV], [m/s], [m²·J⁻¹], [m⁻³·J⁻¹]

clear; clc; close all;

%% Physical Constants
e     = 1.602176634e-19;      % Elementary charge [C]
me    = 9.10938356e-31;       % Electron mass [kg]
c     = 2.998e8;              % Speed of light [m/s]
eps0  = 8.854187817e-12;      % Vacuum permittivity [F/m]
alpha = 1/137.035999084;      % Fine-structure constant [unitless]
r_e   = e^2 / (4*pi*eps0*me*c^2);  % Classical electron radius [m]
Ry_J  = 13.605693 * e;        % Rydberg energy [J]
h     = 6.62607015e-34;       % Planck constant [J·s]
kB    = 1.380649e-23;         % Boltzmann constant [J/K]

%% MPEX Plasma & Geometry
P_input    = 400e3;           % Input power [W]
P_absorbed = 190e3;           % Absorbed by plasma [W]
P_lost     = P_input - P_absorbed;

r_cyl = 0.06; z1 = 3; z2 = 8; % Plasma column: radius and length [m]
V_cyl = pi*r_cyl^2*(z2 - z1); % Volume [m³]

P_vol_input    = P_input / V_cyl;    % Power density input [W/m³]
P_vol_absorbed = P_absorbed / V_cyl;
P_vol_lost     = P_lost / V_cyl;

ne_total = 1e19;   % Total electron density [m⁻³]
Z        = 1;      % Ion charge (hydrogen plasma assumed)

%% Distribution Parameters
Te_eV = 50;              % Effective electron temperature [eV]
Te_J  = Te_eV * e;       % [J]
kappa = 2;               % Kappa distribution parameter

%% Velocity grid (1D-parallel, 1D-perpendicular)
Nv     = 500;
v_th   = sqrt(2*Te_J/me);       % Thermal speed [m/s]
vmax   = 100 * v_th;            % Max velocity for grid [m/s]
v_par  = linspace(-vmax, vmax, Nv);   % [m/s]
v_perp = linspace(0, vmax, Nv);       % [m/s]
[Vpar, Vperp] = meshgrid(v_par, v_perp);
dv_par  = v_par(2)-v_par(1);    % [m/s]
dv_perp = v_perp(2)-v_perp(1);  % [m/s]

V2     = Vpar.^2 + Vperp.^2;     % Total speed squared [m²/s²]
E_grid = 0.5 * me * V2;          % Kinetic energy grid [J]
jac    = 2*pi .* Vperp;          % Cylindrical Jacobian factor

%% κ-distribution and normalization (units: [m⁻⁶·s³], normalized to [m⁻³])
f_unnorm = (1 + V2./(kappa*v_th^2)).^(-(kappa+1));  % unnormalized f(v)
normF    = sum(f_unnorm .* jac,'all') * dv_par * dv_perp;
f_v      = f_unnorm * (ne_total / normF);          % normalized f(v) [m⁻6·s³]

%% Convert f(v) to f(E) [m⁻³·J⁻¹]
E_bins = linspace(0, 1000e3*e, 1000);    % Energy bins [J]
E_cent = 0.5 * (E_bins(1:end-1) + E_bins(2:end)); % Bin centers [J]
dE     = diff(E_bins);
E_flat = E_grid(:);      % Flattened energy grid [J]
f_flat = f_v(:) .* jac(:) * dv_par * dv_perp;  % [m⁻³]

fE = zeros(size(E_cent));   % Energy distribution [m⁻³·J⁻¹]
for i = 1:numel(E_cent)
    mask = E_flat >= E_bins(i) & E_flat < E_bins(i+1);
    fE(i) = sum(f_flat(mask)) / dE(i);
end

%% Photon energy range
K_eV = linspace(10, 400, 1000);  % Photon energy [keV]
K_J  = K_eV * e;                 % Photon energy [J]

%% Velocity and Cross Section Terms
U  = @(E) c * sqrt(E .* (E + 2*me*c^2)) ./ (E + me*c^2);     % Relativistic speed [m/s]
eta  = @(v) alpha * c ./ v;                                  % Sommerfeld parameter [unitless]
F_EW = @(v) (1 - exp(-2*pi*eta(v))) ./ (2*pi*eta(v));        % Elwert correction [unitless]

% Born approximation σB(E,K) [m²/J], from Kunze Eq. 6.80
sigmaB = @(E,K) (16/3) * alpha * r_e^2 .* log(max(2*E./K,1.01)) ./ K;

% Total corrected free-free cross section σ_ff(E,K) [m²/J]
sigma_ff = @(E,K) Z^2 * sigmaB(E,K) .* F_EW(U(E));

%% ε_kappa (numerical emissivity) [W/m³/J]
eps_k = zeros(size(K_J));
for i = 1:numel(K_J)
    m = E_cent > K_J(i);
    eps_k(i) = trapz(E_cent(m), fE(m) .* sigma_ff(E_cent(m),K_J(i)) .* U(E_cent(m)));
end
eps_k = eps_k / e;  % Convert [W/m³/J] → [W/m³/keV]

%% ε_maxwellian (reference) [W/m³/keV]
fM = @(E) ne_total * 2/sqrt(pi) * (1/Te_J)^(3/2) .* sqrt(E) .* exp(-E/Te_J);  % [m⁻³·J⁻¹]
eps_M = zeros(size(K_J));
for i = 1:numel(K_J)
    m = E_cent > K_J(i);
    eps_M(i) = trapz(E_cent(m), fM(E_cent(m)) .* sigma_ff(E_cent(m),K_J(i)) .* U(E_cent(m)));
end
eps_M = eps_M / e;  % Convert to [W/m³/keV]

%% ε_Kunze (analytic, Eq. 6.81) [W/m³/keV]
g_Kunze = 1.2;       % Approximate Gaunt factor
A_Kunze = 1.4e-40;   % Pre-factor from Kunze Eq. 6.81 [W·m³/Hz]
eps_Kunze = A_Kunze * Z^2 * ne_total^2 * (1/sqrt(Te_eV)) .* (1 ./ K_J) * g_Kunze;  % [W/m³/J]
eps_Kunze = eps_Kunze / e;  % Convert to [W/m³/keV]

%% Plot spectral emissivity
figure;
semilogy(K_eV, eps_k, 'b-', 'LineWidth', 1.5); hold on;
semilogy(K_eV, eps_M, 'r--', 'LineWidth', 1.5);
semilogy(K_eV, eps_Kunze, 'k-.', 'LineWidth', 1.5);
legend('\kappa (numeric)','Maxwellian (numeric)','Kunze Eq.6.81');
xlabel('Photon Energy [keV]');
ylabel('Spectral Emissivity \epsilon [W/m³/keV]');
title('Bremsstrahlung Emissivity: κ vs Maxwellian vs Kunze');
grid on;

%% Integrated emissivity (Total Radiated Power)
P_rad_vol_kappa = trapz(K_J, eps_k * e);      % [W/m³]
P_rad_vol_max   = trapz(K_J, eps_M * e);
P_rad_vol_kunze = trapz(K_J, eps_Kunze * e);

P_rad_kappa = P_rad_vol_kappa * V_cyl;   % [W]
P_rad_max   = P_rad_vol_max * V_cyl;
P_rad_kunze = P_rad_vol_kunze * V_cyl;

%% Output
fprintf('--- X-ray Power Summary ---\n');
fprintf('κ  : Vol emissivity = %.3e W/m³ | Total = %.3e W\n', P_rad_vol_kappa, P_rad_kappa);
fprintf('Max: Vol emissivity = %.3e W/m³ | Total = %.3e W\n', P_rad_vol_max, P_rad_max);
fprintf('Kunze Eq.6.81       = %.3e W/m³ | Total = %.3e W\n', P_rad_vol_kunze, P_rad_kunze);