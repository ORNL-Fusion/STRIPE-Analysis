%% Synthetic Spray Diagnostics: X-ray Production in a Uniform Plasma Cylinder
% Computes bremsstrahlung X-ray emission from a 2D kappa-distributed electron population
% contained within a uniform cylindrical plasma (r=0.06 m, z∈[3,8] m). We derive Tₑ self-consistently,
% and use both Kunze's free–free formula and Born approximation (log-based Gaunt factor).

% References:
% [1] H. A. Bethe & W. Heitler, Proc. Roy. Soc. A146 (1934) 83.
% [2] H. J. Kunze, Introduction to Plasma Spectroscopy, Springer (2009), Ch. 6.

% Physical Constants
e    = 1.602e-19;             % Elementary charge [C]
hbar = 1.055e-34;             % Reduced Planck constant [J·s]
c    = 3e8;                   % Speed of light [m/s]
me   = 9.109e-31;             % Electron mass [kg]
alpha= 1/137;                 % Fine-structure constant
eps0 = 8.854e-12;             % Vacuum permittivity [F/m]
Z    = 1;                     % Ion charge state
r_e  = e^2 / (4 * pi * eps0 * me * c^2);  % Classical electron radius [m]
ER   = 13.605693 * e;         % Rydberg energy [J]
kB   = 1.380649e-23;          % Boltzmann constant [J/K]

% Geometry
r_max   = 0.06;               % Cylinder radius [m]
z_start = 3;                  % Start z position [m]
z_end   = 8;                  % End z position [m]
L       = z_end - z_start;   % Length [m]
V_cyl   = pi * r_max^2 * L;  % Volume [m^3]

% Temperature and kappa
kappa = 3;
T_e_eV = 100;                % Temperature in eV
T_e_k  = T_e_eV * e;         % Convert to joules
v_th   = sqrt(2 * T_e_k / me);  % Thermal velocity [m/s] Nv = 200;
v_par = linspace(-5*v_th,5*v_th,Nv);
v_perp = linspace(0,5*v_th,Nv);
[Vpar,Vperp] = meshgrid(v_par,v_perp);
V2 = Vpar.^2 + Vperp.^2;
f2D = (1 + V2./(kappav_th^2)).^(-(kappa+1));
dv_par = v_par(2)-v_par(1);
dv_perp = v_perp(2)-v_perp(1);
norm2D = sum(sum(f2D.(2*pi.*Vperp))) * dv_par * dv_perp;
f2D = f2D / norm2D;

% Particle Initialization
N_particles = 1e6;
r_particles = r_max * sqrt(rand(N_particles,1));
theta_particles = 2*pi * rand(N_particles,1);
z_particles = z_start + L * rand(N_particles,1);
v_par_particles = (2*rand(N_particles,1)-1) * 5*v_th;
v_perp_particles = rand(N_particles,1) * 5*v_th;
V2_particles = v_par_particles.^2 + v_perp_particles.^2;

% Density
ni = N_particles / V_cyl;

% Electron Temperature
mean_V2 = sum(sum(V2 .* f2D .* (2*pi.*Vperp))) * dv_par * dv_perp;
Te_K = me * mean_V2 / (3 * kB);
Te_eV = Te_K * kB / e;

% Emissivity Calculation (Born Approx.)
omega_keV = linspace(0.1,20,200);
omega_J = omega_keV1e3e;
lambda_m = (hbar * c) ./ omega_J * 2*pi;
eps_omega = zeros(size(omega_J));
Vmag = sqrt(V2_particles);
gamma = 0.577215664901532;

for j = 1:numel(omega_J)
omega = omega_J(j);
eta_Z  = Z * e^2 ./ (hbar * Vmag);
eta_nu = hbar * omega ./ (2 * me .* Vmag.^2);
g_ff = (sqrt(3)/pi) .* log( 1 ./ (eta_nu .* max(1, exp(gamma) .* eta_Z)) );
g_ff(g_ff < 1) = 1;
E_part = 0.5 * me .* V2_particles;
log_term = log(2*E_part/(me*c^2) + 1);
sigma = (16/3) * Z^2 * alpha * r_e^2 .* (me*c^2./omega) .* log_term;
emissivity = sigma .* Vmag .* g_ff;
eps_omega(j) = ni * mean(emissivity);
end

% Kunze Analytic
g_ff_kunze = 1;
ne = ni;
eps_kunze = 4.108e-46 * ni * ne * Z^2 * sqrt(ER / (kB * Te_K))...
.* (1 ./ lambda_m.^2) .* exp(-hbar * c * 2 * pi ./ (lambda_m * kB * Te_K))...
* g_ff_kunze;
eps_kunze = eps_kunze / max(eps_kunze);
eps_born = eps_omega / max(eps_omega);

% Plot Comparison
figure;
plot(omega_keV, eps_born, 'k-', 'LineWidth', 2); hold on;
plot(omega_keV, eps_kunze, 'r–', 'LineWidth', 2);
xlabel('Photon Energy (keV)'); ylabel('Normalized \epsilon_X(\omega)');
title('Bremsstrahlung Emissivity: Born vs Kunze'); legend('Born Approx.', 'Kunze Eq. 6.81'); grid on;

% Integration
power_density = trapz(omega_J, eps_omega .* omega_J);
phi_omega = eps_omega ./ (hbar .* omega_J);
photon_density = trapz(omega_J, phi_omega);
total_power = power_density * V_cyl;
total_photons = photon_density * V_cyl;

% Convert to eV
power_density_eV = power_density / e;
total_power_eV_per_s = total_power / e;

% Display
fprintf('\n=== X-ray Production in Cylinder ===\n');
fprintf('Cylinder volume: %.3f m^3\n', V_cyl);
fprintf('Computed ni: %.3e m^-3\n', ni);
fprintf('Power density: %.3e J/m^3/s (%.3e eV/m^3/s)\n', power_density, power_density_eV);
fprintf('Total power:    %.3e J/s (%.3e eV/s)\n', total_power, total_power_eV_per_s);
fprintf('Photon density: %.3e photons/m^3/s\n', photon_density);
fprintf('Total photons:  %.3e photons/s\n', total_photons);

% Plot Emissivity
figure;
plot(omega_keV, eps_omega, 'k', 'LineWidth', 1.5);
xlabel('Photon Energy (keV)'); ylabel('\epsilon_X(\omega) [W\cdot m^{-3}\cdot Hz^{-1}]');
title('Bremsstrahlung Spectral Emissivity'); grid on;