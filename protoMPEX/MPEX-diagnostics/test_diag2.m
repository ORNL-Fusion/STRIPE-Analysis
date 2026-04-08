%% Bremsstrahlung from Kappa-Distributed Electrons in MPEX-Like Plasma Cylinder

clear; clc;

%% Physical Constants
e     = 1.602e-19;         % Elementary charge [C]
me    = 9.109e-31;         % Electron mass [kg]
c     = 3e8;               % Speed of light [m/s]
hbar  = 1.055e-34;         % Reduced Planck constant [J·s]
eps0  = 8.854e-12;         % Vacuum permittivity [F/m]
alpha = 1/137;             % Fine-structure constant
r_e   = e^2 / (4*pi*eps0*me*c^2);  % Classical electron radius [m]

%% Cylinder Geometry (MPEX-like)
r_max   = 0.06;          % Radius [m]
z_start = 3;             % Start z [m]
z_end   = 8;             % End z [m]
L       = z_end - z_start;
V_cyl   = pi * r_max^2 * L;   % Volume [m^3]

%% Kappa Velocity Distribution
kappa = 3;
Te_eV = 10;                % Electron temperature [eV]
Te_J  = Te_eV * e;             % Convert to Joules
v_th  = sqrt(2*Te_J/me);       % Thermal speed [m/s]

Nv = 200;
v_par  = linspace(-5*v_th, 5*v_th, Nv);
v_perp = linspace(0, 5*v_th, Nv);
[Vpar, Vperp] = meshgrid(v_par, v_perp);
V2 = Vpar.^2 + Vperp.^2;
f2D = (1 + V2/(kappa*v_th^2)).^(-(kappa+1));

dv_par  = v_par(2) - v_par(1);
dv_perp = v_perp(2) - v_perp(1);
norm2D  = sum(sum(f2D .* (2*pi.*Vperp))) * dv_par * dv_perp;
f2D     = f2D / norm2D;

%% Particle Sampling
N_particles = 1e6;
v_par_particles  = (2*rand(N_particles,1) - 1) * 5 * v_th;
v_perp_particles = rand(N_particles,1) * 5 * v_th;
v_mag = sqrt(v_par_particles.^2 + v_perp_particles.^2);
E_particles = 0.5 * me * v_mag.^2;

%% Plasma Parameters (Fixed realistic density)
ni = 1e19;   % Ion density [m^-3]

%% Photon Energy Grid [1–100 keV]
K_eV = linspace(0.01, 100, 2000);
K = K_eV * e;  % Convert to Joules

%% Bremsstrahlung Emissivity Loop
epsilon_K = zeros(size(K));
Z = 1;  % Deuterium

for i = 1:length(K)
    Ki = K(i);
    valid = E_particles > Ki;
    if any(valid)
        E_valid = E_particles(valid);
        v_valid = v_mag(valid);
        sigma = (16/3) * alpha * Z^2 * r_e^2 .* log(2*E_valid/Ki) ./ Ki;  % Born approx
        epsilon_K(i) = ni * mean(sigma .* v_valid);  % treat particles as normalized f(v)
    end
end

%% Total Bremsstrahlung Power Density
total_power = trapz(K, epsilon_K);  % [W/m^3]
fprintf('→ Total Bremsstrahlung Power Density: %.3e W/m^3\n', total_power);

% Integrated total bremsstrahlung power [W]
P_total = total_power * V_cyl;
fprintf('→ Integrated Bremsstrahlung Power: %.3f W\n', P_total);

%% Plot: Emissivity Spectrum
figure;
semilogy(K_eV, epsilon_K, 'b', 'LineWidth', 1.5);
xlabel('Photon Energy [keV]');
ylabel('\epsilon_X(K) [W/m^3/keV]');
title('Bremsstrahlung Spectrum from Kappa-Distributed Plasma (n_i = 10^{19} m^{-3})');
grid on; box on;