
% Hybrid MATLAB Script: X-ray Production from 2D Kappa Distribution in a Uniform Plasma Cylinder
% Combines:
% 1) Grid-based 2D kappa velocity distribution for diagnostics
% 2) Particle-based initialization in cylinder
% 3) ni computed from particle count and volume
% 4) Spectral emissivity ε_X(ω), total power, photon rates
% 5) Full visualizations and unit conversion

clear; clc; close all;

%% Physical Constants
e    = 1.602e-19;
hbar = 1.055e-34;
c    = 3e8;
me   = 9.109e-31;
alpha= 1/137;
eps0 = 8.854e-12;
Z    = 1;
r_e  = e^2/(4*pi*eps0*me*c^2);

%% Cylinder Geometry
r_max   = 0.06;
z_start = 3;
z_end   = 8;
L       = z_end - z_start;
V_cyl   = pi * r_max^2 * L;

% Plot Cylinder Geometry
figure('Name','Cylinder Geometry','Position',[100 100 600 400]);
[Xc,Yc,Zc] = cylinder(r_max,100);
Zc = Zc * L + z_start;
surf(Xc,Yc,Zc,'FaceAlpha',0.2,'EdgeColor','none'); hold on;
theta = linspace(0,2*pi,200);
x_c = r_max*cos(theta); y_c = r_max*sin(theta);
plot3(x_c,y_c,z_start*ones(size(theta)),'b','LineWidth',1);
plot3(x_c,y_c,z_end*ones(size(theta)),'b','LineWidth',1);
xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');
title('Uniform Plasma Cylinder (r=0.06 m, z=3–8 m)');
axis equal; grid on; view(3);

%% 2D Kappa Velocity Distribution (Grid-Based)
kappa = 3;
T_e_k = 100 * e;
v_th  = sqrt(2*T_e_k/me);
Nv    = 200;
v_par  = linspace(-5*v_th,5*v_th,Nv);
v_perp = linspace(0,5*v_th,Nv);
[Vpar,Vperp] = meshgrid(v_par,v_perp);
V2 = Vpar.^2 + Vperp.^2;
f2D = (1 + V2./(kappa*v_th^2)).^(-(kappa+1));
dv_par  = v_par(2)-v_par(1);
dv_perp = v_perp(2)-v_perp(1);
norm2D  = sum(sum(f2D.*(2*pi.*Vperp))) * dv_par * dv_perp;
f2D     = f2D / norm2D;

% Plot 2D Distribution
figure('Name','2D Kappa Velocity Distribution','Position',[100 100 600 400]);
contourf(v_par/1e6, v_perp/1e6, f2D,50,'LineColor','none');
xlabel('v_{||} (10^6 m/s)'); ylabel('v_{\perp} (10^6 m/s)');
title(sprintf('2D Kappa Distribution (\kappa=%.1f)',kappa));
colorbar; grid on;

%% Sample Particles in Cylinder
N_particles = 1e6;
r_particles = r_max * sqrt(rand(N_particles,1));
theta_particles = 2*pi * rand(N_particles,1);
z_particles = z_start + L * rand(N_particles,1);
x_particles = r_particles .* cos(theta_particles);
y_particles = r_particles .* sin(theta_particles);

figure('Name','Particle Positions');
scatter3(x_particles, y_particles, z_particles, 1, '.'); axis equal;
xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]');
title('Uniform Particle Distribution in Cylinder'); grid on;

%% Assign Velocities from 2D Kappa Distribution
v_par_particles  = (2*rand(N_particles,1)-1) * 5*v_th;
v_perp_particles = rand(N_particles,1) * 5*v_th;
V2_particles = v_par_particles.^2 + v_perp_particles.^2;
f_kappa_particles = (1 + V2_particles./(kappa*v_th^2)).^(-(kappa+1));

% Plot sampled velocities
figure('Name','Sampled Velocities');
scatter(v_par_particles/1e6, v_perp_particles/1e6, 1, f_kappa_particles, '.');
xlabel('v_{||} (10^6 m/s)'); ylabel('v_{\perp} (10^6 m/s)');
title('Sampled Velocities from Kappa Distribution'); colorbar;

%% Compute ni from Particle Count
ni = N_particles / V_cyl;

% %% Bremsstrahlung Emissivity (Born Approximation)
% % We compute εₓ(ω) = nᵢ ∫ f(v) σ(ω,v) v d³v using
% %
% %   σ(ω,v) = (16/3)·Z²·α·rₑ²·(mₑc²/ω)·ln[2E/(mₑc²)+1],
% %
% % where E = ½ mₑ v².
% %
% % References:
% %  [1] H. A. Bethe & W. Heitler, Proc. Roy. Soc. A146 (1934) 83.
% %  [2] H. J. Kunze, Introduction to Plasma Spectroscopy, Springer (2009), Ch. 4.
% 
omega_keV = linspace(0.1,20,200);
omega_J   = omega_keV*1e3*e;
eps_omega = zeros(size(omega_J));
Vmag = sqrt(V2_particles);
% 
% for j = 1:numel(omega_J)
%     omega = omega_J(j);
%     E_particles = 0.5 * me * V2_particles;
%     log_term = log(2*E_particles/(me*c^2) + 1);
%     sigma = (16/3) * Z^2 * alpha * r_e^2 .* (me*c^2/omega) .* log_term;
%     emissivity = f_kappa_particles .* sigma .* Vmag;
%     eps_omega(j) = ni * mean(emissivity);
% end

%% ——————————————————————————————
%% Compute Free–Free Spectral Emissivity (Kunze Eq. 6.81)
%   ε_ff,λ(λ) = 4.108e-46 · n_i·n_e·Z^2 · sqrt(E_R/(kB·T_e)) · (1/λ^2)
%               · exp[–h·c/(λ·kB·T_e)] · G_ff(T_e,λ)
%
%   Units: ε_ff,λ [W·m⁻³·m⁻¹], λ in metres, n in m⁻³, T_e in K
%
% Constants:
E_R   = 13.605693 * e;           % Rydberg energy [J]
kB    = 1.380649e-23;            % Boltzmann constant [J/K]
% T_e   = 100 * e / kB;            % electron temperature [K] (100 eV)

% ———————————————————————————————————
%  Compute Te from the grid-based 2D kappa PDF

% second moment: <v^2> = ∫(V2) f2D * 2π v_perp dv_perp dv_par
mean_V2 = sum( sum( V2 .* f2D .* (2*pi.*Vperp) ) ) * dv_par * dv_perp;

% invert to get T_e in Kelvin, then convert to eV
Te_K    = me * mean_V2 / (3*kB);
Te_eV   = Te_K * kB / e;

% fprintf('→ Derived electron temperature: %.2f eV (%.0f K)\\n', Te_eV, Te_K);

% now redefine your thermal speed from this Te:
v_th = sqrt(2 * Te_eV * e / me);
fprintf('→ Thermal speed v_{th} = %.2e m/s\\n', v_th);
% ———————————————————————————————————
T_e   = Te_K;
n_i   = ni;                      % ion density [m^-3]
n_e   = ni;                      % electron density ≈ ion density
% G_ff  = ones(size(omega_keV));   % placeholder Gaunt factor (assume 1)

% ——————————————————————————————
% Compute Gaunt factor via Sutherland analytic approximation
%   G_ff(u,γ²) ≈ (√3/π)·ln[4/(γ·u)], clamped to [0.8,1.5]
%
% where γ² = Z²·E_R/(kB·T_e) and u = hν/(kB·T_e)

% Define wavelength grid (convert keV → λ)
omega_J    = omega_keV*1e3*e;      % photon energy [J]
lambda_m   = (hbar .* c) ./ omega_J * 2*pi;  % λ = hc/E

%% ——————————————————————————————
%% Compute Bremsstrahlung Spectral Emissivity with Gaunt Factor
gamma = 0.577215664901532;   % Euler–Mascheroni constant

eps_omega = zeros(size(omega_J));
Vmag      = sqrt(V2_particles);  % particle speeds

for j = 1:numel(omega_J)
    omega = omega_J(j);             % photon energy [J]
    
    % 1) dimensionless parameters for Gaunt factor
    eta_Z  = Z .* e^2 ./ (hbar .* Vmag);
    eta_nu = hbar .* omega ./ (2 * me .* Vmag.^2);
    
    % 2) rough Gaunt factor (Wikipedia approximation)
    g_ff = (sqrt(3)/pi) .* log( 1 ./ (eta_nu .* max(1, exp(gamma) .* eta_Z)) );
    g_ff(g_ff < 1) = 1;              % enforce g_ff ≥ 1
    
    % 3) Born‐approximation cross section term
    E_part   = 0.5 * me .* V2_particles;                 % kinetic energy
    log_term = log(2*E_part/(me*c^2) + 1);               % ln[2E/(m_ec^2)+1]
    sigma    = (16/3) * Z^2 * alpha * r_e^2 ...          % [m^2/Hz]
               .* (me*c^2./omega) .* log_term;
    
    % 4) per‐particle emissivity σ·v·g_ff
    emissivity = sigma .* Vmag .* g_ff;
    
    % 5) ensemble‐average and multiply by density
    eps_omega(j) = ni * mean(emissivity);
end

%% — integrate for power & photon rates as before
power_density  = trapz(omega_J, eps_omega .* omega_J);   % J/m^3/s
phi_omega      = eps_omega ./ (hbar .* omega_J);        % ph/m^3/s/Hz
photon_density = trapz(omega_J, phi_omega);              % ph/m^3/s
total_power    = power_density * V_cyl;                  % J/s
total_photons  = photon_density * V_cyl;                 % ph/s

%% — convert & display
power_density_eV    = power_density / e;                  
total_power_eV_per_s= total_power   / e;

fprintf('\\n=== X-ray Emission (with Gaunt Factor) ===\\n');
fprintf('Total power: %.3e J/s (%.3e eV/s)\\n', total_power, total_power_eV_per_s);
fprintf('Total photons: %.3e photons/s\\n', total_photons);

%% Total Power and Photon Rate
power_density  = trapz(omega_J, eps_omega);                 
phi_omega      = eps_omega ./ (hbar .* omega_J);            
photon_density = trapz(omega_J, phi_omega);                 
total_power    = power_density * V_cyl;
total_photons  = photon_density * V_cyl;

%% Convert to eV
power_density_eV = power_density / e;
total_power_eV_per_s = total_power / e;

%% Display Results
fprintf('\n=== X-ray Production in Cylinder ===\n');
fprintf('Cylinder volume: %.3f m^3\n', V_cyl);
fprintf('Computed ni: %.3e m^-3\n', ni);
fprintf('Power density: %.3e J/m^3/s (%.3e eV/m^3/s)\n', power_density, power_density_eV);
fprintf('Total power:    %.3e J/s (%.3e eV/s)\n', total_power, total_power_eV_per_s);
fprintf('Photon density: %.3e photons/m^3/s\n', photon_density);
fprintf('Total photons:  %.3e photons/s\n', total_photons);

%% Plot Spectral Emissivity
figure('Name','Spectral Emissivity');
plot(omega_keV, eps_omega,'k','LineWidth',1.5);
xlabel('Photon Energy (keV)');
ylabel('\epsilon_X(\omega)');
title('Bremsstrahlung Spectral Emissivity');
grid on;
