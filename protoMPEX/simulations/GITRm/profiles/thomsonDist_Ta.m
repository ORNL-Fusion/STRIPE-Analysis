%% Thomson half-range velocity distribution for Ta_z in MPEX
% Upstream boundary condition with flux profile, Tz = 5 eV, v_par,z = 1000 m/s

clear; clc;

%% Physical constants
eV   = 1.602176634e-19;      % [J/eV]
kB   = 1.380649e-23;         % [J/K]
amu  = 1.66053906660e-27;    % [kg]

% Tantalum ion mass (A ~ 181)
A_Ta = 181;
mTa  = A_Ta * amu;           % [kg]

%% Parameters from slide
Tz_eV  = 5;                  % [eV]
Vpar   = 1000;               % [m/s]
Gamma0 = 1e20;               % [m^-2 s^-1] central flux
r0     = 0.02;               % 2 cm in meters

% Charge state fractions (Ta2+, Ta3+, Ta4+)
Z_list   = [2 3 4];
frac_Z   = [0.08 0.62 0.30]; % must sum to 1

%% Radial grid and flux profile
Nr = 100;
r  = linspace(0, 0.04, Nr);  % 0 to 4 cm
Gamma_tot = Gamma0 * exp( - (r./r0).^12 );   % total Ta flux vs r

% Per-charge-state fluxes
Gamma_z = Gamma_tot(:) .* frac_Z;   % Nr x 3 (implicit expansion in newer MATLAB)

%% Thermal speed (for setting v-grid)
Tz_J = Tz_eV * eV;                 % [J]
v_th = sqrt(2 * Tz_J / mTa);       % [m/s]

%% Velocity grids
% Full symmetric grid for temperature check
Nv_full = 2001;
vmax    = Vpar + 5*v_th;
v_full  = linspace(-vmax, vmax, Nv_full);

% Positive velocities only for boundary distribution
Nv  = 800;
v   = linspace(0, vmax, Nv);       % v >= 0 only

%% --- 1) Build underlying full drifting Maxwellian (for T check) ---

f_full = (1/(sqrt(pi)*v_th)) * exp( -((v_full - Vpar)/v_th).^2 );  % 1D f(v)

% Normalization check (should be ~1)
n_full = trapz(v_full, f_full);

% Moments -> recover T
v_mean  = trapz(v_full, v_full .* f_full) / n_full;
v2_mean = trapz(v_full, (v_full.^2) .* f_full) / n_full;
var_v   = v2_mean - v_mean^2;      % <(v - <v>)^2>

T_rec_eV = mTa * var_v / (eV);  % from <(v - <v>)^2> = kT/m

fprintf('--- Temperature consistency check ---\n');
fprintf('Input   T_z       = %.3f eV\n', Tz_eV);
fprintf('Recovered T_z,1D  = %.3f eV (from full f(v))\n', T_rec_eV);
fprintf('Input   V_par,z   = %.3f m/s\n', Vpar);
fprintf('Mean v from f(v)  = %.3f m/s\n\n', v_mean);

%% --- 2) Half-range Thomson distribution with given flux profile ---

% Shape in v for v >= 0: same drifting Maxwellian, truncated
f_shape = (1/(sqrt(pi)*v_th)) * exp( -((v - Vpar)/v_th).^2 );  % unscaled shape

% Flux carried by this shape (without amplitude) in +v
J_shape = trapz(v, v .* f_shape);   % [m/s * (1/v)] ~ 1

% For each radius, choose amplitude A(r) s.t.
%   Gamma_tot(r) = ∫_0^∞ v * f_T(v,r) dv = A(r) * J_shape
A_r = Gamma_tot / J_shape;          % [m^-2 s^-1] / [dimensionless] -> scaling

% Build f_T(v,r) for total Ta flux (half-range Thomson)
f_T = A_r(:) .* f_shape;            % Nr x Nv (implicit expansion)

%% Example: choose Ta3+ only (dominant charge state)
Gamma_Ta3 = Gamma_tot * frac_Z(2);  % 62% of total
A_r_Ta3   = Gamma_Ta3 / J_shape;
f_T_Ta3   = A_r_Ta3(:) .* f_shape;  % Nr x Nv

%% --- 3) Plots ---

% (a) Total Ta Thomson distribution f_T(v,r)
figure;
imagesc(v, r*100, f_T);  % r in cm
set(gca, 'YDir', 'normal');
xlabel('v_{||,z} [m/s]');
ylabel('r [cm]');
title(sprintf('Thomson half-range f_T(v_{||,z}, r) for Ta (T_z = %.1f eV, V_{||,z} = %g m/s)', ...
              Tz_eV, Vpar));
cb = colorbar;
ylabel(cb, 'f_T(v_{||,z}, r)  [m^{-4} s]');

% (b) Ta3+ distribution at a few radii
figure;
idx_center = round(Nr/2);   % mid-radius
idx_edge   = Nr;            % outer edge

plot(v, f_T_Ta3(1,:),      'LineWidth', 2); hold on;
plot(v, f_T_Ta3(idx_center,:), 'LineWidth', 2);
plot(v, f_T_Ta3(idx_edge,:),   'LineWidth', 2);
grid on;
xlabel('v_{||,3+} [m/s]');
ylabel('f_{T,3+}(v_{||})  [m^{-4} s]');
title('Ta^{3+} Thomson half-range distribution at selected radii');
legend('r = 0', ...
       sprintf('r = %.1f cm', r(idx_center)*100), ...
       sprintf('r = %.1f cm', r(idx_edge)*100), ...
       'Location', 'best');

% (c) Plot the imposed flux profile to verify
figure;
plot(r*100, Gamma_tot, 'LineWidth', 2);
grid on;
xlabel('r [cm]');
ylabel('\Gamma_{Ta} [m^{-2} s^{-1}]');
title('\Gamma_z(r) = 10^{20} exp[-(r / 2cm)^{12}]');

%% ------------------------------------------------------------
%  1D Plots: Full drifting Maxwellian + Half-range Thomson
% ------------------------------------------------------------

% --- Full symmetric drifting Maxwellian (for reference)
f_full = (1/(sqrt(pi)*v_th)) * exp( -((v_full - Vpar)/v_th).^2 );

% --- Half-range Thomson f_T at the center radius (most relevant)
i_mid = round(Nr/2);
f_half = f_T(i_mid,:);     % f(v >= 0, r_mid)


%% (1) Plot: FULL Maxwellian including drift (negative + positive v)
figure;
plot(v_full, f_full, 'LineWidth', 2); hold on;
xline(Vpar, '--k', 'LineWidth', 2);
grid on;
xlabel('v_{||} [m/s]');
ylabel('f(v_{||})');
title(sprintf('Full drifting Maxwellian (Ta_z), T_z = %.1f eV, V_{||} = %g m/s', ...
               Tz_eV, Vpar));
legend('Full Maxwellian', 'Drift velocity V_{||}');


%% (2) Plot: HALF-RANGE Thomson distribution (v ≥ 0) at r_mid
figure;
plot(v, f_half, 'LineWidth', 2); hold on;
xline(Vpar, '--k', 'LineWidth', 2);
grid on;
xlabel('v_{||} [m/s]');
ylabel('f_T(v_{||})');
title(sprintf('Thomson Half-Range Distribution (Ta_z) at r = %.2f cm', ...
               r(i_mid)*100));
legend('Half-range Thomson f_T(v)', 'Drift velocity V_{||}', ...
        'Location', 'best');