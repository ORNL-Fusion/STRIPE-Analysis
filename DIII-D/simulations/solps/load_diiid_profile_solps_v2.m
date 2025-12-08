% Compute normalized flux (\psi_N) and plot SOLPS & experimental n_e profiles at Z=0.27

%% Clear and Load
close all; clear all; clc;

% Load SOLPS extrapolated data and EFIT g
load('extrapolated_data_200882.mat');  % Should contain: X, Y, val_ne, val_Te, etc.
% read_efit_data;                        % Load EFIT structure 'g'

% Replace NaNs in SOLPS ne
val_ne(isnan(val_ne)) = 0;

%% Define Z = 0.27 plane and get corresponding R vector
Z_target = 0.27;
Z_vec    = Y(:,1);
[~, iZ]   = min(abs(Z_vec - Z_target));

R_slice = X(iZ, :);
Z_slice = Z_target * ones(size(R_slice));

%% Plot limiter and SOLPS radial line & query points
figure;
plot(g.lim(1,:), g.lim(2,:), 'k', 'LineWidth', 1.5); hold on;
plot(R_slice, Z_slice, 'ro-', 'LineWidth', 1.2, 'DisplayName', 'SOLPS Slice');

% Define query points
R_query = [2.3, 2.296, 2.293, 2.289, 2.285, 2.281, 2.278, ...
           2.274, 2.27, 2.267, 2.263, 2.259, 2.256, 2.252, ...
           2.248, 2.244, 2.241, 2.237];
Z_query = Z_target * ones(size(R_query));

% Overlay query points on limiter cross-section
plot(R_query, Z_query, 'gs', 'MarkerSize', 8, 'MarkerFaceColor', 'g', ...
     'DisplayName', 'Selected R Points');

xlabel('R [m]'); ylabel('Z [m]');
title('Limiter & SOLPS Slice at Z = 0.27 m');
legend('Location','best');
axis equal; grid on;

%% Compute psi_N from EFIT for full slice and query points
[psi_vals, ierr] = calc_psiN(g, R_slice, Z_slice, []);
psiN_full = psi_vals;

% Interpolate psi_N at the query R positions
psiN_query = interp1(R_slice, psiN_full, R_query, 'linear');

%% Interpolate n_e from SOLPS along full slice and at query points
ne_full = interp2(X, Y, val_ne, R_slice, Z_slice, 'linear', NaN);
ne_query = interp2(X, Y, val_ne, R_query, Z_query, 'linear', NaN);

%% Load experimental He-beam data
data_He = load('he_beam_extracted_ne.mat');  % Contains fields: psiN_exp, ne_exp, ne_err, R_exp
psiN_vec = [1.14985542, 1.13730433, 1.12475158, 1.11219152, 1.09961843, ...
            1.08702662, 1.07441039, 1.06176404, 1.04908357, 1.03637852, ...
            1.02366495, 1.01095895, 0.99827663, 0.98563409, 0.97304743, ...
            0.96053146, 0.94809242, 0.93573293];

%% Plot n_e vs psi_N (SOLPS, He-beam, and query markers)
figure;
% plot(psiN_full, ne_full, 'b-', 'LineWidth', 1.5, 'DisplayName', 'SOLPS'); 
hold on;
errorbar(psiN_query, data_He.ne_exp./1e20, data_He.ne_err./1e20, 'rs--', ...
         'LineWidth', 1.5, 'DisplayName', 'He Beam Exp');

% Overlay query-point values
plot(psiN_query, ne_query./1e20, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'y', ...
     'DisplayName', 'SOLPS at Query R');

xlabel('\\psi_N', 'Interpreter','tex');
ylabel('n_e [10^{20} m^{-3}]', 'Interpreter','tex');
title('n_e vs. \\psi_N at Z = 0.27 m');
legend('Location','best');
grid on;
