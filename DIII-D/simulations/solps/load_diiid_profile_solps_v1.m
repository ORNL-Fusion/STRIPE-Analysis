%% Clear and Load
close all; clear all; clc;

% Load SOLPS extrapolated data and EFIT
load("extrapolated_data_200882.mat");  % Should contain: X, Y, val_ne, val_Te, etc.
% Ensure EFIT g structure is in workspace
% e.g., read_efit_data;

% Replace NaNs in val_ne
val_ne(isnan(val_ne)) = 0;

%% Define Z = 0.27 plane and get corresponding R vector
Z_target = 0.27;
Z_vec = Y(:,1);
[~, iZ] = min(abs(Z_vec - Z_target));

R_slice = X(iZ, :);
Z_slice = Z_target * ones(size(R_slice));

% Optional: check that these points lie inside the limiter
figure;
plot(g.lim(1,:), g.lim(2,:), 'k', 'LineWidth', 1.5); hold on;
plot(R_slice, Z_slice, 'ro');
xlabel('R [m]'); ylabel('Z [m]');
title('Z = 0.27 m Radial Line vs Limiter');
axis equal; grid on;

%% Compute psi_N from EFIT for this line
[psi_vals, ierrs] = calc_psiN(g, R_slice, Z_slice, []);
psiN_vals = psi_vals;

%% Interpolate n_e from SOLPS along this line
ne_SOLPS = interp2(X, Y, val_ne, R_slice, Z_slice, 'linear', NaN);

%% Load experimental He beam data (psiN_exp, ne_exp)
data_He=load('he_beam_extracted_ne.mat');  % should load psiN_exp, ne_exp

% %% Save result (optional)
% save('psiN_comparison_profile.mat', 'R_slice', 'Z_slice', ...
%      'psiN_vals', 'ne_SOLPS', 'psiN_exp', 'ne_exp');

psiN_vec = [1.14985542, 1.13730433, 1.12475158, 1.11219152, 1.09961843, ...
            1.08702662, 1.07441039, 1.06176404, 1.04908357, 1.03637852, ...
            1.02366495, 1.01095895, 0.99827663, 0.98563409, 0.97304743, ...
            0.96053146, 0.94809242, 0.93573293];

%% Plot comparison
figure;
plot(psiN_vals, ne_SOLPS, 'b-o', 'LineWidth', 1.5); hold on;
errorbar(psiN_vec, data_He.ne_exp, data_He.ne_err, 'rs--', ...
    'LineWidth', 1.5, 'DisplayName', 'He Beam Exp');  % 5% error example
xlabel('\psi_N', 'Interpreter', 'tex');
ylabel('n_e [10^{20} m^{-3}]', 'Interpreter', 'tex');
title('n_e vs. \psi_N at Z = 0.27 m');
legend('SOLPS', 'He Beam Exp');
grid on;

%%

%% Define Target R and Z = 0.27
R_query = [2.3, 2.296, 2.293, 2.289, 2.285, 2.281, 2.278, ...
           2.274, 2.27, 2.267, 2.263, 2.259, 2.256, 2.252, ...
           2.248, 2.244, 2.241, 2.237];
Z_query = 0.27 * ones(size(R_query));

%% Interpolate ne from SOLPS at these R-Z points
ne_interp = interp2(X, Y, val_ne, R_query, Z_query, 'linear', NaN);

figure;
plot(R_query, ne_interp, 'b-o', 'LineWidth', 1.5); hold on;
errorbar(data_He.R_exp, data_He.ne_exp, data_He.ne_err, 'rs--', ...
    'LineWidth', 1.5, 'DisplayName', 'He Beam Exp');  % 5% error example
xlabel('\psi_N', 'Interpreter', 'tex');
ylabel('n_e [10^{20} m^{-3}]', 'Interpreter', 'tex');
title('n_e vs. \psi_N at Z = 0.27 m');
legend('SOLPS', 'He Beam Exp');
grid on; 