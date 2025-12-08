clc; clear; close all;

%% Parameters (Modify as needed)
num_samples = 1000; % Number of Monte Carlo samples
sigma_noise = 0.1;  % Standard deviation of noise (fraction of value)

% Load SOLPS Data
load("extrapolated_data.mat");

% Identify OMP index (middle of zS)
[~, ompIndex] = min(abs(zS));

% Extract relevant profiles
R = rS;
ne_obs = neS(:,ompIndex);

% Apply factor of 10 uncertainty bounds
ne_lower = ne_obs / 10;
ne_upper = ne_obs * 10;

% Select range R = 8 to 8.3
fit_indices = (R >= 8) & (R <= 8.3);
R_fit = R(fit_indices);
ne_fit = ne_obs(fit_indices);
ne_upper_fit = ne_upper(fit_indices);
ne_lower_fit = ne_lower(fit_indices);

% Ensure column vectors
R_fit = R_fit(:);
ne_fit = ne_fit(:);
ne_upper_fit = ne_upper_fit(:);
ne_lower_fit = ne_lower_fit(:);

% Add Gaussian noise to simulate measurement uncertainty
ne_noisy = ne_fit .* (1 + sigma_noise * randn(size(ne_fit)));

%% Curve fitting (Exponential Decay Model)
fit_func = @(b, x) b(1) * exp(- (x - R_fit(1)) / b(2));
beta0 = [max(ne_fit), 0.02]; % Initial guess for fitting
beta_hat = nlinfit(R_fit, ne_noisy, fit_func, beta0);

%% Monte Carlo Sampling for Uncertainty Quantification
beta_samples = zeros(num_samples, 2);
for i = 1:num_samples
    ne_perturbed = ne_noisy + sigma_noise * ne_noisy .* randn(size(ne_noisy));
    beta_samples(i, :) = nlinfit(R_fit, ne_perturbed, fit_func, beta0);
end

% Compute mean and standard deviation of fitted parameters
beta_mean = mean(beta_samples, 1);
beta_std = std(beta_samples, [], 1);

% Generate density profiles from sampled parameters
ne_profiles = arrayfun(@(i) fit_func(beta_samples(i, :), R_fit), 1:num_samples, 'UniformOutput', false);
ne_profiles = cell2mat(ne_profiles');

% Compute confidence intervals (2 sigma ~95%)
ne_mean = mean(ne_profiles, 1);
ne_std = std(ne_profiles, 0, 1);
ne_upper_fit = ne_mean(:) + 2 * ne_std(:);
ne_lower_fit = ne_mean(:) - 2 * ne_std(:);

%% Ensure vectors have the same size for fill()
if length(R_fit) ~= length(ne_upper_fit) || length(R_fit) ~= length(ne_lower_fit)
    error('Mismatch in vector sizes: Check dimensions of R_fit, ne_upper_fit, and ne_lower_fit.');
end

%% Plot Results with Shaded Uncertainty Region
figure; hold on;
fill([R_fit; flipud(R_fit)], [ne_upper_fit; flipud(ne_lower_fit)], [0.8 0.8 1], ...
    'EdgeColor', 'none', 'FaceAlpha', 0.5); % Confidence band
plot(R_fit, ne_noisy, 'ko', 'MarkerFaceColor', 'k', 'DisplayName', 'Observed Data'); % Data points
plot(R_fit, fit_func(beta_hat, R_fit), 'r-', 'LineWidth', 2, 'DisplayName', 'Best Fit'); % Best fit
plot(R_fit, ne_mean, 'b--', 'LineWidth', 1.5, 'DisplayName', 'Mean Profile'); % Mean profile
xlabel('Radial Position R [m]'); ylabel('Electron Density [m^{-3}]');
legend('Uncertainty Bound (95%)', 'Observed Data', 'Best Fit', 'Mean Profile');
title('SOL Electron Density Profile (R=8 to 8.3) with Uncertainty Quantification');
grid on;

%% Display Results
disp('Best-fit Parameters:');
disp(['n0 = ', num2str(beta_mean(1)), ' ± ', num2str(beta_std(1))]);
disp(['λ_n = ', num2str(beta_mean(2)), ' ± ', num2str(beta_std(2))]);