clc; clear; close all;

% Load input data (Ensure this file contains Te_CX, ne_CX, CX_rates_full)
load('ChargeExchangeData.mat');

% Ensure column vector format
Te_CX = Te_CX(:);  % Convert to column vector
ne_CX = ne_CX(:);  % Convert to column vector
CX_rates_full = CX_rates_full';  % Ensure proper orientation

% Sorting the Data (Critical for Interpolation)
[Te_CX, Te_sortIdx] = sort(Te_CX, 'ascend');
[ne_CX, ne_sortIdx] = sort(ne_CX, 'ascend');
CX_rates_full = CX_rates_full(Te_sortIdx, ne_sortIdx); % Align matrix

% Debug: Check uniqueness
if numel(unique(Te_CX)) < numel(Te_CX)
    error("🚨 Duplicate values detected in Te_CX. Ensure unique, sorted data.");
end
if numel(unique(ne_CX)) < numel(ne_CX)
    error("🚨 Duplicate values detected in ne_CX. Ensure unique, sorted data.");
end

% Define time span and initial conditions
tspan = [0 10];  % Example time span
f0 = [0; 0.0026; 0.2251; 0.7420; 0.0304];  % Initial values

% ODE solver options
opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);

% Run ODE solver
[tSol, fSol] = ode15s(@(t, f) ionBalanceODE_Velocity_CX(t, f, Te_CX, ne_CX, CX_rates_full), tspan, f0, opts);

% Plot results
figure;
plot(tSol, fSol);
xlabel('Time');
ylabel('Solution Variables');
title('Charge Balance Simulation');
legend('f_1', 'f_2', 'f_3', 'f_4', 'f_5');
grid on;

disp("✅ Simulation Complete.");

disp("🔍 Checking Data Integrity...");
fprintf("Te_CX range: [%.2f, %.2f]\n", min(Te_CX), max(Te_CX));
fprintf("ne_CX range: [%.2f, %.2f]\n", min(ne_CX), max(ne_CX));

% Example test values (Update these if needed)
Te_test = 10.00; 
log_ne_test = 3.10; 

% Find closest match in dataset
[~, idx_Te] = min(abs(Te_CX - Te_test));
[~, idx_ne] = min(abs(ne_CX - log_ne_test));
fprintf("Nearest data points: Te_CX(%.2f) log_ne(%.2f)\n", Te_CX(idx_Te), ne_CX(idx_ne));

% Display corresponding CX_rate value
fprintf("CX_rate at nearest point: %.4e\n", CX_rates_full(idx_Te, idx_ne));

%% Function for ODE system with CX interpolation fix
function dfdt = ionBalanceODE_Velocity_CX(t, f, Te_CX, ne_CX, CX_rates_full)
    % Define the electron temperature and density (example functions)
    Te = 10 + 0.5 * t;  % Example temperature variation
    log_ne = log10(1e13 + 1e12 * sin(t));  % Example density variation
    
    % Clamp values to avoid extrapolation issues
    Te = max(min(Te, max(Te_CX)), min(Te_CX));
    log_ne = max(min(log_ne, max(ne_CX)), min(ne_CX));

    % Debug: Print selected interpolation points
    fprintf("🔎 Interpolating at Te = %.2f, log_ne = %.2f\n", Te, log_ne);
    
    % Perform interpolation (Fixing missing value issues)
    CX_rate = interp2(Te_CX, ne_CX', CX_rates_full, Te, log_ne, 'linear', NaN);

    % Debug: Check if interpolation failed
    if isnan(CX_rate) || isempty(CX_rate)
        fprintf('⚠️ NaN in CX_rate at Te = %.2f, log_ne = %.2f\n', Te, log_ne);
        
        % Nearest-neighbor fallback
        CX_rate = interp2(Te_CX, ne_CX', CX_rates_full, Te, log_ne, 'nearest', 0);
        if isnan(CX_rate)
            CX_rate = 0;  % Final fallback
        end
    end

    % Ensure dfdt has correct size
    dfdt = zeros(size(f));

    % Example ODE system (Replace with actual physics model)
    dfdt(1) = -CX_rate * f(1);
    dfdt(2) = CX_rate * (f(1) - f(2));
    dfdt(3) = CX_rate * (f(2) - f(3));
    dfdt(4) = CX_rate * (f(3) - f(4));
    dfdt(5) = CX_rate * f(4) - 0.1 * f(5);
    
    % Ensure output is a COLUMN vector
    dfdt = dfdt(:);
end