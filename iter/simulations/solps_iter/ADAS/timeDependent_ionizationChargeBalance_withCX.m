clear; clc; close all;

%=====================================================================
% Neon Ionization Balance with Charge Exchange
% (Initial condition: coronal (equilibrium) distribution fSS)
%---------------------------------------------------------------------
% This script:
%  1) Loads ADAS rate coefficients from a NetCDF file
%  2) Interpolates S, R, CX at the given Te, ne
%  3) Computes a self‐consistent "coronal" (equilibrium) fSS
%  4) Solves the time‐dependent ODE for f_z including charge exchange (CX),
%     using fSS as the initial condition
%  5) Plots the time evolution, and shows an inset comparing f_final vs fSS
%=====================================================================

%% 1) Load ADAS Data from NetCDF File
filename = 'ADAS_Rates_Ne.nc';
% filename = 'ADAS_Rates_O.nc';  % (uncomment if you want oxygen instead)

% Ionization data
IonizationData.Temp      = 10.^ncread(filename, 'gridTemperature_Ionization');
IonizationData.Density   = 10.^ncread(filename, 'gridDensity_Ionization');
IonizationData.RateCoeff = 10.^ncread(filename, 'IonizationRateCoeff');

% Recombination data
RecombinationData.Temp      = 10.^ncread(filename, 'gridTemperature_Recombination');
RecombinationData.Density   = 10.^ncread(filename, 'gridDensity_Recombine');
RecombinationData.RateCoeff = 10.^ncread(filename, 'RecombinationRateCoeff');

% Charge‐Exchange data
CXData.Temp      = 10.^ncread(filename, 'gridTemperature_CX');
CXData.Density   = 10.^ncread(filename, 'gridDensity_CX');
CXData.RateCoeff = 10.^ncread(filename, 'ChargeExchangeRateCoeff');

%% 2) Simulation and Plasma Parameters
Zmax = 10;           % Maximum neon charge state (0..10+)
Te   = 10;            % Electron temperature [eV]
ne   = 5e16;         % Electron density [m^-3]

n_n  = 0.02 * ne;    % Neutral density for CX (1% of n_e)

%% 3) Compute S, R, and CX Arrays via Interpolation
S  = zeros(Zmax+1, 1);   % S(z+1) = ionization rate out of charge state z
R  = zeros(Zmax+1, 1);   % R(z+1) = recombination rate out of charge state z
CX = zeros(Zmax+1, 1);   % CX(z+1) = charge‐exchange rate out of charge state z

% Start parallel pool if none exists
if isempty(gcp('nocreate'))
    parpool;
end

parfor j = 1:Zmax
    % Physical charge state z = j-1
    % Ionization rate S_{z} (for z=0..Zmax-1):
    S(j) = interpn( IonizationData.Density, IonizationData.Temp, ...
                    IonizationData.RateCoeff(:, :, j), ...
                    ne, Te, 'linear', 0 );
    % Recombination rate R_{z} (for z=1..Zmax):
    R(j+1) = interpn( RecombinationData.Density, RecombinationData.Temp, ...
                      RecombinationData.RateCoeff(:, :, j), ...
                      ne, Te, 'linear', 0 );
    % CX rate CX_{z} (for z=1..Zmax):
    if size(CXData.RateCoeff, 3) >= j
        CX(j+1) = interpn( CXData.Density, CXData.Temp, ...
                          CXData.RateCoeff(:, :, j), ...
                          ne, Te, 'linear', 0 );
    else
        CX(j+1) = 0;
    end

    % Prevent zeros (avoid divide-by-zero)
    S(j)    = max(S(j),    eps);
    R(j+1)  = max(R(j+1),  eps);
    CX(j+1) = max(CX(j+1), eps);
end

% By construction:
%   S(1..Zmax)  = S₀..S_{Zmax-1},  S(Zmax+1) = 0
%   R(1) = 0,  R(2..Zmax+1) = R₁..R_{Zmax}
%   CX(1) = 0, CX(2..Zmax+1) = CX₁..CX_{Zmax}

%% 4) Compute Self‐Consistent "Coronal" (Equilibrium) Distribution fSS
fSS = zeros(Zmax+1, 1);
fSS(1) = 1;  % Start recursion at f₀ = 1 (will renormalize)

for z = 1:Zmax
    % Ratio f_{z} / f_{z-1} = (S_{z-1} / R_{z}) * exp(-z/3)
    numerator   = S(z);        % S_{z-1}
    denominator = R(z+1);      % R_{z}
    if denominator > 0
        fSS(z+1) = fSS(z) * (numerator / denominator) * exp(-z/3);
    else
        fSS(z+1) = 0;
    end
end
% Normalize so that sum(fSS) = 1
fSS = fSS ./ sum(fSS);

disp('✅ Computed equilibrium (coronal) charge distribution fSS:');
disp(fSS);

%% 5) Solve the ODE System for Ionization Balance with CX
% Initial condition: coronal distribution fSS
tspan = [logspace(-8, -4, 10000), linspace(1e-4, 1e1, 100000)];  % [s]
opts  = odeset('RelTol', 1e-6, 'AbsTol', 1e-10);

[tSol, fSol] = ode23s(@(t, f) ionBalanceODE_CX(t, f, S, R, CX, ne, n_n, Zmax), ...
                      tspan, fSS, opts);

% Normalize each row so that sum(f)=1
fSol = fSol ./ sum(fSol, 2);

%% 6) Extract Final Charge State Distribution
f_final = fSol(end, :);

%% 7) Plot Time Evolution & Inset: Final vs fSS
figure;
hold on;
colorMap = jet(Zmax+1);
for i = 1:(Zmax+1)
    plot(tSol*1e3, fSol(:, i), 'LineWidth', 2, 'Color', colorMap(i,:));
end
xlabel('Time [ms]');
ylabel('Fractional Abundance');
title('Neon Ionization Balance with CX (Initial = Coronal fSS)');
legend(arrayfun(@(z) sprintf('Ne^{%d+}', z), 0:Zmax, 'UniformOutput', false), 'Location', 'best');
set(gca, 'XScale', 'log');
grid on;
hold off;

% Inset axes (compare final distribution vs fSS)
insetAxes = axes('Position', [0.55, 0.55, 0.35, 0.35]);
hold(insetAxes, 'on');

% Bar plot: final distribution
barHandle = bar(insetAxes, 0:Zmax, f_final, 'FaceColor', 'flat', 'FaceAlpha', 0.6);

% Line plot: coronal initial distribution fSS
plot(insetAxes, 0:Zmax, fSS, 'ko-', 'MarkerFaceColor', 'k', 'LineWidth', 2);

xticks(insetAxes, 0:Zmax);
xticklabels(insetAxes, arrayfun(@(z) sprintf('Ne^{%d+}', z), 0:Zmax, 'UniformOutput', false));
xlabel(insetAxes, 'Charge State');
ylabel(insetAxes, 'Fractional Abundance');
title(insetAxes, 'Final vs Coronal fSS');
legend(insetAxes, 'Final Ionization State', 'Coronal fSS', 'Location', 'best');
grid(insetAxes, 'on');

% Color each bar by charge state
for i = 1:(Zmax+1)
    barHandle.CData(i, :) = colorMap(i, :);
end

%% ODE FUNCTION: Ionization Balance with CX (No Velocity)
function dfdt = ionBalanceODE_CX(~, f, S, R, CX, ne, n_n, Zmax)
    % f: (Zmax+1)x1 vector, f(z+1) = fraction at charge state z
    % S, R, CX: length (Zmax+1) arrays with S(z+1)=S_z, R(z+1)=R_z, CX(z+1)=CX_z
    % ne: electron density; n_n: neutral density

    dfdt = zeros(Zmax+1, 1);

    % z = 0 (i = 1)
    dfdt(1) = - ne * S(1) * f(1) + ne * R(2) * f(2) ...   % Ionization loss + Recomb gain
              - n_n * CX(1) * f(1) + n_n * CX(2) * f(2);  % CX loss + CX gain

    % z = 1 .. Zmax-1 (i = 2 .. Zmax)
    for i = 2:Zmax
        dfdt(i) = ne * S(i) * f(i-1) ...                          % Ionization gain
                  - ne * ( S(i+1) + R(i) ) * f(i) ...              % Ionization + Recomb loss
                  + ne * R(i+1) * f(i+1) ...                       % Recomb gain
                  - n_n * CX(i) * f(i) + n_n * CX(i+1) * f(i+1);   % CX loss + CX gain
    end

    % z = Zmax (i = Zmax+1)
    i = Zmax + 1;
    dfdt(i) = ne * S(i-1) * f(i-1) ...        % Ionization gain from Zmax−1 → Zmax
              - ne * R(i)   * f(i) ...        % Recomb loss from Zmax → Zmax−1
              - n_n * CX(i) * f(i);           % CX loss from Zmax → Zmax−1
end