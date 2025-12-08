clear; clc; close all;

%% ✅ Load ADAS Data from NetCDF File
filename = 'ADAS_Rates_O.nc';

% Ionization
IonizationData.Temp = 10.^ncread(filename, 'gridTemperature_Ionization');
IonizationData.Density = 10.^ncread(filename, 'gridDensity_Ionization');
IonizationData.RateCoeff = 10.^ncread(filename, 'IonizationRateCoeff');

% Recombination
RecombinationData.Temp = 10.^ncread(filename, 'gridTemperature_Recombination');
RecombinationData.Density = 10.^ncread(filename, 'gridDensity_Recombine'); 
RecombinationData.RateCoeff = 10.^ncread(filename, 'RecombinationRateCoeff');

%% ✅ Simulation Parameters
Zmax = 8;            % O⁰ to O⁸⁺
Te = 6;              % Electron temperature [eV]
ne = 5e16;           % Electron density [m^-3]

% Time span [s]
tspan = [logspace(-8, -4, 10000), linspace(1e-4, 1e1, 100000)];

%% ✅ Interpolate Rate Coefficients (Ionization and Recombination Only)
if isempty(gcp('nocreate'))
    parpool;
end

S = zeros(Zmax+1, 1);  
R = zeros(Zmax+1, 1);
R(1) = 0;  % No recombination into neutrals

parfor j = 1:Zmax
    S(j) = interpn(IonizationData.Density, IonizationData.Temp, IonizationData.RateCoeff(:,:,j), ne, Te, 'linear', 0);
    R(j) = interpn(RecombinationData.Density, RecombinationData.Temp, RecombinationData.RateCoeff(:,:,j), ne, Te, 'linear', 0);

    % Avoid zero rates
    S(j) = max(S(j), eps);
    R(j) = max(R(j), eps);
end

%% **4) Compute Self-Consistent Initial Charge Distribution (`fSS`)**
% Initialize fSS with neutral dominance
fSS = zeros(Zmax+1, 1);
fSS(1) = 1 ;  % Start with all neutrals

% Compute recursive charge distribution with an enforced decay factor
for z = 1:Zmax
    if S(z) > 0 && R(z+1) > 0
        fSS(z+1) = fSS(z) * (S(z) / R(z+1)) * exp(-z/3);  % Exponential decay factor
    else
        fSS(z+1) = 0;  % Avoid division by zero
    end
end

% Normalize to ensure total sum = 1
fSS = fSS / sum(fSS);
disp('✅ Computed Equilibrium Initial Charge Distribution:');
disp(fSS);

% %% **4) Fully Stripped Initial Charge Distribution (`fSS`)**
% fSS = zeros(Zmax+1, 1);
% fSS(end) = 1;  % All ions initially in fully stripped state (O⁸⁺)
% disp('✅ Fully Stripped Initial Condition (O⁸⁺ only):');
% disp(fSS);
%% ✅ Solve ODE: Ionization + Recombination Only
opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-10);
[tSol, fSol] = ode23s(@(t, f) ionBalanceODE_no_CX(t, f, S, R, ne, Zmax), tspan, fSS, opts);
fSol = fSol ./ sum(fSol, 2);  % Normalize

%% ✅ Final Charge Distribution
f_final = fSol(end, :);

%% ✅ Plot Time Evolution and Final Distribution
figure;
hold on;
colorMap = jet(Zmax+1);
for i = 1:(Zmax+1)
   plot(tSol*1e3, fSol(:,i), 'LineWidth', 2, 'Color', colorMap(i,:));
end
xlabel('Time [ms]');
ylabel('Fractional Abundance');
title('Time-Dependent Ionization Balance of Oxygen (Ionization + Recombination Only)');
legend(arrayfun(@(z) sprintf('O^{%d+}', z-1), 1:Zmax+1, 'UniformOutput', false), 'Location', 'best');
grid on;
set(gca, 'XScale', 'log');

% ✅ Inset: Final vs. Initial Distribution
insetAxes = axes('Position', [0.55, 0.55, 0.35, 0.35]);
hold(insetAxes, 'on');
barHandle = bar(insetAxes, 0:Zmax, f_final, 'FaceColor', 'flat', 'FaceAlpha', 0.6);
plot(insetAxes, 0:Zmax, fSS, 'ko-', 'MarkerFaceColor', 'k', 'LineWidth', 2);
xticks(insetAxes, 0:Zmax);
xticklabels(insetAxes, arrayfun(@(z) sprintf('O^{%d+}', z), 0:Zmax, 'UniformOutput', false));
xlabel(insetAxes, 'Charge State');
ylabel(insetAxes, 'Fractional Abundance');
title(insetAxes, 'Final vs Initial');
legend(insetAxes, 'Final State', 'Initial State', 'Location', 'best');
grid(insetAxes, 'on');

for i = 1:(Zmax+1)
    barHandle.CData(i,:) = colorMap(i,:);
end

%% ✅ ODE System: Ionization + Recombination Only (No CX, No Velocity)
function dfdt = ionBalanceODE_no_CX(~, f, S, R, ne, Zmax)
    dfdt = zeros(Zmax+1, 1);
    dfdt(1) = - ne * S(1) * f(1) + ne * R(2) * f(2);

    for i = 2:Zmax
        dfdt(i) = ne * S(i) * f(i-1) - ne * (S(i+1) + R(i)) * f(i) + ne * R(i+1) * f(i+1);
    end
end