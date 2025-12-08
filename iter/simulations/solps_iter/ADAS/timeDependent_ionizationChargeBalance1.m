clear; clc; close all;

% ✅ Load ADAS Data from NetCDF File
filename = 'ADAS_Rates_O.nc';

% ✅ Ionization Data
IonizationData.Temp = 10.^ncread(filename, 'gridTemperature_Ionization');
IonizationData.Density = 10.^ncread(filename, 'gridDensity_Ionization');
IonizationData.RateCoeff = 10.^ncread(filename, 'IonizationRateCoeff');

% ✅ Recombination Data
RecombinationData.Temp = 10.^ncread(filename, 'gridTemperature_Recombination');
RecombinationData.Density = 10.^ncread(filename, 'gridDensity_Recombine'); 
RecombinationData.RateCoeff = 10.^ncread(filename, 'RecombinationRateCoeff');

% ✅ Charge Exchange Data
CXData.Temp = 10.^ncread(filename, 'gridTemperature_CX');
CXData.Density = 10.^ncread(filename, 'gridDensity_CX'); 
CXData.RateCoeff = 10.^ncread(filename, 'ChargeExchangeRateCoeff'); % Cross-section * velocity

%% **2) Simulation and Plasma Parameters**
Zmax = 8;        % Oxygen charge states: O^0 to O^8+
Te = 6;          % Electron temperature [eV]
ne = 8.87072e17; 

m_p = 1.67e-27;   % Proton mass [kg]
m_O = 16 * m_p;   % Oxygen ion mass [kg]
kB = 1.38e-23;    % Boltzmann constant [J/K]
eV_to_J = 1.602e-19; % eV to Joules

nO = 0.02 * ne;   % Oxygen density (2% of n_e)
n_n = 0.01 * ne;  % Neutral density (1% of n_e) for CX

% Ion sound speed Cs
Cs = sqrt((5/3) * kB * (Te * eV_to_J) / m_O);  % [m/s]

% Velocity parameters
v_parallel_max = Cs;
tau_v = 1e-4;     % Velocity rise time
L = 0.27;         % SOL length [m]
tspan = [logspace(-8, -4, 10000) linspace(1e-4, 1e1, 100000)];

%% **3) Interpolate Rate Coefficients**
if isempty(gcp('nocreate'))
    parpool;
end

S = zeros(Zmax+1, 1);  
R = zeros(Zmax+1, 1);
CX = zeros(Zmax+1, 1);
R(1) = 0;  % No recombination from neutral

parfor j = 1:Zmax
    S(j) = interpn(IonizationData.Density, IonizationData.Temp, IonizationData.RateCoeff(:,:,j), ne, Te, 'linear', 0);
    R(j) = interpn(RecombinationData.Density, RecombinationData.Temp, RecombinationData.RateCoeff(:,:,j), ne, Te, 'linear', 0);
    
    S(j) = max(S(j), eps);
    R(j) = max(R(j), eps);
    
    if size(CXData.RateCoeff, 3) >= j
        CX(j) = interpn(CXData.Density, CXData.Temp, CXData.RateCoeff(:,:,j), ne, Te, 'linear', 0);
    else
        CX(j) = 0;
    end
end

% Last index for CX
if size(CXData.RateCoeff, 3) >= Zmax+1
    CX(Zmax+1) = interpn(CXData.Density, CXData.Temp, CXData.RateCoeff(:,:,Zmax+1), ne, Te, 'linear', 0);
else
    CX(Zmax+1) = 0;
end

%% **4) Initialize with Fully Stripped Oxygen**
fSS = zeros(Zmax+1, 1);
fSS(end) = 1;  % O^8+ only
disp('✅ Initialized with Fully Stripped Oxygen (O^{8+})');
disp(fSS);

%% **5) Solve ODE**
opts = odeset('RelTol',1e-6, 'AbsTol',1e-10);
[tSol, fSol] = ode23s(@(t, f) ionBalanceODE_Velocity_CX(t, f, S, R, CX, ne, n_n, Zmax, v_parallel_max, tau_v, L), tspan, fSS, opts);
fSol = fSol ./ sum(fSol, 2);  % Normalize at each time step

%% **6) Final Distribution**
f_final = fSol(end, :);

%% **7) Plot Time Evolution**
figure;
hold on;
colorMap = jet(Zmax+1);
for i = 1:(Zmax+1)
   plot(tSol*1e3, fSol(:,i), 'LineWidth', 2, 'Color', colorMap(i,:));
end
xlabel('Time (ms)');
ylabel('Fractional Abundance');
title('Time-Dependent Ionization Balance of Oxygen');
legend(arrayfun(@(z) sprintf('O^{%d+}', z-1), 0:Zmax, 'UniformOutput', false), 'Location', 'best');
grid on;
set(gca, 'XScale', 'log');

% ✅ Inset for Final vs Initial
insetAxes = axes('Position', [0.55, 0.55, 0.35, 0.35]); 
hold(insetAxes, 'on');

barHandle = bar(insetAxes, 0:Zmax, f_final, 'FaceColor', 'flat', 'FaceAlpha', 0.6);
plot(insetAxes, 0:Zmax, fSS, 'ko-', 'MarkerFaceColor', 'k', 'LineWidth', 2);

xticks(insetAxes, 0:Zmax);
xticklabels(insetAxes, arrayfun(@(z) sprintf('O^{%d+}', z), 0:Zmax, 'UniformOutput', false));
xlabel(insetAxes, 'Charge State');
ylabel(insetAxes, 'Fractional Abundance');
title(insetAxes, 'Final vs Initial State');
legend(insetAxes, 'Final Ionization State', 'Initial State (O^{8+})', 'Location', 'best');
grid(insetAxes, 'on');

for i = 1:(Zmax+1)
    barHandle.CData(i,:) = colorMap(i,:);
end

%% **8) ODE Function**
function dfdt = ionBalanceODE_Velocity_CX(t, f, S, R, CX, ne, n_n, Zmax, v_parallel_max, tau_v, L)
    dfdt = zeros(Zmax+1, 1);
    v_parallel = v_parallel_max * (1 - exp(-t / tau_v));
    v_effect = v_parallel / L;

    % O^0
    dfdt(1) = -ne * S(1) * f(1) + ne * R(2) * f(2) ...
              - v_effect * (f(2) - f(1)) ...
              - n_n * CX(1) * f(1) + n_n * CX(2) * f(2);

    % O^1+ to O^7+
    for i = 2:Zmax
        dfdt(i) = ne * S(i) * f(i-1) ...
                - ne * (S(i+1) + R(i)) * f(i) ...
                + ne * R(i+1) * f(i+1);
        dfdt(i) = dfdt(i) ...
                - v_effect * (f(i+1) - f(i-1)) / 2 ...
                - n_n * CX(i) * f(i) + n_n * CX(i+1) * f(i+1);
    end

    % O^8+
    dfdt(Zmax+1) = ne * S(Zmax) * f(Zmax) - ne * R(Zmax+1) * f(Zmax+1);
    dfdt(Zmax+1) = dfdt(Zmax+1) ...
                 - v_effect * (f(Zmax+1) - f(Zmax-1)) / 2 ...
                 - n_n * CX(Zmax+1) * f(Zmax+1);
end