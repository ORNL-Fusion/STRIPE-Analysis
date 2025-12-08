clear; clc; close all;

% ✅ Load ADAS Data from NetCDF File
filename = 'ADAS_Rates_Ne.nc';
% filename = 'ADAS_Rates_O.nc';


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
Zmax = 10;        % Maximum charge state (Neon: from 0 to 10+)
Te = 7;          % Electron temperature [eV]
ne = 5e16; 

m_p = 1.67e-27;   % Proton mass [kg]
m_Ne = 16 * m_p;  % Neon ion mass [kg]
kB = 1.38e-23;    % Boltzmann constant [J/K]
eV_to_J = 1.602e-19; % Conversion from eV to Joules

nNe = 0.02 * ne;   % ✅ Neon density (2% of n_e)
n_n = 0.02 * ne;   % ✅ Neutral density (1% of n_e) for Charge Exchange

% Compute ion sound speed Cs
Cs = sqrt((5/3) * kB * (Te * eV_to_J) / m_Ne);  % [m/s]

% Use Cs as parallel velocity
v_parallel_max = Cs;     % ✅ Max plasma velocity (m/s)
tau_v = 1e-4;            % ✅ Velocity rise time (100 µs)
L = 0.27;                % ✅ SOL characteristic length (m)
tspan = [logspace(-8, -4, 10000) linspace(1e-4, 1e1, 100000)]; % Improved resolution
nP = 1;                 % Scaling factor

%% **3) Compute S, R, and CX Arrays via Interpolation**
if isempty(gcp('nocreate'))
    parpool;
end

% ✅ Fix index bounds for Ionization, Recombination, and Charge Exchange Rates
S = zeros(Zmax+1, 1);  
R = zeros(Zmax+1, 1);
CX = zeros(Zmax+1, 1);
R(1) = 0;  % No recombination for neutrals

parfor j = 1:Zmax
    S(j) = interpn(IonizationData.Density, IonizationData.Temp, IonizationData.RateCoeff(:,:,j), ne, Te, 'linear', 0);
    R(j) = interpn(RecombinationData.Density, RecombinationData.Temp, RecombinationData.RateCoeff(:,:,j), ne, Te, 'linear', 0);

    % Prevent zero values
    S(j) = max(S(j), eps);
    R(j) = max(R(j), eps);
    1
    % Charge Exchange only if valid dataset index
    if size(CXData.RateCoeff, 3) >= j
        CX(j) = interpn(CXData.Density, CXData.Temp, CXData.RateCoeff(:,:,j), ne, Te, 'linear', 0);
    else
        CX(j) = 0;
    end
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

%% **5) Solve the ODE System for Ionization Balance with Charge Exchange**
opts = odeset('RelTol',1e-6, 'AbsTol',1e-10); % ✅ Improved solver accuracy
[tSol, fSol] = ode23s(@(t, f) ionBalanceODE_Velocity_CX(t, f, S, R, CX, ne, n_n, Zmax, v_parallel_max, tau_v, L), tspan, fSS, opts);

% Normalize solution at each time step
fSol = fSol ./ sum(fSol, 2);

%% **6) Extract Final Charge State Distribution**
f_final = fSol(end, :);

%% **7) Plot Time Evolution with Inset for Final vs Equilibrium Initial Distribution**
figure;
hold on;
colorMap = jet(Zmax+1);
for i = 1:(Zmax+1)
   plot(tSol*1e3, fSol(:,i), 'LineWidth', 2, 'Color', colorMap(i,:));
end
xlabel('Time (ms)');
ylabel('Fractional Abundance');
title('Time-Dependent Ionization Balance of Neon with Charge Exchange & Velocity');
legend(arrayfun(@(z) sprintf('Ne^{%d+}', z-1), 1:Zmax+1, 'UniformOutput', false), 'Location', 'best');
grid on;
set(gca, 'XScale', 'log');

% ✅ Add inset for Final vs. Equilibrium Initial Charge Distribution
insetAxes = axes('Position', [0.55, 0.55, 0.35, 0.35]); % Adjust position and size
hold(insetAxes, 'on');

% ✅ Bar plot for final charge distribution
barHandle = bar(insetAxes, 0:Zmax, f_final, 'FaceColor', 'flat', 'FaceAlpha', 0.6);

% ✅ Line plot for equilibrium initial charge distribution
plot(insetAxes, 0:Zmax, fSS, 'ko-', 'MarkerFaceColor', 'k', 'LineWidth', 2);

xticks(insetAxes, 0:Zmax);
xticklabels(insetAxes, arrayfun(@(z) sprintf('Ne^{%d+}', z), 0:Zmax, 'UniformOutput', false));
xlabel(insetAxes, 'Charge State');
ylabel(insetAxes, 'Fractional Abundance');
title(insetAxes, 'Final vs Equilibrium Initial');
legend(insetAxes, 'Final Ionization State', 'Equilibrium Initial State', 'Location', 'best');
grid(insetAxes, 'on');

% Apply distinct colors to each bar
for i = 1:(Zmax+1)
    barHandle.CData(i,:) = colorMap(i,:);
end
%% **ODE Function: Ionization Balance with Charge Exchange**
function dfdt = ionBalanceODE_Velocity_CX(t, f, S, R, CX, ne, n_n, Zmax, v_parallel_max, tau_v, L)
    dfdt = zeros(Zmax+1, 1);
    
    v_parallel = v_parallel_max * (1 - exp(-t / tau_v));
    v_effect = v_parallel / L; 

    dfdt(1) = - ne * S(1) * f(1) + ne * R(2) * f(2) - v_effect * (f(2) - f(1)) - n_n * CX(1) * f(1) + n_n * CX(2) * f(2);

    for i = 2:Zmax
        dfdt(i) = ne * S(i) * f(i-1) - ne * (S(i+1) + R(i)) * f(i) + ne * R(i+1) * f(i+1);
        dfdt(i) = dfdt(i) - v_effect * (f(i+1) - f(i-1)) / 2 - n_n * CX(i) * f(i) + n_n * CX(i+1) * f(i+1);
    end
end