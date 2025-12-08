%% timeDependentIonization_ADAS_Ne_Velocity.m
% Computes the time-dependent ionization balance for Neon with velocity dependence.

clear; clc; close all;

%% 1) Load ADAS Data from NetCDF File
filename = 'ADAS_Rates_Ne.nc';

% --- Ionization Data ---
IonizationData.Temp = 10.^ncread(filename, 'gridTemperature_Ionization');
IonizationData.Density = 10.^ncread(filename, 'gridDensity_Ionization');
IonizationData.RateCoeff = 10.^ncread(filename, 'IonizationRateCoeff');

% --- Recombination Data ---
RecombinationData.Temp = 10.^ncread(filename, 'gridTemperature_Recombination');
RecombinationData.Density = 10.^ncread(filename, 'gridDensity_Recombination');
RecombinationData.RateCoeff = 10.^ncread(filename, 'RecombinationRateCoeff');

%% 2) Simulation and Plasma Parameters
Zmax = 10;                % Maximum charge state (Neon: from 0 to 10+)
Te = 10;              % Electron temperature [eV]
ne = 8.87072e17; 

m_p = 1.67e-27;   % Proton mass [kg]
m_Ne = 20 * m_p;  % Neon ion mass [kg]
kB = 1.38e-23;    % Boltzmann constant [J/K]
eV_to_J = 1.602e-19; % Conversion from eV to Joules

nNe = 0.02 * ne;          % ✅ Neon density (2% of n_e)

% Compute ion sound speed Cs
Cs = sqrt((5/3) * kB * (Te * eV_to_J) / m_Ne);  % [m/s]

% Use Cs as parallel velocity
v_parallel_max = Cs;     % ✅ Max plasma velocity (m/s)
tau_v = 1e-4;             % ✅ Velocity rise time (100 µs)
L = 0.27;                  % ✅ SOL characteristic length (m)
tspan = [logspace(-8, -4, 10000) linspace(1e-4, 1e1, 100000)]; % Improved resolution
nP = 1;                 % Scaling factor

%% 3) Compute S and R (α) Arrays via Interpolation
if isempty(gcp('nocreate'))
    parpool;
end

% ✅ Fix index bounds for Ionization Rates
S = zeros(Zmax+1, 1);  
parfor j = 0:(Zmax-1)  
    S(j+1) = nP * interpn(IonizationData.Density, IonizationData.Temp, IonizationData.RateCoeff(:,:,j+1), ne, Te, 'linear', 0);
end
S(Zmax+1) = 0;  

% ✅ Fix index bounds for Recombination Rates
R = zeros(Zmax+1, 1);
R(1) = 0;
parfor j = 1:Zmax
    R(j+1) = nP * interpn(RecombinationData.Density, RecombinationData.Temp, RecombinationData.RateCoeff(:,:,j), ne, Te, 'linear', 0);
end

%% 4) Set Initial Condition (Start Fully Neutral)
f0 = zeros(Zmax+1, 1); % ✅ Include neutrals in initial conditions
f0(1) = 1.0;           % ✅ Start with fully neutral Neon (100% Ne⁰)

%% 5) Solve the ODE System for Ionization Balance with Velocity (Using ode23s)
opts = odeset('RelTol',1e-6, 'AbsTol',1e-10); % ✅ Improved solver accuracy
[tSol, fSol] = ode23s(@(t, f) ionBalanceODE_Velocity(t, f, S, R, ne, Zmax, v_parallel_max, tau_v, L), tspan, f0, opts);

% Normalize solution at each time step
fSol = fSol ./ sum(fSol, 2);

%% 6) Compute Coronal Equilibrium Charge Distribution
fCE = zeros(Zmax+1, 1);
prodRatios = ones(Zmax, 1);
parfor z = 1:Zmax
    prodRatios(z) = prod(S(1:z) ./ R(2:z+1)); % ✅ Compute ratio of ionization/recombination rates
end
fCE(1) = 1 / (1 + sum(prodRatios)); % ✅ Normalize to ensure total sum is 1
for z = 1:Zmax
    fCE(z+1) = fCE(1) * prodRatios(z);
end

% Extract final charge state distribution
f_final = fSol(end, :);

%% 7) Plot Time Evolution with Inset for Final vs Coronal Equilibrium
figure;
hold on;

colorMap = jet(Zmax+1);
for i = 1:(Zmax+1)
   plot(tSol*1e3, fSol(:,i), 'LineWidth', 2, 'Color', colorMap(i,:));
end
xlabel('Time (ms)');
ylabel('Fractional Abundance');
title('Time-Dependent Ionization Balance of Neon with Velocity');
legend(arrayfun(@(z) sprintf('Ne^{%d+}', z-1), 1:Zmax+1, 'UniformOutput', false), 'Location', 'best');
grid on;
set(gca, 'XScale', 'log');
xlim([1e-2 1e4]);

% ✅ Add inset for Final vs. Coronal Equilibrium
insetAxes = axes('Position', [0.55, 0.55, 0.35, 0.35]); % Adjust position and size
hold(insetAxes, 'on');

% ✅ Bar plot for final charge distribution
barHandle = bar(insetAxes, 0:Zmax, f_final, 'FaceColor', 'flat', 'FaceAlpha', 0.6);

% ✅ Line plot for coronal equilibrium
plot(insetAxes, 0:Zmax, fCE, 'ko-', 'MarkerFaceColor', 'k', 'LineWidth', 2);

xticks(insetAxes, 0:Zmax);
xticklabels(insetAxes, arrayfun(@(z) sprintf('Ne^{%d+}', z), 0:Zmax, 'UniformOutput', false));
xlabel(insetAxes, 'Charge State');
ylabel(insetAxes, 'Fractional Abundance');
title(insetAxes, 'Final vs Coronal Equilibrium');
legend(insetAxes, 'Final Ionization State', 'Coronal Equilibrium', 'Location', 'best');
grid(insetAxes, 'on');

% Apply distinct colors to each bar
for i = 1:(Zmax+1)
    barHandle.CData(i,:) = colorMap(i,:);
end

%% ODE Function: Ionization Balance with Time-Dependent Velocity
function dfdt = ionBalanceODE_Velocity(t, f, S, R, ne, Zmax, v_parallel_max, tau_v, L)
    dfdt = zeros(Zmax+1, 1);
    
    % ✅ Time-dependent velocity to reduce stiffness
    v_parallel = v_parallel_max * (1 - exp(-t / tau_v));
    v_effect = v_parallel / L; % Plasma transport term

    % ✅ Neutral Neon (Ne⁰)
    dfdt(1) = - ne * S(1) * f(1) + ne * R(2) * f(2) - v_effect * (f(2) - f(1));

    % ✅ Intermediate charge states with improved transport effects
    for i = 2:Zmax
        dfdt(i) = ne * S(i) * f(i-1) - ne * (S(i+1) + R(i)) * f(i) + ne * R(i+1) * f(i+1);
        dfdt(i) = dfdt(i) - v_effect * (f(i+1) - f(i-1)) / 2;
    end
    
    dfdt(Zmax+1) = ne * S(Zmax) * f(Zmax) - ne * R(Zmax+1) * f(Zmax+1) - v_effect * (f(Zmax+1) - f(Zmax));

end