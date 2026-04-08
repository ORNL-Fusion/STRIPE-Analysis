%% Cleaned / corrected oxygen charge-state equilibrium script
clear all
close all
clc

filename = 'ADAS_Rates_O.nc';

% --- Read ADAS grids & rates (log10 stored -> convert to linear) ---
IonizationData.Temp    = 10.^ncread(filename,'gridTemperature_Ionization');   % [eV]
IonizationData.Density = 10.^ncread(filename,'gridDensity_Ionization');       % [m^-3]
IonizationData.RateCoeff= 10.^ncread(filename, pick_var(filename, ...
    {'IonizationRateCoeff_log10_SI','IonizationRateCoeff'}));                 % [m^3/s]

RecombinationData.Temp    = 10.^ncread(filename,'gridTemperature_Recombination');
RecombinationData.Density = 10.^ncread(filename, pick_var(filename, ...
    {'gridDensity_Recombine','gridDensity_Recombination'}));
RecombinationData.RateCoeff = 10.^ncread(filename, pick_var(filename, ...
    {'RecombinationRateCoeff_log10_SI','RecombinationRateCoeff'}));

% --- user / physics inputs ---
% Temperatures to evaluate (eV). You can supply a vector, e.g. [10 15 20].
T_eV = 15;          % eV (scalar or vector)
n_e  = 6e18;        % electron density [m^-3]
Zmax = 8;           % maximum nuclear charge for O -> 0..8

% Ensure T_eV is a column vector for indexing
T_eV = T_eV(:);
nT = numel(T_eV);
assert_in_range(n_e, 'n_e', IonizationData.Density);
assert_in_range(n_e, 'n_e', RecombinationData.Density);
for it = 1:nT
    assert_in_range(T_eV(it), sprintf('T_eV(%d)', it), IonizationData.Temp);
    assert_in_range(T_eV(it), sprintf('T_eV(%d)', it), RecombinationData.Temp);
end

% Pre-allocate S (ionization) and a (recombination)
% S: nT x (Zmax+1) rates (for ionization from charge Z -> Z+1)
% a: nT x (Zmax+1) rates (for recombination from Z -> Z-1)
S = zeros(nT, Zmax+1);
a = zeros(nT, Zmax+1);

% Interpolate ADAS tables for each Z
% Note: ADAS indexing: RateCoeff(:,:,k) should correspond to process for charge index k (check your file)
for ii = 1:(Zmax+1)
    Z = ii-1;
    % Ionization rate S for Z -> Z+1 is stored in IonizationData.RateCoeff(:,:,Z+1)
    if Z <= Zmax-1
        S(:,ii) = interpn(IonizationData.Density, IonizationData.Temp, ...
                          squeeze(IonizationData.RateCoeff(:,:,ii)), n_e, T_eV, 'linear', 0);
    end
    % Recombination rate a for Z -> Z-1 is stored in RecombinationData.RateCoeff(:,:,Z)
    if Z > 0
        a(:,ii) = interpn(RecombinationData.Density, RecombinationData.Temp, ...
                          squeeze(RecombinationData.RateCoeff(:,:,ii)), n_e, T_eV, 'linear', 0);
    end
end

% --- Eq. (3)-(4) coronal equilibrium populations for each temperature ---
conc = zeros(Zmax+1, nT);
N_total = 1.0;
for it = 1:nT
    % Eq. (3): n_z = (prod_{i=0}^{z-1} S_i/alpha_{i+1}) * n_0
    prod_terms = zeros(Zmax,1);
    running_prod = 1;
    for z = 1:Zmax
        running_prod = running_prod * (S(it,z) / max(a(it,z+1), realmin));
        prod_terms(z) = running_prod;
    end

    % Eq. (4): n_0 = N / (1 + sum_{z=1}^Z prod_{i=0}^{z-1} S_i/alpha_{i+1})
    n0 = N_total / (1 + sum(prod_terms));
    p = zeros(Zmax+1,1);
    p(1) = n0;
    for z = 1:Zmax
        p(z+1) = prod_terms(z) * n0;
    end

    % Normalize fractional distribution
    p = p / sum(p);
    conc(:,it) = p;
end

% --- Plot results for last temperature (or all if vector) ---
figure(1)
hold on
colors = lines(nT);
for it = 1:nT
    plot(0:Zmax, conc(:,it), '-o', 'LineWidth', 2, 'Color', colors(it,:));
end
title({'Equilibrium Charge State Distribution','Oxygen (ADAS)'} ,'Interpreter','none')
xlabel('Charge State [Z]')
ylabel('Fraction')
set(gca,'FontSize',14)
axis([0 Zmax 0 1])
legendStrings = arrayfun(@(x) sprintf('T=%.1f eV', T_eV(x)), 1:nT, 'UniformOutput', false);
legend(legendStrings,'Location','best')

% --- Example: compute neutral ionization time and mean free path for chosen T ---
% Use neutral ionization rate (Z=0 -> 1) stored in RateCoeff(:,:,1)
% pick first temperature index for demonstration
it0 = 1;
S0 = S(it0,1);   % [m^3/s]
if S0 > 0
    tion = 1 / ( n_e * S0 );       % ionization time [s]
else
    tion = Inf;
end

% Thermal velocity of oxygen atoms (use atomic mass A=16)
A_amu = 16;                         % oxygen atomic mass
m_kg = A_amu * 1.66053906660e-27;   % kg
ti0_eV = 4;                         % example neutral temperature (eV) -> set as needed
vTh = sqrt(2 * ti0_eV * 1.602176634e-19 / m_kg);  % m/s

mfp = vTh * tion;   % mean free path [m]
fprintf('For T=%g eV: S0=%.3e m^3/s, tion=%.3e s, vTh=%.3e m/s, mfp=%.3e m\n', ...
        T_eV(it0), S0, tion, vTh, mfp);

function v = pick_var(ncfile, candidates)
    info = ncinfo(ncfile);
    names = {info.Variables.Name};
    for i = 1:numel(candidates)
        if any(strcmp(names, candidates{i}))
            v = candidates{i};
            return;
        end
    end
    error('None of the candidate variables were found in %s.', ncfile);
end

function assert_in_range(x, name, grid)
    gmin = min(grid);
    gmax = max(grid);
    assert(x >= gmin && x <= gmax, ...
        '%s=%.6g is outside data range [%.6g, %.6g].', name, x, gmin, gmax);
end
