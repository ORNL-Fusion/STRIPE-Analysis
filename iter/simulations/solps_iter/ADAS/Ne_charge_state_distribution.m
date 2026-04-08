clear all
close all
clc

% Data source: 'adf11' (direct scd/acd files) or 'nc' (prebuilt NetCDF).
source_mode = 'adf11';

scd_file = 'scd89_ne.dat';
acd_file = 'acd89_ne.dat';
filename = 'ADAS_Rates_Ne.nc';

if strcmpi(source_mode, 'adf11')
    [Te_s, ne_s, scd_log10_cgs] = ADF11s(scd_file);
    [Te_a, ne_a, acd_log10_cgs] = ADF11a(acd_file);
    IonizationData.Temp = 10.^Te_s(:);
    IonizationData.Density = 10.^ne_s(:) * 1e6; % cm^-3 -> m^-3
    RecombinationData.Temp = 10.^Te_a(:);
    RecombinationData.Density = 10.^ne_a(:) * 1e6; % cm^-3 -> m^-3
    IonizationData.RateCoeff = 10.^scd_log10_cgs ./ 1e6;
    RecombinationData.RateCoeff = 10.^acd_log10_cgs ./ 1e6;
else
    % Read temperature/density grids (stored as log10 values in the NetCDF).
    IonizationData.Temp = 10.^ncread(filename, 'gridTemperature_Ionization');
    IonizationData.Density = 10.^ncread(filename, 'gridDensity_Ionization');
    RecombinationData.Temp = 10.^ncread(filename, 'gridTemperature_Recombination');
    RecombinationData.Density = 10.^ncread(filename, pick_var(filename, ...
        {'gridDensity_Recombine', 'gridDensity_Recombination'}));
    % Support both legacy and newer variable names.
    IonizationData.RateCoeff = 10.^ncread(filename, pick_var(filename, ...
        {'IonizationRateCoeff_log10_SI', 'IonizationRateCoeff'}));
    RecombinationData.RateCoeff = 10.^ncread(filename, pick_var(filename, ...
        {'RecombinationRateCoeff_log10_SI', 'RecombinationRateCoeff'}));
end

% Plasma conditions for equilibrium charge-state distribution.
Te_eV = 50;
ne_m3 = 1e18;
assert_in_range(Te_eV, 'Te_eV', IonizationData.Temp);
assert_in_range(ne_m3, 'ne_m3', IonizationData.Density);
assert_in_range(Te_eV, 'Te_eV', RecombinationData.Temp);
assert_in_range(ne_m3, 'ne_m3', RecombinationData.Density);

% Ne charge states: q = 0..10 (11 states), ionization transitions q=0..9.
Zmax = 10;
nStates = Zmax + 1;

% Eq. (3)-(4) coronal solution using rate coefficients at fixed (ne, Te).
% S_i: ionization i -> i+1 for i = 0..Zmax-1
% alpha_{i+1}: recombination i+1 -> i for i = 0..Zmax-1
S_coeff = zeros(Zmax, 1);
alpha_coeff = zeros(Zmax, 1);
for i = 0:(Zmax - 1)
    S_coeff(i + 1) = interpn(IonizationData.Density, IonizationData.Temp, ...
        IonizationData.RateCoeff(:, :, i + 1), ne_m3, Te_eV, 'linear', 0);
    alpha_coeff(i + 1) = interpn(RecombinationData.Density, RecombinationData.Temp, ...
        RecombinationData.RateCoeff(:, :, i + 1), ne_m3, Te_eV, 'linear', 0);
end

% Eq. (3): n_z = (prod_{i=0}^{z-1} S_i/alpha_{i+1}) * n_0
prod_terms = zeros(Zmax, 1);
running_prod = 1;
for z = 1:Zmax
    running_prod = running_prod * (S_coeff(z) / max(alpha_coeff(z), realmin));
    prod_terms(z) = running_prod;
end

% Eq. (4): n_0 = N / (1 + sum_{z=1}^{Z} prod_{i=0}^{z-1} S_i/alpha_{i+1})
N_total = 1.0;
n0 = N_total / (1 + sum(prod_terms));
conc = zeros(nStates, 1);
conc(1) = n0;
for z = 1:Zmax
    conc(z + 1) = prod_terms(z) * n0;
end
conc = conc / sum(conc);

figure('Color', 'w');
plot(0:Zmax, conc, '--o', 'LineWidth', 1.8, 'MarkerSize', 4);
title({ ...
    sprintf('Equilibrium Charge State Distribution of Ne at T_e = %.2f eV', Te_eV), ...
    sprintf('n_e = %.2e m^{-3}', ne_m3) ...
    });
xlabel('Charge State q');
ylabel('Distribution Fraction');
set(gca, 'FontSize', 14);
grid on;
box on;
xlim([0 Zmax]);

% Single-energy ionization mean free path example.
mass_amu = 20.1797;
Ti_imp_eV = 36;
Te_mfp_eV = 4;
ne_mfp_m3 = 1e17;
impurity_charge = 0; % Ne^0 -> Ne^1+
assert_in_range(Te_mfp_eV, 'Te_mfp_eV', IonizationData.Temp);
assert_in_range(ne_mfp_m3, 'ne_mfp_m3', IonizationData.Density);

v_imp = sqrt(2 * Ti_imp_eV * 1.602176634e-19 / (mass_amu * 1.66053906660e-27));
k_ion_mfp = interpn(IonizationData.Density, IonizationData.Temp, ...
    IonizationData.RateCoeff(:, :, impurity_charge + 1), ne_mfp_m3, Te_mfp_eV, 'linear', 0);
tion = 1 / (ne_mfp_m3 * k_ion_mfp);
mfp = v_imp * tion;

fprintf('Ne^%d ionization at Te=%.2f eV, ne=%.2e m^-3\n', impurity_charge, Te_mfp_eV, ne_mfp_m3);
fprintf('k_ion = %.6e m^3/s, tion = %.6e s, mfp = %.6e m\n', k_ion_mfp, tion, mfp);

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
