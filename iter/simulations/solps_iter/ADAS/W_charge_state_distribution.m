clear all
close all
clc

% Data source: 'adf11' (direct scd/acd files) or 'nc' (prebuilt NetCDF).
source_mode = 'adf11';

scd_file = 'scd89_w.dat';
acd_file = 'acd89_w.dat';
filename = 'ADAS_Rates_W.nc';

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
    RecombinationData.Density = 10.^ncread(filename, 'gridDensity_Recombine');
    % Support both legacy and new variable names.
    IonizationData.RateCoeff = 10.^ncread(filename, pick_var(filename, ...
        {'IonizationRateCoeff_log10_SI', 'IonizationRateCoeff'}));
    RecombinationData.RateCoeff = 10.^ncread(filename, pick_var(filename, ...
        {'RecombinationRateCoeff_log10_SI', 'RecombinationRateCoeff'}));
end

% Plasma conditions for equilibrium charge-state distribution.
Te_eV = 15;
ne_m3 = 6e18;
assert_in_range(Te_eV, 'Te_eV', IonizationData.Temp);
assert_in_range(ne_m3, 'ne_m3', IonizationData.Density);
assert_in_range(Te_eV, 'Te_eV', RecombinationData.Temp);
assert_in_range(ne_m3, 'ne_m3', RecombinationData.Density);

% W charge states: q = 0..74 (75 states), ionization transitions q=0..73.
Zmax = 74;
nStates = Zmax + 1;

% Transition rates [1/s] at fixed (ne, Te).
S_rate = zeros(nStates, 1);  % ionization q -> q+1
a_rate = zeros(nStates, 1);  % recombination q -> q-1

for q = 0:(Zmax - 1)
    k_ion = interpn(IonizationData.Density, IonizationData.Temp, ...
        IonizationData.RateCoeff(:, :, q + 1), ne_m3, Te_eV, 'linear', 0);
    S_rate(q + 1) = ne_m3 * k_ion;
end

for q = 1:Zmax
    k_rec = interpn(RecombinationData.Density, RecombinationData.Temp, ...
        RecombinationData.RateCoeff(:, :, q), ne_m3, Te_eV, 'linear', 0);
    a_rate(q + 1) = ne_m3 * k_rec;
end

% Solve steady-state balance A*c = b with normalization sum(c)=1.
A = zeros(nStates, nStates);
b = zeros(nStates, 1);
A(1, :) = 1;
b(1) = 1;

for q = 1:Zmax
    idx = q + 1; % row/column index for charge state q

    % Inflow from q-1 by ionization.
    A(idx, idx - 1) = S_rate(idx - 1);

    % Loss from q by ionization and recombination.
    A(idx, idx) = -(S_rate(idx) + a_rate(idx));

    % Inflow from q+1 by recombination (not present at fully stripped state).
    if idx < nStates
        A(idx, idx + 1) = a_rate(idx + 1);
    end
end

conc = A \ b;
conc = max(conc, 0);
conc = conc / sum(conc);

figure('Color', 'w');
plot(0:Zmax, conc, '--o', 'LineWidth', 1.8, 'MarkerSize', 4);
title({ ...
    sprintf('Equilibrium Charge State Distribution of W at T_e = %.2f eV', Te_eV), ...
    sprintf('n_e = %.2e m^{-3}', ne_m3) ...
    });
xlabel('Charge State q');
ylabel('Distribution Fraction');
set(gca, 'FontSize', 14);
grid on;
box on;
xlim([0 Zmax]);

% Single-energy ionization mean free path example.
mass_amu = 183.84;
Ti_imp_eV = 4;
Te_mfp_eV = 20;
ne_mfp_m3 = 1e19;
impurity_charge = 0; % W^0 -> W^1+
assert_in_range(Te_mfp_eV, 'Te_mfp_eV', IonizationData.Temp);
assert_in_range(ne_mfp_m3, 'ne_mfp_m3', IonizationData.Density);

v_imp = sqrt(2 * Ti_imp_eV * 1.602176634e-19 / (mass_amu * 1.66053906660e-27));
k_ion_mfp = interpn(IonizationData.Density, IonizationData.Temp, ...
    IonizationData.RateCoeff(:, :, impurity_charge + 1), ne_mfp_m3, Te_mfp_eV, 'linear', 0);
tion = 1 / (ne_mfp_m3 * k_ion_mfp);
mfp = v_imp * tion;

fprintf('W^%d ionization at Te=%.2f eV, ne=%.2e m^-3\n', impurity_charge, Te_mfp_eV, ne_mfp_m3);
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
