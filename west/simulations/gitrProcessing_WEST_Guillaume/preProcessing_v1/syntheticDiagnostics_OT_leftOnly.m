close all;
clear;
clc;

%% === Load Species Data for Thermal and Ohmic ===
species = {'o8', 'o7', 'o6', 'o5'};
species_list = {'o8', 'o7', 'o6', 'o5'};
colors = {[0 0 1], [1 0.5 0], [1 0 1], [0 0 0]};
markers = {'s', 'v', 'o', '+'};

for i = 1:numel(species)
    dataTh{i} = readmatrix(['Targets_thermal_', species{i}, '.txt']);
    dataOh{i} = readmatrix(['Targets_thermal_', species{i},'_ohmic', '.txt']);
end

%% === Define Limiter Indices ===
z_vertL = [0.4; 0.3473; 0.1689; 0.0749; 0.0106; -0.0851; -0.3042; -0.4];
thermalL_indices = {
    [3863, 3858, 3859, 3860, 3870, 3476, 3476, 3449, 3450]; 
    [595, 633, 633, 634, 636, 460, 460, 447, 446]; 
    [550, 550, 562, 564, 492, 488, 485, 426, 2660];
    [508, 516, 3918, 3920, 3919, 480, 471, 3563, 2715]; 
    [4004, 4003, 12589, 12591, 12592, 3600, 3576, 3575, 2746]; 
    [926, 940, 12621, 12621, 12623, 3631, 3632, 3621, 2779];
    [4091, 4110, 4112, 12664, 12662, 12662, 3727, 3725, 3702]; 
    [4141, 4145, 4145, 4146, 4146, 4147, 3752, 3751, 3749]
};

%% === Preallocate Arrays ===
nRegions = length(z_vertL);
nSpecies = length(species);

niL_Th = zeros(nRegions, nSpecies); TiL_Th = niL_Th;
VL_Th  = niL_Th; fluxL_Th = niL_Th;

niL_Oh = niL_Th; TiL_Oh = niL_Th;
VL_Oh  = niL_Th; fluxL_Oh = niL_Th;

neL_Th = niL_Th; TeL_Th = niL_Th;
neL_Oh = niL_Th; TeL_Oh = niL_Th;

%% === Compute Region-Averaged Quantities (Left Only) ===
for s = 1:nSpecies
    dTh = dataTh{s}; dOh = dataOh{s};
    for r = 1:nRegions
        niL_Th(r,s) = mean(dTh(thermalL_indices{r}, 11));
        TiL_Th(r,s) = mean(dTh(thermalL_indices{r}, 4));
        VL_Th(r,s)  = max(dTh(thermalL_indices{r}, 1));
        fluxL_Th(r,s) = mean(dTh(thermalL_indices{r},11) .* dTh(thermalL_indices{r},5));
        neL_Th(r,s) = mean(dTh(thermalL_indices{r}, 2));
        TeL_Th(r,s) = mean(dTh(thermalL_indices{r}, 3));

        niL_Oh(r,s) = mean(dOh(thermalL_indices{r}, 11));
        TiL_Oh(r,s) = mean(dOh(thermalL_indices{r}, 4));
        VL_Oh(r,s)  = max(dOh(thermalL_indices{r}, 1));
        fluxL_Oh(r,s) = mean(dOh(thermalL_indices{r},11) .* dOh(thermalL_indices{r},5));
        neL_Oh(r,s) = mean(dOh(thermalL_indices{r}, 2));
        TeL_Oh(r,s) = mean(dOh(thermalL_indices{r}, 3));
    end
end

%% === Plot Thermal Case (Left Limiter Only) ===
figure('Name', 'Thermal Case: Left Limiter Profiles', 'Position', [100, 100, 1000, 800]);
titles = {'(a)', '(b)', '(c)', '(d)'};
subplot_titles = {'n_O [m^{-3}]', 'T_O [eV]', '\Gamma_O [m^{-2}·s^{-1}]', 'V_{Sheath}^{Thermal} [V]'};
varsL_Th = {niL_Th, TiL_Th, fluxL_Th, VL_Th};

for v = 1:4
    subplot(4,1,v); hold on;
    for i = 1:nSpecies
        plot(z_vertL, varsL_Th{v}(:,i), [markers{i} '-'], 'Color', colors{i}, ...
            'DisplayName', [subplot_titles{v}(1:2), ' O^{', species_list{i}(2), '+}']);
    end
    ylabel(subplot_titles{v});
    if v == 4, xlabel('z [m]'); end
    title(titles{v});
    legend('Location', 'best', 'NumColumns', 2); grid on;
end

%% === Plot Ohmic Case (Left Limiter Only) ===
figure('Name', 'Ohmic Case: Left Limiter Profiles', 'Position', [100, 100, 1000, 800]);
varsL_Oh = {niL_Oh, TiL_Oh, fluxL_Oh, VL_Oh};

for v = 1:4
    subplot(4,1,v); hold on;
    for i = 1:nSpecies
        plot(z_vertL, varsL_Oh{v}(:,i), [markers{i} '-'], 'Color', colors{i}, ...
            'DisplayName', [subplot_titles{v}(1:2), ' O^{', species_list{i}(2), '+}']);
    end
    ylabel(subplot_titles{v});
    if v == 4, xlabel('z [m]'); end
    title(titles{v});
    legend('Location', 'best', 'NumColumns', 2); grid on;
end

%% === Plot Electron Profiles (Thermal) ===
figure('Name', 'Thermal Case: Electron Profiles (Left Limiter)', 'Position', [100, 100, 1000, 600]);
subplot(2,1,1); hold on;
for i = 1:nSpecies
    plot(z_vertL, neL_Th(:,i), [markers{i} '-'], 'Color', colors{i}, ...
        'DisplayName', ['n_e O^{', species_list{i}(2), '+}']);
end
ylabel('n_e [m^{-3}]'); xlabel('z [m]');
title('(e) Thermal Case: Electron Density');
legend('Location', 'best'); grid on;

subplot(2,1,2); hold on;
for i = 1:nSpecies
    plot(z_vertL, TeL_Th(:,i), [markers{i} '-'], 'Color', colors{i}, ...
        'DisplayName', ['T_e O^{', species_list{i}(2), '+}']);
end
ylabel('T_e [eV]'); xlabel('z [m]');
title('(f) Thermal Case: Electron Temperature');
legend('Location', 'best'); grid on;

%% === Plot Electron Profiles (Ohmic) ===
figure('Name', 'Ohmic Case: Electron Profiles (Left Limiter)', 'Position', [100, 100, 1000, 600]);
subplot(2,1,1); hold on;
for i = 1:nSpecies
    plot(z_vertL, neL_Oh(:,i), [markers{i} '-'], 'Color', colors{i}, ...
        'DisplayName', ['n_e O^{', species_list{i}(2), '+}']);
end
ylabel('n_e [m^{-3}]'); xlabel('z [m]');
title('(e) Ohmic Case: Electron Density');
legend('Location', 'best'); grid on;

subplot(2,1,2); hold on;
for i = 1:nSpecies
    plot(z_vertL, TeL_Oh(:,i), [markers{i} '-'], 'Color', colors{i}, ...
        'DisplayName', ['T_e O^{', species_list{i}(2), '+}']);
end
ylabel('T_e [eV]'); xlabel('z [m]');
title('(f) Ohmic Case: Electron Temperature');
legend('Location', 'best'); grid on;

%% === Erosion Flux Calculations ===
yield_files_th = {'yields_thermal_o5.csv', 'yields_thermal_o6.csv', ...
                  'yields_thermal_o7.csv', 'yields_thermal_o8.csv'};
yield_files_oh = {'yields_thermal_o5_ohmic.csv', 'yields_thermal_o6.csv', ...
                  'yields_thermal_o7_ohmic.csv', 'yields_thermal_o8_ohmic.csv'};

yield_data_th = cell(1,4); yield_data_oh = cell(1,4);
for i = 1:4
    yield_data_th{i} = readmatrix(yield_files_th{i});
    yield_data_oh{i} = readmatrix(yield_files_oh{i});
end

erosionL_Th = zeros(nRegions, 4); erosionL_Oh = zeros(nRegions, 4);

for s = 1:4
    dTh = dataTh{s}; dOh = dataOh{s};
    if s == 1
        Y_th = [0; yield_data_th{4}(:,8)];
        Y_oh = [0; yield_data_oh{4}(:,8)];
    elseif s == 2
        Y_th = [0; yield_data_th{3}(:,7)];
        Y_oh = [0; yield_data_oh{3}(:,7)];
    elseif s == 3
        Y_th = [0; yield_data_th{2}(:,6)];
        Y_oh = [0; yield_data_oh{2}(:,6)];
    else
        Y_th = [0; yield_data_th{1}(:,5)];
        Y_oh = [0; yield_data_oh{1}(:,5)];
    end

    erosion_full_T = dTh(:,11) .* abs(dTh(:,5)) .* Y_th;
    erosion_full_O = dOh(:,11) .* abs(dOh(:,5)) .* Y_oh;

    for r = 1:nRegions
        erosionL_Th(r,s) = mean(erosion_full_T(thermalL_indices{r}));
        erosionL_Oh(r,s) = mean(erosion_full_O(thermalL_indices{r}));
    end
end

sxbValueL = [9.58779; 10.8139; 12.1153; 12.8246; 13.4241; 13.242; 11.1808; 9.70449];
flux_thL = cell(1,4); flux_ohL = cell(1,4);
for s = 1:4
    flux_thL{s} = erosionL_Th(:,s) ./ (4 * pi * sxbValueL);
    flux_ohL{s} = erosionL_Oh(:,s) ./ (4 * pi * sxbValueL);
end
flux_thL_total = flux_thL{1} + flux_thL{2} + flux_thL{3} + flux_thL{4};
flux_ohL_total = flux_ohL{1} + flux_ohL{2} + flux_ohL{3} + flux_ohL{4};

%% === Plot Erosion Flux (Thermal) ===
figure('Name', 'Thermal Erosion Flux (Left Limiter)');
for i = 1:4
    plot(z_vertL, 10*pi.*flux_thL{i}, 'LineStyle', '-', 'Marker', markers{i}, ...
        'Color', colors{i}, 'DisplayName', ['O^{', species_list{i}(2), '+} (STRIPE)']); hold on;
end
plot(z_vertL, 10*pi.*flux_thL_total, 'g*-', 'DisplayName', 'I_{\phi, W-I} (STRIPE)');
xlabel('z [m]');
ylabel('Brightness I_\phi(W-I) [m^{-2}s^{-1}St^{-1}]');
title('(b) ICRH Case w/ Thermal Sheath (Left Limiter)');
legend show; grid on;

%% === Plot Erosion Flux (Ohmic) ===
figure('Name', 'Ohmic Erosion Flux (Left Limiter)');
for i = 1:4
    plot(z_vertL, flux_ohL{i}, 'LineStyle', '-', 'Marker', markers{i}, ...
        'Color', colors{i}, 'DisplayName', ['O^{', species_list{i}(2), '+}']); hold on;
end
plot(z_vertL, flux_ohL_total, 'g*-', 'DisplayName', 'Total');
xlabel('z [m]');
ylabel('Brightness I_\phi(W-I) [m^{-2}s^{-1}St^{-1}]');
title('Ohmic Case: (Left Limiter)');
legend show; grid on;
