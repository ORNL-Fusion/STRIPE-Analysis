close all;
clear;
clc;

%% === Load Species Data for Thermal and Ohmic ===
species = {'o8', 'o7', 'o6', 'o5'};
species_list = {'o8', 'o7', 'o6', 'o5'};

% Custom colors and markers
colors = {[0 0 1], [1 0.5 0], [1 0 1], [0 0 0]};
markers = {'s', 'v', 'o', '+'};

% Data loading
for i = 1:numel(species)
    dataTh{i} = readmatrix(['Targets_thermal_', species{i}, '.txt']);
    dataOh{i} = readmatrix(['Targets_thermal_', species{i},'_ohmic', '.txt']);
end


%% === Define Limiter Regions ===
z_vertL = [0.4; 0.3473; 0.1689; 0.0749; 0.0106; -0.0851; -0.3042; -0.4];
z_vertR = [0.4; 0.3713; 0.1828; 0.0839; 0.0034; -0.0913; -0.2782; -0.4];

thermalL_indices = {
    [3863, 3858, 3859, 3860, 3870, 3476, 3476, 3449, 3450]; 
    [595, 633, 633, 634, 636, 460, 460, 447, 446 ]; 
    [550, 550, 562, 564, 492, 488, 485, 426, 2660];
    [508, 516, 3918, 3920, 3919, 480, 471, 3563, 2715]; 
    [4004, 4003, 12589, 12591, 12592, 3600, 3576, 3575, 2746]; 
    [926, 940, 12621, 12621, 12623, 3631, 3632, 3621, 2779 ];
    [4091, 4110, 4112, 12664, 12662, 12662 3727, 3725, 3702 ]; 
    [4141, 4145, 4145, 4146, 4146, 4147, 3752, 3751, 3749]
};
thermalR_indices = {
    [8045, 8048, 13332, 13330, 13329, 13329, 8680, 8680, 8679 ]; 
    [8043, 8044, 13340, 13328, 13327, 13326, 8640, 8616, 8613 ]; 
    [8801, 8802, 8802, 8575,8575, 8552, 8550, 8550, 11113]; 
    [1549, 1552, 8444, 8528, 8526, 8525, 8519, 8517, 11163];
    [8447, 8447, 8458, 8460, 8494, 8494, 8491, 8489, 11195]; 
    [1269, 1270, 8473, 8474, 8468, 8468, 8328, 8321, 11218];
    [1772, 1770, 1783, 13303, 13304, 8264, 8261, 8238, 8258]; 
    [1745, 1747, 1755, 1754, 8217, 1795,  1794, 1793, 1805] 
};

%% === Preallocate Arrays ===
nRegions = length(z_vertL);
nSpecies = length(species);

% Quantities to compute
niL_Th = zeros(nRegions, nSpecies); niR_Th = niL_Th;
TiL_Th = niL_Th; TiR_Th = niL_Th;
VL_Th  = niL_Th; VR_Th = niL_Th;
fluxL_Th = niL_Th; fluxR_Th = niL_Th;

niL_Oh = niL_Th; niR_Oh = niL_Th;
TiL_Oh = niL_Th; TiR_Oh = niL_Th;
VL_Oh  = niL_Th; VR_Oh = niL_Th;
fluxL_Oh = niL_Th; fluxR_Oh = niL_Th;

%% === Compute Region-Averaged Quantities ===
for s = 1:nSpecies
    dTh = dataTh{s}; dOh = dataOh{s};
    for r = 1:nRegions
        % Thermal
        niL_Th(r,s) = mean(dTh(thermalL_indices{r}, 3));
        niR_Th(r,s) = mean(dTh(thermalR_indices{r}, 3));
        TiL_Th(r,s) = mean(dTh(thermalL_indices{r}, 4));
        TiR_Th(r,s) = mean(dTh(thermalR_indices{r}, 4));
        VL_Th(r,s)  = max(dTh(thermalL_indices{r}, 1));
        VR_Th(r,s)  = max(dTh(thermalR_indices{r}, 1));
        fluxL_Th(r,s) = mean(dTh(thermalL_indices{r},11) .* dTh(thermalL_indices{r},5));
        fluxR_Th(r,s) = mean(dTh(thermalR_indices{r},11) .* dTh(thermalR_indices{r},5));
        % Ohmic
        niL_Oh(r,s) = mean(dOh(thermalL_indices{r}, 11));
        niR_Oh(r,s) = mean(dOh(thermalR_indices{r}, 11));
        TiL_Oh(r,s) = mean(dOh(thermalL_indices{r}, 4));
        TiR_Oh(r,s) = mean(dOh(thermalR_indices{r}, 4));
        VL_Oh(r,s)  = max(dOh(thermalL_indices{r}, 1));
        VR_Oh(r,s)  = max(dOh(thermalR_indices{r}, 1));
        fluxL_Oh(r,s) = mean(dOh(thermalL_indices{r},11) .* dOh(thermalL_indices{r},5));
        fluxR_Oh(r,s) = mean(dOh(thermalR_indices{r},11) .* dOh(thermalR_indices{r},5));
    end
end



%% === Plot Thermal Case ===
figure('Name', 'Thermal: n_i, T_i, Flux, V (Dual-Axis, SI Units)', ...
       'Position', [100, 100, 1000, 800]);

titles = {'(a)', '(b)', '(c)', '(d)'};
subplot_titles = {
    'n_O [m^{-3}]', ...
    'T_O [eV]', ...
    '\Gamma_O (n_i·V) [m^{-2}·s^{-1}]', ...
    'V_{Sheath}^{Thermal} [V]'};

varsL_Th = {niL_Th, TiL_Th, fluxL_Th, VL_Th};
varsR_Th = {niR_Th, TiR_Th, fluxR_Th, VR_Th};

for v = 1:3
    subplot(3,1,v); hold on;

    % Left axis
    yyaxis left;
    for i = 1:length(species_list)
        plot(z_vertL, varsL_Th{v}(:,i), [markers{i} '-'], 'Color', colors{i}, ...
            'DisplayName', ['O^{', species_list{i}(2), '+}']);
    end
    ylabel([subplot_titles{v}, ' (L)']);

    % Right axis
    yyaxis right;
    for i = 1:length(species_list)
        plot(z_vertR, varsR_Th{v}(:,i), [markers{i} '--'], 'Color', colors{i}, ...
            'DisplayName', ['O^{', species_list{i}(2), '+} (R)']);
    end
    ylabel([subplot_titles{v}, ' (R)']);

    if v == 3, xlabel('Vertical position along Q2 limiters: z [m]');
    end
    % title('ICRH Case w/ Thermal Sheath: Plasma Profiles at APL Surface');
     title([titles{v}]);
    legend('Location', 'best', 'NumColumns', 2); grid on; box on; 
end
figure('Name', 'Ohmic: n_i, T_i, Flux, V (Dual-Axis, SI Units)', ...
       'Position', [100, 100, 1000, 800]);

titles = {'(a)', '(b)', '(c)', '(d)'};
subplot_titles = {
    'n_O [m^{-3}]', ...
    'T_O [eV]', ...
    '\Gamma_O [m^{-2}·s^{-1}]', ...
    'V_{Sheath}^{Thermal} [V]'};

varsL_Oh = {niL_Oh, TiL_Oh, fluxL_Oh, VL_Oh};
varsR_Oh = {niR_Oh, TiR_Oh, fluxR_Oh, VR_Oh};

for v = 1:3
    subplot(3,1,v); hold on;

    % Left axis
    yyaxis left;
    for i = 1:length(species_list)
        plot(z_vertL, varsL_Oh{v}(:,i), [markers{i} '-'], 'Color', colors{i}, ...
            'DisplayName', ['O^{', species_list{i}(2), '+}']);
    end
    ylabel([subplot_titles{v}, ' (L)']);

    % Right axis
    yyaxis right;
    for i = 1:length(species_list)
        plot(z_vertR, varsR_Oh{v}(:,i), [markers{i} '--'], 'Color', colors{i}, ...
            'DisplayName', ['O^{', species_list{i}(2), '+} (R)']);
    end
    ylabel([subplot_titles{v}, ' (R)']);

    if v == 3, xlabel('Vertical position along Q2 limiters: z [m]'); 
    end
    title([titles{v}]);
    legend('Location', 'best', 'NumColumns', 2); grid on; box on; 
end

%% Erosion flux calculations

% === Load Thermal and Ohmic Yields ===
yield_files_th = {'yields_thermal_o5.csv', 'yields_thermal_o6.csv', ...
                  'yields_thermal_o7.csv', 'yields_thermal_o8.csv'};
yield_files_oh = {'yields_thermal_o5_ohmic.csv', 'yields_thermal_o6.csv', ...
                  'yields_thermal_o7_ohmic.csv', 'yields_thermal_o8_ohmic.csv'};

yield_data_th = cell(1,4);
yield_data_oh = cell(1,4);
for i = 1:4
    yield_data_th{i} = readmatrix(yield_files_th{i});
    yield_data_oh{i} = readmatrix(yield_files_oh{i});
end

% === Preallocate ===
erosionL_Th = zeros(8, 4);  erosionR_Th = zeros(8, 4);
erosionL_Oh = zeros(8, 4);  erosionR_Oh = zeros(8, 4);

for s = 1:4
    data_T = dataTh{s};
    data_O = dataOh{s};

    % Extract yield columns for each species
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

    % Compute erosion = n_i · v · Y
    erosion_full_T = data_T(:,11) .* abs(data_T(:,5)) .* Y_th;
    erosion_full_O = data_O(:,11) .* abs(data_O(:,5)) .* Y_oh;

    for r = 1:8
        erosionL_Th(r,s) = mean(erosion_full_T(thermalL_indices{r}));
        erosionR_Th(r,s) = mean(erosion_full_T(thermalR_indices{r}));
        erosionL_Oh(r,s) = mean(erosion_full_O(thermalL_indices{r}));
        erosionR_Oh(r,s) = mean(erosion_full_O(thermalR_indices{r}));
    end
end

% === Region-Averaged Erosion Flux (SXB Normalized) ===
sxbValueL = [9.58779; 10.8139; 12.1153; 12.8246; 13.4241; 13.242; 11.1808; 9.70449];

flux_thL = cell(1,4); flux_thR = cell(1,4);
flux_ohL = cell(1,4); flux_ohR = cell(1,4);

for s = 1:4
    flux_thL{s} = erosionL_Th(:,s) ./ (4 * pi * sxbValueL);
    flux_thR{s} = erosionR_Th(:,s) ./ (4 * pi * sxbValueL);
    flux_ohL{s} = erosionL_Oh(:,s) ./ (4 * pi * sxbValueL);
    flux_ohR{s} = 10.*erosionR_Oh(:,s) ./ (4 * pi * sxbValueL);
end

flux_thL_total = flux_thL{1} + flux_thL{2} + flux_thL{3} + flux_thL{4};
flux_thR_total = flux_thR{1} + flux_thR{2} + flux_thR{3} + flux_thR{4};
flux_ohL_total = flux_ohL{1} + flux_ohL{2} + flux_ohL{3} + flux_ohL{4};
flux_ohR_total = flux_ohR{1} + flux_ohR{2} + flux_ohR{3} + flux_ohR{4};


% Custom colors and markers
colors = {[0 0 1], [1 0.5 0], [1 0 1], [0 0 0]};
markers = {'s', 'v', 'o', '+'};

%% === Thermal Erosion Flux ===
figure('Name', 'Thermal Erosion Flux');

% Left Limiter
subplot(2,1,1);
plot(z_vertL, 10*pi.*flux_thL{1}, 'LineStyle', '-', 'Marker', markers{1}, 'Color', colors{1}, 'DisplayName', 'O^{8+} (STRIPE)'); hold on;
plot(z_vertL, 10*pi.*flux_thL{2}, 'LineStyle', '-', 'Marker', markers{2}, 'Color', colors{2}, 'DisplayName', 'O^{7+} (STRIPE)');
plot(z_vertL, 10*pi.*flux_thL{3}, 'LineStyle', '-', 'Marker', markers{3}, 'Color', colors{3}, 'DisplayName', 'O^{6+} (STRIPE)');
plot(z_vertL, 10*pi.*flux_thL{4}, 'LineStyle', '-', 'Marker', markers{4}, 'Color', colors{4}, 'DisplayName', 'O^{5+} (STRIPE)');
plot(z_vertL, 10*pi.*flux_thL_total, 'g*-', 'DisplayName', 'I_{\phi, W-I} (STRIPE)');
xlabel('Vertical position along Q2 limiters: z [m]'); 
ylabel('Brightness I_\phi(W-I) [m^{-2}s^{-1}St^{-1}]');
title('(a) ICRH Case w/ Thermal Sheath: Left Limiter');
legend show; grid on; box on; 

% Right Limiter
subplot(2,1,2);
plot(z_vertR, 100*pi.*flux_thR{1}, 'LineStyle', '-', 'Marker', markers{1}, 'Color', colors{1}, 'DisplayName', 'O^{8+} (STRIPE)'); hold on;
plot(z_vertR, 100*pi.*flux_thR{2}, 'LineStyle', '-', 'Marker', markers{2}, 'Color', colors{2}, 'DisplayName', 'O^{7+} (STRIPE)');
plot(z_vertR, 100*pi.*flux_thR{3}, 'LineStyle', '-', 'Marker', markers{3}, 'Color', colors{3}, 'DisplayName', 'O^{6+} (STRIPE)');
plot(z_vertR, 100*pi.*flux_thR{4}, 'LineStyle', '-', 'Marker', markers{4}, 'Color', colors{4}, 'DisplayName', 'O^{5+} (STRIPE)');
plot(z_vertR, 100*pi.*flux_thR_total, 'g*-', 'DisplayName', 'I_{\phi, W-I} (STRIPE)');
xlabel('Vertical position along Q2 limiters: z [m]'); 
ylabel('Brightness I_{\phi, W-I} [m^{-2}s^{-1}St^{-1}]');
title('(b) ICRH Case w/ Thermal Sheath: Right Limiter');
legend show; grid on; box on; 

%% === Ohmic Erosion Flux ===
figure('Name', 'Ohmic Erosion Flux');

% Left Limiter
subplot(2,1,1);
plot(z_vertL, flux_ohL{1}, 'LineStyle', '-', 'Marker', markers{1}, 'Color', colors{1}, 'DisplayName', 'O^{8+} (STRIPE)'); hold on;
plot(z_vertL, flux_ohL{2}, 'LineStyle', '-', 'Marker', markers{2}, 'Color', colors{2}, 'DisplayName', 'O^{7+} (STRIPE)');
plot(z_vertL, flux_ohL{3}, 'LineStyle', '-', 'Marker', markers{3}, 'Color', colors{3}, 'DisplayName', 'O^{6+} (STRIPE)');
plot(z_vertL, flux_ohL{4}, 'LineStyle', '-', 'Marker', markers{4}, 'Color', colors{4}, 'DisplayName', 'O^{5+} (STRIPE)');
plot(z_vertL, flux_ohL_total, 'g*-', 'DisplayName', 'I_{\phi, W-I} (STRIPE)');
xlabel('Vertical position along Q2 limiters: z [m]');
ylabel('Brightness I_{\phi, W-I} [m^{-2}s^{-1}St^{-1}]');
title('(a) Ohmic Case: Left Limiter');
legend show; grid on; box on; 

% Right Limiter
subplot(2,1,2);
plot(z_vertR, flux_ohR{1}, 'LineStyle', '-', 'Marker', markers{1}, 'Color', colors{1}, 'DisplayName', 'O^{8+} (STRIPE)'); hold on;
plot(z_vertR, flux_ohR{2}, 'LineStyle', '-', 'Marker', markers{2}, 'Color', colors{2}, 'DisplayName', 'O^{7+} (STRIPE)');
plot(z_vertR, flux_ohR{3}, 'LineStyle', '-', 'Marker', markers{3}, 'Color', colors{3}, 'DisplayName', 'O^{6+} (STRIPE)');
plot(z_vertR, flux_ohR{4}, 'LineStyle', '-', 'Marker', markers{4}, 'Color', colors{4}, 'DisplayName', 'O^{5+} (STRIPE)');
plot(z_vertR, flux_ohR_total, 'g*-', 'DisplayName', 'I_{\phi, W-I} (STRIPE)');
xlabel('Vertical position along Q2 limiters: z [m]');
ylabel('Brightness I_{\ph, W-I} [m^{-2}s^{-1}St^{-1}]');
title('(b) Ohmic Case: Right Limiter');
legend show; grid on; box on; 

%% === Raw Erosion Fluxes & LOS-Averaged Yields (No SXB Normalization) ===

% LOS-Averaged Yields
avgY_thL = zeros(8,4); avgY_thR = zeros(8,4);
avgY_ohL = zeros(8,4); avgY_ohR = zeros(8,4);

for s = 1:4
    data_T = dataTh{s};
    data_O = dataOh{s};

    % Extract same Y arrays as used before
    if s == 1
        Y_th = [0; yield_data_th{4}(:,8)];
        Y_oh = [0; yield_data_oh{4}(:,8)];
    elseif s == 2
        Y_th = [0; yield_data_th{3}(:,7)];
        Y_oh = [0; yield_data_oh{3}(:,7)];
    elseif s == 3
        Y_th = 0.5.*[0; yield_data_th{2}(:,6)];
        Y_oh = 0.5.*[0; yield_data_oh{2}(:,6)];
    else
        Y_th = 0.5.*[0; yield_data_th{1}(:,5)];
        Y_oh = 0.5.*[0; yield_data_oh{1}(:,5)];
    end

    for r = 1:8
        avgY_thL(r,s) = mean(Y_th(thermalL_indices{r}));
        avgY_thR(r,s) = mean(Y_th(thermalR_indices{r}));
        avgY_ohL(r,s) = mean(Y_oh(thermalL_indices{r}));
        avgY_ohR(r,s) = mean(Y_oh(thermalR_indices{r}));
    end
end

%% === Thermal: Yields and Erosion Fluxes (Raw, No S/XB) ===
figure('Name', 'Thermal: Yields and Erosion Fluxes (Left vs Right Limiters)', ...
       'Position', [100, 100, 1000, 800]);

%% === (a) Sputtering Yields ===
subplot(2,1,1); hold on;

yyaxis left
for i = 1:4
    Z = species_list{i}(2);
    plot(z_vertL, avgY_thL(:,i), [markers{i} '-'], ...
        'Color', colors{i}, 'LineWidth', 1.5, ...
        'DisplayName', ['$\rm O^{', Z, '+}$ (L)']);
end
ylabel('Y_{eff} (Left)');
ylim([0, 5e-3]);

yyaxis right
for i = 1:4
    Z = species_list{i}(2);
    plot(z_vertR, avgY_thR(:,i), [markers{i} '--'], ...
        'Color', colors{i}, 'LineWidth', 1.5, ...
        'DisplayName', ['$\rm O^{', Z, '+}$ (R)']);
end
ylabel('Y_{eff} (Right)');
xlabel('z [m]');
ylim([0, 5e-3]);
title('(a)');
legend('Location', 'northeastoutside', 'NumColumns', 2, 'Interpreter', 'latex');
grid on; box on;

%% === (b) Erosion Fluxes (Raw, No S/XB) - Log Scale ===
flux_eroded_raw_thL = cell(1,4); flux_eroded_raw_thR = cell(1,4);
total_flux_thL = zeros(8,1); total_flux_thR = zeros(8,1);
for i = 1:4
    flux_eroded_raw_thL{i} = erosionL_Th(:,i);
    flux_eroded_raw_thR{i} = erosionR_Th(:,i);
    total_flux_thL = total_flux_thL + flux_eroded_raw_thL{i};
    total_flux_thR = total_flux_thR + flux_eroded_raw_thR{i};
end

subplot(2,1,2); hold on;

yyaxis left
for i = 1:4
    Z = species_list{i}(2);
    plot(z_vertL, flux_eroded_raw_thL{i}, [markers{i} '-'], ...
        'Color', colors{i}, 'LineWidth', 1.5, ...
        'DisplayName', ['$\rm O^{', Z, '+}$ (L)']);
end
plot(z_vertL, total_flux_thL, 'g-*', 'LineWidth', 1.8, ...
    'DisplayName', 'Total (L)');
ylabel('\Gamma_{gross, W} [m^{-2}s^{-1}] (Left)');
set(gca, 'YScale', 'log');
ylim([1e10, 1e18]);

yyaxis right
for i = 1:4
    Z = species_list{i}(2);
    plot(z_vertR, flux_eroded_raw_thR{i}, [markers{i} '--'], ...
        'Color', colors{i}, 'LineWidth', 1.5, ...
        'DisplayName', ['$\rm O^{', Z, '+}$ (R)']);
end
plot(z_vertR, total_flux_thR, 'g--*', 'LineWidth', 1.8, ...
    'DisplayName', 'Total (R)');
ylabel('\Gamma_{gross, W} [m^{-2}s^{-1}] (Right)');
set(gca, 'YScale', 'log');
ylim([1e10, 1e18]);

xlabel('z [m]');
title('(b)');
legend('Location', 'northeastoutside', 'NumColumns', 2, 'Interpreter', 'latex');
grid on; box on;

%% === Ohmic: Yields and Erosion Fluxes (Raw, No S/XB) ===
figure('Name', 'Ohmic: Yields and Erosion Fluxes (Left vs Right Limiters)', ...
       'Position', [100, 100, 1000, 800]);

%% === (a) Sputtering Yields ===
subplot(2,1,1); hold on;

yyaxis left
for i = 1:4
    Z = species_list{i}(2);
    plot(z_vertL, avgY_ohL(:,i), [markers{i} '-'], ...
        'Color', colors{i}, 'LineWidth', 1.5, ...
        'DisplayName', ['$\rm O^{', Z, '+}$ (L)']);
end
ylabel('Y_{eff} (Left)');
ylim([0, 5e-3]);

yyaxis right
for i = 1:4
    Z = species_list{i}(2);
    plot(z_vertR, avgY_ohR(:,i), [markers{i} '--'], ...
        'Color', colors{i}, 'LineWidth', 1.5, ...
        'DisplayName', ['$\rm O^{', Z, '+}$ (R)']);
end
ylabel('Y_{eff} (Right)');
xlabel('z [m]');
ylim([0, 5e-3]);
title('(a)');
legend('Location', 'northeastoutside', 'NumColumns', 2, 'Interpreter', 'latex');
grid on; box on;

%% === (b) Erosion Fluxes (Raw, No S/XB) - Log Scale ===
flux_eroded_raw_ohL = cell(1,4); flux_eroded_raw_ohR = cell(1,4);
total_flux_ohL = zeros(8,1); total_flux_ohR = zeros(8,1);
for i = 1:4
    flux_eroded_raw_ohL{i} = erosionL_Oh(:,i);
    flux_eroded_raw_ohR{i} = erosionR_Oh(:,i);
    total_flux_ohL = total_flux_ohL + flux_eroded_raw_ohL{i};
    total_flux_ohR = total_flux_ohR + flux_eroded_raw_ohR{i};
end

subplot(2,1,2); hold on;

yyaxis left
for i = 1:4
    Z = species_list{i}(2);
    plot(z_vertL, flux_eroded_raw_ohL{i}, [markers{i} '-'], ...
        'Color', colors{i}, 'LineWidth', 1.5, ...
        'DisplayName', ['$\rm O^{', Z, '+}$ (L)']);
end
plot(z_vertL, total_flux_ohL, 'g-*', 'LineWidth', 1.8, ...
    'DisplayName', 'Total (L)');
ylabel('\Gamma_{gross, W} [m^{-2}s^{-1}] (Left)');
set(gca, 'YScale', 'log');
ylim([1e10, 1e18]);

yyaxis right
for i = 1:4
    Z = species_list{i}(2);
    plot(z_vertR, flux_eroded_raw_ohR{i}, [markers{i} '--'], ...
        'Color', colors{i}, 'LineWidth', 1.5, ...
        'DisplayName', ['$\rm O^{', Z, '+}$ (R)']);
end
plot(z_vertR, total_flux_ohR, 'g--*', 'LineWidth', 1.8, ...
    'DisplayName', 'Total (R)');
ylabel('\Gamma_{gross, W} [m^{-2}s^{-1}] (Right)');
set(gca, 'YScale', 'log');
ylim([1e10, 1e18]);

xlabel('z [m]');
title('(b)');
legend('Location', 'northeastoutside', 'NumColumns', 2, 'Interpreter', 'latex');
grid on; box on;