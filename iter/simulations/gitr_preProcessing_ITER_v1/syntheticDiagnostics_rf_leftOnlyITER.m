close all; clear; clc;

%% === Load Input Data ===
dataRF = readmatrix('Targets.txt');

% Load sputtering yields for Neon charge states
yields_Ne1 = readmatrix('yields_Ne+.csv');
yields_Ne2 = readmatrix('yields_Ne2+.csv');
yields_Ne3 = readmatrix('yields_Ne3+.csv');
yields_Ne4 = readmatrix('yields_Ne3+.csv');
yields_Ne5 = readmatrix('yields_Ne5+.csv');
yields_Ne6 = readmatrix('yields_Ne6+.csv');
yields_Ne7 = readmatrix('yields_Ne7+.csv');

% Neon fractional contributions (example values, update if needed)
ni_dataRF_Ne1 = 0.02 * 0.23 * dataRF(:,11);
ni_dataRF_Ne2 = 0.02 * 0.26 * dataRF(:,11);
ni_dataRF_Ne3 = 0.02 * 0.205 * dataRF(:,11);
ni_dataRF_Ne4 = 0.02 * 0.12 * dataRF(:,11);
ni_dataRF_Ne5 = 0.02 * 0.067 * dataRF(:,11);
ni_dataRF_Ne6 = 0.02 * 0.031 * dataRF(:,11);
ni_dataRF_Ne7 = 0.02 * 0.013 * dataRF(:,11);

% Velocity
v_dataRF_Ne = abs(dataRF(:,5));

%% === Erosion Flux Computation ===
ero_dataRF_Ne1 = ni_dataRF_Ne1 .* v_dataRF_Ne .* [0;yields_Ne1(:,1)];
ero_dataRF_Ne2 = ni_dataRF_Ne2 .* v_dataRF_Ne .* [0;yields_Ne2(:,2)];
ero_dataRF_Ne3 = ni_dataRF_Ne3 .* v_dataRF_Ne .* [0;yields_Ne3(:,3)];
ero_dataRF_Ne4 = ni_dataRF_Ne4 .* v_dataRF_Ne .* [0;yields_Ne4(:,3)];
ero_dataRF_Ne5 = ni_dataRF_Ne5 .* v_dataRF_Ne .* [0;yields_Ne5(:,5)];
ero_dataRF_Ne6 = ni_dataRF_Ne6 .* v_dataRF_Ne .* [0;yields_Ne6(:,6)];
ero_dataRF_Ne7 = ni_dataRF_Ne7 .* v_dataRF_Ne .* [0;yields_Ne7(:,7)];

%% === Left Limiter Indices (10-segment version) ===
rfL_indices = {
    [9284,9220,9151,9087,19541,19202,19201,19166,19165,19162];
    [9295,9231,9168,9104,19556,19551,19505,19227,19221,19277];
    [9311,9247,9182,9120,19582,19567,19562,19569,19310,19311];
    [9308,9244,9180,9116,19599,19593,19431,19422,19419,19401];
    [9330,9266,9194,9130,19622,19490,19488,19455,19444,19441];
    [8729,8665,8601,8537,18470,18654,18560,18587,18570,18571];
    [8720,8656,8590,8523,18453,18451,18496,18502,18498,18483];
    [8712,8647,8586,8522,18425,18411,18415,18419,18216,18197];
    [8705,8641,8576,8511,18370,18367,18379,18308,18287,18125];
    [8736,8637,8573,8508,18350,18336,18331,18328,18245,18244];
};

z_vertL = linspace(-0.3, 1.6, 10);  % Adjust as needed
sxbValueL = [9.58, 10.81, 12.11, 12.82, 13.42, 13.24, 11.18, 9.70, 9.10, 8.80];  % Extend to 10 regions

%% === Average Brightness per Region ===
data_rfL = {ero_dataRF_Ne1, ero_dataRF_Ne2, ero_dataRF_Ne3, ...
            ero_dataRF_Ne4, ero_dataRF_Ne5, ero_dataRF_Ne6, ero_dataRF_Ne7};

flux_rfL = cell(1,7);
for i = 1:7
    temp = zeros(10,1);
    for j = 1:10
        temp(j) = mean(data_rfL{i}(rfL_indices{j}));
    end
    % flux_rfL{i} = temp ./ (4 * pi * sxbValueL(:));
    flux_rfL{i} = temp;
end

flux_total_L = zeros(10,1);
for i = 1:7
    flux_total_L = flux_total_L + flux_rfL{i};
end

%% === Plot Neon Erosion Flux ===
species_labels = {'Ne^+', 'Ne^{2+}', 'Ne^{3+}', 'Ne^{4+}', ...
                  'Ne^{5+}', 'Ne^{6+}', 'Ne^{7+}'};
colors = lines(7);
markers = {'o','s','^','v','d','+','*'};

figure('Name','Neon RF-Induced Erosion Flux (Left Limiter, 10 Zones)');
hold on;
for i = 1:7
    plot(z_vertL, flux_rfL{i}, ...
        'Color', colors(i,:), ...
        'Marker', markers{i}, ...
        'LineStyle', '-', ...
        'LineWidth', 1.5, ...
        'DisplayName', species_labels{i});
end
plot(z_vertL, flux_total_L, 'k*-', 'LineWidth', 2, 'DisplayName', 'Total Neon Erosion');

xlabel('z [m]');
ylabel('I_\phi (Ne) [m^{-2}s^{-1}St^{-1}]');
title('Neon-Induced W Erosion Flux (RF Sheath, Left Limiter, 10 Segments)');
legend('Location','northeastoutside');
grid on; box on;

%% === Plasma LOS-Averaged Quantities ===
ne_rfL = zeros(10, 1);
Te_rfL = zeros(10, 1);
nev_rfL = zeros(10, 1);
Vsh_rfL = zeros(10, 1);  % RF sheath voltage estimate

for j = 1:10
    indices = rfL_indices{j};
    ne_rfL(j) = mean(dataRF(indices, 2));        % [m^-3], factor 50 from scaling
    Te_rfL(j) = mean(dataRF(indices, 3));              % [eV]
    flow = mean(dataRF(indices, 5));                   % [m/s]
    nev_rfL(j) = ne_rfL(j) * flow;                     % [m^-2 s^-1]
    Vsh_rfL(j) = max(dataRF(indices, 1));                        % Simple Bohm sheath approx
end

%% === Plot: Plasma Profiles Along Left Limiter ===
figure('Name', 'LOS-Averaged Plasma Parameters (Left Limiter)', ...
       'Position', [100, 100, 1000, 800]);

subplot(4,1,1);
plot(z_vertL, ne_rfL, 'b-o', 'LineWidth', 1.5);
ylabel('n_e [m^{-3}]'); grid on; title('LOS-Averaged Plasma Profiles (Left Limiter)');

subplot(4,1,2);
plot(z_vertL, Te_rfL, 'r-s', 'LineWidth', 1.5);
ylabel('T_e [eV]'); grid on;
ylim([0 15])

subplot(4,1,3);
plot(z_vertL, nev_rfL, 'm-^', 'LineWidth', 1.5);
ylabel('n_ev [m^{-2}s^{-1}]'); grid on;

subplot(4,1,4);
plot(z_vertL, Vsh_rfL, 'k-d', 'LineWidth', 1.5);
ylabel('V_{sheath}^{RF} [V]');
xlabel('z [m]');
ylim([0 3000])
grid on;