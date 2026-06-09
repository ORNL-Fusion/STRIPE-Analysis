close all;
clear;
clc;

%% === Load Data ===
dataRF = readmatrix('Targets.txt');

yields_Ne1 = readmatrix('yields_Ne+.csv');
yields_Ne2 = readmatrix('yields_Ne2+.csv');
yields_Ne3 = readmatrix('yields_Ne3+.csv');
yields_Ne4 = readmatrix('yields_Ne3+.csv');
yields_Ne5 = readmatrix('yields_Ne5+.csv');
yields_Ne6 = readmatrix('yields_Ne6+.csv');
yields_Ne7 = readmatrix('yields_Ne7+.csv');

%% === Extract and Compute ===



ni_dataRF_Ne1 = 0.02*0.23.*dataRF(:,11);
ni_dataRF_Ne2 = 0.02*0.26.*dataRF(:,11);
ni_dataRF_Ne3 = 0.02*0.205378.*dataRF(:,11);
ni_dataRF_Ne4 = 0.02*0.205378.*dataRF(:,11);
ni_dataRF_Ne5 = 0.02*0.067.*dataRF(:,11);
ni_dataRF_Ne6 = 0.02*0.031.*dataRF(:,11);
ni_dataRF_Ne7 = 0.02*0.013.*dataRF(:,11);

v_dataRF_Ne = abs(dataRF(:,5));




ero_dataRF_Ne1 = ni_dataRF_Ne1 .* v_dataRF_Ne .* [0;yields_Ne1(:,1)];
ero_dataRF_Ne2 = ni_dataRF_Ne2 .* v_dataRF_Ne .* [0;yields_Ne2(:,2)];
ero_dataRF_Ne3 = ni_dataRF_Ne3 .* v_dataRF_Ne .* [0;yields_Ne3(:,3)];
ero_dataRF_ne4 = ni_dataRF_Ne4 .* v_dataRF_Ne .* [0;yields_Ne4(:,4)];
ero_dataRF_Ne5 = ni_dataRF_Ne5 .* v_dataRF_Ne .* [0;yields_Ne5(:,5)];
ero_dataRF_Ne6 = ni_dataRF_Ne6 .* v_dataRF_Ne .* [0;yields_Ne6(:,6)];
ero_dataRF_Ne7 = ni_dataRF_Ne7 .* v_dataRF_Ne .* [0;yields_Ne7(:,7)];


ero_dataRF_total = ero_dataRF_Ne1 + ero_dataRF_Ne2 + ero_dataRF_Ne3 + ero_dataRF_Ne4+ ero_dataRF_Ne5+ero_dataRF_Ne6+ero_dataRF_Ne7;

%% === Load Geometry ===
if ~exist('x1', 'var')
    fid = fopen('gitrGeometryPointPlane3d_comsol.cfg');
    for i = 1:20
        tline = fgetl(fid);
        if i > 2
            evalc(tline);
        end
    end
    fclose(fid);
end

subset = 1:length(x1);
X = [x1(subset), x2(subset), x3(subset)]';
Y = [y1(subset), y2(subset), y3(subset)]';
Z = [z1(subset), z2(subset), z3(subset)]';

ni_data = dataRF_o8(:,11);
te_data = dataRF_o8(:,3);
flux = v_dataRF_o8;

%% === Plot Figures ===
% load('coolwarm.mat', 'coolwarm_rgb')
% 
% figure(1);
% patch(X, Y, Z, ni_data, 'FaceAlpha', 1, 'EdgeColor', 'none');
% title('Density'); xlabel('X'); ylabel('Y'); zlabel('Z'); colormap(coolwarm_rgb); colorbar;
% 
% figure(2);
% patch(X, Y, Z, te_data, 'FaceAlpha', 0.3, 'EdgeAlpha', 0.3);
% title('Temperature'); xlabel('X'); ylabel('Y'); zlabel('Z'); colorbar;
% 
% figure(3);
% patch(X, Y, Z, ero_dataRF_total, 'FaceAlpha', 1, 'EdgeAlpha', 0.05);
% title('Erosion Flux'); xlabel('X'); ylabel('Y'); zlabel('Z'); colormap(coolwarm_rgb); colorbar;

%% === Average Over Regions ===
z_vertL = [-0.3000,   -0.0889,    0.1222,    0.3333,    0.5444,    0.7556,    0.9667,    1.1778,    1.3889,    1.6000];
z_vertR = [0.4; 0.3713; 0.1828; 0.0839; 0.0034; -0.0913; -0.2782; -0.4];
sxbValueL = [9.58779;10.8139;12.1153;12.8246;13.4241;13.242;11.1808;9.70449];

rfL_indices = {
    [9284,9220,9151, 9087, 19541, 19202, 19201, 19166, 19165, 19162];
    [9295,9231, 9168, 9104, 19556, 19551, 19505, 19227, 19221, 19277];
    [9311, 9247, 9182, 9120, 19582, 19567, 19562, 19569, 19310, 19311];
    [9308, 9244, 9180, 9116, 19599, 19593, 19431, 19422, 19419, 19401];
    [9330, 9266, 9194, 9130, 19622, 19490, 19488, 19455, 19444, 19441];
    [8729, 8665, 8601, 8537, 18470, 18654, 18560, 18587, 18570, 18571];
    [8720, 8656, 8590, 8523, 18453, 18451, 18496, 18502, 18498, 18483];
    [8712, 8647, 8586, 8522, 18425, 18411, 18415, 18419, 18216, 18197];
    [8705, 8641, 8576, 8511, 18370, 18367, 18379, 18308, 18287, 18125];
    [8736, 8637, 8573, 8508, 18350, 18336, 18331, 18328,18245, 18244 ];
};

thermalL_indices = {
    [9284,9220,9151, 9087, 19541, 19202, 19201, 19166, 19165, 19162];
    [9295,9231, 9168, 9104, 19556, 19551, 19505, 19227, 19221, 19277];
    [9311, 9247, 9182, 9120, 19582, 19567, 19562, 19569, 19310, 19311];
    [9308, 9244, 9180, 9116, 19599, 19593, 19431, 19422, 19419, 19401];
    [9330, 9266, 9194, 9130, 19622, 19490, 19488, 19455, 19444, 19441];
    [8729, 8665, 8601, 8537, 18470, 18654, 18560, 18587, 18570, 18571];
    [8720, 8656, 8590, 8523, 18453, 18451, 18496, 18502, 18498, 18483];
    [8712, 8647, 8586, 8522, 18425, 18411, 18415, 18419, 18216, 18197];
    [8705, 8641, 8576, 8511, 18370, 18367, 18379, 18308, 18287, 18125];
    [8736, 8637, 8573, 8508, 18350, 18336, 18331, 18328,18245, 18244 ];
};

data_rfL = {ero_dataRF_o8, ero_dataRF_o7, ero_dataRF_o6, ero_dataRF_o5};
flux_rfL = cell(1,4);

for i = 1:4
    temp = zeros(8,1);
    for j = 1:8
        temp(j) = mean(data_rfL{i}(rfL_indices{j}));
    end
    flux_rfL{i} = temp ./ (4 * pi * sxbValueL);
end

flux_total_L = abs(flux_rfL{1} + flux_rfL{2} + flux_rfL{3} + flux_rfL{4});


data_rfR = {ero_dataRF_o8, ero_dataRF_o7, ero_dataRF_o6, ero_dataRF_o5};
flux_rfR = cell(1,4);

for i = 1:4
    temp = zeros(8,1);
    for j = 1:8
        temp(j) = mean(data_rfR{i}(rfR_indices{j}));
    end
    flux_rfR{i} = temp ./ (4 * pi * sxbValueL);
end

flux_total_R = abs(flux_rfR{1} + flux_rfR{2} + flux_rfR{3} + flux_rfR{4});



%% === Average Over Left Limiter Regions Only ===
species_list = {'o8', 'o7', 'o6', 'o5'};
dataRF_list = {dataRF_o8, dataRF_o7, dataRF_o6, dataRF_o5};
colors = {[0 0 1], [1 0.5 0], [1 0 1], [0 0 0]};
markers = {'s', 'v', 'o', '+'};

ni_L_all = zeros(8, 4);
Ti_L_all = zeros(8, 4);
Te_L_all = zeros(8, 4);
V_L_all  = zeros(8, 4);
flux_L_all = zeros(8, 4);
V_th_sheath_L = zeros(8, 1);  % 3Te ref from O8+

for i = 1:4
    data = dataRF_list{i};
    for j = 1:8
        ni_L_all(j,i) = 50.*mean(data(rfL_indices{j},11));
        Ti_L_all(j,i) = 3.*mean(data(rfL_indices{j},4));
        Te_L_all(j,i) = mean(data(rfL_indices{j},3));
        V_L_all(j,i)  = max(data(rfL_indices{j},1));
        flow_L        = mean(data(rfL_indices{j},5));
        flux_L_all(j,i) = ni_L_all(j,i) * flow_L;
    end
end

v_data = readmatrix('Targets_thermal_o8.txt');

V_thermal_L = zeros(8,1); V_thermal_R = zeros(8,1);
 
for j = 1:8
    V_thermal_L(j) = 3*max(v_data(thermalL_indices{j}, 1));   % Plasma potential [V]
    V_thermal_R(j) = 3*max(v_data(thermalR_indices{j}, 1));
end

%% === Erosion Flux Calculation ===
sxbValueL = [9.58779;10.8139;12.1153;12.8246;13.4241;13.242;11.1808;9.70449];
data_rfL = {ero_dataRF_o8, ero_dataRF_o7, ero_dataRF_o6, ero_dataRF_o5};
flux_rfL = cell(1,4);

for i = 1:4
    temp = zeros(8,1);
    for j = 1:8
        temp(j) = mean(data_rfL{i}(rfL_indices{j}));
    end
    flux_rfL{i} = temp ./ (4 * pi * sxbValueL);
end
flux_total_L = abs(flux_rfL{1} + flux_rfL{2} + flux_rfL{3} + flux_rfL{4});

%% === Plot Erosion Flux (Left Limiter Only) ===
% Custom colors and markers
colors = {[0 0 1], [1 0.5 0], [1 0 1], [0 0 0]};
markers = {'s', 'v', 'o', '+'};

figure('Name', 'RF-Induced Erosion Flux (Left Limiter)');

plot(z_vertL, flux_rfL{1}, 'LineStyle', '-', 'Marker', markers{1}, 'Color', colors{1}, 'DisplayName', 'O^{8+} (STRIPE)'); hold on;
plot(z_vertL, flux_rfL{2}, 'LineStyle', '-', 'Marker', markers{2}, 'Color', colors{2}, 'DisplayName', 'O^{7+} (STRIPE)');
plot(z_vertL, flux_rfL{3}, 'LineStyle', '-', 'Marker', markers{3}, 'Color', colors{3}, 'DisplayName', 'O^{6+} (STRIPE)');
plot(z_vertL, flux_rfL{4}, 'LineStyle', '-', 'Marker', markers{4}, 'Color', colors{4}, 'DisplayName', 'O^{5+} (STRIPE)');
plot(z_vertL, flux_total_L, 'g*-', 'DisplayName', 'I_{\phi, W-I} (STRIPE)');

xlabel('Vertical position along Q2 Left Limiter: z [m]');
ylabel('Brightness I_\phi(W-I) [m^{-2}s^{-1}St^{-1}]');
title('(b) ICRH Case w/ RF Sheath (Left Limiter)');
legend show;
grid on; box on; 

%% === Plot Plasma Profiles (Left Limiter Only) ===
figure('Name', 'RF Plasma Profiles (Left Limiter Only)', ...
       'Position', [100, 100, 1000, 800]);

titles = {'(a)', '(b)', '(c)', '(d)'};
subplot_titles = {
    'n_O [m^{-3}]', ...
    'T_O [eV]', ...
    '\Gamma_O [m^{-2}·s^{-1}]', ...
    'V_{Sheath}^{RF} [V]'};

% ---- Ion Density ----
subplot(4,1,1); hold on;
for i = 1:4
    plot(z_vertL, ni_L_all(:,i), [markers{i} '-'], 'Color', colors{i}, ...
        'DisplayName', ['O^{', species_list{i}(2), '+}']);
end
ylabel(subplot_titles{1}); title([titles{1}]);
% title(['RF Case: ', subplot_titles{1}]);
legend('Location', 'best'); grid on; box on; 
% ---- Ion Temperature ----
subplot(4,1,2); hold on;
for i = 1:4
    plot(z_vertL, Ti_L_all(:,i), [markers{i} '-'], 'Color', colors{i}, ...
        'DisplayName', ['O^{', species_list{i}(2), '+}']);
end
ylabel(subplot_titles{2}); title([titles{2}]);
% title(['RF Case: ', subplot_titles{2}]);
legend('Location', 'best'); grid on; box on; 

% ---- Flux ----
subplot(4,1,3); hold on;
for i = 1:4
    plot(z_vertL, flux_L_all(:,i), [markers{i} '-'], 'Color', colors{i}, ...
        'DisplayName', ['O^{', species_list{i}(2), '+}']);
end
ylabel(subplot_titles{3}); title([titles{3}]);
% title(['RF Case: ', subplot_titles{3}]);
legend('Location', 'best'); grid on; box on; 

% ---- Plasma Potential + Thermal Sheath ----
subplot(4,1,4); hold on;
plot(z_vertL, V_L_all(:,1), [markers{1} '-'], 'Color', colors{1}, ...
    'DisplayName', 'RF Sheath');
plot(z_vertL, V_thermal_L, 'k--', 'LineWidth', 1.2, ...
    'DisplayName', 'Thermal Sheath');
ylabel(subplot_titles{4}); xlabel('Vertical position along Q2 Left Limiter: z [m]');
% title('RF Case: Plasma Potential + Thermal Sheath Reference');
    title([titles{4}]);
legend('Location', 'best'); grid on; box on; 
%% === Plot Yields at RF Indices (Left and Right Limiters) ===

% Define yield data per species
yield_data_list = {
    yields_dataRF_o8(:,8), ...
    yields_dataRF_o7(:,7), ...
    0.5 .* yields_dataRF_o6(:,6), ...
    0.5 .* yields_dataRF_o5(:,5)
};

% Left Limiter
yields_rfL = cell(1, 4);
for i = 1:4
    data_shifted = [0; 0; yield_data_list{i}];
    temp_yield = zeros(8, 1);
    for j = 1:8
        temp_yield(j) = mean(data_shifted(rfL_indices{j}));
    end
    yields_rfL{i} = temp_yield;
end

% Right Limiter
yields_rfR = cell(1, 4);
for i = 1:4
    data_shifted = [0; 0; yield_data_list{i}];
    temp_yield = zeros(8, 1);
    for j = 1:8
        temp_yield(j) = mean(data_shifted(rfR_indices{j}));
    end
    yields_rfR{i} = temp_yield;
end

%% === Raw Erosion Fluxes (No S/XB) for Both Limiters ===

% Left
flux_eroded_raw_rfL = cell(1, 4);
for i = 1:4
    temp = zeros(8, 1);
    for j = 1:8
        temp(j) = mean(data_rfL{i}(rfL_indices{j}));
    end
    flux_eroded_raw_rfL{i} = temp;
end

% Right
flux_eroded_raw_rfR = cell(1, 4);
for i = 1:4
    temp = zeros(8, 1);
    for j = 1:8
        temp(j) = mean(data_rfR{i}(rfR_indices{j}));
    end
    flux_eroded_raw_rfR{i} = temp;
end

%% === Plot Yields and Erosion Fluxes: Left vs Right Limiters on Dual Y-Axes ===

figure('Name', 'Yields and Erosion Fluxes (Dual Axis: Left vs Right Limiters)', ...
       'Position', [100, 100, 1000, 800]);

% === (a) Yields: Left vs Right ===
subplot(2,1,1); hold on;

yyaxis left
for i = 1:4
    plot(z_vertL, yields_rfL{i}, [markers{i} '-'], 'Color', colors{i}, ...
        'LineWidth', 1.5, 'DisplayName', ['Left: O^{', species_list{i}(2), '+}']);
end
ylabel('Yield (Left Limiter)');
ylim([0, max(cellfun(@(x) max(x), yields_rfL)) * 1.2]);

yyaxis right
for i = 1:4
    plot(z_vertL, yields_rfR{i}, [markers{i} '--'], 'Color', colors{i}, ...
        'LineWidth', 1.5, 'DisplayName', ['Right: O^{', species_list{i}(2), '+}']);
end
ylabel('Yield (Right Limiter)');
xlabel('z [m]');
title('(a) Sputtering Yields (Left vs Right Limiter)');
legend('Location', 'northeastoutside'); grid on; box on;

% === (b) Erosion Flux: Left vs Right ===
subplot(2,1,2); hold on;

yyaxis left
for i = 1:4
    plot(z_vertL, flux_eroded_raw_rfL{i}, [markers{i} '-'], 'Color', colors{i}, ...
        'LineWidth', 1.5, 'DisplayName', ['Left: O^{', species_list{i}(2), '+}']);
end
ylabel('I_\phi (Left Limiter) [m^{-2}s^{-1}]');
ylim([0, max(cellfun(@(x) max(x), flux_eroded_raw_rfL)) * 1.2]);

yyaxis right
for i = 1:4
    plot(z_vertL, flux_eroded_raw_rfR{i}, [markers{i} '--'], 'Color', colors{i}, ...
        'LineWidth', 1.5, 'DisplayName', ['Right: O^{', species_list{i}(2), '+}']);
end
ylabel('I_\phi (Right Limiter) [m^{-2}s^{-1}]');
xlabel('z [m]');
title('(b) Erosion Fluxes (Left vs Right Limiter)');
legend('Location', 'northeastoutside'); grid on; box on;