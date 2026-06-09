close all;
clear;
clc;

%% === Load Data ===
dataRF_o8 = readmatrix('Targets_comsol_o8_new.txt');
dataRF_o7 = readmatrix('Targets_comsol_o7_new.txt');
dataRF_o6 = readmatrix('Targets_comsol_o6_new.txt');
dataRF_o5 = readmatrix('Targets_comsol_o5_new.txt');

yields_dataRF_o8 = readmatrix('yields_comsol_o8_new.csv');
yields_dataRF_o7 = readmatrix('yields_comsol_o7_new.csv');
yields_dataRF_o6 = readmatrix('yields_comsol_o7_new.csv');
yields_dataRF_o5 = readmatrix('yields_comsol_o5_new.csv');

%% === Extract and Compute ===

ni_dataRF_o8 = dataRF_o8(:,11);
ni_dataRF_o7 = dataRF_o7(:,11);
ni_dataRF_o6 = dataRF_o6(:,11);
ni_dataRF_o5 = dataRF_o5(:,11);

v_dataRF_o8 = abs(dataRF_o8(:,5));
v_dataRF_o7 = abs(dataRF_o7(:,5));
v_dataRF_o6 = abs(dataRF_o6(:,5));
v_dataRF_o5 = abs(dataRF_o5(:,5));



% ero_dataRF_o8 = ni_dataRF_o8 .* v_dataRF_o8;
% ero_dataRF_o7 = ni_dataRF_o7 .* v_dataRF_o7;
% ero_dataRF_o6 = ni_dataRF_o6 .* v_dataRF_o6;
% ero_dataRF_o5 = ni_dataRF_o5 .* v_dataRF_o5;

% ero_dataRF_o8 = ni_dataRF_o8 ;
% ero_dataRF_o7 = ni_dataRF_o7 ;
% ero_dataRF_o6 = ni_dataRF_o6 ;
% ero_dataRF_o5 = ni_dataRF_o5 ;
% 
% ero_dataRF_o8 = v_dataRF_o8;
% ero_dataRF_o7 = v_dataRF_o7;
% ero_dataRF_o6 = v_dataRF_o6;
% ero_dataRF_o5 = v_dataRF_o5;

ero_dataRF_o8 = v_dataRF_o8 .* [0;yields_dataRF_o8(:,8)];
ero_dataRF_o7 = v_dataRF_o7 .* [0;yields_dataRF_o7(:,7)];
ero_dataRF_o6 = v_dataRF_o6 .* [0;yields_dataRF_o7(:,7)];
ero_dataRF_o5 = v_dataRF_o5 .* [0;yields_dataRF_o5(:,5)];

ero_dataRF_total = ero_dataRF_o8 + ero_dataRF_o7 + ero_dataRF_o6 + ero_dataRF_o5;

%% === Load Geometry ===
if ~exist('x1', 'var')
    fid = fopen('gitrGeometryPointPlane3d.cfg');
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
z_vertL = [0.4; 0.3473; 0.1689; 0.0749; 0.0106; -0.0851; -0.3042; -0.4];
z_vertR = [0.4; 0.3713; 0.1828; 0.0839; 0.0034; -0.0913; -0.2782; -0.4];
sxbValueL = [9.58779;10.8139;12.1153;12.8246;13.4241;13.242;11.1808;9.70449];

rfL_indices = {
    [259, 260, 6186, 6193, 5696, 5696, 5686, 5685, 5564];
    [479, 478, 19161, 19163, 5692, 5676, 5675, 5673,994];
    [5025, 5028, 19122, 19121, 19110, 5638, 5637, 5626, 5224];
    [5078, 5598, 19171, 19169, 5612, 5602, 5604, 5603, 552];
    [434, 435, 448, 5606, 5605, 5608, 19000, 5304, 5296];
    [4683, 4710, 18923, 18904, 18911, 18912, 18912, 5352, 5352];
    [427, 5505, 5508, 18972, 5484, 5483, 5481, 525, 517];
    [400, 19002, 19008, 19064, 19062, 19058, 5780, 512, 510];
};

rfR_indices = {
    [12766, 21433, 21438, 21428, 21380, 21377, 21373, 12032, 12028];
    [12783, 21425, 21430, 21431, 21386, 21385, 21381, 12038,12040];
    [1697, 1714, 1715, 21254, 21204, 21202, 21201, 12342, 12332];
    [12379, 12390, 21181, 21182, 21251, 21251, 21218, 12320, 12312];
    [12509, 1690, 21303, 21304, 12404, 12404, 12397, 12397, 12308];
    [12557, 12566, 21278, 21280, 21316, 21308, 21307, 21305, 12406];
    [12453, 1652, 1640, 1639, 12437, 1563, 1565, 1557, 1552];
    [12496, 12488, 21176, 21173, 21170, 21171, 21164, 11916, 11903];
};

thermalL_indices = {
    [259, 260, 6186, 6193, 5696, 5696, 5686, 5685, 5564];
    [479, 478, 19161, 19163, 5692, 5676, 5675, 5673,994];
    [5025, 5028, 19122, 19121, 19110, 5638, 5637, 5626, 5224];
    [5078, 5598, 19171, 19169, 5612, 5602, 5604, 5603, 552];
    [434, 435, 448, 5606, 5605, 5608, 19000, 5304, 5296];
    [4683, 4710, 18923, 18904, 18911, 18912, 18912, 5352, 5352];
    [427, 5505, 5508, 18972, 5484, 5483, 5481, 525, 517];
    [400, 19002, 19008, 19064, 19062, 19058, 5780, 512, 510];
};
thermalR_indices = {
    [12766, 21433, 21438, 21428, 21380, 21377, 21373, 12032, 12028];
    [12783, 21425, 21430, 21431, 21386, 21385, 21381, 12038,12040];
    [1697, 1714, 1715, 21254, 21204, 21202, 21201, 12342, 12332];
    [12379, 12390, 21181, 21182, 21251, 21251, 21218, 12320, 12312];
    [12509, 1690, 21303, 21304, 12404, 12404, 12397, 12397, 12308];
    [12557, 12566, 21278, 21280, 21316, 21308, 21307, 21305, 12406];
    [12453, 1652, 1640, 1639, 12437, 1563, 1565, 1557, 1552];
    [12496, 12488, 21176, 21173, 21170, 21171, 21164, 11916, 11903];
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



%% === Plot Final Flux ===
% Custom colors and markers
colors = {[0 0 1], [1 0.5 0], [1 0 1], [0 0 0]};
markers = {'s', 'v', 'o', '+'};

figure;

% === Left Limiter Plot ===
subplot(2,1,1)
plot(z_vertL, flux_rfL{1}, 'LineStyle', '-', 'Marker', markers{1}, 'Color', colors{1}, 'DisplayName', 'O^{8+} (STRIPE)'); hold on;
plot(z_vertL, flux_rfL{2}, 'LineStyle', '-', 'Marker', markers{2}, 'Color', colors{2}, 'DisplayName', 'O^{7+} (STRIPE)');
plot(z_vertL, flux_rfL{3}, 'LineStyle', '-', 'Marker', markers{3}, 'Color', colors{3}, 'DisplayName', 'O^{6+} (STRIPE)');
plot(z_vertL, flux_rfL{4}, 'LineStyle', '-', 'Marker', markers{4}, 'Color', colors{4}, 'DisplayName', 'O^{5+} (STRIPE)');
plot(z_vertL, flux_total_L, 'g*-', 'DisplayName', 'I_{\phi, W-I} (STRIPE)');

xlabel('Vertical position along Q2 limiters: z [m]');
ylabel('Brightness I_\phi(W-I) [m^{-2}s^{-1}St^{-1}]');
title('(a) ICRH Case w/RF Sheath (Left Limiter)');
legend show;
grid on;

% === Right Limiter Plot ===
subplot(2,1,2)
plot(z_vertR, flux_rfR{1}, 'LineStyle', '-', 'Marker', markers{1}, 'Color', colors{1}, 'DisplayName', 'O^{8+} (STRIPE)'); hold on;
plot(z_vertR, flux_rfR{2}, 'LineStyle', '-', 'Marker', markers{2}, 'Color', colors{2}, 'DisplayName', 'O^{7+} (STRIPE)');
plot(z_vertR, flux_rfR{3}, 'LineStyle', '-', 'Marker', markers{3}, 'Color', colors{3}, 'DisplayName', 'O^{6+} (STRIPE)');
plot(z_vertR, flux_rfR{4}, 'LineStyle', '-', 'Marker', markers{4}, 'Color', colors{4}, 'DisplayName', 'O^{5+} (STRIPE)');
plot(z_vertR, flux_total_R, 'g*-', 'DisplayName', 'I_{\phi, W-I} (STRIPE)');

xlabel('Vertical position along Q2 limiters: z [m]');
ylabel('Brightness I_\phi(W-I) [m^{-2}s^{-1}St^{-1}]');
title('(b) ICRH Case w/RF Sheath (Right Limiter)');
legend show;
grid on;

%% === Plot n_i, T_i, Flux = n_i·V, and V (Multi-Species, Dual-Axis) with SI Units ===
species_list = {'o8', 'o7', 'o6', 'o5'};
dataRF_list = {dataRF_o8, dataRF_o7, dataRF_o6, dataRF_o5};

% Custom colors and markers: blue, orange, magenta, black
colors = {[0 0 1], [1 0.5 0], [1 0 1], [0 0 0]};
markers = {'s', 'v', 'o', '+'};

ni_L_all = zeros(8, length(species_list));
ni_R_all = zeros(8, length(species_list));
Ti_L_all = zeros(8, length(species_list));
Ti_R_all = zeros(8, length(species_list));
V_L_all  = zeros(8, length(species_list));
V_R_all  = zeros(8, length(species_list));
flux_L_all = zeros(8, length(species_list));
flux_R_all = zeros(8, length(species_list));

% Compute all species data
for i = 1:length(species_list)
    data = dataRF_list{i};

    for j = 1:8
        ne_L = mean(data(rfL_indices{j},2));  % [m^-3]
        ne_R = mean(data(rfR_indices{j},2));
        ni_L = mean(data(rfL_indices{j},11));  % [m^-3]
        ni_R = mean(data(rfR_indices{j},11));

        Te_L = mean(data(rfL_indices{j},3));  % [m^-3]
        Te_R = mean(data(rfR_indices{j},3));

        V_L  = max(data(rfL_indices{j},1));   % [V]
        V_R  = max(data(rfR_indices{j},1));

        flow_L  = mean(data(rfL_indices{j},5));   % [m/s]
        flow_R  = mean(data(rfR_indices{j},5));


        ni_L_all(j,i) = ni_L;
        ni_R_all(j,i) = ni_R;

        ne_L_all(j,i) = ne_L;
        ne_R_all(j,i) = ne_R;


  

        Ti_L_all(j,i) = mean(data(rfL_indices{j},4));  % [eV]
        Ti_R_all(j,i) = mean(data(rfR_indices{j},4));
        Te_L_all(j,i) = mean(data(rfL_indices{j},3));  % [eV]
        Te_R_all(j,i) = mean(data(rfR_indices{j},3));

        V_L_all(j,i) = V_L;
        V_R_all(j,i) = V_R;

        flux_L_all(j,i) = abs(flow_L);  % [m^-2·s^-1] if V is interpreted as velocity
        flux_R_all(j,i) = abs(flow_R);
    end
end


%% === Plot All Quantities: n_i/n_e, T_i/T_e, Flux = n_i·V, V ===
figure('Name', 'n_i, n_e, T_i, T_e, Flux, V (Dual-Axis, SI Units)', ...
       'Position', [100, 100, 1000, 800]);

titles = {'(a)', '(b)', '(c)', '(d)'};
subplot_titles = {
    'n_O [m^{-3}]', ...
    'T_O [eV]', ...
    '\Gamma_O m^{-2}·s^{-1}]', ...
    'V_{Sheath}^{RF} [V]'};

% ---- Ion Density ----
subplot(4,1,1); hold on;
yyaxis left;
for i = 1:length(species_list)
    plot(z_vertL, ni_L_all(:,i), [markers{i} '-'], 'Color', colors{i}, ...
         'DisplayName', ['O^{', species_list{i}(2), '+} (L)']);
end
ylabel([subplot_titles{1}, ' (L)']);

yyaxis right;
for i = 1:length(species_list)
    plot(z_vertR, ni_R_all(:,i), [markers{i} '--'], 'Color', colors{i}, ...
         'DisplayName', ['O^{', species_list{i}(2), '+} (R)']);
end
ylabel([subplot_titles{1}, ' (R)']);

% title(['RF Case: ', subplot_titles{1}]);
title([titles{1}]);
legend('Location', 'best', 'NumColumns', 2); grid on; box on; 

% ---- Ion Temperature ----
subplot(4,1,2); hold on;
yyaxis left;
for i = 1:length(species_list)
    plot(z_vertL, Ti_L_all(:,i), [markers{i} '-'], 'Color', colors{i}, ...
         'DisplayName', ['O^{', species_list{i}(2), '+} (L)']);
end
ylabel([subplot_titles{2}, ' (L)']);

yyaxis right;
for i = 1:length(species_list)
    plot(z_vertR, Ti_R_all(:,i), [markers{i} '--'], 'Color', colors{i}, ...
         'DisplayName', ['O^{', species_list{i}(2), '+} (R)']);
end
ylabel([subplot_titles{2}, ' (R)']);

% title(['RF Case: ', subplot_titles{2}]);
title([titles{2}]);
legend('Location', 'best', 'NumColumns', 2); grid on; box on; 

% ---- Flux = n_i · V ----
subplot(4,1,3); hold on;
yyaxis left;
for i = 1:length(species_list)
    plot(z_vertL, flux_L_all(:,i), [markers{i} '-'], 'Color', colors{i}, ...
         'DisplayName', ['O^{', species_list{i}(2), '+} (L)']);
end
ylabel([subplot_titles{3}, ' (L)']);

yyaxis right;
for i = 1:length(species_list)
    plot(z_vertR, flux_R_all(:,i), [markers{i} '--'], 'Color', colors{i}, ...
         'DisplayName', ['O^{', species_list{i}(2), '+} (R)']);
end
ylabel([subplot_titles{3}, ' (R)']);

% title(['RF Case: ', subplot_titles{3}]); 
title([titles{3}]);
legend('Location', 'best', 'NumColumns', 2); grid on; box on; 

% ---- Plasma Potential ----
%% === Load O8+ Data ===
v_data = readmatrix('Targets_comsol_o8_new.txt');
V_thermal_L = zeros(8,1); V_thermal_R = zeros(8,1);

for j = 1:8
    V_thermal_L(j) = 3.*max(v_data(thermalL_indices{j}, 3));   % Plasma potential [V]
    V_thermal_R(j) = 3.*max(v_data(thermalR_indices{j}, 3));
end
subplot(4,1,4); 
figure;
hold on;
yyaxis left;
for i = 1:1
     plot(z_vertL, V_L_all(:,i), [markers{i} '-'], 'Color', colors{i}, ...
         'DisplayName', 'RF Sheath (L)');
    plot(z_vertL, V_thermal_L, 'k-', 'LineWidth', 1.2, ...
     'DisplayName', 'Thermal Sheath (L)');
  
    


end
ylabel([subplot_titles{4}, ' (L)']);

yyaxis right;
for i = 1:1
    plot(z_vertR, V_R_all(:,i), [markers{i} '--'], 'Color', colors{i}, ...
         'DisplayName', 'RF Sheath (R)');
    plot(z_vertR, V_thermal_R, 'k--', 'LineWidth', 1.2, ...
     'DisplayName', 'Thermal Sheath (R)');
end
ylabel([subplot_titles{4}, ' (R)']);

xlabel('Vertical position along Q2 limiters: z [m]');
% title(['RF Case: ', subplot_titles{4}]);
title([titles{4}]);
legend('Location', 'best', 'NumColumns', 2); grid on; box on; 


%% === Plot Yields at RF Indices (Left and Right Limiters) ===

% Define yield data per species
yield_data_list = {
    yields_dataRF_o8(:,8), ...
    yields_dataRF_o7(:,7), ...
    yields_dataRF_o6(:,6), ...
    yields_dataRF_o5(:,5)
};

% Left Limiter
yields_rfL = cell(1, 4);
for i = 1:4
    data_shifted = [0; yield_data_list{i}];
    temp_yield = zeros(8, 1);
    for j = 1:8
        temp_yield(j) = mean(data_shifted(rfL_indices{j}));
    end
    yields_rfL{i} = temp_yield;
end

% Right Limiter
yields_rfR = cell(1, 4);
for i = 1:4
    data_shifted = [0; yield_data_list{i}];
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

%% === Plot Yields and Erosion Fluxes with Totals (Dual Axis: Left vs Right Limiters) ===

% Custom colors and markers
colors = {[0 0 1], [1 0.5 0], [1 0 1], [0 0 0]};
markers = {'s', 'v', 'o', '+'};

figure('Name', 'Yields and Erosion Fluxes (Left vs Right Limiters)', ...
       'Position', [100, 100, 1000, 800]);

%% === (a) Sputtering Yields ===
subplot(2,1,1); hold on;

yyaxis left
for i = 1:4
    plot(z_vertL, yields_rfL{i}, [markers{i} '-'], ...
        'Color', colors{i}, 'LineWidth', 1.5, ...
        'DisplayName', ['Left: O^{', species_list{i}(2), '+}']);
end
ylabel('Sputtering Yield (Left)');
ylim([0, max(cellfun(@(x) max(x), yields_rfL)) * 1.2]);

yyaxis right
for i = 1:4
    plot(z_vertL, yields_rfR{i}, [markers{i} '--'], ...
        'Color', colors{i}, 'LineWidth', 1.5, ...
        'DisplayName', ['Right: O^{', species_list{i}(2), '+}']);
end
ylabel('Sputtering Yield (Right)');
xlabel('z [m]');
title('(a) Sputtering Yields at RF Indices (Left vs Right Limiter)');
legend('Location', 'northeastoutside'); grid on; box on;

%% === (b) Erosion Fluxes (Raw, No S/XB) ===

% Compute total erosion flux for each limiter
total_flux_L = zeros(8,1);
total_flux_R = zeros(8,1);
for i = 1:4
    total_flux_L = total_flux_L + flux_eroded_raw_rfL{i};
    total_flux_R = total_flux_R + flux_eroded_raw_rfR{i};
end

subplot(2,1,2); hold on;

yyaxis left
for i = 1:4
    plot(z_vertL, flux_eroded_raw_rfL{i}, [markers{i} '-'], ...
        'Color', colors{i}, 'LineWidth', 1.5, ...
        'DisplayName', ['Left: O^{', species_list{i}(2), '+}']);
end
plot(z_vertL, total_flux_L, 'g-*', 'LineWidth', 1.8, ...
    'DisplayName', 'Left: Total I_\phi');
ylabel('Erosion Flux I_\phi [m^{-2}s^{-1}] (Left)');
ylim([0, max(total_flux_L)*1.2]);

yyaxis right
for i = 1:4
    plot(z_vertL, flux_eroded_raw_rfR{i}, [markers{i} '--'], ...
        'Color', colors{i}, 'LineWidth', 1.5, ...
        'DisplayName', ['Right: O^{', species_list{i}(2), '+}']);
end
plot(z_vertL, total_flux_R, 'g--*', 'LineWidth', 1.8, ...
    'DisplayName', 'Right: Total I_\phi');
ylabel('Erosion Flux I_\phi [m^{-2}s^{-1}] (Right)');
xlabel('z [m]');
title('(b) RF-Induced Erosion Fluxes (Left vs Right Limiter)');
legend('Location', 'northeastoutside'); grid on; box on;