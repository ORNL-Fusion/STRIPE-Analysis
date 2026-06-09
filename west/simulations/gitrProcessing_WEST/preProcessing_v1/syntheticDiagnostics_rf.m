close all;
clear;
clc;

%% === Load Data ===
dataRF_o8 = readmatrix('Targets_comsol_o8.txt');
dataRF_o7 = readmatrix('Targets_comsol_o7.txt');
dataRF_o6 = readmatrix('Targets_comsol_o6.txt');
dataRF_o5 = readmatrix('Targets_comsol_o5.txt');

yields_dataRF_o8 = readmatrix('yields_comsol_o8.csv');
yields_dataRF_o7 = readmatrix('yields_comsol_o7.csv');
yields_dataRF_o6 = readmatrix('yields_comsol_o6.csv');
yields_dataRF_o5 = readmatrix('yields_comsol_o5.csv');

%% === Extract and Compute ===

ni_dataRF_o8 = dataRF_o8(:,11);
ni_dataRF_o7 = dataRF_o7(:,11);
ni_dataRF_o6 = dataRF_o6(:,11);
ni_dataRF_o5 = dataRF_o5(:,11);

v_dataRF_o8 = abs(dataRF_o8(:,5));
v_dataRF_o7 = abs(dataRF_o7(:,5));
v_dataRF_o6 = abs(dataRF_o6(:,5));
v_dataRF_o5 = abs(dataRF_o5(:,5));

ne_dataRF_total = dataRF_o5(:,2)+ dataRF_o6(:,2)+ dataRF_o7(:,2)+dataRF_o8(:,2);
v_dataRF_total = abs(dataRF_o5(:,5)+dataRF_o6(:,5)+ dataRF_o7(:,5)+dataRF_o8(:,5));




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

ero_dataRF_o8 = ni_dataRF_o8 .* v_dataRF_o8 .* [0;0;yields_dataRF_o8(:,8)];
ero_dataRF_o7 = ni_dataRF_o7 .* v_dataRF_o7 .* [0;0;yields_dataRF_o7(:,7)];
ero_dataRF_o6 = 0.1.*ni_dataRF_o6 .* v_dataRF_o6 .* [0;0;yields_dataRF_o7(:,7)];
ero_dataRF_o5 = ni_dataRF_o5 .* v_dataRF_o5 .* [0;0;yields_dataRF_o5(:,5)];

ero_dataRF_total = ero_dataRF_o8 + ero_dataRF_o7 + ero_dataRF_o6 + ero_dataRF_o5;
yields_dataRF_total = [0;0;yields_dataRF_o8(:,8)]+[0;0;yields_dataRF_o8(:,7)]+[0;0;yields_dataRF_o8(:,6)]+[0;0;yields_dataRF_o8(:,5)];
writematrix(ero_dataRF_total,'ero_dataRF_total.txt');
writematrix(yields_dataRF_total,'yields_dataRF_total.txt');
writematrix(ne_dataRF_total,'density_dataRF_total.txt');
writematrix(v_dataRF_total,'flow_dataRF_total.txt');

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
z_vertL = [0.4; 0.3473; 0.1689; 0.0749; 0.0106; -0.0851; -0.3042; -0.4];
z_vertR = [0.4; 0.3713; 0.1828; 0.0839; 0.0034; -0.0913; -0.2782; -0.4];
sxbValueL = [9.58779;10.8139;12.1153;12.8246;13.4241;13.242;11.1808;9.70449];

rfL_indices = {
    [3755, 3773, 3777, 3778, 3778, 3798, 3797, 3797, 3790];
    [479, 509, 509, 511, 3812, 3804, 3801, 3801, 438];
    [3963, 3963, 12089, 12092, 12092, 3600, 3600, 3600, 2419];
    [3984, 12106, 12117, 12119, 12120, 3627, 3625, 3625, 2473];
    [4035, 12122, 12123, 12124, 12124, 3632, 3632, 3632, 2483];
    [812, 816, 12151, 12151, 12151, 3659, 3657, 3657, 2537];
    [4100, 4123, 12067, 12067, 12062, 12061, 3491, 3489, 3466];
    [4134, 4134, 12046, 12044, 12043, 12043, 3523, 3522, 3513];
};

rfR_indices = {
    [7714, 7732, 12760, 12758, 12757, 12757, 8248, 8246, 8245];
    [7727, 7728, 12766, 12756, 12755, 12754, 8256, 8232, 8231];
    [8541, 8542, 12799, 12795, 12793, 8174, 8173, 8173, 10855];
    [8556, 12820, 12818, 12811, 12810, 8144, 8119, 8118, 10909];
    [8570, 8581, 8584, 8583, 8111, 8110, 8109, 8086, 10941];
    [1068, 8598, 8600, 8599, 8592, 8604, 8602, 8607, 10967];
    [7909, 7911, 7912, 12736, 12736, 8007, 8006, 1559, 1581];
    [1323, 1521, 1525, 1526, 7971, 7972, 7970, 7969, 7950];
};
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

%% === Plot Combined RF Sheath Flux: Left vs Right Limiter (Log Y) ===
figure('Name','Combined RF Sheath Flux (Left vs Right Limiters, Log Scale)');
colors = {[0 0 1], [1 0.5 0], [1 0 1], [0 0 0]};
markers = {'s', 'v', 'o', '+'};

yyaxis left;
for i = 1:4
    plot(z_vertL, flux_rfL{i}, 'LineStyle', '-', 'Marker', markers{i}, ...
         'Color', colors{i}, 'DisplayName', ['O^{', num2str(9-i), '+} (L)']); hold on;
end
plot(z_vertL, flux_total_L, 'g*-', 'DisplayName', 'I_{\phi, W-I} (L)');
ylabel('Brightness I_{\phi, W-I} Left [m^{-2}s^{-1}St^{-1}]');
set(gca, 'YScale', 'log');

yyaxis right;
for i = 1:4
    plot(z_vertR, flux_rfR{i}, 'LineStyle', '--', 'Marker', markers{i}, ...
         'Color', colors{i}, 'DisplayName', ['O^{', num2str(9-i), '+} (R)']); hold on;
end
plot(z_vertR, flux_total_R, 'g--*', 'DisplayName', 'I_{\phi, W-I} (R)');
ylabel('Brightness I_{\phi, W-I} Right [m^{-2}s^{-1}St^{-1}]');
set(gca, 'YScale', 'log');

xlabel('Vertical position along Q2 limiters: z [m]');
title('Combined ICRH RF Sheath: Left vs Right Limiter Flux (Log Scale)');
legend('Location', 'bestoutside');
grid on; box on;

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

        flux_L_all(j,i) = ni_L * flow_L;  % [m^-2·s^-1] if V is interpreted as velocity
        flux_R_all(j,i) = ni_R * flow_R;
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
subplot(3,1,1); hold on;
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
subplot(3,1,2); hold on;
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
subplot(3,1,3); hold on;
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
v_data = readmatrix('Targets_thermal_o8.txt');

V_thermal_L = zeros(8,1); V_thermal_R = zeros(8,1);
 
for j = 1:8
    V_thermal_L(j) = max(v_data(thermalL_indices{j}, 1));   % Plasma potential [V]
    V_thermal_R(j) = max(v_data(thermalR_indices{j}, 1));
end

figure; 
% subplot(1,1,4); 
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

% Custom colors and markers
colors = {[0 0 1], [1 0.5 0], [1 0 1], [0 0 0]};
markers = {'s', 'v', 'o', '+'};

figure('Name', 'Yields and Erosion Fluxes (Left vs Right Limiters)', ...
       'Position', [100, 100, 1000, 800]);

%% === (a) Sputtering Yields ===
% Custom colors and markers
colors = {[0 0 1], [1 0.5 0], [1 0 1], [0 0 0]};
markers = {'s', 'v', 'o', '+'};

figure('Name', 'Yields and Erosion Fluxes (Left vs Right Limiters)', ...
       'Position', [100, 100, 1000, 800]);

%% === (a) Sputtering Yields ===
% Custom colors and markers
colors = {[0 0 1], [1 0.5 0], [1 0 1], [0 0 0]};
markers = {'s', 'v', 'o', '+'};

figure('Name', 'Yields and Erosion Fluxes (Left vs Right Limiters)', ...
       'Position', [100, 100, 1000, 800]);

%% === (a) Sputtering Yields ===
subplot(2,1,1); hold on;

yyaxis left
for i = 1:4
    Z = species_list{i}(2);  % e.g., '8'
    plot(z_vertL, yields_rfL{i}, [markers{i} '-'], ...
        'Color', colors{i}, 'LineWidth', 1.5, ...
        'DisplayName', ['$\rm O^{', Z, '+}$ (L)']);
end
ylabel('Y_{eff} (Left)');
ylim([0, 0.3]);

yyaxis right
for i = 1:4
    Z = species_list{i}(2);
    plot(z_vertL, yields_rfR{i}, [markers{i} '--'], ...
        'Color', colors{i}, 'LineWidth', 1.5, ...
        'DisplayName', ['$\rm O^{', Z, '+}$ (R)']);
end
ylabel('Y_{eff} (Right)');
xlabel('z [m]');
title('(a)');
legend('Location', 'northeastoutside', 'NumColumns', 2, 'Interpreter', 'latex');
grid on; box on;

%% === (b) Erosion Fluxes (Raw, No S/XB) - Log Scale ===
total_flux_L = zeros(8,1);
total_flux_R = zeros(8,1);
for i = 1:4
    total_flux_L = total_flux_L + flux_eroded_raw_rfL{i};
    total_flux_R = total_flux_R + flux_eroded_raw_rfR{i};
end

subplot(2,1,2); hold on;

yyaxis left
for i = 1:4
    Z = species_list{i}(2);
    plot(z_vertL, flux_eroded_raw_rfL{i}, [markers{i} '-'], ...
        'Color', colors{i}, 'LineWidth', 1.5, ...
        'DisplayName', ['$\rm O^{', Z, '+}$ (L)']);
end
plot(z_vertL, total_flux_L, 'g-*', 'LineWidth', 1.8, ...
    'DisplayName', 'Total (L)');
ylabel('\Gamma_{gross, W} [m^{-2}s^{-1}] (Left)');
set(gca, 'YScale', 'log');
ylim([1e13, 1e20]);

yyaxis right
for i = 1:4
    Z = species_list{i}(2);
    plot(z_vertL, flux_eroded_raw_rfR{i}, [markers{i} '--'], ...
        'Color', colors{i}, 'LineWidth', 1.5, ...
        'DisplayName', ['$\rm O^{', Z, '+}$ (R)']);
end
plot(z_vertL, total_flux_R, 'g--*', 'LineWidth', 1.8, ...
    'DisplayName', 'Total (R)');
ylabel('\Gamma_{gross, W} [m^{-2}s^{-1}] (Right)');
set(gca, 'YScale', 'log');
ylim([1e13, 1e20]);

xlabel('z [m]');
title('(b)');
legend('Location', 'northeastoutside', 'NumColumns', 2, 'Interpreter', 'latex');
grid on; box on;