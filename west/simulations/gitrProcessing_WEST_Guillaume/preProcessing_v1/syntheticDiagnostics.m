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
figure;
plot(z_vertL, flux_rfL{1}, 'ob-', 'DisplayName', 'o8'); hold on;
plot(z_vertL, flux_rfL{2}, 'sr-', 'DisplayName', 'o7');
plot(z_vertL, flux_rfL{3}, '^g-', 'DisplayName', 'o6');
plot(z_vertL, flux_rfL{4}, 'dk-', 'DisplayName', 'o5');
plot(z_vertL, flux_total_L, 'm*-', 'DisplayName', 'Total');

xlabel('z position (z_{vertL})');
ylabel('Flux (normalized by 4\pi·SXB)');
title('RF-induced Flux for All Species and Total');
legend show;
grid on;


%% === Plot Final Flux ===
figure;
plot(z_vertR, flux_rfR{1}, 'ob-', 'DisplayName', 'o8'); hold on;
plot(z_vertR, flux_rfR{2}, 'sr-', 'DisplayName', 'o7');
plot(z_vertR, flux_rfR{3}, '^g-', 'DisplayName', 'o6');
plot(z_vertR, flux_rfR{4}, 'dk-', 'DisplayName', 'o5');
plot(z_vertR, flux_total_R, 'm*-', 'DisplayName', 'Total');

xlabel('z position (z_{vertL})');
ylabel('Flux (normalized by 4\pi·SXB)');
title('RF-induced Flux for All Species and Total');
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

        V_L  = mean(data(rfL_indices{j},1));   % [V]
        V_R  = mean(data(rfR_indices{j},1));

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
figure('Name', 'n_i, n_e, T_i, T_e, Flux, V (Dual-Axis, SI Units)');

% ---- n_i / n_e [m^-3] ----
subplot(4,1,1); hold on;
for i = 1:length(species_list)
    yyaxis left;
    plot(z_vertL, ni_L_all(:,i), [markers{i} '-'], 'Color', colors{i}, ...
        'DisplayName', ['n_i O^{', species_list{i}(2), '+} (L)']);
    yyaxis right;
    plot(z_vertR, ni_R_all(:,i), [markers{i} '--'], 'Color', colors{i}, ...
        'DisplayName', ['n_i O^{', species_list{i}(2), '+} (R)']);
end
% % Electron density (species-independent)
% yyaxis left;
% plot(z_vertL, mean(ne_L_all,2), 'k-.', 'LineWidth', 1.5, 'DisplayName', 'n_e (L)');
% yyaxis right;
% plot(z_vertR, mean(ne_R_all,2), 'k--', 'LineWidth', 1.5, 'DisplayName', 'n_e (R)');

set(gca, 'YScale', 'linear');
ylabel('Density [m^{-3}]');
title('Ion and Electron Density (log scale)');
legend('Location', 'best', 'NumColumns', 2); grid on;

% ---- T_i / T_e [eV] ----
subplot(4,1,2); hold on;
for i = 1:length(species_list)
    yyaxis left;
    plot(z_vertL, Ti_L_all(:,i), [markers{i} '-'], 'Color', colors{i}, ...
        'DisplayName', ['T_i O^{', species_list{i}(2), '+} (L)']);
    yyaxis right;
    plot(z_vertR, Ti_R_all(:,i), [markers{i} '--'], 'Color', colors{i}, ...
        'DisplayName', ['T_i O^{', species_list{i}(2), '+} (R)']);
end
% Electron temperature
% yyaxis left;
% plot(z_vertL, mean(Te_L_all,2), 'k-.', 'LineWidth', 1.5, 'DisplayName', 'T_e (L)');
% yyaxis right;
% plot(z_vertR, mean(Te_R_all,2), 'k--', 'LineWidth', 1.5, 'DisplayName', 'T_e (R)');

set(gca, 'YScale', 'linear');
ylabel('Temperature [eV]');
title('Ion and Electron Temperature (log scale)');
legend('Location', 'best', 'NumColumns', 2); grid on;

% ---- Flux = n_i · V [m^-2·s^-1] ----
subplot(4,1,3); hold on;
for i = 1:length(species_list)
    yyaxis left;
    plot(z_vertL, flux_L_all(:,i), [markers{i} '-'], 'Color', colors{i}, ...
        'DisplayName', ['Flux O^{', species_list{i}(2), '+} (L)']);
    yyaxis right;
    plot(z_vertR, flux_R_all(:,i), [markers{i} '--'], 'Color', colors{i}, ...
        'DisplayName', ['Flux O^{', species_list{i}(2), '+} (R)']);
end
ylabel('Flux [m^{-2}·s^{-1}]');
title('Ion Particle Flux (n_i·V)'); legend('Location', 'best', 'NumColumns', 2); grid on;

% ---- Plasma Potential V [V] ----
subplot(4,1,4); hold on;
for i = 1:length(species_list)
    yyaxis left;
    plot(z_vertL, V_L_all(:,i), [markers{i} '-'], 'Color', colors{i}, ...
        'DisplayName', ['V O^{', species_list{i}(2), '+} (L)']);
    yyaxis right;
    plot(z_vertR, V_R_all(:,i), [markers{i} '--'], 'Color', colors{i}, ...
        'DisplayName', ['V O^{', species_list{i}(2), '+} (R)']);
end
xlabel('z position [m]');
ylabel('Plasma Potential [V]');
title('Plasma Potential'); legend('Location', 'best', 'NumColumns', 2); grid on;