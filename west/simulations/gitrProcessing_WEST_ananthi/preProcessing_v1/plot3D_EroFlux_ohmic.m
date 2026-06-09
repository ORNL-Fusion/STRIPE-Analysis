% Clear and set up
clear; clc; close all;

%% === Load Geometry ===
fid = fopen('gitrGeometryPointPlane3d_thermal.cfg');
for i = 1:20
    tline = fgetl(fid);
    if i > 2 && ischar(tline)
        evalc(tline);
    end
end
fclose(fid);

% Define triangle vertex arrays
subset = 1:length(x1);
X = [transpose(x1(subset)), transpose(x2(subset)), transpose(x3(subset))];
Y = [transpose(y1(subset)), transpose(y2(subset)), transpose(y3(subset))];
Z = [transpose(z1(subset)), transpose(z2(subset)), transpose(z3(subset))];

% Cylindrical coordinates
R = sqrt(X.^2 + Y.^2);
tor_x0 = 1.49345169307177;
tor_y0 = 0.011285617584263547;
phi = atan2(X, Y);
theta = atan2(Z - tor_y0, R - tor_x0);

% Load colormap
load('coolwarm.mat', 'coolwarm_rgb');

% Plotting function
plot_linear_patch = @(data) patch(transpose(phi), transpose(theta), ...
    0 .* transpose(theta), data, 'FaceAlpha', 1, 'EdgeAlpha', 0.3);
%% === Load Ohmic Data for O6 to O8 ===
data_o6 = readmatrix('Targets_thermal_o6_ohmic.txt');
data_o7 = readmatrix('Targets_thermal_o7_ohmic.txt');
data_o8 = readmatrix('Targets_thermal_o8_ohmic.txt');

yields_o6 = readmatrix('yields_thermal_o6_ohmic.csv');
yields_o7 = readmatrix('yields_thermal_o7_ohmic.csv');
yields_o8 = readmatrix('yields_thermal_o8_ohmic.csv');


% Case labels (used for plot titles)
cases = {'O^{6+}', 'O^{7+}', 'O^{8+}'};

% Organize data into cell arrays
data_all = {data_o6, data_o7, data_o8};
yields_all = {yields_o6, yields_o7, yields_o8};

% Storage for processed yield and erosion flux
yield_all = cell(1, 3);
erosion_all = cell(1, 3);

% Number of triangles
nFaces = length(x1);

%% === Compute yields and erosion fluxes ===
for i = 1:3
    data = data_all{i};
    yields = yields_all{i};

    ni = data(:,11);
    v = abs(data(:,5));
    yield_raw = yields(:, end);  % Last column = sputtering yield

    yield_full = zeros(nFaces, 1);

    if length(yield_raw) == nFaces
        yield_full = yield_raw;
    elseif length(yield_raw) == nFaces - 2
        yield_full(3:end) = yield_raw;
        fprintf('Adjusted yield length for %s by padding 2 zeros.\n', cases{i});
    elseif length(yield_raw) < nFaces
        yield_full(end - length(yield_raw) + 1:end) = yield_raw;
        warning('Padded yield for %s: expected %d, got %d\n', cases{i}, nFaces, length(yield_raw));
    else
        yield_full = yield_raw(1:nFaces);
        warning('Truncated yield for %s: expected %d, got %d\n', cases{i}, nFaces, length(yield_raw));
    end

    erosion_flux = yield_full .* ni .* v;

    yield_all{i} = yield_full;
    erosion_all{i} = erosion_flux;
end

%% === Plotting: 2x3 layout (O8 → O6, left to right) ===
figure('Position', [100, 100, 1600, 800]);
t = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

panel_labels = {'(a)', '(b)', '(c)', '(d)', '(e)', '(f)'};

% === Top row: Sputtering Yield ===
for i = 1:3
    nexttile(i);  % Tiles 1, 2, 3
    plot_linear_patch(yield_all{4 - i});  % Reverse order: O8 → O6
    colormap("parula");
    colorbar('eastoutside');
    clim([0 0.5]);
    axis equal;
    % xlim([1.3 1.6]); ylim([-0.4 0.4]);
    title(['Y_{eff} (' cases{4 - i} ') [Ohmic]'], 'FontSize', 14, 'Interpreter', 'tex');
    xlabel('R [m]'); ylabel('Z [m]');
    set(gca, 'FontSize', 12);
    text(1.31, 0.35, panel_labels{i}, 'FontSize', 14, 'FontWeight', 'bold');
end

% === Bottom row: Gross Erosion Flux ===
for i = 1:3
    nexttile(i + 3);  % Tiles 4, 5, 6
    plot_linear_patch(erosion_all{4 - i});
    colormap("parula");
    colorbar('eastoutside');
    clim([0 5e14]);
    axis equal;
    % xlim([1.3 1.6]); ylim([-0.4 0.4]);
    title(['\Gamma_{gross, W} (' cases{4 - i} ') [Ohmic]'], 'FontSize', 14, 'Interpreter', 'tex');
    xlabel('R [m]'); ylabel('Z [m]');
    set(gca, 'FontSize', 12);
    text(1.31, 0.35, panel_labels{i+3}, 'FontSize', 14, 'FontWeight', 'bold');
end

% Overall title
title(t, 'Sputtering Yield and Gross Erosion Flux (Ohmic, O8 to O6)', 'FontSize', 20);