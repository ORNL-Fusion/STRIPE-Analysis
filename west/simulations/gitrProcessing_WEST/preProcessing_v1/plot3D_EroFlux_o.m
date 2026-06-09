% Clear and set up
clear; clc; close all;

%% === Load Geometry ===
fid = fopen('gitrGeometryPointPlane3d_comsol.cfg');
for i = 1:20
    tline = fgetl(fid);
    if i > 2 && ischar(tline)
        evalc(tline);
    end
end
fclose(fid);

% Geometry
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

%% === Load Data for O6 to O8 Only ===
dataRF_o6 = readmatrix('Targets_comsol_o6.txt');
dataRF_o7 = readmatrix('Targets_comsol_o7.txt');
dataRF_o8 = readmatrix('Targets_comsol_o8.txt');

yields_dataRF_o6 = readmatrix('yields_comsol_o6.csv');
yields_dataRF_o7 = readmatrix('yields_comsol_o7.csv');
yields_dataRF_o8 = readmatrix('yields_comsol_o8.csv');

% Case labels
cases = {'o6', 'o7', 'o8'};

% Assign data into cell arrays
data_all = {dataRF_o6, dataRF_o7, dataRF_o8};
yields_all = {yields_dataRF_o6, yields_dataRF_o7, yields_dataRF_o8};

% Output arrays to store for plotting
yield_all = cell(1, 3);
erosion_all = cell(1, 3);

% Loop to compute yield and erosion flux per case
for i = 1:3
    data = data_all{i};
    yields = yields_all{i};

    ni = data(:,11);             % Ion density
    v = abs(data(:,5));          % Ion velocity magnitude

    % Use last column of yield data
    yield_col = size(yields, 2);
    yield_raw = yields(:, yield_col);

    % Match to geometry face count
    num_faces = size(phi, 1);
    yield_full = zeros(num_faces, 1);
    if length(yield_raw) == num_faces
        yield_full = yield_raw;
    elseif length(yield_raw) == num_faces - 2
        yield_full(3:end) = yield_raw;
    else
        warning('Yield length mismatch in %s', cases{i});
        yield_full(end - length(yield_raw) + 1:end) = yield_raw;
    end

    % Compute erosion flux
    erosion_flux = yield_full .* ni .* v;

    % Store for plotting
    yield_all{i} = yield_full;
    erosion_all{i} = erosion_flux;
end

%% === Plotting: 3 rows × 2 columns ===
% figure('Position', [100, 100, 1600, 800]);
% t = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% Latex-style panel labels
panel_labels = {'(a)', '(b)', '(c)', '(d)', '(e)', '(f)'};
cases = {'O^{6+}', 'O^{7+}', 'O^{8+}'};

% === Top row: Sputtering Yield ===
for i = 1:3
    nexttile(i);
    plot_linear_patch(yield_all{4 - i});  % Reverse order: O8 → O6
    colormap("parula");
    colorbar('eastoutside');
    clim([0 0.5]);
    axis equal;
    xlim([-0.2 0.2]); ylim([-0.4 0.4]);
   title(['Y_{eff} (' cases{4 - i} ') [ICRH]'], 'FontSize', 14, 'Interpreter', 'tex');

    xlabel('\phi [rad]'); ylabel('\theta [rad]');
    set(gca, 'FontSize', 12);
    text(-0.18, 0.35, panel_labels{i}, 'FontSize', 14, 'FontWeight', 'bold');
end

% === Bottom row: Gross Erosion Flux ===
for i = 1:3
    nexttile(i + 3);
    plot_linear_patch(erosion_all{4 - i});
    colormap("parula");
    colorbar('eastoutside');
    clim([0 5e18]);
    axis equal;
    xlim([-0.2 0.2]); ylim([-0.4 0.4]);
    title(['\Gamma_{gross, W} (' cases{4 - i} ') [ICRH]'], 'FontSize', 14, 'Interpreter', 'tex');
    xlabel('\phi [rad]'); ylabel('\theta [rad]');
    set(gca, 'FontSize', 12);
    text(-0.18, 0.35, panel_labels{i+3}, 'FontSize', 14, 'FontWeight', 'bold');
end

% Overall figure title
% title(t, 'Sputtering Yield and Gross Erosion Flux (O8 to O6)', 'FontSize', 20);