
% Clear workspace and close all figures
% clear; clc; close all;

% Load geometry if not already in workspace
if ~exist('x1', 'var')
    fid = fopen('gitrGeometryPointPlane3d.cfg');
    if fid == -1
        error('Failed to open gitrGeometryPointPlane3d.cfg');
    end
    for i = 1:20
        tline = fgetl(fid);
        if i > 2 && ischar(tline)
            evalc(tline);
        end
    end
    fclose(fid);
end

% Ensure X, Y, Z are defined
if ~exist('X', 'var') || ~exist('Y', 'var') || ~exist('Z', 'var')
    error('Geometry variables X, Y, and Z are not defined. Check gitrGeometryPointPlane3d.cfg.');
end

% Define subset of data points
subset = 1:length(x1);
X = [transpose(x1(subset)), transpose(x2(subset)), transpose(x3(subset))];
Y = [transpose(y1(subset)), transpose(y2(subset)), transpose(y3(subset))];
Z = [transpose(z1(subset)), transpose(z2(subset)), transpose(z3(subset))];

% Convert to cylindrical coordinates
R = sqrt(X.^2 + Y.^2);
tor_x0 = 1.49345169307177;
tor_y0 = 0.011285617584263547;
phi = atan2(X, Y);
theta = atan2(Z - tor_y0, R - tor_x0);

% Function to plot with log-scale color using patch
plot_log_patch = @(data, title_text) patch(transpose(X), transpose(Y), ...
    transpose(Z), log10(abs(data)), 'FaceAlpha', 1, 'EdgeAlpha', 0.3);

plot_linear_patch = @(data, title_text) patch(transpose(X), transpose(Y), ...
    transpose(Z), (abs(data)), 'FaceAlpha', 1, 'EdgeAlpha', 0.3);


% Define species, yield files, target files, ni fractions
species = {'b1+', 'b2+', 'b3+', 'c2+', 'c3+', 'c4+'};
target_files = {'Targets_boron.txt', 'Targets_boron.txt', 'Targets_boron.txt', ...
                'Targets_carbon.txt', 'Targets_carbon.txt', 'Targets_carbon.txt'};
yield_files = {'yields_b1+.csv', 'yields_b2+.csv', 'yields_b3+.csv', ...
               'yields_c2+.csv', 'yields_c3+.csv', 'yields_c4+.csv'};
yield_column_indices = [1, 2, 3, 2, 3, 4];
ni_fractions = [0.12, 0.12, 0.16, 0.06, 0.06, 0.08];

% Loop over each species and create plots
for i = 1:length(species)
    data = readmatrix(target_files{i});
    yields_data = readmatrix(yield_files{i});

    % Extract plasma parameters
    potential_data = data(:,1);
    ne_data = data(:,2);
    te_data = data(:,3);
    v_data = data(:,5);
    flux_data = ne_data .* v_data;

    % Apply yield and ion density fraction
    Y_eff = [0; 0; yields_data(:, yield_column_indices(i))];
    ni_data = ni_fractions(i) * ne_data;

    % Calculate sputtering yield and erosion flux
    sputtering_yield = Y_eff;
    ero_data = sputtering_yield .* ni_data .* v_data;

    % Create figure
    figure('Name', ['Species: ', species{i}]);

    % Subplot 1: Sputtering Yield
    subplot(2,1,1);
    plot_linear_patch(sputtering_yield, ['Sputtering Yield for ', species{i}]);
    title(['Sputtering Yield for ', species{i}], 'FontSize', 16);
    colorbar('eastoutside');
   xlabel('X [m]', 'FontSize', 10);
    ylabel('Y [m]', 'FontSize', 10);
    zlabel('Z [m]', 'FontSize', 10);
    set(gca, 'FontSize', 10);
    axis equal tight;
    view(30, 30);

    % Subplot 2: Gross Erosion Flux (Log Scale)
    subplot(2,1,2);
    plot_log_patch(ero_data, ['Gross Erosion Flux for ', species{i}]);
    title(['Gross Erosion Flux for ', species{i}], 'FontSize', 16);
    colorbar('eastoutside');
    xlabel('X [m]', 'FontSize', 10);
    ylabel('Y [m]', 'FontSize', 10);
    zlabel('Z [m]', 'FontSize', 10);
    set(gca, 'FontSize', 10);
      axis equal tight;
       view(30, 30);
    
end

% Create a single figure for all subplots
figure('Name', 'Sputtering Yields and Gross Erosion Fluxes for All Species');

for i = 1:length(species)
    data = readmatrix(target_files{i});
    yields_data = readmatrix(yield_files{i});

    % Plasma params
    potential_data = data(:,1);
    ne_data = data(:,2);
    te_data = data(:,3);
    v_data = data(:,5);
    ni_data = ni_fractions(i) * ne_data;

    % Yield and erosion
    Y_eff = [0; 0; yields_data(:, yield_column_indices(i))];
    sputtering_yield = Y_eff;
    ero_data = sputtering_yield .* ni_data .* v_data;

    % Subplot for sputtering yield
    subplot(2, length(species), i);
    plot_linear_patch(sputtering_yield, ['Sputtering Yield for ', species{i}]);
    title(['Sputtering Yield for ', species{i}], 'FontSize', 16);
   colorbar('southoutside');
    view(30, 30);
    axis equal tight;
    xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
    set(gca, 'FontSize', 10);

    % Subplot for erosion flux
    subplot(2, length(species), i + length(species));
    plot_log_patch(ero_data, ['Gross Erosion Flux for ', species{i}]);
    title(['Gross Erosion Flux for ', species{i}], 'FontSize', 16);
    colorbar('southoutside');
    view(30, 30);
    axis equal tight;
    xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
    set(gca, 'FontSize', 10);
end

sgtitle('Top: Sputtering Yields | Bottom: Gross Erosion Fluxes', 'FontSize', 16);

% Create a tiled layout for compact display
figure('Name', 'Sputtering Yields and Gross Erosion Fluxes for All Species');
tiledlayout(2, 6, 'TileSpacing', 'compact', 'Padding', 'compact');

for i = 1:length(species)
    data = readmatrix(target_files{i});
    yields_data = readmatrix(yield_files{i});

    % Plasma params
    potential_data = data(:,1);
    ne_data = data(:,2);
    te_data = data(:,3);
    v_data = data(:,5);
    ni_data = ni_fractions(i) * ne_data;

    % Yield and erosion
    Y_eff = [0; 0; yields_data(:, yield_column_indices(i))];
    sputtering_yield = Y_eff;
    ero_data = sputtering_yield .* ni_data .* v_data;

    % Top row: sputtering yield
    nexttile(i);
     plot_linear_patch(sputtering_yield, ['Sputtering Yield for ', species{i}]);
    title(['Y_{eff} (', species{i}, ')'], 'FontSize', 12);
    load('coolwarm.mat', 'coolwarm_rgb')
    colormap(coolwarm_rgb)
    colorbar('southoutside');
    view(30, 30); axis equal;
    xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
    set(gca, 'FontSize', 10);clim([0 0.35])

    % Bottom row: erosion flux
    nexttile(i + length(species));
     plot_log_patch(ero_data, ['Gross Erosion Flux for ', species{i}]);
    title(['\Gamma_{gross, W} (', species{i}, ')'], 'FontSize', 12);
    load('coolwarm.mat', 'coolwarm_rgb')
    colormap(coolwarm_rgb)
    colorbar('southoutside');
    view(30, 30); axis equal;
    xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
    set(gca, 'FontSize', 10);clim([16 19])
end

sgtitle('Top: Sputtering Yields | Bottom: Gross Erosion Fluxes', 'FontSize', 16);

writematrix(ero_data, 'ero_data.txt');