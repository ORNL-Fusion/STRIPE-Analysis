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

% Define subset of data points and build X,Y,Z from geometry
subset = 1:length(x1);
X = [transpose(x1(subset)), transpose(x2(subset)), transpose(x3(subset))];
Y = [transpose(y1(subset)), transpose(y2(subset)), transpose(y3(subset))];
Z = [transpose(z1(subset)), transpose(z2(subset)), transpose(z3(subset))];

% Ensure X, Y, Z are defined (after building them)
if ~exist('X', 'var') || ~exist('Y', 'var') || ~exist('Z', 'var')
    error('Geometry variables X, Y, and Z are not defined. Check gitrGeometryPointPlane3d.cfg.');
end

% Plotting functions
plot_log_patch = @(data, title_text) patch(transpose(X), transpose(Y), ...
    transpose(Z), log10(abs(data(1:size(X,1)))), 'FaceAlpha', 1, 'EdgeAlpha', 0.3);
plot_linear_patch = @(data, title_text) patch(transpose(X), transpose(Y), ...
    transpose(Z), abs(data(1:size(X,1))), 'FaceAlpha', 1, 'EdgeAlpha', 0.3);

% Define species and files (Neon + D+)
species = {'D^+', 'Ne1^+', 'Ne2^+', 'Ne3^+', 'Ne4^+', 'Ne5^+', 'Ne6^+', 'Ne7^+', 'Ne8^+', 'Ne9^+', 'Ne10^+'};
target_files = {
    '../ieads_D+/Targets_Ne1+.txt', ...         % FIX: correct D+ targets path
    '../ieads_Ne1+/Targets_Ne1+.txt', ...
    '../ieads_Ne2+/Targets_Ne2+.txt', ...
    '../ieads_Ne3+/Targets_Ne3+.txt', ...
    '../ieads_Ne4+/Targets_Ne4+.txt', ...
    '../ieads_Ne5+/Targets_Ne5+.txt', ...
    '../ieads_Ne6+/Targets_Ne6+.txt', ...
    '../ieads_Ne7+/Targets_Ne7+.txt', ...
    '../ieads_Ne8+/Targets_Ne8+.txt', ...
    '../ieads_Ne9+/Targets_Ne9+.txt', ...
    '../ieads_Ne10+/Targets_Ne10+.txt'
};
yield_files = {
    '../ieads_D+/yields_D+.csv', ...
    '../ieads_Ne1+/yields_Ne1+.csv', ...
    '../ieads_Ne2+/yields_Ne2+.csv', ...
    '../ieads_Ne3+/yields_Ne3+.csv', ...
    '../ieads_Ne4+/yields_Ne4+.csv', ...
    '../ieads_Ne5+/yields_Ne5+.csv', ...
    '../ieads_Ne6+/yields_Ne6+.csv', ...
    '../ieads_Ne7+/yields_Ne7+.csv', ...
    '../ieads_Ne8+/yields_Ne8+.csv', ...
    '../ieads_Ne9+/yields_Ne9+.csv', ...
    '../ieads_Ne10+/yields_Ne10+.csv'
};

% Yield columns per file (Neon/D+: use column 1 everywhere)
yield_column_indices = ones(1, numel(yield_files));

% === Plot each species individually ===
for i = 1:length(species)
    data = readmatrix(target_files{i});
    yields_data = readmatrix(yield_files{i});
    % yields_data_D = readmatrix(yield_files{1});
    Y_eff = yields_data(:, yield_column_indices(i));   % FIX: use col 1, no padding

    % Plasma parameters (per file definitions: col2=ne, col4=ni, col6=v)
    ne_data = data(:,2);
    v_data  = data(:,5);                                % FIX: use column 6 for flow
    ni_data = data(:,11);                                % FIX: use column 4 for ion density

    sputtering_yield = [0;Y_eff];
    ero_data = sputtering_yield .* ni_data .* v_data;

    % Plot per species
    figure('Name', ['Species: ', species{i}]);

    % Subplot 1: Sputtering Yield
    subplot(2,1,1);
    plot_linear_patch(sputtering_yield, ['Sputtering Yield for ', species{i}]);
    title(['Sputtering Yield for ', species{i}], 'FontSize', 16);
    colorbar('eastoutside');
    xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
    set(gca, 'FontSize', 10); axis equal tight; view(30, 30);

    % Subplot 2: Gross Erosion Flux
    subplot(2,1,2);
    plot_log_patch(ero_data, ['Gross Erosion Flux for ', species{i}]);
    title(['Gross Erosion Flux for ', species{i}], 'FontSize', 16);
    colorbar('eastoutside');
    xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
    set(gca, 'FontSize', 10); axis equal tight; view(30, 30);
end

