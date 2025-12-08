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

% Plotting functions
plot_log_patch = @(data, title_text) patch(transpose(X), transpose(Y), ...
    transpose(Z), log10(abs(data(1:size(X,1)))), 'FaceAlpha', 1, 'EdgeAlpha', 0.3);
plot_linear_patch = @(data, title_text) patch(transpose(X), transpose(Y), ...
    transpose(Z), abs(data(1:size(X,1))), 'FaceAlpha', 1, 'EdgeAlpha', 0.3);

% Define species and files
species = {'D', 'C1+', 'C2+', 'C3+', 'C4+', 'C5+', 'C6+'};
target_files = {
    '../ieads_D/Targets_D.txt', ...
    '../ieads_c1+/Targets_c1.txt', ...
    '../ieads_c2+/Targets_c2.txt', ...
    '../ieads_c3+/Targets_c3.txt', ...
    '../ieads_c4+/Targets_c4.txt', ...
    '../ieads_c5+/Targets_c5.txt', ...
    '../ieads_c6+/Targets_c6.txt'
};
yield_files = {
    '../ieads_D/yields_D.csv', ...
    '../ieads_c1+/yields_c1.csv', ...
    '../ieads_c2+/yields_c2.csv', ...
    '../ieads_c3+/yields_c3.csv', ...
    '../ieads_c4+/yields_c4.csv', ...
    '../ieads_c5+/yields_c5.csv', ...
    '../ieads_c6+/yields_c6.csv'
};

% Yield columns per file
yield_column_indices = [1, 1, 2, 3, 4, 5, 6];

% === Plot each species individually ===
for i = 1:length(species)
    data = readmatrix(target_files{i});
    yields_data = readmatrix(yield_files{i});
    Y_eff = yields_data(:, yield_column_indices(i));

    % Plasma parameters
    ne_data = data(:,2);
    v_data = data(:,6);
    ni_data = data(:,4);  % Ion density
    sputtering_yield = Y_eff;
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

%% === Combined subplot layout ===
figure('Name', 'Sputtering Yields and Gross Erosion Fluxes for All Species');

for i = 1:length(species)
    data = readmatrix(target_files{i});
    yields_data = readmatrix(yield_files{i});
    Y_eff = yields_data(:, yield_column_indices(i));
    ne_data = data(:,2);
    v_data = data(:,6);
    ni_data = data(:,4);
    sputtering_yield = Y_eff;
    ero_data = sputtering_yield .* ni_data .* v_data;

    % Yield subplot
    subplot(2, length(species), i);
    plot_linear_patch(sputtering_yield, ['Sputtering Yield for ', species{i}]);
    title(['Sputtering Yield (', species{i}, ')'], 'FontSize', 12);
    colorbar('southoutside'); view(30, 30); axis equal tight;
    xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]'); set(gca, 'FontSize', 10);

    % Erosion subplot
    subplot(2, length(species), i + length(species));
    plot_log_patch(ero_data, ['Gross Erosion Flux for ', species{i}]);
    title(['Erosion Flux (', species{i}, ')'], 'FontSize', 12);
    colorbar('southoutside'); view(30, 30); axis equal tight;
    xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]'); set(gca, 'FontSize', 10);
end

sgtitle('Top: Sputtering Yields | Bottom: Gross Erosion Fluxes', 'FontSize', 16);



%% === Carbon species only ===
carbon_species = {'C1+', 'C2+', 'C3+', 'C4+', 'C5+', 'C6+'};
carbon_target_files = {
    '../ieads_c1+/Targets_c1.txt', ...
    '../ieads_c2+/Targets_c2.txt', ...
    '../ieads_c3+/Targets_c3.txt', ...
    '../ieads_c4+/Targets_c4.txt', ...
    '../ieads_c5+/Targets_c5.txt', ...
    '../ieads_c6+/Targets_c6.txt'
};
carbon_yield_files = {
    '../ieads_c1+/yields_c1.csv', ...
    '../ieads_c2+/yields_c2.csv', ...
    '../ieads_c3+/yields_c3.csv', ...
    '../ieads_c4+/yields_c4.csv', ...
    '../ieads_c5+/yields_c5.csv', ...
    '../ieads_c6+/yields_c6.csv'
};
carbon_yield_cols = [1, 2, 3, 4, 5, 6];  % Column index for each species

%% === Figure 1: Carbon Sputtering Yields ===
figure('Name', 'Carbon Sputtering Yields');
tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

for i = 1:length(carbon_species)
    data = readmatrix(carbon_target_files{i});
    yields_data = readmatrix(carbon_yield_files{i});
    col_idx = carbon_yield_cols(i);

    if size(yields_data, 2) < col_idx
        error('Yield file %s does not have column %d.', carbon_yield_files{i}, col_idx);
    end

    Y_eff = yields_data(:, col_idx);

    % Match geometry
    Y_eff = Y_eff(1:size(X,1));

    fprintf('Species: %s | Yield min: %.2e, max: %.2e\n', ...
        carbon_species{i}, min(Y_eff), max(Y_eff));

    subplot(3, 2, i);
    plot_linear_patch(Y_eff);
    title(['Y_{eff} (', carbon_species{i}, ')'], 'FontSize', 12);

    % if exist('coolwarm.mat', 'file')
    %     load('coolwarm.mat', 'coolwarm_rgb');
    %     colormap(coolwarm_rgb);
    % end

    colormap("hsv");colorbar; view(30, 30); axis tight equal;
    xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
    set(gca, 'FontSize', 10);
    clim([0 1.6]);  % Uncomment to enforce limits
end

sgtitle('Sputtering Yields for Carbon Species (C^{1+}–C^{6+})', 'FontSize', 16);


% Save erosion flux from last species for inspection
writematrix(ero_data, 'ero_data.txt');

%% === Updated plot_log_patch function ===
plot_log_patch = @(data) patch( ...
    'XData', transpose(X), ...
    'YData', transpose(Y), ...
    'ZData', transpose(Z), ...
    'CData', repmat((abs(data(:))), 1, 3)', ...
    'FaceAlpha', 1, ...
    'EdgeAlpha', 0.3, ...
    'FaceColor', 'interp');

%% === Figure: Carbon Gross Erosion Fluxes ===
figure('Name', 'Carbon Gross Erosion Fluxes');
tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

for i = 1:length(carbon_species)
    % Load species-specific data
    data = readmatrix(carbon_target_files{i});
    yields_data = readmatrix(carbon_yield_files{i});
    col_idx = carbon_yield_cols(i);

    % Safety check for column existence
    if size(yields_data, 2) < col_idx
        error('Yield file %s does not have column %d.', carbon_yield_files{i}, col_idx);
    end

    % Extract & trim to match geometry
    Y_eff = yields_data(:, col_idx);
    Y_eff = Y_eff(1:size(X,1));

    ni_data = data(:, 4); ni_data = ni_data(1:size(X,1));
    v_data  = data(:, 6); v_data  = v_data(1:size(X,1));

    % Compute erosion flux
    ero_data = Y_eff .* ni_data .* v_data;

    % Sanitize flux data to avoid log10(0) or negative issues
    ero_data(ero_data <= 0 | isnan(ero_data) | isinf(ero_data)) = 1e-30;
    ero_data = ero_data(1:size(X,1));

    % Debug info
    fprintf('Species: %s | ero_data min: %.2e, max: %.2e\n', ...
        carbon_species{i}, min(ero_data), max(ero_data));

    % Plot in tiled subplot
    subplot(3, 2, i);
    plot_log_patch(ero_data);
    title(['\Gamma_{gross} (', carbon_species{i}, ')'], 'FontSize', 12);

    if exist('coolwarm.mat', 'file')
        load('coolwarm.mat', 'coolwarm_rgb');
        colormap(coolwarm_rgb);
    end

    colorbar('southoutside');
    view(30, 30); axis equal;
    xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
    set(gca, 'FontSize', 10);
    
    % Optional: enforce clim once plots are verified
    clim([18 20]);
end

sgtitle('Gross Erosion Fluxes for Carbon Species (C^{1+}–C^{6+})', 'FontSize', 16);

%% === Figure: Total Yields and Gross Erosion Summary ===
figure('Name', 'Total Carbon Yields and Erosion Summary');
tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

% --- Compute total yield and total erosion flux across species ---
total_Y = zeros(size(X,1), 1);
total_ero = zeros(size(X,1), 1);

for i = 1:length(carbon_species)
    data = readmatrix(carbon_target_files{i});
    yields_data = readmatrix(carbon_yield_files{i});
    col_idx = carbon_yield_cols(i);

    % Trim to geometry
    Y_eff = yields_data(:, col_idx);
    Y_eff = Y_eff(1:size(X,1));

    ni_data = data(:, 4); ni_data = ni_data(1:size(X,1));
    v_data  = data(:, 6); v_data  = v_data(1:size(X,1));

    ero_data = Y_eff .* ni_data .* v_data;
    ero_data(ero_data <= 0 | isnan(ero_data) | isinf(ero_data)) = 1e-30;

    % Accumulate totals
    total_Y   = total_Y + Y_eff;
    total_ero = total_ero + ero_data;
end

% --- Subplot 1: Total Sputtering Yield ---
nexttile;
plot_linear_patch(total_Y);
title('Total Effective Yield (C^{1+}–C^{6+})', 'FontSize', 12);
colormap("hsv");
% colorbar('southoutside');
view(30, 30); axis equal tight;
xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
set(gca, 'FontSize', 10);
clim([0 1.6]);

% --- Subplot 2: Total Gross Erosion Flux ---
nexttile;
plot_log_patch(total_ero);
title('Total Gross Erosion Flux (C^{1+}–C^{6+})', 'FontSize', 12);
colormap("hsv");
% colorbar('southoutside');
view(30, 30); axis equal tight;
xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
set(gca, 'FontSize', 10);

sgtitle('Total Carbon Yields and Gross Erosion Flux', 'FontSize', 16);

% Optionally write totals to file
writematrix([total_Y total_ero], 'carbon_total_yields_erosion.txt');

%% === Updated plot_linear_patch function ===
plot_linear_patch = @(data) patch( ...
    'XData', transpose(X), ...
    'YData', transpose(Y), ...
    'ZData', transpose(Z), ...
    'CData', repmat(abs(data(:)), 1, 3)', ...
    'FaceAlpha', 1, ...
    'EdgeAlpha', 0.3, ...
    'FaceColor', 'interp');

% === Figure: Carbon Gross Erosion Fluxes (Linear Scale) ===
figure('Name', 'Carbon Gross Erosion Fluxes (Linear Scale)');
tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

for i = 1:length(carbon_species)
    % Load target and yield data
    data = readmatrix(carbon_target_files{i});
    yields_data = readmatrix(carbon_yield_files{i});
    col_idx = carbon_yield_cols(i);

    if size(yields_data, 2) < col_idx
        error('Yield file %s does not have column %d.', carbon_yield_files{i}, col_idx);
    end

    Y_eff = yields_data(:, col_idx);
    Y_eff = Y_eff(1:size(X,1));

    ni_data = data(:, 4); ni_data = ni_data(1:size(X,1));
    v_data  = data(:, 6); v_data  = v_data(1:size(X,1));

    % Erosion flux calculation
    ero_data = Y_eff .* ni_data .* v_data;
    ero_data(isnan(ero_data) | isinf(ero_data)) = 0;

    % Truncate to match mesh size
    ero_data = ero_data(1:size(X,1));

    % Debug print
    fprintf('Species: %s | Linear Flux min: %.2e, max: %.2e\n', ...
        carbon_species{i}, min(ero_data), max(ero_data));

    % Plot
    subplot(3, 2, i);
    plot_linear_patch(ero_data);
    title(['\Gamma_{gross} (', carbon_species{i}, ')'], 'FontSize', 12);

    % if exist('coolwarm.mat', 'file')
    %     load('coolwarm.mat', 'coolwarm_rgb');
    %     colormap(coolwarm_rgb);
    % end

    colorbar;
    colormap("hsv")
    view(30, 30); axis equal tight;
    xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
    set(gca, 'FontSize', 10);

    % Optional: enforce clim once verified (adjust based on your data)
    % clim([0 2e20]);
end

sgtitle('Gross Erosion Fluxes for Carbon Species (Linear Scale)', 'FontSize', 16);

%% === Update: 2D Patch Function in (phi, theta) ===

% Geometry triangles: X, Y, Z as N×3 matrices

R = sqrt(X.^2 + Y.^2);
tor_x0 = 0; %1.49345169307177;
tor_y0 = 2.21%0.011285617584263547;
phi = atan2(X, Y);  % N×3
theta = atan2(Z - tor_y0, R - tor_x0);  % N×3

%% === Plotting Functions in Cylindrical Coordinates (phi–theta plane) ===
plot_linear_patch_phi_theta = @(data) patch( ...
    'XData', transpose(phi), ...
    'YData', transpose(theta), ...
    'ZData', 0 .* transpose(theta), ...
    'CData', repmat(abs(data(:)), 1, 3)', ...
    'FaceAlpha', 1, ...
    'EdgeAlpha', 0.3, ...
    'FaceColor', 'interp');

%% === Figure 1: Carbon Sputtering Yields (φ–θ) ===
figure('Name', 'Carbon Yields in Cylindrical Coordinates');
tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

for i = 1:length(carbon_species)
    % Load data
    data = readmatrix(carbon_target_files{i});
    yields_data = readmatrix(carbon_yield_files{i});
    col_idx = carbon_yield_cols(i);

    if size(yields_data, 2) < col_idx
        error('Yield file %s does not have column %d.', carbon_yield_files{i}, col_idx);
    end

    % Extract and truncate
    Y_eff = yields_data(:, col_idx);
    Y_eff = Y_eff(1:size(X,1));

    subplot(3, 2, i);
    plot_linear_patch_phi_theta(Y_eff);
    title(['Y_{eff} (', carbon_species{i}, ')'], 'FontSize', 12);
    colormap("jet");colorbar;
    xlabel('\phi [rad]'); ylabel('\theta [rad]');
    view(2); %axis equal tight;

    % if exist('coolwarm.mat', 'file')
    %     load('coolwarm.mat', 'coolwarm_rgb');
    %     colormap(coolwarm_rgb);
    % end
    clim([0.9 1.6])
     xlim([-0.484 0.484]);
    ylim([-0.7032 -0.635])
    box on;
end

% sgtitle('Sputtering Yields in (φ, θ) Cylindrical Coordinates', 'FontSize', 16);

%% === Figure 2: Carbon Erosion Fluxes (φ–θ) ===
figure('Name', 'Carbon Erosion Fluxes in Cylindrical Coordinates');
tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

for i = 1:length(carbon_species)
    % Load data
    data = readmatrix(carbon_target_files{i});
    yields_data = readmatrix(carbon_yield_files{i});
    col_idx = carbon_yield_cols(i);

    if size(yields_data, 2) < col_idx
        error('Yield file %s does not have column %d.', carbon_yield_files{i}, col_idx);
    end

    Y_eff = yields_data(:, col_idx);
    Y_eff = Y_eff(1:size(X,1));
    ni_data = data(:, 4); ni_data = ni_data(1:size(X,1));
    v_data  = data(:, 6); v_data  = v_data(1:size(X,1));

    ero_data = Y_eff .* ni_data .* v_data;
    ero_data(isnan(ero_data) | isinf(ero_data)) = 0;
    ero_data = ero_data(1:size(X,1));

    subplot(3, 2, i);
    plot_linear_patch_phi_theta(ero_data);
    title(['\Gamma_{gross} (', carbon_species{i}, ')'], 'FontSize', 12);
    colorbar;
    colormap('jet')
     xlim([-0.484 0.484]);
    ylim([-0.7032 -0.635])
    xlabel('\phi [rad]'); ylabel('\theta [rad]');
    view(2); %axis equal tight;

    % if exist('coolwarm.mat', 'file')
    %     load('coolwarm.mat', 'coolwarm_rgb');
    %     colormap(coolwarm_rgb);
    % end
    clim([0 1e20])
    box on;
end

% sgtitle('Gross Erosion Fluxes in (φ, θ) Cylindrical Coordinates', 'FontSize', 16);

%% === Combine D⁺ and Carbon for Total Yields & Erosion Flux ===

% --- Load Deuterium data ---
data_D = readmatrix('../ieads_D/Targets_D.txt');
yields_D = readmatrix('../ieads_D/yields_D.csv');
Y_D = yields_D(:,1);
Y_D = Y_D(1:size(X,1));
ni_D = data_D(:,4); ni_D = ni_D(1:size(X,1));
v_D  = data_D(:,6); v_D  = v_D(1:size(X,1));
flux_D = Y_D .* ni_D .* v_D;
flux_D(isnan(flux_D) | isinf(flux_D)) = 0;

% --- Initialize carbon sums ---
Y_C_total = zeros(size(X,1), 1);
flux_C_total = zeros(size(X,1), 1);

% --- Sum over carbon species ---
for i = 1:length(carbon_species)
    data = readmatrix(carbon_target_files{i});
    yields_data = readmatrix(carbon_yield_files{i});
    col_idx = carbon_yield_cols(i);

    Y_eff = yields_data(:, col_idx);
    Y_eff = Y_eff(1:size(X,1));

    ni = data(:, 4); ni = ni(1:size(X,1));
    v  = data(:, 6); v  = v(1:size(X,1));

    flux = Y_eff .* ni .* v;
    flux(isnan(flux) | isinf(flux)) = 0;

    Y_C_total = Y_C_total + Y_eff;
    flux_C_total = flux_C_total + flux;
end

% === Plot: Combined Yields and Fluxes in Cylindrical Coordinates ===
figure('Name', 'Yields and Erosion Fluxes for D⁺ and C Species');
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% --- Subplot 1: D⁺ Yields ---
nexttile(1);
plot_linear_patch_phi_theta(Y_D);
title('Y_{eff} (D^+)', 'FontSize', 12);
colormap("jet"); colorbar;
xlabel('\phi [rad]'); ylabel('\theta [rad]'); view(2);
 xlim([-0.484 0.484]);
    ylim([-0.7032 -0.635])
     ylabel(colorbar, '[Particles/m^2/s]');
box on;

% --- Subplot 2: D⁺ Erosion Flux ---
nexttile(2);
plot_linear_patch_phi_theta(flux_D);
title('\Gamma_{gross} (D^+)', 'FontSize', 12);
colormap("jet"); colorbar;
xlabel('\phi [rad]'); ylabel('\theta [rad]'); view(2);
 xlim([-0.484 0.484]);
    ylim([-0.7032 -0.635])
     ylabel(colorbar, '[Particles/m^2/s]');
box on;

% --- Subplot 3: Carbon Yields Total ---
nexttile(3);
plot_linear_patch_phi_theta(Y_C_total);
title('Y_{eff} (C^{1+}–C^{6+})', 'FontSize', 12);
colormap("jet"); colorbar;
xlabel('\phi [rad]'); ylabel('\theta [rad]'); view(2);
 xlim([-0.484 0.484]);
    ylim([-0.7032 -0.635])
box on;
 ylabel(colorbar, '[Particles/m^2/s]');

% --- Subplot 4: Carbon Erosion Flux Total ---
nexttile(4);
plot_linear_patch_phi_theta(flux_C_total);
title('\Gamma_{gross} (C^{1+}–C^{6+})', 'FontSize', 12);
colormap("jet"); colorbar;
xlabel('\phi [rad]'); ylabel('\theta [rad]'); view(2);
% clim([0 7e18]); 
box on;
 xlim([-0.484 0.484]);
    ylim([-0.7032 -0.635])
   ylabel(colorbar, '[Particles/m^2/s]');

sgtitle('Case # 200882', 'FontSize', 16);

%%
 figure;
theta=57.296.*data(:,11);
    plot_linear_patch_phi_theta(theta);
    title(['\Gamma_{gross} (', carbon_species{i}, ')'], 'FontSize', 12);
    colorbar;
    colormap('jet')
     xlim([-0.484 0.484]);
    ylim([-0.7032 -0.635])
    xlabel('\phi [rad]'); ylabel('\theta [rad]');
    view(2); %axis equal tight;


  
%% === Compute Total Densities and Flows (D + Carbon) ===

% --- Deuterium (D⁺) ---
ni_D = data_D(:,4); ni_D = ni_D(1:size(X,1));
v_D  = data_D(:,6); v_D  = v_D(1:size(X,1));

% initialize carbon totals
ni_C_total = zeros(size(X,1),1);
nv_C_total = zeros(size(X,1),1); % for density-weighted velocity

for i = 1:length(carbon_species)
    data = readmatrix(carbon_target_files{i});

    ni = data(:,4); ni = ni(1:size(X,1));
    v  = data(:,6); v  = v(1:size(X,1));

    % accumulate
    ni_C_total = ni_C_total + ni;
    nv_C_total = nv_C_total + ni .* v;
end

% total density = D + all carbons
density_total = ni_D + ni_C_total;

% total flow velocity = density-weighted average
flow_total = v_D;

% === Compute Total Sputtering Yield (D⁺ + C¹⁺–C⁶⁺) ===

% Deuterium yield (already loaded earlier)
Y_D = yields_D(:,1);
Y_D = Y_D(1:size(X,1));

% Initialize and sum carbon yields
Y_C_total = zeros(size(X,1),1);
for i = 1:length(carbon_species)
    yields_data = readmatrix(carbon_yield_files{i});
    col_idx = carbon_yield_cols(i);
    Y_eff = yields_data(:, col_idx);
    Y_eff = Y_eff(1:size(X,1));
    Y_C_total = Y_C_total + Y_eff;
end

% Total yield
Y_total = Y_D + Y_C_total;

%% === Save to .txt file ===
writematrix(Y_total, 'yields_total_200882.txt');
writematrix(density_total, 'density_total_200882.txt');
writematrix(flow_total,    'flow_total_200882.txt');