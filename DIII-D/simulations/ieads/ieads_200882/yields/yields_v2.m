%% Plot 3D sputtering yields for D+ and total Cx+

% Clear workspace and close all figures
% clear; clc; close all;

%% === Load geometry if not already in workspace ===
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

% Ensure geometry variables exist
if ~exist('x1', 'var') || ~exist('x2', 'var') || ~exist('x3', 'var') || ...
   ~exist('y1', 'var') || ~exist('y2', 'var') || ~exist('y3', 'var') || ...
   ~exist('z1', 'var') || ~exist('z2', 'var') || ~exist('z3', 'var')
    error('Geometry variables are not fully defined in gitrGeometryPointPlane3d.cfg.');
end

%% === Build triangle mesh ===
subset = 1:length(x1);

X = [transpose(x1(subset)), transpose(x2(subset)), transpose(x3(subset))];
Y = [transpose(y1(subset)), transpose(y2(subset)), transpose(y3(subset))];
Z = [transpose(z1(subset)), transpose(z2(subset)), transpose(z3(subset))];

nFaces = size(X,1);

%% === Patch plotting function for 3D yield ===
plot_linear_patch_3d = @(data) patch( ...
    'XData', transpose(X), ...
    'YData', transpose(Y), ...
    'ZData', transpose(Z), ...
    'CData', repmat(abs(data(:)), 1, 3)', ...
    'FaceColor', 'interp', ...
    'FaceAlpha', 1.0, ...
    'EdgeAlpha', 0.2);

%% === File definitions ===

% D+ files
target_file_D = '../ieads_D/Targets_D.txt';
yield_file_D  = '../ieads_D/yields_D.csv';

% Carbon charge-state files
carbon_species = {'C1+', 'C2+', 'C3+', 'C4+', 'C5+', 'C6+'};
carbon_yield_files = {
    '../ieads_c1+/yields_c1.csv', ...
    '../ieads_c2+/yields_c2.csv', ...
    '../ieads_c3+/yields_c3.csv', ...
    '../ieads_c4+/yields_c4.csv', ...
    '../ieads_c5+/yields_c5.csv', ...
    '../ieads_c6+/yields_c6.csv'
};

% Column indices in each carbon yield file
carbon_yield_cols = [1, 2, 3, 4, 5, 6];

%% === Load D+ yield ===
yields_D = readmatrix(yield_file_D);

if size(yields_D,2) < 1
    error('D yield file does not contain column 1: %s', yield_file_D);
end

Y_D = yields_D(:,1);
Y_D = Y_D(1:nFaces);

fprintf('D+ yield min = %.3e, max = %.3e\n', min(Y_D), max(Y_D));

%% === Compute total Cx+ yield ===
Y_C_total = zeros(nFaces,1);

for i = 1:length(carbon_species)
    yields_C = readmatrix(carbon_yield_files{i});
    col_idx = carbon_yield_cols(i);

    if size(yields_C,2) < col_idx
        error('Yield file %s does not have column %d.', carbon_yield_files{i}, col_idx);
    end

    Y_eff = yields_C(:,col_idx);
    Y_eff = Y_eff(1:nFaces);

    Y_C_total = Y_C_total + Y_eff;

    fprintf('%s yield min = %.3e, max = %.3e\n', ...
        carbon_species{i}, min(Y_eff), max(Y_eff));
end

fprintf('Total Cx+ yield min = %.3e, max = %.3e\n', ...
    min(Y_C_total), max(Y_C_total));

%% === Plot 3D yields: D+ and total Cx+ ===
figure('Name', '3D Yields for D+ and Total Cx+', 'Color', 'w');
tiledlayout(1,2, 'TileSpacing', 'compact', 'Padding', 'compact');

% --- D+ yield ---
nexttile;
plot_linear_patch_3d(Y_D);
title('Y_{eff} (D^+)', 'FontSize', 14);
xlabel('X [m]');
ylabel('Y [m]');
zlabel('Z [m]');
axis equal tight;
view(30,30);
colormap('jet');
cb1 = colorbar;
ylabel(cb1, 'Yield');
set(gca, 'FontSize', 11);
box on;

% --- Total Cx+ yield ---
nexttile;
plot_linear_patch_3d(Y_C_total);
title('Y_{eff} Total (C^{1+}–C^{6+})', 'FontSize', 14);
xlabel('X [m]');
ylabel('Y [m]');
zlabel('Z [m]');
axis equal tight;
view(30,30);
colormap('jet');
cb2 = colorbar;
ylabel(cb2, 'Yield');
set(gca, 'FontSize', 11);
box on;

sgtitle('3D Sputtering Yields: D^+ and Total Cx^+', 'FontSize', 16);

%% === Optional: save totals ===
writematrix(Y_D,       'yield_D_3d.txt');
writematrix(Y_C_total, 'yield_total_Cx_3d.txt');