% Clear workspace and close all figures
% clear; clc; close all;

% Load data
% load('diiid_geom_5cmFar.mat')
data = readmatrix('Targets.txt');
yields_data = readmatrix('yields_newRef_360kw.csv');

% Extract relevant data
potential_data = data(:,1);
ne_data = data(:,2);
te_data = data(:,3);
v_data = data(:,5);
flux_data = ne_data .* v_data;
ero_data = yields_data(:,1) .* ne_data .* v_data;
sputtering_yield = yields_data(:,1);

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
plot_log_patch = @(data, title_text) patch(transpose(phi), transpose(theta), ...
    0 .* transpose(theta), log10(abs(data)), 'FaceAlpha', 1, 'EdgeAlpha', 0.3);

% Create figure with three subplots
figure;

% Subplot 1: Potential (Log Scale)
subplot(4,1,1);
plot_log_patch(potential_data, 'Rectified DC Voltage (Log Scale)');
title('Rectified DC Voltage (Log Scale)', 'FontSize', 18);
colorbar('eastoutside');
xlabel('Phi [rad]', 'FontSize', 18);
ylabel('Theta [rad]', 'FontSize', 18);
zlabel('Z [m]', 'FontSize', 18);
set(gca, 'FontSize', 18);

% Subplot 2: Sputtering Yield
subplot(4,1,2);
patch(transpose(phi), transpose(theta), 0 .* transpose(theta), sputtering_yield, ...
    'FaceAlpha', 1, 'EdgeAlpha', 0.3);
title('Effective Sputtering Yield (D→C)', 'FontSize', 18);
colorbar('eastoutside');
xlabel('Phi [rad]', 'FontSize', 18);
ylabel('Theta [rad]', 'FontSize', 18);
zlabel('Z [m]', 'FontSize', 18);
set(gca, 'FontSize', 18);


% Subplot 3: Gross Erosion Flux (Log Scale)
subplot(4,1,3);
plot_log_patch(ero_data, 'Gross Erosion Flux (D→C) (Log Scale)');
title('Gross Erosion Flux (D→C) (Log Scale)', 'FontSize', 18);
colorbar('eastoutside');
xlabel('Phi [rad]', 'FontSize', 18);
ylabel('Theta [rad]', 'FontSize', 18);
zlabel('Z [m]', 'FontSize', 18);
set(gca, 'FontSize', 18);

% Subplot 3: Gross Erosion Flux/D-flux
subplot(4,1,4);
patch(transpose(phi), transpose(theta), 0 .* transpose(theta), ero_data./flux_data, ...
    'FaceAlpha', 1, 'EdgeAlpha', 0.3);
title('Gross Erosion Flux/D-flux', 'FontSize', 18);
colorbar('eastoutside');
xlabel('Phi [rad]', 'FontSize', 18);
ylabel('Theta [rad]', 'FontSize', 18);
zlabel('Z [m]', 'FontSize', 18);
set(gca, 'FontSize', 18);


% Overall title for the figure
sgtitle('Plasma-Wall Interaction on DIII-D Helicon Antenna', 'FontSize', 20);

% Subplot 3: Gross Erosion Flux/D-flux
figure; 
patch(transpose(phi), transpose(theta), 0 .* transpose(theta), potential_data, ...
    'FaceAlpha', 1, 'EdgeAlpha', 0.3);
title('Gross Erosion Flux/D-flux', 'FontSize', 18);
colorbar('eastoutside');
xlabel('Phi [rad]', 'FontSize', 18);
ylabel('Theta [rad]', 'FontSize', 18);
zlabel('Z [m]', 'FontSize', 18);
set(gca, 'FontSize', 18);


figure;
subplot(2,1,1)
plot_log_patch(potential_data, 'Rectified DC Voltage (Log Scale)');
title('Rectified DC Voltage (Log Scale)', 'FontSize', 18);
colorbar('eastoutside');
xlabel('Phi [rad]', 'FontSize', 18);
ylabel('Theta [rad]', 'FontSize', 18);
zlabel('Z [m]', 'FontSize', 18);
set(gca, 'FontSize', 18);

% Subplot 2: Sputtering Yield
subplot(2,1,2);
patch(transpose(phi), transpose(theta), 0 .* transpose(theta), sputtering_yield, ...
    'FaceAlpha', 1, 'EdgeAlpha', 0.3);
title('Effective Sputtering Yield (D→C)', 'FontSize', 18);
colorbar('eastoutside');
xlabel('Phi [rad]', 'FontSize', 18);
ylabel('Theta [rad]', 'FontSize', 18);
zlabel('Z [m]', 'FontSize', 18);
set(gca, 'FontSize', 18);