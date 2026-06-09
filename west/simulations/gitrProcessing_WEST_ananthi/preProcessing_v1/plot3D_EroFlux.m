% Clear workspace and close all figures
clear; clc; close all;

% Load data
% load('diiid_geom_5cmFar.mat')
data = readmatrix('Targets_comsol_o7.txt');
yields_data = readmatrix('yields_comsol_o7.csv');

% Extract relevant data
potential_data = data(:,1);
ne_data = data(:,2);
ni_data = data(:,11);
te_data = data(:,3);
v_data = abs(data(:,5));
flux_data = [0;0;yields_data(:,7)].*ni_data .* v_data;
ero_data =  flux_data;
sputtering_yield = yields_data(:,7);

% % Load geometry if not already in workspace
% if ~exist('x1', 'var')
%     fid = fopen('gitrGeometryPointPlane3d_comsol.cfg');
%     if fid == -1
%         error('Failed to open gitrGeometryPointPlane3d.cfg');
%     end
%     for i = 1:20
%         tline = fgetl(fid);
%         if i > 2 && ischar(tline)
%             evalc(tline);
%         end
%     end
%     fclose(fid);
% end
% 
% % Ensure X, Y, Z are defined
% if ~exist('X', 'var') || ~exist('Y', 'var') || ~exist('Z', 'var')
%     error('Geometry variables X, Y, and Z are not defined. Check gitrGeometryPointPlane3d.cfg.');
% end

if (exist('x1') == 0)
fid = fopen(strcat('gitrGeometryPointPlane3d_comsol.cfg'));
% fid = fopen(strcat('gitrGeometryPointPlane3d_thermal.cfg'));
tline = fgetl(fid);
tline = fgetl(fid);
for i=1:18
tline = fgetl(fid);
evalc(tline);
end
Zsurface = Z;
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

plot_linear_patch = @(data, title_text) patch(transpose(phi), transpose(theta), ...
    0 .* transpose(theta), ((data)), 'FaceAlpha', 1, 'EdgeAlpha', 0.3);


% Create figure with four subplots

% % Subplot 1: Potential (Log Scale)
% subplot(4,1,1);
% figure
% % plot_log_patch(potential_data, 'Rectified DC Voltage (Log Scale)');
% % title('Rectified DC Voltage (Log Scale)', 'FontSize', 18);
% % colorbar('eastoutside');
% % xlabel('Phi [rad]', 'FontSize', 18);
% % ylabel('Theta [rad]', 'FontSize', 18);
% % zlabel('Z [m]', 'FontSize', 18);
% % set(gca, 'FontSize', 18);
% 
% % Subplot 2: Sputtering Yield
% subplot(2,1,1);
% % figure;
% plot_log_patch([0;0;yields_data(:,7)], 'Gross Erosion Flux (D→C) (Log Scale)');
% title(' Flux (Log Scale)', 'FontSize', 18);
% colorbar('eastoutside');
% xlabel('Phi [rad]', 'FontSize', 18);
% ylabel('Theta [rad]', 'FontSize', 18);
% zlabel('Z [m]', 'FontSize', 18);
% set(gca, 'FontSize', 18);
% axis equal;
% xlim([-0.2 0.2])
% ylim([-0.4 0.4])
% % clim([-2 0])
% 
% 
% % Subplot 3: Gross Erosion Flux (Log Scale)
% subplot(2,1,2);
% % figure;
% plot_log_patch(ero_data, 'Gross Erosion Flux (D→C) (Log Scale)');
% title('n_e (Log Scale)', 'FontSize', 18);
% colorbar('eastoutside');
% xlabel('Phi [rad]', 'FontSize', 18);
% ylabel('Theta [rad]', 'FontSize', 18);
% zlabel('Z [m]', 'FontSize', 18);
% set(gca, 'FontSize', 18);
% axis equal;
% xlim([-0.2 0.2])
% ylim([-0.4 0.4])
% % clim([14 18])

figure
% plot_log_patch(potential_data, 'Rectified DC Voltage (Log Scale)');
% title('Rectified DC Voltage (Log Scale)', 'FontSize', 18);
% colorbar('eastoutside');
% xlabel('Phi [rad]', 'FontSize', 18);
% ylabel('Theta [rad]', 'FontSize', 18);
% zlabel('Z [m]', 'FontSize', 18);
% set(gca, 'FontSize', 18);

% Subplot 2: Sputtering Yield
subplot(2,1,1);
% figure;
plot_linear_patch([0;0;yields_data(:,7)], 'Gross Erosion Flux (D→C) (Log Scale)');
title(' Flux (Log Scale)', 'FontSize', 18);
colorbar('eastoutside');
xlabel('Phi [rad]', 'FontSize', 18);
ylabel('Theta [rad]', 'FontSize', 18);
zlabel('Z [m]', 'FontSize', 18);
set(gca, 'FontSize', 18);
axis equal;
xlim([-0.2 0.2])
ylim([-0.4 0.4])
% clim([-2 0])


% Subplot 3: Gross Erosion Flux (Log Scale)
subplot(2,1,2);
% figure;
plot_linear_patch(ero_data, 'Gross Erosion Flux (D→C) (Log Scale)');
title('n_e (Log Scale)', 'FontSize', 18);
colorbar('eastoutside');
xlabel('Phi [rad]', 'FontSize', 18);
ylabel('Theta [rad]', 'FontSize', 18);
zlabel('Z [m]', 'FontSize', 18);
set(gca, 'FontSize', 18);
axis equal;
xlim([-0.2 0.2])
ylim([-0.4 0.4])
clim([1e14 5e18])
% Load the coolwarm RGB values
% Load the custom RGB colormap
load('coolwarm.mat', 'coolwarm_rgb');  % Ensure coolwarm_rgb is Nx3 matrix in [0,1] range
colormap(coolwarm_rgb);

% Add colorbar after colormap to reflect it
colorbar('eastoutside');

% % Subplot 3: Gross Erosion Flux/D-flux
% subplot(4,1,4);
% % figure;
% patch(transpose(phi), transpose(theta), 0 .* transpose(theta), ero_data, ...
%     'FaceAlpha', 1, 'EdgeAlpha', 0.3);
% title('Gross Erosion Flux/D-flux', 'FontSize', 18);
% colorbar('eastoutside');
% xlabel('Phi [rad]', 'FontSize', 18);
% ylabel('Theta [rad]', 'FontSize', 18);
% zlabel('Z [m]', 'FontSize', 18);
% set(gca, 'FontSize', 18);
% 

% Overall title for the figure
sgtitle('Plasma-Wall Interaction on DIII-D Helicon Antenna', 'FontSize', 20);


% figure; 
% patch(transpose(X), transpose(Y), transpose(Z), v_data, ...
%     'FaceAlpha', 1, 'EdgeAlpha', 0.3);
% title('Flow on surface', 'FontSize', 18);
% colorbar('eastoutside');
% xlabel('X', 'FontSize', 18);
% ylabel('Y', 'FontSize', 18);
% zlabel('Z', 'FontSize', 18);
% set(gca, 'FontSize', 18);
% 
% figure; 
% patch(transpose(X), transpose(Y), transpose(Z),flux_data, ...
%     'FaceAlpha', 1, 'EdgeAlpha', 0.3);
% title('Gross Erosion Flux/D-flux', 'FontSize', 18);
% colorbar('eastoutside');
% xlabel('X', 'FontSize', 18);
% ylabel('Y', 'FontSize', 18);
% zlabel('Z', 'FontSize', 18);
% set(gca, 'FontSize', 18);