% Clear workspace and close all figures
% clear; clc; close all;

% Load data
% load('diiid_geom_5cmFar.mat')
data = readmatrix('Targets.txt');
yields_Ne1 = readmatrix('yields_Ne+.csv');
yields_Ne2 = readmatrix('yields_Ne2+.csv');
yields_Ne3 = readmatrix('yields_Ne3+.csv');
yields_Ne4 = readmatrix('yields_Ne3+.csv');
yields_Ne5 = readmatrix('yields_Ne5+.csv');
yields_Ne6 = readmatrix('yields_Ne6+.csv');
yields_Ne7 = readmatrix('yields_Ne7+.csv');

yields=yields_Ne1(:,1)+ yields_Ne2(:,2)+yields_Ne3(:,3)+yields_Ne4(:,3)...
    +yields_Ne5(:,5)+yields_Ne6(:,6)+yields_Ne7(:,7);

% Extract relevant data
potential_data = data(:,1);
ne_data = data(:,2);
ne_data_ele= data(:,2);
ni_data= data(:,11);
te_data = data(:,3);
v_data = data(:,5);
flux_data =abs([0;yields].*ni_data .* v_data);
ero_data =  flux_data;
test = data(:,10);
% sputtering_yield = yields_data(:,2);

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
angle = atan2(Z - tor_y0, R - tor_x0);

% Function to plot with log-scale color using patch
plot_log_patch = @(data, title_text) patch(transpose(phi), transpose(Z), ...
    0 .* transpose(Z), log10(abs(data)), 'FaceAlpha', 1, 'EdgeAlpha', 0.1);

% Function to plot with linear-scale color using patch
plot_linear_patch = @(data, title_text) patch(transpose(phi), transpose(Z), ...
    0 .* transpose(Z), (data), 'FaceAlpha', 1, 'EdgeAlpha', 0.1);

% Create figure with three subplots
figure;

%% Subplot 1: Rectified DC Voltage
subplot(4,1,1);
plot_log_patch(potential_data, 'Rectified DC Voltage (Log Scale)');
% title('Rectified DC Voltage (Log Scale)', 'FontSize', 18);
% xlabel('Phi [rad]', 'FontSize', 18);
ylabel('$\theta$ [rad]','Interpreter', 'latex', 'FontSize', 18);
set(gca, 'FontSize', 18);
xlim([1.42 1.72]); ylim([-0.42 1.7]);
load('coolwarm.mat', 'coolwarm_rgb');
colormap(coolwarm_rgb);
clim([0 4]); % log10(V)
cb1 = colorbar('eastoutside');
cb1.Ticks = 0:1:4;
cb1.TickLabels = arrayfun(@(x) sprintf('10^{%d}', x), cb1.Ticks, 'UniformOutput', false);
ylabel(cb1, '$\rm V_{\mathrm{sheath}}^{\mathrm{RF}}$ [V]', 'Interpreter', 'latex');box on; 

%% Subplot 2: Sputtering Yield
subplot(4,1,2);
% plot_log_patch([0; yields], 'Gross Erosion Flux (D→C) (Log Scale)');
plot_linear_patch(test, 'Gross Erosion Flux (D→C) (Log Scale)');
% title('Gross Erosion Flux (Log Scale)', 'FontSize', 18);
% xlabel('Phi [rad]', 'FontSize', 18);
ylabel('$\theta$ [rad]','Interpreter', 'latex', 'FontSize', 18);
set(gca, 'FontSize', 18);
xlim([1.42 1.72]); ylim([-0.42 1.7]);
colormap(coolwarm_rgb);

clim([-1.5 1.5])

cb2 = colorbar('eastoutside');
% cb2.Ticks = log_min:1:log_max;
% cb2.TickLabels = arrayfun(@(x) sprintf('10^{%d}', x), cb2.Ticks, 'UniformOutput', false);
ylabel(cb2, '$\rm Y_{\mathrm{eff}} (Ne^{x+} \rightarrow W)$', 'Interpreter', 'latex');box on; 
%% Subplot 3: Sputtering Yield
subplot(4,1,3);
plot_log_patch([0; yields], 'Gross Erosion Flux (D→C) (Log Scale)');

% title('Gross Erosion Flux (Log Scale)', 'FontSize', 18);
% xlabel('Phi [rad]', 'FontSize', 18);
ylabel('$\theta$ [rad]','Interpreter', 'latex', 'FontSize', 18);
set(gca, 'FontSize', 18);
xlim([1.42 1.72]); ylim([-0.42 1.7]);
colormap(coolwarm_rgb);
colorbar
% Yield range: 1e-4 to 5e0 → log10 range: -4 to log10(5) ≈ 0.7
log_min = 0;
log_max = ceil(log10(5));  % = 1
clim([log_min log_max]);

cb2 = colorbar('eastoutside');
cb2.Ticks = log_min:1:log_max;
cb2.TickLabels = arrayfun(@(x) sprintf('10^{%d}', x), cb2.Ticks, 'UniformOutput', false);
ylabel(cb2, '$\rm Y_{\mathrm{eff}} (Ne^{x+} \rightarrow W)$', 'Interpreter', 'latex');box on; 
%% Subplot 3: Gross Erosion Flux (Log Scale)
subplot(4,1,4);
plot_log_patch(flux_data, 'n_e (Log Scale)');
% title('n_e (Log Scale)', 'FontSize', 18);
xlabel('$\phi$ [rad]','Interpreter', 'latex', 'FontSize', 18);
ylabel('$\theta$ [rad]','Interpreter', 'latex', 'FontSize', 18);
set(gca, 'FontSize', 18);
xlim([1.42 1.72]); ylim([-0.42 1.7]);
colormap(coolwarm_rgb);
clim([14 18]); % log10(n_e)
cb3 = colorbar('eastoutside');
cb3.Ticks = 14:1:18;
cb3.TickLabels = arrayfun(@(x) sprintf('10^{%d}', x), cb3.Ticks, 'UniformOutput', false);
ylabel(cb3, '$\rm \Gamma_{gross, W}$ [m$^{-3}$]', 'Interpreter', 'latex');box on; 
% Subplot 3: Gross Erosion Flux/D-flux
% subplot(4,1,4);
% figure;
% patch(transpose(phi), transpose(theta), 0 .* transpose(theta), ero_data./flux_data, ...
%     'FaceAlpha', 1, 'EdgeAlpha', 0.3);
% title('Gross Erosion Flux/D-flux', 'FontSize', 18);
% colorbar('eastoutside');
% xlabel('Phi [rad]', 'FontSize', 18);
% ylabel('Theta [rad]', 'FontSize', 18);
% zlabel('Z [m]', 'FontSize', 18);
% set(gca, 'FontSize', 18);


% Overall title for the figure
sgtitle('Plasma-Wall Interaction on DIII-D Helicon Antenna', 'FontSize', 20);


figure; 
patch(transpose(X), transpose(Y), transpose(Z), ni_data, ...
    'FaceAlpha', 1, 'EdgeAlpha', 0.3);
title('Neon Density on surface', 'FontSize', 18);
colorbar('eastoutside');
xlabel('X', 'FontSize', 18);
ylabel('Y', 'FontSize', 18);
zlabel('Z', 'FontSize', 18);
set(gca, 'FontSize', 18);

figure; 
patch(transpose(X), transpose(Y), transpose(Z),test, ...
    'FaceAlpha', 1, 'EdgeAlpha', 0.3);
title('Gross Erosion Flux/D-flux', 'FontSize', 18);
colorbar('eastoutside');
xlabel('X', 'FontSize', 18);
ylabel('Y', 'FontSize', 18);
zlabel('Z', 'FontSize', 18);
set(gca, 'FontSize', 18);


figure;

%% Subplot 1: Rectified DC Voltage
subplot(1,2,1);
plot_log_patch(potential_data, 'Rectified DC Voltage (Log Scale)');
title('$\rm V_{\mathrm{sheath}}^{\mathrm{RF}}$ [V]', 'Interpreter', 'latex');box on; 
xlabel('$\phi$ [rad]','Interpreter', 'latex', 'FontSize', 18);
ylabel('$\theta$ [rad]','Interpreter', 'latex', 'FontSize', 18);
set(gca, 'FontSize', 18);
xlim([1.42 1.72]); ylim([-0.42 1.7]);
load('coolwarm.mat', 'coolwarm_rgb');
colormap(coolwarm_rgb);
clim([0 4]); % log10(V)
cb1 = colorbar('eastoutside');
cb1.Ticks = 0:1:4;
cb1.TickLabels = arrayfun(@(x) sprintf('10^{%d}', x), cb1.Ticks, 'UniformOutput', false);
% ylabel(cb1, '$\rm V_{\mathrm{sheath}}^{\mathrm{RF}}$ [V]', 'Interpreter', 'latex');box on; 

subplot(1,2,2);
plot_log_patch(3.*te_data, 'Rectified DC Voltage (Log Scale)');
title('$\rm V_{\mathrm{sheath}}^{\mathrm{Thermal}}$ [V]', 'Interpreter', 'latex');box on; 
xlabel('$\phi$ [rad]','Interpreter', 'latex', 'FontSize', 18);
ylabel('$\theta$ [rad]','Interpreter', 'latex', 'FontSize', 18);
set(gca, 'FontSize', 18);
xlim([1.42 1.72]); ylim([-0.42 1.7]);
load('coolwarm.mat', 'coolwarm_rgb');
colormap(coolwarm_rgb);
clim([0 2]); % log10(V)
cb1 = colorbar('eastoutside');
cb1.Ticks = 0:1:4;
cb1.TickLabels = arrayfun(@(x) sprintf('10^{%d}', x), cb1.Ticks, 'UniformOutput', false);
% ylabel(cb1, '$\rm V_{\mathrm{sheath}}^{\mathrm{Thermal}}$ [V]', 'Interpreter', 'latex');box on; 