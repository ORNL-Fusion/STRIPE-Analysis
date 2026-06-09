clc;
clear;
close all;
format long;

%% =========================
% 0) Settings
%% =========================
dt = 1e-9;                % [s]
nP = 5e6;                 % test particles used to build histogram
source_strength = 1.165799141267602e+16 ;   % [particles/s]

%%  Skimmer Info
skimmer_total_area=0.0080365;
skimmer_flux_per_meter=source_strength.*skimmer_total_area
% volume = 0.
% weight_factor= source_strength .* dt ./ (nP .* volume);

weight_factor= source_strength .* dt ./ (nP)


% --- IMPORTANT (Cartesian 3D) ---
% Out-of-plane thickness [m] for converting 2D (x-z) bins into 3D volumes.
% If you don't know it, set dy=1 and interpret density as per-meter in y.
dy = 1.0;  % [m]  <-- SET THIS to your box thickness in the 3rd direction

% ---- File paths ----
path1 = "../tilted_targets/0_degrees/MPEX_runs/";
file1 = "gitrm-spec.nc";

%% =========================
% 1) Read GITRm spec data
%% =========================
gridr_2d = ncread(path1 + file1, 'gridr_2d');  % x-like axis [m]
gridz_2d = ncread(path1 + file1, 'gridz_2d');  % z axis [m]
[x, z2d] = meshgrid(gridr_2d, gridz_2d);

data_1    = ncread(path1 + file1, 'n_2d');     % typically [nR x nZ x nCharge x nSpec]

% --- Show raw histogram for first channel ---
data_trim = squeeze(data_1(:,:,1,1)).';        % -> [nZ x nR]
data_trim(data_trim <= 0) = NaN;

figure('Name','n_2d (raw histogram)');
imagesc(gridr_2d, gridz_2d, data_trim);
set(gca,'YDir','normal');
xlabel('x [m]'); ylabel('z [m]');
title('n\_2d histogram (<=0 -> NaN)');
colorbar;
xlim([0 0.05]);

%% =========================
% 2) Cartesian cell volumes from edges (NO azimuthal symmetry)
%    volume = dx * dz * dy
%% =========================
xr = gridr_2d(:);
z  = gridz_2d(:);

if numel(xr) < 2 || numel(z) < 2
    error('gridr_2d and gridz_2d must each have at least 2 points.');
end

dx = diff(xr);
dz = diff(z);

x_edges = [xr(1) - dx(1)/2; (xr(1:end-1)+xr(2:end))/2; xr(end) + dx(end)/2];
z_edges = [z(1)  - dz(1)/2; (z(1:end-1)+z(2:end))/2; z(end)  + dz(end)/2];

dX = (x_edges(2:end) - x_edges(1:end-1));   % [nR x 1]
dZ = (z_edges(2:end) - z_edges(1:end-1));   % [nZ x 1]

% Cartesian volume per (z,x) bin:
volume = (dZ) .* (dX.') .* dy;              % [nZ x nR]

%% =========================
% 3) Convert histogram -> density [m^-3]
%% =========================
density = source_strength .* data_trim .* dt ./ (nP .* volume);
density(~isfinite(density)) = 0;

%% =========================
% 4) Plot density map (x-z)
%% =========================
figure('Name','Ta Density Map (x-z)');
h = pcolor(x, z2d, density);
set(h, 'EdgeColor', 'none');
xlabel('x [m]'); ylabel('z [m]');
axis tight;
c = colorbar;
ylabel(c, 'Ta Density [m^{-3}]');
title('45° tilt (Cartesian volume)');
set(gca, 'ColorScale', 'log', 'FontSize', 14);
xlim([0 0.05]);
clim([1e12 1e17])

%% =========================
% 5) Plot Ta density maps by charge state (resolved only)
%     index = 1 -> Sum, 2..6 -> Ta^0 .. Ta^4+
%% =========================
figure('Name','Ta Density Maps by Charge State','Position',[100 100 1200 800]);

charge_labels = { ...
    '\Sigma Ta', ...
    'Ta^0', ...
    'Ta^+', ...
    'Ta^{2+}', ...
    'Ta^{3+}', ...
    'Ta^{4+}'};

for i = 1:6
    data_trim = squeeze(data_1(:,:,i,1)).';  % [nZ x nR]
    data_trim(data_trim <= 0) = NaN;

    density_i = source_strength .* data_trim .* dt ./ (nP .* volume);
    density_i(~isfinite(density_i)) = 0;

    subplot(2,3,i)
    h = pcolor(x, z2d, density_i);
    set(h,'EdgeColor','none');
    set(gca,'YDir','normal','ColorScale','log','FontSize',16);
    axis tight;
    xlim([0 0.05]);
    clim([1e12 1e17])

    xlabel('x [m]');
    ylabel('z [m]');
    title(charge_labels{i});
end

cb = colorbar('Position',[0.92 0.11 0.02 0.78]);
ylabel(cb,'Ta Density [m^{-3}]');
colormap(jet)

sgtitle('Ta Density Maps: \Sigma Ta and Charge-Resolved States (45° Tilt)', ...
        'FontSize',16);

%% =========================
% 6) Consistency check: sum(2..6) vs channel 1
%% =========================
density_sum = 0;

for i = 2:6
    tmp = squeeze(data_1(:,:,i,1)).';
    tmp(tmp <= 0) = NaN;
    density_sum = density_sum + source_strength .* tmp .* dt ./ (nP .* volume);
end
density_sum(~isfinite(density_sum)) = 0;

density_all = squeeze(data_1(:,:,1,1)).';
density_all(density_all <= 0) = NaN;
density_all = source_strength .* density_all .* dt ./ (nP .* volume);
density_all(~isfinite(density_all)) = 0;

fprintf('Max relative diff (sum vs index=1): %.2e\n', ...
    max(abs(density_sum(:)-density_all(:)) ./ max(density_all(:),1)));