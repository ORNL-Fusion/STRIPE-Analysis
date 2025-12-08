%% Initialization
close all;
clear all;
clc;

%% Load SOLPS Data
fileSOLPS = 'interpolated_values.nc';
rS = ncread(fileSOLPS, 'gridr');    % Radial grid (Nr)
zS = ncread(fileSOLPS, 'gridz');    % Vertical grid (Nz)
neS = ncread(fileSOLPS, 'ne');      % Electron density (Nr x Nz)

%% Data Cleaning
n_min = 1e10;                      % Minimum density floor
neS(neS <= 0 | ~isfinite(neS)) = NaN;  % Replace invalid values with NaN

% Read EFIT data
read_efit_data;

% Additional Cleaning
neS(neS == 0) = NaN; % Ensure NaNs instead of zeros

%% Grid Preparation
rkron = kron(rS, ones(size(zS))')';
zkron = kron(zS, ones(size(rS))');

% Increase resolution of the mesh grid
num_points = 500; % Adjust resolution as needed
[X, Y] = meshgrid(linspace(min(g.lim(1,:)), max(g.lim(1,:)), num_points), ...
                   linspace(min(g.lim(2,:)), max(g.lim(2,:)), num_points));

% Identify valid data points
idx = find(~isnan(neS'));
val = zeros(size(X));

%% Midplane Density Extrapolation (Polynomial Fit)
figure;
mpfx = linspace(8.2, 8.25, 1000);
mpfy = mpfx * 0;  % Z = 0 midplane

fitDensityAtMidplane = interp2(rS, zS, neS', mpfx, mpfy, 'linear', NaN);
p = polyfit(mpfx, log(fitDensityAtMidplane), 1);

% Define full radial range
mpx = linspace(8.1, 9, 1000);
mpy = mpx * 0;

densityAtMidplane = interp2(rS, zS, neS', mpx, mpy, 'linear', NaN);
densityAtMidplane(isnan(densityAtMidplane)) = 0;

% Plot density profile
semilogy(mpx, densityAtMidplane);
writematrix([mpx; densityAtMidplane], 'nemidplane.csv');
hold on;

% Extrapolation using polynomial fit
interpfn = (mpx - 8.2) / (8.25 - 8.2);
interpfn(interpfn < 0) = 0;
interpfn(interpfn > 1) = 1;
extrapolatedne1d = interpfn .* exp(p(2) + mpx * p(1)) + (1 - interpfn) .* densityAtMidplane;

semilogy(mpx, extrapolatedne1d, 'k.');

%% Extrapolation on 2D Grid Using Nearest Valid Points
rcoords = rkron(idx);
zcoords = zkron(idx);
coords = [rcoords'; zcoords'];
okValues = neS';
okValues = okValues(idx);

% Compute extrapolated values
for i = 1:size(X, 1)
    for j = 1:size(X, 2)
        x = X(i, j);
        y = Y(i, j);
        [~, idx_min] = min(vecnorm(coords - [x; y])); % Find nearest valid data point
        val(i, j) = okValues(idx_min) * exp(p(1) * vecnorm(coords(:, idx_min) - [x; y]));
    end
end

%% Visualization of Extrapolated Density (Using imagesc)
figure;
% imagesc([min(X(:)) max(X(:))], [min(Y(:)) max(Y(:))], val);  
pcolor(X,Y,val);
shading interp;
set(gca, 'YDir', 'normal'); % Ensure correct axis orientation
colorbar;
set(gca, 'ColorScale', 'log');
title('Extrapolated Electron Density');
clim([10^16, 10^21]);
hold on;
plot(g.lim(1,:), g.lim(2,:), 'r');

%% Save Extrapolated Data as .mat File
save('extrapolated_data.mat');

% % Also save CSV for external usage
% writematrix(X, 'extrapolatedR.csv');
% writematrix(Y, 'extrapolatedZ.csv');
% writematrix(val, 'extrapolatedne.csv');