%% Initialization
close all;
clear all;
clc;

%% Load SOLPS Data
fileSOLPS = 'interpolated_values_200882.nc';
rS = ncread(fileSOLPS, 'gridr');    % Radial grid (Nr)
zS = ncread(fileSOLPS, 'gridz');    % Vertical grid (Nz)
neS = ncread(fileSOLPS, 'ne');      % Electron density (Nr x Nz)
TeS = ncread(fileSOLPS, 'te');      % Electron temperature (Nr x Nz)
TiS = ncread(fileSOLPS, 'ti'); 
gradTiS=-ncread(fileSOLPS,'gradTi');
gradTirS=ncread(fileSOLPS,'gradTir');
vrS = ncread(fileSOLPS,'vr');
vtS = ncread(fileSOLPS,'vt');
vzS = ncread(fileSOLPS,'vz');

% Read EFIT data
read_efit_data;

%% Extract Wall Information from g.lim
r_wall = max(g.lim(1,:));  % Maximum radial position of the wall

%% Preview SOLPS Data at Outer Midplane (OMP)
figure;

% Identify OMP index (middle of zS)
[~, ompIndex] = min(abs(zS)); % Assuming OMP is at the midplane

% Plot Electron Density (ne) at OMP (Log Scale)
subplot(2,3,1);
semilogy(rS, neS(:,ompIndex), 'b', 'LineWidth', 2);
xlabel('Radial Position (m)');
ylabel('n_e (m^{-3})');
title('Electron Density at OMP');
grid on;
xline(r_wall, '--k', 'Wall Boundary', 'LineWidth', 2);
% xlim([8 8.6]);  % Zoom into 8 to 8.6m

% Plot Electron Temperature (Te) at OMP (Log Scale)
subplot(2,3,2);
semilogy(rS, TeS(:,ompIndex), 'r', 'LineWidth', 2);
xlabel('Radial Position (m)');
ylabel('T_e (eV)');
title('Electron Temperature at OMP');
grid on;
xline(r_wall, '--k', 'Wall Boundary', 'LineWidth', 2);
% xlim([8 8.6]);

% Plot Temperature Gradient (gradTi) at OMP (Linear Scale)
subplot(2,3,3);
plot(rS, gradTiS(:,ompIndex), 'g', 'LineWidth', 2);
xlabel('Radial Position (m)');
ylabel('\nabla T_i');
title('Ion Temperature Gradient at OMP');
grid on;
xline(r_wall, '--k', 'Wall Boundary', 'LineWidth', 2);
% xlim([8 8.6]);

% Plot Radial Velocity (vr) at OMP (Linear Scale)
subplot(2,3,4);
plot(rS, vrS(:,ompIndex), 'm', 'LineWidth', 2);
xlabel('Radial Position (m)');
ylabel('v_r (m/s)');
title('Radial Velocity at OMP');
grid on;
xline(r_wall, '--k', 'Wall Boundary', 'LineWidth', 2);
% xlim([8 8.6]);

% Plot Poloidal Velocity (vt) at OMP (Linear Scale)
subplot(2,3,5);
plot(rS, vtS(:,ompIndex), 'c', 'LineWidth', 2);
xlabel('Radial Position (m)');
ylabel('v_t (m/s)');
title('Poloidal Velocity at OMP');
grid on;
xline(r_wall, '--k', 'Wall Boundary', 'LineWidth', 2);
xlim([8 8.6]);

% Plot Parallel Velocity (vz) at OMP (Linear Scale)
subplot(2,3,6);
plot(rS, vzS(:,ompIndex), 'k', 'LineWidth', 2);
xlabel('Radial Position (m)');
ylabel('v_z (m/s)');
title('Parallel Velocity at OMP');
grid on;
xline(r_wall, '--k', 'Wall Boundary', 'LineWidth', 2);
% xlim([8 8.6]);

sgtitle('SOLPS Data Preview at Outer Midplane (OMP) - Zoomed (8m to 8.6m)');
%% Data Cleaning
n_min = 1e10;  % Minimum density floor
T_min = 10;    % Minimum temperature floor

% Replace invalid or out-of-range values with NaN
neS(neS <= 0 | ~isfinite(neS)) = NaN;
TeS(TeS <= 0 | ~isfinite(TeS)) = NaN;
vrS(~isfinite(vrS)) = NaN;
vzS(~isfinite(vzS)) = NaN;
vtS(~isfinite(vtS)) = NaN;
gradTiS(~isfinite(gradTiS)) = NaN;

% Enforce minimum Te
TeS(TeS < T_min) = T_min;


%% Grid Preparation
rkron = kron(rS, ones(size(zS))')';
zkron = kron(zS, ones(size(rS))');

% Increase resolution of the mesh grid
num_points_X = 1000;
num_points_Y = 1000;
[X, Y] = meshgrid(linspace(min(g.lim(1,:)), max(g.lim(1,:)), num_points_X), ...
                   linspace(min(g.lim(2,:)), max(g.lim(2,:)), num_points_Y));

% Identify valid data points
idx = find(~isnan(neS'));
rcoords = rkron(idx);
zcoords = zkron(idx);
coords = [rcoords, zcoords];  % Ensure `coords` is N x 2 (for knnsearch)

% okValues_ne = neS';
% okValues_ne = okValues_ne(idx);
% 
% okValues_Te = TeS';
% okValues_Te = okValues_Te(idx);

okValues_ne = neS';
okValues_ne = okValues_ne(idx);
okValues_Te = TeS';
okValues_Te = okValues_Te(idx);
okValues_vr = vrS';
okValues_vr = okValues_vr(idx);
okValues_vz = vzS';
okValues_vz = okValues_vz(idx);
okValues_vt = vtS';
okValues_vt = okValues_vt(idx);
okValues_gradTi = gradTiS';
okValues_gradTi = okValues_gradTi(idx);

%% Midplane Density & Temperature Extrapolation (Polynomial Fit)
mpfx = linspace(2.28, 2.31, 1000);
mpfy = mpfx * 0;  % Midplane Z = 0
% 
% % Interpolate electron density and temperature along midplane
% fitDensityAtMidplane = interp2(rS, zS, neS', mpfx, mpfy, 'linear', NaN);
% p_ne = polyfit(mpfx, log(fitDensityAtMidplane), 1);
% 
% fitTeAtMidplane = interp2(rS, zS, TeS', mpfx, mpfy, 'linear', NaN);
% p_Te = polyfit(mpfx, log(fitTeAtMidplane), 1);
% 
% Define full radial range
mpx = linspace(2.28, 2.5, 1000);
% 
% densityAtMidplane = interp2(rS, zS, neS', mpx, mpfx * 0, 'linear', NaN);
% densityAtMidplane(isnan(densityAtMidplane)) = 0;
% 
% TeAtMidplane = interp2(rS, zS, TeS', mpx, mpfx * 0, 'linear', NaN);
% TeAtMidplane(isnan(TeAtMidplane)) = T_min; % Ensure minimum Te is 1
% 
% % **Interpolate velocities and gradTi using linear fit**
% fitVrAtMidplane = interp2(rS, zS, vrS', mpfx, mpfy, 'linear', NaN);
% p_vr = polyfit(mpfx, fitVrAtMidplane, 1);
% 
% fitVzAtMidplane = interp2(rS, zS, vzS', mpfx, mpfy, 'linear', NaN);
% p_vz = polyfit(mpfx, fitVzAtMidplane, 1);
% 
% fitVtAtMidplane = interp2(rS, zS, vtS', mpfx, mpfy, 'linear', NaN);
% p_vt = polyfit(mpfx, fitVtAtMidplane, 1);
% 
% fitGradTiAtMidplane = interp2(rS, zS, gradTiS', mpfx, mpfy, 'linear', NaN);
% p_gradTi = polyfit(mpfx, fitGradTiAtMidplane, 1);
% 
% % Extrapolation function
% interpfn = (mpx - 2.28) / (2.5 - 2.28);
% interpfn = min(max(interpfn, 0), 1); % Ensure 0 ≤ interpfn ≤ 1
% 
% % Extrapolated density & temperature
% extrapolatedne1d = interpfn .* exp(p_ne(2) + mpx * p_ne(1)) + (1 - interpfn) .* densityAtMidplane;
% extrapolatedTe1d = max(interpfn .* exp(p_Te(2) + mpx * p_Te(1)) + (1 - interpfn) .* TeAtMidplane, T_min);
% 
% % **Extrapolated Velocities & gradTi (Linear Decay)**
% extrapolatedVr1d = interpfn .* (p_vr(2) + mpx * p_vr(1));
% extrapolatedVz1d = interpfn .* (p_vz(2) + mpx * p_vz(1));
% extrapolatedVt1d = interpfn .* (p_vt(2) + mpx * p_vt(1));
% extrapolatedGradTi1d = interpfn .* (p_gradTi(2) + mpx * p_gradTi(1));
% 

%% --- Midplane Polynomial Fit with Filtering ---
fitDensityAtMidplane = interp2(rS, zS, neS', mpfx, mpfy, 'linear', NaN);
valid_ne = isfinite(fitDensityAtMidplane) & fitDensityAtMidplane > 0;
p_ne = polyfit(mpfx(valid_ne), log(fitDensityAtMidplane(valid_ne)), 1);

fitTeAtMidplane = interp2(rS, zS, TeS', mpfx, mpfy, 'linear', NaN);
valid_Te = isfinite(fitTeAtMidplane) & fitTeAtMidplane > 0;
p_Te = polyfit(mpfx(valid_Te), log(fitTeAtMidplane(valid_Te)), 1);

% **Interpolate velocities and gradTi using linear fit**
fitVrAtMidplane = interp2(rS, zS, vrS', mpfx, mpfy, 'linear', NaN);
p_vr = polyfit(mpfx, fitVrAtMidplane, 1);

fitVzAtMidplane = interp2(rS, zS, vzS', mpfx, mpfy, 'linear', NaN);
p_vz = polyfit(mpfx, fitVzAtMidplane, 1);

fitVtAtMidplane = interp2(rS, zS, vtS', mpfx, mpfy, 'linear', NaN);
p_vt = polyfit(mpfx, fitVtAtMidplane, 1);

fitGradTiAtMidplane = interp2(rS, zS, gradTiS', mpfx, mpfy, 'linear', NaN);
p_gradTi = polyfit(mpfx, fitGradTiAtMidplane, 1);

%% --- 1D Extrapolation with Flooring ---
densityAtMidplane = interp2(rS, zS, neS', mpx, mpfx * 0, 'linear', NaN);
densityAtMidplane(isnan(densityAtMidplane)) = 0;

TeAtMidplane = interp2(rS, zS, TeS', mpx, mpfx * 0, 'linear', NaN);
TeAtMidplane(isnan(TeAtMidplane)) = T_min;

interpfn = (mpx - 2.28) / (2.31 - 2.28);
interpfn = min(max(interpfn, 0), 1);

extrapolatedne1d = max(interpfn .* exp(p_ne(2) + mpx * p_ne(1)) + ...
                       (1 - interpfn) .* densityAtMidplane, n_min);
extrapolatedTe1d = max(interpfn .* exp(p_Te(2) + mpx * p_Te(1)) + ...
                       (1 - interpfn) .* TeAtMidplane, T_min);
%% **Plot Midplane Extrapolations for All Variables**
% figure;
% semilogy(mpx, extrapolatedne1d, 'b', 'DisplayName', 'Extrapolated Density');
% hold on;
% semilogy(mpx, extrapolatedTe1d, 'r', 'DisplayName', 'Extrapolated Temperature');
% plot(mpx, extrapolatedVr1d, 'g', 'DisplayName', 'Extrapolated Vr');
% plot(mpx, extrapolatedVz1d, 'm', 'DisplayName', 'Extrapolated Vz');
% plot(mpx, extrapolatedVt1d, 'c', 'DisplayName', 'Extrapolated Vt');
% plot(mpx, extrapolatedGradTi1d, 'k', 'DisplayName', 'Extrapolated Grad-Ti');
% 
% legend;
% xlabel('R [m]');
% ylabel('Value');
% title('Midplane Extrapolations for All Variables');
% hold off;
% 
% % Save midplane data
% writematrix([mpx; densityAtMidplane], 'nemidplane.csv');
% writematrix([mpx; TeAtMidplane], 'Temidplane.csv');
% 
%% Efficient Extrapolation on 2D Grid Using `parfor` and Cell Arrays
batch_size = 5000; % Process data in batches to reduce memory load
num_points = numel(X);
% 
% % Convert `X, Y` into an N x 2 matrix
extrap_coords = [X(:), Y(:)];  
% 
% % Create cell arrays for parallel storage
% val_ne_cell = cell(ceil(num_points / batch_size), 1);
% val_Te_cell = cell(ceil(num_points / batch_size), 1);
% 
% parfor batch_idx = 1:ceil(num_points / batch_size)
%     batch_start = (batch_idx - 1) * batch_size + 1;
%     batch_end = min(batch_idx * batch_size, num_points);
% 
%     % Get subset of extrapolation points (Ensure format is N x 2)
%     extrap_subset = extrap_coords(batch_start:batch_end, :);  
% 
%     % Find nearest points efficiently (Ensure inputs are both N x 2)
%     idx_nearest = knnsearch(coords, extrap_subset, 'K', 1);
% 
%     % Compute distances
%     distances = vecnorm(coords(idx_nearest, :) - extrap_subset, 2, 2);
% 
%     % Compute extrapolated values
%     val_ne_cell{batch_idx} = okValues_ne(idx_nearest) .* exp(p_ne(1) * distances);
%     val_Te_cell{batch_idx} = max(okValues_Te(idx_nearest) .* exp(p_Te(1) * distances), T_min); % Ensure Te ≥ 1
% 
%     val_vr_cell{batch_idx} = okValues_vr(idx_nearest) .* (1 - distances/max(distances));
%     val_vz_cell{batch_idx} = okValues_vz(idx_nearest) .* (1 - distances/max(distances));
%     val_vt_cell{batch_idx} = okValues_vt(idx_nearest) .* (1 - distances/max(distances));
%     val_gradTi_cell{batch_idx} = okValues_gradTi(idx_nearest) .* (1 - distances/max(distances));
% end

%% --- 2D Extrapolation with Controlled Decay ---
L_decay = 0.02;  % decay length in meters

parfor batch_idx = 1:ceil(num_points / batch_size)
    batch_start = (batch_idx - 1) * batch_size + 1;
    batch_end = min(batch_idx * batch_size, num_points);
    
    extrap_subset = extrap_coords(batch_start:batch_end, :);
    idx_nearest = knnsearch(coords, extrap_subset, 'K', 1);
    distances = vecnorm(coords(idx_nearest, :) - extrap_subset, 2, 2);

    % Controlled exponential decay
    val_ne_cell{batch_idx} = max(okValues_ne(idx_nearest) .* exp(-distances / L_decay), n_min);
    val_Te_cell{batch_idx} = max(okValues_Te(idx_nearest) .* exp(-distances / L_decay), T_min);

    % Linear decay for others
    norm_factor = max(distances);
    decay_linear = max(1 - distances / norm_factor, 0);
    val_vr_cell{batch_idx} = okValues_vr(idx_nearest) .* decay_linear;
    val_vz_cell{batch_idx} = okValues_vz(idx_nearest) .* decay_linear;
    val_vt_cell{batch_idx} = okValues_vt(idx_nearest) .* decay_linear;
    val_gradTi_cell{batch_idx} = max(okValues_gradTi(idx_nearest) .* decay_linear, 0);
end

% Convert cell arrays back to matrices
val_ne = cell2mat(val_ne_cell);
val_Te = cell2mat(val_Te_cell);

val_vr = cell2mat(val_vr_cell);
val_vz = cell2mat(val_vz_cell);
val_vt = cell2mat(val_vt_cell);
val_gradTi = cell2mat(val_gradTi_cell);


% Reshape to match grid size
val_ne = reshape(val_ne, size(X));
val_Te = reshape(val_Te, size(X));

val_vr = reshape(val_vr, size(X));
val_vz = reshape(val_vz, size(X));

val_vt = reshape(val_vt, size(X));
val_gradTi = reshape(val_gradTi, size(X));


%% **Visualization**
variables = {val_ne, val_Te, val_vr, val_vz, val_vt, val_gradTi};
var_names = {'Electron Density (neS)', 'Electron Temperature (TeS)', ...
             'Radial Velocity (vrS)', 'Poloidal Velocity (vzS)', ...
             'Toroidal Velocity (vtS)', 'Ion Temperature Gradient (gradTiS)'};

for i = 1:length(variables)
    figure;
    imagesc([min(X(:)) max(X(:))], [min(Y(:)) max(Y(:))], variables{i});
    set(gca, 'YDir', 'normal');
    colorbar;
    title(['Extrapolated ', var_names{i}]);
    hold on;
    plot(g.lim(1,:), g.lim(2,:), 'r');
end


%% **Apply Mask AFTER 2D Extrapolation**
% %% **Apply Mask ONLY Beyond Limiter Boundary (on X, Y grid)**
% % Create a mask for points outside the limiter in the extrapolated grid (X, Y)
% outside_limiter_mask_XY = ~inpolygon(X(:), Y(:), g.lim(1,:), g.lim(2,:));  
% outside_limiter_mask_XY = reshape(outside_limiter_mask_XY, size(X));  % Ensure same shape as val_* variables
% 
% % Apply the mask to extrapolated values **only outside the limiter**
% val_ne(outside_limiter_mask_XY) = NaN;
% val_Te(outside_limiter_mask_XY) = NaN;
% val_vr(outside_limiter_mask_XY) = NaN;
% val_vz(outside_limiter_mask_XY) = NaN;
% val_vt(outside_limiter_mask_XY) = NaN;
% val_gradTi(outside_limiter_mask_XY) = NaN;

%% **Visualization**
variables = {val_ne, val_Te, val_vr, val_vz, val_vt, val_gradTi};
var_names = {'Electron Density (neS)', 'Electron Temperature (TeS)', ...
             'Radial Velocity (vrS)', 'Poloidal Velocity (vzS)', ...
             'Toroidal Velocity (vtS)', 'Ion Temperature Gradient (gradTiS)'};

for i = 1:length(variables)
    figure;
    imagesc([min(X(:)) max(X(:))], [min(Y(:)) max(Y(:))], variables{i});
    set(gca, 'YDir', 'normal');
    colorbar;
    title(['Extrapolated ', var_names{i}]);
    hold on;
    plot(g.lim(1,:), g.lim(2,:), 'r');
end


%% **Save Extrapolated Data**
save('extrapolated_data_196154.mat');

disp('Extrapolation complete with mask applied only beyond limiter boundary.');
