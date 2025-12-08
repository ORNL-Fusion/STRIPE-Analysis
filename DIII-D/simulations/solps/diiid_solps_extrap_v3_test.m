%% Initialization
close all;
clear all;
clc;

%% Load SOLPS Data
% fileSOLPS = 'interpolated_values_200882_wallFix.nc';
fileSOLPS = 'interpolated_values_196154.nc';
rS = ncread(fileSOLPS, 'gridr');    % Radial grid (Nr)
zS = ncread(fileSOLPS, 'gridz');    % Vertical grid (Nz)
neS = ncread(fileSOLPS, 'ne');      % Electron density (Nr x Nz)
TeS = ncread(fileSOLPS, 'te');      % Electron temperature (Nr x Nz)
gradTiS=-ncread(fileSOLPS,'gradTi');
gradTirS=ncread(fileSOLPS,'gradTir');
vrS = ncread(fileSOLPS,'vr');
vtS = ncread(fileSOLPS,'vt');
vzS = ncread(fileSOLPS,'vz');

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

% Read EFIT data
read_efit_data;

%%
figure;
imagesc(rS, zS, neS');  
set(gca, 'YDir', 'normal');
colorbar;
set(gca, 'ColorScale', 'linear');
title('Extrapolated Vr');
hold on;
plot(g.lim(1,:), g.lim(2,:), 'r');

figure;
imagesc(rS, zS, vzS');  
set(gca, 'YDir', 'normal');
colorbar;
set(gca, 'ColorScale', 'linear');
title('Extrapolated Vz');
hold on;
plot(g.lim(1,:), g.lim(2,:), 'r');

figure;
imagesc(rS, zS, vtS');  
set(gca, 'YDir', 'normal');
colorbar;
set(gca, 'ColorScale', 'linear');
title('Extrapolated Vt');
hold on;
plot(g.lim(1,:), g.lim(2,:), 'r');

figure;
imagesc(rS, zS, gradTiS');  
set(gca, 'YDir', 'normal');
colorbar;
set(gca, 'ColorScale', 'linear');
title('Extrapolated Vt');
hold on;
plot(g.lim(1,:), g.lim(2,:), 'r');

%% 2) Compute psi_N on SOLPS grid (fixed)
% rS : [Nr×1], zS : [Nz×1]
[RR, ZZ] = meshgrid(rS, zS);        % RR,ZZ are [Nz×Nr]

% Flatten to 1×(Nz*Nr) so calc_psiN sees matching vectors
R_vec = RR(:).';   % 1×(Nz*Nr)
Z_vec = ZZ(:).';

% Compute psi_N
[psiN_vec, ~] = calc_psiN(g, R_vec, Z_vec, 0);

% Reshape back to the original 2D grid
psiN = reshape(psiN_vec, size(RR));   % [Nz×Nr]


%% Grid Preparation
rkron = kron(rS, ones(size(zS))')';
zkron = kron(zS, ones(size(rS))');

% Increase resolution of the mesh grid
num_points = 1000;
[X, Y] = meshgrid(linspace(min(g.lim(1,:)), max(g.lim(1,:)), num_points), ...
                   linspace(min(g.lim(2,:)), max(g.lim(2,:)), num_points));

% Identify valid data points
idx = find(~isnan(neS'));
rcoords = rkron(idx);
zcoords = zkron(idx);
coords = [rcoords, zcoords];  % Ensure `coords` is N x 2 (for knnsearch)

okValues_ne = neS';
okValues_ne = okValues_ne(idx);

okValues_Te = TeS';
okValues_Te = okValues_Te(idx);

%% Midplane Density & Temperature Extrapolation (Polynomial Fit)
% mpfx = linspace(8.2, 8.25, 1000); % ITER
mpfx = linspace(2.14333, 2.20675, 1000); % DIII-D 196154
% mpfx = linspace(2.28122, 2.3156, 1000); % DIII-D # 200882

mpfy = mpfx * 0;  % Midplane Z = 0

% Interpolate electron density and temperature along midplane
fitDensityAtMidplane = interp2(rS, zS, neS', mpfx, mpfy, 'linear', NaN);
p_ne = polyfit(mpfx, log(fitDensityAtMidplane), 1);

fitTeAtMidplane = interp2(rS, zS, TeS', mpfx, mpfy, 'linear', NaN);
p_Te = polyfit(mpfx, log(fitTeAtMidplane), 1);

% Define full radial range
% mpx = linspace(8.1, 9, 1000); % ITER
mpx = linspace(2.14333, 2.5, 1000); % DIII-D # 196154
% mpx = linspace(2.28122, 2.5, 1000); % DIII-D # 200882


densityAtMidplane = interp2(rS, zS, neS', mpx, mpfx * 0, 'linear', NaN);
densityAtMidplane(isnan(densityAtMidplane)) = 0;

TeAtMidplane = interp2(rS, zS, TeS', mpx, mpfx * 0, 'linear', NaN);
TeAtMidplane(isnan(TeAtMidplane)) = T_min; % Ensure minimum Te is 1

% Extrapolation function
% interpfn = (mpx - 8.2) / (8.25 - 8.2); % ITER
interpfn = (mpx - 2.14333) / (2.20675 - 2.14333); % DIII-D # 196154
% interpfn = (mpx - 2.281) / (2.3156 - 2.28122); % DIII-D # 200882
interpfn = min(max(interpfn, 0), 1); % Ensure 0 ≤ interpfn ≤ 1

% Extrapolated density & temperature
extrapolatedne1d = interpfn .* exp(p_ne(2) + mpx * p_ne(1)) + (1 - interpfn) .* densityAtMidplane;
extrapolatedTe1d = max(interpfn .* exp(p_Te(2) + mpx * p_Te(1)) + (1 - interpfn) .* TeAtMidplane, T_min);

% --- Plot Density and Temperature Profiles ---
figure;
semilogy(mpx, densityAtMidplane, 'b');  % Plot ne
hold on;
semilogy(mpx, extrapolatedne1d, 'b.'); % Extrapolated ne
semilogy(mpx, TeAtMidplane, 'r');  % Plot Te
semilogy(mpx, extrapolatedTe1d, 'r.'); % Extrapolated Te

legend('Density Midplane', 'Extrapolated Density', 'Temperature Midplane', 'Extrapolated Temperature');
xlabel('R [m]');
ylabel('Value');
title('Midplane Density and Temperature Extrapolation');
hold off;

% Save midplane data
writematrix([mpx; densityAtMidplane], 'nemidplane.csv');
writematrix([mpx; TeAtMidplane], 'Temidplane.csv');

%% Efficient Extrapolation on 2D Grid Using `parfor` and Cell Arrays
num_points = numel(X);

pool = gcp();  % This will start a pool if none is running
num_workers = pool.NumWorkers;

% Target 10 batches per worker
batch_size = ceil(num_points / (10 * num_workers));
% Convert `X, Y` into an N x 2 matrix
extrap_coords = [X(:), Y(:)];  
% Create cell arrays for parallel storage
val_ne_cell = cell(ceil(num_points / batch_size), 1);
val_Te_cell = cell(ceil(num_points / batch_size), 1);

parfor batch_idx = 1:ceil(num_points / batch_size)
    batch_start = (batch_idx - 1) * batch_size + 1;
    batch_end = min(batch_idx * batch_size, num_points);
    
    % Get subset of extrapolation points (Ensure format is N x 2)
    extrap_subset = extrap_coords(batch_start:batch_end, :);  

    % Find nearest points efficiently (Ensure inputs are both N x 2)
    idx_nearest = knnsearch(coords, extrap_subset, 'K', 1);
    
    % Compute distances
    distances = vecnorm(coords(idx_nearest, :) - extrap_subset, 2, 2);

    % Compute extrapolated values
    val_ne_cell{batch_idx} = okValues_ne(idx_nearest) .* exp(p_ne(1) * distances);
    val_Te_cell{batch_idx} = max(okValues_Te(idx_nearest) .* exp(p_Te(1) * distances), T_min); % Ensure Te ≥ 1
end

% Convert cell arrays back to matrices
val_ne = cell2mat(val_ne_cell);
val_Te = cell2mat(val_Te_cell);

% Reshape to match grid size
val_ne = reshape(val_ne, size(X));
val_Te = reshape(val_Te, size(X));

%% Visualization of Extrapolated Electron Density (Using imagesc)
figure;
imagesc([min(X(:)) max(X(:))], [min(Y(:)) max(Y(:))], val_ne);  
set(gca, 'YDir', 'normal');
colorbar;
set(gca, 'ColorScale', 'log');
title('Extrapolated Electron Density (n_e)');
clim([10^1, 10^20]);
hold on;
plot(g.lim(1,:), g.lim(2,:), 'r');

%% Visualization of Extrapolated Electron Temperature (Using imagesc)
figure;
imagesc([min(X(:)) max(X(:))], [min(Y(:)) max(Y(:))], val_Te);  
set(gca, 'YDir', 'normal');
colorbar;
set(gca, 'ColorScale', 'log');
title('Extrapolated Electron Temperature (Te)');
clim([1, max(val_Te(:))]);
hold on;
plot(g.lim(1,:), g.lim(2,:), 'r');

%% Grid Preparation
[ZZ, RR] = meshgrid(zS, rS);  % Create a grid matching SOLPS data dimensions

% --- Basic plasma mask (for finding edge index only) ---
plasma_mask = neS > n_min;  

% --- Exponential Decay Parameters ---
decay_length = 0.2;  % Adjust decay rate as needed

%% ----- Update Velocities and gradTi Using Exponential Decay (no limiter mask) -----
for i = 1:length(zS)
    % Last valid index for density is proxy for plasma edge
    edge_idx = find(plasma_mask(:, i), 1, 'last');  
    if ~isempty(edge_idx) && edge_idx < length(rS)
        R_edge = rS(edge_idx);

        % Edge values
        edge_vr     = vrS(edge_idx, i);
        edge_vz     = vzS(edge_idx, i);
        edge_vt     = vtS(edge_idx, i);
        edge_gradTi = gradTiS(edge_idx, i);

        % Apply exponential decay beyond edge (no mask condition)
        for j = edge_idx+1:length(rS)
            decay = exp(-(rS(j) - R_edge) / decay_length);
            vrS(j, i)     = edge_vr     * decay;
            vzS(j, i)     = edge_vz     * decay;
            vtS(j, i)     = edge_vt     * decay;
            gradTiS(j, i) = edge_gradTi * decay;
        end
    end
end

%% Plot Updated Velocities and gradTi
figure;
imagesc(rS, zS, vrS');  
set(gca, 'YDir', 'normal', 'ColorScale', 'linear');
colorbar; title('Extrapolated Radial Velocity (vr)');

figure;
imagesc(rS, zS, vzS');  
set(gca, 'YDir', 'normal', 'ColorScale', 'linear');
colorbar; title('Extrapolated Poloidal Velocity (vz)');

figure;
imagesc(rS, zS, vtS');  
set(gca, 'YDir', 'normal', 'ColorScale', 'linear');
colorbar; title('Extrapolated Toroidal Velocity (vt)');

figure;
imagesc(rS, zS, gradTiS');  
set(gca, 'YDir', 'normal', 'ColorScale', 'linear');
colorbar; title('Extrapolated Ion Temperature Gradient (gradTi)');

%%
% Interpolate onto new X-Y grid
val_vr     = interp2(rS, zS, vrS',     X, Y, 'linear');
val_vz     = interp2(rS, zS, vzS',     X, Y, 'linear');
val_vt     = interp2(rS, zS, vtS',     X, Y, 'linear');
val_gradTi = interp2(rS, zS, gradTiS', X, Y, 'linear');

%% 2D Visualization of Extrapolated Data
vars2D     = {val_ne, val_Te, val_vr, val_vz, val_vt, val_gradTi};
varNames2D = {'Electron Density (neS)', 'Electron Temperature (TeS)', ...
              'Radial Velocity (vrS)', 'Poloidal Velocity (vzS)', ...
              'Toroidal Velocity (vtS)', 'Ion Temperature Gradient (gradTiS)'};

for i = 1:length(vars2D)
    figure;
    imagesc([min(X(:)) max(X(:))], [min(Y(:)) max(Y(:))], vars2D{i});
    set(gca, 'YDir', 'normal', 'ColorScale', 'linear');
    colorbar;
    title(['Extrapolated ', varNames2D{i}]);
end

%% Save the data
save('extrapolated_data_196154.mat'); % DIII-D #196154

%% End of Script