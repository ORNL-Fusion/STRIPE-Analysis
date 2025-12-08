%% Initialization
close all;
clear all;
clc;

%% Load SOLPS Data
fileSOLPS = 'interpolated_values.nc';
rS = ncread(fileSOLPS, 'gridr');    % Radial grid (Nr)
zS = ncread(fileSOLPS, 'gridz');    % Vertical grid (Nz)
neS = ncread(fileSOLPS, 'ne');      % Electron density (Nr x Nz)
TeS = ncread(fileSOLPS, 'te');      % Electron temperature (Nr x Nz)
vrS = ncread(fileSOLPS, 'vr');      % Radial velocity
vzS = ncread(fileSOLPS, 'vz');      % Poloidal velocity
vtS = ncread(fileSOLPS, 'vt');      % Toroidal velocity
gradTiS = -ncread(fileSOLPS, 'gradTi'); % Ion temperature gradient

% Read EFIT data
read_efit_data;

%% **Grid Preparation**
num_points = 100; % Define resolution for the extrapolated grid
[X, Y] = meshgrid(linspace(min(rS), max(rS), num_points), ...
                   linspace(min(zS), max(zS), num_points));

% Flatten X, Y for easier computations
extrap_coords = [X(:), Y(:)];  
num_total_points = numel(X);

%% **Data Cleaning**
n_min = 1e10;  % Minimum density floor
T_min = 10;    % Minimum temperature floor

% Replace invalid values with NaN
neS(neS <= 0 | ~isfinite(neS)) = NaN;
TeS(TeS <= 0 | ~isfinite(TeS)) = NaN;
vrS(~isfinite(vrS)) = NaN;
vzS(~isfinite(vzS)) = NaN;
vtS(~isfinite(vtS)) = NaN;
gradTiS(~isfinite(gradTiS)) = NaN;

% Ensure minimum Te value
TeS(TeS < T_min) = T_min;

%% **Midplane Extrapolation**
mpfx = linspace(8.2, 8.25, 1000);
mpfy = mpfx * 0;  % Midplane Z = 0

% **Interpolate ne and Te using log fit**
fitDensityAtMidplane = interp2(rS, zS, neS', mpfx, mpfy, 'linear', NaN);
p_ne = polyfit(mpfx, log(fitDensityAtMidplane), 1);

fitTeAtMidplane = interp2(rS, zS, TeS', mpfx, mpfy, 'linear', NaN);
p_Te = polyfit(mpfx, log(fitTeAtMidplane), 1);

% **Interpolate velocities and gradTi using linear fit**
fitVrAtMidplane = interp2(rS, zS, vrS', mpfx, mpfy, 'linear', NaN);
p_vr = polyfit(mpfx, fitVrAtMidplane, 1);

fitVzAtMidplane = interp2(rS, zS, vzS', mpfx, mpfy, 'linear', NaN);
p_vz = polyfit(mpfx, fitVzAtMidplane, 1);

fitVtAtMidplane = interp2(rS, zS, vtS', mpfx, mpfy, 'linear', NaN);
p_vt = polyfit(mpfx, fitVtAtMidplane, 1);

fitGradTiAtMidplane = interp2(rS, zS, gradTiS', mpfx, mpfy, 'linear', NaN);
p_gradTi = polyfit(mpfx, fitGradTiAtMidplane, 1);

% **Define full radial range**
mpx = linspace(8.1, 9, 1000);
interpfn = (mpx - 8.2) / (8.25 - 8.2);
interpfn = min(max(interpfn, 0), 1); % Ensure 0 ≤ interpfn ≤ 1

% **Extrapolated ne & Te (Exponential Decay)**
extrapolatedne1d = interpfn .* exp(p_ne(2) + mpx * p_ne(1));
extrapolatedTe1d = max(interpfn .* exp(p_Te(2) + mpx * p_Te(1)), T_min);

% **Extrapolated Velocities & gradTi (Linear Decay)**
extrapolatedVr1d = interpfn .* (p_vr(2) + mpx * p_vr(1));
extrapolatedVz1d = interpfn .* (p_vz(2) + mpx * p_vz(1));
extrapolatedVt1d = interpfn .* (p_vt(2) + mpx * p_vt(1));
extrapolatedGradTi1d = interpfn .* (p_gradTi(2) + mpx * p_gradTi(1));

%% **2D Extrapolation Using `parfor`**
batch_size = 5000; 
num_batches = ceil(num_total_points / batch_size);

% **Predefine cell arrays**
val_ne_cell = cell(num_batches, 1);
val_Te_cell = cell(num_batches, 1);
val_vr_cell = cell(num_batches, 1);
val_vz_cell = cell(num_batches, 1);
val_vt_cell = cell(num_batches, 1);
val_gradTi_cell = cell(num_batches, 1);

parfor batch_idx = 1:num_batches
    batch_start = (batch_idx - 1) * batch_size + 1;
    batch_end = min(batch_idx * batch_size, num_total_points);
    
    extrap_subset = extrap_coords(batch_start:batch_end, :);  

    % Find nearest SOLPS points
    idx_nearest = knnsearch([rS(:), zS(:)], extrap_subset, 'K', 1);

    % **Apply Exponential Decay for ne & Te**
    val_ne_cell{batch_idx} = neS(idx_nearest) .* exp(-abs(extrap_subset(:,1) - rS(idx_nearest)));
    val_Te_cell{batch_idx} = max(TeS(idx_nearest) .* exp(-abs(extrap_subset(:,1) - rS(idx_nearest))), T_min);

    % **Keep Velocities & gradTi Constant**
    val_vr_cell{batch_idx} = vrS(idx_nearest);
    val_vz_cell{batch_idx} = vzS(idx_nearest);
    val_vt_cell{batch_idx} = vtS(idx_nearest);
    val_gradTi_cell{batch_idx} = gradTiS(idx_nearest);
end

% **Convert cell arrays back to numeric arrays**
val_ne = cell2mat(val_ne_cell);
val_Te = cell2mat(val_Te_cell);
val_vr = cell2mat(val_vr_cell);
val_vz = cell2mat(val_vz_cell);
val_vt = cell2mat(val_vt_cell);
val_gradTi = cell2mat(val_gradTi_cell);

% **Reshape to match grid size**
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

%% **Save Data**
save('extrapolated_data.mat', 'rS', 'zS', 'X', 'Y', 'val_ne', 'val_Te', 'val_vr', 'val_vz', 'val_vt', 'val_gradTi');