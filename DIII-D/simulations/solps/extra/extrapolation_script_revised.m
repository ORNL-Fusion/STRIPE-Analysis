
%% Initialization
close all;
clear all;
clc;

%% Load SOLPS Data
fileSOLPS = 'interpolated_values_196154.nc';
rS = ncread(fileSOLPS, 'gridr');
zS = ncread(fileSOLPS, 'gridz');
neS = ncread(fileSOLPS, 'ne');
TeS = ncread(fileSOLPS, 'te');
TiS = ncread(fileSOLPS, 'ti');
gradTiS = -ncread(fileSOLPS, 'gradTi');
gradTeS = ncread(fileSOLPS, 'gradTe');
vrS = ncread(fileSOLPS, 'vr');
vtS = ncread(fileSOLPS, 'vt');
vzS = ncread(fileSOLPS, 'vz');

read_efit_data;

%% Wall boundary
r_wall = max(g.lim(1,:));

%% Grid flattening
rkron = kron(rS, ones(size(zS))')';
zkron = kron(zS, ones(size(rS))');
idx = find(~isnan(neS'));
coords = [rkron(idx), zkron(idx)];

okValues_ne = neS'; okValues_ne = okValues_ne(idx);
okValues_Te = TeS'; okValues_Te = okValues_Te(idx);
okValues_vr = vrS'; okValues_vr = okValues_vr(idx);
okValues_vz = vzS'; okValues_vz = okValues_vz(idx);
okValues_vt = vtS'; okValues_vt = okValues_vt(idx);
okValues_gradTi = gradTiS'; okValues_gradTi = okValues_gradTi(idx);
okValues_gradTe = gradTeS'; okValues_gradTe = okValues_gradTe(idx);

%% Midplane fit range
mpfx = linspace(2.14333, 2.20675, 1000);
mpfy = mpfx * 0;

fitDensityAtMidplane = interp2(rS, zS, neS', mpfx, mpfy, 'linear', NaN);
p_ne = polyfit(mpfx(isfinite(fitDensityAtMidplane)), ...
               log(fitDensityAtMidplane(isfinite(fitDensityAtMidplane))), 1);

fitTeAtMidplane = interp2(rS, zS, TeS', mpfx, mpfy, 'linear', NaN);
p_Te = polyfit(mpfx(isfinite(fitTeAtMidplane)), ...
               log(fitTeAtMidplane(isfinite(fitTeAtMidplane))), 1);

fitVrAtMidplane = interp2(rS, zS, vrS', mpfx, mpfy, 'linear', NaN);
valid_vr = isfinite(fitVrAtMidplane);
p_vr = polyfit(mpfx(valid_vr), fitVrAtMidplane(valid_vr), 1);

fitVzAtMidplane = interp2(rS, zS, vzS', mpfx, mpfy, 'linear', NaN);
valid_vz = isfinite(fitVzAtMidplane);
p_vz = polyfit(mpfx(valid_vz), fitVzAtMidplane(valid_vz), 1);

fitVtAtMidplane = interp2(rS, zS, vtS', mpfx, mpfy, 'linear', NaN);
valid_vt = isfinite(fitVtAtMidplane);
p_vt = polyfit(mpfx(valid_vt), fitVtAtMidplane(valid_vt), 1);

fitGradTiAtMidplane = interp2(rS, zS, gradTiS', mpfx, mpfy, 'linear', NaN);
valid_gti = isfinite(fitGradTiAtMidplane);
p_gradTi = polyfit(mpfx(valid_gti), fitGradTiAtMidplane(valid_gti), 1);

fitGradTeAtMidplane = interp2(rS, zS, gradTeS', mpfx, mpfy, 'linear', NaN);
valid_gte = isfinite(fitGradTeAtMidplane);
p_gradTe = polyfit(mpfx(valid_gte), fitGradTeAtMidplane(valid_gte), 1);

%% Extrapolated profile generation
mpx = linspace(2.14333, 2.5, 1000);
interpfn = min(max((mpx - 2.14333) / (2.20675 - 2.14333), 0), 1);

extrapolatedne1d = interpfn .* exp(p_ne(2) + mpx * p_ne(1));
extrapolatedTe1d = max(interpfn .* exp(p_Te(2) + mpx * p_Te(1)), 10);
extrapolatedVr1d = interpfn .* (p_vr(2) + mpx * p_vr(1));
extrapolatedVz1d = interpfn .* (p_vz(2) + mpx * p_vz(1));
extrapolatedVt1d = interpfn .* (p_vt(2) + mpx * p_vt(1));
extrapolatedGradTi1d = interpfn .* (p_gradTi(2) + mpx * p_gradTi(1));
extrapolatedGradTe1d = interpfn .* (p_gradTe(2) + mpx * p_gradTe(1));

%% Grid setup
[X, Y] = meshgrid(linspace(min(g.lim(1,:)), max(g.lim(1,:)), 1000), ...
                  linspace(min(g.lim(2,:)), max(g.lim(2,:)), 1000));
extrap_coords = [X(:), Y(:)];
num_points = numel(X);
batch_size = 5000;

% Initialize storage
val_ne_cell = cell(ceil(num_points / batch_size), 1);
val_Te_cell = cell(ceil(num_points / batch_size), 1);
val_vr_cell = cell(ceil(num_points / batch_size), 1);
val_vz_cell = cell(ceil(num_points / batch_size), 1);
val_vt_cell = cell(ceil(num_points / batch_size), 1);
val_gradTi_cell = cell(ceil(num_points / batch_size), 1);
val_gradTe_cell = cell(ceil(num_points / batch_size), 1);

%% Extrapolation with fixed decay and density floor
parfor batch_idx = 1:ceil(num_points / batch_size)
    batch_start = (batch_idx - 1) * batch_size + 1;
    batch_end = min(batch_idx * batch_size, num_points);
    extrap_subset = extrap_coords(batch_start:batch_end, :);
    idx_nearest = knnsearch(coords, extrap_subset, 'K', 1);
    distances = vecnorm(coords(idx_nearest, :) - extrap_subset, 2, 2);
    distances = min(distances, 0.1);
    L = 0.05;
    decay = max(0, 1 - distances / L);

    val_ne_raw = okValues_ne(idx_nearest);
    val_ne_raw(val_ne_raw < 1e10) = 1e10;
    val_ne_cell{batch_idx} = val_ne_raw .* exp(p_ne(1) * distances);
    val_Te_cell{batch_idx} = max(okValues_Te(idx_nearest) .* exp(p_Te(1) * distances), 10);
    val_vr_cell{batch_idx} = okValues_vr(idx_nearest) .* decay;
    val_vz_cell{batch_idx} = okValues_vz(idx_nearest) .* decay;
    val_vt_cell{batch_idx} = okValues_vt(idx_nearest) .* decay;
    val_gradTi_cell{batch_idx} = okValues_gradTi(idx_nearest) .* decay;
    val_gradTe_cell{batch_idx} = okValues_gradTe(idx_nearest) .* decay;
end

% Reshape
val_ne = reshape(cell2mat(val_ne_cell), size(X));
val_Te = reshape(cell2mat(val_Te_cell), size(X));
val_vr = reshape(cell2mat(val_vr_cell), size(X));
val_vz = reshape(cell2mat(val_vz_cell), size(X));
val_vt = reshape(cell2mat(val_vt_cell), size(X));
val_gradTi = reshape(cell2mat(val_gradTi_cell), size(X));
val_gradTe = reshape(cell2mat(val_gradTe_cell), size(X));

% Mask outside limiter
outside_mask = ~inpolygon(X(:), Y(:), g.lim(1,:), g.lim(2,:));
outside_mask = reshape(outside_mask, size(X));
val_ne(outside_mask) = NaN;
val_Te(outside_mask) = NaN;
val_vr(outside_mask) = NaN;
val_vz(outside_mask) = NaN;
val_vt(outside_mask) = NaN;
val_gradTi(outside_mask) = NaN;
val_gradTe(outside_mask) = NaN;


%% **Visualization**
variables = {val_ne, val_Te, val_vr, val_vz, val_vt, val_gradTi, val_gradTe};
var_names = {'Electron Density (neS)', 'Electron Temperature (TeS)', ...
             'Radial Velocity (vrS)', 'Poloidal Velocity (vzS)', ...
             'Toroidal Velocity (vtS)', 'Ion Temperature Gradient (gradTiS)', 'Electron Temperature Gradient (gradTiS)'};

for i = 1:length(variables)
    figure;
    imagesc([min(X(:)) max(X(:))], [min(Y(:)) max(Y(:))], variables{i});
    set(gca, 'YDir', 'normal');
    colorbar;
    title(['Extrapolated ', var_names{i}]);
    hold on;
    plot(g.lim(1,:), g.lim(2,:), 'r');
end

%% Save result
save('extrapolated_data_196154.mat');
disp('Extrapolation complete and masked beyond limiter boundary.');