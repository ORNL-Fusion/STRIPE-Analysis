%% SOLPS 2D Extrapolation: Local Midplane Fit + Radial Decay
% Electron density & temperature: exponential radial decay beyond SOLPS domain
% Velocities & gradT: nearest‐neighbor + linear decay

%% Initialization
close all; clear all; clc;

%% Load SOLPS Data
fileSOLPS = 'interpolated_values_196154.nc';
rS      = ncread(fileSOLPS, 'gridr');    % (Nr)
zS      = ncread(fileSOLPS, 'gridz');    % (Nz)
neS     = ncread(fileSOLPS, 'ne');      % (Nr x Nz)
TeS     = ncread(fileSOLPS, 'te');      % (Nr x Nz)
gradTiS = -ncread(fileSOLPS,'gradTi');  % (Nr x Nz)
gradTeS =  ncread(fileSOLPS,'gradTe');  % (Nr x Nz)
vrS     = ncread(fileSOLPS, 'vr');
vtS     = ncread(fileSOLPS, 'vt');
vzS     = ncread(fileSOLPS, 'vz');

% Read EFIT geometry (g.lim)
read_efit_data;
r_wall    = max(g.lim(1,:));   % vessel outer radius
rSOLPS_max = max(rS);

%% Data Cleaning & Floors
T_min = 10;
TeS(TeS < T_min) = T_min;

%% 1D Midplane (Z=0) Fit over Last 5 cm Band
[~, idx0] = min(abs(zS));         % midplane index
R_mid    = rS(:);
ne_mid   = neS(:,idx0);
Te_mid   = TeS(:,idx0);

% Band subset: last 5 cm
band = (R_mid >= (rSOLPS_max - 0.05)) & (R_mid <= rSOLPS_max);
% Log‐linear fit of band
p_ne = polyfit(R_mid(band), log(ne_mid(band)), 1);
p_Te = polyfit(R_mid(band), log(Te_mid(band)), 1);
% Boundary values
ne_wall = interp1(R_mid, ne_mid, rSOLPS_max, 'linear');
Te_wall = interp1(R_mid, Te_mid, rSOLPS_max, 'linear');

%% Build High‐Res 2D Grid
N = 1000;
[X, Y] = meshgrid(linspace(min(rS), rSOLPS_max+0.5, N), ...
                  linspace(min(zS), max(zS), N));

%% Electron density & temperature via interp2 + radial decay
ne2D = interp2(rS, zS, neS', X, Y, 'linear', NaN);
Te2D = interp2(rS, zS, TeS', X, Y, 'linear', NaN);

mask_ex = X > rSOLPS_max;

ne2D(mask_ex) = ne_wall .* exp(p_ne(1) * (X(mask_ex) - rSOLPS_max));
Te2D(mask_ex) = max(Te_wall .* exp(p_Te(1) * (X(mask_ex) - rSOLPS_max)), T_min);

%% Prepare for velocities & gradT extrapolation via knn + linear
% Flatten original SOLPS coords & values
rkron = kron(rS, ones(size(zS))')';
zkron = kron(zS, ones(size(rS))');
coords_orig = [rkron(:), zkron(:)];
mask_val    = ~isnan(neS(:));
coords       = coords_orig(mask_val, :);

okValues_vr    = vrS(:);    okValues_vr    = okValues_vr(mask_val);
okValues_vz    = vzS(:);    okValues_vz    = okValues_vz(mask_val);
okValues_vt    = vtS(:);    okValues_vt    = okValues_vt(mask_val);
okValues_gradTi = gradTiS(:); okValues_gradTi = okValues_gradTi(mask_val);
okValues_gradTe = gradTeS(:); okValues_gradTe = okValues_gradTe(mask_val);

%% Parallel knn + linear decay
pts     = numel(X);
pool    = gcp();  % open parallel pool
W       = pool.NumWorkers;
batch_sz = ceil(pts/(10*W));
coords_ex = [X(:), Y(:)];

% Pre‐allocate
val_vr_cell     = cell(ceil(pts/batch_sz),1);
val_vz_cell     = cell(ceil(pts/batch_sz),1);
val_vt_cell     = cell(ceil(pts/batch_sz),1);
val_gradTi_cell = cell(ceil(pts/batch_sz),1);
val_gradTe_cell = cell(ceil(pts/batch_sz),1);

parfor b = 1:ceil(pts/batch_sz)
    i0 = (b-1)*batch_sz + 1;
    i1 = min(b*batch_sz, pts);
    sub = coords_ex(i0:i1,:);
    idxNN = knnsearch(coords, sub);
    d     = vecnorm(coords(idxNN,:) - sub, 2, 2);
    linf  = 1 - d./max(d);

    val_vr_cell{b}     = okValues_vr(idxNN)    .* linf;
    val_vz_cell{b}     = okValues_vz(idxNN)    .* linf;
    val_vt_cell{b}     = okValues_vt(idxNN)    .* linf;
    val_gradTi_cell{b} = okValues_gradTi(idxNN) .* linf;
    val_gradTe_cell{b} = okValues_gradTe(idxNN) .* linf;
end

% Concatenate & reshape
val_vr     = reshape(cell2mat(val_vr_cell),     size(X));
val_vz     = reshape(cell2mat(val_vz_cell),     size(X));
val_vt     = reshape(cell2mat(val_vt_cell),     size(X));
val_gradTi = reshape(cell2mat(val_gradTi_cell), size(X));
val_gradTe = reshape(cell2mat(val_gradTe_cell), size(X));

%% Visualization Example
figure;
imagesc([min(X(:)) max(X(:))], [min(Y(:)) max(Y(:))], ne2D);
set(gca,'YDir','normal'); colorbar;
title('Extended n_e: interp2 + radial decay');
hold on; plot(g.lim(1,:), g.lim(2,:), 'r', 'LineWidth',1.5);

%% Save Extrapolated Data
save('extrapolated_data_196154.mat', ...
     'X','Y','ne2D','Te2D','val_vr','val_vz','val_vt','val_gradTi','val_gradTe');

disp('Extrapolation complete: ne/Te radial decay; v & gradT linear decay.');