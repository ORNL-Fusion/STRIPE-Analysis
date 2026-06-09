clc;
clear;
close all;
format long;

%% 1) Load EFIT grid & limiter geometry
% load('extrapolated_data_196154.mat');  % expects: X, Y, g (includes g.lim)
load('solps_iter.mat');  % expects: X, Y, g (includes g.lim)
read_efit_data;
%% 2) Compute ψ_N on the extrapolated grid
psiN_flat = calc_psiN(g, X(:), Y(:), []);
psiN = reshape(psiN_flat, size(X));  % ψ_N on [X, Y] mesh

%% 3) Settings
psiN_mask = 1.0;                   % We'll evaluate inside ψ_N ≤ 1
dt = 1e-9;                         % Time window (seconds)
nP = 10e6;                          % Number of test particles
% source_strength = 9.442619262858435e+17; % # 196154 % particles/sec from GITRm
source_strength = 5.880221792431569e+19 ; % # 200882
n_Spec = 1;

%% 4) Load GITRm spec data
path1 = "../ITER_antenna_runs/";
% path1 = "../diiid-helicon/DIII-D_helicon_runs_full_196154/";
file1 = "gitrm-spec.nc";

gridr_2d = ncread(path1 + file1, 'gridr_2d');  % [nR]
gridz_2d = ncread(path1 + file1, 'gridz_2d');  % [nZ]
[x, y] = meshgrid(gridr_2d, gridz_2d);         % [nZ × nR]

data_1 = ncread(path1 + file1, 'n_2d');        % [nR × nZ × 1 × nSpec]
data_trim = squeeze(data_1(:,:,1,1))';         % [nZ × nR]
data_trim(data_trim <= 0) = NaN;               % Remove unphysical values

%% 5) Compute cell volumes [m³]
r1 = gridr_2d;
r2 = [gridr_2d(2:end); gridr_2d(end) + diff(gridr_2d(1:2))];
z1 = gridz_2d;
z2 = [gridz_2d(2:end); gridz_2d(end) + diff(gridz_2d(1:2))];
volume = 2 * pi * (z2 - z1) .* (0.5 * (r2.^2 - r1.^2));  % [nZ × 1]
volume = repmat(volume, 1, length(r1));  % [nZ × nR]

%% 6) Compute physical carbon density [particles/m³]
density = source_strength .* data_trim .* dt ./ (nP .* volume);
density(isnan(density)) = 0;
volume(isnan(volume)) = 0;

%% 7) Interpolate ψ_N to (x, y) grid
psiN_interp = interp2(X, Y, psiN, x, y, 'linear', nan);
mask_inside = (psiN_interp <= psiN_mask);
mask_inside(isnan(mask_inside)) = 0;

%% 8) Apply limiter mask (set density=0 outside limiter)
mask_in_limiter = inpolygon(x, y, g.lim(1,:), g.lim(2,:));
density(~mask_in_limiter) = 0;
volume(~mask_in_limiter) = 0;
mask_inside = mask_inside & mask_in_limiter;

%% 9) Integrate carbon atom count and average densities
total_atoms = sum(density(:) .* volume(:), 'omitnan');
total_vol   = sum(volume(:), 'omitnan');
avg_density_total = total_atoms / total_vol;

inside_atoms = sum(density(mask_inside) .* volume(mask_inside), 'omitnan');
inside_vol   = sum(volume(mask_inside), 'omitnan');
avg_density_inside = inside_atoms / inside_vol;

%% 10) Print results
fprintf('=== Carbon Density Integration ===\n');
fprintf('Total integrated atoms        : %.3e particles\n', total_atoms);
fprintf('Total volume                  : %.3e m^3\n', total_vol);
fprintf('Average density (full domain) : %.3e particles/m^3\n\n', avg_density_total);

fprintf('Inside ψ_N ≤ %.2f:\n', psiN_mask);
fprintf('  Integrated atoms            : %.3e particles\n', inside_atoms);
fprintf('  Volume inside ψ_N ≤ %.2f    : %.3e m^3\n', psiN_mask, inside_vol);
fprintf('  Average density (ψ_N ≤ %.2f): %.3e particles/m^3\n', psiN_mask, avg_density_inside);
fprintf('  Fraction of atoms inside    : %.2f %%\n', 100 * inside_atoms / total_atoms);

%% 11) Plotting
figure('Name','Carbon Density Map');
h = pcolor(x, y, density);
set(h, 'EdgeColor', 'none');
xlabel('r [m]'); ylabel('z [m]');
axis equal tight;
% xlim([1 2.5]); ylim([-1.5 1.5]);
c = colorbar;
ylabel(c, 'Carbon Density [m^{-3}]');
title('Carbon Density from GITRm');
set(gca, 'ColorScale', 'log', 'FontSize', 16);
hold on;
contour(X, Y, psiN, [psiN_mask psiN_mask], 'k--', 'LineWidth', 1.5);
plot(g.lim(1,:), g.lim(2,:), 'k', 'LineWidth', 1.0);
legend({'\psi_N = 1 contour', 'Limiter'}, 'Location', 'best');