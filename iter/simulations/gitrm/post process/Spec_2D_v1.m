clc;
clear;
close all;
format long;

%% 1) Load EFIT gfile structure (must produce 'g' with fields used by calc_psiN)
% If your read_efit_data script populates g, call it:
% load('solps_iter.mat');  % optional, only if it contains g
read_efit_data;

assert(exist('g','var')==1, 'g (gfile structure) not found. read_efit_data must create g.');

%% 2) Settings
psiN_mask = 1.0;          % inside psi_N <= 1
dt = 1e-9;                % [s]
nP = 5e6;                % test particles
source_strength = 3.158559659169904e+18 ; % [particles/s]
spec_index = 1;

%% 3) Load GITRm spec data
path1 = "../ITER_antenna_runs/";
file1 = "gitrm-spec.nc";

gridr = double(ncread(fullfile(path1,file1), 'gridr_2d'));  % [nR]
gridz = double(ncread(fullfile(path1,file1), 'gridz_2d'));  % [nZ]
[x, y] = meshgrid(gridr, gridz);                            % [nZ x nR]

n_2d = ncread(fullfile(path1,file1), 'n_2d');               % often [nR x nZ x 1 x nSpec]
n_2d = squeeze(n_2d(:,:,1,spec_index));                     % [nR x nZ]
data = double(n_2d).';                                      % force to [nZ x nR]
data(~isfinite(data) | data <= 0) = 0;

%% 4) Compute psi_N directly on the GITRm grid (no X,Y needed)
psiN_flat = calc_psiN(g, x(:), y(:), []);
psiN_gitr = reshape(psiN_flat, size(x));                    % [nZ x nR]
mask_psin = isfinite(psiN_gitr) & (psiN_gitr <= psiN_mask);

%% 5) Compute 2D toroidal cell volumes [m^3] on the (gridr,gridz) mesh
r = gridr(:); z = gridz(:);
dr = diff(r); dz = diff(z);

r_edges = [r(1)-dr(1)/2; 0.5*(r(1:end-1)+r(2:end)); r(end)+dr(end)/2];  % [nR+1]
z_edges = [z(1)-dz(1)/2; 0.5*(z(1:end-1)+z(2:end)); z(end)+dz(end)/2];  % [nZ+1]
r_edges = max(r_edges, 0);

dZ = diff(z_edges);                       % [nZ]
dA = 0.5 * diff(r_edges.^2);              % [nR]
volume2D = 2*pi * (dZ * dA.');            % [nZ x nR]

%% 6) Physical density [particles/m^3]
density = source_strength .* data .* dt ./ (nP .* volume2D);
density(~isfinite(density)) = 0;

%% 7) Limiter mask (set density=0 outside limiter polygon)
if isfield(g,'lim') && ~isempty(g.lim)
    wall_r = g.lim(1,:); wall_z = g.lim(2,:);
    good = isfinite(wall_r) & isfinite(wall_z);
    wall_r = wall_r(good); wall_z = wall_z(good);

    % Close polygon for safety
    if wall_r(1) ~= wall_r(end) || wall_z(1) ~= wall_z(end)
        wall_r(end+1) = wall_r(1);
        wall_z(end+1) = wall_z(1);
    end

    mask_wall = inpolygon(x, y, wall_r, wall_z);
else
    warning('No g.lim found. Proceeding without limiter mask.');
    mask_wall = true(size(x));
    wall_r = []; wall_z = [];
end

mask_inside = mask_psin & mask_wall;

density(~mask_wall) = 0;
volume2D(~mask_wall) = 0;

%% 8) Integrals and averages
total_atoms = sum(density(:) .* volume2D(:));
total_vol   = sum(volume2D(:));
avg_density_total = total_atoms / max(total_vol, eps);

inside_atoms = sum(density(mask_inside) .* volume2D(mask_inside));
inside_vol   = sum(volume2D(mask_inside));
avg_density_inside = inside_atoms / max(inside_vol, eps);

%% 9) Print results
fprintf('=== Carbon Density Integration ===\n');
fprintf('Total integrated atoms        : %.6e particles\n', total_atoms);
fprintf('Total volume                  : %.6e m^3\n', total_vol);
fprintf('Average density (full domain) : %.6e particles/m^3\n\n', avg_density_total);

fprintf('Inside psi_N <= %.2f:\n', psiN_mask);
fprintf('  Integrated atoms            : %.6e particles\n', inside_atoms);
fprintf('  Volume inside               : %.6e m^3\n', inside_vol);
fprintf('  Average density (inside)    : %.6e particles/m^3\n', avg_density_inside);
fprintf('  Fraction of atoms inside    : %.3f %%\n', 100 * inside_atoms / max(total_atoms, eps));

% --- Remove zeros / non-positive values for log scale ---
density_plot = density;
density_plot(density_plot <= 0) = NaN;

%% 10) Plotting
figure('Name','W Density Map');
h = pcolor(x, y, density_plot);
set(h, 'EdgeColor', 'none');
axis equal tight;
xlabel('r [m]'); ylabel('z [m]');
set(gca, 'ColorScale', 'log', 'FontSize', 16);
c = colorbar; ylabel(c, 'W Density [m^{-3}]');

hold on;
contour(x, y, psiN_gitr, [psiN_mask psiN_mask], 'k--', 'LineWidth', 1.5);

if ~isempty(wall_r)
    plot(wall_r, wall_z, 'k', 'LineWidth', 1.0);
    legend({'Density','\psi_N contour','Limiter'}, 'Location','best');
else
    legend({'Density','\psi_N contour'}, 'Location','best');
end

xlim([7.4 8.4]);
ylim([-0.5 1.5])