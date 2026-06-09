clc;
clear;
close all;
format long;

%% =========================
% 1) LOAD EFIT DATA
% ==========================
read_efit_data;
assert(exist('g','var')==1, 'g (gfile structure) not found. read_efit_data must create g.');

%% =========================
% 2) USER SETTINGS
% ==========================
psiN_mask = 1.0;
dt = 1e-9;
nP = 5e6;
source_strength = 1.232362472103309e+20;
spec_index = 1;
charge_index = 1;   % make consistent with first script

% GITRm file
path1 = "../Guillaume_runs/1p1e17/WEST_runs/";
file1 = "gitrm-spec.nc";

% Geometry files
wall_file = "wall_coordinates.txt";
sep1_file = "separatrix1_coordinates.txt";
sep2_file = "separatrix2_coordinates.txt";

% Optional antenna centroid file
centroid_file = "surf_fields_rot.txt";
antenna_shrink_factor = 0.85;

%% =========================
% 3) LOAD GITR DATA
% ==========================
gridr = double(ncread(fullfile(path1,file1), 'gridr_2d'));   % [nR]
gridz = double(ncread(fullfile(path1,file1), 'gridz_2d'));   % [nZ]
[x, y] = meshgrid(gridr, gridz);                             % [nZ x nR]

n_2d_all = ncread(fullfile(path1,file1), 'n_2d');
n_2d = squeeze(n_2d_all(:,:,charge_index,spec_index));       % [nR x nZ]
data = double(n_2d).';                                       % [nZ x nR]
data(~isfinite(data) | data <= 0) = 0;

%% =========================
% 4) COMPUTE psi_N ON GITR GRID
% ==========================
have_psin = false;
psiN_gitr = nan(size(x));
mask_psin = false(size(x));

if exist('calc_psiN','file') == 2
    try
        psiN_flat = calc_psiN(g, x(:), y(:), []);
        psiN_gitr = reshape(psiN_flat, size(x));
        mask_psin = isfinite(psiN_gitr) & (psiN_gitr <= psiN_mask);
        have_psin = true;
    catch
        warning('calc_psiN evaluation failed. Will fall back to separatrix polygon if available.');
    end
end

%% =========================
% 5) BUILD 2D TOROIDAL CELL VOLUMES
% ==========================
r = gridr(:);
z = gridz(:);

dr = diff(r);
dz = diff(z);

assert(~isempty(dr), 'gridr_2d must contain at least 2 points.');
assert(~isempty(dz), 'gridz_2d must contain at least 2 points.');

r_edges = [r(1)-dr(1)/2; 0.5*(r(1:end-1)+r(2:end)); r(end)+dr(end)/2];
z_edges = [z(1)-dz(1)/2; 0.5*(z(1:end-1)+z(2:end)); z(end)+dz(end)/2];
r_edges = max(r_edges, 0);

dZ = diff(z_edges);                  % [nZ]
dA = 0.5 * diff(r_edges.^2);         % [nR]
volume2D = 2*pi * (dZ * dA.');       % [nZ x nR]

%% =========================
% 6) CONVERT TO PHYSICAL DENSITY
% ==========================
density = source_strength .* data .* dt ./ (nP .* volume2D);
density(~isfinite(density)) = 0;

%% =========================
% 7) LOAD WALL FROM TXT
% ==========================
assert(isfile(wall_file), 'Wall file not found: %s', wall_file);

wall_data = readmatrix(wall_file);
assert(size(wall_data,2) >= 2, 'Wall file must have at least 2 columns [R Z].');

wall_r = wall_data(:,1);
wall_z = wall_data(:,2);

good_wall = isfinite(wall_r) & isfinite(wall_z);
wall_r = wall_r(good_wall);
wall_z = wall_z(good_wall);

assert(~isempty(wall_r), 'Wall polygon is empty after cleaning.');
assert(numel(wall_r) == numel(wall_z), 'Wall R and Z arrays must have same length.');

if wall_r(1) ~= wall_r(end) || wall_z(1) ~= wall_z(end)
    wall_r(end+1) = wall_r(1);
    wall_z(end+1) = wall_z(1);
end

mask_wall = inpolygon(x, y, wall_r, wall_z);

%% =========================
% 8) LOAD SEPARATRIX TXT FALLBACK
% ==========================
sep1_r = [];
sep1_z = [];
sep2_r = [];
sep2_z = [];
mask_sep = false(size(x));
have_separatrix = false;

if isfile(sep1_file)
    sep1_data = readmatrix(sep1_file);
    if size(sep1_data,2) >= 2
        sep1_r = sep1_data(:,1);
        sep1_z = sep1_data(:,2);
        good1 = isfinite(sep1_r) & isfinite(sep1_z);
        sep1_r = sep1_r(good1);
        sep1_z = sep1_z(good1);
    end
end

if isfile(sep2_file)
    sep2_data = readmatrix(sep2_file);
    if size(sep2_data,2) >= 2
        sep2_r = sep2_data(:,1);
        sep2_z = sep2_data(:,2);
        good2 = isfinite(sep2_r) & isfinite(sep2_z);
        sep2_r = sep2_r(good2);
        sep2_z = sep2_z(good2);
    end
end

if ~isempty(sep1_r) && ~isempty(sep2_r)
    sep_r = [sep1_r; flipud(sep2_r)];
    sep_z = [sep1_z; flipud(sep2_z)];

    good_sep = isfinite(sep_r) & isfinite(sep_z);
    sep_r = sep_r(good_sep);
    sep_z = sep_z(good_sep);

    if sep_r(1) ~= sep_r(end) || sep_z(1) ~= sep_z(end)
        sep_r(end+1) = sep_r(1);
        sep_z(end+1) = sep_z(1);
    end

    mask_sep = inpolygon(x, y, sep_r, sep_z);
    have_separatrix = true;
end

%% =========================
% 9) OPTIONAL ANTENNA OUTLINE
% ==========================
r_centroid = [];
z_centroid = [];
antenna_r = [];
antenna_z = [];
have_antenna_outline = false;

if isfile(centroid_file)
    centroid_data = readmatrix(centroid_file);

    if size(centroid_data,2) >= 3
        z_centroid = centroid_data(:,3);
        r_centroid = sqrt(centroid_data(:,1).^2 + centroid_data(:,2).^2);

        good_centroid = isfinite(r_centroid) & isfinite(z_centroid);
        r_centroid = r_centroid(good_centroid);
        z_centroid = z_centroid(good_centroid);

        if numel(r_centroid) >= 3
            pts = unique([r_centroid z_centroid], 'rows');

            if size(pts,1) >= 3
                try
                    k = boundary(pts(:,1), pts(:,2), antenna_shrink_factor);
                catch
                    k = convhull(pts(:,1), pts(:,2));
                end

                antenna_r = pts(k,1);
                antenna_z = pts(k,2);

                if antenna_r(1) ~= antenna_r(end) || antenna_z(1) ~= antenna_z(end)
                    antenna_r(end+1) = antenna_r(1);
                    antenna_z(end+1) = antenna_z(1);
                end

                have_antenna_outline = true;
            end
        end
    end
end

%% =========================
% 10) COMBINED MASKS
% ==========================
density(~mask_wall) = 0;
volume2D(~mask_wall) = 0;

if have_psin
    mask_inside = mask_wall & mask_psin;
elseif have_separatrix
    mask_inside = mask_wall & mask_sep;
else
    mask_inside = mask_wall;
    warning('No psi_N or separatrix available. Inside region defaults to wall interior only.');
end

%% =========================
% 11) INTEGRALS AND AVERAGES
% ==========================
total_atoms = sum(density(:) .* volume2D(:), 'omitnan');
total_vol   = sum(volume2D(:), 'omitnan');
avg_density_total = total_atoms / max(total_vol, eps);

inside_atoms = sum(density(mask_inside) .* volume2D(mask_inside), 'omitnan');
inside_vol   = sum(volume2D(mask_inside), 'omitnan');
avg_density_inside = inside_atoms / max(inside_vol, eps);

%% =========================
% 12) PRINT RESULTS
% ==========================
fprintf('=== Density Integration ===\n');
fprintf('Total integrated atoms        : %.6e particles\n', total_atoms);
fprintf('Total volume                  : %.6e m^3\n', total_vol);
fprintf('Average density (wall domain) : %.6e particles/m^3\n\n', avg_density_total);

if have_psin
    fprintf('Inside psi_N <= %.2f:\n', psiN_mask);
elseif have_separatrix
    fprintf('Inside separatrix polygon:\n');
else
    fprintf('Inside region (wall only fallback):\n');
end

fprintf('  Integrated atoms            : %.6e particles\n', inside_atoms);
fprintf('  Volume inside               : %.6e m^3\n', inside_vol);
fprintf('  Average density (inside)    : %.6e particles/m^3\n', avg_density_inside);
fprintf('  Fraction of atoms inside    : %.3f %%\n', 100 * inside_atoms / max(total_atoms, eps));

%% =========================
% 13) PLOTTING - SAME METHOD AS FIRST SCRIPT
% ==========================
density_plot = density;
density_plot(density_plot <= 0) = NaN;

figure('Name','Density Map');
h = pcolor(x, y, density_plot);
set(h, 'EdgeColor', 'none');
axis equal tight;
xlabel('R [m]');
ylabel('Z [m]');
title('Density Map');
set(gca, 'ColorScale', 'log', 'FontSize', 16);
c = colorbar;
ylabel(c, 'Density [m^{-3}]');
hold on;

legend_entries = {'Density'};

if have_psin
    contour(x, y, psiN_gitr, [psiN_mask psiN_mask], 'k--', 'LineWidth', 1.5);
    legend_entries{end+1} = '\psi_N = 1';
elseif have_separatrix
    plot(sep1_r, sep1_z, 'k--', 'LineWidth', 1.5);
    plot(sep2_r, sep2_z, 'k--', 'LineWidth', 1.5);
    legend_entries{end+1} = 'Separatrix branches';
end

plot(wall_r, wall_z, 'k', 'LineWidth', 1.2);
legend_entries{end+1} = 'Wall';

if ~isempty(r_centroid)
    plot(r_centroid, z_centroid, '.', 'Color', [1 0.7 0.7], 'MarkerSize', 4);
end

if have_antenna_outline
    patch(antenna_r, antenna_z, 'w', ...
        'EdgeColor', 'k', ...
        'LineWidth', 2.0, ...
        'FaceAlpha', 1.0);
    legend_entries{end+1} = 'Antenna';
end

legend(legend_entries, 'Location', 'best');
xlim([min(gridr) max(gridr)]);
ylim([min(gridz) max(gridz)]);