clc;
clear;
close all;
format long;

%% =========================
% 1) USER SETTINGS
% ==========================
psiN_mask = 1.0;          % inside psi_N <= 1
dt = 1e-9;                % [s]
nP = 0.5e6;                 % number of test particles

source_strength = 5.175509e+20 ; % [particles/s] - High Density WEST
% source_strength = 2.021521e+19; % [particles/s] - Low Density WEST
spec_index = 1;           % species index in n_2d
charge_index = 1;         % charge-state slice used in n_2d(:,:,charge_index,spec_index)

% GITR file location
% path1 = "../WEST_runs/";
% path1 = "../Guillaume_runs/1p1e17/WEST_runs/"; % Guillaume Low Density
path1 = "../Guillaume_runs/7p3e18/WEST_runs/"; % Guillaume Low Density
file1 = "gitrm-spec.nc";

% Geometry files
wall_file = "wall_from_mesh.txt";
sep1_file = "separatrix1_from_mesh.txt";
sep2_file = "separatrix2_from_mesh.txt";

% Antenna centroid coordinates from COMSOL/rotated surface file
centroid_file = "surf_fields_rot.txt";

% Boundary tuning for antenna outline
antenna_shrink_factor = 0.85;   % 0-1, smaller = tighter outline

%% =========================
% 2) LOAD GITR DATA
% ==========================
gridr = double(ncread(fullfile(path1,file1), 'gridr_2d'));   % [nR]
gridz = double(ncread(fullfile(path1,file1), 'gridz_2d'));   % [nZ]
[x, y] = meshgrid(gridr, gridz);                             % [nZ x nR]

n_2d_all = ncread(fullfile(path1,file1), 'n_2d');

% Expected common shape: [nR x nZ x nCharge x nSpec]
n_2d = squeeze(n_2d_all(:,:,charge_index,spec_index));       % [nR x nZ]
data = double(n_2d).';                                       % -> [nZ x nR]

data(~isfinite(data) | data <= 0) = 0;

%% =========================
% 3) BUILD 2D TOROIDAL CELL VOLUMES
% ==========================
r = gridr(:);
z = gridz(:);

dr = diff(r);
dz = diff(z);

assert(~isempty(dr), 'gridr_2d must contain at least 2 points.');
assert(~isempty(dz), 'gridz_2d must contain at least 2 points.');

r_edges = [r(1)-dr(1)/2; 0.5*(r(1:end-1)+r(2:end)); r(end)+dr(end)/2];
z_edges = [z(1)-dz(1)/2; 0.5*(z(1:end-1)+z(2:end)); z(end)+dz(end)/2];

r_edges = max(r_edges, 0);   % avoid negative R

dZ = diff(z_edges);                  % [nZ]
dA = 0.5 * diff(r_edges.^2);         % [nR]
volume2D = 2*pi * (dZ * dA.');       % [nZ x nR]

%% =========================
% 4) CONVERT TO PHYSICAL DENSITY
% ==========================
density = source_strength .* data .* dt ./ (nP .* volume2D);
density(~isfinite(density)) = 0;

%% =========================
% 5) LOAD WALL FROM TXT
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

% Close polygon if needed
if wall_r(1) ~= wall_r(end) || wall_z(1) ~= wall_z(end)
    wall_r(end+1) = wall_r(1);
    wall_z(end+1) = wall_z(1);
end

mask_wall = inpolygon(x, y, wall_r, wall_z);

%% =========================
% 6) LOAD SEPARATRIX BRANCHES FROM TXT
% ==========================
sep1_r = [];
sep1_z = [];
sep2_r = [];
sep2_z = [];

mask_sep = false(size(x));
have_separatrix = false;
have_psin = false;
psiN_gitr = nan(size(x));

% ---- separatrix branch 1 ----
if isfile(sep1_file)
    sep1_data = readmatrix(sep1_file);
    assert(size(sep1_data,2) >= 2, 'sep1 file must have at least 2 columns [R Z].');

    sep1_r = sep1_data(:,1);
    sep1_z = sep1_data(:,2);

    good1 = isfinite(sep1_r) & isfinite(sep1_z);
    sep1_r = sep1_r(good1);
    sep1_z = sep1_z(good1);
end

% ---- separatrix branch 2 ----
if isfile(sep2_file)
    sep2_data = readmatrix(sep2_file);
    assert(size(sep2_data,2) >= 2, 'sep2 file must have at least 2 columns [R Z].');

    sep2_r = sep2_data(:,1);
    sep2_z = sep2_data(:,2);

    good2 = isfinite(sep2_r) & isfinite(sep2_z);
    sep2_r = sep2_r(good2);
    sep2_z = sep2_z(good2);
end

% Build closed polygon from lower + upper branch if both are available
if ~isempty(sep1_r) && ~isempty(sep2_r)
    sep_r = [sep1_r; flipud(sep2_r)];
    sep_z = [sep1_z; flipud(sep2_z)];

    good_sep = isfinite(sep_r) & isfinite(sep_z);
    sep_r = sep_r(good_sep);
    sep_z = sep_z(good_sep);

    if ~isempty(sep_r) && numel(sep_r) == numel(sep_z)
        if sep_r(1) ~= sep_r(end) || sep_z(1) ~= sep_z(end)
            sep_r(end+1) = sep_r(1);
            sep_z(end+1) = sep_z(1);
        end
        mask_sep = inpolygon(x, y, sep_r, sep_z);
        have_separatrix = true;
    end
end

% Priority 2: true psi_N if calc_psiN and g are available
if exist('calc_psiN','file') == 2 && exist('g','var') == 1
    try
        psiN_flat = calc_psiN(g, x(:), y(:), []);
        psiN_gitr = reshape(psiN_flat, size(x));
        mask_psin = isfinite(psiN_gitr) & (psiN_gitr <= psiN_mask);
        have_psin = true;
    catch
        warning('calc_psiN was found, but psi_N evaluation failed. Falling back to separatrix polygon if available.');
    end
end

%% =========================
% 7) LOAD ANTENNA CENTROIDS AND BUILD OUTER OUTLINE
% ==========================
r_centroid = [];
z_centroid = [];

antenna_r = [];
antenna_z = [];
have_antenna_outline = false;

if isfile(centroid_file)
    centroid_data = readmatrix(centroid_file);

    assert(size(centroid_data,2) >= 3, ...
        'Centroid file must have at least 3 columns for x, y, z.');

    % Assuming:
    % col 1 = x
    % col 2 = y
    % col 3 = z
    z_centroid = centroid_data(:,3);
    r_centroid = sqrt(centroid_data(:,1).^2 + centroid_data(:,2).^2);

    good_centroid = isfinite(r_centroid) & isfinite(z_centroid);
    r_centroid = r_centroid(good_centroid);
    z_centroid = z_centroid(good_centroid);

    fprintf('Loaded %d antenna centroid points from %s\n', numel(r_centroid), centroid_file);

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
else
    warning('Centroid file not found: %s', centroid_file);
end

%% =========================
% 8) COMBINED MASKS
% ==========================
density(~mask_wall) = 0;
volume2D(~mask_wall) = 0;

if have_psin
    mask_inside = mask_wall & mask_psin;
elseif have_separatrix
    mask_inside = mask_wall & mask_sep;
else
    mask_inside = mask_wall;
    warning('No separatrix or psi_N available. "Inside" region defaults to wall interior only.');
end

%% =========================
% 9) INTEGRALS AND AVERAGES
% ==========================
total_atoms = sum(density(:) .* volume2D(:), 'omitnan');
total_vol   = sum(volume2D(:), 'omitnan');
avg_density_total = total_atoms / max(total_vol, eps);

inside_atoms = sum(density(mask_inside) .* volume2D(mask_inside), 'omitnan');
inside_vol   = sum(volume2D(mask_inside), 'omitnan');
avg_density_inside = inside_atoms / max(inside_vol, eps);

%% =========================
% 10) PRINT RESULTS
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
% 11) PLOTTING
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

% Optional: faint centroid points for debugging
if ~isempty(r_centroid)
    plot(r_centroid, z_centroid, '.', 'Color', [1 0.7 0.7], 'MarkerSize', 4);
end

% Antenna filled region (white mask)
if have_antenna_outline
    patch(antenna_r, antenna_z, 'r', ...
        'EdgeColor', 'k', ...
        'LineWidth', 2.0, ...
        'FaceAlpha', 1.0);   % fully opaque white

    legend_entries{end+1} = 'Antenna';
end

legend(legend_entries, 'Location', 'best');
xlim([min(gridr) max(gridr)]);
ylim([min(gridz) max(gridz)]);