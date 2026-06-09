%% Normalize Te so global Te_min = 1 eV + visualize 2D + 1D (wall, OMP)
clear; clc; close all;

%% ---- User inputs ----
infile  = "profiles_iter_single.nc";
outfile = "profiles_iter_single_TeMin1eV.nc";

% Define OMP location in your grid:
% Default assumes OMP is midplane at z ~ 0
z_omp_target = 0.0;

% Define "wall" side:
% Default: wall at max(x). If your wall is at min(x), set wall_side = "min"
wall_side = "max";  % "max" or "min"

%% ---- Copy original to new file ----
copyfile(infile, outfile);

%% ---- Read variables ----
x  = ncread(outfile, "x");   % nX x 1
z  = ncread(outfile, "z");   % nZ x 1
Te = ncread(outfile, "te");  % expected nX x nZ

% Basic sanity checks
nX = numel(x);
nZ = numel(z);

if ~isequal(size(Te), [nX nZ])
    error("Unexpected Te size. Expected [%d %d], got [%d %d].", nX, nZ, size(Te,1), size(Te,2));
end

%% ---- Global scaling so min(Te) = 1 eV ----
% (Optional safety: ignore non-positive values if present)
Te_safe = Te;
Te_safe(Te_safe <= 0) = NaN;

Te_min = min(Te_safe(:), [], "omitnan");
if ~isfinite(Te_min) || Te_min <= 0
    error("Te_min is invalid (<=0 or NaN). Check the Te field.");
end

scale  = 1.0 / 4;
Te_new = Te * scale;

% Write back to outfile
ncwrite(outfile, "te", Te_new);

fprintf("Wrote normalized Te to: %s\n", outfile);
fprintf("Original Te_min = %.6g eV\n", Te_min);
fprintf("Scale factor     = %.6g\n", scale);
fprintf("New Te_min       = %.6g eV\n", min(Te_new(:), [], "omitnan"));

%% ---- Indices for OMP and wall ----
[~, iz_omp] = min(abs(z - z_omp_target));

switch lower(wall_side)
    case "max"
        ix_wall = nX;
    case "min"
        ix_wall = 1;
    otherwise
        error('wall_side must be "max" or "min".');
end

%% ---- 2D visualizations ----
% For plotting with imagesc/pcolor, it's convenient to use Z as vertical axis
% Te is (nX x nZ), so transpose for display as (z x x)
Te0_plot  = Te';      % nZ x nX
Te1_plot  = Te_new';  % nZ x nX
ratioPlot = Te1_plot ./ Te0_plot;
diffPlot  = Te1_plot - Te0_plot;

% Original Te (2D)
figure("Name","Te original (2D)");
imagesc(x, z, Te0_plot);
set(gca,"YDir","normal");
xlabel("x"); ylabel("z"); title("Original T_e(x,z) [eV]");
colorbar;

% Normalized Te (2D)
figure("Name","Te normalized (2D) - min=1eV");
imagesc(x, z, Te1_plot);
set(gca,"YDir","normal");
xlabel("x"); ylabel("z"); title("Normalized T_e(x,z) [eV], with min=1 eV");
colorbar;

% Ratio
figure("Name","Te ratio (new/orig)");
imagesc(x, z, ratioPlot);
set(gca,"YDir","normal");
xlabel("x"); ylabel("z"); title("T_e^{new} / T_e^{orig}");
colorbar;

% Difference
figure("Name","Te difference (new-orig)");
imagesc(x, z, diffPlot);
set(gca,"YDir","normal");
xlabel("x"); ylabel("z"); title("T_e^{new} - T_e^{orig} [eV]");
colorbar;

%% ---- 1D profiles ----
% 1) Wall profile: Te(z) at x = wall
Te_wall_orig = Te(ix_wall, :);      % 1 x nZ
Te_wall_new  = Te_new(ix_wall, :);  % 1 x nZ

figure("Name","Wall Te(z) profile");
plot(z, Te_wall_orig, "k-", "LineWidth", 1.5); hold on;
plot(z, Te_wall_new,  "r--", "LineWidth", 1.5);
grid on;
xlabel("z"); ylabel("T_e (eV)");
title(sprintf("Wall profile: T_e(z) at x = %s(x) (ix=%d)", wall_side, ix_wall));
legend("Original","Normalized (min=1 eV)", "Location","best");

% 2) OMP profile: Te(x) at z ~ 0
Te_omp_orig = Te(:, iz_omp);      % nX x 1
Te_omp_new  = Te_new(:, iz_omp);  % nX x 1

figure("Name","OMP Te(x) profile");
plot(x, Te_omp_orig, "k-", "LineWidth", 1.5); hold on;
plot(x, Te_omp_new,  "r--", "LineWidth", 1.5);
grid on;
xlabel("x"); ylabel("T_e (eV)");
title(sprintf("OMP profile: T_e(x) at z ≈ %.6g (iz=%d)", z(iz_omp), iz_omp));
legend("Original","Normalized (min=1 eV)", "Location","best");

%% ---- Optional: print key numbers at OMP ----
Te_omp_min_orig = min(Te_omp_orig, [], "omitnan");
Te_omp_min_new  = min(Te_omp_new,  [], "omitnan");
Te_omp_mean_orig = mean(Te_omp_orig, "omitnan");
Te_omp_mean_new  = mean(Te_omp_new,  "omitnan");

fprintf("\nOMP slice stats (z ≈ %.6g):\n", z(iz_omp));
fprintf("  Original: min = %.6g eV, mean = %.6g eV\n", Te_omp_min_orig, Te_omp_mean_orig);
fprintf("  New:      min = %.6g eV, mean = %.6g eV\n", Te_omp_min_new,  Te_omp_mean_new);