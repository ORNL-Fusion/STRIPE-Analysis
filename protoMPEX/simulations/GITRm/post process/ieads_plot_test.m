
clc; clear; close all;

%% ----------------------------
% USER SETTINGS
%% ----------------------------
tilt = 85;
file = sprintf('../tilted_targets/test/%d_degrees/MPEX_runs_v1/gitrm-surface.nc', tilt);
% file = sprintf('../tilted_targets/test/%d_degrees/MPEX_runs_test/energyAngleImpact.nc', tilt);

% Energy / angle definitions (your stated binning)
N_E_BINS      = 5000;
E_MAX_eV      = 500;
N_ANGLE_BINS  = 180;
ANGLE_MAX_DEG = 90;

% Surface selection
SKIP_FIRST_NSURFACES = 6;   % skip first N surfaces
DROP_LAST_NSURFACES  = 3;   % drop last M surfaces

% % Species / charge indices
species_idx = 1;
charge_idx  = 1;

% OPTIONAL: physical normalization (ONLY if you know these)
A_surf_m2 = [];   % e.g. 1.0e-6  (total area of the selected surface patches), leave [] if unknown
T_phys_s  = [];   % e.g. 1.0e-3  (physical time represented by tally), leave [] if unknown

% Plot controls
XMAX_eV = 500;      % x-axis limit for energy
USE_LOG10 = false;  % true to plot log10(A + eps)
%% ----------------------------

%% READ DATA
surfEDist_all = ncread(file,'surfEDist');
% Expected dims: angle x energy x surface x species x charge (or compatible)

%% SURFACE RANGE
nSurfaces  = size(surfEDist_all, 3);
surf_start = SKIP_FIRST_NSURFACES + 1;
surf_end   = nSurfaces - DROP_LAST_NSURFACES;

if surf_start > surf_end
    error('Surface selection invalid: no surfaces left. Check SKIP/DROP values.');
end
surf_idx = surf_start:surf_end;

%% SUM OVER SELECTED SURFACES (weights)
A = squeeze(sum(surfEDist_all(:,:,surf_idx,species_idx,charge_idx), 3));
% A -> [angle x energy] (ideally)

figure; pcolor(A); shading interp;
xlim([0 1000])

% Fix orientation if needed (ensure columns are energy)
if size(A,2) ~= N_E_BINS && size(A,1) == N_E_BINS
    A = A.';
end

[nAngleFile, nEnergyFile] = size(A);

%% BUILD AXES
% Real energy axis (bin centers)
E_eV = linspace(E_MAX_eV/(2*N_E_BINS), ...
                E_MAX_eV - E_MAX_eV/(2*N_E_BINS), ...
                N_E_BINS);

% Angle axis (0..90 deg)
if nAngleFile ~= N_ANGLE_BINS
    % Fall back to file angle count
    angle_deg = linspace(0, ANGLE_MAX_DEG, nAngleFile);
else
    angle_deg = linspace(0, ANGLE_MAX_DEG, N_ANGLE_BINS);
end

% Resample energy dimension if file doesn't match N_E_BINS
if nEnergyFile ~= N_E_BINS
    x_old = linspace(0,1,nEnergyFile);
    x_new = linspace(0,1,N_E_BINS);
    A = interp1(x_old, A.', x_new, 'linear', 0).';
end

%% BIN WIDTHS
dE = E_MAX_eV / N_E_BINS;               % eV/bin  (here 0.1 eV)
dTheta_deg = ANGLE_MAX_DEG / N_ANGLE_BINS; % deg/bin (here 0.5 deg)

%% ----------------------------
% UNITS / NORMALIZATION NOTES
% A is a WEIGHTED TALLY PER BIN:
%   units(A) = "weighted particles per (E-bin, angle-bin)" summed over chosen surfaces
% If weights represent # of real ions, then A is in "real ions per bin" (still not per time/area).
%% ----------------------------

%% SANITY CHECKS (Jeremy-friendly)
total_weighted_particles = sum(A(:));
fprintf('Selected surfaces: %d to %d (count=%d)\n', surf_start, surf_end, numel(surf_idx));
fprintf('Sum of weights over all (E,theta) bins = %.6e (weighted particles)\n', total_weighted_particles);

% PDF / percent per bin (discrete): sums to 100%
percent_per_bin = 100 * (A / total_weighted_particles);
fprintf('Sum of percent_per_bin over all bins = %.12f %% (should be 100%%)\n', sum(percent_per_bin(:)));

%% OPTIONAL: physical IEAD (only if A_surf_m2 and T_phys_s are provided)
IEAD_phys = [];
if ~isempty(A_surf_m2) && ~isempty(T_phys_s)
    % m^-2 s^-1 eV^-1 deg^-1
    IEAD_phys = A ./ (A_surf_m2 * T_phys_s * dE * dTheta_deg);
end

%% ----------------------------
% PLOTS
%% ----------------------------

% 1) RAW WEIGHTED COUNTS PER BIN
figure('Color','w','Position',[120 120 900 520]);
if USE_LOG10
    Z = log10(A + eps);
    p = pcolor(E_eV, angle_deg, Z);
    cb_label = 'log_{10}(weighted particles / bin)';
else
    p = pcolor(E_eV, angle_deg, A);
    cb_label = 'weighted particles / bin';
end
set(p,'EdgeColor','none');
xlabel('Energy (eV)');
ylabel('Angle (deg)');
title(sprintf('GITRm surfEDist (RAW) | tilt=%d° | surfaces %d–%d', tilt, surf_start, surf_end));
cb = colorbar; ylabel(cb, cb_label);
set(gca,'YDir','normal');
xlim([0 XMAX_eV]);
ylim([0 ANGLE_MAX_DEG]);

% 2) PERCENT OF TOTAL (PER BIN) — sums to 100%
figure('Color','w','Position',[150 150 900 520]);
p = pcolor(E_eV, angle_deg, percent_per_bin);
set(p,'EdgeColor','none');
xlabel('Energy (eV)');
ylabel('Angle (deg)');
title(sprintf('IEAD as PDF: %% of total per bin | tilt=%d° | surfaces %d–%d', tilt, surf_start, surf_end));
cb = colorbar; ylabel(cb, '% of total (per bin)');
set(gca,'YDir','normal');
xlim([0 XMAX_eV]);
ylim([0 ANGLE_MAX_DEG]);

% 3) PHYSICAL IEAD (if enabled)
if ~isempty(IEAD_phys)
    figure('Color','w','Position',[180 180 900 520]);
    if USE_LOG10
        Z = log10(IEAD_phys + eps);
        p = pcolor(E_eV, angle_deg, Z);
        cb_label = 'log_{10}(IEAD) [m^{-2} s^{-1} eV^{-1} deg^{-1}]';
    else
        p = pcolor(E_eV, angle_deg, IEAD_phys);
        cb_label = 'IEAD [m^{-2} s^{-1} eV^{-1} deg^{-1}]';
    end
    set(p,'EdgeColor','none');
    xlabel('Energy (eV)');
    ylabel('Angle (deg)');
    title(sprintf('PHYSICAL IEAD | tilt=%d° | surfaces %d–%d', tilt, surf_start, surf_end));
    cb = colorbar; ylabel(cb, cb_label);
    set(gca,'YDir','normal');
    xlim([0 XMAX_eV]);
    ylim([0 ANGLE_MAX_DEG]);
else
    fprintf('\nPhysical IEAD not computed (A_surf_m2 or T_phys_s not provided).\n');
    fprintf('If you provide them, units will be: m^-2 s^-1 eV^-1 deg^-1\n\n');
end