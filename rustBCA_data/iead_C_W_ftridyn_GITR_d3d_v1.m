%% RustBCA: read 2025 schema (energy_eV, angle_deg, Yields) and plot/interp
% Revised:
%  - explicitly supports up to 80 keV (8e4 eV) and clips/plots to that range
%  - uses LOG-SPACED energy grids for interpolation + plotting
%  - uses log x-scale everywhere energy appears
%  - keeps Y (yield) on linear scale unless you enable log color scaling

close all; clear;

% ------------------------ Input file --------------------------
% file = 'RustBCA_NeonW_80keV.nc';
file = 'ftridynSelf_new_higResol.nc';

% file = 'RustBCA_DonW.nc';


% Dimensions & grids
energy = double(ncread(file,'energy_eV'));   % eV, [nE x 1] (or [1 x nE])
angle  = double(ncread(file,'angle_deg'));   % deg, [nA x 1] (or [1 x nA])

energy = energy(:);
angle  = angle(:);

% Core arrays
Y = double(ncread(file,'Yields'));           % expected [nA x nE] = [angle x energy]

% ---- Ensure Y orientation is [angle x energy] ----
if isequal(size(Y), [numel(energy) numel(angle)])
    Y = Y.';  % transpose if stored as [energy x angle]
end

assert(isequal(size(Y), [numel(angle) numel(energy)]), ...
    'Unexpected Yields size. Expected [numel(angle) x numel(energy)].');

% ---- Sort grids if needed (interpn requires monotonic vectors) ----
[energy, iE] = sort(energy, 'ascend');
[angle,  iA] = sort(angle,  'ascend');
Y = Y(iA, iE);

% ---------------------- Limit to <= 80 keV ----------------------
Emax_plot = 8e4; % 80 keV in eV
maskE = (energy > 0) & (energy <= Emax_plot);
energy80 = energy(maskE);
Y80 = Y(:, maskE);

if isempty(energy80)
    error('No energies found in (0, 80 keV]. Check energy_eV units/values.');
end

% ------------------------------- Surface plot (Yield) -------------------------------
figure('Color','w');
h = pcolor(energy80, angle, Y80);   % Z is [angle x energy]
set(h,'EdgeColor','none');
set(gca,'XScale','log');           % log energy axis
% set(gca,'ColorScale','log');      % OPTIONAL: log color scaling (useful if Y spans decades)
xlabel('Ion energy (eV)');
ylabel('Incidence angle (deg)');
title('Sputtering Yield (Ne on W)'); % update label as appropriate
colorbar; box on; grid on;

% Colormap (if available)
if exist('coolwarm.mat','file')
    S = load('coolwarm.mat', 'coolwarm_rgb');
    if isfield(S,'coolwarm_rgb')
        colormap(S.coolwarm_rgb)
    end
end

% ------------------- Interpolate Y(angle,energy) to custom grids --------------------
% Use LOG-SPACED energy sampling for better resolution across decades
nE_dense = 80000;
nA_dense = 180;

Emin = max(min(energy80), 1);      % avoid <=0 for logspace
Emax = max(energy80);

energy_dist = logspace(log10(Emin), log10(Emax), nE_dense);
angle_dist  = linspace(min(angle),  max(angle),  nA_dense);

[AA, EE] = ndgrid(angle_dist, energy_dist);
Y_interp = interpn(angle, energy80, Y80, AA, EE, 'linear', 0); % 0 outside range
Y_interp(~isfinite(Y_interp)) = 0;

% -------------------------- Angle sweeps at chosen energies -------------------------
% energies in eV (include up to 80 keV). Values outside range get dropped.
energies_to_plot = [20 70 100 500 1e3 5e3 1e4 2e4 4e4 8e4];
energies_to_plot = energies_to_plot(energies_to_plot>=min(energy80) & energies_to_plot<=max(energy80));

figure('Color','w');

subplot(2,1,1); hold on;
cols = jet(numel(energies_to_plot));
for k = 1:numel(energies_to_plot)
    yk = interpn(angle, energy80, Y80, angle_dist, energies_to_plot(k)*ones(size(angle_dist)), 'linear', 0);
    plot(angle_dist, yk, 'LineWidth', 1.6, 'Color', cols(k,:));
end
xlabel('Incidence angle (deg)');
ylabel('Sputtering yield (atoms/ion)');
title('Yield vs Angle at selected energies (≤80 keV)');
legend(compose('%.3g eV', energies_to_plot),'Location','northwest'); grid on; box on;

% -------------------------- Energy sweeps at chosen angles --------------------------
angles_to_plot = [5 10 30 50 60 70 80 85 88];
angles_to_plot = angles_to_plot(angles_to_plot>=min(angle) & angles_to_plot<=max(angle));

subplot(2,1,2); hold on;
cols2 = jet(numel(angles_to_plot));
for k = 1:numel(angles_to_plot)
    % Evaluate on LOG-SPACED energy grid and plot vs energy_dist (log x-axis)
    yk = interpn(angle, energy80, Y80, angles_to_plot(k)*ones(size(energy_dist)), energy_dist, 'linear', 0);
    yk = movmean(yk, 9); % light smoothing in log-sampled space
    plot(energy_dist, yk, 'LineWidth', 1.6, 'Color', cols2(k,:));
end
set(gca,'XScale','log');
xlim([Emin Emax]);
xlabel('Ion energy (eV)');
ylabel('Sputtering yield (atoms/ion)');
title('Yield vs Energy at selected angles (energy on log scale)');
legend(compose('%g°', angles_to_plot),'Location','northwest'); grid on; box on;

% ===================== Helpers (index by nearest physical value) ====================
% idxAngle  = @(deg) max(1, min(numel(angle),  round(interp1(angle,  1:numel(angle),  deg, 'nearest','extrap'))));
% idxEnergy = @(eV)  max(1, min(numel(energy80),round(interp1(energy80,1:numel(energy80),eV,  'nearest','extrap'))));