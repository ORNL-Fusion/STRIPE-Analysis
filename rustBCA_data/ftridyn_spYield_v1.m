%% RustBCA: read ftridyn schema (E, A, spyld/rfyld) and plot/interp
% Revised for file with variables:
%   E (energy, eV), A (angle, deg), spyld/rfyld (size nA x nE)

close all; clear;

% ------------------------ Input file --------------------------
file = '/Users/78k/ORNL Dropbox/Atul Kumar/work/STRIPE-Analysis/rustBCA_data/ftridynSelf_new_higResol.nc';
% file = '/Users/78k/ORNL Dropbox/Atul Kumar/work/STRIPE-Analysis/rustBCA_data/ftridyn_DonW.nc';


% Choose yield type: 'spyld', 'rfyld', or 'total'
yield_mode = 'spyld';

% ------------------------ Read grids --------------------------
energy = double(ncread(file, 'E'));   % eV
angle  = double(ncread(file, 'A'));   % deg
energy = energy(:);
angle  = angle(:);

% ------------------------ Read yield --------------------------
switch lower(yield_mode)
    case 'spyld'
        Y = double(ncread(file, 'spyld'));
        yLabel = 'Sputtering yield (spyld)';
    case 'rfyld'
        Y = double(ncread(file, 'rfyld'));
        yLabel = 'Reflection yield (rfyld)';
    case 'total'
        Y = double(ncread(file, 'spyld')) + double(ncread(file, 'rfyld'));
        yLabel = 'Total yield (spyld + rfyld)';
    otherwise
        error('yield_mode must be ''spyld'', ''rfyld'', or ''total''.');
end

% Ensure Y is [angle x energy] = [nA x nE]
if isequal(size(Y), [numel(energy), numel(angle)])
    Y = Y.';
end
assert(isequal(size(Y), [numel(angle), numel(energy)]), ...
    'Unexpected yield size. Expected [numel(angle) x numel(energy)].');

% Sort monotonic for interpolation
[energy, iE] = sort(energy, 'ascend');
[angle,  iA] = sort(angle,  'ascend');
Y = Y(iA, iE);

% ---------------------- Limit to <= 80 keV ----------------------
Emax_plot = 8e4; % eV
maskE = (energy > 0) & (energy <= Emax_plot);
energy80 = energy(maskE);
Y80 = Y(:, maskE);

if isempty(energy80)
    error('No energies found in (0, 80 keV]. Check E units/values.');
end

% -------------------------- Surface plot --------------------------
figure('Color', 'w');
h = pcolor(energy80, angle, Y80); % Z is [angle x energy]
set(h, 'EdgeColor', 'none');
set(gca, 'XScale', 'log');
xlabel('Ion energy (eV)');
ylabel('Incidence angle (deg)');
title([yLabel ' (Ne on W)']);
colorbar; box on; grid on;

if exist('coolwarm.mat','file')
    S = load('coolwarm.mat', 'coolwarm_rgb');
    if isfield(S, 'coolwarm_rgb')
        colormap(S.coolwarm_rgb);
    end
end

% --------------------- Dense interpolation grids ---------------------
nE_dense = 80000;
nA_dense = 180;

Emin = max(min(energy80), 1); % avoid log(<=0)
Emax = max(energy80);

energy_dist = logspace(log10(Emin), log10(Emax), nE_dense);
angle_dist  = linspace(min(angle), max(angle), nA_dense);

[AA, EE] = ndgrid(angle_dist, energy_dist);
Y_interp = interpn(angle, energy80, Y80, AA, EE, 'linear', 0);
Y_interp(~isfinite(Y_interp)) = 0;

% --------------------- Angle sweeps at chosen energies ---------------------
energies_to_plot = [20 70 100 500 1e3 5e3 1e4 2e4 4e4 8e4];
energies_to_plot = energies_to_plot( ...
    energies_to_plot >= min(energy80) & energies_to_plot <= max(energy80));

figure('Color', 'w');

subplot(2,1,1); hold on;
cols = jet(numel(energies_to_plot));
for k = 1:numel(energies_to_plot)
    yk = interpn(angle, energy80, Y80, angle_dist, ...
        energies_to_plot(k) * ones(size(angle_dist)), 'linear', 0);
    plot(angle_dist, yk, 'LineWidth', 1.6, 'Color', cols(k,:));
end
xlabel('Incidence angle (deg)');
ylabel(yLabel);
title('Yield vs Angle at selected energies (<= 80 keV)');
legend(compose('%.3g eV', energies_to_plot), 'Location', 'northwest');
grid on; box on;

% --------------------- Energy sweeps at chosen angles ---------------------
angles_to_plot = [5 10 30 50 60 70 80 85 88];
angles_to_plot = angles_to_plot( ...
    angles_to_plot >= min(angle) & angles_to_plot <= max(angle));

subplot(2,1,2); hold on;
cols2 = jet(numel(angles_to_plot));
for k = 1:numel(angles_to_plot)
    yk = interpn(angle, energy80, Y80, ...
        angles_to_plot(k) * ones(size(energy_dist)), energy_dist, 'linear', 0);
    yk = movmean(yk, 9);
    plot(energy_dist, yk, 'LineWidth', 1.6, 'Color', cols2(k,:));
end
set(gca, 'XScale', 'log');
xlim([Emin Emax]);
xlabel('Ion energy (eV)');
ylabel(yLabel);
title('Yield vs Energy at selected angles (log energy axis)');
legend(compose('%g°', angles_to_plot), 'Location', 'northwest');
grid on; box on;
