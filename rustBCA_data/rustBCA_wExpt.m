%% RustBCA + EXPT overlay (AUTO schema detect: old E/A/spyld OR new energy_eV/angle_deg/Yields)
% For D/W and Ne/W:
%  - Expt data exist mainly at 0° => overlays appear on the energy sweep at 0°
%  - Script ALWAYS makes a dedicated 0° comparison figure

close all; clearvars;

%% ------------------------ Inputs --------------------------
file     = 'RustBCA_DonW.nc';     % or 'RustBCA_NeonW_80keV.nc'
% file   = 'RustBCA_NeonW_80keV.nc';

expt_csv = 'sputtering_expt_only_IPP9_82_W.csv';

Ecap_eV = [];      % optional cap, e.g. 8e4; [] = use full RustBCA max energy
use_log_colorscale = false;

% Lineout controls
energies_to_plot = [20 70 100 200 500 1e3 5e3 1e4 2e4 4e4 8e4];
angles_to_plot   = [0 5 10 30 50 60 70 80 85 88];     % includes 0°

% Grids for smooth curves
nE_dense = 2500;
nA_dense = 181;

% Expt overlay tolerance for "0°"
expt_ang_tol_deg = 0.5;

%% ------------------------ Infer projectile/target from filename ------------------------
[proj, targ, ttl] = inferPairFromFilename(file);
if proj == "" || targ == ""
    error('Cannot infer projectile/target from filename. Rename like RustBCA_DonW.nc or RustBCA_NeonW_80keV.nc');
end

%% ------------------------ AUTO-detect schema & read RustBCA ------------------------
info = ncinfo(file);
vars = string({info.Variables.Name});

hasOld = all(ismember(["E","A","spyld"], vars));
hasNew = all(ismember(["energy_eV","angle_deg","Yields"], vars));

if hasOld
    fprintf('Detected OLD schema (E, A, spyld)\n');
    E = double(ncread(file,'E')); E = E(:);
    A = double(ncread(file,'A')); A = A(:);
    Y = double(ncread(file,'spyld'));   % [nA x nE]
elseif hasNew
    fprintf('Detected NEW schema (energy_eV, angle_deg, Yields)\n');
    E = double(ncread(file,'energy_eV')); E = E(:);
    A = double(ncread(file,'angle_deg')); A = A(:);
    Y = double(ncread(file,'Yields'));    % [nA x nE] (hopefully)
    if isequal(size(Y), [numel(E) numel(A)])
        Y = Y.'; % transpose if [nE x nA]
    end
else
    fprintf('\nVariables in %s:\n', file);
    disp(vars.');
    error('Unknown schema. Need either (E,A,spyld) or (energy_eV,angle_deg,Yields).');
end

% Ensure orientation [angle x energy]
if isequal(size(Y), [numel(E) numel(A)])
    Y = Y.';
end
assert(isequal(size(Y), [numel(A) numel(E)]), 'Yield array must be [numel(A) x numel(E)].');

% Sort monotonic
[E, iE] = sort(E,'ascend');
[A, iA] = sort(A,'ascend');
Y = Y(iA, iE);

% Energy range (auto max + optional cap)
Emax_data = max(E(isfinite(E)));
if isempty(Ecap_eV), Emax_plot = Emax_data; else, Emax_plot = min(Ecap_eV, Emax_data); end

maskE = (E > 0) & (E <= Emax_plot);
Euse = E(maskE);
Yuse = Y(:, maskE);

if isempty(Euse), error('No energies found in (0, Emax_plot].'); end
fprintf('Using RustBCA E range: [%.3g, %.3g] eV\n', min(Euse), max(Euse));

%% ------------------------ Load EXPT-only CSV and filter to pair ------------------------
if exist(expt_csv,'file') ~= 2
    error('Experimental CSV not found: %s', expt_csv);
end

Texpt = readtable(expt_csv);
if iscell(Texpt.projectile), Texpt.projectile = string(Texpt.projectile); end
if iscell(Texpt.target),     Texpt.target     = string(Texpt.target);     end

if any(strcmpi(Texpt.Properties.VariableNames,'E_keV'))
    Texpt.E_eV = 1e3 * Texpt.E_keV;
elseif any(strcmpi(Texpt.Properties.VariableNames,'E_eV'))
    % ok
else
    error('CSV missing E_keV or E_eV column.');
end

Texpt = Texpt(Texpt.projectile==proj & Texpt.target==targ, :);
Texpt = Texpt(isfinite(Texpt.E_eV) & Texpt.E_eV>0 & isfinite(Texpt.angle_deg) & isfinite(Texpt.yield_atoms_per_ion), :);
Texpt = Texpt(Texpt.E_eV <= Emax_plot, :);

%% ------------------------ Interpolant ------------------------
Ysim_at = @(ang, en) interpn(A, Euse, Yuse, ang, en, 'linear', 0);

%% ------------------------ Yield map ------------------------
figure('Color','w');
h = pcolor(Euse, A, Yuse);
set(h,'EdgeColor','none');
set(gca,'XScale','log');
if use_log_colorscale, set(gca,'ColorScale','log'); end
xlabel('Ion energy (eV)');
ylabel('Incidence angle (deg)');
title(sprintf('Sputtering Yield (%s) - RustBCA', ttl));
colorbar; grid on; box on;

%% ------------------------ Grids for smooth lineouts ------------------------
Emin = max(min(Euse), 1);
Emax = max(Euse);
Egrid = logspace(log10(Emin), log10(Emax), nE_dense);
Agrid = linspace(min(A), max(A), nA_dense);

energies_to_plot = energies_to_plot(energies_to_plot>=Emin & energies_to_plot<=Emax);
angles_to_plot   = angles_to_plot(angles_to_plot>=min(A) & angles_to_plot<=max(A));

%% ------------------------ Two-panel lineouts ------------------------
figure('Color','w');

% ---- Panel 1: Y vs Angle at selected energies ----
subplot(2,1,1); hold on;
colsE = jet(max(1,numel(energies_to_plot)));
for k = 1:numel(energies_to_plot)
    Ek = energies_to_plot(k);
    yk = Ysim_at(Agrid, Ek*ones(size(Agrid)));
    plot(Agrid, yk, 'LineWidth', 1.6, 'Color', colsE(k,:));
end
set(gca,'YScale','log');
xlabel('Incidence angle (deg)');
ylabel('Yield (atoms/ion)');
title('Yield vs Angle (RustBCA)');
legend(compose('%g eV', energies_to_plot),'Location','northwest');
grid on; box on;

% ---- Panel 2: Y vs Energy at selected angles (includes 0°) ----
subplot(2,1,2); hold on;
colsA = jet(max(1,numel(angles_to_plot)));
for k = 1:numel(angles_to_plot)
    Ak = angles_to_plot(k);
    yk = Ysim_at(Ak*ones(size(Egrid)), Egrid);
    yk = movmean(yk, 9);
    plot(Egrid, yk, 'LineWidth', 1.6, 'Color', colsA(k,:));
end

% Overlay EXPT ~0° points (D/Ne/He tables are at 0°)
if ~isempty(Texpt)
    m0 = abs(Texpt.angle_deg - 0) <= expt_ang_tol_deg;
    if any(m0)
        [Eex, is] = sort(Texpt.E_eV(m0));
        Yex = Texpt.yield_atoms_per_ion(m0); Yex = Yex(is);
        scatter(Eex, Yex, 46, 'k', 'filled', 'MarkerFaceAlpha', 0.65);
        legend([compose('%g°', angles_to_plot), "EXPT ~0°"], 'Location','northwest');
    else
        legend(compose('%g°', angles_to_plot),'Location','northwest');
    end
else
    legend(compose('%g°', angles_to_plot),'Location','northwest');
end

set(gca,'XScale','log'); set(gca,'YScale','log');
xlim([Emin Emax]);
xlabel('Ion energy (eV)');
ylabel('Yield (atoms/ion)');
title('Yield vs Energy (RustBCA; includes EXPT at ~0°)');
grid on; box on;

%% ------------------------ Dedicated 0° comparison figure ------------------------
figure('Color','w'); hold on;
Ysim0 = Ysim_at(0*ones(size(Egrid)), Egrid);
plot(Egrid, Ysim0, 'LineWidth', 2.4);

if ~isempty(Texpt)
    m0 = abs(Texpt.angle_deg - 0) <= expt_ang_tol_deg;
    if any(m0)
        [Eex, is] = sort(Texpt.E_eV(m0));
        Yex = Texpt.yield_atoms_per_ion(m0); Yex = Yex(is);
        scatter(Eex, Yex, 55, 'k', 'filled', 'MarkerFaceAlpha', 0.7);
        legend('RustBCA @ 0°','EXPT ~0°','Location','northwest');
    else
        legend('RustBCA @ 0°','Location','northwest');
    end
else
    legend('RustBCA @ 0°','Location','northwest');
end

set(gca,'XScale','log'); set(gca,'YScale','log');
xlim([Emin Emax]);
xlabel('Ion energy (eV)');
ylabel('Yield (atoms/ion)');
title(sprintf('0° comparison: %s', ttl));
grid on; box on;

%% ------------------------ helpers ------------------------
function [proj, targ, ttl] = inferPairFromFilename(file)
    fname = lower(file);
    if contains(fname,'donw')
        proj = "D";  targ = "W"; ttl = "D on W";
    elseif contains(fname,'neonw')
        proj = "Ne"; targ = "W"; ttl = "Ne on W";
    elseif contains(fname,'heonw')
        proj = "He"; targ = "W"; ttl = "He on W";
    elseif contains(fname,'wonw')
        proj = "W";  targ = "W"; ttl = "W on W";
    else
        proj = ""; targ = ""; ttl = "Unknown";
    end
end