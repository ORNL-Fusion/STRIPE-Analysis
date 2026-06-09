% view_sputtering_yields.m
% Quick viewer for surface_model_GITRm_rustbca_C_W_d.nc
% Reads: Angles (40), Energies (50), Projectiles (2), Targets (2)
% Variables: Physical_Sputtering (E x A x target x projectile)
%           Reflection_Yield, Energy_Reflection_Yield
%
% Usage: edit the selection variables below, then run the script.

clearvars; close all; clc;

% -------- User selections (change these) ----------
file = 'surface_model_GITRm_rustbca_C_W_d.nc';  % path to file

% Choose indices (1-based) or set to [] to prompt
proj_idx   = [];   % e.g. 1 or 2; set [] to prompt
target_idx = [];   % e.g. 1 or 2; set [] to prompt

% Plotting options
do_log_color = true;   % use log color scale for the map (useful for yields)
interp_fine = false;   % set true to interpolate onto a finer grid for plots
energy_fine_pts = 300; % only used if interp_fine = true
angle_fine_pts  = 300;
% --------------------------------------------------

% --- Read file ---
info = ncinfo(file);
% display brief summary
fprintf('File: %s\n', file);
fprintf('Variables present: ');
disp({info.Variables.Name});

Angles    = ncread(file, 'Angles');      % [angle_dim x 1]
Energies  = ncread(file, 'Energies');    % [energy_dim x 1]
Proj_raw  = ncread(file, 'Projectiles'); % numeric labels
Targ_raw  = ncread(file, 'Targets');     % numeric labels

% read yields (sizes: energy x angle x target x projectile)
PhysSput  = ncread(file, 'Physical_Sputtering');      % [nE x nA x nT x nP]
ReflY     = ncread(file, 'Reflection_Yield');         % same dims
EReflY    = ncread(file, 'Energy_Reflection_Yield');  % same dims

nE = numel(Energies);
nA = numel(Angles);
nP = numel(Proj_raw);
nT = numel(Targ_raw);

fprintf('Grid sizes: nE=%d, nA=%d, nProj=%d, nTarg=%d\n', nE, nA, nP, nT);

% Prompt selection if needed
if isempty(proj_idx)
    fprintf('Available projectiles (raw values):\n');
    disp(Proj_raw(:)');
    proj_idx = input(sprintf('Choose projectile index (1..%d): ', nP));
end
if isempty(target_idx)
    fprintf('Available targets (raw values):\n');
    disp(Targ_raw(:)');
    target_idx = input(sprintf('Choose target    index (1..%d): ', nT));
end

% Safety checks
if proj_idx < 1 || proj_idx > nP || target_idx < 1 || target_idx > nT
    error('proj_idx or target_idx out of range.');
end

% Extract 2D yield matrix for chosen combo
% Note: stored as (energy x angle) => size [nE x nA]
Y_phys = squeeze(PhysSput(:,:, target_idx, proj_idx));
Y_refl = squeeze(ReflY(:,:,   target_idx, proj_idx));
Y_erefl= squeeze(EReflY(:,:,  target_idx, proj_idx));

% Confirm orientation
assert(isequal(size(Y_phys), [nE, nA]), 'Unexpected matrix shape for yields');

% For plotting with pcolor / imagesc we want matrix as [nA x nE] (rows=angle, cols=energy)
Z_phys = Y_phys.';    % now size [nA x nE]
Z_refl = Y_refl.';    % same
Z_erefl= Y_erefl.';

% Optional interpolation to fine grid for nicer maps
if interp_fine
    E_fine = linspace(min(Energies), max(Energies), energy_fine_pts);
    A_fine = linspace(min(Angles),   max(Angles),   angle_fine_pts);
    % interp2 expects Z sized [length(Angles) x length(Energies)] => Z_phys is correct
    [E_q, A_q] = meshgrid(E_fine, A_fine);  % E_q columns vary over energies, rows over angles
    Z_phys_fine = interp2(Energies, Angles, Z_phys, E_q, A_q, 'linear', 0);
    Z_refl_fine = interp2(Energies, Angles, Z_refl, E_q, A_q, 'linear', 0);
    Z_erefl_fine= interp2(Energies, Angles, Z_erefl,E_q,A_q, 'linear', 0);

    E_plot = E_fine;
    A_plot = A_fine;
    Zp = Z_phys_fine;
    Zr = Z_refl_fine;
else
    E_plot = Energies;
    A_plot = Angles;
    Zp = Z_phys;
    Zr = Z_refl;
end

% Labels for legend from raw arrays (convert numeric labels to strings if needed)
proj_label = num2str(Proj_raw(proj_idx));
targ_label = num2str(Targ_raw(target_idx));

% ----------------- Figure 1 : Physical sputtering map ------------------
figure('Name','Physical Sputtering (map)','NumberTitle','off','Units','normalized','Position',[0.1 0.1 0.7 0.6]);
h = pcolor(E_plot, A_plot, Zp); 
h.EdgeColor = 'none';
xlabel('Energy [eV]');
ylabel('Incidence angle [deg]');
title(sprintf('Physical sputtering yield: target=%s, projectile=%s', targ_label, proj_label));
if do_log_color
    set(gca,'ColorScale','log');
end
set(gca,'XScale','log');   % often useful because energies span decades
cb = colorbar;
cb.Label.String = 'Yield (atoms/ion)';
colormap(parula);
axis tight; box on;

% ----------------- Figure 2 : Reflection map (optional) --------------
figure('Name','Reflection Yield (map)','NumberTitle','off');
h2 = pcolor(E_plot, A_plot, Zr);
h2.EdgeColor = 'none';
xlabel('Energy [eV]'); ylabel('Angle [deg]');
title(sprintf('Reflection yield: target=%s, projectile=%s', targ_label, proj_label));
if do_log_color
    set(gca,'ColorScale','log');
end
set(gca,'XScale','log');
colorbar; colormap(parula); axis tight; box on;

% -------------- Figure 3 : Line-cuts (Yield vs Energy at chosen angles) ------------
angles_to_plot = [0, 15, 30, 45, 60, 75];   % degrees; will choose nearest existing angles
% Find nearest angle indices in A_plot
[~, ai] = arrayfun(@(a) min(abs(A_plot - a)), angles_to_plot);
ai = unique(ai);

figure('Name','Yield vs Energy at selected angles','NumberTitle','off');
hold on;
for k = 1:numel(ai)
    plot(E_plot, Zp(ai(k),:), 'LineWidth', 1.6);
end
set(gca,'XScale','log');
xlabel('Energy [eV]'); ylabel('Yield (atoms/ion)');
legend(compose('Angle = %.1f°', A_plot(ai)));
title('Physical sputtering: yield vs energy for selected incidence angles');
grid on; box on; hold off;

% -------------- Figure 4 : Line-cuts (Yield vs Angle at chosen energies) ------------
energies_to_plot = [20, 100, 500, 1000]; % eV
[~, ei] = arrayfun(@(ev) min(abs(E_plot - ev)), energies_to_plot);
ei = unique(ei);

figure('Name','Yield vs Angle at selected energies','NumberTitle','off');
hold on;
for k = 1:numel(ei)
    plot(A_plot, Zp(:, ei(k)), 'LineWidth', 1.6);
end
xlabel('Incidence angle [deg]'); ylabel('Yield (atoms/ion)');
legend(compose('Energy = %.1f eV', E_plot(ei)));
title('Physical sputtering: yield vs angle for selected energies');
grid on; box on; hold off;

% -------------- Save slices or data optionally ----------------
save_slices = false;
if save_slices
    out_pref = sprintf('yields_target%s_proj%s', targ_label, proj_label);
    writematrix(Zp, [out_pref '_physical.csv']);
    writematrix(Zr, [out_pref '_reflection.csv']);
    fprintf('Saved CSV slices with prefix %s_\n', out_pref);
end

% -------------- End --------------
fprintf('Done. Use the top-of-script variables to change projectile/target or plotting options.\n');