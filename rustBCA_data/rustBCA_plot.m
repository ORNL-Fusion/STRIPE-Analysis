% rustBCA_plot.m
% Visualize sputtering & reflection yields from Ta2O5_results.nc
% Robust to older var names (E/A/spyld) and new ones (energy_eV/angle_deg/Yields)

close all; clear;

% --------- Settings ---------
ncfile = 'Ta2O5_results.nc';   % change if needed
use_log_colors = true;         % log color scale for heatmaps
energies_to_plot = [50,100, 200, 500, 1000, 5000];  % eV
angles_to_plot   = [0, 10, 30, 50, 60, 70, 80];         % degrees

% --------- Read coordinates & arrays (with fallbacks) ---------
E      = read_var(ncfile, {'energy_eV','E'});          % nE
A      = read_var(ncfile, {'angle_deg','A'});          % nA
Yields = read_var(ncfile, {'Yields','spyld'});         % typically [nE x nA]

% Reflection yield is optional
try
    Yref = read_var(ncfile, {'Yields_ref','rfyld'});
catch
    Yref = [];
end

% Ensure column vectors
E = E(:); A = A(:);

% Ensure Yields is 2D and oriented as [nA x nE] for plotting vs (A,E)
sz = size(Yields);
if numel(sz) ~= 2
    Yields = squeeze(Yields);
    sz = size(Yields);
end
if isequal(sz, [numel(E) numel(A)])
    spyld = Yields.';     % [nA x nE]
elseif isequal(sz, [numel(A) numel(E)])
    spyld = Yields;       % already [nA x nE]
else
    error('Unexpected Yields size %s for E=%d, A=%d', mat2str(sz), numel(E), numel(A));
end

% Reflection yield (optional)
rfyld = [];
if ~isempty(Yref)
    szr = size(Yref);
    if numel(szr) ~= 2
        Yref = squeeze(Yref);
        szr = size(Yref);
    end
    if isequal(szr, [numel(E) numel(A)])
        rfyld = Yref.'; 
    elseif isequal(szr, [numel(A) numel(E)])
        rfyld = Yref;
    else
        warning('Unexpected Yields_ref size %s; skipping reflection heatmap.', mat2str(szr));
    end
end

% Build plotting grid
[AA, EE] = ndgrid(A, E);  % AA,EE: [nA x nE]
if ~isequal(size(AA), size(EE), size(spyld))
    error('Grid/data size mismatch: size(AA)=%s, size(EE)=%s, size(spyld)=%s', ...
          mat2str(size(AA)), mat2str(size(EE)), mat2str(size(spyld)));
end

% --------- Heatmap: Sputtering Yield ---------
figure('Name','Sputtering Yield (Ta2O5)','Color','w');
p = pcolor(EE, AA, spyld);
set(p,'EdgeColor','none');
set(gca,'XScale','log');
if use_log_colors, set(gca,'ColorScale','log'); end
xlabel('Energy [eV]'); ylabel('Angle [deg]');
title({'Sputtering Yield','Projectile on Ta_2O_5'});
colorbar;

% --------- Heatmap: Reflection Yield (if present) ---------
if ~isempty(rfyld)
    figure('Name','Reflection Yield (Ta2O5)','Color','w');
    p2 = pcolor(EE, AA, rfyld);
    set(p2,'EdgeColor','none');
    set(gca,'XScale','log');
    if use_log_colors, set(gca,'ColorScale','log'); end
    xlabel('Energy [eV]'); ylabel('Angle [deg]');
    title({'Reflection Yield','Projectile on Ta_2O_5'});
    colorbar;
end

% --------- Smooth slices using griddedInterpolant ---------
F = griddedInterpolant({A, E}, spyld, 'linear', 'none');

energy_dist = linspace(max(min(E(E>0)), 1e-3), max(E), 10000);
angle_dist  = linspace(min(A), max(A),  90);
[AAq,EEq]   = ndgrid(angle_dist, energy_dist);
Yq          = F(AAq, EEq);  % interpolated yield

% Yield vs Angle (for several energies)
[~, idxE] = arrayfun(@(v) min(abs(energy_dist - v)), energies_to_plot);
figure('Name','Yield vs Angle (Ta2O5)','Color','w'); hold on; grid on; box on;
for k = 1:numel(idxE)
    plot(angle_dist, Yq(:, idxE(k)), 'LineWidth', 1.5, ...
        'DisplayName', sprintf('%g eV', energies_to_plot(k)));
end
xlabel('Incidence angle (deg)'); ylabel('Sputtering yield (atoms/ion)');
title('Yield vs Angle for Different Energies');
legend('Location','northwest');

% Yield vs Energy (for several angles)
[~, idxA] = arrayfun(@(v) min(abs(angle_dist - v)), angles_to_plot);
figure('Name','Yield vs Energy (Ta2O5)','Color','w'); hold on; grid on; box on;
for k = 1:numel(idxA)
    plot(energy_dist, Yq(idxA(k),:), 'LineWidth', 1.5, ...
        'DisplayName', sprintf('%g°', angles_to_plot(k)));
end
set(gca,'XScale','log');
xlabel('Ion energy (eV)'); ylabel('Sputtering yield (atoms/ion)');
title('Yield vs Energy for Different Angles');
legend('Location','northwest');

% ==========================
% Local helper (for scripts)
% ==========================
function data = read_var(fname, candidates)
    % Try a list of candidate variable names; return the first found.
    lastErr = [];
    for k = 1:numel(candidates)
        try
            data = ncread(fname, candidates{k});
            return;
        catch err
            lastErr = err; %#ok<NASGU>
        end
    end
    error('None of the variables exist in %s: %s', fname, strjoin(candidates, ', '));
end