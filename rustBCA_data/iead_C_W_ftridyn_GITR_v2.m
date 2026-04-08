%% RustBCA + PURE-EXPT (IPP 9/82) : match expt & sim grids for lineouts
% - Works with RustBCA schema: E, A, spyld
% - Uses ONLY expt energies/angles (from CSV) for plotting lineouts
% - For each expt E (angle-sweep) and each expt angle (energy-sweep):
%       sim is interpolated onto the expt grid (no "nearest index")
% - Automatically works for: D/W, He/W, Ne/W, W/W (based on filename OR you can override)
%
% CSV columns:
% projectile,target,E_keV,angle_deg,yield_atoms_per_ion,source
%
% NOTES:
% - D/He/Ne on W: expt data are mainly 0° only => energy-sweep panel will show 0°; angle-sweep will be empty
% - W on W: expt includes angle-dependent digitized points => both panels populate

close all; clearvars;

%% ------------------------ Inputs ------------------------
file     = 'RustBCA_WonW.nc';  % <- change to RustBCA_DonW.nc / RustBCA_HeonW.nc / RustBCA_NeonW_80keV.nc
expt_csv = 'sputtering_expt_digitized_IPP9_82_W.csv';

Emax_plot = 8e4;   % eV (80 keV)
do_log_color = true;   % yield-map colormap scaling

% Energy/angle matching tolerances (only used to group expt points by "same energy" etc.)
% We will still PLOT at the raw expt points; these are only for grouping.
energy_reltol   = 0.06;  % 6% grouping window for digitized expt energies
energy_abstol_eV = 10;   % absolute grouping window at low E (eV)
ang_group_tol   = 0.5;   % deg grouping window (digitized angles)

%% ------------------------ Infer projectile/target from filename (or override) ------------------------
[proj, targ, ttl] = inferPairFromFilename(file);
if proj == "" || targ == ""
    error('Could not infer projectile/target from filename. Rename file like RustBCA_DonW.nc, ... or edit inferPairFromFilename().');
end

%% ------------------------ Read RustBCA netCDF ------------------------
E = double(ncread(file,'E')); E = E(:);      % energy [eV]
A = double(ncread(file,'A')); A = A(:);      % angle [deg]
spyld = double(ncread(file,'spyld'));        % [nA x nE]

if ~isequal(size(spyld), [numel(A) numel(E)])
    error('spyld size is %s but expected [%d %d] = [numel(A) numel(E)].', ...
        mat2str(size(spyld)), numel(A), numel(E));
end

% Sort monotonic
[E, iE] = sort(E,'ascend');
[A, iA] = sort(A,'ascend');
spyld = spyld(iA, iE);

% Energy limit
maskE = (E > 0) & (E <= Emax_plot);
E80 = E(maskE);
spyld80 = spyld(:, maskE);

if isempty(E80)
    error('No energies found in (0, 80 keV]. Check E units/values.');
end

%% ------------------------ Load & filter experimental CSV ------------------------
if exist(expt_csv,'file') ~= 2
    error('Experimental CSV not found: %s', expt_csv);
end

expt = readtable(expt_csv);
if iscell(expt.projectile), expt.projectile = string(expt.projectile); end
if iscell(expt.target),     expt.target     = string(expt.target);     end
if iscell(expt.source),     expt.source     = string(expt.source);     end

% energy to eV
if any(strcmpi(expt.Properties.VariableNames,'E_keV'))
    expt.E_eV = 1e3 * expt.E_keV;
elseif any(strcmpi(expt.Properties.VariableNames,'E_eV'))
    expt.E_eV = expt.E_eV;
else
    error('CSV missing E_keV or E_eV column.');
end

% filter to pair
expt = expt(expt.projectile==proj & expt.target==targ, :);
expt = expt(isfinite(expt.E_eV) & expt.E_eV>0 & isfinite(expt.angle_deg) & isfinite(expt.yield_atoms_per_ion), :);
expt = expt(expt.E_eV <= Emax_plot, :);

if isempty(expt)
    error('No experimental rows found in CSV for %s on %s.', proj, targ);
end

%% ------------------------ Yield map (simulation) ------------------------
figure(1)
h = pcolor(E80, A, spyld80);
h.EdgeColor = 'none';
colorbar
set(gca,'XScale','log');
if do_log_color, set(gca,'ColorScale','log'); end
xlabel('Energy [eV]'); ylabel('Angle [deg]');
title({'Sputtering Yield', sprintf('%s (RustBCA)', ttl)});
grid on; box on;

%% ------------------------ Build interpolation handle for simulation ------------------------
% Interpolate sim yield onto arbitrary expt (A,E) points:
% spyld80 is defined on grid (A, E80)
Ysim_at = @(angles_deg, energies_eV) interpn(A, E80, spyld80, angles_deg, energies_eV, 'linear', 0);

%% ======================== MATCHED LINEOUTS (EXPT GRID) ========================

% -------- Panel A: Yield vs Energy at expt angles (EXPT angles) --------
% Group expt by (approximately) the same angle.
% We'll pick representative unique angles using rounding/grouping.
exptAngles = expt.angle_deg;
angKey = round(exptAngles/ang_group_tol)*ang_group_tol;   % clustered angles
uAng = unique(angKey);
uAng = sort(uAng);

% Keep angles within sim A-range
uAng = uAng(uAng>=min(A) & uAng<=max(A));

% -------- Panel B: Yield vs Angle at expt energies (EXPT energies) --------
exptE = expt.E_eV;
% cluster energies on log scale-ish using relative tolerance:
% Create bins by sorting and accumulating within tolerance.
[uE] = clusterEnergies(exptE, energy_reltol, energy_abstol_eV);
uE = sort(uE);
uE = uE(uE>=min(E80) & uE<=max(E80));

figure(2);

% -------------------- Subplot 1: Y vs Angle at expt energies --------------------
subplot(2,1,1); hold on;

if isempty(uE)
    text(0.1,0.5,'No angle-dependent experimental energies for this pair (likely only 0° table).','Units','normalized');
    axis off;
else
    colsE = jet(numel(uE));
    for k = 1:numel(uE)
        Ek = uE(k);

        % EXPT points near this Ek (grouping only)
        mE = abs(expt.E_eV - Ek) <= max(energy_abstol_eV, energy_reltol*Ek);
        exA = expt.angle_deg(mE);
        exY = expt.yield_atoms_per_ion(mE);

        if isempty(exA), continue; end

        % Sort by angle
        [exA, is] = sort(exA);
        exY = exY(is);

        % SIM evaluated exactly at those expt angles (same Ek)
        simY = Ysim_at(exA, Ek*ones(size(exA)));

        % Plot sim as line through expt angles, expt as markers
        plot(exA, simY, '-', 'LineWidth', 1.6, 'Color', colsE(k,:));
        scatter(exA, exY, 34, 'MarkerEdgeColor', colsE(k,:), 'MarkerFaceColor','none', 'LineWidth', 1.2);
    end

    set(gca,'YScale','log');
    xlabel('Incidence angle [deg]');
    ylabel('Yield [atoms/ion]');
    title('Y vs Angle at experimental energies (line=RustBCA @ expt angles, markers=expt)');
    legend(compose('%g eV',uE), 'Location','NorthWest');
    grid on; box on;
end

% -------------------- Subplot 2: Y vs Energy at expt angles --------------------
subplot(2,1,2); hold on;

if isempty(uAng)
    text(0.1,0.5,'No experimental angles found in-range.','Units','normalized');
    axis off;
else
    colsA = jet(numel(uAng));
    for k = 1:numel(uAng)
        Ak = uAng(k);

        % EXPT points near this angle (grouping only)
        mA = abs(expt.angle_deg - Ak) <= ang_group_tol;
        exE = expt.E_eV(mA);
        exY = expt.yield_atoms_per_ion(mA);

        if isempty(exE), continue; end

        % Sort by energy
        [exE, is] = sort(exE);
        exY = exY(is);

        % SIM evaluated exactly at those expt energies (same Ak)
        simY = Ysim_at(Ak*ones(size(exE)), exE);

        % Plot sim as line vs energy (log x), expt as markers
        plot(exE, simY, '-', 'LineWidth', 1.6, 'Color', colsA(k,:));
        scatter(exE, exY, 34, 'MarkerEdgeColor', colsA(k,:), 'MarkerFaceColor','none', 'LineWidth', 1.2);
    end

    set(gca,'XScale','log'); set(gca,'YScale','log');
    xlim([max(min(expt.E_eV),1) min(max(expt.E_eV),Emax_plot)]);
    xlabel('Ion energy [eV]');
    ylabel('Yield [atoms/ion]');
    title('Y vs Energy at experimental angles (line=RustBCA @ expt energies, markers=expt)');
    legend(compose('%.2g°',uAng), 'Location','NorthWest');
    grid on; box on;
end

%% ------------------------ Quick console report ------------------------
fprintf('Pair: %s on %s | expt points: %d\n', proj, targ, height(expt));
fprintf('Expt angle span: [%.2f, %.2f] deg\n', min(expt.angle_deg), max(expt.angle_deg));
fprintf('Expt energy span: [%.3g, %.3g] eV\n', min(expt.E_eV), max(expt.E_eV));
fprintf('Unique clustered expt angles: %d | unique clustered expt energies: %d\n', numel(uAng), numel(uE));

return

%% ======================== Local helper functions ========================
function [proj, targ, ttl] = inferPairFromFilename(file)
    fname = lower(file);
    if contains(fname,'donw')
        proj = "D";  targ = "W"; ttl = "D on W";
    elseif contains(fname,'heonw')
        proj = "He"; targ = "W"; ttl = "He on W";
    elseif contains(fname,'neonw')
        proj = "Ne"; targ = "W"; ttl = "Ne on W";
    elseif contains(fname,'wonw')
        proj = "W";  targ = "W"; ttl = "W on W";
    else
        proj = ""; targ = ""; ttl = "Unknown";
    end
end

function uE = clusterEnergies(E, reltol, abstol)
    % Cluster energies into representative values:
    % Two energies belong together if |Ei-Ej| <= max(abstol, reltol*Ecenter)
    E = E(:);
    E = E(isfinite(E) & E>0);
    if isempty(E), uE = []; return; end
    E = sort(E,'ascend');
    clusters = {};
    cur = E(1);
    curList = cur;

    for i = 2:numel(E)
        Ei = E(i);
        tol = max(abstol, reltol*cur);
        if abs(Ei - cur) <= tol
            curList(end+1,1) = Ei; %#ok<AGROW>
            cur = median(curList);
        else
            clusters{end+1} = curList; %#ok<AGROW>
            curList = Ei;
            cur = Ei;
        end
    end
    clusters{end+1} = curList;

    uE = zeros(numel(clusters),1);
    for k = 1:numel(clusters)
        uE(k) = median(clusters{k});
    end
end