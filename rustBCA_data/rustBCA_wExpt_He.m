%% He_on_W_full_pipeline_2025schema_FIXED.m
% He+ on W:
% - reads NetCDF schemas with either energy_eV/angle_deg/Yields or E/A/spyld
% - cleans Yields (removes NaN/Inf + clips absurd fill values)
% - overlays EXPT 0° points (your screenshot table)
% - labels "0°" curve as the *minimum RustBCA angle* if 0° is not present
% - optional: fits Eckstein/Yamamura angular factor to RustBCA (normalized at Amin)
%   and anchors modeled yields to EXPT 0° -> writes CSV

close all; clear; clc;

%% ------------------------ Inputs ------------------------
file    = 'ftridyn_HeonW.nc';
out_csv = 'He_on_W_EXPT0_plus_EcksteinAngleFit.csv';

Emax_plot = 8e4;                 % eV max for plots

theta_pts = [0 15 30 45 60 75 80];  % angles (deg) to output in CSV (must include 0)

do_fit_angle_model = true;
fit_theta_max = 80;              % fit only up to this angle
min_sim_yield = 1e-14;           % ignore if Y(Amin,E) below this

use_fmincon_if_available = true;
lb = [0.0, 0.0];                 % bounds for [f,g]
ub = [6.0, 30.0];
p0 = [1.0, 1.0];

%% ------------------------ Read NetCDF schema (robust names/layouts) ------------------------
info = ncinfo(file);
allPaths = strings(0,1);
for iv = 1:numel(info.Variables)
    allPaths(end+1,1) = string(info.Variables(iv).Name); %#ok<AGROW>
end
for ig = 1:numel(info.Groups)
    g = info.Groups(ig);
    gname = string(g.Name);
    for jv = 1:numel(g.Variables)
        allPaths(end+1,1) = "/" + gname + "/" + string(g.Variables(jv).Name); %#ok<AGROW>
    end
end

ePath = "";
energyCandidates = {'energy_eV','E'};
for k = 1:numel(energyCandidates)
    n = string(energyCandidates{k});
    hit = endsWith(allPaths, "/" + n) | (allPaths == n);
    if any(hit), ePath = allPaths(find(hit,1,'first')); break; end
end

aPath = "";
angleCandidates = {'angle_deg','A'};
for k = 1:numel(angleCandidates)
    n = string(angleCandidates{k});
    hit = endsWith(allPaths, "/" + n) | (allPaths == n);
    if any(hit), aPath = allPaths(find(hit,1,'first')); break; end
end

yPath = "";
yieldCandidates = {'Yields','spyld'};
for k = 1:numel(yieldCandidates)
    n = string(yieldCandidates{k});
    hit = endsWith(allPaths, "/" + n) | (allPaths == n);
    if any(hit), yPath = allPaths(find(hit,1,'first')); break; end
end

if strlength(ePath)==0 || strlength(aPath)==0 || strlength(yPath)==0
    error(['Unable to find required variables in %s. ', ...
           'Need energy (energy_eV/E), angle (angle_deg/A), yield (Yields/spyld).'], file);
end

energy = double(ncread(file, ePath)); energy = energy(:);
angle  = double(ncread(file, aPath)); angle  = angle(:);
Yraw   = double(ncread(file, yPath));

% Normalize yield orientation to [angle x energy]
sz = size(Yraw);
if isequal(sz, [numel(angle), numel(energy)])
    Y = Yraw;
elseif isequal(sz, [numel(energy), numel(angle)])
    Y = Yraw.';
else
    Ys = squeeze(Yraw); % handles nS=1 leading dim in ftridyn files
    if isequal(size(Ys), [numel(angle), numel(energy)])
        Y = Ys;
    elseif isequal(size(Ys), [numel(energy), numel(angle)])
        Y = Ys.';
    else
        error('Unexpected yield array size %s. Expected [%d x %d] or [%d x %d].', ...
            mat2str(sz), numel(angle), numel(energy), numel(energy), numel(angle));
    end
end

% Sort monotonic grids
[energy, iE] = sort(energy,'ascend');
[angle,  iA] = sort(angle,'ascend');
Y = Y(iA, iE);

% Limit energy range
maskE = (energy > 0) & isfinite(energy) & (energy <= Emax_plot);
Euse = energy(maskE);
Yuse = Y(:, maskE);

if isempty(Euse)
    error('No energies found in (0, Emax_plot]. Check energy_eV values/units.');
end

Emin = max(min(Euse), 1);
Emax = max(Euse);

% Determine min-angle (effective "0°")
Amin = min(angle);
Amax = max(angle);
has_true_zero = any(abs(angle - 0) < 1e-9);

% Build label for the "0°" curve in legend:
if has_true_zero
    A0_for_sim = 0;
    label0 = "0°";
else
    A0_for_sim = Amin;
    label0 = sprintf('%.1f° (min)', Amin);
end

fprintf('He file: E range [%.3g, %.3g] eV | angle range [%.3g, %.3g] deg\n', Emin, Emax, Amin, Amax);
fprintf('Using %s as the "0°" simulation reference curve.\n', label0);

%% ------------------------ CLEAN yields (CRITICAL) ------------------------
% Remove NaN/Inf
Yuse(~isfinite(Yuse)) = 0;

% Clip negatives
Yuse(Yuse < 0) = 0;

% Clip absurd fill/sentinel values (unphysical for sputtering yields here)
% If your file legitimately has yields > 1e3, increase this.
Yuse(Yuse > 1e3) = 0;

fprintf('Yield map stats after cleaning: min=%g max=%g\n', min(Yuse(:)), max(Yuse(:)));

% Interpolant: Y(angle,energy)
Ysim = @(ang_deg, en_eV) interpn(angle, Euse, Yuse, ang_deg, en_eV, 'linear', 0);

%% ------------------------ Experimental He+ -> W (0°) from your screenshot ------------------------
Eexp_keV = [0.250 0.500 0.500 1.000 2.000 2.000 4.000 8.000 20.000]; % keV
Eexp = 1e3 * Eexp_keV;  % eV

Yexp = [0.00378 0.00883 0.01230 0.02310 0.02876 0.02857 0.03367 0.04272 0.03700]; % atoms/ion

Eexp = Eexp(:); Yexp = Yexp(:);

% Keep only energies in sim range
okE = (Eexp >= Emin*0.999) & (Eexp <= Emax*1.001) & isfinite(Yexp) & (Yexp>=0);
Eexp = Eexp(okE); Yexp = Yexp(okE);
nEexp = numel(Eexp);

%% ------------------------ Plot 1: Yield map (RustBCA) ------------------------
figure('Color','w');
h = pcolor(Euse, angle, Yuse);
set(h,'EdgeColor','none');
set(gca,'XScale','log');
xlabel('Ion energy (eV)');
ylabel('Incidence angle (deg)');
title('He^+ on W: Sputtering Yield map (RustBCA)');
colorbar; grid on; box on;

%% ------------------------ Plot 2: Yield vs Energy (RustBCA + EXPT 0°) ------------------------
figure('Color','w'); hold on;

Egrid = logspace(log10(Emin), log10(Emax), 2000);

A_list_plot = [0 10 30 50 60 70 80 85];
A_list_plot = A_list_plot(A_list_plot>=Amin & A_list_plot<=Amax);

% Replace requested 0° with Amin if true 0 isn't present
A_list_sim = A_list_plot;
A_list_sim(abs(A_list_plot-0)<1e-12) = A0_for_sim;

% Build legend labels (fix zero label to whatever min angle is)
labels = strings(size(A_list_plot));
for k = 1:numel(A_list_plot)
    if abs(A_list_plot(k)-0)<1e-12
        labels(k) = string(label0);
    else
        labels(k) = sprintf('%g°', A_list_plot(k));
    end
end

cols = jet(numel(A_list_sim));
for k = 1:numel(A_list_sim)
    yk = Ysim(A_list_sim(k)*ones(size(Egrid)), Egrid);
    yk = movmean(yk, 7);
    yk(yk<=0) = NaN; % avoid log artifacts
    plot(Egrid, yk, 'LineWidth', 1.6, 'Color', cols(k,:));
end

% EXPT 0° points
scatter(Eexp, Yexp, 90, 'k', 'filled', 'MarkerFaceAlpha', 0.9);

set(gca,'XScale','log'); set(gca,'YScale','log');
xlabel('Ion energy (eV)');
ylabel('Yield (atoms/ion)');
title('He^+ on W: Yield vs Energy (RustBCA + EXPT 0°)');
legend([cellstr(labels) {'EXPT 0°'}], 'Location','northwest');
grid on; box on;

%% ------------------------ Plot 3: Yield vs Angle (RustBCA) ------------------------
figure('Color','w'); hold on;

angle_grid = linspace(Amin, Amax, 241);
energies_to_plot = [200 500 1000 2000 4000 8000 20000];
energies_to_plot = energies_to_plot(energies_to_plot>=Emin & energies_to_plot<=Emax);

cols2 = jet(numel(energies_to_plot));
for k = 1:numel(energies_to_plot)
    Ek = energies_to_plot(k);
    yk = Ysim(angle_grid, Ek*ones(size(angle_grid)));
    yk(yk<=0) = NaN;
    plot(angle_grid, yk, 'LineWidth', 1.6, 'Color', cols2(k,:));
end
set(gca,'YScale','log');
xlabel('Incidence angle (deg)');
ylabel('Yield (atoms/ion)');
title('He^+ on W: Yield vs Angle (RustBCA)');
legend(compose('%g eV', energies_to_plot), 'Location','northwest');
grid on; box on;

%% ------------------------ OPTIONAL: Fit Eckstein/Yamamura angular factor to RustBCA ------------------------
if do_fit_angle_model
    % Angular factor:
    % F(theta) = cos(theta)^(-f) * exp(-g*(sec(theta)-1))
    Fang = @(th_rad, f, g) (cos(th_rad)).^(-f) .* exp(-g.*(1./cos(th_rad) - 1));

    theta_fit = linspace(0, fit_theta_max, 161)';   % deg
    theta_fit_rad = deg2rad(theta_fit);

    fit_f = nan(nEexp,1);
    fit_g = nan(nEexp,1);
    Y_model = nan(nEexp, numel(theta_pts));

    has_fmincon = exist('fmincon','file')==2 && use_fmincon_if_available;

    for i = 1:nEexp
        Ei = Eexp(i);

        % normalize RustBCA ratio at A0_for_sim (0 if present, else Amin)
        Y0_sim = Ysim(A0_for_sim, Ei);
        if ~(isfinite(Y0_sim) && Y0_sim > min_sim_yield)
            warning('Skipping E=%.1f eV due to invalid Ysim(ref)=%g (ref=%s)', Ei, Y0_sim, label0);
            continue;
        end

        Ycurve = Ysim(theta_fit, Ei*ones(size(theta_fit)));
        R = Ycurve ./ Y0_sim;

        mask = isfinite(R) & (R>0);
        if nnz(mask) < 10
            warning('Not enough valid points to fit at E=%.1f eV', Ei);
            continue;
        end

        th = theta_fit_rad(mask);
        Rv = R(mask);

        obj = @(p) sum( (log(max(Rv,1e-30)) - log(max(Fang(th,p(1),p(2)),1e-30))).^2 );

        if has_fmincon
            try
                opts = optimoptions('fmincon','Display','off','Algorithm','sqp');
                psol = fmincon(obj, p0, [], [], [], [], lb, ub, [], opts);
            catch
                has_fmincon = false;
            end
        end

        if ~has_fmincon
            invsig = @(t,lo,hi) lo + (hi-lo) ./ (1 + exp(-t));
            obj_t = @(t) obj([invsig(t(1),lb(1),ub(1)), invsig(t(2),lb(2),ub(2))]);
            t0 = [0 0];
            tsol = fminsearch(obj_t, t0, optimset('Display','off','MaxIter',2000,'MaxFunEvals',2000));
            psol = [invsig(tsol(1),lb(1),ub(1)), invsig(tsol(2),lb(2),ub(2))];
        end

        fit_f(i) = psol(1);
        fit_g(i) = psol(2);

        % modeled yields at theta_pts anchored to EXPT 0° (note: theta_pts includes 0)
        th_pts = deg2rad(theta_pts);
        Fpts = Fang(th_pts, fit_f(i), fit_g(i));
        Y_model(i,:) = Yexp(i) .* Fpts;
    end

    % Write CSV
    angle_col_names = compose('Y_%gdeg', theta_pts);
    T = table(Eexp(:), Yexp(:), fit_f(:), fit_g(:), 'VariableNames', {'E_eV','Y0_expt','f_fit','g_fit'});
    for j = 1:numel(theta_pts)
        T.(angle_col_names{j}) = Y_model(:,j);
    end
    writetable(T, out_csv);
    fprintf('Wrote CSV: %s\n', out_csv);

    % Diagnostic plot at ~1 keV
    [~, idx1keV] = min(abs(Eexp-1000));
    if ~isnan(fit_f(idx1keV))
        figure('Color','w'); hold on;
        thd = linspace(0,fit_theta_max,201)';
        Rsim = Ysim(thd, Eexp(idx1keV)*ones(size(thd))) ./ max(Ysim(A0_for_sim,Eexp(idx1keV)), eps);
        Rfit = Fang(deg2rad(thd), fit_f(idx1keV), fit_g(idx1keV));
        plot(thd, Rsim, 'LineWidth',2);
        plot(thd, Rfit, '--', 'LineWidth',2);
        set(gca,'YScale','log');
        xlabel('Incidence angle (deg)');
        ylabel(sprintf('Y(\\theta)/Y(%s)', label0));
        title(sprintf('Fit check at E=%.0f eV | ref=%s | f=%.3g, g=%.3g', ...
            Eexp(idx1keV), label0, fit_f(idx1keV), fit_g(idx1keV)));
        legend('RustBCA ratio','Eckstein fit','Location','best');
        grid on; box on;
    end
end
