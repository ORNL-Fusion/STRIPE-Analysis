%% W_on_W_full_pipeline_with_HECHTL90_EXPT.m
% Full pipeline for W on W:
% - Reads RustBCA_WonW.nc (auto-detects schema)
% - Uses experimental 0° table from HECHTL90 (screenshot you provided)
% - Fits Eckstein/Yamamura angular factor to RustBCA angular ratios (normalized at A0)
% - Anchors fitted angular model to the EXPT 0° yields
% - Produces comparison plots and writes CSV with modeled angle points
%
% Drop this script in the folder with RustBCA_WonW.nc and run.

close all; clear; clc;

%% ------------------------ User inputs ------------------------
% file = 'rustbca_WonW.nc';     % RustBCA NetCDF file for W on W
file = 'ftridynSelf_new_higResol.nc'; 
out_csv = 'W_on_W_HECHTL90_EXPT_plus_AngleModel.csv';

% choose angles (deg) to output in table (must include 0)
theta_pts = [0 15 30 45 60 75 80];

% fitting controls
fit_theta_max = 80;           % max angle (deg) used for fitting to avoid grazing noise
min_sim_yield  = 1e-14;       % ignore sim yields smaller than this when fitting

% optimization controls
use_fmincon_if_available = true;   % prefer fmincon if Optimization Toolbox present
lb = [0.0, 0.0];   % bounds for [f,g]
ub = [6.0, 30.0];  % generous upper bounds
p0 = [1.0, 1.0];   % initial guess for [f,g]

%% ------------------------ Read RustBCA (auto-detect schema) ------------------------
info = ncinfo(file);
varlist = string({info.Variables.Name});

hasOld = all(ismember(["E","A","spyld"], varlist));
hasNew = all(ismember(["energy_eV","angle_deg","Yields"], varlist));

if hasOld
    fprintf('Detected OLD schema (E,A,spyld)\n');
    E = double(ncread(file,'E'));   E = E(:);
    A = double(ncread(file,'A'));   A = A(:);
    Y = double(ncread(file,'spyld'));    % expected [angle x energy]
elseif hasNew
    fprintf('Detected NEW schema (energy_eV,angle_deg,Yields)\n');
    E = double(ncread(file,'energy_eV')); E = E(:);
    A = double(ncread(file,'angle_deg')); A = A(:);
    Y = double(ncread(file,'Yields'));
    if isequal(size(Y), [numel(E) numel(A)])
        Y = Y.'; % transpose if raw is [energy x angle]
    end
else
    fprintf('Variables in file:\n'); disp(varlist.');
    error('Unknown netCDF schema. Expected either (E,A,spyld) or (energy_eV,angle_deg,Yields).');
end

% ensure orientation [angle x energy]
if isequal(size(Y), [numel(E) numel(A)])
    Y = Y.';
end
assert(isequal(size(Y), [numel(A) numel(E)]), 'Yield array must be [numel(angle) x numel(energy)].');

% sort monotonic
[E, iE] = sort(E, 'ascend');
[A, iA] = sort(A, 'ascend');
Y = Y(iA, iE);

% positive energies only
maskE = (E > 0) & isfinite(E);
Euse = E(maskE);
Yuse = Y(:, maskE);

Emin = max(min(Euse), 1);
Emax = max(Euse);
fprintf('Using RustBCA energy range: [%.3g, %.3g] eV\n', Emin, Emax);

% effective 0° from RustBCA grid (some sims may not include exact 0)
A0 = min(A);
fprintf('Using A0 = %.3f° as effective 0° for RustBCA normalization/plots.\n', A0);

% create interpolation handle: Y(angle,energy)
Ysim = @(ang, en) interpn(A, Euse, Yuse, ang, en, 'linear', 0);

%% ------------------------ Experimental W -> W (0°) - HECHTL90 screenshot ------------------------
% Energies (keV) and yields at theta=0 from your screenshot (HECHTL90)
Eexp_keV = [0.100 0.150 0.300 1.000 3.000 10.000];  % keV
Eexp = 1e3 * Eexp_keV;                              % eV

Yexp = [0.03300 0.05900 0.22700 1.03000 2.33000 3.82000];  % yields (atoms/ion)

Eexp = Eexp(:);
Yexp = Yexp(:);

% filter to available Euse range
okE = (Eexp >= Emin*0.999) & (Eexp <= Emax*1.001) & isfinite(Yexp) & (Yexp>=0);
Eexp = Eexp(okE); Yexp = Yexp(okE);
if isempty(Eexp)
    error('No experimental energies fall within RustBCA energy range.');
end

%% ------------------------ Fit Eckstein/Yamamura angular factor per energy ------------------------
% Angular factor:
% F(theta) = cos(theta)^(-f) * exp(-g*(sec(theta)-1))
Fang = @(theta_rad, f, g) (cos(theta_rad)).^(-f) .* exp(-g.*(1./cos(theta_rad) - 1));

theta_fit = linspace(0, fit_theta_max, 161)';   % deg
theta_fit_rad = deg2rad(theta_fit);

nEexp = numel(Eexp);
fit_f = nan(nEexp,1);
fit_g = nan(nEexp,1);
Y_model = nan(nEexp, numel(theta_pts));

has_fmincon = exist('fmincon','file')==2 && use_fmincon_if_available;

for i = 1:nEexp
    Ei = Eexp(i);

    % RustBCA angle curve (dense)
    Ysim_curve = Ysim(theta_fit, Ei*ones(size(theta_fit)));

    % normalize to A0 instead of 0 (robust if sim has no exact 0°)
    Y0_sim = Ysim(A0, Ei);
    if ~(isfinite(Y0_sim) && Y0_sim > min_sim_yield)
        warning('Skipping E=%.1f eV due to invalid Ysim(A0)=%.3g (A0=%.3f°)', Ei, Y0_sim, A0);
        continue;
    end

    % relative RustBCA angular factor (normalized to A0)
    R = Ysim_curve ./ Y0_sim;
    mask_valid = isfinite(R) & (R>0);
    if nnz(mask_valid) < 10
        warning('Not enough valid RustBCA points to fit at E=%.1f eV', Ei);
        continue;
    end

    th = theta_fit_rad(mask_valid);   % radians
    Rv = R(mask_valid);

    % objective: least-squares in log-space between Rv and Fang(th;f,g)
    obj = @(p) sum( (log(max(Rv,1e-30)) - log(max(Fang(th,p(1),p(2)),1e-30))).^2 );

    % try fmincon with bounds if available
    if has_fmincon
        try
            opts = optimoptions('fmincon','Display','off','Algorithm','sqp');
            psol = fmincon(obj, p0, [], [], [], [], lb, ub, [], opts);
        catch
            has_fmincon = false;
        end
    end

    if ~has_fmincon
        % bounded search via transform + fminsearch
        invsig = @(t,lo,hi) lo + (hi-lo) ./ (1 + exp(-t));
        obj_t = @(t) obj([invsig(t(1),lb(1),ub(1)), invsig(t(2),lb(2),ub(2))]);
        t0 = [log((p0(1)-lb(1)+1e-6)/(ub(1)-p0(1)+1e-6)), log((p0(2)-lb(2)+1e-6)/(ub(2)-p0(2)+1e-6))];
        tsol = fminsearch(obj_t, t0, optimset('Display','off','MaxIter',2000,'MaxFunEvals',2000));
        psol = [invsig(tsol(1),lb(1),ub(1)), invsig(tsol(2),lb(2),ub(2))];
    end

    fit_f(i) = psol(1);
    fit_g(i) = psol(2);

    % produce modeled yields at theta_pts anchored to experimental Y0
    th_pts_rad = deg2rad(theta_pts);
    Fpts = Fang(th_pts_rad, fit_f(i), fit_g(i));
    Y_model(i,:) = Yexp(i) .* Fpts;
end

%% ------------------------ Build output table and write CSV ------------------------
angle_col_names = compose('Y_%gdeg', theta_pts);
T = table(Eexp(:), Yexp(:), fit_f(:), fit_g(:), 'VariableNames', {'E_eV','Y0_expt','f_fit','g_fit'});
for j = 1:numel(theta_pts)
    T.(angle_col_names{j}) = Y_model(:,j);
end
writetable(T, out_csv);
fprintf('Wrote output table: %s (%d rows)\n', out_csv, height(T));

%% ------------------------ Plots ------------------------
Egrid = logspace(log10(max(Emin,1)), log10(Emax), 2000);
Agrid = linspace(min(A), max(A), 181);

% 1) Yield map (RustBCA) + EXPT 0° overlay (plotted at A0 so it appears)
figure('Color','w');
h = pcolor(Euse, A, Yuse);
set(h,'EdgeColor','none'); set(gca,'XScale','log');
xlabel('Ion energy (eV)'); ylabel('Incidence angle (deg)');
title('W on W: Sputtering yield map (RustBCA) + HECHTL90 EXPT @0°');
colorbar; grid on; box on;

hold on;
scatter(Eexp, A0*ones(size(Eexp)), 96, 'k', 'filled', 'MarkerFaceAlpha', 0.9);
% annotate each point with its yield (optional)
for ii=1:numel(Eexp)
    text(Eexp(ii)*1.05, A0, sprintf(' %.3g', Yexp(ii)), 'VerticalAlignment','bottom', 'FontSize',8);
end
hold off;

% 2) Angle dependence (sim) at representative energies (angle axis linear)
figure('Color','w'); hold on;
% choose energies for plotting (in eV) within RustBCA range
Elist = unique( max(min([100 300 800 1000 2000 2500 5000 10000], Emax), Emin) );
cols = jet(numel(Elist));
for k = 1:numel(Elist)
    yk = Ysim(Agrid, Elist(k)*ones(size(Agrid)));
    plot(Agrid, yk, 'LineWidth',1.6,'Color',cols(k,:));
end
set(gca,'YScale','log');
xlabel('Incidence angle (deg)'); ylabel('Yield (atoms/ion)');
title('W on W: Yield vs Angle (RustBCA)');
legend(compose('%g eV', Elist),'Location','northwest'); grid on; box on;

% 3) Energy dependence (sim curves) + EXPT 0°
figure('Color','w'); hold on;
A_list_plot = [0 10 30 50 60 70 80 85];
A_list_plot_sim = A_list_plot;
A_list_plot_sim(A_list_plot==0) = A0; % effective 0° for RustBCA

maskA = (A_list_plot_sim>=min(A) & A_list_plot_sim<=max(A));
A_list_plot = A_list_plot(maskA);
A_list_plot_sim = A_list_plot_sim(maskA);

cols2 = jet(numel(A_list_plot_sim));
for k = 1:numel(A_list_plot_sim)
    yk = Ysim(A_list_plot_sim(k)*ones(size(Egrid)), Egrid);
    yk = movmean(yk,9);
    plot(Egrid, yk, 'LineWidth',1.6,'Color',cols2(k,:));
end

% overlay EXPT 0° points
scatter(Eexp, Yexp, 100, 'k', 'filled', 'MarkerFaceAlpha', 0.9);
legend_entries = [compose('%g°', A_list_plot), "HECHTL90 EXPT 0°"];
legend(legend_entries, 'Location','northwest');
set(gca,'XScale','log'); set(gca,'YScale','log');
xlabel('Ion energy (eV)'); ylabel('Yield (atoms/ion)');
title('W on W: Yield vs Energy (RustBCA + HECHTL90 EXPT 0°)'); grid on; box on;

% 4) Dedicated 0° comparison: RustBCA(A0) vs EXPT, plus modeled angular points
figure('Color','w'); hold on;
plot(Egrid, Ysim(A0*ones(size(Egrid)), Egrid), 'LineWidth', 2.5, 'Color',[0 0.4470 0.7410]);
scatter(Eexp, Yexp, 120, 'k', 'filled');

% also plot modeled Y at selected theta_pts (for each Eexp we have modeled yields)
for i = 1:nEexp
    for j = 1:numel(theta_pts)
        plot(Eexp(i), Y_model(i,j), 'o', 'MarkerEdgeColor','none', 'MarkerFaceColor', [0.85 0.33 0.10], 'MarkerSize',7);
    end
end

set(gca,'XScale','log'); set(gca,'YScale','log');
xlabel('Ion energy (eV)'); ylabel('Yield (atoms/ion)');
title(sprintf('W on W: RustBCA(A0=%.2f°) vs HECHTL90 EXPT 0° (and modeled angular points)', A0));
legend('RustBCA (A0)','HECHTL90 EXPT 0°','Modeled angle points','Location','northwest');
grid on; box on;

% 5) Quick diagnostic: show one fit vs RustBCA ratio for representative energy (closest to 1 keV)
[~, idx1keV] = min(abs(Eexp - 1000));
if ~isnan(fit_f(idx1keV)) && ~isnan(fit_g(idx1keV))
    figure('Color','w'); hold on;
    thd = linspace(0,fit_theta_max,201)';
    Rsim = Ysim(thd, Eexp(idx1keV)*ones(size(thd))) ./ max(Ysim(A0, Eexp(idx1keV)), eps);
    Rfit = Fang(deg2rad(thd), fit_f(idx1keV), fit_g(idx1keV));
    plot(thd, Rsim, 'LineWidth',2);
    plot(thd, Rfit, '--', 'LineWidth',2);
    set(gca,'YScale','log');
    xlabel('Incidence angle (deg)'); ylabel('Y(\theta)/Y(0)');
    title(sprintf('Angular fit check at E = %.0f eV (normalized at %.2f°)', Eexp(idx1keV), A0));
    legend('RustBCA ratio','Eckstein fit','Location','best'); grid on; box on;
end

fprintf('Done. Outputs:\n - CSV: %s\n', out_csv);