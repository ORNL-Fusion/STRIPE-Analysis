%% D_on_W_full_pipeline.m
% Full pipeline: RustBCA read, EXPT 0° overlay, Eckstein/Yamamura angle-fit to RustBCA,
% produce CSV with modeled angle points and plots.

close all; clear; clc;

%% ------------------------ User inputs ------------------------
file = 'ftridyn_DonW.nc';                    % RustBCA NetCDF file for D on W
out_csv = 'D_on_W_EXPT0_plus_EcksteinAngleFit.csv';

% Choose angles (deg) to output in table (must include 0)
theta_pts = [0 30 45 60 75 80];

% Fit controls
fit_theta_max = 80;           % max angle (deg) used for fitting to avoid grazing noise
min_sim_yield  = 1e-14;       % ignore sim yields smaller than this when fitting

% Optimization controls
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

% ------------------------ KEY FIX: effective 0° from RustBCA grid ------------------------
A0 = min(A);   % RustBCA may not include 0 exactly; use smallest available angle
fprintf('Using A0 = %.3f° as effective 0° for RustBCA normalization/plots.\n', A0);

% positive energies only
maskE = (E > 0) & isfinite(E);
Euse = E(maskE);
Yuse = Y(:, maskE);

Emin = max(min(Euse), 1);
Emax = max(Euse);
fprintf('Using RustBCA energy range: [%.3g, %.3g] eV\n', Emin, Emax);

% create interpolation handle: Y(angle,energy)
Ysim = @(ang, en) interpn(A, Euse, Yuse, ang, en, 'linear', 0);

%% ------------------------ Experimental D+ -> W (0°) - IPP 9/82 ------------------------
% ---------------- Experimental D+ → W data (θ = 0° only) ----------------
% IPP 9/82, pages 116–117
Eexp = [ ...
   0.25 0.35 0.50 1.0 1.5 2.0 2.0 2.0 ...
   4.0 8.0 15.0 30.0 50.0 50.0 100.0 ] * 1e3;   % eV

Yexp = [ ...
   1.5e-4 3.1e-4 7.0e-4 2.11e-3 4.7e-3 ...
   4.56e-3 4.83e-3 4.52e-3 7.0e-3 ...
   6.38e-3 4.7e-3 4.15e-3 2.4e-3 4.34e-3 2.7e-3 ];

Eexp = Eexp(:); Yexp = Yexp(:);

% filter to available Euse range
okE = (Eexp >= Emin*0.999) & (Eexp <= Emax*1.001) & isfinite(Yexp) & (Yexp>=0);
Eexp = Eexp(okE); Yexp = Yexp(okE);

if isempty(Eexp)
    error('No experimental energies fall within RustBCA energy range.');
end

%% ------------------------ Fitting Eckstein/Yamamura angular factor per energy ------------------------
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

    % get rustBCA angle curve (dense)
    Ysim_curve = Ysim(theta_fit, Ei*ones(size(theta_fit)));

    % ---- KEY FIX: normalize to A0 instead of 0 ----
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

% 1) Yield map (RustBCA)
figure('Color','w');
h = pcolor(Euse, A, Yuse);
set(h,'EdgeColor','none'); set(gca,'XScale','log');
xlabel('Ion energy (eV)'); ylabel('Incidence angle (deg)');
title('D^+ on W: Sputtering yield map (RustBCA)');
colorbar; grid on; box on;

% 2) Angle dependence (sim) at few energies (angle axis linear)
figure('Color','w'); hold on;
Elist = unique( max(min([100 200 500 1000 5000 10000 50000], Emax), Emin) );
cols = jet(numel(Elist));
for k = 1:numel(Elist)
    yk = Ysim(Agrid, Elist(k)*ones(size(Agrid)));
    plot(Agrid, yk, 'LineWidth',1.6,'Color',cols(k,:));
end
set(gca,'YScale','log');
xlabel('Incidence angle (deg)'); ylabel('Yield (atoms/ion)');
title('D on W: Yield (RustBCA)');
legend(compose('%g eV', Elist),'Location','northwest'); grid on; box on;
ylim([1e-6 1e0])

% 3) Energy dependence (sim curves) + EXPT 0°
figure('Color','w'); hold on;
A_list_plot = [0 10 30 50 60 70 80 85];

A_list_plot_sim = A_list_plot;
A_list_plot_sim(A_list_plot==0) = A0; % effective 0° for sim

maskA = (A_list_plot_sim>=min(A) & A_list_plot_sim<=max(A));
A_list_plot = A_list_plot(maskA);
A_list_plot_sim = A_list_plot_sim(maskA);

cols2 = jet(numel(A_list_plot_sim));
for k = 1:numel(A_list_plot_sim)
    yk = Ysim(A_list_plot_sim(k)*ones(size(Egrid)), Egrid);
    yk = movmean(yk,9);
    plot(Egrid, yk, 'LineWidth',1.6,'Color',cols2(k,:));
end

scatter(Eexp, Yexp, 80, 'k', 'filled', 'MarkerFaceAlpha', 0.8);

set(gca,'XScale','log'); set(gca,'YScale','log');
xlabel('Ion energy (eV)'); ylabel('Yield (atoms/ion)');
title('D on W: Yield (RustBCA)'); grid on; box on;
legend([compose('%g°', A_list_plot) "EXPT 0°"], 'Location','northwest');ylim([1e-6 1e0])

% 4) Dedicated 0° comparison (RustBCA uses A0)
figure('Color','w'); hold on;
plot(Egrid, Ysim(A0*ones(size(Egrid)), Egrid), 'LineWidth', 2.5);
scatter(Eexp, Yexp, 100, 'k', 'filled');

for i = 1:nEexp
    for j = 1:numel(theta_pts)
        plot(Eexp(i), Y_model(i,j), 'o', 'MarkerEdgeColor','none', ...
            'MarkerFaceColor', [0.85 0.33 0.10], 'MarkerSize',6);
    end
end

set(gca,'XScale','log'); set(gca,'YScale','log');
xlabel('Ion energy (eV)'); ylabel('Yield (atoms/ion)');
title(sprintf('D^+ on W: RustBCA(A0=%.2f°) vs EXPT 0° (and modeled angular points)', A0));
legend(sprintf('RustBCA %.2f° (effective 0°)',A0),'EXPT 0°','Modeled angle points','Location','northwest');
grid on; box on;

% 5) Quick diagnostic: fit vs RustBCA ratio at representative energy (~1 keV)
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