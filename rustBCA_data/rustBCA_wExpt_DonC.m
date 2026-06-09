%% C_sputtering_full_pipeline.m
% Full pipeline for:
%   1) D on C
%   2) C on C
%
% Reads RustBCA NetCDF, overlays IPP 9/82 report data at normal incidence,
% fits Eckstein/Yamamura angular factor to RustBCA angular dependence,
% anchors modeled angular points to report Y(0), writes CSV, and plots.
%
% Plot style:
%   - faint black points: all report experimental rows
%   - bold black points: average report yield at each unique energy
%   - gray diamonds: report calculated values
%
% Report source:
%   IPP 9/82 "Sputtering Data"
%   - D on C experimental normal-incidence data: pages 43-44
%   - D on C calculated normal-incidence data: page 45
%   - C on C report/self-sputtering data: pages 48-49

close all; clear; clc;

%% ------------------------ User inputs ------------------------
case_name = 'C_on_C';   % 'D_on_C' or 'C_on_C'

switch case_name
    case 'D_on_C'
        file    = 'ftridyn_DonC.nc';
        % file    = 'rustBCA_DonC.nc';
        out_csv = 'D_on_C_report_plus_EcksteinAngleFit.csv';
    case 'C_on_C'
        file    = 'ftridynConC.nc';
        out_csv = 'C_on_C_report_plus_EcksteinAngleFit.csv';
    otherwise
        error('Unknown case_name. Use ''D_on_C'' or ''C_on_C''.');
end

% Choose angles (deg) to output in table (must include 0)
theta_pts = [0 30 45 60 75 80];

% Fit controls
fit_theta_max = 80;     % max angle (deg) used in angular fit
min_sim_yield = 1e-14;  % ignore smaller sim yields during fit

% Optimization controls
use_fmincon_if_available = true;
lb = [0.0, 0.0];   % bounds for [f,g]
ub = [6.0, 30.0];
p0 = [1.0, 1.0];

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
assert(isequal(size(Y), [numel(A) numel(E)]), ...
    'Yield array must be [numel(angle) x numel(energy)].');

% sort monotonic
[E, iE] = sort(E, 'ascend');
[A, iA] = sort(A, 'ascend');
Y = Y(iA, iE);

% effective 0 degree from available grid
A0 = min(A);
fprintf('Using A0 = %.3f deg as effective 0 deg for RustBCA normalization.\n', A0);

% positive energies only
maskE = (E > 0) & isfinite(E);
Euse = E(maskE);
Yuse = Y(:, maskE);

Emin = max(min(Euse), 1);
Emax = max(Euse);
fprintf('Using RustBCA energy range: [%.3g, %.3g] eV\n', Emin, Emax);

% interpolation handle: Y(angle,energy)
Ysim = @(ang, en) interpn(A, Euse, Yuse, ang, en, 'linear', 0);

%% ------------------------ Load report data ------------------------
[Eexp_all, Yexp_all, Ecalc_rep, Ycalc_rep] = get_report_data(case_name);

assert(numel(Eexp_all) == numel(Yexp_all), 'Eexp and Yexp must have the same length.');
assert(numel(Ecalc_rep) == numel(Ycalc_rep), 'Ecalc_rep and Ycalc_rep must have the same length.');

% filter to available sim energy range
okE = (Eexp_all >= Emin*0.999) & (Eexp_all <= Emax*1.001) & isfinite(Yexp_all) & (Yexp_all >= 0);
Eexp_all = Eexp_all(okE);
Yexp_all = Yexp_all(okE);

okC = (Ecalc_rep >= Emin*0.999) & (Ecalc_rep <= Emax*1.001) & isfinite(Ycalc_rep) & (Ycalc_rep >= 0);
Ecalc_rep = Ecalc_rep(okC);
Ycalc_rep = Ycalc_rep(okC);

if isempty(Eexp_all)
    error('No report normal-incidence data fall within RustBCA energy range.');
end

% averaged experimental points by unique energy (Option 3)
[Eexp, ~, ic] = unique(Eexp_all);
Yexp = accumarray(ic, Yexp_all, [], @mean);
Yexp_std = accumarray(ic, Yexp_all, [], @std);
Nexp = accumarray(ic, Yexp_all, [], @numel);

% replace NaN std for singleton groups
Yexp_std(~isfinite(Yexp_std)) = 0;

%% ------------------------ Fit Eckstein/Yamamura angular factor ------------------------
Fang = @(theta_rad, f, g) ...
    (cos(theta_rad)).^(-f) .* exp(-g .* (1 ./ cos(theta_rad) - 1));

theta_fit = linspace(0, fit_theta_max, 161)';    % deg
theta_fit_rad = deg2rad(theta_fit);

nEexp = numel(Eexp);
fit_f = nan(nEexp,1);
fit_g = nan(nEexp,1);
Y_model = nan(nEexp, numel(theta_pts));

has_fmincon = exist('fmincon','file') == 2 && use_fmincon_if_available;

for i = 1:nEexp
    Ei = Eexp(i);

    % RustBCA angle curve
    Ysim_curve = Ysim(theta_fit, Ei * ones(size(theta_fit)));

    % normalize to smallest available angle
    Y0_sim = Ysim(A0, Ei);
    if ~(isfinite(Y0_sim) && Y0_sim > min_sim_yield)
        warning('Skipping E = %.1f eV due to invalid Ysim(A0) = %.3g', Ei, Y0_sim);
        continue;
    end

    R = Ysim_curve ./ Y0_sim;
    mask_valid = isfinite(R) & (R > 0);
    if nnz(mask_valid) < 10
        warning('Not enough valid RustBCA points to fit at E = %.1f eV', Ei);
        continue;
    end

    th = theta_fit_rad(mask_valid);
    Rv = R(mask_valid);

    obj = @(p) sum((log(max(Rv,1e-30)) - log(max(Fang(th,p(1),p(2)),1e-30))).^2);

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
        t0 = [ ...
            log((p0(1)-lb(1)+1e-6)/(ub(1)-p0(1)+1e-6)), ...
            log((p0(2)-lb(2)+1e-6)/(ub(2)-p0(2)+1e-6)) ...
        ];
        tsol = fminsearch(obj_t, t0, ...
            optimset('Display','off','MaxIter',2000,'MaxFunEvals',2000));
        psol = [invsig(tsol(1),lb(1),ub(1)), invsig(tsol(2),lb(2),ub(2))];
    end

    fit_f(i) = psol(1);
    fit_g(i) = psol(2);

    % modeled yields at requested angles, anchored to averaged report Y(0)
    th_pts_rad = deg2rad(theta_pts);
    Fpts = Fang(th_pts_rad, fit_f(i), fit_g(i));
    Y_model(i,:) = Yexp(i) .* Fpts;
end

%% ------------------------ Build output table ------------------------
angle_col_names = compose('Y_%gdeg', theta_pts);

T = table(Eexp(:), Yexp(:), Yexp_std(:), Nexp(:), fit_f(:), fit_g(:), ...
    'VariableNames', {'E_eV','Y0_report_avg','Y0_report_std','N_report','f_fit','g_fit'});

for j = 1:numel(theta_pts)
    T.(angle_col_names{j}) = Y_model(:,j);
end

writetable(T, out_csv);
fprintf('Wrote output table: %s (%d rows)\n', out_csv, height(T));

%% ------------------------ Plots ------------------------
Egrid = logspace(log10(max(Emin,1)), log10(Emax), 2000);
Agrid = linspace(min(A), max(A), 181);

% 1) Yield map
figure('Color','w');
h = pcolor(Euse, A, Yuse);
set(h,'EdgeColor','none');
set(gca,'XScale','log');
xlabel('Ion energy (eV)');
ylabel('Incidence angle (deg)');
title(sprintf('%s: Sputtering yield map (RustBCA)', strrep(case_name,'_',' ')));
colorbar; grid on; box on;

% 2) Angle dependence at selected energies
figure('Color','w'); hold on;
Elist = unique(max(min([100 200 500 1000 5000 10000 50000], Emax), Emin));
cols = jet(numel(Elist));
for k = 1:numel(Elist)
    yk = Ysim(Agrid, Elist(k) * ones(size(Agrid)));
    plot(Agrid, yk, 'LineWidth', 1.6, 'Color', cols(k,:));
end
set(gca,'YScale','log');
xlabel('Incidence angle (deg)');
ylabel('Yield (atoms/ion)');
title(sprintf('%s: Angular dependence (RustBCA)', strrep(case_name,'_',' ')));
legend(compose('%g eV', Elist), 'Location', 'northwest');
grid on; box on;

% 3) Energy dependence: RustBCA curves + report data
figure('Color','w'); hold on;
A_list_plot = [0 10 30 50 60 70 80 85];
A_list_plot_sim = A_list_plot;
A_list_plot_sim(A_list_plot == 0) = A0;

maskA = (A_list_plot_sim >= min(A) & A_list_plot_sim <= max(A));
A_list_plot = A_list_plot(maskA);
A_list_plot_sim = A_list_plot_sim(maskA);

cols2 = jet(numel(A_list_plot_sim));
for k = 1:numel(A_list_plot_sim)
    yk = Ysim(A_list_plot_sim(k) * ones(size(Egrid)), Egrid);
    yk = movmean(yk, 9);
    plot(Egrid, yk, 'LineWidth', 1.6, 'Color', cols2(k,:));
end

% All experimental rows faintly
scatter(Eexp_all, Yexp_all, 20, 'k', 'filled', ...
    'MarkerFaceAlpha', 0.22, ...
    'MarkerEdgeAlpha', 0.22);

% Averaged experimental points boldly
scatter(Eexp, Yexp, 70, 'k', 'filled', ...
    'MarkerFaceAlpha', 0.95);

% Calculated report points as gray diamonds
if ~isempty(Ecalc_rep)
    scatter(Ecalc_rep, Ycalc_rep, 70, 'd', ...
        'MarkerEdgeColor', [0.35 0.35 0.35], ...
        'MarkerFaceColor', [0.85 0.85 0.85], ...
        'LineWidth', 1.1, ...
        'MarkerFaceAlpha', 0.95);
end

set(gca,'XScale','log');
set(gca,'YScale','log');
xlabel('Ion energy (eV)');
ylabel('Yield (atoms/ion)');
title(sprintf('%s: Energy dependence', strrep(case_name,'_',' ')));

lgd = [compose('%g° RustBCA', A_list_plot), ...
       "0° expt. (IPP 9/82, all data)", ...
       "0° expt. (IPP 9/82, energy-averaged)"];
if ~isempty(Ecalc_rep)
    lgd = [lgd, "0° Bohdansky fit (IPP 9/82)"];
end
legend(lgd, 'Location', 'northwest');

grid on; box on;

% 4) Dedicated 0 degree comparison
figure('Color','w'); hold on;

plot(Egrid, Ysim(A0 * ones(size(Egrid)), Egrid), 'LineWidth', 2.5);

% All experimental rows faintly
scatter(Eexp_all, Yexp_all, 20, 'k', 'filled', ...
    'MarkerFaceAlpha', 0.22, ...
    'MarkerEdgeAlpha', 0.22);

% Averaged experimental points boldly
scatter(Eexp, Yexp, 85, 'k', 'filled', ...
    'MarkerFaceAlpha', 0.95);

% Calculated report points
if ~isempty(Ecalc_rep)
    scatter(Ecalc_rep, Ycalc_rep, 70, 'd', ...
        'MarkerEdgeColor', [0.35 0.35 0.35], ...
        'MarkerFaceColor', [0.85 0.85 0.85], ...
        'LineWidth', 1.1, ...
        'MarkerFaceAlpha', 0.95);
end

% Modeled angle points anchored to averaged Y(0)
for i = 1:nEexp
    for j = 1:numel(theta_pts)
        plot(Eexp(i), Y_model(i,j), 'o', ...
            'MarkerEdgeColor', 'none', ...
            'MarkerFaceColor', [0.85 0.33 0.10], ...
            'MarkerSize', 6);
    end
end

set(gca,'XScale','log');
set(gca,'YScale','log');
xlabel('Ion energy (eV)');
ylabel('Yield (atoms/ion)');
title(sprintf('%s: RustBCA(A0=%.2f°) vs report 0° + modeled angular points', ...
    strrep(case_name,'_',' '), A0));

legend_entries = { ...
    sprintf('0° RustBCA (effective %.2f° grid)', A0), ...
    '0° expt. (IPP 9/82, all data)', ...
    '0° expt. (IPP 9/82, energy-averaged)'};

if ~isempty(Ecalc_rep)
    legend_entries{end+1} = '0° Bohdansky fit (IPP 9/82)';
end

legend_entries{end+1} = 'Eckstein/Yamamura angular scaling (anchored to expt. avg)';

legend(legend_entries, 'Location', 'northwest');
grid on; box on;

% 5) Angular fit diagnostic near 1 keV
[~, idx1keV] = min(abs(Eexp - 1000));
if ~isnan(fit_f(idx1keV)) && ~isnan(fit_g(idx1keV))
    figure('Color','w'); hold on;
    thd = linspace(0, fit_theta_max, 201)';
    Rsim = Ysim(thd, Eexp(idx1keV) * ones(size(thd))) ./ max(Ysim(A0, Eexp(idx1keV)), eps);
    Rfit = Fang(deg2rad(thd), fit_f(idx1keV), fit_g(idx1keV));
    plot(thd, Rsim, 'LineWidth', 2);
    plot(thd, Rfit, '--', 'LineWidth', 2);
    set(gca,'YScale','log');
    xlabel('Incidence angle (deg)');
    ylabel('Y(\theta)/Y(0)');
    title(sprintf('Angular fit check at E = %.0f eV (normalized at %.2f°)', ...
        Eexp(idx1keV), A0));
    legend('RustBCA ratio', 'Eckstein/Yamamura fit', 'Location', 'best');
    grid on; box on;
end

fprintf('Done. Output CSV: %s\n', out_csv);

%% ========================================================================
function [Eexp, Yexp, Ecalc, Ycalc] = get_report_data(case_name)
% Returns report data in eV and atoms/ion

switch case_name

    case 'D_on_C'
        % Experimental D -> C from IPP 9/82 pages 43-44
        D_on_C_exp = [ ...
            0.020  0.03560
            0.020  0.03880
            0.030  0.03620
            0.050  0.03770
            0.050  0.03570
            0.050  0.04810
            0.050  0.03820
            0.050  0.03880
            0.050  0.04150
            0.050  0.03830
            0.050  0.04290
            0.060  0.03370
            0.060  0.01760
            0.070  0.01380
            0.080  0.03170
            0.080  0.03800
            0.080  0.02070
            0.100  0.02390
            0.100  0.02320
            0.100  0.03170
            0.100  0.02720
            0.100  0.02960
            0.100  0.02720
            0.100  0.02950
            0.100  0.01760
            0.120  0.02980
            0.120  0.03270
            0.120  0.00596
            0.150  0.02670
            0.250  0.02390
            0.300  0.02280
            0.300  0.03320
            0.350  0.02120
            0.500  0.02220
            0.500  0.02860
            0.500  0.03630
            0.500  0.02630
            0.500  0.04980
            0.500  0.05900
            1.000  0.01800
            1.000  0.01670
            1.000  0.02030
            1.000  0.02050
            1.000  0.02310
            1.000  0.02070
            1.000  0.01710
            1.000  0.01780
            1.000  0.03250
            1.000  0.03210
            1.000  0.02290
            1.000  0.02050
            1.000  0.01550
            1.000  0.01830
            1.000  0.03620
            1.000  0.03140
            1.000  0.02960
            1.000  0.03210
            1.000  0.02780
            1.000  0.04500
            2.000  0.03360
            2.000  0.01770
            2.000  0.01140
            2.000  0.01320
            2.000  0.01150
            2.000  0.01170
            2.000  0.02880
            2.000  0.02300
            2.000  0.01860
            2.000  0.01740
            2.000  0.01600
            2.000  0.02360
            4.000  0.02290
            4.000  0.00770
            8.000  0.00742
            8.000  0.00495
        ];

        Eexp = D_on_C_exp(:,1) * 1e3;
        Yexp = D_on_C_exp(:,2);

        % Calculated D -> C from IPP 9/82 page 45
        D_on_C_calc = [ ...
            0.050  0.001910
            0.050  0.001860
            0.070  0.005780
            0.100  0.009570
            0.100  0.009800
            0.140  0.012200
            0.200  0.014400
            0.300  0.016100
            0.400  0.015300
            0.500  0.012900
            0.700  0.014300
            1.000  0.012300
            2.000  0.009650
            5.000  0.007400
        ];

        Ecalc = D_on_C_calc(:,1) * 1e3;
        Ycalc = D_on_C_calc(:,2);

    case 'C_on_C'
        % Report/self-sputtering C -> C data from IPP 9/82 pages 48-49
        C_on_C_exp = [ ...
            0.075  0.01800
            0.100  0.12800
            0.150  0.12500
            0.300  0.23300
            0.600  0.34100
            1.000  0.37000
            3.000  0.44300
            3.000  0.32900
            3.000  0.47400
            10.000 0.50700
        ];

        Eexp = C_on_C_exp(:,1) * 1e3;
        Yexp = C_on_C_exp(:,2);

        % Calculated C -> C from IPP 9/82 page 49
        C_on_C_calc = [ ...
            0.050  0.000740
            0.070  0.003700
            0.100  0.012000
            0.100  0.012000
            0.150  0.031000
            0.200  0.056000
            0.300  0.089000
            0.500  0.140000
            0.700  0.180000
            1.000  0.195000
            1.000  0.195000
            1.500  0.215000
            2.000  0.225000
            3.000  0.240000
            5.000  0.230000
            10.000 0.210000
        ];

        Ecalc = C_on_C_calc(:,1) * 1e3;
        Ycalc = C_on_C_calc(:,2);

    otherwise
        error('Unknown case_name.');
end

Eexp  = Eexp(:);
Yexp  = Yexp(:);
Ecalc = Ecalc(:);
Ycalc = Ycalc(:);

end