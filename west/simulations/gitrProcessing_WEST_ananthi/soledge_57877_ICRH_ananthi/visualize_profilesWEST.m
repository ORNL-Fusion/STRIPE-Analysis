%% Visualize GITR/SOLEDGE profiles from profilesWEST.nc
% Standalone visualization script.
%
% Reads from:
%   profilesWEST.nc
%
% Plots:
%   - B-field: br, bt, bz
%   - electrons: ne, te
%   - D+: ni, ti, vr, vt, vz, vp
%   - oxygen O1+--O8+: no1--no8, to1--to8, velocities if present
%   - geometry stored in .nc:
%       wall_r, wall_z
%       sep1_r, sep1_z
%       sep2_r, sep2_z
%       antenna_centroid_r, antenna_centroid_z
%   - mid-Z and mid-R lineouts for ne and Te
%   - optional OMP experimental reflectometry comparison
%
% This script does NOT write or modify profilesWEST.nc.

clc;
close all;

%% ---------------------------
%  User settings
%  ---------------------------
ncfile = 'profilesWEST.nc';

save_plots = false;
plot_dir = 'profilesWEST_visualization';

% Plot switches
plot_B_fields       = true;
plot_electrons      = true;
plot_Dplus          = true;
plot_oxygen_density = true;
plot_oxygen_temp    = true;
plot_oxygen_velocity = false;   % many figures; set true if needed
plot_lineouts       = true;

% Oxygen charge states
nO = 8;

% Experimental OMP ne comparison
compare_exp_omp_ne = true;
exp_ne_file = 'ICRH_57877_8s.mat';
exp_t_min = 5;
exp_t_max = 10;
omp_z_target = 0.0;

% Plot style
fontSize = 18;
use_log_density = true;

if save_plots && ~exist(plot_dir, 'dir')
    mkdir(plot_dir);
end

assert(isfile(ncfile), 'NetCDF file not found: %s', ncfile);

%% ---------------------------
%  Read coordinates and geometry
%  ---------------------------
R = double(ncread(ncfile, 'x'));
Z = double(ncread(ncfile, 'z'));

geom = readGeometryFromNc(ncfile);

fprintf('Read grid from %s: nR=%d, nZ=%d\n', ncfile, numel(R), numel(Z));

%% ---------------------------
%  Read main profiles
%  ---------------------------
profiles = struct();

profiles.br = readNcIfExists(ncfile, 'br');
profiles.bt = readNcIfExists(ncfile, 'bt');
profiles.bz = readNcIfExists(ncfile, 'bz');

profiles.ne = readNcIfExists(ncfile, 'ne');
profiles.te = readNcIfExists(ncfile, 'te');

profiles.ni = readNcIfExists(ncfile, 'ni');
profiles.ti = readNcIfExists(ncfile, 'ti');

profiles.vr = readNcIfExists(ncfile, 'vr');
profiles.vt = readNcIfExists(ncfile, 'vt');
profiles.vz = readNcIfExists(ncfile, 'vz');
profiles.vp = readNcIfExists(ncfile, 'vp');

%% ---------------------------
%  Plot B fields
%  ---------------------------
if plot_B_fields
    plotIfAvailable(R, Z, profiles.br, 'B_R', '$B_R$ [T]', false, geom, save_plots, plot_dir, 'Br', fontSize);
    plotIfAvailable(R, Z, profiles.bt, 'B_t', '$B_t$ [T]', false, geom, save_plots, plot_dir, 'Bt', fontSize);
    plotIfAvailable(R, Z, profiles.bz, 'B_Z', '$B_Z$ [T]', false, geom, save_plots, plot_dir, 'Bz', fontSize);
end

%% ---------------------------
%  Plot electrons
%  ---------------------------
if plot_electrons
    plotIfAvailable(R, Z, profiles.ne, 'Electron density $n_e$', '$n_e$ [m$^{-3}$]', use_log_density, geom, save_plots, plot_dir, 'ne', fontSize);
    plotIfAvailable(R, Z, profiles.te, 'Electron temperature $T_e$', '$T_e$ [eV]', false, geom, save_plots, plot_dir, 'Te', fontSize);
end

%% ---------------------------
%  Plot D+
%  ---------------------------
if plot_Dplus
    plotIfAvailable(R, Z, profiles.ni, 'D$^+$ density $n_i$', '$n_i$ [m$^{-3}$]', use_log_density, geom, save_plots, plot_dir, 'ni_Dplus', fontSize);
    plotIfAvailable(R, Z, profiles.ti, 'D$^+$ temperature $T_i$', '$T_i$ [eV]', false, geom, save_plots, plot_dir, 'Ti_Dplus', fontSize);

    plotIfAvailable(R, Z, profiles.vr, 'D$^+$ $V_R$', '$V_R$ [m/s]', false, geom, save_plots, plot_dir, 'Vr_Dplus', fontSize);
    plotIfAvailable(R, Z, profiles.vt, 'D$^+$ $V_t$', '$V_t$ [m/s]', false, geom, save_plots, plot_dir, 'Vt_Dplus', fontSize);
    plotIfAvailable(R, Z, profiles.vz, 'D$^+$ $V_Z$', '$V_Z$ [m/s]', false, geom, save_plots, plot_dir, 'Vz_Dplus', fontSize);
    plotIfAvailable(R, Z, profiles.vp, 'D$^+$ $|V|$', '$|V|$ [m/s]', false, geom, save_plots, plot_dir, 'Vp_Dplus', fontSize);
end

%% ---------------------------
%  Plot oxygen O1+ to O8+
%  ---------------------------
for q = 1:nO
    no_name = sprintf('no%d', q);
    to_name = sprintf('to%d', q);

    no_q = readNcIfExists(ncfile, no_name);
    to_q = readNcIfExists(ncfile, to_name);

    if plot_oxygen_density
        plotIfAvailable(R, Z, no_q, sprintf('O$^{%d+}$ density', q), ...
            sprintf('$n_{O^{%d+}}$ [m$^{-3}$]', q), use_log_density, ...
            geom, save_plots, plot_dir, sprintf('no%d', q), fontSize);
    end

    if plot_oxygen_temp
        plotIfAvailable(R, Z, to_q, sprintf('O$^{%d+}$ temperature', q), ...
            sprintf('$T_{O^{%d+}}$ [eV]', q), false, ...
            geom, save_plots, plot_dir, sprintf('to%d', q), fontSize);
    end

    if plot_oxygen_velocity
        vopar = readNcIfExists(ncfile, sprintf('vopar_o%d', q));
        vro   = readNcIfExists(ncfile, sprintf('vro%d', q));
        vto   = readNcIfExists(ncfile, sprintf('vto%d', q));
        vzo   = readNcIfExists(ncfile, sprintf('vzo%d', q));
        vpo   = readNcIfExists(ncfile, sprintf('vpo%d', q));

        plotIfAvailable(R, Z, vopar, sprintf('O$^{%d+}$ $V_{||}$', q), '$V_{||}$ [m/s]', false, geom, save_plots, plot_dir, sprintf('vopar_o%d', q), fontSize);
        plotIfAvailable(R, Z, vro,   sprintf('O$^{%d+}$ $V_R$', q), '$V_R$ [m/s]', false, geom, save_plots, plot_dir, sprintf('vro%d', q), fontSize);
        plotIfAvailable(R, Z, vto,   sprintf('O$^{%d+}$ $V_t$', q), '$V_t$ [m/s]', false, geom, save_plots, plot_dir, sprintf('vto%d', q), fontSize);
        plotIfAvailable(R, Z, vzo,   sprintf('O$^{%d+}$ $V_Z$', q), '$V_Z$ [m/s]', false, geom, save_plots, plot_dir, sprintf('vzo%d', q), fontSize);
        plotIfAvailable(R, Z, vpo,   sprintf('O$^{%d+}$ $|V|$', q), '$|V|$ [m/s]', false, geom, save_plots, plot_dir, sprintf('vpo%d', q), fontSize);
    end
end

%% ---------------------------
%  Lineouts for ne and Te
%  ---------------------------
if plot_lineouts && ~isempty(profiles.ne) && ~isempty(profiles.te)
    plotNeTeLineouts(R, Z, profiles.ne, profiles.te, save_plots, plot_dir, fontSize);
end

%% ---------------------------
%  Experimental OMP reflectometry comparison
%  ---------------------------
if compare_exp_omp_ne && ~isempty(profiles.ne)
    plotExperimentVsSimOMPne(exp_ne_file, exp_t_min, exp_t_max, ...
        R, Z, profiles.ne, omp_z_target, geom, save_plots, plot_dir, fontSize);
end

disp('>>>> Visualization complete')

%% ========================================================================
%  Local functions
% ========================================================================

function A = readNcIfExists(ncfile, varname)
    info = ncinfo(ncfile);
    names = string({info.Variables.Name});

    if any(names == string(varname))
        A = double(ncread(ncfile, varname));
        A(~isfinite(A)) = 0;
    else
        A = [];
        fprintf('Variable not found, skipping: %s\n', varname);
    end
end

function geom = readGeometryFromNc(ncfile)
    geom = struct();

    geom.wall_r = readVectorIfExists(ncfile, 'wall_r');
    geom.wall_z = readVectorIfExists(ncfile, 'wall_z');
    geom.have_wall = ~isempty(geom.wall_r) && ~isempty(geom.wall_z);

    geom.sep1_r = readVectorIfExists(ncfile, 'sep1_r');
    geom.sep1_z = readVectorIfExists(ncfile, 'sep1_z');
    geom.have_sep1 = ~isempty(geom.sep1_r) && ~isempty(geom.sep1_z);

    geom.sep2_r = readVectorIfExists(ncfile, 'sep2_r');
    geom.sep2_z = readVectorIfExists(ncfile, 'sep2_z');
    geom.have_sep2 = ~isempty(geom.sep2_r) && ~isempty(geom.sep2_z);

    geom.antenna_centroid_r = readVectorIfExists(ncfile, 'antenna_centroid_r');
    geom.antenna_centroid_z = readVectorIfExists(ncfile, 'antenna_centroid_z');
    geom.have_antenna_centroid = ~isempty(geom.antenna_centroid_r) && ~isempty(geom.antenna_centroid_z);
end

function v = readVectorIfExists(ncfile, varname)
    info = ncinfo(ncfile);
    names = string({info.Variables.Name});

    if any(names == string(varname))
        v = double(ncread(ncfile, varname));
        v = v(:);
        v = v(isfinite(v));
    else
        v = [];
    end
end

function plotIfAvailable(R, Z, A, titleStr, cbarLabel, useLog, geom, save_plots, plot_dir, fname, fontSize)
    if isempty(A)
        return;
    end

    plot2DProfile(R, Z, A, titleStr, cbarLabel, useLog, geom, save_plots, plot_dir, fname, fontSize);
end

function plot2DProfile(R, Z, A, titleStr, cbarLabel, useLog, geom, save_plots, plot_dir, fname, fontSize)
    figure('Color','w', 'Name', titleStr);

    Aplot = A.';

    if useLog
        Aplot(Aplot <= 0) = NaN;
        h = pcolor(R, Z, Aplot);
        set(h, 'EdgeColor','none');
        set(gca, 'ColorScale','log');
    else
        h = pcolor(R, Z, Aplot);
        set(h, 'EdgeColor','none');
    end

    set(gca, 'YDir','normal');
    set(gca, 'FontName','times', 'FontSize', fontSize);

    xlabel('$R$ [m]', 'Interpreter','latex', 'FontSize', fontSize);
    ylabel('$Z$ [m]', 'Interpreter','latex', 'FontSize', fontSize);
    title(titleStr, 'Interpreter','latex');

    axis equal tight;
    c = colorbar;
    ylabel(c, cbarLabel, 'Interpreter','latex');

    hold on;

    if isfield(geom,'have_sep1') && geom.have_sep1
        plot(geom.sep1_r, geom.sep1_z, 'k--', 'LineWidth', 1.5);
    end

    if isfield(geom,'have_sep2') && geom.have_sep2
        plot(geom.sep2_r, geom.sep2_z, 'k--', 'LineWidth', 1.5);
    end

    if isfield(geom,'have_wall') && geom.have_wall
        plot(geom.wall_r, geom.wall_z, 'k', 'LineWidth', 1.6);
    end

    if isfield(geom,'have_antenna_centroid') && geom.have_antenna_centroid
        plot(geom.antenna_centroid_r, geom.antenna_centroid_z, '.', ...
            'Color', [1 0.55 0.55], 'MarkerSize', 5);
    end

    xlim([min(R) max(R)]);
    ylim([min(Z) max(Z)]);
    box on;

    if save_plots
        saveas(gcf, fullfile(plot_dir, [fname '.png']));
    end
end

function plotNeTeLineouts(R, Z, ne, te, save_plots, plot_dir, fontSize)
    z_mid_value = 0.5 * (min(Z) + max(Z));
    R_mid_value = 0.5 * (min(R) + max(R));

    [~, iz_mid] = min(abs(Z - z_mid_value));
    [~, iR_mid] = min(abs(R - R_mid_value));

    fprintf('Lineout ne/Te vs R at Z = %.6g m, index %d\n', Z(iz_mid), iz_mid);
    fprintf('Lineout ne/Te vs Z at R = %.6g m, index %d\n', R(iR_mid), iR_mid);

    figure('Color','w', 'Name','Mid-Z lineout ne and Te vs R');
    yyaxis left
    plot(R, ne(:,iz_mid), 'LineWidth', 2.0);
    ylabel('$n_e$ [m$^{-3}$]', 'Interpreter','latex');

    yyaxis right
    plot(R, te(:,iz_mid), 'LineWidth', 2.0);
    ylabel('$T_e$ [eV]', 'Interpreter','latex');

    xlabel('$R$ [m]', 'Interpreter','latex');
    title(sprintf('Mid-Z lineout at Z = %.4f m', Z(iz_mid)), 'Interpreter','latex');
    set(gca,'FontName','times','FontSize',fontSize);
    grid on;
    box on;

    if save_plots
        saveas(gcf, fullfile(plot_dir, 'lineout_midZ_ne_Te_vs_R.png'));
    end

    figure('Color','w', 'Name','Mid-R lineout ne and Te vs Z');
    yyaxis left
    plot(Z, ne(iR_mid,:), 'LineWidth', 2.0);
    ylabel('$n_e$ [m$^{-3}$]', 'Interpreter','latex');

    yyaxis right
    plot(Z, te(iR_mid,:), 'LineWidth', 2.0);
    ylabel('$T_e$ [eV]', 'Interpreter','latex');

    xlabel('$Z$ [m]', 'Interpreter','latex');
    title(sprintf('Mid-R lineout at R = %.4f m', R(iR_mid)), 'Interpreter','latex');
    set(gca,'FontName','times','FontSize',fontSize);
    grid on;
    box on;

    if save_plots
        saveas(gcf, fullfile(plot_dir, 'lineout_midR_ne_Te_vs_Z.png'));
    end
end

function plotExperimentVsSimOMPne(exp_ne_file, exp_t_min, exp_t_max, R, Z, ne, omp_z_target, geom, save_plots, plot_dir, fontSize)
    disp('>>>> Plotting experimental vs simulation OMP ne comparison')

    if ~isfile(exp_ne_file)
        warning('Experimental file not found: %s', exp_ne_file);
        return;
    end

    S = load(exp_ne_file);

    try
        t_exp = S.data.WDP.S57877.reflec.t;
        idx_t = find(t_exp >= exp_t_min & t_exp <= exp_t_max);

        if isempty(idx_t)
            warning('No experimental points found in selected time window %.3f -- %.3f s.', exp_t_min, exp_t_max);
            return;
        end

        r_exp  = S.data.WDP.S57877.reflec.position.r(:,idx_t);
        ne_exp = S.data.WDP.S57877.reflec.ne(:,idx_t);

        r_exp(~isfinite(r_exp)) = NaN;
        ne_exp(~isfinite(ne_exp) | ne_exp <= 0) = NaN;

        r_mean  = mean(r_exp , 2, 'omitnan');
        ne_mean = mean(ne_exp, 2, 'omitnan');
        ne_std  = std(ne_exp , 0, 2, 'omitnan');

        good_exp = isfinite(r_mean) & isfinite(ne_mean) & ne_mean > 0;

        [~, iz_omp] = min(abs(Z - omp_z_target));
        ne_sim = ne(:,iz_omp);

        good_sim = isfinite(R) & isfinite(ne_sim) & ne_sim > 0;

        fprintf('Using simulation OMP lineout at Z = %.6f m, index %d\n', Z(iz_omp), iz_omp);

        figure('Color','w', 'Name','Experimental vs Simulation OMP ne');
        hold on;

        semilogy(r_exp, ne_exp, '-', ...
            'Color', [0.80 0.80 0.80], ...
            'LineWidth', 0.8);

        errorbar(r_mean(good_exp), ne_mean(good_exp), ne_std(good_exp), ...
            'ko', ...
            'MarkerFaceColor','k', ...
            'LineWidth',1.5, ...
            'DisplayName','Reflectometry average');

        semilogy(R(good_sim), ne_sim(good_sim), ...
            'r-', ...
            'LineWidth',3.0, ...
            'DisplayName', sprintf('SOLEDGE/GITR input, Z = %.4f m', Z(iz_omp)));

        set(gca,'YScale','log');
        xlim([2.9 3.1]);
        ylim([1e16 1e20]);

        xlabel('$R$ [m]', 'Interpreter','latex');
        ylabel('$n_e$ [m$^{-3}$]', 'Interpreter','latex');
        title(sprintf('OMP density comparison, %.1f -- %.1f s', exp_t_min, exp_t_max), ...
            'Interpreter','latex');

        set(gca, 'FontName','times', 'FontSize', fontSize);
        grid on;
        box on;
        legend('Location','best');

        if isfield(geom,'have_wall') && geom.have_wall
            xline(min(geom.wall_r), 'k:', 'LineWidth',1.2, ...
                'Label','Wall min R', 'LabelOrientation','horizontal');
            xline(max(geom.wall_r), 'k:', 'LineWidth',1.2, ...
                'Label','Wall max R', 'LabelOrientation','horizontal');
        end

        if save_plots
            saveas(gcf, fullfile(plot_dir, 'comparison_OMP_ne_expt_vs_sim_zoom_log.png'));
        end

    catch ME
        warning('Could not process experimental reflectometry density data: %s', ME.message);
    end
end
