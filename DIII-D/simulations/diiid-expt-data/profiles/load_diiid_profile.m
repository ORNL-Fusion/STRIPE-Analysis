% plot_profiles_with_tanh_and_errorbars.m
% -------------------------------------------------------------------------
% Loads data, extracts tanh coefficients for n_e and T_e (handling both stored 'coeffs' or fallback y-array),
% re-evaluates both fits at raw psi points, and plots raw data with error bars,
% the original fitted profiles, and the re-evaluated fits in two stacked subplots.
% -------------------------------------------------------------------------

clear; close all;

%% 1) Load your workspace
% S    = load('profs_196154_1600_py3.mat');
S    = load('profs_200882_1800_py3.mat');
data = S.data;

%% ===== Electron density (n_e) =====
tanh_ne   = data.netanhpsi;
psi_ne_f  = tanh_ne.x(:);
ne_f      = tanh_ne.y(:);  % y-array as coefficients if 'coeffs' missing

if isfield(tanh_ne,'coeffs')
    coeffs_ne = tanh_ne.coeffs(:);
    disp('Using stored tanh coefficients for n_e');
else
    coeffs_ne = ne_f;
    disp('Using netanhpsi.y as tanh coefficients for n_e');
end

% raw density + error
if isfield(data,'nedatpsi')
    raw_ne   = data.nedatpsi;
    psi_ne_r = raw_ne.x(:);
    ne_raw   = raw_ne.y(:);
    if isfield(raw_ne,'yerr')
        ne_err = raw_ne.yerr(:);
    else
        ne_err = [];
    end
else
    psi_ne_r = [];
    ne_raw   = [];
    ne_err   = [];
end

% re-evaluate density fit at raw psi
ne_model = evaluate_tanh_fit_update(coeffs_ne, psi_ne_r, 'tanh');

%% ===== Electron temperature (T_e) =====
tanh_te   = data.tetanhpsi;
psi_te_f  = tanh_te.x(:);
te_f      = tanh_te.y(:);  % y-array as coefficients if 'coeffs' missing

if isfield(tanh_te,'coeffs')
    coeffs_te = tanh_te.coeffs(:);
    disp('Using stored tanh coefficients for T_e');
else
    coeffs_te = te_f;
    disp('Using tetanhpsi.y as tanh coefficients for T_e');
end

% raw temperature + error
if isfield(data,'tedatpsi')
    raw_te   = data.tedatpsi;
    psi_te_r = raw_te.x(:);
    te_raw   = raw_te.y(:);
    if isfield(raw_te,'yerr')
        te_err = raw_te.yerr(:);
    else
        te_err = [];
    end
else
    psi_te_r = [];
    te_raw   = [];
    te_err   = [];
end

% re-evaluate temperature fit at raw psi
te_model = evaluate_tanh_fit_update(coeffs_te, psi_te_r, 'tanh');

%% 2) Plot both in one figure with subplots
figure('Name','n_e & T_e tanh Fits','Color','w','Position',[200 200 700 600]);

% --- Subplot 1: n_e ---
subplot(2,1,1);
hold on; box on;
if ~isempty(ne_raw)
    if ~isempty(ne_err)
        errorbar(psi_ne_r, ne_raw, ne_err, 'k.', ...
                 'MarkerSize',8,'LineWidth',1,'DisplayName','n_e raw ±σ');
    else
        plot(psi_ne_r, 1e20.*ne_raw, 'k.', 'MarkerSize',10, 'DisplayName','n_e raw');
    end
end

plot(psi_ne_r, ne_model, 'b--','LineWidth',1.5, 'DisplayName','n_e re-eval');
xlabel('\psi','Interpreter','tex');
ylabel('n_e','Interpreter','tex');
title('Electron Density & tanh Fit','Interpreter','none');
legend('Location','best');
grid on;

% --- Subplot 2: T_e ---
subplot(2,1,2);
hold on; box on;
if ~isempty(te_raw)
    if ~isempty(te_err)
        errorbar(psi_te_r, te_raw, te_err, 'k.', ...
                 'MarkerSize',8,'LineWidth',1,'DisplayName','T_e raw ±σ');
    else
        plot(psi_te_r, te_raw, 'k.', 'MarkerSize',10, 'DisplayName','T_e raw');
    end
end

plot(psi_te_r, te_model, 'b--','LineWidth',1.5, 'DisplayName','T_e re-eval');
xlabel('\psi','Interpreter','tex');
ylabel('T_e','Interpreter','tex');
title('Electron Temperature & tanh Fit','Interpreter','none');
legend('Location','best');
grid on;


