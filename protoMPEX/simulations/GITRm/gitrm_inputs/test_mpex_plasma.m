%% plot_tantalum_fluxes_1D_2D.m
% Plots total and per-charge-state Tantalum particle fluxes:
%   1D radial:   Gamma_Z(r)
%   2D (z,r):    Gamma_Z(z,r)
% Uses workspace vars if available; otherwise reconstructs doc defaults.

clear; clc;

%% -------- Load from workspace if available --------
vars_ok = evalin('base','exist(''profiles'',''var'')');
if vars_ok
    P = evalin('base','profiles');
    r        = P.r_m;                   % [m]
    a        = P.meta.a_m;              % [m]
    mOrder   = P.meta.mOrder;
    SG       = P.SG;
    % Tantalum fractions (normalize for safety)
    f2 = P.tantalum.fractions.Ta2p; 
    f3 = P.tantalum.fractions.Ta3p; 
    f4 = P.tantalum.fractions.Ta4p; S=f2+f3+f4; f2=f2/S; f3=f3/S; f4=f4/S;
else
    % Rebuild doc defaults
    a = 0.02; mOrder = 12;
    r_max = 0.05; Nr = 501; r = linspace(0, r_max, Nr);
    SG = exp(-(r./a).^mOrder);
    f2=0.08; f3=0.62; f4=0.30; S=f2+f3+f4; f2=f2/S; f3=f3/S; f4=f4/S;
end

%% -------- Tantalum fluxes (doc) --------
% Total Ta flux: Gamma_Z(r) = 1e20 * SG  [m^-2 s^-1]
GammaZ_1D  = 1e20 * SG;                
% Split by charge state (same velocity, same fractions at sheath edge)
GammaZ2_1D = f2 * GammaZ_1D;
GammaZ3_1D = f3 * GammaZ_1D;
GammaZ4_1D = f4 * GammaZ_1D;

%% -------- 1D plots (radial) --------
r_cm = r*100;
figure(41); clf;
plot(r_cm, GammaZ_1D, 'k-', 'LineWidth', 1.7); hold on;
plot(r_cm, GammaZ2_1D, '--', 'LineWidth', 1.4);
plot(r_cm, GammaZ3_1D, '-.', 'LineWidth', 1.4);
plot(r_cm, GammaZ4_1D, ':',  'LineWidth', 1.8);
grid on; xlabel('r [cm]'); ylabel('\Gamma_Z [m^{-2}\,s^{-1}]');
title('Tantalum particle flux (1D radial)');
legend('\Gamma_{Z,tot}','\Gamma_{Ta^{2+}}','\Gamma_{Ta^{3+}}','\Gamma_{Ta^{4+}}','Location','northeast');

%% -------- 2D (z,r) maps with no parallel variation --------
Lz = 1.0; Nz = 401; Nr2 = 301;
z = linspace(0, Lz, Nz);
r2 = linspace(0, max(r), Nr2);
[Z,R] = meshgrid(z, r2);
SG2 = exp(-(R./a).^mOrder);

% Total and per-charge flux maps
GammaZ_2D  = 1e20 * SG2;         % [m^-2 s^-1]
GammaZ2_2D = f2 * GammaZ_2D;
GammaZ3_2D = f3 * GammaZ_2D;
GammaZ4_2D = f4 * GammaZ_2D;

% 2D plots
z_cm = z*100; r2_cm = r2*100;
figure(42); clf;
subplot(2,2,1);
imagesc(z_cm, r2_cm, GammaZ_2D); set(gca,'YDir','normal'); axis tight;
xlabel('z [cm]'); ylabel('r [cm]'); title('\Gamma_{Z,tot}(z,r) [m^{-2}\,s^{-1}]'); colorbar;

subplot(2,2,2);
imagesc(z_cm, r2_cm, GammaZ2_2D); set(gca,'YDir','normal'); axis tight;
xlabel('z [cm]'); ylabel('r [cm]'); title('\Gamma_{Ta^{2+}}(z,r)'); colorbar;

subplot(2,2,3);
imagesc(z_cm, r2_cm, GammaZ3_2D); set(gca,'YDir','normal'); axis tight;
xlabel('z [cm]'); ylabel('r [cm]'); title('\Gamma_{Ta^{3+}}(z,r)'); colorbar;

subplot(2,2,4);
imagesc(z_cm, r2_cm, GammaZ4_2D); set(gca,'YDir','normal'); axis tight;
xlabel('z [cm]'); ylabel('r [cm]'); title('\Gamma_{Ta^{4+}}(z,r)'); colorbar;
colormap(turbo);
sgtitle('Tantalum particle fluxes (2D z–r)');

%% -------- (Optional) consistency check against density & velocity --------
% If density and Ta velocity exist in the base workspace, verify Gamma = n*V
if vars_ok && isfield(P,'tantalum') && isfield(P.tantalum,'Vpar')
    try
        % Build 1D Gamma from stored n_tot and Vpar; compare to spec
        nTa_1D = P.tantalum.n_tot; VTa_1D = P.tantalum.Vpar;
        Gamma_from_nV_1D = nTa_1D .* VTa_1D;
        figure(43); clf;
        plot(r_cm, GammaZ_1D, 'k-', 'LineWidth',1.6); hold on;
        plot(r_cm, Gamma_from_nV_1D, 'r--', 'LineWidth',1.4);
        grid on; xlabel('r [cm]'); ylabel('\Gamma_Z [m^{-2}\,s^{-1}]');
        title('Consistency: \Gamma_Z(r) vs n_Z(r) V_{\parallel,Z}');
        legend('Spec: 1e20·SG', 'Computed: n_Z V_{\parallel,Z}', 'Location','northeast');
    catch
    end
end

%% -------- 1D profiles: electrons, temperatures, and Ta fluxes (separate figures) --------
r_cm = r*100;

% ---------- Helper: try common field locations ----------
getField = @(S, path) getfield(S, path{:}); %#ok<GFLD>  (used only if isfield checks pass)

% ---------- Electron density ne(r) ----------
ne_ok = false;
if vars_ok
    if isfield(P,'ne')
        ne_1D = P.ne; ne_ok = true;
    elseif isfield(P,'electron') && isfield(P.electron,'n')
        ne_1D = P.electron.n; ne_ok = true;
    elseif isfield(P,'plasma') && isfield(P.plasma,'ne')
        ne_1D = P.plasma.ne; ne_ok = true;
    elseif isfield(P,'profiles') && isfield(P.profiles,'ne')
        ne_1D = P.profiles.ne; ne_ok = true;
    end
end

if ~ne_ok
    ne0 = 1e19;            % [m^-3] default peak
    ne_1D = ne0 * SG;       % shaped like SG
end

% ---------- Electron temperature Te(r) ----------
Te_ok = false;
if vars_ok
    if isfield(P,'Te')
        Te_1D = P.Te; Te_ok = true;
    elseif isfield(P,'electron') && isfield(P.electron,'T')
        Te_1D = P.electron.T; Te_ok = true;
    elseif isfield(P,'plasma') && isfield(P.plasma,'Te')
        Te_1D = P.plasma.Te; Te_ok = true;
    end
end

if ~Te_ok
    Te0 = 30;              % [eV] default peak
    Te_1D = Te0 * (0.2 + 0.8*SG);  % keep finite edge value
end

% ---------- Ion temperature Ti(r) ----------
Ti_ok = false;
if vars_ok
    if isfield(P,'Ti')
        Ti_1D = P.Ti; Ti_ok = true;
    elseif isfield(P,'ion') && isfield(P.ion,'T')
        Ti_1D = P.ion.T; Ti_ok = true;
    elseif isfield(P,'plasma') && isfield(P.plasma,'Ti')
        Ti_1D = P.plasma.Ti; Ti_ok = true;
    end
end

if ~Ti_ok
    Ti0 = 20;              % [eV] default peak
    Ti_1D = Ti0 * (0.2 + 0.8*SG);
end

% ---------- Figure A: ne(r) ----------
figure(50); clf;
plot(r_cm, ne_1D, 'LineWidth', 1.8);
grid on;
xlabel('r [cm]');
ylabel('n_e [m^{-3}]');
title('Electron density profile (1D)');
set(gca,'FontSize',12);

% ---------- Figure B: Te(r) and Ti(r) ----------
figure(51); clf;
plot(r_cm, Te_1D, 'LineWidth', 1.8); hold on;
plot(r_cm, Ti_1D, '--', 'LineWidth', 1.8);
grid on;
xlabel('r [cm]');
ylabel('T [eV]');
title('Temperature profiles (1D)');
legend('T_e','T_i','Location','northeast');
set(gca,'FontSize',12);

% ---------- Figure C: Ta fluxes (total + charge states) ----------
figure(52); clf;
plot(r_cm, GammaZ_1D,  'k-', 'LineWidth', 1.8); hold on;
plot(r_cm, GammaZ2_1D, '--', 'LineWidth', 1.6);
plot(r_cm, GammaZ3_1D, '-.', 'LineWidth', 1.6);
plot(r_cm, GammaZ4_1D, ':',  'LineWidth', 2.0);
grid on;
xlabel('r [cm]');
ylabel('\Gamma_Z [m^{-2}\,s^{-1}]');
title('Tantalum particle flux (1D radial)');
legend('\Gamma_{Z,tot}','\Gamma_{Ta^{2+}}','\Gamma_{Ta^{3+}}','\Gamma_{Ta^{4+}}', ...
       'Location','northeast');
set(gca,'FontSize',12);

%% -------- Save (optional) --------
save('tantalum_fluxes_1D_2D.mat', 'r','r_cm','z','z_cm','GammaZ_1D','GammaZ2_1D','GammaZ3_1D','GammaZ4_1D', ...
     'GammaZ_2D','GammaZ2_2D','GammaZ3_2D','GammaZ4_2D','a','mOrder');