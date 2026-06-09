% Q4a_GITRm_sheathModel1_FULL_signedPhi0_inset_logx_dOverRhoL.m
% Full MATLAB script to plot GITRm sheath model #1 potential
% (normalized to electron temperature) vs distance from surface normal.
%
% UPDATED per request:
%   - Plot distance normalized to DEUTERIUM gyroradius: x = d / rho_{L,D}
%   - Solve/plot out to 10 deuterium gyro-radii  -> x_max = 1e1
%   - X-axis is LOG scale from 1e-3 to 1e1
%
% IMPORTANT SIGN CONVENTION:
%   Physics: phi0 (wall wrt plasma) is NEGATIVE for electron sheath:
%       phi0 = -(Te/e)*ln(sqrt(mi/(2*pi*me)))
%   GITRm model #1 uses a POSITIVE "pot" = |phi0| as amplitude of the drop.
%   You can plot:
%       (A) signed:   e*phi/Te (negative in sheath)
%       (B) magnitude: |e*phi|/Te
%   For Q4a, typically plot signed potential (negative). Toggle below.

clear; close all; clc;

%% ------------------------- USER INPUTS -------------------------
Te_eV = 5;             % Electron temperature [eV]
ne    = 1.0e19;        % Electron density [m^-3]
B     = 0.5;           % Magnetic field magnitude [T]

ion = "D";             % "H", "D", "T", or "He4"
Zi  = 1;               % ion charge state for rho_L (usually 1)

angles_deg = [0 45 85];  % angle wrt surface normal [deg]

% Requested x-range in d/rhoL and log-scale limits
x_min = 1e-3;          % d/rhoL min
x_max = 1e1;           % d/rhoL max  (== 10 gyro-radii)

Npts  = 6000;          % number of log-spaced points (increase if you want)

make_E_plot = true;

% Plot signed potential (negative) or magnitude (positive)?
plot_signed_potential = true;

%% ------------------------- CONSTANTS -------------------------
e    = 1.602176634e-19;       % C
eps0 = 8.8541878128e-12;      % F/m
me   = 9.1093837015e-31;      % kg
mp   = 1.67262192369e-27;     % kg

%% ------------------------- ION MASS -------------------------
switch ion
    case "H"
        mi = 1.0 * mp;
    case "D"
        mi = 2.0 * mp;
    case "T"
        mi = 3.0 * mp;
    case "He4"
        mi = 4.0 * mp;
    otherwise
        error('Unknown ion option. Use "H","D","T","He4".');
end

%% ------------------------- DERIVED SCALES -------------------------
Te_J = Te_eV * e;

% Debye length (electron Debye length)
lambdaD = sqrt(eps0 * Te_J / (ne * e^2));

% Ion sound speed (simple): cs = sqrt(Te/mi)
cs = sqrt(Te_J / mi);

% Ion Larmor radius: rho_L = mi*cs / (Zi*e*B)
rhoL = mi * cs / (Zi * e * B)

%% ------------------------- SELF-CONSISTENT FLOATING POTENTIAL -------------------------
% Current balance gives:
%   e*phi0/Te = - ln( sqrt(mi/(2*pi*me)) )   (negative)
phi0_norm_signed = -log( sqrt(mi/(2*pi*me)) );  % e*phi0/Te  (negative)
phi0_V_signed    = phi0_norm_signed * Te_eV;    % Volts (negative)

% GITRm "pot" is the magnitude of the potential drop
pot_V = abs(phi0_V_signed);                     % Volts (positive amplitude)

%% ------------------------- GITRm fd POLYNOMIAL -------------------------
% EXACT polynomial from your code (angle must be in DEGREES)
fd_poly = @(a) ( ...
    0.98992 + 5.1220E-03.*a - 7.0040E-04.*a.^2 + ...
    3.3591E-05.*a.^3 - 8.2917E-07.*a.^4 + ...
    9.5856E-09.*a.^5 - 4.2682E-11.*a.^6 );

%% ------------------------- GRID (LOG in d/rhoL) -------------------------
% Requested: x-axis log scale from 1e-3 to 1e1, and solve to 10 rhoL
x_rho = logspace(log10(x_min), log10(x_max), Npts).';  % column vector

d  = x_rho * rhoL;         % meters
xi = d / lambdaD;          % d/lambdaD (not plotted, but sometimes useful)

%% ------------------------- COMPUTE phi(d), E(d) -------------------------
phi_V     = zeros(numel(d), numel(angles_deg));   % Volts (signed or mag depends below)
phi_norm  = zeros(size(phi_V));                   % e*phi/Te (dimensionless)
E_Vpm     = zeros(size(phi_V));                   % V/m (magnitude from expression)

for k = 1:numel(angles_deg)
    th = angles_deg(k);
    fd = fd_poly(th);

    % Clamp to [0,1] (polynomial fit can slightly overshoot)
    fd = max(0.0, min(1.0, fd));

    % GITRm model #1 potential magnitude (positive drop)
    phi_mag = pot_V .* ( fd .* exp(-d./(2.0*lambdaD)) + (1.0-fd) .* exp(-d./rhoL) );

    % Signed potential: wall is negative relative to plasma
    if plot_signed_potential
        phi_V(:,k) = -phi_mag;
    else
        phi_V(:,k) = phi_mag;
    end

    % Corresponding E-field magnitude from GITRm (as implemented)
    E_Vpm(:,k) = pot_V .* ( fd./(2.0*lambdaD) .* exp(-d./(2.0*lambdaD)) + ...
                            (1.0-fd)./rhoL      .* exp(-d./rhoL) );

    % Normalize: e*phi/Te = phi(Volts)/Te(Volts) = phi/Te_eV
    phi_norm(:,k) = phi_V(:,k) ./ Te_eV;
end

%% ------------------------- PRINT SUMMARY -------------------------
fprintf('--- GITRm sheath model #1 (inputs/derived) ---\n');
fprintf('Te = %.3g eV, ne = %.3g m^-3, B = %.3g T, ion = %s (Zi=%d)\n', Te_eV, ne, B, ion, Zi);
fprintf('lambdaD = %.6g m\n', lambdaD);
fprintf('rhoL    = %.6g m\n', rhoL);
fprintf('rhoL/lambdaD = %.6g\n', rhoL/lambdaD);

fprintf('\n--- Floating potential (self-consistent) ---\n');
fprintf('Signed:  e*phi0/Te = %.6g  (should be negative)\n', phi0_norm_signed);
fprintf('Signed:  phi0 = %.6g V\n', phi0_V_signed);
fprintf('GITRm uses pot = |phi0| = %.6g V  (positive amplitude)\n', pot_V);

fprintf('\n--- Plot domain ---\n');
fprintf('x = d/rhoL from %.3e to %.3e (i.e., up to %.3g gyro-radii)\n', x_min, x_max, x_max);
fprintf('Corresponding d from %.3e m to %.3e m\n', min(d), max(d));

fprintf('\n--- fd(theta) (theta in degrees, wrt surface normal) ---\n');
for k = 1:numel(angles_deg)
    th = angles_deg(k);
    fprintf('theta = %g deg -> fd = %.6f\n', th, fd_poly(th));
end

%% ------------------------- FIGURE 1: Q4a plot (log x) -------------------------
hFig = figure('Color','w','Units','inches','Position',[1 1 7.6 4.9]);
semilogx(x_rho, phi_norm, 'LineWidth', 1.8);
grid on; box on;

xlim([x_min x_max]);
xlabel('Distance from surface (normal),  d / \rho_{L,D}','FontSize',11);

if plot_signed_potential
    ylabel('Electric potential,  e\phi / T_e  (signed)','FontSize',11);
else
    ylabel('Electric potential magnitude,  |e\phi| / T_e','FontSize',11);
end

legend(arrayfun(@(t)sprintf('\\theta = %g^\\circ', t), angles_deg, 'UniformOutput', false), ...
       'Location','northeast');
title('GITRm sheath model','FontSize',12);

% Annotation strings (equations + values)
lambdaD_str = sprintf('$$\\lambda_D = %.3e\\ \\mathrm{m}$$', lambdaD);
rhoL_str    = sprintf('$$\\rho_{L,D} = %.3e\\ \\mathrm{m}$$', rhoL);

phi0_eqn_str = '$$\phi_0=-\frac{T_e}{e}\ln\!\left(\sqrt{\frac{m_i}{2\pi m_e}}\right),\ \ \frac{e\phi_0}{T_e}=-\ln\!\left(\sqrt{\frac{m_i}{2\pi m_e}}\right)$$';
phi0_val_str = sprintf('$$\\frac{e\\phi_0}{T_e}=%.4g,\\quad \\phi_0=%.4g\\ \\mathrm{V},\\quad \\mathrm{pot}=|\\phi_0|=%.4g\\ \\mathrm{V}$$', ...
                       phi0_norm_signed, phi0_V_signed, pot_V);

phi_model_str = '$$\phi(d)=\phi_0\Big[f_d(\theta)e^{-d/(2\lambda_D)}+(1-f_d(\theta))e^{-d/\rho_L}\Big]$$';

if plot_signed_potential
    sign_note = '$$\mathrm{(signed\ sheath:\ wall\ negative)}$$';
else
    sign_note = '$$\mathrm{(magnitude\ only)}$$';
end

fd_poly_str1 = '$$f_d(\theta)=0.98992+5.122{\times}10^{-3}\theta-7.004{\times}10^{-4}\theta^2+3.3591{\times}10^{-5}\theta^3$$';
fd_poly_str2 = '$$\qquad\qquad\ -8.2917{\times}10^{-7}\theta^4+9.5856{\times}10^{-9}\theta^5-4.2682{\times}10^{-11}\theta^6\ \ (\theta\ \mathrm{deg})$$';

eqn = {phi_model_str, sign_note, fd_poly_str1, fd_poly_str2, phi0_eqn_str, phi0_val_str, lambdaD_str, rhoL_str};

annotation('textbox',[0.05 0.44 0.90 0.46], ...
    'String', eqn, ...
    'Interpreter','latex', ...
    'FontSize',10, ...
    'EdgeColor','none', ...
    'HorizontalAlignment','center');

%% ------------------ INSET: fd polynomial vs cos(theta) ------------------
insetPos = [0.62 0.62 0.28 0.28]; % [left bottom width height] in normalized units
axInset = axes('Position', insetPos);

theta = linspace(0,90,500);
fd_git = fd_poly(theta);
fd_cos = cosd(theta);

plot(axInset, theta, fd_git, 'LineWidth', 1.6); hold(axInset,'on');
plot(axInset, theta, fd_cos, '--', 'LineWidth', 1.3); hold(axInset,'off');
grid(axInset,'on');
xlabel(axInset, '\theta (deg)', 'FontSize',8);
ylabel(axInset, 'f_D', 'FontSize',8);
xlim(axInset, [0 90]);
ylim(axInset, [0 1.05]);
set(axInset, 'FontSize',8, 'Box','on', 'Color','white');
legend(axInset, {'GITRm poly','cos(\theta)'}, 'FontSize',7, 'Location','southwest');

saveas(hFig, 'Figure1_Q4a_GITRm_model1_phiNorm_vs_dOverRhoL_logx_with_inset_signedPhi0.png');
fprintf('\nSaved: Figure1_Q4a_GITRm_model1_phiNorm_vs_dOverRhoL_logx_with_inset_signedPhi0.png\n');

%% ------------------------- OPTIONAL: Plot E-field magnitude -------------------------
if make_E_plot
    hFig2 = figure('Color','w','Units','inches','Position',[1 1 7.6 4.6]);
    semilogx(x_rho, E_Vpm, 'LineWidth', 1.8);
    grid on; box on;
    xlim([x_min x_max]);
    xlabel('Distance from surface (normal),  d / \rho_{L,D}','FontSize',11);
    ylabel('Sheath electric field magnitude |E| [V/m]','FontSize',11);
    legend(arrayfun(@(t)sprintf('\\theta = %g^\\circ', t), angles_deg, 'UniformOutput', false), ...
           'Location','northeast');
    title('GITRm sheath model #1 electric field magnitude','FontSize',12);

    saveas(hFig2, 'Figure2_GITRm_model1_Efield_vs_dOverRhoL_logx.png');
    fprintf('Saved: Figure2_GITRm_model1_Efield_vs_dOverRhoL_logx.png\n');
end