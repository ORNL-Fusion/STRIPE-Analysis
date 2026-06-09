close all; clear; clc;

outnc = fullfile(pwd, 'profiles_iter_single.nc');

%% === MIDPLANE PROFILE CHECK (Z ~= Z0) from ITER single-species NetCDF ===
fprintf('\n=== Visualizing MIDPLANE profiles from %s ===\n', outnc);

% --- Read back main fields ---
x = ncread(outnc,'x');   % [nX]
z = ncread(outnc,'z');   % [nZ]

ne  = ncread(outnc,'ne');    % [nX x nZ]
ni  = ncread(outnc,'ni');    % [nX x nZ]
te  = ncread(outnc,'te');    % [nX x nZ]
ti  = ncread(outnc,'ti');    % [nX x nZ]

vr  = ncread(outnc,'vr');    % [nX x nZ]
vt  = ncread(outnc,'vt');    % [nX x nZ]
vz  = ncread(outnc,'vz');    % [nX x nZ]

% Optional extras if you want them later:
% psiN = ncread(outnc,'psiN');

[nX, nZ] = size(ne);
fprintf('File grid: nX=%d, nZ=%d\n', nX, nZ);

% --- Choose midplane target ---
Z0 = 0.0;   % midplane target [m] (adjust if needed)

% --- Find closest Z index to Z0 ---
[~, iz0] = min(abs(z - Z0));
fprintf('Midplane cut: Z0=%.6f m, using nearest grid z(iz0)=%.6f m at iz0=%d\n', ...
        Z0, z(iz0), iz0);

% --- Midplane (z=const) 1D profiles vs x ---
ne_mid = ne(:, iz0);
ni_mid = ni(:, iz0);
te_mid = te(:, iz0);
ti_mid = ti(:, iz0);

vr_mid = vr(:, iz0);
vt_mid = vt(:, iz0);
vz_mid = vz(:, iz0);

Vmag_mid = sqrt(vr_mid.^2 + vt_mid.^2 + vz_mid.^2);

%% === Plot midplane profiles ===
figure('Name','Midplane profiles: ne, ni, Te, Ti');

subplot(4,1,1);
plot(x, ne_mid, 'LineWidth', 2); grid on;
ylabel('n_e [m^{-3}]');
title(sprintf('Midplane cut at z = %.6f m', z(iz0)));

subplot(4,1,2);
plot(x, ni_mid, 'LineWidth', 2); grid on;
ylabel('n_i [m^{-3}]');

subplot(4,1,3);
plot(x, te_mid, 'LineWidth', 2); grid on;
ylabel('T_e [eV]');

subplot(4,1,4);
plot(x, ti_mid, 'LineWidth', 2); grid on;
ylabel('T_i [eV]');
xlabel('x [m]');

figure('Name','Midplane flow profiles (vr, vt, vz, |v|)');

subplot(4,1,1);
plot(x, vr_mid, 'LineWidth', 2); grid on; ylabel('v_r [m/s]');
title(sprintf('Midplane flows at z = %.6f m', z(iz0)));

subplot(4,1,2);
plot(x, vt_mid, 'LineWidth', 2); grid on; ylabel('v_t [m/s]');

subplot(4,1,3);
plot(x, vz_mid, 'LineWidth', 2); grid on; ylabel('v_z [m/s]');

subplot(4,1,4);
plot(x, Vmag_mid, 'LineWidth', 2); grid on; ylabel('|v| [m/s]'); xlabel('x [m]');

fprintf('✅ Midplane profile visualization complete.\n');