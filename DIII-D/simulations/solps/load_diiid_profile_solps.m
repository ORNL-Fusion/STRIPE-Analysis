%% Clear and Load
% close all;
% clear all;

% Load your extrapolated SOLPS data on GITR grid
load("extrapolated_data_200882.mat");

% If you have EFIT data, ensure you have 'g' in workspace
% e.g.,
% read_efit_data;

%% Check grid sizes
disp('Grid sizes:');
disp(size(X));
disp(size(Y));
disp(size(val_ne));

%% Replace NaNs with zeros
val_ne(isnan(val_ne)) = 0;
val_Te(isnan(val_Te)) = 0;
val_gradTi(isnan(val_gradTi)) = 0;
% val_gradTe(isnan(val_gradTe)) = 0;
val_vr(isnan(val_vr)) = 0;
val_vz(isnan(val_vz)) = 0;
val_vt(isnan(val_vt)) = 0;

%% Compute psiN safely
[psiN_flat, ~] = calc_psiN(g, X(:), Y(:), []);
psiN = reshape(psiN_flat, size(X));

%% 1) 2D Map of psiN
figure;
imagesc([min(X(:)) max(X(:))], [min(Y(:)) max(Y(:))], val_ne);
set(gca, 'YDir', 'normal');
colorbar;
colormap jet;
title('2D Map of Normalized Poloidal Flux (\psi_N)');
xlabel('R [m]');
ylabel('Z [m]');
hold on;
plot(g.lim(1,:), g.lim(2,:), 'k', 'LineWidth', 1);

%% 2) 2D Map of ne with psiN contours
figure;
imagesc([min(X(:)) max(X(:))], [min(Y(:)) max(Y(:))], val_ne);
set(gca, 'YDir', 'normal');
colorbar;
colormap turbo;
title('Electron Density with \psi_N contours');
xlabel('R [m]');
ylabel('Z [m]');
hold on;
plot(g.lim(1,:), g.lim(2,:), 'k', 'LineWidth', 1);
contour(X, Y, psiN', [0.5 0.7 0.9 1.0 1.2], 'w', 'LineWidth', 1);

%% 3) Scatter plot ne vs psiN (all grid points)
ne_flat = val_ne(:);
figure;
scatter(psiN_flat, ne_flat, 5, '.');
xlabel('\(\psi_N\)','Interpreter','latex','FontSize',18);
ylabel('\(n_e\) [m^{-3}]','Interpreter','latex','FontSize',18);
title('Electron density vs. \psi_N (all grid points)','FontSize',16);
grid on;

%% 4) Extract OMP midplane profiles
% Midplane index
Zvec = Y(:,1);
[~, iMid] = min(abs(Zvec));

% Extract along midplane
psiN_mp = psiN(iMid, :);
ne_mp   = val_ne(iMid, :);
Te_mp   = val_Te(iMid, :);
R_mp    = X(iMid, :);

% Use rmaxis as reference Rgeo
Rgeo = g.rmaxis;

% Select OMP points (R > Rgeo)
OMP_idx = R_mp > Rgeo;
psiN_OMP = psiN_mp(OMP_idx);
ne_OMP   = ne_mp(OMP_idx);
Te_OMP   = Te_mp(OMP_idx);

% Sort by psiN
[psiN_OMP_sorted, sortIdx] = sort(psiN_OMP);
ne_OMP_sorted = ne_OMP(sortIdx);
Te_OMP_sorted = Te_OMP(sortIdx);

%% 5) Plot ne vs psiN at OMP
figure;
plot(psiN_OMP_sorted, ne_OMP_sorted./1e20, 'o-', 'LineWidth', 1.5);hold on;
plot(psiN_OMP, ne_OMP./1e20, 'ro-', 'LineWidth', 1.5);
grid on;
xlabel('\(\psi_N\)','Interpreter','latex','FontSize',18);
ylabel('\(n_e\) [m^{-3}]','Interpreter','latex','FontSize',18);
title('Electron density vs. \psi_N at Outboard Midplane (OMP)','FontSize',16);

%% 6) Plot Te vs psiN at OMP
figure;
plot(psiN_OMP_sorted, Te_OMP_sorted/1000, 'o-', 'LineWidth', 1.5, 'Color', [0.85 0.33 0.1]);
grid on;
xlabel('\(\psi_N\)','Interpreter','latex','FontSize',18);
ylabel('\(T_e\) [eV]','Interpreter','latex','FontSize',18);
title('Electron temperature vs. \psi_N at Outboard Midplane (OMP)','FontSize',16);

disp('Script completed successfully.');