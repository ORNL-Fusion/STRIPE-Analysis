close all;
clearvars -global;
load("extrapolated_data.mat");

%% Preview SOLPS Data at Outer Midplane (OMP)
figure;

% Identify OMP index (middle of zS)
[~, ompIndex] = min(abs(zS)); % Assuming OMP is at the midplane

% Plot Electron Density (ne) at OMP (Log Scale)
subplot(2,3,1);
semilogy(rS, neS(:,ompIndex), 'b', 'LineWidth', 2);
xlabel('Radial Position (m)');
ylabel('n_e (m^{-3})');
title('Electron Density at OMP');
grid on;

% Plot Electron Temperature (Te) at OMP (Log Scale)
subplot(2,3,2);
semilogy(rS, TeS(:,ompIndex), 'r', 'LineWidth', 2);
xlabel('Radial Position (m)');
ylabel('T_e (eV)');
title('Electron Temperature at OMP');
grid on;

% Plot Temperature Gradient (gradTi) at OMP (Linear Scale)
subplot(2,3,3);
plot(rS, gradTiS(:,ompIndex), 'g', 'LineWidth', 2);
xlabel('Radial Position (m)');
ylabel('\nabla T_i');
title('Ion Temperature Gradient at OMP');
grid on;

% Plot Radial Velocity (vr) at OMP (Linear Scale)
subplot(2,3,4);
plot(rS, vrS(:,ompIndex), 'm', 'LineWidth', 2);
xlabel('Radial Position (m)');
ylabel('v_r (m/s)');
title('Radial Velocity at OMP');
grid on;

% Plot Poloidal Velocity (vt) at OMP (Linear Scale)
subplot(2,3,5);
plot(rS, vtS(:,ompIndex), 'c', 'LineWidth', 2);
xlabel('Radial Position (m)');
ylabel('v_t (m/s)');
title('Poloidal Velocity at OMP');
grid on;

% Plot Parallel Velocity (vz) at OMP (Linear Scale)
subplot(2,3,6);
plot(rS, vzS(:,ompIndex), 'k', 'LineWidth', 2);
xlabel('Radial Position (m)');
ylabel('v_z (m/s)');
title('Parallel Velocity at OMP');
grid on;

sgtitle('SOLPS Data Preview at Outer Midplane (OMP)');

figure;
semilogy(mpx, densityAtMidplane, 'b');  % Plot ne
hold on;
semilogy(mpx, extrapolatedne1d, 'b.'); % Extrapolated ne
semilogy(mpx, TeAtMidplane, 'r');  % Plot Te
semilogy(mpx, extrapolatedTe1d, 'r.'); % Extrapolated Te

% num_points_fine = 1000;
% [X_fine, Y_fine] = meshgrid(linspace(min(X(:)), max(X(:)), num_points_fine), ...
%                             linspace(min(Y(:)), max(Y(:)), num_points_fine));

% % Interpolate the velocities/gradTi onto the new X-Y grid
% val_vr = interp2(rS, zS, vrS', X, Y, 'linear');
% val_vz= interp2(rS, zS, vzS', X, Y, 'linear');
% val_vt= interp2(rS, zS, vtS', X, Y, 'linear');
% val_gradTi = interp2(rS, zS, gradTiS', X, Y, 'linear');

% Interpolate the velocities/gradTi onto the new X-Y grid

val_br = interp2(rS, zS, Br, X, Y, 'linear');
val_bz= interp2(rS, zS, Bz, X, Y, 'linear');
val_bt= interp2(rS, zS, Bt, X, Y, 'linear');

%% 5) Visualization of 2D Extrapolated Data
vars2D = {val_ne, val_Te, val_vr, val_vz, val_vt, val_gradTi, val_br, val_bz, val_bt };
varNames2D = {'Electron Density (neS)', 'Electron Temperature (TeS)', ...
              'Radial Velocity (vrS)', 'Poloidal Velocity (vzS)', ...
              'Toroidal Velocity (vtS)', 'Ion Temperature Gradient (gradTiS)',...
              'Radial B-field (Br)', 'Poloidal B-field (Bz)','Toroidal B-field (Bt)'};

for i = 1:length(vars2D)
    figure;
    imagesc([min(X(:)) max(X(:))], [min(Y(:)) max(Y(:))], vars2D{i});
    set(gca, 'YDir', 'normal');
    colorbar;
    set(gca, 'ColorScale', 'linear');
    title(['Extrapolated ', varNames2D{i}]);
    hold on;
    plot(g.lim(1,:), g.lim(2,:), 'r');
end

%% 6) Save Extrapolated Data
save('extrapolatedData_solps_ITER.mat');

% fnames = {'extrapolatedne.csv','extrapolatedTe.csv','extrapolatedVr.csv',...
%           'extrapolatedVz.csv','extrapolatedVt.csv','extrapolatedGradTi.csv'};
% for i = 1:length(vars2D)
%     writematrix(vars2D{i}, fnames{i});
% end
