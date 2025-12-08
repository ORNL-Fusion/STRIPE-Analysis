close all;
clearvars -global;
load("extrapolated_data.mat");


%% SOLPS ITER Data on GITR Grid

%% Updating SOLPS profileS in GITR
% ================================

%%Modify variable here
% -------------------
% % 
% val_br=Br;
% val_bz=Bz;
% val_bt=Bt;


vars2D = {val_ne, val_Te, val_vr, val_vz, val_vt, val_gradTi};
varNames2D = {'Electron Density (neS)', 'Electron Temperature (TeS)', ...
              'Radial Velocity (vrS)', 'Poloidal Velocity (vzS)', ...
              'Toroidal Velocity (vtS)', 'Ion Temperature Gradient (gradTiS)'};

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

val_ne(isnan(val_ne))=0;
val_Te(isnan(val_Te))=0;
val_gradTi(isnan(val_gradTi))=0;

val_vr(isnan(val_vr))=0;
val_vz(isnan(val_vz))=0;
val_vt(isnan(val_vt))=0;

