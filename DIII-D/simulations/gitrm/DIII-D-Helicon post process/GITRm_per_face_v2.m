clc; clear; close all; format long

%% --- Load GITRm per-face data ---
path1 = '../diiid-helicon/DIII-D_helicon_runs_196154/gitrm-surface2.nc';
ncid1 = netcdf.open(path1,'NOWRITE');

coords      = netcdf.getVar(ncid1, 2);
grossPlasma = netcdf.getVar(ncid1, 3);
grossDepo   = netcdf.getVar(ncid1, 4);
grossEroRaw = netcdf.getVar(ncid1, 5);

netcdf.close(ncid1);

sz = size(grossEroRaw,1);

%% --- Unpack triangle vertices ---
coords_mat = reshape(coords, [9, sz]).';

x1 = coords_mat(:,1); y1 = coords_mat(:,2); z1 = coords_mat(:,3);
x2 = coords_mat(:,4); y2 = coords_mat(:,5); z2 = coords_mat(:,6);
x3 = coords_mat(:,7); y3 = coords_mat(:,8); z3 = coords_mat(:,9);

%% --- Triangle area ---
P1 = [x1 y1 z1];
P2 = [x2 y2 z2];
P3 = [x3 y3 z3];

v1 = P2 - P1;
v2 = P3 - P1;

cp = cross(v1, v2, 2);
ElementArea = 0.5 * sqrt(sum(cp.^2,2));
ElementArea(ElementArea <= 0) = eps;

%% --- Physics scaling ---
erosion_rate=4.090184441246757e+17;
nP=5e6;
erosionPP=erosion_rate/nP;
erosionPP = erosion_rate / nP;

grossPlasmaA = erosionPP .* grossPlasma ./ ElementArea;
grossDepoA   = erosionPP .* grossDepo   ./ ElementArea;
grossEroRawA = erosionPP .* grossEroRaw ./ ElementArea;

grossEroA = grossPlasmaA + grossEroRawA;
netEroA   = grossEroA - grossDepoA;
netDepoA  = grossDepoA - grossEroRawA;

%% --- Plot settings ---
vars = { ...
    grossEroA,  'Gross Erosion'; ...
    grossDepoA, 'Gross Deposition'; ...
    netEroA,    'Net Erosion'; ...
    };

Xtri = [x1 x2 x3]';
Ytri = [y1 y2 y3]';
Ztri = [z1 z2 z3]';

%% === Fixed color limits (LINEAR) ===
cmin = 1e16;
cmax = 1e18;

figure('Name','3D Surface Fluxes (Linear Scale)', ...
       'Color','w', ...
       'Position',[100 100 1400 900]);

tiledlayout(3,1,'TileSpacing','compact','Padding','compact');

for k = 1:size(vars,1)

    nexttile

    data = vars{k,1};

    % Use magnitude for visualization (avoid negatives breaking colormap)
    data_plot = abs(data);

    % Optional clipping for stable color scaling
    data_plot(data_plot < cmin) = cmin;
    data_plot(data_plot > cmax) = cmax;

    patch(Xtri, Ytri, Ztri, data_plot', ...
        'FaceColor','flat', ...
        'EdgeColor','none', ...
        'FaceAlpha',1);

    axis equal tight
    view(35,28)
    grid on
    box on

    xlabel('X [m]')
    ylabel('Y [m]')
    zlabel('Z [m]')

    title(vars{k,2}, 'FontSize', 14, 'FontWeight','bold')

    colormap(turbo)
    clim([cmin cmax])

    cb = colorbar;
    ylabel(cb, 'particles m^{-2} s^{-1}')

    set(gca,'FontSize',14,'LineWidth',1)
    axis equal tight;

    camlight headlight
    lighting gouraud
end

sgtitle('Case 196154: Surface Fluxes', ...
        'FontSize',18,'FontWeight','bold');

%% --- Save outputs ---
writematrix(grossEroA,  'grossErosion_linear.txt');
writematrix(grossDepoA, 'grossDeposition_linear.txt');
writematrix(netEroA,    'netErosion_linear.txt');
writematrix(netDepoA,   'netDeposition_linear.txt');