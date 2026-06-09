%% plot_VDC_with_antenna_overlay.m
% Plot V_DC(phi,z) with antenna layout overlay

clc; clear; close all;

%% INPUT
csvFile = 'sheath_voltage_phi_z.csv';

% Antenna/strap layout in [phi_min phi_max z_min z_max]
% EDIT THESE VALUES to match your MPEX antenna geometry
antennaBoxes = [
    -135  -105   0.035   0.090;   % left upper strap
    45    85   0.025   0.095;   % right upper strap
    -90   -55  -0.090  -0.025;   % left lower strap
    70   115  -0.045   0.010;   % right lower strap
    ];

antennaLabels = {'Strap 1','Strap 2','Strap 3','Strap 4'};

%% READ DATA
T = readtable(csvFile);

phi = T.phi;
z   = T.z;
VDC = T.V_DC;

good = isfinite(phi) & isfinite(z) & isfinite(VDC);
phi = phi(good);
z   = z(good);
VDC = VDC(good);

%% INTERPOLATE V_DC TO PHI-Z GRID
nPhi = 350;
nZ   = 300;

phiGrid = linspace(min(phi),max(phi),nPhi);
zGrid   = linspace(min(z),max(z),nZ);

[PHI,Z] = meshgrid(phiGrid,zGrid);

VDC_grid = griddata(phi,z,VDC,PHI,Z,'natural');

%% PLOT V_DC WITH ANTENNA OVERLAY
figure('Color','w','Position',[100 100 1050 760])

imagesc(phiGrid,zGrid,VDC_grid)
set(gca,'YDir','normal')
colormap(turbo)

hcb = colorbar;
ylabel(hcb,'V_{DC} [V]','FontSize',14,'FontWeight','bold')

xlabel('Azimuthal coordinate \phi [deg]','FontSize',18,'FontWeight','bold')
ylabel('Axial coordinate z [m]','FontSize',18,'FontWeight','bold')
title('DC Sheath Voltage V_{DC}(\phi,z) with Antenna Layout', ...
    'FontSize',20,'FontWeight','bold')

hold on

for i = 1:size(antennaBoxes,1)

    ph1 = antennaBoxes(i,1);
    ph2 = antennaBoxes(i,2);
    z1  = antennaBoxes(i,3);
    z2  = antennaBoxes(i,4);

    rectangle('Position',[ph1,z1,ph2-ph1,z2-z1], ...
        'EdgeColor','w', ...
        'LineWidth',2.5, ...
        'LineStyle','--');

    text(mean([ph1 ph2]),mean([z1 z2]),antennaLabels{i}, ...
        'Color','w', ...
        'FontSize',13, ...
        'FontWeight','bold', ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','middle');
end

xlim([-180 180])
ylim([min(zGrid) max(zGrid)])

set(gca,'FontSize',16,'LineWidth',1.4)
box on
grid on

%% OPTIONAL: SAVE FIGURE
exportgraphics(gcf,'VDC_with_antenna_overlay.png','Resolution',300);