%% get genray.nc file
clear
close all
clc

saveFig = 1;

%% Load LUFS data:
% =========================================================================
LUFS = load('step_0a_Profile_A_LUFS.mat');
for nn = 1:numel(LUFS.z_resLayer)
    rng = find(LUFS.r_resLayer{nn}{1} > 1);    
    LUFS.r_resLayer{nn}{1}(rng) = [];
    LUFS.z_resLayer{nn}{1}(rng) = [];
end
%% Load GENRAY-C data:
% =========================================================================
% Select dataset:
dataset = 'genray_2200.0_1630.0_1.0_0.8_-30.0_shifted.nc';

% File name:
filename = ['genray-c/',dataset];

center_distance = 0; %3.22 for Proto-MPEX, 3.4 for ICH
center_distance = 3.27; %3.22 for Proto-MPEX, 3.4 for ICH

% Extract data:
% 0-d data:
% -------------------------------------------------------------------------
freq =  ncread(filename, 'freqcy');
B_cyl = freq/6.5e9/2; %second harmonic

% 1-D data:
% -------------------------------------------------------------------------
B_x = ncread(filename, 'eqdsk_x');
B_z = ncread(filename, 'eqdsk_z') + center_distance;
B = ncread(filename, 'bmodprofxz');
ne = ncread(filename, 'densprofxz')*1e6;
ray_x = ncread(filename, 'wx') /100.; %m
ray_y = ncread(filename, 'wy') /100.; %m
ray_z = ncread(filename, 'wz') /100. + center_distance; %m
ray_dist = ncread(filename, 'ws')/100; %m 
ray_pwr = ncread(filename, 'delpwr')*1e-7; %W
ray_Ex = ncread(filename, 'cwexde'); %normalized polarization
ray_Ex = ray_Ex(:,:,1).^2 + ray_Ex(:,:,2).^2;
ray_Ey = ncread(filename, 'cweyde');
ray_Ey = ray_Ey(:,:,1).^2 + ray_Ey(:,:,2).^2;
ray_Ez = ncread(filename, 'cwezde');
ray_Ez = ray_Ez(:,:,1).^2 + ray_Ez(:,:,2).^2;
ray_Bx = ncread(filename, 'sb_x')/1e4;
ray_By = ncread(filename, 'sb_y')/1e4;
ray_Bz = ncread(filename, 'sb_z')/1e4;
ray_B = ncread(filename, 'sbtot')/1e4;
ray_ne = ncread(filename, 'sene')*1e6;
ray_Te = ncread(filename, 'ste')*1e3;
ray_npar = ncread(filename, 'wnpar');
ray_nperp = ncread(filename, 'wnper');
ray_n = sqrt(ray_npar.^2 + ray_nperp.^2);
ray_num = numel(ray_x(1,:));
rho_power = ncread(filename, 'rho_bin_center');
power = ncread(filename, 'powden')*1e-11;
power_e = ncread(filename, 'powden_e')*1e-11;
power_col = ncread(filename, 'powden_cl')*1e-11;
powertot_e = ncread(filename, 'powtot_e')*1e-7; %W
powertol_col = ncread(filename, 'powtot_cl')*1e-7; %W
den_x = ncread(filename,'w_x_densprof_nc');
den_x = den_x(ceil(numel(den_x)/2):end);
den = ncread(filename, 'w_dens_vs_x_nc')*1e6/1e19;
den = den(ceil(numel(den)/2):end);
temp = ncread(filename, 'w_temp_vs_x_nc')*1e3;
temp = temp(ceil(numel(temp)/2):end);
xscan = ncread(filename, 'xscan');
rhoscan = ncread(filename, 'rhoscan');
nscan = numel(rhoscan);
power_x = interp1(rhoscan(1:nscan/2), xscan(1:nscan/2), rho_power);
power_x(1)= 0;

% 2-D data:
% -------------------------------------------------------------------------
rgrid = ncread(filename, 'Rgrid');
zgrid = ncread(filename, 'Zgrid') + center_distance;
power_2D = ncread(filename, 'spwr_rz_e');
B_2D = ncread(filename, 'bmodprofxz');
den_1D_l = zeros(530,1);
den_1D_h = den_1D_l;
UH_1D_l = den_1D_h;
UH_1D_h = den_1D_h;
for ii = 1:numel(den_1D_l)
    den_1D_l(ii) = interp1(ne(1:530,ii), B_x(1:530), .98e19);
    den_1D_h(ii) = interp1(ne(531:1060,ii), B_x(531:1060), .98e19);
    UH_1D_l(ii) = interp1(sqrt(81.* ne(1:530,ii) + 784e18* B(1:530,ii).^2), ...
        B_x(1:530), 28e9);
    UH_1D_h(ii) = interp1(sqrt(81.*ne(531:1060,ii) + 784e18* B(531:1060,ii).^2), ...
        B_x(531:1060), 28e9);
end

% write to file:
% -------------------------------------------------------------------------
csvwrite(strcat(erase(filename, '.nc'), '.csv'), [ray_x, ray_y, ray_z, ...
    ray_dist, ray_ne, ray_Te ...
    ,ray_pwr, ray_Bx,ray_By,ray_Bz, ray_npar, ray_nperp, ray_Ex, ...
    ray_Ey, ray_Ez]);

%% Plot: Schematic with ray and B profile:
% =========================================================================
close all

% On-axis magnetic field:
B_arr = B(530,:);
index = find(B_arr > B_cyl);

% Switches:
select_Genray_B = 0;

% Formatting:
fontSize.label = 14;
fontSize.axes = 13;
zMax_plot = 3.66;
zMin_plot = 2.88;

figure('color','w')
box on
set(gcf,'Position',[20    89   682   534])

% Magnetic field:
% -------------------------------------------------------------------------
ax1 = subplot(10,1,[1:3]);
hold on

if select_Genray_B
    plot(ax1, B_z, B_arr, 'k', 'LineWidth', 2); 
else
    plot(LUFS.z1Db,LUFS.B_onAxis, 'k', 'LineWidth', 2);
end

% Formatting:
box on
set(gca,'fontSize',fontSize.axes,'fontName','Times')
% title('Magnetic field [T]','interpreter','latex','fontSize',fontSize.label);
xlabel(ax1, 'z [m]','interpreter','latex','fontSize',fontSize.label);
ylabel(ax1, 'B [T]','interpreter','latex','fontSize',fontSize.label);
xlim([zMin_plot,zMax_plot]);
ylim([0, max(B_arr) * 1.2]); 

% Geometry and ray:
% -------------------------------------------------------------------------
ax2 = subplot(10,1,[5:10]);
hold on;

try 
    clear hLayer
end
hLayer(1) = plot(B_z, den_1D_l, 'k');
hLayer(2) = plot(B_z, den_1D_h, 'k');
hLayer(3) = plot(B_z, UH_1D_l, 'b');
hLayer(4) = plot(B_z, UH_1D_h, 'b');
set(hLayer,'LineWidth',2)

for ii = 1:ray_num
    rayindex = max(find(ray_x(:,ii) ~= 0));
    surface(ax2, [ray_z(1:rayindex,ii)'; ray_z(1:rayindex,ii)'],  [-ray_x(1:rayindex,ii)'; -ray_x(1:rayindex,ii)'], ...
    [ray_z(1:rayindex,ii)'; ray_z(1:rayindex,ii)'], [ray_pwr(1:rayindex,ii)'./1e3; ray_pwr(1:rayindex,ii)'./1e3], ...
    'facecol', 'no', 'edgecol', 'interp', 'linew', 3);
end
colormap(flipud(hot));

% Draw vacuum vessel:
z_offset = 0*3.29;
plot(LUFS.vessel_1_U.z-z_offset,+LUFS.vessel_1_U.r,'k','lineWidth',1)
plot(LUFS.vessel_1_U.z-z_offset,-LUFS.vessel_1_U.r,'k','lineWidth',1)

% Draw ICH window:
w = diff(LUFS.ichSleeve.z(2:3));
h = diff(LUFS.heliconWindow.r(3:end));
x = LUFS.ichSleeve.z(1);
y = LUFS.ichSleeve.r;
pos = [x +y(2) w h];
hrect = rectangle('Position',pos);
hrect.FaceColor = 'c';
hrect.EdgeColor = 'c';
pos = [x -y(2)-h w h];
hrect = rectangle('Position',pos);
hrect.FaceColor = 'c';
hrect.EdgeColor = 'c';

% Electron cyclotron resonant surfaces:
lineColor = {'r.','g.'};
for nn = 1:numel(LUFS.z_resLayer)
    hResU(nn) = plot(LUFS.z_resLayer{nn}{1}-z_offset,+LUFS.r_resLayer{nn}{1},lineColor{nn});
    hResL(nn) = plot(LUFS.z_resLayer{nn}{1}-z_offset,-LUFS.r_resLayer{nn}{1},lineColor{nn});
end
set([hResU,hResL],'MarkerSize',9)

% Formatting:
box on
set(gca,'fontSize',fontSize.axes,'fontName','Times')
% caxis(ax2, [0,1e2]);
ylabel(ax2, 'r [m]','interpreter','latex','fontSize',fontSize.label);
xlabel(ax2, 'z [m]','interpreter','latex','fontSize',fontSize.label); 
ylim([-.11, .05]);
xlim([zMin_plot,zMax_plot]);
ylim([-1,1]*0.12)

if 0
    ylabel(colorbar, 'Power (%) left in Ray');    
end

% Arrows:
% Text arrows fields:
fields.HorizontalAlignment = 'center';
fields.Color = 'k';
fields.FontSize = 12;
fields.Interpreter = 'Latex';
fields.HeadLength = 5;
fields.HeadWidth = 5;
fields.HeadStyle = 'vback2';

% UH resonance layer:
fields.String = 'UHR';
x =      [3    ,3    ];
y =      [-0.06,-0.025];
hta = myTextArrow(ax2,x,y,fields);

% 2nd harmonic resonance layer:
fields.String = '2$^{nd}$';
fields.FontSize = 14;
x =      [3.1 ,3.15 ];
y =      [-0.08,-0.08];
fields.Interpreter = 'latex';
hta = myTextArrow(ax2,x,y,fields);

% 3rd harmonic resonance layer:
fields.String = '3$^{rd}$';
x =      [3.25 ,3.21 ];
y =      [0.08,0.067];
fields.Interpreter = 'latex';
hta = myTextArrow(ax2,x,y,fields);

fields.String = [];
x =      [3.275,3.318 ];
y =      [0.08,0.067];
fields.Interpreter = 'latex';
hta = myTextArrow(ax2,x,y,fields);

% ICH window:
fields.String = 'ICH window';
fields.FontSize = 12;
x =      [3.5  ,3.5   ];
y =      [-0.075,-0.05];
hta = myTextArrow(ax2,x,y,fields);

% Zoom-in region for ray;
% z_zoom = [3.24,3.255];
% r_zoom = [-0.045,-0.024];
z_zoom = [+3.229,+3.252] - 0.06;
r_zoom = [-0.040,-0.022] + 0.006;

hZ(1) = line(z_zoom         ,[1,1]*r_zoom(1));
hZ(2) = line(z_zoom         ,[1,1]*r_zoom(2));
hZ(3) = line([1,1]*z_zoom(1),r_zoom         );
hZ(4) = line([1,1]*z_zoom(2),r_zoom         );

set(hZ,'Color','k','lineStyle','-','LineWidth',1)

if 1
    % Beam direction:
    x =      [0.490713587487781,0.470185728250244];
    y =      [0.113607990012484,0.176029962546817];
    hArr = annotation('textarrow',x,y);
    hArr.LineWidth = 4;
    hArr.Color = 'k';
    hArr.HeadLength = 15;
    hArr.HeadWidth = 15;
end

% Save figure:
% =========================================================================
if saveFig
    figureName = ['Step_14a_RayTracing_Ray2D_Bprofile'];
    % PDF figure:
    exportgraphics(gcf,[figureName,'.pdf'],'Resolution',600,'ContentType', 'vector') 

    % TIFF figure:
    exportgraphics(gcf,[figureName,'.tiff'],'Resolution',600) 
end

%% Plot: Zoom-in ray:
% =========================================================================

fontSize.legend = 11;
fontSize.label  = 14;
fontSize.axes   = 11;

figure('color','w')
% set(gcf,'Position',[710    89   400   534])
set(gcf,'Position',[356   216   441   371])

% Ray:
% -------------------------------------------------------------------------
ax1 = gca;
%subplot(1,2,1);
hold on;

clear hLayer
hLayer(1) = plot(B_z, den_1D_l, 'k');
hLayer(2) = plot(B_z, den_1D_h, 'k');
hLayer(3) = plot(B_z, UH_1D_l, 'b');
hLayer(4) = plot(B_z, UH_1D_h, 'b');
set(hLayer,'LineWidth',2)

for ii = 1:ray_num
    rayindex = max(find(ray_x(:,ii) ~= 0));
    term.z = [ray_z(1:rayindex,ii)'; ray_z(1:rayindex,ii)'];
    term.x = [-ray_x(1:rayindex,ii)'; -ray_x(1:rayindex,ii)'];
    term.z = [ray_z(1:rayindex,ii)'; ray_z(1:rayindex,ii)'];
    term.u = [ray_pwr(1:rayindex,ii)'./1e3; ray_pwr(1:rayindex,ii)'./1e3];
    surface(ax1, term.z , term.x,term.z , term.u, ...
    'facecol', 'no', 'edgecol', 'interp', 'linew', 3);
end
colormap(flipud(hot));

% Formatting:
box on
set(gca,'fontSize',fontSize.axes,'fontName','Times')
caxis(ax1, [0,1e2]);
ylabel(ax1, 'r [m]','interpreter','latex','fontSize',fontSize.label);
xlabel(ax1, 'z [m]','interpreter','latex','fontSize',fontSize.label);
axis equal
xlim(z_zoom);
ylim(r_zoom);

if  0
    lineColor = {'r.','g.'};
    for nn = 1:numel(LUFS.z_resLayer)
        hResU(nn) = plot(LUFS.z_resLayer{nn}{1}-z_offset,+LUFS.r_resLayer{nn}{1},lineColor{nn});
        hResL(nn) = plot(LUFS.z_resLayer{nn}{1}-z_offset,-LUFS.r_resLayer{nn}{1},lineColor{nn});
    end
    set([hResU,hResL],'MarkerSize',1,'LineStyle','-','LineWidth',2)
end

if 1
    ylabel(colorbar, 'Power (%) left in Ray','interpreter','latex','fontSize',fontSize.label);    
end

% Arrows:
% Text arrows fields:
fields.HorizontalAlignment = 'center';
fields.Color = 'k';
fields.FontSize = 12;
fields.Interpreter = 'Latex';
fields.HeadLength = 5;
fields.HeadWidth = 5;
fields.HeadStyle = 'vback2';

if 0
    % Beam direction:
    x =      [0.718820861678004,0.63718820861678];
    y =      [0.22811051212938,0.383647798742138];
    hArr = annotation('textarrow',x,y);
    hArr.LineWidth = 5;
    hArr.Color = 'k';
    hArr.HeadLength = 20;
    hArr.HeadWidth = 20;
end

hLeg = legend([hLayer([1,3])],'O-mode','UHR');
set(hLeg,'Location','southwest','interpreter','latex','fontsize',fontSize.legend)

% Save figure:
% =========================================================================
if saveFig
    figureName = ['Step_14a_RayTracing_Zpom'];
    % PDF figure:
    exportgraphics(gcf,[figureName,'.pdf'],'Resolution',600,'ContentType', 'vector') 

    % TIFF figure:
    exportgraphics(gcf,[figureName,'.tiff'],'Resolution',600) 
end

%% Plot: Density profile:
% =========================================================================

fontSize.legend = 11;
fontSize.label  = 14;

yMax_plot = 3e19;

figure('color','w')
set(gcf,'Position',[423   300   377   292])

ax2 = gca;
%subplot(1,2,2);
box on
hold on
set(gca,'XDir','reverse')

% Absorption profile:
hP(1) = plot(-power_x,power_e*0.8*yMax_plot/max(power_col), 'k','lineWidth',3);
hP(2) = plot(-power_x,power_col*0.8*yMax_plot/max(power_col),'r','lineWidth',3);

% layers:
rng_layer = find(B_z > 3.24,1);
hLine(1) = line([1,1]*den_1D_l(rng_layer),[0,5e19],'lineWidth',2,'color','k');
hLine(2) = line([1,1]*UH_1D_l(rng_layer),[0,5e19],'lineWidth',2,'color','b');

hR(1) = plot(-den_x,den*1e19, 'k--', 'LineWidth', 2);

set(gca,'fontSize',fontSize.axes,'fontName','Times')
ylabel(ax2, '$n_e$ [m$^{-3}$]','interpreter','latex','fontSize',fontSize.label);
xlabel(ax2, 'r [m]','interpreter','latex','fontSize',fontSize.label);
xlim(-[0.035,0])
ylim([0,yMax_plot])

clear hLeg
hLeg = legend([hP(1),hP(2),hR(1),hLine],'Resonant','Collisional','$n_e$','O-mode','UHR');
set(hLeg,'Location','southwest','interpreter','latex','fontsize',fontSize.legend,'Box','off')

% Save figure:
% =========================================================================
if saveFig
    figureName = ['Step_14a_RayTracing_densityProfile'];
    % PDF figure:
    exportgraphics(gcf,[figureName,'.pdf'],'Resolution',600,'ContentType', 'vector') 

    % TIFF figure:
    exportgraphics(gcf,[figureName,'.tiff'],'Resolution',600) 
end

%% Plot: dispersion relation

fontSize.label(1) = 18;
fontSize.label(2) = 16;
fontSize.title = 15;
fontSize.axes  = 13;

fontWeight = 'normal';

% Angular frequency:
w = 2*pi*28e9; 

% X term:
X = w_pe(ray_ne)/w;

% Y term:
Y = w_ce(ray_B)/w;

% Refractive index:
n_perp = ray_nperp;

figure('color','w')
hold on
box on
set(gca,'fontSize',fontSize.axes,'fontName','Times')

for ii = 1:ray_num
    rayindex = max(find(ray_x(:,ii) ~= 0));
    term.z = [X(1:rayindex,ii)'; X(1:rayindex,ii)'];
    term.x = [n_perp(1:rayindex,ii)'; n_perp(1:rayindex,ii)'];
    term.z = [X(1:rayindex,ii)'; X(1:rayindex,ii)'];
    term.u = [ray_pwr(1:rayindex,ii)'./1e3; ray_pwr(1:rayindex,ii)'./1e3];
       
        surface(term.z , term.x,term.z , term.u, ...
    'facecol', 'no', 'edgecol', 'interp', 'linew', 5);
end
colormap(flipud(hot));

% UHR layer:
X_UH = sqrt( w_pe(ray_ne).^2 + w_ce(ray_B).^2 )./w;
i_UH = find(X_UH >= 1,1);
X_UH = X(i_UH,1);
h_UH = line([1,1]*X_UH,[1e-2,1e3]);

% O-mode:
X_O = sqrt( w_pe(ray_ne).^2 )./w;
i_O = find(X_O >= 1,1);
X_O = X(i_O,1);
h_O = line([1,1]*X_O,[1e-2,1e3]);

% Label:
title(['Incidence angle $\approx$ 30','$^{\circ}$, ','$ \omega_{ce}/\omega$ $\approx$  0.4'],'fontSize',fontSize.title,'interpreter','latex');
xlabel('$ \omega_{pe}/\omega$','interpreter','latex','fontSize',fontSize.label(1))
ylabel('$ n_{\perp} $','interpreter','latex','fontSize',fontSize.label(1))
hT{1} = text(0.88,0.012,'UHR','fontSize',fontSize.label(2),'interpreter','latex','color','bl','rotation',90,'FontWeight',fontWeight);
hT{2} = text(1.05,5.9401,'O-cutoff','fontSize',fontSize.label(2),'interpreter','latex','color','k','rotation',90,'FontWeight',fontWeight);
ylabel(colorbar, 'Power (%) left in Ray','interpreter','latex','fontSize',fontSize.label(2));    

% Formatting:
set(gca,'Yscale','log')
set(h_UH,'lineWidth',2,'LineStyle','--','Color','b')
set(h_O ,'lineWidth',2,'LineStyle','--','Color','k') 

ylim([1e-2,1e3])
xlim([0.3,1.4])

% Arrows:
% Text arrows fields:
fields.HorizontalAlignment = 'center';
fields.Color = 'k';
fields.FontSize = 14;
fields.Interpreter = 'Latex';
fields.HeadLength = 5;
fields.HeadWidth = 5;
fields.HeadStyle = 'vback2';

% O-mode branch:
fields.String = 'O-mode';
x =      [0.6,0.6  ];
y =      [450,375];
hta = myTextArrow(gca,x,y,fields);

% SX-mode branch:
fields.String = 'SX-mode';
x =      [1.1,1.05];
y =      [300,300];
hta = myTextArrow(gca,x,y,fields);

% B-mode branch:
fields.String = 'B-mode';
x =      [0.7,0.8];
y =      [0.8e3,0.8e3];
hta = myTextArrow(gca,x,y,fields);

% Save figure:
% =========================================================================
if saveFig
    figureName = ['Step_14a_RayTracing_Absorption_Dispersion'];
    % PDF figure:
    exportgraphics(gcf,[figureName,'.pdf'],'Resolution',600,'ContentType', 'vector') 

    % TIFF figure:
    exportgraphics(gcf,[figureName,'.tiff'],'Resolution',600) 
end