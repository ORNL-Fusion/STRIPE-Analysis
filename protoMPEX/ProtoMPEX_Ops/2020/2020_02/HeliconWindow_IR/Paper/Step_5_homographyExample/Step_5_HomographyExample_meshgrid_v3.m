% Homography on old helicon window:
% Introduction:
% Pictures have been taken of the helicon window innner surface. For the
% setup used in this script, the datum for the window has been chosen to be
% on-axis at the axial location where the duct tape has been positioned.
% the duct tape has been used to indicate the approximate location where
% the antenna helical straps are connected to the transverse current rings

% In the current view, we are looking from the Target view towards the dump
% side

% Object's coordinate system:
% Right handed coordinate system
% Datum is on-axis at the upstream edge of window
% "x" is vertical
% "y" is horizontal
% "z" into the page and pointing towards Target

clc
close all
clearvars

figureName{1} = 'Step_5_Window_ImagePlane'; 
figureName{2} = 'Step_5_Window_UnwrappedSurface';

% Object's dimensions in physical space:
% =========================================================================
% Radius of window:
Ra = 12.4/2;
% Lenght of window:
La = 12.11*2.54;

% Load image:
% =========================================================================
t0 = tic;
d = imread('AlN_OldHeliconWindow.jpg');
rng_1 = 2700:-1:500;
rng_2 = 1300:3500;
im = double(d(rng_1,rng_2,1));
dd = double(d(rng_1,rng_2,:));
t0 = toc(t0);
disp(['Loading image: ',num2str(t0),' [s]'])

% For testing the effect of image size:
if 0
    rng = 1:600;
    im = im(rng,rng);
end

% Image dimensions:
rows = size(im,1);
cols = size(im,2);

% Assign coordinates to image plane:
% =========================================================================
x_im = linspace(-0.5,+0.5,rows); % "x" is vertical
y_im = linspace(-0.5,+0.5,cols); % "y" is vertical

%% Plot raw data:
close all

figure('color','w')
subplot(1,2,1)
ax(1) = gca;
hold on
plotType = 2;
switch plotType
    case 1
        surf(y_im,x_im,im,'LineStyle','none');
    case 2
        image(y_im,x_im,d(rng_1,rng_2,:));
end
view([0,90])
xlabel('y [mm]','Interpreter','latex','FontSize',12)
ylabel('x [mm]','Interpreter','latex','FontSize',12)
axis image
% Draw coordinate system:
line([0,0.1],[0,0  ],[300,300],'color','k','LineWidth',2)
text(0.1,0,300,'$y$','FontSize',12,'Interpreter','latex')
line([0,0  ],[0,0.1],[300,300],'color','k','LineWidth',2)
text(0,0.1,300,'$x$','FontSize',12,'Interpreter','latex')

subplot(1,2,2)
ax(2) = gca;
hold on
plotType = 1;
switch plotType
    case 1
        surf(y_im,x_im,im,'LineStyle','none');
    case 2
        image(y_im,x_im,d(rng_1,rng_2,:));
end
view([0,90])
axis image

% Draw coordinate system:
unitVectorColor = 'k';
line([0,0.1],[0,0  ],[300,300],'color',unitVectorColor,'LineWidth',2)
text(0.1,0,300,'$y$','FontSize',12,'Interpreter','latex')
line([0,0  ],[0,0.1],[300,300],'color',unitVectorColor,'LineWidth',2)
text(0,0.1,300,'$x$','FontSize',12,'Interpreter','latex')

% Calculte grid and homographic reconstruction:

% Define a 2D grid on the object's inner surface and with respect to
% object's datum:
% =========================================================================
% Azimuthal angle:
t_deg = (90 + 35 + linspace(+368,+0,4000));
t     = t_deg*pi/180;
% Along axis of window:
zo = linspace(+0,+La,1000)';
% Circumference:
yo = Ra*t;

% Create 3D grid of window's inner surface:
% =========================================================================
% Define and create vector q_s:
[tt,qz] = meshgrid(t,zo);
qs = Ra*tt;
qx = Ra*cos(tt);
qy = Ra*sin(tt);

% Offset vector L_s:
% =========================================================================
% Defines the location of the camera relative to the object's datum:
Lx = +0;
Ly = -0.45;
Lz = 22 + La;
L = [Lx, Ly, Lz]';

% Rotation angles:
% =========================================================================
% Rotation angles needed to rotate the objects's frame to the camera's frame:
% Rotation about yy:
a = +0*pi/180;
% Rotation about xx:
b = +(1.0 + 180)*pi/180;

% Create rotation matrix:
% =========================================================================
% Rotation matrix defined as: [e] = Ryx*[s]
% where [e] is camera's referene frame
% [s] is object's reference frame
Ryx = [+cos(a), +sin(a)*sin(b), +sin(a)*cos(b) ;...
       +0     , +cos(b)       , -sin(b)        ;...
       -sin(a), +cos(a)*sin(b), +cos(a)*cos(b)];
   
% Define unit vectors:
% =========================================================================
ex = [1 0 0]';
ey = [0 1 0]';
ez = [0 0 1]';

% Transform object vector q to camera frame (q -> Q subject to L and Ryx):
% =========================================================================
g = @(qx,qy,qz) Ryx'*( [qx,qy,qz]' - L );

 % Project object's grid to image plane:
 % =========================================================================
t0 = tic;
f =1.1*1.45;
xi = zeros(numel(zo),numel(yo));
yi = zeros(numel(zo),numel(yo));
zi = zeros(numel(zo),numel(yo));
for jj = 1:numel(yo)
     for kk = 1:numel(zo)
         Q = g(qx(kk,jj),qy(kk,jj),qz(kk,jj));
         xi(kk,jj) = -f*(Q'*ex)./(Q'*ez);
         yi(kk,jj) = -f*(Q'*ey)./(Q'*ez);
         zi(kk,jj) = -f*(Q'*ez)./(Q'*ez);
     end
end
t0 = toc(t0);
disp(['Transforming 3D grid into image plane grid: ',num2str(t0),' [s]'])

% World frame:
if 0
    figure; 
    hold on
    plot3(qx,qy,qz,'k.-')
    plot3(Lx,Ly,Lz,'ro')
    axis image
    xlim([-3,3]*Ra)
    ylim([-3,3]*Ra)
    grid on
    xlabel('x')
    ylabel('y')
    zlabel('z')
end 

% plot 3D grid to image plane:
% =========================================================================
gridColor = 'k';
g_dz = 90;
g_dt = 90;

figure(1)
subplot(1,2,2);
xx = -xi(1:g_dz:end,1:g_dt:end);
yy = -yi(1:g_dz:end,1:g_dt:end);
zz = -zi(1:g_dz:end,1:g_dt:end);
plot3(yy ,xx ,zz*500,gridColor)
plot3(yy',xx',zz'*500,gridColor)
line([0,1],[0,0],[0,0])
line([0,0],[0,1],[0,0])
axis image
xlabel('y [mm]','Interpreter','latex','FontSize',12)
ylabel('x [mm]','Interpreter','latex','FontSize',12)
ylim([-0.5,0.5])
xlim([-0.5,0.5])

% Final formatting:
% =========================================================================
set(ax,'FontName','times','FontSize',11)
text(ax(1),-0.5,0.6,'(a)','Interpreter','latex','FontSize',13)
text(ax(2),-0.5,0.6,'(b)','Interpreter','latex','FontSize',13)

% Save figure:
% =========================================================================
saveas(gcf,figureName{1},'tiffn')

%% Use cylindrical grid to interpolate image:
% Inputs need to be multiplied by -1 to produce reflection of image plane.
% Recall that pin-hole camera images are mirror reflected.

t0 = tic;
[Y,X] = meshgrid(y_im,x_im);
v = interp2(Y,X,im,-yi,-xi);
disp('Using interp2 function ...')
t0 = toc(t0);
disp(['Interpolating data: ',num2str(t0),' [s]'])

%% Plot data:
figure('color','w')
hold on
surf((tt-min(t))*180/pi,qz,v,'LineStyle','none')
set(gca,'XTick',[0:45:360],'XDir','reverse')
set(gca,'YTick',[0:5:La])
xlabel('Angle [deg]','Interpreter','Latex','FontSize',14)
ylabel('z [cm]','Interpreter','Latex','FontSize',14)
view([0,90])
axis tight
colormap('bone')
xlim([0,360])

% Set the aspect ratio based on geometry of window:
set(gca,'PlotBoxAspectRatio',[1 La/(2*pi*Ra) 1])

% Draw antenna straps:
% Disjointed helical strap
L(1) = line(    + [90,00] - min(t) - 35       ,[5,25] ,[300,300]);
L(2) = line(    + [90,00] - min(t) - 35  +360 ,[5,25] ,[300,300]);
% Full helical strap
L(3) = line(180 + [90,00] - min(t) - 35      ,[5,25] ,[300,300]);
L(4) = line(180 + [90,00] - min(t) - 35 +360 ,[5,25] ,[300,300]);
% Transerve strap - Dump side:
L(5) = line([0,360] - min(t)        ,[5,5]  ,[300,300]);
% Transerve strap - Target side:
L(6) = line([0,360] - min(t)        ,[25,25],[300,300]);
% HV gap:
% L(7) = line(-10 + [240,230] - min(t),[14,16],[300,300]);

% Add color and thickness:
lineWidth = 10;
set(L,'color','r','LineWidth',lineWidth)
% set(L(7),'color','w','LineWidth',lineWidth)

% Label HV and GND side:
ht(1) = text(82 ,7.5 ,300,'HV','FontSize',14,'Interpreter','latex','Color','r','EdgeColor','r');
ht(2) = text(319,23,300,'GND','FontSize',14,'Interpreter','latex','Color','r');

% Label DUMP and TARGET side:
ht(3) = text(225,2,300,'DUMP','FontSize',13,'Interpreter','latex','Color','w');
ht(4) = text(225,28,300,'TARGET','FontSize',13,'Interpreter','latex','Color','w');

% Annotation for spot near HV:
[a0,z0] = dsxy2figxy(250,11);
ha = annotation('ellipse',[a0,z0,0.1,0.13]);
set(ha,'LineStyle',':','Color','w','LineWidth',2)

% Annotation for FP locations:
% FP#1:
plot3(125,10,300,'wo','MarkerSize',5,'MarkerFaceColor','w')
hfp = text(130,10,300,'FP $\#$ 1','FontSize',11,'Interpreter','latex','Color','w');
hfp.HorizontalAlignment = 'right';
% FP#3:
plot3(85,19,300,'wo','MarkerSize',5,'MarkerFaceColor','w')
hfp = text(90,19,300,'FP $\#$ 3','FontSize',11,'Interpreter','latex','Color','w');
hfp.HorizontalAlignment = 'right';
% FP#2:
plot3(125 + 180,10,300,'wo','MarkerSize',5,'MarkerFaceColor','w')
hfp = text(130 + 180,10,300,'FP $\#$ 2','FontSize',11,'Interpreter','latex','Color','w');
hfp.HorizontalAlignment = 'right';
% FP#4:
plot3(85 + 180,19,300,'wo','MarkerSize',5,'MarkerFaceColor','w')
hfp = text(90 + 180,19,300,'FP $\#$ 4','FontSize',11,'Interpreter','latex','Color','w');
hfp.HorizontalAlignment = 'right';

% Final formatting:
% =========================================================================
set(gca,'FontName','times','FontSize',12)

% Saving figure:
% =========================================================================
saveas(gcf,figureName{2},'tiffn')