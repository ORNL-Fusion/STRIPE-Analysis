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
im = double(d(2700:-1:500,1300:3500,1));
dd = double(d(2700:-1:500,1300:3500,:));
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

figure
hold on
surf(y_im,x_im,im,'LineStyle','none')
view([0,90])
xlabel('y')
ylabel('x')
axis image
% Draw coordinate system:
line([0,0.1],[0,0  ],[300,300],'color','k','LineWidth',2)
text(0.1,0,300,'$y$','FontSize',12,'Interpreter','latex')
line([0,0  ],[0,0.1],[300,300],'color','k','LineWidth',2)
text(0,0.1,300,'$x$','FontSize',12,'Interpreter','latex')

figure(2)
hold on
surf(y_im,x_im,im,'LineStyle','none')
view([0,90])
axis image
% Draw coordinate system:
line([0,0.1],[0,0  ],[300,300],'color','k','LineWidth',2)
text(0.1,0,300,'$y$','FontSize',12,'Interpreter','latex')
line([0,0  ],[0,0.1],[300,300],'color','k','LineWidth',2)
text(0,0.1,300,'$x$','FontSize',12,'Interpreter','latex')

%% Calculte grid and homographic reconstruction:

% Define a 2D grid on the object's inner surface and with respect to
% object's datum:
% =========================================================================
% Azimuthal angle:
t_deg = (180 + 50 + linspace(+0,-360,360));
t     = t_deg*pi/180;
% Along axis of window:
zo = linspace(+0,+La,100)';
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
Ly = +0.55;
Lz = -22;
L = [Lx, Ly, Lz]';

% Rotation angles:
% =========================================================================
% Rotation angles needed to rotate the objects's frame to the camera's frame:
% Rotation about yy:
a = +0*pi/180;
% Rotation about xx:
b = +1.0*pi/180;

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
figure(2)
xx = -xi(1:10:end,1:10:end);
yy = -yi(1:10:end,1:10:end);
zz = -zi(1:10:end,1:10:end);
plot3(yy ,xx ,zz*500,'k-')
plot3(yy',xx',zz'*500,'k-')
line([0,1],[0,0],[0,0])
line([0,0],[0,1],[0,0])
axis image
xlabel('y')
ylabel('x')
ylim([-0.5,0.5])
xlim([-0.5,0.5])

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
surf(fliplr(tt-min(t))*180/pi,qz,flipud(v),'LineStyle','none')
set(gca,'XTick',[0:45:360],'XDir','reverse')
set(gca,'YTick',[0:5:La])
xlabel('Angle [deg]','Interpreter','Latex','FontSize',14)
ylabel('z [cm]','Interpreter','Latex','FontSize',14)
view([0,90])
axis tight
colormap('bone')
xlim([0,360])

% Draw antenna straps:
L(1) = line(    + [90,00] - min(t) - 35  ,[5,25] ,[300,300]);
L(2) = line(180 + [90,00] - min(t) - 35   ,[5,25] ,[300,300]);
L(3) = line([0,360] - min(t)        ,[5,5]  ,[300,300]);
L(4) = line([0,360] - min(t)        ,[25,25],[300,300]);
L(5) = line(-10 + [240,230] - min(t),[14,16],[300,300]);

% Add color and thickness:
set(L,'color','r','LineWidth',10)
set(L(5),'color','w','LineWidth',10)

% Label HV and GND side:
text(230,11,300,'$GND$','FontSize',13,'Interpreter','latex','Color','r')
text(255,18,300,'$HV$','FontSize',13,'Interpreter','latex','Color','r')

% Label DUMP and TARGET side:
text(225,2,300,'$TARGET$','FontSize',13,'Interpreter','latex','Color','w')
text(225,28,300,'$DUMP$','FontSize',13,'Interpreter','latex','Color','w')

% Saving figure:
InputStructure.prompt = {['Would you like to save figure? Yes [1], No [0]']};
InputStructure.option.WindowStyle = 'normal';
saveFig = GetUserInput(InputStructure);

if saveFig
    figureName = 'figure4_Reconstruction';
    saveas(gcf,figureName,'tiffn')
end