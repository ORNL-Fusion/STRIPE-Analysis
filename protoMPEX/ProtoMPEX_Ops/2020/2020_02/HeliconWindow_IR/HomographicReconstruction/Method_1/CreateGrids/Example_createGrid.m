% Objective:
% Apply homographic reconstruction to the IR data taken in Proto-MPEX
% Stitch raw images together 

%% Process: 
% 1- Read the dataset spreadsheet
% 2- Assemble the shot series
% 3- Gather data into structures
% 4- Calculate temperature difference and mirror image

% Notes:
% In this code, we do not extract the data directly from the .seq files; we
% obtained the data from previosly extracted data stored in .mat files.
% The effect of the mirror is applied in step 4

clc
clear all
close all

%% 1- Read the dataset spreadsheet
% =========================================================================
spreadsheetName = '\HeliconWindowIR_2020_02_XPs.xlsx';
home = cd;
cd ..
cd ..
addr = pwd;
cd(home)
T = readtable([addr,spreadsheetName],'Sheet',1);
groups = unique(T.Group);

%% 2- Assemble the shot series:
% =========================================================================

% MPEX-Limit, 140 kW:
shotSeries{1} = [29077,29082,29049,29117,29120,29128];
rfPower{1}    = [142  ,142  ,142  ,135  ,150  ,141  ];

shotSeries{2} = [29076,29081,29047,29115,29122,29127];
rfPower{2}    = [126  ,126  ,119  ,120  ,116  ,118  ];

% MPEX-Limit, 100 kW:
shotSeries{3} = [29075,29080,29066,29114,29123,29126];
rfPower{3}    = [107  ,102  ,105  ,100  ,97   ,96   ];

% MPEX-Limit, 80 kW:
shotSeries{4} = [29070,29079,29067,29113,29124,29125];
rfPower{4}    = [83   ,88   ,83   ,78   ,74   ,72   ];

% Window-Limit, 130 kW:
shotSeries{5} = [29101,29100,29106,29145,29136,29132];
rfPower{5}    = [130  ,131  ,132  ,132  ,133  ,133  ];

% Window-Limit, 115 kW:
shotSeries{6} = [29102,29099,29105,29144,29137,29131];
rfPower{6}    = [117  ,116  ,116  ,119  ,118  ,119  ];

% Window-Limit, 100 kW:
shotSeries{7} = [29103,29098,29108,29143,29138,29130];
rfPower{7}    = [103  ,104  ,103  ,103  ,103  ,103  ];

% Window-Limit, 80 kW:
shotSeries{8} = [29104,29097,29107,29141,29139,29129];
rfPower{8}    = [88   ,86   ,88   ,88   ,88   ,87  ];

if 0
% Enable if one needs to search for shots:
% Select shots corresponding to 6 views at the same RF power:
% =========================================================================
    viewType = 'bottom';
    limitMode = 'window';
    viewSide = 'non-pit';
    rfPwr = [70,90];
    n = find(strcmpi(T.limitMode,limitMode) & strcmpi(T.viewType,viewType) & strcmpi(T.viewSide,viewSide) & T.X>rfPwr(1) & T.X<rfPwr(2))
    disp( {'ii: ', num2str(T.Group(n(end))), 'jj: ', num2str((n(end))) } )
    T.rfPwrNet(n(end))
    T.shot(n(end))
end
%% 3- Gather data into structures:
% =========================================================================
cmptType = 2;
switch cmptType
    case 1
        seriesToAnalyze = 1:numel(shotSeries);
    case 2
        seriesToAnalyze = 4;
end

for kk = seriesToAnalyze
    rootAddress = 'C:\Users\nfc\Documents\ProtoMPEX_Ops\2020\';
    for jj = 1:numel(shotSeries{kk})
        % shot index:
        n = find(T.shot == shotSeries{kk}(jj));

        % Define address and file name:
        relPath = [T.date{n}]; relPath = [relPath(1:7),'\',relPath];
        fileName = ['\dataset_',num2str(T.dataset(n)),'_IRdata.mat'];
        filePath = [rootAddress,relPath,'\HeliconWindowIR',fileName];

        % load structure of .mat file:
        d = matfile(filePath);

        % Indentify shot needed within .mat file:
        m = find(d.shot == T.shot(n));

        % Extract required data from .mat file:
        f{kk}{jj}.shot = d.shot(1,m);
        f{kk}{jj}.limitMode = d.limitMode(1,m);
        f{kk}{jj}.viewType  = d.viewType(1,m);
        f{kk}{jj}.viewSide  = d.viewSide(1,m);
        f{kk}{jj}.temperature   = cell2mat(d.temperature(1,m));
        f{kk}{jj}.t_temperature = cell2mat(d.t_temperature(1,m));
        f{kk}{jj}.rfPwr = d.X(1,m);
        f{kk}{jj}.t0Plasma = d.t0Plasma(1,m);
        f{kk}{jj}.intf = cell2mat(d.intf(1,m));
        f{kk}{jj}.rngPlasma = cell2mat(d.rngPlasma(1,m));
    end
end

%% 4- Calculate temperature difference and mirror image:
% =========================================================================
% Dimensions of microbolometer chip
pixelsize = 17e-6; %[m]
% "x" is the vertical axis:
xchip = 240*pixelsize;
% "y" is the horizontal axis
ychip = 640*pixelsize;

for kk = seriesToAnalyze;
    for jj = 1:numel(shotSeries{kk})
        % Frame prior to plasma:
        n0 = 1;

        % Calculate dT and mirror image:
        for rr = 1:size(f{kk}{jj}.temperature,3)
            f{kk}{jj}.dT(:,:,rr)   = f{kk}{jj}.temperature(:,end:-1:1,rr) - f{kk}{jj}.temperature(:,end:-1:1,n0);
        end
        
        % Calculate the time base:
        f{kk}{jj}.t_dT = f{kk}{jj}.t_temperature - f{kk}{jj}.t_temperature(1);

        % Clear memory:
        f{kk}{jj}.temperature = [];
    end
end

%% 5- Get hottest frames to prepare for stitching:
% =========================================================================
for kk = seriesToAnalyze;
    for jj = 1:6;
        [~,fmax] = max(max(max(f{kk}{jj}.dT)));
        im{kk}{jj} = f{kk}{jj}.dT(:,:,fmax-1);
    end
end

% Assign coordinates:
[row,col,~] = size(f{kk}{jj}.dT);
ar1 = col/row;
x_im = linspace(-xchip/2,+xchip/2,row);
y_im = linspace(-ychip/2,+ychip/2,col);

%% 6- Stitch images together:
% =========================================================================
close all

% Description of shot series:
% kk = 1; % MPEX-Limit  , 140 kW:
% kk = 2; % MPEX-Limit  , 120 kW:
% kk = 3; % MPEX-Limit  , 100 kW:
% kk = 4; % MPEX-Limit  , 80  kW:
% kk = 5; % Window-Limit, 130 kW:
% kk = 6; % Window-Limit, 115 kW:
% kk = 7; % Window-Limit, 100 kW:
% kk = 8; % Window-Limit, 80  kW:

kk = seriesToAnalyze;
limitMode = cell2mat(f{kk}{1}.limitMode);
if 0
    viewSide = 'Pit';
else
    viewSide = 'Non-Pit';
end

% Plot raw data:
figure
ii = 1;
for jj = 1:6;
    if strcmpi(viewSide,f{kk}{jj}.viewSide)
        subplot(3,1,ii)
        surf(im{kk}{jj},'LineStyle','none')
        box on
        view([0,90])
        set(gca,'PlotBoxAspectRatio',[ar1 1 1])
        xlim([0,size(im{kk}{jj},2)])
        ylim([0,size(im{kk}{jj},1)])
        if ii == 1
            title([limitMode,' , ',viewSide, ' , RF: ', num2str(mean(rfPower{kk}),4),' kW'],'FontSize',11)
        end
        caxis([0,15])
        ii = ii + 1;
    end
end

% Shift and crop raw data:
switch kk
    case {1,2,3,4}
       switch viewSide
           case 'Pit'
               crop_rng{kk}{1}   = [077:240];
               crop_shift{kk}(1) = 4;
               crop_rng{kk}{2}   = [035:172];
               crop_shift{kk}(2) = 0;
               crop_rng{kk}{3}   = [001:162];
               crop_shift{kk}(3) = 5;
               compositeImage = [circshift(im{kk}{3}(crop_rng{kk}{3},:),crop_shift{kk}(3),2);...
                                 circshift(im{kk}{2}(crop_rng{kk}{2},:),crop_shift{kk}(2),2);...
                                 circshift(im{kk}{1}(crop_rng{kk}{1},:),crop_shift{kk}(1),2)];
           case 'Non-Pit'
               crop_rng{kk}{4}   = [063:240];
               crop_shift{kk}(4) = -8;
               crop_rng{kk}{5}   = [065:168];
               crop_shift{kk}(5) = 0;
               crop_rng{kk}{6}   = [001:163];
               crop_shift{kk}(6) = 0;
               compositeImage = [circshift(im{kk}{6}(crop_rng{kk}{6},:),crop_shift{kk}(6),2);...
                                 circshift(im{kk}{5}(crop_rng{kk}{5},:),crop_shift{kk}(5),2);...
                                 circshift(im{kk}{4}(crop_rng{kk}{4},:),crop_shift{kk}(4),2)];
       end
    case {5,6,7,8}
        switch viewSide
           case 'Pit'
               crop_rng{kk}{1}   = [077:240];
               crop_shift{kk}(1) = 4;
               crop_rng{kk}{2}   = [037:172];
               crop_shift{kk}(2) = 0;
               crop_rng{kk}{3}   = [001:181];
               crop_shift{kk}(3) = 8;
               compositeImage = [circshift(im{kk}{3}(crop_rng{kk}{3},:),crop_shift{kk}(3),2);...
                                 circshift(im{kk}{2}(crop_rng{kk}{2},:),crop_shift{kk}(2),2);...
                                 circshift(im{kk}{1}(crop_rng{kk}{1},:),crop_shift{kk}(1),2)];
           case 'Non-Pit'
               crop_rng{kk}{4}   = [095:240];
               crop_shift{kk}(4) = -10;
               crop_rng{kk}{5}   = [061:195];
               crop_shift{kk}(5) = 0;
               crop_rng{kk}{6}   = [001:159];
               crop_shift{kk}(6) = 0;
               compositeImage = [circshift(im{kk}{6}(crop_rng{kk}{6},:),crop_shift{kk}(6),2);...
                                 circshift(im{kk}{5}(crop_rng{kk}{5},:),crop_shift{kk}(5),2);...
                                 circshift(im{kk}{4}(crop_rng{kk}{4},:),crop_shift{kk}(4),2)];
        end
end

[row_cI,col_cI] = size(compositeImage); 
 
figure
surf(compositeImage,'LineStyle','none')
box on
view([0,90])
set(gca,'PlotBoxAspectRatio',[col_cI/row_cI 1 1])
xlim([0,col_cI])
ylim([0,row_cI])
title([limitMode,' , ',viewSide, ' , RF: ', num2str(mean(rfPower{kk}),4),' kW'],'FontSize',11)
caxis([0,15])
colormap('hot')

%% 6 Create grid on image plane:
% =========================================================================
% close all

% Create x and y coordinates for composite image:
% =========================================================================
[row,col,~] = size(compositeImage);
ar2 = col/row;
% "x" is the vertical axis:
xchip = row*pixelsize;
x_cI = linspace(-xchip/2,+xchip/2,row);
% "y" is the horizontal axis
ychip = col*pixelsize;
y_cI = linspace(-ychip/2,+ychip/2,col);

% Object's dimensions in physical space:
% =========================================================================
% Radius of window:
Ra = 0.90*0.5*12.4*1e-2; % [m]
% Length of window:
bb = 0.64;
La = bb*12.11*2.54*1e-2; % [m]

% Ra = 0.5*12.4*1e-2; % [m]
% % Lenght of window:
% La = 12.11*2.54*1e-2; % [m]

% Define a 2D grid on the object's inner surface and wrt object's datum:
% =========================================================================
% Azimuthal angle:
t_1D  = linspace(+180+25,+360-25,180)*pi/180;
% Along axis of window:
z_1D = linspace(+0,+La,100)';
% 2D grid 
[t_2D,z_2D] = meshgrid(t_1D,z_1D);

% Create 3D grid of window's inner surface wrt to objects datum:
% =========================================================================
qx = Ra*cos(t_2D);
qy = Ra*sin(t_2D);
qz = z_2D;

% Initial guess:
% =========================================================================
% Axial location of camera relative to object's datum:
L(3) = -34.5*1e-2; % [m]
% Rotation about yy to rotate the objects's frame to the camera's frame:
a = -0*0.2*pi/180;
% Rotation about xx to rotate the objects's frame to the camera's frame:
b = +13*pi/180;

% Pin hole camera factors:
% =========================================================================
% Focal length:
C.f = 24.6*1e-3;
% Reflection factor
C.Omega = -1;

% Create rotation matrix:
% =========================================================================
% Rotation matrix defined as: [e] = Ryx*[s]
% where [e] is camera's referene frame
% [s] is object's reference frame
R1 = [+cos(a)       ; +0     ; -sin(a)       ];
R2 = [+sin(a)*sin(b); +cos(b); +cos(a)*sin(b)];
R3 = [+sin(a)*cos(b); -sin(b); +cos(a)*cos(b)];
R  = [R1,R2,R3];

% Reference points:
% =========================================================================
refType = 1;
switch refType
    case 1
        % Ref point on cylinder:
        q_star = [+0         ,-Ra        ,0  ]';
        % Ref point on image plane:
        r_star = [-0.102*1e-3,-1.92*1e-3 ,C.f ]';
    case 2
        % Ref point on cylinder:
        q_star = [+0         ,-Ra         ,La  ]';
        % Ref point on image plane:
        r_star = [-0.008*1e-3,+0.825*1e-3 ,C.f ]';
end


% Calibration matrix:
% =========================================================================
M = r_star*R3' + C.f*C.Omega*R';
G = M(1:2,1:2);
H = -M(1:2,3);
P = inv(G)*H;

% Define offset vector:
% =========================================================================
L(1) = q_star(1) + P(1)*(L(3) - q_star(3));
L(2) = q_star(2) + P(2)*(L(3) - q_star(3));
L = L'

 % Project object's grid to image plane:
 % =========================================================================
t0 = tic;
for jj = 1:numel(t_1D)
     for kk = 1:numel(z_1D)
         qq = [qx(kk,jj),qy(kk,jj),qz(kk,jj)]';
         rr = PinHoleCamera(C,R,L,qq);
         xi(kk,jj) = rr(1);
         yi(kk,jj) = rr(2);
         zi(kk,jj) = rr(3);
     end
end
t0 = toc(t0);
disp(['Transforming 3D grid into image plane grid: ',num2str(t0),' [s]'])

%%
% Plot Composite image:
% =========================================================================
figure('color','w')
hold on
surf(y_cI*1e3,x_cI*1e3,compositeImage,'LineStyle','none')
box on
view([0,90])
xlabel('y [mm]')
ylabel('x [mm]')
set(gca,'PlotBoxAspectRatio',[ar2 1 1])
xlim([-1,1]*max(y_cI*1e3))
ylim([-1,1]*max(x_cI*1e3))
colormap('hot')
% Draw coordinate system:
line([0,0.1],[0,0  ],[1,1]*1e2,'color','k','LineWidth',2)
text(0.1,0,300,'$e_2$','FontSize',12,'Interpreter','latex')
line([0,0  ],[0,0.1],[1,1]*1e2,'color','k','LineWidth',2)
text(0,0.1,300,'$e_1$','FontSize',12,'Interpreter','latex')

% Plot 3D grid to image plane:
% =========================================================================
nn = 5;
xx = xi(1:nn:end,1:nn:end);
yy = yi(1:nn:end,1:nn:end);
zz = zi(1:nn:end,1:nn:end);
plot3(yy*1e3  ,xx*1e3 ,zz*500 ,'w-')
plot3(yy'*1e3 ,xx'*1e3,zz'*500,'w-')
caxis([0,50])
%
% Use cylindrical 3D grid to interpolate image:
% =========================================================================
t0 = tic;
[Y,X] = meshgrid(y_cI,x_cI);
v = interp2(Y,X,compositeImage,yi,xi);
disp('Using interp2 function ...')
t0 = toc(t0);
disp(['Interpolating data: ',num2str(t0),' [s]'])

% Plot data:
% =========================================================================
figure('color','w')
hold on
surf(t_2D*180/pi,(qz/bb)*1e2,v,'LineStyle','none')
set(gca,'XTick',[0:45:360],'XDir','reverse')
set(gca,'YTick',[0:5:(La/bb)*1e2])
xlabel('Angle [deg]','Interpreter','Latex','FontSize',14)
ylabel('z [cm]','Interpreter','Latex','FontSize',14)
view([0,90])
axis tight
colormap('bone')
colormap('hot')
xlim([0,360])
caxis([0,20])
colorbar

% Draw antenna:
% Transverse straps:
line([360,000],[26,26],[50,50],'color','k','LineWidth',10)
line([360,000],[04,04],[50,50],'color','k','LineWidth',10)
% Bottom helical strap:
line([225,135],[04,26],[50,50],'color','k','LineWidth',10)
% HV top side helical strap:
line([045,000],[04,15],[50,50],'color','k','LineWidth',10)
% GND top side helical strap:
line([360,315],[15,26],[50,50],'color','k','LineWidth',10)

% World frame:
if 0
    figure; 
    hold on
    plot3(qy,qx,-qz,'k.-')
    plot3(Ly,Lx,-Lz,'ro')
    axis image
    xlim([-3,3]*Ra)
    ylim([-3,3]*Ra)
    grid on
    xlabel('x')
    ylabel('y')
    zlabel('z')
end 
