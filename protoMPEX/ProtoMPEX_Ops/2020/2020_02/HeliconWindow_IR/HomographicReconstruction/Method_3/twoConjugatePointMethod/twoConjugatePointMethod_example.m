% Objective:
% Constrain the image plane grid using the two-conjugate points method

% Context:
% We have developed a method that relies on two conjugate points in order
% to fully constrain the location of the object's grid

%% Process: 
% 1- Load shot series.mat files:

clc
clear all
close all

% Select shotSeries and view:
% =========================================================================
kk = 5;
jj = 3;

% Load data:
% =========================================================================
homeAddress = cd;
cd ..
cd ..
cd Method_2
fileName = ['ShotSeriesData_',num2str(kk)];
dat = load(fileName);
cd(homeAddress)

%%
close all

% Object's dimensions in physical space:
% =========================================================================
% Radius of window:
aa = 1;
Ra = aa*0.5*12.6*1e-2; % [m]
% Length of window:
bb = 1.;
La = bb*12.11*2.54*1e-2 - 2.0*1e-2; % [m]

% Pin hole camera factors:
% =========================================================================
% Focal length:
C.f = 24.6*1e-3;
% Reflection factor
C.Omega = -1;
% focal length and reflection factor lumped together:
m = C.f*C.Omega;
% "z" offset of camera:
L3 = -0.39;
    
% Creat chip coordinates:
% =========================================================================
xI = dat.u{jj}.xI;
yI = dat.u{jj}.yI;

% Create image:
% =========================================================================
fr = 55;
im = dat.u{jj}.dT(:,:,fr);

%%
close all

% Plot data:
% =========================================================================
figure('color','w')
surf(yI,xI,im,'LineStyle','none')
view([0,90])
axis image

% Define conjugate points:
% =========================================================================
% Recall that "x" axis is vertical and "y" axis is horizontal
% Use [vertical,horizontal,z]
% Point A:
r_star.A = [-1.75e-3,+1.1e-3,C.f]';
q_star.A = [0      ,Ra     ,0  ]';

% r_star.A = [-1.699e-3,+1.03e-3,C.f]';
% q_star.A = [0      ,Ra     ,0  ]';

% top
r_star.A = [+1.95e-3,-2.02e-3,C.f]';
q_star.A = [Ra      ,0       ,0  ]';

r_star.A = [-1.565e-3,+1.098e-3,C.f]';
q_star.A = [0        ,+Ra      ,0  ]';

% Point B:
r_star.B = [-1.938e-3,-1.65e-3,C.f]';
q_star.B = [0        ,Ra      ,La ]';

r_star.B = [-1.87e-3,-1.675e-3,C.f]';
q_star.B = [0        ,Ra      ,La ]';

% r_star.B = [-1.972e-3,-1.677e-3,C.f]';
% q_star.B = [0        ,Ra      ,La ]';


r_star.B = [-1.87e-3,-1.645e-3,C.f]';
q_star.B = [0        ,Ra      ,La ]';


% pA = ginput(1);
% r_star.A = [pA(2),pA(1),C.f]';
% q_star.A = [0    ,+Ra  ,0  ]';
% r_star.A = [pA(2),pA(1),C.f]';
% q_star.A = [-Ra  ,0    ,0  ]';

% pB = ginput(1);
% r_star.B = [pB(2),pB(1),C.f]';
% q_star.B = [0    ,+Ra  ,La ]';
% 

homeAddress = cd;
dataAddress  = 'C:\Users\nfc\Documents\ProtoMPEX_Ops\2020\2020_02\HeliconWindow_IR\HomographicReconstruction\Method_3\shotSeriesMatFiles';
cd(dataAddress)
fileName = ['shotSeries_',num2str(kk),'_Params_Pit_Bottom.mat'];
load(fileName)
cd(homeAddress)

r_star.A = [params.rxA*1e-3      ,params.ryA*1e-3      ,C.f             ]';
q_star.A = [eval(params.qxA),eval(params.qyA),eval(params.qzA)]';

r_star.B = [params.rxB*1e-3      ,params.ryB*1e-3      ,C.f             ]';
q_star.B = [eval(params.qxB),eval(params.qyB),eval(params.qzB)]';

% Reduced conjugate points:
% =========================================================================
% Point A:
u.A  = q_star.A(1:2);
v.A  = r_star.A(1:2);
Q3.A = q_star.A(3) - L3;

% Point B:
u.B  = q_star.B(1:2);
v.B  = r_star.B(1:2);
Q3.B = q_star.B(3) - L3;

% Create rotation matrix:
% =========================================================================
% Rotation matrix defined as: [e] = Ryx*[s]
% where [e] is camera's referene frame
% [s] is object's reference frame

% a: Rotation about yy to rotate the objects's frame to the camera's frame
% b: Rotation about xx to rotate the objects's frame to the camera's frame

R{1,1} = @(a,b) +cos(a);
R{1,2} = @(a,b) +sin(a).*sin(b);
R{1,3} = @(a,b) +sin(a).*cos(b);
R{2,1} = @(a,b) 0;
R{2,2} = @(a,b) +cos(b);
R{2,3} = @(a,b) -sin(b);
R{3,1} = @(a,b) -sin(a);
R{3,2} = @(a,b) +cos(a).*sin(b);
R{3,3} = @(a,b) +cos(a).*cos(b);

R1 = @(a,b) [R{1,1}(a,b); R{2,1}(a,b); R{3,1}(a,b)];
R2 = @(a,b) [R{1,2}(a,b); R{2,2}(a,b); R{3,2}(a,b)];
R3 = @(a,b) [R{1,3}(a,b); R{2,3}(a,b); R{3,3}(a,b)];
RR = @(a,b) [R1(a,b), R2(a,b), R3(a,b)];

% M matrix:
% =========================================================================
fld = ['A','B'];
for pp  = 1:2 
    for rr = 1:2
        for cc = 1:3
            k = fld(pp);
            r_i = r_star.(k)(rr);
            % M_ij       =        r_i*R_j3         + m*R_ji
            M{rr,cc}.(k) = @(a,b) r_i*R{cc,3}(a,b) + m*R{cc,rr}(a,b);
        end
    end
end

% L1 and L2 linear equation:
% =========================================================================
for pp = 1:2
    % Select conjugate point:
    k = fld(pp);
    % Adjugate matrix of G:
    adjG.(k) = @(a,b) [+M{2,2}.(k)(a,b), -M{1,2}.(k)(a,b);...
                       -M{2,1}.(k)(a,b), +M{1,1}.(k)(a,b)];   
    % Determinant of G:
    detG.(k) = @(a,b)  M{1,1}.(k)(a,b)*M{2,2}.(k)(a,b) -  M{2,1}.(k)(a,b)*M{1,2}.(k)(a,b);                    
    % Factor K
    K.(k) = @(a,b) Q3.(k)./detG.(k)(a,b);
    % Equation for offset vector:
    t.(k) = @(a,b) u.(k) + K.(k)(a,b)*adjG.(k)(a,b)*[M{1,3}.(k)(a,b),M{2,3}.(k)(a,b)]';
end

% Assemble F:
% =========================================================================
F = @(a,b) (u.A-u.B) + ...
    K.A(a,b)*adjG.A(a,b)*[M{1,3}.A(a,b);M{2,3}.A(a,b)] - ...
    K.B(a,b)*adjG.B(a,b)*[M{1,3}.B(a,b);M{2,3}.B(a,b)];

% Put function in form needed for non-linear solver:
Fsolve = @(x) F(x(1),x(2));

% Solve non-linear equation:
% =========================================================================
tic
[x, ~, ~, exitflag, output, ~] = newtonraphson(Fsolve, [4,-13]*pi/180);
toc

% Assign solution to angles and offset vector:
% =========================================================================
% View angles:
a = x(1);
b = x(2);

% Offset vector:
L1= [1,0]*t.A(a,b);
L2 = [0,1]*t.A(a,b);
L = [L1,L2,L3]';

    % Define 2D grid on the object's inner surface and wrt object's datum:
    % =====================================================================
    % Azimuthal angle:
    t_1D{jj} = linspace(0,2*pi,360);
    % Along axis of window:
    z_1D{jj} = linspace(+0,+La,50)';
    % 2D grid 
    [t_2D{jj},z_2D{jj}] = meshgrid(t_1D{jj},z_1D{jj});

    % Create 3D grid of window's inner surface wrt to objects datum:
    % =====================================================================
    qx{jj} = Ra*cos(t_2D{jj});
    qy{jj} = Ra*sin(t_2D{jj});
    qz{jj} = z_2D{jj};
    
    % Create 2D grid on image plane:
    % =====================================================================
    t0 = tic;
        for cc = 1:numel(t_1D{jj})
             for rr = 1:numel(z_1D{jj})
                 qq = [qx{jj}(rr,cc),qy{jj}(rr,cc),qz{jj}(rr,cc)]';
                 ri = PinHoleCamera(C,RR(a,b),L,qq);
                 xi(rr,cc) = ri(1);
                 yi(rr,cc) = ri(2);
                 zi(rr,cc) = ri(3);
             end
        end
    t0 = toc(t0);
    disp(['Transforming 3D grid into image plane grid: ',num2str(t0),' [s]'])
    
% Plot image plane with grid:
% =========================================================================

figure('color','w')
ii = 1;
pos = [2,4,6,1,3,5];
    if numel(1) == 1;
    else
            subplot(3,2,pos(jj))
    end
    hold on
    surf(yI*1e3,xI*1e3,im,'LineStyle','none')
    title(['jj: ',num2str(jj)])
    [rows,cols,~] = size(im);
    box on
    view([0,90])
    xlabel('y [mm]')
    ylabel('x [mm]')
    set(gca,'PlotBoxAspectRatio',[cols/rows 1 1])
    xlim([-1,1]*max(yI*1e3))
    ylim([-1,1]*max(xI*1e3))
    colormap('hot')
    % Draw coordinate system:
    line([0,0.1],[0,0  ],[1,1]*1e2,'color','k','LineWidth',2)
    text(0.1,0,300,'$e_2$','FontSize',12,'Interpreter','latex')
    line([0,0  ],[0,0.1],[1,1]*1e2,'color','k','LineWidth',2)
    text(0,0.1,300,'$e_1$','FontSize',12,'Interpreter','latex')
    caxis([0,20]) 

    % Plot 3D grid to image plane:
    % =========================================================================
    nn = 2;
    xx = xi(1:nn:end,1:nn:end);
    yy = yi(1:nn:end,1:nn:end);
    zz = zi(1:nn:end,1:nn:end);
    plot3(yy*1e3  ,xx*1e3 ,zz*500 ,'w-')
    plot3(yy'*1e3 ,xx'*1e3,zz'*500,'w-')
    ii = ii + 1;
    
    % Plot conjugate points:
    % =====================================================================
    plot3(r_star.A(2)*1e3,r_star.A(1)*1e3,100,'w.','MarkerSize',18)
    plot3(r_star.B(2)*1e3,r_star.B(1)*1e3,500,'w.','MarkerSize',18)


%% 3- Interpolate image plane with 2D grid:
% =========================================================================
t0 = tic;
[Y,X] = meshgrid(yI,xI);
for fr = 1:size(dat.u{jj}.dT,3)
    vv{jj}(:,:,fr) = interp2(Y,X,dat.u{jj}.dT(:,:,fr),yi,xi,'*linear');
end
disp('Using interp2 function ...')
t0 = toc(t0);
disp(['Interpolating data: ',num2str(t0),' [s]'])

% Plot data:
% =========================================================================
fr = 34;
figure('color','w')
for fr = 5:5:55
    surf(t_2D{jj}*180/pi,(z_2D{jj}/bb)*1e2,vv{jj}(:,:,fr),'LineStyle','none')
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
    title(['frame: ',num2str(fr)])


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

    drawnow
end
title(['jj: ',num2str(jj)])

return


%% 4- Stitch images into a single temperature distribution map:

% Find the number of frames of the shortest video in the shot series:
for jj = 1:6
    [rw,cl,fr(jj)] = size(v{jj});
end
fr = min(fr);

% Create time base for composite image:
t_dT = u{1}.t_dT(1:fr);

% Create coordinates of the composite image:
phi_1D = linspace(0,360,2*numel(t_1D{1}));
[phi_2D,s_2D] = meshgrid(phi_1D,z_1D{1});

% Allocate memory:
z1 = zeros(rw,cl);
z2 = zeros(rw,cl);
z  = zeros(rw,2*cl,fr);

% Create composite frames:
tic
for ff = 1:size(z,3)
    for rr = 1:size(z1,1)
        for cc = 1:size(z1,2)
            z1(rr,cc) = max([v{4}(rr,cc,ff),v{5}(rr,cc,ff),v{6}(rr,cc,ff)],[],'omitnan');
            z2(rr,cc) = max([v{1}(rr,cc,ff),v{2}(rr,cc,ff),v{3}(rr,cc,ff)],[],'omitnan');
        end
    end
    z(:,:,ff) = [z2,z1];
end
toc

% Plot and animate composite frames:
figure('color','w')
for fr = 5:5:55
        surf(phi_2D,(s_2D/bb)*1e2,z(:,:,fr),'LineStyle','none')
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
        title([cell2mat(u{jj}.limitMode),', RF: ',num2str(u{jj}.rfPwr),' [kW] , frame: ',num2str(fr)])

        % Draw antenna:
        w = 13;
        % Transverse straps:
        line([360,000],[26,26],[50,50],'color','k','LineWidth',w)
        line([360,000],[04,04],[50,50],'color','k','LineWidth',w)
        % Bottom helical strap:
        line([225,135],[04,26],[50,50],'color','k','LineWidth',w)
        % HV top side helical strap:
        line([045,000],[04,15],[50,50],'color','k','LineWidth',w)
        % GND top side helical strap:
        line([360,315],[15,26],[50,50],'color','k','LineWidth',w)
        drawnow
end

%% 5- Remove areas affected by vignetting:

if 0
    angleRemove = 1+[0:25,[180-25:180+25],[359-25:359]];
    cl_rng{1} = find(phi_1D <= 30);
    cl_rng{2} = find(phi_1D >= (180-20) & phi_1D <= (180+35) );
    cl_rng{3} = find(phi_1D >= (360-30) & phi_1D <= 360 );

    for ff = 1:size(z,3)
        for ss = 1:numel(cl_rng)
               z(:,cl_rng{ss},ff) = [30];
        end
    end
end

figure('color','w')
for fr = 5:5:55
        surf(phi_2D,(s_2D/bb)*1e2,z(:,:,fr),'LineStyle','none')
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
        title([cell2mat(u{jj}.limitMode),', RF: ',num2str(u{jj}.rfPwr),' [kW] , frame: ',num2str(fr)])

        % Draw antenna:
        w = 13;
        % Transverse straps:
        line([360,000],[26,26],[50,50],'color','k','LineWidth',w)
        line([360,000],[04,04],[50,50],'color','k','LineWidth',w)
        % Bottom helical strap:
        line([225,135],[04,26],[50,50],'color','k','LineWidth',w)
        % HV top side helical strap:
        line([045,000],[04,15],[50,50],'color','k','LineWidth',w)
        % GND top side helical strap:
        line([360,315],[15,26],[50,50],'color','k','LineWidth',w)
        drawnow
end

if saveFig
    figureName = [cell2mat(u{jj}.limitMode),'_RF_',num2str(u{jj}.rfPwr),'kW'];
    saveas(gcf,figureName,'tiffn')
end

%% 6- Save composite data into .mat files:
% =========================================================================

% Organizing output data into a structure:
% =========================================================================
f.dT = z;
f.t_dT = t_dT;
f.phi_2D = phi_2D;
f.s_2D = s_2D;
f.limitMode = u{1}.limitMode;
for jj = 1:6
    f.shots(jj) = u{jj}.shot;
    f.rfPwr(jj) = u{jj}.rfPwr;
    f.thermalParam(jj) = u{jj}.thermalParam;
end

% Saving data:
% =========================================================================
if saveData
    t1 = tic;
    disp('Saving data ...')
    fileName = ['CompositeData_ShotSeries_',num2str(kk),'.mat'];
    save(fileName,'f')
    t1 = toc(t1);
    disp(['Data Saved!! took ',num2str(t1),' seconds'])
    beep
end

% notes:
% The next step is to apply the inverse method to the temperature
% calculation
% Calculate power integrated over the entire area
% extrapolate to 200 kW
% Reconstruct the surface deposition image of the old window with the same
% cooridnate system as akk the thermal data:
% then we need to fill in the gaps where we do not have data
