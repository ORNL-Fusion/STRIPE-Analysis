% Objective:
% Load shot series .mat files and apply grids to all views: 

%% Process: 
% 1- Load shot series.mat files:
% 6- Save composite data into .mat files:
clc
clear all
close all

saveData = 1;
saveFig  = 1;

for kk = 1:8
kk

%% 1- Load shot series.mat files:
% =========================================================================
% kk = 8;
% kk = 1;
fileName = ['ShotSeriesData_',num2str(kk)];
load(fileName)

%% 2- Calculate grids and project object to image plane:
% =========================================================================
close all

% Object's dimensions in physical space:
% =========================================================================
% Radius of window:
aa = 1;
Ra = aa*0.5*12.4*1e-2; % [m]
% Length of window:
bb = 0.64;
bb = 0.70;
% bb = 1;
La = bb*12.11*2.54*1e-2; % [m]

% Pin hole camera factors:
% =========================================================================
% Focal length:
C.f = 24.6*1e-3;
% Reflection factor
C.Omega = -1;

viewsToAnalyze = [1:6];
% viewsToAnalyze = [4:6];
% viewsToAnalyze = [1:3];
% viewsToAnalyze = 6;

% Define independent variables on pin-hole camera model:
% =========================================================================
% L3: Axial location of camera relative to object's datum
% a: Rotation about yy to rotate the objects's frame to the camera's frame
% b: Rotation about xx to rotate the objects's frame to the camera's frame
switch kk
    case {1,2,3,4}
        % Define independent variables:
        a  = [+3   ,0    ,-1.85  ,+3   ,0    ,-3   ]*pi/180; % [Rad]
        b  = [-14  ,-14  ,-14    ,12   ,12   ,12 ]*pi/180; % [Rad]
        L3 = [-34.5,-34.5,-34.5  ,-34.5,-34.5,-34.5]*1e-2; % [m]

        a  = [+4.25,0    ,-1.85  ,+3   ,0    ,-3   ]*pi/180; % [Rad]
        b  = [-11.51  ,-14  ,-14    ,12   ,12   ,12 ]*pi/180; % [Rad]
        L3 = [-34.5,-34.5,-34.5  ,-34.5,-34.5,-34.5]*1e-2; % [m]
        
        % Define angular range of object's grid:
        t1 = [0    ,00   ,00   ,180  ,180  ,180  ]*pi/180;
        t2 = [180  ,180  ,180  ,360  ,360  ,360  ]*pi/180;         
    case {5,6,7,8}
        % Define independent variables:
        a  = [+3   ,0    ,-3     ,+3   ,0    ,-2.65]*pi/180; % [Rad]
        b  = [-14  ,-14  ,-14    ,12   ,12   ,12   ]*pi/180; % [Rad]
        L3 = [-34.5,-34.5,-34.5  ,-34.5,-34.5,-34.5]*1e-2; % [m]

        % Define angular range of object's grid:
        t1 = [0    ,00   ,00   ,360  ,225  ,180  ]*pi/180;
        t2 = [180  ,180  ,180  ,180  ,315  ,315  ]*pi/180;
end

% Define reference points between image plane and object' inner surface:
% =========================================================================
switch kk
    case {1,2,3,4}
        %           [vertical   ,horizontal ,z   ]
        q_star{1} = [+0         ,+Ra         ,La  ]';
        r_star{1} = [-1.938*1e-3,-1.65*1e-3 ,C.f ]';
        
        q_star{2} = [+0         ,+Ra        ,La  ]';
        r_star{2} = [-0.300*1e-3,-1.60*1e-3 ,C.f ]';
         
        q_star{3} = [+0         ,+Ra        ,La  ]';
        r_star{3} = [+1.69*1e-3,-1.73*1e-3,C.f ]';
        
        q_star{4} = [+0         ,-Ra        ,La   ]';
        r_star{4} = [-1.700*1e-3,+0.945*1e-3 ,C.f ]';
        
        q_star{5} = [+0         ,-Ra         ,La  ]';
        r_star{5} = [-0.025*1e-3,+0.785*1e-3 ,C.f ]';
         
        q_star{6} = [+0         ,-Ra         ,La  ]';
        r_star{6} = [+1.65*1e-3,+0.785*1e-3 ,C.f ]';
    case {5,6,7,8}
        q_star{1} = [+0         ,+Ra         ,La  ]';
        r_star{1} = [-1.938*1e-3,-1.65*1e-3 ,C.f ]';
        
        q_star{2} = [+0         ,+Ra        ,La  ]';
        r_star{2} = [-0.300*1e-3,-1.60*1e-3 ,C.f ]';
         
        q_star{3} = [+0         ,+Ra        ,La  ]';
        r_star{3} = [+2.11*1e-3,-1.80*1e-3,C.f ]';
        
        q_star{4} = [+0         ,-Ra        ,La   ]';
        r_star{4} = [-1.700*1e-3,+0.945*1e-3 ,C.f ]';
        
        q_star{5} = [+0         ,-Ra         ,La  ]';
        r_star{5} = [-0.025*1e-3,+0.785*1e-3 ,C.f ]';
         
        q_star{6} = [+0         ,-Ra         ,La  ]';
        r_star{6} = [+1.68*1e-3,+0.81*1e-3 ,C.f ]';
end

for jj = viewsToAnalyze
    % Get hottest frame:
    % =====================================================================
    fr = 55;
    im{jj} = u{jj}.dT(:,:,fr);
    
    % Creat chip coordinates:
    % =====================================================================
    xI{jj} = u{jj}.xI;
    yI{jj} = u{jj}.yI;
    
    % Define 2D grid on the object's inner surface and wrt object's datum:
    % =====================================================================
    % Azimuthal angle:
    t_1D{jj} = linspace(t1(jj),t2(jj),50);
    % Along axis of window:
    z_1D{jj} = linspace(+0,+La,50)';
    % 2D grid 
    [t_2D{jj},z_2D{jj}] = meshgrid(t_1D{jj},z_1D{jj});

    % Create 3D grid of window's inner surface wrt to objects datum:
    % =====================================================================
    qx{jj} = Ra*cos(t_2D{jj});
    qy{jj} = Ra*sin(t_2D{jj});
    qz{jj} = z_2D{jj};
    
    % Create rotation matrix:
    % =========================================================================
    % Rotation matrix defined as: [e] = Ryx*[s]
    % where [e] is camera's referene frame
    % [s] is object's reference frame
    R1     = [+cos(a(jj))           ; +0         ; -sin(a(jj))           ];
    R2     = [+sin(a(jj))*sin(b(jj)); +cos(b(jj)); +cos(a(jj))*sin(b(jj))];
    R3     = [+sin(a(jj))*cos(b(jj)); -sin(b(jj)); +cos(a(jj))*cos(b(jj))];
    R{jj}  = [R1,R2,R3];

    % Calibration matrix:
    % =========================================================================
    M = r_star{jj}*R3' + C.f*C.Omega*R{jj}';
    G = M(1:2,1:2);
    H = -M(1:2,3);
    P = inv(G)*H;

    % Define offset vector:
    % =========================================================================
    L{jj}(1) = q_star{jj}(1) + P(1)*(L3(jj) - q_star{jj}(3));
    L{jj}(2) = q_star{jj}(2) + P(2)*(L3(jj) - q_star{jj}(3));
    L{jj}(3) = L3(jj);
    if size(L{jj},2)>1
        L{jj}    = L{jj}';
    end
    
    % Create 2D grid on image plane:
    % =====================================================================
    t0 = tic;
        for cc = 1:numel(t_1D{jj})
             for rr = 1:numel(z_1D{jj})
                 qq = [qx{jj}(rr,cc),qy{jj}(rr,cc),qz{jj}(rr,cc)]';
                 ri = PinHoleCamera(C,R{jj},L{jj},qq);
                 xi{jj}(rr,cc) = ri(1);
                 yi{jj}(rr,cc) = ri(2);
                 zi{jj}(rr,cc) = ri(3);
             end
        end
    t0 = toc(t0);
    disp(['Transforming 3D grid into image plane grid: ',num2str(t0),' [s]'])
end

    
% Plot image plane with grid:
% =========================================================================

figure('color','w')
ii = 1;
pos = [2,4,6,1,3,5];
for jj = viewsToAnalyze
    if numel(viewsToAnalyze) == 1;
    else
            subplot(3,2,pos(jj))
    end
    hold on
    surf(yI{jj}*1e3,xI{jj}*1e3,im{jj},'LineStyle','none')
    title(['jj: ',num2str(jj)])
    [rows,cols,~] = size(im{jj});
    box on
    view([0,90])
    xlabel('y [mm]')
    ylabel('x [mm]')
    set(gca,'PlotBoxAspectRatio',[cols/rows 1 1])
    xlim([-1,1]*max(yI{jj}*1e3))
    ylim([-1,1]*max(xI{jj}*1e3))
    colormap('hot')
    % Draw coordinate system:
    line([0,0.1],[0,0  ],[1,1]*1e2,'color','k','LineWidth',2)
    text(0.1,0,300,'$e_2$','FontSize',12,'Interpreter','latex')
    line([0,0  ],[0,0.1],[1,1]*1e2,'color','k','LineWidth',2)
    text(0,0.1,300,'$e_1$','FontSize',12,'Interpreter','latex')
    caxis([0,25]) 

    % Plot 3D grid to image plane:
    % =========================================================================
    nn = 2;
    xx = xi{jj}(1:nn:end,1:nn:end);
    yy = yi{jj}(1:nn:end,1:nn:end);
    zz = zi{jj}(1:nn:end,1:nn:end);
    plot3(yy*1e3  ,xx*1e3 ,zz*500 ,'w-')
    plot3(yy'*1e3 ,xx'*1e3,zz'*500,'w-')
    ii = ii + 1;
end

%% 3- Interpolate image plane with 2D grid:
% =========================================================================
for jj = viewsToAnalyze;
    t0 = tic;
    [Y,X] = meshgrid(yI{jj},xI{jj});
    for fr = 1:size(u{jj}.dT,3)
        v{jj}(:,:,fr) = interp2(Y,X,u{jj}.dT(:,:,fr),yi{jj},xi{jj},'*linear');
    end
    disp('Using interp2 function ...')
    t0 = toc(t0);
    disp(['Interpolating data: ',num2str(t0),' [s]'])
end

% Plot data:
% =========================================================================
fr = 34;
for jj = viewsToAnalyze;
figure('color','w')
for fr = 5:5:55
    surf(t_2D{jj}*180/pi,(z_2D{jj}/bb)*1e2,v{jj}(:,:,fr),'LineStyle','none')
    set(gca,'XTick',[0:45:360],'XDir','reverse')
    set(gca,'YTick',[0:5:(La/bb)*1e2])
    xlabel('Angle [deg]','Interpreter','Latex','FontSize',14)
    ylabel('z [cm]','Interpreter','Latex','FontSize',14)
    view([0,90])
    axis tight
    colormap('bone')
    colormap('hot')
    xlim([0,360])
    caxis([0,25])
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
end

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
        caxis([0,25])
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
        caxis([0,25])
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
end
disp('End of Script!!')