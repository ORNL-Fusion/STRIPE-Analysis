% Objective:
% Using the constrained view parameters created on step 3, we now need to
% stitch 2D maps from a given shot series into a single 2D map spanning 0
% to 360 degrees
%% Process: 
% ####################################
% Process needs to be described here
% ####################################

clc
clear all
close all

saveData = 0;
saveFig  = 0;

% Address where IR .mat files are stored:
% D = dir('Step_2_*');
% dirFlag = [D.isdir];
% n = find(dirFlag == 1)
% dataAddress_step_2 = D(n).name;

clc
clear
close all

saveData = 0;
saveFig  = 0;

% Folder where this script lives
scriptDir = fileparts(mfilename('fullpath'));
parentDir = fileparts(scriptDir);

% ---- Find Step_2_* ----
D2 = dir(fullfile(scriptDir, 'Step_2_*'));
D2 = D2([D2.isdir]);

if isempty(D2)
    D2 = dir(fullfile(parentDir, 'Step_2_*'));
    D2 = D2([D2.isdir]);
end

if isempty(D2)
    error('No Step_2_* directory found in script folder or parent folder.');
else
    dataAddress_step_2 = fullfile(D2(1).folder, D2(1).name);
end

% ---- Find Step_3_* ----
D3 = dir(fullfile(scriptDir, 'Step_3_*'));
D3 = D3([D3.isdir]);

if isempty(D3)
    D3 = dir(fullfile(parentDir, 'Step_3_*'));
    D3 = D3([D3.isdir]);
end

if isempty(D3)
    error('No Step_3_* directory found in script folder or parent folder.');
else
    dataAddress_step_3 = fullfile(D3(1).folder, D3(1).name);
end

disp(['Step 2 folder found: ', dataAddress_step_2])
disp(['Step 3 folder found: ', dataAddress_step_3])

% kk represents the shot series:
for kk = 1


%% 1- Load shot series.mat files:
% =========================================================================
homeAddress = cd;
cd(dataAddress_step_2)
fileName = ['ShotSeries_',num2str(kk),'_IRdata.mat'];
thermaldata = load(fileName);
cd(homeAddress)

% Assign variables:
dT    = thermaldata.dT;
t_dT  = thermaldata.t_dT;
xI    = thermaldata.xI;
yI    = thermaldata.yI;
rfPwr = thermaldata.X;
shot  = thermaldata.shot;

%% 2- Calculate grids and project object to image plane:
% =========================================================================
close all

% Object's dimensions in physical space:
% =========================================================================
% Radius of window:
Ra = 0.5*12.4*1e-2; % [m]
% Length of window:
La = 12.11*2.54*1e-2; % [m]

% Pin hole camera factors:
% =========================================================================
% Focal length:
focalLength = thermaldata.focalLength;
C.f = focalLength;
% Reflection factor
Omega = -1;
C.Omega = Omega;

viewsToAnalyze = [1:6];

% Define independent variables on pin-hole camera model:
% =========================================================================
% L3: Axial location of camera relative to object's datum
% a: Rotation about yy to rotate the objects's frame to the camera's frame
% b: Rotation about xx to rotate the objects's frame to the camera's frame
% switch kk
%     case {1,2,3,4}
%         % For shotsSeries 1 to 4, we need to use the same parameters
%         % We Just need to identify which viewType and viewSide:
%         
% %         gg = '1';
% %         % Pit side:
% %         % -----------------------------------------------------------------
% %         % Top view:
% %         viewSide = 'Pit';
% %         viewType = 'Top';
% %         fileName = ['shotSeries_',gg,'_Params_',viewSide,'_',viewType,'.mat'];  
% %         cd(dataAddress_step_3)
% %         load(fileName)
% %         cd(homeAddress)
%         
%         
%         
%         % Define independent variables:
%         a  = [+3   ,0    ,-1.85  ,+3   ,0    ,-3   ]*pi/180; % [Rad]
%         b  = [-14  ,-14  ,-14    ,12   ,12   ,12 ]*pi/180; % [Rad]
%         L3 = [-34.5,-34.5,-34.5  ,-34.5,-34.5,-34.5]*1e-2; % [m]
% 
%         a  = [+4.25,0    ,-1.85  ,+3   ,0    ,-3   ]*pi/180; % [Rad]
%         b  = [-11.51  ,-14  ,-14    ,12   ,12   ,12 ]*pi/180; % [Rad]
%         L3 = [-34.5,-34.5,-34.5  ,-34.5,-34.5,-34.5]*1e-2; % [m]
%         
%         % Define angular range of object's grid:
%         t1 = [0    ,00   ,00   ,180  ,180  ,180  ]*pi/180;
%         t2 = [180  ,180  ,180  ,360  ,360  ,360  ]*pi/180;         
%     case {5,6,7,8}
%         % Define independent variables:
%         a  = [+3   ,0    ,-3     ,+3   ,0    ,-2.65]*pi/180; % [Rad]
%         b  = [-14  ,-14  ,-14    ,12   ,12   ,12   ]*pi/180; % [Rad]
%         L3 = [-34.5,-34.5,-34.5  ,-34.5,-34.5,-34.5]*1e-2; % [m]
% 
%         % Define angular range of object's grid:
%         t1 = [0    ,00   ,00   ,360  ,225  ,180  ]*pi/180;
%         t2 = [180  ,180  ,180  ,180  ,315  ,315  ]*pi/180;
% end

% Define reference points between image plane and object' inner surface:
% =========================================================================
% switch kk
%     case {1,2,3,4}
%         %           [vertical   ,horizontal ,z   ]
%         q_star{1} = [+0         ,+Ra         ,La  ]';
%         r_star{1} = [-1.938*1e-3,-1.65*1e-3 ,C.f ]';
%         
%         q_star{2} = [+0         ,+Ra        ,La  ]';
%         r_star{2} = [-0.300*1e-3,-1.60*1e-3 ,C.f ]';
%          
%         q_star{3} = [+0         ,+Ra        ,La  ]';
%         r_star{3} = [+1.69*1e-3,-1.73*1e-3,C.f ]';
%         
%         q_star{4} = [+0         ,-Ra        ,La   ]';
%         r_star{4} = [-1.700*1e-3,+0.945*1e-3 ,C.f ]';
%         
%         q_star{5} = [+0         ,-Ra         ,La  ]';
%         r_star{5} = [-0.025*1e-3,+0.785*1e-3 ,C.f ]';
%          
%         q_star{6} = [+0         ,-Ra         ,La  ]';
%         r_star{6} = [+1.65*1e-3,+0.785*1e-3 ,C.f ]';
%     case {5,6,7,8}
%         q_star{1} = [+0         ,+Ra         ,La  ]';
%         r_star{1} = [-1.938*1e-3,-1.65*1e-3 ,C.f ]';
%         
%         q_star{2} = [+0         ,+Ra        ,La  ]';
%         r_star{2} = [-0.300*1e-3,-1.60*1e-3 ,C.f ]';
%          
%         q_star{3} = [+0         ,+Ra        ,La  ]';
%         r_star{3} = [+2.11*1e-3,-1.80*1e-3,C.f ]';
%         
%         q_star{4} = [+0         ,-Ra        ,La   ]';
%         r_star{4} = [-1.700*1e-3,+0.945*1e-3 ,C.f ]';
%         
%         q_star{5} = [+0         ,-Ra         ,La  ]';
%         r_star{5} = [-0.025*1e-3,+0.785*1e-3 ,C.f ]';
%          
%         q_star{6} = [+0         ,-Ra         ,La  ]';
%         r_star{6} = [+1.68*1e-3,+0.81*1e-3 ,C.f ]';
% end

for jj = viewsToAnalyze
    % Identify which limit mode:
    % =====================================================================
    limitMode{jj} = thermaldata.limitMode{jj};
    
    % Determine view type and side:
    % =====================================================================
    viewType{jj} = thermaldata.viewType{jj};
    viewSide{jj} = thermaldata.viewSide{jj};
    
    % Identify the file with the contrained view paramters:
    % =====================================================================    
    % Select which parameter file to use:
    switch kk
        case {1,2,3,4}
            gg = '1';
        case {5,6,7,8}
            gg = '5';
    end
    fileName = ['shotSeries_',gg,'_Params_',viewSide{jj},'_',viewType{jj},'.mat'];  
    cd(dataAddress_step_3)
    load(fileName)
    cd(homeAddress)
    
    % offset length vector:
    % =====================================================================
    L{jj} = [params.L1,params.L2,params.L3]';
    
    % Angular range:
    % =====================================================================
    t1(jj) = params.angle1*pi/180;
    t2(jj) = params.angle2*pi/180;
    
    % Assemble conjugate points:
    % =====================================================================
    % Point A:
    rA = [params.rxA*1e-3 ,params.ryA*1e-3 ]'  ; % [m]
    qA = [eval(params.qxA),eval(params.qyA),0]'; % [m]
    % Point B:
    rB = [params.rxB*1e-3  ,params.ryB*1e-3   ]'; % [m]
    qB = [eval(params.qxB),eval(params.qyB),La]'; % [m]
              
    % Get hottest frame:
    % =====================================================================
    fr = 55;
    im{jj} = dT{jj}(:,:,fr);
        
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
    
    % Constrain view:
    % =====================================================================
    m = focalLength*Omega;
    [a,b,L{jj},rotMatrix] = ConstrainView(m,rA,qA,rB,qB,L{jj});

    % Create 2D grid on image plane:
    % =====================================================================
    t0 = tic;
        [xi{jj},yi{jj},zi{jj}] = CreateImagePlaneGrid(qx{jj},qy{jj},qz{jj},rotMatrix,L{jj},focalLength,Omega);
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
    surf(yI*1e3,xI*1e3,im{jj},'LineStyle','none')
    title(['jj: ',num2str(jj)])
    [rows,cols,~] = size(im{jj});
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
    xx = xi{jj}(1:nn:end,1:nn:end);
    yy = yi{jj}(1:nn:end,1:nn:end);
    zz = zi{jj}(1:nn:end,1:nn:end);
    plot3(yy*1e3  ,xx*1e3 ,zz*500 ,'w-')
    plot3(yy'*1e3 ,xx'*1e3,zz'*500,'w-')
    ii = ii + 1;
end

%% 3- Interpolate image plane with 2D grid:
% =========================================================================

for jj = 1:6
    [rw,cl,fr(jj)] = size(dT{jj});
end
frNum = min(fr);

for jj = viewsToAnalyze;
    t0 = tic;
    [Y,X] = meshgrid(yI,xI);
    for fr = 1:frNum
        v{jj}(:,:,fr) = interp2(Y,X,dT{jj}(:,:,fr),yi{jj},xi{jj},'*linear');
    end
    disp('Using interp2 function ...')
    t0 = toc(t0);
    disp(['Interpolating data: ',num2str(t0),' [s]'])
end

% Plot data:
% =========================================================================
fr = 34;
bb = 1;
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
end

%% 4- Stitch images into a single temperature distribution map:

% Find the number of frames of the shortest video in the shot series:
for jj = 1:6
    [rw,cl,fr(jj)] = size(v{jj});
end
fr = min(fr);

% Create time base for composite image:
t_dT = t_dT{1}(1:fr);

% Create coordinates of the composite image:
phi_1D = linspace(0,360,2*numel(t_1D{1}));
[phi_2D,s_2D] = meshgrid(phi_1D,z_1D{1});

% Allocate memory:
z1 = zeros(rw,cl);
z2 = zeros(rw,cl);
z  = zeros(rw,2*cl,fr);

% Create composite frames:
tic
switch kk
    case {1,3,4,5,6,7,8}
    for ff = 1:size(z,3)
        for rr = 1:size(z1,1)
            for cc = 1:size(z1,2)
                z1(rr,cc) = max([v{4}(rr,cc,ff),v{5}(rr,cc,ff),v{6}(rr,cc,ff)],[],'omitnan');
                z2(rr,cc) = max([v{1}(rr,cc,ff),v{2}(rr,cc,ff),v{3}(rr,cc,ff)],[],'omitnan');
            end
        end
        z(:,:,ff) = [z2,z1];
    end
    case {2}
    for ff = 1:size(z,3)
        for rr = 1:size(z1,1)
            for cc = 1:size(z1,2)
                z1(rr,cc) = max([v{1}(rr,cc,ff),v{2}(rr,cc,ff),v{3}(rr,cc,ff)],[],'omitnan');
                z2(rr,cc) = max([v{4}(rr,cc,ff),v{5}(rr,cc,ff),v{6}(rr,cc,ff)],[],'omitnan');
            end
        end
        z(:,:,ff) = [z2,z1];
    end
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
        title([cell2mat(limitMode(jj)),', RF: ',num2str(round(mean(rfPwr))),' [kW] , frame: ',num2str(fr)])

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
        title([(limitMode{jj}),', RF: ',num2str(round(mean(rfPwr))),' [kW] , frame: ',num2str(fr)])

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
    figureName = [(limitMode{jj}),'_RF_',num2str(round(mean(rfPwr))),'kW'];
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
f.limitMode = limitMode{1};
for jj = 1:6
    f.shots(jj) = shot(jj);
    f.rfPwr(jj) = rfPwr(jj);
    f.thermalParam(jj) = thermaldata.thermalParam;
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

%% 4.3- Write stitched theta-z coordinates and data to CSV:
% =========================================================================
% Export directly from stitched composite plot coordinates:
%
% phi_2D -> theta [deg]
% s_2D   -> z [m]
% z      -> stitched data

csvFileName = ['StitchedMap_ShotSeries_',num2str(kk),'.csv'];

Frame      = [];
Time_s     = [];
Theta_deg  = [];
Z_m        = [];
Value      = [];

for ff = 1:size(z,3)

    % Current stitched frame
    valueFrame = z(:,:,ff);

    % Flatten stitched coordinates
    thetaVec = phi_2D(:);
    zVec     = s_2D(:);
    valVec   = valueFrame(:);

    % Remove invalid points
    validIdx = ~isnan(valVec);

    thetaVec = thetaVec(validIdx);
    zVec     = zVec(validIdx);
    valVec   = valVec(validIdx);

    % Store data
    Frame      = [Frame; ff*ones(numel(valVec),1)];
    Time_s     = [Time_s; t_dT(ff)*ones(numel(valVec),1)];
    Theta_deg  = [Theta_deg; thetaVec];
    Z_m        = [Z_m; zVec];
    Value      = [Value; valVec];

end

% Create output table
T = table(Frame,Time_s,Theta_deg,Z_m,Value);

% Write CSV file
writetable(T,csvFileName);

disp(['CSV file written: ',csvFileName])

%% 4.4- Read CSV file and plot for cross-check:
% =========================================================================
% This verifies that the exported CSV reproduces the stitched theta-z map.

frCheck = min(55,size(z,3));

% Read CSV
Tcheck = readtable(csvFileName);

% Select one frame
idx = Tcheck.Frame == frCheck;

thetaCSV = Tcheck.Theta_deg(idx);
zCSV     = Tcheck.Z_m(idx);
valueCSV = Tcheck.Value(idx);

% Rebuild 2D grid from CSV points
thetaUnique = unique(thetaCSV);
zUnique     = unique(zCSV);

[thetaCSV_2D,zCSV_2D] = meshgrid(thetaUnique,zUnique);

valueCSV_2D = griddata(thetaCSV,zCSV,valueCSV, ...
    thetaCSV_2D,zCSV_2D,'linear');

% Plot original stitched map and CSV-reconstructed map
figure('color','w')

subplot(1,2,1)
surf(phi_2D,s_2D,z(:,:,frCheck),'LineStyle','none')
set(gca,'XTick',0:45:360,'XDir','reverse')
xlabel('Angle [deg]','Interpreter','Latex','FontSize',14)
ylabel('z [m]','Interpreter','Latex','FontSize',14)
title(['Original stitched map, frame ',num2str(frCheck)])
view([0,90])
axis tight
xlim([0,360])
caxis([0,20])
colormap('hot')
colorbar

subplot(1,2,2)
surf(thetaCSV_2D,zCSV_2D,valueCSV_2D,'LineStyle','none')
set(gca,'XTick',0:45:360,'XDir','reverse')
xlabel('Angle [deg]','Interpreter','Latex','FontSize',14)
ylabel('z [m]','Interpreter','Latex','FontSize',14)
title(['CSV reconstructed map, frame ',num2str(frCheck)])
view([0,90])
axis tight
xlim([0,360])
caxis([0,20])
colormap('hot')
colorbar