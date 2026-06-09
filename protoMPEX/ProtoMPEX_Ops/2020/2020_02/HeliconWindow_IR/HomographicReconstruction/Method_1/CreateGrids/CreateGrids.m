% Objective:
% Take the previosly made composite IR videos and construct the associated
% cylindrical grids and then produce the 2D time-dependent temperature maps

clear all
close all

% Select shot series to analyze:
% =========================================================================
shotsToAnalyze = 8;
saveData = 1;
saveFig  = 1;

for kk = shotsToAnalyze
% Load data:
% =========================================================================
homeAddress = cd;
cd .. 
cd CompositeIRdata\
fileName = ['dT_comp_shotSeries_',num2str(kk)];
load(fileName)
cd(homeAddress)

% Object's dimensions in physical space:
% =========================================================================
% Radius of window:
Ra = 0.90*0.5*12.4*1e-2; % [m]
% Length of window:
bb = 0.64;
La = bb*12.11*2.54*1e-2; % [m]

t_rng{1} = [180+25,360-25];
t_rng{2} = [25    ,180-25];

t_rng{1} = [180,360];
t_rng{2} = [0  ,180];

for n = 1:2
    % Define a 2D grid on the object's inner surface and wrt object's datum:
    % =========================================================================
    % Azimuthal angle:
    t_1D{n} = linspace(t_rng{n}(1),t_rng{n}(2),180)*pi/180;
    % Along axis of window:
    z_1D{n} = linspace(+0,+La,100)';
    % 2D grid 
    [t_2D{n},z_2D{n}] = meshgrid(t_1D{n},z_1D{n});

    % Create 3D grid of window's inner surface wrt to objects datum:
    % =========================================================================
    qx{n} = Ra*cos(t_2D{n});
    qy{n} = Ra*sin(t_2D{n});
    qz{n} = z_2D{n};
end

% Pin hole camera factors:
% =========================================================================
% Focal length:
C.f = 24.6*1e-3;
% Reflection factor
C.Omega = -1;

%% Create 2D grid:
% =========================================================================

n = 1;
% Initial guess:
% =========================================================================
% Axial location of camera relative to object's datum:
L{n}(3) = -34.5*1e-2; % [m]
% Rotation about yy to rotate the objects's frame to the camera's frame:
a{n} = -0*0.2*pi/180;
% Rotation about xx to rotate the objects's frame to the camera's frame:
b{n} = +13*pi/180;

% Reference points:
% =========================================================================
refType = 1;
switch refType
    case 1
        % Ref point on cylinder:
        q_star{n} = [+0         ,-Ra        ,0  ]';
        % Ref point on image plane:
        r_star{n} = [-0.102*1e-3,-1.92*1e-3 ,C.f ]';
    case 2
        % Ref point on cylinder:
        q_star{n} = [+0         ,-Ra         ,La  ]';
        % Ref point on image plane:
        r_star{n} = [-0.008*1e-3,+0.825*1e-3 ,C.f ]';
end

n = 2;
% Initial guess:
% =========================================================================
% Axial location of camera relative to object's datum:
L{n}(3) = -34.5*1e-2; % [m]
% Rotation about yy to rotate the objects's frame to the camera's frame:
a{n}    = +0.7*pi/180;
% Rotation about xx to rotate the objects's frame to the camera's frame:
b{n}    = -15*pi/180;

% Reference points:
% =========================================================================
refType = 2;
switch refType
case 1
    % Ref point on cylinder:
    q_star{n} = [+0         ,-Ra        ,0  ]';
    % Ref point on image plane:
    r_star{n} = [-0.102*1e-3,-1.92*1e-3 ,C.f ]';
case 2
    % Ref point on cylinder:
    q_star{n} = [+0         ,+Ra         ,0   ]';
    % Ref point on image plane:
    r_star{n} = [-0.008*1e-3,+1.115*1e-3 ,C.f ]';
end

for n = 1:2
    % Create rotation matrix:
    % =========================================================================
    % Rotation matrix defined as: [e] = Ryx*[s]
    % where [e] is camera's referene frame
    % [s] is object's reference frame
    R1    = [+cos(a{n})          ; +0        ; -sin(a{n})          ];
    R2    = [+sin(a{n})*sin(b{n}); +cos(b{n}); +cos(a{n})*sin(b{n})];
    R3    = [+sin(a{n})*cos(b{n}); -sin(b{n}); +cos(a{n})*cos(b{n})];
    R{n}  = [R1,R2,R3];

    % Calibration matrix:
    % =========================================================================
    M = r_star{n}*R3' + C.f*C.Omega*R{n}';
    G = M(1:2,1:2);
    H = -M(1:2,3);
    P = inv(G)*H;

    % Define offset vector:
    % =========================================================================
    L{n}(1) = q_star{n}(1) + P(1)*(L{n}(3) - q_star{n}(3));
    L{n}(2) = q_star{n}(2) + P(2)*(L{n}(3) - q_star{n}(3));
    L{n} = L{n}';
end

% Project object's grid to image plane:
% =========================================================================
t0 = tic;
for n = 1:2
    for jj = 1:numel(t_1D{n})
         for kk = 1:numel(z_1D{n})
             qq = [qx{n}(kk,jj),qy{n}(kk,jj),qz{n}(kk,jj)]';
             rr = PinHoleCamera(C,R{n},L{n},qq);
             xi{n}(kk,jj) = rr(1);
             yi{n}(kk,jj) = rr(2);
             zi{n}(kk,jj) = rr(3);
         end
    end
end
t0 = toc(t0);
disp(['Transforming 3D grid into image plane grid: ',num2str(t0),' [s]'])

% Plot Composite image:
% =========================================================================

figure('color','w')
fr = 55;
for n = 1:2
    subplot(1,2,n)
    hold on
    surf(y_cI{n}*1e3,x_cI{n}*1e3,dT_comp{n}(:,:,fr),'LineStyle','none')
    [rows,cols,~] = size(dT_comp{n});
    box on
    view([0,90])
    xlabel('y [mm]')
    ylabel('x [mm]')
    set(gca,'PlotBoxAspectRatio',[cols/rows 1 1])
    xlim([-1,1]*max(y_cI{n}*1e3))
    ylim([-1,1]*max(x_cI{n}*1e3))
    colormap('hot')
    % Draw coordinate system:
    line([0,0.1],[0,0  ],[1,1]*1e2,'color','k','LineWidth',2)
    text(0.1,0,300,'$e_2$','FontSize',12,'Interpreter','latex')
    line([0,0  ],[0,0.1],[1,1]*1e2,'color','k','LineWidth',2)
    text(0,0.1,300,'$e_1$','FontSize',12,'Interpreter','latex')
    
    % Plot 3D grid to image plane:
    % =========================================================================
    nn = 5;
    xx = xi{n}(1:nn:end,1:nn:end);
    yy = yi{n}(1:nn:end,1:nn:end);
    zz = zi{n}(1:nn:end,1:nn:end);
    plot3(yy*1e3  ,xx*1e3 ,zz*500 ,'w-')
    plot3(yy'*1e3 ,xx'*1e3,zz'*500,'w-')
    caxis([0,30]) 
end

if saveFig
    meanRF = mean(metadata.rfPwr);
    figureName = ['Grid_',metadata.limitMode,'_Limit_',num2str(round(meanRF),3),'_kW'];
    saveas(gcf,figureName,'tiffn')
end


%%
% Use cylindrical 3D grid to interpolate image:
% =========================================================================
fr = 55;
frMin = min([size(dT_comp{1},3),size(dT_comp{2},3)]);
for n = 1:2;
    t0 = tic;
    [Y,X] = meshgrid(y_cI{n},x_cI{n});
    for fr = 1:frMin
        v{n}(:,:,fr) = interp2(Y,X,dT_comp{n}(:,:,fr),yi{n},xi{n});
    end
    disp('Using interp2 function ...')
    t0 = toc(t0);
    disp(['Interpolating data: ',num2str(t0),' [s]'])
end
%%

% Plot data:
% =========================================================================
fr = 34;
figure('color','w')
hold on
for fr = 2:4:80
    surf(t_2D{1}*180/pi,(z_2D{1}/bb)*1e2,v{1}(:,:,fr),'LineStyle','none')
    surf(t_2D{2}*180/pi,(z_2D{2}/bb)*1e2,v{2}(:,:,fr),'LineStyle','none')
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

%%
if saveFig
    meanRF = mean(metadata.rfPwr);
    figureName = ['TemperatureDist_',metadata.limitMode,'_Limit_',num2str(round(meanRF),3),'_kW'];
    saveas(gcf,figureName,'tiffn')
end
disp('End of script!!')
end