% Objective:
% Constrain the image plane grid using the two-conjugate points method

% Context:
% We have developed a method that relies on two conjugate points in order
% to fully constrain the locatio of the object's grid
% In this script, we plot the solution for the viewing angles alpha and
% beta but we do not use a root-finder to extract the solution. The
% extractio of the roots is done in the next version of this script test3

% For all practical purposes, this is not the latest developement of the
% code!


%% Process: 
% 1- Load shot series.mat files:

clc
clear all
close all

% Select shotSeries and view:
% =========================================================================
kk = 1;
jj = 1;

% Load data:
% =========================================================================
homeAddress = cd;
cd ..
cd Method_2
fileName = ['ShotSeriesData_',num2str(kk)];
dat = load(fileName);
cd(homeAddress)

%%
% Object's dimensions in physical space:
% =========================================================================
% Radius of window:
aa = 1;
Ra = aa*0.5*12.4*1e-2; % [m]
% Length of window:
bb = 1;
La = bb*12.11*2.54*1e-2; % [m]

% Pin hole camera factors:
% =========================================================================
% Focal length:
C.f = 24.6*1e-3;
% Reflection factor
C.Omega = -1;
% focal length and reflection factor lumped together:
m = C.f*C.Omega;
% "z" offset of camera:
L3 = -34.5*1e-2;
    
% Creat chip coordinates:
% =========================================================================
xI = dat.u{jj}.xI;
yI = dat.u{jj}.yI;

% Create image:
% =========================================================================
fr = 55;
im = dat.u{jj}.dT(:,:,fr);

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
r_star.A = [-2.0e-3,+1.2e-3,C.f]';
q_star.A = [0      ,Ra     ,0  ]';

% Point B:
r_star.B = [-1.938e-3,-1.65e-3,C.f]';
q_star.B = [0        ,Ra      ,La ]';

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

% Assemble G:
% =========================================================================
for pp = 1:2
    k = fld(pp);
    G.(k) = @(a,b) [M{1,1}.(k)(a,b), M{1,2}.(k)(a,b);...
                    M{2,1}.(k)(a,b), M{2,2}.(k)(a,b)];
                
    detG.(k) = @(a,b)  M{1,1}.(k)(a,b)*M{2,2}.(k)(a,b) -  M{2,1}.(k)(a,b)*M{1,2}.(k)(a,b);                    
    
    K.(k) = @(a,b) Q3.(k)./detG.(k)(a,b);
    
    g.(k) = @(a,b) inv(G.(k)(a,b));             
end

% Assemble H:
% =========================================================================
for pp = 1:2
    k = fld(pp);
    H.(k) = @(a,b) -Q3.(k)*[M{1,3}.(k)(a,b); M{2,3}.(k)(a,b)];
end

% Assemble F1 and F2:
% =========================================================================
F1A = @(a,b) K.A(a,b)*(+M{2,2}.A(a,b)*M{1,3}.A(a,b) - M{1,2}.A(a,b)*M{2,3}.A(a,b));
F1B = @(a,b) K.B(a,b)*(+M{2,2}.B(a,b)*M{1,3}.B(a,b) - M{1,2}.B(a,b)*M{2,3}.B(a,b));
F1  = @(a,b) F1A(a,b) - F1B(a,b);

F2A = @(a,b) K.A(a,b)*(-M{2,1}.A(a,b)*M{1,3}.A(a,b) + M{1,1}.A(a,b)*M{2,3}.A(a,b));
F2B = @(a,b) K.B(a,b)*(-M{2,1}.B(a,b)*M{1,3}.B(a,b) + M{1,1}.B(a,b)*M{2,3}.B(a,b));
F2  = @(a,b) F2A(a,b) - F2B(a,b);

F = @(x) [F1(x(1),x(2)),F2(x(1),x(2))];

tic
[x, ~, ~, exitflag, output, ~] = newtonraphson(F, [3,-10]*pi/180);
toc

x*180/pi

for pp = 1:2
    k = fld(pp);
    l{pp} = u.(k) - inv(G.(k)(x(1),x(2)))*H.(k)(x(1),x(2));
end

l{1}*100

% 
% % Assemble F:
% % =========================================================================
% ex = [1,0]';
% ey = [0,1]';
% F  = @(a,b) g.A(a,b)*H.A(a,b) - g.B(a,b)*H.B(a,b) ;
% F1 = @(a,b) g.A(a,b)*Q3.A*M{1,3}.A(a,b) - g.B(a,b)*Q3.B*M{1,3}.B(a,b);
% F2 = @(a,b) g.A(a,b)*Q3.A*M{2,3}.A(a,b) - g.B(a,b)*Q3.B*M{2,3}.B(a,b);
% 

a1  = 4.25*pi/180; % [Rad]
b1  = -11.51*pi/180; % [Rad]

a = linspace(-20,20,100)*pi/180;
b = linspace(-20,20,100)*pi/180;

tic
for rr = 1:numel(a)
    for cc = 1:numel(b)
        f1(rr,cc) = F1(a(rr),b(cc));
        f2(rr,cc) = F2(a(rr),b(cc));
    end
end
toc

figure;
subplot(1,2,1)
hold on
plot(b*180/pi,abs(f1(:,1)))
plot(b*180/pi,abs(f1(:,50)))
plot(b*180/pi,abs(f1(:,100)))
xlabel('\beta')
subplot(1,2,2)
hold on
plot(a*180/pi,abs(f2(1,:)))
plot(a*180/pi,abs(f2(50,:)))
plot(a*180/pi,abs(f2(100,:)))
xlabel('\alpha')

figure;
subplot(1,2,1)
surf(a*180/pi,b*180/pi,abs(f1),'LineStyle','none')
xlabel('\alpha')
ylabel('\beta')
view([0,90])
axis square
title('f1')
subplot(1,2,2)
surf(a*180/pi,b*180/pi,abs(f2),'LineStyle','none')
xlabel('\alpha')
ylabel('\beta')
view([0,90])
axis square
title('f2')

figure;
hold on
surf(a*180/pi,b*180/pi,abs((f1-f2)).*abs((f1+f2)),'LineStyle','none')
% surf(a*180/pi,b*180/pi,abs((f1+f2)),'LineStyle','none')

figure
surf(a*180/pi,b*180/pi,abs(((f1+f2).^0.5)).*abs(((f1-f2).^0.5)),'LineStyle','none')
xlabel('\alpha')
ylabel('\beta')

[aa,bb] = meshgrid(a,b);

figure; hold on
surf(aa*180/pi,bb*180/pi,abs(F1(aa,bb)),'LineStyle','none')
surf(aa*180/pi,bb*180/pi,abs(F2(aa,bb)),'LineStyle','none')

xlabel('\alpha')
ylabel('\beta')

return

for pp  = 1:2
    M11{pp} = @(a,b) r_star{pp}(1)*R13(a,b) + m*R11(a,b);
    M12{pp} = @(a,b) r_star{pp}(1)*R23(a,b) + m*R21(a,b);
    M13{pp} = @(a,b) r_star{pp}(1)*R33(a,b) + m*R31(a,b);
    
    M21{pp} = @(a,b) r_star{pp}(2)*R13(a,b) + m*R12(a,b);
    M22{pp} = @(a,b) r_star{pp}(2)*R23(a,b) + m*R22(a,b);
    M23{pp} = @(a,b) r_star{pp}(2)*R33(a,b) + m*R32(a,b);
end

% K matrix:
% =========================================================================
K11 = @(a,b) M11{2}(a,b).*M22{1}(a,b) - M12{2}(a,b).*M21{1}(a,b);
K21 = @(a,b) M21{2}(a,b).*M22{1}(a,b) - M22{2}(a,b).*M21{1}(a,b);
K12 = @(a,b) M12{2}(a,b).*M11{1}(a,b) - M11{2}(a,b).*M12{1}(a,b);
K22 = @(a,b) M22{2}(a,b).*M11{1}(a,b) - M21{2}(a,b).*M12{1}(a,b);

% Determinant term:
% =========================================================================
A = @(a,b) M11{1}(a,b).*M22{1}(a,b) - M21{1}(a,b).*M12{1}(a,b);

% Trascendental equations:
% =========================================================================
F1 = @(a,b) A(a,b).*M13{2}(a,b).*Q3{2} - (K11(a,b).*M13{1}(a,b) + K12(a,b).*M23{1}(a,b)).*Q3{1};
F2 = @(a,b) A(a,b).*M23{2}(a,b).*Q3{2} - (K21(a,b).*M13{1}(a,b) + K22(a,b).*M23{1}(a,b)).*Q3{1};

%%
fun = @(x) [F1(x(1),x(2)),F2(x(1),x(2))];
options = optimset('TolX',1e-12); 
[x, ~, ~, exitflag, output, ~] = newtonraphson(fun, [11.75,-4.15]*pi/180,options);


%%
% a: Rotation about yy to rotate the objects's frame to the camera's frame
% b: Rotation about xx to rotate the objects's frame to the camera's frame

a = linspace(-0,16,1e3)*pi/180;
b = linspace(-6,0,1e3)*pi/180;

a = linspace(-45,45,1e3)*pi/180;
b = linspace(-45,45,1e3)*pi/180;

[aa,bb] = meshgrid(a,b);

figure; hold on
surf(aa*180/pi,bb*180/pi,abs(F1(aa,bb)),'LineStyle','none')
surf(aa*180/pi,bb*180/pi,abs(F2(aa,bb)),'LineStyle','none')

xlabel('\alpha')
ylabel('\beta')
% zlim([0,0.5e-6])
% caxis([0,0.5e-6])


%%
a = 11.99;
b = -3.68;

for pp = 1:2;
    G{pp} = [+M11{pp}(a,b),+M12{pp}(a,b);...
         +M21{pp}(a,b),+M22{pp}(a,b)];
    l{pp} = v{pp} + inv(G{pp})*[M13{pp}(a,b);M23{pp}(a,b)]*Q3{pp};
end

l{1}
l{2}

%%
return

for kk = 1


%% 1- Load shot series.mat files:
% =========================================================================
cd ..
cd Method_2
fileName = ['ShotSeriesData_',num2str(kk)];
load(fileName)
cd(homeAddress)

%% 2- Calculate grids and project object to image plane:
% =========================================================================
close all

% Object's dimensions in physical space:
% =========================================================================
% Radius of window:
aa = 1;
Ra = aa*0.5*12.4*1e-2; % [m]
% Length of window:
bb = 1;
La = bb*12.11*2.54*1e-2; % [m]

% Pin hole camera factors:
% =========================================================================
% Focal length:
C.f = 24.6*1e-3;
% Reflection factor
C.Omega = -1;
% "z" offset of camera:
L3 = 34.5*1e-2;


end


return

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
end