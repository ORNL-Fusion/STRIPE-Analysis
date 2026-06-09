% Objective:
% Preview the IR data as fast as practically possible for use during
% experiments

clear all
close all
clc

% Select shot series to plot:
% =========================================================================
xp = 1;
rotate_flag = 1;
switch xp
    case 1 
        shot       = 30000 + [956];
        addressLoc =         [ 1];        
end

%% Load seq files and get Temperature data:
% =========================================================================
a{1} = ['\\mpexserver\ProtoMPEX_Data\IR_Camera\2020_11_24'];
a{2} = ['\\mpexserver\ProtoMPEX_Data\IR_Camera\2020_11_17'];
a{3} = ['\\mpexserver\ProtoMPEX_Data\IR_Camera\2020_11_20'];

disp('Extracting data ...')
dum1 = tic;
for s = 1:length(shot)
    % Specify IR file name
    pathName = [a{addressLoc(s)},'\'];
    fileName = ['Shot ',num2str(shot(s)),'.seq'];
    fileName2 = ['Shot-0',num2str(shot(s)),'.seq'];

    % Set thermal parameters:
    thermalParam.ExternalOpticsTransmission = 0.7;
    thermalParam.AtmosphericTemperature = 24;
    thermalParam.Distance = 1;
    thermalParam.ExternalOpticsTemperature = 24;
    thermalParam.ReferenceTemperature = 24;
    thermalParam.Transmission = 1;
    thermalParam.RelativeHumidity = 0;
    thermalParam.ReflectedTemperature = 24;
    thermalParam.emissivity = 0.3;

    % Define extraction options:
    extractOptions.frames = 40:1:250;
    extractOptions.frameRate = 100;
    
    % Get the data:
    disp(['Extracting data from ',fileName,' ...'])
    try
        [g{s},t_f{s},seq{s}] = ExtractDataFromSeqFile(fileName,pathName,thermalParam,extractOptions);
    catch
       [g{s},t_f{s},seq{s}] = ExtractDataFromSeqFile(fileName2,pathName,thermalParam,extractOptions);
        warning('Filename2 was used for the IR file...')
    end
    
    % Display temperature range selected:
    tempRange = seq{s}.ThermalImage.CameraInformation.Range;
    disp(['Temperature range of shot #',num2str(shot(s)),' is: ',num2str(tempRange.Minimum),' to ',num2str(tempRange.Maximum),' deg'])
       
    % Crop data:
    for ii = 1:size(g{s},3)
        f{s}(:,:,ii) = g{s}(:,130:400,ii);
    end
        
    % Extract data during plasma pulse:
    options.n_Before = 10;
    options.n_After = 1;
    [intf{s},rng{s},t0(s)] = GetDataDuringPlasma(f{s},t_f{s},options);
    rng{s} = 1:size(f{s},3);
          
    % Convert intensity to temperature:
    p = IntensityTempConv(thermalParam.emissivity,f{s}(:,:,rng{s}),seq{s});
    t_temp{s} = t_f{s}(rng{s});
    
    % Rotate image:
    if rotate_flag
        angle = -46;
        [r] = rotate_IR_data(p,angle);
        temp{s} = r;
    else
        temp{s} = p;
    end   
end

disp('Data extraction completed!')
disp(['Time taken: ',num2str(toc(dum1))])

% Plot the integral of the intensity over time:
if 1
    figure('color','w')
    grid on
    hold on
    lineColor = {'k','bl','r','g','m','c','k','bl','r','g','m','c'};
    for s = 1:length(shot)
               plot(intf{s}         - min(intf{s}),lineColor{s})
        h(s) = plot(rng{s},intf{s}(rng{s}) - min(intf{s}),lineColor{s},'LineWidth',3);
        legendText{s} = num2str(shot(s));
    end
    legend(h,legendText)
    box on
end

%% Select Crop area:
close all
s = 1;
[~,frame] = max(intf{1}(rng{1}));
frame = 100:2:175

% Define range of pixels to extract
col_rng = [125:245];
row_rng = [106:249];

col_rng = [145:282] + 10 ;
row_rng = [128:275] - 10;

disp(['col_rng size: ', num2str(numel(col_rng))])
disp(['row_rng size: ', num2str(numel(row_rng))])

% Plot region to be cropped:
options.colorbar = 1;
options.magnitudePlotMode = 1;
options.removeAxesTicks = 0;
options.shot = shot(s);
options.mirrorImage = 0;
options.frames = frame;
options.zlim = [0,10];
PlayMovieFromArray(temp{s},options)
hold on

% Define area to be cropped:
x1 = col_rng(1);
y1 = row_rng(1);
drow = max(row_rng) - min(row_rng);
dcol = max(col_rng) - min(col_rng);
x2 = x1 + dcol;
y2 = y1 + drow;

% Draw area to be cropped:
hdum = line([x1,x1,x2,x2,x1],[y1,y2,y2,y1,y1],40*ones(1,5));
hdum.Color = 'k';
hdum.LineWidth = 2;
hdum.LineStyle = '--';
%% Crop data
try 
    clear temp_c
end
for s = 1:length(shot)
        for ii = 1:size(temp{s},3)
            temp_c{s}(:,:,ii) = temp{s}(row_rng,col_rng,ii);
            dT{s}(:,:,ii) = temp_c{s}(:,:,ii) - temp_c{s}(:,:,1);
        end    
        t_dT{s} = t_temp{s} - t_temp{s}(1);
end

if 1
    s = 1;
    [~,frame] = max(intf{s}(rng{1}));
    options.colorbar = 1;
    options.magnitudePlotMode = 1;
    options.removeAxesTicks = 0;
    options.shot = shot(s);
    options.mirrorImage = 0;
    options.frames = frame;
    options.zlim = [0,150];
    PlayMovieFromArray(temp_c{s},options)
end

%% Plot data
% close all
% Select shot:
% =========================================================================
s = 1;
% Movie rendering options:
% =========================================================================
try 
    clear options
end
options.colorbar = 1;
options.magnitudePlotMode = 1;
options.removeAxesTicks = 0;
options.shot = shot(s);
options.mirrorImage = 0;

% Select quantity to plot:
% =========================================================================
if 1
    plotType = 2;
    switch plotType
        case 1
            options.frames = 95:1:size(temp{s},3)-25;
            options.zlim = [0,150];
            PlayMovieFromArray(temp_c{s},options)
        case 2
            options.frames = 95:2:size(temp{s},3)-2;
            options.zlim = [0,max(max(dT{s}(:,:,frame)))];
            PlayMovieFromArray(dT{s},options)
    end
end
%% Normalize data:
s = 1;
% W plate Material properties: ------------------------------------
rho = 19300; 
kt   = 173;
cp  = 134;
% Thermal diffusivity
a = kt/(rho*cp);
% Order of magnitude heat flux
q0 = 1e6;
% Plate thickness
Lz = 0.75e-3;
% Plate width:
Lx = 90e-3;
% Plate height:
Ly = 90e-3;

% Characteristic time scale
t_star = Lz*Lz/a;
% Characteristic temperature
T_star = q0*Lz/kt;

% Calculate normalized time:
t_R = t_dT{1}'/t_star;
dt_R = mean(diff(t_R));

% Initial temperature distribution:
T0 = dT{s}(:,:,1);

% Define size of data:
 [I,J,R] = size(dT{s});

% Calculate normalized temperature
for r = 1:R
    u_IJR(:,:,r) = (dT{s}(:,:,r) - T0)./T_star;
end

%% Compute 2D Fourier series representation:
x_J = linspace(0,1,J)';
y_I = linspace(0,1,I)';

% Define the number of Fourier modes:
% "x":
M = 60;
% "y":
N = 60;
% "z":
L = 70;

% Define the number of "m" and "n" terms to include in the Fourier cosine
% expansion:
m_M = [0:1:(M-1)]';
n_N = [0:1:(N-1)]';
l_L = [0:1:L-1]';

tic
for r = 1:length(t_R)
    [a1,a2,a3,a4] = FourierCosine_2D_v3(x_J,y_I,u_IJR(:,:,r),M,N);
    % "fu_IJR" is the Fourier-cosine approximation of "u_IJR"  
    fu_IJR(:,:,r) = a1;
    % "U_NM" are the time dependent Fourier cosine coefficients of "u"
    U_NM{r} = a2;
end
% Phi_JM is the eigenfunction matrices for all "x" and all "m"
Phi_JM = a3;
% Phi_IN is the eigenfunction matrices for all "y" and all "n"
Phi_IN = a4;
toc

% Compare input and Fourier representation:
s = 1;
[~,fr] = max(intf{s}(rng{s}));
xx = Lx*(x_J - 0.5);
yy = Ly*(y_I - 0.5);

figure
subplot(1,2,1)
surf(xx,yy,u_IJR(:,:,fr)*T_star,'LineStyle','none')
view([0,90])
axis('square')
title('Input')

subplot(1,2,2)
surf(xx,yy,fu_IJR(:,:,fr)*T_star,'LineStyle','none')
view([0,90])
axis('square')
title('Fourier representation')

%% Test Fourier representation:
if 1
    figure; 
    set(gcf,'color','w')
    maxU = max(max(U_NM{fr}));
    for r = 90:2:length(t_R)
        surf(m_M,n_N,abs(U_NM{r}),'LineStyle','none')
        caxis([0,maxU/10])
        zlim([0,maxU/10])
        view([0,90])
        title(['Temporal evolution of Fourier-cosine coefficients, fr: ',num2str(r)])
        xlabel('m mode')
        ylabel('n mode')
        drawnow
        pause(0.1)
    end

    % Compare the measured data "u"  with the MxN term Fourier-cosine represenation
    % of "u"
    s = 1;
    [~,r] = max(intf{s}(rng{s}));   
    figure; 
    set(gcf,'color','w')
    subplot(2,2,1);
    surf(x_J,y_I,u_IJR(:,:,r)*T_star,'LineStyle','none')
    view([0,90])
    axis('square')
    zlim([0,300*T_star])
    title('Measured Temp data')

    subplot(2,2,2);
    surf(x_J,y_I,fu_IJR(:,:,r)*T_star,'LineStyle','none')
    view([0,90])
    axis('square')
    zlim([0,300*T_star])
    title('Fourier-cosine representation')

    subplot(2,2,[3:4])
    MaxValue = max(max(u_IJR(:,:,r)));
    surf(x_J,y_I,100*((fu_IJR(:,:,r)-u_IJR(:,:,r))./MaxValue),'LineStyle','none')
    view([0,90])
    colorbar
    zlim([1e-4,50])
    caxis([0,50])
    axis('square')
    title('% Error in approximation')
end

%% 
%==========================================================================
% Assemble U_NM 
%==========================================================================

% Assemble U_NM in a manner that is suitable for the next steps
for r = 1:R
    for n = 1:N
        for m = 1:M
            U_R{n,m}(r,1) = U_NM{r}(n,m);
        end
    end
end

%% 
%==========================================================================
% Calculate Toeplitz matrix
%==========================================================================
% The purpose of the Toeplitz matrix is to perform convolution on a vector.
% Its the finite-dimensional representation of the convolution integral
% operator

% P_RR{N,M,K}
% This cell has NxMxK matrices.
% This consumes the largest memory in the calculation
% By evaluating only one "z_K" location we can reduce the size to NxM
% P_RR_new = P_RR{:,:,k};

% For Back-side  imaging use z_K = 1;
% For Front-side imaging use z_K = 0;
z_K = 1;
% Calculate Toeplitz matrix for a single value of z_K
[P_RR] = Toeplitz_NM(m_M,n_N,l_L,Lx,Ly,Lz,z_K,t_R);

%%
%==========================================================================
% Conjugate gradient method to find surface heat flux that produced the
% surface temperature distribution that was measured with IR camera
%==========================================================================
% In this section we use the time dependenent 2D Fourier cosine
% coefficients "U_NM" of the measured temperature "u" and solve the inverse
% heat conduction problem for each "m" and "n" mode. The outcome of this is
% to produce the time-dependent 2DFourier cosine coefficients of the
% surface heat flux Q_NM. We then use Q_MN in a 2D Fourier cosine series to
% compute the time dependent surface heat flux q_IJR

% This process is essentially a minimization/optimization problem.
% we use the Conjugate gradient method:
% Input data:
% P_RR: convolution operator
% Q_iteration{1}: initial guess at "m,n" fourier cosine heat flux coefficient, column vector
% U_data: "m,n" Fourier-cosine coeficient of experimental data, column vector
% Ni: number of iteration
% dtt: time step, need to make = 1 since we have accounted for this in the P_RR operator
% ExitCond: Exit condition, 1 or 2
% Val: Value associated with ExitCond

% Output data:
% U_dataInferred: calculated "m,n" Fourier-cosine coeficient of experimental data, column vector
% Residual: "m,n" residual, column vector
% Q_iteration: minimizing "m,n" fourier cosine heat flux coefficient, structure, column vector

% Form:
% [Q_iteration,U_dataInferred,Residual] = IHCP_ConjugateGradient(P_RR{m,n},Q_iteration,U_data,dtt,Ni,ExitCond,Val)


% Apply the conjugate gradient method for all "m" and "n" modes:
% Time step
dtt = 1; % Always equal to 1
% Set the maximum # of iterations for the minimization search
Ni_a = 4; 
Ni_b = 2;

tic
for m = 1:M
    for n = 1:N
        % "m,n" Fourier-cosine coeficient of experimental data:
        U_data = U_R{n,m};
        
        if sqrt(m*m + n*n) < 120
            % Initialize "n,m" heat flux fourier-cosine coeffiecient guess:
            Q_iteration{1} = (diff(U_R{n,m})./dt_R);
            Q_iteration{1} = [Q_iteration{1};Q_iteration{1}(end)];
  
            % Inverse method:
            [Q_iteration,U_dataInferred,Residual] = ...
            IHCP_ConjugateGradient(P_RR{n,m},Q_iteration,U_data,dtt,Ni_a,3,[]);
        else
            Q_iteration{1} = zeros(size(t_R));
            [Q_iteration,U_dataInferred,Residual] = ...
            IHCP_ConjugateGradient(P_RR{n,m},Q_iteration,U_data,dtt,Ni_b,3,[]);
        end
        
        % Assign last iteration as solution:
        Q_R{n,m} = Q_iteration{end};
        UdInfer_R{n,m} = U_dataInferred{end};
        Res{n,m} = Residual;
            
        % Clear dummy variables before proceeding to the next "m,n" coeficient
        clear Q_iteration U_dataInferred Residual
    end
    disp(['m: ',num2str(m)])
end

% At this point, we have calculated all the "m,n" Fourier-cosine
% coefficients of the surface heat flux
toc

%%
%==========================================================================
% Plot results:
%==========================================================================
% Q_R is the fourier coeficient as function of time for all "M" and "N" modes
% Plot Q_R{1,1} and compare it to d/dt(U_R{1,1})
m = 1;
n = 1;

figure; 
hold on
hQ(1) = plot(Q_R{n,m},'r.-','Markersize',8)
hQ(2) = plot((diff(U_R{n,m})./diff(t_R)),'g.-')
legend(hQ,'Q_R','\propto dU_R/dt')
titleText = ['Q_R{',num2str(n),',',num2str(m),'}'];
title(titleText)
box on
set(gcf,'color','w')

% Compare the measured T (fourier component) with the reconstructed one
figure; 
subplot(2,1,1)
hold on
hU(1) = plot(UdInfer_R{n,m},'LineWidth',2)
hU(2) = plot(U_R{n,m},'r.-','Markersize',8)
legend(hU,'Reconstructed u_R','U_R')
titleText = ['U_R{',num2str(n),',',num2str(m),'}'];
title(titleText)
box on

m = 1;
n = 2;
subplot(2,1,2)
hold on
hU(1) = plot(UdInfer_R{n,m},'LineWidth',2)
hU(2) = plot(U_R{n,m},'r.-','Markersize',8)
legend(hU,'Reconstructed u_R','U_R')
titleText = ['U_R{',num2str(n),',',num2str(m),'}'];
title(titleText)
box on
set(gcf,'color','w')

%%
% =========================================================================
% Arrange Q_R in a manner that can be used in a Fourier-cosine series
% =========================================================================
for m = 1:M
    for n = 1:N
            for r = 1:R
                if sqrt(m*m + n*n) > 120 && sqrt(m*m + n*n) < 130
                    Q_NM{r}(n,m) = 0.01*Q_R{n,m}(r);
                else
                    Q_NM{r}(n,m) = Q_R{n,m}(r);
                end
            end
    end
end

figure; 
set(gcf,'color','w')
Qmax = (max(max(Q_NM{fr})));
for r = 90:2:160
    surf(abs(Q_NM{r}),'LineStyle','none')
    zlim([0,Qmax*1.3])
    view([0,90])
    view([60,30])
    title(['Temporal evolution of Fourier-cosine coefficients, fr: ',num2str(r)])
    xlabel('m mode')
    ylabel('n mode')
    drawnow
    pause(0.1)
end

%==========================================================================
% Assemble heat flux function: --------------------------------------------
% Use a Fourier-cosine series to compute the surface heat flux
%==========================================================================
dx = Lx*(x_J(2)-x_J(1));
dy = Ly*(y_I(2)-y_I(1));
rng1 = 10:120;
for r = 1:R
        f = Phi_IN*Q_NM{r}*Phi_JM';
        q_IJR(:,:,r) = f;
        pwr(r) = sum(sum(f(:,:)))*dx*dy*q0;
end

% ########################################################################
% ########################################################################
% q_IJR is the solution to the entire problem
% This is the time-dependent surface heat flux that produces the measured 2D
% time-dependent surface temperatures measured with IR camera.
% ########################################################################
% ########################################################################

%%
%==========================================================================
% Plot results
%==========================================================================
s = 1;
zMax = 1*max(max(max(q_IJR)));
figure
set(gcf,'color','w')
[~,fr] = max(pwr);
% for r = [1:5:205,150];
% for r = 120:1:144;
for r = [90:2:175,140]
    zz = (q_IJR(:,:,r) + q_IJR(:,:,r+1) + q_IJR(:,:,r+3))/3; 
%     zz = q_IJR(:,:,r); 

    surf(xx*1e3,yy*1e3,zz,'LineStyle','none')
    colormap(flipud(hot))
    caxis([-0,zMax])
    zlim([-0,zMax])
    view([-45,45])
    view([00,90])
    line([-45,+45],[+45,+45],[1,1],'color','k','LineWidth',2)
    line([-45,+45],[-45,-45],[1,1],'color','k','LineWidth',2)
    line([+45,+45],[-45,+45],[1,1],'color','k','LineWidth',2)
    line([-45,-45],[-45,+45],[1,1],'color','k','LineWidth',2)
    
    colorbar
    set(gca,'XTick',[-50:10:50])
    title(['shot: ',num2str(shot(s)),' , frame:',num2str(r)])
    drawnow
    pause(0.01)
end

s = 1;
figure 
plot(pwr*1e-3,'k','LineWidth',2)
box on
ylabel('[kW]')
title(['shot: ',num2str(shot(s))])
ylim([0,1])
xlim([0,200])