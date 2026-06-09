% Step_B:
% Load temperature data
% Calculate heat flux

clear all
close all
clc

% Select experiment to load:
% =========================================================================
xp = 5;
rotate_flag = 1;

% Shot cases:
% =========================================================================
switch xp
    case 1
        shot       = 29000 + [777];
        date       = '2020_04_08';
        % Define range of pixels to extract:
        col_rng = [112:256];
        row_rng = [106:249];
        % Notes:
        % Thermal range: 100 to 650 deg     
    case 2
        shot       = 29000 + [778];
        date       = '2020_04_08';
        % Define range of pixels to extract:
        col_rng = [112:256];
        row_rng = [106:249];
        % Notes:
        % Thermal range: 100 to 650 deg      
    case 3
        shot       = 29000 + [824];
        date       = '2020_04_15';
        % Define range of pixels to extract:
        col_rng = [112:256];
        row_rng = [106:249];
        % Notes:
        % Thermal range: -40 to 150 deg
    case 4
        shot       = 29000 + [854];
        date       = '2020_04_15';
        % Define range of pixels to extract:
        col_rng = [112:256];
        row_rng = [106:249];
        % Notes:
        % Thermal range: -40 to 150 deg
    case 5
        shot       = 29000 + [855];
        date       = '2020_04_15';
        % Define range of pixels to extract:
        col_rng = [112:256];
        row_rng = [106:249];
        fr_substract = 160;
        % Notes:
        % Thermal range: -40 to 150 deg
    case 6
        shot       = 29000 + [691];
        date       = '2020_03_31';
        % Define range of pixels to extract:
        col_rng = [112:256];
        row_rng = [106:249];
        % Notes:
        % Thermal range: -40 to 150 deg
    case 7
        shot       = 30000 + [854];
        date       = '2020_11_17';
        % Define range of pixels to extract:
        col_rng = [112:256] + 40;
        row_rng = [106:249] + 12;
        fr_substract = 160;
        % Notes:
        % Thermal range: -40 to 150 deg
        % Target is slightly tilted so target appears narrower
    case 8
        shot       = 30000 + [698];
        date       = '2020_11_06';
        % Define range of pixels to extract:
        col_rng = [112:256] + 40;
        row_rng = [106:249] + 12;
        % Notes:
        % Thermal range: -40 to 150 deg
        % Target is slightly tilted so target appears narrower
end

% Check that col_rng and row_rng are square
disp(['Column dimension: ',num2str(numel(col_rng))]);
disp(['Row dimension: ',num2str(numel(row_rng))]);

% Create figure names:


% Load seq files and get Temperature data:
% =========================================================================
dataAddress = ['\\mpexserver\ProtoMPEX_Data\IR_Camera\',date];

disp('Extracting data ...')
dum1 = tic;
    % Specify IR file name:
    pathName = [dataAddress,'\'];
    fileName = ['Shot ',num2str(shot),'.seq'];

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
    [g,t_f,seq] = ExtractDataFromSeqFile(fileName,pathName,thermalParam,extractOptions);
    
    % Display temperature range selected:
    tempRange = seq.ThermalImage.CameraInformation.Range;
    disp(['Temperature range of shot #',num2str(shot),' is: ',num2str(tempRange.Minimum),' to ',num2str(tempRange.Maximum),' deg'])
    
    % Initial crop:
    for ii = 1:size(g,3)
        f(:,:,ii) = g(:,130:400,ii);
    end
        
    % Extract data during plasma pulse:
    options.n_Before = 10;
    options.n_After = 20;
    [intf,rngPlasma,t0] = GetDataDuringPlasma(f,t_f,options);
    dintf = diff(intf);
    rngPlasma = find(dintf>2,1) - options.n_Before;
    rngPlasma = [rngPlasma:size(f,3)- options.n_After];
    %     rngPlasma = 1:size(f,3);
          
    % Convert intensity to temperature:
    p = IntensityTempConv(thermalParam.emissivity,f(:,:,rngPlasma),seq);
    t_temp = t_f(rngPlasma);
    
    % Rotate image:
    if rotate_flag
        angle = -45;
        [r] = rotate_IR_data(p,angle);
        temp = r;
    else
        temp = p;
    end   


disp('Data extraction completed!')
disp(['Time taken: ',num2str(toc(dum1))])

% Plot the integral of the intensity over time:
% =========================================================================
if 1
    figure('color','w')
    grid on
    hold on
    lineColor = 'k';
    plot(intf         - min(intf),lineColor)
    h = plot(rngPlasma,intf(rngPlasma) - min(intf),lineColor,'LineWidth',3);
    legendText = num2str(shot);
    legend(h,legendText)
    box on  
end

%% Select Crop area:
% =========================================================================

% Select hottest frame:
[~,frame] = max(intf(rngPlasma));

% Select area:
if isempty(col_rng) && isempty(row_rng)
    [row_rng,col_rng,~] = size(temp);
    row_rng = 1:row_rng;
    col_rng = 1:col_rng;
end

% col_rng = [112:256] + 40;
% row_rng = [106:249] + 12;


% Plot region to be cropped:
options.colorbar = 1;
options.magnitudePlotMode = 1;
options.removeAxesTicks = 0;
options.shot = shot;
options.mirrorImage = 0;
options.frames = frame;
options.zlim = [0,10];
PlayMovieFromArray(temp,options)
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

%% Crop data:
% =========================================================================
try 
    clear temp_c
end


% Create cropped frames:
        for ii = 1:size(temp,3)
            temp_c(:,:,ii) = temp(row_rng,col_rng,ii);
            dT(:,:,ii) = temp_c(:,:,ii) - temp_c(:,:,1);
        end    
        t_dT = t_temp- t_temp(1);


% View cropped frame:
if 0
    s = 1;
    [~,frame] = max(intf(rngPlasma));
    options.colorbar = 1;
    options.magnitudePlotMode = 1;
    options.removeAxesTicks = 0;
    options.shot = shot;
    options.mirrorImage = 0;
    options.frames = frame;
    options.zlim = [0,150];
    PlayMovieFromArray(temp_c,options)
end

%% Plot data:
% =========================================================================

% Movie rendering options:
% =========================================================================
try 
    clear options
end
options.colorbar = 1;
options.magnitudePlotMode = 1;
options.removeAxesTicks = 0;
options.shot = shot;
options.mirrorImage = 0;

% Select quantity to plot:
% =========================================================================
if 1
    plotType = 2;
    switch plotType
        case 1
            options.frames = 95:1:size(temp,3)-25;
            options.zlim = [0,150];
            PlayMovieFromArray(temp_c,options)
        case 2
            options.frames = 100:2:size(temp,3)-2;
            options.zlim = [0,100];
            PlayMovieFromArray(dT,options)
    end
end

%% Normalize data:

% W plate Material properties:
plateMaterial = 'W';
rho = 19300; 
kt   = 173;
cp  = 134;

% Thermal diffusivity
a = kt/(rho*cp);

% Plate thickness
Lz = 0.75e-3;
% Plate width:
Lx = 90e-3;
% Plate height:
Ly = 90e-3;

% Characterstic heat flux:
q_star = 1e6;
% Characteristic time scale
t_star = Lz*Lz/a;
% Characteristic temperature
T_star = q_star*Lz/kt;

% Calculate normalized time:
t_R = t_dT'/t_star;
dt_R = mean(diff(t_R));

% Initial temperature distribution:
T0 = dT(:,:,1);

% Define size of data:
[I,J,R] = size(dT);

% Calculate normalized temperature
for r = 1:R
    u_IJR(:,:,r) = (dT(:,:,r) - T0)./T_star;
end

%% Compute 2D Fourier series representation:
% =========================================================================

% Create normalized coordinates for image:
% -------------------------------------------------------------------------
x_J = linspace(0,1,J)';
y_I = linspace(0,1,I)';

% Define the number of Fourier modes:
% -------------------------------------------------------------------------
% "x":
M = 80;
% "y":
N = 80;
% "z":
L = 100;

% Define the number of "m" and "n" terms to include in the Fourier cosine
% expansion:
% -------------------------------------------------------------------------
m_M = [0:1:(M-1)]';
n_N = [0:1:(N-1)]';
l_L = [0:1:L-1]';

% Calculate Fourier cosine coefficients:
% -------------------------------------------------------------------------
tic
for r = 1:length(t_R)
    [a1,a2,a3,a4] = FourierCosine_2D_v3(x_J,y_I,u_IJR(:,:,r),M,N);
    % "fu_IJR" is the Fourier-cosine approximation of "u_IJR"  
    fu_IJR(:,:,r) = a1;
    % "U_NM" are the time dependent Fourier cosine coefficients of "u"
    U_NM{r} = a2;
end

% Assign eigenfunction:
% -------------------------------------------------------------------------
% Phi_JM is the eigenfunction matrices for all "x" and all "m"
Phi_JM = a3;
% Phi_IN is the eigenfunction matrices for all "y" and all "n"
Phi_IN = a4;
toc

% Plot input temperature and fourier represetation:
% =========================================================================
figure('color','w')

% Font sizes:
fontSize.axes = 11;
fontSize.title = 13;
fontSize.legend = 11;
fontSize.label = 11;
fontSize.colorbar = 10;

% Plasma radius based on flux mapping:
Rp = 3.4;

% Create image coordinates:
xx = Lx*(x_J - 0.5);
yy = Ly*(y_I - 0.5);

% Choose the brigthest frame:
[~,fr] = max(intf(rngPlasma));

% Input temperature:
subplot(1,2,1)
hold on;
ax(1) = gca;
rngPlot = find(xx*1e2 > -4.4 & xx*1e2 <+4.4);
surf(xx(rngPlot)*1e2,yy(rngPlot)*1e2,u_IJR(rngPlot,rngPlot,fr)*T_star,'LineStyle','none')
view([0,90])
axis('image')
title('Measured temperature','interpreter','latex','fontSize',fontSize.title)
xlabel('x [cm]','interpreter','latex','fontSize',fontSize.label)
ylabel('y [cm]','interpreter','latex','fontSize',fontSize.label)
hC(1) = colorbar;
box on
% Add LUFS:
xCirc = Rp*cos(linspace(0,2*pi));
yCirc = Rp*sin(linspace(0,2*pi));
hCirc(1) = plot3(xCirc,yCirc,ones(size(xCirc))*1e3,'k','LineStyle',':','LineWidth',2);

% Fourier representation:
subplot(1,2,2)
hold on
ax(2) = gca;
surf(xx(rngPlot)*1e2,yy(rngPlot)*1e2,fu_IJR(rngPlot,rngPlot,fr)*T_star,'LineStyle','none')
view([0,90])
axis('square')
title('Fourier representation','interpreter','latex','fontSize',fontSize.title)
xlabel('x [cm]','interpreter','latex','fontSize',fontSize.label)
ylabel('y [cm]','interpreter','latex','fontSize',fontSize.label)
hC(2) = colorbar;
colormap((flipud(hot)));
box on
% Add LUFS:
hCirc(2) = plot3(xCirc,yCirc,ones(size(xCirc))*1e3,'k','LineStyle',':','LineWidth',2);

% Formatting:
set(ax,'fontName','Times','fontSize',fontSize.axes);
set(hC,'fontSize',fontSize.colorbar,'Location','northoutside')
xlim(ax,0.5*Lx*[-1,1]*1e2)
ylim(ax,0.5*Ly*[-1,1]*1e2)

% Visualize fourier components:
% =========================================================================
figure('color','w')
fontSize.axes = 12;
fontSize.title = 13;
fontSize.legend = 11;
fontSize.label = 12;
fontSize.colorbar = 12;

% Fourier spectrum:
Umn_abs = abs(U_NM{fr});
Umn_abs_max = max(max(Umn_abs));
hFC = pcolor(m_M,n_N,Umn_abs/Umn_abs_max);
xlim([0,20])
ylim([0,20])
view([0,90])
axis('square')
title('Fourier-cosine spectrum','interpreter','latex','fontSize',fontSize.title)
xlabel('m','interpreter','latex','fontSize',fontSize.label)
ylabel('n','interpreter','latex','fontSize',fontSize.label)
colormap((flipud(hot)));
hFCC = colorbar;
box on

% Formatting:=
set(gca,'fontName','Times','fontSize',fontSize.axes);
set(hFCC,'fontSize',fontSize.colorbar,'Location','Eastoutside')
hFC.LineStyle = ':';
%% Test Fourier representation:
% =========================================================================

% Compare the measured data "u"  with the MxN term Fourier-cosine represenation
% of "u"

figure('color','w'); 

subplot(2,2,1);
surf(x_J,y_I,u_IJR(:,:,fr)*T_star,'LineStyle','none')
view([0,90])
axis('square')
zlim([0,300*T_star])
title('Measured Temp data')

subplot(2,2,2);
surf(x_J,y_I,fu_IJR(:,:,fr)*T_star,'LineStyle','none')
view([0,90])
axis('square')
zlim([0,300*T_star])
title('Fourier-cosine representation')

subplot(2,2,[3:4])
MaxValue = max(max(u_IJR(:,:,fr)));
surf(x_J,y_I,100*((fu_IJR(:,:,fr)-u_IJR(:,:,fr))./MaxValue),'LineStyle','none')
view([0,90])
colorbar
zlim([1e-4,50])
caxis([0,50])
axis('square')
title('% Error in approximation')

%% Assemble U_NM 
%==========================================================================

% Assemble U_NM in a manner that is suitable for the next steps:
for r = 1:R
    for n = 1:N
        for m = 1:M
            U_R{n,m}(r,1) = U_NM{r}(n,m);
        end
    end
end

%% Calculate Toeplitz matrix
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

%% Conjugate gradient:
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

%% Compare temperature fourier component with reconstructed temperature:
%==========================================================================

% Q_R is the fourier coeficient as function of time for all "M" and "N" modes
% Plot Q_R{1,1} and compare it to d/dt(U_R{1,1})

% Choose M and N mode number to compare:
m = 1;
n = 1;

% Inferred fourier component of heat flux:
% =========================================================================
figure('color','w'); 
hold on
hQ(1) = plot(Q_R{n,m},'r.-','Markersize',8);
hQ(2) = plot((diff(U_R{n,m})./diff(t_R)),'g.-');
legend(hQ,'Q_R','\propto dU_R/dt')
titleText = ['Q_R{',num2str(n),',',num2str(m),'}'];
title(titleText)
box on

% Compare the measured T (fourier component) with the reconstructed one:
% =========================================================================
figure('color','w'); 

% Font sizes:
fontSize.axes = 11;
fontSize.title = 13;
fontSize.legend = 11;
fontSize.label = 12;
fontSize.colorbar = 12;

% First comparison:
% -------------------------------------------------------------------------
% Choose M and N mode number to compare:
m = [1,1];
n = [1,3];

try
    clear legendText
end

for ii = 1:numel(m)
    subplot(2,1,ii)
    hold on

    % Plot data:
    hU(1) = plot(t_dT,UdInfer_R{n(ii),m(ii)},'k','LineWidth',2);
    legendText{1} = ['Reconstructed ','$U_{',num2str(n(ii)),',',num2str(m(ii)),'}$'];
    hU(2) = plot(t_dT,U_R{n(ii),m(ii)},'rsq-','Markersize',5);
    legendText{2} = ['$U_{',num2str(n(ii)),',',num2str(m(ii)),'}$'];

    % Labels:
    ylabel('$\Delta$T [K]','interpreter','Latex','fontSize',fontSize.label)
    xlabel('time [s]','interpreter','Latex','fontSize',fontSize.label)

    % Legend:
    hLeg = legend(hU,legendText);
    set(hLeg,'interpreter','Latex','fontSize',fontSize.legend,'Location','best')

    % Formatting:
%     xlim([0.7,2])
    box on
    set(gca,'fontName','Times','fontSize',fontSize.axes)
end
%% Assemble Fourier Cosine solution for the heat flux:
% =========================================================================

% Arrange Q_R in a manner that can be used in a Fourier-cosine series:
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

% Assemble heat flux function:
%==========================================================================
% Use a Fourier-cosine series to compute the surface heat flux
dx = Lx*(x_J(2)-x_J(1));
dy = Ly*(y_I(2)-y_I(1));
rng1 = 10:120;
for r = 1:R
        f = Phi_IN*Q_NM{r}*Phi_JM';
        q_IJR(:,:,r) = f;
        pwr(r) = sum(sum(f(5:end-5,5:end-5)))*dx*dy*q_star;
end

% ########################################################################
% ########################################################################
% q_IJR is the solution to the entire problem
% This is the time-dependent surface heat flux that produces the measured 2D
% time-dependent surface temperatures measured with IR camera.
% ########################################################################
% ########################################################################

%% Plot results
%==========================================================================

% note:
% Plot 2D distribution of highest heat flux
% plot temporal evolution of total power

zMax = max(max(max(q_IJR)));
zMax = 1;

figure('color','w')

% Font sizes:
fontSize.axes = 11;
fontSize.title = 13;
fontSize.legend = 11;
fontSize.label = 12;
fontSize.colorbar = 12;

[~,fr] = max(pwr);
fr_substract = 80; 
for r = [1:3:length(rngPlasma),40]
    zz = q_IJR(:,:,r) - 0*q_IJR(:,:,fr_substract); 

    % Plot heat flux:
    surf(xx*1e3,yy*1e3,zz,'LineStyle','none')
    
    % Plot LUFS:
    hold on
%     hCirc(2) = plot3(xCirc*10,yCirc*10,ones(size(xCirc))*zMax,'k','LineStyle',':','LineWidth',2);
    
    % Draw plate boundaries:
    line([-45,+45],[+45,+45],[1,1],'color','k','LineWidth',2)
    line([-45,+45],[-45,-45],[1,1],'color','k','LineWidth',2)
    line([+45,+45],[-45,+45],[1,1],'color','k','LineWidth',2)
    line([-45,-45],[-45,+45],[1,1],'color','k','LineWidth',2)
        
    % Labels:
    set(gca,'XTick',[-50:10:50])
    title(['shot: ',num2str(shot),' , frame:',num2str(r),', [MWm$^{-2}$]'],'interpreter','latex','fontSize',fontSize.title)
    
    % Formatting:
    colormap(flipud(hot))
    caxis([-0,zMax])
    zlim([-0,zMax])
    view([-45,45])
    view([00,90])
%     axis image
    colorbar
    hold off    
    set(gca,'fontName','Times','fontSize',fontSize.axes)

    drawnow
    pause(0.01)
end

% Save figure:
figureName = ['shot_',num2str(shot),'_HeatFluxMap'];
saveas(gcf, figureName,'tiff')

figure('color','w')
plot(pwr*1e-3,'k','LineWidth',2)

% Labels:
ylabel('[kW]','interpreter','latex','fontSize',fontSize.label)
title(['shot: ',num2str(shot),' , 100 Hz frame rate'],'interpreter','latex','fontSize',fontSize.title)

% Formatting:
ylim([-0.01,2])
xlim([0,70])
xlabel('frame')
grid on
box on
set(gca,'fontName','Times','fontSize',fontSize.axes)

% Save figure:
figureName = ['shot_',num2str(shot),'_IntegratedPower'];
saveas(gcf, figureName,'tiff')

return

%% Save data:
% =========================================================================

% Data information:
% -------------------------------------------------------------------------
data.shot = shot;
data.thermalParams = d.thermalParams;
data.tempRange = d.tempRange;

% Plate dimensions:
% -------------------------------------------------------------------------
data.Lx = Lx;
data.Ly = Ly;
data.Lz = Lz;

% Plate material properties:
% -------------------------------------------------------------------------
data.material = plateMaterial;
data.rho = rho;
data.kt = kt;
data.cp = cp;
data.a = a;

% Characterstic quantities:
% -------------------------------------------------------------------------
data.q_star = q_star;
data.t_star = t_star;
data.T_star = T_star;

% Initial temperature distribution:
% -------------------------------------------------------------------------
data.T0 = T0;

% Fourier cosine data:
% -------------------------------------------------------------------------
% Fourier modes:
data.FourierCosine.M = M; % "x"
data.FourierCosine.N = N; % "y"
data.FourierCosine.L = L; % "z"
% Normalized temperature fourier coefficients:
data.U_R = U_R;
data.U_NM = U_NM;
% Normalized heat flux fourier coefficients:
data.Q_R = Q_R;
data.Q_NM = Q_NM;

% Normalized variables:
% -------------------------------------------------------------------------
% Normalized Input temperature data:
data.u_IJR = u_IJR;
% Normalized Fourier representation of u_IJR:
data.fu_IJR = fu_IJR;
% Normalized Heat flux:
data.q_IJR = q_IJR(:,:,r);
% Normalized coordinates:
data.t_R = t_R;
data.x_J = linspace(0,1,J)';
data.y_I = linspace(0,1,I)';

% Plate coordinates, physical units:
% -------------------------------------------------------------------------
data.xx = Lx*(x_J - 0.5);
data.yy = Ly*(y_I - 0.5);

% Integrated IR emission:
% -------------------------------------------------------------------------
data.intf = intf;

% Save data:
% -------------------------------------------------------------------------
fileName = ['Step_B_heatFlux2D_',num2str(shot),'.mat'];
save(fileName,'data');
