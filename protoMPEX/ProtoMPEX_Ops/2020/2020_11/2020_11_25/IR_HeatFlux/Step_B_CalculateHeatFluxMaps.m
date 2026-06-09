% Step_B:
% Load temperature data
% Calculate heat flux

clear all
close all
clc

saveData = 1;

% Load data:
% =========================================================================
dataToLoad = 'Step_A_temp2D_30698.mat';
load(dataToLoad);

% Create variables:
% =========================================================================
t_dT = d.t_dT2D;
dT = d.dT2D;
intf = d.intf;
shot = d.shot;

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
M = 70;
% "y":
N = 50;
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
Rp = 2.06;

% Create image coordinates:
xx = Lx*(x_J - 0.5);
yy = Ly*(y_I - 0.5);

% Choose the brigthest frame:
[~,fr] = max(intf);

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
xCirc = Rp*cos(linspace(0,2*pi))+0.2;
yCirc = Rp*sin(linspace(0,2*pi))-0.2;
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

% Save figure:
saveas(gcf,'Step_B_dT_2D','tiffn')

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

% Save figure:
saveas(gcf,'Step_B_CosineFourierContent','tiffn')

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

% Save figure:
saveas(gcf,'Step_B_ErrorFourierCosineRepresentation','tiffn')

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
    xlim([0.7,2])
    box on
    set(gca,'fontName','Times','fontSize',fontSize.axes)
end

% Save figure:
saveas(gcf,'Step_B_TemperatureReconstruction','tiffn')

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
        pwr(r) = sum(sum(f(:,:)))*dx*dy*q_star;
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
% zMax = 2;

figure('color','w')

% Font sizes:
fontSize.axes = 11;
fontSize.title = 13;
fontSize.legend = 11;
fontSize.label = 12;
fontSize.colorbar = 12;

[~,fr] = max(pwr);
for r = [90:1:150,fr]
    zz = q_IJR(:,:,r); 

    % Plot heat flux:
    surf(xx*1e3,yy*1e3,zz,'LineStyle','none')
    
    % Plot LUFS:
    hold on
    hCirc(2) = plot3(xCirc*10,yCirc*10,ones(size(xCirc))*zMax,'k','LineStyle',':','LineWidth',2);
    
    % Draw plate boundaries:
    line([-45,+45],[+45,+45],[1,1],'color','k','LineWidth',2)
    line([-45,+45],[-45,-45],[1,1],'color','k','LineWidth',2)
    line([+45,+45],[-45,+45],[1,1],'color','k','LineWidth',2)
    line([-45,-45],[-45,+45],[1,1],'color','k','LineWidth',2)
        
    % Labels:
    set(gca,'XTick',[-50:10:50])
    title(['shot: ',num2str(shot),' , frame:',num2str(r)],'interpreter','latex','fontSize',fontSize.title)
    
    % Formatting:
    colormap(flipud(hot))
    caxis([-0,zMax])
    zlim([-0,zMax])
    view([-45,45])
    view([00,90])
    axis image
    colorbar
    hold off    
    set(gca,'fontName','Times','fontSize',fontSize.axes)

    drawnow
    pause(0.01)
end

% Save figure:
saveas(gcf,'Step_B_HeatFlux_2DProfile','tiffn')

figure('color','w')
plot(t_dT,pwr*1e-3,'k','LineWidth',2)

% Labels:
ylabel('[kW]','interpreter','latex','fontSize',fontSize.label)
title(['shot: ',num2str(shot)],'interpreter','latex','fontSize',fontSize.title)

% Formatting:
ylim([-0.05,18])
xlim([0.5,1.9])
grid on
box on
set(gca,'fontName','Times','fontSize',fontSize.axes)

% Save figure:
saveas(gcf,'Step_B_TotalHeatTarget','tiffn')

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
