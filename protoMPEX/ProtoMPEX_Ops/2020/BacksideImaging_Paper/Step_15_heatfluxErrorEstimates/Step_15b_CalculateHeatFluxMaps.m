% Step_15b:
% Load temperature data
% Calculate heat flux

clear all
close all
clc

computeFlag = 0;

saveData_1 = 0;
saveData_2 = 0;

% Load data:
inputfileName = dir('*plateTemperature.mat');
load(inputfileName.name);

if computeFlag
    
% Select shots to analyse:
shotsToAnalyze = [1];

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
    
% Normalize data:
% =========================================================================
for kk = shotsToAnalyze       
    for ee = 1:numel(plate{kk})

        % Extract data:
        t_dT{kk}{ee} = plate{kk}{ee}.t;
        dT{kk}{ee}   = plate{kk}{ee}.dT;

        % Calculate normalized time:
        t_R{kk}{ee} = t_dT{kk}{ee}'/t_star;
        dt_R(kk,ee) = mean(diff(t_R{kk}{ee}));

        % Initial temperature distribution:
        T0{kk}{ee} = dT{kk}{ee}(:,:,1);

        % Define size of data:
        [I,J,R] = size(dT{kk}{ee});

        % Calculate normalized temperature
        for r = 1:R
            u_IJR{kk}{ee}(:,:,r) = (dT{kk}{ee}(:,:,r) - T0{kk}{ee})./T_star;
        end

    end
end

try 
    clear dT
end

% Compute 2D Fourier series representation:
% =========================================================================
for kk = shotsToAnalyze
    % Input:
    % I
    % J
    % t_R
    % u_IJR
    % Tstar

    % Output:
    % fu_IJR
    % U_MN
    % Phi_JM
    % Phi_IN

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
    for ee = 1:numel(plate{kk})

        tic
        for r = 1:length(t_R{kk}{ee})
            [a1,a2,a3,a4] = FourierCosine_2D_v3(x_J,y_I,u_IJR{kk}{ee}(:,:,r),M,N);
            % "fu_IJR" is the Fourier-cosine approximation of "u_IJR"  
            fu_IJR{kk}{ee}(:,:,r) = a1;
            % "U_NM" are the time dependent Fourier cosine coefficients of "u"
            U_NM{kk}{ee}{r} = a2;
        end

        % Assign eigenfunction:
        % -------------------------------------------------------------------------
        % Phi_JM is the eigenfunction matrices for all "x" and all "m"
        Phi_JM{kk}{ee} = a3;
        % Phi_IN is the eigenfunction matrices for all "y" and all "n"
        Phi_IN{kk}{ee} = a4;
        toc

    end
end

% Calculate heat flux:
%==========================================================================
for kk = shotsToAnalyze
    
    for ee = 1:numel(plate{kk})

        % Assemble U_NM in a manner that is suitable for the next steps:
        for r = 1:R
            for n = 1:N
                for m = 1:M
                    U_R{kk}{ee}{n,m}(r,1) = U_NM{kk}{ee}{r}(n,m);
                end
            end
        end

        % Calculate Toeplitz matrix
        % The purpose of the Toeplitz matrix is to perform convolution on a vector.
        % Its the finite-dimensional representation of the convolution integral
        % operator

        % P_RR{N,M,K}
        % This cell has NxMxK matrices.
        % This consumes the largest memory in the calculation
        % By evaluating only one "z_K" location we can reduce the size to NxM
        % P_RR_new = P_RR{:,:,k};

        % Calculate Toeplitz matrix for a single value of z_K:
        % ====================================================
        % For Back-side  imaging use z_K = 1;
        % For Front-side imaging use z_K = 0;
        z_K = 1;
        [P_RR] = Toeplitz_NM(m_M,n_N,l_L,Lx,Ly,Lz,z_K,t_R{kk}{ee});


        % Conjugate gradient:
        % ===================
        % Conjugate gradient method to find surface heat flux that produced the
        % surface temperature distribution that was measured with IR camera

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
                U_data = U_R{kk}{ee}{n,m};

                if sqrt(m*m + n*n) < 120
                    % Initialize "n,m" heat flux fourier-cosine coeffiecient guess:
                    Q_iteration{1} = (diff(U_R{kk}{ee}{n,m})./dt_R(kk,ee));
                    Q_iteration{1} = [Q_iteration{1};Q_iteration{1}(end)];

                    % Inverse method:
                    [Q_iteration,U_dataInferred,Residual] = ...
                    IHCP_ConjugateGradient(P_RR{n,m},Q_iteration,U_data,dtt,Ni_a,3,[]);
                else
                    Q_iteration{1} = zeros(size(t_R{kk}{ee}));
                    [Q_iteration,U_dataInferred,Residual] = ...
                    IHCP_ConjugateGradient(P_RR{n,m},Q_iteration,U_data,dtt,Ni_b,3,[]);
                end

                % Assign last iteration as solution:
                Q_R{kk}{ee}{n,m} = Q_iteration{end};
                UdInfer_R{kk}{ee}{n,m} = U_dataInferred{end};
                Res{ee}{n,m} = Residual;

                % Clear dummy variables before proceeding to the next "m,n" coeficient
                clear Q_iteration U_dataInferred Residual
            end
            disp(['m: ',num2str(m)])
        end

        % At this point, we have calculated all the "m,n" Fourier-cosine
        % coefficients of the surface heat flux
        toc
    end
end

% Assemble Fourier Cosine solution for the heat flux:
% =========================================================================
for kk = shotsToAnalyze
    for ee = 1:numel(plate{kk})

        % Arrange Q_R in a manner that can be used in a Fourier-cosine series:
        for m = 1:M
            for n = 1:N
                    for r = 1:R
                        if sqrt(m*m + n*n) > 120 && sqrt(m*m + n*n) < 130
                            Q_NM{r}(n,m) = 0.01*Q_R{kk}{ee}{n,m}(r);
                        else
                            Q_NM{r}(n,m) = Q_R{kk}{ee}{n,m}(r);
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
            f = Phi_IN{kk}{ee}*Q_NM{r}*Phi_JM{kk}{ee}';
            q_IJR{kk}{ee}(:,:,r) = f;
            pwr{kk}{ee}(r) = sum(sum(f(5:end-5,5:end-5)))*dx*dy*q_star;
    end

    % ########################################################################
    % ########################################################################
    % q_IJR is the solution to the entire problem
    % This is the time-dependent surface heat flux that produces the measured 2D
    % time-dependent surface temperatures measured with IR camera.
    % ########################################################################
    % ########################################################################
    end
end
    
if saveData_1
    fileName = ['Step_15b_heatFluxVariousEmissivities.mat'];
    save(fileName)
end

else
    fileName = ['Step_15b_heatFluxVariousEmissivities.mat'];
    load(fileName)
end

saveFig    = 1;

%% Plot hottest frame:
% =========================================================================
for kk = shotsToAnalyze
    for ee = 1:numel(plate{kk})

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
        [~,fr] = max(plate{kk}{ee}.integratedIntensity);

        % Temperature range:
        Tmax = 600;

        % Input temperature:
        % ---------------------------------------------------------------------
        subplot(1,2,1)
        hold on;
        ax(1) = gca;
        rngPlot = find(xx*1e2 > -4.4 & xx*1e2 <+4.4);
        surf(xx(rngPlot)*1e2,yy(rngPlot)*1e2,u_IJR{kk}{ee}(rngPlot,rngPlot,fr)*T_star,'LineStyle','none');

        % Formating:
        view([0,90])
        axis('image')
        box on
        caxis([0,Tmax])
        zlim([0,Tmax])

        % Labels:
        title('Measured temperature','interpreter','latex','fontSize',fontSize.title)
        xlabel('x [cm]','interpreter','latex','fontSize',fontSize.label)
        ylabel('y [cm]','interpreter','latex','fontSize',fontSize.label)
        legendText = ['Shot: ',num2str(plate{kk}{ee}.shot),', emiss: ',num2str(plate{kk}{ee}.emissivity)];
        hL = text(-4,4,600,legendText);
        set(hL,'interpreter','Latex','fontSize',fontSize.legend)
        hC(1) = colorbar;
        hT = text(-1.5,8.5,'${\Delta}T$ [K]');
        set(hT,'interpreter','Latex','fontSize',fontSize.label)

        % Add LUFS:
        xCirc = Rp*cos(linspace(0,2*pi));
        yCirc = Rp*sin(linspace(0,2*pi));
        hCirc(1) = plot3(xCirc,yCirc,ones(size(xCirc))*Tmax,'k','LineStyle',':','LineWidth',2);

        % Fourier representation:
        % ---------------------------------------------------------------------
        subplot(1,2,2)
        hold on
        ax(2) = gca;
        surf(xx(rngPlot)*1e2,yy(rngPlot)*1e2,fu_IJR{kk}{ee}(rngPlot,rngPlot,fr)*T_star,'LineStyle','none')

        % Formating:
        view([0,90])
        axis('image')
        box on
        caxis([0,Tmax])
        zlim([0,Tmax])

        % Labels:    
        title('Fourier representation','interpreter','latex','fontSize',fontSize.title)
        xlabel('x [cm]','interpreter','latex','fontSize',fontSize.label)
        ylabel('y [cm]','interpreter','latex','fontSize',fontSize.label)
        legendText = ['Shot: ',num2str(plate{kk}{ee}.shot),', emiss: ',num2str(plate{kk}{ee}.emissivity)];
        hL = text(-4,4,600,legendText);
        set(hL,'interpreter','Latex','fontSize',fontSize.legend)
        hC(2) = colorbar;
        colormap((flipud(hot)));
        hT = text(-1.5,8.5,'${\Delta}T$ [K]');
        set(hT,'interpreter','Latex','fontSize',fontSize.label)

        % Add LUFS:
        hCirc(2) = plot3(xCirc,yCirc,ones(size(xCirc))*Tmax,'k','LineStyle',':','LineWidth',2);

        % Global formatting:
        set(ax,'fontName','Times','fontSize',fontSize.axes);
        set(hC,'fontSize',fontSize.colorbar,'Location','northoutside')
        xlim(ax,0.5*Lx*[-1,1]*1e2)
        ylim(ax,0.5*Ly*[-1,1]*1e2)

        % Save figure:
        if saveFig
            figureName = ['Step_15b_FourierRepresentation_',num2str(round(plate{kk}{ee}.pwr28GHz)),'_kW_',num2str(plate{kk}{ee}.shot)];
            saveas(gcf,figureName,'tiffn')
        end

    end
end   


%% Plot time evolution of power incident on target:
% =========================================================================
kk = 1;

try
    clear legendText
end

% Font sizes:
fontSize.axes = 11;
fontSize.title = 13;
fontSize.legend = 11;
fontSize.label = 13;
fontSize.colorbar = 12;

figure('color','w')
set(gcf,'Position',[438   305   482   313])
hold on
box on

hA(1) = area(t_dT{kk}{1} + 3.9,pwr{kk}{1}*1e-3);
hA(2) = area(t_dT{kk}{3} + 3.9,pwr{kk}{3}*1e-3);
hL(1) = plot(t_dT{kk}{2} + 3.9,pwr{kk}{2}*1e-3,'k','LineWidth',2);

% Formatting:
set(hA(1),'FaceColor',[1 1 1]*0.8,'EdgeColor',[0 0 0])
set(hA(2),'FaceColor',[1 1 1],'EdgeColor',[0 0 0])
set(gca,'fontName','Times','fontSize',fontSize.axes)
xlim([4.0,4.9])
ylim([0,16])

% Labels:
title(['shot: ',num2str(plate{kk}{ee}.shot)],'interpreter','latex','fontSize',fontSize.title)
ylabel('Power to Target [kW]','interpreter','latex','fontSize',fontSize.label)
xlabel('time [s]','interpreter','latex','fontSize',fontSize.label)
hLeg = legend([hL(1),hA(1)],'$\epsilon$ = 0.3','$\epsilon$ = 0.3 $\pm$ 0.03');

% Formatting:
set(hLeg,'interpreter','latex','fontSize',fontSize.legend,'Location','northwest')

% Save figure:
if saveFig
    figureName = 'step_15b_Uncertainty_power';

    % PDF figure:
    exportgraphics(gcf,[figureName,'.pdf'],'Resolution',300) 

    % TIFF figure:
    exportgraphics(gcf,[figureName,'.tiff'],'Resolution',600) 
end

%% Plot error bars vrs microwave power:
% =========================================================================

rng_pwr = find(pwr{1}{2} > 0.5e3);
err_pwr = (pwr{1}{1} - pwr{1}{3})./pwr{1}{2};

figure('color','w');
hold on
box on
hE(1) = plot(t_dT{kk}{2}          + 3.9,pwr{1}{2}/max(pwr{1}{2})        ,'k','LineWidth',2);
hE(2) = plot(t_dT{kk}{2}(rng_pwr) + 3.9,err_pwr(rng_pwr),'r','LineWidth',2);

hLeg = legend(hE(2),'error');
set(hLeg,'interpreter','latex','fontSize',fontSize.legend,'Location','northwest')

set(gca,'fontName','Times','fontSize',fontSize.axes)
xlim([4.0,4.9])
ylim([0,1])

title(['shot: ',num2str(plate{kk}{ee}.shot),', mean error = ',num2str(mean(err_pwr(rng_pwr)))],'interpreter','latex','fontSize',fontSize.title)
ylabel('Error','interpreter','latex','fontSize',fontSize.label)
xlabel('time [s]','interpreter','latex','fontSize',fontSize.label)

% The data indicates that the mean error is about 17% accross the power
% range: error = P/dP = 0.17 or P = Pmean +- error*Pmean/2

% Save figure:
if saveFig
    figureName = 'step_15b_Uncertainty_value';

    % PDF figure:
    exportgraphics(gcf,[figureName,'.pdf'],'Resolution',300) 

    % TIFF figure:
    exportgraphics(gcf,[figureName,'.tiff'],'Resolution',600) 
end

%% Plot time evolution of hot spot heat flux:

% Extract the hot spot heat flux over time:
% Select averaging region:
X = 66;
Y = 113;
dX = 1;
dY = 1;
rngX = X-dX:X+dX;
rngY = Y-dY:Y+dY;

for jj = 1:numel(t_dT{kk}{ee})
    for ee = 1:numel(plate{kk})
        q_peaked{ee}(jj) = mean(q_IJR{1}{ee}(rngX,rngY,jj) - q_IJR{1}{ee}(rngX,rngY,88),'all')*q_star;
    end
end

figure('color','w'); 
set(gcf,'Position',[438   305   482   313])
hold on
box on
% for jj = 1:numel(t_dT{kk}{ee})
%     for ee = 1:numel(plate{kk})
%         plot(t_dT{1}{1} + 3.9,q_peaked{ee}*1e-6,'k')
%     end
% end

hA(1) = area(t_dT{kk}{1} + 3.9,q_peaked{1}*1e-6);
hA(2) = area(t_dT{kk}{3} + 3.9,q_peaked{3}*1e-6);
hL(1) = plot(t_dT{kk}{2} + 3.9,q_peaked{2}*1e-6,'k','LineWidth',2);

% Formatting:
set(hA(1),'FaceColor',[1 1 1]*0.8,'EdgeColor',[0 0 0])
set(hA(2),'FaceColor',[1 1 1],'EdgeColor',[0 0 0])
set(gca,'fontName','Times','fontSize',fontSize.axes)

set(gca,'fontName','Times','fontSize',fontSize.axes)
xlim([4.0,4.9])
ylim([0,25])

title(['shot: ',num2str(plate{kk}{ee}.shot)],'interpreter','latex','fontSize',fontSize.title)
ylabel('Peak heat flux [MWm$^{-2}$]','interpreter','latex','fontSize',fontSize.label)
xlabel('time [s]','interpreter','latex','fontSize',fontSize.label)
hLeg = legend([hL(1),hA(1)],'$\epsilon$ = 0.3','$\epsilon$ = 0.3 $\pm$ 0.03');

% Formatting:
set(hLeg,'interpreter','latex','fontSize',fontSize.legend,'Location','northwest')

% Save figure:
if saveFig
    figureName = 'step_15b_Uncertainty_peakHeatflux';

    % PDF figure:
    exportgraphics(gcf,[figureName,'.pdf'],'Resolution',300) 

    % TIFF figure:
    exportgraphics(gcf,[figureName,'.tiff'],'Resolution',600) 
end