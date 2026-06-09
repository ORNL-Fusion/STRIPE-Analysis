% Step_1b:
% Load temperature data
% Calculate heat flux

clear all
close all
clc

% Load data:
inputfileName = dir('*Cropped.mat');
load(inputfileName.name);

% Plasma radius based on flux mapping:
Rp          = [3.448       ,2.876       ,2.422      ];
profileType = {'profile A' ,'profile B' ,'profile C'};

computeFlag = 1;
saveData = 0;
saveFig = 1;

if computeFlag
    
    % Select shots to analyse:
    shotsToAnalyze = [1:3];

    %% Normalize data:
    for kk = shotsToAnalyze
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
        t_R{kk} = t_dT{kk}'/t_star;
        dt_R(kk) = mean(diff(t_R{kk}));

        % Initial temperature distribution:
        T0{kk} = dT{kk}(:,:,1);

        % Define size of data:
        [I,J,R] = size(dT{kk});

        % Calculate normalized temperature
        for r = 1:R
            u_IJR{kk}(:,:,r) = (dT{kk}(:,:,r) - T0{kk})./T_star;
        end    
    end

    try 
        clear dT
    end

    %% Compute 2D Fourier series representation:
    % =========================================================================

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

    for kk = shotsToAnalyze
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
        for r = 1:length(t_R{kk})
            [a1,a2,a3,a4] = FourierCosine_2D_v3(x_J,y_I,u_IJR{kk}(:,:,r),M,N);
            % "fu_IJR" is the Fourier-cosine approximation of "u_IJR"  
            fu_IJR{kk}(:,:,r) = a1;
            % "U_NM" are the time dependent Fourier cosine coefficients of "u"
            U_NM{kk}{r} = a2;
        end

        % Assign eigenfunction:
        % -------------------------------------------------------------------------
        % Phi_JM is the eigenfunction matrices for all "x" and all "m"
        Phi_JM{kk} = a3;
        % Phi_IN is the eigenfunction matrices for all "y" and all "n"
        Phi_IN{kk} = a4;
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

        % Create image coordinates:
        xx = Lx*(x_J - 0.5);
        yy = Ly*(y_I - 0.5);

        % Choose the brigthest frame:
        [~,fr] = max(intf_c{kk});

        % Temperature range:
        Tmax = 400;

        % Input temperature:
        % ---------------------------------------------------------------------
        subplot(1,2,1)
        hold on;
        ax(1) = gca;
        rngPlot = find(xx*1e2 > -4.4 & xx*1e2 <+4.4);
        surf(xx(rngPlot)*1e2,yy(rngPlot)*1e2,u_IJR{kk}(rngPlot,rngPlot,fr)*T_star,'LineStyle','none');

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
        legendText = ['Shot: ',num2str(shot(kk))];
        hL = text(-4,4,600,legendText);
        set(hL,'interpreter','Latex','fontSize',fontSize.legend)
        hC(1) = colorbar;
        hT = text(-1.5,8.5,'${\Delta}T$ [K]');
        set(hT,'interpreter','Latex','fontSize',fontSize.label)

        % Add LUFS:
        xCirc{kk} = Rp(kk)*cos(linspace(0,2*pi));
        yCirc{kk} = Rp(kk)*sin(linspace(0,2*pi));
        hCirc(1) = plot3(xCirc{kk},yCirc{kk},ones(size(xCirc{kk}))*Tmax,'k','LineStyle',':','LineWidth',2);

        % Fourier representation:
        % ---------------------------------------------------------------------
        subplot(1,2,2)
        hold on
        ax(2) = gca;
        surf(xx(rngPlot)*1e2,yy(rngPlot)*1e2,fu_IJR{kk}(rngPlot,rngPlot,fr)*T_star,'LineStyle','none')

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
        legendText = ['Shot: ',num2str(shot(kk))];
        hL = text(-4,4,600,legendText);
        set(hL,'interpreter','Latex','fontSize',fontSize.legend)
        hC(2) = colorbar;
        colormap((flipud(hot)));
        hT = text(-1.5,8.5,'${\Delta}T$ [K]');
        set(hT,'interpreter','Latex','fontSize',fontSize.label)

        % Add LUFS:
        hCirc(2) = plot3(xCirc{kk},yCirc{kk},ones(size(xCirc{kk}))*Tmax,'k','LineStyle',':','LineWidth',2);

        % Global formatting:
        set(ax,'fontName','Times','fontSize',fontSize.axes);
        set(hC,'fontSize',fontSize.colorbar,'Location','northoutside')
        xlim(ax,0.5*Lx*[-1,1]*1e2)
        ylim(ax,0.5*Ly*[-1,1]*1e2)

        % Save figure:
        figureName = ['Step_2c_FourierRepresentation_',num2str(round(PS2_current(kk))),'_kW_',num2str(shot(kk))];
        saveas(gcf,figureName,'tiffn')

        % Visualize fourier components:
        % =========================================================================
        figure('color','w')
        fontSize.axes = 12;
        fontSize.title = 13;
        fontSize.legend = 11;
        fontSize.label = 12;
        fontSize.colorbar = 12;

        % Fourier spectrum:
        Umn_abs = abs(U_NM{kk}{fr});
        Umn_abs_max = max(max(Umn_abs));
        hFC = pcolor(m_M,n_N,Umn_abs/Umn_abs_max);
        xlim([0,20])
        ylim([0,20])
        view([0,90])
        axis('square')
        title(['Shot ',num2str(shot(kk)),', Fourier-cosine spectrum'],'interpreter','latex','fontSize',fontSize.title)
        xlabel('m','interpreter','latex','fontSize',fontSize.label)
        ylabel('n','interpreter','latex','fontSize',fontSize.label)
        colormap((flipud(hot)));
        hFCC = colorbar;
        box on

        % Formatting:=
        set(gca,'fontName','Times','fontSize',fontSize.axes);
        set(hFCC,'fontSize',fontSize.colorbar,'Location','Eastoutside')
        hFC.LineStyle = ':';

        if kk == 10
            % Save figure:
            figureName = ['Step_2c_FourierContent_',num2str(shot(kk))];
            saveas(gcf,figureName,'tiffn')
        end


    end

    %% Test Fourier representation:
    % =========================================================================
    % Compare the measured data "u"  with the MxN term Fourier-cosine
    % represenation of "u"

    % Input:
    % y_I
    % x_J
    % u_IJR
    % fu_IJR
    % Tstar

    kk = 1;

    figure('color','w'); 

    subplot(2,2,1);
    surf(x_J,y_I,u_IJR{kk}(:,:,fr)*T_star,'LineStyle','none')
    view([0,90])
    axis('square')
    zlim([0,300*T_star])
    title('Measured Temp data')

    subplot(2,2,2);
    surf(x_J,y_I,fu_IJR{kk}(:,:,fr)*T_star,'LineStyle','none')
    view([0,90])
    axis('square')
    zlim([0,300*T_star])
    title('Fourier-cosine representation')

    subplot(2,2,[3:4])
    MaxValue = max(max(u_IJR{kk}(:,:,fr)));
    surf(x_J,y_I,100*((fu_IJR{kk}(:,:,fr)-u_IJR{kk}(:,:,fr))./MaxValue),'LineStyle','none')
    view([0,90])
    colorbar
    zlim([1e-4,50])
    caxis([0,50])
    axis('square')
    title('% Error in approximation')


    %% Extract heat flux:
    %==========================================================================

    for kk = shotsToAnalyze
        % Assemble U_NM in a manner that is suitable for the next steps:
        for r = 1:R
            for n = 1:N
                for m = 1:M
                    U_R{kk}{n,m}(r,1) = U_NM{kk}{r}(n,m);
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

        % For Back-side  imaging use z_K = 1;
        % For Front-side imaging use z_K = 0;
        z_K = 1;
        % Calculate Toeplitz matrix for a single value of z_K
        [P_RR] = Toeplitz_NM(m_M,n_N,l_L,Lx,Ly,Lz,z_K,t_R{kk});


        % Conjugate gradient:

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
                U_data = U_R{kk}{n,m};

                if sqrt(m*m + n*n) < 120
                    % Initialize "n,m" heat flux fourier-cosine coeffiecient guess:
                    Q_iteration{1} = (diff(U_R{kk}{n,m})./dt_R(kk));
                    Q_iteration{1} = [Q_iteration{1};Q_iteration{1}(end)];

                    % Inverse method:
                    [Q_iteration,U_dataInferred,Residual] = ...
                    IHCP_ConjugateGradient(P_RR{n,m},Q_iteration,U_data,dtt,Ni_a,3,[]);
                else
                    Q_iteration{1} = zeros(size(t_R{kk}));
                    [Q_iteration,U_dataInferred,Residual] = ...
                    IHCP_ConjugateGradient(P_RR{n,m},Q_iteration,U_data,dtt,Ni_b,3,[]);
                end

                % Assign last iteration as solution:
                Q_R{kk}{n,m} = Q_iteration{end};
                UdInfer_R{kk}{n,m} = U_dataInferred{end};
                Res{n,m} = Residual;

                % Clear dummy variables before proceeding to the next "m,n" coeficient
                clear Q_iteration U_dataInferred Residual
            end
            disp(['m: ',num2str(m)])
        end

        % At this point, we have calculated all the "m,n" Fourier-cosine
        % coefficients of the surface heat flux
        toc
    end
    %% Compare temperature Fourier component with reconstructed temperature:
    %==========================================================================

    % Q_R is the fourier coeficient as function of time for all "M" and "N" modes
    % Plot Q_R{1,1} and compare it to d/dt(U_R{1,1})

    try
        kk = 10;
        Q_R{kk}{n,m};
    catch
        kk = 1;
    end
            
    % Choose M and N mode number to compare:
    m = 1;
    n = 1;

    % Inferred fourier component of heat flux:
    % =========================================================================
    figure('color','w'); 
    hold on
    hQ(1) = plot(Q_R{kk}{n,m},'r.-','Markersize',8);
    hQ(2) = plot((diff(U_R{kk}{n,m})./diff(t_R{kk})),'g.-');
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
        hU(1) = plot(t_dT{kk},UdInfer_R{kk}{n(ii),m(ii)},'k','LineWidth',2);
        legendText{1} = ['Reconstructed ','$U_{',num2str(n(ii)),',',num2str(m(ii)),'}$'];
        hU(2) = plot(t_dT{kk},U_R{kk}{n(ii),m(ii)},'rsq-','Markersize',5);
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

    for kk = shotsToAnalyze
        % Arrange Q_R in a manner that can be used in a Fourier-cosine series:
        for m = 1:M
            for n = 1:N
                    for r = 1:R
                        if sqrt(m*m + n*n) > 120 && sqrt(m*m + n*n) < 130
                            Q_NM{r}(n,m) = 0.01*Q_R{kk}{n,m}(r);
                        else
                            Q_NM{r}(n,m) = Q_R{kk}{n,m}(r);
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
                f = Phi_IN{kk}*Q_NM{r}*Phi_JM{kk}';
                q_IJR{kk}(:,:,r) = f;
                pwr{kk}(r) = sum(sum(f(5:end-5,5:end-5)))*dx*dy*q_star;
        end

        % ########################################################################
        % ########################################################################
        % q_IJR is the solution to the entire problem
        % This is the time-dependent surface heat flux that produces the measured 2D
        % time-dependent surface temperatures measured with IR camera.
        % ########################################################################
        % ########################################################################
    end
    %% Plot results
    %==========================================================================
    close all

    % Font sizes:
    fontSize.axes = 16;
    fontSize.title = 16;
    fontSize.legend = 11;
    fontSize.label = 16;
    fontSize.colorbar = 13;
    fontSize.profileName= 19;
    
    % Select highest heat flux value:
    zMax = 0.8*max(max(max(q_IJR{1})));
    zMin = 1;
    zMin_caxis = -0.5;
    
    % Select frame to substract:
    fr_substract = 85; 
    substractFlag = 1;

    % Time averaging
    timeAverage = 1;
    
    % Surface of contour:
    useSurf    =  0;
    useContour1 = 1;
    useContour2 = 0;
    
    % Heat flux 2D map:
    % -------------------------------------------------------------------------
    for kk = shotsToAnalyze
        figure('color','w')

        % Select hottest frame:
        [~,fr] = max(pwr{kk})

        if timeAverage
            ff = movmean(q_IJR{kk},3,1);
        else
            ff = q_IJR{kk};
        end
        
        for r = [fr]%,20:5:length(intf_c{kk})-30,fr]

            % Create data:
            zz = ff(:,:,r) - substractFlag*ff(:,:,fr_substract);

            % Plot heat flux:
            if useSurf == 1
                surf(xx*1e3,yy*1e3,zz,'LineStyle','none')
            elseif useContour1 == 1
                zLevels = linspace(zMin,zMax,15);
                [cCon{kk},hCon{kk}] = contourf(xx*1e3,yy*1e3,zz,zLevels);
                set(hCon{kk},'LineStyle','none');
            end

            hold on
            
            % Contours lines:
            if useContour2 == 1
                zLevels = [1,1];
                [cCon2{kk},hCon2{kk}] = contour(xx*1e3,yy*1e3,zz,zLevels);
                set(hCon2{kk},'LineStyle','-','lineColor','k');
            end
 
            % Plot LUFS:
            hCirc(2) = plot3(xCirc{kk}*10 - 3 ,yCirc{kk}*10,ones(size(xCirc{kk}))*zMax,'k','LineStyle',':','LineWidth',3);

            % Draw plate boundaries:
            line([-45,+45],[+45,+45],[1,1],'color','k','LineWidth',2)
            line([-45,+45],[-45,-45],[1,1],'color','k','LineWidth',2)
            line([+45,+45],[-45,+45],[1,1],'color','k','LineWidth',2)
            line([-45,-45],[-45,+45],[1,1],'color','k','LineWidth',2)

            % Labels:
            set(gca,'XTick',[-50:10:50])
            title(['shot: ',num2str(shot(kk)),' $q(x_*,y_*,t)$ [MWm$^{-2}$]'],'interpreter','latex','fontSize',fontSize.title)
            ylabel('$y_*$ [mm]','interpreter','latex','fontSize',fontSize.label)
            xlabel('$x_*$ [mm]','interpreter','latex','fontSize',fontSize.label)
            
            % Profile type:
            hText = text(-40,40,10,profileType{kk});
            set(hText,'interpreter','latex','fontSize',fontSize.profileName)

            % Formatting:
            colormap(flipud(hot))
            caxis([zMin_caxis,zMax])
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
        
        % Arrows indicating 28 GHz injection angle:
        % -----------------------------------------------------------------
        % Text arrows fields:
        fields.HorizontalAlignment = 'center';
        fields.Color = 'g';
        fields.FontSize = 11;
        fields.Interpreter = 'Latex';
        fields.HeadLength = 10;
        fields.HeadWidth = 10;
        fields.HeadStyle = 'vback2';
        fields.LineWidth = 3;

        fields.String = '';
        switch kk
            case 1
                x =  [38  ,29];
            case 2
                x =  [38  ,25];                
            case 3
                x =  [38  ,21];
        end
        y =      [6   ,6 ];
        hta = myTextArrow(gca,x,y,fields);
        
        y =      [3   ,3 ];
        hta = myTextArrow(gca,x,y,fields);
        
        y =      [0   ,0 ];
        hta = myTextArrow(gca,x,y,fields);
        
        y =      [-3  ,-3 ];
        hta = myTextArrow(gca,x,y,fields);
       
        y =      [-6  ,-6 ];
        hta = myTextArrow(gca,x,y,fields);
                
        % Save figure:
        % =========================================================================
        if saveFig
            figureName = ['step_2c_HeatFluxMap_',num2str(round(PS2_current(kk))),'_A_shot_',num2str(shot(kk))];

            % PDF figure:
            exportgraphics(gcf,[figureName,'.pdf'],'Resolution',600,'ContentType', 'vector') 

            % TIFF figure:
            exportgraphics(gcf,[figureName,'.tiff'],'Resolution',600) 
        end
        
        if 0
            % Temporal evolution of heating rate:
            % -------------------------------------------------------------------------
            figure('color','w')
            plot(pwr{kk}*1e-3,'k','LineWidth',2)

            % Labels:
            ylabel('[kW]','interpreter','latex','fontSize',fontSize.label)
            title(['shot: ',num2str(shot(kk)),' , 100 Hz frame rate'],'interpreter','latex','fontSize',fontSize.title)

            % Formatting:
            ylim([-0.01,20])
            % xlim([0,70])
            xlabel('frame')
            grid on
            box on
            set(gca,'fontName','Times','fontSize',fontSize.axes)
            
            % Save figure:
            % =========================================================================
            if saveFig
                figureName = ['step_2c_IntegratedPower_',num2str(round(PS2_current(kk))),'_A_shot_',num2str(shot(kk))];

                % PDF figure:
                exportgraphics(gcf,[figureName,'.pdf'],'Resolution',600,'ContentType', 'vector') 

                % TIFF figure:
                exportgraphics(gcf,[figureName,'.tiff'],'Resolution',600) 
            end               
        end     
        
    end

    % Save data:
    % =====================================================================
    if saveData
        fileName = ['Step_2c_heatFlux2D.mat'];
        save(fileName)
    end

else
    fileName = ['Step_2c_heatFlux2D.mat'];
    load(fileName)
end