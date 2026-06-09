% Step_2d_PostProcessHeatFlux
clear all
close all
clc

% Load data:
% =========================================================================
fileName = ['Step_2c_heatFlux2D.mat'];
load(fileName)

% Flags:
% =========================================================================
saveFig = 1;
addUncertainty = 1;

% Rename and rescale some variables:
% =========================================================================
pwr_Target = pwr;    % Power coupled to Target in [W]

%% Extract magnetic field data:
Bdata = load('Step_1_Magneticfield_Bscan.mat');

% Extract magnetic field strength at specified location:
rng_z1 = find(Bdata.z1D > 3.5 & Bdata.z1D < 4);
rng_z2 = find(Bdata.z1D >= 4.2,1);
rng_z3 = find(Bdata.z1D > 2 & Bdata.z1D < 2.25);

for jj = 1:numel(Bdata.B2D)
    B_transport(jj) = max(Bdata.B2D{jj}(rng_z1,1));
    B_target(jj)    = Bdata.B2D{jj}(rng_z2,1);
    B_limiter(jj)   = min(Bdata.B2D{jj}(rng_z1,1));
end

%% Extract power during ECH
% =========================================================================
for ii = 1:numel(pwr_Target)
    % ECH power value:
    [pwr_ech_peak(ii),i_pwr_ech_peak(ii)] = max(pwr_Target{ii});
    di = 1;
    rng = [-di,+di] + i_pwr_ech_peak(ii);
    d_peak_ech(ii) = std(pwr_Target{ii}(rng),1);
    
    % Helicon power prior to ECH:
    rng = [36:62];
    pwr_helicon_mean(ii) = mean(pwr_Target{ii}(rng));
    
    % Error in measurements due emissivity uncertainty:
    d_pwr_helicon_mean(ii) = 0.17*pwr_helicon_mean(ii)/2;
    d_pwr_ech_peak(ii)     = 0.17*pwr_ech_peak(ii)/2;
end

%% Extract 2D heat flux profiles:

% Select frame to substract:
fr_substract = 85; 
substractFlag = 1;

% Time averaging
timeAverage = 1;

for kk = 1:numel(q_IJR)
    
     % Select hottest frame:
    [~,fr28] = max(pwr{kk});

    if timeAverage
        zz = movmean(q_IJR{kk},3,1);
    else
        zz = q_IJR{kk};
    end
    
    heatFlux{kk} = zz(:,:,fr28)- substractFlag*zz(:,:,fr_substract);
end

%% Clean up region outside plasma:

% Threshold radius:
rThreshold = [3.6,3.2,2.8]*1e-2; % [m]
    
for kk = 1:numel(q_IJR)

    % Initialize data:
    heatFlux_clean{kk} = heatFlux{kk};
    
    % Produce data:
    for ii = 1:size(heatFlux{kk},1)
        for jj = 1:size(heatFlux{kk},2)
            % Create radius:
            rr = sqrt(xx(jj)^2 + yy(ii)^2);
            
            if (heatFlux_clean{kk}(ii,jj) > 1 && rr >= rThreshold(kk))
                heatFlux_clean{kk}(ii,jj) = heatFlux_clean{kk}(ii,jj)*exp(-0.5*rr/rThreshold(kk));
            end
        end
    end
end

for kk = 1:3
    figure
    mesh(xx,yy,heatFlux_clean{kk})
end
%% Extract 50% power contour lines:
close all

% Hot spot:
pwr_fraction = [0.2,0.32,0.45];
pwr_fraction = [0.2,0.28,0.45]*0.9;

% 50% power:
% pwr_fraction = [0.17,0.19,0.22];

for kk = 1:numel(q_IJR)

    % Define threshold heat flux:
    qThreshold(kk) = max(max(heatFlux{kk}))*pwr_fraction(kk);

    % Define spatial increaments:
    dxx = xx(2) - xx(1);
    dyy = yy(2) - yy(1);

    % Initialize data:
    heatFlux_peaked{kk} = heatFlux_clean{kk};
    heatFlux_halo{kk}   = heatFlux_clean{kk};
    
    % Produce data:
    for ii = 1:size(heatFlux{kk},1)
        for jj = 1:size(heatFlux{kk},2)
            if heatFlux_peaked{kk}(ii,jj) < qThreshold(kk)
                heatFlux_peaked{kk}(ii,jj) = 0;
            else
                heatFlux_halo{kk}(ii,jj) = 0;
            end
            
            if heatFlux_halo{kk}(ii,jj) < 0.1
                heatFlux_halo{kk}(ii,jj) = 0;
            end
        end
    end
    
    % Smooth data:
    frame_length = 10;
    heatFlux_peaked{kk} = movmean(heatFlux_peaked{kk},frame_length,1);
    heatFlux_peaked{kk} = movmean(heatFlux_peaked{kk},frame_length,2);
    
    % Contour line:
    C = contour(xx,yy,heatFlux_peaked{kk},[1,1]*qThreshold(kk));
    xx_peaked{kk} = C(1,2:end);
    yy_peaked{kk} = C(2,2:end);
    
    % Remove large numners:
    rng_rmv = find(abs(xx_peaked{kk}) > 1);
    xx_peaked{kk}(rng_rmv) = [];
    yy_peaked{kk}(rng_rmv) = [];
    
    rng_rmv = find(abs(yy_peaked{kk}) > 1);
    xx_peaked{kk}(rng_rmv) = [];
    yy_peaked{kk}(rng_rmv) = [];
    
    % Intergrate power [W]:
    pwr_peaked(kk) = sum(sum(heatFlux_peaked{kk}*1e6))*dxx*dxx;
    pwr_halo(kk)   = sum(sum(heatFlux_halo{kk}*1e6))*dxx*dxx;

end

if 1
    for kk = 1:numel(q_IJR)
            % Plot data:
            figure; 
            subplot(1,2,1)
            mesh(xx*1e3,yy*1e3,heatFlux_peaked{kk})
            zlim([0,22])
            view([0,90])
            axis square

            subplot(1,2,2)
            mesh(xx*1e3,yy*1e3,heatFlux_halo{kk})
            zlim([0,22])
            view([0,90])
            axis square
    end
end


figure;
hold on
hPwr(1) = plot(B_transport,pwr_ech_peak*1e-3,'ksq');
hPwr(2) = plot(B_transport,pwr_peaked*1e-3,'ro');
hPwr(3) = plot(B_transport,pwr_halo*1e-3,'go');
plot(B_transport,(pwr_halo + pwr_peaked)*1e-3,'c');
 
legend(hPwr,'Total','Peaked','Halo')
ylim([0,15])
xlim([0,1.2])

%% Plot: 2D heat flux

% Font sizes:
fontSize.axes = 16;
fontSize.title = 16;
fontSize.legend = 11;
fontSize.label = 16;
fontSize.colorbar = 13;
fontSize.profileName= 19;

% Select highest heat flux value:
zMax = 0.8*max(max(max(q_IJR{1})));
zMax = 13;
zMin = 1;
zMin_caxis = -0.5;

% Surface of contour:
useSurf    =  0;
useContour1 = 1;
useContour2 = 0;

% Heat flux 2D map:
% -------------------------------------------------------------------------
for kk = shotsToAnalyze
    figure('color','w')

    % Plot heat flux:
    if useSurf == 1
        surf(xx*1e3,yy*1e3,heatFlux{kk},'LineStyle','none')
    elseif useContour1 == 1
        zLevels = linspace(zMin,zMax,15);
        [cCon{kk},hCon{kk}] = contourf(xx*1e3,yy*1e3,heatFlux{kk},zLevels);
        set(hCon{kk},'LineStyle','none');
    end

    % Contours lines:
    if useContour2 == 1
        zLevels = [1,1];
        [cCon2{kk},hCon2{kk}] = contour(xx*1e3,yy*1e3,heatFlux{kk},zLevels);
        set(hCon2{kk},'LineStyle','-','lineColor','k');
    end

   % Plot hot spot contour line:
    hold on
    plot(xx_peaked{kk}*1e3,yy_peaked{kk}*1e3,'g','LineWidth',3)
    
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
    
    % Rectangle were halo heat flux is calculated:
    switch kk
        case 1          
            dx = 3.5;
            dy = 3.5;
            mean_x = -2;
            mean_y = 29.5;
            rng_x = find(xx*1e3 > mean_x - dx & xx*1e3 < mean_x + dx);
            rng_y = find(yy*1e3 > mean_y - dy & yy*1e3 < mean_y + dy);            
        case 2            
            dx = 3.5;
            dy = 3.5;
            mean_x = -2;
            mean_y = 23.5;
            rng_x = find(xx*1e3 > mean_x - dx & xx*1e3 < mean_x + dx);
            rng_y = find(yy*1e3 > mean_y - dy & yy*1e3 < mean_y + dy);  
        case 3           
            dx = 3.5;
            dy = 3.5;
            mean_x = -2;
            mean_y = 18.5;
            rng_x = find(xx*1e3 > mean_x - dx & xx*1e3 < mean_x + dx);
            rng_y = find(yy*1e3 > mean_y - dy & yy*1e3 < mean_y + dy);            
    end
    
    % Plot halo heat flux sampling region:
    hold on
    [xxg,yyg] = meshgrid(xx(rng_x),yy(rng_y));
    hR = rectangle('Position',[ mean_x - dx, mean_y - dy, 2*dx, 2*dy]);
    hR.EdgeColor = 'k';
    hR.LineWidth = 2;
    

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

%% Halo and peak heat flux scaling:

% Select area and extract data from:
for kk = shotsToAnalyze

%     zz13 = heatFlux{kk} - substractFlag*heatFlux{kk};

    % Plot data:
    figure('color','w')
    [cCon3{kk},hCon3{kk}] = contourf(xx*1e3,yy*1e3,heatFlux{kk},zLevels);
    
    % Calculate halo heat flux sampling region:
    switch kk
        case 1          
            dx = 3.5;
            dy = 3.5;
            mean_x = -2;
            mean_y = 29.5;
            rng_x = find(xx*1e3 > mean_x - dx & xx*1e3 < mean_x + dx);
            rng_y = find(yy*1e3 > mean_y - dy & yy*1e3 < mean_y + dy);            
        case 2            
            dx = 3.5;
            dy = 3.5;
            mean_x = -2;
            mean_y = 23.5;
            rng_x = find(xx*1e3 > mean_x - dx & xx*1e3 < mean_x + dx);
            rng_y = find(yy*1e3 > mean_y - dy & yy*1e3 < mean_y + dy);
            
        case 3
            rng_x = find(xx*1e3 > -7 & xx*1e3 < 3);
            rng_y = find(yy*1e3 > 15 & yy*1e3 < 22);
            
            dx = 3.5;
            dy = 3.5;
            mean_x = -2;
            mean_y = 18.5;
            rng_x = find(xx*1e3 > mean_x - dx & xx*1e3 < mean_x + dx);
            rng_y = find(yy*1e3 > mean_y - dy & yy*1e3 < mean_y + dy);            
    end
    
    % Calculate halo heat flux:
    [xxg,yyg]   = meshgrid(xx(rng_x),yy(rng_y));
    q_halo(kk)  = mean(heatFlux{kk}(rng_y,rng_x),'all');
    dq_halo(kk) = q_halo(kk)*0.17/2;
    
    % Plot halo heat flux sampling region:
    hold on
    plot(xxg*1e3,yyg*1e3,'k.','MarkerSize',4)
    
    % Calculate halo heat flux sampling region:
    switch kk
        case 1
            di = 1;
            mean_x = 23.75;
            mean_y = -4.09;
            rng_x = find(xx*1e3 > mean_x - di & xx*1e3 < mean_x + di);
            rng_y = find(yy*1e3 > mean_y - di & yy*1e3 < mean_y + di);
        case 2          
            di = 2;
            mean_x = 18.7;
            mean_y = 0.9;
            rng_x = find(xx*1e3 > mean_x - di & xx*1e3 < mean_x + di);
            rng_y = find(yy*1e3 > mean_y - di & yy*1e3 < mean_y + di);
        case 3
            di = 1;
            mean_x = 16.8;
            mean_y = 2.2;
            rng_x = find(xx*1e3 > mean_x - di & xx*1e3 < mean_x + di);
            rng_y = find(yy*1e3 > mean_y - di & yy*1e3 < mean_y + di);
        end

    % Calculate peak heat flux:
    [xxg,yyg] = meshgrid(xx(rng_x),yy(rng_y));
    q_peak(kk)  = mean(heatFlux{kk}(rng_y,rng_x),'all');
    dq_peak(kk) = q_peak(kk)*0.17/2;
    
    % Plot peak heat flux sampling region:
    plot(xxg*1e3,yyg*1e3,'k.','MarkerSize',4)
    
    % Format data:
    set(hCon3{kk},'LineStyle','none');
    colorbar
end

% Plot magnetic field scaling of heat fluxes:
% =========================================================================
fontSize.label  = 17;
fontSize.legend = 16;
fontSize.axes   = 16;
fontSize.text   = 17;

figure('color','w')
hold on

if addUncertainty
    hB(1) = errorbar(B_transport,q_peak,dq_peak,'rsq','MarkerFaceColor','r');
    hB(2) = errorbar(B_transport,q_halo,dq_halo,'ko','MarkerFaceColor','k');
else
    hB(1) = plot(B_transport,q_peak,'rsq','MarkerFaceColor','r');
    hB(2) = plot(B_transport,q_halo,'ko','MarkerFaceColor','k');
end

xlim([0,1.2])
ylim([0,15])
box on

% Linear trend:
P28 = polyfit(B_transport,q_peak,1);
P13 = polyfit(B_transport,q_halo,1);
X = linspace(0,1.2);
y28 = polyval(P28,X);
y13 = polyval(P13,X);
plot(X,y28,'r','lineStyle','--','LineWidth',1);
plot(X,y13,'k','lineStyle','--','LineWidth',1);

% Text A:
hText = text(0.19*2,10-0.5,'A');
set(hText,'FontSize',fontSize.text,'interpreter','latex')
hText.EdgeColor = 'r';
hText.Color = 'r';

% Text B:
hText = text(0.29*2,7.7-0.5,'B');
set(hText,'FontSize',fontSize.text,'interpreter','latex')
hText.EdgeColor = 'bl';
hText.Color = 'bl';

% Text C:
hText = text(0.42*2,4-0.5,'C');
set(hText,'FontSize',fontSize.text,'interpreter','latex')
hText.EdgeColor = 'k';
hText.Color = 'k';

% Formatting:
grid on
set(gca,'FontName','Times','FontSize',fontSize.axes);
set(hB,'MarkerSize',9,'LineWidth',1.5)
hLeg = legend(hB,'Peaked w/ 28 GHz','Halo w/ 28 GHz');
set(hLeg,'interpreter','Latex','fontSize',fontSize.legend);
ylabel('Heat flux [MWm$^{-2}$]','interpreter','Latex','fontSize',fontSize.label);
% xlabel('Target $|B|$ [T]','interpreter','Latex','fontSize',fontSize.label);
xlabel('peak $|B|$ near z = 3.8 m, [T]','interpreter','Latex','fontSize',fontSize.label);

% Save figure:
% =========================================================================
if saveFig
    figureName = 'Step_2d_Peak_Halo_heatflux_MagneticScaling';

    % PDF figure:
    exportgraphics(gcf,[figureName,'.pdf'],'Resolution',300) 

    % TIFF figure:
    exportgraphics(gcf,[figureName,'.tiff'],'Resolution',600) 
end

%% Plot: Total power to target:

fontSize.legend = 11;

% All time traces of pwr to Target:
% =========================================================================
figure('color','w')
hold on
for ii = 1:numel(pwr_Target)
    hP(ii) = plot(pwr_Target{ii}*1e-3);
end
set(hP,'LineWidth',1.5)
% xlim([0,15])
ylim([0,15])
box on

% Peak power traces as a function of target magnetic field:
% =========================================================================
fontSize.label  = 17;
fontSize.legend = 14;
fontSize.axes   = 16;
fontSize.text   = 17;

figure('color','w')
hold on

if addUncertainty 
    hB(1) = errorbar(B_transport,pwr_ech_peak*1e-3,d_pwr_ech_peak*1e-3,'rsq','MarkerFaceColor','r');
    hB(2) = errorbar(B_transport,pwr_helicon_mean*1e-3,d_pwr_helicon_mean*1e-3,'k^','MarkerFaceColor','k');
else
    hB(1) = plot(B_transport,pwr_ech_peak*1e-3,'rsq','MarkerFaceColor','r');
    hB(2) = plot(B_transport,pwr_helicon_mean*1e-3,'ko','MarkerFaceColor','k');
end

% Add Peaked power contribution:
y = pwr_peaked*1e-3;
hB(3) = errorbar(B_transport,y,y*0.17/2,'gsq');

% Add Halo power contribution:
y = pwr_halo*1e-3;
hB(4) = errorbar(B_transport,y,y*0.17/2,'ksq');

set(hB,'MarkerSize',6,'LineWidth',1.5)
% xlim([0,0.6])
xlim([0,1.2])
ylim([0,20])
box on

% Linear trend:
P28 = polyfit(B_transport,pwr_ech_peak*1e-3    ,1);
P13 = polyfit(B_transport,pwr_helicon_mean*1e-3,1);
X = linspace(0,1.2);
y28 = polyval(P28,X);
y13 = polyval(P13,X);
plot(X,y28,'r','lineStyle','--','LineWidth',1);
plot(X,y13,'k','lineStyle','--','LineWidth',1);

% Text A:
hText = text(0.19*2,10.5-0.5 + 4.3,'A');
set(hText,'FontSize',fontSize.text,'interpreter','latex')
hText.EdgeColor = 'r';
hText.Color = 'r';

% Text B:
hText = text(0.29*2,8-0.5 + 3.8,'B');
set(hText,'FontSize',fontSize.text,'interpreter','latex')
hText.EdgeColor = 'bl';
hText.Color = 'bl';

% Text C:
hText = text(0.42*2,4.65-0.5 + 3.5,'C');
set(hText,'FontSize',fontSize.text,'interpreter','latex')
hText.EdgeColor = 'k';
hText.Color = 'k';

% Formatting:
grid on
set(gca,'FontName','Times','FontSize',fontSize.axes);
set(hB,'MarkerSize',9)
hLeg = legend(hB,'Total power w/ 28 GHz','Helicon only','Peaked','halo');
set(hLeg,'interpreter','Latex','fontSize',fontSize.legend,'Box','off');
ylabel('Power on Target [kW]','interpreter','Latex','fontSize',fontSize.label);
xlabel('peak $|B|$ near z = 3.8 m, [T]','interpreter','Latex','fontSize',fontSize.label);


% Save figure:
% =========================================================================
if saveFig
    figureName = 'Step_2d_HeatFlux_MagneticScaling';

    % PDF figure:
    exportgraphics(gcf,[figureName,'.pdf'],'Resolution',300) 

    % TIFF figure:
    exportgraphics(gcf,[figureName,'.tiff'],'Resolution',600) 
end

%% Plotting efficiency:
% -------------------------------------------------------------------------
fontSize.axes = 17;
fontSize.label = 17;

% Calculate efficiency:
eta28 = (pwr_ech_peak - pwr_helicon_mean)/68e3; 

% Uncertainty in efficiency based on uncertainty in emissivity:
deta28 = 0.17*eta28/2;

figure('color','w')
hold on

if addUncertainty
    hEff = errorbar(B_transport,100*eta28,100*deta28,'ksq');
else
    hEff = plot(B_transport,100*eta28,'ksq');
end

% Labels:
ylabel('Efficiency [\%]','Interpreter','Latex','FontSize',fontSize.label)
% xlabel('Target $|B|$ [T]','interpreter','Latex','fontSize',fontSize.label);
xlabel('peak $|B|$ near z = 3.8 m, [T]','interpreter','Latex','fontSize',fontSize.label);

% Text A:
hText = text(0.19*2,19+1,'A');
set(hText,'FontSize',fontSize.text,'interpreter','latex')
hText.EdgeColor = 'r';
hText.Color = 'r';

% Text B:
hText = text(0.29*2,14+1.5,'B');
set(hText,'FontSize',fontSize.text,'interpreter','latex')
hText.EdgeColor = 'bl';
hText.Color = 'bl';

% Text C:
hText = text(0.42*2,7+3,'C');
set(hText,'FontSize',fontSize.text,'interpreter','latex')
hText.EdgeColor = 'k';
hText.Color = 'k';

% Formatting:
grid on
set(gca,'FontName','Times','FontSize',fontSize.axes);
set(hEff,'Markersize',9,'LineWidth',2)
xlim([0,1.2])
ylim([0,25])
box on
set(gca,'fontName','Times','fontSize',fontSize.axes)

% Save figure:
% =========================================================================
if saveFig
    figureName = 'step_2d_HeatingEfficiency';

    % PDF figure:
    exportgraphics(gcf,[figureName,'.pdf'],'Resolution',300) 

    % TIFF figure:
    exportgraphics(gcf,[figureName,'.tiff'],'Resolution',600) 
end