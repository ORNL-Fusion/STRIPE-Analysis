% Step 2: 
% Post process heat flux data

clear all
close all

% Load data:
% =========================================================================
load('Step_2c_heatFlux2D.mat');

saveFig = 1;
addUncertainty = 1;

% Rename and rescale some variables:
% =========================================================================
pwr_Target   = pwr;    % Power coupled to Target in [W]
pwr28GHz     = pwr28GHz*1e3; % 28 GHz power in [W]

%% Extract power during ECH:
% =========================================================================
close all

% Time traces of coupled power:
% -------------------------------------------------------------------------
figure('color','w')
hold on

fontSize.label = 16;
fontSize.axes = 16;

colorRng = {'k','bl','c','m','g','c','k','bl','c','m','g','c'};

t0 = 0.31;
tRF_start = 4.15;

for ii = numel(shot):-1:1
    % Time trace:
    t_pwr_Target{ii} = t_R{ii}*t_star - t0 + tRF_start;
    
    % Time trace of power coupled to Target:
    hPwrTarget(ii) = plot(t_pwr_Target{ii},pwr_Target{ii}*1e-3,'LineWidth',2,'Color',colorRng{ii});
    
    % Find peak power for ECH:
    rng = find(t_pwr_Target{ii} > 4.45 & t_pwr_Target{ii} < 4.58);
    [pwr_Target_ech(ii),ipeak_ech] = max(pwr_Target{ii}(rng));
    ipeak_ech = rng(ipeak_ech);
    
    % Find Helicon power prior to ECH pulse:
    rng = find(t_pwr_Target{ii} > 4.40 & t_pwr_Target{ii} < 4.45);
    [pwr_Target_helicon(ii),ipeak_helicon] = max(pwr_Target{ii}(rng));
    ipeak_helicon = rng(ipeak_helicon);
    
    % Plot peak:
    hPwrPeak = plot(t_pwr_Target{ii}(ipeak_ech),pwr_Target_ech(ii)*1e-3);
    set(hPwrPeak,'Marker','o','MarkerFaceColor',colorRng{ii})
    
    % ECH power contribution
    pwr_ech = pwr_Target_ech - pwr_Target_helicon;
end
ylim([0,16])
ylabel('Power to target [kW]','Interpreter','Latex','FontSize',fontSize.label)
xlabel('time [s]','Interpreter','Latex','FontSize',fontSize.label)
grid on

box on
set(gca,'fontName','Times','fontSize',fontSize.axes)

% Save figure:
figureName = ['step_2d_timeEvolutionPowerToTarget'];
saveas(gcf, figureName,'tiff');

%% Extract 2D heat flux profiles:

% Select frame to substract:
fr_substract = 90; 
substractFlag = 1;

% Time averaging
timeAverage = 0;

for kk = 1:10
    % Select hottest frame:
    [~,fr] = max(pwr{kk});

    for r = fr

        if timeAverage
                heatFlux{kk} = ( q_IJR{kk}(:,:,r) + q_IJR{kk}(:,:,r+1) + q_IJR{kk}(:,:,r+2) )/3 - substractFlag*q_IJR{kk}(:,:,fr_substract); 
        else
                heatFlux{kk} = q_IJR{kk}(:,:,r) - substractFlag*q_IJR{kk}(:,:,fr_substract); 
        end
    end
end

%% Extract 50% power countour lines:

for kk = 1:10

    % Define threshold heat flux:
    qThreshold(kk) = max(max(heatFlux{kk}))*0.14;

    % Define spatial increaments:
    dxx = xx(2) - xx(1);
    dyy = yy(2) - yy(1);

    % Initialize data:
    heatFlux_peaked{kk} = heatFlux{kk};
    heatFlux_halo{kk}   = heatFlux{kk};
    
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
    frame_length = 5;
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

    if 0
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
    
    % Intergrate power [W]:
    pwr_peaked(kk) = sum(sum(heatFlux_peaked{kk}*1e6))*dxx*dxx;
    pwr_halo(kk)   = sum(sum(heatFlux_halo{kk}*1e6))*dxx*dxx;

end

figure;
hold on
hPwr(1) = plot(pwr28GHz*1e-3,pwr_Target_ech*1e-3,'ksq');
hPwr(2) = plot(pwr28GHz*1e-3,pwr_peaked*1e-3,'ro');
hPwr(3) = plot(pwr28GHz*1e-3,pwr_halo*1e-3,'go');
plot(pwr28GHz*1e-3,(pwr_halo + pwr_peaked)*1e-3,'c');

legend(hPwr,'Total','Peaked','Halo')

%% Plot: 2D heat flux

shotsToAnalyze = 1:10;

% Font sizes:
fontSize.axes = 16;
fontSize.title = 16;
fontSize.legend = 11;
fontSize.label = 16;
fontSize.colorbar = 13;

% Select highest heat flux value:
try
    kk = 10;
    pwr{kk};
catch
    kk = 1;
end
zMax = 22;
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

    % Select hottest frame:
    [~,fr] = max(pwr{kk});

    for r = fr
        
        if useSurf == 1
            surf(xx*1e3,yy*1e3,movmean(movmean(heatFlux{kk},3,1),3,2),'LineStyle','none')
        elseif useContour1 == 1
            zLevels = linspace(zMin,zMax,15);
            [cCon{kk},hCon{kk}] = contourf(xx*1e3,yy*1e3,heatFlux{kk},zLevels);
            set(hCon{kk},'LineStyle','none');
        end

        % Plot 50% power line:
        hold on
        plot(xx_peaked{kk}*1e3,yy_peaked{kk}*1e3,'g','LineWidth',3)
        
        % Plot LUFS:
        hCirc(2) = plot3(xCirc*10 - 3 ,yCirc*10,ones(size(xCirc))*zMax,'k','LineStyle',':','LineWidth',3);

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

        % DLP path:
        if kk == 9
            hDLP = line([-35,40],[-20,35]-13,[20,20]);
            set(hDLP,'lineStyle','--','Color','k','LineWidth',2)
            hSpot = plot3(24,10.1,20,'sq');
            set(hSpot,'Color','k','MarkerFaceColor','k','MarkerSize',10)

        end
        
        % Formatting:
        colormap(flipud(hot))
        caxis([zMin_caxis,zMax])
        zlim([zMin,zMax])
%         caxis([0,zMax])
%         zlim([0,zMax])
        xlim([-1,1]*45)
        ylim([-1,1]*45)
        view([-45,45])
        view([00,90])
        axis image
        hCB = colorbar;
        hCB.Ticks = [0:5:20];
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
    
    x =  [38  ,29];
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
    dx = 3.5;
    dy = 3.5;
    mean_x = -15;
    mean_y = 27.5;
    rng_x = find(xx*1e3 > mean_x - dx & xx*1e3 < mean_x + dx);
    rng_y = find(yy*1e3 > mean_y - dy & yy*1e3 < mean_y + dy);
    
    % Plot halo heat flux sampling region:
    hold on
    [xxg,yyg] = meshgrid(xx(rng_x),yy(rng_y));
    q_halo(kk) = mean(heatFlux{kk}(rng_y,rng_x),'all');
    hR = rectangle('Position',[ mean_x - dx, mean_y - dy, 2*dx, 2*dy]);
    hR.EdgeColor = 'k';
    hR.LineWidth = 2;
    
    q_peak(kk) = max(heatFlux{kk},[],'all');
    
    % Uncertainty in heat flux due to uncertainty in emissivity:
    dq_halo(kk) = 0.17*q_halo(kk)/2;
    dq_peak(kk) = 0.17*q_peak(kk)/2;
    
    % Save figure:
    % =========================================================================
    if saveFig
        figureName = ['step_2c_HeatFluxMap_',num2str(round(pwr28GHz(kk))),'_kW_shot_',num2str(shot(kk))];

        % PDF figure:
        exportgraphics(gcf,[figureName,'.pdf'],'Resolution',300) 

        % TIFF figure:
        exportgraphics(gcf,[figureName,'.tiff'],'Resolution',600) 
    end
    
    % Temporal evolution of heating rate:
    % -------------------------------------------------------------------------    
    if 0
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
            figureName = ['step_2c_IntegratedPower_',num2str(round(pwr28GHz(kk))),'_kW_shot_',num2str(shot(kk))];

            % PDF figure:
            exportgraphics(gcf,[figureName,'.pdf'],'Resolution',300) 

            % TIFF figure:
            exportgraphics(gcf,[figureName,'.tiff'],'Resolution',600) 
        end
    end
end

%% Calculate heat flux along DLP path

if 1

    % Define path:
    path{1}.xq = linspace(-35,24)+2;
    path{1}.yq = linspace(-25,10)+1;

    path{2}.xq = linspace(-35,24)+2;
    path{2}.yq = linspace(-25,10)+3;

    % Extract heat flux:
    zq_mean = zeros(size(path{1}.xq));
    for ii = 1:numel(path)
        path{ii}.zq = interp2(xx*1e3,yy*1e3,heatFlux{9},path{ii}.xq,path{ii}.yq);
        path{ii}.sq = linspace(0,7,numel(path{ii}.zq));
        zq_mean = zq_mean + path{ii}.zq;
    end
    zq_mean = zq_mean/ii;

    % Plot data:
    figure
    hold on
    surf(xx*1e3,yy*1e3,heatFlux{9},'LineStyle','none')
    for ii = 1:numel(path)
        plot3(path{ii}.xq,path{ii}.yq,path{ii}.zq,'k.')
    end
    hold on
    hCirc(2) = plot3(xCirc*10 - 3 ,yCirc*10,ones(size(xCirc))*zMax,'k','LineStyle',':','LineWidth',3);

    colormap(flipud(hot))
    caxis([zMin_caxis,zMax])
    zlim([zMin,zMax])
    view([-45,45])
    view([00,90])
    axis image
    hCB = colorbar;
    hCB.Ticks = [0:5:20];
    hold off    
    set(gca,'fontName','Times','fontSize',fontSize.axes)

    % Plot heat flux along path:
    figure('color','w');
    hold on
    for ii = 1:numel(path)
        plot(path{ii}.sq - 3.5,path{ii}.zq)
    end
    sq_mean = path{ii}.sq - 3.5;
    plot(sq_mean,zq_mean,'k','lineWidth',2)
    ylim([0,12])
    box on

    if 0
        % Save data:
        fileName = 'Step_2d_HeatFluxAlongDlpPath';
        varList = {'sq_mean','zq_mean'};
        save(fileName,varList{:})
    end
end

%% Plot: Scaling of halo heat flux coupled to Target with ECH power:
% -------------------------------------------------------------------------
fontSize.label  = 16;
fontSize.legend = 16;
fontSize.axes   = 15;
fontSize.text   = 17;

% Linear fit:
x = linspace(0,100);

figure('color','w')
hold on

if addUncertainty
    hH(1) = errorbar(pwr28GHz*1e-3,q_halo,dq_halo,'ko','MarkerFaceColor','k');
    hH(2) = errorbar(pwr28GHz*1e-3,q_peak,dq_peak,'rsq','MarkerFaceColor','r');
else
    hH(1) = plot(pwr28GHz*1e-3,q_halo,'ko','MarkerFaceColor','k');
    hH(2) = plot(pwr28GHz*1e-3,q_peak,'rsq','MarkerFaceColor','r');
end

% Linear fits:
% Halo:
hHFit(1) = plot(x,polyval(polyfit(pwr28GHz*1e-3,q_halo,1),x),'k','LineWidth',1,'LineStyle','--');
% Peak:
hHFit(2) = plot(x,polyval(polyfit(pwr28GHz*1e-3,q_peak,1),x),'r','LineWidth',1,'LineStyle','--');

% Labels:
ylabel('Heat flux [MWm$^{-2}$]','Interpreter','Latex','FontSize',fontSize.label)
xlabel('28 GHz microwave power [kW]','Interpreter','Latex','FontSize',fontSize.label)

% Formatting:
set(gca,'FontName','Times','FontSize',fontSize.axes);
set(hH,'MarkerSize',9,'LineWidth',1.5)
hLeg = legend(hH,'Halo w/ 28 GHz','Peaked w/ 28 GHz');
set(hLeg,'interpreter','Latex','fontSize',fontSize.legend,'Location','NorthWest');
xlim([0,80])
ylim([0,25])
box on
set(gca,'fontName','Times','fontSize',fontSize.axes)
grid on

% Save figure:
% =========================================================================
if saveFig
    figureName = 'step_2d_HaloHeatFlux_vs_28GHz';

    % PDF figure:
    exportgraphics(gcf,[figureName,'.pdf'],'Resolution',300) 

    % TIFF figure:
    exportgraphics(gcf,[figureName,'.tiff'],'Resolution',600) 
end

%% Plot: Scaling of power coupled to Target with ECH power:
% -------------------------------------------------------------------------
% Linear fit:
x = linspace(0,100);

fontSize.legend = 13;

figure('color','w')
hold on

if addUncertainty
    d_pwr_Target_ech = 0.17*pwr_Target_ech/2;
    hPwr(1) = errorbar(pwr28GHz*1e-3,pwr_Target_ech*1e-3,d_pwr_Target_ech*1e-3,'rsq');
else
    hPwr(1) = plot(pwr28GHz*1e-3,pwr_Target_ech*1e-3,'rsq');    
end

hlinFit = plot(x,polyval(polyfit(pwr28GHz*1e-3,pwr_Target_ech*1e-3,1),x),'r--','LineWidth',2);

% Add Peaked power contribution:
y = pwr_peaked(2:end)*1e-3;
hPwr(2) = errorbar(pwr28GHz(2:end)*1e-3,y,y*0.17/2,'gsq');

% Add Halo power contribution:
y = pwr_halo(2:end)*1e-3;
hPwr(3) = errorbar(pwr28GHz(2:end)*1e-3,y,y*0.17/2,'ksq');

% Labels:
ylabel('Power to Target [kW]','Interpreter','Latex','FontSize',fontSize.label)
xlabel('28 GHz microwave power [kW]','Interpreter','Latex','FontSize',fontSize.label)
hLeg = legend(hPwr,'Total power','Peaked','Halo');
set(hLeg,'Interpreter','Latex','FontSize',fontSize.legend,'Location','NorthWest')

% Formatting:
set(hPwr,'Markersize',9,'LineWidth',2)
xlim([0,80])
ylim([0,16])
box on
set(gca,'fontName','Times','fontSize',fontSize.axes)
grid on

% Save figure:
% =========================================================================
if saveFig
    figureName = 'step_2d_PeakPowerToTarget';

    % PDF figure:
    exportgraphics(gcf,[figureName,'.pdf'],'Resolution',300) 

    % TIFF figure:
    exportgraphics(gcf,[figureName,'.tiff'],'Resolution',600) 
end

%% Plot: efficiency
% -------------------------------------------------------------------------
figure('color','w')

if addUncertainty
    eta28   = pwr_ech./pwr28GHz;
    deta28 = 0.17*eta28/2; 
    hEff = errorbar(pwr28GHz*1e-3,100*eta28,100*deta28,'ksq');
else
    hEff = plot(pwr28GHz*1e-3,100*pwr_ech./pwr28GHz,'ksq');
end

% Labels:
ylabel('Efficiency [\%]','Interpreter','Latex','FontSize',fontSize.label)
xlabel('28 GHz microwave power [kW]','Interpreter','Latex','FontSize',fontSize.label)

% Formatting:
grid on
set(hEff,'Markersize',9,'LineWidth',2)
xlim([0,80])
ylim([0,30])
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
