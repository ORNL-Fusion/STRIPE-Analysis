% Postprocess heat flux calculations:
% Load heat flux 2D maps for MPEX limit case power
% Integrate the total power on window and extrapolate to 200 kW
% Scale highest heat flux map to 200 kW

clear all
close all
clc

%% Coupled power:
% Rick gave me values for the vacuum loading and plasma loading typical of
% Proto-MPEX:
Rv = 0.19; % Ohms
Rp = 2; % Ohms, this ranges between 2 - 6 Ohms

%% 1- Load all data:
% =========================================================================

ShotSeriesToAnalyze = 1:8;

% Location where the data is located:
dataAddress = 'C:\Users\nfc\Documents\ProtoMPEX_Ops\2020\2020_02\HeliconWindow_IR\HomographicReconstruction\Method_2\';

disp('Loading data...')

t0 = tic;
for kk = ShotSeriesToAnalyze
    fileName{kk} = ['HeatfluxData_ShotSeries_',num2str(kk),'.mat'];
    filePath = [dataAddress,fileName{kk}];
    load(filePath);
    S{kk} = w;
    clear w
end
t0 = toc(t0);

disp(['Time to load data: ',num2str(t0)])

%% 2- Animate heat flux data:
% =========================================================================

animateData = 1;

if animateData
    % Select shot series:
    kk = 1;
    
    figure('color','w')
    
    % Window length:
    La = 30/100;
    
    % Window length scaling factor:
    bb = S{kk}.bb;
    
    % Mean RF power:
    mean_RF = round(mean(S{kk}.rfPwr));
    
    % Animation:
    videoFrames = 5:1:40;
    for fr = videoFrames
            % time-averaged heat flux:
            heatflux = (S{kk}.qnorm(:,:,fr) + S{kk}.qnorm(:,:,fr+1) + S{kk}.qnorm(:,:,fr+2))/3;
            
            % Plot surface:
            surf(S{kk}.phi_2D,S{kk}.z_2D*1e2,heatflux*1e3,'LineStyle','none')
            
            % Formatting:
            set(gca,'XTick',[0:45:360],'XDir','reverse')
            set(gca,'YTick',[0:5:(La/bb)*1e2])
            view([0,90])
            axis tight
            colormap('bone')
            colormap('hot')
            xlim([0,360])
            caxis([0,800])
            colorbar
            
            % Labels:
            title(['Heat flux [kW/m^2] , ',cell2mat(S{kk}.limitMode),' at ',num2str(mean_RF),' kW , frame: ',num2str(fr)])
            xlabel('Angle [deg]','Interpreter','Latex','FontSize',14)
            ylabel('z [cm]','Interpreter','Latex','FontSize',14)
            
            % Draw antenna:
            lw = 13;
            % Transverse straps:
            line([360,000],[26,26],[500,500],'color','k','LineWidth',lw)
            line([360,000],[04,04],[500,500],'color','k','LineWidth',lw)
            % Bottom helical strap:
            line([225,135],[04,26],[500,500],'color','k','LineWidth',lw)
            % HV top side helical strap:
            line([045,000],[04,15],[500,500],'color','k','LineWidth',lw)
            % GND top side helical strap:
            line([360,315],[15,26],[500,500],'color','k','LineWidth',lw)
            
            
            drawnow
            pause(0.1)
    end
end
%% 3- Integration of the surface heat flux::
% =========================================================================

% Window internal radius:
Ra = 6.25/100; %[m]

% Factor to account for missing areas of image:
% IR data has a 300 deg view
ff = 1.2; % 360 deg /300 deg = 1.2

% Calculate Steady-state total heat, ProtoMPEX:
for kk = ShotSeriesToAnalyze
    % Increments:
    dz     = S{kk}.z_2D(2,1) - S{kk}.z_2D(1,1);
    dtheta = (S{kk}.phi_2D(1,2) - S{kk}.phi_2D(1,1))*pi/180;
    ds     = dtheta*Ra; % arc length
    
    % Total number of frames:
    maxfr(kk)  = size(S{kk}.qnorm,3);
    
    % Characteristic heat flux:
    q0 = S{kk}.q0;
    
    % Loop over frames:
    for fr = 1:maxfr(kk)
        % 2D surface heat flux:
        q = S{kk}.qnorm(:,:,fr)*q0;
        
        % Heating rate per unit length:
        % Integrate q over angle
        dPdz{kk}(:,fr) = ff*sum(q,2,'omitnan' )*ds;
        
        % Total heating rate Q [W]:
        % Integrate q over angle and axial distance
        Q{kk}(fr) = ff*sum(sum(q),'omitnan' )*dz*ds;
    end
    
    % Time vector:
    t_star = S{kk}.t_star;
    t_Q{kk} = S{kk}.t_qnorm*t_star;
end

% Time-averaged heating rate:
t1 = 0.085;
t2 = 0.53;
for kk = ShotSeriesToAnalyze
    rng = find(t_Q{kk} >= t1 & t_Q{kk} <=t2);
    Q_timeAverage(kk) = mean(Q{kk}(rng));
    dQ_timeAverage(kk) = std(Q{kk}(rng),1);
end

% Steady-state heating rate:
t1 = 0.45;
t2 = 0.5;
for kk = ShotSeriesToAnalyze
    rng = find(t_Q{kk} >= t1 & t_Q{kk} <=t2);
    Q_steadyState(kk) = mean(Q{kk}(rng));
    dQ_steadyState(kk) = std(Q{kk}(rng),1);
    RF(kk)   = mean(S{kk}.rfPwr);
end

% Steady-state heating rate data from UCSD in hydrogen with window limit:
RF_UCSD            = [2   , 4  , 6  , 8  , 10 , 12 ,  14,  16,  18,  20];
Q_SteadyState_UCSD = [0.58,1.33,1.87,2.34,2.81,3.38,3.91,4.41,4.93,5.24];

%% Plotting data:
% =========================================================================

% -------------------------------------------------------------------------
% Heating rate per unit length:
% -------------------------------------------------------------------------
% Select shot series to compare:
seriesA = 1;
seriesB = seriesA+4;

% Select frames for animation:
frEnd =  50;
animationFrames = 1:4:frEnd;

figure('color','w')
for fr = animationFrames
    
    % Plot data:
    hdP(1) = plot(dPdz{seriesA}(:,fr)*1e-3,S{seriesA}.z_2D(:,1)*1e2,'k','LineWidth',2);
    hold on
    hdP(2) = plot(dPdz{seriesB}(:,fr)*1e-3,S{seriesB}.z_2D(:,1)*1e2,'r','LineWidth',2);
    
    % Draw Transverse traps:
    line([0,300],[ 4, 4],'Color','k','LineWidth',10)
    line([0,300],[26,26],'Color','k','LineWidth',10)

    hold off
    
    % Legend:
    lP = legend(hdP,{'MPEX-limiter','Window-limiter'});
    set(lP,'interpreter','Latex','FontSize',12,'box','off','Location','best')

    % Labels:
    ylabel('z [cm]','interpreter','Latex','FontSize',14)
    xlabel('[kW m$^{-1}$]','interpreter','Latex','FontSize',14)
    
    % Limits:
    xlim([0,300])
    ylim([0,30])
    
    drawnow
    pause(0.1)
end

% -------------------------------------------------------------------------
% Heating rate as a function of time:
% -------------------------------------------------------------------------
figure('color','w')
hold on
C = {'k','k','k','k','r','r','r','r'};
for kk = ShotSeriesToAnalyze
    t_star = S{kk}.t_star;
    t_Q{kk} = S{kk}.t_qnorm*t_star;
    hQ(kk) = plot(t_Q{kk},Q{kk}*1e-3,C{kk},'LineWidth',2);
end
legendText{2} = 'Window-limiter';
legendText{1} = 'MPEX-limiter';
lQ = legend([hQ(1),hQ(5)],legendText);
set(lQ,'interpreter','Latex','FontSize',14)
set(gca,'fontName','Times','FontSize',12)
xlim([0,0.65])
ylim([0,100])
xlabel('t [s]','Interpreter','Latex','FontSize',13)
ylabel('[kW]','Interpreter','Latex','FontSize',13)
box on
grid on

% -------------------------------------------------------------------------
% Steady-state heating rate as a function RF power:
% -------------------------------------------------------------------------
figure('color','w')
hold on

% Linear fits:
P{1} = polyfit(RF(1:4),Q_steadyState(1:4),1);
P{2} = polyfit(RF(5:8),Q_steadyState(5:8),1);
xRF = linspace(60,200); % [W]

% MPEX-limiter:
hQss(1) = errorbar(RF(1:4),Q_steadyState(1:4)*1e-3,dQ_steadyState(1:4)*1e-3,'ko-','lineWidth',2);
plot(xRF,polyval(P{1},xRF)*1e-3,'k--','lineWidth',2)

% Window-limiter:
hQss(2) = errorbar(RF(5:8),Q_steadyState(5:8)*1e-3,dQ_steadyState(5:8)*1e-3,'ro-','lineWidth',2);
plot(xRF,polyval(P{2},xRF)*1e-3,'r--','lineWidth',2)

% UCSD data in hydrogen:
hQss(3) = plot(RF_UCSD,Q_SteadyState_UCSD,'blsq','MarkerFaceColor','bl','MarkerSize',7);

% Formatting:
set(gca,'Fontname','Times','FontSize',12)
xlim([0,200])
ylim([0,90])
box on
grid on

% Legend:
legendText{1} = ['MPEX-limiter, Proto-MPEX'];
legendText{2} = ['Window-limiter, Proto-MPEX'];
legendText{3} = ['Window-limiter, RF-PISCES'];
lQss = legend([hQss,hQss],legendText);
set(lQss,'interpreter','Latex','FontSize',12,'Location','NorthWest')

% Labels:
xlabel('RF power [kW]','interpreter','LAtex','FontSize',13)
ylabel('Total Heat [kW]','interpreter','LAtex','FontSize',13)
title('Steady-state heating rate','interpreter','LAtex','FontSize',14)

% Save figure:
saveFig = 1;
figureName = 'SteadyStateHeatingRate_vs_RfPwr';

if saveFig
    saveas(gcf,figureName,'tiffn')
end

% Plot data relative to coupled power:
a1 = [gca,lQss];
f2 = figure('color','w');
a2 = copyobj(a1,f2);
% Convert X data to the coupled power:
for ii = 1:numel(a2(1).Children)
    Pnet = a2(1).Children(ii).XData;
    Pcoupled = Pnet*(Rp/(Rv+Rp));
    a2(1).Children(ii).XData = Pcoupled;
end
set(gca,'XTick',[0:25:200])
xlabel('Coupled RF power [kW]')
% Save figure:
saveFig = 1;
figureName = 'SteadyStateHeatingRate_vs_CoupledRfPwr';

if saveFig
    saveas(gcf,figureName,'tiffn')
end

% Note:
% The linear scaling shows that at 200 kW the total power to the window in
% the MPEX limiting case reaches 29 kW while at 141 kW it is 20 kW;
% Hence the fractional increase from 141 to 200 kW is 1.45, hence it is 45%
% higher
% For the Window-limiter, the factor from 132 to 200 kW is 1.6


% -------------------------------------------------------------------------
% Time-averaged heating rate as a function RF power:
% -------------------------------------------------------------------------
% This data takes into account the transient response of the heat flux and
% so this quantity is directly comparable to the FP data:

figure('color','w')
hold on

% Linear fits:
P{1} = polyfit(RF(1:4),Q_timeAverage(1:4),1);
P{2} = polyfit(RF(5:8),Q_timeAverage(5:8),1);
xRF = linspace(60,200); % [W]

% MPEX-limiter:
hQss(1) = plot(RF(1:4),Q_timeAverage(1:4)*1e-3,'ko-','lineWidth',2);
plot(xRF,polyval(P{1},xRF)*1e-3,'k--','lineWidth',2)

% Window-limiter:
hQss(2) = plot(RF(5:8),Q_timeAverage(5:8)*1e-3,'ro-','lineWidth',2);
plot(xRF,polyval(P{2},xRF)*1e-3,'r--','lineWidth',2)

% UCSD data in hydrogen:
hQss(3) = plot(RF_UCSD,Q_SteadyState_UCSD,'blsq','MarkerFaceColor','bl','MarkerSize',7);

% Formatting:
set(gca,'Fontname','Times','FontSize',12)
xlim([0,200])
ylim([0,90])
box on
grid on

% Legend:
legendText{1} = ['MPEX-limiter, Proto-MPEX'];
legendText{2} = ['Window-limiter, Proto-MPEX'];
legendText{3} = ['Window-limiter, CSDX'];
lQss = legend([hQss,hQss],legendText);
set(lQss,'interpreter','Latex','FontSize',14,'Location','NorthWest')

% Labels:
xlabel('RF power [kW]','interpreter','LAtex','FontSize',13)
ylabel('Total Heat [kW]','interpreter','LAtex','FontSize',13)
title('Time-averaged heating rate','interpreter','LAtex','FontSize',14)

% Save figure:
saveFig = 1;
figureName = 'timeAveragedHeatingRate_vs_RfPwr';

if saveFig
    saveas(gcf,figureName,'tiffn')
end

%% Create and save heat flux movie:
% =========================================================================

% Choose shotseries:
shotsToAnalyze = [1:8];

saveVideo = 1;

% Loop for all shot series:
for kk = shotsToAnalyze
    
    % Write to CLI:
    disp(['Video for shot series: ', num2str(kk)])

    % Create a video writer object:
    videoName = ['video_HeatFluxShotSeries_',num2str(kk),'.avi'];
    v = VideoWriter(videoName);
    open(v);

    % Start of animation:
    figure('color','w')
    mean_RF = RF(kk);
    La = 30/100;
    bb = S{kk}.bb;
    extrapolFactor(1:4) = 1;
    extrapolFactor(5:8) = 1;
    maxFr = size(S{kk}.qnorm,3);

    % Loop for all frames during the RF pulse:
    for fr = 1:1:60
            switch kk
                case {1,2,3,4}
                    y{1}(:,:,fr) = q0*(S{kk}.qnorm(:,:,fr) + S{kk}.qnorm(:,:,fr+1) + S{kk}.qnorm(:,:,fr+2)+ S{kk}.qnorm(:,:,fr+3))/4;
                    surf(S{kk}.phi_2D,S{kk}.z_2D*1e2,y{1}(:,:,fr)*1e-3,'LineStyle','none')
                case {5,6,7,8}
                    y{2}(:,:,fr) = q0*(S{kk}.qnorm(:,:,fr) + S{kk}.qnorm(:,:,fr+1) + S{kk}.qnorm(:,:,fr+2))/3;
                    surf(S{kk}.phi_2D,S{kk}.z_2D*1e2,y{2}(:,:,fr)*1e-3,'LineStyle','none')
            end
            
            % Formatting:
            set(gca,'XTick',[0:45:360],'XDir','reverse')
            set(gca,'YTick',[0:5:(La/bb)*1e2])
            xlabel('Angle [deg]','Interpreter','Latex','FontSize',14)
            ylabel('z [cm]','Interpreter','Latex','FontSize',14)
            view([0,90])
            axis tight
            colormap('bone')
            colormap('hot')
            xlim([0,360])
            caxis([0,800])
            colorbar
            title(['Heat flux [kW/m^2] , ',cell2mat(S{kk}.limitMode),' at ',num2str(mean_RF),' kW , frame: ',num2str(fr)])

            % Draw antenna:
            lw = 13;
            % Transverse straps:
            line([360,000],[26,26],[500,500],'color','k','LineWidth',lw)
            line([360,000],[04,04],[500,500],'color','k','LineWidth',lw)
            % Bottom helical strap:
            line([225,135],[04,26],[500,500],'color','k','LineWidth',lw)
            % HV top side helical strap:
            line([045,000],[04,15],[500,500],'color','k','LineWidth',lw)
            % GND top side helical strap:
            line([360,315],[15,26],[500,500],'color','k','LineWidth',lw)
            set(gca,'PlotBoxAspectRatio',[1 0.5 1])
            drawnow

            % Save frame and write to video object:
            if saveVideo
                 frame = getframe(gcf);
                 writeVideo(v,frame);
            end
    end
    
    close(v)
end

%% 200 kW heat fluxes:

% Choose shotseries:
kk = 1;

% Create a video writer object:
saveVideo = 1;
videoName = ['video_HeatFluxExtrapolation_MPEX_200kW.avi'];
v = VideoWriter(videoName);
open(v);

figure('color','w')
mean_RF = 200;
La = 30/100;
bb = S{kk}.bb;
% Extrapolation factor for MPEX-limiter:
extrapolFactor(1) = 1.45;
% Extrapolation factor for Window-limiter:
extrapolFactor(5) = 1.60;

% Maximum number of frames:
maxFr = size(S{kk}.qnorm,3);

% Plot data extrapolated to 200 kW:
for fr = 1:1:60
        switch kk
            case 1 % 142 kW with MPEX-limiter scenario:
                y{1}(:,:,fr) = q0*extrapolFactor(kk)*(S{kk}.qnorm(:,:,fr) + S{kk}.qnorm(:,:,fr+1) + S{kk}.qnorm(:,:,fr+2)+ S{kk}.qnorm(:,:,fr+3))/4;
                surf(S{kk}.phi_2D,S{kk}.z_2D*1e2,y{1}(:,:,fr)*1e-3,'LineStyle','none')
            case 5 % 132 kW with Window-limiter scenario:
                y{2}(:,:,fr) = q0*extrapolFactor(kk)*(S{kk}.qnorm(:,:,fr) + S{kk}.qnorm(:,:,fr+1) + S{kk}.qnorm(:,:,fr+2))/3;
                surf(S{kk}.phi_2D,S{kk}.z_2D*1e2,y{2}(:,:,fr)*1e-3,'LineStyle','none')
        end
        
        % Formatting:
        set(gca,'XTick',[0:45:360],'XDir','reverse')
        set(gca,'YTick',[0:5:(La/bb)*1e2])
        xlabel('Angle [deg]','Interpreter','Latex','FontSize',14)
        ylabel('z [cm]','Interpreter','Latex','FontSize',14)
        view([0,90])
        axis tight
        colormap('bone')
        colormap('hot')
        xlim([0,360])
        caxis([0,1200])
        colorbar
        title(['Heat flux [kW/m^2] , ',cell2mat(S{kk}.limitMode),' at ',num2str(mean_RF),' kW , frame: ',num2str(fr)])

        % Draw antenna:
        lw = 13;
        % Transverse straps:
        line([360,000],[26,26],[500,500],'color','k','LineWidth',lw)
        line([360,000],[04,04],[500,500],'color','k','LineWidth',lw)
        % Bottom helical strap:
        line([225,135],[04,26],[500,500],'color','k','LineWidth',lw)
        % HV top side helical strap:
        line([045,000],[04,15],[500,500],'color','k','LineWidth',lw)
        % GND top side helical strap:
        line([360,315],[15,26],[500,500],'color','k','LineWidth',lw)
        set(gca,'PlotBoxAspectRatio',[1 0.5 1])
        drawnow

        if saveVideo
             frame = getframe(gcf);
             writeVideo(v,frame);
        end
end

close(v)

%% Filling in the gaps for the MPEX-limiter 200k W heat flux distribution
close all

% Select the shot series:
kk = 1;

% Plot figure:
figure;
% Select frame at the end of the RF pulse where the 
fr = 40;
surf(y{1}(:,:,fr)*1e-3,'LineStyle','none')
view([0,90])
axis tight
colormap('bone')
colormap('hot')
xlim([0,100])
caxis([0,1000])
colorbar
title(['Heat flux [kW/m^2] , ',cell2mat(S{kk}.limitMode),' at ',num2str(mean_RF),' kW , frame: ',num2str(fr)])

% Start the "gap filling" process:
% -------------------------------------------------------------------------

% Create the variable to modify:
% rr is in Wm^-2
rr = y{1};

rngrow = [1:50]';
rngcol = [15:30];
rngcolnew = [85:100]';
%     rndrng1 = randi(5,5,1);
% rndrng2 = randi(9,9,1);

% Apply modifications:
for fr = 1:size(rr,3)
    rr(:,rngcolnew,fr) = rr(:,rngcol(randperm(16)),fr);
    rr(:,1:16,fr) = rr(:,rngcol(randperm(16)),fr);

    if 0
        dum1 = rr(:,65:70,fr);
        rr(:,60:65,fr) = circshift(dum1,3,1);
        dum1 = rr(:,65:70,fr);
        rr(:,55:60,fr) = circshift(dum1,6,1);
    elseif 0
        dum1 = rr(:,[68:76]-3,fr);
        rr(:,[59:67],fr) = fliplr(flipud(dum1));
        
        dum1 = rr(:,[73:80],fr);
        rr(:,[52:59],fr) = fliplr(flipud(dum1));
        
        dum1 = rr(:,[80:92],fr);
        rr(:,[40:52],fr) = fliplr(flipud(dum1));
     elseif 1
        dum1 = rr(:,[65:76]+1,fr);
        rr(:,[53:64],fr) = circshift(fliplr(flipud(dum1)),1,1);
        
        rng = [71:80];
        dum1 = rr(:,rng,fr);
        rr(:,[50:59],fr) = fliplr(flipud(dum1));
        
        rng = [80:90];
        dum1 = rr(:,rng,fr);
        rng = [40:50]+2;
        rr(:,rng,fr) = fliplr(flipud(dum1));
        
        dum1 = rr([1:5],[1:31],fr);
        rr([1:5],[40:70],fr) = fliplr(flipud(dum1));
        
        dum1 = rr([1:4],[1:100],fr);
        rr([47:50],[1:100],fr) = fliplr(flipud(dum1));
     end
end

for fr = 1:size(rr,3)-3
            rr(:,:,fr) = (rr(:,:,fr)+rr(:,:,fr+1)+rr(:,:,fr+2)+rr(:,:,fr+3))/4;
end


% Plot intermediate results of the "filling" process:
figure;
for fr = 40;
    surf((1/3)*(rr(:,:,fr)+rr(:,:,fr+1)+rr(:,:,fr+2))*1e-3,'LineStyle','none')
    view([0,90])
    axis tight
    colormap('bone')
    colormap('hot')
    xlim([0,100])
    caxis([0,1200])
    colorbar
    title(['Heat flux [kW/m^2] , ',cell2mat(S{kk}.limitMode),' at ',num2str(mean_RF),' kW , frame: ',num2str(fr)])
    drawnow
end

% Apply smoothing to data:
% -------------------------------------------------------------------------
% Create a new variable called rr_smooth which contains the heat flux with
% the filled-in maps:

for fr = 40
    for ii = 1:50
            rr_smooth(ii,:,fr) = sgolay_t(rr(ii,:,fr),3,7); % in units of Wm-2
    end
end

for fr = 40
    for jj = 1:100
            rr_smooth(:,jj,fr) = sgolay_t(rr_smooth(:,jj,fr),3,7);
    end
end


%% Plot heat flux maps extrapolated to 200 kW:
% This includes both "filled-in" and "as-is" 2D heat flux maps from the
% MPEX-limiter configuration:

% not filled-in MPEX-limiter heat flux map extrapolated to 200 kW:
% This represents the data as it is obtained from the IR analsys but only
% extrapolated to 200 kW by applying a factor of 1.45
% -------------------------------------------------------------------------
figure('color','w')
mean_RF = 200;
La = 30/100;
bb = S{1}.bb;
kk = 1;
for fr = 40
        surf(S{kk}.phi_2D,S{kk}.z_2D*1e2,y{1}(:,:,fr)*1e-3,'LineStyle','none')
        set(gca,'XTick',[0:45:360],'XDir','reverse')
        set(gca,'YTick',[0:5:(La/bb)*1e2])
        xlabel('Angle [deg]','Interpreter','Latex','FontSize',14)
        ylabel('z [cm]','Interpreter','Latex','FontSize',14)
        view([0,90])
        axis tight
        colormap('bone')
        colormap('hot')
        xlim([0,360])
        caxis([0,1200])
        zlim([0,1200])
        colorbar
        title(['Heat flux [kW/m^2] , ',cell2mat(S{kk}.limitMode),' at ',num2str(mean_RF),' kW , frame: ',num2str(fr)])

        % Draw antenna:
        lw = 13;
        % Transverse straps:
        line([360,000],[26,26],[500,500],'color','k','LineWidth',lw)
        line([360,000],[04,04],[500,500],'color','k','LineWidth',lw)
        % Bottom helical strap:
        line([225,135],[04,26],[500,500],'color','k','LineWidth',lw)
        % HV top side helical strap:
        line([045,000],[04,15],[500,500],'color','k','LineWidth',lw)
        % GND top side helical strap:
        line([360,315],[15,26],[500,500],'color','k','LineWidth',lw)
        drawnow
%         pause(0.1)
end

if 1
    figureName = ['Heatflux_Partial2D_Extrapolated_200kW'];
    saveas(gcf,figureName,'tiffn')
end

% "Filled-in" MPEX-limiter heat flux map extrapolated to 200 kW:
% -------------------------------------------------------------------------
figure('color','w')
mean_RF = 200;
La = 30/100;
bb = S{1}.bb;
kk = 1;
for fr = 40
        surf(S{kk}.phi_2D,S{kk}.z_2D*1e2,rr_smooth(:,:,fr)*1e-3,'LineStyle','none')
        set(gca,'XTick',[0:45:360],'XDir','reverse')
        set(gca,'YTick',[0:5:(La/bb)*1e2])
        xlabel('Angle [deg]','Interpreter','Latex','FontSize',14)
        ylabel('z [cm]','Interpreter','Latex','FontSize',14)
        view([0,90])
        axis tight
        colormap('bone')
        colormap('hot')
        xlim([0,360])
        caxis([0,1200])
                zlim([0,1200])
        colorbar
        title(['Heat flux [kW/m^2] , ',cell2mat(S{kk}.limitMode),' at ',num2str(mean_RF),' kW , frame: ',num2str(fr)])

        % Draw antenna:
        lw = 13;
        % Transverse straps:
        line([360,000],[26,26],[500,500],'color','k','LineWidth',lw)
        line([360,000],[04,04],[500,500],'color','k','LineWidth',lw)
        % Bottom helical strap:
        line([225,135],[04,26],[500,500],'color','k','LineWidth',lw)
        % HV top side helical strap:
        line([045,000],[04,15],[500,500],'color','k','LineWidth',lw)
        % GND top side helical strap:
        line([360,315],[15,26],[500,500],'color','k','LineWidth',lw)
        drawnow
%         pause(0.1)
end

if 1
    figureName = ['Heatflux_complete2D_Extrapolated_200kW'];
    saveas(gcf,figureName,'tiffn')
end

%% Save steady-state frame for MPEX-limiter extrapolation to 200 kW

% Create a data set for the thermo-mechanical analsys:

% 2D heat flux maps in Wm^-2:
d.q2D = rr_smooth(:,:,fr);
% Coordinates of the data arrays:
% Angular coordinate:
d.phi_q2D = S{1}.phi_2D*pi/180;
% Axial coordinate:
d.z_q2D   = S{1}.z_2D;

% RF power value:
d.R_vacuum = Rv;
d.R_plasma = Rp;
d.RFpower_NET = 200*1e3;
d.RFpower_Coupled = d.RFpower_NET*(d.R_plasma/(d.R_vacuum+d.R_plasma))*1e3;

% Geometry of window:
d.RadiusWindow = 0.0620; 
d.LengthWindow = 0.3;

% Thermal parameters used for the data analysis:
% This includes the directional emissivity of the AlN window at 12 deg
% (0.55) and optics transmission factors due to ZnSe IR viewing window
% (0.7)
d.thermalParams = S{1}.thermalParam;

% Comments and metdata:
d.comment{1} = 'This heat flux 2D data is based on IR-based measurements performed in Proto-MPEX';
d.comment{2} = 'Heat flux has been linearly extrapolated from experiments performed at 80, 120 and 140 kW in Proto-MPEX';
d.comment{3} = 'q2D is the heat flux map in [W]';
d.comment{4} = 'phi_q2D is the azimuthal angle coordinate in radians. Zero radians corresponds to 12 o clock in physical space';
d.comment{5} = 'z_q2D is axial coordinate (along window) in [m]';
d.comment{6} = 'Ra is the inner radius and La is the length of the helicon window in Proto-MPEX';

fileName = ['HeatFlux2D_ExtrapolatedTo_200kW_created_',date];

% Save data:
saveData = 1;
if saveData
    disp('Saving 200 kW data ...')
    save(fileName,'d')
    disp('Saving complete!')
end
