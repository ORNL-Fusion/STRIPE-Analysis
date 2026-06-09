% Postprocess heat flux calculations:
% Load heat flux 2D maps for MPEX limit case power
% Integrate the total power on window and extrapolate to 175 kW coupled
% power
% Scale highest heat flux map to 175 kW coupled
% Scale heat flux to account for larger surface area of MPEX window

clear all
close all
clc

%% Coupled power:
% Rick gave me values for the vacuum loading and plasma loading typical of
% Proto-MPEX:
Rv = 0.19; % Ohms
Rp = 2; % Ohms, this ranges between 2 - 6 Ohms

% These values are required to extract the coupled RF power from the RF
% injected power.

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

% The data is stored in the structure "S". this contains the heat flux and
% the coordinates.

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
    mean_RF_NET = round(mean(S{kk}.rfPwr));
    
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
            title(['Heat flux [kW/m^2] , ',cell2mat(S{kk}.limitMode),' at ',num2str(mean_RF_NET),' kW , frame: ',num2str(fr)])
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
P_steadyState{1} = polyfit(RF(1:4),Q_steadyState(1:4),1);
P_steadyState{2} = polyfit(RF(5:8),Q_steadyState(5:8),1);
xRF = linspace(60,200); % [W]

% MPEX-limiter:
hQss(1) = errorbar(RF(1:4),Q_steadyState(1:4)*1e-3,dQ_steadyState(1:4)*1e-3,'ko-','lineWidth',2);
plot(xRF,polyval(P_steadyState{1},xRF)*1e-3,'k--','lineWidth',2)

% Window-limiter:
hQss(2) = errorbar(RF(5:8),Q_steadyState(5:8)*1e-3,dQ_steadyState(5:8)*1e-3,'ro-','lineWidth',2);
plot(xRF,polyval(P_steadyState{2},xRF)*1e-3,'r--','lineWidth',2)

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
figureName = 'Step_4c_SteadyStateHeatingRate_vs_RfPwr';

if saveFig
    saveas(gcf,figureName,'tiffn')
end

% -------------------------------------------------------------------------
% Steady-state heating rate as a function of COUPLED RF power:
% -------------------------------------------------------------------------
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
figureName = 'Step_4c_SteadyStateHeatingRate_vs_CoupledRfPwr';

if saveFig
    saveas(gcf,figureName,'tiffn')
end

% Note 2020_11_10:
% The linear scaling shows that at 200 kW the total power to the window in
% the MPEX limiting case reaches 40 kW while at 141 kW it is 29 kW;
% Hence the fractional increase from 141 to 200 kW is 1.45, hence it is 45%
% higher
% For the Window-limiter, the factor from 132 to 200 kW is 1.61


% -------------------------------------------------------------------------
% Time-averaged heating rate as a function RF power:
% -------------------------------------------------------------------------
% This data takes into account the transient response of the heat flux and
% so this quantity is directly comparable to the FP data:

figure('color','w')
hold on

% Linear fits:
% MPEX-limiter:
P_timeAverage{1} = polyfit(RF(1:4),Q_timeAverage(1:4),1);
% Window-limiter:
P_timeAverage{2} = polyfit(RF(5:8),Q_timeAverage(5:8),1);
xRF = linspace(60,200); % [W]

% MPEX-limiter:
hQss(1) = plot(RF(1:4),Q_timeAverage(1:4)*1e-3,'ko-','lineWidth',2);
plot(xRF,polyval(P_timeAverage{1},xRF)*1e-3,'k--','lineWidth',2)

% Window-limiter:
hQss(2) = plot(RF(5:8),Q_timeAverage(5:8)*1e-3,'ro-','lineWidth',2);
plot(xRF,polyval(P_timeAverage{2},xRF)*1e-3,'r--','lineWidth',2)

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
figureName = 'Step_4c_timeAveragedHeatingRate_vs_RfPwr';

if saveFig
    saveas(gcf,figureName,'tiffn')
end

%% Create and save heat flux movie:
% =========================================================================

% Choose shotseries:
shotsToAnalyze = [];%[1:8];

saveVideo = 0;

% Loop for all shot series:
for kk = shotsToAnalyze
    
    % Write to CLI:
    disp(['Video for shot series: ', num2str(kk)])

    % Create a video writer object:
    videoName = ['Step_4c_video_HeatFluxShotSeries_',num2str(kk),'.avi'];
    v = VideoWriter(videoName);
    open(v);

    % Start of animation:
    figure('color','w')
    mean_RF_NET = RF(kk);
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
            title(['Heat flux [kW/m^2] , ',cell2mat(S{kk}.limitMode),' at ',num2str(mean_RF_NET),' kW , frame: ',num2str(fr)])

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

%% 175 kW coupled power heat fluxes:

% =========================================================================
% To extrapolate to 175 kW coupled power, we need to first convert the net
% power to coupled power. This is on by using the vacuum loading and plasma
% loading
% =========================================================================

format short

% The coupled power to extrapolate to is:
mean_RF_Coupled = 175;

% Conversion factor:
net2coupled = Rp/(Rp+Rv);

% The corresponding net power is"
mean_RF_NET = mean_RF_Coupled/net2coupled

% Choose shotseries:
kk = 1;

% Create a video writer object:
saveVideo = 1;
videoName = ['Step_4c_video_HeatFluxExtrapolation_MPEX_200kW.avi'];
v = VideoWriter(videoName);
open(v);

figure('color','w')

La = 30/100;
bb = S{kk}.bb;

% Calculate extrapolation factor to apply to data:
% -------------------------------------------------------------------------
% Extrapolation factor for MPEX-limiter:
extrapolFactor(1) = polyval(P_steadyState{1},mean_RF_NET)/polyval(P_steadyState{1},RF(1))
 
% Extrapolation factor for Window-limiter:
extrapolFactor(5) = polyval(P_steadyState{2},mean_RF_NET)/polyval(P_steadyState{2},RF(5));

% Maximum number of frames:
maxFr = size(S{kk}.qnorm,3);

% Plot data extrapolated to 175 kW coupled:
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
        title(['Heat flux [kW/m^2] , ',cell2mat(S{kk}.limitMode),' at ',num2str(mean_RF_Coupled),' kW Coupled, frame: ',num2str(fr)])

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

% ##########################################################################################################################################################################################################
% ##########################################################################################################################################################################################################
% ##########################################################################################################################################################################################################
% ##########################################################################################################################################################################################################
% ##########################################################################################################################################################################################################
% ##########################################################################################################################################################################################################
% ##########################################################################################################################################################################################################

% The following material takes the Proto-MPEX data and uses it to produce
% heat flux extrapolations for the MPEX helicon window and ICH window

%% Filling in the gaps for the MPEX-limiter 175 kW COUPLED power heat flux distribution
close all

% Select the shot series:
kk = 1;

% Plot figure:
figure;

% Select frame at the end of the RF pulse:
% This represents the steady-state part of the pulse
fr = 40;

surf(y{1}(:,:,fr)*1e-3,'LineStyle','none')
view([0,90])
axis tight
colormap('bone')
colormap('hot')
xlim([0,100])
caxis([0,1000])
colorbar
title(['Heat flux [kW/m^2] , ',cell2mat(S{kk}.limitMode),' at ',num2str(mean_RF_NET),' kW , frame: ',num2str(fr)])

% Start the "gap filling" process:
% -------------------------------------------------------------------------

% Create the variable to modify:
% rr is in Wm^-2
rr = y{1};

rngrow = [1:50]';
rngcol = [15:30];
rngcolnew = [85:100]';

% Apply modifications to fill in the gaps:
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
    title(['Heat flux [kW/m^2] , ',cell2mat(S{kk}.limitMode),' at ',num2str(mean_RF_Coupled),' kW Coupled , frame: ',num2str(fr)])
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


%% Plot heat flux maps extrapolated to 175 kW COUPLED power:
% This includes both "filled-in" and "as-is" 2D heat flux maps from the
% MPEX-limiter configuration:

% not filled-in MPEX-limiter heat flux map extrapolated to 175 kW coupled:
% This represents the data as it is obtained from the IR analsys but only
% extrapolated to 175 kW coupled by applying a factor of 1.36
% -------------------------------------------------------------------------
figure('color','w')
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
        title(['Heat flux [kW/m^2] , ',cell2mat(S{kk}.limitMode),' at ',num2str(mean_RF_Coupled),' kW Coupled, frame: ',num2str(fr)])

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
    figureName = ['Step_4c_Heatflux_Partial2D_Extrapolated_175kW_Coupled'];
    saveas(gcf,figureName,'tiffn')
end

% "Filled-in" MPEX-limiter heat flux map extrapolated to 175 kW coupled:
% -------------------------------------------------------------------------
figure('color','w')
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
        title(['Heat flux [kW/m^2] , ',cell2mat(S{kk}.limitMode),' at ',num2str(mean_RF_Coupled),' kW Coupled, frame: ',num2str(fr)])

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
    figureName = ['Step_4c_Heatflux_complete2D_Extrapolated_175kW_Coupled'];
    saveas(gcf,figureName,'tiffn')
end

% ##########################################################################################################################################################################################################
% ##########################################################################################################################################################################################################
% ##########################################################################################################################################################################################################

%% Scale heat flux to account for larger surface area of the MPEX window:
close all

% Set dimensions:
% =========================================================================
% Proto-MPEX window and antenna dimensions:
% -------------------------------------------------------------------------
Proto.L1 = 4.75*1e-2; %[m] Distance from antenna to edge facing Target
Proto.L2 = 3.25*1e-2; %[m] Distance from antenna to edge facing Dump
Proto.La = 22*1e-2; %[m] antenna length
Proto.R  = 6.256*1e-2; %[m] window inner radius
Proto.Lw = 30*1e-2; % [m] Antenna length

% MPEX window and antenna dimensions:
% -------------------------------------------------------------------------
MPEX.windowType = {'DES' ,'MCS'};
MPEX.L1    = [8.636 ,7.216]*1e-2; %[m] 
MPEX.L2    = [9.016 ,7.597]*1e-2; %[m]
MPEX.La    = [25    ,25   ]*1e-2; %[m] Antenna length
MPEX.R     = [6.35  ,6.35 ]*1e-2; %[m] Window inner radius
MPEX.Lw    = [42.641,39.8 ]*1e-2; %[m] Antenna length

% Select Proto-MPEX data set to analyze:
% =========================================================================
% Select frame with steady-state heat flux:
fr = 40;
ii = 1;
Proto.q2D = rr_smooth(:,:,fr); % [Wm^-2]

% Coordinates associated with Proto.q2D
Proto.z_q2D = S{ii}.z_2D; %[m]
Proto.dz = Proto.z_q2D(2,1)-Proto.z_q2D(1,1);
Proto.phi_q2D = S{ii}.phi_2D; % [Rad]

% Identify the Proto-MPEX regions 1,2 and 3:
% =========================================================================
% The "z" coordinate associated with the Proto-MPEX window is S{ii}.z_2D
% It has 100 elements in angle and 50 elements along "z" and dim is 50*100
rng1 = find( Proto.z_q2D(:,1) < Proto.L2);
rng2 = find( Proto.z_q2D(:,1) > (Proto.Lw-Proto.L1));
rng3 = find( Proto.z_q2D(:,1) >= Proto.L2 & Proto.z_q2D(:,1) <= (Proto.Lw-Proto.L1) );

% Select the heat flux for regions 1,2 and 3:
% Region 1:
Proto.q1     = Proto.q2D(rng1,:);
Proto.z_q1   = Proto.z_q2D(rng1,:);
Proto.phi_q1 = Proto.phi_q2D(rng1,:);

% Region 2:
Proto.q2     = Proto.q2D(rng2,:);
Proto.z_q2   = Proto.z_q2D(rng2,:);
Proto.phi_q2 = Proto.phi_q2D(rng2,:);

% Region 3:
Proto.q3     = Proto.q2D(rng3,:);
Proto.z_q3   = Proto.z_q2D(rng3,:);
Proto.phi_q3 = Proto.phi_q2D(rng3,:);

% Calculate the mean heat flux in regions 1 and 2:
Proto.q1_mean = mean(Proto.q1,'all');
Proto.q2_mean = mean(Proto.q2,'all');

% MPEX heat flux calculations:
% =========================================================================
% Scale Proto.q3 to MPEX dimensions:
% -------------------------------------------------------------------------
for jj = 1:numel(MPEX.windowType)
    scaleFactor(jj) = (Proto.R/MPEX.R(jj))*(Proto.La/MPEX.La(jj));
    MPEX.q3{jj} = Proto.q3*scaleFactor(jj);
end

% Create MPEX "z" coordinate for region 3:
% -------------------------------------------------------------------------
for jj = 1:numel(MPEX.windowType)
    MPEX.dz(jj) = Proto.dz*MPEX.La(jj)/Proto.La;
    dum_z   = [MPEX.L1(jj) + MPEX.dz(jj) :MPEX.dz(jj):(MPEX.L1(jj) + MPEX.La(jj))]';
    dum_phi = [Proto.phi_q3(1,:)]';
    [MPEX.phi_q3{jj},MPEX.z_q3{jj}] = meshgrid(dum_phi,dum_z);
end

% L1 and L2 extension:
% -------------------------------------------------------------------------
extensionMethod = 1;
switch extensionMethod
    case 1
        % Calculate average heat flux on regions 1 and 2 based on Proto-MPEX data:
        for jj = 1:numel(MPEX.windowType)
            MPEX.q1_mean(jj) = Proto.q1_mean*scaleFactor(jj);
            MPEX.q2_mean(jj) = Proto.q2_mean*scaleFactor(jj);
        end

        % Create heat flux for region 1 and 2:
        for jj = 1:numel(MPEX.windowType)
            % Region 1:
            dum_z   = [0 :MPEX.dz(jj):(MPEX.L1(jj))]';
            dum_phi = [Proto.phi_q1(1,:)]';
            [MPEX.phi_q1{jj},MPEX.z_q1{jj}] = meshgrid(dum_phi,dum_z);
            % Heat flux:
            MPEX.q1{jj} = MPEX.q1_mean(jj)*ones(size(MPEX.z_q1{jj}));

            % Region 2:
            dum_z   = [(MPEX.L1(jj) + MPEX.La(jj) + MPEX.dz(jj)) :MPEX.dz(jj):(MPEX.Lw(jj))]';
            dum_phi = [Proto.phi_q2(1,:)]';
            [MPEX.phi_q2{jj},MPEX.z_q2{jj}] = meshgrid(dum_phi,dum_z);
            % Heat flux:
            MPEX.q2{jj} = MPEX.q2_mean(jj)*ones(size(MPEX.z_q2{jj}));
        end
    case 2
        % No other method has been developed yet
        
end

% Assemble the MPEX heat flux map:
% -------------------------------------------------------------------------
for jj = 1:numel(MPEX.windowType)
    MPEX.q2D{jj}     = [MPEX.q1{jj}     ;MPEX.q3{jj}     ;MPEX.q2{jj}     ];
    MPEX.z_q2D{jj}   = [MPEX.z_q1{jj}   ;MPEX.z_q3{jj}   ;MPEX.z_q2{jj}   ];
    MPEX.phi_q2D{jj} = [MPEX.phi_q1{jj} ;MPEX.phi_q3{jj} ;MPEX.phi_q2{jj} ];
end

% Check data and heat flux map stretching:
% =========================================================================
% Proto-MPEX data:
% -------------------------------------------------------------------------
figure('color','w')
subplot(1,2,2)
hold on
ax(1) = gca;
ht = title('Proto-MPEX data and regions 1,2 and 3');
    set(ht,'interpreter','latex','fontSize',12)
surf(Proto.phi_q1,Proto.z_q1,Proto.q1*1e-3,'LineStyle','none')
surf(Proto.phi_q2,Proto.z_q2,Proto.q2*1e-3,'LineStyle','none')
surf(Proto.phi_q3,Proto.z_q3,Proto.q3*1e-3,'LineStyle','none')
set(gca,'XTick',[0:45:360],'XDir','reverse')
xlim(ax,[0,360])
xlabel('\theta [deg]')
ylim(ax,[0,max(Proto.Lw)])
caxis([0,1200])
colormap('hot')
colorbar

% Annotations:
% Region 1:
htext = text(110,0.28,1200,'Region 1');
set(htext,'Interpreter','Latex','FontSize',14,'EdgeColor','w','Color','w')
% Region 2:
htext = text(110,0.015,1200,'Region 2');
set(htext,'Interpreter','Latex','FontSize',14,'EdgeColor','w','Color','w')
% Region 3:
htext = text(110,0.15,1200,'Region 3');
set(htext,'Interpreter','Latex','FontSize',14,'EdgeColor','w','Color','w')

subplot(1,2,1)
hold on
ax(2) = gca;
ht = title('Proto-MPEX data [kWm$^{-2}$]');
set(ht,'interpreter','latex','fontSize',12)
surf(Proto.phi_q2D,Proto.z_q2D,Proto.q2D*1e-3,'LineStyle','none')
set(gca,'XTick',[0:45:360],'XDir','reverse')
xlim(ax,[0,360])
xlabel('\theta [deg]')
ylim(ax,[0,max(Proto.Lw)])
ylabel('z [m]')
clear ax      
caxis([0,1200])
colormap('hot')
colorbar
set(gcf,'position',[172   171   965   420])

if 1
    figureName = ['Step_4c_ProtoMPEX_HeatFlux_175kW_CoupledPower'];
    saveas(gcf,figureName,'tiffn')
end

% MPEX data:
% -------------------------------------------------------------------------
jj = 1;
figure('color','w')
hold on
title(['MPEX data and region 3, ',num2str(MPEX.windowType{jj})])
surf(MPEX.phi_q1{jj},MPEX.z_q1{jj},MPEX.q1{jj}*1e-3,'LineStyle','none')
surf(MPEX.phi_q2{jj},MPEX.z_q2{jj},MPEX.q2{jj}*1e-3,'LineStyle','none')
surf(MPEX.phi_q3{jj},MPEX.z_q3{jj},MPEX.q3{jj}*1e-3,'LineStyle','none')

set(gca,'XTick',[0:45:360],'XDir','reverse')
xlim(gca,[0,360])
ylim(gca,[0,max(MPEX.Lw)])
caxis([0,1200])
colormap('hot')
colorbar

% Completed MPEX data:
% -------------------------------------------------------------------------
figure('color','w')
for jj = 1:numel(MPEX.windowType)
    subplot(1,2,jj)
    hold on
    ht = title(['Completed MPEX 2D heat flux [kWm$^{-2}$], ',num2str(MPEX.windowType{jj})]);
    set(ht,'interpreter','latex','fontSize',12)
    surf(MPEX.phi_q2D{jj},MPEX.z_q2D{jj},MPEX.q2D{jj}*1e-3,'LineStyle','none')
    set(gca,'XTick',[0:45:360],'XDir','reverse')
    xlim(gca,[0,360])
    xlabel('\theta [deg]')
    ylabel('z [m]')
    ylim(gca,[0,max(MPEX.Lw)])
    caxis([min(min(MPEX.q2D{jj})),max(max(MPEX.q2D{jj}))]*1e-3)     
    caxis([0,1200])
    colormap('hot')
    colorbar
end
set(gcf,'position',[172   171   965   420])

if 1
    figureName = ['Step_4c_MPEX_HeatFlux_DES_MCS_175kW_CoupledPower'];
    saveas(gcf,figureName,'tiffn')
end

%% Save steady-state frame for MPEX-limiter extrapolation to 175 kW Coupled

% Create a data set for the thermo-mechanical analsys:
% =========================================================================

for jj = 1:numel(MPEX.q2D)
    % 2D heat flux maps in Wm^-2:
    dataSet{jj}.q2D = MPEX.q2D{jj};

    % Coordinates of the data arrays:
    % Angular coordinate:
    dataSet{jj}.phi_q2D = MPEX.phi_q2D{jj}*pi/180;
    % Axial coordinate:
    dataSet{jj}.z_q2D   = MPEX.z_q2D{jj};

    % RF power value:
    dataSet{jj}.R_vacuum = Rv;
    dataSet{jj}.R_plasma = Rp;
    dataSet{jj}.RFpower_NET = mean_RF_NET
    dataSet{jj}.RFpower_Coupled = mean_RF_Coupled

    % Geometry of window:
    dataSet{jj}.RadiusWindow = MPEX.R(jj)
    dataSet{jj}.LengthWindow = MPEX.Lw(jj)
    dataSet{jj}.AntennaLength = MPEX.La(jj);
    dataSet{jj}.L1 = MPEX.L1(jj);
    dataSet{jj}.L2 = MPEX.L2(jj);

    dataSet{jj}.MPEX_sealType = MPEX.windowType{jj}
    dataSet{jj}.MagneticFluxMapping = 'MPEX-limiter';
   
    % Thermal parameters used for the data analysis:
    % This includes the directional emissivity of the AlN window at 12 deg
    % (0.55) and optics transmission factors due to ZnSe IR viewing window
    % (0.7)
    dataSet{jj}.thermalParams = S{1}.thermalParam;

    % Comments and metadata:
    dataSet{jj}.comment{1} = 'This heat flux 2D data is based on IR-based measurements performed in Proto-MPEX';
    dataSet{jj}.comment{2} = 'Heat flux has been linearly extrapolated from experiments performed at 80, 120 and 140 kW NET power in Proto-MPEX';
    dataSet{jj}.comment{3} = 'q2D is the heat flux map in [W]';
    dataSet{jj}.comment{4} = 'phi_q2D is the azimuthal angle coordinate in radians. Zero radians corresponds to 12 o clock in physical space';
    dataSet{jj}.comment{5} = 'z_q2D is axial coordinate (along window) in [m]';
end

fileName = ['Step_4c_MpexHeatFlux2D_175kW_CoupledPwr_created_',date];

% Save data:
saveData = 1;
if saveData
    disp('Saving 175 kW coupled power data ...')
    save(fileName,'dataSet')
    disp('Saving complete!')
end


% ##########################################################################################################################################################################################################
% ##########################################################################################################################################################################################################
% ##########################################################################################################################################################################################################

%% Use Proto-MPEX data to create MPEX ICH window heat flux
% Created 2020_12_18

close all

% Set dimensions:
% =========================================================================
% Proto-MPEX window and antenna dimensions:
% -------------------------------------------------------------------------
Proto_Hel.L1           = 4.75*1e-2; %[m] Distance from antenna to edge facing Target
Proto_Hel.L2           = 3.25*1e-2; %[m] Distance from antenna to edge facing Dump
Proto_Hel.La           = 22*1e-2; %[m] antenna length
Proto_Hel.R            = 6.256*1e-2; %[m] window inner radius
Proto_Hel.Lw           = 30*1e-2; % [m] Antenna length
Proto_Hel.PwrRfCoupled = 175*1e3; % [W] Coupled RF power
Proto_Hel.PwrWindow    = 40e3; % [W] Power coupled to window from 175 kW coupled RF power

% MPEX ICH window and antenna dimensions:
% -------------------------------------------------------------------------
MPEX_ICH.windowType   = {'DES Inlet' ,'DES Oulet'};
MPEX_ICH.L1           = [45.5        ,71.1       ]*1e-3; %[m] L1
MPEX_ICH.La           = [249.885     ,249.885    ]*1e-3; %[m] L3 Antenna length
MPEX_ICH.L2           = [68.9        ,61.8       ]*1e-3; %[m] L2
MPEX_ICH.Lw           = [364.5       ,382.9      ]*1e-3; %[m] Antenna length
MPEX_ICH.R            = [56          ,56         ]*0.5*1e-3; %[m] Window inner radius
MPEX_ICH.PwrRfCoupled = [155         ,155        ]*1e3; %[W] ICH RF coupled power
MPEX_ICH.PwrWindow    = MPEX_ICH.PwrRfCoupled*(Proto_Hel.PwrWindow/Proto_Hel.PwrRfCoupled);

% Select Proto-MPEX data set to analyze:
% =========================================================================
% Select frame with steady-state heat flux:
fr = 40;
ii = 1;
Proto_Hel.q2D = rr_smooth(:,:,fr); % [Wm^-2]

% Coordinates associated with Proto.q2D
Proto_Hel.z_q2D   = S{ii}.z_2D; %[m]
Proto_Hel.dz      = Proto_Hel.z_q2D(2,1)-Proto_Hel.z_q2D(1,1);
Proto_Hel.phi_q2D = S{ii}.phi_2D; % [Rad]

% Identify the Proto-MPEX regions 1,2 and 3:
% =========================================================================
% The "z" coordinate associated with the Proto-MPEX window is S{ii}.z_2D
% It has 100 elements in angle and 50 elements along "z" and dim is 50*100
rng2 = find( Proto_Hel.z_q2D(:,1) < Proto_Hel.L2);
rng1 = find( Proto_Hel.z_q2D(:,1) > (Proto_Hel.Lw-Proto_Hel.L1));
rng3 = find( Proto_Hel.z_q2D(:,1) >= Proto_Hel.L2 & Proto_Hel.z_q2D(:,1) <= (Proto_Hel.Lw-Proto_Hel.L1) );

% Select the heat flux for regions 1,2 and 3:
% Region 1:
Proto_Hel.q1     = Proto_Hel.q2D(rng1,:);
Proto_Hel.z_q1   = Proto_Hel.z_q2D(rng1,:);
Proto_Hel.phi_q1 = Proto_Hel.phi_q2D(rng1,:);

% Region 2:
Proto_Hel.q2     = Proto_Hel.q2D(rng2,:);
Proto_Hel.z_q2   = Proto_Hel.z_q2D(rng2,:);
Proto_Hel.phi_q2 = Proto_Hel.phi_q2D(rng2,:);

% Region 3:
Proto_Hel.q3     = Proto_Hel.q2D(rng3,:);
Proto_Hel.z_q3   = Proto_Hel.z_q2D(rng3,:);
Proto_Hel.phi_q3 = Proto_Hel.phi_q2D(rng3,:);

% Calculate the mean heat flux in regions 1 and 2:
Proto_Hel.q1_mean = mean(Proto_Hel.q1,'all');
Proto_Hel.q2_mean = mean(Proto_Hel.q2,'all');

% MPEX ICH window heat flux calculations:
% =========================================================================
% Scale Proto_Hel.q3 to MPEX ICH window dimensions:
% -------------------------------------------------------------------------
for jj = 1:numel(MPEX_ICH.windowType)
    
    % Calculate the scaling factor for the heat flux:
    % Change in length:
    lengthFactor    = Proto_Hel.Lw/MPEX_ICH.Lw(jj);
    % Change in radius:
    radiusFactor    = Proto_Hel.R/MPEX_ICH.R(jj);
    % Change in window heating rate:
    powerFactor     = MPEX_ICH.PwrWindow(jj)/Proto_Hel.PwrWindow;
    % Scaling factor that includes all changes in geometry and power:
    scaleFactor(jj) = radiusFactor*lengthFactor*powerFactor;
    
    % Apply scaling factor:
    MPEX_ICH.q3{jj} = Proto_Hel.q3*scaleFactor(jj);
    
    % Apply the effect of left-handed ICH antenna:
    % Mirror image of heat flux
    MPEX_ICH.q3{jj} = fliplr(MPEX_ICH.q3{jj});
    
    if jj == 2
        % Rotate heat flux since antenna is azimutally displaced by 90 deg:
        dphi = mean(diff(Proto_Hel.phi_q3(1,:)'));
        phi_offset = round(90/dphi);
        
        % Shift azimuthally:
        MPEX_ICH.q3{jj} = circshift(MPEX_ICH.q3{jj},[0,phi_offset]);
    end
    
end

% Create MPEX ICH "z" coordinate for region 3:
% -------------------------------------------------------------------------
for jj = 1:numel(MPEX_ICH.windowType)
    MPEX_ICH.dz(jj) = Proto_Hel.dz*(MPEX_ICH.La(jj)/Proto_Hel.La);
    dum_z   = [MPEX_ICH.L2(jj) + MPEX_ICH.dz(jj) :MPEX_ICH.dz(jj):(MPEX_ICH.L2(jj) + MPEX_ICH.La(jj))]';
    dum_phi = [Proto_Hel.phi_q3(1,:)]';
    [MPEX_ICH.phi_q3{jj},MPEX_ICH.z_q3{jj}] = meshgrid(dum_phi,dum_z);
end

% L1 and L2 extension:
% -------------------------------------------------------------------------
extensionMethod = 1;
switch extensionMethod
    case 1
              
        % Calculate average heat flux on regions 1 and 2 based on Proto-MPEX data:
        for jj = 1:numel(MPEX.windowType)
            % Region 1:
            MPEX_ICH.q1_mean(jj) = Proto_Hel.q1_mean*scaleFactor(jj);
            
            % Region 2:         
            MPEX_ICH.q2_mean(jj) = Proto_Hel.q2_mean*scaleFactor(jj);
        end

        % Create heat flux for region 1 and 2:
        for jj = 1:numel(MPEX_ICH.windowType)
            % Region 2:
            dum_z   = [0 :MPEX_ICH.dz(jj):(MPEX_ICH.L2(jj))]';
            dum_phi = [Proto_Hel.phi_q2(1,:)]';
            [MPEX_ICH.phi_q2{jj},MPEX_ICH.z_q2{jj}] = meshgrid(dum_phi,dum_z);
            % Heat flux:
            MPEX_ICH.q2{jj} = MPEX_ICH.q2_mean(jj)*ones(size(MPEX_ICH.z_q2{jj}));

            % Region 1:
            dum_z   = [(MPEX_ICH.L2(jj) + MPEX_ICH.La(jj) + MPEX_ICH.dz(jj)) :MPEX_ICH.dz(jj):(MPEX_ICH.Lw(jj))]';
            dum_phi = [Proto_Hel.phi_q1(1,:)]';
            [MPEX_ICH.phi_q1{jj},MPEX_ICH.z_q1{jj}] = meshgrid(dum_phi,dum_z);
            % Heat flux:
            MPEX_ICH.q1{jj} = MPEX_ICH.q1_mean(jj)*ones(size(MPEX_ICH.z_q1{jj}));
        end
    case 2
        % No other method has been developed yet
        
end

% Assemble the MPEX heat flux map:
% -------------------------------------------------------------------------
for jj = 1:numel(MPEX_ICH.windowType)
    MPEX_ICH.q2D{jj}     = [MPEX_ICH.q2{jj}     ;MPEX_ICH.q3{jj}     ;MPEX_ICH.q1{jj}     ];
    MPEX_ICH.z_q2D{jj}   = [MPEX_ICH.z_q2{jj}   ;MPEX_ICH.z_q3{jj}   ;MPEX_ICH.z_q1{jj}   ];
    MPEX_ICH.phi_q2D{jj} = [MPEX_ICH.phi_q2{jj} ;MPEX_ICH.phi_q3{jj} ;MPEX_ICH.phi_q1{jj} ];
end

% Check data and heat flux map stretching:
% =========================================================================
% Proto-MPEX data:
% -------------------------------------------------------------------------
figure('color','w')
subplot(1,2,2)
hold on
ax(1) = gca;
ht = title('Proto-MPEX data and regions 1,2 and 3');
    set(ht,'interpreter','latex','fontSize',12)
surf(Proto_Hel.phi_q1,Proto_Hel.z_q1,Proto_Hel.q1*1e-3,'LineStyle','none')
surf(Proto_Hel.phi_q2,Proto_Hel.z_q2,Proto_Hel.q2*1e-3,'LineStyle','none')
surf(Proto_Hel.phi_q3,Proto_Hel.z_q3,Proto_Hel.q3*1e-3,'LineStyle','none')
set(gca,'XTick',[0:45:360],'XDir','reverse')
xlim(ax,[0,360])
xlabel('\theta [deg]')
ylim(ax,[0,max(Proto_Hel.Lw)])
caxis([0,1200])
colormap('hot')
colorbar

% Annotations:
% Region 1:
htext = text(110,0.28,1200,'Region 1');
set(htext,'Interpreter','Latex','FontSize',14,'EdgeColor','w','Color','w')
% Region 2:
htext = text(110,0.015,1200,'Region 2');
set(htext,'Interpreter','Latex','FontSize',14,'EdgeColor','w','Color','w')
% Region 3:
htext = text(110,0.15,1200,'Region 3');
set(htext,'Interpreter','Latex','FontSize',14,'EdgeColor','w','Color','w')

subplot(1,2,1)
hold on
ax(2) = gca;
ht = title('Proto-MPEX data [kWm$^{-2}$]');
set(ht,'interpreter','latex','fontSize',12)
surf(Proto_Hel.phi_q2D,Proto_Hel.z_q2D,Proto_Hel.q2D*1e-3,'LineStyle','none')
set(gca,'XTick',[0:45:360],'XDir','reverse')
xlim(ax,[0,360])
xlabel('\theta [deg]')
ylim(ax,[0,max(Proto_Hel.Lw)])
ylabel('z [m]')
clear ax      
caxis([0,1200])
colormap('hot')
colorbar
set(gcf,'position',[172   171   965   420])

if 1
    figureName = ['Step_4c_ProtoMPEX_Helicon_HeatFlux_175kW_CoupledPower'];
    saveas(gcf,figureName,'tiffn')
end

% MPEX data:
% -------------------------------------------------------------------------
jj = 1;
figure('color','w')
hold on
title(['MPEX ICH data and region 3, ',num2str(MPEX_ICH.windowType{jj})])
surf(MPEX_ICH.phi_q1{jj},MPEX_ICH.z_q1{jj},MPEX_ICH.q1{jj}*1e-3,'LineStyle','none')
surf(MPEX_ICH.phi_q2{jj},MPEX_ICH.z_q2{jj},MPEX_ICH.q2{jj}*1e-3,'LineStyle','none')
surf(MPEX_ICH.phi_q3{jj},MPEX_ICH.z_q3{jj},MPEX_ICH.q3{jj}*1e-3,'LineStyle','none')

set(gca,'XTick',[0:45:360],'XDir','reverse')
xlim(gca,[0,360])
ylim(gca,[0,max(MPEX_ICH.Lw)])
caxis([0,1200])
colormap('hot')
colorbar

% Completed MPEX data:
% -------------------------------------------------------------------------
figure('color','w')
for jj = 1:numel(MPEX_ICH.windowType)
    subplot(1,2,jj)
    hold on
    ht = title(['MPEX ICH heat flux [kWm$^{-2}$], ',num2str(MPEX_ICH.windowType{jj})]);
    set(ht,'interpreter','latex','fontSize',12)
    surf(MPEX_ICH.phi_q2D{jj},MPEX_ICH.z_q2D{jj},MPEX_ICH.q2D{jj}*1e-3,'LineStyle','none')
    set(gca,'XTick',[0:45:360],'XDir','reverse')
    xlim(gca,[0,360])
    xlabel('\theta [deg]')
    ylabel('z [m]')
    ylim(gca,[0,max(MPEX_ICH.Lw)])
    caxis([min(min(MPEX_ICH.q2D{jj})),max(max(MPEX_ICH.q2D{jj}))]*1e-3)     
    caxis([0,1800])
    colormap('hot')
    colorbar
end
set(gcf,'position',[172   171   965   420])

if 1
    figureName = ['Step_4c_MPEX_ICH HeatFlux_DES_MCS_175kW_CoupledPower'];
    saveas(gcf,figureName,'tiffn')
end

%% Save steady-state frame for MPEX ICH window extrapolation to 155 kW Coupled

% Create a data set for the thermo-mechanical analsys:
% =========================================================================

for jj = 1:numel(MPEX_ICH.q2D)
    % ICH 2D heat flux maps in Wm^-2:
    dataSet_ICH{jj}.q2D = MPEX_ICH.q2D{jj};

    % Coordinates of the data arrays:
    % Angular coordinate:
    dataSet_ICH{jj}.phi_q2D = MPEX_ICH.phi_q2D{jj}*pi/180;
    % Axial coordinate:
    dataSet_ICH{jj}.z_q2D   = MPEX_ICH.z_q2D{jj};

    % RF power value:
    dataSet_ICH{jj}.RFpower_Coupled = MPEX_ICH.PwrRfCoupled(jj);

    % Geometry of window:
    dataSet_ICH{jj}.RadiusWindow = MPEX_ICH.R(jj)
    dataSet_ICH{jj}.LengthWindow = MPEX_ICH.Lw(jj)
    dataSet_ICH{jj}.AntennaLength = MPEX_ICH.La(jj);
    dataSet_ICH{jj}.L1 = MPEX_ICH.L1(jj);
    dataSet_ICH{jj}.L2 = MPEX_ICH.L2(jj);

    dataSet_ICH{jj}.MPEX_sealType = MPEX_ICH.windowType{jj}
    dataSet_ICH{jj}.MagneticFluxMapping = 'MPEX-limiter';
   
    % Thermal parameters used for the data analysis:
    % This includes the directional emissivity of the AlN window at 12 deg
    % (0.55) and optics transmission factors due to ZnSe IR viewing window
    % (0.7)
    
    dataSet_ICH{jj}.thermalParams = S{1}.thermalParam;

    % Comments and metadata:
    dataSet_ICH{jj}.comment{1} = 'This heat flux 2D data is based on Helicon window IR-based measurements performed in Proto-MPEX';
    dataSet_ICH{jj}.comment{2} = 'Heat flux has been linearly extrapolated from experiments performed at 80, 120 and 140 kW Helicon NET power in Proto-MPEX';
    dataSet_ICH{jj}.comment{3} = 'q2D is the ICH window heat flux map in [W]';
    dataSet_ICH{jj}.comment{4} = 'phi_q2D is the azimuthal angle coordinate in radians. Zero radians corresponds to 12 o clock in physical space';
    dataSet_ICH{jj}.comment{5} = 'z_q2D is axial coordinate (along window) in [m]';
end

fileName = ['Step_4c_MPEX_ICH_HeatFlux2D_155kW_CoupledPwr_created_',date];

% Save data:
saveData = 1;
if saveData
    disp('Saving 155 kW coupled power ICH heat flux data ...')
    save(fileName,'dataSet_ICH')
    disp('Saving complete!')
end
