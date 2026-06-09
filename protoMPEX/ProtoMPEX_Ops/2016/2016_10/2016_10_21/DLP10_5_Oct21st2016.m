% In this script I run two DLP fitting codes, one that I have written and
% Rick Goulding's code as of Aug 22nd 2016.

% Script written by Juan Caneses Oct 27th 2016 unless otherwise stated

close all
clear all

% ########################################################################
%  PROBE VALUES:
% ########################################################################
% On shot 11098 and 99, DLP 10.5 was inserted. the day before, we had ran the 2s
% 100 kW pulses and damaged DLP 10.5. upon inspection by John C, the probe
% tip length is about 0.5 mm in length with a 25% uncertainty

% DLP 6.5 1.01 mm
% DLP 9.5 1.2 mm

shotlist = 11000 + [98 ,  99];
%shotlist = 11000 + [94 , 123];
%shotlist = 11000 + [95,97]; % before TS

% shotlist = 11000 + [123,125]; After TS
%shotlist = 11000 + [100,101];


L_tip = 0.5/1000; % probe tip length in mm, 
%L_tip = 1.2/1000; % probe tip length in mm, 
%L_tip = 1.01/1000; % probe tip length in mm, 

fctr  =  1.2./(L_tip*1000);

AMU = 2; % D2 gas
D_tip = 0.254/1000; % [m]
rp = D_tip/2; % [m]

% ########################################################################
% MDS+ ADDRESS:
% ########################################################################\
Stem = '\MPEX::TOP.';
Branch = 'MACHOPS1:';
RootAddress = [Stem,Branch];
DataAddress{1} = [RootAddress,'LP_V_RAMP']; % V
DataAddress{2} = [RootAddress,'TARGET_LP']; % I

% ########################################################################
% DLP fitting procedure condition configuration variables
% ########################################################################
Config.tStart = 4.15; % [s]
Config.tEnd = 4.32;
Config.V_Att = 2;  % Output voltage of DLP box (Voltage) = V_Att*Digitized data 
Config.I_Att = 2;  % Output voltage of DLP box (Current) = I_att*Digitized data
Config.V_cal = [12.05,0.205];   % Voltage output of DLP = V_cal(1)*Output voltage of DLP box + V_cal(2)
Config.I_cal = [-142.5, 0.015]; % Current output of DLP = (I_cal(2) + Output voltage of DLP box)/Ical(1)
Config.FilterDataInput = 1; % Filter input data with savitsky Golay filter order 3 frame size 7
Config.TimeMode = 2; % Effective time of sweep, (1) start of ramp or (2) mean time of ramp
Config.AMU = 2; % Ion mass in AMU
Config.L_tip = L_tip;
Config.D_tip = D_tip;
Config.FitFunction = 2; 

% ########################################################################
% JUAN CANESES,DATA ANALYSIS:
% ########################################################################
if 1
Config.AreaType = 1; % Cylindrical + cap
[Ni_a,Te_a,time_a,Ifit,I_dlp,V_dlp,t_dlp,Vsweep,Isweep]...
    = DLP_fit_V5(Config,shotlist,DataAddress);
else
Config.AreaType = 2; % Cylindrical no cap
[Ni_a,Te_a,time_a,Ifit,I_dlp,V_dlp,t_dlp,Vsweep,Isweep]...
    = DLP_fit_V5(Config,shotlist,DataAddress);
end

% ########################################################################
% RICK GOULDING, DATA ANALYSIS:
% ########################################################################
% the output of Rick's code is ne [cm-3], te [V], tim(:,1) [s] 
plot_ne_te_spool9p5_050616

% To correct for the smaller probe area (0.5 mm and not 0.12) we use a
% multiplication factor
ne = ne*fctr*1e6; % Convert from cm-3 to m-3


%% Plot data:
close all

C = {'Black','Red','Blue','Green'};
PlotTimeStart = 4.18; 
PlotTimeEnd = 4.32;

figure; 
k = 1;
for s = 1:length(shotlist)
    
subplot(2,2,k); hold on
hRG(s) = plot(tim(:,s),ne(:,s),'r','LineWidth',2,'LineStyle','-')

plot(time_a{s},Ni_a{s}{1},'k:','linewidth',1);
hJC_a(s) = plot(time_a{s},0.5*(Ni_a{s}{1}+Ni_a{s}{2}),'k','linewidth',2);
plot(time_a{s},Ni_a{s}{2},'k:','linewidth',1);

L = legend([hJC_a(s),hRG(s)],{[num2str(shotlist(s)'),...
    ' Caneses '],[num2str(shotlist(s)'),' Goulding']},...
    'location','NorthWest')
set(L,'Fontsize',7)
     
title('$ n_e $ $ [m^{-3}] $','interpreter','Latex','FontSize',13,'Rotation',0)
xlim([PlotTimeStart,PlotTimeEnd])
ylim([0,8e19])
ax(1) = gca;

k = k+1;
subplot(2,2,k); hold on
plot(time_a{s},Te_a{s},'k','linewidth',2)
plot(tim(:,s),te(:,s),'r','LineWidth',2,'LineStyle','-')
title('$ T_e $ $ [eV] $','interpreter','Latex','FontSize',13,'Rotation',0)
ylim([0,4])
xlim([PlotTimeStart,PlotTimeEnd])
ax(2) = gca;

set(gcf,'color','w')
set(ax,'Box','on')
k = k + 1;
end

%% Fit check:
if 1
    figure;     
    hold on
    for s = 1:length(shotlist);
        plot(t_dlp{s},V_dlp{s})
        plot(t_dlp{s},I_dlp{s}*1000)
    end
    grid on
    ylim([-100,100])
    set(gcf,'color','w')
    xlim([PlotTimeStart,PlotTimeEnd])
    set(gca,'PlotBoxAspectRatio',[1 0.8 1])
    
        figure;     
    hold on
    for s = 1:length(shotlist);
        plot(t_dlp{s},I_dlp{s}*1000)
    end
    grid on
    ylim([-25,25])
    set(gcf,'color','w')
    xlim([PlotTimeStart,PlotTimeEnd])
    set(gca,'PlotBoxAspectRatio',[1 0.8 1])
end

if 1
    
     figure;
    for c = 1:25;
        subplot(5,5,c); hold on
        plot(Vsweep{s}{c},Isweep{s}{c}*1e3,'k')
        plot(Vsweep{s}{c},Ifit{s}{c}*1e3,'r')
        ht = title(['T_e: ',num2str(Te_a{s}(c)),'t: ',num2str(time_a{s}(c))]);
        set(ht,'FontSize',5)
    end
    
    figure;
        for c = 26:50;
        subplot(5,5,c-25); hold on
        plot(Vsweep{s}{c},Isweep{s}{c}*1e3,'k')
        plot(Vsweep{s}{c},Ifit{s}{c}*1e3,'r')
        ht = title(['T_e: ',num2str(Te_a{s}(c)),'t: ',num2str(time_a{s}(c))]);
        set(ht,'FontSize',5)
        end
        
    figure;
        for c = 51:75;
        subplot(5,5,c-50); hold on
        plot(Vsweep{s}{c},Isweep{s}{c}*1e3,'k')
        plot(Vsweep{s}{c},Ifit{s}{c}*1e3,'r')
        ht = title(['T_e: ',num2str(Te_a{s}(c)),'t: ',num2str(time_a{s}(c))]);
        set(ht,'FontSize',5)
        end  
    
end
