% In this script I run two DLP fitting codes, one that I have written and
% Rick Goulding's code as of Aug 22nd 2016.

% Script written by Juan Caneses Aug 26th 2016 unless otherwise stated

close all
clear all

shotlist = 10000 + [101];

% ########################################################################
%  PROBE VALUES:
% ########################################################################
Spool = 9.5;
AMU = 2; % D2 gas
D_tip = 0.254/1000; % [m]
rp = D_tip/2; % [m]

switch Spool
    case 6.5
        L_tip = 1.3/1000;  % [m]
    case 9.5
        L_tip = 1.2/1000;  % [m]
    case 4.5
        L_tip = 1.2/1000;  % [m]
    case 10.5
        L_tip = 1.2/1000;  % [m]
end

% ########################################################################
% MDS+ ADDRESS:
% ########################################################################\
Address = 1;
switch Address
    case 1
        Stem = '\MPEX::TOP.';
        Branch = 'MACHOPS1:';
        RootAddress = [Stem,Branch];
        DataAddress{1} = [RootAddress,'LP_V_RAMP']; % V
        DataAddress{2} = [RootAddress,'TARGET_LP']; % I
        Data = [RootAddress,'PWR_28GHz']; % ECH
             
        mdsconnect('mpexserver');
        [SHOT,~]=mdsopen('MPEX',shotlist); % see 8636
        [Gyrotron,~]  = mdsvalue(Data);
        [t_Gyrotron,~] = mdsvalue(['DIM_OF(',Data,')']);

    case 2
        Stem = '\MPEX::TOP.';
        Branch = 'MPEX1:';
        RootAddress = [Stem,Branch];
        DataAddress{1} = [RootAddress,'RF_REF_PWR']; % V
        DataAddress{2} = [RootAddress,'TARGET_LP']; % I
end
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
Config.AreaType = 1; % Cylindrical + cap
[Ni_a,Te_a,time_a,Ifit,I_dlp,V_dlp,t_dlp,Vsweep,Isweep] = DLP_fit_V1(Config,shotlist,DataAddress);

Config.AreaType = 2; % Cylindrical no cap
[Ni_b,Te_b,time_b,~,~,~,~,~,~] = DLP_fit_V1(Config,shotlist,DataAddress);

% ########################################################################
% RICK GOULDING, DATA ANALYSIS:
% ########################################################################
% the output of Rick's code is ne [cm-3], te [V], tim(:,1) [s] 
plot_ne_te_spool9p5_050616
ne = ne*1e6; % Convert from cm-3 to m-3

%% Plot data:
%close all

C = {'Black','Red','Blue','Green'};
PlotTimeStart = 4.16; 
PlotTimeEnd = 4.32;


for s = 1:length(shotlist)
    
figure; 
subplot(2,2,1); hold on

hJC_a(s) = plot(time_a{s},Ni_a{s},'k');
hJC_b(s) = plot(time_b{s},Ni_b{s},'k','LineStyle',':');
hRG(s) = plot(tim(:,1),ne(:,s),'r','LineWidth',1,'LineStyle','-')
legend([hJC_a(s),hJC_b(s),hRG(s)],{[num2str(shotlist(s)'),...
    ' Caneses "a"'],[num2str(shotlist(s)'),' Caneses "b"'],[num2str(shotlist(s)'),' Goulding']})
     
title('$ n_e $ $ [m^{-3}] $','interpreter','Latex','FontSize',13,'Rotation',0)
xlim([PlotTimeStart,PlotTimeEnd])
ylim([0,6e19])
ax(1) = gca;

subplot(2,2,2); hold on
plot(time_a{s},Te_a{s},'k')
plot(time_b{s},Te_b{s},'k','LineStyle',':');
plot(tim(:,1),te(:,s),'r','LineWidth',1,'LineStyle','-')
plot(t_Gyrotron,Gyrotron,'g');
title('$ T_e $ $ [eV] $','interpreter','Latex','FontSize',13,'Rotation',0)
ylim([0,8])
xlim([PlotTimeStart,PlotTimeEnd])
ax(2) = gca;

subplot(2,2,3); hold on
plot(time_a{s},Ni_a{s}.*Te_a{s}*e_c,'k')
plot(time_b{s},Ni_b{s}.*Te_b{s}*e_c,'k','LineStyle',':')
plot(tim(:,1),ne(:,s).*te(:,s)*e_c,'r','LineWidth',1,'LineStyle','-')
title('$ P_e $ $ [Pa] $','interpreter','Latex','FontSize',13,'Rotation',0)
set(gcf,'position',[417 106 574 503])
ylim([0,20])
xlim([PlotTimeStart,PlotTimeEnd])
ax(3) = gca;

set(gcf,'color','w')
set(ax,'Box','on')

end

%% Fit check:
if 1
    figure;     hold on
    for s = 1:length(shotlist);
    plot(t_dlp{s},V_dlp{s})
    plot(t_dlp{s},I_dlp{s}*1000)
    end
    grid on
    ylim([-100,100])
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
