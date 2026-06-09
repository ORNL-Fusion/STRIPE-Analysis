% The following script collects data from July to September 2016. It
% focuses on on-axis data (y = 0,+1) in spool 9.5.

%==========================================================================
% Changes relative to V1:
% Sept 21st 2016
% We have incorporated the use of Rick's DLP code in addition to my
% DLP_fit_v2 code
% Just DLP_fit_v2 takes 130 seconds to compute
% Just RIck's code takes 
%==========================================================================

close all
clear all

AnalysisType = 1;
% 1: Calculates DLP characteristics one by on XP runs
% 2: Calculates all shots at once with my DLP code and Ricks's
% 4: Loads the results from case 2

switch AnalysisType
    case 1
% Use Segmented analysis:
%% April 27th:
% shotlist = [8411,8412]; % x2 attenuators on V and I
% shotlist = [8413,8414,8415]; % x1 attenuators on V and I

%% May 6th:
% shotlist = [8515,8516,8517,8522,8523];
%% May 11th:
%  shotlist = [8644,8645,8646,8647,8648,8649,8650,8651,8652,8653,8654,8655,8656];
%% May 19th:
% shotlist = [8846,8847,8848,8851,8856,8857,8883,8884,8885]; % 21 sec
%% June 1st:
% shotlist = [8937,8938,8946];
%% June 15th:
%shotlist = 8982;
%% June 17th:
% shotlist = 9044;
%% June 29th:
% shotlist = [9141,9142,9143];
%% July 1st:
%shotlist = [9179,9180,9181,9182,9183,9184,9185,9186,9187,9202,9203,9204];
%% July 6th:
% shotlist = [9247,9248];
%% July 13th:
% shotlist = [9437,9438,9442,9443,9444];
%% July 14th:
% shotlist = [9447,9449,9450,9452];
%% July 26th:
% shotlist = 9500 + [15:22,24:28,30:35];
%% July 28th:
% shotlist = 9500 + [55,60,65:68,70,74:80,86:87,89:91,94,102];
%% Aug 2nd:
% shotlist = 9600 + [8:10,37,39,43,46,47];
% 9608, there is a 20% difference between Rick and my code. it appears to
% be related to the probe characteristic been asymmetric
%% Aug 3rd:
shotlist = 9600 + [79:82,85,86];
%% Aug 12th:
% shotlist = 9700 + [73,77:86,94,95,101:103,105,106,110,117,121,133,142,151:153];
%% Sept 16th:
% shotlist = 1e4 + 300 + [46:50,52:56,70,72,79,81,82,85,86,92,94,95,99,100,101,103:114,116:119];
% shotlist  = 1e4 + 300 + [103:114,116:119]; 

% Some of these shots display a great difference between Rick's and my
% code. the main reason seems to be the asymmetry in the DLP characteristic
    case 2
% Analyze everything at once:
% April 27th:
sn{1} = [8411,8412]; % x2 attenuators on V and I
sn1b = [8413,8414,8415]; % x1 attenuators on V and I
% May 6th:
sn{2} = [8515,8516,8517,8522,8523];
% May 11th:
sn{3} = [8644,8645,8646,8647,8648,8649,8650,8651,8652,8653,8654,8655,8656];
% May 19th:
sn{4} = [8846,8847,8848,8851,8856,8857,8883,8884,8885];
% June 1st:
sn{5} = [8937,8938,8946];
% June 15th:
sn{6} = 8982;
% June 17th:
sn{7} = 9044;
% June 29th:
sn{8} = [9141,9142,9143];
% July 1st:
sn{9} = [9179,9180,9181,9182,9183,9184,9185,9186,9187,9202,9203,9204];
% July 6th:
sn{10} = [9247,9248];
% July 13th:
sn{11} = [9437,9438,9442,9443,9444];
% July 14th:
sn{12} = [9447,9449,9450,9452];
% July 26th:
sn{13} = 9500 + [15:22,24:28,30:35];
% July 28th:
sn{14} = 9500 + [55,60,65:68,70,74:80,86:87,89:91,94,102];
% Aug 2nd:
sn{15} = 9600 + [8:10,37,39,43,46,47];
% Aug 3rd:
sn{16} = 9600 + [79:82,85,86];
% Aug 12th:
sn{17} = 9700 + [73,77:86,94,95,101:103,105,106,110,117,121,133,142,151:153];
% Sept 16th:
sn{18} = 1e4 + 300 + [46:50,52:56,70,72,79,81,82,85,86,92,94,95,99,100,101,103:114,116:119];

shotlist = [sn{1},sn{2},sn{3},sn{4},sn{5},sn{6},sn{7},sn{8},sn{9},sn{10},sn{11},sn{12},...
    sn{13},sn{14},sn{15},sn{16},sn{17},sn{18}];
    case 3
% May 6th 2016, Spool 6.5 radial scan, use 4:1 attenuator on Current         
        shotlist = [8602:1:8622];
    case 4
        load('file2.mat')
end

if AnalysisType ~= 4
% Acquiring Ne and Te data
Stem = '\MPEX::TOP.';
Branch = 'MACHOPS1:';
RootAddress = [Stem,Branch];
DataAddress{1} = [RootAddress,'LP_V_RAMP']; % V
DataAddress{2} = [RootAddress,'TARGET_LP']; % I

Config.tStart = 4.16; % [s]
Config.tEnd = 4.32;
Config.V_Att = 2;  % Output voltage of DLP box (Voltage) = V_Att*Digitized data 
if AnalysisType == 3
    Config.I_Att = 4;  % Output voltage of DLP box (Current) = I_att*Digitized data
else
    Config.I_Att = 2;  % Output voltage of DLP box (Current) = I_att*Digitized data

end
Config.V_cal = [12.05,0.205];   % Voltage output of DLP = V_cal(1)*Output voltage of DLP box + V_cal(2)
Config.I_cal = [-142.5, 0.015]; % Current output of DLP = (I_cal(2) + Output voltage of DLP box)/Ical(1)
Config.FilterDataInput = 1; % Filter input data with savitsky Golay filter order 3 frame size 7
Config.TimeMode = 2; % Effective time of sweep, (1) start of ramp or (2) mean time of ramp
Config.AMU = 2; % Ion mass in AMU
Config.L_tip = 1.2/1000;
Config.D_tip = 0.254/1000; % [m]
Config.FitFunction = 2; 
Config.AreaType = 2; 
% 1 - Cylindrical + cap
% 2 - Cylindrical only
% 3 - x2 Projected only

shotlist_Main = shotlist;

% ########################################################################
% JUAN CANESES,DATA ANALYSIS:
% ########################################################################
[Ni,Te,time,Ifit,Ip,Vp,tm,Vsweep,Isweep] = DLP_fit_V2(Config,shotlist_Main,DataAddress);

% ########################################################################
% RICK GOULDING, DATA ANALYSIS:
% ########################################################################
% the output of Rick's code is ne [cm-3], te [V], tim(:,1) [s] 
BadShots = [8937,9522,9518,9521,9528,9780,9803,10394,10411];
for s = 1:length(shotlist)
        shotlist = shotlist_Main(s);
    if ~sum(shotlist == BadShots)
        plot_ne_te_spool9p5_050616
        Te_rg{s} = te;
        Ne_rg{s} = ne*1e6; % Convert from cm-3 to m-3
        t_rg{s} = tim;
    else
        Ne_rg{s} = NaN;
        Te_rg{s} = NaN;
        t_rg{s} = NaN;
    end
end
end
%%
t = t_zero(shotlist_Main);
tGather = 4.28;

figure; 
subplot(2,2,1); hold on
for s = 1:length(shotlist_Main)
    h(s) = plot(time{s},Ni{s},'lineWidth',1);
    c = get(h(s),'color');
    plot(t_rg{s},Ne_rg{s},'color',c,'LineWidth',1,'LineStyle',':')
    
    L{s} = [num2str(shotlist_Main(s)),' ,t=',num2str(t{s}(10:14))];
    
    if AnalysisType ~= 4
    % Gather data into table
    nt = find(time{s}>=tGather); nt = nt(1:3);
    T(s,2) = mean(Ni{s}(nt)); 
    T(s,3) = mean(Te{s}(nt));
    T(s,1) = shotlist_Main(s);
    
    if isnan(t_rg{s})
        T_rg(s,2) = NaN; 
        T_rg(s,3) = NaN;
        T_rg(s,1) = shotlist_Main(s);
    else
        nt = find(t_rg{s}>=tGather); nt = nt(1:3);
        T_rg(s,2) = mean(Ne_rg{s}(nt)); 
        T_rg(s,3) = mean(Te_rg{s}(nt));
        T_rg(s,1) = shotlist_Main(s);
    end
    end
end
legend(h,L,'location','NorthWest')
title('$ n_e $ $ [m^{-3}] $','interpreter','Latex','FontSize',13,'Rotation',0)
%set(gca,'PlotBoxAspectRatio',[1 1 1])
ylim([0,6e19])
xlim([4.15,4.32])

subplot(2,2,2); hold on
for s = 1:length(shotlist_Main)
    plot(time{s},Te{s},'lineWidth',1)
    c = get(h(s),'color');
    plot(t_rg{s},Te_rg{s},'color',c,'LineWidth',1,'LineStyle',':')    
end
title('$ T_e $ $ [eV] $','interpreter','Latex','FontSize',13,'Rotation',0)
ylim([0,10])
xlim([4.15,4.32])

subplot(2,2,3); hold on
for s = 1:length(shotlist_Main)
    plot(time{s},e_c.*Ni{s}.*Te{s},'lineWidth',1)
    c = get(h(s),'color');
    plot(t_rg{s},e_c*Te_rg{s}.*Ne_rg{s},'color',c,'LineWidth',1,'LineStyle',':')
end
ylim([0,20])
xlim([4.15,4.32])
title('$ P_e $ $ [Pa] $','interpreter','Latex','FontSize',13,'Rotation',0)

set(findobj('-Property','YTick'),'box','on')
set(gcf,'color','w')

figure; plot(T(:,2),T(:,3),'ko'); hold on; plot(T_rg(:,2),T_rg(:,3),'ro')
ylim([0,10]);
xlim([0,7e19]);

if 1
%     s = 2;
    figure; hold on
    plot(tm{s},Vp{s})
    plot(tm{s},Ip{s}*1000)
    ylim([-100,100])
    xlim([4.15,4.32])
    grid on
end

if 1
    
     figure;
    for c = 1:25;
        subplot(5,5,c); hold on
        plot(Vsweep{s}{c},Isweep{s}{c}*1e3,'k')
        plot(Vsweep{s}{c},Ifit{s}{c}*1e3,'r')
        ht = title(['T_e: ',num2str(Te{s}(c)),'t: ',num2str(time{s}(c))]);
        set(ht,'FontSize',5)
    end
    
    figure;
        for c = 26:50;
        subplot(5,5,c-25); hold on
        plot(Vsweep{s}{c},Isweep{s}{c}*1e3,'k')
        plot(Vsweep{s}{c},Ifit{s}{c}*1e3,'r')
        ht = title(['T_e: ',num2str(Te{s}(c)),'t: ',num2str(time{s}(c))]);
        set(ht,'FontSize',5)
        end
        
    figure;
        for c = 51:75;
        subplot(5,5,c-50); hold on
        plot(Vsweep{s}{c},Isweep{s}{c}*1e3,'k')
        plot(Vsweep{s}{c},Ifit{s}{c}*1e3,'r')
        ht = title(['T_e: ',num2str(Te{s}(c)),'t: ',num2str(time{s}(c))]);
        set(ht,'FontSize',5)
        end  
    
end