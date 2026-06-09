% The following script collects data from July to September 2016. It
% focuses on on-axis data (y = 0,+1) in spool 9.5.

close all
clear all

%% April 27th:
%shotlist = [8411,8412]; % x2 attenuators on V and I
%shotlist = [8413,8414,8415]; % x1 attenuators on V and I

%% May 6th:
%shotlist = [8515,8516,8517,8522,8523];
%% May 11th:
% shotlist = [8644,8645,8646,8647,8648,8649,8650,8651,8652,8653,8654,8655,8656];
%% May 19th:
%shotlist = [8846,8847,8848,8851,8856,8857,8883,8884,8885];
%% June 1st:
% shotlist = [8937,8938,8946];
%% June 15th:
%shotlist = 8982;
%% June 17th:
%shotlist = 9044;
%% June 29th:
%shotlist = [9141,9142,9143];
%% July 1st:
%shotlist = [9179,9180,9181,9182,9183,9184,9185,9186,9187,9202,9203,9204];
%% July 6th:
%shotlist = [9247,9248];
%% July 13th:
%shotlist = [9437,9438,9442,9443,9444];
%% July 14th:
%shotlist = [9447,9449,9450,9452];
%% July 26th:
%shotlist = 9500 + [15:22,24:28,30:35];
%% July 28th:
% shotlist = 9500 + [55,60,65:68,70,74:80,86:87,89:91,94,102];
%% Aug 2nd:
% shotlist = 9600 + [8:10,37,39,43,46,47];
%% Aug 3rd:
% shotlist = 9600 + [79:82,85,86];
%% Aug 12th:
% shotlist = 9700 + [73,77:86,94,95,101:103,105,106,110,117,121,133,142,151:153];
%% Sept 16th:
%shotlist = 1e4 + 300 + [46:50,52:56,70,72,79,81,82,85,86,92,94,95,99,100,101,103:114,116:119];
% shotlist  = 1e4 + 300 + [103:114,116:119]; 


% Acquiring Ne and Te data
Stem = '\MPEX::TOP.';
Branch = 'MACHOPS1:';
RootAddress = [Stem,Branch];
DataAddress{1} = [RootAddress,'LP_V_RAMP']; % V
DataAddress{2} = [RootAddress,'TARGET_LP']; % I

Config.tStart = 4.16; % [s]
Config.tEnd = 4.32;
Config.V_Att = 2;  % Output voltage of DLP box (Voltage) = V_Att*Digitized data 
Config.I_Att = 4;  % Output voltage of DLP box (Current) = I_att*Digitized data
Config.V_cal = [12.05,0.205];   % Voltage output of DLP = V_cal(1)*Output voltage of DLP box + V_cal(2)
Config.I_cal = [-142.5, 0.015]; % Current output of DLP = (I_cal(2) + Output voltage of DLP box)/Ical(1)
Config.FilterDataInput = 1; % Filter input data with savitsky Golay filter order 3 frame size 7
Config.TimeMode = 2; % Effective time of sweep, (1) start of ramp or (2) mean time of ramp
Config.AMU = 2; % Ion mass in AMU
Config.L_tip = 1.2/1000;
Config.D_tip = 0.254/1000; % [m]
Config.FitFunction = 2; 
Config.AreaType = 1; % Cylindrical + cap
[Ni,Te,time,Ifit,Ip,Vp,tm,Vsweep,Isweep] = DLP_fit_V2(Config,shotlist,DataAddress);

%%
C = {'k','r','bl','g','m','k:','r:','bl:','g:','m:',...
    'k','r','bl','g','m','k:','r:','bl:','g:','m:','k',...
    'r','bl','g','m','k:','r:','bl:','g:','m:','k','r',...
    'bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:'};

t = t_zero(shotlist);
tGather = 4.28;

figure; 
subplot(2,2,1); hold on
for s = 1:length(shotlist)
    h(s) = plot(time{s},Ni{s},C{s},'lineWidth',1);
    L{s} = [num2str(shotlist(s)),' ,t=',num2str(t{s}(10:14))];
    % Gather data into table
    nt = find(time{s}>=tGather); nt = nt(1:3);
    T(s,2) = mean(Ni{s}(nt)); 
    T(s,3) = mean(Te{s}(nt));
    T(s,1) = shotlist(s); 
end
legend(h,L,'location','NorthWest')
title('$ n_e $ $ [m^{-3}] $','interpreter','Latex','FontSize',13,'Rotation',0)
%set(gca,'PlotBoxAspectRatio',[1 1 1])
ylim([0,6e19])
xlim([4.15,4.32])

subplot(2,2,2); hold on
for s = 1:length(shotlist)
    plot(time{s},Te{s},C{s},'lineWidth',1)
end
title('$ T_e $ $ [eV] $','interpreter','Latex','FontSize',13,'Rotation',0)
ylim([0,10])
xlim([4.15,4.32])

subplot(2,2,3); hold on
for s = 1:length(shotlist)
    plot(time{s},e_c.*Ni{s}.*Te{s},C{s},'lineWidth',1)
end
ylim([0,20])
xlim([4.15,4.32])
title('$ P_e $ $ [Pa] $','interpreter','Latex','FontSize',13,'Rotation',0)

set(findobj('-Property','YTick'),'box','on')
set(gcf,'color','w')


T


if 1
    figure; hold on
    plot(tm{s},Vp{s})
    plot(tm{s},Ip{s}*1000)
    ylim([-100,100])
    xlim([4.15,4.32])
end

if 0
    
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