% Sept 16th Bdot probe XP
close all
clear all

shotlist = 1e4 + 300 + [85];
shotlist = 1e4 + 800 + [64:66]
shotlist = 1e4 + 800 + [66,70,71,72,73];
 shotlist = 1e4 + 800 + [71:75,78];
 shotlist = 1e4 + 800 + [64,78,79];

 shotlist = 1e4 + 800 + [79,83,84,85];
 shotlist = 1e4 + 800 + [85,86,88,92,93,94,95];
 shotlist = 1e4 + 800 + [93,94,95,96,99,100,102];

 shotlist = 1e4 + 800 + [99,100,102,104,106];

  shotlist = 1e4 + 800 + [18,64,79,104,106];
  
% October 13th
shotlist = 1e4 + 900 + [22,23,27,28,29,30];
shotlist = 1e4 + 900 + [33,34,35,36,37,38,39,40,41];
shotlist = 1e4 + 900 + [36:39,40,41,42,43,44,45,46,47,49,50];

shotlist = 1e4 + 900 + [52,53,54,55,56];
shotlist = 1e4 + 900 + [9,14,18];
shotlist = 8600 + [5:10]

shotlist = 1e4 + 300 + [72,78]
shotlist = 8883

% shotlist = 9983;
shotlist = 1e4 + [904,1072,1073,1074,1075];

shotlist = 1e4 + [1076,1077,1078,1082,1091,1093,1094];
shotlist = 1e4 + 1e3 + [94,95,97,98,99]
shotlist = 1e4 + 1e3 + [100,101]
shotlist = 1e4 + 1e3 + [97,123,125];

% 88, LP_V_RAMP bad

% shotlist = 1e4 + 300 + [];

% Acquiring Ne and Te data
Stem = '\MPEX::TOP.';
Branch = 'MACHOPS1:';
RootAddress = [Stem,Branch];
DataAddress{1} = [RootAddress,'LP_V_RAMP']; % V
DataAddress{2} = [RootAddress,'TARGET_LP']; % I
Data = [RootAddress,'PWR_28GHz'];
% Data = [RootAddress,'PG2'];

[G,t_G] = my_mdsvalue_v3(shotlist,Data);

Config.tStart = 4.15; % [s]
Config.tEnd = 4.32;
Config.V_Att = 2;  % Output voltage of DLP box (Voltage) = V_Att*Digitized data 
Config.I_Att = 2;  % Output voltage of DLP box (Current) = I_att*Digitized data
Config.V_cal = [12.05,0.205];   % Voltage output of DLP = V_cal(1)*Output voltage of DLP box + V_cal(2)
Config.I_cal = [-142.5, 0.015]; % Current output of DLP = (I_cal(2) + Output voltage of DLP box)/Ical(1)
Config.FilterDataInput = 1; % Filter input data with savitsky Golay filter order 3 frame size 7
Config.TimeMode = 2; % Effective time of sweep, (1) start of ramp or (2) mean time of ramp
Config.AMU = 2; % Ion mass in AMU
Config.L_tip = 1.2/1000;
Config.D_tip = 0.254/1000; % [m]
Config.FitFunction = 2; 
Config.AreaType = 1; % Cylindrical + cap
[Ni,Te,time,Ifit,Ip,Vp,tm,Vsweep,Isweep] = DLP_fit_V4(Config,shotlist,DataAddress);

%%
close all
C = {'k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:'};
t = t_zero(shotlist);

figure; 
subplot(2,2,1); hold on
for s = 1:length(shotlist)
    h(s) = plot(time{s},Ni{s},C{s},'lineWidth',2);
    L{s} = [num2str(shotlist(s)),' ,t=',num2str(t{s}(10:14))];
end
legend(h,L,'location','NorthWest')
title('$ n_e $ $ [m^{-3}] $','interpreter','Latex','FontSize',13,'Rotation',0)
%set(gca,'PlotBoxAspectRatio',[1 1 1])
ylim([0,8e19])
xlim([4.15,4.32])

subplot(2,2,2); hold on
for s = 1:length(shotlist)
    plot(time{s},Te{s},C{s},'lineWidth',2)
    plot(t_G{s},0.5*G{s},C{s},'lineWidth',1)
end
title('$ T_e $ $ [eV] $','interpreter','Latex','FontSize',13,'Rotation',0)
ylim([0,10])
xlim([4.15,4.32])

subplot(2,2,3); hold on
for s = 1:length(shotlist)
    plot(time{s},e_c.*Ni{s}.*Te{s},C{s},'lineWidth',2)
   plot(t_G{s},G{s}*3,C{s},'lineWidth',1)
end
ylim([0,10])
xlim([4.15,4.32])
title('$ P_e $ $ [Pa] $','interpreter','Latex','FontSize',13,'Rotation',0)

set(findobj('-Property','YTick'),'box','on')
set(gcf,'color','w')


if 1
    figure; hold on
    for s = 1:length(shotlist)
    plot(tm{s},Vp{s})
    h(s) = plot(tm{s},Ip{s}*1000);
    ylim([-100,100])
    xlim([4.15,4.32])
    grid on
    end
    legend(h,num2str(shotlist'))
end

if 1
    
     figure;
    for c = 1:25;
        subplot(5,5,c); hold on
        plot(Vsweep{s}{c},Isweep{s}{c}*1e3,'k')
        plot(Vsweep{s}{c},Ifit{s}{c}*1e3,'r')
        ht = title(['T_e: ',num2str(Te{s}(c))]);
        set(ht,'FontSize',5)
    end
    
    figure;
        for c = 26:50;
        subplot(5,5,c-25); hold on
        plot(Vsweep{s}{c},Isweep{s}{c}*1e3,'k')
        plot(Vsweep{s}{c},Ifit{s}{c}*1e3,'r')
        ht = title(['T_e: ',num2str(Te{s}(c))]);
        set(ht,'FontSize',5)
        end
        
    figure;
        for c = 51:75;
        subplot(5,5,c-50); hold on
        plot(Vsweep{s}{c},Isweep{s}{c}*1e3,'k')
        plot(Vsweep{s}{c},Ifit{s}{c}*1e3,'r')
        ht = title(['T_e: ',num2str(Te{s}(c))]);
        set(ht,'FontSize',5)
        end  
    
end