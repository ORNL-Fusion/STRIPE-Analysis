close all
clear all

shotlist = 13700 + [81,98,99,100,101];
shotlist = 13700 + [81,103,104,105];

shotlist = 13700 + [81,105,106,107,108]; % 81 reference shot from April 7th, 105 Today's shot

% shotlist = 13700 + [81,96,102];

shotlist = 13700 + [108,109]; % 112 kW ECH OFF and ON
shotlist = 13700 + [111,112]; % 100 kW ECH OFF and ON
DLPType='10MP';

shotlist = 13700 + [113,114,115,116,117]; % 100 and 110 kW ECH OFF, conditioning DLP6.5
DLPType='6';

shotlist = 13700 + [120,121,122,124]; % 110 kW ECH ON and OFF
DLPType='6';

shotlist = 13700 + [125,126]; % 110 kW ECH ON, DLP 6.5 and 10.5 are on-axis
DLPType='10MP';

% Begin radial scan
shotlist = 13800 + [27,28,29,30,31,32,33];
r        =         [0 ,1 ,1 ,2 ,2 ,2 ,2 ];
echPulse =         [1 ,0 ,1 ,0 ,1 ,1 ,1 ];

shotlist = 13800 + [32,33,36,37,38,40,44 ,45 ];
r        =         [2 ,2 ,3 ,3 ,4 ,4 ,3.5,3.5];
echPulse =         [0 ,1 ,1 ,0 ,0 ,1 ,1  ,0  ];

shotlist = 13800 + [46,47,48,49,50,51,52,53]
shotlist = 13800 + [52,53]
shotlist = 13800 + [32,33,54,55,56,57]
shotlist = 13800 + [56,57,58,59,60,61,62]

shotlist = 13800 + [63,64,65,66,67]
shotlist = 13800 + [66,67,68,69,70]
DLPType = '6';

shotlist = 13800 + [8,9]
shotlist = 13900 + [28,29]
DLPType = '10MP';

shotlist = 13850;

% Finding when the DLP 6.5 broke and got damaged
% shotlist = 13800 + [13:17,19]; % First try at DLP 6.5, on-axis
% shotlist = 13800 + [20]; % First try at DLP 6.5 + ECH, on-axis, doesnt work
shotlist = 13800 + [13:17,19,21:23]; % First try at DLP 6.5 + ECH, on-axis, doesnt work
% shotlist = 13800 + [22,24]; % DLP 6.5 + ECH, on-axis
% DLPType = '6'; Config.tStart = 4.13; Config.tEnd = 4.35;
% shotlist = 13800 + [25,26]; % leave DLP 6.5 on-axis, sweep DLP 10.5 on-axis
% DLPType = '10MP'; Config.tStart = 4.13; Config.tEnd = 4.35;
% shotlist = 13800 + [22,27,28,29]; 
% shotlist = 13800 + [22,28,30]; % DLP 6.5 y = 0,1,2, ECH OFF
% shotlist = 13800 + [30:33]; % DLP 6.5 y = 2,2,2
% shotlist = 13800 + [36,37,38]; % DLP 6.5 y = 3,3,4   
% shotlist = 13800 + [38,40,41]; % DLP 6.5 y = 4,4   
% shotlist = 13800 + [44,45]; % DLP 6.5 y = 3.5   
shotlist = 13800 + [46,47,48,49,50,51]; % DLP 6.5 y = 3.5. Probe shatter based on video   
% shotlist = 13800 + [61,62]; % DLP 6.5 y = 0 
DLPType = '6'; Config.tStart = 4.13; Config.tEnd = 4.35;


AddressType='s'; % s for standard
CalType='niso'; % niso for not isolated, is for isolated

% Acquiring Ne and Te data
Stem = '\MPEX::TOP.';
Branch = 'MACHOPS1:';
RootAddress = [Stem,Branch];

switch AddressType
    case 's'
DataAddress{1} = [RootAddress,'LP_V_RAMP']; % V
DataAddress{2} = [RootAddress,'TARGET_LP']; % I
Config.V_Att = 2;  % Output voltage of DLP box (Voltage) = V_Att*Digitized data 
Config.I_Att = 5;  % Output voltage of DLP box (Current) = I_att*Digitized data
    case 'n'
DataAddress{1} = [RootAddress,'INT_4MM_1']; % V
DataAddress{2} = [RootAddress,'INT_4MM_2']; % I
Config.V_Att = 1;  % Output voltage of DLP box (Voltage) = V_Att*Digitized data 
Config.I_Att = 1;  % Output voltage of DLP box (Current) = I_att*Digitized data
end

DA{1} = [RootAddress,'PWR_28GHZ']; % Isx
[ECH,t_ech]   = my_mdsvalue_v2(shotlist,DA(1))

DA{1} = [RootAddress,'RF_FWD_PWR'];
[RF,t_rf]   = my_mdsvalue_v2(shotlist,DA(1))

switch CalType 
    case 'iso'
Config.V_cal = [(0.46e-3)^-1,0]; % Voltage output of DLP = V_cal(1)*Output voltage of DLP box + V_cal(2) %  2.1739e+03
Config.I_cal = [-1,0]; % Current output of DLP = (I_cal(2) + Output voltage of DLP box)/Ical(1)
    case 'niso'
Config.V_cal = [12.05,0.205];   % Voltage output of DLP = V_cal(1)*Output voltage of DLP box + V_cal(2)
Config.I_cal = [-142.5, -0.6 + 0.015]; % Current output of DLP = (I_cal(2) + Output voltage of DLP box)/Ical(1)
Config.I_cal = [-142.5, 0.015]; % Current output of DLP = (I_cal(2) + Output voltage of DLP box)/Ical(1)
end

switch DLPType 
    case '4'
DLP = 4.5;
    case '9'
DLP = 9.5;
    case '6'
DLP = 6.5;
Config.L_tip = 1.0/1000; % 1.0 as of April 11th 2017
    case '10'
DLP = 10.5;
Config.L_tip = 1.2/1000;
    case '10MP'
DLP = 10.5;
Config.L_tip = 1.8/1000;
end

Config.FilterDataInput = 1; % Filter input data with savitsky Golay filter order 3 frame size 7
Config.SGF = 7;
Config.TimeMode = 2; % Effective time of sweep, (1) start of ramp or (2) mean time of ramp
Config.AMU = 2; % Ion mass in AMU
Config.D_tip = 0.254/1000; % [m]
Config.FitFunction = 2; 
Config.AreaType = 1; % 1: Cylindrical + cap

[ni,Te,time,Ifit,Ip,Vp,tm,Vsweep,Isweep] = DLP_fit_V5(Config,shotlist,DataAddress);

for s = 1:length(shotlist)
    Ni{s} = 0.5*(ni{s}{1} + ni{s}{2});
end

% note, 


%%
close all
C = {'k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:'};
t = t_zero(shotlist);

TimePlotEnd = 4.38;

figure; 
subplot(2,2,1); hold on
for s = 1:length(shotlist)
    plot(time{s},ni{s}{1},C{s},'lineWidth',1);
    plot(time{s},ni{s}{2},C{s},'lineWidth',1);
    h(s) = plot(time{s},Ni{s},C{s},'lineWidth',2);
    hold on
    plot(t_ech{s}(1:length(ECH{s})),ECH{s}*0.2e19,C{s},'lineWidth',0.5)
    plot(t_rf{s}(1:length(RF{s})),(RF{s}.^2)*6e19,C{s},'lineWidth',0.5)

end
title('$ n_e $ $ [m^{-3}] $','interpreter','Latex','FontSize',13,'Rotation',0)
legend(h,['DLP ',num2str(DLP)],'location','SouthEast')
%set(gca,'PlotBoxAspectRatio',[1 1 1])
ylim([0,6e19])
xlim([4.15,TimePlotEnd])

subplot(2,2,2); hold on
for s = 1:length(shotlist)
    h(s) = plot(time{s},Te{s},C{s},'lineWidth',2)
    L{s} = [num2str(shotlist(s)),' ,t=',num2str(t{s}(10:14))];
end
legend(h,L,'location','NorthEast')
title('$ T_e $ $ [eV] $','interpreter','Latex','FontSize',13,'Rotation',0)
ylim([0,6])
xlim([4.15,TimePlotEnd])

subplot(2,2,3); hold on
for s = 1:length(shotlist)
    plot(time{s},e_c.*Ni{s}.*Te{s},C{s},'lineWidth',2)
end
ylim([0,30])
xlim([4.15,TimePlotEnd])
title('$ P_e $ $ [Pa] $','interpreter','Latex','FontSize',13,'Rotation',0)

set(findobj('-Property','YTick'),'box','on')
set(gcf,'color','w')

GasCompare_2017_04_11

if 1
    figure; hold on
    for s = 1:length(shotlist)
    plot(tm{s},Vp{s})
    h(s) = plot(tm{s},Ip{s}*1000);
    ylim([-100,100])
    xlim([4.15,TimePlotEnd])
    grid on
    end
    legend(h,num2str(shotlist'))
    set(gcf,'color','w')
    box on
end

if 1
    
     figure;
    for c = 1:25;
        subplot(5,5,c); hold on
        plot(Vsweep{s}{c},Isweep{s}{c}*1e3,'k')
        plot(Vsweep{s}{c},Ifit{s}{c}*1e3,'r')
        ht = title(['T_e: ',num2str(Te{s}(c))]);
        set(ht,'FontSize',5); grid on
    end
    
    figure;
        for c = 26:50;
        subplot(5,5,c-25); hold on
        plot(Vsweep{s}{c},Isweep{s}{c}*1e3,'k')
        plot(Vsweep{s}{c},Ifit{s}{c}*1e3,'r')
        ht = title(['T_e: ',num2str(Te{s}(c))]);
        set(ht,'FontSize',5); grid on
        end
        
    figure;
        for c = 51:75;
        subplot(5,5,c-50); hold on
        plot(Vsweep{s}{c},Isweep{s}{c}*1e3,'k')
        plot(Vsweep{s}{c},Ifit{s}{c}*1e3,'r')
        ht = title(['T_e: ',num2str(Te{s}(c))]);
        set(ht,'FontSize',5); grid on
        end  
    
end

% % #########################################################################
% % Steady state values
% % #########################################################################
% 
% % Before ECH ==============================================================
% for s = 1:length(shotlist)
%     if shotlist(s) == 13611 || shotlist(s) == 13615
%         rng0 = find(time{s}>=4.25 & time{s}<=4.26);
%     else
%         rng0 = find(time{s}>=4.23 & time{s}<=4.24);
%     end
%         ne0(s)  = mean(Ni{s}(rng0));
%         dne0(s) = std(Ni{s}(rng0),1,2);
%         Te0(s)  = mean(Te{s}(rng0));
%         dTe0(s) = std(Te{s}(rng0),1,2);
% end
% 
% % During ECH: =============================================================
% for s = 1:length(shotlist)
% %     if shotlist(s) == 13611 || shotlist(s) == 13615
% %         rng0 = find(time{s}>=4.25 & time{s}<=4.26);
%     if shotlist(s) == 13599
%         rng1 = find(time{s}>=4.25 & time{s}<=4.29);
%     else
%         rng1 = find(time{s}>=4.27 & time{s}<=4.30);
%     end
%         ne1(s)  = mean(Ni{s}(rng1));
%         dne1(s) = std(Ni{s}(rng1),1,2);
%         Te1(s)  = mean(Te{s}(rng1));
%         dTe1(s) = std(Te{s}(rng1),1,2);
% end
% 
% figure; 
% subplot(1,2,1); hold on
% h = errorbar(PS1,ne1,dne1,'ro-','LineWidth',1)
% h = errorbar(PS1,ne0,dne0,'ko-','LineWidth',1); 
% ylim([0,3]*1e19)
% xlim([2000,5500]); 
% box on
% t = title('$n_e$ $[m^{-3}]$ '); t.Interpreter = 'latex'; t.FontSize = 13;
% 
% subplot(1,2,2); hold on
% h = errorbar(PS1,Te1,dTe1,'ro-','LineWidth',1)
% h = errorbar(PS1,Te0,dTe0,'ko-','LineWidth',1)
% ylim([0,3])
% xlim([2000,5500]); 
% box on
% t = title('$T_e$ $[eV]$ '); t.Interpreter = 'latex'; t.FontSize = 13;
% 
% set(gcf,'color','w')