% October 6th Bdot probe XP
close all
clear all

shotlist = 1e4 + 700 + [35,38,39,41,42,43,44,45,46:48];

% shotlist = 1e4 + 700 + [47,48,49,50,51];
shotlist = 1e4 + 700 + [35,38,39,41,42];
shotlist = 1e4 + 700 + [54,55,56,57,59,60]

shotlist = 1e4 + 700 + [87,88,90,91];


% Acquiring Ne and Te data
Stem = '\MPEX::TOP.';
Branch = 'MACHOPS1:';
RootAddress = [Stem,Branch];
DataAddress{1} = [RootAddress,'LP_V_RAMP']; % V
DataAddress{2} = [RootAddress,'TARGET_LP']; % I
Data{5} = [RootAddress,'PWR_28GHz'];
[G,t_G] = my_mdsvalue_v3(shotlist,Data{5});

Data{1} = ['\MPEX::TOP.MACHOPS1:INT_4MM_1'];
Data{2} = ['\MPEX::TOP.MACHOPS1:INT_4MM_2'];
Data{3} = ['\MPEX::TOP.MACHOPS1:INT_2MM_1'];
Data{4} = ['\MPEX::TOP.MACHOPS1:INT_2MM_2'];

[VmagRatio,P,tBdot] = BdotPhaseDetect_v1(shotlist,Data);

Config.tStart = 4.10; % [s]
Config.tEnd = 4.38;
Config.V_Att = 2;  % Output voltage of DLP box (Voltage) = V_Att*Digitized data 
Config.I_Att = 2;  % Output voltage of DLP box (Current) = I_att*Digitized data
Config.V_cal = [12.05,0.205];   % Voltage output of DLP = V_cal(1)*Output voltage of DLP box + V_cal(2)
Config.I_cal = [-142.5, 1*0.015]; % Current output of DLP = (I_cal(2) + Output voltage of DLP box)/Ical(1)
Config.FilterDataInput = 1; % Filter input data with savitsky Golay filter order 3 frame size 7
Config.TimeMode = 2; % Effective time of sweep, (1) start of ramp or (2) mean time of ramp
Config.AMU = 2; % Ion mass in AMU
Config.L_tip = 1.2/1000;
Config.D_tip = 0.254/1000; % [m]
Config.FitFunction = 2; 
Config.AreaType = 1; % Cylindrical + cap
[Ni,Te,time,Ifit,Ip,Vp,tm,Vsweep,Isweep] = DLP_fit_V2(Config,shotlist,DataAddress);

%%
close all

figure; subplot(2,1,1); hold on
for s = 1:length(shotlist)
    hP(s) = plot(P{s});
end
title('Phase')
legend(hP,num2str(shotlist'))
ylim([-180,180])

subplot(2,1,2); hold on
for s = 1:length(shotlist)
    hM(s) = plot(VmagRatio{s});
end
title('Magnitude')
legend(hM,num2str(shotlist'))
%xlim([4.15,4.32])
ylim([0,10])

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
ylim([0,6e19])
xlim([4.15,4.32])

subplot(2,2,2); hold on
for s = 1:length(shotlist)
    plot(time{s},Te{s},C{s},'lineWidth',2)
    plot(t_G{s},G{s},C{s},'lineWidth',1)
end
title('$ T_e $ $ [eV] $','interpreter','Latex','FontSize',13,'Rotation',0)
ylim([0,10])
xlim([4.15,4.32])

subplot(2,2,3); hold on
for s = 1:length(shotlist)
    plot(time{s},e_c.*Ni{s}.*Te{s},C{s},'lineWidth',2)
    plot(t_G{s},G{s}*3,C{s},'lineWidth',1)
end
ylim([0,40])
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

if 0
    
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