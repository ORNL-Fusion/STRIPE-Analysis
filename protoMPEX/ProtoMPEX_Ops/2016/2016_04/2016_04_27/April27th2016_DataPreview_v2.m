close all
clear all

Extract = 0;
% NO PG4 data is available for this XP

if Extract
shotlist   = 8400 + [13 ,14 ,15 ,16 ,17 ,18,22,23,24 ,25 , 26, 29, 30, 32, 34, 35, 36];
TR2        =        [250,210,180,150,120,90,60,30,280,280,310,350,400,500,600,330,310];

AttType='a'; 
Config.L_tip = 1.2/1000;

Config.tStart = 4.15; % [s]
Config.tEnd = 4.32;

% Acquiring Ne and Te data
Stem = '\MPEX::TOP.';
Branch = 'MACHOPS1:';
RootAddress = [Stem,Branch];
DataAddress{1} = [RootAddress,'LP_V_RAMP']; % V
DataAddress{2} = [RootAddress,'TARGET_LP']; % I
Config.V_cal = [12.05,0.205];   % Voltage output of DLP = V_cal(1)*Output voltage of DLP box + V_cal(2)
Config.I_cal = [-142.5, 0.015]; % Current output of DLP = (I_cal(2) + Output voltage of DLP box)/Ical(1)

switch AttType
    case 'a'
Config.V_Att = 1;  % Output voltage of DLP box (Voltage) = V_Att*Digitized data 
Config.I_Att = 1;  % Output voltage of DLP box (Current) = I_att*Digitized data
    case 'b'
Config.V_Att = 1;  % Output voltage of DLP box (Voltage) = V_Att*Digitized data 
Config.I_Att = 1;  % Output voltage of DLP box (Current) = I_att*Digitized data
    end

Config.FilterDataInput = 1; % Filter input data with savitsky Golay filter order 3 frame size 7
Config.SGF = 5;
Config.TimeMode = 2; % Effective time of sweep, (1) start of ramp or (2) mean time of ramp
Config.AMU = 2; % Ion mass in AMU
Config.D_tip = 0.254/1000; % [m]
Config.FitFunction = 2; 
Config.AreaType = 1; % Cylindrical + cap

[ni,Te,time,Ifit,Ip,Vp,tm,Vsweep,Isweep] = DLP_fit_V5(Config,shotlist,DataAddress);

for s = 1:length(shotlist)
    Ni{s} = 0.5*(ni{s}{1} + ni{s}{2});
end

% Extract Gas data for specfic shots
Address = {'\MPEX::TOP.MACHOPS1:PG1'};
[PG1,t_pg1] = my_mdsvalue_v2(shotlist,Address);
Address = {'\MPEX::TOP.MACHOPS1:PG2'};
[PG2,t_pg2] = my_mdsvalue_v2(shotlist,Address);
Address = {'\MPEX::TOP.MACHOPS1:PG3'};
[PG3,t_pg3] = my_mdsvalue_v2(shotlist,Address);
Address = {'\MPEX::TOP.MACHOPS1:PG4'};
[PG4,t_pg4] = my_mdsvalue_v2(shotlist,Address);

    save('B_Helicon_Scan')
else
    load('B_Helicon_Scan')
end

%%
close all
C = {'k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:'};
t = t_zero(shotlist);

figure; hold on
for s = 1:length(shotlist)
    if 1
    plot3(time{s},TR2(s)*ones(size(ni{s}{1})),ni{s}{1},C{s},'lineWidth',1);
    plot3(time{s},TR2(s)*ones(size(ni{s}{2})),ni{s}{2},C{s},'lineWidth',1);
    end
    h(s) = plot3(time{s},TR2(s)*ones(size(ni{s}{2})),Ni{s},C{s},'lineWidth',2);
end
title('$ n_e $ $ [m^{-3}] $','interpreter','Latex','FontSize',13,'Rotation',0)
view([105,28])
zlim([0,8e19])

% GAS PRESSURE PG4
figure; hold on
for s = 1:length(shotlist)
    h(s) = plot3(t_pg3{s},TR2(s)*ones(size(PG3{s})),PG3{s}*2,C{s},'lineWidth',2);
end
title('$ P_n $ $ [mTorr] $','interpreter','Latex','FontSize',13,'Rotation',0)
view([105,28])

figure; hold on
for s = 1:length(shotlist)
    h(s) = plot3(time{s},TR2(s)*ones(size(Te{s})),Te{s},C{s},'lineWidth',2)
    L{s} = [num2str(shotlist(s)),' ,t=',num2str(t{s}(10:14))];
end
legend(h,L,'location','NorthEast')
title('$ T_e $ $ [eV] $','interpreter','Latex','FontSize',13,'Rotation',0)
view([105,28])
zlim([0,10])

figure; hold on
for s = 1:length(shotlist)
    plot3(time{s},TR2(s)*ones(size(Te{s})),e_c.*Ni{s}.*Te{s},C{s},'lineWidth',2)
end
title('$ P_e $ $ [Pa] $','interpreter','Latex','FontSize',13,'Rotation',0)
view([105,28])
zlim([0,10])

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