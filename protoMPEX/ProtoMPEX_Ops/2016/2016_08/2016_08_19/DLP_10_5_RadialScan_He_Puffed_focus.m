close all
clear all

CMPT = 0;

if CMPT
    shotlist = 9990 + [5,7   ,9   ,10,11  ,12  ,13  ,14 ,15  ,16, 17 ,18  ,19 ,20   ,  21,   22,23,   24,25  ,18,19  ,2    ];
    % 8 was a trip
    x        =        [0,0.25,0.25,1 ,1.25,1.25,1.25,1.5,1.75,2 ,2.25,2.25,2.5,-0.25,-0.5,-0.75,-1,-1.25,-1.5,-2,-2.25,-2.5];
  
    Config.tStart = 4.15; % [s]
    Config.tEnd = 4.285;

    % Acquiring Ne and Te data
    Stem = '\MPEX::TOP.';
    Branch = 'MACHOPS1:';
    RootAddress = [Stem,Branch];
    DataAddress{1} = [RootAddress,'LP_V_RAMP']; % V
    DataAddress{2} = [RootAddress,'TARGET_LP']; % I
    Config.V_Att = 2;  % Output voltage of DLP box (Voltage) = V_Att*Digitized data 
    Config.I_Att = 2;  % Output voltage of DLP box (Current) = I_att*Digitized data

    Config.V_cal = [12.05,0.205];   % Voltage output of DLP = V_cal(1)*Output voltage of DLP box + V_cal(2)
    Config.I_cal = [-142.5, -0.6 + 0.015]; % Current output of DLP = (I_cal(2) + Output voltage of DLP box)/Ical(1)
  
    DLP = 10.5;
    Config.L_tip = 1.03/1000;

    Config.FilterDataInput = 1; % Filter input data with savitsky Golay filter order 3 frame size 7
    Config.SGF = 7;
    Config.TimeMode = 2; % Effective time of sweep, (1) start of ramp or (2) mean time of ramp
    Config.AMU = 2; % Ion mass in AMU
    Config.D_tip = 0.254/1000; % [m]
    Config.FitFunction = 2; 
    Config.AreaType = 1; % Cylindrical + cap

    [ni,Te,time,Ifit,Ip,Vp,tm,Vsweep,Isweep] = DLP_fit_V5(Config,shotlist,DataAddress);

    for s = 1:length(shotlist)
        Ni{s} = 0.5*(ni{s}{1} + ni{s}{2});
    end
    
    address{1} = '\MPEX::TOP.MACHOPS1:PG4';
    [PG4,t_pg4] = my_mdsvalue_v2(shotlist,address(1));
    
    address{1} = '\MPEX::TOP.MACHOPS1:PG1';
    [PG9,t_pg9] = my_mdsvalue_v2(shotlist,address(1));

    address{1} = '\MPEX::TOP.MACHOPS1:RF_FWD_PWR';
    [Prf,t_prf] = my_mdsvalue_v2(shotlist,address(1));
    
    DA{1} = [RootAddress,'PWR_28GHZ'];
    [ECH,t_ech]   = my_mdsvalue_v2(shotlist,DA(1));

    save('Aug19th_2016_DLP105_RadialScan')
else
    load('Aug19th_2016_DLP105_RadialScan')

end

%%
close all
C = {'k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:'};
t = t_zero(shotlist);

figure; 
% subplot(2,2,1); 
hold on
for s = 1:length(shotlist)
%     plot(time{s},ni{s}{1},C{s},'lineWidth',1);
%     plot(time{s},ni{s}{2},C{s},'lineWidth',1);
    rng = find(time{s}>=4.22 & time{s}<=4.27);
    h(s) = plot3(time{s}(rng),x(s)*ones(size(Ni{s}(rng))),Ni{s}(rng),C{s},'lineWidth',2);
    hold on
end
title('$ n_e $ $ [m^{-3}] $','interpreter','Latex','FontSize',13,'Rotation',0)
zlim([0,2e19])
view([90,0])

figure; 
hold on
for s = 1:length(shotlist)
    rng = find(time{s}>=4.20 & time{s}<=4.27);
    h(s) = plot3(time{s}(rng),x(s)*ones(size(Te{s}(rng))),Te{s}(rng),C{s},'lineWidth',2);
    hold on
end
title('$ T_e $ $ [eV] $','interpreter','Latex','FontSize',13,'Rotation',0)
zlim([0,10])
view([90,0])

figure; 
hold on
for s = 1:length(shotlist)
    rng = find(time{s}>=4.20 & time{s}<=4.27);
    h(s) = plot3(time{s}(rng),x(s)*ones(size(Te{s}(rng))),e_c*Te{s}(rng).*Ni{s}(rng),C{s},'lineWidth',2);
    hold on
end
title('$ T_e $ $ [eV] $','interpreter','Latex','FontSize',13,'Rotation',0)
zlim([0,15])
view([90,0])
% #########################################################################
%% Gas pressure
address{1} = '\MPEX::TOP.MACHOPS1:PG4';
[PG4_0,t_pg4_0] = my_mdsvalue_v2(9998,address(1));
address{1} = '\MPEX::TOP.MACHOPS1:PG1';
[PG9_0,t_pg9_0] = my_mdsvalue_v2(9998,address(1));

figure
hold on
plot(t_pg4_0{1},PG4_0{1}*10,'k','lineWidth',3);
% plot(t_pg9_0{1},PG9_0{1}*10,'k','lineWidth',3);

plot(t_prf{8},70*Prf{8}.^2,'g','lineWidth',3);

for s = 1:length(shotlist)
    h4(s) = plot(t_pg4{s},PG4{s}*10,C{s},'lineWidth',2);
    h9(s) = plot(t_pg9{s},PG9{s}*10,C{s},'lineWidth',1);

    
    hold on
end
ylim([0,20])
xlim([4.14,4.6])

%%
% ########################################################################
figure; 
subplot(2,2,1); hold on
for s = 1:length(shotlist)
%     plot(time{s},ni{s}{1},C{s},'lineWidth',1);
%     plot(time{s},ni{s}{2},C{s},'lineWidth',1);
    h(s) = plot(time{s},Ni{s},C{s},'lineWidth',2);
    hold on
end
title('$ n_e $ $ [m^{-3}] $','interpreter','Latex','FontSize',13,'Rotation',0)
legend(h,['DLP ',num2str(DLP)],'location','SouthEast')
%set(gca,'PlotBoxAspectRatio',[1 1 1])
ylim([0,8e19])
xlim([4.15,4.32])

subplot(2,2,2); hold on
for s = 1:length(shotlist)
    h(s) = plot(time{s},Te{s},C{s},'lineWidth',2)
    L{s} = [num2str(shotlist(s)),' ,t=',num2str(t{s}(10:14))];
end
legend(h,L,'location','NorthEast')
title('$ T_e $ $ [eV] $','interpreter','Latex','FontSize',13,'Rotation',0)
ylim([0,12])
xlim([4.15,4.32])

subplot(2,2,3); hold on
for s = 1:length(shotlist)
    plot(time{s},e_c.*Ni{s}.*Te{s},C{s},'lineWidth',2)
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
    hp(s) = plot(t_pg4{s},PG4{s}*40,'r');
    hrf(s) = plot(t_prf{s},Prf{s}*100,'g');
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