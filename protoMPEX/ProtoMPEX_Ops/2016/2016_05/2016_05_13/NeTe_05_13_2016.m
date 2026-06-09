close all
clear all

CMPT = 0;

if CMPT
shotlist = 8800 + [10 , 11, 12, 13, 15, 17, 18, 19, 20, 21, 23, 25, 27, 29, 31, 32, 33, 34, 83, 84, 85, 87,88];
TR2      =        [210,270,290,310,330,350,370,390,410,430,450,470,500,600,230,200,170,140,170,140,110, 60, 0];

AddressType='s'; % s for standard
CalType='niso'; % niso for not isolated, iso for isolated
DLPType='9';

Config.tStart = 4.15; % [s]
Config.tEnd = 4.32;

% Acquiring Ne and Te data
Stem = '\MPEX::TOP.';
Branch = 'MACHOPS1:';
RootAddress = [Stem,Branch];

switch AddressType
    case 's'
DataAddress{1} = [RootAddress,'LP_V_RAMP']; % V
DataAddress{2} = [RootAddress,'TARGET_LP']; % I
Config.V_Att = 2;  % Output voltage of DLP box (Voltage) = V_Att*Digitized data 
Config.I_Att = 2;  % Output voltage of DLP box (Current) = I_att*Digitized data
    case 'n'
DataAddress{1} = [RootAddress,'INT_4MM_1']; % V
DataAddress{2} = [RootAddress,'INT_4MM_2']; % I
Config.V_Att = 1;  % Output voltage of DLP box (Voltage) = V_Att*Digitized data 
Config.I_Att = 1;  % Output voltage of DLP box (Current) = I_att*Digitized data
end

DA{1} = [RootAddress,'PWR_28GHZ']; % Isx
[ECH,t_ech]   = my_mdsvalue_v2(shotlist,DA(1));
DA{1} = [RootAddress,'RF_FWD_PWR'];
[RF,t_rf]   = my_mdsvalue_v2(shotlist,DA(1));

DA{1} = '\MPEX::TOP.MACHOPS1:PG2'; % 2.5
[PG2,t_pg2]   = my_mdsvalue_v2(shotlist,DA(1));
% PG4 was not installed for the May 13th 2016 experiment
% All other PGs are strongly affected by the magnetic field

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
end

Config.L_tip = 1.1/1000;

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
    save('NeTe_2016_05_13')
else
    load('NeTe_2016_05_13')    
end
% note, 


%%
close all
C = {'k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:'};
t = t_zero(shotlist);

TimePlotEnd = 4.33;

figure; 
subplot(2,2,1); hold on
for s = 1:length(shotlist)
%     plot(time{s},ni{s}{1},C{s},'lineWidth',1);
%     plot(time{s},ni{s}{2},C{s},'lineWidth',1);
    h(s) = plot(time{s},Ni{s},C{s},'lineWidth',2);
    hold on
%     plot(t_ech{s}(1:length(ECH{1})),ECH{s}*0.5e19,C{s},'lineWidth',0.5)
end
title('$ n_e $ $ [m^{-3}] $','interpreter','Latex','FontSize',13,'Rotation',0)
legend(h,['DLP ',num2str(DLP)],'location','SouthEast')
%set(gca,'PlotBoxAspectRatio',[1 1 1])
ylim([0,7e19])
xlim([4.15,TimePlotEnd])

subplot(2,2,2); hold on
for s = 1:length(shotlist)
    h(s) = plot(time{s},Te{s},C{s},'lineWidth',2)
    L{s} = [num2str(shotlist(s)),' ,t=',num2str(t{s}(10:14))];
end
legend(h,L,'location','NorthEast')
title('$ T_e $ $ [eV] $','interpreter','Latex','FontSize',13,'Rotation',0)
ylim([0,12])
xlim([4.15,TimePlotEnd])

subplot(2,2,3); hold on
for s = 1:length(shotlist)
    plot(time{s},e_c.*Ni{s}.*Te{s},C{s},'lineWidth',2)
end
ylim([0,40])
xlim([4.15,TimePlotEnd])
title('$ P_e $ $ [Pa] $','interpreter','Latex','FontSize',13,'Rotation',0)

set(findobj('-Property','YTick'),'box','on')
set(gcf,'color','w')


if 1
    figure; hold on
    [a,b] = sort(TR2,2,'ascend');
    seq = shotlist(b);
    for s = 1:length(b)
    subplot(5,5,s); hold on
    %plot(tm{s},Vp{s})
    h(b(s)) = plot(tm{b(s)}   ,abs(Ip{b(s)})*1000);
%               plot(tm{b(s)}   ,abs(Ifit{b(s)})*1000,'r')
              plot(t_rf{b(s)} ,150*RF{b(s)}.^2       ,'lineWidth',1,'color','r')
              plot(t_pg2{b(s)},15*PG2{b(s)}           ,'lineWidth',1,'color','g')
              ylim([-10,60])
    xlim([4.15,TimePlotEnd])
%     legend(h(b(s)),num2str(shotlist(b(s))'))
    legend(h(b(s)),num2str(TR2(b(s))'),'location','NorthWest')
    grid on

    end
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
        for c = 51:67;
        subplot(5,5,c-50); hold on
        plot(Vsweep{s}{c},Isweep{s}{c}*1e3,'k')
        plot(Vsweep{s}{c},Ifit{s}{c}*1e3,'r')
        ht = title(['T_e: ',num2str(Te{s}(c))]);
        set(ht,'FontSize',5); grid on
        end  
    
end

%% Ploting density as a function of time
close all
ta{1} = 4.25;
tb{1} = 4.27;

ta{2} = 4.20;
tb{2} = 4.21;

ta{3} = 4.18;
tb{3} = 4.19;


figure; hold on
for s = 1:length(shotlist)
    h(s) = plot3(time{s},TR2(s)*ones(size(Ni{s})),Ni{s});
    rng = find(time{s}>= ta{1} & time{s}<=tb{1});
    plot3(time{s}(rng),TR2(s)*ones(size(Ni{s}(rng))),Ni{s}(rng),'k','lineWidth',1);
    ne_m1(s) = mean(Ni{s}(rng));
    dne_m1(s) = std(Ni{s}(rng),1,2);
    
    rng = find(time{s}>= ta{2} & time{s}<=tb{2});
    plot3(time{s}(rng),TR2(s)*ones(size(Ni{s}(rng))),Ni{s}(rng),'r','lineWidth',1);
    ne_m2(s) = mean(Ni{s}(rng));
    dne_m2(s) = std(Ni{s}(rng),1,2);
    
    
    rng = find(time{s}>= ta{3} & time{s}<=tb{3});
    plot3(time{s}(rng),TR2(s)*ones(size(Ni{s}(rng))),Ni{s}(rng),'g','lineWidth',1);
    ne_m3(s) = mean(Ni{s}(rng));
    dne_m3(s) = std(Ni{s}(rng),1,2);
end
view([30,30])
zlim([0,1e20])

figure; hold on
[a,b] = sort(TR2);
errorbar(TR2(b),ne_m1(b),dne_m1(b),'ko-')
errorbar(TR2(b),ne_m2(b),dne_m2(b),'r.-')
errorbar(TR2(b),ne_m3(b),dne_m3(b),'g.-')

xlim([0,600])
ylabel('n_e [m^{-3}]')
xlabel('TR2 [A]')
set(gcf,'color','w')
grid on

for s = 1:length(shotlist)
    rng = find(t_rf{s}>4.29 & t_rf{s}<4.30);
    if RF{s}(rng)<0.2
        RF_T (:,s) =  RF{s-1};
    else
        RF_T (:,s) =  RF{s};
    end
end

RF_m = mean(RF_T,2)*26;
dRF_m = std(RF_T,1,2)*26;
RF_m_u = RF_m + 0.5*dRF_m;
RF_m_l = RF_m - 0.5*dRF_m;

figure; hold on 
plot(t_rf{1},RF_m.^2,'bl','LineWidth',2)
plot(t_rf{1},RF_m_u.^2,'bl')
plot(t_rf{1},RF_m_l.^2,'bl')

rng = find(t_rf{s}>= ta{1} & t_rf{s}<=tb{1});
plot(t_rf{1}(rng),RF_m(rng).^2,'k','LineWidth',5)

rng = find(t_rf{s}>= ta{2} & t_rf{s}<=tb{2});
plot(t_rf{1}(rng),RF_m(rng).^2,'r','LineWidth',5)

rng = find(t_rf{s}>= ta{3} & t_rf{s}<=tb{3});
plot(t_rf{1}(rng),RF_m(rng).^2,'g','LineWidth',5)
xlim([4.13,4.33])
xlabel('Time [s]')
ylabel('RF power [kW]')

set(gcf,'color','w')
grid on; box on

save('PawelData')

% figure; hold on
% for s = 1:length(shotlist)
%     h(s) = plot3(time{b(s)},TR2(b(s))*ones(size(Te{b(s)})),Te{b(s)});
% end
% view([30,30])

