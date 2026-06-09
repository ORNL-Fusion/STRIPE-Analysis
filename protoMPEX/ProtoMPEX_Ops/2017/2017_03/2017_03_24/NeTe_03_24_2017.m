close all
clear all

shotlist = 13000 + [494:495];

shotlist = 13000 + [495,499,500,501,502];
shotlist = 13000 + [509,510]; % ECH ON, x = 0.5
shotlist = 13000 + [509,511]; % ECH OFF
shotlist = 13000 + [511,512]; % ECH ON  x = 1.0
shotlist = 13000 + [511,514,515,516,517]; % ECH OFF

shotlist = 13000 + [517,526]; % ECH ON, x = 1.5 
shotlist = 13000 + [517,527]; % ECH OFF
shotlist = 13000 + [527,528]; % ECH ON, x = 2.0
shotlist = 13000 + [527,529]; % ECH OFF
shotlist = 13000 + [529,530]; % ECH ON, x = 2.5
shotlist = 13000 + [529,531]; % ECH OFF
shotlist = 13000 + [531,532]; % ECH ON, x = 3.0
shotlist = 13000 + [531,533]; % ECH OFF
shotlist = 13000 + [533,534]; % ECH ON, x = 3.5
shotlist = 13000 + [508,535]; % ECH OFF
shotlist = 13000 + [535,536,537]; % ECH ON, x = 0
shotlist = 13000 + [535,538]; % ECH OFF
shotlist = 13000 + [538,539,540]; % ECH ON, x = -0.5
shotlist = 13000 + [538,541]; % ECH OFF
shotlist = 13000 + [541,542]; % ECH ON, x = -1
shotlist = 13000 + [541,543]; % ECH OFF
shotlist = 13000 + [543,544]; % ECH ON, x = -1.5
shotlist = 13000 + [543,546]; % ECH OFF
shotlist = 13000 + [546,545]; % ECH ON, x = -2.0
shotlist = 13000 + [546,548]; % ECH OFF
shotlist = 13000 + [547,548]; % ECH ON, x = -2.5

% Begin Global field scan
shotlist = 13000 + [549 ,550 ,551 ,552 ,553 ,554 ,555 ,556 ,557 ,558 ,559 ,560 ];
GF       =       + [4000,3900,3800,3800,3700,3600,3500,3400,4100,4200,4300,4400];

shotlist = 13000 + [561:567];
% % shotlist = 13000 + [568,569];
% 
% shotlist = 13000 + [551]; % 3800 A
% shotlist = 13000 + [560]; % 4400 A
% shotlist = 13000 + [563]; % 4700 A
% shotlist = 13000 + [564]; % 4900 A
% shotlist = 13000 + [556]; % 3400 A
% shotlist = 13000 + [556,564];
% 
% % Radial scan
% shotlist = 13000 + [531,532]; % x = 3 cm
% shotlist = 13000 + [511,512]; % x = 1 cm
% 
% shotlist = 13000 + [565:567]; % x = 1 cm


AddressType='s'; % s for standard
CalType='niso'; % niso for not isolated, iso for isolated
DLPType='9';

Config.tStart = 4.18; % [s]
Config.tEnd = 4.4;

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

Config.L_tip = 1.2/1000;
Config.L_tip = 1.9/1000;
% Config.L_tip = 1.2/1000;


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

% note, 


%%
close all
C = {'k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:'};
t = t_zero(shotlist);

TimePlotEnd = 4.38;

figure; 
subplot(2,2,1); hold on
for s = 1:length(shotlist)
%     plot(time{s},ni{s}{1},C{s},'lineWidth',1);
%     plot(time{s},ni{s}{2},C{s},'lineWidth',1);
    h(s) = plot(time{s},Ni{s},C{s},'lineWidth',2);
    hold on
    plot(t_ech{s}(1:length(ECH{1})),ECH{s}*0.5e19,C{s},'lineWidth',0.5)
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
    for s = 1:length(shotlist)
    %plot(tm{s},Vp{s})
    h(s) = plot(tm{s},Ip{s}*1000);
    ylim([-100,100])
    xlim([4.15,TimePlotEnd])
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

%GasCompare_03_24_17