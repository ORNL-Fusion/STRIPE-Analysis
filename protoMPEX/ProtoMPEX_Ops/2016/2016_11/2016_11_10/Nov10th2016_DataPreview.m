close all
clear all

shot = 67; 

switch shot
    case 65; shotlist = 11300 + [65]; AddressType='s';CalType='niso';DLPType='9'; % FP
    case 66; shotlist = 11300 + [66]; AddressType='n';CalType='iso' ;DLPType='9'; % T
    case 67; shotlist = 11300 + [67]; AddressType='n';CalType='iso' ;DLPType='9'; % FP
    case 68; shotlist = 11300 + [68]; AddressType='n';CalType='iso' ;DLPType='9'; % FP
    case 70; shotlist = 11300 + [70]; AddressType='n';CalType='iso' ;DLPType='9'; % FP
    case 71; shotlist = 11300 + [71]; AddressType='n';CalType='iso' ;DLPType='4'; % FP low Ne, noisy
    case 72; shotlist = 11300 + [72]; AddressType='n';CalType='iso' ;DLPType='4'; % FP low Ne, noisy
    case 73; shotlist = 11300 + [73]; AddressType='n';CalType='iso' ;DLPType='4'; % FP low Ne, noisy        
    case 74; shotlist = 11300 + [74]; AddressType='n';CalType='iso' ;DLPType='4'; % FP 4e19 Ne, noisy        
    case 75; shotlist = 11300 + [75]; AddressType='n';CalType='iso' ;DLPType='4'; % FP 4e19 Ne, noisy        
    case 76; shotlist = 11300 + [76]; AddressType='n';CalType='iso' ;DLPType='4'; % FP low Ne, noisy, changed match     
    case 77; shotlist = 11300 + [77]; AddressType='n';CalType='iso' ;DLPType='4'; % FP 4e19 Ne, noisy        
    case 78; shotlist = 11300 + [78]; AddressType='n';CalType='iso' ;DLPType='4'; % FP 4e19 Ne, noisy        
    case 80; shotlist = 11300 + [80]; AddressType='n';CalType='iso' ;DLPType='4'; % FP low Ne, noisy        
    case 81; shotlist = 11300 + [81]; AddressType='n';CalType='iso' ;DLPType='4'; % FP low Ne, noisy        
    case 82; shotlist = 11300 + [82]; AddressType='n';CalType='iso' ;DLPType='4'; % T        
    case 84; shotlist = 11300 + [84]; AddressType='n';CalType='iso' ;DLPType='4'; % FP low Ne, noisy        
    case 85; shotlist = 11300 + [85]; AddressType='n';CalType='iso' ;DLPType='4'; % NO sweep, I offset       
    case 87; shotlist = 11300 + [87]; AddressType='n';CalType='iso' ;DLPType='4'; % FP 4e19 Ne, less noisy       
    case 88; shotlist = 11300 + [88]; AddressType='n';CalType='iso' ;DLPType='4'; % FP low Ne, noisy        
    case 89; shotlist = 11300 + [89]; AddressType='n';CalType='iso' ;DLPType='4'; % FP low Ne, noisy, t<=4.3     
    case 90; shotlist = 11300 + [90]; AddressType='n';CalType='iso' ;DLPType='4'; % no sweep      
    case 91; shotlist = 11300 + [91]; AddressType='n';CalType='iso' ;DLPType='4'; % no sweep, 4.22<=t<=4.32
    case 92; shotlist = 11300 + [92]; AddressType='n';CalType='iso' ;DLPType='4'; % no sweep       
    case 93; shotlist = 11300 + [93]; AddressType='n';CalType='iso' ;DLPType='4'; % no sweep       
    case 94; shotlist = 11300 + [94]; AddressType='n';CalType='iso' ;DLPType='4'; % FP low Ne, noisy        
    case 95; shotlist = 11300 + [95]; AddressType='n';CalType='iso' ;DLPType='4'; % no sweep, swap cables     
    case 98; shotlist = 11300 + [98]; AddressType='n';CalType='iso' ;DLPType='4'; % FP, 4.19<=t<=4.32   
    case 98.5; shotlist = 11300 + [98]; AddressType='s';CalType='niso' ;DLPType='9'; % FP, 4.18<=t<=4.32   
    case 99; shotlist = 11300 + [99]; AddressType='n';CalType='iso' ;DLPType='4'; % FP,  
    case 99.5; shotlist = 11300 + [99]; AddressType='s';CalType='niso' ;DLPType='9'; % FP, data not good 
    case 100; shotlist = 11300 + [100]; AddressType='n';CalType='iso' ;DLPType='4'; % FP,              
    case 100.5; shotlist = 11300 + [100]; AddressType='s';CalType='niso' ;DLPType='9'; % FP,              
    case 101; shotlist = 11300 + [101]; AddressType='n';CalType='iso' ;DLPType='4'; % FP,              
    case 101.5; shotlist = 11300 + [101]; AddressType='s';CalType='niso' ;DLPType='9'; % FP,       
end

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

switch CalType 
    case 'iso'
Config.V_cal = [(0.46e-3)^-1,0]; % Voltage output of DLP = V_cal(1)*Output voltage of DLP box + V_cal(2)
Config.I_cal = [-1,0]; % Current output of DLP = (I_cal(2) + Output voltage of DLP box)/Ical(1)
    case 'niso'
Config.V_cal = [12.05,0.205];   % Voltage output of DLP = V_cal(1)*Output voltage of DLP box + V_cal(2)
Config.I_cal = [-142.5, 0.015]; % Current output of DLP = (I_cal(2) + Output voltage of DLP box)/Ical(1)
end

switch DLPType 
    case '4'
DLP = 4.5;
    case '9'
DLP = 9.5;
end

Config.FilterDataInput = 1; % Filter input data with savitsky Golay filter order 3 frame size 7
Config.SGF = 7;
Config.TimeMode = 2; % Effective time of sweep, (1) start of ramp or (2) mean time of ramp
Config.AMU = 2; % Ion mass in AMU
Config.L_tip = 1.2/1000;
Config.D_tip = 0.254/1000; % [m]
Config.FitFunction = 2; 
Config.AreaType = 1; % Cylindrical + cap

[ni,Te,time,Ifit,Ip,Vp,tm,Vsweep,Isweep] = DLP_fit_V5(Config,shotlist,DataAddress);

for s = 1:length(shotlist)
    Ni{s} = 0.5*(ni{s}{1} + ni{s}{2});
end



%%
close all
C = {'k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:'};
t = t_zero(shotlist);

figure; 
subplot(2,2,1); hold on
for s = 1:length(shotlist)
    plot(time{s},ni{s}{1},C{s},'lineWidth',1);
    plot(time{s},ni{s}{2},C{s},'lineWidth',1);
    h(s) = plot(time{s},Ni{s},C{s},'lineWidth',2);
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