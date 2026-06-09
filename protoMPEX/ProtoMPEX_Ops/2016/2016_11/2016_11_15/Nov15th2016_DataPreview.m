close all
clear all

shot = 33.4; 

switch shot
    case 3; shotlist = 11400 + [3]; AddressType='s';CalType='niso';DLPType='9'; % Trip
    case 4; shotlist = 11400 + [4]; AddressType='s';CalType='niso';DLPType='9'; % Trip
    case 5; shotlist = 11400 + [5]; AddressType='s';CalType='niso';DLPType='9'; % Trip     
    case 6; shotlist = 11400 + [6]; AddressType='s';CalType='niso';DLPType='9'; % FP    
    case 19; shotlist = 11400 + [19]; AddressType='s';CalType='niso';DLPType='9'; % FP    
    case 21; shotlist = 11400 + [21]; AddressType='s';CalType='niso';DLPType='9'; % FP    
    case 26; shotlist = 11400 + [26]; AddressType='s';CalType='niso';DLPType='9'; % FP       
    case 30.4; shotlist = 11400 + [30]; AddressType='n';CalType='iso';DLPType='4'; % FP       
    case 30.9; shotlist = 11400 + [30]; AddressType='s';CalType='niso';DLPType='9'; % FP       
    case 32.4; shotlist = 11400 + [32]; AddressType='n';CalType='iso';DLPType='4'; % FP       
    case 32.9; shotlist = 11400 + [32]; AddressType='s';CalType='niso';DLPType='9'; % FP 

    case 33.4; shotlist = 11400 + [33]; AddressType='n';CalType='iso';DLPType='4'; % FP       
    case 33.9; shotlist = 11400 + [33]; AddressType='s';CalType='niso';DLPType='9'; % FP  
    case 34.4; shotlist = 11400 + [34]; AddressType='n';CalType='iso';DLPType='4'; % FP       
    case 34.9; shotlist = 11400 + [34]; AddressType='s';CalType='niso';DLPType='9'; % FP     
    case 35.4; shotlist = 11400 + [35]; AddressType='n';CalType='iso';DLPType='4'; % FP       
    case 35.9; shotlist = 11400 + [35]; AddressType='s';CalType='niso';DLPType='9'; % FP 
    case 36.4; shotlist = 11400 + [36]; AddressType='n';CalType='iso';DLPType='4'; % FP       
    case 36.9; shotlist = 11400 + [36]; AddressType='s';CalType='niso';DLPType='9'; % FP  
    case 38.4; shotlist = 11400 + [38]; AddressType='n';CalType='iso';DLPType='4'; % FP       
    case 38.9; shotlist = 11400 + [38]; AddressType='s';CalType='niso';DLPType='9'; % FP 
    case 39.4; shotlist = 11400 + [39]; AddressType='n';CalType='iso';DLPType='4'; % FP       
    case 39.9; shotlist = 11400 + [39]; AddressType='s';CalType='niso';DLPType='9'; % FP  
    case 41.4; shotlist = 11400 + [41]; AddressType='n';CalType='iso';DLPType='4'; % FP       
    case 41.9; shotlist = 11400 + [41]; AddressType='s';CalType='niso';DLPType='9'; % FP 
    case 42.4; shotlist = 11400 + [42]; AddressType='n';CalType='iso';DLPType='4'; % FP       
    case 42.9; shotlist = 11400 + [42]; AddressType='s';CalType='niso';DLPType='9'; % FP 
    case 45.4; shotlist = 11400 + [45]; AddressType='n';CalType='iso';DLPType='4'; % FP       
    case 45.9; shotlist = 11400 + [45]; AddressType='s';CalType='niso';DLPType='9'; % FP  
        
        
    case 47.4; shotlist = 11400 + [47]; AddressType='n';CalType='iso';DLPType='4'; % FP       
    case 47.9; shotlist = 11400 + [47]; AddressType='s';CalType='niso';DLPType='9'; % FP 
    case 48.4; shotlist = 11400 + [48]; AddressType='n';CalType='iso';DLPType='4'; % FP       
    case 48.9; shotlist = 11400 + [48]; AddressType='s';CalType='niso';DLPType='9'; % FP 
    case 49.4; shotlist = 11400 + [49]; AddressType='n';CalType='iso';DLPType='4'; % FP       
    case 49.9; shotlist = 11400 + [49]; AddressType='s';CalType='niso';DLPType='9'; % FP 
    case 50.4; shotlist = 11400 + [50]; AddressType='n';CalType='iso';DLPType='4'; % FP       
    case 50.9; shotlist = 11400 + [50]; AddressType='s';CalType='niso';DLPType='9'; % FP   
    case 51.4; shotlist = 11400 + [51]; AddressType='n';CalType='iso';DLPType='4'; % FP       
    case 51.9; shotlist = 11400 + [51]; AddressType='s';CalType='niso';DLPType='9'; % FP  
    case 53.4; shotlist = 11400 + [53]; AddressType='n';CalType='iso';DLPType='4'; % FP       
    case 53.9; shotlist = 11400 + [53]; AddressType='s';CalType='niso';DLPType='9'; % FP       
    case 54.4; shotlist = 11400 + [54]; AddressType='n';CalType='iso';DLPType='4'; % FP       
    case 54.9; shotlist = 11400 + [54]; AddressType='s';CalType='niso';DLPType='9'; % FP  
    case 55.4; shotlist = 11400 + [55]; AddressType='n';CalType='iso';DLPType='4'; % FP       
    case 55.9; shotlist = 11400 + [55]; AddressType='s';CalType='niso';DLPType='9'; % FP    
    case 57.4; shotlist = 11400 + [63]; AddressType='n';CalType='iso';DLPType='4'; % FP       
    case 57.9; shotlist = 11400 + [63]; AddressType='s';CalType='niso';DLPType='9'; % FP            
end

Config.tStart = 4.20; % [s]
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
Config.V_cal = [(0.46e-3)^-1,0]; % Voltage output of DLP = V_cal(1)*Output voltage of DLP box + V_cal(2) %  2.1739e+03
Config.I_cal = [-1,0]; % Current output of DLP = (I_cal(2) + Output voltage of DLP box)/Ical(1)
    case 'niso'
Config.V_cal = [12.05,0.205];   % Voltage output of DLP = V_cal(1)*Output voltage of DLP box + V_cal(2)
Config.I_cal = [-142.5, 0.015]; % Current output of DLP = (I_cal(2) + Output voltage of DLP box)/Ical(1)
end

switch DLPType 
    case '4'
DLP = 4.5;
Config.L_tip = 1.2/1000;
Config.D_tip = 0.5/1000; % [m]
    case '9'
DLP = 9.5;
Config.L_tip = 1.2/1000;
Config.D_tip = 0.254/1000; % [m]
end

Config.FilterDataInput = 1; % Filter input data with savitsky Golay filter order 3 frame size 7
Config.SGF = 7;
Config.TimeMode = 2; % Effective time of sweep, (1) start of ramp or (2) mean time of ramp
Config.AMU = 2; % Ion mass in AMU
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