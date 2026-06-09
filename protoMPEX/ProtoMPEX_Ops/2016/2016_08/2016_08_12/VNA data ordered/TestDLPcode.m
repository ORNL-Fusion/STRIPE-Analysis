% Test
close all
clear all

% Acquiring Ne and Te data
Stem = '\MPEX::TOP.';
Branch = 'MACHOPS1:';
RootAddress = [Stem,Branch];
DataAddress{1} = [RootAddress,'LP_V_RAMP']; % V
DataAddress{2} = [RootAddress,'TARGET_LP']; % I

Config.tStart = 4.16; % [s]
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
shotlist = 9700 + [95,96];
[Ni,Te,time,Ifit,Ip,Vp,tm,Vsweep,Isweep] = DLP_fit_V1(Config,shotlist,DataAddress);

%%
C = {'k','r','bl','g','m','k:','r:','bl:','g:','m:'};
t = t_zero(shotlist);

figure; 
subplot(2,1,1); hold on
for s = 1:length(shotlist)
    h(s) = plot(time{s},Ni{s},C{s});
    L{s} = [num2str(shotlist(s)),' ,t=',num2str(t{s}(10:14))];
end
%L = num2str(shotlist');
%leg(1) = [L(1,:),' ,t = ',num2str(t{1}(10:14))];
% legend(h,[num2str(shotlist')],'location','NorthWest')
legend(h,L,'location','NorthWest')

ylim([0,6e19])

subplot(2,1,2); hold on
for s = 1:length(shotlist)
    plot(time{s},Te{s},C{s})
end
ylim([0,10])


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
        ht = title(['T_e: ',num2str(Te{s}(c)),'t: ',num2str(time_a{s}(c))]);
        set(ht,'FontSize',5)
        end
        
    figure;
        for c = 51:75;
        subplot(5,5,c-50); hold on
        plot(Vsweep{s}{c},Isweep{s}{c}*1e3,'k')
        plot(Vsweep{s}{c},Ifit{s}{c}*1e3,'r')
        ht = title(['T_e: ',num2str(Te_a{s}(c)),'t: ',num2str(time_a{s}(c))]);
        set(ht,'FontSize',5)
        end  
    
end