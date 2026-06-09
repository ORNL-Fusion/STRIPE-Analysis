% May 3rd 2017, 
% DLP 1.5 and DLP 10.5 are both on-axis in the plasma during ECH pulses, we
% have measured the traces simultaneously

close all
clear all

sweep = 'niso'; % niso: not isolated (DLP 10.5), % iso: isolated (DLP 1.5)

% DLP 1.5 and DLP 10.5 are used at the same time in the following shots:
shotlist = 14200 + [65,66,69,97];

switch sweep
    case 'niso'
DLPType = '10MP'; % DLP 10.5 vertical
ChannelType = '1' ; % "1" for TARGET_LP,LP_V_RAMP, "2" for other options
Config.tStart = 4.15; % [s]
Config.tEnd = 4.45;
SweepType = 'niso'; % "iso" for isolated sweep, "niso" for non-isolated sweep
AttType = 'Vx2,Ix5'; % Attenuation on the digitized signals
Config.Center_V = 1; % Remove offset on V: 1 (yes) 0(no)
Config.Center_I = 0; % Remove offset on I: 1 (yes) 0(no)

    case 'iso'
% Isolated sweep drive
DLPType = '1MP'; % DLP 1.5 horizontal
ChannelType = '4';%2 ; % "1" for TARGET_LP,LP_V_RAMP, "2" for other options
Config.tStart = 4.18; % [s]
Config.tEnd = 4.45;
SweepType = 'iso'; % "iso" for isolated sweep, "niso" for non-isolated sweep
AttType = 'Vx1,Ix1'; % Attenuation on the digitized signals
Config.Center_V = 1; % Remove offset on V: 1 (yes) 0(no)
Config.Center_I = 1; % Remove offset on I: 1 (yes) 0(no)
end

%%
% #########################################################################
% DLP FITTING ROUTINE
% #########################################################################
switch DLPType 
    case '4'
        DLP = 4.5;
        Config.L_tip = 1.2/1000; % 1.0 as of April 11th 2017
        Config.D_tip = 2*0.254/1000; % [m]
    case '9'
        DLP = 9.5;
    case '6'
        DLP = 6.5;
        Config.L_tip = 1.0/1000; % 1.0 as of April 11th 2017
        Config.D_tip = 0.254/1000; % [m]
    case '10'
        DLP = 10.5;
        Config.L_tip = 1.2/1000;
        Config.D_tip = 0.254/1000; % [m]
    case '10MP'
        DLP = 10.5;
        Config.L_tip = 1.8/1000;
        Config.D_tip = 0.254/1000; % [m]
    case '4MP'
        DLP = 4.5;
        Config.L_tip = 1.7/1000;
        Config.D_tip = 0.254/1000; % [m]
    case '1MP'
        DLP = 1.5;
        Config.L_tip = 1.7/1000;
        Config.D_tip = 0.254/1000; % [m]
end
% #########################################################################
switch AttType
    case 'Vx2,Ix5'
        Config.V_Att = 2;  % Output voltage of DLP box (Voltage) = V_Att*Digitized data 
        Config.I_Att = 5;  % Output voltage of DLP box (Current) = I_att*Digitized data
    case 'Vx1,Ix1'
        Config.V_Att = 1;  % Output voltage of DLP box (Voltage) = V_Att*Digitized data 
        Config.I_Att = 1;  % Output voltage of DLP box (Current) = I_att*Digitized data
end
% #########################################################################
Stem = '\MPEX::TOP.'; Branch = 'MACHOPS1:'; RootAddress = [Stem,Branch];
switch ChannelType
    case '1'
        DataAddress{1} = [RootAddress,'LP_V_RAMP']; % V
        DataAddress{2} = [RootAddress,'TARGET_LP']; % I
    case '2'
        DataAddress{1} = [RootAddress,'INT_4MM_1']; % V
        DataAddress{2} = [RootAddress,'INT_4MM_2']; % I
    case '3'
        DataAddress{1} = [RootAddress,'GEN_RF_PWR']; % V
        DataAddress{2} = [RootAddress,'INT_4MM_2'];  % I
    case '4'
        DataAddress{1} = [RootAddress,'ICH_LP']; % V
        DataAddress{2} = [RootAddress,'INT_4MM_2'];  % I
end
% #########################################################################
switch SweepType
    case 'iso'
Config.V_cal = [(0.46e-3)^-1,0]; % Voltage output of DLP = V_cal(1)*Output voltage of DLP box + V_cal(2) %  2.1739e+03
Config.I_cal = [-1,0]; % Current output of DLP = (I_cal(2) + Output voltage of DLP box)/Ical(1)
    case 'niso'
Config.V_cal = [12.05,0.205];   % Voltage output of DLP = V_cal(1)*Output voltage of DLP box + V_cal(2)
Config.I_cal = [-142.5, -0.6 + 0.015]; % Current output of DLP = (I_cal(2) + Output voltage of DLP box)/Ical(1)
Config.I_cal = [-142.5, 0.015]; % Current output of DLP = (I_cal(2) + Output voltage of DLP box)/Ical(1)
end
% #########################################################################
Config.FilterDataInput = 1; % Filter input data with savitsky Golay filter order 3 frame size 7
Config.SGF = 7;
Config.TimeMode = 2; % Effective time of sweep, (1) start of ramp or (2) mean time of ramp
Config.AMU = 2; % Ion mass in AMU
Config.FitFunction = 2; 
Config.AreaType = 1; % 1: Cylindrical + cap
% #########################################################################
[ni,Te,time,Ifit,Ip,Vp,tm,Vsweep,Isweep,GlitchFlag,SSQres] = DLP_fit_V5_2(Config,shotlist,DataAddress);
for s = 1:length(shotlist)
    Ni{s} = 0.5*(ni{s}{1} + ni{s}{2});
end
DA{1} = [RootAddress,'RF_FWD_PWR'];
[RF,t_rf]   = my_mdsvalue_v2(shotlist,DA(1))
DA{1} = [RootAddress,'PWR_28GHZ']; % Isx
[ECH,t_ech]   = my_mdsvalue_v2(shotlist,DA(1))
% #########################################################################
%%
% #########################################################################
% PLOT DATA:
% #########################################################################

close all
GasCompareALL_2017_04_21
C = {'k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:'};
t = t_zero(shotlist);

TimePlotEnd = 4.5;
ssqLevel = 0.001;

figure; 
subplot(2,2,1); hold on
for s = 1:length(shotlist)
    rng = find(SSQres{s}<ssqLevel);
    rng = find(GlitchFlag{s}==0);
%      plot(time{s},ni{s}{1},C{s},'lineWidth',1);
%      plot(time{s},ni{s}{2},C{s},'lineWidth',1);
    h(s) = plot(time{s}(rng),Ni{s}(rng),C{s},'lineWidth',2);

    hold on
    plot(t_rf{s}(1:length(RF{s})),(RF{s}.^2)*6e19,C{s},'lineWidth',0.5)
    plot(t_ech{s}(1:length(ECH{s})),ECH{s}*0.5e19,C{s},'lineWidth',0.5)

end
title('$ n_e $ $ [m^{-3}] $','interpreter','Latex','FontSize',13,'Rotation',0)
legend(h,['DLP ',num2str(DLP)],'location','SouthEast')
%set(gca,'PlotBoxAspectRatio',[1 1 1])
ylim([0,8e19])
xlim([4.15,TimePlotEnd])

subplot(2,2,2); hold on
for s = 1:length(shotlist)
    rng = find(SSQres{s}<ssqLevel);
    rng = find(GlitchFlag{s}==0);
    h(s) = plot(time{s}(rng),Te{s}(rng),C{s},'lineWidth',2)
    L{s} = [num2str(shotlist(s)),' ,t=',num2str(t{s}(10:14))];
end
legend(h,L,'location','NorthEast')
title('$ T_e $ $ [eV] $','interpreter','Latex','FontSize',13,'Rotation',0)
ylim([0,6])
xlim([4.15,TimePlotEnd])

subplot(2,2,3); hold on
for s = 1:length(shotlist)
    rng = find(SSQres{s}<ssqLevel);
    rng = find(GlitchFlag{s}==0);
    plot(time{s}(rng),e_c.*Ni{s}(rng).*Te{s}(rng),C{s},'lineWidth',2)
end
ylim([0,30])
xlim([4.15,TimePlotEnd])
title('$ P_e $ $ [Pa] $','interpreter','Latex','FontSize',13,'Rotation',0)

set(findobj('-Property','YTick'),'box','on')
set(gcf,'color','w')

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
        plot(Vsweep{s}{c},Isweep{s}{c}*1e3,'k');
        hvi(c) = plot(Vsweep{s}{c},Ifit{s}{c}*1e3,'r')
        ht = title(['T_e: ',num2str(Te{s}(c))]);
        set(ht,'FontSize',5); grid on
    end
    
    figure;
        for c = 26:50;
        subplot(5,5,c-25); hold on
        plot(Vsweep{s}{c},Isweep{s}{c}*1e3,'k')
        hvi(c) = plot(Vsweep{s}{c},Ifit{s}{c}*1e3,'r')
        ht = title(['T_e: ',num2str(Te{s}(c))]);
        set(ht,'FontSize',5); grid on
        end
        
    figure;
        for c = 51:75;
        subplot(5,5,c-50); hold on
        plot(Vsweep{s}{c},Isweep{s}{c}*1e3,'k')
        hvi(c) = plot(Vsweep{s}{c},Ifit{s}{c}*1e3,'r')
        ht = title(['T_e: ',num2str(Te{s}(c))]);
        set(ht,'FontSize',5); grid on
        end  
        
        figure;
        for c = 76:100;
        subplot(5,5,c-75); hold on
        plot(Vsweep{s}{c},Isweep{s}{c}*1e3,'k')
        hvi(c) = plot(Vsweep{s}{c},Ifit{s}{c}*1e3,'r')
        ht = title(['T_e: ',num2str(Te{s}(c))]);
        set(ht,'FontSize',5); grid on
        end      
        
        rng = find(SSQres{s}>ssqLevel);
        rng = find(GlitchFlag{s}==1);

        for s = 1:length(rng)
        set(hvi(rng(s)),'LineWidth',2)
        end
end