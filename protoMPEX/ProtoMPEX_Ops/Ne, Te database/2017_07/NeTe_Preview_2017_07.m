% This code was created on 2017_07_23 for the purpose of helping calculate
% the data needed for the Ne Te database to be created on July 2017

clear all
close all

ProbeLoc = 'A';
FitShow = 1; 
% Since we now have several probe drives, we can run up to 4 probes
% simultanously.

switch ProbeLoc
    case 'A'

shotlist = 12000 + [193,201,205,207];
shotlist = 12000 + [207:211];
shotlist = 12000 + [222,229];
% shotlist = [12193,12201,12205,12207,12208,12209,12210,12211,12222];
% shotlist = [12205];

DLPType = '9';
Config.FitFunction = 2; 

shotlist = 15658;

DLPType = '10MP';
Config.FitFunction = 2; 

ChannelType = '1' ; % "1 for TARGET_LP,LP_V_RAMP, "2" for other options
Config.tStart = 4.18; % [s]
Config.tEnd = 4.5;
SweepType = 'niso'; % "iso" for isolated sweep, "niso" for non-isolated sweep
AttType = 'Vx2,Ix5'; % Attenuation on the digitized signals
Config.Center_V = 1; % Remove offset on V: 1 (yes) 0(no)
Config.Center_I = 0; % Remove offset on I: 1 (yes) 0(no)

shotlist = 15658;
shotlist = 14178;

DLPType = '10MP';
Config.FitFunction = 2; 

ChannelType = '1' ; % "1 for TARGET_LP,LP_V_RAMP, "2" for other options
Config.tStart = 4.18; % [s]
Config.tEnd = 4.5;
SweepType = 'niso'; % "iso" for isolated sweep, "niso" for non-isolated sweep
AttType = 'Vx2,Ix5'; % Attenuation on the digitized signals
Config.Center_V = 1; % Remove offset on V: 1 (yes) 0(no)
Config.Center_I = 0; % Remove offset on I: 1 (yes) 0(no)


% shotlist = [12474,12477];
% DLPType = '9';
% Config.FitFunction = 2; 
% 
% ChannelType = '1' ; % "1 for TARGET_LP,LP_V_RAMP, "2" for other options
% Config.tStart = 4.18; % [s]
% Config.tEnd = 4.32;
% SweepType = 'niso'; % "iso" for isolated sweep, "niso" for non-isolated sweep
% AttType = 'Vx2,Ix5'; % Attenuation on the digitized signals
% Config.Center_V = 1; % Remove offset on V: 1 (yes) 0(no)
% Config.Center_I = 0; % Remove offset on I: 1 (yes) 0(no)

    case 'B'
        % No voltage data
        shotlist = 14600 + [64];

DLPType = '1MP'; % DLP 1.5 horizontal
ChannelType = '4' ; % "1" for TARGET_LP,LP_V_RAMP, "2" for other options
Config.FitFunction = 2; 

Config.tStart = 4.15; % [s]
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
        Config.L_tip = 1.2/1000;
        Config.D_tip = 0.254/1000; % [m]        
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
    case 'TDLP'
        DLP = 11.5;
        Config.L_tip = 0/1000;
        Config.D_tip = 0.8/1000; % [m]
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
Config.AreaType = 1; % 1: Cylindrical + cap
% #########################################################################
[ni,Te,time,Ifit,Ip,Vp,tm,Vsweep,Isweep,GlitchFlag,SSQres,StdRes,StdResNorm]...
    = DLP_fit_V5_3(Config,shotlist,DataAddress);
% [ni,Te,time,Ifit,Ip,Vp,tm,Vsweep,Isweep,GlitchFlag,SSQres] = DLP_fit_V5_2(Config,shotlist,DataAddress);
% [ni,Te,time,Ifit,Ip,Vp,tm,Vsweep,Isweep] = DLP_fit_V5(Config,shotlist,DataAddress);

%Ni,Te,time,Ifit,Ip,Vp,tm,Vsweep,Isweep,GlitchFlag,SSQres
for s = 1:length(shotlist)
    Ni{s} = 0.5*(ni{s}{1} + ni{s}{2});
end
% #########################################################################
DA{1} = [RootAddress,'PWR_28GHZ']; % Isx
[ECH,t_ech]   = my_mdsvalue_v2(shotlist,DA(1))
% #########################################################################
DA{1} = [RootAddress,'RF_FWD_PWR'];
[RF,t_rf]   = my_mdsvalue_v2(shotlist,DA(1))
% #########################################################################


%%
% #########################################################################
% PLOT DATA:
% #########################################################################

close all
% C = {'k.','r.','bl.','g.','m.','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:'};
C = {'k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:'};
t = t_zero(shotlist);

TimePlotEnd = 4.5;
try
    clear h
end

figure; 
subplot(2,2,1); hold on
for s = 1:length(shotlist)
    condition = abs(StdResNorm{s})<10/100 & StdRes{s} == real(StdRes{s});
    PlotRng{s} = find(condition);
    noPlotRng{s} = find(~condition);

    if ~isempty(PlotRng{s})
            h(s) = plot(time{s}(PlotRng{s}),Ni{s}(PlotRng{s}),C{s},'lineWidth',1);
            set(h(s),'Marker','o','MarkerEdgeColor',C{s}(1),'MarkerSize',2)
%             plot(time{s}(PlotRng{s}),Ni{s}(PlotRng{s}),'k-','lineWidth',2);
            %      plot(time{s},ni{s}{1},C{s},'lineWidth',1);
            %      plot(time{s},ni{s}{2},C{s},'lineWidth',1);
            legend(h,['DLP ',num2str(DLP)],'location','SouthEast')
    end
    hold on
    plot(t_ech{s}(1:length(ECH{s})),ECH{s}*0.5e19,C{s},'lineWidth',0.5)
    plot(t_rf{s}(1:length(RF{s})),(RF{s}.^2)*6e19,C{s},'lineWidth',0.5)

end
title('$ n_e $ $ [m^{-3}] $','interpreter','Latex','FontSize',13,'Rotation',0)

%set(gca,'PlotBoxAspectRatio',[1 1 1])
ylim([0,8e19])
xlim([4.15,TimePlotEnd])

subplot(2,2,2); hold on
for s = 1:length(shotlist)
    if ~isempty(PlotRng{s})
            h(s) = plot(time{s}(PlotRng{s}),Te{s}(PlotRng{s}),C{s},'lineWidth',2)
            L{s} = [num2str(shotlist(s)),' ,t=',num2str(t{s}(10:14))];
            legend(h,L,'location','NorthEast')
    end

end
title('$ T_e $ $ [eV] $','interpreter','Latex','FontSize',13,'Rotation',0)
ylim([0,6])
xlim([4.15,TimePlotEnd])

subplot(2,2,3); hold on
for s = 1:length(shotlist)
        if ~isempty(PlotRng{s})
            Pe = e_c.*Ni{s}.*Te{s};
            plot(time{s}(PlotRng{s}),Pe(PlotRng{s}),C{s},'lineWidth',2)
        end
end
ylim([0,30])
xlim([4.15,TimePlotEnd])
title('$ P_e $ $ [Pa] $','interpreter','Latex','FontSize',13,'Rotation',0)

set(findobj('-Property','YTick'),'box','on')
set(gcf,'color','w')

if FitShow
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

if FitShow
    
     figure;
    for c = 1:25;
        subplot(5,5,c); hold on
        plot(Vsweep{s}{c},Isweep{s}{c}*1e3,'k')
        plot(Vsweep{s}{c},Ifit{s}{c}*1e3,'r')
        ht = title(['T_e: ',num2str(Te{s}(c),2),' , t : ',num2str(time{s}(c),3)]);
        set(ht,'FontSize',5); grid on
    end
    
    figure;
        for c = 26:50;
        subplot(5,5,c-25); hold on
        plot(Vsweep{s}{c},Isweep{s}{c}*1e3,'k')
        plot(Vsweep{s}{c},Ifit{s}{c}*1e3,'r')
        ht = title(['T_e: ',num2str(Te{s}(c),2),' , t : ',num2str(time{s}(c),3)]);
        set(ht,'FontSize',5); grid on
        end
        
    figure;
        for c = 51:75;
        subplot(5,5,c-50); hold on
        plot(Vsweep{s}{c},Isweep{s}{c}*1e3,'k')
        plot(Vsweep{s}{c},Ifit{s}{c}*1e3,'r')
        ht = title(['T_e: ',num2str(Te{s}(c),2),' , t : ',num2str(time{s}(c),3)]);
        set(ht,'FontSize',5); grid on
        end  
        
        figure;
        for c = 76:100;
        subplot(5,5,c-75); hold on
        plot(Vsweep{s}{c},Isweep{s}{c}*1e3,'k')
        plot(Vsweep{s}{c},Ifit{s}{c}*1e3,'r')
        ht = title(['T_e: ',num2str(Te{s}(c),2),' , t : ',num2str(time{s}(c),3)]);
        set(ht,'FontSize',5); grid on
        end      
    
end