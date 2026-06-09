clear all
close all

ProbeLoc = 'B';
ProbeLoc = 'A';
FitShow = 1; 
% Since we now have several probe drives, we can run up to 4 probes
% simultanously.

switch ProbeLoc
    case 'A'
% =========================================================================
% First few shots at the MAB conditions
shotlist = [17900 + [92,93,95,97,98,100]];
shotlist = [17900 + [100,101,103]];
SweepType = 'niso'; % "iso" for isolated sweep, "niso" for non-isolated sweep
AttType = 'Vx2,Ix5'; % Attenuation on the digitized signals
% =========================================================================
shotlist = [17900 + [104,105,107]];
SweepType = 'iso'; % "iso" for isolated sweep, "niso" for non-isolated sweep
AttType = 'Vx1,Ix1'; % Attenuation on the digitized signals
% =========================================================================
% shotlist = [17900 + [106]];
% SweepType = 'niso'; % "iso" for isolated sweep, "niso" for non-isolated sweep
% AttType = 'Vx2,Ix5'; % Attenuation on the digitized signals
% =========================================================================
% BEGIN DLP 10.5 RADIAL SCAN, NO ECH
shotlist = [17900 + [107,108,109,110,111,112,113]];
shotlist = [17900 + [ 110, 111, 112, 113, 114, 116, 117, 118, 122, 124, 126]];
r =                 [ 0.5, 0.0,-0.5,-1.0,-1.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5] ;
SweepType = 'iso'; % "iso" for isolated sweep, "niso" for non-isolated sweep
AttType = 'Vx1,Ix1'; % Attenuation on the digitized signals
% =========================================================================
% % % BEGIN DLP 10.5 RADIAL SCAN, with ECH
 shotlist = [17900 + [121, 123, 125, 127, 128, 130, 131, 132, 133, 134, 135, 136]];
r =                  [2.0, 2.5, 3.0, 3.5, 1.5, 1.0, 0.5, 0.0,-0.5,-1.0,-1.5,-2.0] ;
SweepType = 'iso'; % "iso" for isolated sweep, "niso" for non-isolated sweep
AttType = 'Vx1,Ix1'; % Attenuation on the digitized signals
% % =========================================================================
% % =========================================================================
% % BEGIN DLP 10.5 RADIAL SCAN, with w/ and w/o ECH
 shotlist = [17900 + [121, 123, 125, 127, 128, 130, 131, 132, 133, 134, 135, 136, 110, 111, 112, 113, 114, 116, 117, 118, 122, 124, 126]];
r =                  [2.0, 2.5, 3.0, 3.5, 1.5, 1.0, 0.5, 0.0,-0.5,-1.0,-1.5,-2.0, 0.5, 0.0,-0.5,-1.0,-1.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5] ;
ECH_status =         [  1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,    0,  0,   0,   0,   0,   0,   0,   0,   0,   0,   0] ;
SweepType = 'iso'; % "iso" for isolated sweep, "niso" for non-isolated sweep
AttType = 'Vx1,Ix1'; % Attenuation on the digitized signals
% % % =========================================================================
DLPType = '10MP';
Config.FitFunction = 2; 
ChannelType = '1' ;

Config.tStart = 4.16; % [s]
Config.tEnd = 4.53;
Config.Center_V = 1; % Remove offset on V: 1 (yes) 0(no)
Config.Center_I = 1; % Remove offset on I: 1 (yes) 0(no)
    
case 'B'
% =========================================================================
% =========================================================================
        
Config.FitFunction = 2; 

DLPType = '1MP'; % DLP 1.5 horizontal
ChannelType = '2' ; % "1" for TARGET_LP,LP_V_RAMP, "2" for other options

Config.tStart = 4.16; % [s]
Config.tEnd = 4.45;
SweepType = 'iso'; % "iso" for isolated sweep, "niso" for non-isolated sweep
AttType = 'Vx1,Ix1'; % Attenuation on the digitized signals
Config.Center_V = 0; % Remove offset on V: 1 (yes) 0(no)
Config.Center_I = 0; % Remove offset on I: 1 (yes) 0(no)
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
        Config.L_tip = 1.01/1000; % 
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
        Config.L_tip = 1.6/1000;% according to Nischal 2017_11_21
        Config.D_tip = 0.254/1000; % [m]
    case '4MP'
        DLP = 4.5;
        Config.L_tip = 1.7/1000;
        Config.D_tip = 0.254/1000; % [m]
    case '1MP'
        DLP = 2.5;
        Config.L_tip = 1.7/1000;
        Config.D_tip = 0.254/1000; % [m]
    case 'TDLP'
        DLP = 11.5;
        Config.L_tip = 0/1000;
        Config.D_tip = 0.8/1000; % [m]
    case 'IF'
        DLP = 11.5;
        Config.L_tip = 0/1000;
        Config.D_tip = 1/1000; % [m]
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
        DataAddress{1} = [RootAddress,'INT_4MM_2']; % V
        DataAddress{2} = [RootAddress,'INT_4MM_1']; % I
    case '3'
        DataAddress{1} = [RootAddress,'GEN_RF_PWR']; % V
        DataAddress{2} = [RootAddress,'INT_4MM_2'];  % I
    case '4'
        DataAddress{1} = [RootAddress,'LP_V_RAMP']; % V
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
% [ni,Te,time,Ifit,Ip,Vp,tm,Vsweep,Isweep,GlitchFlag,SSQres] = DLP_fit_V5_2(Config,shotlist,DataAddress);
[ni,Te,time,Ifit,Ip,Vp,tm,Vsweep,Isweep,GlitchFlag,SSQres,StdRes,StdResNorm] = DLP_fit_V5_4(Config,shotlist,DataAddress);

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

% close all
GasCompareALL_2017_11_29
C = {'k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:'};
t = t_zero(shotlist);

TimePlotEnd = 4.7;

figure; 
subplot(2,2,1); hold on
for s = 1:length(shotlist)
    GoodFits{s} = GlitchFlag{s} == 0 & StdResNorm{s}<=0.15 & Ni{s}>0 & Ni{s}<1e21 & Te{s}<=15;
%     GoodFits{s} = GlitchFlag{s} == 0 & StdResNorm{s}<=0.15 & Ni{s}>0 & Ni{s}<1e21 & Te{s}<=10;

%      plot(time{s},ni{s}{1},C{s},'lineWidth',1);
%      plot(time{s},ni{s}{2},C{s},'lineWidth',1);
    h(s) = plot(time{s}(GoodFits{s}),Ni{s}(GoodFits{s}),C{s},'lineWidth',2);

    hold on
    plot(t_ech{s}(1:length(ECH{s})),ECH{s}*1e19,C{s},'lineWidth',0.5)
    plot(t_rf{s}(1:length(RF{s})),(RF{s}.^2)*15e19,C{s},'lineWidth',0.5)

end
title('$ n_e $ $ [m^{-3}] $','interpreter','Latex','FontSize',13,'Rotation',0)
legend(h,['DLP ',num2str(DLP)],'location','SouthEast')
%set(gca,'PlotBoxAspectRatio',[1 1 1])
ylim([0,7e19])
xlim([4.15,TimePlotEnd])

subplot(2,2,2); hold on
for s = 1:length(shotlist)
    h(s) = plot(time{s}(GoodFits{s}),Te{s}(GoodFits{s}),C{s},'lineWidth',2)
    L{s} = [num2str(shotlist(s)),' ,t=',num2str(t{s}(10:14))];
end
legend(h,L,'location','NorthEast')
title('$ T_e $ $ [eV] $','interpreter','Latex','FontSize',13,'Rotation',0)
ylim([0,6])
xlim([4.15,TimePlotEnd])

subplot(2,2,3); hold on
for s = 1:length(shotlist)
    Pe{s} = e_c.*Ni{s}.*Te{s};
    plot(time{s}(GoodFits{s}),Pe{s}(GoodFits{s}),C{s},'lineWidth',2)
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
    h(s) = plot(tm{s},(Ip{s})*1000);
    ylim([-100,100])
    xlim([4.15,TimePlotEnd])
    grid on
    end
    legend(h,num2str(shotlist'))
    set(gcf,'color','w')
    box on
end

%%
% close all
Ymax = 150; % mA
if FitShow
    for k = 1:8
         figure;
    for c = (1 + (k-1)*25):(25 + (k-1)*25);
        subplot(5,5,c- (k-1)*25 ); hold on
        try
        plot(Vsweep{s}{c},Isweep{s}{c}*1e3,'k')
        if GoodFits{s}(c) == 1
            plot(Vsweep{s}{c},Ifit{s}{c}*1e3,'r')
        else
            plot(Vsweep{s}{c},Ifit{s}{c}*1e3,'g')
        end
        ht = title(['T_e: ',num2str(Te{s}(c),2),' ,t = ',num2str(time{s}(c),5),' s']);
        set(ht,'FontSize',5); grid on
%         ylim(Ymax*[-1,1])
        catch
            warning('error')
            continue
        end
    end
    end
end
%%
% #########################################################################
% Steady state values
% #########################################################################

% close all

% Before ECH ==============================================================
for s = 1:length(shotlist)
    if shotlist(s) == 13611 || shotlist(s) == 13615
        rng0 = time{s}>=4.3 & time{s}<=4.35;
    else
        rng0 = time{s}>=4.3 & time{s}<=4.35;
    end
        ne0(s)  = mean(Ni{s}(rng0 & GoodFits{s}));
        dne0(s) = std(Ni{s}(rng0 & GoodFits{s}),1,2);
        Te0(s)  = mean(Te{s}(rng0 & GoodFits{s}));
        dTe0(s) = std(Te{s}(rng0 & GoodFits{s}),1,2);
end

% During ECH: =============================================================
for s = 1:length(shotlist)
%     if shotlist(s) == 13611 || shotlist(s) == 13615
%         rng0 = find(time{s}>=4.25 & time{s}<=4.26);
    if shotlist(s) == 13891
        rng1 = time{s}>=4.45 & time{s}<=4.53;
    else
        rng1 = time{s}>=4.45 & time{s}<=4.53;
    end
        ne1(s)  = mean(Ni{s}(rng1 & GoodFits{s}));
        dne1(s) = std(Ni{s}(rng1 & GoodFits{s}),1,2);
        Te1(s)  = mean(Te{s}(rng1 & GoodFits{s}));
        dTe1(s) = std(Te{s}(rng1 & GoodFits{s}),1,2);
end

[a,b] = sort(r);

figure; 
subplot(2,1,1); hold on
h = errorbar(r(b),ne1(b),dne1(b),'ro-','LineWidth',1)
h = errorbar(r(b),ne0(b),dne0(b),'ko-','LineWidth',1); 
ylim([0,8]*1e19)
xlim([-4,4]); 
box on
t = title('$n_e$ $[m^{-3}]$ '); t.Interpreter = 'latex'; t.FontSize = 13;

subplot(2,1,2); hold on
h = errorbar(r(b),Te1(b),dTe1(b),'ro-','LineWidth',1)
h = errorbar(r(b),Te0(b),dTe0(b),'ko-','LineWidth',1)
ylim([0,10])
xlim([-4,4]); 
box on
t = title('$T_e$ $[eV]$ '); t.Interpreter = 'latex'; t.FontSize = 13;

set(gcf,'color','w')

% Save data to Excel spreadsheet:
% Create table to export to EXCEL
if 0
G = [shotlist(b)',r(b)'-offset,ne0(b)',dne0(b)',Te0(b)',dTe0(b)'];
F = {'Shot','R [cm]','ne[m^-3]','dne','Te[eV]','dTe'};
FileName = 'ECH_OFF_NeTe_10_5_2017_11_29.xlsx';
xlswrite(FileName,[F;num2cell(G)]);
end