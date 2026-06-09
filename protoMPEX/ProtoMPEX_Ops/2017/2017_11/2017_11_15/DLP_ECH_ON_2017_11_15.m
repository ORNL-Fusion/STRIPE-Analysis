close all 
clear all

shotlist = 17600 + [69 , 70,  73, 104, 105, 106, 108, 109, 110, 167, 168, 169, 172, 174, 175, 176];
r        =         [0.5,0.0,-0.5,-1.0,-1.5,-2.0,-2.5,-3.0,-3.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 4.5];

% shotlist = 17600 + [ 172, 174, 175, 176];

Config.tStart = 4.3;
Config.tEnd = 4.50;
Config.tEnd = 4.46;

% shotlist = 17774
% Config.tStart = 4.15;
% Config.tEnd = 4.67;

DLPType = '10MP'; % DLP-MP 10.5 vertical
SweepType = 'niso'; % "iso" for isolated sweep, "niso" for non-isolated sweep
AttType = 'Vx2,Ix5'; % Attenuation on the digitized signals
ChannelType = '1' ; % "1" for TARGET_LP,LP_V_RAMP, "2" for other options
Config.Center_V = 1; % Remove offset on V: 1 (yes) 0(no)
Config.Center_I = 1; % Remove offset on I: 1 (yes) 0(no)

%%
% #########################################################################
% DLP FITTING ROUTINE
% #########################################################################
switch DLPType 
    case '4'
        DLP = 4.5;
        Config.L_tip = 1.1/1000; % is it 1.0 or 1.2 mm ?
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
    case '2MP'
        DLP = 2.5;
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
        DataAddress{1} = [RootAddress,'INT_4MM_2']; % V
        DataAddress{2} = [RootAddress,'INT_4MM_1']; % I
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
% #########################################################################
% DA{1} = [RootAddress,'PWR_28GHZ']; 
% [ECH,t_ech]   = my_mdsvalue_v2(shotlist,DA(1))
% #########################################################################
DA{1} = [RootAddress,'RF_FWD_PWR'];
[RF,t_rf]   = my_mdsvalue_v2(shotlist,DA(1));
% #########################################################################
%%
% #########################################################################
% PLOT DATA:
% #########################################################################

close all
% GasCompareALL_2017_11_09
C = {'k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:'};
t = t_zero(shotlist);

TimePlotEnd = 4.7;

figure; 
subplot(2,2,1); hold on
for s = 1:length(shotlist)
    nGoodFit{s} = find(SSQres{s}<0.03 & ~GlitchFlag{s});
    
    plot(time{s}(nGoodFit{s}),ni{s}{1}(nGoodFit{s}),C{s},'lineWidth',1);
    plot(time{s}(nGoodFit{s}),ni{s}{2}(nGoodFit{s}),C{s},'lineWidth',1);
    h(s) = plot(time{s}(nGoodFit{s}),Ni{s}(nGoodFit{s}),C{s},'lineWidth',2);
    hold on
    plot(t_rf{s}(1:length(RF{s})),(RF{s}.^2)*6e19,C{s},'lineWidth',0.5)

end
title('$ n_e $ $ [m^{-3}] $','interpreter','Latex','FontSize',13,'Rotation',0)
legend(h,['DLP ',num2str(DLP)],'location','SouthEast')
%set(gca,'PlotBoxAspectRatio',[1 1 1])
ylim([0,6e19])
xlim([4.15,TimePlotEnd])

subplot(2,2,2); hold on
for s = 1:length(shotlist)
    h(s) = plot(time{s}(nGoodFit{s}),Te{s}(nGoodFit{s}),C{s},'lineWidth',2);
    
    L{s} = [num2str(shotlist(s)),' ,t=',num2str(t{s}(10:14))];
end
legend(h,L,'location','NorthEast')
title('$ T_e $ $ [eV] $','interpreter','Latex','FontSize',13,'Rotation',0)
ylim([0,10])
xlim([4.15,TimePlotEnd])

subplot(2,2,3); hold on
for s = 1:length(shotlist)
        plot(time{s}(nGoodFit{s}),e_c.*Ni{s}(nGoodFit{s}).*Te{s}(nGoodFit{s}),C{s},'lineWidth',2)
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
    
end

%%
% #########################################################################
% Steady state values
% #########################################################################

% close all
for s = 1:length(shotlist)
        rng0 = find(time{s}(nGoodFit{s})>=4.3 & time{s}(nGoodFit{s})<=4.4);
        ne0(s)  = mean(Ni{s}(nGoodFit{s}(rng0)));
        dne0(s) = std(Ni{s}(nGoodFit{s}(rng0)),1,2);
        Te0(s)  = mean(Te{s}(nGoodFit{s}(rng0)));
        dTe0(s) = std(Te{s}(nGoodFit{s}(rng0)),1,2);
end

[a,b] = sort(r);
offset = 0.5;

figure; 
subplot(2,1,1); hold on
h = errorbar(r(b)-offset,ne0(b),dne0(b),'ko-','LineWidth',1); 
ylim([0,3]*1e19)
xlim([-6,6]); 
box on
t = title('$n_e$ $[m^{-3}]$ '); t.Interpreter = 'latex'; t.FontSize = 13;

subplot(2,1,2); hold on
h = errorbar(r(b)-offset,Te0(b),dTe0(b),'ko-','LineWidth',1);
ylim([0,8]);
xlim([-6,6]); 
box on
t = title('$T_e$ $[eV]$ '); t.Interpreter = 'latex'; t.FontSize = 13;

set(gcf,'color','w')


figure; 
plot(r(b)-offset,Te0(b).*ne0(b),'ko-')


% Save data to Excel spreadsheet:
% Create table to export to EXCEL
if 0
G = [shotlist(b)',r(b)'-offset,ne0(b)',dne0(b)',Te0(b)',dTe0(b)'];
F = {'Shot','R [cm]','ne[m^-3]','dne','Te[eV]','dTe'};
FileName = 'ECH_ON_NeTe_10_5_2017_11_15.xlsx';
xlswrite(FileName,[F;num2cell(G)]);
end

