clear all
close all

ProbeLoc = 'A';
FitShow = 0; 
RawDataShow = 1;
% 
shotlist = 24000 + [526,528,529,530,532,533];
PS1      =         [2.0,3.0,4.0,5.0,6.0,7.0];
PS2      =         [2.0,2.0,2.0,2.0,2.0,2.0];
DLPType = '12MP';

% PS1 SCAN
% 1st scan with PS2 = 3500 A
shotlist = 24000 + [504,503,505,506,502,508];
PS1      =         [2.0,3.0,4.0,5.0,6.0,7.0];
DLPType = '1MP';

% #########################################################################

SweepType = 'iso'; % "iso" for isolated sweep, "niso" for non-isolated sweep
AttType = 'Vx1,Ix1'; % Attenuation on the digitized signals
Config.FitFunction = 2; 
ChannelType = '6' ;
Config.tStart = 4.15; % [s]
Config.tEnd = 4.9;
Config.Center_V = 1; % Remove offset on V: 1 (yes) 0(no)
Config.Center_I = 1; % Remove offset on I: 1 (yes) 0(no)
    

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
%         Config.L_tip = 1.6/1000;% according to Nischal 2017_11_21
        Config.L_tip = 0.75/1000; % As of 2018_08_27 after Magnetic field upgrade
        Config.D_tip = 0.254/1000; % [m]
    case '4MP'
        DLP = 4.5;
        Config.L_tip = 1.7/1000;
        Config.D_tip = 0.254/1000; % [m]
    case '1MP'
        DLP = 1.5;
        Config.L_tip = 1.7/1000;
        Config.L_tip = 1.1/1000; % As of 2018_08_31
        Config.L_tip = 3.0/1000; % As of 2018_09_19
        Config.D_tip = 0.254/1000; % [m]
    case 'TDLP'
        DLP = 11.5;
        Config.L_tip = 0/1000;
        Config.D_tip = 0.8/1000; % [m]
    case 'IF'
        DLP = 11.5;
        Config.L_tip = 0/1000;
        Config.D_tip = 1/1000; % [m]
    case '10FluxProbe';
        DLP = 10.5;
        Config.L_tip = 0/1000;
        Config.D_tip = ((2^(1/4))*1.22)/1000; % [m]
    case '12MP';
        DLP = 12.5;
        Config.L_tip = 1.1/1000; 
        Config.L_tip = 2.4/1000; % As of 2018_10_24 
        Config.D_tip = 0.254/1000; % [m]
   case '8MP';
        DLP = 8.5;
        Config.L_tip = 2/1000;
        Config.L_tip = 2.3/1000; % 2018_11_06
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
    case '3'
        DataAddress{1} = [RootAddress,'GEN_RF_PWR']; % V
        DataAddress{2} = [RootAddress,'INT_4MM_2'];  % I
    case '4'
        DataAddress{1} = [RootAddress,'LP_V_RAMP']; % V
        DataAddress{2} = [RootAddress,'INT_4MM_2'];  % I
    case '5'
        DataAddress{1} = [RootAddress,'LP_V_RAMP']; % V
        DataAddress{2} = [RootAddress,'EA_CURRENT'];  % I
    case '6'
        DataAddress{1} = [RootAddress,'LP_V_RAMP']; % V
        DataAddress{2} = [RootAddress,'LP_1'];  % I
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
Config.SGF = 11;
Config.TimeMode = 2; % Effective time of sweep, (1) start of ramp or (2) mean time of ramp
Config.AMU = 2; % Ion mass in AMU
Config.AreaType = 1; % 1: Cylindrical + cap
% #########################################################################
% [ni,Te,time,Ifit,Ip,Vp,tm,Vsweep,Isweep,GlitchFlag,SSQres] = DLP_fit_V5_2(Config,shotlist,DataAddress);
[ni,Te,Isat,time,Ifit,Ip,Vp,tm,Vsweep,Isweep,GlitchFlag,SSQres,StdRes,StdResNorm] = DLP_fit_V5_7(Config,shotlist,DataAddress);

%Ni,Te,time,Ifit,Ip,Vp,tm,Vsweep,Isweep,GlitchFlag,SSQres
for s = 1:length(shotlist)
    Ni{s} = 0.5*(ni{s}{1} + ni{s}{2});
    Is{s} = 0.5*(Isat{s}{1} + Isat{s}{2});
end
% #########################################################################
DA{1} = [RootAddress,'PWR_105GHZ']; % Isx
[ECH,t_ech]   = my_mdsvalue_v2(shotlist,DA(1))
% % #########################################################################
DA{1} = [RootAddress,'RF_FWD_PWR'];
[RF,t_rf]   = my_mdsvalue_v2(shotlist,DA(1))
% #########################################################################

%%
% #########################################################################
% PLOT DATA:
% #########################################################################

close all
% Add code that looks at the date and then looks for a GasCompare file with
% that name

% -------------------------------------------------------------------------
% The following script checks the existence of the "GasCompare" script and
% runs if it exists
ScriptList = dir('*.m');
for n = 1:length(ScriptList)
    if strcmp(ScriptList(n).name(1:13),'GasCompareALL');
        run(ScriptList(n).name)
        break % the script was been found
    else
        error('Could not find GasCompare script')
    end
end
% -------------------------------------------------------------------------

C = {'k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:'};
t = t_zero(shotlist);

TimePlotStart = 4.1;
TimePlotEnd = 5.0;

figure; 
subplot(2,2,1); hold on
for s = 1:length(shotlist)
    GoodFits{s} = GlitchFlag{s} == 0 & StdResNorm{s}<=0.6 & Ni{s}>0 & Ni{s}<30e19 & Te{s}<=20 & Is{s}>0.01e-3;
      
%     plot(time{s}(GoodFits{s}),ni{s}{1}(GoodFits{s}),C{s},'lineWidth',1);
%     plot(time{s}(GoodFits{s}),ni{s}{2}(GoodFits{s}),C{s},'lineWidth',1);
      
    t_Ni_gf{s} = time{s}(GoodFits{s});
    Ni_gf{s} = Ni{s}(GoodFits{s});
    
    h(s) = plot(t_Ni_gf{s},Ni_gf{s},C{s},'lineWidth',2);

    hold on
    plot(t_ech{s}(1:length(ECH{s})),abs(ECH{s})*1e19,C{s},'lineWidth',0.5)
    
    if s == 2
        continue
    end
     plot(t_rf{s}(1:length(RF{s})),(RF{s}.^2)*10e19,C{s},'lineWidth',0.5)

end
title('$ n_e $ $ [m^{-3}] $','interpreter','Latex','FontSize',13,'Rotation',0)
legend(h,['DLP ',num2str(DLP)],'location','SouthEast')
%set(gca,'PlotBoxAspectRatio',[1 1 1])
% ylim([0,7e18])
% ylim([0,1e19])
ylim([0,1.0e20])
xlim([TimePlotStart,TimePlotEnd])

subplot(2,2,2); hold on
for s = 1:length(shotlist)
        t_Te_gf{s} = time{s}(GoodFits{s});
        Te_gf{s} = sgolay_t(Te{s}(GoodFits{s}),3,7);
        
        h(s) = plot(t_Te_gf{s},Te_gf{s},C{s},'lineWidth',2);
        
        L{s} = [num2str(shotlist(s)),' ,t=',num2str(t{s}(10:14))];
       plot(t_ech{s}(1:length(ECH{s})),2*abs(ECH{s}),C{s},'lineWidth',0.5)

end
legend(h,L,'location','NorthEast')
title('$ T_e $ $ [eV] $','interpreter','Latex','FontSize',13,'Rotation',0)
ylim([0,17])
xlim([TimePlotStart,TimePlotEnd])

subplot(2,2,3); hold on
for s = 1:length(shotlist)
    Pe_gf{s} = e_c.*Ni_gf{s}.*Te_gf{s};
    plot(t_Ni_gf{s},Pe_gf{s},C{s},'lineWidth',2)
end
ylim([0,30])
xlim([TimePlotStart,TimePlotEnd])
title('$ P_e $ $ [Pa] $','interpreter','Latex','FontSize',13,'Rotation',0)

set(findobj('-Property','YTick'),'box','on')
set(gcf,'color','w')

if RawDataShow
    figure; hold on
    for s = 1:length(shotlist)
    plot(tm{s},Vp{s})
    h(s) = plot(tm{s},(Ip{s})*1000);
    ylim([-300,300])
    xlim([4.15,TimePlotEnd])
    grid on
    end
    legend(h,num2str(shotlist'))
    set(gcf,'color','w')
    box on
end

%%
% #########################################################################
% Steady state values
% #########################################################################
close all

for s = 1:length(shotlist) 
    rngT = [time{s}>=4.635 & time{s}<=4.640]; 
    nes(s)    = mean(Ni{s}(rngT & GoodFits{s}));
    dnes(s)   = std(Ni{s}(rngT & GoodFits{s}),1,2);
    Tes(s)    = mean(Te{s}(rngT & GoodFits{s}));
    dTes(s)   = std(Te{s}(rngT & GoodFits{s}),1,2);
    [M,I] = max(GasPressure{1}{s});
    MaxPressure(s)= M;
end

R = PS1*30/7;

figure; 
subplot(1,2,2);
hold on
plot(R,nes,'ko-')
switch DLPType
    case '12MP'
        ylim([0,4e19])
    case '1MP'
        ylim([0,10e19])
end
xlim([0,35])
axis('square')
box on
set(gca,'FontName','Times','FontSize',12)
ylabel('$n_e$ $[m^-3]$','Interpreter','latex')
xlabel('R','Interpreter','latex')

subplot(1,2,1);
plot(R,MaxPressure,'ksq-')
ylim([0,0.7])
xlim([0,35])
axis('square')
box on
set(gcf,'color','w')
set(gca,'FontName','Times','FontSize',12)
ylabel('$P_{D_2}$ [Pa]','Interpreter','latex')
xlabel('R','Interpreter','latex')

figure; 
hold on
plot(R,Tes,'ko-')
switch DLPType
    case '12MP'
        ylim([0,5])
    case '1MP'
        ylim([0,5])
end
xlim([0,35])
axis('square')
box on
set(gca,'FontName','Times','FontSize',12)
ylabel('$T_e$ $[eV]$','Interpreter','latex')
xlabel('R','Interpreter','latex')

%% Convert data into Excel
% Save data to Excel spreadsheet:
% Create table to export to EXCEL
if 0
G = [shotlist(b)',PS1(b)',ne{p}(b)',dne{p}(b)',Te{p}(b)',dTe{p}(b)'];
F = {'Shot','R [cm]','ne[m^-3]','dne','Te[eV]','dTe'};
FileName = 'NeTe_Spool_10_5_2018_01_16.xlsx';
xlswrite(FileName,[F;num2cell(G)]);
end