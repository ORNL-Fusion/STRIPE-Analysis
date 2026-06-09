clear all
close all

ProbeLoc = 'A';
FitShow = 0; 
RawDataShow = 1;

% #########################################################################
% 2019_08_08
% DLP 12.5
shotlist = [27000 + [356,357,358]];
shotlist = [27000 + [357,361,362,363,364,365,366]];
shotlist = [27000 + [366,367,368],26827];
% Best shots of the day:
shotlist = [27000 + [357,366,369]];

% Start with IR helicon window experimetns
% Comparing two cases with different PS1 current 
shotlist = [27000 + [368,371]];

% reference conditions for IR helicon window experimetns
% Bottom view
shotlist = [27000 + [372,373,374,375,376,377,380]];
Coil2    =          [520,420,320,220,120, 50,  0];

% shotlist = [27000 + [377,379,380]];
% Coil2    =          [50 ,50 ,0  ];

% Top view
shotlist = [27000 + [381,382,384,383,385,386,387,389]];
Coil2    =          [0  ,120,220,320,420,520,720,920];

% Increasing PS3
shotlist = [27000 + [389,391,393,394]];
PS3      =          [220,250,280,220];


% Start woth downstream  puffing as per shot 15947 on 2017_08_02
shotlist = [27000 + [398,399,400,401,404,405]];
PS3      =          [220,220,220,250,280,310];

shotlist = [27000 + [401,410,411,412,413]];
PS3      =          [250,310,250,250,250];

Config.AMU = 2; % Ion mass in AMU
DLPType = '12MP';

SweepType = 'iso'  ;AttType = 'Vx1,Ix1';
% #########################################################################

Config.FitFunction = 2; 
ChannelType = '6' ;
Config.tStart = 4.15; % [s]
Config.tEnd = 4.8;
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
    case '6MP'
        DLP = 6.5;
        Config.L_tip = 1.0/1000; % 1.0 as of April 11th 2017
        Config.L_tip = 1.4/1000; % 1.4 as of Nov 29th 2018
        Config.L_tip = 1.7/1000; % 1.7 as of Dec 13th 2018
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
        Config.L_tip = 2.8/1000; % As of 2018_09_19
        Config.L_tip = 2.4/1000; % As of 2018_12_21
        Config.L_tip = 2.25/1000; % As of 2019_02_07
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
        Config.L_tip = 2/1000; % As of 2019_08_08 
        Config.D_tip = 0.254/1000; % [m]
   case '8MP';
        DLP = 8.5;
        Config.L_tip = 2/1000;
        Config.L_tip = 2.2/1000; % 2018_11_06
        Config.L_tip = 1.5/1000; % 2018_12_04
        Config.D_tip = 0.254/1000; % [m]
end
% #########################################################################
switch AttType
    case 'Vx1,Ix5'
        Config.V_Att = 1;  % Output voltage of DLP box (Voltage) = V_Att*Digitized data 
        Config.I_Att = 5;  % Output voltage of DLP box (Current) = I_att*Digitized data 
    case 'Vx1,Ix2'
        Config.V_Att = 1;  % Output voltage of DLP box (Voltage) = V_Att*Digitized data 
        Config.I_Att = 2;  % Output voltage of DLP box (Current) = I_att*Digitized data    
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
DA{1} = [RootAddress,'RF_CURRENT'];
[RF,t_rf]   = my_mdsvalue_v2(shotlist,DA(1))
% #########################################################################
DA{1} = [RootAddress,'ICH_FWD_PWR'];
[ICH_F,t_ICH_F]   = my_mdsvalue_v2(shotlist,DA(1))
DA{1} = [RootAddress,'ICH_REF_PWR'];
[ICH_R,t_ICH_R]   = my_mdsvalue_v2(shotlist,DA(1))
%%
% #########################################################################
% PLOT DATA:
% #########################################################################

% close all
% Add code that looks at the date and then looks for a GasCompare file with
% that name

% -------------------------------------------------------------------------
% The following script checks the existence of the "GasCompare" script and
% runs if it exists
close all

ScriptList = dir('*.m');
for n = 1:length(ScriptList)
    if strcmp(ScriptList(n).name(1:13),'GasCompareALL');
        run(ScriptList(n).name)
        break % the script was been found
    else
%         error('Could not find GasCompare script')
    end
end
% -------------------------------------------------------------------------

C = {'k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:'};
t = t_zero(shotlist);

TimePlotStart = 4.14;
TimePlotEnd = 5.3;

figure; hold on
for s = 1:length(shotlist)
    Vpwr = RF{s};
    RFPWR = 10.^( ((27*Vpwr)-32.54)/10 + 7.2);
    plot(t_rf{s}(1:end-1),RFPWR*1e-3,C{s})
end

figure; 
subplot(2,2,1); hold on
for s = 1:length(shotlist)
    GoodFits{s} = GlitchFlag{s} == 0 & StdResNorm{s}<=0.5 & Ni{s}>0 & Ni{s}<30e19 & Te{s}<=30 & Is{s}>0.01e-3;
      
    t_Ni_gf{s} = time{s}(GoodFits{s});
    Ni_gf{s} = Ni{s}(GoodFits{s});
    
    h(s) = plot(t_Ni_gf{s},Ni_gf{s},C{s},'lineWidth',2);

    hold on
%     plot(t_ech{s}(1:length(ECH{s})),abs(ECH{s})*1e19,C{s},'lineWidth',0.5)
    
     plot(t_rf{s}(1:length(RF{s})),(RF{s}.^2)*3e19,C{s},'lineWidth',1.5)
     plot(t_ICH_F{s}(1:length(ICH_F{s})),(ICH_F{s}.^2)*3e19,C{s},'lineWidth',1.5)
%      plot(t_ICH_R{s}(1:length(ICH_R{s})),(ICH_R{s}.^2)*3e19,C{s},'lineWidth',1.5)

end
title('$ n_e $ $ [m^{-3}] $','interpreter','Latex','FontSize',13,'Rotation',0)
legend(h,['DLP ',num2str(DLP)],'location','SouthEast')
%set(gca,'PlotBoxAspectRatio',[1 1 1])
% ylim([0,7e18])
% ylim([0,1e19])
ylim([0,2e20])
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
ylim([0,6])
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

%% Convert data into Excel
% Save data to Excel spreadsheet:
% Create table to export to EXCEL
if 0
    for s = 1:length(shotlist)
        G = [t_Ni_gf{s}',Ni_gf{s}',Te_gf{s}'];
        F = {'time [s]','ne[m^-3]','Te[eV]'};
        FileName = ['Shot_',num2str(shotlist(s)),'_NeTe_Spool_10_5.xlsx']
        xlswrite(FileName,[F;num2cell(G)]);
    end
end

%%
% close all

Ymax = 150; % mA

if FitShow
    for k = 1:10
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
%         ht = title(['I_s: ',num2str(Is{s}(c),2),' ,t = ',num2str(time{s}(c),5),' s']);
        set(ht,'FontSize',6); grid on
        set(gcf,'color','w')
        set(gca,'FontSize',5)
        box on
%         ylim(Ymax*[-1,1])
        catch
            warning('error')
            continue
        end
    end
    end
end

return
%%
% #########################################################################
% Steady state values
% #########################################################################
% close all

for s = 1:length(shotlist)
    rngT{1} = [time{s}>=4.20 & time{s}<=4.23]; 
    rngT{2} = [time{s}>=4.44 & time{s}<=4.48];      
    rngT{2} = [time{s}>=4.4 & time{s}<=4.5]; 
    rngT{3} = [time{s}>=4.6 & time{s}<=4.63]; 
for p = 1:3
        nes{p}(s)    = mean(Ni{s}(rngT{p} & GoodFits{s}));
        dnes{p}(s)   = std(Ni{s}(rngT{p} & GoodFits{s}),1,2);
        Tes{p}(s)    = mean(Te{s}(rngT{p} & GoodFits{s}));
        dTes{p}(s)   = std(Te{s}(rngT{p} & GoodFits{s}),1,2);
    end
end

if 1
    figure; 
    subplot(2,1,1); hold on
    errorbar(R,nes{3},dnes{3},'ko')
    errorbar(R,nes{2},dnes{2},'ro')

    ylim([0,5e19])
    
    subplot(2,1,2); hold on
    errorbar(R,Tes{3},dTes{3},'ko')
    errorbar(R,Tes{2},dTes{2},'ro')
    ylim([0,6])
end

rOffset = 6.5;
rngECH_ON  = find(ECHpwr == 1);
rngECH_OFF = find(ECHpwr == 0);
[a1,b1] = sort(R(rngECH_ON ));
[a0,b0] = sort(R(rngECH_OFF));

figure; 
subplot(2,1,1); hold on
for p = 2:3
y(1) = errorbar(R(rngECH_ON(b1))-rOffset,nes{p}(rngECH_ON(b1)),dnes{p}(rngECH_ON(b1)),'LineWidth',1,'color','r','marker','o')
% y(2) = errorbar(R(rngECH_OFF(b0))-rOffset,nes{p}(rngECH_OFF(b0)),dnes{p}(rngECH_OFF(b0)),'LineWidth',1,'color','k')
end
set(y(1),'color','k')
% legend([y(1),y(2)],'28 GHz ON','28 GHz OFF')

% h(1).Color = 'k';
% h(2).Color = 'r';
% h(3).Color = 'g';
% h(4).Color = 'bl';
% h(1).Marker = 'sq';
% h(2).Marker = 'o';
% h(3).Marker = '^';
% h(4).Marker = '*';

ylim([0,0.7]*1e19); ylabel('$n_e$ $[m^{-3}]$','Interpreter','Latex')
xlim([-4,4]); xlabel('$R$ $[cm]$','Interpreter','Latex')
box on
grid on
t = title('$n_e$ $[m^{-3}]$ '); t.Interpreter = 'latex'; t.FontSize = 13;

subplot(2,1,2); hold on
for p = 2:3
k(1) = errorbar(R(rngECH_ON(b1))-rOffset,Tes{p}(rngECH_ON(b1)),dTes{p}(rngECH_ON(b1)),'LineWidth',1,'color','r','marker','o')
% k(2) = errorbar(R(rngECH_OFF(b0))-rOffset,Tes{p}(rngECH_OFF(b0)),dTes{p}(rngECH_OFF(b0)),'LineWidth',1,'color','k')
end
set(k(1),'color','k')
% legend([k(1),k(2)],'28 GHz ON','28 GHz OFF','location','south')
ylim([0,7]); xlabel('$R$ $[cm]$','Interpreter','Latex')
xlim([-4,4]); ylabel('$T_e$ $[eV]$','Interpreter','Latex')
box on
grid on
t = title('$T_e$ $[eV]$ '); t.Interpreter = 'latex'; t.FontSize = 13;

set(gcf,'color','w','position',[360.3333  184.3333  399.3333  433.3333])


figure

HeatFlux_ECH_ON  = e_c*4*nes{3}(rngECH_ON(b1)).*Tes{3}(rngECH_ON(b1)).*C_s(Tes{3}(rngECH_ON(b1)),2);
HeatFlux_ECH_OFF = e_c*4*nes{2}(rngECH_ON(b1)).*Tes{2}(rngECH_ON(b1)).*C_s(Tes{2}(rngECH_ON(b1)),2);

plot(R(rngECH_ON(b1))-rOffset, (HeatFlux_ECH_ON- HeatFlux_ECH_OFF)*1e-6,'ko-')
xlim([-4,4]); ylabel('$[MWm^{-2}]$','Interpreter','Latex')


return
% Plot ne and Te on-axis
n1 = find(shotlist == 20822); % ECH ON
n0 = find(shotlist == 20808); % ECH OFF

figure
subplot(2,1,1)
title('DLP 10.5 on-axis')
hold on
k(2) = plot(t_Ni_gf{n1},Ni_gf{n1},'r','lineWidth',2);
k(1) = plot(t_Ni_gf{n0},Ni_gf{n0},'k','lineWidth',2);
rngRF = find(t_rf{n0}>=4.14 & t_rf{n0}<=5.2);
plot(t_rf{n0}(rngRF),(RF{n0}(rngRF).^2)*2e19,'k','lineWidth',0.5)
plot(t_rf{n1}(rngRF),(RF{n1}(rngRF).^2)*2e19,'r','lineWidth',0.5)
plot(t_ech{n1}(1:length(ECH{n1})),ECH{n1}*1e18,'r','lineWidth',0.5)
ylim([0,10e18])
xlim([3.9,5.4])
box on

subplot(2,1,2)
title('DLP 10.5 on-axis')
hold on
k(2) = plot(t_Te_gf{n1},Te_gf{n1},'r','lineWidth',2);
k(1) = plot(t_Te_gf{n0},Te_gf{n0},'k','lineWidth',2);
legend([k(1),k(2)],'28 GHz OFF','28 GHz ON','location','northeast')
plot(t_rf{n0}(rngRF),(RF{n0}(rngRF).^2)*10,'k','lineWidth',0.5)
plot(t_rf{n1}(rngRF),(RF{n1}(rngRF).^2)*10,'r','lineWidth',0.5)
plot(t_ech{n1}(1:length(ECH{n1})),ECH{n1}*2,'r','lineWidth',0.5)
ylim([0,10])
xlim([3.9,5.4])
set(gcf,'color','w','position',[360.3333  184.3333  399.3333  433.3333])
box on

%% Convert data into Excel
% Save data to Excel spreadsheet:
% Create table to export to EXCEL
if 0
G = [shotlist(b)',PS1(b)',ne{p}(b)',dne{p}(b)',Te{p}(b)',dTe{p}(b)'];
F = {'Shot','R [cm]','ne[m^-3]','dne','Te[eV]','dTe'};
FileName = 'NeTe_Spool_10_5_2018_01_16.xlsx';
xlswrite(FileName,[F;num2cell(G)]);
end