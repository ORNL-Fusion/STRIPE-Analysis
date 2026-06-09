clear all
close all

ProbeLoc = 'A'; % DLP 10.5
% ProbeLoc = 'B'; % DLP 1.5
FitShow = 0; 
RawDataShow = 1;

switch ProbeLoc
    case 'A'
% #########################################################################
% =========================================================================
% 2018_04_05:
% First few shots:
% =========================================================================
shotlist = [20000 + [884,885]];
shotlist = [20000 + [889]];

% % Compare the 2.5 and 2.0 kA cases
shotlist = [20000 + [891,892,893,896]];

% Test if TS measurements of the plasma can be done
% 100 kW, TR2 = 100 A, PS1/2 = 2kA
shotlist = [20000 + [897]];
% 70 kW, TR2 = 120 A, PS1/2 = 4.5kA
shotlist = [20000 + [899]];
% 100 kW, TR2 = 160 A, PS1/2 = 4.5kA
shotlist = [20000 + [900]];
% 100 kW, TR2 = 160 A, PS1/2 = 2.5kA
shotlist = [20000 + [901]];
% 80 kW, TR2 = 160 A, PS1/2 = 2.5kA
shotlist = [20000 + [902]];
% Results show that at the high field, high RF power the TS system was able
% to see scattered photons but at the required low densities it was not
% =========================================================================

% 40 kW, TR2 = 100 A, PS1/2 = 2.5kA
shotlist = [20000 + [903]];
% 40 kW, TR2 = 100 A, compare 2.5 and 2.0 kA cases
shotlist = [20000 + [903,904]];

% =========================================================================
% 1st radial scan at DLP 10.5 
% 40 kW, TR2 = 100 A, PS1/2 = 2.0kA
shotlist = [20000 + [904,905,906,907,909,910,911,914,915,916,917,918]];
R        =          [6.0,5.5,5.0,4.5,4.0,3.5,6.5,7.0,7.5,8.0,8.5,9.0] ;
ECHpwr   =          [0  ,0  ,0  ,0  ,0  ,0  ,0  ,1  ,1  ,1  ,1  ,1  ] ;

% DLP 10.5 on-axis with ECH
shotlist = [20000 + [904,919,920]];

shotlist = [20000 + [914,915,916,917,918,919]];
R        =          [7.0,7.5,8.0,8.5,9.0,6.5] ;
ECHpwr   =          [1  ,1  ,1  ,1  ,1  ,1  ] ;
% It was found that plasma density was oscillating in time which lead to
% very noisy DLP measurements, we needed to find better base plasma
% conditions.
% To solve this, we need to increase PS2 from 2 kA to 2.5 kW (coils 1 to 6)
% while keeping PS1 = 2.0 kA to remove the effect of the 2nd harmonic
% resonance
% =========================================================================
% Test and optimize the PS1 = 2.0 kA, PS2 = 2.5 kA case
shotlist = [20000 + [923,924,925,926,927]];
% First direct comparison of the effect of the 2nd harmonic
shotlist = [20000 + [927,928]];

% shotlist = [20000 + [927,928,929,930,931,932]];
% shotlist = [20000 + [931,933,934,935]];
% 
% Compare response between Target and dump sides during 28 GHz injection
shotlist = [20000 + [935,936]];
% 
% Begin radial scan at dlp 10.5
% use GoodFits{s} = GlitchFlag{s} == 0 & StdResNorm{s}<=0.13 & Ni{s}>0 & Ni{s}<6e18 & Te{s}<=8 & Is{s}>0.1e-3;
shotlist = [20000 + [935,939,940,941,942,943,944,945,946,948,952,957,955,956]];
R        =          [6.0,5.5,5.0,4.5,4.0,3.5,3.0,6.5,7.0,7.5,8.0,8.5,9.0,9.5] ;
ECHpwr   =          [1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ] ;

% % 
% % 28 GHz power scan
% shotlist = [20000 + [961,962,963,964,965,966,969]];
% shotlist = [20000 + [966,969,971,972]];

% HELIOS relevant shots:
% shotlist = [20000 + [961]];

SweepType = 'iso'; % "iso" for isolated sweep, "niso" for non-isolated sweep
AttType = 'Vx1,Ix1'; % Attenuation on the digitized signals
DLPType = '10MP';
Config.FitFunction = 2; 
ChannelType = '5' ;
Config.tStart = 4.15; % [s]
Config.tEnd = 5.2;
Config.Center_V = 0; % Remove offset on V: 1 (yes) 0(no)
Config.Center_I = 0; % Remove offset on I: 1 (yes) 0(no)
    
case 'B'
% =========================================================================
% =========================================================================
        
Config.FitFunction = 2; 

DLPType = '1MP'; % DLP 1.5 horizontal
ChannelType = '5' ; % "1" for TARGET_LP,LP_V_RAMP, "2" for other options

Config.tStart = 4.25; % [s]
Config.tEnd = 4.7;
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
%         Config.L_tip = 1.6/1000;% according to Nischal 2017_11_21
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
    case 'IF'
        DLP = 11.5;
        Config.L_tip = 0/1000;
        Config.D_tip = 1/1000; % [m]
    case '10FluxProbe';
        DLP = 10.5;
        Config.L_tip = 0/1000;
        Config.D_tip = ((2^(1/4))*1.22)/1000; % [m]
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
[ni,Te,Isat,time,Ifit,Ip,Vp,tm,Vsweep,Isweep,GlitchFlag,SSQres,StdRes,StdResNorm] = DLP_fit_V5_6(Config,shotlist,DataAddress);

%Ni,Te,time,Ifit,Ip,Vp,tm,Vsweep,Isweep,GlitchFlag,SSQres
for s = 1:length(shotlist)
    Ni{s} = 0.5*(ni{s}{1} + ni{s}{2});
    Is{s} = 0.5*(Isat{s}{1} + Isat{s}{2});
end
% #########################################################################
DA{1} = [RootAddress,'PWR_28GHZ']; % Isx
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
GasCompareALL_2018_04_05
C = {'k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:'};
t = t_zero(shotlist);

TimePlotStart = 4.1;
TimePlotEnd = 6;

figure; 
subplot(2,2,1); hold on
for s = 1:length(shotlist)
    GoodFits{s} = GlitchFlag{s} == 0 & StdResNorm{s}<=0.13 & Ni{s}>0 & Ni{s}<6e18 & Te{s}<=8 & Is{s}>0.1e-3;
      
%     plot(time{s}(GoodFits{s}),ni{s}{1}(GoodFits{s}),C{s},'lineWidth',1);
%     plot(time{s}(GoodFits{s}),ni{s}{2}(GoodFits{s}),C{s},'lineWidth',1);
      
    t_Ni_gf{s} = time{s}(GoodFits{s});
    Ni_gf{s} = Ni{s}(GoodFits{s});
    
    h(s) = plot(t_Ni_gf{s},Ni_gf{s},C{s},'lineWidth',2);

    hold on
    plot(t_ech{s}(1:length(ECH{s})),ECH{s}*0.1e19,C{s},'lineWidth',0.5)
    
    if s == 2
        continue
    end
    plot(t_rf{s}(1:length(RF{s})),(RF{s}.^2)*1e19,C{s},'lineWidth',0.5)

end
title('$ n_e $ $ [m^{-3}] $','interpreter','Latex','FontSize',13,'Rotation',0)
legend(h,['DLP ',num2str(DLP)],'location','SouthEast')
%set(gca,'PlotBoxAspectRatio',[1 1 1])
ylim([0,7e18])
% ylim([0,10e19])
xlim([TimePlotStart,TimePlotEnd])

subplot(2,2,2); hold on
for s = 1:length(shotlist)
        t_Te_gf{s} = time{s}(GoodFits{s});
        Te_gf{s} = sgolay_t(Te{s}(GoodFits{s}),3,7);
        
        h(s) = plot(t_Te_gf{s},Te_gf{s},C{s},'lineWidth',2);
        
        L{s} = [num2str(shotlist(s)),' ,t=',num2str(t{s}(10:14))];
        plot(t_ech{s}(1:length(ECH{s})),2*ECH{s},C{s},'lineWidth',0.5)

end
legend(h,L,'location','NorthEast')
title('$ T_e $ $ [eV] $','interpreter','Latex','FontSize',13,'Rotation',0)
ylim([0,10])
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
    for k = 1:16
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

% return
%%
% #########################################################################
% Steady state values
% #########################################################################
% close all

for s = 1:length(shotlist)
    rngT{1} = [time{s}>=4.20 & time{s}<=4.23]; 
    rngT{3} = [time{s}>=4.44 & time{s}<=4.48];  % ECH ONF
    rngT{2} = [time{s}>=4.69 & time{s}<=4.75];  % ECH OFF
for p = 1:3
        nes{p}(s)    = mean(Ni{s}(rngT{p} & GoodFits{s}));
        dnes{p}(s)   = std(Ni{s}(rngT{p} & GoodFits{s}),1,2);
        Tes{p}(s)    = mean(Te{s}(rngT{p} & GoodFits{s}));
        dTes{p}(s)   = std(Te{s}(rngT{p} & GoodFits{s}),1,2);
    end
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
end
set(y(1),'color','k')

ylim([0,0.7]*1e19); ylabel('$n_e$ $[m^{-3}]$','Interpreter','Latex')
xlim([-4,4]); xlabel('$R$ $[cm]$','Interpreter','Latex')
box on
grid on
t = title('$n_e$ $[m^{-3}]$ '); t.Interpreter = 'latex'; t.FontSize = 13;

subplot(2,1,2); hold on
for p = 2:3
k(1) = errorbar(R(rngECH_ON(b1))-rOffset,Tes{p}(rngECH_ON(b1)),dTes{p}(rngECH_ON(b1)),'LineWidth',1,'color','r','marker','o')
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

HeatFlux_ECH_OFF  = e_c*4*nes{3}(rngECH_ON(b1)).*Tes{3}(rngECH_ON(b1)).*C_s(Tes{3}(rngECH_ON(b1)),2);
HeatFlux_ECH_ON = e_c*4*nes{2}(rngECH_ON(b1)).*Tes{2}(rngECH_ON(b1)).*C_s(Tes{2}(rngECH_ON(b1)),2);

plot(R(rngECH_ON(b1))-rOffset, (HeatFlux_ECH_ON- HeatFlux_ECH_OFF)*1e-6,'ko-')
xlim([-4,4]); ylabel('$[MWm^{-2}]$','Interpreter','Latex')
set(gcf,'position',[463.6667  395.6667  456.6667  222.0000],'color','w')
box on
ylim([0,0.25])
xlabel('$R$ $[cm]$','Interpreter','Latex')
% Plot ne and Te on-axis
n1 = find(shotlist == 20935); % ECH ON

figure
subplot(2,1,1)
title('DLP 10.5 on-axis')
hold on
k(2) = plot(t_Ni_gf{n1},Ni_gf{n1},'r','lineWidth',2);
rngRF = find(t_rf{n1}>=4.14 & t_rf{n1}<=5.2);
plot(t_ech{n1}(1:length(ECH{n1})),ECH{n1}*1e18,'r','lineWidth',0.5)
ylim([0,10e18])
xlim([3.9,5.4])
box on
title('Shot: 20935')

subplot(2,1,2)
title('DLP 10.5 on-axis')
hold on
k(2) = plot(t_Te_gf{n1},Te_gf{n1},'r','lineWidth',2);
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