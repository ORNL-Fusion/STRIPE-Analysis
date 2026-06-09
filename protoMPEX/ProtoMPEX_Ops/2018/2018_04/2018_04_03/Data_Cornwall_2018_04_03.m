clear all
close all

ProbeLoc = 'A'; % DLP 10.5
FitShow = 0; 
RawDataShow = 1;

switch ProbeLoc
    case 'A'
shotlist = [20000 + [808,809,810,812,813,814 ,815,816,817,818,819,820,822,823,824,827,828,829,830,831,832,834,835]];
R        =          [6.0,5.5,5.0,4.0,4.5,4.75,6.5,7.0,7.5,8.0,8.5,9.0,6.5,6.0,5.5,4.0,6.5,7.0,7.5,8.0,8.5,5.0,4.5];
ECHpwr   =          [0  ,0  ,0  ,0  ,0  ,0   ,0  ,0  ,0  ,0  ,0  ,0  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ]; 

SweepType = 'iso'; % "iso" for isolated sweep, "niso" for non-isolated sweep
AttType = 'Vx1,Ix1'; % Attenuation on the digitized signals
DLPType = '10MP';
Config.FitFunction = 2; 
ChannelType = '5' ;
Config.tStart = 4.15; % [s]
Config.tEnd = 5.5;
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
% #########################################################################
DA{1} = [RootAddress,'RF_FWD_PWR'];
[RF,t_rf]   = my_mdsvalue_v2(shotlist,DA(1))
% #########################################################################

A_fp = 0.25*sqrt(2)*pi*(1.22/1000)^2;

%%
% #########################################################################
% PLOT DATA:
% #########################################################################

close all
GasCompareALL_2018_04_03
C = {'k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:'};
t = t_zero(shotlist);

TimePlotStart = 4.1;
TimePlotEnd = 5.5;


figure; 
subplot(2,2,1); hold on
for s = 1:length(shotlist)
    GoodFits{s} = GlitchFlag{s} == 0 & StdResNorm{s}<=0.14 & Ni{s}>0 & Ni{s}<0.8e19 & Te{s}<=8 & Is{s}>0.3e-3;
            
    t_Ni_gf{s} = time{s}(GoodFits{s});
    Ni_gf{s} = sgolay_t(Ni{s}(GoodFits{s}),3,11);
    
    h(s) = plot(t_Ni_gf{s},Ni_gf{s},C{s},'lineWidth',2);

    hold on
    plot(t_ech{s}(1:length(ECH{s})),ECH{s}*1e18,C{s},'lineWidth',0.5)
    
    if s == 2
        continue
    end
    plot(t_rf{s}(1:length(RF{s})),(RF{s}.^2)*0.25e19,C{s},'lineWidth',0.5)

end
title('$ n_e $ $ [m^{-3}] $','interpreter','Latex','FontSize',13,'Rotation',0)
legend(h,['DLP ',num2str(DLP)],'location','SouthEast')
ylim([0,1.5e19])
xlim([TimePlotStart,TimePlotEnd])

subplot(2,2,2); hold on
for s = 1:length(shotlist)
        t_Te_gf{s} = time{s}(GoodFits{s});
        Te_gf{s} = sgolay_t(Te{s}(GoodFits{s}),3,11);
        
        h(s) = plot(t_Te_gf{s},Te_gf{s},C{s},'lineWidth',2);
        
        L{s} = [num2str(shotlist(s)),' ,t=',num2str(t{s}(10:14))];
        plot(t_ech{s}(1:length(ECH{s})),2*ECH{s},C{s},'lineWidth',0.5)

end
legend(h,L,'location','NorthEast')
title('$ T_e $ $ [eV] $','interpreter','Latex','FontSize',13,'Rotation',0)
ylim([0,20])
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
close all

for s = 1:length(shotlist)
    rngT{1} = [t_Ni_gf{s}>=4.20 & t_Ni_gf{s}<=4.23]; 
    rngT{2} = [t_Ni_gf{s}>=4.7 & t_Ni_gf{s}<=4.8];  
    rngT{3} = [t_Ni_gf{s}>=4.46 & t_Ni_gf{s}<=4.51]; 
    for p = 1:3
        nes{p}(s)    = mean(Ni_gf{s}(rngT{p}));
        dnes{p}(s)   = std(Ni_gf{s}(rngT{p}),1,2);
        Tes{p}(s)    = mean(Te_gf{s}(rngT{p}));
        dTes{p}(s)   = std(Te_gf{s}(rngT{p}),1,2);
    end
end

rOffset = 6.5;
rngECH_ON  = find(ECHpwr == 1);
rngECH_OFF = find(ECHpwr == 0);
[a1,b1] = sort(R(rngECH_ON ));
[a0,b0] = sort(R(rngECH_OFF));

figure; 
subplot(2,1,1); hold on
for p = 3
k(1) = errorbar(R(rngECH_ON(b1))-rOffset,nes{p}(rngECH_ON(b1)),dnes{p}(rngECH_ON(b1)),'LineWidth',1,'color','r','marker','o')
k(2) = errorbar(R(rngECH_OFF(b0))-rOffset,nes{p}(rngECH_OFF(b0)),dnes{p}(rngECH_OFF(b0)),'LineWidth',1,'color','k')
end
legend([k(1),k(2)],'28 GHz ON','28 GHz OFF')


ylim([0,0.7]*1e19); ylabel('$n_e$ $[m^{-3}]$','Interpreter','Latex')
xlim([-4,4]); xlabel('$R$ $[cm]$','Interpreter','Latex')
box on
grid on
t = title('$n_e$ $[m^{-3}]$ '); t.Interpreter = 'latex'; t.FontSize = 13;

subplot(2,1,2); hold on
for p = 3
k(1) = errorbar(R(rngECH_ON(b1))-rOffset,Tes{p}(rngECH_ON(b1)),dTes{p}(rngECH_ON(b1)),'LineWidth',1,'color','r','marker','o')
k(2) = errorbar(R(rngECH_OFF(b0))-rOffset,Tes{p}(rngECH_OFF(b0)),dTes{p}(rngECH_OFF(b0)),'LineWidth',1,'color','k')
end
legend([k(1),k(2)],'28 GHz ON','28 GHz OFF','location','south')
ylim([0,7]); xlabel('$R$ $[cm]$','Interpreter','Latex')
xlim([-4,4]); ylabel('$T_e$ $[eV]$','Interpreter','Latex')
box on
grid on
t = title('$T_e$ $[eV]$ '); t.Interpreter = 'latex'; t.FontSize = 13;

set(gcf,'color','w','position',[360.3333  184.3333  399.3333  433.3333])

% Plot ne and Te on-axis
n1 = find(shotlist == 20822); % ECH ON
n0 = find(shotlist == 20808); % ECH OFF
rng_t = find(t_Ni_gf{n1}<5.135);

figure
subplot(2,1,1)
title('DLP 10.5 on-axis')
hold on
k(2) = plot(t_Ni_gf{n1}(rng_t),Ni_gf{n1}(rng_t),'r','lineWidth',2);
k(1) = plot(t_Ni_gf{n0}(rng_t),Ni_gf{n0}(rng_t),'k','lineWidth',2);
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
k(2) = plot(t_Te_gf{n1}(rng_t),Te_gf{n1}(rng_t),'r','lineWidth',2);
k(1) = plot(t_Te_gf{n0}(rng_t),Te_gf{n0}(rng_t),'k','lineWidth',2);
legend([k(1),k(2)],'28 GHz OFF','28 GHz ON','location','northeast')
plot(t_rf{n0}(rngRF),(RF{n0}(rngRF).^2)*10,'k','lineWidth',0.5)
plot(t_rf{n1}(rngRF),(RF{n1}(rngRF).^2)*10,'r','lineWidth',0.5)
plot(t_ech{n1}(1:length(ECH{n1})),ECH{n1}*2,'r','lineWidth',0.5)
ylim([0,10])
xlim([3.9,5.4])
set(gcf,'color','w','position',[360.3333  184.3333  399.3333  433.3333])
box on

figure

HeatFlux_ECH_ON  = e_c*4*nes{3}(rngECH_ON(b1)).*Tes{3}(rngECH_ON(b1)).*C_s(Tes{3}(rngECH_ON(b1)),2);
HeatFlux_ECH_OFF = e_c*4*nes{2}(rngECH_ON(b1)).*Tes{2}(rngECH_ON(b1)).*C_s(Tes{2}(rngECH_ON(b1)),2);

plot(R(rngECH_ON(b1))-rOffset, (HeatFlux_ECH_ON- HeatFlux_ECH_OFF)*1e-6,'ko-')
xlim([-4,4]); ylabel('$[MWm^{-2}]$','Interpreter','Latex')
set(gcf,'position',[463.6667  395.6667  456.6667  222.0000],'color','w')
box on
xlabel('$R$ $[cm]$','Interpreter','Latex')

%% 
if 0
    % ECH OFF radial scan
DataToWrite.FileName = 'DLP_radialScan_ECH_OFF.xlsx';
DataToWrite.Column{1}.Heading = 'Shot';
DataToWrite.Column{1}.Data = shotlist(rngECH_OFF(b0));
DataToWrite.Column{2}.Heading = 'R [cm]';
DataToWrite.Column{2}.Data = R(rngECH_OFF(b0))-Offset;
DataToWrite.Column{3}.Heading = 'Ne [m^-3]';
DataToWrite.Column{3}.Data = nes{3}(rngECH_OFF(b0));
DataToWrite.Column{4}.Heading = 'dNe [m^-3]';
DataToWrite.Column{4}.Data = dnes{3}(rngECH_OFF(b0));
DataToWrite.Column{5}.Heading = 'Te [eV]';
DataToWrite.Column{5}.Data = Tes{3}(rngECH_OFF(b0));
DataToWrite.Column{6}.Heading = 'dTe [eV]';
DataToWrite.Column{6}.Data = dTes{3}(rngECH_OFF(b0));
DataToExcel(DataToWrite)

    % ECH ON radial scan
clear DataToWrite
DataToWrite.FileName = 'DLP_radialScan_ECH_ON.xlsx';
DataToWrite.Column{1}.Heading = 'Shot';
DataToWrite.Column{1}.Data = shotlist(rngECH_ON(b1));
DataToWrite.Column{2}.Heading = 'R [cm]';
DataToWrite.Column{2}.Data = R(rngECH_ON(b1))-Offset;
DataToWrite.Column{3}.Heading = 'Ne [m^-3]';
DataToWrite.Column{3}.Data = nes{3}(rngECH_ON(b1));
DataToWrite.Column{4}.Heading = 'dNe [m^-3]';
DataToWrite.Column{4}.Data = dnes{3}(rngECH_ON(b1));
DataToWrite.Column{5}.Heading = 'Te [eV]';
DataToWrite.Column{5}.Data = Tes{3}(rngECH_ON(b1));
DataToWrite.Column{6}.Heading = 'dTe [eV]';
DataToWrite.Column{6}.Data = dTes{3}(rngECH_ON(b1));
DataToExcel(DataToWrite)

    % Time traces, DLP with and without ECH
clear DataToWrite
DataToWrite.FileName = 'DLP_TimeTrace_ECH_OFF_ON_Shot_20808_20822.xlsx';
DataToWrite.Column{1}.Heading = 't [s]';
DataToWrite.Column{1}.Data = t_Ni_gf{n0}(rng_t);
DataToWrite.Column{2}.Heading = 'n_e [m^-3], ECH OFF';
DataToWrite.Column{2}.Data = Ni_gf{n0}(rng_t);
DataToWrite.Column{3}.Heading = 'T_e [eV], ECH OFF';
DataToWrite.Column{3}.Data = Te_gf{n0}(rng_t);
DataToWrite.Column{4}.Heading = 'n_e [m^-3], ECH ON';
DataToWrite.Column{4}.Data = Ni_gf{n1}(rng_t);
DataToWrite.Column{5}.Heading = 'T_e [eV], ECH ON';
DataToWrite.Column{5}.Data = Te_gf{n1}(rng_t);
DataToExcel(DataToWrite)

    % Time traces, ECH and RF power
clear DataToWrite
rng_t_rf = find(t_rf{n0}>4.149 & t_rf{n0}<5.23);
DataToWrite.FileName = 'Pwr_TimeTrace_ECH_OFF_ON_Shot_20808_20822.xlsx';
DataToWrite.Column{1}.Heading = 't [s]';
DataToWrite.Column{1}.Data = t_rf{n0}(rng_t_rf);
DataToWrite.Column{2}.Heading = 'RF power [A.U], ECH OFF';
DataToWrite.Column{2}.Data = RF{n0}(rng_t_rf).^2;
DataToWrite.Column{3}.Heading = 't [s]';
DataToWrite.Column{3}.Data = t_rf{n1}(rng_t_rf);
DataToWrite.Column{4}.Heading = 'RF power [A.U], ECH ON';
DataToWrite.Column{4}.Data = RF{n1}(rng_t_rf).^2;
DataToWrite.Column{5}.Heading = 't [s]';
DataToWrite.Column{5}.Data = t_ech{n0}(rng_t_rf);
DataToWrite.Column{6}.Heading = 'ECH power [A.U], ECH OFF';
DataToWrite.Column{6}.Data = ECH{n0}(rng_t_rf);
DataToWrite.Column{7}.Heading = 't [s]';
DataToWrite.Column{7}.Data = t_ech{n1}(rng_t_rf);
DataToWrite.Column{8}.Heading = 'ECH power [A.U], ECH ON';
DataToWrite.Column{8}.Data = ECH{n1}(rng_t_rf);
DataToWrite.Column{9}.Heading = 'Filtered ECH power [A.U], ECH ON';
DataToWrite.Column{9}.Data = sgolay_t(ECH{n1}(rng_t_rf),3,57);
DataToExcel(DataToWrite)
end