clear all
close all

ProbeLoc = 'B';
ProbeLoc = 'A';
FitShow = 0; 
RawDataShow = 1;

XP = 4
% 1 - Finding plasma core
% 2 - Reducing first puff from 4.5 to 3.5 v
% 3 - Radial scan at 100 kW, TR2 =  180 A, PS1/2 = 4500 A
% 4 - RF power scan
% 5 - Radial scan at 70 kW , TR2 =  180 A, PS1/2 = 4500 A
switch ProbeLoc
    case 'A'
        switch XP
            case 1
% #########################################################################
% 2018_04_13:
% finding plasma core:
shotlist = [21000 + [186,188,189,190,192]];
shotlist = [21000 + [192,195,196,197,198]];
% x = 5.5 cm
shotlist = [21000 + [198 ,199,200 ,202 ,203 ,204 ,205 ]];
side     =          ['R1','0','L1','R1','R2','R3','R1'] ;
x        =          [5.5 ,5.5,5.5 ,5.5 ,5.5 ,5.5 ,5.5 ] ;

shotlist = [21000 + [205 ,206 ,207 ,208 ,209 ,210 ,212 ,213 ]];
side     =          ['R1','R1','0' ,'L1','R2','R3','R1','R1'] ;
x        =          [5.5 ,5.75,5.75,5.75,5.75,5.75,5.75,5.75] ;

shotlist = [21000 + [213 ,214 ,215 ,216 ,217 ,218 ,219   ,220 ]];
side     =          ['R1','R1','R1','R1','R1','R1','R0.5','R0'] ;
x        =          [5.75,6.0 ,6.0 ,6.25,6.5 ,6.25,6.25  ,6.25] ;

shotlist = [21000 + [219   ,220 ,222   ]];
side     =          ['R0.5','R0','R0.5'] ;
x        =          [6.25  ,6.25,6.0   ] ;
% #########################################################################
            case 2
% #########################################################################
% Reduce 1st puff from 4.5 v to 3.5 v
shotlist = [21000 + [222   ,223   ,224   ]];
side     =          ['R0.5','R0.5','R0.5'] ;
x        =          [6.0   ,6.0   ,6.0   ] ;
% #########################################################################
            case 3
% #########################################################################
% Radial scan at 100 kW, 3.5 v gas flow
shotlist = [21000 + [224,225,227,228,230,231,232,233,234,235,236,237]];
x        =          [6.0,5.5,5.0,4.5,4.0,3.5,6.5,7.0,7.5,8.0,8.5,6.0] ;
shotlist = [21000 + [225,227,228,230,231,232,233,234,235,236,237]];
x        =          [5.5,5.0,4.5,4.0,3.5,6.5,7.0,7.5,8.0,8.5,6.0] ;
% #########################################################################
            case 4
% #########################################################################
% RF power scan
shotlist = [21000 + [238,239,241,245,246,247]];
RFpwr =             [10 ,9  ,8  ,7  ,6  ,5  ] ;

% shotlist = [21000 + [238,245,247]];
% RFpwr =             [10 ,7  ,5  ] ;
% #########################################################################
            case 5
% #########################################################################
% Radial scan at 70 kW
shotlist = [21000 + [248,249,251,252,253,254]];
x        =          [6.0,6.5,7.0,7.5,8.0,8.5] ;

shotlist = [21000 + [248,249,255,256,257,258,259]];
x        =          [6.0,6.5,5.5,5.0,4.5,4.0,6.0] ;

shotlist = [21000 + [248,249,251,252,253,254,255,256,257,258,259]];
x        =          [6.0,6.5,7.0,7.5,8.0,8.5,5.5,5.0,4.5,4.0,6.0] ;

shotlist = [21000 + [248,249,251,252,253,254,255,256,257,258]];
x        =          [6.0,6.5,7.0,7.5,8.0,8.5,5.5,5.0,4.5,4.0] ;

% #########################################################################
            case 6
% #########################################################################
% TR2 scan, on -axis flux probe, 100 kW
shotlist = [21000 + [260,261,262,263,264]];
TR2 =               [180,160,140,120,100] ;
RFpwr =               [180,160,140,120,100] ;

% #########################################################################
            case 7
% #########################################################################
% PS1 scan, on -axis flux probe, 100 kW
shotlist = [21000 + [260,265,266,267,268,269,270]];
PS1 =               [4.5,4.0,3.5,3.0,2.5,2.0,2.0] ;
% #########################################################################
            case 8
% #########################################################################
% Radial scan, PS1 = 2kA 100 kW
% shotlist = [21000 + [269,270,271,272,273,274,275]];
% x =                 [6.0,6.5,7.0,8.0,9.0,5.5,5.0] ;
% #########################################################################
        end
        
SweepType = 'iso'; % "iso" for isolated sweep, "niso" for non-isolated sweep
AttType = 'Vx1,Ix1'; % Attenuation on the digitized signals
DLPType = '10FluxProbe';
Config.FitFunction = 1; 

if 0
% #########################################################################
% =========================================================================
shotlist = [20000 + [572,573,574]];
DLPType = '1MP';
% =========================================================================
% #########################################################################
end

Config.FitFunction = 2; 
ChannelType = '5' ;

Config.tStart = 4.15; % [s]
Config.tEnd = 4.7;
Config.Center_V = 0; % Remove offset on V: 1 (yes) 0(no)
Config.Center_I = 0; % Remove offset on I: 1 (yes) 0(no)
    
case 'B'
% =========================================================================
% =========================================================================
        
Config.FitFunction = 2; 

DLPType = '1MP'; % DLP 1.5 horizontal
ChannelType = '2' ; % "1" for TARGET_LP,LP_V_RAMP, "2" for other options

Config.tStart = 4.16; % [s]
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
Config.FilterDataInput = 0; % Filter input data with savitsky Golay filter order 3 frame size 7
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
GasCompareALL_2018_04_13
C = {'k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:'};
t = t_zero(shotlist);

TimePlotStart = 4.1;
TimePlotEnd = 4.7;

figure; 
subplot(2,2,1); hold on
for s = 1:length(shotlist)
    GoodFits{s} = GlitchFlag{s} == 0 & StdResNorm{s}<=0.2 & Ni{s}>0 & Ni{s}<1e21 & Te{s}<=30;
      
%     plot(time{s}(GoodFits{s}),ni{s}{1}(GoodFits{s}),C{s},'lineWidth',1);
%     plot(time{s}(GoodFits{s}),ni{s}{2}(GoodFits{s}),C{s},'lineWidth',1);
      
    t_Ni_gf{s} = time{s}(GoodFits{s});
    Ni_gf{s} = Ni{s}(GoodFits{s});
    
    h(s) = plot(t_Ni_gf{s},Ni_gf{s},C{s},'lineWidth',2);

    hold on
    plot(t_ech{s}(1:length(ECH{s})),ECH{s}*3e19,C{s},'lineWidth',0.5)
    
    if s == 2
        continue
    end
    plot(t_rf{s}(1:length(RF{s})),(RF{s}.^2)*15e19,C{s},'lineWidth',0.5)

end
title('$ n_e $ $ [m^{-3}] $','interpreter','Latex','FontSize',13,'Rotation',0)
legend(h,['DLP ',num2str(DLP)],'location','SouthEast')
%set(gca,'PlotBoxAspectRatio',[1 1 1])
ylim([0,14e19])
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
%     plot(tm{s},Vp{s})
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
    for k = 1:12
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
%         ht = title(['T_e: ',num2str(Te{s}(c),2),' ,t = ',num2str(time{s}(c),5),' s']);
        ht = title(['I_s: ',num2str(Is{s}(c),2),' ,t = ',num2str(time{s}(c),5),' s']);
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
%% Steady state values

A_fp = 0.25*sqrt(2)*pi*(1.2/1000)^2;
Ac = pi*(2.0/100)^2;

for s = 1:length(shotlist)
    rngT{1} = [time{s}>=4.20 & time{s}<=4.23]; 
    rngT{2} = [time{s}>=4.25 & time{s}<=4.35];  
    rngT{3} = [time{s}>=4.6 & time{s}<=4.62]; 
    for p = 1:3
        Ifp{p}(s)    = mean(Is{s}(rngT{p} & GoodFits{s}));
        dIfp{p}(s)   = std(Is{s}(rngT{p} & GoodFits{s}),1,2);
        flux{p}(s)  = Ifp{p}(s)/(e_c*A_fp);
        dflux{p}(s) = dIfp{p}(s)/(e_c*A_fp);
        Tefp{p}(s)    = mean(Te{s}(rngT{p} & GoodFits{s}));
        dTefp{p}(s)   = std(Te{s}(rngT{p} & GoodFits{s}),1,2);
        
        
        pg1_ss{p}(s)     = mean(interp1(Time,P{1}{s} - Offset{1}{s},time{s}(rngT{p})));
        dpg1_ss{p}(s)    = std(interp1(Time,P{1}{s} - Offset{1}{s},time{s}(rngT{p})),1);
        pg2_ss{p}(s)     = mean(interp1(Time,P{2}{s} - Offset{2}{s},time{s}(rngT{p})));
        pg3_ss{p}(s)     = mean(interp1(Time,P{3}{s} - Offset{3}{s},time{s}(rngT{p})));
        pg4_ss{p}(s)     = mean(interp1(Time,P{4}{s} - Offset{4}{s},time{s}(rngT{p})));

    end
end

%% Power scan
if XP == 4
figure
[a,b] = sort(RFpwr);
subplot(2,2,1); hold on
for p = 3
h(p) = errorbar(RFpwr(b),flux{p}(b),dflux{p}(b));
end
xlim([4,11])
t = title('${\Gamma}_+$ $[m^{-2}s^{-1}]$ '); t.Interpreter = 'latex'; t.FontSize = 13;
box on

subplot(2,2,2); hold on
for p = 3
h(p) = errorbar(RFpwr(b),pg1_ss{p}(b),dpg1_ss{p}(b),'ko-')
end
xlim([4,11])
ylim([0,0.3])
t = title('$P_{Target}$ $[Pa]$ '); t.Interpreter = 'latex'; t.FontSize = 13;
box on

subplot(2,2,3); hold on
for p = 3
h2 = plot(RFpwr(b),pg2_ss{p}(b),'ko-')
h4 = plot(RFpwr(b),pg4_ss{p}(b),'ro-')
end
legend([h2,h4],'PG 2.5','PG 4.5')
xlim([4,11])
ylim([0,20])
t = title('$P_{Source}$ $[Pa]$ '); t.Interpreter = 'latex'; t.FontSize = 13;
box on

subplot(2,2,4); hold on
for p = 3
h(p) = plot(RFpwr(b),pg3_ss{p}(b),'ko-')
end
xlim([4,11])
ylim([-0.5,2.5])
t = title('$P_{CC}$ $[Pa]$ '); t.Interpreter = 'latex'; t.FontSize = 13;
box on

set(gcf,'color','w')
end
%% TR2 scan
if XP == 6
figure
[a,b] = sort(RFpwr);
subplot(2,2,1); hold on
for p = 3
h(p) = errorbar(RFpwr(b),flux{p}(b),dflux{p}(b));
end
xlim([80,200])
ylim([0,15]*1e23)
t = title('${\Gamma}_+$ $[m^{-2}s^{-1}]$ '); t.Interpreter = 'latex'; t.FontSize = 13;
box on

subplot(2,2,2); hold on
for p = 3
h(p) = plot(RFpwr(b),pg1_ss{p}(b),'ko-')
end
xlim([80,200])
ylim([0,2.5])
t = title('$P_{Target}$ $[Pa]$ '); t.Interpreter = 'latex'; t.FontSize = 13;
box on

subplot(2,2,3); hold on
for p = 3
h2 = plot(RFpwr(b),pg2_ss{p}(b),'ko-')
h4 = plot(RFpwr(b),pg4_ss{p}(b),'ro-')
end
legend([h2,h4],'PG 2.5','PG 4.5')
xlim([80,200])
ylim([0,10])
t = title('$P_{Source}$ $[Pa]$ '); t.Interpreter = 'latex'; t.FontSize = 13;
box on

subplot(2,2,4); hold on
for p = 3
h(p) = plot(RFpwr(b),pg3_ss{p}(b),'ko-')
end
xlim([80,200])
ylim([-0.5,2.5])
t = title('$P_{CC}$ $[Pa]$ '); t.Interpreter = 'latex'; t.FontSize = 13;
box on

set(gcf,'color','w')
end
%% PS1 scan
if XP == 7
figure
[a,b] = sort(PS1);
subplot(2,2,1); hold on
for p = 3
h(p) = errorbar(PS1(b),flux{p}(b),dflux{p}(b));
end
% xlim([80,200])
ylim([0,15]*1e23)
t = title('${\Gamma}_+$ $[m^{-2}s^{-1}]$ '); t.Interpreter = 'latex'; t.FontSize = 13;
box on

subplot(2,2,2); hold on
for p = 3
h(p) = plot(PS1(b),pg1_ss{p}(b),'ko-')
end
% xlim([80,200])
ylim([0,2.5])
t = title('$P_{Target}$ $[Pa]$ '); t.Interpreter = 'latex'; t.FontSize = 13;
box on

subplot(2,2,3); hold on
for p = 3
h2 = plot(PS1(b),pg2_ss{p}(b),'ko-')
h4 = plot(PS1(b),pg4_ss{p}(b),'ro-')
end
legend([h2,h4],'PG 2.5','PG 4.5')
% xlim([80,200])
ylim([0,10])
t = title('$P_{Source}$ $[Pa]$ '); t.Interpreter = 'latex'; t.FontSize = 13;
box on

subplot(2,2,4); hold on
for p = 3
h(p) = plot(PS1(b),pg3_ss{p}(b),'ko-')
end
% xlim([80,200])
ylim([-0.5,2.5])
t = title('$P_{CC}$ $[Pa]$ '); t.Interpreter = 'latex'; t.FontSize = 13;
box on

set(gcf,'color','w')
end
%% Radial scans
if XP == 3 || XP == 5
    close all
switch XP
    case 5
        rOffset = 6.25;
    case 3
        rOffset = 6.0;
end
R = x-rOffset;
[a,b] = sort(R);

figure; 
subplot(2,1,1); hold on
for p = 3
h(p) = errorbar(R(b),Ifp{p}(b),dIfp{p}(b))
end
% h(1).Color = 'k';
% h(2).Color = 'r';
% h(3).Color = 'g';
% h(4).Color = 'bl';
% h(1).Marker = 'sq';
% h(2).Marker = 'o';
% h(3).Marker = '^';
% h(4).Marker = '*';

ylim([0,0.5])
xlim([-7.5,7.5]); 
box on
grid on
t = title('$I_+$ $[A]$ '); t.Interpreter = 'latex'; t.FontSize = 13;

subplot(2,1,2); hold on
for p = 3
h(p) = errorbar(R(b),flux{p}(b),dflux{p}(b),'LineWidth',1)
end
% h(1).Color = 'k';
% h(2).Color = 'r';
% h(3).Color = 'g';
% h(4).Color = 'm';
% h(1).Marker = 'sq';
% h(2).Marker = 'o';
% h(3).Marker = '^';
% h(4).Marker = '*';

ylim([0,1.5]*1e24)
xlim([-7.5,7.5]); 
box on
grid on
t = title('${\Gamma}_+$ $[m^{-2}s^{-1}]$ '); t.Interpreter = 'latex'; t.FontSize = 13;

set(gcf,'color','w')

figure; 
subplot(2,1,1); hold on
for p = 3
h(p) = errorbar(R(b),Tefp{p}(b),dTefp{p}(b))
end
% h(1).Color = 'k';
% h(2).Color = 'r';
% h(3).Color = 'g';
% h(4).Color = 'bl';
% h(1).Marker = 'sq';
% h(2).Marker = 'o';
% h(3).Marker = '^';
% h(4).Marker = '*';

ylim([0,10])
xlim([-7.5,7.5]); 
box on
grid on
t = title('$T_e$ $[eV]$ '); t.Interpreter = 'latex'; t.FontSize = 13;
end
%% Convert data into Excel
% Save data to Excel spreadsheet:
% Create table to export to EXCEL
if 0
G = [shotlist(b)',PS1(b)',ne{p}(b)',dne{p}(b)',Te{p}(b)',dTe{p}(b)'];
F = {'Shot','R [cm]','ne[m^-3]','dne','Te[eV]','dTe'};
FileName = 'NeTe_Spool_10_5_2018_01_16.xlsx';
xlswrite(FileName,[F;num2cell(G)]);
end
%% Calculating the total ion flux across the surface

if XP == 3 || XP == 5
    r = R(b)/100;
    f = flux{p}(b);
    switch XP
        case 5
            rq = linspace(-2.25,2.25,500)/100;
        case 3
            rq = linspace(-2.59,2.59,500)/100;
    end
    drq = rq(2)-rq(1);
%     fq = interp1(r,f,rq,'pchip','extrapolation');
    fq = interp1(r,f,rq,'cubic');


    % method 1:
    fq_p = fq(rq>0);
    rq_p = rq(rq>0);

    fq_n = fq(rq<0);
    rq_n = rq(rq<0);

    Fp = fq_p.*rq_p;
    Fn = fq_n.*rq_n;

    figure;
    subplot(2,1,1);
    hold on
    for p = 3
        h(p) = errorbar(R(b),flux{p}(b),dflux{p}(b),'LineWidth',1)
    end
    plot(rq*1e2,fq)
    grid on
    box on
    
    subplot(2,1,2);
    hold on
    plot(rq_p*1e2,Fp,'k.-')
    plot(rq_n*1e2,-Fn,'r.-')
    grid on
    box on

    Sp = +2*pi*trapz(Fp)*drq;
    Sn = -2*pi*trapz(Fn)*drq;

    S = (Sp+Sn)/2;
    title([num2str(S*1e-20),'e20 s^{-1} ions'])
end





