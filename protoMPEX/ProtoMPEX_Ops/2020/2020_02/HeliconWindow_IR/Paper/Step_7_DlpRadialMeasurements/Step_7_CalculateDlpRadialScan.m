% Create the DLP 1.5 radial scan for the IR paper
% DLP 1.5 with window limiting configuration showing the decay length of
% the plasma relative to the Last Uninterrupted Flux Surface (LUFS)
% =========================================================================
clear all
close all

figureName{1} = 'Step_7_DlpRadialScan';
figureName{2} = 'Step_7_ElectronPressure';

% Data plotting options:
% =========================================================================
FitShow = 0; 
RawDataShow = 0;

% 2018_03_07 DLP 1.5 radial scan:
% =========================================================================
shotlist = [20000 + [194,195,196,197 ,198 ,199 ,200 ,201 ,202 ,203 , 204 ,205 ,206,207 ,208,210,211]];
R        =          [8.5,9.0,9.5,10.0,10.5,11.0,11.5,12.0,12.5,13.0,13.5,14.0,14.5,15.0,8.0,7.5,8.5];

% Calculation setup:
% =========================================================================
SweepType = 'iso'; % "iso" for isolated sweep, "niso" for non-isolated sweep
AttType = 'Vx1,Ix1'; % Attenuation on the digitized signals
DLPType = '1MP';
Config.FitFunction = 2; 
ChannelType = '5' ;
Config.tStart = 4.15; % [s]
Config.tEnd = 4.7;
Config.Center_V = 0; % Remove offset on V: 1 (yes) 0(no)
Config.Center_I = 0; % Remove offset on I: 1 (yes) 0(no)
    

%% DLP FITTING ROUTINE
% Select DLP type:
% =========================================================================
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
        DLP = 2.5;
        Config.L_tip = 1.8/1000;
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

% Select attenuator types:
% =========================================================================
switch AttType
    case 'Vx2,Ix5'
        Config.V_Att = 2;  % Output voltage of DLP box (Voltage) = V_Att*Digitized data 
        Config.I_Att = 5;  % Output voltage of DLP box (Current) = I_att*Digitized data
    case 'Vx1,Ix1'
        Config.V_Att = 1;  % Output voltage of DLP box (Voltage) = V_Att*Digitized data 
        Config.I_Att = 1;  % Output voltage of DLP box (Current) = I_att*Digitized data
end

% Select address to fetch data:
% =========================================================================
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

% Select the amplification factors of the Voltage sweep:
% =========================================================================
switch SweepType
    case 'iso'
Config.V_cal = [(0.46e-3)^-1,0]; % Voltage output of DLP = V_cal(1)*Output voltage of DLP box + V_cal(2) %  2.1739e+03
Config.I_cal = [-1,0]; % Current output of DLP = (I_cal(2) + Output voltage of DLP box)/Ical(1)
    case 'niso'
Config.V_cal = [12.05,0.205];   % Voltage output of DLP = V_cal(1)*Output voltage of DLP box + V_cal(2)
Config.I_cal = [-142.5, -0.6 + 0.015]; % Current output of DLP = (I_cal(2) + Output voltage of DLP box)/Ical(1)
Config.I_cal = [-142.5, 0.015]; % Current output of DLP = (I_cal(2) + Output voltage of DLP box)/Ical(1)
end

% Select additional DLP analysis options:
% =========================================================================
Config.FilterDataInput = 1; % Filter input data with savitsky Golay filter order 3 frame size 7
Config.SGF = 11;
Config.TimeMode = 2; % Effective time of sweep, (1) start of ramp or (2) mean time of ramp
Config.AMU = 2; % Ion mass in AMU
Config.AreaType = 1; % 1: Cylindrical + cap

% Apply the DLP fitting routine:
% =========================================================================
[ni,Te,time,Ifit,Ip,Vp,tm,Vsweep,Isweep,GlitchFlag,SSQres,StdRes,StdResNorm] = DLP_fit_V5_5(Config,shotlist,DataAddress);

% Calculate the mean of the plasma density:
% =========================================================================
for s = 1:length(shotlist)
    Ni{s} = 0.5*(ni{s}{1} + ni{s}{2});
end

% Extract the RF power trace:
% =========================================================================
DA{1} = [RootAddress,'RF_FWD_PWR'];
[RF,t_rf]   = my_mdsvalue_v2(shotlist,DA(1))

%% Select good fits:
% =========================================================================
for s = 1:length(shotlist)
    % Select good fits based on certain conditions:
    GoodFits{s} = GlitchFlag{s} == 0 & StdResNorm{s}<=0.09 & Ni{s}>0 & Ni{s}<1e21 & Te{s}<=30;
            
    % Extract good fits:
    % Plasma density:
    t_Ni_gf{s} = time{s}(GoodFits{s});
    Ni_gf{s} = Ni{s}(GoodFits{s});
    % Electron temperature:
    t_Te_gf{s} = time{s}(GoodFits{s});
    Te_gf{s} = Te{s}(GoodFits{s});
    
    % Derived quantities:
    Pe_gf{s} = e_c.*Ni_gf{s}.*Te_gf{s};
end

%% Plot time dependent data:
% =========================================================================
close all
C = {'k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl'};
t = t_zero(shotlist);

TimePlotStart = 4.1;
TimePlotEnd = 4.7;

% Plot ni, Te and Pe:
% =========================================================================
figure('color','w'); 
% Plasma density:
subplot(2,2,1); hold on
for s = 1:length(shotlist)
    h(s) = plot(t_Ni_gf{s},Ni_gf{s},C{s},'lineWidth',2);
    plot(t_rf{s}(1:length(RF{s})),(RF{s}.^2)*15e19,C{s},'lineWidth',0.5)
end
title('$ n_e $ $ [m^{-3}] $','interpreter','Latex','FontSize',13,'Rotation',0)
legend(h,['DLP ',num2str(DLP)],'location','SouthEast')
ylim([0,8e19])
xlim([TimePlotStart,TimePlotEnd])

% Electron temperature:
subplot(2,2,2); hold on
for s = 1:length(shotlist)
        h(s) = plot(t_Te_gf{s},Te_gf{s},C{s},'lineWidth',2);
        L{s} = [num2str(shotlist(s)),' ,t=',num2str(t{s}(10:14))];
end
legend(h,L,'location','NorthEast')
title('$ T_e $ $ [eV] $','interpreter','Latex','FontSize',13,'Rotation',0)
ylim([0,6])
xlim([TimePlotStart,TimePlotEnd])

% Electron pressure:
subplot(2,2,3); hold on
for s = 1:length(shotlist)
    plot(t_Ni_gf{s},Pe_gf{s},C{s},'lineWidth',2)
end
ylim([0,30])
xlim([TimePlotStart,TimePlotEnd])
title('$ P_e $ $ [Pa] $','interpreter','Latex','FontSize',13,'Rotation',0)

set(findobj('-Property','YTick'),'box','on')
set(gcf,'color','w')

% Plot raw voltage and current data:
% =========================================================================
if RawDataShow
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

%% Calculate steady state values:
% =========================================================================
close all

% Extract steady-state values:
% =========================================================================
for s = 1:length(shotlist)
    rngT{1} = [time{s}>=4.20 & time{s}<=4.23]; 
    rngT{2} = [time{s}>=4.25 & time{s}<=4.35];  
    rngT{3} = [time{s}>=4.50 & time{s}<=4.6]; 
    for p = 1:3
        ness(s,p)  = mean(Ni{s}(rngT{p} & GoodFits{s}))';
        dness(s,p) = std(Ni{s}(rngT{p} & GoodFits{s}),1,2)';
        Tess(s,p)  = mean(Te{s}(rngT{p} & GoodFits{s}))';
        dTess(s,p) = std(Te{s}(rngT{p} & GoodFits{s}),1,2)';
    end
end

% Select the data to plot:
% =========================================================================
rng = 1:3;
ne_SS  = mean(ness(:,rng) ,2,'omitnan');
dne_SS = mean(dness(:,rng),2,'omitnan');
Te_SS  = mean(Tess(:,rng) ,2,'omitnan');
dTe_SS = mean(dTess(:,rng),2,'omitnan');

% Sort radial position:
% =========================================================================
[a,b] = sort(R);

% Plot data:
% =========================================================================

% Plasma density and temperature:
% -------------------------------------------------------------------------
figure('color','w')
hold on

% Draw LUFS:
rDLP = 2.73;
fillColor = [0.9 0.9 0.9];
fillColor = [0.8 1.0 0.8];
hPatch = patch(rDLP*[-1 -1 1 1],12*[0 1 1 0],fillColor);
hPatch.EdgeColor = fillColor;
hL(1) = line(+[1,1]*rDLP,[0,1e20],'LineWidth',1);
hL(2) = line(-[1,1]*rDLP,[0,1e20],'LineWidth',1);
set(hL,'LineStyle','--','Color','k')

% Plasma density:
hP(1) = errorbar(R(b)-11.25,ne_SS(b)*1e-19,dne_SS(b)*1e-19,'LineWidth',1);
hP(1).Marker = 'o'; 
hP(1).MarkerEdgeColor = 'k';
hP(1).Color = 'k';
legendText{1} = '$n_e$ $[\times 10^{19}$ m$^{-3}]$ ';

% Electron temperature:
hP(2) = errorbar(R(b)-11.25,Te_SS(b),dTe_SS(b),'LineWidth',1);
hP(2).Marker = 'sq'; 
hP(2).MarkerEdgeColor = 'r';
hP(2).Color = 'r';
legendText{2} = '$T_e$ [eV] ';

% Limits:
ylim([0,12])
xlim(5*[-1,1]);
xlabel('r [cm]','interpreter','Latex','FontSize',13)
set(gca,'XTick',[-5:1:5])

% Legend:
hLeg = legend(hP,legendText);
hLeg.Interpreter = 'latex'; hLeg.FontSize = 10;

% Final formatting:
box on
grid on
set(gcf,'position',[387.0000  356.3333  300  261.6667])
set(gca,'FontName','times','FontSize',11)

% Save figure:
% =========================================================================
saveas(gcf,figureName{1},'tiffn')

% Electron pressure profile:
% -------------------------------------------------------------------------
figure('color','w')
hold on

% Draw LUFS:
rDLP = 2.73;
fillColor = [0.9 0.9 0.9];
fillColor = [0.8 1.0 0.8];
hPatch = patch(rDLP*[-1 -1 1 1],50*[0 1 1 0],fillColor);
hPatch.EdgeColor = fillColor;
hL(1) = line(+[1,1]*rDLP,[0,50],'LineWidth',1);
hL(2) = line(-[1,1]*rDLP,[0,50],'LineWidth',1);
set(hL,'LineStyle','--','Color','k')

Pe_SS  = e_c*ne_SS.*Te_SS;
dPe_SS = Pe_SS.*sqrt( (dne_SS./ne_SS).^2   +  (dTe_SS./Te_SS).^2  );

hPe(1) = errorbar(R(b)-11.25,Pe_SS(b),dPe_SS(b),'LineWidth',1);
hPe(1).Marker = 'o'; 
hPe(1).MarkerEdgeColor = 'k';
hPe(1).Color = 'k';
legendText{1} = '$P_e$ [Pa]';

if 0
    % Draw exponential decay
    r_SOL = rDLP + linspace(0,2);
    decayLength = 0.5;
    SOL = 25.71*exp(-(r_SOL-rDLP)/decayLength);
    hSOL = plot(r_SOL,SOL,'r','LineWidth',2);
end

% Limits:
ylim([0,50])
xlim(5*[-1,1]);
xlabel('r [cm]','interpreter','Latex','FontSize',13)
set(gca,'XTick',[-5:1:5])

% Legend:
hLeg = legend(hPe,legendText);
hLeg.Interpreter = 'latex'; hLeg.FontSize = 10;

% Final formatting:
box on
grid on
set(gcf,'position',[387.0000  356.3333  300  261.6667])
set(gca,'FontName','times','FontSize',11)

% Save figure:
% =========================================================================
saveas(gcf,figureName{2},'tiffn')

