% Objective:
% Gather Fluoroptic probe data and plot it. 

%% Process: 
% 1- Read the dataset spreadsheet
% 2- Assemble the shot series
% 3- Gather data into structures
% 4- Plot data
% 5- Save figure

% Notes:
% In this code, we do not extract the data directly from the server; we
% obtain the data from previosly extracted data stored in .mat files.

clc
clear all
close all

saveFig  = 1;

figureName{1} = 'Step_8_FluoropticProbe';
figureName{2} = 'Step_8_TotalHeatLoad_FP';
figureName{3} = 'Step_8_Compare_IR_FP';


%% 1- Read the dataset spreadsheet
% =========================================================================
spreadsheetName = '\HeliconWindowIR_2020_02_XPs.xlsx';
targetLocation = 'C:\Users\nfc\Documents\ProtoMPEX_Ops\2020\2020_02\HeliconWindow_IR';
homeLocation = cd;

T = readtable([targetLocation,spreadsheetName],'Sheet',1);
groups = unique(T.Group);

%% 2- Assemble the shot series:
% =========================================================================

% MPEX-Limit, 140 kW:
% Shot 29049 corresponds to the case with PS1 = 3500 A which showed
% significant reduction in heat flux on window. This has approximately 1 cm
% gap between window and LCFS. All other shots that PS1 = 4500 A
% The change in PS1 has only observable on shot the bottom pit view
shotSeries{1} = [29077,29082,29049,29117,29120,29128];
rfPower{1}    = [142  ,142  ,142  ,135  ,150  ,141  ];

% MPEX-Limit, 120 kW:
shotSeries{2} = [29076,29081,29047,29115,29122,29127];
rfPower{2}    = [126  ,126  ,119  ,120  ,116  ,118  ];

% MPEX-Limit, 100 kW:
shotSeries{3} = [29075,29080,29066,29114,29123,29126];
rfPower{3}    = [107  ,102  ,105  ,100  ,97   ,96   ];

% MPEX-Limit, 80 kW:
shotSeries{4} = [29070,29079,29067,29113,29124,29125];
rfPower{4}    = [83   ,88   ,83   ,78   ,74   ,72   ];

% Window-Limit, 130 kW:
shotSeries{5} = [29101,29100,29106,29145,29136,29132];
rfPower{5}    = [130  ,131  ,132  ,132  ,133  ,133  ];

% Window-Limit, 115 kW:
shotSeries{6} = [29102,29099,29105,29144,29137,29131];
rfPower{6}    = [117  ,116  ,116  ,119  ,118  ,119  ];

% Window-Limit, 100 kW:
shotSeries{7} = [29103,29098,29108,29143,29138,29130];
rfPower{7}    = [103  ,104  ,103  ,103  ,103  ,103  ];

% Window-Limit, 80 kW:
shotSeries{8} = [29104,29097,29107,29141,29139,29129];
rfPower{8}    = [88   ,86   ,88   ,88   ,88   ,87  ];

%% 3- Gather data into structures:
% =========================================================================
seriesToAnalyze = 1:8;

for kk = seriesToAnalyze
    rootAddress = 'C:\Users\nfc\Documents\ProtoMPEX_Ops\2020\';
    for jj = 1:numel(shotSeries{kk})
        % shot index:
        n = find(T.shot == shotSeries{kk}(jj));

        % Define address and file name:
        relPath = [T.date{n}]; relPath = [relPath(1:7),'\',relPath];
        fileName = ['\dataset_',num2str(T.dataset(n)),'_DlpGasFp.mat'];
        filePath = [rootAddress,relPath,'\HeliconWindowIR',fileName];
        
        % load .mat file:
        d = load(filePath);

        % Indentify shot needed within .mat file:
        m = find(d.shot == T.shot(n));
        
         f.shot(kk,jj) = d.shot(m);
         f.rfPwr(kk,jj) = d.X(m);
         for ii = 1:4
             f.dataFP_ss{ii}(kk,jj) = d.dataFP_ss{ii}(m);
         end

    end
end

%% 4- Calculate volume-averaged heat load:
% Window dimensions:
Rw = 0.5*12/100; %[m]
tw = 6/1000; % [m]
Lw  = 30/100; %[m]
dt = 0.5;

% Peaking Factor:
% The 8 degree variation on the FP even after 10 seconds after the
% discharge indicates a significant variation on the temperature
% distribution along the window's outer-surface. To account for this we
% introduce a peaking factor.
% Give the variation in the IR derived dQ/dz, a peaking factor of 0.7
% appears to be adequate.
peakFactor = 0.70;

% Volume:
V = Lw*2*pi*Rw*tw;

% Material properties:
rho = 3300; 
cp = 740;

% Mean temperatures:
% MPEX-Limit:
shotSeriesRange = 1:4;
% Loop over all FPs:
for kk = 1:4
    dum_rfpwr(:,kk) = mean(f.rfPwr(shotSeriesRange,:),2);
    dum_dT(:,kk) = mean(f.dataFP_ss{kk}(shotSeriesRange,:),2);
end
rfPwr.MpexLimit = mean(dum_rfpwr,2);
dT.MpexLimit = mean(dum_dT,2);

% Window-Limit:
shotSeriesRange = 5:8;
% Loop over all FPs:
for kk = 1:4
    dum_rfpwr(:,kk) = mean(f.rfPwr(shotSeriesRange,:),2);
    dum_dT(:,kk) = mean(f.dataFP_ss{kk}(shotSeriesRange,:),2);
end
rfPwr.WindowLimit = mean(dum_rfpwr,2);
dT.WindowLimit = mean(dum_dT,2);

% Total heat load:
% Window Limit:
P.WindowLimit = peakFactor*rho*V*cp*dT.WindowLimit/dt;
dP.WindowLimit = 0.15*P.WindowLimit;
% MPEX Limit:
P.MpexLimit = peakFactor*rho*V*cp*dT.MpexLimit/dt;
dP.MpexLimit = 0.15*P.MpexLimit;

%% 4- Plot data:
% =========================================================================
close all


% Resolve all 4 FP:
% =========================================================================
figure('color','w')

fontSize1 = 11;
fontSize2 = 11;
fontSizeLeg = 10;
fontSizeLet = 13;

xvec = linspace(60,200);

ht = tiledlayout(2,2);
ht.TileSpacing = 'compact';
ht.Padding = 'compact';
letters = {'(a)','(b)','(c)','(d)'};

for ii = 1:4
    nexttile 
    hold on
    ax(ii) = gca;

    % MPEX-Limit:
    rng = 1:4;
    hp1{ii} = plot(f.rfPwr(rng,:),f.dataFP_ss{ii}(rng,:),'ksq','MarkerFaceColor','k');
    plot(xvec,polyval(polyfit(f.rfPwr(rng,:),f.dataFP_ss{ii}(rng,:),1),xvec),'k')
    
    % Window-limit:
    rng = 5:8;
    hp2{ii} = plot(f.rfPwr(rng,:),f.dataFP_ss{ii}(rng,:),'ro','MarkerFaceColor','r');
    plot(xvec,polyval(polyfit(f.rfPwr(rng,:),f.dataFP_ss{ii}(rng,:),1),xvec),'r')

    % labels:
    if ii == 3 || ii == 4
        xlabel('[kW]','Interpreter','latex','FontSize',fontSize1)
    end
    if ii == 1 || ii == 3
        ylabel('$\Delta$T [K]','Interpreter','latex','FontSize',fontSize1)
    end
    
    % Annotations:
    titleText{ii} = ['FP $\#$ ',num2str(ii)];
    titlePos{ii} = [140,4.25];
    hText = text(titlePos{ii}(1),titlePos{ii}(2),titleText{ii},'interpreter','Latex','fontSize',12)
    hText.EdgeColor = 'k';
    hText.Position = titlePos{ii};
    % Letters:
    text(10,35,letters{ii},'FontSize',fontSizeLet,'Interpreter','latex')
    
    % Formatting:
    box on
    ylim([0,40])
    xlim([0,200])
end

% Final formatting:
hL = legend([hp1{4}(1),hp2{4}(1)],'MPEX-Limiter','Window-Limiter');
hL.Location = 'northeast';
hL.Interpreter = 'latex';
hL.FontSize = fontSizeLeg; 
hL.Box = 'off';
set(ax,'FontName','times','FontSize',fontSize2)

% Save figure:
% =========================================================================
saveas(gcf,figureName{1},'tiffn')

% Total heat load:
% =========================================================================
figure('color','w')

fontSizeLabel = 12;
fontSizeAx = 12;
fontSizeLeg = 10;

% Window-Limit:
hold on
hQ_fp(1) = errorbar(rfPwr.WindowLimit,P.WindowLimit*1e-3,dP.WindowLimit*1e-3,'ro','MarkerFaceColor','r');
plot(xvec,polyval(polyfit(rfPwr.WindowLimit,P.WindowLimit*1e-3,1),xvec),'r')

% MPEX-Limit:
hQ_fp(2) = errorbar(rfPwr.MpexLimit,P.MpexLimit*1e-3,dP.MpexLimit*1e-3,'ksq','MarkerFaceColor','k');
plot(xvec,polyval(polyfit(rfPwr.MpexLimit,P.MpexLimit*1e-3,1),xvec),'k')

% Labels:
xlabel('RF [kW]','Interpreter','latex','FontSize',fontSizeLabel)
ylabel(' Total Heat [kW]','Interpreter','latex','FontSize',fontSizeLabel)

% Final formatting:
set(hQ_fp,'MarkerSize',8)
set(gcf,'Position',[ 562   360   358   258])
set(gcf,'Position',[485   281   435   337])
xlim([0,200])
ylim([0,90])
grid on
box on

% Legend:
hL = legend([hQ_fp(1),hQ_fp(2)],'Window-Limiter','MPEX-Limiter');
hL.Location = 'best';
hL.Interpreter = 'latex';
hL.FontSize = fontSizeLeg; 
hL.Box = 'off';
set(gca,'FontName','times','FontSize',fontSizeAx)

% Save figure:
% =========================================================================
saveas(gcf,figureName{2},'tiffn')

% Compare IR time-integrated data with FP heating rate:
% =========================================================================
% Load data:
targetLocation = 'C:\Users\nfc\Documents\ProtoMPEX_Ops\2020\2020_02\HeliconWindow_IR\Paper\Step_10_HeatFluxMaps\';
fileName = 'step_10_timeAveragedHeatingRate_IRdata.mat';
d = load([targetLocation,fileName]);

figure('color','w')

fontSizeLabel = 12;
fontSizeAx = 12;
fontSizeLeg = 11;

% Window-Limit (FP):
hold on
hQ_fp_IR(1) = errorbar(rfPwr.WindowLimit,P.WindowLimit*1e-3,dP.WindowLimit*1e-3,'ro','MarkerFaceColor','r');
plot(xvec,polyval(polyfit(rfPwr.WindowLimit,P.WindowLimit*1e-3,1),xvec),'r')

% Window-Limit (IR):
rng = 5:8;
hQ_fp_IR(2) = errorbar(d.Q_RF(rng),d.Q_timeAverage(rng)*1e-3,d.dQ_timeAverage(rng)*1e-3,'ro','MarkerFaceColor','w');
plot(xvec,polyval(polyfit(d.Q_RF(rng),d.Q_timeAverage(rng)*1e-3,1),xvec),'r--')

% MPEX-Limit:
hQ_fp_IR(3) = errorbar(rfPwr.MpexLimit,P.MpexLimit*1e-3,dP.MpexLimit*1e-3,'ksq','MarkerFaceColor','k');
plot(xvec,polyval(polyfit(rfPwr.MpexLimit,P.MpexLimit*1e-3,1),xvec),'k')

% Window-Limit (IR):
rng = 1:4;
hQ_fp_IR(4) = errorbar(d.Q_RF(rng),d.Q_timeAverage(rng)*1e-3,d.dQ_timeAverage(rng)*1e-3,'ksq','MarkerFaceColor','w');
plot(xvec,polyval(polyfit(d.Q_RF(rng),d.Q_timeAverage(rng)*1e-3,1),xvec),'k--')

% Labels:
xlabel('RF [kW]','Interpreter','latex','FontSize',fontSizeLabel)
ylabel('$\dot{Q}$ [kW]','Interpreter','Latex','FontSize',fontSizeLabel)
title('Time-averaged heating rate $\dot{Q}$','interpreter','Latex','FontSize',fontSizeLabel)

% Final formatting:
set(hQ_fp_IR,'MarkerSize',6)
set(gcf,'Position',[ 562   360   358   258])
set(gcf,'Position',[485   281   435   337])
xlim([0,200])
ylim([0,90])
grid on
box on

% Legend:
hL = legend([hQ_fp_IR],'Window-Limiter (FP)','Window-Limiter (IR)','MPEX-Limiter (FP)','MPEX-Limiter (IR)');
hL.Location = 'best';
hL.Interpreter = 'latex';
hL.FontSize = fontSizeLeg; 
hL.Box = 'off';
set(gca,'FontName','times','FontSize',fontSizeAx)

% Save figure:
% =========================================================================
saveas(gcf,figureName{3},'tiffn')