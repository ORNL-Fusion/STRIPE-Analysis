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

figureName = 'Step_6_FpData';

saveData = 1;
saveFig  = 1;

%% 1- Read the dataset spreadsheet
% =========================================================================
spreadsheetName = '\HeliconWindowIR_2020_02_XPs.xlsx';
targetLocation = 'C:\Users\nfc\Documents\ProtoMPEX_Ops\2020\2020_02\HeliconWindow_IR';
homeLocation = cd;

% cd(targetLocation);
% addr = pwd;
% cd(homeLocation)

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

%% 4- Plot data:
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

    % Window-Limit:
    rng = 1:4;
    hp1{ii} = plot(f.rfPwr(rng,:),f.dataFP_ss{ii}(rng,:),'ksq-');
    plot(xvec,polyval(polyfit(f.rfPwr(rng,:),f.dataFP_ss{ii}(rng,:),1),xvec),'k')
    
    % MPEX-limit:
    rng = 5:8;
    hp2{ii} = plot(f.rfPwr(rng,:),f.dataFP_ss{ii}(rng,:),'ro-');
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
hL = legend([hp1{4}(1),hp2{4}(1)],'MPEX-Limiter','Window-Limiter');
hL.Location = 'northeast';
hL.Interpreter = 'latex';
hL.FontSize = fontSizeLeg; 
hL.Box = 'off';
set(ax,'FontName','times','FontSize',fontSize2)

% Save figure:
% =========================================================================
saveas(gcf,figureName,'tiffn')

%% 4- Calculate volume-averaged heat load:
% Window dimensions:
Rw = 0.5*12/100; %[m]
tw = 6/1000; % [m]
Lw  = 30/100; %[m]
dt = 0.5;

% Peaking Factor:
peakFactor = 0.8;

% Volume:
V = Lw*pi*( (Rw+tw)^2 - (Rw)^2 );
V = Lw*2*pi*Rw*tw;

% Material properties:
rho = 3300; 
cp = 740;

% Mean temperatures:
ii = 4;
% Window-Limit:
rng = 1:4;
rfPwr.WindowLimit = mean(f.rfPwr(rng,:),2);
dT.WindowLimit = mean(f.dataFP_ss{ii}(rng,:),2);

% MPEX-Limit:
rng = 5:8;
rfPwr.MpexLimit = mean(f.rfPwr(rng,:),2);
dT.MpexLimit = mean(f.dataFP_ss{ii}(rng,:),2);

% Total heat load:
P.WindowLimit = peakFactor*rho*V*cp*dT.WindowLimit/dt;
P.MpexLimit = peakFactor*rho*V*cp*dT.MpexLimit/dt;

% Total heat load:
% =========================================================================
figure('color','w')

fontSizeLabel = 12;
fontSizeAx = 12;
fontSizeLeg = 10;

% Window-Limit:
hold on
hQ_fp(1) = plot(rfPwr.WindowLimit,P.WindowLimit*1e-3,'ksq','MarkerFaceColor','k');
plot(xvec,polyval(polyfit(rfPwr.WindowLimit,P.WindowLimit*1e-3,1),xvec),'k')

% MPEX-Limit:
hQ_fp(2) = plot(rfPwr.MpexLimit,P.MpexLimit*1e-3,'ro','MarkerFaceColor','r');
plot(xvec,polyval(polyfit(rfPwr.MpexLimit,P.MpexLimit*1e-3,1),xvec),'r')

% Labels:
xlabel('RF [kW]','Interpreter','latex','FontSize',fontSizeLabel)
ylabel(' Total Heat [kW]','Interpreter','latex','FontSize',fontSizeLabel)

% Final formatting:
set(hQ_fp,'MarkerSize',8)
set(gcf,'Position',[ 562   360   358   258])
xlim([0,200])
ylim([0,90])
grid on
box on

% Legend:
hL = legend([hQ_fp(2),hQ_fp(1)],'Window-Limiter','MPEX-Limiter');
hL.Location = 'best';
hL.Interpreter = 'latex';
hL.FontSize = fontSizeLeg; 
hL.Box = 'off';
set(gca,'FontName','times','FontSize',fontSizeAx)
