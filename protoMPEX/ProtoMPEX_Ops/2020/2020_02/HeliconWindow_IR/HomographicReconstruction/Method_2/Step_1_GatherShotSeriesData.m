% Objective:
% Load IR data and create shot series each with associated 6 views 

%% Process: 
% 1- Read the dataset spreadsheet
% 2- Assemble the shot series
% 3- Gather data into structures
% 4- Calculate temperature difference and mirror image
% 5- Save shot series data into .mat files

% Notes:
% In this code, we do not extract the data directly from the .seq files; we
% obtained the data from previosly extracted data stored in .mat files.
% The effect of the mirror is applied in step 4

% clc
clear all
close all

saveData = 1;
saveFig  = 1;

%% 1- Read the dataset spreadsheet
% =========================================================================
spreadsheetName = '\HeliconWindowIR_2020_02_XPs.xlsx';
home = cd;
cd ..
cd ..
addr = pwd;
cd(home)
T = readtable([addr,spreadsheetName],'Sheet',1);
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
for seriesToAnalyze = 1:8
seriesToAnalyze

if numel(seriesToAnalyze) > 1;
    error('Do not analyze more than one shot series at once due to memory restrictions..')
end

for kk = seriesToAnalyze
    rootAddress = 'C:\Users\nfc\Documents\ProtoMPEX_Ops\2020\';
    for jj = 1:numel(shotSeries{kk})
        % shot index:
        n = find(T.shot == shotSeries{kk}(jj));

        % Define address and file name:
        relPath = [T.date{n}]; relPath = [relPath(1:7),'\',relPath];
        fileName = ['\dataset_',num2str(T.dataset(n)),'_IRdata.mat'];
        filePath = [rootAddress,relPath,'\HeliconWindowIR',fileName];

        % load structure of .mat file:
        d = matfile(filePath);

        % Indentify shot needed within .mat file:
        m = find(d.shot == T.shot(n));

        % Extract required data from .mat file:
        f{kk}{jj}.shot = d.shot(1,m);
        f{kk}{jj}.limitMode = d.limitMode(1,m);
        f{kk}{jj}.viewType  = d.viewType(1,m);
        f{kk}{jj}.viewSide  = d.viewSide(1,m);
        f{kk}{jj}.temperature   = cell2mat(d.temperature(1,m));
        f{kk}{jj}.t_temperature = cell2mat(d.t_temperature(1,m));
        f{kk}{jj}.rfPwr = d.X(1,m);
        f{kk}{jj}.t0Plasma = d.t0Plasma(1,m);
        f{kk}{jj}.intf = cell2mat(d.intf(1,m));
        f{kk}{jj}.rngPlasma = cell2mat(d.rngPlasma(1,m));
        f{kk}{jj}.thermalParam = d.thermalParam;
    end
end

%% 4- Calculate temperature difference and mirror image:
% =========================================================================
% Dimensions of microbolometer chip in FLIR A655sc:
pixelsize = 17e-6; %[m]
% Number of pixels:
Nx = 240;
Ny = 640;
% "x" is the vertical axis:
xchip = Nx*pixelsize;
% "y" is the horizontal axis
ychip = Ny*pixelsize;

for kk = seriesToAnalyze;
    for jj = 1:numel(shotSeries{kk})
        % Frame prior to plasma:
        n0 = 1;

        % Calculate dT and mirror image:
        for rr = 1:size(f{kk}{jj}.temperature,3)
            f{kk}{jj}.dT(:,:,rr)   = f{kk}{jj}.temperature(:,end:-1:1,rr) - f{kk}{jj}.temperature(:,end:-1:1,n0);
        end
        
        % Calculate the time base:
        f{kk}{jj}.t_dT = f{kk}{jj}.t_temperature - f{kk}{jj}.t_temperature(1);

        % Clear memory:
        f{kk}{jj}.temperature = [];
        
        % Create chip coordinates:
        % "x" is the vertical axis:
        f{kk}{jj}.xI = linspace(-xchip/2,+xchip/2,Nx)';
        % "y" is the horizontal axis:
        f{kk}{jj}.yI = linspace(-ychip/2,+ychip/2,Ny)';
    end
end

% Plot raw data:
figure
pos = [2,4,6,1,3,5];
for jj = [1:6];
        subplot(3,2,pos(jj))
        surf(f{kk}{jj}.dT(:,:,40),'LineStyle','none')
        box on
        view([0,90])
        set(gca,'PlotBoxAspectRatio',[640/240 1 1])
        xlim([0,640])
        ylim([0,240])
        title([cell2mat(f{kk}{jj}.limitMode),' , '...
            ,cell2mat(f{kk}{jj}.viewSide), ' , RF: ',...
            num2str(rfPower{kk}(jj),4),' kW'],'FontSize',11)
        caxis([0,15])
end

if saveFig
    % figureName = [cell2mat(f{kk}{1}.limitMode),'_Limit_',num2str(round(mean(rfPower{kk})),3),'_kW'];
    figureName = ['ShotSeriesData_',num2str(kk)];
    saveas(gcf,figureName,'tiffn')
end

% return
%% 5- Save data into .mat files:
% =========================================================================

% Description of shot series:
% kk = 1; % MPEX-Limit  , 140 kW:
% kk = 2; % MPEX-Limit  , 120 kW:
% kk = 3; % MPEX-Limit  , 100 kW:
% kk = 4; % MPEX-Limit  , 80  kW:
% kk = 5; % Window-Limit, 130 kW:
% kk = 6; % Window-Limit, 115 kW:
% kk = 7; % Window-Limit, 100 kW:
% kk = 8; % Window-Limit, 80  kW:

if saveData
    t1 = tic;
    disp('Saving data ...')
    fileName = ['ShotSeriesData_',num2str(kk),'.mat'];
    u = f{kk};
    save(fileName,'u')
    t1 = toc(t1);
    disp(['Data Saved!! took ',num2str(t1),' seconds'])
    beep
end
end
disp('End of script!!!')