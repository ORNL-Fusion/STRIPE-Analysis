% Objective:
% Load IR data and create composite IR videos 

%% Process: 
% 1- Read the dataset spreadsheet
% 2- Assemble the shot series
% 3- Gather data into structures
% 4- Calculate temperature difference and mirror image
% 5- Stitch images together

% Notes:
% In this code, we do not extract the data directly from the .seq files; we
% obtained the data from previosly extracted data stored in .mat files.
% The effect of the mirror is applied in step 4

clc
clear all
close all
saveData = 1;
saveFig = 1;

%% 1- Read the dataset spreadsheet
% =========================================================================
spreadsheetName = '\HeliconWindowIR_2020_02_XPs.xlsx';
homeAddress = cd;
dataAddress = 'C:\Users\nfc\Documents\ProtoMPEX_Ops\2020\2020_02\HeliconWindow_IR\HomographicReconstruction\Method_1\CompositeIRdata';
tableAddress = 'C:\Users\nfc\Documents\ProtoMPEX_Ops\2020\2020_02\HeliconWindow_IR';
% cd ..
% cd ..
% cd ..
% addr = pwd;
% cd(homeAddress)
T = readtable([tableAddress,spreadsheetName],'Sheet',1);
groups = unique(T.Group);

%% 2- Assemble the shot series:
% =========================================================================

% MPEX-Limit, 140 kW:
shotSeries{1} = [29077,29082,29049,29117,29120,29128];
rfPower{1}    = [142  ,142  ,142  ,135  ,150  ,141  ];

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

if 0
% Enable if one needs to search for shots:
% Select shots corresponding to 6 views at the same RF power:
% =========================================================================
    viewType = 'bottom';
    limitMode = 'window';
    viewSide = 'non-pit';
    rfPwr = [70,90];
    n = find(strcmpi(T.limitMode,limitMode) & strcmpi(T.viewType,viewType) & strcmpi(T.viewSide,viewSide) & T.X>rfPwr(1) & T.X<rfPwr(2))
    disp( {'ii: ', num2str(T.Group(n(end))), 'jj: ', num2str((n(end))) } )
    T.rfPwrNet(n(end))
    T.shot(n(end))
end
%% 3- Gather data into structures:
% =========================================================================
seriesToAnalyze = 8

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
    end
end

%% 4- Calculate temperature difference and mirror image:
% =========================================================================
% Dimensions of microbolometer chip in FLIR A655sc:
pixelsize = 17e-6; %[m]
% "x" is the vertical axis:
xchip = 240*pixelsize;
% "y" is the horizontal axis
ychip = 640*pixelsize;

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
    end
end

%% 5- Stitch images together:
% =========================================================================
% close all

% Description of shot series:
% kk = 1; % MPEX-Limit  , 140 kW:
% kk = 2; % MPEX-Limit  , 120 kW:
% kk = 3; % MPEX-Limit  , 100 kW:
% kk = 4; % MPEX-Limit  , 80  kW:
% kk = 5; % Window-Limit, 130 kW:
% kk = 6; % Window-Limit, 115 kW:
% kk = 7; % Window-Limit, 100 kW:
% kk = 8; % Window-Limit, 80  kW:

for kk = seriesToAnalyze;
% Plot raw data:
figure
pos = [2,4,6,1,3,5];
for jj = [1:6];
        subplot(3,2,pos(jj))
        surf(f{kk}{jj}.dT(:,:,55),'LineStyle','none')
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

figureName = [cell2mat(f{kk}{1}.limitMode),'_Limit_',num2str(round(mean(rfPower{kk})),3),'_kW'];
saveas(gcf,figureName,'tiffn')

% return

% Shift and crop raw data:
% Data needs to be cropped in the vertical direction and shifted in the
% horizontal direction.
% The amount of vertical cropping and horizonal shift are recorded in the
% following variables:
switch kk
    case {1,2,3,4}
        % Pit side:
        crop_rng{kk}{1}   = [077:240];
        crop_shift{kk}(1) = 4;
        crop_rng{kk}{2}   = [035:172];
        crop_shift{kk}(2) = 0;
        crop_rng{kk}{3}   = [001:162];
        crop_shift{kk}(3) = 5;

        % Non-Pit side:
        viewSide{2} = f{kk}{4}.viewSide;
        crop_rng{kk}{4}   = [063:240];
        crop_shift{kk}(4) = -8;
        crop_rng{kk}{5}   = [065:168];
        crop_shift{kk}(5) = 0;
        crop_rng{kk}{6}   = [001:163];
        crop_shift{kk}(6) = 0; 
                     
    case {5,6,7,8}
        % Pit side:
        viewSide{1} = f{kk}{1}.viewSide;
        crop_rng{kk}{1}   = [077:240];
        crop_shift{kk}(1) = 4;
        crop_rng{kk}{2}   = [037:172];
        crop_shift{kk}(2) = 0;
        crop_rng{kk}{3}   = [001:181];
        crop_shift{kk}(3) = 8;
        
        % Non-Pit side:
       viewSide{2} = f{kk}{4}.viewSide;
       crop_rng{kk}{4}   = [095:240];
       crop_shift{kk}(4) = -10;
       crop_rng{kk}{5}   = [061:195];
       crop_shift{kk}(5) = 0;
       crop_rng{kk}{6}   = [001:159];
       crop_shift{kk}(6) = 0;           
end

% Pit side:
viewSide_comp{kk}{2} = cell2mat(f{kk}{1}.viewSide);
frameSize = [size(f{kk}{1}.dT,3),size(f{kk}{2}.dT,3),size(f{kk}{3}.dT,3)];
for fr = 1:min(frameSize)
       frame{1} = f{kk}{1}.dT(:,:,fr);                   
       frame{2} = f{kk}{2}.dT(:,:,fr);             
       frame{3} = f{kk}{3}.dT(:,:,fr);
       dT_comp{2}(:,:,fr) = [circshift(frame{3}(crop_rng{kk}{3},:),crop_shift{kk}(3),2);...
                             circshift(frame{2}(crop_rng{kk}{2},:),crop_shift{kk}(2),2);...
                             circshift(frame{1}(crop_rng{kk}{1},:),crop_shift{kk}(1),2)];
end
% Non-Pit side:
viewSide_comp{kk}{1} = cell2mat(f{kk}{4}.viewSide);
frameSize = [size(f{kk}{4}.dT,3),size(f{kk}{5}.dT,3),size(f{kk}{6}.dT,3)];
for fr = 1:min(frameSize)
       frame{4} = f{kk}{4}.dT(:,:,fr);                   
       frame{5} = f{kk}{5}.dT(:,:,fr);             
       frame{6} = f{kk}{6}.dT(:,:,fr);
       dT_comp{1}(:,:,fr) = [circshift(frame{6}(crop_rng{kk}{6},:),crop_shift{kk}(6),2);...
                             circshift(frame{5}(crop_rng{kk}{5},:),crop_shift{kk}(5),2);...
                             circshift(frame{4}(crop_rng{kk}{4},:),crop_shift{kk}(4),2)];
end  

if 0
    n = 1;
    options.colorbar = 1;
    options.magnitudePlotMode = 1;
    options.removeAxesTicks = 0;
    options.mirrorImage = 0;
    options.frames = 15:70;
    options.zlim = [0,15];
    PlayMovieFromArray(dT_comp{n},options)
end
 
%% Plot near 360-degree view of thermograph:
figure
for n = 1:2
    subplot(1,2,n)
    % Create new coordinates for the "composite" chip
    [row_cI,col_cI,~] = size(dT_comp{n}); 
    % "x" is the vertical axis:
    xchip = row_cI*pixelsize;
    x_cI{n} = linspace(-xchip/2,+xchip/2,row_cI)';
    % "y" is the horizontal axis
    ychip = col_cI*pixelsize;
    y_cI{n} = linspace(-ychip/2,+ychip/2,col_cI)';
    
    % Select subset of data:
    rng1 = find(x_cI{n}*1e3 > -5 & x_cI{n}*1e3 < +5);
    rng2 = find(y_cI{n}*1e3 > -3 & y_cI{n}*1e3 < +3);

    % Plot surface of the temperature map:
    surf(y_cI{n}(rng2)*1e3,x_cI{n}(rng1)*1e3,dT_comp{n}(rng1,rng2,55),'LineStyle','none')
    box on
    view([0,90])
%     set(gca,'PlotBoxAspectRatio',[col_cI/row_cI 1 1])
    xlim([-ychip/2,+ychip/2]*1e3)
    ylim([-xchip/2,+xchip/2]*1e3)
    title([cell2mat(f{kk}{1}.limitMode),' , ',viewSide_comp{kk}{n}, ' , mean RF: ', num2str(mean(rfPower{kk}),4),' kW'],'FontSize',11)
    caxis([0,25])
    colormap('hot')
    colorbar        
    axis image
end

if saveFig
    % Save figures:
    figureName = ['dT_comp',cell2mat(f{kk}{1}.limitMode),'_Limit_',num2str(round(mean(rfPower{kk})),3),'_kW_',viewSide_comp{kk}{n}];
    saveas(gcf,figureName,'tiffn')
end

if saveData
    % Save data to .mat file:
    fileName = ['dT_comp_shotSeries_',num2str(kk)];
    metadata.shotSeries = kk;
    metadata.rfPwr = rfPower{kk};
    metadata.shots = shotSeries{kk};
    metadata.limitMode = cell2mat(f{kk}{1}.limitMode);
    tic
    disp('Saving data ...')
    save(fileName,'dT_comp','x_cI','y_cI','metadata')
    disp('Saving complete!')
    beep
    toc
end
end