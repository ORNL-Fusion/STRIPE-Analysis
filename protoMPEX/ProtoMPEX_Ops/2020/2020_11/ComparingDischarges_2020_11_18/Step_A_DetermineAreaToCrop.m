% Step_A:
% Load data from seq. file
% Rotate frames
% Determine the area to crop

clear all
close all
clc

% Select experiment to load:
% =========================================================================
xp = 7;
rotate_flag = 1;

% Shot cases:
% =========================================================================
switch xp
    case 1
        shot       = 29000 + [777];
        date       = '2020_04_08';
        % Define range of pixels to extract:
        col_rng = [112:256];
        row_rng = [106:249];
        % Notes:
        % Thermal range: 100 to 650 deg     
    case 2
        shot       = 29000 + [778];
        date       = '2020_04_08';
        % Define range of pixels to extract:
        col_rng = [112:256];
        row_rng = [106:249];
        % Notes:
        % Thermal range: 100 to 650 deg      
    case 3
        shot       = 29000 + [824];
        date       = '2020_04_15';
        % Define range of pixels to extract:
        col_rng = [112:256];
        row_rng = [106:249];
        % Notes:
        % Thermal range: -40 to 150 deg
    case 4
        shot       = 29000 + [854];
        date       = '2020_04_15';
        % Define range of pixels to extract:
        col_rng = [112:256];
        row_rng = [106:249];
        % Notes:
        % Thermal range: -40 to 150 deg
    case 5
        shot       = 29000 + [855];
        date       = '2020_04_15';
        % Define range of pixels to extract:
        col_rng = [112:256];
        row_rng = [106:249];
        % Notes:
        % Thermal range: -40 to 150 deg
    case 6
        shot       = 29000 + [691];
        date       = '2020_03_31';
        % Define range of pixels to extract:
        col_rng = [112:256];
        row_rng = [106:249];
        % Notes:
        % Thermal range: -40 to 150 deg
    case 7
        shot       = 30000 + [854];
        date       = '2020_11_17';
        % Define range of pixels to extract:
        col_rng = [112:256] + 40;
        row_rng = [106:249] + 12;
        % Notes:
        % Thermal range: -40 to 150 deg
        % Target is slightly tilted so target appears narrower
end

% Check that col_rng and row_rng are square
disp(['Column dimension: ',num2str(numel(col_rng))]);
disp(['Row dimension: ',num2str(numel(row_rng))]);

% Load seq files and get Temperature data:
% =========================================================================
dataAddress = ['\\mpexserver\ProtoMPEX_Data\IR_Camera\',date];

disp('Extracting data ...')
dum1 = tic;
for s = 1:length(shot)
    % Specify IR file name:
    pathName = [dataAddress,'\'];
    fileName = ['Shot ',num2str(shot(s)),'.seq'];

    % Set thermal parameters:
    thermalParam.ExternalOpticsTransmission = 0.7;
    thermalParam.AtmosphericTemperature = 24;
    thermalParam.Distance = 1;
    thermalParam.ExternalOpticsTemperature = 24;
    thermalParam.ReferenceTemperature = 24;
    thermalParam.Transmission = 1;
    thermalParam.RelativeHumidity = 0;
    thermalParam.ReflectedTemperature = 24;
    thermalParam.emissivity = 0.3;

    % Define extraction options:
    extractOptions.frames = 40:1:250;
    extractOptions.frameRate = 100;
    
    % Get the data:
    disp(['Extracting data from ',fileName,' ...'])
    [g{s},t_f{s},seq{s}] = ExtractDataFromSeqFile(fileName,pathName,thermalParam,extractOptions);
    
    % Display temperature range selected:
    tempRange = seq{1}.ThermalImage.CameraInformation.Range;
    disp(['Temperature range of shot #',num2str(shot(s)),' is: ',num2str(tempRange.Minimum),' to ',num2str(tempRange.Maximum),' deg'])
    
    % Crop data:
    for ii = 1:size(g{s},3)
        f{s}(:,:,ii) = g{s}(:,130:400,ii);
    end
        
    % Extract data during plasma pulse:
    options.n_Before = 10;
    options.n_After = 10;
    [intf{s},~,t0(s)] = GetDataDuringPlasma(f{s},t_f{s},options);
    rngPlasma{s} = 1:size(f{s},3);
          
    % Convert intensity to temperature:
    p = IntensityTempConv(thermalParam.emissivity,f{s}(:,:,rngPlasma{s}),seq{s});
    t_temp{s} = t_f{s}(rngPlasma{s});
    
    % Rotate image:
    if rotate_flag
        angle = -45;
        [r] = rotate_IR_data(p,angle);
        temp{s} = r;
    else
        temp{s} = p;
    end   
end

disp('Data extraction completed!')
disp(['Time taken: ',num2str(toc(dum1))])

% Plot the integral of the intensity over time:
if 1
    figure('color','w')
    grid on
    hold on
    lineColor = {'k','bl','r','g','m','c','k','bl','r','g','m','c'};
    for s = 1:length(shot)
               plot(intf{s}         - min(intf{s}),lineColor{s})
        h(s) = plot(rngPlasma{s},intf{s}(rngPlasma{s}) - min(intf{s}),lineColor{s},'LineWidth',3);
        legendText{s} = num2str(shot(s));
    end
    legend(h,legendText)
    box on
end

%% Select Crop area:
% =========================================================================
% Select shot:
s = 1;

% Select hottest frame:
[~,frame] = max(intf{s}(rngPlasma{s}));

% Select area:
if isempty(col_rng) && isempty(row_rng)
    [row_rng,col_rng,~] = size(temp{s});
    row_rng = 1:row_rng;
    col_rng = 1:col_rng;
end

col_rng = [112:256] + 40;
row_rng = [106:249] + 12;


% Plot region to be cropped:
options.colorbar = 1;
options.magnitudePlotMode = 1;
options.removeAxesTicks = 0;
options.shot = shot(s);
options.mirrorImage = 0;
options.frames = frame;
options.zlim = [0,10];
PlayMovieFromArray(temp{s},options)
hold on

% Define area to be cropped:
x1 = col_rng(1);
y1 = row_rng(1);
drow = max(row_rng) - min(row_rng);
dcol = max(col_rng) - min(col_rng);
x2 = x1 + dcol;
y2 = y1 + drow;

% Draw area to be cropped:
hdum = line([x1,x1,x2,x2,x1],[y1,y2,y2,y1,y1],40*ones(1,5));
hdum.Color = 'k';
hdum.LineWidth = 2;
hdum.LineStyle = '--';

%% Crop data:
% =========================================================================
try 
    clear temp_c
end


% Create cropped frames:
for s = 1:length(shot)
        for ii = 1:size(temp{s},3)
            temp_c{s}(:,:,ii) = temp{s}(row_rng,col_rng,ii);
            dT{s}(:,:,ii) = temp_c{s}(:,:,ii) - temp_c{s}(:,:,1);
        end    
        t_dT{s} = t_temp{s} - t_temp{s}(1);
end

% View cropped frame:
if 0
    s = 1;
    [~,frame] = max(intf{s}(rngPlasma{1}));
    options.colorbar = 1;
    options.magnitudePlotMode = 1;
    options.removeAxesTicks = 0;
    options.shot = shot(s);
    options.mirrorImage = 0;
    options.frames = frame;
    options.zlim = [0,150];
    PlayMovieFromArray(temp_c{s},options)
end

%% Plot data:
% =========================================================================

% Movie rendering options:
% =========================================================================
try 
    clear options
end
options.colorbar = 1;
options.magnitudePlotMode = 1;
options.removeAxesTicks = 0;
options.shot = shot(s);
options.mirrorImage = 0;

% Select quantity to plot:
% =========================================================================
if 1
    plotType = 2;
    switch plotType
        case 1
            options.frames = 95:1:size(temp{s},3)-25;
            options.zlim = [0,150];
            PlayMovieFromArray(temp_c{s},options)
        case 2
            options.frames = 100:2:size(temp{s},3)-2;
            options.zlim = [0,100];
            PlayMovieFromArray(dT{s},options)
    end
end