% Step_A:
% Load data from seq. file
% Rotate frames
% Convert to temperature
% Crop the frame to capture W plate only
% Save temperature data

clear all
close all
clc

% Select experiment to load:
% =========================================================================
xp = 1;
rotate_flag = 1;

% Shot cases:
% =========================================================================
switch xp
    case 1 
       shot       =  30698;
       addressLoc =  1;        
end

% Load seq files and get Temperature data:
% =========================================================================
a{1} = cd;

disp('Extracting data ...')
dum1 = tic;
for s = 1:length(shot)
    % Specify IR file name:
    pathName = [a{addressLoc(s)},'\'];
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
        f{s}(:,:,ii) = g{s}(:,120:420,ii);
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
        angle = -47;
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
frame = frame + 40;
frame = 80:3:211;

% Define range of pixels to extract:
col_rng = [157:299];
row_rng = [138:282];

% Plot region to be cropped:
options.colorbar = 1;
options.magnitudePlotMode = 1;
options.removeAxesTicks = 0;
options.shot = shot(s);
options.mirrorImage = 0;
options.frames = frame;
options.zlim = [0,3];
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


% Creat cropped frames:
for s = 1:length(shot)
        for ii = 1:size(temp{s},3)
            temp_c{s}(:,:,ii) = temp{s}(row_rng,col_rng,ii);
            dT{s}(:,:,ii) = temp_c{s}(:,:,ii) - temp_c{s}(:,:,1);
        end    
        t_dT{s} = t_temp{s} - t_temp{s}(1);
end

if 1
    s = 1;
    [~,frame] = max(intf{s}(rngPlasma{1}));
    frame = 80:3:211;
    options.colorbar = 1;
    options.magnitudePlotMode = 1;
    options.removeAxesTicks = 0;
    options.shot = shot(s);
    options.mirrorImage = 0;
    options.frames = frame;
    options.zlim = [0,20];
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
           [~,frame] = max(intf{s}(rngPlasma{1}));
            options.frames = 80:3:211;
            options.zlim = [0,max(max(dT{s}(:,:,frame)))];
            PlayMovieFromArray(dT{s},options)
    end
end

%% Save data:
% =========================================================================

% shot:
d.shot = shot(s);

% thermal parameters:
d.thermalParams = seq{s}.ThermalImage.ThermalParameters;
d.tempRange = [tempRange.Minimum,tempRange.Maximum];

% data:
% Absolute temperature:
d.temp2D = temp{s};
d.t_temp2D = t_temp{s};
% Change in temperature:
d.dT2D = dT{s};
d.t_dT2D = t_dT{s};
% Integrated emission:
d.intf = intf{s};

% Save data:
fileName = ['Step_A_temp2D_',num2str(shot(s)),'.mat'];
save(fileName,'d');