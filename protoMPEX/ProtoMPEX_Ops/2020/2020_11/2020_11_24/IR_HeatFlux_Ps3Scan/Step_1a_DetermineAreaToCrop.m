% Step_1a:
% Load data from seq. file
% Rotate frames
% Determine the area to crop

clear all
close all
clc

% All shots in the dataset:
shot     = 30000 + [976,977,979,980,981,982,983,984,985];
PS3      =         [220,250,300,325,350,375,400,425,450];

% Select shots to analyse:
shotlist = 1:length(shot);

% Plot flags:
plotFig_1 = 1;
plotFig_2 = 1;

% Save flags:
saveFig_1 = 1;
saveFig_2= 1;
saveData = 1;

% Load IR data:
% =========================================================================
for kk = shotlist
    
% switch shot(kk)
%     case {shot}
% end

% Apply rotation:
rotate_flag = 1;

% Define initial column crop range:
colCropRng{kk} = [130:400];

% Define range of pixels to extract:
col_rng{kk} = [112:256] + 40;
row_rng{kk} = [106:249] + 10;

% Check that col_rng and row_rng are square
disp(['Column dimension: ',num2str(numel(col_rng))]);
disp(['Row dimension: ',num2str(numel(row_rng))]);

% Load seq files and get Temperature data:
dataAddress = [cd,'\IR_rawData'];

disp('Extracting data ...')
dum1 = tic;
% Specify IR file name:
pathName = [dataAddress,'\'];
fileName = ['Shot ',num2str(shot(kk)),'.seq'];

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
[g,t_f,seq{kk}] = ExtractDataFromSeqFile(fileName,pathName,thermalParam,extractOptions);

% Display temperature range selected:
tempRange = seq{kk}.ThermalImage.CameraInformation.Range;
disp(['Temperature range of shot #',num2str(shot(kk)),' is: ',num2str(tempRange.Minimum),' to ',num2str(tempRange.Maximum),' deg'])

% Crop data:
for ii = 1:size(g,3)
    f(:,:,ii) = g(:,colCropRng{kk},ii);
end

% Extract data during plasma pulse:
options.n_Before = 10;
options.n_After = 10;
[intf] = CalculateIntegratedIntensity(f,t_f,options);
ii0 = find(movmean(diff(intf),10)>2,1);
rngPlasma{kk} = (ii0-30):(ii0+80);

% Integrated intensity durnig plasma shot:
intf_c{kk} = intf(rngPlasma{kk});

% Creating the time vector of the cropped intf_c data:
dt = t_f(2) - t_f(1);
N = numel(intf_c{kk});
t_temp{kk} = 0:dt:dt*(N-1);

% Convert intensity to temperature:
p = IntensityTempConv(thermalParam.emissivity,f(:,:,rngPlasma{kk}),seq{kk});

% Rotate image:
if rotate_flag
    angle = -45;
    [r] = rotate_IR_data(p,angle);
    temp{kk} = r;
else
    temp{kk} = p;
end   

disp('Data extraction completed!')
disp(['Time taken: ',num2str(toc(dum1))])

end

% Plot the integral of the intensity over time:
% =========================================================================
if plotFig_1
    figure('color','w')
    lineColor = {'k','bl','r','g','m','c','k','bl','r','g','m','c'};
    grid on
    hold on
    for kk = shotlist
        h(kk) = plot(intf_c{kk}  - min(intf_c{kk}),lineColor{kk},'LineWidth',3);
        legendText{kk} = ['kk = ',num2str(kk),', ',num2str(shot(kk)),', ECH power: ',num2str(PS3(kk))];
    end
    hL = legend(h,legendText);
    set(hL,'Location','best')
    box on
    
    % Save figure:
    if saveFig_1
        figureName = 'Step_1a_IntegratedIntensityVsTime';
        saveas(gcf,figureName,'tiffn')
    end
end

%% Select Crop area:
% =========================================================================
% Select hottest shot:
[~,kkmax] = max(PS3);

for kk = shotlist
    % Select hottest frame:
    [~,frame] = max(intf_c{kkmax});
    
    % Plot region to be cropped:
    options.colorbar = 1;
    options.magnitudePlotMode = 1;
    options.removeAxesTicks = 0;
    options.shot = shot(kk);
    options.mirrorImage = 0;
    options.frames = frame;
    options.zlim = [0,40];
    PlayMovieFromArray(temp{kk},options)
    hold on

    % Define area to be cropped:
    x1 = col_rng{kk}(1);
    y1 = row_rng{kk}(1);
    drow = max(row_rng{kk}) - min(row_rng{kk});
    dcol = max(col_rng{kk}) - min(col_rng{kk});
    x2 = x1 + dcol;
    y2 = y1 + drow;

    % Draw area to be cropped:
    hdum = line([x1,x1,x2,x2,x1],[y1,y2,y2,y1,y1],40*ones(1,5));
    if options.zlim > 10
        hdum.Color = 'w';
    else
        hdum.Color = 'k';
    end
    hdum.LineWidth = 2;
    hdum.LineStyle = '--';
end

%% Crop data:
figure('color','w')

zlimRng  = [0,60];
caxisRng = [0,60];

for kk = shotlist
    subplot(3,4,kk)
    hold on
    
    % Clear previous cropped data if any:
    try 
        clear temp_c
    end

    % Create cropped frames:
    for ii = 1:size(temp{kk},3)
        temp_c{kk}(:,:,ii) = temp{kk}(row_rng{kk},col_rng{kk},ii);
        dT{kk}(:,:,ii) = temp_c{kk}(:,:,ii) - temp_c{kk}(:,:,1);
    end    
    t_dT{kk} = t_temp{kk} - t_temp{kk}(1);

    % Plot data:
    surf(dT{kk}(:,:,ii)-dT{kk}(:,:,1),'LineStyle','none')
    
    % Formatting:
    colormap(flipud(hot))
    caxis(caxisRng)
    zlim(zlimRng)
    axis('square')
    box on
   
    % Labels:
    title(['shot: ',num2str(shot(kk))])
    
    if kk == numel(shotlist)
        hC = colorbar;
        set(hC,'Position',[0.5382    0.1111    0.0213    0.2087])
        hT = text(250,75,'${\Delta}T$ [K]');
        set(hT,'Interpreter','Latex','fontSize',12)
    end
    
    % Need to add the flux mapping to the data
end

% Save figure:
if saveFig_2
    figureName = 'Step_1a_TemperatureMaps';
    saveas(gcf,figureName,'tiffn')
end

%% Save data:
if saveData
    varList = {'shot','PS3','seq','dT','t_dT','intf_c'};
    fileName = 'step_1a_dT_Cropped.mat';
    save(fileName',varList{:})
end