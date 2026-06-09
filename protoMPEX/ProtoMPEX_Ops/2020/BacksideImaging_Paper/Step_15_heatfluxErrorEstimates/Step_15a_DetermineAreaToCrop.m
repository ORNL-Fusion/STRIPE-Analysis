% Step_15a:
% We are going to use the data from the microwave power scan. This data set
% has the same camera view throughout the entire scan and thus thus the
% image cropping becomes trivial: one set of conditions for all shots.

% The main steps are as follows:
% Extract data from seq file
% Raw data is cropped to reduce file size
% Cropped image is rotated 45 degrees
% 

clear all
close all
clc

% Define shots in dataset:
% =========================================================================
shot        = [29655];
pwr28GHz    = [70.07];

% Select shots to analyse:
shotlist = 1:length(shot);

% Save flags:
% =========================================================================
saveFig = 1;
saveData = 1;

% Define location of IR data:
% =========================================================================
dataAddress = [cd,'\IR_rawData'];
pathName = [dataAddress,'\'];

% Define cropping information:
% =========================================================================
% Define initial column crop range:
columnsToCrop{1} = [130:400];

% Define range of pixels to extract:
col_rng = [112:256];
row_rng = [106:249];

% Check that col_rng and row_rng are square
disp(['Column dimension: ',num2str(numel(col_rng))]);
disp(['Row dimension: ',num2str(numel(row_rng))]);

% Set thermal parameters:
% =========================================================================
thermalParam.ExternalOpticsTransmission = 0.7;
thermalParam.AtmosphericTemperature = 24;
thermalParam.Distance = 1;
thermalParam.ExternalOpticsTemperature = 24;
thermalParam.ReferenceTemperature = 24;
thermalParam.Transmission = 1;
thermalParam.RelativeHumidity = 0;
thermalParam.ReflectedTemperature = 24;

% Set emissivity range:
emissivity = [0.27,0.3,0.33];

% Define frame extraction options:
% =========================================================================
extractOptions.frames = 40:1:250;
extractOptions.frameRate = 100;
options.n_Before = 10;
options.n_After = 10;

% Get IR data:
% =========================================================================
for kk = shotlist
    
    disp('Extracting data ...')
    dum1 = tic;
    
    % Specify IR file name:
    fileName = ['Shot ',num2str(shot(kk)),'.seq'];

    % Get the data:
    disp(['Extracting data from ',fileName,' ...'])
    [IRdata.Raw.intensity,IRdata.Raw.t,seq{kk}] = ExtractDataFromSeqFile(fileName,pathName,thermalParam,extractOptions);

    % Display temperature range selected:
    tempRange = seq{kk}.ThermalImage.CameraInformation.Range;
    disp(['Temperature range of shot #',num2str(shot(kk)),' is: ',num2str(tempRange.Minimum),' to ',num2str(tempRange.Maximum),' deg'])

    % Crop data:
    for ii = 1:size(IRdata.Raw.intensity,3)
        IRdata.Cropped.intensity(:,:,ii) = IRdata.Raw.intensity(:,columnsToCrop{kk},ii);
        IRdata.Cropped.t = IRdata.Raw.t;
    end

    % Calculate integrated intensity vrs time:
    [IRdata.Cropped.integratedIntensity,~,~] = GetDataDuringPlasma(IRdata.Cropped.intensity,IRdata.Cropped.t,options);

    % Find start of plasma pulse based on intensity:
    IRdata.Cropped.framePlasmaStart = find(movmean(diff(IRdata.Cropped.integratedIntensity),10)>2,1);
    
    % Get the frames during plasma pulse:
    IRdata.Cropped.framesDuringPlasma = (IRdata.Cropped.framePlasmaStart-30):(IRdata.Cropped.framePlasmaStart+80);

    % Calculate temperature based on emissivity range:
    for ee = 1:numel(emissivity)
        
        % Integrated IR intensity durnig plasma shot:
        temperature{kk}{ee}.integratedIntensity = IRdata.Cropped.integratedIntensity(IRdata.Cropped.framesDuringPlasma);
        
        % time vector:
        temperature{kk}{ee}.t = IRdata.Cropped.t(IRdata.Cropped.framesDuringPlasma)-IRdata.Cropped.t(IRdata.Cropped.framesDuringPlasma(1));

        % Select emissivity;
        temperature{kk}{ee}.emissivity = emissivity(ee);

        % Display emissivity under use:
        disp(['Using emissivity value of ',num2str(temperature{kk}{ee}.emissivity),' ...'])
    
        % Convert intensity to temperature:
        p = IntensityTempConv(temperature{kk}{ee}.emissivity,IRdata.Cropped.intensity(:,:,IRdata.Cropped.framesDuringPlasma),seq{kk});

        % Rotate image:
        angle = -45;
        temperature{kk}{ee}.value = rotate_IR_data(p,angle);
    end

    disp('Data extraction completed!')
    disp(['Time taken: ',num2str(toc(dum1))])

end

%% Select Crop area:
% =========================================================================
% Select hottest shot:
[~,kkmax] = max(pwr28GHz);

for kk = shotlist
    % Select hottest frame:
    [~,frame] = max(temperature{kk}{1}.integratedIntensity);

    % Plot region to be cropped:
    options.colorbar = 1;
    options.magnitudePlotMode = 1;
    options.removeAxesTicks = 0;
    options.shot = shot(kk);
    options.mirrorImage = 0;
    options.frames = frame;
    options.zlim = [0,550];
    PlayMovieFromArray(temperature{kk}{1}.value,options)
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
    if options.zlim > 10
        hdum.Color = 'w';
    else
        hdum.Color = 'k';
    end
    hdum.LineWidth = 2;
    hdum.LineStyle = '--';
end

%% Crop data and get plate temperature distribution:
% =========================================================================

% Formmatting:
zlimRng  = [0,450];
caxisRng = [0,450];

% Calculate dT:
% =============
for kk = shotlist
    
    figure('color','w')
    
    % Clear previous cropped data if any:
    try 
        clear y
    end

    % Create cropped frames:
    for ee = 1:numel(temperature{kk})
        
        subplot(1,numel(emissivity),ee)
        hold on
        
        for ii = 1:size(temperature{kk}{ee}.value,3)
            
            % Extract plate temperature profile:
            y(:,:,ii) = temperature{kk}{ee}.value(row_rng,col_rng,ii);
            
            % Calculate temperature difference:
            plate{kk}{ee}.dT(:,:,ii) = y(:,:,ii) - y(:,:,1);
        end    
        
        % Time vector:
        plate{kk}{ee}.t = temperature{kk}{ee}.t;
        
        % Emissivity:
        plate{kk}{ee}.emissivity = temperature{kk}{ee}.emissivity;
        
        % Shot:
        plate{kk}{ee}.shot = shot(kk);
        
        % Power:
        plate{kk}{ee}.pwr28GHz = pwr28GHz(kk);
        
        % Integrated IR intensity
        plate{kk}{ee}.integratedIntensity = temperature{kk}{ee}.integratedIntensity;
        
        % Plot data:
        surf(plate{kk}{ee}.dT(:,:,ii) - plate{kk}{ee}.dT(:,:,1),'LineStyle','none')
        
        % Formatting:
        colormap(flipud(hot))
        caxis(caxisRng)
        zlim(zlimRng)
        axis('square')
        box on
        
        % Labels:
        title(['shot: ',num2str(shot(kk)),' emiss: ',num2str(plate{kk}{ee}.emissivity)])
    end
        
    if kk == numel(shotlist)
        hC = colorbar;
        set(hC,'Position',[0.5382    0.1111    0.0213    0.2087])
        hT = text(250,75,'${\Delta}T$ [K]');
        set(hT,'Interpreter','Latex','fontSize',12)
    end
    
    % Save figure:
    % ============
    if saveFig
        figureName = ['Step_15a_plateTemperature_',num2str(shot(kk))];
        saveas(gcf,figureName,'tiffn')
    end
    
end

%% Save data:
if saveData
    varList = {'plate','seq'};
    fileName = 'step_15a_plateTemperature.mat';
    save(fileName',varList{:})
end