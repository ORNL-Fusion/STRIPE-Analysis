% Objective:
% Preview the IR data as fast as practically possible for use during
% experiments

clear all
close all
clc

% Select shot series to plot:
% =========================================================================
xp = 3;
switch xp
    case 1 
        shot       = 29000 + [499];
        rfPwrNet   =         [100];
        addressLoc =         [  3];
        viewType   = 'WindowLim';         
    case 2 % Helicon window IR with ECH
        shot       = 29000 + [520,521];
        rfPwrNet   =         [100,100];
        addressLoc =         [3  ,  3];
        viewType   = 'WindowLim';    
    case 3
        shot       = 29000 + [520,521,522,524,518];
        rfPwrNet   =         [100,100,100,100,100];
        addressLoc =         [3  ,  3,  3,  3,  3];
        viewType   = 'WindowLim'; 
        
                shot       = 29000 + [524,518];
        rfPwrNet   =         [100,100];
        addressLoc =         [  3,  3];
        viewType   = 'WindowLim'; 
        
        shot       = 29000 + [517,518,520,521,524,525];
        rfPwrNet   =         [100,100,100,100,100,100];
        addressLoc =         [  3,  3,3  ,3  ,3  ,3  ];
        viewType   = 'WindowLim'; 
end

%% Seq. files
% Load seq files:
% =========================================================================
a{1} = [cd,'\IR_Data'];
a{2} = ['\\mpexserver\ProtoMPEX_Data\IR_Camera\2020_02_11'];
a{3} = ['\\mpexserver\ProtoMPEX_Data\IR_Camera\2020_02_28'];

disp('Extracting data ...')
dum1 = tic;
for s = 1:length(shot)
    % Specify IR file name
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
    extractOptions.frames = 50:1:300;
    extractOptions.frameRate = 100;
    
    % Get the data:
    disp(['Extracting data from ',fileName,' ...'])
    [f{s},t_f{s},seq{s}] = ExtractDataFromSeqFile(fileName,pathName,thermalParam,extractOptions);
    
    % Extract data during plasma pulse:
    options.n_Before = 10;
    options.n_After = 50;
    [intf{s},rng{s},t0(s)] = GetDataDuringPlasma(f{s},t_f{s},options);
    rng{s} = 1:size(f{s},3);
      
    % Convert intensity to temperature:
    temperature{s} = IntensityTempConv(thermalParam.emissivity,f{s}(:,:,rng{s}),seq{s});
    t_temperature{s} = t_f{s}(rng{s});
end
nshots = length(shot);
disp('Data extraction completed!')
disp(['Time taken: ',num2str(toc(dum1))])

figure('color','w')
grid on
hold on
lineColor = {'k','bl','r','g','m','c','k','bl','r','g','m','c'};
for s = 1:length(shot)
           plot(t_f{s}      - t0(s)   ,intf{s}         - min(intf{s}),lineColor{s})
    h(s) = plot(t_f{s}(rng{s}) - t0(s),intf{s}(rng{s}) - min(intf{s}),lineColor{s},'LineWidth',3);
    legendText{s} = num2str(shot(s));
end
legend(h,legendText)
box on

% Gather Fluroptic probe data

%% Plot data
% close all
% Select shot:
% =========================================================================
s = 2;
% Movie rendering options:
% =========================================================================
try 
    clear options
end
options.colorbar = 1;
options.magnitudePlotMode = 1;
options.removeAxesTicks = 0;
options.shot = shot(s);
options.mirrorImage = 1;

% Select quantity to plot:
% =========================================================================
plotType = 1;
switch plotType
    case 1
        options.frames = 80:2:size(temperature{s},3);
        options.zlim = [0,20];
        PlayMovieFromArray(temperature{s},options)
    case 2
        options.frames = 0:2:size(f{s},3);
        options.zlim = [0,2000];
        PlayMovieFromArray(f{s},options)
end