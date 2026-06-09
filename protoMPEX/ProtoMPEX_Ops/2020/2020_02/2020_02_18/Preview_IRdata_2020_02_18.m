% Objective:
% Preview the IR data as fast as practically possible for use during
% experiments

clear all
close all
clc

% Select shot series to plot:
% =========================================================================
xp = 2;
switch xp
    case 1 
        shot       = 29000 + [172];
        rfPwrNet   =         [100];
        addressLoc =         [  3];
        viewType   = 'WindowLim'; 
    case 2 
        shot       = 29000 + [184];
        rfPwrNet   =         [100];
        addressLoc =         [  3];
        viewType   = 'MPEXLim'; 
    case 3

end

%% Seq. files
% Load seq files:
% =========================================================================
a{1} = [cd,'\IR_Data'];
a{2} = ['\\mpexserver\ProtoMPEX_Data\IR_Camera\2020_02_11'];
a{3} = ['\\mpexserver\ProtoMPEX_Data\IR_Camera\2020_02_18'];

disp('Extracting data ...')
dum1 = tic;
for s = 1:length(shot)
    % Specify IR file name
    pathName = [a{addressLoc(s)},'\'];
    fileName = ['Shot-0',num2str(shot(s)),'.seq'];
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
s = 1;
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
plotType = 1;
switch plotType
    case 1
        options.frames = 15:2:size(temperature{s},3);
        options.zlim = [0,145];
        PlayMovieFromArray(temperature{s},options)
    case 2
        options.frames = 40:2:size(f{s},3);
        options.zlim = [0,2000];
        PlayMovieFromArray(f{s},options)
end