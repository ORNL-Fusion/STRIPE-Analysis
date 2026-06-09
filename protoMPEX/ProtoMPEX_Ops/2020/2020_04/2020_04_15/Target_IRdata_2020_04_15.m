% Objective:
% Preview the IR data as fast as practically possible for use during
% experiments

clear all
close all
clc

% Select shot series to plot:
% =========================================================================
xp = 1;
rotate_flag = 1;
switch xp
    case 1 
        shot       = 29000 + [822];
        rfPwrNet   =         [100];
        addressLoc =         [  4];
        viewType   = 'WindowLim'; 
        
    case 2 
    case 3

end

%% Seq. files
% Load seq files:
% =========================================================================
a{1} = [cd,'\IR_TargetData'];
a{2} = ['\\mpexserver\ProtoMPEX_Data\IR_Camera\2020_03_31'];
a{3} = ['\\mpexserver\ProtoMPEX_Data\IR_Camera\2020_04_08'];
a{4} = ['\\mpexserver\ProtoMPEX_Data\IR_Camera\2020_04_15'];

disp('Extracting data ...')
dum1 = tic;
for s = 1:length(shot)
    % Specify IR file name
    pathName = [a{addressLoc(s)},'\'];
    fileName = ['Shot-0',num2str(shot(s)),'.seq'];
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
    
    % Crop data:
    for ii = 1:size(g{s},3)
        f{s}(:,:,ii) = g{s}(:,130:400,ii);
    end
        
    % Extract data during plasma pulse:
    options.n_Before = 10;
    options.n_After = 50;
    [intf{s},rng{s},t0(s)] = GetDataDuringPlasma(f{s},t_f{s},options);
    rng{s} = 1:size(f{s},3);
          
    % Convert intensity to temperature:
    p = IntensityTempConv(thermalParam.emissivity,f{s}(:,:,rng{s}),seq{s});
    t_temperature{s} = t_f{s}(rng{s});
    
    % Rotate image:
    if rotate_flag
        angle = -45;
        [r] = rotate_IR_data(p,angle);
            % Crop data:
        for ii = 1:size(r,3)
            q(:,:,ii) = r(50:300,98:272,ii);
        end
        temperature{s} = q;
    else
        temperature{s} = p;
    end

    
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
        options.frames = 95:1:size(temperature{s},3)-25;
        options.zlim = [0,120];
        PlayMovieFromArray(temperature{s},options)
    case 2
        options.frames = 0:2:size(f{s},3);
        options.zlim = [0,2000];
        PlayMovieFromArray(f{s},options)
end