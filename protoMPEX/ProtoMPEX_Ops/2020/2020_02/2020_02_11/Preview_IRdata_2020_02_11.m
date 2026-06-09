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
    case 1 % PREVIEW DATA:
        shot       = 28000 + [1070,1074];
        rfPwrNet   =         [80  ,80];
        addressLoc =         [1   ,1];
        viewType   = 'Top'; 
    case 2 % Bottom view with PS1 = 3500 A
        shot       = 29000 + [47 ,48 ,49 ,50 ,51 ,52 ,53 ,54 ];
        rfPwrNet   =         [119,128,139,153,162,102,102,102];
        addressLoc =         [1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ];
        viewType   = 'Bottom';   
   case 3 % Bottom view with PS1 = 4500 A:
        shot       = 29000 + [63 ,64 ,65 ,66 ,67 ,68 ];
        rfPwrNet   =         [141,149,158,105,83 ,75 ];
        addressLoc =         [1  ,1  ,1  ,1  ,1  ,1  ];
        viewType   = 'Bottom';      
    case 4 % TOP VIEW DATA:
        shot       = 28000 + [];
        rfPwrNet   =         [];
        addressLoc =         [];
        viewType   = 'Top';
    case 5 % Data comparison:
        shot       = 28000 + [927,1047];
        rfPwrNet   =         [121,119 ];
        addressLoc =         [3  ,2   ];
        viewType   = 'Bottom';
    case 6 % Data comparison:
        shot       = 28000 + [928,1048];
        rfPwrNet   =         [129,128 ];
        addressLoc =         [3  ,2   ];
        viewType   = 'Bottom';
    case 7 % Data comparison:
        shot       = 28000 + [1054,1055,1056];
        rfPwrNet   =         [102 ,102 ,102 ];
        addressLoc =         [1   ,1   ,1   ];
        viewType   = 'Bottom';
    case 8 % Data comparison:
        shot       = 28000 + [1049,1063];
        rfPwrNet   =         [128 ,128  ];
        addressLoc =         [2   ,2    ];
        viewType   = 'Bottom';
    case 9 % Data comparison:
        shot       = 28000 + [1051,1065,928];
        rfPwrNet   =         [128 ,128 ,128];
        addressLoc =         [2   ,2   ,3  ];
        viewType   = 'Bottom';
    case 10 % Data comparison:
        shot       = 28000 + [1047,1057];
        rfPwrNet   =         [119 ,124 ];
        addressLoc =         [1   ,1   ];
        viewType   = 'Bottom';
    case 11 % Data comparison:
        shot       = 28000 + [1057,1061];
        rfPwrNet   =         [124 ,124 ];
        addressLoc =         [1   ,1   ];
        viewType   = 'Bottom';
end

%% Seq. files
% Load seq files:
% =========================================================================
a{1} = [cd,'\IR_Data'];
a{2} = ['\\mpexserver\ProtoMPEX_Data\IR_Camera\2020_02_11'];
a{3} = ['\\mpexserver\ProtoMPEX_Data\IR_Camera\2020_02_05'];

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
    thermalParam.emissivity = 0.9;

    % Define extraction options:
    extractOptions.frames = 100:1:250;
    extractOptions.frameRate = 100;
    
    % Get the data:
    disp(['Extracting data from ',fileName,' ...'])
    [f{s},t_f{s},seq{s}] = ExtractDataFromSeqFile(fileName,pathName,thermalParam,extractOptions);
    
    % Extract data during plasma pulse:
    options.n_Before = 10;
    options.n_After = 50;
    [intf{s},rng{s},t0(s)] = GetDataDuringPlasma(f{s},t_f{s},options);
      
    % Convert intensity to temperature:
    temperature{s} = IntensityTempConv(thermalParam.emissivity,f{s}(rng{s}),seq{s});
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
s = 4;
s = 5;
% s = find(shot == 29065);
s = 5;

% if strcmpi(viewType,'bottom')
% elseif strcmpi(viewType,'top')
% end

% Movie rendering options:
% =========================================================================
try 
    clear options
end
s = 1;
options.colorbar = 1;
options.magnitudePlotMode = 1;
options.removeAxesTicks = 0;
% options.shot = shot(s);
options.mirrorImage = 1;

% Select quantity to plot:
% =========================================================================
plotType = 2;
switch plotType
    case 1
        options.frames = 1:2:size(temperature{s},3);
        options.zlim = [0,17];
        PlayMovieFromArray(temperature{s},options)
    case 2
        options.frames = 40:2:size(f{s},3);
        options.zlim = [0,2000];
        PlayMovieFromArray(f{s},options)
end