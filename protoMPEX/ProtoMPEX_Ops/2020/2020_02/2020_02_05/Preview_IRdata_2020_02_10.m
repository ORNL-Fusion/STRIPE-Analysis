% Objective:
% Preview the IR data as fast as practically possible for use during
% experiments

clear all
close all
clc

xp = 2;
switch xp
    case 1 % PREVIEW DATA:
        shot       = 28000 + [928];
        rfPwrNet   =         [80];
        addressLoc =         [1  ];
        viewType   = 'Bottom'; 
    case 2 % BOTTOM VIEW DATA:
        shot       = 28000 + [923,924,926,927,928,929];
        rfPwrNet   =         [73 ,80 ,104,120,129,90 ];
        addressLoc =         [1  ,1  ,1  ,1  ,1  ,1  ];
        viewType   = 'Bottom';   
   case 3 % MIDDLE VIEW DATA:
        shot       = 28000 + [930,931,932,933,934];
        rfPwrNet   =         [80 ,90 ,106,121,129];
        addressLoc =         [1  ,1  ,1  ];
        viewType   = 'Middle';      
    case 4 % TOP VIEW DATA:
        shot       = 28000 + [915,916,917,918,919,920,921];
        rfPwrNet   =         [104,107,119,128,89 ,67 ,76 ];
        addressLoc =         [1  ,1  ,1  ,1  ,1  ,1  ,1  ];
        viewType   = 'Top'; 
   case 5 % Bottom view with gas scan:
        shot       = 28000 + [926,935];
        rfPwrNet   =         [104,106];
        addressLoc =         [1  ,1  ];
        viewType   = 'Bottom view gas scan'; 
   case 6 % Upstream limit:
        shot       = 28000 + [960,962,963,964,967,968,969,971];
        rfPwrNet   =         [67,86  ,102,120,120,120,117,117];
        addressLoc =         [1 ,1   ,1  ,1  ,1  ,1  ,1  ,1  ];
        viewType   = 'Bottom view Upstream limit';         
end

%% 
% =========================================================================
% Load seq files and create RawData
% =========================================================================
a{1} = [cd,'\IR_Data'];
% a{1} = ['\\mpexserver\ProtoMPEX_Data\IR_Camera\2020_02_05'];

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
    
    % Define extraction options:
    extractOptions.frames = 150:1:230;
    
    % Get the data:
    disp(['Extracting data from ',fileName,' ...'])
    [intensityRaw{s},seq{s}] = ExtractDataFromSeqFile(fileName,pathName,thermalParam,extractOptions);
    
    % Convert intensity to temperature:
    temperature{s} = IntensityTempConv(0.9,intensityRaw{s},seq{s});
end
nshots = length(shot);
disp('Data extraction completed!')
disp(['Time taken: ',num2str(toc(dum1))])
%% Preview data
% Select shot:
s = 4;
s = 5;

% Movie rendering options:
try 
    clear options
end
options.frames = 1:2:size(intensityRaw{s},3);
options.colorbar = 1;
options.magnitudePlotMode = 1;
options.removeAxesTicks = 0;
options.shot = shot(s);
options.zlim = [0,20];
options.mirrorImage = 1;

% Select quantity to plot:
plotType = 1;
switch plotType
    case 1
        PlayMovieFromArray(temperature{s},options)
    case 2
        PlayMovieFromArray(intensityRaw{s},options)
end