% Objective:
% =========================================================================
% Using the shotSeries spreadsheet, load .seq files for each shot series
% and save the temperature and MIRRORED dT data and metadata into a .mat file

% Context:
% =========================================================================
% For each power level ~ 80, 100, 120 and 140, there are 6 different views
% of the helicon window. Each group of 6 views is called a "shot series"
% and need to be grouped into one single file

clc
clear all
close all

t0 = tic; % Time total process
disp('Start of task 1 ##################################################');

% Extract data from .seq file or read from .mat file:
% =========================================================================
cmpt = 0;
saveFigure = 0;

% Read dataset spreadsheet:
% =========================================================================
homeAddress = cd;
spreadsheetName = 'Step_1_ShotSeries_HeliconWindowIR_2020_02.xlsx';
T = readtable(spreadsheetName,'Sheet',1);
shotSeries = unique(T.shotSeries);
numberOfShotSeries = numel(shotSeries);

% Select shotSeries to analyse with the script:
% =========================================================================
shotSeriesToAnalyze = 1:8;
% shotSeriesToAnalyze = 7:8;

% Root address:
% =========================================================================
rootAddress = 'C:\Users\nfc\Documents\ProtoMPEX_Ops\2020\2020_02\';

% Destination address:
% =========================================================================
cd Step_2_shotSeriesMatFiles
destAddress = pwd;
cd(homeAddress)

% Get data:
% =========================================================================
% ii is index for shotSeries
% jj is index for shot within a shotSeries

if cmpt
% Variables to save in .mat file:
% -------------------------------------------------------------------------
        varList = {...
            'shot',...
            'scanType',...
            'viewType',...
            'viewSide',...
            'limitMode',...
            'X',...
            't_temperature',...
            'temperature',...
            'thermalParam',...
            'extractOptions',...
            'intensity_mean',...
            't_intensity',...
            't0Plasma',...
            'pulseLength',...
            'xI',...
            'yI',...
            'dT',...
            't_dT',...
            'focalLength'...
            };
          
% Define the name of .mat file:
% -------------------------------------------------------------------------
    for ii = shotSeriesToAnalyze
        % Define the name of the .mat file:
        matFileName{ii}= ['shotSeries_',num2str(shotSeries(ii)),'_IRdata.mat'];

        % Each shot has a unique index relative to the entire spreadsheet:
        [shotIndex{ii},~] = find(T.shotSeries == shotSeries(ii));
        numberOfShotsInShotSeries(ii) = numel(shotIndex{ii});
    end
    
% Define address of each shot:
% -------------------------------------------------------------------------
   for ii = shotSeriesToAnalyze
        % For each shot, get the address where the .seq file is located:
        for jj = 1:numberOfShotsInShotSeries(ii)
            dateOfShot = T.date{shotIndex{ii}(jj)};
            address{ii}{jj} = [rootAddress,dateOfShot,'\HeliconWindowIR\IR_Data'];
        end
   end
   
% Extract the data:
% -------------------------------------------------------------------------
    for ii = shotSeriesToAnalyze
        disp(['shotSeries ',num2str(ii),' -------------------------------------'])
        
        % Start timer to record load time for each shotSeries:
        tii = tic;
        clearvars(varList{:})
        
        for jj = 1:6
            % Start timer to record load time for each shot:
            % -------------------------------------------------------------
            tjj = tic; 
            
            % Metadata:
            % -------------------------------------------------------------
            shot(jj)        = T.shot(shotIndex{ii}(jj));
            scanType(jj)    = T.scanType(shotIndex{ii}(jj));
            viewType(jj)    = T.viewType(shotIndex{ii}(jj));
            viewSide(jj)    = T.viewSide(shotIndex{ii}(jj));
            limitMode(jj)   = T.limitMode(shotIndex{ii}(jj));
            pulseLength(jj) = T.pulseLength(shotIndex{ii}(jj));
            X(jj)           = T.X(shotIndex{ii}(jj));
            
            seqFileName = ['Shot ',num2str( T.shot(shotIndex{ii}(jj)) ),'.seq'];
            pathName = [address{ii}{jj},'\'];

            % Set thermal parameters:
            % -------------------------------------------------------------
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
            % -------------------------------------------------------------
            extractOptions.frames = 100:1:250;
            extractOptions.frameRate = 100;
            
            % Get the intensity data:
            % -------------------------------------------------------------            
            disp(['Extracting data from ',seqFileName,' ...'])
            [intensity{jj},t_intensity{jj},seq{jj}] = ExtractDataFromSeqFile(seqFileName,pathName,thermalParam,extractOptions);
            
            % Extract data during plasma pulse and compute integrated intensity:
            % -------------------------------------------------------------           
            options.n_Before = 10;
            options.n_After = 50;
            [intensity_mean{jj},rngPlasma{jj},t0Plasma(jj)] = GetDataDuringPlasma(intensity{jj},t_intensity{jj},options);
                            
            % Convert intensity to temperature:
            % -------------------------------------------------------------
            dum1 = intensity{jj}(:,:,rngPlasma{jj});
            temperature{jj} = IntensityTempConv(thermalParam.emissivity,dum1,seq{jj});
            t_temperature{jj} = t_intensity{jj}(rngPlasma{jj});
            
            % Calculate dT and mirror image:
            % -------------------------------------------------------------            
            % Frame prior to plasma:
            n0 = 1;
            for kk = 1:size(temperature{jj},3)
                dT{jj}(:,:,kk)   = temperature{jj}(:,end:-1:1,kk) - temperature{jj}(:,end:-1:1,n0);
            end
        
            % Calculate the time base:
            % -------------------------------------------------------------            
            t_dT{jj} = t_temperature{jj} - t_temperature{jj}(n0);
            
            % Dimensions of microbolometer chip in FLIR A655sc:
            % -------------------------------------------------------------            
            pixelsize = 17e-6; %[m]
            focalLength = 24.6*1e-3; % [m]
            % Number of pixels:
            Nx = 240;
            Ny = 640;
            % "x" is the vertical axis:
            xchip = Nx*pixelsize;
            % "y" is the horizontal axis
            ychip = Ny*pixelsize;

            % Create chip coordinates:
            % "x" is the vertical axis:
            xI = linspace(-xchip/2,+xchip/2,Nx)';
            % "y" is the horizontal axis:
            yI = linspace(-ychip/2,+ychip/2,Ny)';
                        
            % Stop timer:
            % -------------------------------------------------------------
            tjj = toc(tjj);
            disp(['Shot ',num2str(shot(jj)),' loaded from .seq file in ',num2str(tjj),' sec']);
        end
        
        tii = toc(tii);
        disp(['Total time to load .seq files from shotSeries: ',num2str(tii),' sec']);
        
% Plot the integrated intensity (intensity_mean):
% -------------------------------------------------------------------------
        figure('color','w')
        grid on
        hold on
        lineColor = {'k','bl','r','g','m','c','k','bl','r','g','m','c'};
        for jj = 1:numberOfShotsInShotSeries(ii)
                        plot(t_intensity{jj}                - t0Plasma(jj),intensity_mean{jj}                - min(intensity_mean{jj}),lineColor{jj})
                h(jj) = plot(t_intensity{jj}(rngPlasma{jj}) - t0Plasma(jj),intensity_mean{jj}(rngPlasma{jj}) - min(intensity_mean{jj}),lineColor{jj},'LineWidth',3);
                legendText{jj} = num2str(shot(jj));
        end
        set(gca,'FontName','Times','FontSize',11)
        title(['ShotSeries: ',num2str(ii)],'interpreter','Latex','FontSize',14)
        xlabel('time [s]','interpreter','Latex','FontSize',14)
        ylabel('$\int{intensity(t,x)}dx^2$','interpreter','Latex','FontSize',14)
        legend(h,legendText)
        xlim([-0.2,1])
        box on
        
 % Saving data plot for the iith shotSeries:
% -------------------------------------------------------------------------
        figureName = ['IntegratedIntensity_ShotSeries_',num2str(ii)];
        if saveFigure
            cd(destAddress)
            saveas(gcf,figureName,'tiffn')
            cd(homeAddress)
        end
        
% Saving the data for the iith shotSeries:
% -------------------------------------------------------------------------
        t1 = tic;
        disp(['Saving data into .mat file for shotSeries ',num2str(ii),' ...'])
        cd(destAddress)
        save(matFileName{ii},varList{:},'-v7.3')
        cd(homeAddress)
        t1 = toc(t1);

        disp(['Data saved for ShotSeries ',num2str(ii),' ... took ',num2str(t1),' sec'])
    end
end

t0 = toc(t0);
disp(['Total calculation time: ',num2str(t0),' sec'])

% Total time to run entire script: ~400 sec
% For each shotSeries, time to load x6 seq files: X sec

% Once .mat files have been created:
% Time to load pre-saved .mat file containing x6 shots: 9 sec

% This means that by using the .mat files instead of using the .seq files
% directly we speed up the process by X times

return
%% Test datasets
close all
clc

% The function matfile allows you to access data inside a .mat file without
% loading the entire .mat file. 
% For example, if only one shot is required from a dataset, we can use this
% functionality to extract only the shot that is required.

% Define shot series and shot to plot:
ii = 5;
jj = 1;

% Assembly fileName:
fileName = ['shotSeries_',num2str(ii),'_IRdata.mat'];

% Extract data:
d = matfile(['shotSeriesMatfiles\',fileName]);
A = cell2mat(d.dT(:,jj));

% Limit mode:
disp(['limitMode: ',cell2mat(d.limitMode(:,jj))])

% View type and side:
disp(['viewType: ',cell2mat(d.viewType(:,jj))])
disp(['viewSide: ',cell2mat(d.viewSide(:,jj))])

% RF power level:
disp(['RF power: ',num2str(d.X(:,jj)),' [kW]'])

% Define options structure:
options.frames = 1:3:size(A,3)-10;
options.colorbar = 1;
options.magnitudePlotMode = 2;
options.removeAxesTicks = 0;
options.shot = d.shot(:,jj);
options.zlim = [0,12];
options.mirrorImage = 0;

% Play movie:
PlayMovieFromArray(A,options)