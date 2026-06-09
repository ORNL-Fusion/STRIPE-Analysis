% Task 1: Load .seq files and create temperature data

% note:
% The data extracted and saved into the .mat files has not been corrected
% for the effect of the mirror periscope. 
% The best way to resolve this is to extract the temperature data and then
% calculate the dT data and applying the mirroring (column-wise flip) in
% that operation,

clc
clear all
close all

t0 = tic; % Time total process
disp('Start of task 1 ##################################################');

% Extract data from .seq file or read from .mat file:
% =========================================================================
cmpt = 1;
saveFigure = 1;

% Read dataset spreadsheet:
% =========================================================================
spreadsheetName = 'HeliconWindowIR_2020_02_11.xlsx';
T = readtable(spreadsheetName,'Sheet',1);
datasets = unique(T.dataset);
numberOfDatasets = numel(datasets);

datasetToAnalyze = 7;
datasetToAnalyze = 1:9;

% Root address:
% =========================================================================
currentFolder = pwd;
dummy1 = findstr(currentFolder,'\');
rootAddress = currentFolder(1:dummy1(end));

% ii is index for dataset
% jj is index for shot within a dataset

% Get data:
% =========================================================================
if cmpt
    % Variables to save in .mat file:
    % ---------------------------------------------------------------------
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
            'intf',...
            't_f',...
            'rngPlasma',...
            't0Plasma',...
            'pulseLength'...
            };
          
    % Define the name of .mat file:
    % ---------------------------------------------------------------------
    for ii = datasetToAnalyze
        % Define the name of the .mat file:
        matFileName{ii}= ['dataset_',num2str(datasets(ii)),'_IRdata.mat'];

        % For each data set, assign an index to each of the shots:
        [shotIndex{ii},~] = find(T.dataset == datasets(ii));
        numberOfShotsInDataset(ii) = numel(shotIndex{ii});
    end
    
    % Define address of each shot:
    % ---------------------------------------------------------------------
   for ii = datasetToAnalyze
        % For each shot, get the address where the .seq file is located:
        getAddressMethod = 3;
        switch getAddressMethod
            case 1 
                    for jj = 1:numberOfShotsInDataset(ii)
                        dateOfShot = T.date{shotIndex{ii}(jj)};
                        address{ii}{jj} = [rootAddress,dateOfShot(1:end-3),'\',dateOfShot,'\IR_RawData'];
                    end
            case 2
                    for jj = 1:numberOfShotsInDataset(ii)
                        dateOfShot = T.date{shotIndex{ii}(jj)};
                        address{ii}{jj} = [rootAddress,dateOfShot,'\IR_Data'];
                    end
            case 3
                    for jj = 1:numberOfShotsInDataset(ii)
                        dateOfShot = T.date{shotIndex{ii}(jj)};
                        address{ii}{jj} = [pwd,'\IR_Data'];
                    end
        end
   end

    % Extract the data:
    % ---------------------------------------------------------------------
    for ii = datasetToAnalyze
        disp(['Dataset ',num2str(ii),' -------------------------------------'])
        tii = tic;
        clearvars(varList{:})
        
        % Extract raw intensity data, subset and convert to temperature:
        for jj = 1:numberOfShotsInDataset(ii)
            shot(jj) = T.shot(shotIndex{ii}(jj));
            scanType(jj) = T.scanType(shotIndex{ii}(jj));
            viewType(jj) = T.viewType(shotIndex{ii}(jj));
            viewSide(jj) = T.viewSide(shotIndex{ii}(jj));
            limitMode(jj) = T.limitMode(shotIndex{ii}(jj));
            pulseLength(jj) = T.pulseLength(shotIndex{ii}(jj));

            X(jj) = T.X(shotIndex{ii}(jj));
            
            seqFileName = ['Shot ',num2str( T.shot(shotIndex{ii}(jj)) ),'.seq'];
            pathName = [address{ii}{jj},'\'];

            % -------------------------------------------------------------
            % -------------------------------------------------------------
            tjj = tic;
            
            % Set parameters for the extraction of data from seq files:
            % ---------------------------------------------------------------------
            % Set thermal parameters:
            thermalParam.ExternalOpticsTransmission = 0.7;
            thermalParam.AtmosphericTemperature = 24;
            thermalParam.Distance = 1;
            thermalParam.ExternalOpticsTemperature = 24;
            thermalParam.ReferenceTemperature = 24;
            thermalParam.Transmission = 1;
            thermalParam.RelativeHumidity = 0;
            thermalParam.ReflectedTemperature = 24;
            thermalParam.emissivity = 0.55;
            % Define extraction options:
            extractOptions.frames = 100:1:250;
            if ii >= 7
                extractOptions.frames = 100:1:300;
            end
            extractOptions.frameRate = 100;
            
            % Get the intensity data:
            disp(['Extracting data from ',seqFileName,' ...'])
            [f{jj},t_f{jj},seq{jj}] = ExtractDataFromSeqFile(seqFileName,pathName,thermalParam,extractOptions);
            
            % Extract data during plasma pulse and compute integrated intensity:
            options.n_Before = 10;
            options.n_After = 50;
            [intf{jj},rngPlasma{jj},t0Plasma(jj)] = GetDataDuringPlasma(f{jj},t_f{jj},options);
                            
            % Convert intensity to temperature:
            dum1 = f{jj}(:,:,rngPlasma{jj});
            temperature{jj} = IntensityTempConv(thermalParam.emissivity,dum1,seq{jj});
            t_temperature{jj} = t_f{jj}(rngPlasma{jj});
                        
            tjj = toc(tjj);
            % -------------------------------------------------------------
            % -------------------------------------------------------------

            disp(['Shot ',num2str(shot(jj)),' loaded from .seq file in ',num2str(tjj),' sec']);
        end
        
        tii = toc(tii);
        disp(['Total time to load .seq files from dataset: ',num2str(tii),' sec']);
        
        % Plot the integrated intensity (intf):
        figure('color','w')
        grid on
        hold on
        lineColor = {'k','bl','r','g','m','c','k','bl','r','g','m','c'};
        for jj = 1:numberOfShotsInDataset(ii)
                        plot(t_f{jj}                - t0Plasma(jj),intf{jj}                - min(intf{jj}),lineColor{jj})
                h(jj) = plot(t_f{jj}(rngPlasma{jj}) - t0Plasma(jj),intf{jj}(rngPlasma{jj}) - min(intf{jj}),lineColor{jj},'LineWidth',3);
                legendText{jj} = num2str(shot(jj));
        end
        set(gca,'FontName','Times','FontSize',11)
        title(['Dataset: ',num2str(ii)],'interpreter','Latex','FontSize',14)
        xlabel('time [s]','interpreter','Latex','FontSize',14)
        ylabel('$\int{f(t,x)}dx^2$','interpreter','Latex','FontSize',14)
        legend(h,legendText)
        xlim([-0.2,1])
        box on
        figureName = ['IntegratedIntensity_Dataset_',num2str(ii)];
        if saveFigure
            saveas(gcf,figureName,'tiffn')
        end
        
        % Saving the data for the iith dataset:
        t1 = tic;
        disp(['Saving data into .mat file for dataset ',num2str(ii),' ...'])
        save(matFileName{ii},varList{:},'-v7.3')
        t1 = toc(t1);

        disp(['Data saved for dataset ',num2str(ii),' ... took ',num2str(t1),' sec'])
    end

end

t0 = toc(t0);
disp(['Total calculation time: ',num2str(t0),' sec'])

% Notes:
% Total calculation time: 225.0814 sec
% It takes 89 sec to extract all data from .mat files

return
%% Test datasets
close all

% The function matfile allows you to access data inside a .mat file without
% loading the entire .mat file. 
% For example, if only one shot is required from a dataset, we can use this
% functionality to extract only the shot that is required.

ii = 8;
jj = 1;
d = matfile(['dataset_',num2str(ii),'_IRdata.mat']);
A = cell2mat(d.temperature(:,jj));

options.frames = 1:3:size(A,3)-10;
options.colorbar = 1;
options.magnitudePlotMode = 1;
options.removeAxesTicks = 0;
options.shot = d.shot(:,jj);
options.zlim = [0,12];
options.mirrorImage = 1;

PlayMovieFromArray(A,options)


ii = 5;
jj = 6;
d = matfile(['dataset_',num2str(ii),'_IRdata.mat']);
A = cell2mat(d.temperature(:,jj));

options.frames = 1:3:size(A,3)-10;
options.colorbar = 1;
options.magnitudePlotMode = 1;
options.removeAxesTicks = 0;
options.shot = d.shot(:,jj);
options.zlim = [0,12];
options.mirrorImage = 1;

PlayMovieFromArray(A,options)

%%
ii = 9;
jj = 2;
d = matfile(['dataset_',num2str(ii),'_IRdata.mat']);
A = cell2mat(d.temperature(:,jj));

options.frames = 1:3:size(A,3)-10;
options.colorbar = 1;
options.magnitudePlotMode = 1;
options.removeAxesTicks = 0;
options.shot = d.shot(:,jj);
options.zlim = [0,12];
options.mirrorImage = 1;

PlayMovieFromArray(A,options)


ii = 4;
jj = 1;
d = matfile(['dataset_',num2str(ii),'_IRdata.mat']);
A = cell2mat(d.temperature(:,jj));

options.frames = 1:3:size(A,3)-10;
options.colorbar = 1;
options.magnitudePlotMode = 1;
options.removeAxesTicks = 0;
options.shot = d.shot(:,jj);
options.zlim = [0,12];
options.mirrorImage = 1;

PlayMovieFromArray(A,options)

%%
close all

ii = 7;
jj = 4;
d = matfile(['dataset_',num2str(ii),'_IRdata.mat']);
A = cell2mat(d.temperature(:,jj));

options.frames = 1:3:size(A,3)-10;
options.colorbar = 1;
options.magnitudePlotMode = 1;
options.removeAxesTicks = 0;
options.shot = d.shot(:,jj);
options.zlim = [0,12];
options.mirrorImage = 1;

PlayMovieFromArray(A,options)

ii = 8;
jj = 1;
d = matfile(['dataset_',num2str(ii),'_IRdata.mat']);
A = cell2mat(d.temperature(:,jj));

options.frames = 1:3:size(A,3)-10;
options.colorbar = 1;
options.magnitudePlotMode = 1;
options.removeAxesTicks = 0;
options.shot = d.shot(:,jj);
options.zlim = [0,12];
options.mirrorImage = 1;

PlayMovieFromArray(A,options)

ii = 9;
jj = 2;
d = matfile(['dataset_',num2str(ii),'_IRdata.mat']);
A = cell2mat(d.temperature(:,jj));

options.frames = 1:3:size(A,3)-10;
options.colorbar = 1;
options.magnitudePlotMode = 1;
options.removeAxesTicks = 0;
options.shot = d.shot(:,jj);
options.zlim = [0,12];
options.mirrorImage = 1;

PlayMovieFromArray(A,options)

%%
close all

ii = 5;
jj = 7;
d = matfile(['dataset_',num2str(ii),'_IRdata.mat']);
A = cell2mat(d.temperature(:,jj));

options.frames = 1:3:size(A,3)-10;
options.colorbar = 1;
options.magnitudePlotMode = 1;
options.removeAxesTicks = 0;
options.shot = d.shot(:,jj);
options.zlim = [0,12];
options.mirrorImage = 1;

PlayMovieFromArray(A,options)

ii = 6;
jj = 5;
d = matfile(['dataset_',num2str(ii),'_IRdata.mat']);
A = cell2mat(d.temperature(:,jj));

options.frames = 1:3:size(A,3)-10;
options.colorbar = 1;
options.magnitudePlotMode = 1;
options.removeAxesTicks = 0;
options.shot = d.shot(:,jj);
options.zlim = [0,12];
options.mirrorImage = 1;

PlayMovieFromArray(A,options)

ii = 4;
jj = 3;
d = matfile(['dataset_',num2str(ii),'_IRdata.mat']);
A = cell2mat(d.temperature(:,jj));

options.frames = 1:3:size(A,3)-10;
options.colorbar = 1;
options.magnitudePlotMode = 1;
options.removeAxesTicks = 0;
options.shot = d.shot(:,jj);
options.zlim = [0,12];
options.mirrorImage = 1;

PlayMovieFromArray(A,options)


ii = 1;
jj = 4;
d = matfile(['dataset_',num2str(ii),'_IRdata.mat']);
A = cell2mat(d.temperature(:,jj));

options.frames = 1:3:size(A,3)-10;
options.colorbar = 1;
options.magnitudePlotMode = 1;
options.removeAxesTicks = 0;
options.shot = d.shot(:,jj);
options.zlim = [0,12];
options.mirrorImage = 1;

PlayMovieFromArray(A,options)

% Notes:
% - use matfile() to extract individual shots from .mat files without
% needing to load the entire file into memory.
% - Before any analysis can be done, we need to determine the start of the RF
% for every shot and produce a new time stamp relative to the start of the
% RF

