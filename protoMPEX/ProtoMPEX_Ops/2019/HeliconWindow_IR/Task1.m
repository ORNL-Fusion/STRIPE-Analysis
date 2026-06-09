% Task 1: Load .seq files

clc
clear all
close all

totalElapsedTime = tic;
disp('Start of task 1 ##################################################');

% Read dataset spreadsheet
spreadsheetName = 'HeliconWindowIR_XPs.xlsx';
T = readtable(spreadsheetName,'Sheet',2);
datasets = unique(T.dataset);
numberOfDatasets = numel(datasets);

currentFolder = pwd;
dummy1 = findstr(currentFolder,'\');
rootAddress = currentFolder(1:dummy1(end));

% ii is index for number of shots
% jj is index for shot within a dataset

cmpt = 0;
if cmpt
for ii = 1:numberOfDatasets
    % Define the name of the .mat file
    rawDatafileName{ii}= ['dataset_',num2str(datasets(ii)),'_raw.mat'];
    % For each data set, assign an index to each of the shots
    [shotIndex{ii},~] = find(T.dataset == datasets(ii));
    numberOfShotsInDataset(ii) = numel(shotIndex{ii});
    
    % For each shot, get the address where the .seq file is located
    for jj = 1:numberOfShotsInDataset(ii)
        dateOfShot = T.date{shotIndex{ii}(jj)};
        address{ii}{jj} = [rootAddress,dateOfShot(1:end-3),'\',dateOfShot,'\IR_RawData'];
    end
end

% Load the Atlats SDK
atPath = getenv('FLIR_Atlas_MATLAB');
atImage = strcat(atPath,'Flir.Atlas.Image.dll');
asmInfo = NET.addAssembly(atImage);

for ii = 1:numberOfDatasets
    disp(['Dataset ',num2str(ii),' -------------------------------------'])
    loadElapsedTime_dataset = tic;
    for jj = 1:numberOfShotsInDataset(ii)
        loadElapsedTime_singleShot = tic;
        PATHNAME = [address{ii}{jj},'\'];
        seqFileName = ['Shot ',num2str( T.shot(shotIndex{ii}(jj)) ),'.seq'];
        videoFileName=[PATHNAME seqFileName];
        file = Flir.Atlas.Image.ThermalImageFile(videoFileName);
        
        % Define the seq file and associated thermal parameters
        seq = file.ThermalSequencePlayer();
        seq.ThermalImage.ThermalParameters.ExternalOpticsTransmission = 0.7;
        seq.ThermalImage.ThermalParameters.AtmosphericTemperature = 24;
        seq.ThermalImage.ThermalParameters.Distance = 1;
        seq.ThermalImage.ThermalParameters.ExternalOpticsTemperature = 24;
        seq.ThermalImage.ThermalParameters.ReferenceTemperature = 24;
        seq.ThermalImage.ThermalParameters.Transmission = 1;
        seq.ThermalImage.ThermalParameters.RelativeHumidity = 0;
        seq.ThermalImage.ThermalParameters.ReflectedTemperature = 24;

        %Get the pixels
        img = seq.ThermalImage.ImageProcessing.GetPixelsArray;
        im = double(img);

        intensityRaw{jj}(:,:,1) = im;
        shotNumber(jj) = T.shot(shotIndex{ii}(jj));
        fr = 1;
        if(seq.Count > 1)
            while(seq.Next())
                img = seq.ThermalImage.ImageProcessing.GetPixelsArray;
                im = double(img);
                intensityRaw{jj}(:,:,fr) = im(end:-1:1,:);         
                fr = fr + 1;
            end
        end
        loadElapsedTime_singleShot = toc(loadElapsedTime_singleShot);
        disp(['Shot ',num2str(shotNumber(jj)),' loaded from .seq file in ',num2str(loadElapsedTime_singleShot),' sec']);
    end
    loadElapsedTime_dataset = toc(loadElapsedTime_dataset);
    disp(['Total time to load .seq files from dataset: ',num2str(loadElapsedTime_dataset),' sec']);
    
    saveElapsedtime = tic;
    disp(['Saving data into .mat file for dataset ',num2str(ii),' ...'])
    save(rawDatafileName{ii},'intensityRaw','shotNumber','loadElapsedTime_dataset','-v7.3')
    clear intensityRaw im shotNumber
    saveElapsedtime = toc(saveElapsedtime);
    save(rawDatafileName{ii},'saveElapsedtime','-append')

    disp(['Data saved for dataset ',num2str(ii),' ... took ',num2str(saveElapsedtime),' sec'])
end
end

totalElapsedTime = toc(totalElapsedTime);
disp(['Total calculation time: ',num2str(totalElapsedTime),' sec'])

% Notes:
% The entire process took 1600 seconds. the largest time sink is the saving
% process with 1200 sec and producing 9.3 GB in total
% Use the following to read and plot the table on the command line
% T = readtable('HeliconWindowIR_XPs.xlsx','Sheet',2)
% here T is a structure where each field corresponds to a column from the
% spreadsheet
% It appears that loading .mat files of raw data takes 66% of the time is
% requires to load them direclty from the .seq files
% dataset 1: 101 sec compared to 55 sec, which leads to 7.2 sec per .seq
% load and to 3.92 sec per .mat load
% dataset 2:  85 sec compared to 54 sec

% How do we prove that we have saved the correct files?

