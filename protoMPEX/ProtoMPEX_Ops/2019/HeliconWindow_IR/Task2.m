% Task 2: Define start of RF pulse and image masks
clc
close all
disp('Start of task 2 ##################################################');

%% SECTION 1
% =========================================================================
% Compute prompt:
InputStructure.prompt = {['Would you like to reload raw data from .mat files? Yes [1], No [0]']};
beep
cmpt = GetUserInput(InputStructure);

% =========================================================================
% Clear memory
if cmpt
    clearvars('-except','cmpt');
end

% #########################################################################
% SECTION 1:
sectionName = 'Load .mat file';
disp(['Section: ',sectionName])
% #########################################################################

% =========================================================================
% Load raw data from .mat files:
if cmpt
    % ---------------------------------------------------------------------
    % Input required from user:
    datasetToAnalyze = 11;
    
    % ---------------------------------------------------------------------
    % Load raw data:
    disp('Loading raw data from .mat file...')
    loadMatFile_ElapsedTime = tic;
    rawDatafileName= ['dataset_',num2str(datasetToAnalyze),'_raw.mat'];
    load(rawDatafileName)
    loadMatFile_ElapsedTime = toc(loadMatFile_ElapsedTime);
    disp(['Time required to load raw data: ',num2str(loadMatFile_ElapsedTime),' s'])
end

% =========================================================================
% Define variables specific to this section:
if cmpt
    vars_section1 = '';
    vars_section1 = who;
    vars_section1 = setdiff(vars_section1,[]);
end

% =========================================================================
disp('***********************Section 1 completed*************************')
%% SECTION 2
% =========================================================================
% Clear memory and figures:
try
    close(findobj('Tag','section2'))
    clearvars('-except',vars_section1{:});
end

% #########################################################################
% SECTION 2:
sectionName = 'Start of RF pulse';
disp(['Section: ',sectionName])
% #########################################################################

% =========================================================================
% Compute prompt:
InputStructure.prompt = {['Would you like to recompute the start of the RF? Yes [1], No [0]']};
beep
cmpt = GetUserInput(InputStructure);

% =========================================================================
% Define the start of the RF:
if cmpt
    disp('Computing...')
    
    % ---------------------------------------------------------------------
    % Read dataset spreadsheet:
    spreadsheetName = 'HeliconWindowIR_XPs.xlsx';
    datasetTable = readtable(spreadsheetName,'Sheet',2);
    
    % ---------------------------------------------------------------------
    % Identify shots from the dataset table:
    [shotIndex,~] = find(datasetTable.dataset == datasetToAnalyze);
    shots = datasetTable.shot(shotIndex);
    nshots = numel(shots);
    
    tic
    % Offsets to reject data from the entire time trace
    n1_offset = 10;
    n2_offset = 40;
    % Offsets to widen the RF pulse window
    n_Before = 10;
    n_After = 60;
    frameRate = 100;
    dt = 1/frameRate;
    
    % ---------------------------------------------------------------------
    % Define start and end of pulse and define data range
    for s = 1:nshots
        [Nx,Ny,Nz] = size(intensityRaw{s}); 

        for ii = 1:Nz
           intensityRaw_mean{s}(ii) = mean(mean(intensityRaw{s}(1:40:Nx,1:40:Ny,ii))); 
        end
        t_intensityRaw_mean{s} = 0:dt:(Nz-1)*dt;

        % Find the start and end of the RF pulse:
        [~,n1(s)] = max(diff(intensityRaw_mean{s}(n1_offset:end),1));
        n1(s) = n1(s) + n1_offset - 1;
        [~,n2(s)] = min(diff(intensityRaw_mean{s}(n2_offset:end),1));
        n2(s) = n2(s) + n2_offset - 1;

        % Start of the RF pulse relative to the raw data's time trace
        t0_intensityRaw_mean(s) = t_intensityRaw_mean{s}(n1(s));

        % Define the time window over which the RF is on:
        % Include some points before and after the RF
        rng{s} = [n1(s)-n_Before:n2(s)+n_After];
    end
    
    % ---------------------------------------------------------------------
    % Extract the relevent data from "intensityRaw":
    for s = 1:nshots
        intensity{s} = intensityRaw{s}(:,:,rng{s});
        t_intensity{s} = 0:dt:(length(rng{s})-1)*dt;
    end    
    toc 
else
    % ---------------------------------------------------------------------
    % Load existing data
    disp('Loading precalculated data from section 2...')
    fileName = ['dataset_',num2str(datasetToAnalyze),'_postprocess.mat'];
    load(fileName)
    disp('Data loaded succesfully')
end

% =========================================================================
% Plot data to confirm correctness of calculation
X = datasetTable.X(shotIndex);
scanType = datasetTable.scanType(shotIndex);

figure; 
set(gcf,'Tag','section2')
hold on
for s = 1:nshots
    intensityRaw_mean_noOffset{s} = intensityRaw_mean{s}-min(intensityRaw_mean{s});
    plot(t_intensityRaw_mean{s}-t0_intensityRaw_mean(s),intensityRaw_mean_noOffset{s}.^1,'LineWidth',0.5);
    hIR(s) = plot(t_intensityRaw_mean{s}(rng{s})-t0_intensityRaw_mean(s),intensityRaw_mean_noOffset{s}(rng{s}).^1,'LineWidth',3);
end
title('intensityRaw_mean_noOffset','Interpreter','none')
ylabel('{\Delta}intensity')
xlabel('t [s]')
xlim([-0.5,2])
legend(hIR,num2str(X))
box on; set(gcf,'color','w')

% =========================================================================
% Define variables specific to this section:
if cmpt
    vars_section2 = '';
    vars_section2 = who;
    vars_section2 = setdiff(vars_section2,vars_section1);
end

% =========================================================================
% Save prompt:
if cmpt
    InputStructure.prompt = {['Would you like to save? Yes [1], No [0]']};
    InputStructure.option.WindowStyle = 'normal';
    svdt = GetUserInput(InputStructure);
    if svdt        
        variableNames = vars_section2;
        fileName = ['dataset_',num2str(datasetToAnalyze),'_postprocess.mat'];
        SaveData
    end
end

% =========================================================================
disp('***********************Section 2 completed*************************')

%% SECTION 3
% =========================================================================
% Clear memory and figures:
try
    close(findobj('Tag','section3'))
    clearvars('-except',vars_section1{:},vars_section2{:});
end

% #########################################################################
% SECTION 3:
sectionName = 'Calculate gradients';
disp(['Section: ',sectionName])
% #########################################################################

% =========================================================================
% Compute prompt:
InputStructure.prompt = {['Would you like to recompute image Mask? Yes [1], No [0]']};
beep
cmpt = GetUserInput(InputStructure);

% =========================================================================
% Compute:
if cmpt
    disp('Computing...')
    tic
    % Extract the gradient
    for s = 1:nshots
        % intensity contains a subset from the initial Raw intensity
        % The first element of intensity is a few frames prior to the RF

        % Frame just prior to the RF:
        intensity0{s} = intensity{s}(:,:,1);
        % Frame a few frames after the RF:
        intensity1{s} = intensity{s}(:,:,end);

        % Gradient of the reference frame:
        [intensity1_gx,intensity1_gy] = gradient(intensity0{s});
        intensity1_gabs{s} = sqrt(intensity1_gy.^2 + intensity1_gx.^2);

        % Need to remove effect from edges from intensity1_gabs:
        intensity1_gabs{s}(1,:) = 0;
        intensity1_gabs{s}(end,:) = 0;
        intensity1_gabs{s}(:,1) = 0;
        intensity1_gabs{s}(:,end) = 0;
    end
    
    edges = logspace(-1,4,1e2);
    % Statistics on intensity0
    % Determine the most likely peaks in the intensity in the presence of spikes in
    % some pixels
    for s = 1:nshots
        intensity0_mean = mean(mean(intensity0{s}));
        intensity0_noOffset{s} = abs(intensity0{s} - intensity0_mean);
        % Produce statistics of pixel intensities:
        [intensity0_hist{s},~] = histcounts(intensity0_noOffset{s},edges,'normalization','cdf');
        % Determine the fraction of outliers to exclude
        minPixCountFrac(s) = 1500/numel(intensity0{s});
        % Based on this fraction, determine the most probably peak (mpp)
        intensity0_mpp(s) = edges(find((1-intensity0_hist{s})<minPixCountFrac(s),1));
    end
    
   % Statistics on "intensity1_gabs"
    for s = 1:nshots
        % Perform a smoothing operation of data
        intensity1_gabs_smooth{s} = Blur(intensity1_gabs{s},1)/3;
        [intensity1_hist{s},~] = histcounts(abs(intensity1_gabs_smooth{s}),edges,'normalization','cdf');
        % Define the most propable gradient (mpg)
        intensity1_mpg(s) = edges(find((1-intensity1_hist{s})<minPixCountFrac(s),1));
        % Find locations of most probable gradients
        [nx_mpg{s},ny_mpg{s}] = find(abs(intensity1_gabs_smooth{s})>intensity1_mpg(s));
    end  
    toc
end
 
% =========================================================================
% Get input from user
if cmpt
    disp('Selecting edges...')
    
    instancesType = 1;
    switch instancesType
        case 1
            instances = 1:nshots;
        case 2
            % Use this case if one needs to correct a few shots and does
            % not wish to recalculate all shots again
            fileName = ['dataset_',num2str(datasetToAnalyze),'_postprocess.mat'];
            load(fileName)
            instances = 4;
    end
    
    figure
    for s = instances
        exitFlag = 0;
        while exitFlag == 0
            clf
            disp(['Shot #',num2str(s),' of ',num2str(nshots)])
            subplot(2,2,[1,2])
            surf(abs(intensity1_gabs_smooth{s}),'LineStyle','none')
            set(gcf,'Tag','section3','color','w')
            title(['shot: ',num2str(shots(s)),', #',num2str(s),' of ',num2str(nshots)])

            set(gca,'XTick',[],'YTick',[])
            axis('image')
            view([0,90])
            colormap(flipud(hot))
            caxis([0,intensity1_mpg(s)])
            zlim([0,intensity1_mpg(s)*2])

            hold on
            for ii = 1:length(ny_mpg{s})-1
                mpg{s}(ii) = intensity1_gabs_smooth{s}(nx_mpg{s}(ii),ny_mpg{s}(ii));
                plot3(ny_mpg{s}(ii),nx_mpg{s}(ii),abs(mpg{s}(ii)),'g.')
            end
            hold off
        
           % --------------------------------------------------------------
           % Get user to select edges:
           disp('Waiting for user input...')
           [nx_mpg_ginput(s),ny_mpg_ginput(s)] = ginput(1);
           nx_mpg_ginput(s) = floor(nx_mpg_ginput(s));
           ny_mpg_ginput(s) = floor(ny_mpg_ginput(s));
           
           mpg_ginput(s) = intensity0_noOffset{s}(ny_mpg_ginput(s),nx_mpg_ginput(s));
           windowArea{s} = double(intensity0_noOffset{s}>mpg_ginput(s));

           % --------------------------------------------------------------
           % Show effect of selection:
           subplot(2,2,3)
           surf(intensity0_noOffset{s},'LineStyle','none')
           view([0,90])
           axis('image')
           set(gca,'XTick',[],'YTick',[])
           title(['shot: ',num2str(shots(s)),', Initial'])
           subplot(2,2,4)
           surf(intensity0_noOffset{s}.*windowArea{s},'LineStyle','none')
           view([0,90])
           axis('image')
           set(gca,'XTick',[],'YTick',[])
           title(['shot: ',num2str(shots(s)),', Masked'])

           % --------------------------------------------------------------
           % Query user to select or reject section:
           disp('Waiting for user to accept or reject selection...')
            InputStructure.prompt = {['Keep selection? Yes [1], No [0]']};
            exitFlag = GetUserInput(InputStructure);         
            switch exitFlag
                case 0
                    disp('Repeat selection')
                case 1
                    disp('Selection accepted')                   
            end % switch           
        end % while
    end % for
    
    % =====================================================================
    % Define variables specific to this section:
    vars_section3 = '';
    vars_section3 = who;
    vars_section3 = setdiff(vars_section3,[vars_section2;vars_section1]);
    
    % =====================================================================
    % Save prompt:
    disp('Save prompt...')
    InputStructure.prompt = {['Would you like to save? Yes [1], No [0]']};
    InputStructure.option.WindowStyle = 'normal';
    beep
    svdt = GetUserInput(InputStructure);
    if svdt        
        variableNames = vars_section3;
        fileName = ['dataset_',num2str(datasetToAnalyze),'_postprocess.mat'];
        SaveData
    else
        disp('User choose not to save data')
    end
else
    % ---------------------------------------------------------------------
    % Load existing data
    disp('Loading existing data...')
    fileName = ['dataset_',num2str(datasetToAnalyze),'_postprocess.mat'];
    load(fileName)
    disp('Data loaded succesfully')
end

% =========================================================================
% Plot data:
disp('Plotting main results from section 3...')
switch nshots
    case 1
        nCol = 1;
        nRow = 1;
    case 2
        nCol = 2;
        nRow = 1;
    case {3,4}
        nCol = 2;
        nRow = 2;
    case {5,6}
        nCol = 3;
        nRow = 2;
    case {7,8,9}
        nCol = 3;
        nRow = 3;
    case {10,11,12}
        nCol = 4;
        nRow = 3;
    case {13,14,15,16}
        nCol = 4;
        nRow = 4;
end

figure
set(gcf,'Tag','section3')
for s = 1:nshots
    subplot(nRow,nCol,s)
    surf(intensity0_noOffset{s}.*windowArea{s},'LineStyle','none')
    colormap(flipud(hot))
    view([0,90])
    axis('image')
    set(gca,'XTick',[],'YTick',[])
    title(['shot: ',num2str(shots(s)),', #',num2str(s),' of ',num2str(nshots)],'FontSize',7)
end

% =========================================================================
disp('***********************Section 3 completed*************************')

%% SECTION 4
% =========================================================================
% Clear memory and figures:
try
    close(findobj('Tag','section4'))
    clearvars('-except',[vars_section1{:};vars_section2{:};vars_section3{:}]);
end

% #########################################################################
% SECTION 4:
sectionName = 'Convert intensity to temperature';
disp(['Section: ',sectionName])
% #########################################################################

% =========================================================================
% Compute prompt:
InputStructure.prompt = {['Would you like to recalculate temperature? Yes [1], No [0]']};
beep
cmpt = GetUserInput(InputStructure);
if cmpt
    tic
    disp('Loading .seq files')
    % ---------------------------------------------------------------------
    % Load .seq files
    currentFolder = pwd;
    dummy1 = findstr(currentFolder,'\');
    rootAddress = currentFolder(1:dummy1(end));
    
    % ---------------------------------------------------------------------
    % Load the Atlats SDK
    atPath = getenv('FLIR_Atlas_MATLAB');
    atImage = strcat(atPath,'Flir.Atlas.Image.dll');
    asmInfo = NET.addAssembly(atImage);
    
    % ---------------------------------------------------------------------
    % For each shot, get the address where the .seq file is located
    for s = 1:nshots
        dateOfShot = datasetTable.date{shotIndex(s)};
        address{s} = [rootAddress,dateOfShot(1:end-3),'\',dateOfShot,'\IR_RawData'];
        
        PATHNAME = [address{s},'\'];
        seqFileName = ['Shot ',num2str( shots(s) ),'.seq'];
        videoFileName=[PATHNAME seqFileName];
        file = Flir.Atlas.Image.ThermalImageFile(videoFileName);
        
        % Define the seq file and associated thermal parameters
        seq{s} = file.ThermalSequencePlayer();
        seq{s}.ThermalImage.ThermalParameters.ExternalOpticsTransmission = 0.7;
        seq{s}.ThermalImage.ThermalParameters.AtmosphericTemperature = 24;
        seq{s}.ThermalImage.ThermalParameters.Distance = 1;
        seq{s}.ThermalImage.ThermalParameters.ExternalOpticsTemperature = 24;
        seq{s}.ThermalImage.ThermalParameters.ReferenceTemperature = 24;
        seq{s}.ThermalImage.ThermalParameters.Transmission = 1;
        seq{s}.ThermalImage.ThermalParameters.RelativeHumidity = 0;
        seq{s}.ThermalImage.ThermalParameters.ReflectedTemperature = 24;
    end
    disp('Succesfully loaded .seq files')
    
    % ---------------------------------------------------------------------
    % Calculate surface temperature
    disp('Calculating surface temperature')
    emissivity = 0.97;
    for s = 1:nshots
        [temperature{s}] = IntensityTempConv(emissivity,intensity{s},seq{s});
        for fr = 1:size(temperature{s},3)
                deltaT{s}(:,:,fr) = (temperature{s}(:,:,fr)-temperature{s}(:,:,1)).*windowArea{s};
                deltaT_peak{s}(fr) = mean(mean(deltaT{s}(1:40:Nx,1:40:Ny,fr)));
        end
        [~,index_deltaT_max(s)] = max(deltaT_peak{s}); 
        t_deltaT{s} = t_intensity{s};
    end
    
    % ---------------------------------------------------------------------
    % clear large variables to reduce saving and loading times:
    clearvars temperature
    toc
 end % cmpt

% =========================================================================
% Select temperature probes
for s = 1:nshots
    figure('Tag','section4')  
    
    % -----------------------------------------------------------------
    if cmpt
        subplot(2,1,1)
        for fr = 1:2:size(deltaT{s},3)
            surf(deltaT{s}(:,:,fr),'LineStyle','none')
            view([0,90])
            caxis([0,16])
            set(gca,'XTick',[],'YTick',[])
            axis('image')
            colormap(flipud(hot))
            title(['frame: ',num2str(fr),', time: ',num2str(t_deltaT{s}(fr)),' sec'])
            drawnow
        end
    end
    
    % ---------------------------------------------------------------------
    exitFlag = 0;
    while exitFlag == 0 
        subplot(2,1,1)
        deltaT_max = deltaT{s}(:,:,index_deltaT_max(s));
        surf(deltaT_max,'LineStyle','none')
        view([0,90])
        caxis([0,16])
        set(gca,'XTick',[],'YTick',[])
        axis('image')
        colormap(flipud(hot))
        title(['Hottest frame, shot: ',num2str(shots(s)),', #',num2str(s),' of ',num2str(nshots)])
        colorbar
        pause(3)
        
        % -----------------------------------------------------------------
        if cmpt
            % Get input from user:
            disp('Waiting for user input...')
            [yCenter(s),xCenter(s)] = ginput(1);
            xCenter(s) = round(xCenter(s));
            yCenter(s) = round(yCenter(s));
        else
            % Load precalculate data:
            fileName = ['dataset_',num2str(datasetToAnalyze),'_postprocess.mat'];
            load(fileName,'xCenter','yCenter')
            exitFlag = 1;
        end % cmpt
        
        % -------------------------------------------------------------
        % Plot location
        rngy_tempProbe = [(yCenter(s)-10):(yCenter(s)+10)];
        rngx_tempProbe = [(xCenter(s)-10):(xCenter(s)+10)];
        [x_tempProbe,y_tempProbe] = meshgrid(rngy_tempProbe,rngx_tempProbe);
        z_tempProbe = 20*ones(size(y_tempProbe));
        hold on
        plot3(x_tempProbe(:,end)  ,y_tempProbe(:,end)  ,z_tempProbe(1,:),'k.','MarkerSize',0.5)
        plot3(x_tempProbe(1,:  )  ,y_tempProbe(1,:  )  ,z_tempProbe(1,:),'k.','MarkerSize',0.5)
        plot3(x_tempProbe(:,1  )  ,y_tempProbe(:,1  )  ,z_tempProbe(1,:),'k.','MarkerSize',0.5)
        plot3(x_tempProbe(end,:  ),y_tempProbe(end,:  ),z_tempProbe(1,:),'k.','MarkerSize',0.5)
        hold off
        drawnow
        
        % -------------------------------------------------------------
        % Plot temperature probe:
        for fr = 1:size(deltaT{s},3)
            tempProbe_mean{s}(fr) = mean(mean(deltaT{s}(rngx_tempProbe,rngy_tempProbe,fr)));
            t_tempProbe_mean{s} = t_deltaT{s};
        end
        subplot(2,1,2)
        plot(t_tempProbe_mean{s},tempProbe_mean{s},'LineWidth',2)
        box on
        grid on
        ylim([0,20])
        xlim([0,1])
        ylabel('$\Delta T$ [C]','Interpreter','latex','FontSize',12)
        xlabel('time [s]','Interpreter','latex','FontSize',12)
        set(gca,'FontName','times','FontSize',11)
        
        % -------------------------------------------------------------
        if cmpt          
            % Query user to accept or reject selection:
            disp('Waiting for user to accept or reject selection...')
            InputStructure.prompt = {['Keep selection? Yes [1], No [0]']};
            exitFlag = GetUserInput(InputStructure);         
            switch exitFlag
                case 0
                    disp('Repeat selection')
                case 1
                    disp('Selection accepted')                   
            end % switch    
        end % cmpt 
        
    end % while  
    
end % for 

% =========================================================================
% Review data and save:
if cmpt   
    % ---------------------------------------------------------------------
    % Review selections
    hDummy = findobj('Tag','section4');
    try
        for s = 1:length(hDummy)
            figure(hDummy(s))
            pause(1)
        end
    end
    clearvars hDummy
    
    % ---------------------------------------------------------------------
    % Define variables specific to this section:
    vars_section4 = '';
    vars_section4 = who;
    vars_section4 = setdiff(vars_section4,[vars_section3;vars_section2;vars_section1]);
    
    % ---------------------------------------------------------------------
    % Save data
    disp('Save prompt...')
    InputStructure.prompt = {['Would you like to save? Yes [1], No [0]']};
    InputStructure.option.WindowStyle = 'normal';
    beep
    svdt = GetUserInput(InputStructure);
    if svdt        
        variableNames = vars_section4;
        fileName = ['dataset_',num2str(datasetToAnalyze),'_postprocess.mat'];
        SaveData
    else
        disp('Data for section 4 not saved')
    end % if
end % cmpt

% =========================================================================
% Save figures from section4
disp('Save prompt...')
InputStructure.prompt = {['Would you like to save figures from section 4? Yes [1], No [0]']};
InputStructure.option.WindowStyle = 'normal';
beep
svfig = GetUserInput(InputStructure);
if svfig
    folderName = ['dataset_',num2str(datasetToAnalyze),'_figures'];
    if exist(folderName) == 0
        mkdir(folderName);
    end
    h_dummy = findobj('Tag','section4');
    address_home = cd;
    for ii = 1:length(h_dummy)
        figure(h_dummy(ii))
        figureName = ['ds_',num2str(datasetToAnalyze),'_section4_shot_',num2str(shots(ii))];
        cd([address_home,'\',folderName])
        saveas(gcf,figureName,'tiffn')
        cd(address_home)
    end
end

% =========================================================================
disp('***********************Section 4 completed*************************')

% =========================================================================
% End of script
disp('End of script!!')

return

%% COMMENTS
% need to convert deltaT into heat flux
% The final task is to plot the various data sets and combinations of

% To mirror flip the images use
% surf(640:-1:1,1:240,deltaT{1}(:,:,80),'LineStyle','none'); view([0,90]); axis('image')
% or just flip the image data
% figure; surf(deltaT{1}(:,640:-1:1,80),'LineStyle','none'); view([0,90]); axis('image')

% #########################################################################
% Notes:
% figures can be saved into a .mat file by their handle
% figure('Tag','figureTag'); plot(x,sin(x));
% figureObjectHandle = findobj('Tag','figureTag')
% load('dataset_12_postprocess.mat','figureObjectHandle')
% This object includes the entire figure with the associated lines and axes
% axesHandle = figureObjectHandle.Children
% lineHandle = figureObjectHandle.Children.Children

% In this way, one can create a figure then create a figure handle to 
% encapsulate all the data in it.
% If one saves the figure to a .mat file, one can "replot" the data just by
% calling the figure handle as follows:
% figure(figureObjectHandle)

% More information for this can be found under "Graphics Objects"


