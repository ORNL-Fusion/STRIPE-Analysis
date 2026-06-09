% Preprocessor 1
% This script was created on 2019_11_22
% This script loads the IR data from experiments performed on the
% 2019_06_19 and 2019_06_27
% This script deals only with the raw data from the IR camera which we
% refer to as intensityRaw and intensity

% The objective of preprocessor 1 is to:
% 1 - Identify the start of the RF pulse
% 2 - Identify the window area from the intensity data
% 3 - Extract intensity data during RF pulse
% 4 - Save the data needed for post-processing

clc % Clear command window
clear all
close all
 
%%
sectionName = 'Load IR video files';
disp(sectionName)

% Input required from user: -----------------------------------------------
xp = 1;
% -------------------------------------------------------------------------

% Compute prompt ----------------------------------------------------------
InputStructure.prompt = {['",',sectionName,'"? Yes [1], No [0]']};
cmpt = GetUserInput(InputStructure);
% -------------------------------------------------------------------------

if cmpt
    switch xp
        case 1                
            shot       = 26000 + [658,659,660];
            rfPwrNet   =         [104,122,135];
            addressLoc =         [1  ,1  ,1  ];
            viewType = 'Window limit, bottom';

        case 2        
            shot       = 26000 + [665,666,667,668,669,670,674,735,737,738,739,741,742];
            rfPwrNet   =         [169,161,154,135,114,96 ,144,68 ,94 ,87 ,78 ,58 ,47 ];
            addressLoc =         [1  ,1  ,1  ,1  ,1  ,1  ,1  ,2  ,2  ,2  ,2  ,2  ,2  ];
            viewType = 'Window limit, middle, ';

        case 3
            shot       = 26000 + [646,647,648,649,650,653,655,656,743,745,746,747,748,750];
            rfPwrNet   =         [56 ,89 ,138,153,150,134,165,104,95 ,88 ,78 ,68 ,58 ,48 ];
            addressLoc =         [1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,2  ,2  ,2  ,2  ,2  ,2  ];
            viewType = 'Window limit, top';
    end
    
    % Define initial variables:
    initialVars = {'xp','shot','rfPwrNet','addressLoc','viewType'};
    initialVars{end+1} = 'initialVars';
   
    % Define address where data is to be found
    RootAddress = 'C:\Users\nfc\Documents\ProtoMPEX_Ops\2019\2019_06';
    a{1} = [RootAddress,'\2019_06_19','\IR_RawData'];
    a{2} = [RootAddress,'\2019_06_27','\IR_RawData'];
    
    % Define file where data is to be stored
%     fileName = ['preprocessData_xp_',num2str(xp)];  
end

% Load IR files and create intensityRaw -----------------------------------
if cmpt 
    Task_1
    % Save prompt:
    clearvars InputStructure
    InputStructure.prompt = {'Save data? Yes [1], No [0]'};
    InputStructure.defaultAnswer = {'0'};
    svdata = GetUserInput(InputStructure);
end

% Save data ---------------------------------------------------------------
if svdata
    variableNames = initialVars;
    fileName = ['preprocessData_xp_',num2str(xp)];  
    SaveData
else
     tic
     disp('Loading previously saved data...');
     load(fileName);
     toc
end
% -------------------------------------------------------------------------

disp(viewType)

return
%%
% =========================================================================
% 1 - Determine the start of the RF pulse:
% =========================================================================
% To determine the start and end of the RF pulse we use the following
% method:
% For each frame, we calculate the mean pixel intensity for the entire image
% This operation results in a 1D vector called "intensityRaw_mean"
% We use a diff operation to find the starting and ending edges of
% "intensityRaw_mean", these edges are called n1 and n2
% We apply an offset to n1 and n2 to widen the time window
% From these edges we define a range which we use to extract a subset of
% data from the main raw data

% Section name ------------------------------------------------------------
sectionName = 'CalculateStartOfRf';
disp(sectionName)
try; close(findobj('Tag',sectionName)); end
% -------------------------------------------------------------------------
    

% Compute prompt ----------------------------------------------------------
prompt = {['Compute "',sectionName,'"? Yes [1], No [0]']};
% Add third option to kill program
name = '';
defaultAnswer = {'1'};
options.WindowStyle = 'normal';
answer = inputdlg(prompt,name,[1,length(prompt{:})+1],defaultAnswer,options);
if isempty(answer); 
    answer{1} = '0';
end
cmpt = str2num(answer{:});
% -------------------------------------------------------------------------

if cmpt
    disp('computing...')
    tic
    % Taken from AnalysisMethod_1
    % Offsets to reject data from the entire time trace
    n1_offset = 10;
    n2_offset = 40;
    % Offsets to widen the RF pulse window
    n_Before = 3;
    n_After = 60;
    frameRate = 100;
    dt = 1/frameRate;
    
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
    
    % Extract the relevent data from "intensityRaw":
    for s = 1:nshots
        intensity{s} = intensityRaw{s}(:,:,rng{s});
        t_intensity{s} = 0:dt:(length(rng{s})-1)*dt;
    end    
    toc 
    
     % Plot data to confirm correctness of calculation
    if 1
        figure; 
        set(gcf,'Tag',sectionName)
        hold on
        for s = 1:nshots
            intensityRaw_mean_noOffset{s} = intensityRaw_mean{s}-min(intensityRaw_mean{s});
            plot(t_intensityRaw_mean{s}-t0_intensityRaw_mean(s),intensityRaw_mean_noOffset{s}.^1,'LineWidth',0.5);
            hIR(s) = plot(t_intensityRaw_mean{s}(rng{s})-t0_intensityRaw_mean(s),intensityRaw_mean_noOffset{s}(rng{s}).^1,'LineWidth',2);
        end
        title('intensityRaw_mean_noOffset','Interpreter','none')
        ylabel('{\Delta}intensity')
        xlabel('t [s]')
        xlim([-0.5,2])
        legend(hIR,num2str(rfPwrNet'))
        box on; set(gcf,'color','w')
    end
end

% Save prompt ----------------------------------------------------------
fileName = ['preprocessData_xp_',num2str(xp)];    
variableNames = {'t_intensity','intensity'};

if cmpt    
    prompt = {'Save data? Yes [1], No [0]'};
    name = ''; defaultAnswer = {'0'}; options.WindowStyle = 'normal';
    answer = inputdlg(prompt,name,[1,length(prompt{:})+1],defaultAnswer,options);
    if isempty(answer); 
        answer{1} = '0';
    end
    saveData = str2num(answer{:});
    
     if saveData
        disp('Saving data...')
        tic
        for ii = 1:length(variableNames)
            try
                save(fileName,variableNames{ii},'-append');
            catch
                save(fileName,variableNames{ii})
            end
        end
        toc
        disp('Save complete!!'); msgbox('Save complete!!');
        whos ('-file',fileName')
     else
         disp('Not saving data')
     end
end
% -------------------------------------------------------------------------

%%
% =========================================================================
% 2 - Determine the window area
% =========================================================================

sectionName = 'DetermineWindowArea';

home
disp(sectionName)
try; close(findobj('Tag',sectionName)); end

tic

prompt = {['Compute "',sectionName,'"? Yes [1], No [0]']};
name = ''; defaultAnswer = {'0'}; options.WindowStyle = 'normal';

answer = inputdlg(prompt,name,[1,length(prompt{:})+1],defaultAnswer,options);
if isempty(answer); 
    answer{1} = '0';
end
cmptThresholds = str2num(answer{:});

return

% Extract the gradient
for s = 1:nshots
    % intensity contains a subset from the initial Raw intensity
    % The first element of intensity is a few frames prior to the RF
    
    % Frame just prior to the RF:
    intensity0{s} = intensity{s}(:,:,1);
    % Frame a few frames after the RF:
    intensity1{s} = intensity{s}(:,:,end);
    
    % Gradient of the reference frame:
    [intensity1_gx,intensity1_gy] = gradient(intensity1{s});
    intensity1_gabs{s} = sqrt(intensity1_gy.^2 + intensity1_gx.^2);
    
    % Need to remove effect from edges from intensity1_gabs:
    intensity1_gabs{s}(1,:) = 0;
    intensity1_gabs{s}(end,:) = 0;
    intensity1_gabs{s}(:,1) = 0;
    intensity1_gabs{s}(:,end) = 0;
end
toc

try
    hf = findobj('Tag','hist_intensity0');
    close(hf)
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

% Plot the reference frame with colormap based on "intensity0_mpp"
figure
for s = 1:nshots
    subplot(ceil(sqrt(nshots)),ceil(sqrt(nshots)),s)
    surf(intensity0_noOffset{s},'LineStyle','none'); set(gcf,'Tag','intensity0_hist'); caxis([0,3000]); view([0,90])
    axis('equal')
    set(gca,'XTick',[],'YTick',[])
    hTitle = title([num2str(shot(s)),' , RF: ',num2str(rfPwrNet(s)),' kW']);
    hTitle.FontSize = 8;
    clear hTitle
    caxis([0,intensity0_mpp(s)*1.2])
    colormap('hot')
        colormap(flipud(hot))
    xlim([0,size(intensity0_noOffset{s},2)])
    ylim([0,size(intensity0_noOffset{s},1)])
end

% Statistics on "intensity1_gabs"
for s = 1:nshots
    % Perform a smoothing operation of data
    intensity1_gabs_smooth{s} = Blur(intensity1_gabs{s},1)/3;
    [intensity1_hist{s},~] = histcounts(abs(intensity1_gabs_smooth{s}),edges,'normalization','cdf');
    % Define the most propable gradient (mpg)
    intensity1_mpg(s) = edges(find((1-intensity1_hist{s})<minPixCountFrac(s),1));
end

edgeSelectMode = 2;
switch edgeSelectMode
    case 1 % Only one instance
        instances = 1;
    case 2 % Define all instances
        instances = 1:nshots;
end

for s = instances
    figure; 
    set(gcf,'Tag','userSelectEdges')
    surf(abs(intensity1_gabs_smooth{s}),'LineStyle','none')
    axis('equal')
    view([0,90])
    colormap(flipud(hot))
    caxis([0,intensity1_mpg(s)])
    zlim([0,intensity1_mpg(s)*2])

    [nx_mpg{s},ny_mpg{s}] = find(abs(intensity1_gabs_smooth{s})>intensity1_mpg(s));
    hold on
    for ii = 1:length(ny_mpg{s})-1
        mpg{s}(ii) = intensity1_gabs_smooth{s}(nx_mpg{s}(ii),ny_mpg{s}(ii));
        plot3(ny_mpg{s}(ii),nx_mpg{s}(ii),abs(mpg{s}(ii)),'g.')
    end

    if cmptThresholds
          [nx_mpg_ginput(s),ny_mpg_ginput(s)] = ginput(1);
          nx_mpg_ginput(s) = floor(nx_mpg_ginput(s));
          ny_mpg_ginput(s) = floor(ny_mpg_ginput(s));
    end
    close(findobj('Tag','userSelectEdges'))
end
   
if cmptThresholds
    if edgeSelectMode == 1;
        nx_mpg_ginput = ones(nshots)*nx_mpg_ginput(1);
        ny_mpg_ginput = ones(nshots)*ny_mpg_ginput(1);
    end
    save(['preprocessData_xp_',num2str(xp)],'nx_mpg_ginput','ny_mpg_ginput');
else
    load(['preprocessData_xp_',num2str(xp)],'nx_mpg_ginput','ny_mpg_ginput');
end

for s = 1:nshots
    mpg_ginput(s) = intensity0_noOffset{s}(ny_mpg_ginput(s),nx_mpg_ginput(s));
    windowArea{s} = double(intensity0_noOffset{s}>mpg_ginput(s));
end

if cmptThresholds
   save(['preprocessData_xp_',num2str(xp)],'mpg_ginput','windowArea','-append');
   save(['preprocessData_xp_',num2str(xp)],'t_intensity','intensity','-append');
   save(['preprocessData_xp_',num2str(xp)],'shot','rfPwrNet','-append'); 
   save(['preprocessData_xp_',num2str(xp)],'viewType','-append'); 
   save(['preprocessData_xp_',num2str(xp)],'seq','-append'); 
end

% Plot the reference frame with colormap based on "mostProbablePeak"
figure
for s = 1:nshots
    subplot(ceil(sqrt(nshots)),ceil(sqrt(nshots)),s)
    surf(intensity0_noOffset{s}.*windowArea{s},'LineStyle','none'); set(gcf,'Tag','intensity0_hist');
%     zlim([intensity0_noOffset{s}(nx_mpg_ginput(s),nx_mpg_ginput(s)),5000])
    caxis([0,3000]); view([0,90])
    axis('equal')
    set(gca,'XTick',[],'YTick',[])
    hTitle = title([num2str(shot(s)),' , RF: ',num2str(rfPwrNet(s)),' kW']);
    hTitle.FontSize = 8;
    clear hTitle
    caxis([0,intensity0_mpp(s)*1.2])
    colormap(flipud(hot))
    xlim([0,size(intensity0_noOffset{s},2)])
    ylim([0,size(intensity0_noOffset{s},1)])
end

% =========================================================================
% Main outcome of this script is to produce the masks that we apply to the
% intensity to extract the window area
% =========================================================================

%% Window heating with RF

magConfig{1} = 'Window limit, bottom';
magConfig{2} = 'Window limit, middle';
magConfig{3} = 'Window limit, top';


try
    hf = findobj('Tag','intensityRawVid');
    close(hf)
end

s = 9;

switch xp
    case 1
        xCenter = 113;
        yCenter = 376;

        xCenter = 125;
        yCenter = 340;
    case 2
        xCenter = 64;
        yCenter = 367;
    case 3
        xCenter = 113;
        yCenter = 257;
        [~,s] = max(rfPwrNet);
end

rngy_intensity_probe = [(yCenter-10):(yCenter+10)];
rngx_intensity_probe = [(xCenter-10):(xCenter+10)];

figure;
set(gcf,'color','w')
set(gcf,'Tag','intensityRawVid')
for fr = 1:1:size(intensity{s},3)
    intensity_windowArea = intensity{s}(:,:,fr).*windowArea{s}  -  intensity{s}(:,:,1).*windowArea{s};
    surf(intensity_windowArea,'LineStyle','none')
    view([0,90])
    zlim([0,intensity0_mpp(s)*6])
    caxis([0,intensity0_mpp(s)*6])
    colormap(flipud(hot))
    title(['shot: ',num2str(shot(s)),', frame: ',num2str(fr),', P_{RF}: ',num2str(rfPwrNet(s)),' kW, ',magConfig{xp}])
    axis('equal')
    xlim([0,size(intensity0_noOffset{s},2)])
    ylim([0,size(intensity0_noOffset{s},1)])
    hold on
    plot3(rngy_intensity_probe,rngx_intensity_probe,intensity_windowArea(rngx_intensity_probe,rngy_intensity_probe),'k.')
    plot3(flip(rngy_intensity_probe),rngx_intensity_probe,intensity_windowArea(rngx_intensity_probe,rngy_intensity_probe),'k.')
    hold off
    drawnow
end

for s = 1:nshots
    for fr = 1:1:size(intensity{s},3)
        intensity_windowArea = intensity{s}(:,:,fr).*windowArea{s}  -  intensity{s}(:,:,1).*windowArea{s};
        intensity_probe{s}(fr) = mean(mean(intensity_windowArea(rngx_intensity_probe,rngy_intensity_probe)));
    end
end

figure
set(gcf,'Tag','intensityRawVid')
hold on
[~,b] = sort(rfPwrNet);
for s = 1:nshots
    h_probe(s) = plot(t_intensity{b(s)},intensity_probe{b(s)},'LineWidth',2);
end
title([magConfig{xp},' xp ',num2str(xp)])
box on
grid on
set(gcf,'color','w')
ylim([0,2.5e3])
l_intensity_probe = legend([h_probe],num2str(rfPwrNet(b)'));

SaveFig = 1;
if SaveFig
    saveas(gcf,[magConfig{xp}],'tiffn')
end

% The change in intensity appears to be in the linear regime, since for the
% worst case it spans less than 2e3 units