% Preprocessor 1
% This script loads the IR data from experiments performed on the
% 2020_02_04
% This script deals only with the raw data from the IR camera which we
% refer to as intensityRaw and intensity

% The objective of preprocessor 1 is to:
% 1 - Identify the start of the RF pulse
% 2 - Identify the window area from the intensity data
% 3 - Extract intensity data during RF pulse
% 4 - Save the data needed for post-processing

clear all
close all

xp = 1;
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
a{1} = ['\\mpexserver\ProtoMPEX_Data\IR_Camera\2020_02_05'];
a{2} = ['\\mpexserver\ProtoMPEX_Data\IR_Camera\2019_12_19'];

thermalParam.ExternalOpticsTransmission = 0.7;
thermalParam.AtmosphericTemperature = 24;
thermalParam.Distance = 1;
thermalParam.ExternalOpticsTemperature = 24;
thermalParam.ReferenceTemperature = 24;
thermalParam.Transmission = 1;
thermalParam.RelativeHumidity = 0;
thermalParam.ReflectedTemperature = 24;

for s = 1:length(shot)
    % Specify IR file name
    pathName = [a{addressLoc(s)},'\'];
    fileName = ['Shot-0',num2str(shot(s)),'.seq'];

    % Get the data
    [intensityRaw{s},seq{s}] = ExtractDataFromSeqFile(fileName,pathName,thermalParam);
%     [intensityRaw{s},seq{s}] = ExtractDataFromSeqFile(fileName,pathName,[]);
end

nshots = length(shot);


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
% 
'Data subset'
% Offsets to reject data from the entire time trace
n1_offset = 10;
n2_offset = 10;
% Offsets to widen the RF pulse window
n_Before = 10;
n_After = 20;
frameRate = 100;
dt = 1/frameRate;

try
    hf = findobj('Tag','DefineDataRng'); close(hf);
end

tic
% Define start and end of pulse and define data range
for s = 1:nshots    
    [Nx,Ny,Nz] = size(intensityRaw{s}); 

    for ii = 1:Nz
       intensityRaw_mean{s}(ii) = mean(mean(intensityRaw{s}(1:20:Nx,1:20:Ny,ii))); 
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

% Plot data to confirm correctness of calculation
if 1
    figure; 
    set(gcf,'Tag','DefineDataRng')
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

% Extract the relevent data from "intensityRaw":
for s = 1:nshots
    intensity{s} = intensityRaw{s}(:,:,rng{s});
    t_intensity{s} = 0:dt:(length(rng{s})-1)*dt;
end
toc

%%
% =========================================================================
% 2 - Determine the window area
% =========================================================================
'Determine the window area'
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
toc

%% 
% =========================================================================
% =========================================================================

cmptThresholds = 1;

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

edgeSelectMode = 1;
switch edgeSelectMode
    case 1 % Only one instance
        instances = 1;
    case 2 % Define all instances
        instances = 1:nshots;
end

for s = instances
    MySurf1(abs(intensity1_gabs_smooth{s}))
    axis('equal')
    view([0,90])
    set(gcf,'Tag','intensity0_hist')
    colormap(flipud(hot))
    caxis([0,intensity1_mpg(s)])
    zlim([0,intensity1_mpg(s)*2])

    [nx_mpg{s},ny_mpg{s}] = find(abs(intensity1_gabs_smooth{s})>intensity1_mpg(s));
    hold on
    for ii = 1:length(ny_mpg{s})-1
        mpg{s}(ii) = intensity1_gabs_smooth{s}(nx_mpg{s}(ii),ny_mpg{s}(ii));
        plot3(ny_mpg{s}(ii),nx_mpg{s}(ii),abs(mpg{s}(ii)),'g.')
    end

    switch cmptThresholds
    case 0
          load(['preprocessData_xp_',num2str(xp)],'nx_mpg_ginput','ny_mpg_ginput');
    case 1
          [nx_mpg_ginput(s),ny_mpg_ginput(s)] = ginput(1);
          nx_mpg_ginput(s) = floor(nx_mpg_ginput(s));
          ny_mpg_ginput(s) = floor(ny_mpg_ginput(s));
          save(['preprocessData_xp_',num2str(xp)],'nx_mpg_ginput','ny_mpg_ginput');
    end
end

if edgeSelectMode == 1;
    nx_mpg_ginput = ones(nshots)*nx_mpg_ginput(1);
    ny_mpg_ginput = ones(nshots)*ny_mpg_ginput(1);
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

magConfig{1} = 'Window limit';
magConfig{2} = 'Window limit';
magConfig{3} = 'Window limit';
magConfig{4} = 'Window limit';

try
    hf = findobj('Tag','intensityRawVid');
    close(hf)
end

s = 1;

switch viewType
    case {'Middle'}
        xCenter = 194;
        yCenter = 340;
    case {'Bottom'}
        xCenter = 103;
        yCenter = 422;
    case {'Top'}
        xCenter = 156;
        yCenter = 430;
    case {'Bottom view gas scan'}
        xCenter = 103;
        yCenter = 422;
    case {'Bottom view Upstream limit'}
        xCenter = 103;
        yCenter = 422;
end

rngy_intensity_probe = [(yCenter-10):(yCenter+10)];
rngx_intensity_probe = [(xCenter-10):(xCenter+10)];

for fr = 1:1:size(intensity{s},3)
    intensity_windowArea(:,:,fr) = intensity{s}(:,:,fr).*windowArea{s}  -  intensity{s}(:,:,1).*windowArea{s};
end

intenMax = max(max(max(intensity_windowArea)));
intenMin = min(min(min(intensity_windowArea)));

figure;
set(gcf,'color','w')
set(gcf,'Tag','intensityRawVid')
for fr = 1:1:size(intensity{s},3)
    surf(intensity_windowArea(:,:,fr),'LineStyle','none')
    view([0,90])
    zlim([0,2000])
    caxis([0,2000])
    colormap(flipud(hot))
    title(['shot: ',num2str(shot(s)),', frame: ',num2str(fr),', P_{RF}: ',num2str(rfPwrNet(s)),' kW, ',magConfig{xp}])
    axis('equal')
    xlim([0,size(intensity0_noOffset{s},2)])
    ylim([0,size(intensity0_noOffset{s},1)])
    hold on
    plot3(rngy_intensity_probe,rngx_intensity_probe,intensity_windowArea(rngx_intensity_probe,rngy_intensity_probe,fr),'k.')
    plot3(flip(rngy_intensity_probe),rngx_intensity_probe,intensity_windowArea(rngx_intensity_probe,rngy_intensity_probe,fr),'k.')
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
title(['Intensity ',viewType,', ',magConfig{xp}])
box on
grid on
set(gcf,'color','w')
ylim([0,2.5e3])
l_intensity_probe = legend([h_probe],num2str(rfPwrNet(b)'));

SaveFig = 0;
if SaveFig
    figureName = ['Intensity ',viewType];
    saveas(gcf,figureName,'tiffn')
end