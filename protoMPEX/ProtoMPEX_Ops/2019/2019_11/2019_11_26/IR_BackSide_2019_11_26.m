clear all
close all

XP = 1; 

% Define address location to get IR data
RootAddress = '\\mpexserver\ProtoMPEX_Data\IR_Camera';
a{1} = [RootAddress,'\2019_11_25'];
AddressLoc =         [1  ];     

switch XP
    case 1
        shot       = 28000 + [379];

        
    case 2
 
    case 3  

end

%% 
% =========================================================================
% Load seq files and create RawData
% =========================================================================

for s = 1:length(shot)
% Load the Atlats SDK
atPath = getenv('FLIR_Atlas_MATLAB');
atImage = strcat(atPath,'Flir.Atlas.Image.dll');
asmInfo = NET.addAssembly(atImage);
%open the IR-file'
PATHNAME = [a{AddressLoc(s)},'\'];
FILENAME = ['Shot ',num2str(shot(s)),'.seq'];
videoFileName=[PATHNAME FILENAME];
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

%Get the pixels
img = seq{s}.ThermalImage.ImageProcessing.GetPixelsArray;
im = double(img);
    
RawData{s}(:,:,1) = im;
fr = 1;
if(seq{s}.Count > 1)
    while(seq{s}.Next())
        img = seq{s}.ThermalImage.ImageProcessing.GetPixelsArray;
        im = double(img);
        RawData{s}(:,:,fr) = im(end:-1:1,:);         
        fr = fr + 1;
    end
end
end

%% Determine start of RF pulse
try 
    figure(09); close
end

rngx = [340:370];
rngy = [140:170];

for s = 1:size(RawData,1)
    N = numel(RawData{1}(:,:,1));
    N_Axis = numel(RawData{s}(rngy,rngx,1));

    for ii = 1:size(RawData{1},3)
        IntensityIntegrated(ii) = sum(sum(RawData{s}(1:20:end,1:20:end,ii)))/N;
        IntensityIntegrated_Axis(ii) = sum(sum(RawData{s}(rngy,rngx,ii)))/N;

    end
end

figure(09); 
hold on
plot(IntensityIntegrated_Axis - mean(IntensityIntegrated_Axis(1:50)) ,'r')
plot(IntensityIntegrated - mean(IntensityIntegrated(1:50)),'k')
yyaxis right
plot(diff(IntensityIntegrated),'g')

%%
try 
    figure(10); close
end

rngx = [340:365];
rngy = [140:165];

figure(10)
for ii = 100:200
    hold on
    surf(RawData{1}(:,:,ii),'LineStyle','none')
    plot3(rngx,rngy,ones(size(rngx))*2e4,'r.')
    hold off
    title(['Frame: ',num2str(ii),' , shot: ',num2str(shot)])
%     set(gca,'PlotBoxAspectRatio',[2.5 1 1])
    axis image
    caxis([1.37,1.7]*1e4)
    colorbar
%     zlim([1.37e4,1.54e4])
    view([0,90])
    drawnow
end

for ii = 1:size(RawData{1}(3,3,:),3)-1
    CentralTemp(ii) = 0.5*(mean(mean(RawData{1}(rngy,rngx,ii))) + mean(mean(RawData{1}(rngy,rngx,ii+1))));
end

figure; 
hold on
plot(CentralTemp-mean(CentralTemp(1:50)))
plot(10*diff(CentralTemp),'r','LineWidth',2)

return

%%
for ii = 2:498
    dI(:,:,ii) = (RawData{1}(:,:,ii+1)- RawData{1}(:,:,ii-1))/2;
end
%%
try 
    figure(11); close
end

figure(11)
for ii = 125:200
    surf(dI(:,:,ii),'LineStyle','none')
    title(['Frame: ',num2str(ii),' , shot: ',num2str(shot)])
    axis image
    caxis([-10,0.8e3])
%     colorbar
    view([0,90])
    drawnow
    pause(0.1)
end