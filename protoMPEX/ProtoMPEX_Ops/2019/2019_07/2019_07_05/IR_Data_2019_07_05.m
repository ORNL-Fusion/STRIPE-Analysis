close all 
clear all


% NOTES:
% - Determine the best shot to do determine the emissivity calibration
% - For a shot, need to create the WinArea variable

 
shot(1) = 26809;  
shot(1) = 26806;  
shot(1) = 26805;  
shot(1) = 26804;  

a{1} = '\\mpexserver\ProtoMPEX_Data\IR_Camera\2019_07_05';


shot(1) = 27435;  

a{1} = '\\mpexserver\ProtoMPEX_Data\IR_Camera\2019_08_09';

mov = 1;

for s = 1:length(shot)
PATHNAME = [a{s},'\'];
try
FILENAME = ['Shot ',num2str(shot(s)),'.seq'];
catch
end
% FILENAME = ['Shot-0',num2str(shot(s)),'.seq'];    
videoFileName=[PATHNAME FILENAME];


% Load the Atlats SDK
atPath = getenv('FLIR_Atlas_MATLAB');
atImage = strcat(atPath,'Flir.Atlas.Image.dll');
asmInfo = NET.addAssembly(atImage);
%open the IR-file'
file = Flir.Atlas.Image.ThermalImageFile(videoFileName);
seq = file.ThermalSequencePlayer();
%Get the pixels
img = seq.ThermalImage.ImageProcessing.GetPixelsArray;
im = double(img);
    
RawData{s}(:,:,1) = im;
fr = 1;
if(seq.Count > 1)
    while(seq.Next())
        img = seq.ThermalImage.ImageProcessing.GetPixelsArray;
        im = double(img);
        RawData{s}(:,:,fr) = im(end:-1:1,:);         
        fr = fr + 1;
    end
end
end

msgbox('Data acquired')

% =========================================================================
%  Calculate surface temperature
% =========================================================================
emissivity = 0.8;
[TempMeasured,dTMeasured] = IntensityTempConv2(emissivity,0.2*emissivity,RawData{1},seq);
% framerate = 100;
framerate = seq.FrameRate;
dt = 1/framerate;
time = 0:dt:(fr-2)*dt;

% =========================================================================
% Determine the start of the RF pulse:
% =========================================================================
Nx = size(TempMeasured(:,:,:),1);
Ny = size(TempMeasured(:,:,:),2);
Nz = size(TempMeasured(:,:,:),3);

PixelCount = Nx*Ny;

for s = 1:Nz
   MeanTempMeasured(s) = sum(sum(TempMeasured(1:20:end,1:20:end,s)))/PixelCount; 
end

n1 = 40;
D0 = TempMeasured(:,:,n1);

figure; 
hold on
plot(time,MeanTempMeasured,'k.')
title('Sum')
ylabel('surface T [C]')
xlabel('t [s]')
xlim([3,7])

%% Movie
figure;
n1 = 100
for ii = n1:4:Nz
    surf(TempMeasured(:,:,ii)-TempMeasured(:,:,n1),'LineStyle','none')
    view([0,90])
    caxis([0,50])
    title(['frame: ',num2str(ii)])
    axis('equal')
    drawnow
end
