%##### Load image #####
%close all 
clear all
[FILENAME, PATHNAME, FILTERINDEX] = uigetfile('*.jpg;*.seq', 'Choose IR file (jpg) or radiometric sequence (seq)');
videoFileName=[PATHNAME FILENAME];

% Load the Atlats SDK
atPath = getenv('FLIR_Atlas_MATLAB');
atImage = strcat(atPath,'Flir.Atlas.Image.dll');
asmInfo = NET.addAssembly(atImage);
%open the IR-file
file = Flir.Atlas.Image.ThermalImageFile(videoFileName);
seq = file.ThermalSequencePlayer();
%Get the pixels
img = seq.ThermalImage.ImageProcessing.GetPixelsArray;
im = double(img);
    
Data(:,:,1) = im;
fr = 1;
if(seq.Count > 1)
    while(seq.Next())
        img = seq.ThermalImage.ImageProcessing.GetPixelsArray;
        im = double(img);
        Data(:,:,fr) = im;
        fr = fr + 1;
    end
end

%%
close all

f0 = 20;
f1 = 36;

rngy = [280:480];
rngx = [120:320];
IR_f10 = Data(rngx,rngy,f1)-Data(rngx,rngy,f0);
IR_f10n = IR_f10/max(max(IR_f10));

% need to create a script that tracks the circumference of the plasma from
% frame 27 using the correlation function 

% how can we integrate and calculate the mean emission at the edge and then
% the mean emission at the core?

% there is also a GUI MATLAB file to obtain and play the seq files

figure;
h(1) = surf(Data(:,:,f0),'LineStyle','none')
zlim([1.3,1.7]*1e4)
view([0,90])
caxis([1.3,1.6]*1e4)
colormap('jet')

figure;
h(2) = surf(Data(:,:,f1),'LineStyle','none')
zlim([1.3,1.7]*1e4)
caxis([1.3,1.6]*1e4)
view([0,90])
colormap('jet')

% Delta emission
% absolute
figure;
contourf(IR_f10,80,'LineStyle','none')
zlim([1.3,1.7]*1e4 - 1.3e4)
caxis([1.3,1.57]*1e4 - 1.3e4)
axis(gca,'equal')
colormap('jet')
title(FILENAME)

% normalized
figure;
contourf(IR_f10n,80,'LineStyle','none')
zlim([0,1])
caxis([0,0.9])
axis(gca,'equal')
colormap('jet')
title(FILENAME)


