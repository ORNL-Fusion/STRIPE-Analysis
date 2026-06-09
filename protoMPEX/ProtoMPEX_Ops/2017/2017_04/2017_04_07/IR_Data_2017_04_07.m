% IR data analysis March 29th 2017
close all 
clear all

shot = 13767; 
shot = 13770; % with ECH

shot = 13781; % PS1/2 = 4000 A, TR2 = 180 A, TR1 =600 A, U puff, ECH ON
% shot = 13782; % PS1/2 = 4000 A, TR2 = 180 A, TR1 =600 A, U puff, ECH OFF

% shot = 13782; % PS1/2 = 2500 A, TR2 = 100 A, TR1 =600 A, U puff, ECH ON, RF trip

% shot = 13788; % PS1/2 = 4000 A, TR2 = 160 A, TR1 =600 A, U puff, ECH ON, RF FP
% shot = 13790; % PS1/2 = 3500 A, TR2 = 140 A, TR1 =600 A, U puff, ECH ON, RF FP
% shot = 13791; % PS1/2 = 3000 A, TR2 = 120 A, TR1 =600 A, U puff, ECH ON, RF FP

% shot = 13792; % PS1/2 = 2500 A, TR2 = 100 A, TR1 =600 A, U puff, ECH ON, RF trip

shot = 13747;
a = '\\mpexserver\ProtoMPEX_Data\IR_Camera\2017_04_07';
PATHNAME = [a,'\'];
FILENAME = ['shot ',num2str(shot),'.seq'];
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
%         Data(:,:,fr) = im;
        Data(:,:,fr) = im(end:-1:1,:);
        fr = fr + 1;
    end
end

%%
close all

f0 = 20;
f1 = 36;
% 25 begins, 28 looks good, 30 helicon has developed
% 36 is fully developed
% frame 32 the ECH is just starting

rngy = [280:480] - 20;
rngx = [100:310]  + 40;
IR_f10 = Data(rngx,rngy,f1)-Data(rngx,rngy,f0);
IR_f10n = IR_f10/max(max(IR_f10));

a = 30e4

figure;
h(1) = surf(Data(:,:,f0),'LineStyle','none')
zlim([1.48,1.52]*a)
view([0,90])
caxis([1.48,1.54]*a)
colormap('jet')

figure;
h(2) = surf(Data(:,:,f1),'LineStyle','none')
zlim([1.48,1.52]*a)
caxis([1.48,1.54]*a)
view([0,90])
colormap('jet')

% Delta emission
% absolute
figure;
contourf(IR_f10,80,'LineStyle','none')
% zlim([1.3,1.7]*1e4 - 1.3e4)
caxis([0,150]*a*1e-4)
axis(gca,'equal')
colormap('jet')
title(FILENAME)

% normalized
% figure;
% contourf(IR_f10n,80,'LineStyle','none')
% zlim([0,1])
% caxis([0,0.9])
% axis(gca,'equal')
% colormap('jet')
% title(FILENAME)


