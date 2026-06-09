% IR data analysis MAY 13th 2016
close all 
clear all
PATHNAME = [cd,'\'];
FILENAME = ['shot ',num2str(8833),'.seq'];
videoFileName=[PATHNAME FILENAME];

% shotlist = 8800 + [10 , 11, 12, 13, 15, 17, 18, 19, 20, 21, 23, 25, 27, 29, 31, 32, 33, 34, 83, 84, 85, 87,88];
% TR2      =        [210,270,290,310,330,350,370,390,410,430,450,470,500,600,230,200,170,140,170,140,110, 60, 0];

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
f1 = 44;
% 38 begins
% 44 fully developed

rngy = [280:480];
rngx = [100:310];
IR_f10 = Data(rngx,rngy,f1)-Data(rngx,rngy,f0);
IR_f10n = IR_f10/max(max(IR_f10));

figure;
h(1) = surf(Data(:,:,f0),'LineStyle','none')
zlim([1.48,1.52]*1e4)
view([0,90])
caxis([1.48,1.54]*1e4)
colormap('jet')

figure;
h(2) = surf(Data(:,:,f1),'LineStyle','none')
zlim([1.48,1.52]*1e4)
caxis([1.48,1.54]*1e4)
view([0,90])
colormap('jet')

% Delta emission
% absolute
figure;
contourf(IR_f10,80,'LineStyle','none')
% zlim([1.3,1.7]*1e4 - 1.3e4)
caxis([0,250])
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

%% 
% close all

mov = 1;
if mov
% figure(1)
for s = 2:70
%    imshow(Data(:,:,s)- 0*Data(:,:,1),[]);
   dTdt_a(:,:,s) = Data(:,:,s)- Data(:,:,s-1);
   dTdt_b(:,:,s) = Data(:,:,s+1)- Data(:,:,s-1);
   dTdt_c(:,:,s) = (-11/6)*Data(:,:,s)+(3)*Data(:,:,s+1)-(3/2)*Data(:,:,s+2)+(1/3)*Data(:,:,s+3);

%    drawnow;
%    title(num2str(s))
%    pause(0.05)
end
end

%%
close all
figure
for s = 40%30:60
   surf(dTdt_b(:,:,s),'LineStyle','none')
   set(gcf,'Position',[641.0000  299.0000  279.3333  218.6667])
   view([0,90])
   title(num2str(s))
   zlim([-20,150])
   xlim([250,500])
   ylim([80,300])
   axis('square')
   caxis([0,50])
   colormap(flipud(hot))
    pause(0.0001)
end




