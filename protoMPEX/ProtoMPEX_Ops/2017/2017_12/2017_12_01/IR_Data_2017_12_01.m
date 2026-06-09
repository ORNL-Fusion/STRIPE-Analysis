close all 
clear all

% ECH OFF, 
shot(1) = 18158;
shot(1) = 18194;
shot(1) = 18190; % TR2 = 260 A
% shot(1) = 18210; % TR2 = 120 A
shot(1) = 18165; % TR2 = 220 A
shot(1) = 18179; % TR2 = 160 A

a{1} = '\\mpexserver\ProtoMPEX_Data\IR_Camera\2017_12_01';

mov = 1;

for s = 1:length(shot)
PATHNAME = [a{s},'\'];
FILENAME = ['Shot ',num2str(shot(s)),'.seq'];
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
    
Data{s}(:,:,1) = im;
fr = 1;
if(seq.Count > 1)
    while(seq.Next())
        img = seq.ThermalImage.ImageProcessing.GetPixelsArray;
        im = double(img);
%         Data(:,:,fr) = im;
        Data{s}(:,:,fr) = im(end:-1:1,:);         
        fr = fr + 1;

    end
end
end
%%
close all
if mov
figure(1)
dT = 0.001;
rngX = [200:-1:1];
rngY = [400:-1:270];
for s = 45:120
   imshow(Data{1}(:,:,s) - 0*Data{1}(:,:,51),[]);
%    imshow(Data{1}(rngX,rngY,s)- Data{1}(rngX,rngY,s-1),[]);
%  y0(s) = Data{1}(round(mean(rngX)),round(mean(rngY)),s);
  y0(s) = Data{1}(102,336,s);
  y1(s) = Data{1}(136,315,s);

   title(num2str(s))
   pause(dT)
end
end
line([336,336],[102,102],'marker','.')
line([315,315],[136,136],'marker','.')

figure; hold on
plot(y0)
plot(y1,'r')

%%
n = 1;
N{n} = length(Data{n}(1,1,:))-1;
for s = 2:N{n}
       dTdt_a{n}(:,:,s) = Data{n}(:,:,s)  - Data{n}(:,:,s-1);
       dTdt_b{n}(:,:,s) = Data{n}(:,:,s+1)- Data{n}(:,:,s-1);
end
%%
rngy = [280:480];
rngx = [100:310];
figure;
n = 1;
for s = 1:N{n}
   surf(dTdt_b{n}(:,:,s),'LineStyle','none')
   set(gcf,'Position',[600 300  230  200])
   view([0,90])
   title(num2str(s))
   zlim([-20,150])
%     xlim([0,210])
%     ylim([0,215])
    axis('square')
   caxis([0,50])
   colormap(flipud(hot))
       set(gca,'XTickLabel',[])
           set(gca,'YTickLabel',[])
           box on
    pause(0.05)

end