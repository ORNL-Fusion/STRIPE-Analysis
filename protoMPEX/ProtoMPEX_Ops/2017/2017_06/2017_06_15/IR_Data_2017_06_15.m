close all 
clear all


shot = 14872;
a = '\\mpexserver\ProtoMPEX_Data\IR_Camera\2017_06_09';

% shot = 14881;
% a = '\\mpexserver\ProtoMPEX_Data\IR_Camera\2017_06_14';

% shot = 14983;
% a = '\\mpexserver\ProtoMPEX_Data\IR_Camera\2017_06_15';

% shot = 14984;
% a = '\\mpexserver\ProtoMPEX_Data\IR_Camera\2017_06_15';

shot = 14987;
a = '\\mpexserver\ProtoMPEX_Data\IR_Camera\2017_06_15';

shot = 15098;
a = '\\mpexserver\ProtoMPEX_Data\IR_Camera\2017_06_16';

shot = 15100;
a = '\\mpexserver\ProtoMPEX_Data\IR_Camera\2017_06_16';

shot = 15101;
a = '\\mpexserver\ProtoMPEX_Data\IR_Camera\2017_06_16';

shot = 15103;
a = '\\mpexserver\ProtoMPEX_Data\IR_Camera\2017_06_16';

shot = 15104;
a = '\\mpexserver\ProtoMPEX_Data\IR_Camera\2017_06_16';

shot = 15111;
a = '\\mpexserver\ProtoMPEX_Data\IR_Camera\2017_06_16';

mov = 1;

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
if mov
figure(1)
for s = 10:40
   imshow(Data(:,:,s),[]);
   drawnow;
   title(num2str(s))
   pause(0.5)
end
end
%%
% figure(2)
% for s = 1:60
%    imshow(Data(:,:,s)-Data(:,:,10),[]);
%    drawnow;
%    title(num2str(s))
%    pause(0.01)
% end
%%
figure(3)
Rloc = 352;
for s = 2:60
   RadialHeatProfile(:,s) = Data(:,Rloc,s)-Data(:,Rloc,20); 
   RadialHeatProfile2D(:,:,s) = Data(:,:,s)-Data(:,:,s-1); 
   plot(RadialHeatProfile(:,s));
   ylim([0,16000])
   drawnow;
   title(num2str(s))
   pause(0.01)
end

% figure(4)
% for s = 1:40
%    surf(RadialHeatProfile2D(:,:,s),'lineStyle','none');
%    drawnow;
%    title(num2str(s))
%    pause(0.01)
% end
%%
close all
figure; 
surf(RadialHeatProfile(:,1:60),'lineStyle','none')

% Data(:,359,25)-Data(:,359,24)
% Data(:,:,25)-Data(:,:,15)
figure; hold on
for s = 1:45
    RPD(:,s) = sgolay_t(RadialHeatProfile(:,s),3,7);
    RPD2D(:,:,s) = sgolay_t(RadialHeatProfile2D(:,:,s),3,27);
    plot(RPD(:,s))
end
xlim([150,350])

% figure; 
% surf(diff(RPD,2),'lineStyle','none')

figure; hold on
dt = 19.8543/1000;
t = [1:size(RPD,2)-1]*dt -0.45 + 4.15 ;
plot(t,diff(RPD(273,:)),'k.-')
plot(t,diff(RPD(207,:)),'r.-')

for s = 1:44
    ff(:,s) = RPD(:,s+1)-RPD(:,s);
    ff2D(:,:,s) = RPD2D(:,:,s+1)-RPD2D(:,:,s);
end

figure;
R = 1:length(ff(:,1));
%contourf(t(20:40),R,ff(:,20:40),20,'lineStyle','none')
%caxis([0,800])
% surf(t(20:40),R,ff(:,20:40),'lineStyle','none')
surf(ff(:,:),'lineStyle','none')
colormap('hot')
%%
figure
for s = 1:39;
   surf(ff2D(:,:,s),'lineStyle','none');
    view([-75,50])
%     view([0,90])
    zlim([0,1100])
    caxis([0,100])
    colormap('hot')
    ylim([100,400])
    xlim([200,500])
    axis('square')
   drawnow;
   title(num2str(s))
   pause(0.1)
end

%%
% close all

f0 = 20;
f1 = 40;
% 26 begins, 28 looks good, 30 helicon has developed
% 36 is fully developed

% for 300 ms use frame 42 for fully developed heat,

rngy = [280:480] - 30;
rngx = [100:310]  + 60;
IR_f10 = Data(rngx,rngy,f1)-Data(rngx,rngy,f0);
IR_f10n = IR_f10/max(max(IR_f10));

a = 30e4;

% Delta emission
% absolute
figure;
contourf(IR_f10,80,'LineStyle','none')
caxis_max = 300;
caxis([0,caxis_max]*a*1e-4)
axis(gca,'equal')
colormap('jet')
th = title([FILENAME,', caxis: ',num2str(caxis_max),', f = ',num2str(f1)])
set(th,'FontSize',14)

% normalized
% figure;
% contourf(IR_f10n,80,'LineStyle','none')
% zlim([0,1])
% caxis([0,0.9])
% axis(gca,'equal')
% colormap('jet')
% title(FILENAME)


