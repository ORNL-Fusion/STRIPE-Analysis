close all 
clear all


shot = 6514;

mov = 1;
UploadType = 2;

switch UploadType
    case 1
        a = '\\mpexserver\ProtoMPEX_Data\IR_Camera\2015_12_11';
    case 2
        a = cd;
end

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
for s = 1:40
   imshow(Data(:,:,s),[]);
%    imshow(Data(190:230,300:345,s),[]);
   drawnow;
   title(num2str(s))
   pause(0.5)
end
end
%%
Rloc = 323;
for s = 2:60
   S_RadialHeatProfile(:,s) = Data(:,Rloc,s)-Data(:,Rloc,12); 
   RadialHeatProfile2D(:,:,s) = Data(:,:,s)-Data(:,:,s-1); 
end

figure
for s = 1:60
   surf(RadialHeatProfile2D(190:230,300:345,s),'lineStyle','none');
    zlim([0,1000])
    xlim([0,50])
    ylim([0,50])
    drawnow;
   title(num2str(s))
   pause(0.05)
end
%%
close all
figure; 
surf(S_RadialHeatProfile(190:230,1:end),'lineStyle','none')
xlim([30,60])

figure; hold on
for s = 1:59
    S_RPD(:,s) = sgolay_t(S_RadialHeatProfile(:,s),3,7);
    RPD2D(:,:,s) = sgolay_t(RadialHeatProfile2D(:,:,s),3,27);
    plot(S_RPD(:,s))
end
% xlim([150,350])

%%
figure; hold on
dt = 19.8543/1000;
t = [1:size(S_RPD,2)-1]*dt -0.45 + 4.15 ;
plot(t,diff(S_RPD(210,:)),'k.-')
plot(t,diff(S_RPD(199,:)),'r.-')

for s = 1:58
    ff(:,s) = S_RPD(:,s+1)-S_RPD(:,s);
    ff2D(:,:,s) = RPD2D(:,:,s+1)-RPD2D(:,:,s);
end

figure;
R = 1:length(ff(:,1));
%contourf(t(20:40),R,ff(:,20:40),20,'lineStyle','none')
%caxis([0,800])
% surf(t(20:40),R,ff(:,20:40),'lineStyle','none')
% surf(ff(:,:),'lineStyle','none')
contourf(ff(:,:),30,'lineStyle','none')

ylim([190,230])
xlim([20,60])


colormap('hot')
%%
figure
for s = 1:24;
   surf(ff2D(:,:,s),'lineStyle','none');
    view([-75,50])
%     view([0,90])
    zlim([0,1000])
    caxis([0,1000])
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


