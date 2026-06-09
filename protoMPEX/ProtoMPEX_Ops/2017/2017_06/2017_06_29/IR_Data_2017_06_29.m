close all 
clear all

% ECH ON
shot = 15260;
% ECH OFF
shot = 15251;
% ECH ON
shot = 15247;
% % ECH OFF
shot = 15208;
a = '\\mpexserver\ProtoMPEX_Data\IR_Camera\2017_06_27';

% June 29th 2017
% First few shots after the N2 prefill
shot = 15399;
a = '\\mpexserver\ProtoMPEX_Data\IR_Camera\2017_06_29';

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
for s = 1:40
   imshow(Data(:,:,s),[]);
   drawnow;
   title(num2str(s))
   pause(0.05)
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
Rloc = 352;
figure; hold on
for s = 2:60
   S_RadialHeatProfile(:,s) = Data(:,Rloc,s)-Data(:,Rloc,12); 
   RadialHeatProfile2D(:,:,s) = Data(:,:,s)-Data(:,:,s-1); 
   plot(S_RadialHeatProfile(:,s))
   ylim([0,1e4])
   drawnow;
   pause(0.05)
end

figure
surf(S_RadialHeatProfile,'LineStyle','none')
%%
close all
figure
for s = 1:29 % use 13 to 22
   surf(RadialHeatProfile2D(150:450,250:500,s),'lineStyle','none');
   view([-70,70])
    zlim([0,4000])
    xlim([0,300])
    ylim([0,300])
    drawnow;
   title(num2str(s))
   pause(0.1)
end

figure;
subplot(1,2,1); hold on
s = 25;
surf(RadialHeatProfile2D(150:450,250:500,s),'lineStyle','none');
axis('square')
    zlim([0,4000])
    xlim([0,300])
    ylim([0,300])
    
subplot(1,2,2); hold on
s = 20;
surf(RadialHeatProfile2D(150:450,250:500,s),'lineStyle','none');
axis('square')
    zlim([0,4000])
    xlim([0,300])
    ylim([0,300])
return
%%
figure; 
contourf(RadialHeatProfile2D(210:300,320:410,21),20,'LineStyle','none'); colormap('hot'); caxis([0,3000])
%%
close all
figure; 
surf(S_RadialHeatProfile(:,1:60),'lineStyle','none')

% Data(:,359,25)-Data(:,359,24)
% Data(:,:,25)-Data(:,:,15)
figure; hold on
for s = 1:45
    S_RPD(:,s) = sgolay_t(S_RadialHeatProfile(:,s),3,7);
    RPD2D(:,:,s) = sgolay_t(RadialHeatProfile2D(:,:,s),3,27);
    plot(S_RPD(:,s))
end
% xlim([150,350])

%%
figure; hold on
dt = 19.8543/1000;
t = [1:size(S_RPD,2)-1]*dt -0.45 + 4.15 ;
plot(t,diff(S_RPD(255,:)),'k.-')
plot(t,diff(S_RPD(228,:)),'r.-')

for s = 1:44
    ff(:,s) = S_RPD(:,s+1)-S_RPD(:,s);
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
close all
figure
for s = 10:28;
   surf(ff2D(:,:,s),'lineStyle','none');
    view([-75,50])
%     view([0,90])
    zlim([0,1000])
    caxis([10,1000])
    colormap('hot')
    ylim([150,450])
    xlim([200,500])
%     axis('square')
   drawnow;
   title(num2str(s))
   pause(0.1)
end

%%
figure;
surf(ff2D(:,:,20),'lineStyle','none');


%%
% close all

f0 = 15;
f1 = 30;

rngy = [280:480] - 30;
rngx = [100:310]  + 60;
IR_f10 = Data(rngx,rngy,f1)-Data(rngx,rngy,f0);
IR_f10n = IR_f10/max(max(IR_f10));

a = 30e4;

% Delta emission
% absolute
figure;
contourf(IR_f10,80,'LineStyle','none')
caxis_max = 350;
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


