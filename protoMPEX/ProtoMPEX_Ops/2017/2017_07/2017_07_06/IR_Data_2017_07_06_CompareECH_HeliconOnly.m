close all 
clear all

% ECH ON, 
shot(1) = 15579;
a{1} = '\\mpexserver\ProtoMPEX_Data\IR_Camera\2017_07_06';
% % ECH ON
shot(2) = 15585;
a{2} = '\\mpexserver\ProtoMPEX_Data\IR_Camera\2017_07_06';

mov = 1;

for s = 1:2
PATHNAME = [a{s},'\'];
FILENAME = ['shot ',num2str(shot(s)),'.seq'];
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
if mov
figure(1)
dT = 0.001;
for s = 1:40
   subplot(1,2,1)
   imshow(Data{1}(:,:,s),[]);
   title(num2str(s))
   pause(dT)

   subplot(1,2,2)
   imshow(Data{2}(:,:,s),[]);
   title(num2str(s))
   pause(dT)
end
end

%%
Rloc = 352;
figure; hold on
for s = 2:60
   S_RadialHeatProfile{1}(:,s) = Data{1}(:,Rloc,s)-Data{1}(:,Rloc,12); 
   RadialHeatProfile2D{1}(:,:,s) = Data{1}(:,:,s)-Data{1}(:,:,s-1); 
   
   S_RadialHeatProfile{2}(:,s) = Data{2}(:,Rloc,s)-Data{2}(:,Rloc,12); 
   RadialHeatProfile2D{2}(:,:,s) = Data{2}(:,:,s)-Data{2}(:,:,s-1); 

   subplot(1,2,1); hold on
   plot(S_RadialHeatProfile{1}(:,s))
   ylim([0,2e4])
   title(num2str(shot(1)))
%    drawnow;
   pause(0.05)
   
      subplot(1,2,2); hold on
   plot(S_RadialHeatProfile{2}(:,s))
   ylim([0,2e4])
   title(num2str(shot(2)))
%    drawnow;
   pause(0.05)
end

%%
figure;
subplot(1,2,1);
surf(S_RadialHeatProfile{1},'LineStyle','none')
set(gca,'PlotBoxAspectRatio',[1 1 1])
xlim([0,50])
zlim([0,15e3])

subplot(1,2,2);
surf(S_RadialHeatProfile{2},'LineStyle','none')
set(gca,'PlotBoxAspectRatio',[1 1 1])
xlim([0,50])
zlim([0,15e3])
%%
% close all
figure
set(gcf,'Position',[  0.0750    0.1977    1.1347    0.4200]*1e3,'color','w')
el = 90;

for s = 8:36 % use 13 to 22
    
           subplot(1,2,2)
   surf(RadialHeatProfile2D{2}(150:450,250:500,s),'lineStyle','none');
   HeatFluxOnAxis{1}(s) = max(max(RadialHeatProfile2D{1}(150:450,250:500,s)));
   view([-70,el])
    zlim([0,5000])
        caxis([0,5000])
    xlim([0,300])
    ylim([0,300])
   title(num2str(s))
   colormap('jet')
%    drawnow
%    pause(0.1)
   
    subplot(1,2,1)
   surf(RadialHeatProfile2D{1}(150:450,250:500,s),'lineStyle','none');
   HeatFluxOnAxis{2}(s) = max(max(RadialHeatProfile2D{2}(150:450,250:500,s)));
   line([125,125],[0,500],[0,0])
   view([-70,el])
    zlim([0,5000])
    caxis([0,5000])
    xlim([0,300])
    ylim([0,300])
   title(num2str(s))
   colormap('jet')

   pause (0.15)
  
end

figure; hold on
h(1) = plot(HeatFluxOnAxis{1},'r')
h(2) = plot(HeatFluxOnAxis{2},'k')
legend(h,num2str(shot'))
%%
% close all
PlotType = 1;

figure; 
subplot(2,1,1)
F = RadialHeatProfile2D{1}(:,370,:);
if PlotType
    contourf(permute(F,[3 1 2])',20,'LineStyle','none')
else
    surf(permute(F,[3 1 2])','LineStyle','none')
end
xlim([18,35])
ylim([150,400])
zlim([0,5000]);caxis([0,5000])
   colormap('jet')
title(num2str(shot(1)))

subplot(2,1,2)
F = RadialHeatProfile2D{2}(:,370,:);
if PlotType
    contourf(permute(F,[3 1 2])',20,'LineStyle','none')
else
    surf(permute(F,[3 1 2])','LineStyle','none')
end
xlim([18,35])
ylim([150,400])
zlim([0,5000]);caxis([0,5000])
   colormap('jet')
title(num2str(shot(2)))
get(gcf,'color','w')

% figure; 
% subplot(2,1,1)
% F = RadialHeatProfile2D{1}(:,370,:);
% % contourf(permute(F,[3 1 2])',20,'LineStyle','none')
% surf(permute(F,[3 1 2])','LineStyle','none')
% xlim([18,35])
% ylim([150,400])
% zlim([0,4000]);caxis([0,3500])
%    colormap('jet')
% title(num2str(shot(1)))
% 
% subplot(2,1,2)
% F = RadialHeatProfile2D{2}(:,370,:);
% % contourf(permute(F,[3 1 2])',20,'LineStyle','none')
% surf(permute(F,[3 1 2])','LineStyle','none')
% xlim([18,35])
% ylim([150,400])
% zlim([0,4000]);caxis([0,3500])
%    colormap('jet')
% title(num2str(shot(2)))
% get(gcf,'color','w')

% figure;
% subplot(1,2,1); hold on
% s = 25;
% surf(RadialHeatProfile2D(150:450,250:500,s),'lineStyle','none');
% axis('square')
%     zlim([0,4000])
%     xlim([0,300])
%     ylim([0,300])
%     
% subplot(1,2,2); hold on
% s = 20;
% surf(RadialHeatProfile2D(150:450,250:500,s),'lineStyle','none');
% axis('square')
%     zlim([0,4000])
%     xlim([0,300])
%     ylim([0,300])
return
%%
figure; 
contourf(RadialHeatProfile2D(210:300,320:410,21),20,'LineStyle','none'); colormap('jet'); caxis([0,3000])
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
colormap('jet')
%%
close all
figure
for s = 10:28;
   surf(ff2D(:,:,s),'lineStyle','none');
    view([-75,50])
%     view([0,90])
    zlim([0,1000])
    caxis([10,1000])
    colormap('jet')
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


