close all 
clear all

shot(1) = 26553;
a{1} = '\\mpexserver\ProtoMPEX_Data\IR_Camera\2018_06_06';

mov = 1;

for s = 1:length(shot)
PATHNAME = [a{s},'\'];
FILENAME = ['Shot ',num2str(shot(s)),'.seq'];
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
rngX = [240:-1:1]+0;
rngY = [450:-1:300]+0;
graph = 2;
tIR = (1:size(Data{1},3))*(1/98) + 4.17 - 32/98;

% Calculate a quantity that approximates the heat flux at the plate:
for s = 1:(size(Data{1},3)-1)
    % Time derivative:
    f = Data{1}(:,:,s+1) - Data{1}(:,:,s);
    dTdt(:,:,s) = f;
    % Integrated heat flux over the entire image
    SdT(s) = sum(sum(dTdt(:,:,s)));
end

% Plot the heat flux proxy
for s = 25:(135-1)
    set(gcf,'Position',[38.3333   62.3333  989.3333  548.6667])
    set(gcf,'Position',[38.3333   62.3333  450 550])

    surf(mean(dTdt(:,:,[s:1:s+4]),3),'LineStyle','none');
    % Set delay between frames
    dT = 0.1;
    
    % Set color and axis limits and view
    caxis([0,350])
    zlim([-10,0.8e3])
%     xlim([270,450])
%     ylim([0,140])
    view([0,90])

    % Color map
    colormap('jet')
    % Axis
    axis('equal')
   
   % Indicate the start and end of RF
   if 1
       if s >= 122 && s <= 127
               title([num2str(tIR(s),3),' ,End of RF'])
       elseif s >= 32 && s <= 35
               title([num2str(tIR(s),3),' ,Start of RF'])
       else
              title(num2str(tIR(s),3))
       end
   end
   pause(dT)
end
end
%%
% Gather the RF and ECH pulses
Stem = '\MPEX::TOP.'; Branch = 'MACHOPS1:'; RootAddress = [Stem,Branch];
DA{1} = [RootAddress,'PWR_28GHZ']; % Isx
[ECH,t_ech]   = my_mdsvalue_v2(shot,DA(1));
DA{1} = [RootAddress,'RF_FWD_PWR'];
[RF,t_rf]   = my_mdsvalue_v2(shot,DA(1));
tSdT = tIR(1:length(SdT));

close all
% Select the data from the RF tarce less than 5.18 to avoid plotting
% regions that may be full of glitches
rngRF = find(t_rf{1}<=5.18);
figure; hold on
h(1) = plot(t_rf{1}(rngRF),-2.5*RF{1}(rngRF),'r')
h(2) = plot(t_ech{1}(1:end-1),sgolay_t(0.5*ECH{1},3,167),'m')
% Plot the normalized total heat flux to the target. this signal is
% normalized to the peak of the 28 GHz effect. Neglect the spike in heat
% flux that occurs at the start of the RF pulse (t<4.3)
h(3) = plot(tSdT,SdT/max(SdT(find(tSdT>4.29))),'k','LineWidth',2)
legend(h,'RF','28 GHz','Target Heat flux')
set(gcf,'color','w'); box on
ylim([-0.1,1])
xlim([4,5.5])
grid on
title(['Shot: ',num2str(shot)])
set(gcf,'position',[395.6667  377.6667  560.0000  224.0000])

% =========================================================================
% Compare all three stages:
% =========================================================================
MaxCaxis = 450;
% Location of maximum heat flux
[~,b1] = max(SdT(find(tSdT>4.29)));
b1 = b1 + find(tSdT>4.29,1)-2
% Heat flux before ECH pulse
figure;
subplot(1,3,1);
rng1 = [b1-3:-1:b1-7];
surf(mean(dTdt(:,:,rng1),3),'LineStyle','none')
zlim([-10,600])
axis('equal')
xlim([260,440])
ylim([0,150])
caxis([0,MaxCaxis])
view([0,90])
colormap('jet')
title('RF only')
set(gca,'XTick',[],'YTick',[])

% Max heat flux
subplot(1,3,2);
rng2 = [b1:1:b1+3];
surf(mean(dTdt(:,:,rng2),3),'LineStyle','none')
zlim([-10,600])
axis('equal')
xlim([260,440])
ylim([0,150])
caxis([0,MaxCaxis])
view([0,90])
colormap('jet')
title('RF+28GHz, Max')
set(gca,'XTick',[],'YTick',[])

% 3nd ECH pulse heat flux
subplot(1,3,3);
Offset3 = 11;
dOffset3 = 7;
if shot == 20506
    Offset3 = 23;
    dOffset3 = 5;
elseif shot == 20504
    Offset3 = 19;
    dOffset3 = 3;
end
rng3 = [b1+Offset3:1:b1+Offset3+dOffset3];
surf(mean(dTdt(:,:,rng3),3),'LineStyle','none')
zlim([-10,600])
axis('equal')
xlim([260,440])
ylim([0,150])
caxis([0,MaxCaxis])
view([0,90])
colormap('jet')
title('RF+28GHz, steady state')
set(gca,'XTick',[],'YTick',[])
set(gcf,'color','w')


figure(1); 
hold on
plot(tSdT(rng1),SdT(rng1)/max(SdT(find(tSdT>4.29))),'g.-','LineWidth',3)
plot(tSdT(rng2),SdT(rng2)/max(SdT(find(tSdT>4.29))),'g.-','LineWidth',3)
plot(tSdT(rng3),SdT(rng3)/max(SdT(find(tSdT>4.29))),'g.-','LineWidth',3)
% =========================================================================
% Compare 1st and 3rd stages:
% =========================================================================
MaxCaxis = 200;
% Heat flux before ECH pulse
figure;
subplot(1,3,1);
surf(mean(dTdt(:,:,rng1),3),'LineStyle','none')
zlim([-10,600])
axis('equal')
xlim([260,440])
ylim([0,150])
caxis([0,MaxCaxis])
view([0,90])
colormap('jet')
title('RF only')
set(gca,'XTick',[],'YTick',[])

% 3nd ECH pulse heat flux
subplot(1,3,2);
surf(mean(dTdt(:,:,rng3),3),'LineStyle','none')
zlim([-10,600])
axis('equal')
xlim([260,440])
ylim([0,150])
caxis([0,MaxCaxis])
view([0,90])
colormap('jet')
title('RF+28 GHz')
set(gca,'XTick',[],'YTick',[])

% difference between ECH heat flux and helicon only heat flux
subplot(1,3,3);
surf(mean(dTdt(:,:,rng3),3)-mean(dTdt(:,:,rng1),3),'LineStyle','none')
zlim([-10,600])
axis('equal')
xlim([260,440])
ylim([0,150])
caxis([0,MaxCaxis])
view([0,90])
colormap('jet')
title('28 GHz effect only')
set(gca,'XTick',[],'YTick',[])

set(gcf,'color','w')




