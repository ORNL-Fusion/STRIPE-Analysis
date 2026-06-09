close all 
clear all

shot(1) = 20837; 
shot(1) = 20822; 
% shot(1) = 20827; 
% shot(1) = 20833; 

a{1} = '\\mpexserver\ProtoMPEX_Data\IR_Camera\2018_04_03';

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
rngX = [150:-1:1]+0;
rngY = [400:-1:270]+30;
graph = 2;
tIR = (1:size(Data{1},3))*(1/98) + 4.17 - 32/98;

for s = 1:(size(Data{1},3)-1)
        f = Data{1}(rngX,rngY,s+1) - Data{1}(rngX,rngY,s);
        dTdt(:,:,s) = f;
        dTdt0(s) = f(75,75);
        SdT(s) = sum(sum(dTdt(:,:,s)));
end

for s = 1:(135-1)
      set(gcf,'Position',[38.3333   62.3333  989.3333  548.6667])

switch graph
    case 1
   surf(Data{1}(rngX,rngY,s) - 1*Data{1}(rngX,rngY,1),'LineStyle','none');
   zlim([0,10e3])
   dT = 0.02;

    case 2

  surf(mean(dTdt(:,:,[s:1:s+5]),3),'LineStyle','none');
   zlim([-10,0.8e3])
  dT = 0.1;
  caxis([0,300])
  colormap('jet')
  end

    
    axis('square')
%     view([20,60])
    view([0,90])
    

   if s >= 122 && s <= 127
           title([num2str(tIR(s),3),' ,End of RF'])
   elseif s >= 32 && s <= 35
           title([num2str(tIR(s),3),' ,Start of RF'])
   else
          title(num2str(tIR(s),3))
   end

      pause(dT)

end
end
%%
Stem = '\MPEX::TOP.'; Branch = 'MACHOPS1:'; RootAddress = [Stem,Branch];
DA{1} = [RootAddress,'PWR_28GHZ']; % Isx
[ECH,t_ech]   = my_mdsvalue_v2(shot,DA(1));
DA{1} = [RootAddress,'RF_FWD_PWR'];
[RF,t_rf]   = my_mdsvalue_v2(shot,DA(1));
%%
close all
tSdT = (1:length(SdT))*(1/98) + 4.17 - 32/98;
rngRF = find(t_rf{1}<=5.18);
figure; hold on
h(1) = plot(t_rf{1}(rngRF),-2.5*RF{1}(rngRF),'r')
% h(2) = plot(t_ech{1}(1:end-1),0.5*ECH{1},'g')
h(2) = plot(t_ech{1}(1:end-1),sgolay_t(0.5*ECH{1},3,167),'m')
h(3) = plot(tSdT,SdT/max(SdT),'k','LineWidth',2)
legend(h,'RF','28 GHz','Target Heat flux')
set(gcf,'color','w'); box on
ylim([-0.1,1])
xlim([4,5.5])
grid on
% figure; surf(Data{1}(rngX,rngY,b1+1),'LineStyle','none')

% =========================================================================
% Compare all three stages:
% =========================================================================
MaxCaxis = 450;
[~,b1] = max(SdT);
% Heat flux before ECH pulse
figure;
subplot(1,3,1);
surf(mean(dTdt(:,:,[b1-3:-1:b1-7]),3),'LineStyle','none')
axis('square')
zlim([-10,600])
xlim([0,130])
ylim([0,150])
caxis([0,MaxCaxis])
view([0,90])
colormap('jet')

% Max heat flux
subplot(1,3,2);
surf(mean(dTdt(:,:,[b1:1:b1+5]),3),'LineStyle','none')
axis('square')
zlim([-10,600])
xlim([0,130])
ylim([0,150])
caxis([0,MaxCaxis])
view([0,90])
colormap('jet')

% 3nd ECH pulse heat flux
subplot(1,3,3);
[~,b1] = find(tSdT>=4.762,1);
surf(mean(dTdt(:,:,[b1:1:b1+9]),3),'LineStyle','none')
axis('square')
zlim([-10,600])
xlim([0,130])
ylim([0,150])
caxis([0,MaxCaxis])
view([0,90])
colormap('jet')

% =========================================================================
% Compare 1st and 3rd stages:
% =========================================================================
MaxCaxis = 150;
[~,b1] = max(SdT);
% Heat flux before ECH pulse
figure;
subplot(1,2,1);
surf(mean(dTdt(:,:,[b1-3:-1:b1-5]),3),'LineStyle','none')
axis('square')
zlim([-10,600])
xlim([0,130])
ylim([0,150])
caxis([0,MaxCaxis])
view([0,90])
colormap('jet')

% 3nd ECH pulse heat flux
subplot(1,2,2);
[~,b1] = find(tSdT>=4.762,1);
[~,b1] = find(tSdT>=4.497,1);
surf(mean(dTdt(:,:,[b1:1:b1+4]),3),'LineStyle','none')
axis('square')
zlim([-10,600])
xlim([0,130])
ylim([0,150])
caxis([0,MaxCaxis])
view([0,90])
colormap('jet')

