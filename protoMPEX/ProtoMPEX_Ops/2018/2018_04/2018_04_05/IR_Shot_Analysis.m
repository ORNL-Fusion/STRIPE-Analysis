close all 
clear all

% shot(1) = 20501;
% shot(1) = 20505;
% shot(1) = 20506; % good core heating, 400 ms ECH pulse
shot(1) = 20504; % very good core heating, 2.5kA
% shot(1) = 20503; % no helicon plasma, 2.5kA
% shot(1) = 20502; % very good core heating, 2.0kA ,but IR is not good near second stage

% shot(1) = 20497; % Good example of core heating
% shot(1) = 20490; % Very small ECH pulse
% shot(1) = 20489; % Helicon only, short pulse
% shot(1) = 20491; % hollow core during ECH

a{1} = '\\mpexserver\ProtoMPEX_Data\IR_Camera\2018_03_26';

% 2018_04_05
shot(1) = 20958; % all probes removed  
a{1} = '\\mpexserver\ProtoMPEX_Data\IR_Camera\2018_04_05';

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

for s = 1:(size(Data{1},3)-1)
    f = Data{1}(:,:,s+1) - Data{1}(:,:,s);
    dTdt(:,:,s) = f;
    dTdt0(s) = f(75,75);
    SdT(s) = sum(sum(dTdt(:,:,s)));
end

for s = 25:(135-1)
    set(gcf,'Position',[38.3333   62.3333  989.3333  548.6667])
    surf(mean(dTdt(:,:,[s:1:s+4]),3),'LineStyle','none');
    zlim([-10,0.8e3])
    dT = 0.1;
    caxis([0,350])
    colormap('jet')
    axis('equal')
    xlim([270,450])
    view([0,90])
   
   if s >= 122 && s <= 127
           title([num2str(tIR(s),3),' ,End of RF'])
   elseif s >= 32 && s <= 35
           title([num2str(tIR(s),3),' ,Start of RF'])
   else
          title(num2str(tIR(s),3))
   end

   [~,bb] = max(SdT); 
   if s == bb
   s
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
tSdT = tIR(1:length(SdT));

close all
rngRF = find(t_rf{1}<=5.18);
figure; hold on
h(1) = plot(t_rf{1}(rngRF),-2.5*RF{1}(rngRF),'r')
h(2) = plot(t_ech{1}(1:end-1),sgolay_t(0.5*ECH{1},3,167),'m')
h(3) = plot(tSdT,SdT/max(SdT),'k','LineWidth',2)
legend(h,'RF','28 GHz','Target Heat flux')
set(gcf,'color','w'); box on
ylim([-0.1,1])
xlim([4,5.5])
grid on
title(['Shot: ',num2str(shot)])

% =========================================================================
% Compare all three stages:
% =========================================================================
MaxCaxis = 450;
% Location of maximum heat flux
[~,b1] = max(SdT);
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

% 3nd ECH pulse heat flux
subplot(1,3,3);
Offset3 = 10;
dOffset3 = 10;
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

figure(1); 
hold on
plot(tSdT(rng1),SdT(rng1)/max(SdT),'g.-','LineWidth',3)
plot(tSdT(rng2),SdT(rng2)/max(SdT),'g.-','LineWidth',3)
plot(tSdT(rng3),SdT(rng3)/max(SdT),'g.-','LineWidth',3)
% =========================================================================
% Compare 1st and 3rd stages:
% =========================================================================
MaxCaxis = 150;
% Heat flux before ECH pulse
figure;
subplot(1,2,1);
surf(mean(dTdt(:,:,rng1),3),'LineStyle','none')
zlim([-10,600])
axis('equal')
xlim([260,440])
ylim([0,150])
caxis([0,MaxCaxis])
view([0,90])
colormap('jet')

% 3nd ECH pulse heat flux
subplot(1,2,2);
surf(mean(dTdt(:,:,rng3),3),'LineStyle','none')
zlim([-10,600])
axis('equal')
xlim([260,440])
ylim([0,150])
caxis([0,MaxCaxis])
view([0,90])
colormap('jet')

%% 2019_01_10
% here we are using the semi-infinite slab solution with surface heat flux
% to figure out what the heat flux is
% this data from the upper hybrid experiments had front imaging

T0 = Data{1}(:,:,24);
for s = 24:size(Data{1},3)
     dTs(:,:,s-23) = (Data{1}(:,:,s) - T0)*1e-3;
     dT2(:,:,s-23) = dTs(:,:,s-23).^2;
end

dt = 1/60;

for s = 1:size(dT2,3)-1
    ddT2(:,:,s) = (dT2(:,:,s+1) - dT2(:,:,s))*dt;
    SqddT2(:,:,s) = sqrt(ddT2(:,:,s));
end
%%
 figure; 
 for s  = 1:100
%surf(ddT2(:,:,3),'LineStyle','none')
set(gcf,'Position',[38.3333   62.3333  989.3333  548.6667])
surf(real(SqddT2(1:150,250:400,s)),'LineStyle','none')
view([0,90])
% view([-20,60])

%     xlim([270,450])
%     ylim([270,450])

%     zlim([0,0.8])
    axis('equal')
    caxis([0,0.5])
    colormap('jet')
    pause(0.1)
 end
%%
for s = 25:(135-1)
    set(gcf,'Position',[38.3333   62.3333  989.3333  548.6667])
    surf(mean(dTdt(:,:,[s:1:s+4]),3),'LineStyle','none');
    zlim([-10,0.8e3])
    dT = 0.1;
    caxis([0,350])
    colormap('jet')
    axis('equal')
    xlim([270,450])
    view([0,90])
   
   if s >= 122 && s <= 127
           title([num2str(tIR(s),3),' ,End of RF'])
   elseif s >= 32 && s <= 35
           title([num2str(tIR(s),3),' ,Start of RF'])
   else
          title(num2str(tIR(s),3))
   end

   [~,bb] = max(SdT); 
   if s == bb
   s
   end
   pause(dT)
end
