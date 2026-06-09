close all 
clear all

% RF power scan
% DLP 10.5 is on-axis
shot(1) = 20998; 
shot(2) = 21002; 

% % RF power scan with ECH, no probes, clean IR
shot(2) = 21003; % 2.8 v
shot(1) = 21006; % 4.0 v

% RF power scan with ECH, no probes
% see the effect of the UHR layer
% shot(1) = 21004; % 3.0 v
% shot(1) = 21005; % 3.5 v
% shot(1) = 21006; % 4.0 v
% shot(2) = 21007; % 4.5 v
% shot(1) = 21008; % 5.0 v
% shot(2) = 21009; % 5.5 v
% shot(2) = 21010; % 6.0 v
% shot(1) = 21011; % 6.5 v
% shot(1) = 21012; % 7.5 v
shot(2) = 21013; % 8.5 v
% shot(1) = 21014; % 9.0 v
shot(1) = 21015; % 9.96 v

% % PS1 field scan:
% shot(1) = 21016; % 2.0 kA
% shot(1) = 21017; % 2.5 kA
% shot(1) = 21018; % 3.0 kA
% shot(1) = 21019; % 3.5 kA
% shot(1) = 21020; % 4.0 kA
% shot(1) = 21021; % 4.5 kA
% 
% % DLP 1.5 radial scan with ECH
% shot(1) = 21022; %
% % ECH power scan , clean IR
% shot(1) = 21041; %Reference
% % shot(1) = 21042; % 68%
% % shot(1) = 21043; % 33%
% % shot(1) = 21046; %
shot(2) = 21046;
shot(1) = 21047;
% shot(1) = 21049;
% 
a{1} = '\\mpexserver\ProtoMPEX_Data\IR_Camera\2018_04_06';
a{2} = '\\mpexserver\ProtoMPEX_Data\IR_Camera\2018_04_06';

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

% Calculate a quantity that approximates the heat flux at the plate:
for n = 1:length(shot)
    for s = 1:(size(Data{n},3)-1)
        % Time derivative:
        tIR{n} = (1:size(Data{n},3))*(1/98) + 4.17 - 32/98;
        f{n} = Data{n}(:,:,s+1) - Data{n}(:,:,s);
        dTdt{n}(:,:,s) = f{n};
        % Integrated heat flux over the entire image
        SdT{n}(s) = sum(sum(dTdt{n}(:,:,s)));
    end
end
tSdT{1} = tIR{1}(1:length(SdT{1}));
tSdT{2} = tIR{2}(1:length(SdT{2}));

% Plot the heat flux proxy
for s = 25:(135-1)
    set(gcf,'Position',[38.3333   62.3333  989.3333  548.6667])
    set(gcf,'Position',[38.3333   62.3333  450 550])
    
    for n = 1:length(shot)
    subplot(1,2,n);
    surf(mean(dTdt{n}(:,:,[s:1:s+4]),3),'LineStyle','none');
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
               title([num2str(tIR{n}(s),3),' ,End of RF'])
       elseif s >= 32 && s <= 35
               title([num2str(tIR{n}(s),3),' ,Start of RF'])
       else
              title(num2str(tIR{n}(s),3))
       end
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

close all
% Select the data from the RF tarce less than 5.18 to avoid plotting
% regions that may be full of glitches
rngRF = find(t_rf{1}<=5.18);
figure;
subplot(2,1,1)
hold on
h(1) = plot(t_rf{1}(rngRF),-2.5*RF{1}(rngRF),'r','LineWidth',2)
h(2) = plot(t_rf{2}(rngRF),-2.5*RF{2}(rngRF),'g','LineWidth',2)
h(3) = plot(t_ech{1}(1:end-1),sgolay_t(0.5*ECH{1},3,167),'r','LineWidth',2)
h(4) = plot(t_ech{2}(1:end-1),sgolay_t(0.5*ECH{2},3,167),'g','LineWidth',2)
% legend(h,'RF','28 GHz','Target Heat flux')
set(gcf,'color','w'); box on
ylim([-0.1,1])
xlim([4,5.5])
grid on
title(['Shot: ',num2str(shot)])

subplot(2,1,2)
hold on
% Plot the normalized total heat flux to the target. this signal is
% normalized to the peak of the 28 GHz effect. Neglect the spike in heat
% flux that occurs at the start of the RF pulse (t<4.3)
h(4) = plot(tSdT{1},SdT{1}/max(SdT{1}(find(tSdT{1}>4.29))),'r','LineWidth',3)
h(5) = plot(tSdT{2},SdT{2}/max(SdT{1}(find(tSdT{1}>4.29))),'g','LineWidth',3)
legend(h,{num2str(shot(1)),num2str(shot(2))})
set(gcf,'color','w'); box on
ylim([-0.1,1])
xlim([4,5.5])
grid on

% =========================================================================
% Compare all three stages:
% =========================================================================
MaxCaxis = 450;

% Location of maximum heat flux
[~,b1] = max(SdT{1}(find(tSdT{1}>4.29)));
b1 = b1 + find(tSdT{1}>4.29,1)-2

for n = 1:length(shot)
% Heat flux before ECH pulse
figure;
subplot(1,3,1);
rng1 = [b1-3:-1:b1-7];
surf(mean(dTdt{n}(:,:,rng1),3),'LineStyle','none')
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
surf(mean(dTdt{n}(:,:,rng2),3),'LineStyle','none')
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
rng3 = [b1+Offset3:1:b1+Offset3+dOffset3];
surf(mean(dTdt{n}(:,:,rng3),3),'LineStyle','none')
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
plot(tSdT{n}(rng1),SdT{n}(rng1)/max(SdT{1}(find(tSdT{1}>4.29))),'k.-','LineWidth',3)
plot(tSdT{n}(rng2),SdT{n}(rng2)/max(SdT{1}(find(tSdT{1}>4.29))),'k.-','LineWidth',3)
plot(tSdT{n}(rng3),SdT{n}(rng3)/max(SdT{1}(find(tSdT{1}>4.29))),'k.-','LineWidth',3)
end
% =========================================================================
% Compare 1st and 3rd stages:
% =========================================================================
MaxCaxis = 150;
% Heat flux before ECH pulse
for n = 1:length(shot)
figure;
subplot(1,3,1);
surf(mean(dTdt{n}(:,:,rng1),3),'LineStyle','none')
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
surf(mean(dTdt{n}(:,:,rng3),3),'LineStyle','none')
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
surf(mean(dTdt{n}(:,:,rng3),3)-mean(dTdt{n}(:,:,rng1),3),'LineStyle','none')
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
end



