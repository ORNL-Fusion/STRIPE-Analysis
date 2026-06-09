% IR data analysis MAY 13th 2016

close all 
clear all
CMPT = 0; % == 1 takes 66 seconds to retrieve data, == 0 takes 1 second

if CMPT == 1
        PATHNAME = [cd,'\'];
        
        shotlist = 8800 + [10 , 11, 12, 13, 15, 17, 18, 19, 20, 21, 23, 25, 27, 29, 31, 32, 33, 34];
        TR2      =        [210,270,290,310,330,350,370,390,410,430,450,470,500,600,230,200,170,140];
        
        shotlist = 8800 + [10 , 11, 12, 13, 15, 17, 18, 19, 20, 21, 23, 25, 27, 29, 31, 32, 33, 34, 83, 84, 85, 87,88];
        TR2      =        [210,270,290,310,330,350,370,390,410,430,450,470,500,600,230,200,170,140,170,140,110, 60, 0];
        [a,b] = sort(TR2);
        
    for s = 1:length(shotlist(b))
        FILENAME{s} = ['shot ',num2str(shotlist(s)),'.seq'];
        videoFileName=[PATHNAME FILENAME{s}];

        % Load the Atlas SDK
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
                irdata{s}(:,:,fr) = im(end:-1:1,:);
                fr = fr + 1;
            end
        end
        % Save only relevant frames in order to keep "Data" to less than 2
        % GB, frame 16 seems to be the hottest for this data set (May 13th)
        Data{s} = irdata{s}(:,:,30:60);
    end
    save('IRData_2016_05_13.mat')
else
    load('IRData_2016_05_13.mat')
end
    
%%
close all
shot = 1;

mov = 1;
if mov
figure(1)
for s = 2:31
   imshow(Data{shot}(:,:,s)- Data{shot}(:,:,1),[]);
   dTdt(:,:,s) = Data{shot}(:,:,s)- Data{shot}(:,:,s-1);
   drawnow;
   title(num2str(s))
   pause(0.05)
end
end

for s = 1:31
   RadialHeatProfile(:,s) = Data{shot}(:,359,s)-Data{shot}(:,359,3); 
end

figure; 
surf(RadialHeatProfile(:,:),'lineStyle','none')

figure; hold on
for s = 1:31
    RPD(:,s) = sgolay_t(RadialHeatProfile(:,s),3,11);
    plot(RPD(:,s))
end
xlim([100,300])


figure; hold on
dt = 19.8543/1000;
t = [1:size(RPD,2)-1]*dt;
plot(t,diff(RPD(218,:)),'k.-')
plot(t,diff(RPD(154,:)),'r.-')


for s = 1:30
    ff(:,s) = RPD(:,s+1)-RPD(:,s);
end

figure;
R = 1:length(ff(:,1));
%contourf(t(20:40),R,ff(:,20:40),20,'lineStyle','none')
%caxis([0,800])
% surf(t(20:40),R,ff(:,20:40),'lineStyle','none')
surf(ff,'lineStyle','none')

colormap('hot')

%%

f0 = 1;
f1 = 16;

rngy = [280:480];
rngx = [100:310];
for s = 1:length(shotlist)
    IR_f10{s} = Data{s}(rngx,rngy,f1)-Data{s}(rngx,rngy,f0);
    IR_f10n{s} = IR_f10{s}./max(max(IR_f10{s}));
end

figure;
h(1) = surf(Data{s}(:,:,f0),'LineStyle','none')
zlim([1.48,1.52]*1e4)
view([0,90])
caxis([1.48,1.54]*1e4)
colormap('jet')

figure;
h(2) = surf(Data{s}(:,:,f1),'LineStyle','none')
zlim([1.48,1.52]*1e4)
caxis([1.48,1.54]*1e4)
view([0,90])
colormap('jet')

% Delta emission
% absolute
figure;
for s = 1:length(shotlist)
subplot(5,5,s); hold on
contourf(IR_f10{b(s)},30,'LineStyle','none')
% zlim([1.3,1.7]*1e4 - 1.3e4)
caxis([0,350])
axis(gca,'equal')
colormap('jet')
title(TR2(b(s)))
%title(FILENAME{b(s)})
end

% normalized
figure;
contourf(IR_f10n{s},80,'LineStyle','none')
zlim([0,1])
caxis([0,0.9])
axis(gca,'equal')
colormap('jet')
title(FILENAME{s})


