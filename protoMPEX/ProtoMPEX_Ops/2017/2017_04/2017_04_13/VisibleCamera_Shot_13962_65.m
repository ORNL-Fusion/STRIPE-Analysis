% This was the first instance of the analysis created on May 15th 2017

clear all
close all

FileName{1} = 'slomo_c_13962';
FileName{2} = 'slomo_bw_13962'
FileName{3} = 'slomo_c_13965';
FileName{4} = 'slomo_bw_13965'

Address = 'C:\Users\nfc\Documents\Proto-MPEX Data Analysis\2017\04 - April\13th - Long RF pulse + ECH\Visible camera selected\';
% Load video data as an array

LoadFile = 0;
if LoadFile == 1 
    for s = 1:length(FileName)
        f{s} = VideoReader([Address,FileName{s},'.mov']);
        d{s} = read(f{s}); % :,:,3,500; 3 colors RGB and 500 frames
    end
    
    Stem = '\MPEX::TOP.'; Branch = 'MACHOPS1:'; RootAddress = [Stem,Branch];
    % DA{1} = [RootAddress,'PWR_28GHZ']; % Isx
    % [ECH,t_ech]   = my_mdsvalue_v2(shotlist,DA(1))
    DA{1} = [RootAddress,'RF_FWD_PWR'];
    [RF,t_rf]   = my_mdsvalue_v2(13970,DA(1));
    save('v'); 
else
    load('v');
     for s = 1:length(FileName)
        d{s} = read(f{s}); % :,:,3,500; 3 colors RGB and 500 frames
     end
    clear rp RadialEmission fr
end

%%
close all
figure; hold on
t = [[1:1:500]-163];
for k = 1:length(FileName)
    for s = 1:1:500
           R{k}(:,:,s) = d{k}(:,:,1,s)/3;
           G{k}(:,:,s) = d{k}(:,:,2,s)/3;
           B{k}(:,:,s) = d{k}(:,:,3,s)/3;
           fr{k}(:,:,s) = (R{k}(:,:,s)+ 1*G{k}(:,:,s) + 1*B{k}(:,:,s));
       c{k}(s) = sum(sum(fr{k}(:,:,s),1),2);
       rp{k}(:,:,s) = fr{k}(1:end,380:400,s); 
       RadialEmission{k}(:,s) = mean(rp{k}(:,:,s),2);
    end
h(k) = plot(t,c{k})
end
plot(t_rf{1}(1:length(RF{1}))*1e3 - 4162,(RF{1}.^2)*1e7)
set(h(1),'color','bl','LineWidth',2)
set(h(2),'color','bl')
set(h(3),'color','r','LineWidth',2)
set(h(4),'color','r')
legend(h,'c no ech','bw no ech','c ech','bw ech')
ylim([0,1e8])
% It looks like the color 
%%
figure;
for k = 1:4
subplot(2,2,k); hold on
surf(RadialEmission{k},'LineStyle','none')
caxis([0,105])
zlim([0,255])
end

figure;
for k = 1:4
subplot(2,1,1); hold on
plot(mean(RadialEmission{1},1))
plot(mean(RadialEmission{3},1))
ylim([0,100])

subplot(2,1,2); hold on
plot(mean(RadialEmission{2},1))
plot(mean(RadialEmission{4},1))
ylim([0,100])
end

k = 1;
figure
hold on;
nf = 150:10:450;
for n = 1:length(nf);
    fdata = RadialEmission{k}(:,nf(n));
plot3(1:length(fdata),nf(n)*ones(size(fdata)),fdata)
end

%%
mov = 1;
k = 3;
r = 1;
if mov
    figure(1)
    for s = 150:2:500
       imshow(fr{k}(:,:,s),[]);
       drawnow;
       title(num2str(s))
       pause(0.0001)
    end
end
