clear all
close all

FileName{1} = 'slomo_c_13962';
FileName{2} = 'slomo_bw_13962'
FileName{3} = 'slomo_c_13966';
FileName{4} = 'slomo_bw_13966'

Address = 'C:\Users\nfc\Documents\Proto-MPEX Data Analysis\2017\04 - April\13th - Long RF pulse + ECH\Visible camera selected\';

LoadFile = 1;
if LoadFile == 1 
    for s = 1:length(FileName)
        f{s} = VideoReader([Address,FileName{s},'.mov']);
    end
    
    Stem = '\MPEX::TOP.'; Branch = 'MACHOPS1:'; RootAddress = [Stem,Branch];
    DA{1} = [RootAddress,'PWR_28GHZ']; % Isx
    [ECH,t_ech]   = my_mdsvalue_v2(13966,DA(1))
    DA{1} = [RootAddress,'RF_FWD_PWR'];
    [RF,t_rf]   = my_mdsvalue_v2(13966,DA(1));
    save('v2'); 
else
    load('v2');
end

for k = 1:length(FileName)
    vidHeight{k} = f{k}.Height;
	vidWidth{k} = f{k}.Width;
    S{k} = struct('cdata',zeros(vidHeight{k},vidWidth{k},3,'uint8'),...
    'colormap',[]);
    n = 1;
    while hasFrame(f{k})
        S{k}(n).cdata = readFrame(f{k});
        n = n+1;
    end
    
end

%% Play movie
close all
k = 1;
figure
set(gcf,'position',[150 150 0.7*vidWidth{k} 0.7*vidHeight{k}]);
set(gca,'units','pixels');
set(gca,'position',[0 0 0.5*vidWidth{k} 0.5*vidHeight{k}]);
movie(S{k},1,180)

%%
VidFrameOffset = 163 - 4165;
t = [[1:1:500]-VidFrameOffset]; % 1 ms per frame, since frame rate is 500 fr/s and 500 ms movie
for k = 1:length(FileName)
    for s = 1:1:500
           R{k}(:,:,s) = S{k}(s).cdata(:,:,1)/3;
           G{k}(:,:,s) = S{k}(s).cdata(:,:,2)/3;
           B{k}(:,:,s) = S{k}(s).cdata(:,:,3)/3;
           fr{k}(:,:,s) = (R{k}(:,:,s)+ 1*G{k}(:,:,s) + 1*B{k}(:,:,s));
       IntensityIntegrated{k}(s) = sum(sum(fr{k}(:,:,s),1),2);
       rp{k}(:,:,s) = fr{k}(1:end,380:400,s); 
       RadialEmission{k}(:,s) = mean(rp{k}(:,:,s),2);
       OnAxisEmission{k}(:,s) = mean(rp{k}(280:320,:,s),2);
    end
end

%%
close all
figure; hold on
for k = 1:length(FileName)
h(k) = plot(t*1e-3,IntensityIntegrated{k})
end
plot(t_rf{1}(1:length(RF{1})),(RF{1}.^2)*1e7)
plot(t_ech{1}(1:length(ECH{1})),(ECH{1})*1e7)

set(h(1),'color','bl','LineWidth',2)
set(h(2),'color','bl')
set(h(3),'color','r','LineWidth',2)
set(h(4),'color','r')
legend(h,'c no ech','bw no ech','c ech','bw ech')
ylim([0,1e8])
% It looks like the color 

address{3} = '\MPEX::TOP.MACHOPS1:PG3'; % PG6.5
[PG3,t3] = my_mdsvalue_v2(13966,address(3));
[PG3_0,t3_0] = my_mdsvalue_v2(13968,address(3));

plot(t3{1}(1:length(PG3{1})),(PG3{1}-PG3_0{1})*2e7)
xlim([4,4.8])
% nta = find(t*1e-3 >= 4.25); nta = nta(1)
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

k = 3;
figure
hold on;
nf = 150:10:450;
for n = 1:length(nf);
    fdata = RadialEmission{k}(:,nf(n));
plot3(1:length(fdata),nf(n)*ones(size(fdata)),fdata)
end

