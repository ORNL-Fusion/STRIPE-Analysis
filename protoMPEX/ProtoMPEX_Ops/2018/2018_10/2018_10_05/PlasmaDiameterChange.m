% Fast visible light camera measurements of the plasma-target region in
% Proto-MPEX

clear all
close all

FileAddress = 'C:\Users\nfc\Documents\ProtoMPEX_Ops\2018\2018_10\2018_10_05';

FileName{1} = 'slomo_1538770650.mov';
fctr(1) = 1;

% Load movie data:
for n = 1;
    
vidObj{n} = VideoReader([FileAddress,'\',FileName{n}]);

vidHeight = vidObj{n}.Height;
vidWidth = vidObj{n}.Width;
vidObj{n}.CurrentTime = 2;

s{n} = struct('cdata',zeros(vidHeight,vidWidth,3,'uint8'),...
        'colormap',[]);

k = 1;
Time(k) = 0;
while hasFrame(vidObj{n})
    s{n}(k).cdata = readFrame(vidObj{n});
    Dalpha{n}{k} = s{n}(k).cdata(:,:,1);
    xn = 550;
    Da_R{n}(k,:) = fctr(n)*s{n}(k).cdata(:,xn,1);
    k = k+1;   
    Time(k) = Time(k-1) + 1e-3; % time in seconds
end
end

Time = Time(1:end-1);

%%
n= 1; k = 140;
figure;
for k = 60:10:450;
F = double(s{n}(k).cdata(:,:,3));
surf(F/max(max(F)),'LineStyle','none')
view([0,90])
colormap('hot')
zlim([0,1])
caxis([0,1])
title(num2str(k))
drawnow
% pause(0.001)
end

%%
close all
figure; 
subplot(1,2,1)
k = 90;
F = double(s{n}(k).cdata(:,:,3));
surf(F/max(max(F)),'LineStyle','none')
view([0,90])
colormap('hot')
zlim([0,1])
caxis([0,1])
title(num2str(k))
xlim([0,500])

subplot(1,2,2)
k = 300;
F = double(s{n}(k).cdata(:,:,3));
surf(F/max(max(F)),'LineStyle','none')
view([0,90])
colormap('hot')
zlim([0,1])
caxis([0,1])
xlim([0,500])

title(num2str(k))