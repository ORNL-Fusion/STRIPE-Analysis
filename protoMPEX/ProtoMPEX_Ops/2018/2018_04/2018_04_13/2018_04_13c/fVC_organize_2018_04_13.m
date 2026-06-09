% 2018_04_17
% The objective of this script is to organize the fVC files and name them
% accordingly

clear all 
close all
PlayMovie = 0;
PlaySurf  = 0;

% Collect the name of all the .mov files
mov = dir('*.mov');
% Collect the name of all the .txt files
movTxt = dir('*.txt');


% The size of mov is 89 elements which is the same number of shots taken
% during the 2018_04_13 experiments

ShotStart = 21187;
% ShotStart = 21238;
% ShotStart = 21241;

ShotEnd   = ShotStart + 88;
Shotlist  = ShotStart:1:ShotEnd;

% Lets us load data
% shot 21247, 5v
% shot 21246, 6v
% shot 21245, 7v
% shot 21241, 8v
% shot 21239, 9v
% shot 21238, 10v

vidObj = VideoReader(mov(Shotlist == 21247).name);

vidHeight = vidObj.Height;
vidWidth = vidObj.Width;
vidObj.CurrentTime = 2;

s = struct('cdata',zeros(vidHeight,vidWidth,3,'uint8'),...
    'colormap',[]);

k = 1;
while hasFrame(vidObj)
    s(k).cdata = readFrame(vidObj);
    Dalpha{k} = s(k).cdata(:,:,1);
    k = k+1;
end
%%
if PlaySurf
    figure
    for k = 80:5:405%length(Dalpha)-60
        % plasma turn off: 
        % plasma turn on : 80:1:105
    surf(Dalpha{k},'lineStyle','none')
    view([0,90])
    caxis([0,255])
    zlim([0,255])
    title(num2str(k))
    pause(0.01)
    end    
end
%%
    figure
    xn = 250;
    xn = 555;
    for k = 80:10:580
%         plot3(110*ones(size(Dalpha{110}(:,xn),1)),k*ones(size(Dalpha{110}(:,xn)),1),Dalpha{110}(:,xn))
        hold on
%         plot(Dalpha{k}(:,xn))
        plot3(k*ones(size(Dalpha{k}(:,xn),1)),1:size(Dalpha{k}(:,xn),1),Dalpha{k}(:,xn))
view([-25,40])
%         hold off
        zlim([0,255])
        xlim([1,500])
        pause(0.0001)
    end
%%
if PlayMovie
    figure
currAxes = axes;
while hasFrame(vidObj)
    vidFrame = readFrame(vidObj);
    image(vidFrame, 'Parent', currAxes);
    currAxes.Visible = 'off';
    pause(1/vidObj.FrameRate);
end
end




