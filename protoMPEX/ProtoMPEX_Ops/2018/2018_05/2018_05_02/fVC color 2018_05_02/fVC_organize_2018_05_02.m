% 2018_05_02
% The objective of this script is to organize the fVC files and name them
% accordingly

clear all 
close all
PlayMovie = 1;
PlaySurf  = 1;

% Collect the name of all the .mov files
mov = dir('*.mov');
% Collect the name of all the .txt files
movTxt = dir('*.txt');

ShotStart = 21838;

% =========================================================================
% Central Chamber pump ON and OFF, optimized fueling
shotlist = 21000 + [845,846,849];
flow     =         [3.5,3.5,3.5];
Shutter  =         [15 ,15 ,15 ];
fctr     =         [1  ,1  ,1  ];
CCTurbo  =         [1  ,0  ,0  ];
% =========================================================================
% Central Chamber pump ON and OFF, optimized fueling
shotlist = 21000 + [850,851,852,853];
flow     =         [3.5,3.5,3.5,3.5];
Shutter  =         [15 ,15 ,15 ,15 ];
fctr     =         [1  ,1  ,1  ,1  ];
CCTurbo  =         [0  ,0  ,0  ,0  ];
% =========================================================================
% effect of gas timing on pumping 
shotlist = 21000 + [845,854,855,857];
% =========================================================================
% Early vrs optimized puffing with ALL pumps ON
shotlist = 21000 + [845,857];
% =========================================================================
% Early vrs optimized puffing, ALL pumps and Ballast turbo only
shotlist = 21000 + [845,857,855];

for n = 1:length(shotlist)
    q = shotlist(n) - ShotStart;
    
    vidObj{n} = VideoReader(mov(q).name);

    vidHeight = vidObj{n}.Height;
    vidWidth = vidObj{n}.Width;
    vidObj{n}.CurrentTime = 2;

    s{n} = struct('cdata',zeros(vidHeight,vidWidth,3,'uint8'),...
        'colormap',[]);

    k = 1;
    while hasFrame(vidObj{n})
        s{n}(k).cdata = readFrame(vidObj{n});
        Dalpha{n}{k} = s{n}(k).cdata(:,:,1);
        xn = 550;
        Da_R{n}(k,:) = fctr(n)*s{n}(k).cdata(:,xn,1);
        k = k+1;
%         
    end
end

%%
figure
for n = 1:4
subplot(2,2,n)
contourf(Da_R{n}',10,'linestyle','none'); xlim([0,700])
title(num2str(shotlist(n)))
end
%%
if PlaySurf
    n = 1;
    figure
    for k = 60:5:450%length(Dalpha)-60
    surf(Dalpha{n}{k},'lineStyle','none')
    view([0,90])
%     view([60,30])
    caxis([0,255])
    zlim([0,255])
    title(num2str(k))
    pause(0.001)
    end    
end

%%
if PlaySurf
    n = 1;
    figure; hold on
    for k = 1:5:450%length(Dalpha)-60
    plot(Da_R{n}(k,:))
    ylim([0,255])
    xlim([0,700])
    title(num2str(k))
    pause(0.001)
    end    
end
%%
if 1
    n = 1;
    figure
    xn = 250;
    xn = 555;
    for k = 50:10:580
%         plot3(110*ones(size(Dalpha{110}(:,xn),1)),k*ones(size(Dalpha{110}(:,xn)),1),Dalpha{110}(:,xn))
        hold on
        plot(Dalpha{n}{k}(:,xn))
%         plot3(k*ones(size(Dalpha{n}{k}(:,xn),1)),1:size(Dalpha{n}{k}(:,xn),1),Dalpha{n}{k}(:,xn))
% view([-25,40])
%         hold off
        ylim([0,255])
        xlim([1,700])
        pause(0.1)
    end
end
%%
if PlayMovie
    figure
    n = 1; 
currAxes = axes;
while hasFrame(vidObj{n})
    vidFrame = readFrame(vidObj{n});
    image(vidFrame, 'Parent', currAxes);
    currAxes.Visible = 'off';
    pause(1./vidObj{n}.FrameRate);
end
end




