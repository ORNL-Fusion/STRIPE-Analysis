% 2018_04_17
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

shotlist = 21000 + [742,743,744,745];
dt       =         [0  ,100,200,300];
Shutter  =         [15 ,15 ,15 ,20 ];
fctr     =         [1  ,1  ,1  ,1.3];

shotlist = 21000 + [746,748];
dt       =         [3.5,5.0];
Shutter  =         [15 ,15 ];
fctr     =         [1  ,1  ];

% =========================================================================
shotlist = 21000 + [742,831,832,833];
flow     =         [3.5,3.5,3.5,3.5];
Shutter  =         [15 ,15 ,15 ,15 ];
fctr     =         [1  ,1  ,1  ,1  ];
CCTurbo  =         [1  ,0  ,0  ,1  ];

shotlist = 21000 + [832,833];
flow     =         [3.5,3.5];
Shutter  =         [15 ,15 ];
fctr     =         [1  ,1  ];
CCTurbo  =         [0  ,1  ];
% =========================================================================

for n = 1:length(shotlist)
    vidObj{n} = VideoReader(mov(n).name);

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
    n = 3;
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




