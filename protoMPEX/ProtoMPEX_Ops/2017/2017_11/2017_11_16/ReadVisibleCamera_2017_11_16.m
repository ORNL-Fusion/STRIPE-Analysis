% 2017_11_16
clear all
close all

shotlist = [17773,17827];

% the index "s" is used if you have several shots
for s = 1:length(shotlist)
    FileName{s} = ['slomo_',num2str(shotlist(s))];
    f{s} = VideoReader([FileName{s},'.mov']);
    % f{s} is video object in matlab
    vidHeight = f{s}.Height;
    vidWidth  = f{s}.Width;
    % Create a structure to store the video frames
    S{s} = struct('cdata',zeros(vidHeight,vidWidth,3,'uint8'),...
    'colormap',[]);
    n = 1; % Initialize the frame index
    while hasFrame(f{s})
        S{s}(n).cdata = readFrame(f{s});
        n = n+1;
    end
end

figure
s = 1; % shot 1
n = 800; % look at frame 800, there are 1500 frames
image(S{s}(n).cdata)

figure
s = 2; % shot 2
n = 800; % look at frame 800, there are 1500 frames
image(S{s}(n).cdata)