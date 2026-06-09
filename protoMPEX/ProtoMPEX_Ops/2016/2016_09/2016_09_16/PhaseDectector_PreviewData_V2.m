% This script reads the data acquired with the RF phase detector installed
% in ProtoMPEX on Sept 16th 2016

% #########################################################################
% THIS SCRIPT USED DURING XP TO MONITOR THE PROGRESS OF THE DATA
% ACQUISITION:
% #########################################################################

close all
clear all

%shot = 1e4 + 300 + [95,92,88,86,85,82,81,80,79,99,105,101,106,108,109,111,113:118];


% Sequence for 
shot = 1e4 + 300 + ... 
    [95   92 88   86 85   82 81  80  79  99  105 100 101 103 104 106 107 108 109 110 111 112 113 114 115 116 117 118];
R = [12.5 12 11.5 11 10.5 10 9.5 9.0 8.5 8.0 7.5 7.5 7.0 7.0 7.0 6.5 6.5 6.0 5.5 5.5 5.0 5.0 4.5 4.0 3.5 3.0 2.5 2.0] ;
Rcenter = 7.5;

Data{1} = ['\MPEX::TOP.MACHOPS1:INT_4MM_1'];
Data{2} = ['\MPEX::TOP.MACHOPS1:INT_4MM_2'];
Data{3} = ['\MPEX::TOP.MACHOPS1:INT_2MM_1'];
Data{4} = ['\MPEX::TOP.MACHOPS1:INT_2MM_2'];

Title{1} = 'Vp0';
Title{2} = 'Vp90';
Title{3} = 'Vmag0';
Title{4} = 'None';

[VmagRatio,P,t] = BdotPhaseDetect_v1(shot,Data);
%%
tStart = 4.29;
tEnd  = 4.30;
close all

for s = 1:length(shot)
    Pa(:,s) = unwrap(P{s}*pi/180)*180/pi;
end
Pa = unwrap(Pa*pi/180,[],1)*180/pi;

figure; hold on
for s = 1:length(shot)
    hM(s) = plot3(t{s},R(s)*ones(size(VmagRatio{s})),sgolay_t(VmagRatio{s},3,21));
end
%legend(hM,num2str(shot'))
xlim([tStart,tEnd])
ylim([0,15])
zlim([0,8])
view([30,30])

Offset = zeros(length(shot));

Offset(find(R == 2.0)) = -2*360;
Offset(find(R == 2.5)) = -2*360;
Offset(find(R == 3.0)) = -1*360;
Offset(find(R == 3.5)) = -1*360;
Offset(find(R == 4.0)) = -1*360;

Offset(find(R == 4.5)) = -1*360;
Offset(find(R == 5.0)) = -360;
Offset([12]) = +360; % R = 7.0
Offset([13]) = +360; % R = 7.5

figure; 
hold on
for s = 1:length(shot)
    F = unwrap(P{s}*pi/180)*180/pi;
    rng = find(t{s}>=4.29 & t{s}<=4.30);
    b(s) = mean(F(rng))/360;
    b(s) = round(b(s));
    plot3(t{s},R(s)*ones(size(F)),(sgolay_t(F,3,5) - b(s)*360 + Offset(s)))
end
xlim([tStart,tEnd])
zlim([-660,360])
ylim([0,15])
view([60,30])