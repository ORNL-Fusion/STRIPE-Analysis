% This script reads the data acquired with the RF phase detector installed
% in ProtoMPEX on Sept 16th 2016

% #########################################################################
% THIS SCRIPT USED DURING XP TO MONITOR THE PROGRESS OF THE DATA
% ACQUISITION:
% #########################################################################

close all
clear all

% Btheta scan
shot = 12900 +  [103,108,109,110,111,112 , 113,114 ,115,116 ,117,118,120,121,122,123,124,125,126,127,128,129];
R =             [8.5,9.5,9  ,9.5,10 ,10.5, 11 ,11.5,12 ,12.5,13,13.5,7.5, 7 ,6.5, 6 ,5.5, 5 ,4.5, 4 ,3.5,3  ];

% % Bz scan
% shot = 12900 +  [130,131,132,133,134,135 ,136,137 ,138,157,158];
% R =             [8  ,8.5,9  ,9.5,10 ,10.5,11 ,11.5,7.5,7.5,8.5];
% % 
% shot = 12900 +  [138,153,154,155,158,159,160,161,162,163,164,165,166,168,169,170,171,172,173];
% R =             [7.5,7.5,7.5,7.5,8.5,8  ,7.5,7  ,6.5,6  ,5.5,5  ,4.5,4  ,4.5,5.5,3.5,3  ,2.5];

shot = 13000 + [38,57,60];
R = [7.5,7.5,7.5]
% 
% 180 degree phase shift check
% shot = 12900 + [158,174];
% R =            [8.5,8.5];
% 
% % x = 11 cm comparison
% shot = 12900 + [136,175];
% R =            [11 ,11 ];
% 
% shot = 12900 + [175,176 ,177,178 ,179];
% R =            [11 ,11.5,12 ,12.5,13 ];
% 
% % entire Bz scan:
% shot = 12900 +  [130,131,132,133,134,135 ,136,137 ,138,157,158,159,160,161,162,163,164,165,166,168,169,170,171,172,173,175,176 ,177,178 ,179];
% R =             [8  ,8.5,9  ,9.5,10 ,10.5,11 ,11.5,7.5,7.5,8.5,8  ,7.5,7  ,6.5,6  ,5.5,5  ,4.5,4  ,4.5,5.5,3.5,3  ,2.5,11 ,11.5,12 ,12.5,13 ];
% % 
% TR1 coil 1, Btheta scan:
% shot = 12900 + [180,183,184,185,186,187,188 ,191,192 ,193,194,195,196,197,198,199,200,201,202,204];
R =            [8  ,8  ,8.5,9  ,9.5,10 ,10.5,11 ,11.5,12 ,13 ,14 ,7.5,7  ,6.5,6  ,5.5,4.5,4.5,3.5];
% 
% % Compare Bt field strength of HMJ and non HMJ at center of discharge
% shot = 12900 + [186,205,206];
% R =            [9.5,9.5,9.5];
% 
% % Reduce TR2 = 260 to 200 A
% shot = 12900 + [186,205,206,208];
% R =            [9.5,9.5,9.5,9.5];
% 
% Bt scan at TR2 = 200 A
shot = 12900 + [208,209,210,211,212,213,214,215,216,217,218,219,220];
R =            [9.5,10 ,11 , 12, 13,9  ,8  ,7  ,6  ,5  ,4  ,8.5,7.5];
% 
% % Bt scan at TR2 = 140 A
% shot = 12900 + [221,222,223,224,225,226,227,228,229,230,231,232,233,234 ,235 ,236 ,237,238 ,239,240,241];
% R =            [7.5,8  ,9  ,10 , 11,12 ,13 ,8.5,7  , 6 ,5  ,4  ,3  ,10.5,10.5,11.5,9.5,12.5,6.5,5.5,4.5];
% 
% % Bz on axis TR2 = 140A
% shot = 12900 + [228,242];
% R =            [8.5,8.5];

% % Bz scan TR2 = 140A
% shot = 12900 + [242,243,244,245,246,247,248,249,250,252,253,254];
% R =            [8.5,9  ,10 ,11 ,12 ,13 ,8  ,7  ,6  ,5  ,4  ,3  ];

Data{1} = ['\MPEX::TOP.MACHOPS1:INT_2MM_1']; % Vp0
Data{2} = ['\MPEX::TOP.MACHOPS1:INT_4MM_2']; % Vp90
Data{3} = ['\MPEX::TOP.MACHOPS1:INT_4MM_1']; % Vmag 
Data{4} = ['\MPEX::TOP.MACHOPS1:INT_2MM_2'];


Title{1} = 'Vp0';
Title{2} = 'Vp90';
Title{3} = 'Vmag0';
Title{4} = 'None';

figure;
for ch = 1:4 % For all channels on digitizer
    [f{ch},t{ch}] = my_mdsvalue_v2(shot,Data(ch));
    % f{channel}{shot}
    subplot(2,2,ch);
    hold on
    for s = 1:length(shot)
            rng = find(t{ch}{s}>= 4.14 & t{ch}{s}<= 4.33);
            t{ch}{s} = t{ch}{s}(rng);
            if isempty(rng)
                rng = find(t{ch}{s+1}>= 4.14 & t{ch}{s+1}<= 4.33);
                t{ch}{s} = t{ch}{s+1}(rng);
            end
            f{ch}{s} = sgolay_t(f{ch}{s}(rng),3,71);
            plot(t{ch}{s},f{ch}{s})
    end
    ylim([0,2]);
    title(Title{ch})
end

%%
[VmagRatio,P,tBdot] = BdotPhaseDetect_v1(shot,Data);

%%
close all
tStart = 4.1;
tEnd  = 4.4;
 
C = {'k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:'};

figure; hold on
for s = 1:length(shot)
hA(s) = plot(tBdot{s},VmagRatio{s},C{s},'marker','none','LineStyle','-');
end
xlim([tStart,tEnd])
legend(hA,num2str(shot'))
ylim([0,10])

figure; hold on
for s = 1:length(shot)
hP(s) = plot(tBdot{s},P{s},C{s},'marker','.','LineStyle','none');
end
xlim([tStart,tEnd])
title('Phase')
ylim([-180,180])
legend(hP,num2str(shot'))

figure; hold on
for s = 1:length(shot)
    rng = find(tBdot{s}>=4.25 & tBdot{s}<=4.31);
    if 0
        plot3(tBdot{s}(rng),R(s)*ones(size(P{s}(rng))),unwrap(P{s}(rng)*pi/180)*180/pi,'marker','.','LineStyle','none')
    else
        plot3(tBdot{s}(rng),R(s)*ones(size(P{s}(rng))),P{s}(rng),'marker','.','LineStyle','none')
    end
end
xlim([tStart,tEnd])
ylim([0,15])
% zlim([-10,10])

figure; hold on
for s = 1:length(shot)
    rng = find(tBdot{s}>=4.24 & tBdot{s}<=4.31);
    plot3(tBdot{s}(rng),R(s)*ones(size(VmagRatio{s}(rng))),VmagRatio{s}(rng))
end
xlim([tStart,tEnd])
ylim([0,15])
zlim([0,8])