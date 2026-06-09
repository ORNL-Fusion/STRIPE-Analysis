% 2017_04_10, ICH Bdot probe measurements, JF Caneses 2018_04_10

% #########################################################################    
% READ ME!!!
% Pawel, since you do not have all the codes to run the "BdotPhaseDetect_v3.m"
% code, just set LoadCalculatedData = 1, GetRawData = 0 and SaveData = 0.
% #########################################################################    

close all
clear all

LoadCalculatedData = 1; % Load previously calculated data

% =========================================================================
% If you want to recalculate your own data enable the follwing variables to
% 1, otherwise just load previously calculate data LoadCalculatedData = 1
GetRawData = 0; % Get raw data
SaveData = 0;   % Enable saving new calculated data
% =========================================================================


if GetRawData
% #########################################################################    
% 180 degree phase test
% shot = 21000 + [82,83];
% #########################################################################    
% Radial scan
shot = [21080:1:21082,21084:1:21085,21087:1:21091];
R     = [-4:0.5:-2,-1:0.5:1];
Atten = [10*ones(3,1)',zeros(7,1)'];
Coil  = (1:1:12);
% #########################################################################    

% coil S1: 9.5 , B_phi_z
Data{1} = ['\MPEX::TOP.MACHOPS1:COIL1_1']; % Vp0
Data{2} = ['\MPEX::TOP.MACHOPS1:COIL1_2']; % Vp90
Data{3} = ['\MPEX::TOP.MACHOPS1:COIL1_3']; % Vmag
Title{1} = 'S1,Vp0';
Title{2} = 'S1,Vp90';
Title{3} = 'S1,Vmag0';

% coil S2: 9.5 , Br
Data{4} = ['\MPEX::TOP.MACHOPS1:COIL1_4']; % Vp0
Data{5} = ['\MPEX::TOP.MACHOPS1:COIL1_5']; % Vp90
Data{6} = ['\MPEX::TOP.MACHOPS1:COIL1_6']; % Vmag
Title{4} = 'S2,Vp0';
Title{5} = 'S2,Vp90';
Title{6} = 'S2,Vmag0';

% coil S3: 10.5, B_phi_z
Data{7} = ['\MPEX::TOP.MACHOPS1:COIL1_7']; % Vp0
Data{8} = ['\MPEX::TOP.MACHOPS1:COIL1_8']; % Vp90
Data{9} = ['\MPEX::TOP.MACHOPS1:COIL1_9']; % Vmag
Title{7} = 'S3,Vp0';
Title{8} = 'S3,Vp90';
Title{9} = 'S3,Vmag0';

% coil S4: 10.5, Br
Data{10} = ['\MPEX::TOP.MACHOPS1:COIL1_10']; % Vp0
Data{11} = ['\MPEX::TOP.MACHOPS1:COIL1_11']; % Vp90
Data{12} = ['\MPEX::TOP.MACHOPS1:COIL1_12']; % Vmag
Title{10} = 'S4,Vp0';
Title{11} = 'S4,Vp90';
Title{12} = 'S4,Vmag0';

Stem = '\MPEX::TOP.'; Branch = 'MACHOPS1:'; RootAddress = [Stem,Branch];

figure;
    for ch = 1:12 % For all channels on digitizer
        [f{ch},t{ch}] = my_mdsvalue_v2(shot,Data(ch));
        % f{channel}{shot}
        subplot(4,3,ch);
        hold on
        for s = 1:length(shot)
                rng = find(t{ch}{s}>= 4.1 & t{ch}{s}<= 4.7);
                t{ch}{s} = t{ch}{s}(rng);
                if isempty(rng)
                    rng = find(t{ch}{s+1}>= 4.1 & t{ch}{s+1}<= 4.7);
                    t{ch}{s} = t{ch}{s+1}(rng);
                end
                f{ch}{s} = f{ch}{s}(rng);%sgolay_t(f{ch}{s}(rng),3,11);
                h(ch,s) = plot(t{ch}{s},f{ch}{s});
        end
        ylim([0,2]);
        title(Title{ch})
        xlim([4,4.8])
    end

[VmagRatio_S1,P_S1,tBdot_S1] = BdotPhaseDetect_v3(shot,Data(1:3));
[VmagRatio_S2,P_S2,tBdot_S2] = BdotPhaseDetect_v3(shot,Data(4:6));
[VmagRatio_S3,P_S3,tBdot_S3] = BdotPhaseDetect_v3(shot,Data(7:9));
[VmagRatio_S4,P_S4,tBdot_S4] = BdotPhaseDetect_v3(shot,Data(10:12));

    if SaveData
    error('Dont recompute the data, talk to Juan')
    save('ICH_Bdot_2017_04_10')
    end
end

if LoadCalculatedData
    load('ICH_Bdot_2017_04_10')
end

% return
%%
tStart = 4.1;
tEnd   = 4.7;
close all
% #########################################################################
% Amplitude data
% The effect of the Attenuators has not been added (JF Caneses 2014_04_13)
% #########################################################################
figure
N = ceil(sqrt(length(shot)));
hold on
for p = 1:length(shot)
plot3(R(p)*ones(size(VmagRatio_S1{p})),tBdot_S1{p}(1:length(VmagRatio_S1{p})),VmagRatio_S1{p})
ylim([tStart,tEnd])
zlim([0,7])
grid on
view([-60,40])
end
title(['S_1, Vmag'])

figure
hold on
for p = 1:length(shot)
h2(p) = plot3(R(p)*ones(size(VmagRatio_S2{p})),tBdot_S2{p}(1:length(VmagRatio_S2{p})),VmagRatio_S2{p})
ylim([tStart,tEnd])
zlim([0,7])
grid on
view([-60,40])
end
title(['S_2, Vmag'])

figure
hold on
for p = 1:length(shot)
h2(p) = plot3(R(p)*ones(size(VmagRatio_S3{p})),tBdot_S4{p}(1:length(VmagRatio_S4{p})),VmagRatio_S4{p})
ylim([tStart,tEnd])
zlim([0,15])
title([num2str(shot(p)),' ,Vmag'])
grid on
view([-60,40])
end
title(['S_4, Vmag'])

% #########################################################################
% Phase data
% #########################################################################
figure
N = ceil(sqrt(length(shot)));
hold on
for p = 1:length(shot)
plot3(R(p)*ones(size(P_S1{p})),tBdot_S1{p}(1:length(P_S1{p})),P_S1{p}-0*mean(P_S1{p}(1:100)))
ylim([tStart,tEnd])
zlim([-250,250])
grid on
view([-60,40])
end
title(['S_1, Phase'])

figure;
hold on
for p = 1:length(shot)
plot3(R(p)*ones(size(P_S2{p})),tBdot_S2{p}(1:length(P_S2{p})),P_S2{p}-0*mean(P_S2{p}(1:100)))
ylim([tStart,tEnd])
zlim([-250,250])
grid on
view([-60,40])
end
title(['S_2, Phase'])

figure; 
hold on
for p = 1:length(shot)
plot3(R(p)*ones(size(P_S4{p})),tBdot_S4{p}(1:length(P_S4{p})),P_S4{p}-0*mean(P_S4{p}(1:100)))
ylim([tStart,tEnd])
zlim([-250,250])
grid on
view([-60,40])
end
title(['S_4, Phase'])