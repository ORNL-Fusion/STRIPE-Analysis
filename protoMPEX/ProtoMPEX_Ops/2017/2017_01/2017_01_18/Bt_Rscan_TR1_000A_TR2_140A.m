% This script reads the data acquired with the RF phase detector installed
% in ProtoMPEX on Sept 16th 2016

% #########################################################################
% THIS SCRIPT USED to post process data:
% It takes  sec to call up and plot all the data N =  shots
% it takes  sec to retrieve the above data if saved in a .mat file
% #########################################################################

close all
clear all

Extract = 0;

if Extract
% Bz scan:
shot = 13000 + [142,143,144 ,145 ,146 ,147 ,148,149,150,152,153,154];
R =            [8.5,9.0,10.0,11.0,12.0,13.0,8.0,7.0,6.0,5.0,4.0,3.0];

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
            frame = 7;
            f{ch}{s} = sgolay_t(f{ch}{s}(rng),3,frame);
            plot(t{ch}{s},f{ch}{s})
    end
    ylim([0,2]);
    title(Title{ch})
end
        [VmagRatio,P,tBdot] = BdotPhaseDetect_v2(shot,Data,7);
        
        % Extract the RF_FWD_PWR trace for all shots:
        Address = {'\MPEX::TOP.MACHOPS1:RF_FWD_PWR'};
        [RF_pwr,t_RF] = my_mdsvalue_v2(shot,Address);
        
        % Extract Gas data for specfic shots
        Address = {'\MPEX::TOP.MACHOPS1:PG4'};
        [PG4,t_pg4] = my_mdsvalue_v2([13147,13154],Address);
        Address = {'\MPEX::TOP.MACHOPS1:PG3'};
        [PG6,t_pg6] = my_mdsvalue_v2([13147,13154],Address);

        
        % Extract Ion sat of DLP6.5 for specific shots
        Address = {'\MPEX::TOP.MACHOPS1:TARGET_LP'}; % Helicon signal
        [Isat,t_isat] = my_mdsvalue_v2(13154,Address);
        
        save('Bt_TR1_000A_TR2_140A_data')
else
        load('Bt_TR1_000A_TR2_140A_data')
end

%% Plot data
close all

tStart = 4.1;
tEnd  = 4.4;
 
C = {'k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:'};

% Preview data
if 0
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
end

% Radial plots:
Tstart = 4.15;
Tend   = 4.31;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Phase---------------------FIGURE 1---------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
POffset = zeros(size(shot));
nOffset2 = find(R<=8.5);
POffset(nOffset2) = 2*pi;
nOffset1 = find(R<=4);
POffset(nOffset1) = 2*pi;

figure; hold on
for s = 1:length(shot)
    rng = find(tBdot{s}>=4.273 & tBdot{s}<=4.295);
    P{s}(rng) = (unwrap(P{s}(rng)*pi/180) + POffset(s))*180/pi;
    plot3(tBdot{s}(rng),R(s)*ones(size(P{s}(rng))),P{s}(rng),'marker','.','LineStyle','none')
end
xlim([tStart,tEnd])
ylim([0,15])
% zlim([-5,5])
view([70,25])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Magnitude-----------------FIGURE 2---------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%==========================================================================
% see entire radial plot in time
if 1
figure; hold on
for s = 1:length(shot)
    hdp = find(VmagRatio{s}>6);
    VmagRatio{s}(hdp) = NaN;
    rng = find(tBdot{s}>=4.15 & tBdot{s}<=4.31);
    if 1
        Amplitude = sgolay_t(VmagRatio{s}(rng),3,61);
    else
        Amplitude = VmagRatio{s}(rng);
    end
    plot3(tBdot{s}(rng),R(s)*ones(size(VmagRatio{s}(rng))),Amplitude)
end
xlim([4.15,4.31])
ylim([0,15])
zlim([0,8])
view([70,25])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RF power trace------------FIGURE 3---------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure; hold on
for s = 1:length(shot)
    plot(t_RF{s},-10*RF_pwr{s})
end
plot(t_pg4{1},PG4{1}*10,'k')
plot(t_pg4{2},PG4{2}*10,'r')
plot(t_pg6{1},PG6{1}*2,'k:')
plot(t_pg6{2},PG6{2}*2,'r:')

    s = 2;
rng = find(tBdot{s}>=4.16 & tBdot{s}<=4.35);
plot(tBdot{s}(rng), 0.7*VmagRatio{s}(rng))

plot(t_isat{1},abs(Isat{1})*2.5,'g')
ylabel('mTorr')
xlim([4.1,4.7])
ylim([0,18])
end
%==========================================================================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% On-axis Bdot in time------FIGURE 4---------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% View on-axis in time
figure; hold on
s = 2;
rng = find(tBdot{s}>=4.16 & tBdot{s}<=4.35);
plot(tBdot{s}(rng), VmagRatio{s}(rng))
plot(tBdot{s}(rng), sgolay_t(VmagRatio{s}(rng),3,371),'k','LineWidth',2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Low density part ------------FIGURE 5------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure; hold on
for s = 1:length(shot)
    rng = find(tBdot{s}>=4.175 & tBdot{s}<=4.18);
    if 0
        Amplitude = sgolay_t(VmagRatio{s}(rng),3,11);
    else
        Amplitude = VmagRatio{s}(rng);
    end
        A_mean(s) = mean(Amplitude);
        A_std(s) = std(Amplitude,1,1);
    plot3(tBdot{s}(rng),R(s)*ones(size(Amplitude)),Amplitude)
end
xlim([4.15,4.20])
ylim([0,15])
zlim([0,8])
view([70,25])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Low density mean ------------FIGURE 6------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure; 
errorbar(R,A_mean,A_std,'ko')
ylim([0,8])
xlim([0,15])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% High density region in time----------FIGURE 7----------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

figure; hold on
for s = 1:length(shot)
    rng = find(tBdot{s}>=4.273 & tBdot{s}<=4.295);
    if 1
        Amplitude = sgolay_t(VmagRatio{s}(rng),3,61);
        Phase     = sgolay_t(P{s}(rng),3,61); 
    else
        Amplitude = VmagRatio{s}(rng);
        Phase     = P{s}(rng);
    end
    A_mean_HD(s) = mean(Amplitude);
    A_std_HD(s) = std(Amplitude,1,1);
    plot3(tBdot{s}(rng),R(s)*ones(size(VmagRatio{s}(rng))),Amplitude)
    
    P_mean_HD(s) = mean(Phase);
    P_std_HD(s) = std(Phase,1,2);
end
xlim([4.23,tEnd])
ylim([0,15])
zlim([0,8])
view([70,25])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% High density region mean-------------FIGURE 8----------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

figure; 
[~,n] = sort(R);
subplot(2,1,1); hold on
errorbar(R(n),A_mean_HD(n),A_std_HD(n),'ko-')
ylim([0,6])
xlim([0,15])

subplot(2,1,2); hold on
errorbar(R(n),P_mean_HD(n)/180,P_std_HD(n)/180,'ko-')
% ylim([0,8])
ylabel('\pi [rad]')
xlim([0,15])
