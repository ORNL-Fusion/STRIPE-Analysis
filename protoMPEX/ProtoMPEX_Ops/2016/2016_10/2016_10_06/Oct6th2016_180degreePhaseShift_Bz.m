% October 6th Bdot probe XP
% Testing the 180 degree phase shift measurement at x = 8 cm

close all
clear all

shotlist = 1e4 + 700 + [35,36,37];
% 35: 90  deg
% 36: 270 deg
% 37: 270 deg

% Shots 35 and 36 are very similar in magnitude from 4.28 to 4.26

% Acquiring Ne and Te data
Stem = '\MPEX::TOP.';
Branch = 'MACHOPS1:';
RootAddress = [Stem,Branch];
DataAddress{1} = [RootAddress,'LP_V_RAMP']; % V
DataAddress{2} = [RootAddress,'TARGET_LP']; % I
Data{5} = [RootAddress,'PWR_28GHz'];
[G,t_G] = my_mdsvalue_v3(shotlist,Data{5});

Data{1} = ['\MPEX::TOP.MACHOPS1:INT_4MM_1'];
Data{2} = ['\MPEX::TOP.MACHOPS1:INT_4MM_2'];
Data{3} = ['\MPEX::TOP.MACHOPS1:INT_2MM_1'];
Data{4} = ['\MPEX::TOP.MACHOPS1:INT_2MM_2'];

[VmagRatio,P,tBdot] = BdotPhaseDetect_v1(shotlist,Data);

% The DLP data for the 180 degree phase shift test has a strong 50 Hz
% component in it making the data analysis difficult.
if 0
    Config.tStart = 4.10; % [s]
    Config.tEnd = 4.38;
    Config.V_Att = 2;  % Output voltage of DLP box (Voltage) = V_Att*Digitized data 
    Config.I_Att = 2;  % Output voltage of DLP box (Current) = I_att*Digitized data
    Config.V_cal = [12.05,0.205];   % Voltage output of DLP = V_cal(1)*Output voltage of DLP box + V_cal(2)
    Config.I_cal = [-142.5, 1*0.015]; % Current output of DLP = (I_cal(2) + Output voltage of DLP box)/Ical(1)
    Config.FilterDataInput = 1; % Filter input data with savitsky Golay filter order 3 frame size 7
    Config.TimeMode = 2; % Effective time of sweep, (1) start of ramp or (2) mean time of ramp
    Config.AMU = 2; % Ion mass in AMU
    Config.L_tip = 1.2/1000;
    Config.D_tip = 0.254/1000; % [m]
    Config.FitFunction = 2; 
    Config.AreaType = 1; % Cylindrical + cap
    [Ni,Te,time,Ifit,Ip,Vp,tm,Vsweep,Isweep] = DLP_fit_V2(Config,shotlist,DataAddress);
end

%%
close all
tStart = 4.28;
tEnd  = 4.30;

% Phase data
figure; 
subplot(2,2,1);
hold on;
for s = 1:length(shotlist)
    hP(s) = plot(tBdot{s},P{s});
end
title('Phase (raw)')
legend(hP,num2str(shotlist'))
ylim([-180,180])
xlim([tStart,tEnd])

subplot(2,2,2);
hold on;
offset = [0,23*360,-5*360];
for s = 1:length(shotlist)
    hP(s) = plot(tBdot{s},unwrap(P{s}*pi/180)*180/pi + offset(s));
end
title('Phase (unwrapped)')
ylim([0,360]); set(gca,'YTick',[0:60:360])
xlim([tStart,tEnd])

subplot(2,2,3);
hold on;
offset = [0,23*360,-5*360];

for s = 2:3
    hP(s) = plot(tBdot{s},unwrap(P{1}*pi/180)*180/pi + offset(1) - (unwrap(P{s}*pi/180)*180/pi + offset(s)));
end
line([0,5],[180,180],'LineStyle','--','LineWidth',2,'color','k')
set(hP(2),'color','red'); set(hP(3),'color','green')
title('Phase difference (unwrapped)')
legend([hP(2:3)],num2str(shotlist(2:3)'),'Location','SouthWest')
ylim([100,200]); set(gca,'YTick',[100:20:200]); grid on
xlim([tStart,tEnd])

subplot(2,2,4); hold on
for s = 1:length(shotlist)
    if s==21
        continue
    end
    hM(s) = plot(tBdot{s},VmagRatio{s});
end
title('Magnitude')
legend(hM,num2str(shotlist'))
ylim([0,10])
xlim([tStart,tEnd])

set(gcf,'color','w')