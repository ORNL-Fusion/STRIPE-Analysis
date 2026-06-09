% Bdot probe data from August 2nd 2016, from VNA.
% script written Aug 29th 2016 by Juan F Caneses

close all
clear all

Scenario = 2; % (1) Compare both 270 degree shots, (2) Compare 90 and 270 degree shots 

% Producing FileName array:
TimeData = importdata('Ch2_Time.csv');
N = size(TimeData(:,1));

for n = 1:N
    a1 = num2str(TimeData(n,1));
    a2 = num2str(TimeData(n,2));
    if length(a2) == 2
            Time{n} = [a1,'_',a2];
    else
            Time{n} = [a1,'_0',a2];
    end
    FileNameCh1{n} = ['Ch1_',Time{n},'.csv'];
    DataCh1{n} = importdata(FileNameCh1{n},',',8);
    
    FileNameCh2{n} = ['Ch2_',Time{n},'.csv'];
    DataCh2{n} = importdata(FileNameCh2{n},',',8);
    
    t{1}{n} = DataCh1{n}. data(:,1) + 4.15;

    tStart = 4.15;
    tEnd   = 4.32;
    rng = find(t{1}{n} >= tStart & t{1}{n} <= tEnd);

    t{1}{n} = DataCh1{n}. data(rng,1) + 4.15;
    Re{1}{n} = DataCh1{n}. data(rng,2);
    Im{1}{n} = DataCh1{n}. data(rng,3);

    t{2}{n} = DataCh2{n}. data(rng,1) + 4.15;
    Re{2}{n} = DataCh2{n}. data(rng,2);
    Im{2}{n} = DataCh2{n}. data(rng,3);
end

%% Acquiring Ne and Te data
if  mdsconnect('mpexserver') == 1
Stem = '\MPEX::TOP.';
Branch = 'MACHOPS1:';
RootAddress = [Stem,Branch];
DataAddress{1} = [RootAddress,'LP_V_RAMP']; % V
DataAddress{2} = [RootAddress,'TARGET_LP']; % I
Config.tStart = 4.15; % [s]
Config.tEnd = 4.32;
Config.V_Att = 2;  % Output voltage of DLP box (Voltage) = V_Att*Digitized data 
Config.I_Att = 2;  % Output voltage of DLP box (Current) = I_att*Digitized data
Config.V_cal = [12.05,0.205];   % Voltage output of DLP = V_cal(1)*Output voltage of DLP box + V_cal(2)
Config.I_cal = [-142.5, 0.015]; % Current output of DLP = (I_cal(2) + Output voltage of DLP box)/Ical(1)
Config.FilterDataInput = 1; % Filter input data with savitsky Golay filter order 3 frame size 7
Config.TimeMode = 2; % Effective time of sweep, (1) start of ramp or (2) mean time of ramp
Config.AMU = 2; % Ion mass in AMU
Config.L_tip = 1.2/1000;
Config.D_tip = 0.254/1000; % [m]
Config.FitFunction = 2; 
Config.AreaType = 1; % Cylindrical + cap

if Scenario == 1 % Compare both 270 degree shots (last two)
    shotlist = 9800 + [52,53]; 
elseif Scenario == 2 % Compare 90 and 270 degree shots 
    shotlist = 9800 + [51,53];
elseif Scenario == 3; % Compare 90 and 90 degree shots 
     shotlist = [9795,9851];
end

[Ni,Te,time,Ifit,Ip,Vp,tm,Vsweep,Isweep] = DLP_fit_V1(Config,shotlist,DataAddress);

C = {'k','k:'};
t_0 = t_zero(shotlist);

figure; 
subplot(4,1,1); hold on
for s = 1:length(shotlist)
    h(s) = plot(time{s},Ni{s},C{s},'LineWidth',1);
    L{s} = [num2str(shotlist(s)),' ,t=',num2str(t_0{s}(10:14))];
end

legend(h,L,'location','NorthWest')

ylim([0,6e19])

subplot(4,1,2); hold on
for s = 1:length(shotlist)
    plot(time{s},Te{s},C{s},'LineWidth',1)
end
ylim([0,10])
end

%% Comparing phase for a 180 degree rotation:

% Shots 28 and 29 are at r = 5 cm and angle 90 degrees
% Shots 30 and 31 same but angle 270 degrees

% Shots 28 data was corrupted

if Scenario == 1; % Compare both 270 degree shots (last two)
    L = [30,31];
elseif Scenario == 2; % Compare 90 and 270 degree shots 
    L = [29,31];
    L = [22,23];
elseif Scenario == 3; % Compare 90 and 90 degree shots 
    L = [13,29];
end

C = {'k','r','k:','r:'};
LW = [1,1,1,1];

% Amplitude
%figure;
subplot(4,1,3); hold on
k = 1;
for s = L;
    Ch = 1; % 'k'
    A{Ch}{s} = sqrt(Re{Ch}{s}.^2  + Im{Ch}{s}.^2);
    plot(t{Ch}{s},A{Ch}{s},C{k},'LineWidth',LW(k));
    k = k + 1;
     
    Ch = 2; % 'r'
    A{Ch}{s} = sqrt(Re{Ch}{s}.^2  + Im{Ch}{s}.^2);
    plot(t{Ch}{s},A{Ch}{s},C{k},'LineWidth',LW(k));
    k = k + 1;
end
ylim([0,1])

% Phase
subplot(4,1,4); hold on
k = 1;
a(L) = [1,1]; % use this to flip the phase and check the 180 degree flip
for s = L;
    Ch = 1;
    P{Ch}{s} = a(s)*(atan2(Im{Ch}{s},Re{Ch}{s}))*180/pi;
    plot(t{Ch}{s},P{Ch}{s},C{k},'LineWidth',LW(k));
    k = k + 1;
    
    Ch = 2;
    P{Ch}{s} = a(s)*(atan2(Im{Ch}{s},Re{Ch}{s}))*180/pi;
    plot(t{Ch}{s},P{Ch}{s},C{k},'LineWidth',LW(k));
    k = k+1;
end
set(gca,'YTick',[-180:60:180])
ylim([-180,180])

set(gcf,'Position',[360 75 560 542],'color','w')
set(findobj('-Property','YTick'),'box','on','xlim',[4.18,4.32])

% RESULTS show that Ch1 ('k') is the signal that flips 180 degrees
% when the probe is rotated 180 degrees. this indicates that the
% probe is responding to the magnetic signal only. we see that the exact
% value is actually 15 degrees short of 180 degrees. this 15 degree angle
% may be due to finite capacitive coupling and/or finite change in the
% plasma conditions.
% The amplitude between the shots oscillates in time but are in average the
% same
% Preliminary results show that the O-scope data Ch3 changes approximately
% 180 degrees phase when the probe is rotated 180 degrees
