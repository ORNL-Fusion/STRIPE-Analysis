% Bdot probe data from August 2nd 2016, from VNA.
% script written Aug 29th 2016 by Juan F Caneses

close all
clear all

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

    tStart = 4.26;
    tEnd   = 4.32;
    rng = find(t{1}{n} >= tStart & t{1}{n} <= tEnd);

    t{1}{n} = DataCh1{n}. data(rng,1) + 4.15;
    Re{1}{n} = DataCh1{n}. data(rng,2);
    Im{1}{n} = DataCh1{n}. data(rng,3);

    t{2}{n} = DataCh2{n}. data(rng,1) + 4.15;
    Re{2}{n} = DataCh2{n}. data(rng,2);
    Im{2}{n} = DataCh2{n}. data(rng,3);
end

% Acquiring Ne and Te data
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
shotlist = 9800 + [50:53];
[Ni,Te,time,Ifit,Ip,Vp,tm,Vsweep,Isweep] = DLP_fit_V1(Config,shotlist,DataAddress);

%%
C = {'k','r','bl','g','m','k:','r:','bl:','g:','m:'};
t_0 = t_zero(shotlist);

figure; 
subplot(2,1,1); hold on
for s = 1:length(shotlist)
    h(s) = plot(time{s},Ni{s},C{s});
    L{s} = [num2str(shotlist(s)),' ,t=',num2str(t_0{s}(10:14))];
end

legend(h,L,'location','NorthWest')

ylim([0,6e19])

subplot(2,1,2); hold on
for s = 1:length(shotlist)
    plot(time{s},Te{s},C{s})
end
ylim([0,10])

%%
% Comparing phase for a 180 degree rotation:
L = [30,31];
L = [29,31];
C = {'k','r','k:','r:'};
LW = [1,1,1,1];

% Amplitude
figure; hold on
k = 1;
for s = L;
    Ch = 1;
    A{Ch}{s} = sqrt(Re{Ch}{s}.^2  + Im{Ch}{s}.^2);
    plot(t{Ch}{s},A{Ch}{s},C{k},'LineWidth',LW(k));
    k = k + 1;
    
    Ch = 2;
    A{Ch}{s} = sqrt(Re{Ch}{s}.^2  + Im{Ch}{s}.^2);
    plot(t{Ch}{s},A{Ch}{s},C{k},'LineWidth',LW(k));
    k = k + 1;
end
ylim([0,4])

% Phase
figure; hold on
k = 1;
for s = L;
    Ch = 1;
    P{Ch}{s} = (atan2(Im{Ch}{s},Re{Ch}{s}))*180/pi;
    plot(t{Ch}{s},P{Ch}{s},C{k},'LineWidth',LW(k));
    k = k + 1;
    
    Ch = 2;
    P{Ch}{s} = (atan2(Im{Ch}{s},Re{Ch}{s}))*180/pi;
    plot(t{Ch}{s},P{Ch}{s},C{k},'LineWidth',LW(k));
    k = k+1;
end
ylim([-pi,pi]*180/pi)
