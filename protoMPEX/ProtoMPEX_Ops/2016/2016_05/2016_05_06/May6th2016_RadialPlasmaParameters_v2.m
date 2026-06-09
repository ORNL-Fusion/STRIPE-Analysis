% Data analysis for May 6th 2016
% See page 245 to 248 Logbook 
% Experiment:
% " D2 plasma, Mode jump plasma test"

close all
clear all

% Connect to Server: 
mdsconnect('mpexserver') 

% Acquiring Ne and Te data
Stem = '\MPEX::TOP.';
Branch = 'MACHOPS1:';
RootAddress = [Stem,Branch];
Data{1} = [RootAddress,'LP_V_RAMP']; % V
Data{2} = [RootAddress,'TARGET_LP']; % I
Data{3} = [Stem,'FSCOPE:','TUBE07:','PMT_VOLT']; % 28 GHz pulse
Data{4} = [Stem,'SHOT_NOTE'];


%% 1 - Spool 6.5, Radial scan
% See Logbook page 247, for shots 8602 and onwards, the current uses a 4:1
% attenuator.

Shots = 8600 + [2:19,22];
R = [0,0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,-0.5,-1,-1.5,-2,-2,-2.5,-3];

% Configuration file:
for s = 1
Config.tStart = 4.15; % [s]
Config.tEnd = 4.32;
Config.V_Att = 2;  % Output voltage of DLP box (Voltage) = V_Att*Digitized data 
Config.I_Att = 4;  % Output voltage of DLP box (Current) = I_att*Digitized data
Config.V_cal = [12.05,0.205];   % Voltage output of DLP = V_cal(1)*Output voltage of DLP box + V_cal(2)
Config.I_cal = [-142.5, 0.015]; % Current output of DLP = (I_cal(2) + Output voltage of DLP box)/Ical(1)
Config.FilterDataInput = 1; % Filter input data with savitsky Golay filter order 3 frame size 7
Config.TimeMode = 2; % Effective time of sweep, (1) start of ramp or (2) mean time of ramp
Config.AMU = 2; % Ion mass in AMU
Config.L_tip = 1.3/1000;
Config.D_tip = 0.254/1000; % [m]
Config.FitFunction = 2; 
Config.AreaType = 1; % Cylindrical + cap
end
[Ni,Te,time,Ifit,Ip,Vp,tm,V,I] = DLP_fit_V2(Config,Shots,Data);

% Plot data:
close all

figure; 
for s = 1:18
    subplot(5,4,s); hold on
    plot(time{s},Ni{s})
    title(['R',num2str(R(s))])
    ylim([0,1e20])
end

t_HD_Gather_Start = 4.28;
t_HD_Gather_End   = 4.31;

t_LD_Gather_Start = 4.16;
t_LD_Gather_End   = 4.19;

for s = 1:length(Shots)
    if s ==16 | s == 19
        Ni_HD(s) = NaN;
        Te_HD(s) = NaN;
        
        Ni_LD(s) = NaN;
        Te_LD(s) = NaN;
        continue
    end
    rng2{s} = find(time{s}>=t_HD_Gather_Start & time{s}<=t_HD_Gather_End);
    Ni_HD(s) = mean(Ni{s}(rng2{s}));
    Te_HD(s) = mean(Te{s}(rng2{s}));
    
    rng3{s} = find(time{s}>=t_LD_Gather_Start & time{s}<=t_LD_Gather_End);
    Ni_LD(s) = mean(Ni{s}(rng3{s}));
    Te_LD(s) = mean(Te{s}(rng3{s}));
end

figure; 
subplot(2,1,1); hold on
plot(R,Ni_HD,'ko')
plot(R,Ni_LD,'ro')
xlim([-5,5])
 
subplot(2,1,2); hold on
plot(R,Te_HD,'ko')
plot(R,Te_LD,'ro')
ylim([0,10]); xlim([-5,5])

% Fit Test:
figure;
% k = length(Shots) - 1;
k = find(R==1.5);
n = 45;
for s = 1:16
subplot(4,4,s); hold on
plot(V{k}{s+n},I{k}{s+n},'k')
plot(V{k}{s+n},Ifit{k}{s+n},'r')
title([num2str(time{k}(s+n)),'r',num2str(R(k))])
end

s = 14;
figure; hold on
plot(V{k}{s+n},I{k}{s+n},'k')
plot(V{k}{s+n},Ifit{k}{s+n},'r')
title(['time : ',num2str(time{k}(s+n)),', R = ',num2str(R(k))])

set(gcf,'color','w')

% Save data to excel spreadsheet:
% remove NaN
rng = ~isnan(Ni_HD);

D = [R(rng)',Ni_HD(rng)',Te_HD(rng)'];
F = {'R [cm]','n_i [m^-3]','T_e [eV]'};
FileName = 'NiTe_Spool_6_5.xlsx';
xlswrite(FileName,[F;num2cell(D)]);

%% 1 - Spool 9.5, Radial scan
% See Logbook page 246, for these shots both the current and voltage
% attenuators are 2:1

Shots = 8500 + [49,51,54,55:63];
R = [0,0.5,1,1.5,2,2.5,3,3.5,-0.5,-1,-1.5,-2]; 

% Configuration file:
for s = 1
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
end
[Ni,Te,time,Ifit,Ip,Vp,tm,V,I] = DLP_fit_V2(Config,Shots,Data);

% Plot data:
clear Ni_HD Ni_LD Te_HD Te_LD

figure; 
for s = 1:length(Ni)
    subplot(5,4,s); hold on
    plot(time{s},Ni{s})
    title(['R',num2str(R(s))])
    ylim([0,1e20])
end

t_HD_Gather_Start = 4.28;
t_HD_Gather_End   = 4.31;

t_LD_Gather_Start = 4.16;
t_LD_Gather_End   = 4.19;

for s = 1:length(Shots)
    if s ==16 | s == 19
        Ni_HD(s) = NaN;
        Te_HD(s) = NaN;
        
        Ni_LD(s) = NaN;
        Te_LD(s) = NaN;
        continue
    end
    rng2{s} = find(time{s}>=t_HD_Gather_Start & time{s}<=t_HD_Gather_End);
    Ni_HD(s) = mean(Ni{s}(rng2{s}));
    Te_HD(s) = mean(Te{s}(rng2{s}));
    
    rng3{s} = find(time{s}>=t_LD_Gather_Start & time{s}<=t_LD_Gather_End);
    Ni_LD(s) = mean(Ni{s}(rng3{s}));
    Te_LD(s) = mean(Te{s}(rng3{s}));
end

figure; 
subplot(2,1,1); hold on
plot(R,Ni_HD,'ko')
plot(R,Ni_LD,'ro')
xlim([-5,5])
 
subplot(2,1,2); hold on
plot(R,Te_HD,'ko')
plot(R,Te_LD,'ro')
ylim([0,10]); xlim([-5,5])

% Fit Test:
figure;
% k = length(Shots) - 1;
k = find(R==1.5);
n = 45;
for s = 1:16
subplot(4,4,s); hold on
plot(V{k}{s+n},I{k}{s+n},'k')
plot(V{k}{s+n},Ifit{k}{s+n},'r')
title([num2str(time{k}(s+n)),'r',num2str(R(k))])
end

s = 14;
figure; hold on
plot(V{k}{s+n},I{k}{s+n},'k')
plot(V{k}{s+n},Ifit{k}{s+n},'r')
title(['time : ',num2str(time{k}(s+n)),', R = ',num2str(R(k))])

set(gcf,'color','w')

% Save data to excel spreadsheet:
% remove NaN
rng = ~isnan(Ni_HD);

D = [R(rng)',Ni_HD(rng)',Te_HD(rng)'];
F = {'R [cm]','n_i [m^-3]','T_e [eV]'};
FileName = 'NiTe_Spool_9_5.xlsx';
xlswrite(FileName,[F;num2cell(D)]);

%% 1 - Spool 10.5, Radial scan
% See Logbook page 246, for these shots both the current and voltage
% attenuators are 2:1

Shots = 8500 + [65:71,75:78];
R = [0,0.5,1,1.5,2,2.5,3,-0.5,-1,-1.5,-2];

% Configuration file:
for s = 1
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
end
[Ni,Te,time,Ifit,Ip,Vp,tm,V,I] = DLP_fit_V2(Config,Shots,Data);

% Plot data:
clear Ni_HD Ni_LD Te_HD Te_LD

figure; 
for s = 1:11
    subplot(3,4,s); hold on
    plot(time{s},Ni{s})
    title(['R',num2str(R(s))])
    ylim([0,1e20])
end

t_HD_Gather_Start = 4.28;
t_HD_Gather_End   = 4.31;

t_LD_Gather_Start = 4.16;
t_LD_Gather_End   = 4.19;

for s = 1:length(Shots)
    if s ==16 | s == 19
        Ni_HD(s) = NaN;
        Te_HD(s) = NaN;
        
        Ni_LD(s) = NaN;
        Te_LD(s) = NaN;
        continue
    end
    rng2{s} = find(time{s}>=t_HD_Gather_Start & time{s}<=t_HD_Gather_End);
    Ni_HD(s) = mean(Ni{s}(rng2{s}));
    Te_HD(s) = mean(Te{s}(rng2{s}));
    
    rng3{s} = find(time{s}>=t_LD_Gather_Start & time{s}<=t_LD_Gather_End);
    Ni_LD(s) = mean(Ni{s}(rng3{s}));
    Te_LD(s) = mean(Te{s}(rng3{s}));
end

figure; 
subplot(2,1,1); hold on
plot(R,Ni_HD,'ko')
plot(R,Ni_LD,'ro')
xlim([-5,5])
 
subplot(2,1,2); hold on
plot(R,Te_HD,'ko')
plot(R,Te_LD,'ro')
ylim([0,10]); xlim([-5,5])

% Save data to excel spreadsheet:
% remove NaN
rng = ~isnan(Ni_HD);
D = [R(rng)',Ni_HD(rng)',Te_HD(rng)'];
F = {'R [cm]','n_i [m^-3]','T_e [eV]'};
FileName = 'NiTe_Spool_10_5.xlsx';
xlswrite(FileName,[F;num2cell(D)]);

% Fit Test:
figure;
% k = length(Shots) - 1;
k = find(R==1.5);
n = 45;
for s = 1:16
subplot(4,4,s); hold on
plot(V{k}{s+n},I{k}{s+n},'k')
plot(V{k}{s+n},Ifit{k}{s+n},'r')
title([num2str(time{k}(s+n)),'r',num2str(R(k))])
end
set(gcf,'color','w')