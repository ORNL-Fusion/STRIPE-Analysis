
% ==================================================
% SUPERSEEDED
% this version was the first attempt to observe the data radially. we had
% only onbserved the data from 4.27 to 4.32. however to understand the
% behaviour of the fields in the low density region we have to include
% ==================================================
% Bdot probe data from August 2nd 2016, from VNA.
% script written Sept 7th 2016 by Juan F Caneses

close all
clear all

% Producing FileName array:
ShotData = importdata('UsableShotMapping_PlasmaVnaOscope.xlsx');
N = length(ShotData.data(:,1));
PlasmaShot = ShotData.data(:,1);
PlasmaShot_tzero = ShotData.data(:,2);
R = ShotData.data(:,3);
Fctr6dB = 9782; % After and including this shot, all VNA data had a 6 dB attenuator in place
% this attenuation in a 50 Ohm system means that the measured amplitude
% must be multiplied by  1.9953 to obtain the signal at the input of the
% attenuator.

% #########################################################################
% Remember that RG58 has about 3 dB attenuation per 100 feet of cable. we
% might need to add this effect to obtain the true voltage at the probe
% #########################################################################

% Oscope shot number is stored in
% ShotData.data(:,6)

for n = 1:N
    a1 = num2str(ShotData.data(n,4));
    a2 = num2str(ShotData.data(n,5));
    if length(a2) == 2
            Time{n} = [a1,'_',a2];
    else
            Time{n} = [a1,'_0',a2];
    end
    FileNameCh1{n} = ['Ch1_',Time{n},'.csv'];
    DataCh1{n} = importdata(FileNameCh1{n},',',8);
    
    FileNameCh2{n} = ['Ch2_',Time{n},'.csv'];
    DataCh2{n} = importdata(FileNameCh2{n},',',8);
    
    if PlasmaShot(n)>= Fctr6dB
        a = 1.9953;
    else
        a = 1;
    end
    
    t{1}{n} = DataCh1{n}. data(:,1) + 4.15;

    tStart = 4.15;
    tEnd   = 4.32;
%     tStart = 4.27;
%     tEnd   = 4.32;
%     tStart = 4.31;
%     tEnd   = 4.32;
    rng = find(t{1}{n} >= tStart & t{1}{n} <= tEnd);

    % #####################################################################
    % Channel 1 from VNA
    t{1}{n} = DataCh1{n}. data(rng,1) + 4.15;
    Re{1}{n} = a*DataCh1{n}. data(rng,2); % Real part
    Im{1}{n} = a*DataCh1{n}. data(rng,3); % Imaginary part
    
    A{1}{n} = sqrt(Re{1}{n}.^2  + Im{1}{n}.^2);  % Amplitude
    P{1}{n} = (atan2(Im{1}{n},Re{1}{n}))*180/pi; % Phase

    % #####################################################################
    % Channel 2 from VNA
    t{2}{n} = DataCh2{n}. data(rng,1) + 4.15;
    Re{2}{n} = a*DataCh2{n}. data(rng,2);
    Im{2}{n} = a*DataCh2{n}. data(rng,3);
    
    A{2}{n} = sqrt(Re{2}{n}.^2  + Im{2}{n}.^2);
    P{2}{n} = (atan2(Im{2}{n},Re{2}{n}))*180/pi;
end

%% Plot all data in subplots

% NOTE:
% By looking at the phase data and then at the amplitude data, one can
% determine if the VNA data cannot be used. 
% For example, for Ch1 and Ch2 shots 1-5,9, 19, 20 and 26 have random phase AND zero amplitude 
% In addition, Ch2, there are some shots (12,14,20)whose phase gitters -180
% degrees in single points, these have to be corrected.

close all
yLimMax_A = 4;
yLimMax_P = 1;
yLimMin_P = -1;
Ch = 1;

% Amplitude and Phase:
figure; 
for s = 1:16;
    subplot(4,4,s); hold on
    yyaxis left
    plot(t{Ch}{s},P{Ch}{s}/180,'-')
    ylim([yLimMin_P,yLimMax_P])
    xlim([tStart,tEnd])
    if R(s) == 4
        s
    end
    yyaxis right
    plot(t{Ch}{s},A{Ch}{s},'lineWidth',1)
    ylim([0,yLimMax_A])
    xlim([tStart,tEnd])
        
    title(['R: ',num2str(R(s)),', s:',num2str(PlasmaShot(s))],'FontSize',8)
end

figure
for s = 1:13;
    subplot(4,4,s); hold on
    yyaxis left
    plot(t{Ch}{s+16},P{Ch}{s+16}/180,'-')
    ylim([yLimMin_P,yLimMax_P])
    xlim([tStart,tEnd])
    
    yyaxis right
    plot(t{Ch}{s+16},A{Ch}{s+16},'lineWidth',1)
    ylim([0,yLimMax_A])
    xlim([tStart,tEnd])
        
    title(['R: ',num2str(R(s+16)),', s:',num2str(PlasmaShot(s+16))],'FontSize',8)
end

%% Plot all data in 3D
close all

BadShots = [1:5,9, 19, 20,26];
% Not all shots are bad, for example, shots 2 to 5 have a finite magnitude
% during the low density mode but the phase is all over the place. this
% means that we can still assess the wavefields during the low density mode
% by looking only at the magnitude.

figure; hold on
Ch = 1;
for s = 1:N
    if sum(s == BadShots)
        continue
    end
    plot3(t{Ch}{s},R(s)*ones(size(A{Ch}{s})),A{Ch}{s})
    %plot3(t{2}{s},R(s)*ones(size(A{2}{s})),A{2}{s},'r')
end
view([60,60])
ylim([0,15])

figure; hold on
for s = 1:N
    if sum(s == BadShots)
        continue
    end
    plot3(t{Ch}{s},R(s)*ones(size(P{Ch}{s})),P{Ch}{s}/180)
end
zlabel('\pi [Rad]')
view([60,60])
ylim([0,15])

%% Phase +-2*pi degree and Gitter corrections
close all
Ch = 1;
for s = 1:length(P{Ch})
    P_cor{Ch}{s} = P{Ch}{s};
end

s = 10; 
figure; hold on
h(1) = plot(P{Ch}{s});
rng = [124,125,129,133:136,156,193,194,196,220:229,308:310,338:341];
P_cor{Ch}{s} = P{Ch}{s};
P_cor{Ch}{s}(rng) = -P{Ch}{s}(rng);
%P_cor{Ch}{s}(171:end) = -P_cor{Ch}{s}(171:end);

h(2) = plot(P_cor{Ch}{s},'r');
title('10')
legend(h,'Raw','Corrected')
ylim([-180,180])

% if phase less than zero, add 360 degrees to it?
% f = (diff(P{Ch}{s})/360);
% nn = peakseek(-f,2,0.7);
% np = peakseek(f,2,0.7);
% for c = 1:length(np)
%     k = np(c)+1;
%     P_cor{Ch}{s}(k:end) = P_cor{Ch}{s}(k:end)-360;
% end
% 
% for c = 1:length(nn)
%     k = nn(c)+1;
%     P_cor{Ch}{s}(k:end) = P_cor{Ch}{s}(k:end) + 360;
% end

%%
s = 11; 
figure; hold on
h(1) = plot(P{Ch}{s});
rng = [156:159,161,162:2:166,171:2:175,179,182,186,196,224,226,239 + [21,92,95]];
P_cor{Ch}{s} = P{Ch}{s};
P_cor{Ch}{s}(rng) = -P{Ch}{s}(rng);
h(2) = plot(P_cor{Ch}{s},'r');
title('11')
legend(h,'Raw','Corrected')
ylim([-180,180])

%%
s = 12; 
figure; hold on
h(1) = plot(P{Ch}{s});
rng = [152,159:165,135,143,145,239 + [2:10,85,92,87:90,102]];
P_cor{Ch}{s} = P{Ch}{s};
P_cor{Ch}{s}(rng) = -P{Ch}{s}(rng);
h(2) = plot(P_cor{Ch}{s},'r');
title(num2str(s))
legend(h,'Raw','Corrected')
ylim([-180,180])

%%
s = 13; 
figure; hold on
h(1) = plot(P{Ch}{s});
rng = [116,123:125,147:152,154,155,157,175];
P_cor{Ch}{s} = P{Ch}{s};
P_cor{Ch}{s}(rng) = -P{Ch}{s}(rng);
%P_cor{Ch}{s}(194:end) = -P_cor{Ch}{s}(194:end);
h(2) = plot(P_cor{Ch}{s},'r');
title(num2str(s))
legend(h,'Raw','Corrected')
ylim([-180,180])

%%
s = 14; 
figure; hold on
h(1) = plot(P{Ch}{s});
rng = [155,158,166,172,175,206,210,219];
P_cor{Ch}{s} = P{Ch}{s};
P_cor{Ch}{s}(rng) = -P{Ch}{s}(rng);
%P_cor{Ch}{s}(194:end) = -P_cor{Ch}{s}(194:end);
h(2) = plot(P_cor{Ch}{s},'r');
title(num2str(s))
legend(h,'Raw','Corrected')
ylim([-180,180])
%%
s = 21; 
figure; hold on
h(1) = plot(P{Ch}{s});
rng = 239 + [1:3,96:102];
P_cor{Ch}{s} = P{Ch}{s};
P_cor{Ch}{s}(rng) = -P{Ch}{s}(rng);
h(2) = plot(P_cor{Ch}{s},'r');
title(num2str(s))
legend(h,'Raw','Corrected')

s = 22; 
figure; hold on
h(1) = plot(P{Ch}{s});
rng = 239 + [83:89,91];
P_cor{Ch}{s} = P{Ch}{s}
P_cor{Ch}{s}(rng) = -P{Ch}{s}(rng);
h(2) = plot(P_cor{Ch}{s},'r');
title(num2str(s))
legend(h,'Raw','Corrected')

PhaseOffsetType = 2;
    PhaseOffset{1} = zeros(size(P{1}));
    PhaseOffset{2} = zeros(size(P{2}));
    
switch PhaseOffsetType
    case 1
        
    case 2
    PhaseOffset{1}(11) = 2;
    PhaseOffset{1}([28,29]) = 1; % these are the shots where the probe was rotated 180 degrees    
    PhaseOffset{1}(find(R == 8.5)) = 2;
    PhaseOffset{1}(find(R == 11)) = 2;
    %PhaseOffset{1}(find(R == 7.5)) = -2;
    %PhaseOffset{1}(find(R == 8)) = -2;

    case 3
        PhaseOffset{1}([10,12]) = -2;
   
end

figure; hold on
for s = 1:N
    if sum(s == BadShots)
        continue
    end
    plot3(t{Ch}{s},R(s)*ones(size(P_cor{Ch}{s})),P_cor{Ch}{s}/180 + PhaseOffset{1}(s))
end
zlabel('\pi [Rad]')
view([60,60])
ylim([0,15])

% Things to do:
% Gather amplitude of wavefields in the low density region
% Correct phase using the if phase is well below 170 degrees then add 360
% degrees to it
% Take statistics of wavefeilds on the high density mode


