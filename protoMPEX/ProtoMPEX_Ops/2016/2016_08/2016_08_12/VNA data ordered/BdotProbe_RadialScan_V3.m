% Analysing Bdot probe data from Aug 12th, 2016
% ==================================================
% Glitches and jumps will be fixed by adding +- 2pi radians 
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
    tEnd   = 4.34;

    rng = find(t{1}{n} >= tStart & t{1}{n} <= tEnd);

    % #####################################################################
    % Channel 1 from VNA
    t{1}{n} = DataCh1{n}. data(rng,1) + 4.15;
    Re{1}{n} = a*DataCh1{n}. data(rng,2); % Real part
    Im{1}{n} = a*DataCh1{n}. data(rng,3); % Imaginary part
    
    A{1}{n} = sqrt(Re{1}{n}.^2  + Im{1}{n}.^2);  % Amplitude [Gain]
    P{1}{n} = unwrap(atan2(Im{1}{n},Re{1}{n})); % Phase [Radians]

    % #####################################################################
    % Channel 2 from VNA
    t{2}{n} = DataCh2{n}. data(rng,1) + 4.15;
    Re{2}{n} = a*DataCh2{n}. data(rng,2);
    Im{2}{n} = a*DataCh2{n}. data(rng,3);
    
    A{2}{n} = sqrt(Re{2}{n}.^2  + Im{2}{n}.^2);
    P{2}{n} = unwrap(atan2(Im{2}{n},Re{2}{n})); % Radians
end

% to acess the Rf fwd power use the following:
% [a,b] = my_mdsvalue_v3(Shot,['\MPEX::TOP.MACHOPS1:RF_FWD_PWR'])

%% Plot all data in subplots

% NOTE:
% By looking at the phase data and then at the amplitude data, one can
% determine if the VNA data cannot be used. 
% For example, for Ch1 and Ch2 shots 1-5,9, 19, 20 and 26 have random phase AND zero amplitude 

close all
yLimMax_A = 4;
yLimMax_P = 4; % pi Radians
yLimMin_P = -8; % pi Radians
Ch = 1;

figure; 
for s = 1:16;
    subplot(4,4,s); hold on
    yyaxis left
    plot(t{Ch}{s},P{Ch}{s}/pi,'-') % in pi Radians
    ylim([yLimMin_P,yLimMax_P])
    xlim([tStart,tEnd])
    if R(s) == 5
        s
    end
    yyaxis right
    plot(t{Ch}{s},A{Ch}{s},'lineWidth',1)
    ylim([0,yLimMax_A])
    xlim([tStart,tEnd])
        
    title(['R: ',num2str(R(s)),', s:',num2str(PlasmaShot(s))],'FontSize',8)
    box on
end
set(gcf,'color','w')

figure
for s = 1:13;
    subplot(4,4,s); hold on
    yyaxis left
    plot(t{Ch}{s+16},P{Ch}{s+16}/pi,'-') % in pi Radians
    ylim([yLimMin_P,yLimMax_P])
    xlim([tStart,tEnd])
    
    yyaxis right
    plot(t{Ch}{s+16},A{Ch}{s+16},'lineWidth',1)
    ylim([0,yLimMax_A])
    xlim([tStart,tEnd])
        
    title(['R: ',num2str(R(s+16)),', s:',num2str(PlasmaShot(s+16))],'FontSize',8)
    box on
end
set(gcf,'color','w')

% Data shows that for ch1, the phase is stable for all times greater and
% including 4.27

%% Extract data from the High density mode region:
close all
ntStable = find(t{1}{s} >= 4.30 & t{1}{s} <= 4.306);
BadShots = [1:5,9, 19, 20,26];

% On October 25th it was found that the phase data obtained from the VNA
% was negative of that from the AD8032 phase detectors, therefore we must
% multiply the phase in this data (VNA) by -1.

figure; hold on
for s = 1:N
    if sum(s == BadShots)
        continue
    end
    hL(s) = plot(t{Ch}{s},P{Ch}{s}/pi);
    % Phase
    P_hd_mean{Ch}(s) = -mean(P{Ch}{s}(ntStable)); % Radians
    dP_hd_mean{Ch}(s) = std(P{Ch}{s}(ntStable),1);
    % Magnitude
    A_hd_mean{Ch}(s) = mean(A{Ch}{s}(ntStable));
    dA_hd_mean{Ch}(s) = std(A{Ch}{s}(ntStable),1);
    
    errorbar(t{Ch}{s}(ntStable),ones(size(ntStable))*P_hd_mean{Ch}(s)/pi,ones(size(ntStable))*dP_hd_mean{Ch}(s)/pi,'LineWidth',2)
    Lshot(s) = PlasmaShot(s);
end
ylabel('\pi [Rad]')
grid on

% 2pi phase correction
Offset = zeros(size(1:N));
OffsetBaseUnit = -2*pi;

% R = 3
Offset(5+2) = 1*OffsetBaseUnit;

% R = 3.5
Offset(5+3) = -1*OffsetBaseUnit;

% R = 4.5
Offset(5+5) = 2*OffsetBaseUnit; % Purple at 4.5
Offset(5+6) = 2*OffsetBaseUnit; % Green at 4.5
Offset(5+7) = 3*OffsetBaseUnit; % Cyan at 4.5

% R = 5
Offset(5+8) = 1*OffsetBaseUnit; % Brown at 5
Offset(27) = -2*OffsetBaseUnit; % Purple at 5
Offset(28) = 0.5*OffsetBaseUnit; % Green at 5, 180 degree flip XP
Offset(29) = 0.5*OffsetBaseUnit; % Cyan at 5,  180 degree flip XP

% R = 5.5
Offset(14) = 0*OffsetBaseUnit; % Blue at 5.5

% R = 6
Offset(15) = -1*OffsetBaseUnit; % Orange at 6

% R = 6.5
Offset(16) = -1*OffsetBaseUnit; % Yellow at 6.5
Offset(17) = 2*OffsetBaseUnit; % Purple at 6.5
Offset(18) = 2*OffsetBaseUnit; % Green at 6.5

% R = 7
Offset(21) = -1*OffsetBaseUnit; % Cyan at 7

% R = 7.5
Offset(22) = 2*OffsetBaseUnit; % Brown at 7.5

% R = 8
Offset(23) = 1*OffsetBaseUnit; % Blue at 8

% R = 8.5
Offset(24) = 0*OffsetBaseUnit; % Orange at 8.5

% R = 11
Offset(25) = 1*OffsetBaseUnit; % Yellow at 11   

figure; hold on
for s = 1:N
    if sum(s == BadShots)
        continue
    end
    hL(s) = plot3(t{Ch}{s}(ntStable),R(s)*ones(size(P{Ch}{s}(ntStable))),(P{Ch}{s}(ntStable) + Offset(s))/pi,'LineWidth',2);
end
grid on
view([80,15])
ylim([0,15])

figure; hold on
for s = 1:N
    if sum(s == BadShots)
        continue
    end
    hL(s) = plot3(t{Ch}{s}(ntStable),R(s)*ones(size(A{Ch}{s}(ntStable))),A{Ch}{s}(ntStable),'LineWidth',2);
end
grid on
view([80,15])
ylim([0,15])

% Plot High density results with errorbars
% R4.5 has shots 10,11,12
% R5 has shots 13,27,28,29
% R6.5 has shots 16,17,18
% GoodShots = [6:8,10:18,21:25,27:29];
GoodShots = [6:8,10,13:15,16,21:25];
xOffset = 0;

figure; 
subplot(2,1,1); hold on
h(1) = errorbar((R(GoodShots)-xOffset),A_hd_mean{Ch}(GoodShots),dA_hd_mean{Ch}(GoodShots),'ko-')
ylabel('$[a.u.]$','interpreter','Latex','Fontsize',13)
xlabel('$ r [cm] $','interpreter','Latex','Fontsize',13)
xlim([0,15])
box on

subplot(2,1,2); hold on
h(1) = errorbar((R(GoodShots)-xOffset),0.4 + (P_hd_mean{Ch}(GoodShots) + Offset(GoodShots))/pi,dP_hd_mean{Ch}(GoodShots)...
    ,'ko-')
ylabel('$\pi [Rad]$','interpreter','Latex','Fontsize',13)
xlabel('$ r [cm] $','interpreter','Latex','Fontsize',13)
xlim([0,15])
box on

set(gcf,'color','w')

% save the data in EXCEL format:
Radius = (R(GoodShots)-xOffset);
PhaseBz = 0.4 + (P_hd_mean{Ch}(GoodShots) + Offset(GoodShots))/pi;
dPhaseBz = dP_hd_mean{Ch}(GoodShots);
ABz = A_hd_mean{Ch}(GoodShots);
dABz = dA_hd_mean{Ch}(GoodShots);

D = [Radius,ABz',dABz',PhaseBz',dPhaseBz'];
F = {'R[cm]','Amp[a.u.]','dAmp [a.u.]','Phase[pi*Rad]','dPhase[pi*Rad]'};
FileName = 'Bz_RadialScan_HeliconModeJump_Aug12th2016.xlsx';
xlswrite(FileName,[F;num2cell(D)]);

%% Extract data from the low density mode region:
close all

figure; 
subplot(2,1,1); hold on
n = 1; % 1 for relative amplitude, 2 for relative power

s = find(R==8);
hld(4) = plot(t{Ch}{s},A{Ch}{s}.^n,'g','LineWidth',2);

s = find(R==5.5); 
hld(3) = plot(t{Ch}{s},A{Ch}{s}.^n,'m','LineWidth',2);

s = find(R==3);
hld(2) = plot(t{Ch}{s},A{Ch}{s}.^n,'r','LineWidth',2);

s = find(R==2.5);
hld(1) = plot(t{Ch}{s},A{Ch}{s}.^n,'k','LineWidth',2);

l = -(R([6,7,14,23]) - 6.5);
legend(hld,{num2str(l)},'location','NorthWest')
ylabel('Wave Amplitude [a.u]')
box on
xlim([tStart,tEnd])
ylim([0,5])

% shot 9827 corresponds to R = 8 cm, this shot cannot be analyzed with DLP
% due to glitches on V but we can read the Isat:
[Isat,tIsat] = my_mdsvalue_v3(9827,'\MPEX::TOP.MACHOPS1:TARGET_LP');

subplot(2,1,2); hold on
[na,nb] = peakseek(abs(Isat{1}),100);
plot(tIsat{1}(na) + 0.02,nb/max(nb),'k-','LineWidth',2)
xlim([tStart,tEnd])
ylim([0,1])
box on
set(gcf,'color','w')
xlabel('time [s]')
ylabel('Norm Isat [a.u]')
% ntStable = find(t{1}{s} >= 4.30 & t{1}{s} <= 4.306);
% BadShots = [1:5,9, 19, 20,26];
% 
% figure; hold on
% for s = 1:N
%     if sum(s == BadShots)
%         continue
%     end
%     hL(s) = plot(t{Ch}{s},P{Ch}{s}/pi);
%     % Phase
%     P_hd_mean{Ch}(s) = mean(P{Ch}{s}(ntStable)); % Radians
%     dP_hd_mean{Ch}(s) = std(P{Ch}{s}(ntStable),1);
%     % Magnitude
%     A_hd_mean{Ch}(s) = mean(A{Ch}{s}(ntStable));
%     dA_hd_mean{Ch}(s) = std(A{Ch}{s}(ntStable),1);
%     
%     errorbar(t{Ch}{s}(ntStable),ones(size(ntStable))*P_hd_mean{Ch}(s)/pi,ones(size(ntStable))*dP_hd_mean{Ch}(s)/pi,'LineWidth',2)
%     Lshot(s) = PlasmaShot(s);
% end
% ylabel('\pi [Rad]')
% grid on



% Things to do:
% Gather amplitude of wavefields in the low density region