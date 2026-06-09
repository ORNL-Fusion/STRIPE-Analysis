% This script reads the data acquired with the RF phase detector installed
% in ProtoMPEX on Sept 16th 2016.
% Looking at the 180 degree phase shift of the data:

% #########################################################################
% VERSION 2 OF CODE:
% This new version uses the "BdotPhaseDetect_v1" code to compute the
% phases.
% #########################################################################

close all
clear all

shot = 1e4 + 300 + [70,72,73,78];
% 70: +90 deg
% 72: +90 deg
% 73: -90 deg
% 78: -90 deg

% Address of data:
Data{1} = ['\MPEX::TOP.MACHOPS1:INT_4MM_1'];
Data{2} = ['\MPEX::TOP.MACHOPS1:INT_4MM_2'];
Data{3} = ['\MPEX::TOP.MACHOPS1:INT_2MM_1'];
Data{4} = ['\MPEX::TOP.MACHOPS1:INT_2MM_2'];

% name of plots:
Title{1} = 'Vp0';
Title{2} = 'Vp90';
Title{3} = 'Vmag0';
Title{4} = 'None';
m_v = 1.8/60; % V/dB slope of the amplitude ratio signal
tStart = 4.23;
tEnd   = 4.33;
ClipInTime = 0;

[VmagRatio,P,t] = BdotPhaseDetect_v1(shot,Data);

%%
close all
for s = 1:length(shot)
    Pa(:,s) = unwrap(P{s}*pi/180)*180/pi;
end

figure; hold on 
for s = 1:length(shot)
    plot(Pa(:,s),'.')
end

Pa = unwrap(Pa*pi/180,[],2)*180/pi;

figure; hold on
for s = 1:length(shot)
    plot(Pa(:,s),'.')
end

%% Plot Amplitude
figure;
hold on
for s = 1:4;    
    plot(t{s},sgolay_t(VmagRatio{s},3,39));
end
ylim([0,10]); 
xlim([4.15,4.33])
title('Amplitude Ratio (Linear)')
ylabel('[Input/reference]')
xlabel('t [s]')                        
box on
set(gcf,'color','w')

%% Plot Phase
N = length(shot); % N = number of shots 
M = length(P{1}); % M = number of points in trace of single shot

for s = 1:4;
    P_u{s}= unwrap(P{s}*pi/180)*180/pi;
    PhaseEnd = mean(P_u{s}((M-100):M));
    P_u{s} = P_u{s} - PhaseEnd;
end

figure; 
subplot(2,1,1); hold on
Offset = [360,0,0,360] + 360 + 80;
for s = 1:N;
h(s) = plot(t{1},P_u{s}-Offset(s));
end
legend(h,num2str(shot'),'location','West')
set(findobj('Type','line'),'LineWidth',1)
ylim([0,400]); xlim([tStart,tEnd])
set(gca,'YTick',[0:60:360]); grid on
ylabel('Angle [deg]')
title('Unwrapped phase'); box on

subplot(2,1,2); hold on
for s = 1:N;
h(s) = plot(t{1},P{s});
end
xlim([tStart,tEnd])
set(findobj('Type','line'),'LineWidth',1)
title('Phase before unwrapping'); box on
ylabel('Angle [deg]')
set(gca,'YTick',[-180:60:180]); grid on
set(gcf,'color','w')

% Results:
% the data clearly indicates a 180 degree phase shift due to the rotation
% of the probe.
% Clearly the unwrapped data is the best method to use