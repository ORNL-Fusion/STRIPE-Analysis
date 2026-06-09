% This code reads data (csv) format from the oscilloscope
% data taken on April 11th 2017
% Code created April 25th 2017

clear all
close all

% Load csv data:
load('RF_data.mat'); % Ch1: fwd power, ch2: bwd power

% Find envelop of fwd power trace
[a,b]=peakseek(CH1,40,0.01);

% A typical shot that had a 11.9 dBm drive is 13850
Stem = '\MPEX::TOP.';
Branch = 'MACHOPS1:';
RootAddress = [Stem,Branch];
DA{1} = [RootAddress,'RF_FWD_PWR'];
[RF,t_rf]   = my_mdsvalue_v2(13850,DA(1))


%%
close all

% Voltage (normalized)
figure; hold on
[n] = max(CH1(100:end));
% plot(TIME,CH1/n)
h(1) = plot(TIME(a),sgolay_t(CH1(a)/n,3,13),'k','lineWidth',2)
[n] = max(-RF{1});
h(2) = plot(t_rf{1}(1:end-1)-4.155,-RF{1}/n,'r','lineWidth',2)
xlim([-0.05,0.25])
ylim([0,1])
legend(h,'Directional coupler','RF power detector','location','NorthWest')
title('Voltage (normalized)')
set(gcf,'color','w')
grid on
box on

% power (normalized)
figure; hold on
[n] = max(CH1(100:end).^2);
% plot(TIME,CH1/n)
h(1) = plot(TIME(a),sgolay_t(CH1(a).^2,3,13)/n,'k','lineWidth',2)
[n] = max(RF{1}.^2);
h(2) = plot(t_rf{1}(1:end-1)-4.155,(RF{1}.^2)/n,'r','lineWidth',2)
xlim([-0.05,0.25])
ylim([0,1])
legend(h,'Directional coupler','RF power detector','location','NorthWest')
title('Power (normalized)')
set(gcf,'color','w')
grid on
box on