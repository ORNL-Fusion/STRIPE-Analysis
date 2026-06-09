% Compute the PEMP 2017 combined power level
clear all
close all

shotlist = 16430;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ECH data:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AcquireData = 0;
if AcquireData
    address{2} = '\MPEX::TOP.MACHOPS1:PWR_28GHZ'; % 28 GHz power trace
    [ECH_28,t_ech28] = my_mdsvalue_v2(shotlist,address(2));
    address{2} = '\MPEX::TOP.MACHOPS1:RF_FWD_PWR'; % RF power trace
    [RF,t_rf] = my_mdsvalue_v2(shotlist,address(2));
    save('RF_ECH_16430.mat')
else
    load('RF_ECH_16430.mat')
end

figure; plot(t_ech28{1}(1:end-1),ECH_28{1}*15.625)
xlim([4.15,4.3])

ECH_kW = ECH_28{1}*15.625;
ECH_t = t_ech28{1}(1:end-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Helicon and ICH power 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('RF_ICH_pwr_Oscilloscope_16430.mat')

figure; 
subplot(1,2,1); hold on
% Plot raw data, 1e6 data points
plot(RF_t,RF_fwd)
% Find envelope:
dnpeak = 50;
% Top peaks:
[locs1,RF_pks1] = peakseek(RF_fwd,dnpeak,0);
% Bottom peaks:
[locs2,RF_pks2] = peakseek(-RF_fwd,dnpeak,0);
% Smoothen data:
frm = 17;
RF_pks1s = sgolay_t(RF_pks1,3,frm);
RF_pks2s = sgolay_t(RF_pks2,3,frm);

% Plot envelope data
plot(RF_t(locs1),RF_pks1s,'g')
plot(RF_t(locs2),-RF_pks2s,'bl')

subplot(1,2,2); hold on
% Plot raw data, 1e6 data points
plot(ICH_t,ICH_fwd)
% Find envelope:
dnpeak = 50;
% Top peaks:
[locs1_ich,ICH_pks1] = peakseek(ICH_fwd,dnpeak,0);
% Bottom peaks:
[locs2_ich,ICH_pks2] = peakseek(-ICH_fwd,dnpeak,0);
% Smoothen data:
frm = 17;
ICH_pks1s = sgolay_t(ICH_pks1,3,frm);
ICH_pks2s = sgolay_t(ICH_pks2,3,frm);

% Plot envelope data
plot(ICH_t(locs1_ich),ICH_pks1s,'g')
plot(ICH_t(locs2_ich),-ICH_pks2s,'bl')

% need to interpolate the data:
t_offset = 4.166;
RF_fwd_1 = interp1(t_offset + RF_t(locs1),RF_pks1s,ECH_t);
RF_fwd_2 = interp1(t_offset + RF_t(locs2),RF_pks2s,ECH_t);

ICH_fwd_1 = interp1(t_offset + ICH_t(locs1_ich),ICH_pks1s,ECH_t);
ICH_fwd_2 = interp1(t_offset + ICH_t(locs2_ich),ICH_pks2s,ECH_t);

figure; hold on
plot(ECH_t,RF_fwd_1);
plot(ECH_t,RF_fwd_2);
plot(ECH_t,ICH_fwd_1);
plot(ECH_t,ICH_fwd_2);

plot(ECH_t,ECH_28{1});

% Plot RF power in kW
figure; hold on
HeliconPower = 100*((RF_fwd_1+RF_fwd_2)/1.28).^2;
ICHPower = 10*((ICH_fwd_1+ICH_fwd_2)/0.194).^2;
ECHPower = ECH_kW;

h(1) = plot(ECH_t,HeliconPower,'g')
h(2) = plot(ECH_t,ICHPower,'r');
h(3) = plot(ECH_t,ECH_kW,'bl')
line([0,5],[200,200],'LineStyle','--','color','k')

h(4) = plot(ECH_t,HeliconPower + ICHPower + ECHPower,'k')
legend(h,'Helicon','ICH','ECH','Combined')
set(h(4),'LineWidth',3)
xlim([4.15,4.32])
box on
ylim([0,250])
ylabel('[kW]')
xlabel('time [s]')
set(gcf,'color','w')
title('Shot: 16430')