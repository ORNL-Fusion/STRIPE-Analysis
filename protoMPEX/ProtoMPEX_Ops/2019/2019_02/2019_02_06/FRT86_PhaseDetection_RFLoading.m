% =========================================================================
%   TESTING PHASE DETECTION IN PROTO-MPEX RF HELICON SYSTEM
%
% We are testing the new 2pi phase detector circuit made for recording the
% RF data into the MDSplus treee
% 2019_02_07
% =========================================================================

clear all
close all

% =========================================================================
% Select shots:
% Only 1 at a time
% =========================================================================

shot = 25600 + [93];

% =========================================================================
% Define Data addresses:
% =========================================================================
Stem = '\MPEX::TOP.'; 
Branch = 'MACHOPS1:'; 
RootAddress = [Stem,Branch];
switch shot
    case {25641,25643}
    % Phase between RF fwd match box and HP signal generator
    % 1 and 7.5 kHz
        Data{1} = ['\MPEX::TOP.MACHOPS1:DET_PHASE2_1']; % Vp0
        Data{2} = ['\MPEX::TOP.MACHOPS1:DET_PHASE2_2']; % Vp90
    case {25644,25646,25648,25658,25671,25674,25682,25687,25693}
        % Phase between RF fwd match box and HP signal generator
        % 50 kHz 
        % 48 is a good shot to show the phase difference
        % 58 low RF power
        Data{1} = ['\MPEX::TOP.MACHOPS1:BDOT_1']; % Vp0
        Data{2} = ['\MPEX::TOP.MACHOPS1:BDOT_2']; % Vp90
end

% Fwd, Ref from XT side and antenna current from RF mathing box
Data{3} = ['\MPEX::TOP.MACHOPS1:RF_FWD_PWR']; % Vp0
Data{4} = ['\MPEX::TOP.MACHOPS1:RF_REF_PWR']; % Vp90
Data{5} = ['\MPEX::TOP.MACHOPS1:RF_CURRENT']; % Vp90

% =========================================================================
% Import data from MDSplus
% =========================================================================
mdsconnect('mpexserver');
[~,~]    = mdsopen('MPEX',shot);
[f00,~]  = mdsvalue(Data(1));
[t00,~]  = mdsvalue(['DIM_OF(',cell2mat(Data(1)),')']);

[f90,~]  = mdsvalue(Data(2));
[t90,~]  = mdsvalue(['DIM_OF(',cell2mat(Data(2)),')']);

[fwd,~]    = mdsvalue(Data(3));
[t_fwd,~]  = mdsvalue(['DIM_OF(',cell2mat(Data(3)),')']);

[ref,~]    = mdsvalue(Data(4));
[t_ref,~]  = mdsvalue(['DIM_OF(',cell2mat(Data(4)),')']);

[I,~]    = mdsvalue(Data(5));
[t_I,~]  = mdsvalue(['DIM_OF(',cell2mat(Data(5)),')']);

a = sgolay_t(f00,3,71);
b = sgolay_t(f90,3,71);

PreviewRawData = 1;
if PreviewRawData
    figure; 
    hold on
    plot(t00(1:end-1),a,'k')
    plot(t90(1:end-1),b,'r')
    ylim([0,2])
end

% =========================================================================
% Extract phase difference from raw data:
% =========================================================================
[P] = RFPhaseDetect(a,f90);
t_P = t00;

%% Plot data
tStartPlot = 4.0;
tEndPlot = 4.8;

figure;
% Plot phase difference
subplot(1,2,1);
hold on
h(1) = plot(t_P(1:end-1),unwrap(P*pi/180)*180/pi,'k.');
xlim([tStartPlot,tEndPlot])
ylim([-360,360])
box on
grid on
legend(h,['shot: ',num2str(shot)])
title('\Delta\phi, fwd and RF drive')
xlabel('time [sec]')
ylabel('[deg]')

% Plot fwd, ref and antenna current
subplot(1,2,2);
hold on
h(1) = plot(t_fwd(1:end-1),fwd-mean(fwd(1:100)),'k');
h(2) = plot(t_ref(1:end-1),ref-mean(ref(1:100)),'r');
h(3) = plot(t_I(1:end-1),5*(I-mean(I(1:100))),'g');
xlim([tStartPlot,tEndPlot])
ylim([0,1.5])
box on
grid on
legend(h,'RF FWD','RF REF','RF CURRENT')
title(['shot: ',num2str(shot)])
xlabel('time [sec]')
ylabel('[a.u.]')
set(gcf,'color','w')