% =========================================================================
%   TESTING PHASE DETECTION IN PROTO-MPEX RF HELICON SYSTEM
%
% We are testing the new 2pi phase detector circuit made for recording the
% RF data into the MDSplus treee
% =========================================================================

clear all
close all

shot = 25600 + [39];
shot = 25600 + [40]; % we are now using detector A to get the phase between
% the RF fwd power from the mathcing box and the wavetek generator

shot = 25600 + [41];



Stem = '\MPEX::TOP.'; Branch = 'MACHOPS1:'; RootAddress = [Stem,Branch];
Data{1} = ['\MPEX::TOP.MACHOPS1:DET_PHASE2_1']; % Vp0
Data{2} = ['\MPEX::TOP.MACHOPS1:DET_PHASE2_2']; % Vp90
Data{3} = ['\MPEX::TOP.MACHOPS1:DET_PHASE2_2']; % Vp90

[f00,t00] = my_mdsvalue_v2(shot,Data(1));
[f90,t90] = my_mdsvalue_v2(shot,Data(2));

figure; 
hold on
plot(t00{1}(1:end-1),f00{1},'k')
plot(t90{1}(1:end-1),f90{1},'r')
ylim([0,2])

[VmagRatio,P,t_P] = BdotPhaseDetect_v3(shot,Data(1:3));

figure;
plot(t_P{1}(1:end-1),unwrap(P{1}*pi/180)*180/pi,'k.')
