% This script reads the data acquired with the RF phase detector installed
% in ProtoMPEX on Sept 16th 2016

% #########################################################################
% THIS SCRIPT USED DURING XP TO MONITOR THE PROGRESS OF THE DATA
% ACQUISITION:
% #########################################################################

close all
clear all

% shot = 1e4 + 300 + [72,73];
% shot = 1e4 + 300 + [114:119];
shot = 1e4 + 700 + [35,36];

Data{1} = ['\MPEX::TOP.MACHOPS1:INT_4MM_1'];
Data{2} = ['\MPEX::TOP.MACHOPS1:INT_4MM_2'];
Data{3} = ['\MPEX::TOP.MACHOPS1:INT_2MM_1'];
Data{4} = ['\MPEX::TOP.MACHOPS1:INT_2MM_2'];

[VmagRatio,P,t] = BdotPhaseDetect_v1(shot,Data);


Title{1} = 'Vp0';
Title{2} = 'Vp90';
Title{3} = 'Vmag0';
Title{4} = 'None';

figure;
for ch = 1:4 % For all channels on digitizer
    [f{ch},t{ch}] = my_mdsvalue_v2(shot,Data(ch));
    % f{channel}{shot}
    subplot(2,2,ch);
    hold on
    for s = 1:length(shot)
            rng = find(t{ch}{s}>= 4.14 & t{ch}{s}<= 4.33);
            t{ch}{s} = t{ch}{s}(rng);
            f{ch}{s} = f{ch}{s}(rng);
            plot(t{ch}{s},f{ch}{s})
    end
    ylim([0,2]);
    title(Title{ch})
end



% Convert phase signals into a single angle
ch = 1; % Vp0
s = 1 ; % Shot 1
for n = 1:length(f{ch}{s}) % For all points in a trace
    % Check the voltages of Vp0 and Vp90
     
end