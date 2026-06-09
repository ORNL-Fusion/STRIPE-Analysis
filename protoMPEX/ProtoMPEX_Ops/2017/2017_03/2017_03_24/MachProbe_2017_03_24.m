% ==================
% This code is designed to retrive data from Mach Probe, plot current vs time
% from MP_n1 and MP_n2 and calculate the current ratio
% Isense A goes to MP tip 1 and Isense B goes to MP tip 2
% First written on Feb 13, 2017
% ==================

close all
clear all

shotlist = 13500+[09,11];

% -------------------------
Config.tStart = 4.18; % [s]
Config.tEnd = 4.36;

% Acquiring Isat current from MP tip 1 and tip 2
Stem = '\MPEX::TOP.';
Branch = 'MACHOPS1:';
RootAddress = [Stem,Branch];

DataAddress{1} = [RootAddress,'INT_2MM_1']; % I_senseA Channel 09
DataAddress{2} = [RootAddress,'INT_2MM_2']; % I_senseA Channel 10

% Mach Probe tip lengths

MP_type = 'MP1' % MP1 or MP2

switch MP_type
    case 'MP1'
        L_mp1 = 5.5; % [mm]
        L_mp2 = 5.5; %[mm]
    case 'MP2'
        L_mp1 = 4.0; % [mm]
        L_mp2 = 5.0; %[mm]
end



% Standard probe tip configuration
% --------------------------------
% MPI_1  = probe tip facing target 
% MPI_2 = probe tip facing dump

[MPI_1] = my_mdsvalue_v2(shotlist,DataAddress(1)); % [V] signal from digitizer from Tip A
[MPI_2] = my_mdsvalue_v2(shotlist,DataAddress(2)); % [V] signal from digitizer from Tip B

% Calibration using 10 Ohm resitor
R = 10.0; % Ohms

%%
figure; hold on
C = {'k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:'};
for s = 1:length(shotlist)
    plot(MPI_1{s}*1000/R, C{s})
    xlabel('time [s]')
    ylabel('Current Isense A [mA]')
    plot(MPI_2{s}*1000/R, C{s+2})
    xlabel('time [s]')
    ylabel('Current Isense A [mA]')
    ylim([-400 100])
    %xlabel('time [s]')
    %ylabel('Current Isense B [mA]')
end
return
ratio = MPI_2./MPI_1
M_number = log(ratio)

return


