% revised_IEAD_sheathEdge_energyAngleImpact.m
clc; clear; close all;

%% ----------------------------
% USER SETTINGS
%% ----------------------------
tilt = 0;
file = sprintf('../tilted_targets/test/%d_degrees/MPEX_runs_test_noCollisions/energyAngleImpact.nc', tilt);



%% ----------------------------
% READ DATA (sheath-edge IEAD)
%% ----------------------------
W = ncread(file,'Weight');  % Size: [nAngles x nEnergies]



figure; pcolor(W); shading interp