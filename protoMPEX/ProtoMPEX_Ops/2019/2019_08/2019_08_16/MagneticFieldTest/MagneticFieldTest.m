% Magnetic field testing on Proto-MPEX
clear all
close all

shotlist = 27000 + [509,510,511,512];
Scale    =         [1e3,1e3,1e3,1e3];

Stem = '\MPEX::TOP.'; 
Branch = 'MACHOPS1:'; 
RootAddress = [Stem,Branch];
% Gaussmeter measurements
DA{1} = [RootAddress,'PS1_I2'];
% TR2 current, COIL 2
DA{2} = [RootAddress,'TR2_I'];
% PS1 current, COIL 1,6,7 and 8
DA{3} = [RootAddress,'PS1_I'];

[RawGaussMeter,t_RawGaussMeter]   = my_mdsvalue_v2(shotlist,DA(1));
[RawTR2,t_RawTR2]   = my_mdsvalue_v2(shotlist,DA(2));
[RawPS1,t_RawPS1]   = my_mdsvalue_v2(shotlist,DA(3));

% Convert raw data into magnetic field in Tesla
for i = 1:numel(shotlist)
    Bmeasured{i}   = -RawGaussMeter{i}*Scale(i);
    t_Bmeasured{i} = t_RawGaussMeter{i}(1:end-1);
    PS1{i}         = RawPS1{i}*1e3;
    t_PS1{i}       = t_RawPS1{i}(1:end-1);
    TR2{i}         = RawTR2{i}*1e3;
    t_TR2{i}       = t_RawTR2{i}(1:end-1);
end

figure
hold on;
for i = 1:numel(shotlist)
    plot(t_Bmeasured{i},Bmeasured{i})
end

figure
hold on;
for i = 1:numel(shotlist)
    plot(t_TR2{i},TR2{i})
end