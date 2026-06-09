clear all
close all

shotlist = 27892;
shotlist = [27890,27893,27895];
% shotlist = [27895,27896];

Stem = '\MPEX::TOP.'; Branch = 'MACHOPS1:'; RootAddress = [Stem,Branch];
DA{1} = [RootAddress,'PS1_I2']; % Isx
[I,t_I]   = my_mdsvalue_v2(shotlist,DA(1));

Stem = '\MPEX::TOP.'; Branch = 'MACHOPS1:'; RootAddress = [Stem,Branch];
DA{1} = [RootAddress,'PS3_V']; % Isx
[V,t_V]   = my_mdsvalue_v2(shotlist,DA(1));

Stem = '\MPEX::TOP.'; Branch = 'MACHOPS1:'; RootAddress = [Stem,Branch];
DA{1} = [RootAddress,'RF_FWD_PWR']; % Isx
[RF,t_RF]   = my_mdsvalue_v2(shotlist,DA(1));

figure;
hold on
for s = 1:length(shotlist)
    plot(t_I{s}(1:end-1),-I{s}*5)
    plot(t_V{s}(1:end-1),V{s})
    plot(t_RF{s}(1:end-1),RF{s})
    grid on
end

box on