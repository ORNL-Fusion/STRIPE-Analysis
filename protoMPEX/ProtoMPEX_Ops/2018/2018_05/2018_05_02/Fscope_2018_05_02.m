% Collect Filtescope data
close all
% clear all

Stem = '\MPEX::TOP.'; Branch = 'FSCOPE:'; RootAddress = [Stem,Branch];

DA{1} = [RootAddress,'TUBE16:PMT_VOLT'];
[f16,t_f16]   = my_mdsvalue_v2(shotlist,DA(1));

DA{1} = [RootAddress,'TUBE24:PMT_VOLT'];
[f24,t_f24]   = my_mdsvalue_v2(shotlist,DA(1));


C = {'k','r','bl','m','g','c'};

figure;
subplot(2,1,1)
hold on
for s = 1:length(shotlist)
    h16(s) = plot(t_f16{s},f16{s},C{s});
end
box on
legend(h16,{num2str(shotlist')},'location','NorthWest')

subplot(2,1,2)
hold on
for s = 1:length(shotlist)
    h24(s) = plot(t_f24{s},f24{s},C{s});
end
box on
legend(h24,{num2str(shotlist')},'location','NorthWest')

set(gcf,'color','w')