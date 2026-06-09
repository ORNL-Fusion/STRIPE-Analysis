% Read and plot the ion saturation current from ion flux probes

close all 
clear all

shotlist = 15000 + [814]; % -10 V

% 2017_08_02
shotlist = 15000 + [944,945,946]; % -5 V
shotlist = 15000 + [947,948]; % -2 V
shotlist = 15000 + [949,950]; % -8 V
shotlist = 15000 + [951,952]; % -12 V
shotlist = 15000 + [953,954]; % -14 V
shotlist = 15000 + [955,956]; % -16 V
shotlist = 15000 + [965,966]; % -10 V
shotlist = 15000 + [967]; % -16 V

R = 17; 

fNames = {'IF','none'}; gNames = {'none','none'};

Stem = '\MPEX::TOP.'; Branch = 'MACHOPS1:'; RootAddress = [Stem,Branch];
DA{1} = [RootAddress,'INT_2MM_1']; 
[fA,t_fA]   = my_mdsvalue_v2(shotlist,DA(1));
DA{1} = [RootAddress,'INT_2MM_2']; 
[fB,t_fB]   = my_mdsvalue_v2(shotlist,DA(1));

DA{1} = [RootAddress,'INT_4MM_1']; 
[gA,t_gA]   = my_mdsvalue_v2(shotlist,DA(1));
DA{1} = [RootAddress,'INT_4MM_2']; 
[gB,t_gB]   = my_mdsvalue_v2(shotlist,DA(1));


%%
C = {'k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:'};
figure; hold on;
for s = 1:length(shotlist)
    IA{s} = fA{s}/R;
    IB{s} = fB{s}/R;

    hA(s) = plot(t_fA{s}(1:end-1),1000*IA{s},'k');
    hB(s) = plot(t_fB{s}(1:end-1),1000*IB{s},'r'); 

    set(hA(s),'LineWidth',1)
end
title(num2str(shotlist))
legend([hA,hB],fNames)
ylabel('[mA]')
xlim([4.12,4.5])
ylim([-10,200])
set(gcf,'color','w')
box on

% figure; hold on;
% for s = 1:length(shotlist)
%     GA{s} = gA{s}/R;
%     GB{s} = gB{s}/R;
%     
%     pA(s) = plot(t_gA{s}(1:end-1),1000*GA{s},'k'); 
%     pB(s) = plot(t_gB{s}(1:end-1),1000*GB{s},'r'); 
%     set(pA(s),'LineWidth',1)
% end
% title(num2str(shotlist))
% legend([pA,pB],gNames)
% ylabel('[mA]')
% xlim([4.12,4.5])
% ylim([-10,50])
% set(gcf,'color','w')
% box on
% 
% % to convert to flux use F = I