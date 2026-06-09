% Read and plot the ion saturation current from ion flux probes

close all 
clear all

shotlist = 15000 + [770];
fNames = {'3A','3B'}; gNames = {'IF','none'};

shotlist = 15000 + [771];
fNames = {'IF','3B'}; gNames = {'3A','none'};

shotlist = 15000 + [776];
fNames = {'3A','3B'}; gNames = {'4A','4B'};

shotlist = 15000 + [777];
fNames = {'3A','3B'}; gNames = {'4B','4A'};

shotlist = 15000 + [778];
fNames = {'3A','3B'}; gNames = {'4A','NONE'};

shotlist = 15000 + [779];
fNames = {'3A','3B'}; gNames = {'none','4B'};

shotlist = 15000 + [782];
fNames = {'3A','3B'}; gNames = {'2A','2B'};

shotlist = 15000 + [784];
fNames = {'3A','3B'}; gNames = {'2A','2B'};

shotlist = 15000 + [786];
fNames = {'2A','2B'}; gNames = {'4A','4B'};


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
    IA{s} = -fA{s}/10;
    IB{s} = -fB{s}/10;

    hA(s) = plot(t_fA{s}(1:end-1),1000*IA{s},'k');
    hB(s) = plot(t_fB{s}(1:end-1),1000*IB{s},'r'); 

    set(hA(s),'LineWidth',1)
end
title(num2str(shotlist))
legend([hA,hB],fNames)
ylabel('[mA]')
xlim([4.12,4.5])
ylim([-10,50])
set(gcf,'color','w')
box on

figure; hold on;
for s = 1:length(shotlist)
    GA{s} = -gA{s}/10;
    GB{s} = -gB{s}/10;
    
    pA(s) = plot(t_gA{s}(1:end-1),1000*GA{s},'k'); 
    pB(s) = plot(t_gB{s}(1:end-1),1000*GB{s},'r'); 
    set(pA(s),'LineWidth',1)
end
title(num2str(shotlist))
legend([pA,pB],gNames)
ylabel('[mA]')
xlim([4.12,4.5])
ylim([-10,50])
set(gcf,'color','w')
box on

% to convert to flux use F = I