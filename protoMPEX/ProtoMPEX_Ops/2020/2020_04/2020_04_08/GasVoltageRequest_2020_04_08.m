% Plot the gas request voltage
% clear all
% close all

% 2020_03_06
% =========================================================================
% Shot to replicate:
shotlist = [27000 + [187]];
shotlist = [29000 + [731,765]];
shotlist = [29000 + [766]]; 

Stem = '\MPEX::TOP.'; Branch = 'MACHOPS1:'; RootAddress = [Stem,Branch];
addressGasFlow{1} = [RootAddress,'MFC_FLOW_D2']; 
addressGasFlow{2} = [RootAddress,'MFC_FLOW_HE']; 

for ii = 1:2
    [dataGasFlow{ii},t_dataGasFlow{ii}]   = my_mdsvalue_v2(shotlist,addressGasFlow(ii));
end

figure;
for s = 1:length(shotlist)
    subplot(length(shotlist),1,s)
    hold on
    h(1) = plot(t_dataGasFlow{1}{s}(1:end-1),dataGasFlow{1}{s},'k')
    h(2) = plot(t_dataGasFlow{2}{s}(1:end-1),dataGasFlow{2}{s},'r')
    legend(h,'Piezo 2.5','Piezo 1.5')
    title(['shot: ',num2str(shotlist(s))])
    xlim([3.5,5.5])
    ylim([0,10])
end