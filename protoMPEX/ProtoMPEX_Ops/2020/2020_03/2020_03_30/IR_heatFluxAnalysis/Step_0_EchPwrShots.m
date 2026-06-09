% Step 0: Determine the ECH power and thus the shot sequence

clear all
close all

% Shot list for the ECH power scan on 2020_03_30:
% =========================================================================
shotlist = [29671,29669,29665,29668,29661,29663,29662,29673,29655,29657];
pwrECH   = [5.92 ,26.12,36.84,38.29,45.41,50.98,52.61,54.07,70.07,70.48];

% MDSplus address:
% =========================================================================
Stem = '\MPEX::TOP.'; Branch = 'MACHOPS1:'; RootAddress = [Stem,Branch];

% ECH power trace:
dataAddress{1} = [RootAddress,'PWR_28GHZ'];
[ech,t_ech]   = my_mdsvalue_v2(shotlist,dataAddress(1));

% RF power trace:
dataAddress{2} = [RootAddress,'RF_FWD_PWR'];
[rf,t_rf]   = my_mdsvalue_v2(shotlist(3),dataAddress(2));
rf_s{1} = movmean(rf{1},300);

% Sort shots:
% =========================================================================
for ii = 1:length(shotlist)
    % Smooth data:
    ech_s{ii} = movmean(ech{ii},300);
    % Find average peak:
    [d1,d2] = max(ech_s{ii});
    ech_peak(ii) = mean(ech_s{ii}(d2-20:d2+20));
end

[~,b] = sort(ech_peak);

figure('color','w')
plot(70*ech_peak(b))
ylabel('ECH power [kW]')
xlabel('index')
box on

% Plot data:
% =========================================================================
figure('color','w')
hold on
for ii = 1:length(shotlist)
    plot3(t_ech{b(ii)}(1:end-1),ii*ones(size(ech_s{(ii)})),70*ech_s{b(ii)})
end
view([60,30])
xlabel('time [s]')
ylabel('index')
zlabel('ECH power [kW]')

% Print order of shots:
% =========================================================================
shotlist_ordered = shotlist(b)
echPwr = round(ech_peak(b)*70,2)

%%
close all

% Plot data:
% =========================================================================
fontSize.label = 12;
fontSize.title = 13;
fontSize.legend = 11;
fontSize.axes = 11;

figure('color','w')
hold on
hP(1) = plot(t_rf{1}(1:end-1),100*rf_s{1},'k','LineWidth',2);
hP(2) = plot(t_ech{end}(1:end-1),70*ech_s{end},'r','LineWidth',2);

% Legend text:
legendText{1} = '13.56 MHz';
legendText{2} = '28 GHz';

% Formatting:
xlim([4.05,4.8])
ylim([0,150])
box on
set(gca,'fontSize',fontSize.axes,'FontName','Times')

% Labels:
ylabel('[kW]','interpreter','Latex','FontSize',fontSize.label)
xlabel('time [s]','interpreter','Latex','FontSize',fontSize.label)

% Legend:
hL = legend(hP,legendText);
set(hL,'fontSize',fontSize.legend)

% Save figure
% =========================================================================
figureName = 'Step_0_RfPowerTraces';
saveas(gcf,figureName,'tiffn')