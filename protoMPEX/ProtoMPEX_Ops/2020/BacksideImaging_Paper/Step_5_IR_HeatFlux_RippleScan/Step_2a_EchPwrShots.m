% Step 2a: Determine the ECH power and thus the shot sequence

clear all
close all

FetchData = 0;
saveData  = 0;

if FetchData
    
    % Shot list for the ECH power scan on 2020_03_30:
    % =========================================================================
    shotlist = [29771,29775,29778];
    PS2      = [3370 ,2360 ,1630 ];

    % MDSplus address:
    % =========================================================================
    Stem = '\MPEX::TOP.'; Branch = 'MACHOPS1:'; RootAddress = [Stem,Branch];

    % Fetch data from server:
    % =====================================================================
    % ECH power trace:
    dataAddress{1} = [RootAddress,'PWR_28GHZ'];
    [ech,t_ech]   = my_mdsvalue_v2(shotlist,dataAddress(1));

    % RF power trace:
    dataAddress{2} = [RootAddress,'RF_FWD_PWR'];
    [rf,t_rf]   = my_mdsvalue_v2(shotlist,dataAddress(2));
    for ii = 1:numel(rf)
        rf_s{ii} = movmean(rf{ii},300);
    end
    
     % PS2 currents:
    dataAddress{3} = [RootAddress,'PS2_I'];
    [ps2,t_ps2]   = my_mdsvalue_v2(shotlist,dataAddress(3));

    % Detect peak ECH power level:
    % =========================================================================
    for ii = 1:length(shotlist)
        % Smooth data:
        ech_s{ii} = movmean(ech{ii},300);
        % Find average peak:
        [d1,d2] = max(ech_s{ii});
        ech_peak(ii) = mean(ech_s{ii}(d2-20:d2+20));
    end
    
    % Plot data:
    % =========================================================================
    figure('color','w')
    hold on
    for ii = 1:length(shotlist)
        plot3(t_ech{ii}(1:end-1),PS2(ii)*ones(size(ech_s{ii})),70*ech_s{ii})
    end
    view([60,30])
    xlabel('time [s]')
    ylabel('index')
    zlabel('28 GHz power [kW]')

    % Save data:
    % =========================================================================
    varList = {'shotlist','PS2',...
               'ech','t_ech','rf','t_rf','ps2','t_ps2',...
               'rf_s','ech_s'};

    if saveData
        fileName = 'Step_2a_EchPwr';
        save(fileName,varList{:})
    end

else
    fileName = 'Step_2a_EchPwr';
    load(fileName)
end

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
hP(1) = plot(t_rf{2}(1:end-1),98*rf_s{2}/max(rf_s{2}),'k','LineWidth',2);
hP(2) = plot(t_ech{1}(1:end-1),70*ech_s{1},'bl','LineWidth',2);
hP(3) = plot(t_ech{2}(1:end-1),70*ech_s{2},'r','LineWidth',2);
hP(4) = plot(t_ech{3}(1:end-1),70*ech_s{3},'g','LineWidth',2);

% Legend text:
legendText{1} = '13.56 MHz';
legendText{2} = '28 GHz, PS2 = 3370 A';
legendText{3} = '28 GHz, PS2 = 2360 A';
legendText{4} = '28 GHz, PS2 = 1630 A';

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
figureName = 'Step_2a_RfPowerTraces';
saveas(gcf,figureName,'tiffn')