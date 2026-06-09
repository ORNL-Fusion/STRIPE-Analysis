% Power scan
clc
clearvars
close all

% Load datasetTable data
load('dataset_12_postprocess.mat','datasetTable')

% =========================================================================
% Gather previosly calculated data:
tic

% Bottom view:
% -------------------------------------------------------------------------
% Group 1: Very first heat flux measurement, window limit
varNames = {'t_heatflux_tempProbe','heatflux_tempProbe','shotIndex'};
group{1} = load('dataset_2_postprocess.mat',varNames{:});

% Group 2: Window limit with MPEX limiter installed
varNames = {'t_heatflux_tempProbe','heatflux_tempProbe','shotIndex'};
group{2} = load('dataset_4_postprocess.mat',varNames{:});

% Group 3: MPEX limit
varNames = {'t_heatflux_tempProbe','heatflux_tempProbe','shotIndex'};
group{3} = load('dataset_5_postprocess.mat',varNames{:});
toc

% Middle view:
% -------------------------------------------------------------------------
% Group 4: Very first heat flux measurement, window limit
varNames = {'t_heatflux_tempProbe','heatflux_tempProbe','shotIndex'};
group{4} = load('dataset_3_postprocess.mat',varNames{:});

% Top view:
% -------------------------------------------------------------------------
% Group 5: Very first heat flux measurement, window limit
varNames = {'t_heatflux_tempProbe','heatflux_tempProbe','shotIndex'};
group{5} = load('dataset_1_postprocess.mat',varNames{:});


for ii = 1:5
    heatflux_peak{ii} = group{ii}.heatflux_tempProbe;
    t_heatflux_peak{ii} = group{ii}.t_heatflux_tempProbe;
    pulseLength{ii} = datasetTable.pulseLength(group{ii}.shotIndex); 
    rfPwrNet{ii} = datasetTable.rfPwrNet(group{ii}.shotIndex); 
    shots{ii} = datasetTable.shot(group{ii}.shotIndex);
    
    for s = 1:length(shots{ii})
    end
    
end

C = {'k','r','bl','g','m','c','k--','r--','bl--','g--','m--','c--','k:','r:','bl:','g:','m:','c:'};

disp('Plot all calculated heat fluxes')
for ii = 1:5
    figure('Tag','powerScan','color','w')
    figureName{ii} = ['heatflux_peak_vs_time_group_',num2str(ii)];
    hold on
    for s = 1:length(shots{ii})
        tdummy = t_heatflux_peak{ii}{s};
        
        ss_type = 2;
        
        switch ss_type
            case 1
                if pulseLength{ii}(s)<300
                    rng = find(tdummy > 0.22 & tdummy <0.265 );
                    fctr = 0.8;
                elseif pulseLength{ii}(s)>300
                    rng = find(tdummy > 0.45 & tdummy <0.52 );
                    fctr = 1;
                end
            case 2
                if pulseLength{ii}(s)<300
                    rng = find(tdummy > 0.21 & tdummy <0.25 );
                    fctr = 1;
                elseif pulseLength{ii}(s)>300
                    rng = find(tdummy > 0.21 & tdummy <0.25 );
                    fctr = 1;
                end
        end
        
        heatflux_peak_ss{ii}(s) = fctr*mean(heatflux_peak{ii}{s}(rng));
        hDummy1(s) = plot(t_heatflux_peak{ii}{s},heatflux_peak{ii}{s}*1e-3,C{s},'LineWidth',2);
        plot(t_heatflux_peak{ii}{s}(rng),heatflux_peak{ii}{s}(rng)*1e-3,C{s},'LineWidth',4);
        xlim([0,1])
        ylim([0,800])
        box on
        grid on
        legendText{s} = [num2str(shots{ii}(s)),', RF: ',num2str(rfPwrNet{ii}(s)),' [kW]'];
    end
    title(['Group ',num2str(ii)])
    xlabel('time [s]','Interpreter','latex','FontSize',12)
    ylabel('[kWm$^{-2}$]','Interpreter','latex','FontSize',12)
    set(gca,'FontName','times','FontSize',11)
    hDummy2 = legend(hDummy1,legendText);
    hDummy2.FontSize = 9;
    set(gca,'PlotBoxAspectRatio',[2 1 1])
    clear hDummy1 hDummy2 
end

figure('Tag','powerScan','color','w')
figureName{end+1} = ['heatflux_peak_vs_rfPower'];
hold on
h_dummy(1) = plot(rfPwrNet{1},heatflux_peak_ss{1}*1e-3,'ko','MarkerFaceColor','k','MarkerSize',9)
h_dummy(2) = plot(rfPwrNet{2},heatflux_peak_ss{2}*1e-3,'rsq','MarkerFaceColor','r','MarkerSize',9)
h_dummy(3) = plot(rfPwrNet{3},heatflux_peak_ss{3}*1e-3,'gsq','MarkerFaceColor','g','MarkerSize',9)
xlim([0,200])
ylim([0,800])
xlabel('RF power [kW]','Interpreter','latex','FontSize',12)
ylabel('[kWm$^{-2}$]','Interpreter','latex','FontSize',12)
box on 
grid on
set(gca,'PlotBoxAspectRatio',[2 1 1])
h_dummy2 = legend(h_dummy,'Group 1, Window limit','Group 2, Window limit','Group 3, MPEX limit');
set(h_dummy2,'interpreter','latex','FontSize',11,'Location','northwest')
clearvars h_dummy


figure('Tag','powerScan','color','w')
figureName{end+1} = ['heatflux_peak_vs_rfPower_differentViews'];
hold on
h_dummy(1) = plot(rfPwrNet{1},heatflux_peak_ss{1}*1e-3,'ko','MarkerFaceColor','k','MarkerSize',9)
h_dummy(2) = plot(rfPwrNet{4},heatflux_peak_ss{4}*1e-3,'rsq','MarkerFaceColor','r','MarkerSize',9)
h_dummy(3) = plot(rfPwrNet{5},heatflux_peak_ss{5}*1e-3,'gsq','MarkerFaceColor','g','MarkerSize',9)
xlim([0,200])
ylim([0,800])
xlabel('RF power [kW]','Interpreter','latex','FontSize',12)
ylabel('[kWm$^{-2}$]','Interpreter','latex','FontSize',12)
box on 
grid on
set(gca,'PlotBoxAspectRatio',[2 1 1])
h_dummy2 = legend(h_dummy,'Group 1, bottom view','Group 4, middle view','Group 5, top view');
set(h_dummy2,'interpreter','latex','FontSize',11,'Location','northwest')
clearvars h_dummy

% =========================================================================
% Save figures from section4
disp('Save prompt...')
InputStructure.prompt = {['Would you like to save figures? Yes [1], No [0]']};
InputStructure.option.WindowStyle = 'normal';
beep
svfig = GetUserInput(InputStructure);
if svfig
    folderName = ['PowerScan_figures'];
    if exist(folderName) == 0
        mkdir(folderName);
    end
    h_dummy = findobj('Tag','powerScan');
    address_home = cd;
    for ii = 1:length(h_dummy)
        figure(h_dummy(ii))
        cd([address_home,'\',folderName])
        saveas(gcf,figureName{ii},'tiffn')
        cd(address_home)
    end
    clearvars h_dummy
end

