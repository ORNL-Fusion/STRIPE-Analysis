% Task4
% load and save the FP, gas and DLP data for all datasets

clc
close all
disp('Start of task 3 ##################################################');

%% SECTION 1
% =========================================================================
% Compute prompt:
InputStructure.prompt = {['Would you like to reload FP data from MPEX server? Yes [1], No [0]']};
beep
cmpt = GetUserInput(InputStructure);

defaultAns = 0;
tic
for datasetToAnalyze = 1:8;
% =========================================================================
% Clear memory
if cmpt
    clearvars('-except','cmpt','datasetToAnalyze','defaultAns');
end

close all

% =========================================================================
% Input required from user:
% datasetToAnalyze = 8;

% #########################################################################
% SECTION 1:
sectionName = 'Load FP data from MPEX server';
disp(['Section: ',sectionName])
% #########################################################################
    
% =========================================================================
% Load FP data MPEX server:
if cmpt
    disp('Loading data from MPEX server...')
    dum1 = who;
    
    % =========================================================================
    % Read dataset spreadsheet:
    spreadsheetName = 'HeliconWindowIR_2020_02_14.xlsx';
    datasetTable = readtable(spreadsheetName,'Sheet',1);

    % =========================================================================
    % Identify shots from the dataset table:
    [shotIndex,~] = find(datasetTable.dataset == datasetToAnalyze);
    shot = datasetTable.shot(shotIndex);
    X    = datasetTable.X(shotIndex);
    nshots = numel(shot);
    
    % ---------------------------------------------------------------------
    % MPEXserver data address:
    Stem = '\MPEX::TOP.'; Branch = 'MACHOPS1:'; RootAddress = [Stem,Branch];
    
    % ---------------------------------------------------------------------
    % FP data:
    removeOffset = 1;
    % Gather data
    addressFP{1} = [RootAddress,'FLUOROPT_1']; 
    addressFP{2} = [RootAddress,'FLUOROPT_2']; 
    addressFP{3} = [RootAddress,'FLUOROPT_3']; 
    addressFP{4} = [RootAddress,'FLUOROPT_4']; 

    for ii = 1:4
        [dataFP_raw{ii},t_dataFP{ii}]   = my_mdsvalue_v2(shot,addressFP(ii));
    end

    for ii = 1:4
        for s = 1:length(shot)
            dataFP{ii}{s} = sgolay_t( ( dataFP_raw{ii}{s}- removeOffset*mean(dataFP_raw{ii}{s}(1:200)) )*20 ,3,51);
            dataFP_ss{ii}(s) = mean(dataFP{ii}{s}( (end-300):end) );
        end
    end
    clearvars ii s
    vars2save = setdiff(who,dum1);
else
    disp('Loading previously saved data ...')
    fileName = ['dataset_',num2str(datasetToAnalyze),'_DlpGasFp.mat'];
    load(fileName)
    disp('Data succesfully loaded!')
end
    
% ---------------------------------------------------------------------
% Plot and check data:
figure
set(gcf,'Tag','FP vs time','color','w')
hold on

for ii = 1:4
    subplot(2,2,ii); hold on
    for s = 1:length(shot)
        h_fp{ii}(s) = plot(t_dataFP{ii}{s}(1:end-1),dataFP{ii}{s},'LineWidth',3);
    end
    legend(h_fp{ii},num2str(shot))
    ylim([0,25])
    title(num2str(addressFP{ii}(numel(RootAddress)+1:end)))
end
clearvars ii s h*

if strcmp(datasetTable.scanType(shotIndex(1)),'RF power')
    figure('Tag','FP vs RfPwr','color','w')
    for ii = 1:4
        subplot(2,2,ii)
        hold on
        pulseLength = datasetTable.pulseLength(shotIndex);
        b = find(pulseLength>300);
        plot(X(b),dataFP_ss{ii}(b),'ro-')
        b = find(pulseLength<300);
        plot(X(b),dataFP_ss{ii}(b),'ko-')
        xlim([0,200])
        ylim([0,35])
        ylabel('$\Delta T$ [C]','Interpreter','latex','FontSize',12)
        box on
        grid on
        title(num2str(addressFP{ii}(numel(RootAddress)+1:end)))
    end
    saveFig = 1;
    if saveFig
        disp('Saving figures...')
        hdum = findobj('Tag','FP vs RfPwr');
        figureName = ['dataset_',num2str(datasetToAnalyze),'_FP_vs_rfPwr'];
        saveas(hdum,figureName,'tiffn')
       disp('Figure saved!')
    end
end
    
if cmpt
    % ---------------------------------------------------------------------
    if defaultAns ~= 1
        % Save prompt:
        disp('Save prompt...')
        InputStructure.prompt = {['Would you like to save? Yes [1], No [0]']};
        InputStructure.option.WindowStyle = 'normal';
        beep
        svdt = GetUserInput(InputStructure);
    else
        svdt = 1;
    end
    if svdt
        disp('Saving FP data ...')
        variableNames = vars2save;
        fileName = ['dataset_',num2str(datasetToAnalyze),'_DlpGasFp.mat'];
        SaveData;
    end
end

% =========================================================================
disp('***********************Section 1 completed*************************')

%% SECTION 2
% =========================================================================
% Compute prompt:
if defaultAns ~= 1
    InputStructure.prompt = {['Would you like to reload Gas data from MPEX server? Yes [1], No [0]']};
    beep
    cmpt = GetUserInput(InputStructure);
else
    cmpt = 1;
end

% #########################################################################
% SECTION 2:
sectionName = 'Load Gas data from MPEX server';
disp(['Section: ',sectionName])
% #########################################################################
    
% =========================================================================
% Load DLP data MPEX server:
if cmpt
    disp('Loading data from MPEX server...')
    dum1 = who;

    % ---------------------------------------------------------------------
    % MPEXserver DLP data address:
    Stem = '\MPEX::TOP.'; Branch = 'MACHOPS1:'; RootAddress = [Stem,Branch];

    % ---------------------------------------------------------------------
    % Gas data:
    addressGF{1} = [RootAddress,'MFC_FLOW_D2'];
    addressGF{2} = [RootAddress,'MFC_FLOW_HE']; 

    for ii = 1:2
        [dataGF{ii},t_dataGF{ii}]   = my_mdsvalue_v2(shot,addressGF(ii));
    end

    addressGP{1} = [RootAddress,'PG1']; 
    addressGP{2} = [RootAddress,'PG2']; 
    addressGP{3} = [RootAddress,'PG3']; 
    addressGP{4} = [RootAddress,'PG4']; 

    shot_cal = 20908; % 2018_04_05, 2.0 kA case, 1 sec RF pulse
    fctr = 0.1333;
    gainPG = [2,2,2,10];
    pressureGaugeLabel = {'12.5','8.5','6.5','2.5'};

    for ii = 1:4
        [dataGP_raw{ii}    ,t_dataGP_raw{ii}    ]   = my_mdsvalue_v2(shot,addressGP(ii));
        [dataGP_cal{ii},t_dataGP_cal{ii}]   = my_mdsvalue_v2(shot_cal,addressGP(ii));

        for s = 1:length(shot)
                L1 = length(dataGP_raw{ii}{s}    );
                L2 = length(dataGP_cal{ii}{1});
                if L1>L2
                    Ldata = L2;
                else
                    Ldata = L1;
                end
            dataGP{ii}{s} = ( dataGP_raw{ii}{s}(1:Ldata)-dataGP_cal{ii}{1}(1:Ldata) )*gainPG(ii)*fctr;
            t_dataGP{ii}{s} = t_dataGP_raw{ii}{s}(1:Ldata);
        end
    end

    clearvars ii s *dataGP_raw *dataGP_cal
    vars2save = setdiff(who,dum1);
else
    disp('Loading previously saved data ...')
    fileName = ['dataset_',num2str(datasetToAnalyze),'_DlpGasFp.mat'];
    load(fileName)
    disp('Data succesfully loaded!')
end

% ---------------------------------------------------------------------
% Plot and check data:

s = 2;
figure;
set(gcf,'Tag','GP','color','w')
subplot(2,1,1)
hold on
for ii = 1:2
    plot(t_dataGF{ii}{s}(1:end-1),dataGF{ii}{s},'LineWidth',3)
end
xlim([3.8,5.5])
ylim([0,8])
ylabel('Piezo voltage [V]','Interpreter','latex','Fontsize',12)
box on
grid on
% h_t = text(5.08,4.1,{magConfig{xp};['shot: ',num2str(shot(s))]});
% set(h_t,'interpreter','Latex','fontsize',13)

subplot(2,1,2)
hold on
for ii = 1:4
    t_rng = find(t_dataGP{ii}{s}>3.9 & t_dataGP{ii}{s}<5.4);
    plot(t_dataGP{ii}{s}(t_rng),dataGP{ii}{s}(t_rng),'LineWidth',3)
end
xlim([3.8,5.5])
ylim([0,1.5])
ylabel('Pressure [Pa]','Interpreter','latex','Fontsize',12)
box on
grid on

% =========================================================================
if cmpt
    if defaultAns ~= 1
        % Save prompt:
        disp('Save prompt...')
        InputStructure.prompt = {['Would you like to save? Yes [1], No [0]']};
        InputStructure.option.WindowStyle = 'normal';
        beep
        svdt = GetUserInput(InputStructure);
    else
        svdt =1;
    end
    if svdt
        disp('Saving FP data ...')
        variableNames = vars2save;
        fileName = ['dataset_',num2str(datasetToAnalyze),'_DlpGasFp.mat'];
        SaveData;
    end
end

% =========================================================================
disp('***********************Section 2 completed*************************')
end
toc
% =========================================================================
% End of script
disp('End of script!!')



% notes:
return

if cmpt
    % var snapshot
    % load data from server
    % var snapshot
    % create vars2save
else
    % load previously saved data
end

% plot data to check

if cmpt
    % save prompt
    if svdt
        % save data from vars2save
    end
end