% Task 3: 
clc
close all
clearvars

disp('Start of task 3 ##################################################');

% #########################################################################
% SECTION 1:
sectionName = 'Load data from Task 2';
disp(['Section: ',sectionName])
% #########################################################################

% =========================================================================
% Load data from .mat files:
% ---------------------------------------------------------------------
% Input required from user:
datasetToAnalyze = 12;
% ---------------------------------------------------------------------
% Load data:
disp('Loading Task 2 data from .mat file...')
fileName = ['dataset_',num2str(datasetToAnalyze),'_postprocess.mat'];
tic
load(fileName,'tempProbe_mean','t_tempProbe_mean','shots','nshots','datasetTable','shotIndex')
disp('Data loaded!')
toc

% =========================================================================
% Define variables specific to this section:
vars_task3.section1 = '';
vars_task3.section1 = who;
vars_task3.section1 = setdiff(vars_task3.section1,[]);

% =========================================================================
disp('***********************Section 1 completed*************************')

%% SECTION 2
% =========================================================================
% Clear memory and figures:
try
    clearvars('-except',vars_task3_section1{:});
    close(findobj('Tag',sectionName))
end

% It occurred to me that instead of using the -except flag, we can instead
% use try-catch to extract vars.task3.section2 from dataset_X_postprocess
% and then use that variable to directly clear all section2 variables

% #########################################################################
% SECTION 2:
sectionName = 'Convert temperature probes into heat flux';
disp(['Section: ',sectionName])
% #########################################################################

% =========================================================================
% Compute prompt:
InputStructure.prompt = {['Would you like to recompute the heat flux? Yes [1], No [0]']};
cmpt = GetUserInput(InputStructure);

if cmpt
    disp('Computing...')
% =========================================================================
    % ALN material properties
    rho = 3300; 
    kt   = 180 ;
    cp  = 740 ;
    a = kt/(rho*cp); % thermal diffusivity

% =========================================================================
    % Geometry
    Lz = 6/1000; % ALN wall thickness in ProtoMPEX

 % =========================================================================
    % Dimensionless and normalized data
    % Order of magnitude heat flux:
    q0 = 1e6; 
    % Characteristic time scale:
    t_star = Lz*Lz/a;
    % Characteristic temperature:
    T_star = q0*Lz/kt;
    % Normalized measured temperature and time:
    for s = 1:nshots
        z{s}    = tempProbe_mean{s}'/T_star; 
        t_z{s}  = t_tempProbe_mean{s}'/t_star;
        dt_z(s) = t_z{s}(2) - t_z{s}(1);
    end

 % =========================================================================
    % Create the convolution operator
    xx0 = 0.0005;
    % To capture the details of the front surface we need to use 1e4
    % partial sum terms and use xx0 = 0.01 to 0.05

    % Partial sum
    % "n" is the index
    disp('Computing Toeplitz matrix...')
    tic
    Ns = 1e4;
    for s = 1:nshots
         for fr = 1:length(t_z{s})
               K0{s}(fr) = G_Impulse_1D(xx0,t_z{s}(fr),Ns);
         end
         P{s} = toeplitz(K0{s},zeros(size(K0{s})));
    end
    toc
    disp('Computation complete')

 % =========================================================================
    %  Conjugate gradient method:
    % Input data:
    % P: convolution operator
    % q{1}: initial guess at heat flux, column vector
    % z: experimental data, column vector
    % Ni: number of iteration
    % Output data:
    % u: calculated temperature, column vector
    % J: residual, column vector
    % q: minimized heat flux, structure, column vector
    disp('Applying inverse method...')
    tic
    for s = 1:nshots
        qdummy{1} = zeros(size(t_z{s}));
        [a1,a2,a3] = IHCP_ConjugateGradient(P{s},qdummy,z{s},dt_z(s),1900,3,[]);
        q{s} = a1{end};
        u{s} = a2{end};
        J{s} = a3;

        heatflux_tempProbe{s} = q0*q{s};
        t_heatflux_tempProbe{s} = t_tempProbe_mean{s};
    end
    toc
    disp('Computation completed!')

% =========================================================================
    % Infinite slab method:
    disp('Applying infinite slab method...')
    tic
    for s = 1:nshots
        dt = t_tempProbe_mean{s}(2) - t_tempProbe_mean{s}(1);
        var_dummy = real(sqrt(0.25*pi*kt*rho*cp*diff(tempProbe_mean{s}.^2)/dt));
        heatflux_tempProbe_infiniteSlab{s} = [var_dummy,var_dummy(end)];
    end
    clear var_dummy
    toc
    disp('Computation complete!')
end

% =========================================================================
if cmpt == 0
    % Load previosly calculated data
    disp('Loading previously calculated data...')
    fileName = ['dataset_',num2str(datasetToAnalyze),'_postprocess.mat'];
    variableNames = {'heatflux_tempProbe','t_heatflux_tempProbe','heatflux_tempProbe_infiniteSlab',...
            't_z','z','u','t_star','T_star'};
    load(fileName,variableNames{:});
end

% =========================================================================
% Plotting solution:
switch nshots
    case 1
        nCol = 1;
        nRow = 1;
    case 2
        nCol = 2;
        nRow = 1;
    case {3,4}
        nCol = 2;
        nRow = 2;
    case {5,6}
        nCol = 3;
        nRow = 2;
    case {7,8,9}
        nCol = 3;
        nRow = 3;
    case {10,11,12}
        nCol = 4;
        nRow = 3;
    case {13,14,15,16}
        nCol = 4;
        nRow = 4;
end
lineColor = {'k','r','bl','g','m','c','k--','r--','bl--','g--','m--','c--','k:','r:','bl:','g:','m:','c:'};

disp('Testing inverse solution...')
figure('Tag','section2','color','w')
for s = 1:nshots
    subplot(nRow,nCol,s)
    hold on
    hDummy1(1) = plot(t_z{s}*t_star,z{s}*T_star,'k.','MarkerSize',12);
    hDummy1(2) = plot(t_z{s}*t_star,u{s}*T_star,'r','LineWidth',1);
    xlim([0,1])
    ylim([0,15])
    box on
    grid on
    title(['Shot: ',num2str(shots(s))])
end
hDummy2 = legend(hDummy1,'Exp.','Reconstruction');
clear hDummy1 hDummy2 

disp('Comparing finite and infinite methods...')
figure('Tag','section2','color','w')
for s = 1:nshots
    subplot(nRow,nCol,s)
    hold on
    hDummy1(1) = plot(t_z{s}*t_star,heatflux_tempProbe{s}*1e-3             ,lineColor{s},'LineWidth',2);
    hDummy1(2) = plot(t_z{s}*t_star,heatflux_tempProbe_infiniteSlab{s}*1e-3,lineColor{s},'LineWidth',1);
    xlim([0,1])
    ylim([0,700])
    box on
    grid on
    title(['Shot: ',num2str(shots(s))])
end
hDummy2 = legend(hDummy1,'Finite','Infinite');
clear hDummy1 hDummy2 

% =========================================================================
if cmpt
    % Define variables specific to this section:
    vars_task3.section2 = '';
    vars_task3.section2 = who;
    vars_task3.section2 = setdiff(vars_task3.section2,vars_task3.section1);

    % Save prompt:
    disp('Save prompt')
    InputStructure.prompt = {['Would you like to save? Yes [1], No [0]']};
    InputStructure.option.WindowStyle = 'normal';
    svdt = GetUserInput(InputStructure);
    if svdt        
        variableNames = {'heatflux_tempProbe','t_heatflux_tempProbe','heatflux_tempProbe_infiniteSlab',...
            't_z','z','u','t_star','T_star'};
        fileName = ['dataset_',num2str(datasetToAnalyze),'_postprocess.mat'];
        SaveData
    end
end

% =========================================================================
% Plot all data:
disp('Plot all calculated heat fluxes')
figure('Tag','section2','color','w')
hold on
X = datasetTable.X(shotIndex);
scanType = datasetTable.scanType(shotIndex);
for s = 1:nshots
    hDummy1(s) = plot(t_z{s}*t_star,heatflux_tempProbe{s}*1e-3,lineColor{s},'LineWidth',2);
    xlim([0,1])
    ylim([0,800])
    box on
    grid on
    legendText{s} = [num2str(shots(s)),', ',scanType{1},': ',num2str(X(s))];
end
title(['Dataset ',num2str(datasetToAnalyze),', ',scanType{1}])
xlabel('time [s]','Interpreter','latex','FontSize',12)
ylabel('[kWm$^{-2}$]','Interpreter','latex','FontSize',12)
set(gca,'FontName','times','FontSize',11)
hDummy2 = legend(hDummy1,legendText);
set(gca,'PlotBoxAspectRatio',[2 1 1])
clear hDummy1 hDummy2 

% =========================================================================
% Save figures from section4
disp('Save prompt...')
InputStructure.prompt = {['Would you like to save figures from section 2? Yes [1], No [0]']};
InputStructure.option.WindowStyle = 'normal';
beep
svfig = GetUserInput(InputStructure);
if svfig
    folderName = ['dataset_',num2str(datasetToAnalyze),'_figures'];
    if exist(folderName) == 0
        mkdir(folderName);
    end
    h_dummy = gcf;
    address_home = cd;
    for ii = 1:length(h_dummy)
        figure(h_dummy(ii))
        figureName = ['ds_',num2str(datasetToAnalyze),'_peakHeatFlux'];
        cd([address_home,'\',folderName])
        saveas(gcf,figureName,'tiffn')
        cd(address_home)
    end
end

% =========================================================================
disp('***********************Section 2 completed*************************')

% =========================================================================
% End of script
disp('End of script!!')

%% SECTION 3

% Here we need plot and analyze data


% Consider saving figures
% Consider saving heat flux data to .mat file
