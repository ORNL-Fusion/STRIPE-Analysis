% Task 1: Load .mat files, compute dT array with mirror flip and compute
% the heat flux at the hot spots using the inverse method

% Process: 
% 1- Read the dataset spreadsheet
% 2- Load .mat files for each dataset
% 3- Calculate dT and apply image mirroring
% 3- Plot temperature fields and select hot spot
% 4- Apply inverse method to extract heat flux at each hot spot
% 5- Plot heat flux data for all datasets

% Note:
% In this code, we do not extract the data directly from the .seq files; we
% obtained the data from previosly extracted data stored in .mat files.

clc
clear all
close all

t0 = tic; % Time total process
disp('Start of task 1 ##################################################');

% Read dataset spreadsheet:
% =========================================================================
spreadsheetName = '\HeliconWindowIR_2020_02_XPs.xlsx';
home = cd;
cd ..
addr = pwd;
cd(home)
T = readtable([addr,spreadsheetName],'Sheet',1);
groups = unique(T.Group);

% Load .mat files:
% =========================================================================
groupsToAnalyze =  1:length(groups);
% groupsToAnalyze = 1:2;

for ii = groupsToAnalyze
    dataset{ii} = T.dataset(find(T.Group == ii,1));
    date{ii}    = T.date{find(T.Group == ii,1)};
    % dataset_1_IRdata
    pathName = ['C:\Users\nfc\Documents\ProtoMPEX_Ops\2020\2020_02\',date{ii},'\HeliconWindowIR\'];
    fileName = ['dataset_',num2str(dataset{ii}),'_IRdata.mat'];
    disp([pathName,fileName])
    d{ii} = load([pathName,fileName]);
end

%% 
% Select hot spot:
% =========================================================================
close all
%groupsToAnalyze = 1:length(groups);
for ii = groupsToAnalyze
    for jj = 1:numel(d{ii}.shot)
        for kk = 1:size(d{ii}.temperature{jj},3)
            deltaT{ii}{jj}(:,:,kk) = d{ii}.temperature{jj}(:,end:-1:1,kk) - d{ii}.temperature{jj}(:,end:-1:1,1);
        end
        t_deltaT{ii}{jj} = d{ii}.t_temperature{jj} - d{ii}.t_temperature{jj}(1);
        % Clear some memory:
        d{ii}.temperature{jj} = [];
    end
    
    [~,dum] = max(d{ii}.X);
    jj = dum;
    options.frames = 1:3:size(deltaT{ii}{jj},3)-30;
    options.frames = 55;
    options.colorbar = 1;
    options.magnitudePlotMode = 1;
    options.removeAxesTicks = 0;
    options.shot = d{ii}.shot(jj);
    options.zlim = [0,16];
    options.mirrorImage = 0;
    PlayMovieFromArray(deltaT{ii}{jj},options)

    % Select temperature probes
    switch ii
        case {1}
            xCenter{ii} = 205;
            yCenter{ii} = 087;
        case {2}
            xCenter{ii} = 205;
            yCenter{ii} = 087;
        case {3}
            xCenter{ii} = 246;
            yCenter{ii} = 162;
        case {4} % Middle
            xCenter{ii} = 249;
            yCenter{ii} = 162;
        case {5} % Bottom
            xCenter{ii} = 329;
            yCenter{ii} = 177;
        case {6} % Bottom
            xCenter{ii} = 276;
            yCenter{ii} = 113;
            xCenter{ii} = 310;
            yCenter{ii} = 045;
        case {7} % Bottom
            xCenter{ii} = 329;
            yCenter{ii} = 115;
        case {8} % Bottom
            xCenter{ii} = 291;
            yCenter{ii} = 175;
        case {9} % Bottom
            xCenter{ii} = 264;
            yCenter{ii} = 115;
        case {10} % Bottom
            xCenter{ii} = 206;
            yCenter{ii} = 107;
        case {11} % Bottom
            xCenter{ii} = 323;
            yCenter{ii} = 116;
        case {12} % Bottom
            xCenter{ii} = 302;
            yCenter{ii} = 058;
        case {13} % Bottom
            xCenter{ii} = 372;
            yCenter{ii} = 159;
    end
     % Plot location
    rngy_tempProbe = [(yCenter{ii}-10):(yCenter{ii}+10)];
    rngx_tempProbe = [(xCenter{ii}-10):(xCenter{ii}+10)];
    [x_tempProbe,y_tempProbe] = meshgrid(rngx_tempProbe,rngy_tempProbe);
    z_tempProbe = 2*deltaT{ii}{jj}(rngy_tempProbe,rngx_tempProbe,60);

    hold on
%     plot3(x_tempProbe,y_tempProbe,z_tempProbe','k.','MarkerSize',1)

    for jj = 1:numel(d{ii}.shot)
        for kk = 1:size(deltaT{ii}{jj},3)
            tempProbe{ii}{jj}(kk) = mean(mean(deltaT{ii}{jj}(rngy_tempProbe,rngx_tempProbe,kk)));
        end
        t_tempProbe{ii}{jj} = t_deltaT{ii}{jj};
    end   

end

%%
% ALN material properties
rho = 3300; 
kt   = 180 ;
cp  = 740 ;
a = kt/(rho*cp); % thermal diffusivity

% =========================================================================
% Geometry
Lz = 6.3500/1000; % ALN wall thickness in ProtoMPEX

% =========================================================================
% Dimensionless and normalized data
% Order of magnitude heat flux:
q0 = 1e6; 
% Characteristic time scale:
t_star = Lz*Lz/a;
% Characteristic temperature:
T_star = q0*Lz/kt;
% Normalized measured temperature and time:
    
for ii = groupsToAnalyze
    for jj = 1:numel(d{ii}.shot)
        z{ii}{jj}    = tempProbe{ii}{jj}'/T_star; 
        t_z{ii}{jj}  = t_tempProbe{ii}{jj}'/t_star;
        dt_z(jj) = t_z{ii}{jj}(2) - t_z{ii}{jj}(1);
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
    for jj = 1:numel(d{ii}.shot)
         for kk = 1:length(t_z{ii}{jj})
               K0{ii}{jj}(kk) = G_Impulse_1D(xx0,t_z{ii}{jj}(kk),Ns);
         end
         P{jj} = toeplitz(K0{ii}{jj},zeros(size(K0{ii}{jj})));
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
    if ii == 2
        disp('ii')
    end
    for jj = 1:numel(d{ii}.shot)
        qdummy{1} = zeros(size(t_z{ii}{jj}));
        [a1,a2,a3] = IHCP_ConjugateGradient(P{jj},qdummy,z{ii}{jj},dt_z(jj),1900,3,[]);
        q{jj} = a1{end};
        u{ii}{jj} = a2{end};
        J{jj} = a3;

        heatflux_tempProbe{ii}{jj} = q0*q{jj};
        t_heatflux_tempProbe{ii}{jj} = t_tempProbe{ii}{jj};
    end
    toc
    disp('Computation completed!') 
end
%%
close all

for ii = groupsToAnalyze
    lineColor = {'k','r','bl','g','m','c','k--','r--','bl--','g--','m--','c--','k:','r:','bl:','g:','m:','c:'};

    if 0
        disp('Testing inverse solution...')
        figure('Tag','section2','color','w')
        for jj = 1:numel(d{ii}.shot)
            subplot(3,3,jj)
            hold on
            hDummy1(1) = plot(t_z{ii}{jj}*t_star,z{ii}{jj}*T_star,'k.','MarkerSize',12);
            hDummy1(2) = plot(t_z{ii}{jj}*t_star,u{ii}{jj}*T_star,'r','LineWidth',1);
            xlim([0,1])
            ylim([0,15])
            box on
            grid on
            title(['Shot: ',num2str(d{ii}.shot(jj))])
        end
        hDummy2 = legend(hDummy1,'Exp.','Reconstruction');
        clear hDummy1 hDummy2 
    end 

    % Plot all data:
    disp('Plot all calculated heat fluxes')
    figure('Tag','section2','color','w')
    hold on
    for jj = 1:numel(d{ii}.shot)
        rng = find(t_tempProbe{ii}{jj} > 0.39 & t_tempProbe{ii}{jj} < 0.45);
        hDummy1(jj) = plot(t_tempProbe{ii}{jj},heatflux_tempProbe{ii}{jj}*1e-3,lineColor{jj},'LineWidth',2);
        
        heatflux_tempProbe_ss{ii} = heatflux_tempProbe{ii}{jj}(rng);
        heatflux_tempProbe_ss_mean{ii}(jj) = mean(heatflux_tempProbe_ss{ii});
        heatflux_tempProbe_ss_std{ii}(jj) = std(heatflux_tempProbe_ss{ii},1,1);

        plot(t_tempProbe{ii}{jj}(rng),heatflux_tempProbe_ss{ii}*1e-3,lineColor{jj},'LineWidth',4);
        legendText{jj} = [num2str(d{ii}.shot(jj)),', ',d{ii}.scanType{1},': ',num2str(d{ii}.X(jj))];
    end
    box on
    grid on
    xlim([0,2*max(d{ii}.pulseLength)*1e-3])
    ylim([0,800])
    title(['Dataset ',num2str(ii),', ',d{ii}.scanType{1},', ',d{ii}.limitMode{1},' limit, ',d{ii}.viewType{1},' view'])
    xlabel('time [s]','Interpreter','latex','FontSize',12)
    ylabel('[kWm$^{-2}$]','Interpreter','latex','FontSize',12)
    set(gca,'FontName','times','FontSize',11)
    hDummy2 = legend(hDummy1,legendText);
    set(gca,'PlotBoxAspectRatio',[2 1 1])
    clear hDummy1 hDummy2 

    saveFigure = 1;
    figureName = ['HeatFlux_Grou_',num2str(ii)];
    if saveFigure
        saveas(gcf,figureName,'tiffn')
    end

end

%%
close all
 figure('color','w'); 
% subplot(2,1,1)
hold on
kk = 1;
C = {'k','r','bl','g'};
legendText = {'E','G','H','F'};
for ii = groupsToAnalyze
    if strcmpi(d{ii}.limitMode,'window') & (strcmpi(d{ii}.viewType,'top') | strcmpi(d{ii}.viewType,'bottom'))
    [~,b] = sort(d{ii}.X);
    hwind(kk) = plot(d{ii}.X(b),heatflux_tempProbe_ss_mean{ii}(b)*1e-3,C{kk},'LineWidth',3);
%     legendText{kk} = [d{ii}.viewSide{1},d{ii}.viewType{1}];
    kk = kk + 1;
    grid on 
    box on
    end
end
set(gca,'FontSize',13,'FontName','Arial monospaced')
axis('square')
L = legend(hwind,legendText);
L.Location = 'northwest';
ylim([0,600])
xlim([0,200])
ylabel('[kWm$^{-2}$]','Interpreter','latex','FontSize',15)
xlabel('RF power [kW]','Interpreter','latex','FontSize',15)
title('Window limiter','Interpreter','latex','FontSize',15)
set(gcf,'Position',[85.0000  240.3333  405.3333  292.6667])

figure('color','w'); 
% subplot(2,1,2)
hold on
kk = 1;
legendText = {'D','A','B','C'};
for ii = groupsToAnalyze
    if strcmpi(d{ii}.limitMode,'mpex') & (strcmpi(d{ii}.viewType,'top') | strcmpi(d{ii}.viewType,'bottom'))
        if ii == 2
            continue
        end
    [~,b] = sort(d{ii}.X);
    hmpex(kk) = plot(d{ii}.X(b),heatflux_tempProbe_ss_mean{ii}(b)*1e-3,C{kk},'LineWidth',3);
%     legendText{kk} = [d{ii}.viewSide{1},d{ii}.viewType{1}];
    kk = kk + 1;
    grid on 
    box on
    end
end
set(gca,'FontSize',13,'FontName','Arial monospaced')
axis('square')
L = legend(hmpex,legendText);
L.Location = 'northwest';
ylim([0,600])
xlim([0,200])
ylabel('[kWm$^{-2}$]','Interpreter','latex','FontSize',15)
xlabel('RF power [kW]','Interpreter','latex','FontSize',15)
title('MPEX-like limiter','Interpreter','latex','FontSize',15)
set(gcf,'Position',[285.0000  240.3333  405.3333  292.6667])

return
t0 = toc(t0);
disp(['Total calculation time: ',num2str(t0),' sec'])





