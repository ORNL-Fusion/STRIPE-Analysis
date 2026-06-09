% Task 2: Convert temperature data into heat flux

clc
clear all
close all

t0 = tic; % Time total process
disp('Start of task 2 ##################################################');

datasetToAnalyze = 1:8;

for ii = datasetToAnalyze
%     clearvars -except t0 ii datasetToAnalyze
    % Load data from .mat file:
    d{ii} = load(['dataset_',num2str(ii),'_IRdata.mat']);

    for jj = 1:numel(d{ii}.shot)
        for kk = 1:size(d{ii}.temperature{jj},3)
            deltaT{ii}{jj}(:,:,kk) = d{ii}.temperature{jj}(:,end:-1:1,kk) - d{ii}.temperature{jj}(:,end:-1:1,1);
        end
        t_deltaT{ii}{jj} = d{ii}.t_temperature{jj} - d{ii}.t_temperature{jj}(1);
    end
    
    jj = 1;
    options.frames = 1:2:size(deltaT{ii}{jj},3)-50;
    options.frames = 60;
    options.colorbar = 1;
    options.magnitudePlotMode = 1;
    options.removeAxesTicks = 0;
    options.shot = d{ii}.shot(jj);
    options.zlim = [0,15];
    options.mirrorImage = 0;
    PlayMovieFromArray(deltaT{ii}{jj},options)

    % Select temperature probes
    switch ii
        case {6} % MPEX Top
            yCenter{ii} = 159;
            xCenter{ii} = 374;
        case {1} % MPEX Top
            yCenter{ii} = 057;
            xCenter{ii} = 293;
            yCenter{ii} = 178;
            xCenter{ii} = 346;
        case {2} % Middle
            yCenter{ii} = 042;
            xCenter{ii} = 314;
        case {3,4,7,8} % Bottom
            yCenter{ii} = 122;
            xCenter{ii} = 323;
        case {5} % Bottom
            yCenter{ii} = 056;
            xCenter{ii} = 299;
       case {7} % Bottom
            yCenter{ii} = 096;
            xCenter{ii} = 328;
       case {8} % Bottom
            yCenter{ii} = 086;
            xCenter{ii} = 341;
    end
     % Plot location
    rngy_tempProbe = [(yCenter{ii}-10):(yCenter{ii}+10)];
    rngx_tempProbe = [(xCenter{ii}-10):(xCenter{ii}+10)];
    [x_tempProbe,y_tempProbe] = meshgrid(rngx_tempProbe,rngy_tempProbe);
    z_tempProbe = 2*deltaT{ii}{jj}(rngy_tempProbe,rngx_tempProbe,60);

    hold on
    plot3(x_tempProbe,y_tempProbe,z_tempProbe','k.','MarkerSize',1)

    for jj = 1:numel(d{ii}.shot)
        for kk = 1:size(d{ii}.temperature{jj},3)
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
    
for ii = datasetToAnalyze
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
for ii = datasetToAnalyze
    lineColor = {'k','r','bl','g','m','c','k--','r--','bl--','g--','m--','c--','k:','r:','bl:','g:','m:','c:'};

    disp('Testing inverse solution...')
    figure('Tag','section2','color','w')
    for jj = 1:numel(d{ii}.shot)
        subplot(3,3,jj)
        hold on
        hDummy1(1) = plot(t_z{ii}{jj}*t_star,z{ii}{jj}*T_star,'k.','MarkerSize',12);
        hDummy1(2) = plot(t_z{ii}{jj}*t_star,u{ii}{jj}*T_star,'r','LineWidth',1);
        xlim([0,1])
        ylim([0,25])
        box on
        grid on
        title(['Shot: ',num2str(d{ii}.shot(jj))])
    end
    hDummy2 = legend(hDummy1,'Exp.','Reconstruction');
    clear hDummy1 hDummy2 


    % Plot all data:
    disp('Plot all calculated heat fluxes')
    figure('Tag','section2','color','w')
    hold on
    for jj = 1:numel(d{ii}.shot)
        hDummy1(jj) = plot(t_z{ii}{jj}*t_star,heatflux_tempProbe{ii}{jj}*1e-3,lineColor{jj},'LineWidth',2);
        legendText{jj} = [num2str(d{ii}.shot(jj)),', ',d{ii}.scanType{1},': ',num2str(d{ii}.X(jj))];
    end
    box on
    grid on
    xlim([0,2*max(d{ii}.pulseLength)*1e-3])
    ylim([0,1200])
    title(['Dataset ',num2str(ii),', ',d{ii}.scanType{1},', ',d{ii}.limitMode{1},' limit, ',d{ii}.viewType{1},' view'])
    xlabel('time [s]','Interpreter','latex','FontSize',12)
    ylabel('[kWm$^{-2}$]','Interpreter','latex','FontSize',12)
    set(gca,'FontName','times','FontSize',11)
    hDummy2 = legend(hDummy1,legendText);
    set(gca,'PlotBoxAspectRatio',[2 1 1])
    clear hDummy1 hDummy2 

    saveFigure = 1;
    figureName = ['HeatFlux_Dataset_',num2str(ii)];
    if saveFigure
        saveas(gcf,figureName,'tiffn')
    end

end

t0 = toc(t0);
disp(['Total calculation time: ',num2str(t0),' sec'])

