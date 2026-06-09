% Objective:
% Extract heat flux at various points (>10) over the 6 views associated
% with a single RF power into order to estimate the total power onto the
% window and the uncertainty:

%% Process: 
% 1- Read the dataset spreadsheet
% 2- Assemble the shot series
% 3- Gather data into structures
% 4- Calculate temperature difference and mirror image
% 5- Create grid points in all views
% 6- Extract time evolution of temperature grid points
% 7- Normalize grid point temperature data
% 8- Create convolution operator
% 9- Calculate heat flux using inverse method
% 10- Extract steady state heat flux from grid points
% 11- Calculate total power to the window
% 12- Plot heat flux data
% 13- Test and inspect solutions

% Notes:
% In this code, we do not extract the data directly from the .seq files; we
% obtained the data from previosly extracted data stored in .mat files.
% The effect of the mirror is applied in step 4

clc
clear all
close all

%% 1- Read the dataset spreadsheet
% =========================================================================
spreadsheetName = '/HeliconWindowIR_2020_02_XPs.xlsx';
home = cd;
cd ..
addr = pwd;
cd(home)
T = readtable([addr,spreadsheetName],'Sheet',1);
groups = unique(T.Group);

%% 2- Assemble the shot series:
% =========================================================================

% MPEX-Limit, 140 kW:
shotSeries{1} = [29077,29082,29049,29117,29120,29128];
rfPower{1}    = [142  ,142  ,142  ,135  ,150  ,141  ];

% MPEX-Limit, 120 kW:
shotSeries{2} = [29115,29122,29127,29076,29081,29047];
rfPower{2}    = [120  ,116  ,118  ,126  ,126  ,119  ];

% MPEX-Limit, 100 kW:
shotSeries{3} = [29075,29080,29066,29114,29123,29126];
rfPower{3}    = [107  ,102  ,105  ,100  ,97   ,96   ];

% MPEX-Limit, 80 kW:
shotSeries{4} = [29070,29079,29067,29113,29124,29125];
rfPower{4}    = [83   ,88   ,83   ,78   ,74   ,72   ];

% Window-Limit, 130 kW:
shotSeries{5} = [29101,29100,29106,29145,29136,29132];
rfPower{5}    = [130  ,131  ,132  ,132  ,133  ,133  ];

% Window-Limit, 115 kW:
shotSeries{6} = [29102,29099,29105,29144,29137,29131];
rfPower{6}    = [117  ,116  ,116  ,119  ,118  ,119  ];

% Window-Limit, 100 kW:
shotSeries{7} = [29103,29098,29108,29143,29138,29130];
rfPower{7}    = [103  ,104  ,103  ,103  ,103  ,103  ];

% Window-Limit, 80 kW:
shotSeries{8} = [29104,29097,29107,29141,29139,29129];
rfPower{8}    = [88   ,86   ,88   ,88   ,88   ,87  ];

% Select shots corresponding to 6 views at the same RF power:
% =========================================================================
viewType = 'bottom';
limitMode = 'window';
viewSide = 'non-pit';
rfPwr = [70,90];
n = find(strcmpi(T.limitMode,limitMode) & strcmpi(T.viewType,viewType) & strcmpi(T.viewSide,viewSide) & T.X>rfPwr(1) & T.X<rfPwr(2))
disp( {'ii: ', num2str(T.Group(n(end))), 'jj: ', num2str((n(end))) } )
T.rfPwrNet(n(end))
T.shot(n(end))

%% 3- Gather data into structures:
% =========================================================================

cmptType = 2;
switch cmptType
    case 1
        seriesToAnalyze = 1:numel(shotSeries);
    case 2
        seriesToAnalyze = 1;
end

for kk = seriesToAnalyze
    rootAddress = 'C:\Users\nfc\Documents\ProtoMPEX_Ops\2020\';
    for jj = 1:numel(shotSeries{kk})
        % shot index:
        n = find(T.shot == shotSeries{kk}(jj));

        % Define address and file name:
        relPath = [T.date{n}]; relPath = [relPath(1:7),'\',relPath];
        fileName = ['\dataset_',num2str(T.dataset(n)),'_IRdata.mat'];
        filePath = [rootAddress,relPath,'\HeliconWindowIR',fileName];

        % load structure of .mat file:
        d = matfile(filePath);

        % Indentify shot needed within .mat file:
        m = find(d.shot == T.shot(n));

        % Extract required data from .mat file:
        f{kk}{jj}.shot = d.shot(1,m);
        f{kk}{jj}.limitMode = d.limitMode(1,m);
        f{kk}{jj}.viewType  = d.viewType(1,m);
        f{kk}{jj}.viewSide  = d.viewSide(1,m);
        f{kk}{jj}.temperature   = cell2mat(d.temperature(1,m));
        f{kk}{jj}.t_temperature = cell2mat(d.t_temperature(1,m));
        f{kk}{jj}.rfPwr = d.X(1,m);
        f{kk}{jj}.t0Plasma = d.t0Plasma(1,m);
        f{kk}{jj}.intf = cell2mat(d.intf(1,m));
        f{kk}{jj}.rngPlasma = cell2mat(d.rngPlasma(1,m));
    end
end

%% 4- Calculate temperature difference and mirror image:
% =========================================================================
for kk = seriesToAnalyze;
    for jj = 1:numel(shotSeries{kk})
        % Frame prior to plasma:
        n0 = 1;

        % Calculate dT and mirror image:
        for rr = 1:size(f{kk}{jj}.temperature,3)
            f{kk}{jj}.dT(:,:,rr)   = f{kk}{jj}.temperature(:,end:-1:1,rr) - f{kk}{jj}.temperature(:,end:-1:1,n0);
        end
        
        % Calculate the time base:
        f{kk}{jj}.t_dT = f{kk}{jj}.t_temperature - f{kk}{jj}.t_temperature(1);

        % Clear memory:
        f{kk}{jj}.temperature = [];
    end
end

%% 5- Create grid points in all views:
% =========================================================================
options.colorbar = 1;
options.magnitudePlotMode = 2;
options.removeAxesTicks = 0;
options.zlim = [0,15];
options.mirrorImage = 0;

cmptflg{1} = 0*[1,1,1,1,1,1];
cmptflg{2} = 0*[0,1,1,1,1,1];
cmptflg{3} = 0*[1,1,1,1,1,1];
cmptflg{4} = 0*[1,1,1,1,1,1];

cmptflg{5} = 0*[0,0,1,0,0,0];
cmptflg{6} = 0*[0,0,1,0,0,0];
cmptflg{7} = 0*[0,0,1,0,0,0];
cmptflg{8} = 0*[0,0,1,0,0,0];

for kk = seriesToAnalyze
    for jj = 1:numel(shotSeries{kk})
        
        % Identify which view side and type:
        % -----------------------------------------------------------------
        if strcmpi(f{kk}{jj}.limitMode,'MPEX')
        switch cell2mat(f{kk}{jj}.viewType)
            case {'Top'}
                viewTypeCase = [1,4];
            case {'Middle'}
                viewTypeCase = [2,5];
            case {'Bottom'}
                viewTypeCase = [3,6];
        end
        switch cell2mat(f{kk}{jj}.viewSide)
            case {'Pit'}
                viewSideCase = [1,2,3];
            case {'non-Pit'}
                viewSideCase = [4,5,6];
        end
        fileid = 'mpexLim';        
        elseif strcmpi(f{kk}{jj}.limitMode,'Window')
            switch cell2mat(f{kk}{jj}.viewType)
            case {'Top'}
                viewTypeCase = [1,4];
            case {'Middle'}
                viewTypeCase = [2,5];
            case {'Bottom'}
                viewTypeCase = [3,6];
        end
        switch cell2mat(f{kk}{jj}.viewSide)
            case {'Pit'}
                viewSideCase = [1,2,3];
            case {'non-Pit'}
                viewSideCase = [4,5,6];
        end
        fileid = 'windLim'; 
        end
         % -----------------------------------------------------------------
             
        % Determine file name based on viewSide and viewType:
        n1 = intersect(viewSideCase,viewTypeCase);
        fileName = [fileid,'_grid_shotSeries_',num2str(1),'_shot_',num2str(n1)];

        % find the hottest frame:
        dum1 = f{kk}{jj}.rngPlasma;
        [~,n0] = max(f{kk}{jj}.intf(dum1));
        f{kk}{jj}.framePlasmaEnd = n0;

        % Plot hottest frame:
        options.frames = n0;
        options.shot = f{kk}{jj}.shot;
        PlayMovieFromArray(f{kk}{jj}.dT,options)
        hold on

        % Select previously defined or produce new grid points:
        if cmptflg{kk}(jj)
            col = 0; row = 0;
            for rr = 1:40
                [c,r] = ginput(1);
                col(rr) = round(c);
                row(rr) = round(r);
                plot3(col(rr),row(rr),f{kk}{jj}.dT(row(rr),col(rr),n0),'ko','MarkerFaceColor','w','MarkerSize',4)
            end         
            save(fileName,'col','row')
        else
            load(fileName)
        end
        f{kk}{jj}.row = row;
        f{kk}{jj}.col = col;
        
        % Plot grid points on dT frame:
        h = plot3(f{kk}{jj}.col,f{kk}{jj}.row,f{kk}{jj}.dT(f{kk}{jj}.row,f{kk}{jj}.col,n0),'ko',...
            'MarkerFaceColor','w','MarkerSize',4);
        pause(1)
    end
    close all
end

%% 6- Extract time evolution of temperature grid points:
% =========================================================================
tic
for kk = seriesToAnalyze;
    for jj = 1:numel(shotSeries{kk})
        disp(['jj: ',num2str(jj)])
        for ss = 1:numel(f{kk}{jj}.row)
            % Define the location of the grid point:
            rw = f{kk}{jj}.row(ss);
            cl = f{kk}{jj}.col(ss);

            % Create range of points for averaging:
            rw_rng = [rw-10:rw+10];
            cl_rng = [cl-10:cl+10];

            % Avoid indices less than 1:
            rw_rng(find(rw_rng<1)) = 1;
                
            % Extract the mean value of dT at grid point:
            for rr = 1:size(f{kk}{jj}.dT,3)
                f{kk}{jj}.dT_grid(rr,ss) = mean(mean(f{kk}{jj}.dT(rw_rng,cl_rng,rr)));
            end
        end
    end
end
toc
%% 7- Normalize grid point temperature data:
% =========================================================================
% ALN material properties:
rho = 3300; 
kt   = 180 ;
cp  = 740 ;
a = kt/(rho*cp); % thermal diffusivity

% Geometry:
Lz = 6.3500/1000; % ALN wall thickness in ProtoMPEX

% Dimensionless and normalized data:
% Order of magnitude heat flux:
q0 = 1e6; 
% Characteristic time scale:
t_star = Lz*Lz/a;
% Characteristic temperature:
T_star = q0*Lz/kt;

% Normalized measured temperature and time:
for kk = seriesToAnalyze;
    for jj = 1:numel(shotSeries{kk})  
        z{kk}{jj}    = f{kk}{jj}.dT_grid/T_star;
        t_z{kk}{jj}  = f{kk}{jj}.t_dT'/t_star;
        dt_z{kk}{jj} = t_z{kk}{jj}(2) - t_z{kk}{jj}(1);
    end
end

%% 8- Create convolution operator:
% =========================================================================
xx0 = 0.0005;
% To capture the details of the front surface we need to use 1e4
% partial sum terms and use xx0 = 0.01 to 0.05

% Partial sum
% "n" is the index
disp('Computing Toeplitz matrix...')
tic
Ns = 1e4;

% Create convolution operator P:
for kk = seriesToAnalyze;
    for jj = 1:numel(shotSeries{kk})
             for rr = 1:length(t_z{kk}{jj})
                   G0{kk}{jj}(rr) = G_Impulse_1D(xx0,t_z{kk}{jj}(rr),Ns);
             end
             P{kk}{jj} = toeplitz(G0{kk}{jj},zeros(size(G0{kk}{jj})));
    end
end
toc
disp('Computation complete')    

%% 9- Calculate heat flux using inverse method:
% =========================================================================
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

% Use conjugate gradient method:
for kk = seriesToAnalyze
    for jj = 1:numel(shotSeries{kk})
        for ss = 1:numel(f{kk}{jj}.row)
            qdummy{1} = zeros(size(t_z{kk}{jj}));
            [a1,a2,a3] = IHCP_ConjugateGradient(P{kk}{jj},qdummy,z{kk}{jj}(:,ss),dt_z{kk}{jj},700,3,[]);
            q{kk}{jj}(:,ss) = a1{end};
            u{kk}{jj}(:,ss) = a2{end};
            J{kk}{jj}(:,ss) = a3;
        end
        heatflux{kk}{jj} = q0*q{kk}{jj};
        t_heatflux{kk}{jj} = f{kk}{jj}.t_dT;
    end
end
toc
disp('Computation completed!') 

%% 10- Extract steady state heat flux from grid points:
% =========================================================================
% Select range of points to time averaging:
rng = 30:45;

% Apply time-averaging operation on each grid point:
for kk = seriesToAnalyze
    for jj = 1:numel(shotSeries{kk})
        for ss = 1:numel(f{kk}{jj}.row)
            q_ss{kk}(jj,ss) = mean(heatflux{kk}{jj}(rng,ss));
        end
    end
end

% Apply averaging over all grid points for each shot:
for kk = seriesToAnalyze
    for jj = 1:numel(shotSeries{kk})
        heatflux_ss{kk}(jj) = mean(q_ss{kk}(jj,:));
        dheatflux_ss{kk}(jj) = std(q_ss{kk}(jj,:),1);
    end
end

%% 11- Calculate total power to the window:

% Dimensions of helicon window:
L = 30/100;
diam = 12.4/100;
A = pi*diam*L;

% Calculate total power on window:
for kk = seriesToAnalyze
    for jj = 1:numel(shotSeries{kk})
        pwr_ss{kk}(jj) = heatflux_ss{kk}(jj)*A;
        dpwr_ss{kk}(jj) = dheatflux_ss{kk}(jj)*A;
    end
    pwr_mean(kk) = mean(pwr_ss{kk});
    dpwr_mean(kk) = std(pwr_ss{kk},1);
    rfPwr_mean(kk) = mean(rfPower{kk});
    drfPwr_mean(kk) = std(rfPower{kk},1);
end

% Apply a straight line fit to the data:
% MPEX-like limiter case:
P1  = polyfit(rfPwr_mean(1:4),pwr_mean(1:4),1);
P1u = polyfit(rfPwr_mean(1:4),pwr_mean(1:4) + 0.66*dpwr_mean(1:4),1);
P1l = polyfit(rfPwr_mean(1:4),pwr_mean(1:4) - 0.66*dpwr_mean(1:4),1);
% Window limiter case:
P2 = polyfit(rfPwr_mean(5:8),pwr_mean(5:8),1);
P2u = polyfit(rfPwr_mean(5:8),pwr_mean(5:8) + 0.66*dpwr_mean(5:8),1);
P2l = polyfit(rfPwr_mean(5:8),pwr_mean(5:8) - 0.66*dpwr_mean(5:8),1);

%% 12- Plot heat flux data:

figure;
% Plot steady-state time-averaged surface heat flux obtained for each shot:
subplot(1,2,1)
hold on
for kk = 1:4
    for jj = 1:numel(shotSeries{kk})
        errorbar(f{kk}{jj}.rfPwr,heatflux_ss{kk}(jj)*1e-3,dheatflux_ss{kk}(jj)*1e-3,'k.')
    end
end
for kk = 5:8
    for jj = 1:numel(shotSeries{kk})
        errorbar(f{kk}{jj}.rfPwr,heatflux_ss{kk}(jj)*1e-3,dheatflux_ss{kk}(jj)*1e-3,'r.')
    end
end
ylim([0,500])
xlim([0,200])
box on
axis('square')

% Total power:
subplot(1,2,2)
hold on
for kk = 1:4
    for jj = 1:numel(shotSeries{kk})
        errorbar(f{kk}{jj}.rfPwr,pwr_ss{kk}(jj)*1e-3,dpwr_ss{kk}(jj)*1e-3,'ko')
    end
end
for kk = 5:8
    for jj = 1:numel(shotSeries{kk})
        errorbar(f{kk}{jj}.rfPwr,pwr_ss{kk}(jj)*1e-3,dpwr_ss{kk}(jj)*1e-3,'ro')
    end
end
ylim([0,50])
xlim([0,200])
box on
axis('square')

% Data from UCSD in hydrogen with window limit:
rf_UCSD  = [2   , 4  , 6  , 8  , 10 , 12 ,  14,  16,  18,  20];
pwr_UCSD = [0.58,1.33,1.87,2.34,2.81,3.38,3.91,4.41,4.93,5.24];

figureName = 'IntegratedWindowPower';
figure;
hold on
h1 = errorbar(rfPwr_mean(1:4),pwr_mean(1:4)*1e-3,dpwr_mean(1:4)*1e-3,'ko','MarkerFaceColor','k','MarkerSize',7);
h2 = errorbar(rfPwr_mean(5:8),pwr_mean(5:8)*1e-3,dpwr_mean(5:8)*1e-3,'ro','MarkerFaceColor','r','MarkerSize',7);
h3 = plot(rf_UCSD,pwr_UCSD,'blsq','MarkerFaceColor','bl','MarkerSize',7);
legendText{1} = ['MPEX-limiter, Proto-MPEX'];
legendText{2} = ['Window-limiter, Proto-MPEX'];
legendText{3} = ['Window-limiter, CSDX'];

x = linspace(60,200);
plot(x,polyval(P1,x)*1e-3,'k--','LineWidth',1)
plot(x,polyval(P1u,x)*1e-3,'k:')
plot(x,polyval(P1l,x)*1e-3,'k:')

plot(x,polyval(P2,x)*1e-3,'r--','LineWidth',1)
plot(x,polyval(P2u,x)*1e-3,'r:')
plot(x,polyval(P2l,x)*1e-3,'r:')

ylim([0,60])
xlim([0,200])
box on
axis('square')
set(gcf,'color','w')
grid on
set(gca,'FontName','console')
ylabel('[kW]','Interpreter','latex','FontSize',13)
xlabel('P$_{RF}$[kW]','Interpreter','latex','FontSize',13)
hl = legend([h1,h2,h3],legendText);
set(hl,'Interpreter','Latex','FontSize',11) 

saveas(gcf,figureName,'tiffn')

return

%% 13- Test and inspect solutions:
% =========================================================================
if 1
    close all
    kk = 2;
    jj = 1;
    f{kk}{jj}.shot
    figure
    for ss = 1:40
        subplot(5,8,ss)
        hold on
        plot(t_z{kk}{jj},z{kk}{jj}(:,ss)*T_star,'k','LineWidth',3)
        plot(t_z{kk}{jj},u{kk}{jj}(:,ss)*T_star,'r.')
        ylim([0,20])
        xlim([0,2])
        titleText = ['r: ',num2str(f{kk}{jj}.row(ss)),', col: ',num2str(f{kk}{jj}.col(ss))];
        title(titleText)
    end

    figure
    for ss = 1:40
        subplot(5,8,ss)
        hold on
        plot(t_heatflux{kk}{jj},heatflux{kk}{jj}(:,ss)*1e-3,'k','LineWidth',3)
        plot(t_heatflux{kk}{jj}(rng),heatflux{kk}{jj}(rng,ss)*1e-3,'r','LineWidth',3)

        ylim([0,800])
        xlim([0,0.6])
        titleText = ['r: ',num2str(f{kk}{jj}.row(ss)),', col: ',num2str(f{kk}{jj}.col(ss))];
        title(titleText)
    end
end
