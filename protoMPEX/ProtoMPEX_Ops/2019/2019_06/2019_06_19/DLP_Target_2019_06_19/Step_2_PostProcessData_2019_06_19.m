% need to complete ne2D plot
% pressure gauges

% Postprocess data:
clear all
close all

% load data:
% =========================================================================
load('Step_1_RawData_2019_06_19.mat')


% Process DLP data:
% =========================================================================

% Assemble input data structure:
% -------------------------------------------------------------------------
dlpRawData = RawData;
dlpRawData.shot(end)          = [];
dlpRawData.MetaData.shot(end) = [];
dlpRawData.V(end)             = [];
dlpRawData.I(end)             = [];

% Define DLP configuration file:
% -------------------------------------------------------------------------
AttType = 'Vx1,Ix1'; % Attenuation on the digitized signals
SweepType = 'iso'; % "iso" for isolated sweep, "niso" for non-isolated sweep

DLPType = '12MP';
Config.AMU = 2; % Ion mass in AMU
Config.tStart = 4.13; % [s]
Config.tEnd = 4.4;

Config.Center_V = 1; % Remove offset on V: 1 (yes) 0(no)
Config.Center_I = 1; % Remove offset on I: 1 (yes) 0(no)
Config.FilterDataInput = 1; % Filter input data with Savitsky-Golay filter order 3
Config.SGF = 7; % Frame of Savitsky-Golay filter 
Config.TimeMode = 2; % Effective time of sweep, (1) start of ramp or (2) mean time of ramp
Config.AreaType = 1; % 1: Cylindrical + cap

% The following script completes the configuration file
ProtoMPEX_DLPInfo_2018_06_04

% Use the DLP fitting routine:
% -------------------------------------------------------------------------
% INPUT: Config (structure), RawData (structure)
% OUTPUT: dlpData (structure)
dlpData = DLP_fit_V7(Config,dlpRawData);

% Extract data from the fitting process
% -------------------------------------------------------------------------
% Assign DLP data to variables
ni = dlpData.Ni;
te = dlpData.Te;
time = dlpData.time;
Ifit = dlpData.Ifit;
Ip = dlpData.Ip;
Vp = dlpData.Vp;
tm = dlpData.tm;
Vsweep = dlpData.Vsweep;
Isweep = dlpData.Isweep;
GlitchFlag = dlpData.GlitchFlag;
StdResNorm = dlpData.StdResNorm;

for s = 1:length(dlpRawData.shot)
    ni_m{s} = 0.5*(ni{s}{1} + ni{s}{2});
end

% Perform check and only select data with "good" fits
% -------------------------------------------------------------------------
NiMax = 1e21;
NiMin = 1e16; 
TeMax = 15; 
TeMin = 0;
NormalizedResidualForGoodFits = 6.4;

for s = 1:length(dlpRawData.shot)
    % For each shot find sweeps with "good" fits
    GoodFits{s} = GlitchFlag{s} == 0 & StdResNorm{s}<= NormalizedResidualForGoodFits ...
        & ni_m{s}>NiMin & ni_m{s}<NiMax & te{s}<=TeMax;
    % Extract "good" fits
    ne{s} = ni_m{s}(GoodFits{s});
    Te{s} = te{s}(GoodFits{s});
    t_ne{s} = time{s}(GoodFits{s});
end

% Plasma radius:
rprobe = RawData.MetaData.rgauge - 9.5;

%% time-dependent 2D array:
% =========================================================================
dii = 5;

[~,b] = sort(rprobe);
for ss = 1:numel(b)     
    for ii = 1:numel(t_ne{b(ss)})-dii         
        ne_m(ii,ss) = mean(ne{b(ss)}(ii:ii+dii));
        Te_m(ii,ss) = mean(Te{b(ss)}(ii:ii+dii));        
        t_ne_m(ii)  = mean(t_ne{b(ss)}(ii:ii+dii));      
    end
end
r_ne_m = rprobe(b);

% Radial location of Last uniterrupted magnetic flux surface at the
% location of DLP 12.5:
% -------------------------------------------------------------------------
rLUFS = 1.6; % [cm]
zLUFS = 4.11; % [m]

% Calculate steady-state values:
% -------------------------------------------------------------------------
r_ne_ss = r_ne_m;
rng = find(t_ne_m > 4.30 & t_ne_m < 4.33);
ne_ss = mean(ne_m(rng,:));
d_ne_ss = std(ne_m(rng,:),1);
Te_ss = mean(Te_m(rng,:));
d_Te_ss = std(Te_m(rng,:),1);

% Plot steady-state values:
figure('color','w')
hold on
hne = errorbar(r_ne_ss,ne_ss*1e-19,d_ne_ss*1e-19,'ksq-');
hTe = errorbar(r_ne_ss,Te_ss      ,d_Te_ss      ,'ro-' );
line(+[1,1]*rLUFS,[0,20],'color','g','LineStyle','--','LineWidth',2)
hLine = line(-[1,1]*rLUFS,[0,20],'color','g','LineStyle','--','LineWidth',2)
xlim([-4,4])
ylim([0,20])
box on
legentText{1} = '$n_e$ $\times$ $10^{19}$ [m$^{-3}$]';
legentText{2} = '$T_e$ [eV]';
legentText{3} = 'LUFS';
hL = legend([hne,hTe,hLine],legentText);
set(hL,'interpreter','Latex','FontSize',13)
set(gca,'FontName','Times','FontSize',12)
set(hne,'MarkerFaceColor','k','MarkerSize',9)
set(hTe,'MarkerFaceColor','r','MarkerSize',9)
xlabel('r [cm]','interpreter','Latex','FontSize',14)

% Save figure:
figureName = 'Step_2_SteadyState_NeTe';
saveas(gcf,figureName,'tiffn')

%% Interpolate data:
% =========================================================================
% Initial grid:
[rr,tt] = meshgrid(r_ne_m,linspace(t_ne_m(1),t_ne_m(end),numel(t_ne_m)));

% Final grid:
t_1D = linspace(4.13,4.4,100);
r_1D = linspace(r_ne_m(1),r_ne_m(end),50);
[r_2D,t_2D] = meshgrid(r_1D,t_1D);

[ne_2D] = interp2(rr,tt,ne_m,r_2D,t_2D,'cubic');
[Te_2D] = interp2(rr,tt,Te_m,r_2D,t_2D,'cubic');

%% Plot processed data:
% =========================================================================
if 0
    figure
    hold on
    [~,b] = sort(rprobe);
    for ss = b
        X = t_ne{ss};
        Z = ne{ss};
        Y = r_ne_m(ss)*ones(size(Z));
        plot3(X,Y,Z)
    end
    view([-130,30])

    figure
    hold on
    for ss = 1:numel(ne)
        X = t_ne{ss};
        Z = Te{ss};
        Y = r_ne_m(ss)*ones(size(Z));
        plot3(X,Y,Z)
    end
    view([-130,30])
end

% Pressure gauges:
% =========================================================================
t_PG3 = RawData.t_PG3{1}(1:end-1);
PG3 = [RawData.PG3{1}]*2;

t_PG4 = RawData.t_PG4{1}(1:end-1);
PG4 = [RawData.PG4{1}]*10;

figure; 
hold on
plot(t_PG3,PG3,'k')
plot(t_PG4,PG4,'r')

%% Plot time-dependent radial scan data:

% Preview data:
% =========================================================================
% RF pulse:
figure('color','w')
t_RF = RawData.t_RF{1}(1:end-1);
RF = RawData.RF{1};
plot(t_RF,RF)
xlim([4,5])
t_rfStart = 4.154;

% Plasma density
figure('color','w');
hold on
surf(r_2D,t_2D,ne_2D,'LineStyle','none')
hLine = line([-3,+3],[1,1]*t_rfStart,[1,1]*1e18,'color','r','LineStyle','-','LineWidth',2);
view([30,60])
box on
xlabel('r [cm]','Interpreter','latex','FontSize',13)
ylabel('t [s]','Interpreter','latex','FontSize',13)
zlabel('$n_e$ [m$^{-3}$]','Interpreter','latex','FontSize',13)
xlim([-2.5,2.5])
ylim([4.13,4.4])

hLeg = legend([hLine],'Start of RF');
hLeg.Location = 'north';
hLeg.Interpreter = 'latex';

colorbar

% Save figure:
figureName = 'Step_2_TimeDependent_NeTe';
saveas(gcf,figureName,'tiffn')