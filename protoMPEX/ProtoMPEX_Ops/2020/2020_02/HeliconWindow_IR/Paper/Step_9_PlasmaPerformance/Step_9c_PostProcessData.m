% need to complete ne2D plot
% pressure gauges

% Postprocess data:
clear all
close all

% load data:
% =========================================================================
load('Step_9b_RawData.mat')

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
Config.tEnd = 4.7;

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
rprobe = RawData.MetaData.rgauge - 10;

% Convert pressure gauge raw data into pressure in Pa:
% =========================================================================
% Convert miliTorr into Pascal:
fctr = 0.13333; 
% Label for each pressure gauge:
T = {'12.5','8.5','6.5','2.5'};

% Convert to Pascal:
for ss = 1:numel(RawData.shot(1:end-1))
    P{1}{ss} = (RawData.PG1{ss})*2*fctr;
    P{2}{ss} = (RawData.PG2{ss})*2*fctr;
    P{3}{ss} = (RawData.PG3{ss})*2*fctr;
    P{4}{ss} = (RawData.PG4{ss})*10*fctr;
    
    t_P{1}{ss} = (RawData.t_PG1{ss}(1:end-1));
    t_P{2}{ss} = (RawData.t_PG2{ss}(1:end-1));
    t_P{3}{ss} = (RawData.t_PG3{ss}(1:end-1));
    t_P{4}{ss} = (RawData.t_PG4{ss}(1:end-1));
end

% Remove offset from data:
% =========================================================================
for ss = 1:numel(RawData.shot(1:end-1))
    for pp = 1:4
        t0 = t_P{pp}{ss};
        rngMean = find(t0>3.6 & t0 <3.8);
        Offset{pp}{ss} = mean(P{pp}{ss}(rngMean));    
        P{pp}{ss} = P{pp}{ss} - Offset{pp}{ss};
    end
end

%% Plot data:
close all
% RF power trace:
% =========================================================================
% check when the RF starts and ends:
t_rfStart = 4.165;
t_rfEnd   = 4.645;
figure
hold on
t_RF = RawData.t_RF{1}(1:end-1);
RF   = RawData.RF{1};
plot(t_RF,RF)
line([1,1]*t_rfStart,[0,1.2],'LineWidth',2,'Color','g')
line([1,1]*t_rfEnd  ,[0,1.2],'LineWidth',2,'Color','g')
xlim([4,5]);

% Plasma density and electron temperature:
% =========================================================================
fontSize_axes = 12;
fontSize_label = 14;
fontSize_title = 13;
fontSize_legend = 11;

figure('color','w')
subplot(2,1,1)
hold on
hLs = line([1,1]*t_rfStart,[0,2]*1e20,'LineWidth',2,'Color','g')
hne(1) = plot(t_ne{1},ne{1},'k','LineWidth',2)
hne(2) = plot(t_ne{2},ne{2},'r','LineWidth',2)
hLe = line([1,1]*t_rfEnd  ,[0,2]*1e20,'LineWidth',2,'Color','g','LineStyle','--')
set(gca,'FontName','Times','FontSize',fontSize_axes)
ylim([0,12]*1e19)
xlim([4,4.9]);
box on
grid on
ylabel('$n_e$ $\times$ $10^{19}$ [m$^{-3}$]','Interpreter','latex','FontSize',fontSize_label);
xlabel('t [s]','Interpreter','latex','FontSize',fontSize_label);
text(4.02,10e19,'Target','Interpreter','latex','FontSize',fontSize_title,'EdgeColor','k')

subplot(2,1,2)
hold on
hLs = line([1,1]*t_rfStart,[0,6],'LineWidth',2,'Color','g')
rng1 = 10:193;
hte(1) = plot(t_ne{1}(rng1),Te{1}(rng1),'k','LineWidth',2)
rng2 = 10:193;
hte(2) = plot(t_ne{2}(rng2),Te{2}(rng2),'r','LineWidth',2)
hLe = line([1,1]*t_rfEnd  ,[0,6],'LineWidth',2,'Color','g','LineStyle','--')
set(gca,'FontName','Times','FontSize',fontSize_axes)
ylim([0,6])
xlim([4,4.9]);
box on
grid on
ylabel('$T_e$ [eV]','Interpreter','latex','FontSize',fontSize_label);
xlabel('t [s]','Interpreter','latex','FontSize',fontSize_label);
legendText{1} = ['Window-limit, shot: ',num2str(dlpData.shot(1))];
legendText{2} = ['MPEX-limit, shot: '  ,num2str(dlpData.shot(2))];
legendText{3} = ['RF start'];
legendText{4} = ['RF end'];
hLeg1 = legend([hte,hLs,hLe],legendText);
set(hLeg1,'Interpreter','Latex','FontSize',fontSize_legend,'Location','best','Box','on')
text(4.02,5.1,'Target','Interpreter','latex','FontSize',fontSize_title,'EdgeColor','k')

% Save figure:
figureName = 'Step_9c_NeTe';
saveas(gcf,figureName,'tiffn')

% Neutral gas pressure:
% =========================================================================
figure('color','w')
subplot(2,1,1)
hold on
hLs = line([1,1]*t_rfStart,[0,2]*1e20,'LineWidth',2,'Color','g')
hLe = line([1,1]*t_rfEnd  ,[0,2]*1e20,'LineWidth',2,'Color','g','LineStyle','--')
hpg4(1) = plot(t_P{4}{1},P{4}{1},'k','LineWidth',2)
hpg4(2) = plot(t_P{4}{2},P{4}{2},'r','LineWidth',2)
set(gca,'FontName','Times','FontSize',fontSize_axes)
ylim([0,2])
xlim([4,4.9]);
box on
grid on
ylabel('[Pa]','Interpreter','latex','FontSize',fontSize_label);
xlabel('t [s]','Interpreter','latex','FontSize',fontSize_label);
legendText{1} = ['Window-limit, shot: ',num2str(dlpData.shot(1))];
legendText{2} = ['MPEX-limit, shot: '  ,num2str(dlpData.shot(2))];
legendText{3} = ['RF start'];
legendText{4} = ['RF end'];
hLeg2 = legend([hpg4,hLs,hLe],legendText);
set(hLeg2,'Interpreter','Latex','FontSize',fontSize_legend,'Location','NorthEast')
text(4.02,1.7,'Plasma source','Interpreter','latex','FontSize',fontSize_title,'EdgeColor','k','BackgroundColor','w')

subplot(2,1,2)
hold on
hLs = line([1,1]*t_rfStart,[0,2]*1e20,'LineWidth',2,'Color','g')
hLe = line([1,1]*t_rfEnd  ,[0,2]*1e20,'LineWidth',2,'Color','g','LineStyle','--')
pp = 3;
hpg3(1) = plot(t_P{pp}{1},P{pp}{1},'k','LineWidth',2)
hpg3(2) = plot(t_P{pp}{2},P{pp}{2},'r','LineWidth',2)
set(gca,'FontName','Times','FontSize',fontSize_axes)
ylim([0,0.3])
xlim([4,4.9]);
box on
grid on
ylabel('[Pa]','Interpreter','latex','FontSize',fontSize_label);
xlabel('t [s]','Interpreter','latex','FontSize',fontSize_label);
text(4.02,0.26,'Electron heating region','Interpreter','latex','FontSize',fontSize_title,'EdgeColor','k','BackgroundColor','w')

% Save figure:
figureName = 'Step_9c_NeutralPressure';
saveas(gcf,figureName,'tiffn')

%% Plot type 2:


%% Functions:

% function [] = GetNeutralGasPressure(shot,shot_cal,gauge,units)
% Inputs:
% shot: double representing shotnumber
% end
