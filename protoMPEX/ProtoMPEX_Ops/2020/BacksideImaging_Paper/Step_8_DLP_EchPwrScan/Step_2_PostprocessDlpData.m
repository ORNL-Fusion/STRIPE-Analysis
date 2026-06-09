% Step 2: Postprocess DLP data

clear all
close all

% Flags
% =========================================================================
saveFig = 0;

% load data:
% =========================================================================
load('Step_1_GetDlpRawData_EchPwrScan_20200330.mat')

% Process DLP data:
% =========================================================================
% Assemble input data structure:
% -------------------------------------------------------------------------
dlpRawData = RawData;

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
Config.SGF = 17; % Frame of Savitsky-Golay filter 
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

% Extract only good fits:
% =========================================================================
NiMax = 1e20;
NiMin = 1e16; 
TeMax = 12; 
TeMin = 0;
NormalizedResidualForGoodFits = 0.3;

for s = 1:length(dlpRawData.shot)
    % For each shot find sweeps with "good" fits
    GoodFits{s} = GlitchFlag{s} == 0 & StdResNorm{s}<= NormalizedResidualForGoodFits ...
        & ni_m{s}>NiMin & ni_m{s}<NiMax & te{s}<=TeMax;
    % Extract "good" fits
    ne{s} = ni_m{s}(GoodFits{s});
    Te{s} = te{s}(GoodFits{s});
    t_ne{s} = time{s}(GoodFits{s});
end

% Assemble radial coordinate:
% =========================================================================
% Plasma pwr28GHz:
pwr28GHz = RawData.MetaData.pwr28GHz';
[~,b] = sort(pwr28GHz);
plasma.pwr28GHz = pwr28GHz(b);

% Interpolate data to a common time base:
% =========================================================================
plasma.t = linspace(4.15,4.65,100)';
for s = 1:numel(b)
    k = b(s);
    plasma.ne(s,:) = interp1(t_ne{k},ne{k},plasma.t);
    plasma.Te(s,:) = interp1(t_ne{k},Te{k},plasma.t);
end

% Plot surface:
% =========================================================================
figure('color','w');
surf(plasma.t,plasma.pwr28GHz,plasma.ne,'LineStyle','none')
view([0,90])

figure('color','w');
surf(plasma.t,plasma.pwr28GHz,plasma.Te,'LineStyle','none')
view([0,90])

figure('color','w');
hold on
for ii = 1:numel(b)
    plot(plasma.t,movmean(plasma.Te(ii,:),3))
end

%%
% Extract helicon and ECH profiles:
% =========================================================================
t1 = 4.45;
t2 = 4.50;

t1 = 4.6;
t2 = 4.63;

t3 = 4.55;
t4 = 4.57;

rng_13MHz = find(plasma.t >= t1 & plasma.t <= t2);
rng_28GHz = find(plasma.t >= t3 & plasma.t <= t4);

ne_13MHz_mean = mean(plasma.ne(:,rng_13MHz),2);
ne_13MHz_std  = std(plasma.ne(:,rng_13MHz),1,2);

ne_28GHz_mean = mean(plasma.ne(:,rng_28GHz),2);
ne_28GHz_std  = std(plasma.ne(:,rng_28GHz),1,2);

Te_13MHz_mean = mean(plasma.Te(:,rng_13MHz),2);
Te_13MHz_std  = std(plasma.Te(:,rng_13MHz),1,2);

Te_28GHz_mean = mean(plasma.Te(:,rng_28GHz),2);
Te_28GHz_std  = std(plasma.Te(:,rng_28GHz),1,2);

Pe_13MHz_mean = e_c*ne_13MHz_mean.*Te_13MHz_mean;
Pe_13MHz_std  = Pe_13MHz_mean.*sqrt( (Te_13MHz_std./Te_13MHz_mean).^2 + (ne_13MHz_std./ne_13MHz_mean).^2 );

Pe_28GHz_mean = e_c*ne_28GHz_mean.*Te_28GHz_mean;
Pe_28GHz_std  = Pe_28GHz_mean.*sqrt( (Te_28GHz_std./Te_28GHz_mean).^2 + (ne_28GHz_std./ne_28GHz_mean).^2 );


%% Plot data: RF traces
% =========================================================================
close all

% 13.56 MHz trace:
RF_13MHz   = RawData.RF_13MHz{5};
RF_13MHz   = 120*RF_13MHz/max(RF_13MHz);
t_RF_13MHz = RawData.t_RF_13MHz{5}(1:end-1);
rng_13MHz  = find(t_RF_13MHz >= t1 & t_RF_13MHz <= t2);

% 28 GHz trace:
RF_28GHz   = movmean(RawData.RF_28GHz{5},10);
RF_28GHz   = 55*RF_28GHz/max(RF_28GHz);
t_RF_28GHz = RawData.t_RF_28GHz{5}(1:end-1);
rng_28GHz  = find(t_RF_28GHz >= t3 & t_RF_28GHz <= t4);

figure('color','w');
hold on
hRF(1) = plot(t_RF_13MHz,RF_13MHz,'k','lineWidth',2);
hRF(2) = plot(t_RF_28GHz,RF_28GHz,'r','lineWidth',2);

% Averaging regions:
plot(t_RF_13MHz(rng_13MHz),RF_13MHz(rng_13MHz),'g','lineWidth',4);
plot(t_RF_28GHz(rng_28GHz),RF_28GHz(rng_28GHz),'g','lineWidth',4);

% Formatting:
set(gca,'FontName','Times','fontSize',11)
box on
xlim([4.1,4.75])
ylim([0,150])

% Labels:
xlabel('time [s]','interpreter','latex','fontSize',13)
ylabel('Power [kW]','interpreter','latex','fontSize',13)
hLeg = legend(hRF,'13.56 MHz','28 GHz');
set(hLeg,'interpreter','latex','FontSize',13)

% Save figure:
% =========================================================================
if saveFig
    figureName = ['Step_2_RfPwrTraces'];

    % PDF figure:
    exportgraphics(gcf,[figureName,'.pdf'],'Resolution',300) 

    % TIFF figure:
    exportgraphics(gcf,[figureName,'.tiff'],'Resolution',600) 
end

%% Plot data:
% =========================================================================

figure('color','w')
set(gcf,'position',[266   154   654   464])

fontSize.label = 10;
fontSize.axes = 10;
fontSize.abc  = 15;
fontSize.legend = 8;

subplot(1,3,1)
hold on
hne(1) = errorbar(plasma.pwr28GHz,ne_13MHz_mean,ne_13MHz_std);
hne(2) = errorbar(plasma.pwr28GHz,ne_28GHz_mean,ne_28GHz_std);
% Formatting:
set(hne(1),'color','k','LineWidth',2)
set(hne(2),'color','r','LineWidth',2)
ylim([0,5E19])
xlim([0,75])
axis square
box on
set(gca,'FontName','Times','fontSize',11)
% Labels:
xlabel('28 GHz [kW]','interpreter','latex','fontSize',fontSize.label)
ylabel('$n_e$ [m$^{-3}$]','interpreter','latex','fontSize',fontSize.label)
text(2.5,4.2E19,'(a)','fontSize',fontSize.abc ,'interpreter','latex')

subplot(1,3,2)
hold on
hTe(1) = errorbar(plasma.pwr28GHz,Te_13MHz_mean,Te_13MHz_std);
hTe(2) = errorbar(plasma.pwr28GHz,Te_28GHz_mean,Te_28GHz_std);
% Formatting:
set(hTe(1),'color','k','LineWidth',2)
set(hTe(2),'color','r','LineWidth',2)
ylim([0,15])
xlim([0,75])
axis square
box on
set(gca,'FontName','Times','fontSize',11)
% Labels:
xlabel('28 GHz [kW]','interpreter','latex','fontSize',fontSize.label)
ylabel('$T_e$ [eV]','interpreter','latex','fontSize',fontSize.label)
text(2.5,12.5,'(b)','fontSize',fontSize.abc ,'interpreter','latex')

subplot(1,3,3)
hold on
hPe(1) = errorbar(plasma.pwr28GHz,Pe_13MHz_mean,Pe_13MHz_std);
hPe(2) = errorbar(plasma.pwr28GHz,Pe_28GHz_mean,Pe_28GHz_std);
% Formatting:
set(hPe(1),'color','k','LineWidth',2)
set(hPe(2),'color','r','LineWidth',2)
ylim([0,60])
xlim([0,75])
axis square
box on
set(gca,'FontName','Times','fontSize',fontSize.axes)
% Labels:
xlabel('28 GHz [kW]','interpreter','latex','fontSize',fontSize.label)
ylabel('$P_e$ [Pa]','interpreter','latex','fontSize',fontSize.label)
hLeg = legend(hPe,'prior to 28GHz (t = 4.45 s)','during 28 GHz (t = 4.55 s)');
set(hLeg,'interpreter','latex','FontSize',fontSize.legend,'Location','north')
hLeg.Position = [0.7077    0.5848    0.1870    0.0717];
text(2,33,'(c)','fontSize',fontSize.abc ,'interpreter','latex')

% htile.TileSpacing = 'compact';

% Save figure:
% =========================================================================
if saveFig
    figureName = ['Step_2_Dlp_28GHzPowerScan'];

    % PDF figure:
    exportgraphics(gcf,[figureName,'.pdf'],'Resolution',300) 

    % TIFF figure:
    exportgraphics(gcf,[figureName,'.tiff'],'Resolution',600) 
end

%% Heat flux based in DLP:

Rie = [0,0.5,0.7];
Eiz = 13.6;
Ediss = 4.5;
W = 0;
Ti = 2;
Ma = 2*m_p;

gamma_sh = @(R,T) 0.5 - 0.5*log(2*pi*(m_e/Ma)*(1 + (Ti./T))).*(1 - R);
gamma_ion = @(R,T) 2*(Ti./T).*(1-R) + 2;
gamma_R = @(R,T) (Eiz - W + 0.5*Ediss.*(1-R))./T;

gamma = @(R,T) gamma_sh(R,T) + gamma_ion(R,T) + gamma_R(R,T);

for ii = 1:numel(Rie)
    pwr_13MHz{ii} = e_c*ne_13MHz_mean.*C_s(Te_13MHz_mean,2).*Te_13MHz_mean.*gamma(Rie(ii),Te_13MHz_mean);
    pwr_28GHz{ii} = e_c*ne_28GHz_mean.*C_s(Te_28GHz_mean,2).*Te_28GHz_mean.*gamma(Rie(ii),Te_28GHz_mean);
end

figure; 
hold on
for ii = 1:numel(Rie)
    plot(pwr28GHz,pwr_28GHz{ii}*1e-6,'r.')
end
ylim([0,12])

if saveFig
    figureName = ['Step_2_DlpHeatFlux'];

    % PDF figure:
    exportgraphics(gcf,[figureName,'.pdf'],'Resolution',300) 

    % TIFF figure:
    exportgraphics(gcf,[figureName,'.tiff'],'Resolution',600) 
end
