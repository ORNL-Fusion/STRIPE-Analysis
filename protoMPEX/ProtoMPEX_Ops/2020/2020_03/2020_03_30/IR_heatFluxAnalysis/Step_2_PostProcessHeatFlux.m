% Step 2: 
% Post process heat flux data

clear all
close all

% Load data:
% =========================================================================
load('Step_1b_heatFlux2D.mat');

% Extract data from variable:
% =========================================================================
shot       = data.shot;
q_star     = data.q_star; % characteristic heat flux in [W]
q_IJR      = data.q_IJR;  % Normalized 2D heat flux
pwr_Target = data.pwr;    % Power coupled to Target in [W]
pwrECH     = data.pwrECH*1e3; % 28 GHz power in [W]
xx         = data.xx; % "x" coordinate in [m]
yy         = data.yy; % "y" coordinate in [m]
t_R        = data.t_R; % Normalized time coordinate in [s]
t_star     = data.t_star; % Characteristic time in [s]

%% Plot data:
% =========================================================================
close all

% Time traces of coupled power:
% -------------------------------------------------------------------------
figure('color','w')
hold on

fontSize.label = 13;
fontSize.axes = 12;

colorRng = {'k','bl','c','m','g','c','k','bl','c','m','g','c'};

t0 = 0.31;
tRF_start = 4.15;

for ii = numel(shot):-1:1
    % Time trace:
    t_pwr_Target{ii} = t_R{ii}*t_star - t0 + tRF_start;
    
    % Time trace of power coupled to Target:
    hPwrTarget(ii) = plot(t_pwr_Target{ii},pwr_Target{ii}*1e-3,'LineWidth',2,'Color',colorRng{ii});
    
    % Find peak power for ECH:
    rng = find(t_pwr_Target{ii} > 4.45 & t_pwr_Target{ii} < 4.58);
    [pwr_Target_ech(ii),ipeak_ech] = max(pwr_Target{ii}(rng));
    ipeak_ech = rng(ipeak_ech);
    
    % Find Helicon power prior to ECH pulse:
    rng = find(t_pwr_Target{ii} > 4.40 & t_pwr_Target{ii} < 4.45);
    [pwr_Target_helicon(ii),ipeak_helicon] = max(pwr_Target{ii}(rng));
    ipeak_helicon = rng(ipeak_helicon);
    
    % Plot peak:
    hPwrPeak = plot(t_pwr_Target{ii}(ipeak_ech),pwr_Target_ech(ii)*1e-3);
    set(hPwrPeak,'Marker','o','MarkerFaceColor',colorRng{ii})
    
    % ECH power contribution
    pwr_ech = pwr_Target_ech - pwr_Target_helicon;
end
ylim([0,16])
ylabel('Power to target [kW]','Interpreter','Latex','FontSize',fontSize.label)
xlabel('time [s]','Interpreter','Latex','FontSize',fontSize.label)

box on
set(gca,'fontName','Times','fontSize',fontSize.axes)

% Save figure:
figureName = ['step_2_timeEvolutionPowerToTarget'];
saveas(gcf, figureName,'tiff');

% Scaling of power coupled to Target with ECH power:
% -------------------------------------------------------------------------
% Linear fit:
x = linspace(0,100);

figure('color','w')
hold on
hPeakPwr = plot(pwrECH*1e-3,pwr_ech*1e-3,'ksq');
hlinFit = plot(x,polyval(polyfit(pwrECH*1e-3,pwr_ech*1e-3,1),x),'r','LineWidth',2);

% Labels:
ylabel('Power to Target [kW]','Interpreter','Latex','FontSize',fontSize.label)
xlabel('28 GHz [kW]','Interpreter','Latex','FontSize',fontSize.label)

% Formatting:
set(hPeakPwr,'Markersize',7,'LineWidth',2)
xlim([0,80])
ylim([0,14])
box on
set(gca,'fontName','Times','fontSize',fontSize.axes)

% Save figure:
figureName = ['step_2_PeakPowerToTarget'];
saveas(gcf, figureName,'tiff');

% Plotting efficiency:
% -------------------------------------------------------------------------
figure('color','w')
hEff = plot(pwrECH*1e-3,100*pwr_ech./pwrECH,'ksq');

% Labels:
ylabel('Efficiency [\%]','Interpreter','Latex','FontSize',fontSize.label)
xlabel('28 GHz [kW]','Interpreter','Latex','FontSize',fontSize.label)

% Formatting:
set(hEff,'Markersize',7,'LineWidth',2)
xlim([0,80])
ylim([0,30])
box on
set(gca,'fontName','Times','fontSize',fontSize.axes)

% Save figure:
figureName = ['step_2_HeatingEfficiency'];
saveas(gcf, figureName,'tiff');

% This represents the fraction of the injected ECH power to the power
% coupled to the target via plasma transport.
% Where does the other 
