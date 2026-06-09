% Step 2: 
% Post process heat flux data

clear all
close all

% Load data:
% =========================================================================
load('Step_2c_heatFlux2D.mat');

% Rename and rescale some variables:
% =========================================================================
pwr_Target = pwr;    % Power coupled to Target in [W]
pwr28GHz   = pwr28GHz*1e3; % 28 GHz power in [W]

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
figureName = ['step_2d_timeEvolutionPowerToTarget'];
saveas(gcf, figureName,'tiff');

% Scaling of power coupled to Target with ECH power:
% -------------------------------------------------------------------------
% Linear fit:
x = linspace(0,100);

figure('color','w')
hold on
hPeakPwr = plot(pwr28GHz*1e-3,pwr_ech*1e-3,'ksq');
hlinFit = plot(x,polyval(polyfit(pwr28GHz*1e-3,pwr_ech*1e-3,1),x),'r','LineWidth',2);

% Labels:
ylabel('Power to Target [kW]','Interpreter','Latex','FontSize',fontSize.label)
xlabel('28 GHz [kW]','Interpreter','Latex','FontSize',fontSize.label)

% Formatting:
set(hPeakPwr,'Markersize',9,'LineWidth',2)
xlim([0,80])
ylim([0,14])
box on
set(gca,'fontName','Times','fontSize',fontSize.axes)

% Save figure:
figureName = ['step_2d_PeakPowerToTarget'];
saveas(gcf, figureName,'tiff');

% Plotting efficiency:
% -------------------------------------------------------------------------
figure('color','w')
hEff = plot(pwr28GHz*1e-3,100*pwr_ech./pwr28GHz,'ksq');

% Labels:
ylabel('Efficiency [\%]','Interpreter','Latex','FontSize',fontSize.label)
xlabel('28 GHz [kW]','Interpreter','Latex','FontSize',fontSize.label)

% Formatting:
set(hEff,'Markersize',9,'LineWidth',2)
xlim([0,80])
ylim([0,30])
box on
set(gca,'fontName','Times','fontSize',fontSize.axes)

% Save figure:
figureName = ['step_2d_HeatingEfficiency'];
saveas(gcf, figureName,'tiff');

% This represents the fraction of the injected ECH power to the power
% coupled to the target via plasma transport.
% Where does the other 
