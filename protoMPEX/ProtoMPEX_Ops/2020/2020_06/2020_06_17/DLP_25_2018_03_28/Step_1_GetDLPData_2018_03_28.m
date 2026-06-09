% DLP 2.5 data radial scan for MAB conditions

% CONTEXT:
% It has been observed that inserting the probe near the helicon antenna
% adversily affects the operation of the plasma. To characterize this we
% have measured the plasma parameters at spool 2.5 and spool 10.5
% simulteaneoulsy. The measurements at spool 10.5 and the gas pressure
% increase at spool 2.5 provide a measure of how much the probe is
% affecting the plasma performance.

% 2020_06_28:
% We need this data for Josh. He is interested in the radial plasma density
% during the first 20-40 ms of the RF pulse:

clear all
close all

fetchDataFromServer = 0;

if fetchDataFromServer
    DLP = 2.5;  % Let Probe "A" be DLP 2.5
    
    shotlist = [20000 + [736,737,738,739,740,742,744,745,746,747,749,751,752,753,754]];
    R        =          [1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5,7.0,7.5,8.0,8.5];
    
    Stem = '\MPEX::TOP.'; Branch = 'MACHOPS1:'; RootAddress = [Stem,Branch];

    DataAddress{1} = [RootAddress,'LP_V_RAMP']; % V
    DataAddress{2} = [RootAddress,'EA_CURRENT']; % I

    % Setting the Configuration file ==========================================
    Config.tStart = 4.13;
    Config.tEnd   = 4.7;
    Config.FitFunction = 2;
    Config.AreaType = 1; % 1: Cylindrical + cap
    Config.V_Att = 1; % Output voltage of DLP box (Voltage) = V_Att*Digitized data 
    Config.I_Att = 1; % Output voltage of DLP box (Current) = I_att*Digitized data
    Config.V_cal = [(0.46e-3)^-1,0]; % Voltage output of DLP = V_cal(1)*Output voltage of DLP box + V_cal(2) %  2.1739e+03
    Config.I_cal = [-1,0]; % Current output of DLP = (I_cal(2) + Output voltage of DLP box)/Ical(1)
    Config.Center_V = 1;
    Config.Center_I = 1; 
    Config.L_tip = 1.7/1000;
    Config.D_tip = 0.254/1000; % [m]
    Config.FilterDataInput = 1; % Filter input data with savitsky Golay filter order 3 frame size 7
    Config.SGF = 11;
    Config.TimeMode = 2; % Effective time of sweep, (1) start of ramp or (2) mean time of ramp
    Config.AMU = 2; % Ion mass in AMU

    % Probe "A"
    [ni,Te,time,Ifit,Ip,Vp,tm,Vsweep,Isweep,GlitchFlag,SSQres,StdRes,StdResNorm] = DLP_fit_V5_5(Config,shotlist,DataAddress);
    for s = 1:length(shotlist)
        Ni{s} = 0.5*(ni{s}{1} + ni{s}{2});
    end  
  
    % Gathering RF power trace data:
    addressRF{1} = '\MPEX::TOP.MACHOPS1:RF_FWD_PWR';
    [RF,t_RF] = my_mdsvalue_v2(shotlist,addressRF(1));
   
             varList = {'DLP',...
               'shotlist','R',...
               'Ni','Te','time',...
               'GlitchFlag','GlitchFlag',...
               'StdResNorm','StdResNorm',...
               'RF','t_RF',...
               'ni','Ip','Vp','tm'};
      
    save('Step_1_GetData_DLP_2018_03_29.mat',varList{:})
else
    load('Step_1_GetData_DLP_2018_03_29.mat')
end


%% Plotting data
close all
C = {'k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:'};
t = t_zero(shotlist);

TimePlotEnd = 4.40;

figure; 
% =========================================================================
subplot(2,1,1); hold on
for s = 1:length(shotlist)
    GoodFits{s} = GlitchFlag{s} == 0 & StdResNorm{s}<=5 & Ni{s}>0 & Ni{s}<1e21 & Te{s}<=20;
%     plot(time{s},Ni{s},C{s},'lineWidth',1);
    h(s) = plot(time{s}(GoodFits{s}),Ni{s}(GoodFits{s}),C{s},'lineWidth',2);

end
title('$ n_e $ $ [m^{-3}] $','interpreter','Latex','FontSize',13,'Rotation',0)
legend(h,['DLP ',num2str(DLP)],'location','NorthEast')
ylim([0,8e19])
xlim([4.15,TimePlotEnd])
% =========================================================================
subplot(2,1,2); hold on
for s = 1:length(shotlist)
    h(s) = plot(time{s}(GoodFits{s}),Te{s}(GoodFits{s}),C{s},'lineWidth',2)
    L{s} = [num2str(shotlist(s)),' ,t=',num2str(t{s}(10:14))];
end
legend(h,L,'location','NorthEast')
title('$ T_e $ $ [eV] $','interpreter','Latex','FontSize',13,'Rotation',0)
ylim([0,6])
xlim([4.15,TimePlotEnd])
% =========================================================================


%% Extracting steady state data:
% #########################################################################
% Steady state values
% #########################################################################

for s = 1:length(shotlist)
    rngne = [time{s}>=4.31 & time{s}<=4.32]; 
    rngte = [time{s}>=4.31 & time{s}<=4.33]; 

    ne(s)  = mean(Ni{s}(rngne & GoodFits{s}));
    dne(s) = std(Ni{s}(rngne & GoodFits{s}),1,2);
    te(s)  = mean(Te{s}(rngte & GoodFits{s}));
    dte(s) = std(Te{s}(rngte & GoodFits{s}),1,2);
end

%%
[a,b] = sort(R);

figure; 
subplot(2,1,1); hold on
hne(1) = errorbar(R(b),ne(b),dne(b),'LineWidth',3,'color','r')

ylim([0,8]*1e19)
xlim([-7.5,7.5]); 
box on
t = title('$n_e$ $[m^{-3}]$ '); t.Interpreter = 'latex'; t.FontSize = 13;
t = xlabel('$R$ $[cm]$ '); t.Interpreter = 'latex'; t.FontSize = 11;

subplot(2,1,2); hold on
hte(1) = errorbar(R(b),te(b),dte(b),'LineWidth',3,'color','r')

ylim([0,8])
xlim([-7.5,7.5]); 
box on
t = title('$T_e$ $[eV]$ '); t.Interpreter = 'latex'; t.FontSize = 13;
t = xlabel('$R$ $[cm]$ '); t.Interpreter = 'latex'; t.FontSize = 11;

set(gcf,'color','w')


%% Time resolved data:
% Since each shot does not have the same numnber of points, we will need to
% divide time in intervals and then interpolate data into the global time
% trace

% Plot data to get a feel for the glitches and radial profiles:
% =========================================================================
figure
hold on
for ss = 1:numel(Ni)
    plot3(time{ss},8.5 - R(ss)*ones(size(Ni{ss})),Ni{ss},'LineWidth',2)
end
zlim([0,1e20])
ylim([0,8])
xlim([4.1,4.2])

figure
hold on
for ss = 1:numel(Ni)
    plot3(tm{ss},8.5 - R(ss)*ones(size(Ip{ss})) - 7.5,Ip{ss},'LineWidth',2)
end
zlim([-0.5,0.5])
ylim([0,8])
xlim([4.1,4.2])

% The data seems to ahve empty spaces, lets look at the data individually:
% =========================================================================
figure
for ss = 1:numel(shotlist)
    subplot(3,5,ss)
    hold on
   
    % Data with glitches:
    plot(time{ss},Ni{ss},'k.')
    
    % Remove some of the glithces:
    GoodFits{ss} = GlitchFlag{ss} == 0 & StdResNorm{ss}<=0.4 & Ni{ss}>0 & Ni{ss}<1e21 & Te{ss}<=20;
    plot(time{ss}(GoodFits{ss}),Ni{ss}(GoodFits{ss}),'ro')
    
    ylim([-1,10]*1e19)
    xlim([4.13,4.19])
end

% Replace all empty locations and NaNs with the mean of the neihbor
% elements:
% =========================================================================
for ss = 1:numel(shotlist)
    for ii = 1:numel(Ni{ss})
         if (Ni{ss}(ii)) < 0 || Ni{ss}(ii) > 1e20 || isempty(Ni{ss}(ii))
            Ni{ss}(ii) = 0;
         end
        if isnan(Ni{ss}(ii))
            Ni{ss}(ii) = 0;
        end
       if (Ni{ss}(ii)) == 0  && ii >1
            Ni{ss}(ii) = Ni{ss}(ii-1);
       end
        
        if (Te{ss}(ii)) < 0 || Te{ss}(ii) > 15 || isempty(Te{ss}(ii))
            Te{ss}(ii) = 0;
         end
        if isnan(Te{ss}(ii))
            Te{ss}(ii) = 0;
        end
       if (Te{ss}(ii)) == 0  && ii >1
            Te{ss}(ii) = Te{ss}(ii-1);
       end
        
    end
end

% Interpolate data to create a 2D array for the density and temperature:
% =========================================================================
t_interp = linspace(4.135,4.24);
for ss = 1:numel(shotlist)
    ne_interp(ss,:) = interp1(time{ss},Ni{ss},t_interp,'cubic');
    Te_interp(ss,:) = interp1(time{ss},Te{ss},t_interp,'cubic');
end

[RR,TT] = meshgrid(8.5 - R,t_interp);

% Determine the start of the RF pulse:
% =========================================================================
figure; hold on; for ii = 1:numel(shotlist); plot(t_RF{ii}(1:end-1),-RF{ii}); end
ylim([-0.1,0.6])
xlim([4.12,4.2])
t_rfStart = 4.1584;

% Plot time-resolved radial data:
% =========================================================================
figure('color','w');
hold on
surf(RR,TT,ne_interp','LineStyle','none')
hL = line([0,7.5],[1,1]*t_rfStart,[1,1]*1e18,'color','r','LineWidth',2);
zlim([-1,10]*1e19)
xlim([0,7.5])
set(gca,'FontName','Times','FontSize',11)
xlabel('r [cm]','Interpreter','latex','FontSize',14)
ylabel('time [s]','Interpreter','latex','FontSize',14)
zlabel('$n_e$ [m$^{-3}$]','Interpreter','latex','FontSize',14)
title('DLP 2.5 time-resolved radial scan, 2018-03-28','Interpreter','latex','FontSize',14)
view([30,60])
colorbar
hLeg = legend(hL,'Start of RF');
hLeg.Location = 'best'

saveFig = 1;

if saveFig
    saveas(gcf,'Step_0_TimeResolvedRadialScan_DLP_2','tiffn')
end

% Reference Magnetic flux at limiter:
B0 = 0.04807;  % [T] 
r0 = 6.256;  %[cm]
phi0 = B0*r0^2;

% Magnetic field at probe locations:
Bprobe = 0.05487; % [T]

% Flux coordinates:
Xi  = Bprobe*(RR.^2)/phi0; 

figure('color','w');
hold on
surf(sqrt(Xi),TT,ne_interp','LineStyle','none')
hL = line([0,7.5],[1,1]*t_rfStart,[1,1]*1e18,'color','r','LineWidth',2);
hXi = line([1,1],[4.15,4.2],[1,1]*1e18,'Color','r','LineWidth',2,'LineStyle','--');
zlim([-1,10]*1e19)
xlim([0,1.3])
ylim([4.15,4.18])
set(gca,'FontName','Times','FontSize',11)
xlabel('$\sqrt{\chi}$','interpreter','latex','Fontsize',12)
ylabel('time [s]','Interpreter','latex','FontSize',14)
zlabel('$n_e$ [m$^{-3}$]','Interpreter','latex','FontSize',14)
title('DLP 2.5 time-resolved radial scan, 2018-03-28','Interpreter','latex','FontSize',14)
view([30,60])
colorbar
hLeg = legend(hL,'Start of RF');
hLeg.Location = 'best';

saveFig = 1;

if saveFig
    saveas(gcf,'Step_1_TimeResolvedRadialScan_DLP_2_Xi','tiffn')
end


%% Electron temperature:
% =========================================================================
figure('color','w');
hold on
surf(RR,TT,Te_interp','LineStyle','none')
hL = line([0,7.5],[1,1]*t_rfStart,[1,1]*1,'color','r','LineWidth',2);
zlim([0,8])
xlim([0,7.5])
set(gca,'FontName','Times','FontSize',11)
xlabel('r [cm]','Interpreter','latex','FontSize',14)
ylabel('time [s]','Interpreter','latex','FontSize',14)
zlabel('$T_e$ [eV]','Interpreter','latex','FontSize',14)
title('DLP 2.5 time-resolved radial scan, 2018-03-28','Interpreter','latex','FontSize',14)
view([30,60])
colorbar
hLeg = legend(hL,'Start of RF');
hLeg.Location = 'best'

saveFig = 1;

if saveFig
    saveas(gcf,'Step_0_TimeResolvedRadialScan_Te_DLP_2','tiffn')
end

% Reference Magnetic flux at limiter:
B0 = 0.04807;  % [T] 
r0 = 6.256;  %[cm]
phi0 = B0*r0^2;

% Magnetic field at probe locations:
Bprobe = 0.05487; % [T]

% Flux coordinates:
Xi  = Bprobe*(RR.^2)/phi0; 

figure('color','w');
hold on
surf(sqrt(Xi),TT,Te_interp','LineStyle','none')
hL = line([0,7.5],[1,1]*t_rfStart,[1,1]*12,'color','r','LineWidth',2);
hXi = line([1,1],[4.15,4.2],[1,1],'Color','r','LineWidth',2,'LineStyle','--');
zlim([0,12])
xlim([0,1.3])
ylim([4.15,4.24])
set(gca,'FontName','Times','FontSize',11)
xlabel('$\sqrt{\chi}$','interpreter','latex','Fontsize',12)
ylabel('time [s]','Interpreter','latex','FontSize',14)
zlabel('$T_e$ [eV]','Interpreter','latex','FontSize',14)
title('DLP 2.5 time-resolved radial scan, 2018-03-28','Interpreter','latex','FontSize',14)
view([30,60])
colorbar
hLeg = legend(hL,'Start of RF');
hLeg.Location = 'best';

saveFig = 1;

if saveFig
    saveas(gcf,'Step_1_TimeResolvedRadialScan_Te_DLP_2_Xi','tiffn')
end


%%  SAVE DATA:
% Assemble data package:
% =========================================================================

% DLP 2.5 data
DLP2.ne = ne_interp';
DLP2.Te = Te_interp';
DLP2.RR = RR;
DLP2.TT = TT;
DLP2.Xi = Xi;
DLP2.shotlist = shotlist;
DLP2.rprobe   = 8.5-R;
DLP2.dateOfExperiment  = '2018_03_28';
DLP2.comment{1} = 'DLP 2.5 data, Window-limiter magnetic configuration, see shot summaries for 2018-03-28';
DLP2.comment{2} = 'Te data is very noisy, I would not use it, ne data I do trust';
DLP2.t_rfStart = t_rfStart;
DLP2.zLoc = 'Spool 2.5';

% Assemble list of variables to save:
varlistSave = {'DLP2'};
save('DLP_2_ProbeData_2018_03_28.mat',varlistSave{:})

