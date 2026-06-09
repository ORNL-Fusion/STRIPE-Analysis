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

fetchDataFromServer = 1;

if fetchDataFromServer
    DLP = 2.5;  % Let Probe "A" be DLP 2.5
    
    shotlist = [20000 + [736,737,738,739,740,742,744,745,746,747,749,751,752,753,754]];
    R        =          [1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5,7.0,7.5,8.0,8.5];
    
    Stem = '\MPEX::TOP.'; Branch = 'MACHOPS1:'; RootAddress = [Stem,Branch];

    DataAddress{1} = [RootAddress,'LP_V_RAMP']; % V
    DataAddress{2} = [RootAddress,'EA_CURRENT']; % I

    % Setting the Configuration file ==========================================
    Config.tStart = 4.13;
    Config.tEnd   = 4.38;
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
    [niA,TeA,timeA,IfitA,IpA,VpA,tmA,VsweepA,IsweepA,GlitchFlagA,SSQresA,StdResA,StdResNormA] = DLP_fit_V5_5(Config,shotlist,DataAddress);
    for s = 1:length(shotlist)
        NiA{s} = 0.5*(niA{s}{1} + niA{s}{2});
    end
    
    % Setting the Configuration file ==========================================
    ConfigB = Config;
    ConfigB.V_Att = 2; % Output voltage of DLP box (Voltage) = V_Att*Digitized data 
    ConfigB.I_Att = 5; % Output voltage of DLP box (Current) = I_att*Digitized data
    ConfigB.V_cal = [12.05,0.205];   % Voltage output of DLP = V_cal(1)*Output voltage of DLP box + V_cal(2)
    ConfigB.I_cal = [-142.5, 0.015]; % Current output of DLP = (I_cal(2) + Output voltage of DLP box)/Ical(1)

    % Probe "B"
    [niB,TeB,timeB,IfitB,IpB,VpB,tmB,VsweepB,IsweepB,GlitchFlagB,SSQresB,StdResB,StdResNormB] = DLP_fit_V5_5(ConfigB,shotlist,DataAddressB);
    for s = 1:length(shotlist)
        NiB{s} = 0.5*(niB{s}{1} + niB{s}{2});
    end

    % Gathering data from Baratron at spool 2.5
    addressPG{1} = '\MPEX::TOP.MACHOPS1:PG2'; % PG2.5
    addressPG{2} = '\MPEX::TOP.MACHOPS1:PG1'; % PG10.5
    Calshot = 17830; % 2017_11_16
    [pgAcal,t_pgAcal] = my_mdsvalue_v2(Calshot,addressPG(1));
    [pgBcal,t_pgBcal] = my_mdsvalue_v2(Calshot,addressPG(2));
    [pgA,t_pgA] = my_mdsvalue_v2(shotlist,addressPG(1));
    [pgB,t_pgB] = my_mdsvalue_v2(shotlist,addressPG(2));
    
    % In cases were both the calibration and raw data have the same starting
    % time and the same sampling rate (1 kHz), a difference in number of
    % elements can be fixed as follows:
    for s = 1:length(shotlist)
        L1 = length(pgA{1}); L2 = length(pgAcal{1});
        if L1>L2; Ldata = L2; else; Ldata = L1; end
        PGA{s} = ( pgA{s}(1:Ldata)-pgAcal{1}(1:Ldata) )*2;
        t_PGA{s} = t_pgA{s}(1:Ldata);
        PGB{s} = ( pgB{s}(1:Ldata)-pgBcal{1}(1:Ldata) )*2;
        t_PGB{s} = t_pgB{s}(1:Ldata);
    end
    
    % Gathering RF power trace data:
    addressRF{1} = '\MPEX::TOP.MACHOPS1:RF_FWD_PWR';
    [RF,t_RF] = my_mdsvalue_v2(shotlist,addressRF(1));
   
             varList = {'DLPA','DLPB',...
               'shotlist','R',...
               'NiA','TeA','timeA',...
               'NiB','TeB','timeB',...
               'PGA','t_PGA',...
               'PGB','t_PGB',...
               'GlitchFlagA','GlitchFlagB',...
               'StdResNormA','StdResNormB',...
               'RF','t_RF',...
               'niA','IpA','VpA','tmA'};
      
    save('Step_1_GetData_DLP_2017_11_21.mat',varList{:})
else
    load('Step_1_GetData_DLP_2017_11_21.mat')
end


%% Plotting data
close all
C = {'k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:'};
t = t_zero(shotlist);

TimePlotEnd = 4.40;

figure; 
% =========================================================================
subplot(2,3,1); hold on
for s = 1:length(shotlist)
    GoodFitsA{s} = GlitchFlagA{s} == 0 & StdResNormA{s}<=0.4 & NiA{s}>0 & NiA{s}<1e21 & TeA{s}<=20;
    plot(timeA{s},NiA{s},C{s},'lineWidth',1);
    h(s) = plot(timeA{s}(GoodFitsA{s}),NiA{s}(GoodFitsA{s}),C{s},'lineWidth',2);

end
title('$ n_e $ $ [m^{-3}] $','interpreter','Latex','FontSize',13,'Rotation',0)
legend(h,['DLP ',num2str(DLP)],'location','NorthEast')
ylim([0,8e19])
xlim([4.15,TimePlotEnd])
% =========================================================================
subplot(2,3,2); hold on
for s = 1:length(shotlist)
    h(s) = plot(timeA{s}(GoodFitsA{s}),TeA{s}(GoodFitsA{s}),C{s},'lineWidth',2)
    L{s} = [num2str(shotlist(s)),' ,t=',num2str(t{s}(10:14))];
end
legend(h,L,'location','NorthEast')
title('$ T_e $ $ [eV] $','interpreter','Latex','FontSize',13,'Rotation',0)
ylim([0,6])
xlim([4.15,TimePlotEnd])
% =========================================================================
subplot(2,3,4); hold on
for s = 1:length(shotlist)
    GoodFitsB{s} = GlitchFlagB{s} == 0 & StdResNormB{s}<=0.4 & NiB{s}>0 & NiB{s}<1e21 & TeB{s}<=20;
    plot(timeB{s},NiB{s},C{s},'lineWidth',1);
    h(s) = plot(timeB{s}(GoodFitsB{s}),NiB{s}(GoodFitsB{s}),C{s},'lineWidth',2);

end
title('$ n_e $ $ [m^{-3}] $','interpreter','Latex','FontSize',13,'Rotation',0)
legend(h,['DLP ',num2str(DLPB)],'location','NorthEast')
ylim([0,7e19])
xlim([4.15,TimePlotEnd])
% =========================================================================
subplot(2,3,5); hold on
for s = 1:length(shotlist)
    h(s) = plot(timeB{s}(GoodFitsB{s}),TeB{s}(GoodFitsB{s}),C{s},'lineWidth',2)
end
title('$ T_e $ $ [eV] $','interpreter','Latex','FontSize',13,'Rotation',0)
ylim([0,6])
xlim([4.15,TimePlotEnd])
% ========================================================================
subplot(2,3,3); hold on
for s = 1:length(shotlist)
    h(s) = plot(t_PGA{s},PGA{s},C{s},'lineWidth',2)
end
title('$ PG2.5 $ $ [Pa] $','interpreter','Latex','FontSize',13,'Rotation',0)
ylim([0,25])
xlim([4,4.5])

subplot(2,3,6); hold on
for s = 1:length(shotlist)
    h(s) = plot(t_PGB{s},PGB{s} - mean(PGB{s}(find(t_PGB{s}>4.1 & t_PGB{s}<4.14))),C{s},'lineWidth',2)
end
title('$ PG10.5 $ $ [Pa] $','interpreter','Latex','FontSize',13,'Rotation',0)
ylim([0,3])
xlim([4,4.5])

set(gcf,'color','w')

%% Extracting steady state data:
% #########################################################################
% Steady state values
% #########################################################################

for s = 1:length(shotlist)
    rngneA = [timeA{s}>=4.31 & timeA{s}<=4.32]; 
    rngteA = [timeA{s}>=4.31 & timeA{s}<=4.33]; 
    rngPGA = [t_PGA{s}>=4.31 & t_PGA{s}<=4.32]; 
    rngneB = [timeB{s}>=4.31 & timeB{s}<=4.32];
    rngPGB = [t_PGB{s}>=4.31 & t_PGB{s}<=4.32]; 

    neA(s)  = mean(NiA{s}(rngneA & GoodFitsA{s}));
    dneA(s) = std(NiA{s}(rngneA & GoodFitsA{s}),1,2);
    teA(s)  = mean(TeA{s}(rngteA & GoodFitsA{s}));
    dteA(s) = std(TeA{s}(rngteA & GoodFitsA{s}),1,2);
    
    neB(s)  = mean(NiB{s}(rngneB & GoodFitsB{s}));
    dneB(s) = std(NiB{s}(rngneB & GoodFitsB{s}),1,2);
    teB(s)  = mean(TeB{s}(rngneB & GoodFitsB{s}));
    dteB(s) = std(TeB{s}(rngneB & GoodFitsB{s}),1,2);
    
    pA(s) = mean(PGA{s}(rngPGA));
    dpA(s) = std(PGA{s}(rngPGA),1,1);
    pB(s) = mean(PGB{s}(rngPGB));
    dpB(s) = std(PGB{s}(rngPGB),1,1);
end

%%
[a,b] = sort(R);

figure; 
subplot(2,1,1); hold on
hne(1) = errorbar(R(b),neA(b),dneA(b),'LineWidth',3,'color','r')
hne(2) = errorbar(R(b),neB(b),dneB(b),'LineWidth',3,'color','k')

ylim([0,8]*1e19)
xlim([-7.5,7.5]); 
box on
t = title('$n_e$ $[m^{-3}]$ '); t.Interpreter = 'latex'; t.FontSize = 13;
t = xlabel('$R$ $[cm]$ '); t.Interpreter = 'latex'; t.FontSize = 11;

subplot(2,1,2); hold on
hte(1) = errorbar(R(b),teA(b),dteA(b),'LineWidth',3,'color','r')
hte(2) = errorbar(R(b),teB(b),dteB(b),'LineWidth',3,'color','k')
hpA(3) = errorbar(R(b),pA(b),dpA(b),'LineWidth',1,'color','r')
hpB(4) = errorbar(R(b),pB(b),dpB(b),'LineWidth',1,'color','k')

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

figure
hold on
for ss = 1:numel(NiA)
    plot3(timeA{ss},R(ss)*ones(size(NiA{ss})),NiA{ss},'LineWidth',2)
end
zlim([0,1e20])
ylim([0,8])
xlim([4.1,4.2])

figure
hold on
for ss = 1:numel(NiA)
    plot3(tmA{ss},R(ss)*ones(size(IpA{ss})),IpA{ss},'LineWidth',2)
end
zlim([-0.5,0.5])
ylim([0,8])
xlim([4.1,4.2])

% The data seems to ahve empty spaces, lets look at the data individually:

figure
for ss = 1:numel(shotlist)
    subplot(3,4,ss)
    hold on
   
    % Data with glitches:
    plot(timeA{ss},NiA{ss},'k.')
    
    % Remove some of the glithces:
    GoodFitsA{ss} = GlitchFlagA{ss} == 0 & StdResNormA{ss}<=0.4 & NiA{ss}>0 & NiA{ss}<1e21 & TeA{ss}<=20;
    plot(timeA{ss}(GoodFitsA{ss}),NiA{ss}(GoodFitsA{ss}),'ro')
    
    ylim([-1,10]*1e19)
    xlim([4.13,4.19])
end

% It seems that for t <= 4.19 sec we can use the data as it comes out of
% the DLP routine:

% Replace all empty locations and NaNs with the mean of the neihbor
% elements:
for ss = 1:numel(shotlist)
    for ii = 1:numel(NiA{ss})
         if (NiA{ss}(ii)) < 0 || NiA{ss}(ii) > 1e22 || isempty(NiA{ss}(ii))
            NiA{ss}(ii) = 0;
         end
        if isnan(NiA{ss}(ii))
            NiA{ss}(ii) = 0;
        end
       if (NiA{ss}(ii)) == 0  && ii >1
            NiA{ss}(ii) = NiA{ss}(ii-1);
       end
        
        if (TeA{ss}(ii)) < 0 || TeA{ss}(ii) > 15 || isempty(TeA{ss}(ii))
            TeA{ss}(ii) = 0;
         end
        if isnan(TeA{ss}(ii))
            TeA{ss}(ii) = 0;
        end
       if (TeA{ss}(ii)) == 0  && ii >1
            TeA{ss}(ii) = TeA{ss}(ii-1);
       end
        
    end
end

t_interp = linspace(4.135,4.195);
for ss = 1:numel(shotlist)
    neA_interp(ss,:) = interp1(timeA{ss},NiA{ss},t_interp,'cubic');
    TeA_interp(ss,:) = interp1(timeA{ss},TeA{ss},t_interp,'cubic');
end

[RR,TT] = meshgrid(R,t_interp);

t_rfStart = 4.1538;

figure('color','w');
hold on
surf(RR,TT,neA_interp','LineStyle','none')
hL = line([0,8],[1,1]*t_rfStart,[1,1]*1e18,'color','r','LineWidth',2);
zlim([-1,10]*1e19)
xlim([0,max(R)])
set(gca,'FontName','Times','FontSize',11)
xlabel('r [cm]','Interpreter','latex','FontSize',14)
ylabel('time [s]','Interpreter','latex','FontSize',14)
zlabel('$n_e$ [m$^{-3}$]','Interpreter','latex','FontSize',14)

% The data from above seems to suggest that the edge density builds up
% first; however when we inspect the ion saturation current traces, it
% appears that data in the core is compromized due to the presence of
% glitches and so we cannot really prove that the edge plasma builds up
% first.
