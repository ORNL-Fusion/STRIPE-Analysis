% DLP 2.5 data radial scan for MAB conditions

% CONTEXT:
% It has been observed that inserting the probe near the helicon antenna
% adversily affects the operation of the plasma. To characterize this we
% have measured the plasma parameters at spool 2.5 and spool 10.5
% simulteaneoulsy. The measurements at spool 10.5 and the gas pressure
% increase at spool 2.5 provide a measure of how much the probe is
% affecting the plasma performance.

clear all
close all

DLPA = 2.5;  % Let Probe "A" be DLP 2.5
DLPB = 10.5; % Let Probe "B" be DLP 10.5

shotlist = 17000 + [952:956];
shotlist = 17000 + [930]; % RF pulse length was 500 ms
% =========================================================================
shotlist = 17000 + [930,931]; % Change RF pulse to 200 ms
% these shots show that when Probe "A" is out or at the edge of the plasma,
% the plasma density at the target is about 4e19 and the gas pressure at
% spool 2.5 shows significant pumping
% =========================================================================
shotlist = 17000 + [931:936]; 
% Here we observe the effect of inserting DLP 2.5 into the axis of the
% plasma while observing the effect it has on DLP 10.5 and baratrons at
% spool 2.5 and 10.5. 
% we observe that inserting DLP 2.5 into the center of the plasma has the
% follwing effects:
% - Decrease the plasma density measured at DLP 10.5 from 4 to 3e19 m-3
% - Decrease the amount of gas pressure produced at the target presumably
% due to plasma-surface recombination
% - Increase in neutral gas pressure at spool 2.5 which corresponds to a
% reduction in the neutral gas pumping by the plasma
% THE MAIN RESULT IS that inserting DLP 2.5 in the center of the plasma
% significantly reduces the amount of plasma that reaches the target and
% significantly modifies the neutral gas pressure at the helicon source.
% these changes indicate the presence of the probe at spool 2.5 affects the
% plasma such that we essentially measure a different discharge
% =========================================================================
shotlist = 17000 + [931,937,938];
% Here we compare having DLP 2.5 fully out and fully in the plasma.
% we the see the formation of very large low freqyency oscillations with
% high levels that exceed 100 mA, this corresponds to very high ne
% - the plasma density at the target is reduced by about 50%
%       From 4 to 2e19 m-3
% - Neutral gas pressure at the target is reduced by 40%
% - Neutral gas pressure at the source is increased by 50%
% =========================================================================
shotlist = 17000 + [ 931, 939, 941, 942, 944, 945, 946, 947,948,949,950,951,952];
R        =         [ 7.5,-3.5,-3.0,-2.5,-2.0,-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.0];
% =========================================================================
shotlist = 17000 + [ 931,952,954,955,957,960];
R        =         [ 7.5,2.0,2.5,3.0,3.5,4.0];
% In this data set we see the first instance when DLP 2.5 only weakly
% affects the plasma as it can be seen from DLP 10.5 response (ne and gas)
% =========================================================================
shotlist = 17000 + [ 931,957,958,960];
R        =         [ 7.5,3.5,3.5,4.0];
% =========================================================================
shotlist = 17000 + [956,958,960,962,963,964,965,967,968,969 ,970 ];
R        =         [3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5,7.0,6.75,7.25];
shotlist = 17000 + [958,960,962,963,964,965,967,968,969 ,970 ];
R        =         [3.5,4.0,4.5,5.0,5.5,6.0,6.5,7.0,6.75,7.25];
% This is the data set required to produce the partial radial scan
% =========================================================================
% shotlist = 17000 + [956,958,960,962,973,974:977];
% R        =         [3.0,3.5,4.0,4.5,2.5,2.0];

FitShow = 1; 
RawDataShow = 1;

Stem = '\MPEX::TOP.'; Branch = 'MACHOPS1:'; RootAddress = [Stem,Branch];

DataAddressA{1} = [RootAddress,'INT_4MM_2']; % V
DataAddressA{2} = [RootAddress,'INT_4MM_1']; % I
DataAddressB{1} = [RootAddress,'LP_V_RAMP']; % V
DataAddressB{2} = [RootAddress,'TARGET_LP']; % I

% Setting the Configuration file ==========================================
ConfigA.tStart = 4.13;
ConfigA.tEnd   = 4.38;
ConfigA.FitFunction = 2;
ConfigA.AreaType = 1; % 1: Cylindrical + cap
ConfigA.V_Att = 1; % Output voltage of DLP box (Voltage) = V_Att*Digitized data 
ConfigA.I_Att = 1; % Output voltage of DLP box (Current) = I_att*Digitized data
ConfigA.V_cal = [(0.46e-3)^-1,0]; % Voltage output of DLP = V_cal(1)*Output voltage of DLP box + V_cal(2) %  2.1739e+03
ConfigA.I_cal = [-1,0]; % Current output of DLP = (I_cal(2) + Output voltage of DLP box)/Ical(1)
ConfigA.Center_V = 1;
ConfigA.Center_I = 1; 
ConfigA.L_tip = 1.7/1000;
ConfigA.D_tip = 0.254/1000; % [m]
ConfigA.FilterDataInput = 1; % Filter input data with savitsky Golay filter order 3 frame size 7
ConfigA.SGF = 11;
ConfigA.TimeMode = 2; % Effective time of sweep, (1) start of ramp or (2) mean time of ramp
ConfigA.AMU = 2; % Ion mass in AMU


% Probe "A"
[niA,TeA,timeA,IfitA,IpA,VpA,tmA,VsweepA,IsweepA,GlitchFlagA,SSQresA,StdResA,StdResNormA] = DLP_fit_V5_5(ConfigA,shotlist,DataAddressA);
for s = 1:length(shotlist)
    NiA{s} = 0.5*(niA{s}{1} + niA{s}{2});
end

% Setting the Configuration file ==========================================
ConfigB = ConfigA;
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

% Gathering RF


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
legend(h,['DLP ',num2str(DLPA)],'location','NorthEast')
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

%% raw data
if RawDataShow
    figure; hold on
    for s = 1:length(shotlist)
    plot(tmA{s},VpA{s})
    h(s) = plot(tmA{s},(IpA{s})*1000);
    ylim([-100,100])
    xlim([4.15,TimePlotEnd])
    grid on
    end
    legend(h,num2str(shotlist'))
    set(gcf,'color','w')
    box on
end
%% Fits on Probe A
if FitShow
    s = 1;
    kTotal = round(length(TeA{s})/25);
    for k = 1:kTotal
         figure;
    for c = (1 + (k-1)*25):(25 + (k-1)*25);
        subplot(5,5,c- (k-1)*25 ); hold on
        try
        plot(VsweepA{s}{c},IsweepA{s}{c}*1e3,'k')
        if GoodFitsA{s}(c) == 1
            plot(VsweepA{s}{c},IfitA{s}{c}*1e3,'r')
        else
            plot(VsweepA{s}{c},IfitA{s}{c}*1e3,'g')
        end
        ht = title(['T_e: ',num2str(TeA{s}(c),2),' ,t = ',num2str(timeA{s}(c),5),' s']);
        set(ht,'FontSize',5); grid on
%         ylim(Ymax*[-1,1])
        catch
            warning('error')
            continue
        end
    end
    end
end
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
return
%%
% close all
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

%% Convert data into Excel
% Save data to Excel spreadsheet:
% Create table to export to EXCEL
if 0
G = [shotlist(b)',R(b)',neA(b)',dneA(b)',teA(b)',dteA(b)'];
F = {'Shot','R [cm]','ne[m^-3]','dne','Te[eV]','dTe'};
FileName = 'NeTe_Spool_2_5_2017_11_21.xlsx';
xlswrite(FileName,[F;num2cell(G)]);
end