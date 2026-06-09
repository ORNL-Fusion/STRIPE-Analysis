clear all
close all

load('NeTe_2016_12_19_to_2017_07_20.mat')

try 
    s_start;
catch
    s_start = 1;
    s_end = length(Shot);
end

%% STEADY STATE DATA: =====================================================
for s = s_start:s_end 
    
    if s == find(Shot == 14178)
        s
    end
    
    % REMOVE BAD SHOTS: ===================================================
    % 1- Reject all shots whose data is corrupted: ============================
    % this is, when the DLP code failed

    try % DataCorruptedShots may not have been saved in the latest run
        if DataCorruptedShots(s) || isnan(Shot(s))
            continue
        end
    end

    % 2 - Reject all shots with bad fits: =====================================
    % DLP code did not fail but produced fits with high residuals. We use the
    % value of StdResNorm as the metric to judge the fit.
    condition = find(StdResNorm{s}<(10/100) & StdResNorm{s} > 0 & StdRes{s} == real(StdRes{s}));
    rngTime = find(time{s}>=Gstart(s) & time{s}<=Gend(s));
    rngMean{s} = intersect(rngTime,condition);
    if length(rngMean{s})<3
%         continue
        HMJ(s) = 0;
    end

    % 3 - Reject all shots with negative values fits: =========================
    % In some cases, the above two filters still let some abnormal results to
    % leak through. We have observed that in some cases we get negative values
    % for Ne and Te. This is clearly unphysical hence we must remove these shots
    % if ~isempty(rng)
    if mean(Ni{s}(rngMean{s}))<0 || mean(Te{s}(rngMean{s}))<0
        NegValueShots(s) = 1;
%         continue
        HMJ(s) = 0;
    end  
    
    % #####################################################################
    % EXCEPTIONS
    Exceptions_2016_12_19_to_2017_07_20
    % #####################################################################

  
     % mean Ne:
     Nea(s) = mean(Ni{s}(rngMean{s}));
     if isnan(Nea(s))
     Nea(s) = 0;
     Tea(s) = 0;
     continue
     end
     % Error bar on Ne:
      % uncertainty due to asymmetry:
     dNea_asym(s) = mean(mean(ni{s}{1}(rngMean{s})) - mean(ni{s}{2}(rngMean{s})));
      % uncertainty due to time variation:
     dNea_var(s) = std(Ni{s}(rngMean{s}),1);
     if dNea_asym(s)>=dNea_var(s)
         dNea(s) = dNea_asym(s);
     else
         dNea(s) = dNea_var(s);
     end

     Tea(s) = mean(Te{s}(rngMean{s}));
     dTea(s) = std(Te{s}(rngMean{s}),1);
     
end
%% Add offsets to radial positions of data
% From the data taken on the XPs on novemeber 2017, it became evident that
% the probes had a position offset relative to the scale that was provided.

% DLP 4.5  vertical  : -4.0  cm, based on 2016_11_17, 2017_05_02
% DLP 4.5  horizontal: +1.0  cm, based on 2017_04_21,
% DLP 6.5  vertical  : +2.5  cm, based on 2016_11_16
% DLP 9.5  vertical  : +0.5  cm, based on 2016_11_16
% DLP 10.5 vertical  : +0.5  cm, based on 2016_11_23

rng4V = find(Date == 20170502 & Spool == 4.5);
R_cm(rng4V) = R_cm(rng4V) - 4; 

% 2017_08_03: The following shotlist has not been included in the current Database
% because we forgot to include the Y_cm column in the logbook
% rng4H = find(Date == 20170421 & Spool == 4.5);
% R_cm(rng4H) = R_cm(rng4H) - 4; 

rng6V = find(Spool == 6.5 & Date>20161116 & Date<20170214);
R_cm(rng6V) = R_cm(rng6V) + 2.5; 
%% Account for DLP 6.5 damage on 2017_04_11
% During an ECH experiment with a radial scan of DLP 6.5, the probe tip was
% damaged while sampling the upper edge of the plasna column where the ECH
% is injected. The probe shattered on Shot 13847 and clearly observed in
% the visible camera videos. Comparing shots 13846 to 13851, the damage
% seems to have caused the area of the probe to approximately double.
% Hence to correct for the increased area we need to mulitply all shots
% with the damaged probe by 0.53 to get the correct plasma density.

rng6Damage = find(Spool == 6.5 & Shot >13847);
Nea(rng6Damage) = Nea(rng6Damage)*0.53;

%% Add TS data: From Nischal
TSD = importdata('TSinfo_2017_06_to_2017_08_v1.xlsx','',1);

[r,c] = find(strcmp(TSD.textdata.Sheet1,'Spool'));
Spool_TS = TSD.data.Sheet1(:,c);
[r,c] = find(strcmp(TSD.textdata.Sheet1,'r [cm]'));
R_TS = TSD.data.Sheet1(:,c);
[r,c] = find(strcmp(TSD.textdata.Sheet1,'Te [eV]'));
Te_TS = TSD.data.Sheet1(:,c);
[r,c] = find(strcmp(TSD.textdata.Sheet1,'ne [m^-3]'));
Ne_TS = TSD.data.Sheet1(:,c);
[r,c] = find(strcmp(TSD.textdata.Sheet1,'ECH '));
ECH_FP_TS = TSD.data.Sheet1(:,c);

% Add TS data from Juergen 2017_08_08:
Ne_TS_PressureScan = [1  ,1.3 ,1.5 ]*1e19; % High to low pressure
Te_TS_PressureScan = [6.9,10.7,19.2];
Spool_TS_PressureScan = 6.5*ones(size(Ne_TS_PressureScan));
R_TS_PressureScan = 1*ones(size(Ne_TS_PressureScan)); 
ECH_FP_TS_PressureScan = 1*ones(size(Ne_TS_PressureScan));

% Add TS data for PMAC 2018, 
% Code added 2018_01_24
TSD_2018_01_24 = importdata('TSinfo_DataFor_PMAC2018.xlsx','',1);
[r,c] = find(strcmp(TSD_2018_01_24.textdata.Sheet1,'Spool'));
Spool_TS_2018_01_24 = TSD_2018_01_24.data.Sheet1(:,c);
[r,c] = find(strcmp(TSD_2018_01_24.textdata.Sheet1,'r [cm]'));
R_TS_2018_01_24 = TSD_2018_01_24.data.Sheet1(:,c);
[r,c] = find(strcmp(TSD_2018_01_24.textdata.Sheet1,'Te [eV]'));
Te_TS_2018_01_24 = TSD_2018_01_24.data.Sheet1(:,c);
[r,c] = find(strcmp(TSD_2018_01_24.textdata.Sheet1,'ne [m^-3]'));
Ne_TS_2018_01_24 = TSD_2018_01_24.data.Sheet1(:,c);
[r,c] = find(strcmp(TSD_2018_01_24.textdata.Sheet1,'ECH '));
ECH_FP_TS_2018_01_24 = TSD_2018_01_24.data.Sheet1(:,c);

% Data send by Nischal on 2018_01_25, code added 2018_01_25
TSD_2018_01_16 = importdata('TSinfo_2018_01_16.xlsx','',1);
[r,c] = find(strcmp(TSD_2018_01_16.textdata.Sheet1,'Spool'));
Spool_TS_2018_01_16 = TSD_2018_01_16.data.Sheet1(:,c);
[r,c] = find(strcmp(TSD_2018_01_16.textdata.Sheet1,'r [cm]'));
R_TS_2018_01_16 = TSD_2018_01_16.data.Sheet1(:,c);
[r,c] = find(strcmp(TSD_2018_01_16.textdata.Sheet1,'Te [eV]'));
Te_TS_2018_01_16 = TSD_2018_01_16.data.Sheet1(:,c);
[r,c] = find(strcmp(TSD_2018_01_16.textdata.Sheet1,'ne [m^-3]'));
Ne_TS_2018_01_16= TSD_2018_01_16.data.Sheet1(:,c);
[r,c] = find(strcmp(TSD_2018_01_16.textdata.Sheet1,'ECH '));
ECH_FP_TS_2018_01_16 = TSD_2018_01_16.data.Sheet1(:,c);

% Append TS data:
Ne_TS =     [Ne_TS     ;Ne_TS_PressureScan'     ;Ne_TS_2018_01_24    ;Ne_TS_2018_01_16    ];
Te_TS =     [Te_TS     ;Te_TS_PressureScan'     ;Te_TS_2018_01_24    ;Te_TS_2018_01_16    ];
Spool_TS =  [Spool_TS  ;Spool_TS_PressureScan'  ;Spool_TS_2018_01_24 ;Spool_TS_2018_01_16 ];
R_TS =      [R_TS      ;R_TS_PressureScan'      ;R_TS_2018_01_24     ;R_TS_2018_01_16     ];
ECH_FP_TS = [ECH_FP_TS ;ECH_FP_TS_PressureScan' ;ECH_FP_TS_2018_01_24;ECH_FP_TS_2018_01_16];

%% Plot #0

srng = (s_start:s_end);
[UniqueDates,nU,nT] = unique(Date(srng)); % UniqueDates => Date(nU,:)  
Ndays = length(nU); 
% Date(find(nT==nd),:)

% Filters:
% On-axis shots:
b = find(R_cm>=-1 & R_cm<=1 & HMJ == 1);
% Shots with no ICH:
c = find(ICH_FP == 0);
% DLP 4.5 V:
% dlp{4} = find(Spool == 4.5 & AxisType == 'V');
dlp{4} = find(Spool == 4.5);
Cdlp{4} = 'ro';
% DLP 10.5 H:
% dlp{10} = find(Spool == 10.5 & AxisType == 'V');
dlp{10} = find(Spool == 10.5);
Cdlp{10} = 'go';
% DLP 6.5 V:
% dlp{6} = find(Spool == 6.5 & AxisType == 'V');
dlp{6} = find(Spool == 6.5);
Cdlp{6} = 'co';
% DLP 9.5 V:
% dlp{9} = find(Spool == 9.5 & AxisType == 'V');
dlp{9} = find(Spool == 9.5);
Cdlp{9} = 'blo';

figure; 
hold on
for np = [4,6,9,10]% for all probe locations
    for nd = 1:Ndays
        % Find all shots that belong to the date Date(nU(nd),:)
        a = find(nT==nd) + s_start - 1;

        % Shots in "a" that satisfy "b" (on-axis
        p = intersect(a,b);
        % Shots that do not have ICH
        q = intersect(p,c);
        % Shots whose ne and Te are no identical to zero
        u = intersect(q,find(Nea ~= 0 | Tea ~= 0));
        % Shots with ECH
        v = intersect(u,find(ECH_FP_XP == 1));
        m{np} = intersect(v,dlp{np});

        if ~isempty(m{np})
            h(nd) = plot3(Nea(m{np}),Tea(m{np}),Shot(m{np}),Cdlp{np}); 
            DateString = datetime(Date(nU(nd)+s_start-1,:),'ConvertFrom','yyyymmdd');
            DateString.Format = 'uuuu-MM-dd';
            set(h(nd),'DisplayName',[char(DateString),', dlp ',num2str(Spool(m{np}(1)))],'MarkerFaceColor',Cdlp{np}(1),...
                'MarkerSize',4)

        else
            continue
        end
    end
end

ylim([0,10]); ylabel('[eV]')
xlim([1e18,1e20]); xlabel('[m^{-3}]')

set(gcf,'color','w')
Scale = 'lin';
set(gca,'Xscale',Scale)
box on
plotbrowser('on')

try 
    clear Shot2Preview
catch
end

%% Plot #1 for Juergen 2017_08_04
% 1.	Te vs ne from DLP and TS, labelled as such:
% DLP open symbols
% TS full symbols
% helicon only circle symbol
% helicon plus ECH triangle symbol
% Probe A (spool 4.5) blue
% Probe B (spool 6.5) red
% Probe C (spool 9.5) green
% Probe D (spool 10.5) yellow
% TS (spool 6.5) red; density
% range 0.5 x 10^19 m-3 to 20 x 10^19 m-3 or logarithmic

DebugMode = 0;

% Filters:
% Shots with no ICH:
c = find(ICH_FP == 0);

% DLP 4.5 V:
Cdlp{4} = {'blo','bl^'};
NameDLP{4} = 'Probe A';
R_max(4) = 3;
% DLP 6.5 V:
Cdlp{6} = {'ro','r^'};
NameDLP{6} = 'Probe B';
R_max(6) = 1.5;
% DLP 9.5 V:
Cdlp{9} = {'go','g^'};
NameDLP{9} = 'Probe C';
R_max(9) = 1;
% DLP 10.5 H:
Cdlp{10} = {'co','c^'};
NameDLP{10} = 'Probe D';
R_max(10) = 1;


PlotECHpulses = [0,1];
Lw = [0.5,1];
Ms = [4,6];

figure; 
hold on
for ne = 1:2
for np = [4,6,9,10]% for all probe locations
    for nd = 1:Ndays
        % Find all shots that belong to the date Date(nU(nd),:)
        a = find(nT==nd) + s_start - 1;
        % On-axis shots:
        b = find(abs(R_cm)<= R_max(np) & HMJ == 1 & Nea' > 0.5e19);
        % Shots in "a" that satisfy "b" (on-axis
        p = intersect(a,b);
        % Shots that do not have ICH
        q = intersect(p,c);
        % Shots whose ne and Te are no identical to zero
        u = intersect(q,find(Nea > 0| Tea ~= 0));
        % Shots with ECH OFF
        v = intersect(u,find(ECH_FP_XP == PlotECHpulses(ne)));
        m{np} = intersect(v,dlp{np});

        if ~isempty(m{np})
            hn{np}{ne}(nd) = plot3(Nea(m{np}),Tea(m{np}),Shot(m{np}),Cdlp{np}{ne}); 
            DateString = datetime(Date(nU(nd)+s_start-1,:),'ConvertFrom','yyyymmdd');
            DateString.Format = 'uuuu-MM-dd';
            if isempty(DateString)
                np
            end
            set(hn{np}{ne}(nd),'DisplayName',[char(DateString),', ',num2str(NameDLP{np})],'MarkerEdgeColor',Cdlp{np}{ne}(1),...
                'MarkerSize',Ms(ne),'LineWidth',Lw(ne))

        else
            continue
        end
    end
end
end
 
ylim([0,25]); ylabel('$T_e$ $[eV]$','Interpreter','latex','FontSize',14)
xlim([6e18,1e20]); xlabel('$n_e$ $[m^{-3}]$','Interpreter','latex','FontSize',14)

set(gcf,'color','w')
Scale = 'log';
set(gca,'Xscale',Scale)
box on
plotbrowser('on')

% #####################################################################
% ADD TS data:

R_TS_Lim = 3;
hold on
% 1 - Spool 6.5 ECH OFF
rngTS{1} = Spool_TS == 6.5 & ECH_FP_TS == 0 & abs(R_TS) <=R_TS_Lim;
hTS(1) = plot3(Ne_TS(rngTS{1}),Te_TS(rngTS{1}),ones(size(Ne_TS(rngTS{1}))));
set(hTS(1),'Marker','o','MarkerFaceColor','r','MarkerEdgeColor','r','LineStyle','none')

% 2 - Spool 11.5  ECH OFF
rngTS{2} = Spool_TS == 11.5 & ECH_FP_TS == 0 & abs(R_TS) <=R_TS_Lim;
hTS(2) = plot3(Ne_TS(rngTS{2}),Te_TS(rngTS{2}),ones(size(Ne_TS(rngTS{2}))));
set(hTS(2),'Marker','o','MarkerFaceColor','c','MarkerEdgeColor','k','LineStyle','none')

% 3 - Spool 6.5 ECH ON
rngTS{3} = Spool_TS == 6.5 & ECH_FP_TS == 1 & abs(R_TS) <=R_TS_Lim;
hTS(3) = plot3(Ne_TS(rngTS{3}),Te_TS(rngTS{3}),ones(size(Ne_TS(rngTS{3}))));
set(hTS(3),'Marker','^','MarkerFaceColor','r','MarkerEdgeColor','r','LineStyle','none')

% 4 - Spool 11.5  ECH ON
rngTS{4} = Spool_TS == 11.5 & ECH_FP_TS == 1 & abs(R_TS) <=R_TS_Lim;
hTS(4) = plot3(Ne_TS(rngTS{4}),Te_TS(rngTS{4}),ones(size(Ne_TS(rngTS{4}))));
set(hTS(4),'Marker','^','MarkerFaceColor','c','MarkerEdgeColor','k','LineStyle','none')

legend([hn{4}{1}(end),hn{6}{1}(end),hn{9}{1}(end),hn{10}{1}(end),...
    hn{6}{2}(end),hn{9}{2}(end),hn{10}{2}(end),hTS(1),hTS(2),hTS(3),hTS(4)],...
    'Probe A ','Probe B','Probe C','Probe D','Probe B + ECH','Probe C + ECH','Probe D + ECH',...
    'TS Central Chamber','TS Target','TS Central Chamber + ECH','TS Target + ECH')


try 
    clear Shot2Preview
catch
end

if DebugMode
    ylim([0,10]); ylabel('[eV]')
    xlim([0,8e19]); xlabel('[m^{-3}]')
    Scale = 'lin';
    set(gca,'Xscale',Scale)
end
%% Plot #2 for Juergen 2017_08_09
% 2.	Te vs ne from DLP and TS, labelled as such:
% DLP open symbols
% helicon only circle symbol
% helicon plus ECH triangle symbol, 
% Probe C (spool 9.5) green
% Probe D (spool 10.5) yellow
% TS (target) black;
% density range 0.5 x 10^19 m-3 to 20 x 10^19 m-3 or logarithmic

DebugMode = 0;

% Filters:
% Shots with no ICH:
c = find(ICH_FP == 0);

% DLP 4.5 V:
Cdlp{4} = {'blo','bl^'};
NameDLP{4} = 'Probe A';
R_max(4) = 3;
% DLP 6.5 V:
Cdlp{6} = {'ro','r^'};
NameDLP{6} = 'Probe B';
R_max(6) = 1.5;
% DLP 9.5 V:
Cdlp{9} = {'go','g^'};
NameDLP{9} = 'Probe C';
R_max(9) = 1;
% DLP 10.5 H:
Cdlp{10} = {'co','c^'};
NameDLP{10} = 'Probe D';
R_max(10) = 1;


PlotECHpulses = [0,1];
Lw = [0.5,1];
Ms = [4,6];

figure; 
hold on
for ne = 1:2
for np = [9,10]% for all probe locations
    for nd = 1:Ndays
        % Find all shots that belong to the date Date(nU(nd),:)
        a = find(nT==nd) + s_start - 1;
        % On-axis shots:
        b = find(abs(R_cm)<= R_max(np) & HMJ == 1 & Nea' > 0.5e19);
        % Shots in "a" that satisfy "b" (on-axis
        p = intersect(a,b);
        % Shots that do not have ICH
        q = intersect(p,c);
        % Shots whose ne and Te are no identical to zero
        u = intersect(q,find(Nea > 0| Tea ~= 0));
        % Shots with ECH OFF
        v = intersect(u,find(ECH_FP_XP == PlotECHpulses(ne)));
        m{np} = intersect(v,dlp{np});

        if ~isempty(m{np})
            hn{np}{ne}(nd) = plot3(Nea(m{np}),Tea(m{np}),Shot(m{np}),Cdlp{np}{ne}); 
            DateString = datetime(Date(nU(nd)+s_start-1,:),'ConvertFrom','yyyymmdd');
            DateString.Format = 'uuuu-MM-dd';
            if isempty(DateString)
                np
            end
            set(hn{np}{ne}(nd),'DisplayName',[char(DateString),', ',num2str(NameDLP{np})],'MarkerEdgeColor',Cdlp{np}{ne}(1),...
                'MarkerSize',Ms(ne),'LineWidth',Lw(ne))

        else
            continue
        end
    end
end
end
 
ylim([0,7]); ylabel('$T_e$ $[eV]$','Interpreter','latex','FontSize',14)
xlim([5e18,1e20]); xlabel('$n_e$ $[m^{-3}]$','Interpreter','latex','FontSize',14)

set(gcf,'color','w')
Scale = 'log';
set(gca,'Xscale',Scale)
box on
plotbrowser('on')

% #####################################################################
% ADD TS data:

R_TS_Lim = 3;
hold on
% 1 - Spool 6.5 ECH OFF
% rngTS{1} = Spool_TS == 6.5 & ECH_FP_TS == 0 & abs(R_TS) <=R_TS_Lim;
% hTS(1) = plot3(Ne_TS(rngTS{1}),Te_TS(rngTS{1}),ones(size(Ne_TS(rngTS{1}))));
% set(hTS(1),'Marker','o','MarkerFaceColor','r','MarkerEdgeColor','r','LineStyle','none')

% 2 - Spool 11.5  ECH OFF
rngTS{2} = Spool_TS == 11.5 & ECH_FP_TS == 0 & abs(R_TS) <=R_TS_Lim;
hTS(2) = plot3(Ne_TS(rngTS{2}),Te_TS(rngTS{2}),ones(size(Ne_TS(rngTS{2}))));
set(hTS(2),'Marker','o','MarkerFaceColor','k','MarkerEdgeColor','k','LineStyle','none')

% 3 - Spool 6.5 ECH ON
% rngTS{3} = Spool_TS == 6.5 & ECH_FP_TS == 1 & abs(R_TS) <=R_TS_Lim;
% hTS(3) = plot3(Ne_TS(rngTS{3}),Te_TS(rngTS{3}),ones(size(Ne_TS(rngTS{3}))));
% set(hTS(3),'Marker','^','MarkerFaceColor','r','MarkerEdgeColor','r','LineStyle','none')

% 4 - Spool 11.5  ECH ON
rngTS{4} = Spool_TS == 11.5 & ECH_FP_TS == 1 & abs(R_TS) <=R_TS_Lim;
hTS(4) = plot3(Ne_TS(rngTS{4}),Te_TS(rngTS{4}),ones(size(Ne_TS(rngTS{4}))));
set(hTS(4),'Marker','^','MarkerFaceColor','k','MarkerEdgeColor','k','LineStyle','none')

legend([hn{9}{1}(end),hn{10}{1}(end),...
    hn{9}{2}(end),hn{10}{2}(end),hTS(2),hTS(4)],...
    'Probe C','Probe D','Probe C + ECH','Probe D + ECH',...
    'TS Target','TS Target + ECH')


try 
    clear Shot2Preview
catch
end

if DebugMode
    ylim([0,7]); ylabel('[eV]')
    xlim([0,8e19]); xlabel('[m^{-3}]')
    Scale = 'lin';
    set(gca,'Xscale',Scale)
end