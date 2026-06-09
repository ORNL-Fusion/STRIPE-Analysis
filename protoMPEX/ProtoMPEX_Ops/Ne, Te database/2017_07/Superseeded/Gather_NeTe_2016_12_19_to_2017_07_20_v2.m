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
        continue
    end

    % 3 - Reject all shots with negative values fits: =========================
    % In some cases, the above two filters still let some abnormal results to
    % leak through. We have observed that in some cases we get negative values
    % for Ne and Te. This is clearly unphysical hence we must remove these shots
    % if ~isempty(rng)
    if mean(Ni{s}(rngMean{s}))<0 || mean(Te{s}(rngMean{s}))<0
        NegValueShots(s) = 1;
        continue
    end  
end
    
for s = s_start:s_end 
    % #####################################################################
    % EXCEPTIONS
    Exceptions_2016_12_19_to_2017_07_20
    % #####################################################################

  
     % mean Ne:
     Nea(s) = mean(Ni{s}(rngMean{s}));

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

%% Plot data

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
        m{np} = intersect(u,dlp{np});

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

%% Comments on data
% DLP 4.5 with NISO
% 12414, 12416, 12415