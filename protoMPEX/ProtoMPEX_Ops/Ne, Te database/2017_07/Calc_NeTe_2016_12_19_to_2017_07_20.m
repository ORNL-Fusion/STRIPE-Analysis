% Created by JF Caneses, 2017_07_28
% Ne and Te database

close all
clear all

% 1 - need to recompute data: DataCorruptedShots did not get saved and
% without it we can not compute the data

state = 1;
% 1: load data
% 2: Compute data
% 3: Update DLPinfo but do NOT compute data
% 4: Compute sections of the data

switch state
    case 1 % Load data, 1531 shots in 15 seconds
        % Load previously calculated Ne and Te
        load('NeTe_2016_12_19_to_2017_07_20.mat')
        % Hence no need to reload DLPinfo, recompute the data or save
        % anything:
        LoadDLPinfo = 0;
        CMPT = 0;
        SaveAllData = 0; 
        SaveDLPinfo = 0;
    case 2 % Compute data
        % 1503 shots in 2253.710 seconds; ~ 1.49 sec per Shot, ethernet
        % 1503 shots in 38 min ;
        
        % Load DLPinfo
        LoadDLPinfo = 1;
        % Compute data
        CMPT = 1;
        % Save computed data and DLPinfo data
        SaveAllData = 1;
    case 3 % DO NOT compute data but just update DLPinfo
        LoadDLPinfo = 1;
        CMPT = 0;
        SaveAllData = 0; 
        SaveDLPinfo = 1;
    case 4 % Compute sections of the data
        % Load DLPinfo
        LoadDLPinfo = 1;
        % Compute data
        CMPT = 1;
        % Do not save computed data and DLPinfo data
        SaveAllData = 0;    
        SaveDLPinfo = 0;  
    case 5
        LoadDLPinfo = 0;
        CMPT = 0;
        SaveAllData = 0; 
        SaveDLPinfo = 0;
end

if LoadDLPinfo
    D = importdata('DLPinfo_2016_12_19_to_2017_07_20.xlsx','',1);
        
        % Year      : 1
        % Month     : 2
        % Day       : 3
        % Shot      : 4
        % spool     : 5
        % V_H       : 6
        % R_cm      : 7
        % Gauge     : 8
        % HMJ       : 9
        % Tstart    : 10
        % Tend      : 11
        % Gstart    : 12
        % Gend      : 13
        % Vcal_1    : 14
        % Vcal_2    : 15
        % Vatt      : 16
        % Ical_1    : 17
        % Ical_2    : 18
        % Iatt      : 19
        % Ltip_mm   : 20
        % Dtip_mm   : 21
        % MDS_V     : 22
        % MDS_I     : 23
        % ICH_FP    : 24
        % ECH_FP    : 25
        
        for s = 1:length(D.colheaders)
            f = D.data(:,s);
            if ~isnan(f)
                assignin('base',D.colheaders{s},f);
            else
                f = char(D.textdata(2:end,s));
                assignin('base',D.colheaders{s},f);
            end
        end
        for s = 1:length(Shot)
                Date(s,:) = yyyymmdd(datetime(Year(s),Month(s),Day(s)));
        end
end


if CMPT == 1
% Address:
Stem = '\MPEX::TOP.';
Branch = 'MACHOPS1:';
RootAddress = [Stem,Branch];
Config.FilterDataInput = 1; % Filter input data with savitsky Golay filter order 3 
Config.SGF = 5; % frame size for the SG filter
Config.TimeMode = 2; % Effective time of sweep, (1) start of ramp or (2) mean time of ramp
Config.AMU = 2; % Ion mass in AMU
Config.FitFunction = 2; 
Config.AreaType = 1; % 1: Cylindrical + cap

% =========================================================================
if state == 4 % This part of the code is done to systematically check the quality of the data
     % Once the entire data set has been checked, we set it to 0 and let
     % the entire sequence compute and then to be saved
%     s_start = find(Shot == 14148);
%     s_end   = find(Shot == 14157);
%     s_start = 1300;
%     s_end   = 1531;
    s_start = find(Shot == 15673);
    s_end   = find(Shot == 15676);
elseif state == 2
    s_start = 1;
    s_end   = length(Shot);
end
% =========================================================================

for s = s_start(1):s_end(end) % For all shots
    
    if s == 228
        s
    end
    if Shot(s) == 12230
        Shot(s)
    end
    
    Shot(s)
    
    DataAddress{1} = [RootAddress,MDS_V(s,:)]; % V
    DataAddress{2} = [RootAddress,MDS_I(s,:)]; % I
    Config.I_Att  = Iatt(s);  % Output voltage of DLP box (Current) = I_att*Digitized data
    Config.tStart = Tstart(s); % [s]
    Config.tEnd   = Tend(s);
    Config.L_tip = Ltip_mm(s)/1000;
    Config.D_tip = Dtip_mm(s)/1000; % [m]
    Config.V_Att = Vatt(s);  % Output voltage of DLP box (Voltage) = V_Att*Digitized data 
    Config.V_cal = [Vcal_1(s),Vcal_2(s)];   % Voltage output of DLP = V_cal(1)*Output voltage of DLP box + V_cal(2)
    Config.I_cal = [Ical_1(s),Ical_2(s)]; % Current output of DLP = (I_cal(2) + Output voltage of DLP box)/Ical(1)

    Config.Center_I = 1;
    
    if Spool(s) == 4.5
            Config.Center_V = 1;
    else
            Config.Center_V = 1;
    end
%     tic
    try
        [a,b,c,d,f,g,j,p,q,w,x,y,z] = DLP_fit_V5_3(Config,Shot(s),DataAddress);
        s
        DataCorruptedShots(s) = 0; 
    catch 
        DataCorruptedShots(s) = 1; 
        continue
    end
%     toc
    ni{s}{1} = a{1}{1}; ni{s}{2} = a{1}{2};
    Ni{s} = 0.5*(ni{s}{1} + ni{s}{2});
    Te{s} = b{1};
    time{s} = c{1};
    Ifit{s} = d{1};
    Ip{s} = f{1};
    Vp{s} = g{1};
    tm{s} = j{1};
    Vsweep{s} = p{1};
    Isweep{s} = q{1};
    GlitchFlag{s} = w{1};
    SSQres{s} = x{1};
    StdRes{s} = y{1};
    StdResNorm{s} = z{1};
end
end

% Save data:
if SaveAllData
   save('NeTe_2016_12_19_to_2017_07_20.mat')
elseif SaveDLPinfo
   load('NeTe_2016_12_19_to_2017_07_20.mat','ni','Ni','Te','time','Ifit','Ip',...
        'Vp','tm','Vsweep','Isweep','GlitchFlag','SSQres','StdRes','StdResNorm','DataCorruptedShots');
   save('NeTe_2016_12_19_to_2017_07_20.mat')
end

%%
% THE FOLLOWING HAS BEEN CREATED INTO AN INDIVIDUAL M-SCRIPT



% try 
%     s_start;
% catch
%     s_start = 1;
%     s_end = length(Shot);
% end
% 
% %% % REMOVE ALL BAD SHOTS ================================================
% for s = s_start:s_end 
% % 1- Reject all shots whose data is corrupted: ============================
% % this is, when the DLP code failed
% 
% try % DataCorruptedShots may not have been saved in the latest run
%     if DataCorruptedShots(s) || isnan(Shot(s))
%         continue
%     end
% end
% 
% % 2 - Reject all shots with bad fits: =====================================
% % DLP code did not fail but produced fits with high residuals. We use the
% % value of StdResNorm as the metric to judge the fit.
% condition = find(StdResNorm{s}<(10/100) & StdResNorm{s} > 0 & StdRes{s} == real(StdRes{s}));
% rngTime = find(time{s}>=Gstart(s) & time{s}<=Gend(s));
% rngMean{s} = intersect(rngTime,condition);
% if length(rngMean{s})<3
%     continue
% end
% 
% % 3 - Reject all shots with negative values fits: =========================
% % In some cases, the above two filters still let some abnormal results to
% % leak through. We have observed that in some cases we get negative values
% % for Ne and Te. This is clearly unphysical hence we must remove these shots
% % if ~isempty(rng)
% if mean(Ni{s}(rngMean{s}))<0 || mean(Te{s}(rngMean{s}))<0
%     NegValueShots(s) = 1;
%     continue
% end
% end
% 
% %% EXCEPTIONS =============================================================
% for s = s_start:s_end 
% % 2017_06_16: The mach probe 10.5 V was swept with the isolation transformer in
% % order to get a the difference in the ion saturation current. this was not
% % recorded on the logbook. shots: 15168, 15170, 15171
% if sum(Shot(s) == [15168, 15170, 15171])
%     rngMean{s} = [];
% end
% % 2017_07_18: In the first 6 shots of the day we had 2 shots that were not
% % in the high density mode but they were recorded as HMJ = 1. 
% if sum(Shot(s) == [15714, 15715])
%     HMJ(s) = 0;
% end
% % Other shots that were not in the high density mode
% if sum(Shot(s) == [15083,15084])
%     HMJ(s) = 0;
% end
% % This shot has unsually high mean Te due to a glitch that could not be
% % detected by the code
% if sum(Shot(s) == [14272,14365])
%     rngMean{s} = rngMean{s}(1:end-2);
% end
% % 2017_07_18: we connected the IF probe with an arbitrary off-axis target
% % DLP and swept it as a DLP. the IV traces are not very good
% if sum(Shot(s) == [15723:15725])
%     HMJ(s) = 0;
% end
% % 2017_01_12: DLP 9.5 on this day have significanlty different measurements
% % when compared to DLP 10.5. Te measured to be 6 eV while DLP 10.5 was
% % measured to be 1.5 eV
% if sum(Shot(s) == [12791,12793])
%     HMJ(s) = 0;
% end
% end
% %% STEADY STATE DATA: =====================================================
% for s = s_start:s_end 
%  % mean Ne:
%  Nea(s) = mean(Ni{s}(rngMean{s}));
%  
%  % Error bar on Ne:
%   % uncertainty due to asymmetry:
%  dNea_asym(s) = mean(mean(ni{s}{1}(rngMean{s})) - mean(ni{s}{2}(rngMean{s})));
%   % uncertainty due to time variation:
%  dNea_var(s) = std(Ni{s}(rngMean{s}),1);
%  if dNea_asym(s)>=dNea_var(s)
%      dNea(s) = dNea_asym(s);
%  else
%      dNea(s) = dNea_var(s);
%  end
%  
%  Tea(s) = mean(Te{s}(rngMean{s}));
%  dTea(s) = std(Te{s}(rngMean{s}),1);
% end
% %% Add offsets to radial positions of data
% % From the data taken on the XPs on novemeber 2017, it became evident that
% % the probes had a position offset relative to the scale that was provided.
% 
% % DLP 4.5  vertical  : -4.0  cm, based on 2016_11_17, 2017_05_02
% % DLP 4.5  horizontal: +1.0  cm, based on 2017_04_21,
% % DLP 6.5  vertical  : +2.5  cm, based on 2016_11_16
% % DLP 9.5  vertical  : +0.5  cm, based on 2016_11_16
% % DLP 10.5 vertical  : +0.5  cm, based on 2016_11_23
% 
% rng4V = find(Date == 20170502 & Spool == 4.5);
% R_cm(rng4V) = R_cm(rng4V) - 4; 
% 
% % 2017_08_03: The following shotlist has not been included in the current Database
% % because we forgot to include the Y_cm column in the logbook
% % rng4H = find(Date == 20170421 & Spool == 4.5);
% % R_cm(rng4H) = R_cm(rng4H) - 4; 
% 
% rng6V = find(Spool == 6.5 & Date>20161116 & Date<20170214);
% R_cm(rng6V) = R_cm(rng6V) + 2.5; 
% %% Account for DLP 6.5 damage on 2017_04_11
% % During an ECH experiment with a radial scan of DLP 6.5, the probe tip was
% % damaged while sampling the upper edge of the plasna column where the ECH
% % is injected. The probe shattered on Shot 13847 and clearly observed in
% % the visible camera videos. Comparing shots 13846 to 13851, the damage
% % seems to have caused the area of the probe to approximately double.
% % Hence to correct for the increased area we need to mulitply all shots
% % with the damaged probe by 0.53 to get the correct plasma density.
% 
% rng6Damage = find(Spool == 6.5 & Shot >13847);
% Nea(rng6Damage) = Nea(rng6Damage)*0.53;
% 
% %% Plot data
% 
% srng = (s_start:s_end);
% [UniqueDates,nU,nT] = unique(Date(srng)); % UniqueDates => Date(nU,:)  
% Ndays = length(nU); 
% % Date(find(nT==nd),:)
% 
% % Filters:
% % On-axis shots:
% b = find(R_cm>=-1 & R_cm<=1 & HMJ == 1);
% % Shots with no ICH:
% c = find(ICH_FP == 0);
% % DLP 4.5 V:
% % dlp{4} = find(Spool == 4.5 & AxisType == 'V');
% dlp{4} = find(Spool == 4.5);
% Cdlp{4} = 'ro';
% % DLP 10.5 H:
% % dlp{10} = find(Spool == 10.5 & AxisType == 'V');
% dlp{10} = find(Spool == 10.5);
% Cdlp{10} = 'go';
% % DLP 6.5 V:
% % dlp{6} = find(Spool == 6.5 & AxisType == 'V');
% dlp{6} = find(Spool == 6.5);
% Cdlp{6} = 'co';
% % DLP 9.5 V:
% % dlp{9} = find(Spool == 9.5 & AxisType == 'V');
% dlp{9} = find(Spool == 9.5);
% Cdlp{9} = 'blo';
% 
% figure; 
% hold on
% for np = [4,6,9,10]% for all probe locations
%     for nd = 1:Ndays
%         % Find all shots that belong to the date Date(nU(nd),:)
%         a = find(nT==nd) + s_start - 1;
% 
%         % Shots in "a" that satisfy "b" (on-axis
%         p = intersect(a,b);
%         % Shots that do not have ICH
%         q = intersect(p,c);
%         % Shots whose ne and Te are no identical to zero
%         u = intersect(q,find(Nea ~= 0 | Tea ~= 0));
%         m{np} = intersect(u,dlp{np});
% 
%         if ~isempty(m{np})
%             h(nd) = plot3(Nea(m{np}),Tea(m{np}),Shot(m{np}),Cdlp{np}); 
%             DateString = datetime(Date(nU(nd)+s_start-1,:),'ConvertFrom','yyyymmdd');
%             DateString.Format = 'uuuu-MM-dd';
%             set(h(nd),'DisplayName',[char(DateString),', dlp ',num2str(Spool(m{np}(1)))],'MarkerFaceColor',Cdlp{np}(1),...
%                 'MarkerSize',4)
% 
%         else
%             continue
%         end
%     end
% end
% 
% ylim([0,10]); ylabel('[eV]')
% xlim([1e18,1e20]); xlabel('[m^{-3}]')
% 
% set(gcf,'color','w')
% Scale = 'lin';
% set(gca,'Xscale',Scale)
% box on
% plotbrowser('on')
% 
% try 
%     clear Shot2Preview
% catch
% end
% 
% %% Comments on data
% % DLP 4.5 with NISO
% % 12414, 12416, 12415
