% Created by JF Caneses, 2017_07_16
% Ne and Te database

close all
clear all

state = 1;
% 1: load data
% 2: Compute data
% 3: Update DLPinfo but do NOT compute data

switch state
    case 1 % Load data
        % Load previously calculated Ne and Te
        load('NeTe_2016_07_28_to_2016_11_15.mat')
        % Hence no need to reload DLPinfo, recompute the data or save
        % anything:
        LoadDLPinfo = 0;
        CMPT = 0;
        SaveAllData = 0; 
        SaveDLPinfo = 0;
    case 2 % Compute data
        % 299 shots in 526 seconds; ~ 1.75 sec per Shot
        % 299 shots in 8:46 min ;
        
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
end

if LoadDLPinfo
    D = importdata('DLPinfo_2016_07_28_to_2016_11_15.xlsx','',1);
        
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

if 0 % This part of the code is done to systematically check the quality of the data
     % Once the entire data set has been checked, we set it to 0 and let
     % the entire sequence compute and then to be saved
    s_start = find(Shot == 9515);
    s_end   = find(Shot == 9535);
else
    s_start = 1;
    s_end   = length(Shot);
end

for s = s_start(1):s_end(end) % For all shots
    if isnan(Shot(s))
        continue
    end
    DataAddress{1} = [RootAddress,MDS_V{s}]; % V
    DataAddress{2} = [RootAddress,MDS_I{s}]; % I
    Config.I_Att  = Iatt(s);  % Output voltage of DLP box (Current) = I_att*Digitized data
    Config.tStart = Tstart(s); % [s]
    Config.tEnd   = Tend(s);
    Config.L_tip = Ltip_mm(s)/1000;
    Config.D_tip = Dtip_mm(s)/1000; % [m]
    Config.V_Att = Vatt(s);  % Output voltage of DLP box (Voltage) = V_Att*Digitized data 
    Config.V_cal = [Vcal_1(s),Vcal_2(s)];   % Voltage output of DLP = V_cal(1)*Output voltage of DLP box + V_cal(2)
    Config.I_cal = [Ical_1(s),Ical_2(s)]; % Current output of DLP = (I_cal(2) + Output voltage of DLP box)/Ical(1)

    Config.Center_I = 1;
    Config.Center_V = 0;
    [a,b,c,d,f,g,j,p,q,x,y] = DLP_fit_V5_2(Config,Shot(s),DataAddress);
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
    GlitchFlag{s} = x{1};
    SSQres{s} = y{1};
end
end

% Save data:
if SaveAllData
   save('NeTe_2016_07_28_to_2016_11_15')
elseif SaveDLPinfo
    load('NeTe_2016_07_28_to_2016_11_15','ni','Ni','Te','time','Ifit','Ip',...
        'Vp','tm','Vsweep','Isweep','GlitchFlag','SSQres');
    save('NeTe_2016_07_28_to_2016_11_15')
end

%% Extract Steady state data
try 
    s_start;
catch
    s_start = 1;
    s_end = length(Shot);
end
for s = s_start:s_end 
    
    if isnan(Shot(s))
        continue
    end
    
 rng = find(time{s}>=Gstart(s) & time{s}<=Gend(s));
 % mean Ne:
 Nea(s) = mean(Ni{s}(rng));
 
 % Error bar on Ne:
  % uncertainty due to asymmetry:
 dNea_asym(s) = mean(mean(ni{s}{1}(rng)) - mean(ni{s}{2}(rng)));
  % uncertainty due to time variation:
 dNea_var(s) = std(Ni{s}(rng),1);
 if dNea_asym(s)>=dNea_var(s)
     dNea(s) = dNea_asym(s);
 else
     dNea(s) = dNea_var(s);
 end
 
 Tea(s) = mean(Te{s}(rng));
 dTea(s) = std(Te{s}(rng),1);
end

%% Plot data
[UniqueDates,nU,nT] = unique(Date); % UniqueDates => Date(nU,:)  
Ndays = length(nU); 
% Date(find(nT==nd),:)

% Filters:
% On-axis shots:
b = find(R_cm>=-3 & R_cm<=3);
% Shots with no ICH:
c = find(ICH_FP == 0);
% DLP 4.5 V:
dlp{4} = find(Spool == 4.5 & AxisType == 'V');
Cdlp{4} = 'ro';
% DLP 10.5 H:
dlp{10} = find(Spool == 10.5 & AxisType == 'V');
Cdlp{10} = 'go';
% DLP 6.5 V:
dlp{6} = find(Spool == 6.5 & AxisType == 'V');
Cdlp{6} = 'co';
% DLP 9.5 V:
dlp{9} = find(Spool == 9.5 & AxisType == 'V');
Cdlp{9} = 'blo';

figure; 
hold on
for np = [4,6,9,10]% for all probe locations
    for nd = 1:Ndays
        % Find all shots that belong to the date Date(nU(nd),:)
        a = find(nT==nd);

        % Shots in "a" that satisfy "b" (on-axis
        p = intersect(a,b);
        % Shots that do not have ICH
        q = intersect(p,c);
        m{np} = intersect(q,dlp{np});

        if ~isempty(m{np})
            h(nd) = plot3(Nea(m{np}),Tea(m{np}),Shot(m{np}),Cdlp{np}); 
            DateString = datetime(Date(nU(nd),:),'ConvertFrom','yyyymmdd');
            DateString.Format = 'uuuu-MM-dd';
            set(h(nd),'DisplayName',[char(DateString),', dlp ',num2str(Spool(m{np}(1)))],'MarkerFaceColor',Cdlp{np}(1))

        else
            continue
        end
    end
end

ylim([0,10]); ylabel('[eV]')
xlim([1e18,1e20]); xlabel('[m^{-3}]')

set(gcf,'color','w')
set(gca,'Xscale','log')
box on
plotbrowser('on')

try 
    clear Shot2Preview
catch
end

