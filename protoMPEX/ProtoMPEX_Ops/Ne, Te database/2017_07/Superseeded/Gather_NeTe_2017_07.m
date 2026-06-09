% Created by JF Caneses, 2017_07_16
% Ne and Te database

close all
clear all

state = 2;
% 1: load data
% 2: Compute data

switch state
    case 1
        load('NeTe_2016_07_28_to_2016_11_15.mat')
        CMPT = 0;
    case 2
        % XX shots in XX seconds; ~ X.X sec per shot
        % XX shots in X:XX min  (Xs);

        
        CMPT = 1; % Compute data
        SaveData = 1;
        D = importdata('DLPinfo_2016_07_28_to_2016_11_15.xlsx','',1);
        
        % Year      : 1
        % Month     : 2
        % Day       : 3
        % shot      : 4
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
        
%         Year     = D.data(:,find(strcmp(D.colheaders,'HMJ'  )==1));
%         Month    = D.data(:,find(strcmp(D.colheaders,'Month')==1));
%         Day      = D.data(:,find(strcmp(D.colheaders,'Day')==1));
%         shot = D.data(:,find(strcmp(D.colheaders,'shot')==1));
%         Spool    = D.data(:,find(strcmp(D.colheaders,'Spool')==1));
%         V_H      = D.data(:,find(strcmp(D.colheaders,'V_H')==1));
%         R_cm        = D.data(:,find(strcmp(D.colheaders,'R_cm')==1));
%         Gauge    = D.data(:,find(strcmp(D.colheaders,'Gauge')==1));
%         HMJ      = D.data(:,find(strcmp(D.colheaders,'HMJ')==1));
%         Tstart   = D.data(:,find(strcmp(D.colheaders,'Tstart')==1));
%         Tend     = D.data(:,find(strcmp(D.colheaders,'Tend')==1));
%         Gstart   = D.data(:,find(strcmp(D.colheaders,'Gstart')==1));
%         Gend     = D.data(:,find(strcmp(D.colheaders,'Gend')==1));
%         Vcal_1   = D.data(:,find(strcmp(D.colheaders,'Vcal_1')==1));
%         Vcal_2   = D.data(:,find(strcmp(D.colheaders,'Vcal_2')==1));
%         Vatt     = D.data(:,find(strcmp(D.colheaders,'Vatt')==1));
%         Ical_1   = D.data(:,find(strcmp(D.colheaders,'Ical_1')==1));
%         Ical_2   = D.data(:,find(strcmp(D.colheaders,'Ical_2')==1));
%         Iatt     = D.data(:,find(strcmp(D.colheaders,'Iatt')==1));
%         Ltip     = D.data(:,find(strcmp(D.colheaders,'Ltip_mm')==1));
%         Dtip     = D.data(:,find(strcmp(D.colheaders,'Dtip_mm')==1));
%         MDS_V    = D.textdata(2:end,find(strcmp(D.colheaders,'MDS_V')==1));
%         MDS_I    = D.textdata(2:end,find(strcmp(D.colheaders,'MDS_I')==1));
%         ICH_FP    = D.data(:,find(strcmp(D.colheaders,'ICP_FP')==1));
%         ECH_FP    = D.data(:,find(strcmp(D.colheaders,'ECH_FP')==1)); % 2017_07_16, need to check the ECH data
for s = 1:length(D.colheaders)
    assignin('base',D.colheaders{s},D.data(:,s));
end
MDS_V    = D.textdata(2:end,find(strcmp(D.colheaders,'MDS_V')==1));
MDS_I    = D.textdata(2:end,find(strcmp(D.colheaders,'MDS_I')==1));
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
    s_start = find(shot == 9515);
    s_end   = find(shot == 9535);
else
    s_start = 1;
    s_end   = length(shot);
end

for s = s_start(1):s_end(end) % For all shots
    if isnan(shot(s))
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

%     [a,b,c,d,f,g,j,p,q] = DLP_fit_V5(Config,shot(s),DataAddress);
    Config.Center_I = 1;
    Config.Center_V = 0;
    [a,b,c,d,f,g,j,p,q,x,y] = DLP_fit_V5_2(Config,shot(s),DataAddress);
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

if SaveData
   save('NeTe_2016_07_28_to_2016_11_15')
end

end

%% Extract Steady state data

for s = s_start:s_end 
    
    if isnan(shot(s))
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
% NeTe = [Month,Day,Year,shot,,,,R_cm,Spool,,,,Nea',Tea'];
fltr{1}  = find(ICH_FP == 0 & Spool==9.5);
fltr{2} = find(R_cm>=-1 & R_cm<=1);
m9 = intersect(fltr{2},fltr{1});

DD = [Month,Day,Year,shot,R_cm,Spool,Nea',Tea'];
NeTe9 = DD(m9,:);

fltr{1}  = find(ICH_FP == 0 & Spool==6.5);
fltr{2} = find(R_cm>=-1.5 & R_cm<=1.5);
m6 = intersect(fltr{2},fltr{1});
NeTe6 = DD(m6,:);

fltr{1}  = find(ICH_FP == 0 & Spool==10.5);
fltr{2} = find(R_cm>=-1 & R_cm<=1);
m10 = intersect(fltr{2},fltr{1});

NeTe10 = DD(m10,:);

f3d = figure; hold on
h(1) = plot3(NeTe9(:,7),NeTe9(:,8),NeTe9(:,4),'ko'); 
set(h(1),'DisplayName',[num2str(DD(1,3)),', 9.5'])

h(2) = plot3(NeTe10(:,7),NeTe10(:,8),NeTe10(:,4),'ro');
set(h(2),'DisplayName',[num2str(DD(1,3)),', 10.5'])

h(3) = plot3(NeTe6(:,7),NeTe6(:,8),NeTe6(:,4),'go');
set(h(3),'DisplayName',[num2str(DD(1,3)),', 6.5'])

ylim([0,10]); ylabel('[eV]')
xlim([1e18,1e20]); xlabel('[m^{-3}]')

plotbrowser('on')
set(gcf,'color','w')
set(gca,'Xscale','log')
box on