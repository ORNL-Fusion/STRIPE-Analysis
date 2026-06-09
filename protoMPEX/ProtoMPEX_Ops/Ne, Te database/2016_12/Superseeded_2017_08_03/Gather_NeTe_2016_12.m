% December 7th 2016
% Ne and Te database

close all
clear all

state = 1;
% 1: load data
% 2: Compute data

switch state
    case 1
        load('NeTe_data.mat')
        CMPT = 0;
    case 2
        %prompt = 'Are you sure you want to compute [Y/N]?';
        %x = input(prompt);
        
        % 95 shots in 245 seconds; ~ 2.6 sec per shot
        % 300 shots in 5:41 min  (341s);

        
        CMPT = 1; % Compute data
        SaveData = 1;
        D = importdata('DLP_info.xlsx','',1);
        
        % Day       : 1
        % Month     : 2
        % Year      : 3
        % Spool     : 4
        % shot      : 5
        % R_cm      : 6
        % Gauge     : 7
        % HMJ       : 8
        % Tstart    : 9
        % Tend      : 10
        % Gstart    : 11
        % Gend      : 12
        % Vcal_1    : 13
        % Vcal_2    : 14
        % Vatt      : 15
        % Ical_1    : 16
        % Ical_2    : 17
        % Iatt      : 18
        % Ltip_mm   : 19
        % Dtip_mm   : 20
        % MDS_V     : 21
        % MDS_I     : 22

        Day      = D.data(:,1);
        Month    = D.data(:,2);
        Year     = D.data(:,3);
        Spool    = D.data(:,4);
        shotlist = D.data(:,5);
        R        = D.data(:,6);
        Gauge    = D.data(:,7);
        HMJ      = D.data(:,8);
        Tstart   = D.data(:,9);
        Tend     = D.data(:,10);
        Gstart   = D.data(:,11);
        Gend     = D.data(:,12);
        Vcal_1   = D.data(:,13);
        Vcal_2   = D.data(:,14);
        Vatt     = D.data(:,15);
        Ical_1   = D.data(:,16);
        Ical_2   = D.data(:,17);
        Iatt     = D.data(:,18);
        Ltip     = D.data(:,19);
        Dtip     = D.data(:,20);
        MDS_V    = D.textdata(2:end,21);
        MDS_I    = D.textdata(2:end,22);
        ICH_F    = D.data(:,23);
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
    s_start = find(shotlist == 9515);
    s_end   = find(shotlist == 9535);
else
    s_start = 1;
    s_end   = length(shotlist);
end

for s = s_start(1):s_end(end) % For all shots
    if isnan(shotlist(s))
        continue
    end
    DataAddress{1} = [RootAddress,MDS_V{s}]; % V
    DataAddress{2} = [RootAddress,MDS_I{s}]; % I
    Config.I_Att  = Iatt(s);  % Output voltage of DLP box (Current) = I_att*Digitized data
    Config.tStart = Tstart(s); % [s]
    Config.tEnd   = Tend(s);
    Config.L_tip = Ltip(s)/1000;
    Config.D_tip = Dtip(s)/1000; % [m]
    Config.V_Att = Vatt(s);  % Output voltage of DLP box (Voltage) = V_Att*Digitized data 
    Config.V_cal = [Vcal_1(s),Vcal_2(s)];   % Voltage output of DLP = V_cal(1)*Output voltage of DLP box + V_cal(2)
    Config.I_cal = [Ical_1(s),Ical_2(s)]; % Current output of DLP = (I_cal(2) + Output voltage of DLP box)/Ical(1)

    [a,b,c,d,f,g,j,p,q] = DLP_fit_V5(Config,shotlist(s),DataAddress);
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
end

if SaveData
   save('NeTe_data')
end

end

%% Extract Steady state data

for s = s_start:s_end 
    
    if isnan(shotlist(s))
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
figure; hold on

fltr{1} = find(Spool == 6.5 & ICH_F == 0);
fltr{2} = find(R>=-1 & R<=1);
m1 = intersect(fltr{2},fltr{1});

fltr{1} = find(Spool == 9.5 & ICH_F == 0);
fltr{2} = find(R>=-1 & R<=1);
m2 = intersect(fltr{2},fltr{1});

fltr{1} = find(Spool == 10.5 & ICH_F == 0);
fltr{2} = find(R>=-1 & R<=1);
m3 = intersect(fltr{2},fltr{1});

plot(Nea(m1),Tea(m1),'ko')
plot(Nea(m2),Tea(m2),'ro')
plot(Nea(m3),Tea(m3),'go')

ylim([0,20])
xlim([0,10e19])

% figure; plot(Vsweep{s}{c},Ifit{s}{c})
% figure; plot(Vsweep{s}{c},Isweep{s}{c})

    if 1
        figure; hold on
        for s = 1:length(Ni)
        plot(time{s},Ni{s})
        end
    end

%%
% NeTe = [Month,Day,Year,shotlist,,,,R,Spool,,,,Nea',Tea'];
fltr{1}  = find(ICH_F == 0 & Spool==9.5);
fltr{2} = find(R>=-1 & R<=1);
m9 = intersect(fltr{2},fltr{1});

DD = [Month,Day,Year,shotlist,R,Spool,Nea',Tea'];
NeTe9 = DD(m9,:);

fltr{1}  = find(ICH_F == 0 & Spool==6.5);
fltr{2} = find(R>=-1.5 & R<=1.5);
m6 = intersect(fltr{2},fltr{1});
NeTe6 = DD(m6,:);

fltr{1}  = find(ICH_F == 0 & Spool==10.5);
fltr{2} = find(R>=-1 & R<=1);
m10 = intersect(fltr{2},fltr{1});

NeTe10 = DD(m10,:);

figure; hold on
h(1) = plot3(NeTe9(:,7),NeTe9(:,8),NeTe9(:,4),'ko'); 
set(h(1),'DisplayName',[num2str(DD(1,3)),', 9.5'])

h(2) = plot3(NeTe10(:,7),NeTe10(:,8),NeTe10(:,4),'ro');
set(h(2),'DisplayName',[num2str(DD(1,3)),', 10.5'])

h(3) = plot3(NeTe6(:,7),NeTe6(:,8),NeTe6(:,4),'go');
set(h(3),'DisplayName',[num2str(DD(1,3)),', 6.5'])

ylim([0,10])
xlim([0,1e20])

plotbrowser('on')
% plottools(gcf,'plotbrowser')