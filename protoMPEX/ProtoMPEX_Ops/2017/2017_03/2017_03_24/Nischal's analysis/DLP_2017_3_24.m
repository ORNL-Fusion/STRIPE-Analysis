
%% Day's experiment focuses on upstream gas fueling with start magnetic
% configuration (as of Jan 13,2016) Gobal field = 6450A, TR2=260A,
% TR1= 700A (but not working for today)
% Gas flow set points and timiing see ShotSummaries
% Reference shot = 12873
%12975 is the long pulse during MAB - Jan 2017


close all
clear all

% Radial Scan

%shotlist = [ 12873, 13400 + [ 76 77 ] ];
%shotlist =  13500 + [ 08 09 11 15 27 29 31 33 38 41 43 ] ; %w/o ECH
%shotlist = 13500 +  [ 37 10 12 26 28 30 32 34 40 42 44 ]; %w/ECH

% EBW coil current scan
%shotlist = 13500 + [ 49 50 52:60 ];
shotlist = 13500+[62:64];

% -------------------------
Config.tStart = 4.2; % [s]
Config.tEnd = 4.36;

% Acquiring Ne and Te data
Stem        = '\MPEX::TOP.';
Branch      = 'MACHOPS1:';
RootAddress = [Stem,Branch];

%% ECH data retrival from MDSplus
DataAddress{3} = [RootAddress, 'PWR_28GHZ'];
[EBW,t_28]     = my_mdsvalue_v2(shotlist,DataAddress(3));

%% Branches for neutral pressure measurements in MDSplus
CalShot = 13300 + [ 56 57 169 ];
CalShot = 13500 + [ 18 19 20 21];
% Cal for 150 ms - 13300+ 56, 57 6450 A
% Cal for 200 ms - 13300+ 136 5500A
% cal for 200 ms - 137 6450 A
% cal for 200 ms - 169 4000 A
% cal for 200 ms - 170 4400 A

DataAddress{4} = [ RootAddress,'PG4' ];
DataAddress{5} = [ RootAddress,'PG3' ];

[PG_4,t4]   = my_mdsvalue_v2( shotlist,DataAddress(4) );
[PG0_4,~]   = my_mdsvalue_v2( CalShot,DataAddress(4) );
[PG_6,t6]   = my_mdsvalue_v2( shotlist,DataAddress(5) );
[PG0_6,~]   = my_mdsvalue_v2( CalShot,DataAddress(5) );

%% Probe choices and variables

% AddressType='s'; % s for standard
% CalType='niso'; % niso for not isolated- "Standard DLP circuit box",
% CalType = 'iso' for isolated - "Transformer box"
% ----------
DLPType='M10';

switch DLPType
    case '4'
        DLP = 4.5;
        AddressType   = 'n';
        CalType       = 'iso';
        Config.L_tip  = 1.2/1000;
        Config.D_tip  = 0.508/1000; % [m]
    case '6'
        DLP = 6.5;
        AddressType  = 's';
        CalType      = 'niso';
        Config.L_tip = 1.2/1000;
        Config.D_tip = 0.254/1000; % [m]
    case '9'
        DLP = 9.5;
        AddressType  = 's';
        CalType      = 'niso';
        Config.L_tip = 1.1/1000;  % 1.1 mm after 11/21 probe change
        Config.D_tip = 0.254/1000; % [m]
    case '10'
        DLP = 10.5;
        AddressType   = 's';
        CalType       = 'niso';
        Config.L_tip  = 1.2/1000;            % John's probe
        Config.D_tip  = 0.254/1000; % [m]
    case 'M10'                              %MDLP at 10.5
        DLP = 10.5;
        AddressType  = 's';
        CalType      = 'niso';
        Config.L_tip = 2./1000;
        Config.D_tip = 0.254/1000; % [m]
    case 'MDLP1'
        DLP = 1.5;
        AddressType  = 'n';
        CalType      = 'iso';
        Config.L_tip = 2./1000;
        Config.D_tip = 0.254/1000; % [m]
    case 'MDLP4'
        DLP = 4.5;
        AddressType  = 'n';
        CalType      = 'iso';
        Config.L_tip = 1.8/1000;
        Config.D_tip = 0.254/1000; % [m]
end

switch AddressType
    case 's'
        DataAddress{1} = [ RootAddress,'LP_V_RAMP' ]; % V
        DataAddress{2} = [ RootAddress,'TARGET_LP' ]; % I
        
        Config.V_Att = 2;  % Output voltage of DLP box (Voltage) = V_Att*Digitized data
        Config.I_Att = 5;  % Output voltage of DLP box (Current) = I_att*Digitized data
    case 'n'
        DataAddress{1} = [ RootAddress,'INT_4MM_1' ]; % V
        DataAddress{2} = [ RootAddress,'INT_4MM_2' ]; % I
        Config.V_Att   = 1;  % Output voltage of DLP box (Voltage) = V_Att*Digitized data
        Config.I_Att   = 1;  % Output voltage of DLP box (Current) = I_att*Digitized data
end

switch CalType
    case 'iso'
        Config.V_cal = [ (0.46e-3)^-1,0 ]; % Voltage output of DLP = V_cal(1)*Output voltage of DLP box + V_cal(2) %  2.1739e+03
        Config.I_cal = [ -1,0 ]; % Current output of DLP = (I_cal(2) + Output voltage of DLP box)/Ical(1)
    case 'niso'
        Config.V_cal = [ 12.05,0.205 ];   % Voltage output of DLP = V_cal(1)*Output voltage of DLP box + V_cal(2)
        Config.I_cal = [ -142.5, 0.015 ]; % Current output of DLP = (I_cal(2) + Output voltage of DLP box)/Ical(1)
end

Config.FilterDataInput = 1; % Filter input data with savitsky Golay filter order 3 frame size 7
Config.SGF             = 7;
Config.TimeMode        = 2; % Effective time of sweep, (1) start of ramp or (2) mean time of ramp
Config.AMU             = 2; % Ion mass in AMU

Config.FitFunction     = 2;
Config.AreaType        = 1; % Cylindrical + cap
Config.e_c             = 1.602e-19;
Config.m_p             = 1.6726e-27;
Config.e_0             = 8.854187817*1e-12;
e_c                    = 1.602e-19;

%% Plasma Parameters calculation funtion

[ni,Te,time,Ifit,Ip,Vp,tm,Vsweep,Isweep] = DLP_fit_V5(Config,shotlist,DataAddress);

for s = 1:length(shotlist)
    Ni{s} = 0.5*(ni{s}{1} + ni{s}{2});
end

%% Plotting Te, ne, neutral pressures and (plasma pressures)
close all
C = {'k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:'};
t = t_zero(shotlist);

figure;
% Electron density
subplot(3,2,1); hold on
for s = 1:length(shotlist)
    % plot(time{s},ni{s}{1},C{s},'lineWidth',1);
    % plot(time{s},ni{s}{2},C{s},'lineWidth',1);
    h(s) = plot(time{s},Ni{s},C{s},'lineWidth',2);
        plot(t_28{s}(1:end-1), EBW{s}*.5e19, C{s}); % plot EBW waveform
end
title('$ n_e $ $ [m^{-3}] $','interpreter','Latex','FontSize',13,'Rotation',0)
legend(h,['DLP ',num2str(DLP)],'location','NorthWest')
set(gca,'Fontsize', 16,'FontWeight','Bold')
ylim([0,5e19])
xlim([4.15,4.45])
grid on


% Electron Temperature
subplot(3,2,2); hold on
for s= 1:length(shotlist)
    h(s) = plot(time{s},Te{s},C{s},'lineWidth',2);
    L{s} = [num2str(shotlist(s)),' ,t=',num2str(t{s}(10:14))];
    hold on
    
    set(gca,'Fontsize',12,'FontWeight','Bold')
end
legend(h,L,'location','NorthEast', 'Fontsize',13)
title('$ T_e $ $ [eV] $','interpreter','Latex','FontSize',13,'Rotation',0)
ylim([0,20])
xlim([4.15,4.45])
set(gca,'Fontsize', 18,'FontWeight','Bold')
grid on

% ===============
% Pressure at 6.5
subplot(3,2,4); hold on
for s = 1 :length(shotlist)
        h(s) = plot(t6{s}(1:end-1),(PG_6{s}-(PG0_6{1}+PG0_6{2}+PG0_6{3}+PG0_6{4})/4)*2, C{s},'lineWidth',2);
end
ylabel('mTorr')
xlabel('time(s)')
title('Baratron at SP 6.5')
set(gca,'Fontsize', 16,'FontWeight','Bold')
xlim([3.5,6.5])
%xlim([4.1,4.5])
grid on

% =================
%Pressure at SP 4.5

subplot(3,2,3); hold on
for s = 1 :length(shotlist)
        h(s) = plot(t4{s}(1:end-1),(PG_4{s}-(PG0_4{1}+PG0_4{2}+PG0_4{3}+PG0_4{4})/4)*2, C{s},'lineWidth',2);
end
%legend(h,{num2str(Shots')}, 'Fontsize', 16)
ylabel('mTorr')
xlabel('time(s)')
title('Baratron at SP 4.5')
xlim([3.5,6.5])
%xlim([4.1,4.5])
set(findobj('-Property','YTick'),'box','on')
set(gcf,'color','w')
set(gca,'Fontsize', 16,'FontWeight','Bold')
grid on

% ===============
% Plasma pressure
subplot(3,2,5); hold on
for s = 1:length(shotlist)
    plot(time{s},e_c.*Ni{s}.*Te{s},C{s},'lineWidth',2)
end
ylim([0,10])
xlim([4.15,4.45])
title('$ P_e $ $ [Pa] $','interpreter','Latex','FontSize',13,'Rotation',0)
grid on
set(gca,'Fontsize', 15,'FontWeight','Bold')

% ========================
% I_sat trace for the shot
if 1
    figure; hold on
    for s = 1:length(shotlist)
        %plot(tm{s},Vp{s})
        h(s) = plot(tm{s},Ip{s}*1000);
        ylim([-100,100])
        xlim([4.15,4.45])
        grid on
    end
    legend(h,num2str(shotlist'))
    set(gca,'Fontsize', 15,'FontWeight','Bold')
end

% ========================

if 0
    figure;
    for c = 1:25;
        subplot(5,5,c); hold on
        plot(Vsweep{s}{c},Isweep{s}{c}*1e3,'k')
        plot(Vsweep{s}{c},Ifit{s}{c}*1e3,'r')
        ht = title(['T_e: ',num2str(Te{s}(c))]);
        set(ht,'FontSize',10); grid on
    end
    
    figure;
    for c = 26:50;
        subplot(5,5,c-25); hold on
        plot(Vsweep{s}{c},Isweep{s}{c}*1e3,'k')
        plot(Vsweep{s}{c},Ifit{s}{c}*1e3,'r')
        ht = title(['T_e: ',num2str(Te{s}(c))]);
        set(ht,'FontSize',10); grid on
    end
    
    figure;
    for c = 51:75;
        subplot(5,5,c-50); hold on
        plot(Vsweep{s}{c},Isweep{s}{c}*1e3,'k')
        plot(Vsweep{s}{c},Ifit{s}{c}*1e3,'r')
        ht = title(['T_e: ',num2str(Te{s}(c))]);
        set(ht,'FontSize',10); grid on
    end
end


%% Selecting portion of the code

s=1
t_LD_Gather_Start = 4.3;
t_LD_Gather_End   = 4.33;
%%
t_LD_Gather_Start = 4.23;
t_LD_Gather_End   = 4.245;
%%

idx = find(time{1} >= t_LD_Gather_Start  & time{1} <= t_LD_Gather_End)

% Standard deviation

% Error analysis within the shot
for s = 1:length(shotlist)
    % Temperature before ECH
    %-----------------------
    x11 = Te{s}(idx); 
    SEM11 = std(x11); % /sqrt(length(x1));              % Standard Error
    e11(s) = SEM11;
    Te11(s) = mean(x11);
    
    % Density
    % -----------
    x21 = Ni{s}(idx);
    SEM21 = std(x21);              % Standard Error
    Ne21(s)= mean(x21);
    e21(s) = SEM21;
    
    summary = [Te11; e11; Ne21; e21]'
end
return
%%
coil = [
    4000
    3900
    3800
    3700
    3600
    3500
    3400
    4100
    4200
    4300
    4400
    4500
    4600
    4700
    4900
    ];

figure %Radial Temperature
plot(coil, Te11, 'bx--', 'MarkerSize',8)
errorbar(coil,Te11,e11, 'bx--', 'LineWidth',1)
hold on;
xlabel('Coil Current [A]')
ylim([1 5])
ylabel('kTe [eV]')
set(gca,'Fontsize', 20,'FontWeight','Bold')
title('Te vs coil current after EBW')
grid on

figure %Radial Density
plot(coil,Ne21, 'bx', 'MarkerSize',8);
errorbar(coil, Ne21, dNeU, 'bx--', 'LineWidth',1)
hold on
xlabel('Coil Current [A]')
ylabel('ne [m^{-3}]')
%legend('kTe Before ECH','kTe During ECH', 'Location','NorthWest')
set(gca,'Fontsize', 20,'FontWeight','Bold')
title('Ne vs coil current after EBW')
grid on
return

