% In this script I run two DLP fitting codes, one that I have written and
% Rick Goulding's code as of Aug 22nd 2016.

% Script written by Juan Caneses Oct 27th 2016 unless otherwise stated

close all
clear all
% ########################################################################
%  PROBE VALUES:
% ########################################################################
% On shot 11098, DLP 10.5 was inserted. the day before, we had ran the 2s
% 100 kW pulses and damaged DLP 10.5. upon inspection by John C, the probe
% tip length is about 0.5 mm in length with a 25% uncertainty


shot = 11000 + [98 ,  99,  123];
L_tip =            [0.5, 0.5, 1.2]/1000; % probe tip length in mm, 
fctr  =  1.2./(L_tip*1000);

% 11098: DLP 10.5
% 11099: DLP 10.5
% 11023: DLP 9.5

AMU = 2; % D2 gas
D_tip = 0.254/1000; % [m]
rp = D_tip/2; % [m]

% ########################################################################
% MDS+ ADDRESS:
% ########################################################################\
Address = 1;
switch Address
    case 1
        Stem = '\MPEX::TOP.';
        Branch = 'MACHOPS1:';
        RootAddress = [Stem,Branch];
        DataAddress{1} = [RootAddress,'LP_V_RAMP']; % V
        DataAddress{2} = [RootAddress,'TARGET_LP']; % I
    case 2
        Stem = '\MPEX::TOP.';
        Branch = 'MPEX1:';
        RootAddress = [Stem,Branch];
        DataAddress{1} = [RootAddress,'RF_REF_PWR']; % V
        DataAddress{2} = [RootAddress,'TARGET_LP']; % I
end

for s = 1:length(shot)
    shotlist = shot(s);
    % ########################################################################
    % DLP fitting procedure condition configuration variables
    % ########################################################################
    Config.tStart = 4.23; % [s]
    Config.tEnd = 4.32;
    Config.V_Att = 2;  % Output voltage of DLP box (Voltage) = V_Att*Digitized data 
    Config.I_Att = 2;  % Output voltage of DLP box (Current) = I_att*Digitized data
    Config.V_cal = [12.05,0.205];   % Voltage output of DLP = V_cal(1)*Output voltage of DLP box + V_cal(2)
    Config.I_cal = [-142.5, 0.015]; % Current output of DLP = (I_cal(2) + Output voltage of DLP box)/Ical(1)
    Config.FilterDataInput = 1; % Filter input data with savitsky Golay filter order 3 frame size 7
    Config.TimeMode = 2; % Effective time of sweep, (1) start of ramp or (2) mean time of ramp
    Config.AMU = 2; % Ion mass in AMU
    Config.L_tip = L_tip(s);
    Config.D_tip = D_tip;
    Config.FitFunction = 2; 

    % ########################################################################
    % JUAN CANESES,DATA ANALYSIS:
    % ########################################################################
    Config.AreaType = 1; % Cylindrical + cap
    [Ni_a{s},Te_a{s},time_a{s},Ifit{s},I_dlp{s},V_dlp{s},t_dlp{s},Vsweep{s},Isweep{s}]...
        = DLP_fit_V4(Config,shotlist,DataAddress);

    Config.AreaType = 2; % Cylindrical no cap
    [Ni_b{s},Te_b{s},time_b{s},~,~,~,~,~,~] = DLP_fit_V4(Config,shotlist,DataAddress);

    % ########################################################################
    % RICK GOULDING, DATA ANALYSIS:
    % ########################################################################
    % the output of Rick's code is ne [cm-3], te [V], tim(:,1) [s] 
    plot_ne_te_spool9p5_050616
    % To correct for the smaller probe area (0.5 mm and not 0.12) we use a
    % multiplication factor
    Time{s} = tim;
    Te{s} = te;
    Ne{s} = ne*fctr(s)*1e6; % Convert from cm-3 to m-3
end


%% Plot data:
%close all

C = {'Black','Red','Blue','Green'};
PlotTimeStart = 4.16; 
PlotTimeEnd = 4.32;


for s = 1:length(shot)
    
figure; 
subplot(2,2,1); hold on

hJC_a(s) = plot(time_a{s}{1},Ni_a{s}{1},'k');
hJC_b(s) = plot(time_b{s}{1},Ni_b{s}{1},'g','LineStyle','-');
hRG(s) = plot(Time{s},Ne{s},'r','LineWidth',1,'LineStyle','-')
legend([hJC_a(s),hJC_b(s),hRG(s)],{[num2str(shot(s)'),...
    ' Caneses "a" (w/ cap)'],[num2str(shot(s)'),' Caneses "b" (no cap)'],[num2str(shot(s)'),' Goulding']})
     
title('$ n_e $ $ [m^{-3}] $','interpreter','Latex','FontSize',13,'Rotation',0)
xlim([PlotTimeStart,PlotTimeEnd])
ylim([0,8e19])
ax(1) = gca;

subplot(2,2,2); hold on
plot(time_a{s}{1},Te_a{s}{1},'k')
plot(time_b{s}{1},Te_b{s}{1},'g','LineStyle','-');
plot(Time{s},Te{s},'r','LineWidth',1,'LineStyle','-')
title('$ T_e $ $ [eV] $','interpreter','Latex','FontSize',13,'Rotation',0)
ylim([0,8])
xlim([PlotTimeStart,PlotTimeEnd])
ax(2) = gca;

subplot(2,2,3); hold on
plot(time_a{s}{1},Ni_a{s}{1}.*Te_a{s}{1}*e_c,'k')
plot(time_b{s}{1},Ni_b{s}{1}.*Te_b{s}{1}*e_c,'g','LineStyle','-')
plot(Time{s},Ne{s}.*Te{s}*e_c,'r','LineWidth',1,'LineStyle','-')
title('$ P_e $ $ [Pa] $','interpreter','Latex','FontSize',13,'Rotation',0)
set(gcf,'position',[417 106 574 503])
ylim([0,15])
xlim([PlotTimeStart,PlotTimeEnd])
ax(3) = gca;

set(gcf,'color','w')
set(ax,'Box','on')

end

%% Fit check:
if 1
    figure;     hold on
    for s = 1:length(shotlist);
    plot(t_dlp{s}{1},V_dlp{s}{1})
    plot(t_dlp{s}{1},I_dlp{s}{1}*1000)
    end
    grid on
    ylim([-100,100])
end

if 1
    
     figure;
    for c = 1:25;
        subplot(5,5,c); hold on
        plot(Vsweep{s}{c}{1},Isweep{s}{c}{1}*1e3,'k')
        plot(Vsweep{s}{c}{1},Ifit{s}{c}{1}*1e3,'r')
        ht = title(['T_e: ',num2str(Te_a{s}{1}(c)),'t: ',num2str(time_a{s}{1}(c))]);
        set(ht,'FontSize',5)
    end
    
    figure;
        for c = 26:50;
        subplot(5,5,c-25); hold on
        plot(Vsweep{s}{c},Isweep{s}{c}*1e3,'k')
        plot(Vsweep{s}{c},Ifit{s}{c}*1e3,'r')
        ht = title(['T_e: ',num2str(Te_a{s}(c)),'t: ',num2str(time_a{s}(c))]);
        set(ht,'FontSize',5)
        end
        
    figure;
        for c = 51:75;
        subplot(5,5,c-50); hold on
        plot(Vsweep{s}{c},Isweep{s}{c}*1e3,'k')
        plot(Vsweep{s}{c},Ifit{s}{c}*1e3,'r')
        ht = title(['T_e: ',num2str(Te_a{s}(c)),'t: ',num2str(time_a{s}(c))]);
        set(ht,'FontSize',5)
        end  
    
end
