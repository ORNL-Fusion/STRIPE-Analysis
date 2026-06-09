clear all
close all

ProbeLoc = 'A';
FitShow = 0; 
RawDataShow = 1;

% #########################################################################
% 2019_11_01
ShotType = 3;

switch ShotType
    case 1
    case 2
    case 3 % RAdial scan
        shotlist = [28000 + [96,  97,98,  99,100,101 ,102 ,103 ,104,105,106,108,109,110,111]];
        R        =          [10,10.5,11,11.5,12 ,12.5,13.0,13.5,10 ,9.5,9.0,8.5,8.0,7.5,7.0] ;
        shotlist = [28000 + [146 ,147,148,149,150,151,152,153,154,155 ,156 ,158 ,159 ]];
        R        =          [10.5,9.5,8.5,7.5,6.5,6.0,7.0,8.0,9.0,11.5,11.0,12.5,13.5];
    case 4 
    case 5
    case 6
    case 7
end
 
Config.AMU = 2; % Ion mass in AMU
DLPType = '12MP';

SweepType = 'iso'  ;AttType = 'Vx1,Ix1';
% #########################################################################

Config.FitFunction = 2; 
ChannelType = '6' ;
Config.tStart = 4.15; % [s]
Config.tEnd = Config.tStart + 0.5;
Config.Center_V = 1; % Remove offset on V: 1 (yes) 0(no)
Config.Center_I = 1; % Remove offset on I: 1 (yes) 0(no)
    

%%
% #########################################################################
% DLP FITTING ROUTINE
% #########################################################################
switch DLPType 
    case '4'
        DLP = 4.5;
        Config.L_tip = 1.2/1000; % 1.0 as of April 11th 2017
        Config.D_tip = 2*0.254/1000; % [m]
    case '9'
        DLP = 9.5;
        Config.L_tip = 1.01/1000; % 
        Config.D_tip = 0.254/1000; % [m]
    case '6MP'
        DLP = 6.5;
        Config.L_tip = 1.0/1000; % 1.0 as of April 11th 2017
        Config.L_tip = 1.4/1000; % 1.4 as of Nov 29th 2018
        Config.L_tip = 1.7/1000; % 1.7 as of Dec 13th 2018
        Config.D_tip = 0.254/1000; % [m]
    case '10'
        DLP = 10.5;
        Config.L_tip = 1.2/1000;
        Config.D_tip = 0.254/1000; % [m]
    case '10MP'
        DLP = 10.5;
        Config.L_tip = 1.8/1000;
%         Config.L_tip = 1.6/1000;% according to Nischal 2017_11_21
        Config.L_tip = 0.75/1000; % As of 2018_08_27 after Magnetic field upgrade
        Config.D_tip = 0.254/1000; % [m]
    case '4MP'
        DLP = 4.5;
        Config.L_tip = 1.7/1000;
        Config.D_tip = 0.254/1000; % [m]
    case '1MP'
        DLP = 1.5;
        Config.L_tip = 1.7/1000;
        Config.L_tip = 1.1/1000; % As of 2018_08_31
        Config.L_tip = 2.8/1000; % As of 2018_09_19
        Config.L_tip = 2.4/1000; % As of 2018_12_21
        Config.L_tip = 2.25/1000; % As of 2019_02_07
        Config.D_tip = 0.254/1000; % [m]
    case 'TDLP'
        DLP = 11.5;
        Config.L_tip = 0/1000;
        Config.D_tip = 0.8/1000; % [m]
    case 'IF'
        DLP = 11.5;
        Config.L_tip = 0/1000;
        Config.D_tip = 1/1000; % [m]
    case '10FluxProbe';
        DLP = 10.5;
        Config.L_tip = 0/1000;
        Config.D_tip = ((2^(1/4))*1.22)/1000; % [m]
    case '12MP';
        DLP = 12.5;
        Config.L_tip = 1.1/1000; 
        Config.L_tip = 2.4/1000; % As of 2018_10_24 
        Config.L_tip = 2/1000; % As of 2019_08_08 
        Config.L_tip = 2.1/1000; % As of 2019_10_09 
        Config.D_tip = 0.254/1000; % [m]
   case '8MP';
        DLP = 8.5;
        Config.L_tip = 2/1000;
        Config.L_tip = 2.2/1000; % 2018_11_06
        Config.L_tip = 1.5/1000; % 2018_12_04
        Config.D_tip = 0.254/1000; % [m]
end
% #########################################################################
switch AttType
    case 'Vx1,Ix5'
        Config.V_Att = 1;  % Output voltage of DLP box (Voltage) = V_Att*Digitized data 
        Config.I_Att = 5;  % Output voltage of DLP box (Current) = I_att*Digitized data 
    case 'Vx1,Ix2'
        Config.V_Att = 1;  % Output voltage of DLP box (Voltage) = V_Att*Digitized data 
        Config.I_Att = 2;  % Output voltage of DLP box (Current) = I_att*Digitized data    
    case 'Vx2,Ix5'
        Config.V_Att = 2;  % Output voltage of DLP box (Voltage) = V_Att*Digitized data 
        Config.I_Att = 5;  % Output voltage of DLP box (Current) = I_att*Digitized data
    case 'Vx1,Ix1'
        Config.V_Att = 1;  % Output voltage of DLP box (Voltage) = V_Att*Digitized data 
        Config.I_Att = 1;  % Output voltage of DLP box (Current) = I_att*Digitized data
end
% #########################################################################
Stem = '\MPEX::TOP.'; Branch = 'MACHOPS1:'; RootAddress = [Stem,Branch];
switch ChannelType
    case '1'
        DataAddress{1} = [RootAddress,'LP_V_RAMP']; % V
        DataAddress{2} = [RootAddress,'TARGET_LP']; % I
    case '2'
        DataAddress{1} = [RootAddress,'INT_4MM_2']; % V
        DataAddress{2} = [RootAddress,'INT_4MM_1']; % I
    case '3'
        DataAddress{1} = [RootAddress,'GEN_RF_PWR']; % V
        DataAddress{2} = [RootAddress,'INT_4MM_2'];  % I
    case '4'
        DataAddress{1} = [RootAddress,'LP_V_RAMP']; % V
        DataAddress{2} = [RootAddress,'INT_4MM_2'];  % I
    case '5'
        DataAddress{1} = [RootAddress,'LP_V_RAMP']; % V
        DataAddress{2} = [RootAddress,'EA_CURRENT'];  % I
    case '6'
        DataAddress{1} = [RootAddress,'LP_V_RAMP']; % V
        DataAddress{2} = [RootAddress,'LP_1'];  % I
end
% #########################################################################
switch SweepType
    case 'iso'
Config.V_cal = [(0.46e-3)^-1,0]; % Voltage output of DLP = V_cal(1)*Output voltage of DLP box + V_cal(2) %  2.1739e+03
Config.I_cal = [-1,0]; % Current output of DLP = (I_cal(2) + Output voltage of DLP box)/Ical(1)
    case 'niso'
Config.V_cal = [12.05,0.205];   % Voltage output of DLP = V_cal(1)*Output voltage of DLP box + V_cal(2)
Config.I_cal = [-142.5, -0.6 + 0.015]; % Current output of DLP = (I_cal(2) + Output voltage of DLP box)/Ical(1)
Config.I_cal = [-142.5, 0.015]; % Current output of DLP = (I_cal(2) + Output voltage of DLP box)/Ical(1)
end
% #########################################################################
Config.FilterDataInput = 1; % Filter input data with savitsky Golay filter order 3 frame size 7
Config.SGF = 11;
Config.TimeMode = 2; % Effective time of sweep, (1) start of ramp or (2) mean time of ramp
Config.AreaType = 1; % 1: Cylindrical + cap
% #########################################################################
% [ni,Te,time,Ifit,Ip,Vp,tm,Vsweep,Isweep,GlitchFlag,SSQres] = DLP_fit_V5_2(Config,shotlist,DataAddress);
[ni,Te,Isat,time,Ifit,Ip,Vp,tm,Vsweep,Isweep,GlitchFlag,SSQres,StdRes,StdResNorm] = DLP_fit_V5_7(Config,shotlist,DataAddress);

%Ni,Te,time,Ifit,Ip,Vp,tm,Vsweep,Isweep,GlitchFlag,SSQres
for s = 1:length(shotlist)
    Ni{s} = 0.5*(ni{s}{1} + ni{s}{2});
    Is{s} = 0.5*(Isat{s}{1} + Isat{s}{2});
end
% #########################################################################
DA{1} = [RootAddress,'PWR_28GHZ']; % Isx
DA{1} = [RootAddress,'PWR_28GHZ']; % Isx
[ECH,t_ech]   = my_mdsvalue_v2(shotlist,DA(1))
% % #########################################################################
DA{1} = [RootAddress,'RF_FWD_PWR'];
DA{1} = [RootAddress,'RF_CURRENT'];
[RF,t_rf]   = my_mdsvalue_v2(shotlist,DA(1))
% #########################################################################
DA{1} = [RootAddress,'ICH_FWD_PWR'];
[ICH_F,t_ICH_F]   = my_mdsvalue_v2(shotlist,DA(1))
DA{1} = [RootAddress,'ICH_REF_PWR'];
[ICH_R,t_ICH_R]   = my_mdsvalue_v2(shotlist,DA(1))
% #########################################################################
DA{1} = [RootAddress,'PS3_V'];
[PS3_V,t_PS3_V]   = my_mdsvalue_v2(shotlist,DA(1))
%%
% #########################################################################
% PLOT DATA:
% #########################################################################

% close all
% Add code that looks at the date and then looks for a GasCompare file with
% that name

% -------------------------------------------------------------------------
% The following script checks the existence of the "GasCompare" script and
% runs if it exists
close all

Yscale = 'lin';
ScriptList = dir('*.m');
for n = 1:length(ScriptList)
    if strcmp(ScriptList(n).name(1:13),'GasCompareALL');
        run(ScriptList(n).name)
        break % the script was been found
    else
%         error('Could not find GasCompare script')
    end
end
% -------------------------------------------------------------------------

C = {'k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:'};
t = t_zero(shotlist);

TimePlotStart = 4.14;
TimePlotEnd = TimePlotStart + 0.6;

figure; hold on
for s = 1:length(shotlist)
    Vpwr = RF{s};
    RFPWR = 10.^( ((27*Vpwr)-32.54)/10 + 7.2);
    plot(t_rf{s}(1:end-1),(RFPWR-mean(RFPWR(1:15)))*1e-3,C{s})
end
set(gcf,'color','w')
box on
grid on
title('Raw RF detector output ')
xlim([TimePlotStart,TimePlotEnd])

figure; 
subplot(2,2,1); hold on
for s = 1:length(shotlist)
    GoodFits{s} = GlitchFlag{s} == 0 & StdResNorm{s}<=0.5 & Ni{s}>0 & Ni{s}<50e19 & Te{s}<=30 & Is{s}>0.01e-3;
      
    t_Ni_gf{s} = time{s}(GoodFits{s} & Ni{s}>1e16);
    Ni_gf{s} = Ni{s}(GoodFits{s} & Ni{s}>1e16);
    
    h(s) = plot(t_Ni_gf{s},Ni_gf{s},C{s},'lineWidth',2);

    hold on
    plot(t_ech{s}(1:length(ECH{s})),abs(ECH{s})*8e19,C{s},'lineWidth',0.5)
    
     plot(t_rf{s}(1:length(RF{s})),(RF{s}.^2)*3e19,C{s},'lineWidth',1.5)
     plot(t_ICH_F{s}(1:length(ICH_F{s})),(ICH_F{s}.^2)*1e19,C{s},'lineWidth',1.5)
%      plot(t_ICH_R{s}(1:length(ICH_R{s})),(ICH_R{s}.^2)*3e19,C{s},'lineWidth',1.5)
end
title('$ n_e $ $ [m^{-3}] $','interpreter','Latex','FontSize',13,'Rotation',0)
legend(h,['DLP ',num2str(DLP)],'location','SouthEast')
%set(gca,'PlotBoxAspectRatio',[1 1 1])
% ylim([0,7e18])
ylim([0,10e19])
ylim([0,2e20])
% ylim([0,5e20])
xlim([TimePlotStart,TimePlotEnd])

yyaxis right
for s = 1:length(shotlist)
        plot(t_PS3_V{s}(1:end-1),PS3_V{s}*382.6,C{s},'lineWidth',1.5,'LineStyle',':')
end
ylim([0,1000])
ylabel('PS3 [A]')
xlim([TimePlotStart,TimePlotEnd])

subplot(2,2,2); hold on
for s = 1:length(shotlist)
        t_Te_gf{s} = time{s}(GoodFits{s} & Ni{s}>1e16);
        Te_gf{s} = sgolay_t(Te{s}(GoodFits{s} & Ni{s}>1e16),3,7);
        
        h(s) = plot(t_Te_gf{s},Te_gf{s},C{s},'lineWidth',2);
        
        L{s} = [num2str(shotlist(s)),' ,t=',num2str(t{s}(10:14))];
       plot(t_ech{s}(1:length(ECH{s})),2*abs(ECH{s}),C{s},'lineWidth',0.5)

end
legend(h,L,'location','NorthEast')
title('$ T_e $ $ [eV] $','interpreter','Latex','FontSize',13,'Rotation',0)
ylim([0,6])
xlim([TimePlotStart,TimePlotEnd])

subplot(2,2,3); hold on
for s = 1:length(shotlist)
    Pe_gf{s} = e_c.*Ni_gf{s}.*Te_gf{s};
    plot(t_Ni_gf{s},Pe_gf{s},C{s},'lineWidth',2)
end
ylim([0,60])
xlim([TimePlotStart,TimePlotEnd])
title('$ P_e $ $ [Pa] $','interpreter','Latex','FontSize',13,'Rotation',0)

set(findobj('-Property','YTick'),'box','on')
set(gcf,'color','w')

if RawDataShow
    figure; hold on
    for s = 1:length(shotlist)
    plot(tm{s},Vp{s})
    h(s) = plot(tm{s},(Ip{s})*1000);
    ylim([-500,500])
    xlim([4.15,TimePlotEnd])
    grid on
    end
    legend(h,num2str(shotlist'))
    set(gcf,'color','w')
    box on
end

%% Convert data into Excel
% Save data to Excel spreadsheet:
% Create table to export to EXCEL
if 0
    for s = 1:length(shotlist)
        G = [t_Ni_gf{s}',Ni_gf{s}',Te_gf{s}'];
        F = {'time [s]','ne[m^-3]','Te[eV]'};
        FileName = ['Shot_',num2str(shotlist(s)),'_NeTe_Spool_10_5.xlsx']
        xlswrite(FileName,[F;num2cell(G)]);
    end
end

%%
% close all

Ymax = 150; % mA

if FitShow
    for k = 1:15
         figure;
    for c = (1 + (k-1)*25):(25 + (k-1)*25);
        subplot(5,5,c- (k-1)*25 ); hold on
        try
        plot(Vsweep{s}{c},Isweep{s}{c}*1e3,'k')
        if GoodFits{s}(c) == 1
            plot(Vsweep{s}{c},Ifit{s}{c}*1e3,'r')
        else
            plot(Vsweep{s}{c},Ifit{s}{c}*1e3,'g')
        end
        ht = title(['T_e: ',num2str(Te{s}(c),2),' ,t = ',num2str(time{s}(c),5),' s']);
%         ht = title(['I_s: ',num2str(Is{s}(c),2),' ,t = ',num2str(time{s}(c),5),' s']);
        set(ht,'FontSize',6); grid on
        set(gcf,'color','w')
        set(gca,'FontSize',5)
        box on
%         ylim(Ymax*[-1,1])
        catch
            warning('error')
            continue
        end
    end
    end
end

% return
%%
% #########################################################################
% Steady state values
% #########################################################################
close all

for s = 1:length(shotlist)
    rngT{1} = [time{s}>=4.28 & time{s}<=4.38]; 
    rngT{2} = [time{s}>=4.42 & time{s}<=4.44];      
for p = 1:length(rngT)
        nes{p}(s)    = mean(Ni{s}(rngT{p} & GoodFits{s}));
        dnes{p}(s)   = std(Ni{s}(rngT{p} & GoodFits{s}),1,2);
        Tes{p}(s)    = mean(Te{s}(rngT{p} & GoodFits{s}));
        dTes{p}(s)   = std(Te{s}(rngT{p} & GoodFits{s}),1,2);
    end
end

R_offset = 10.5

if 1
    figure; 
    subplot(3,1,1); hold on
    [~,b] = sort(R);
    errorbar(R(b)-R_offset,nes{1}(b),dnes{1}(b),'ko-')
    errorbar(R(b)-R_offset,nes{2}(b),dnes{2}(b),'ro-')
    ylim([0,1.5e20])
    ylim([0,1.0e20])
    xlim([-5,5])
    ylabel('$$n_e$$ $$[m^{-3}]$$','Interpreter','latex','FontSize',12)
    grid on
    box on
    
    subplot(3,1,2); hold on
    errorbar(R(b)-R_offset,Tes{1}(b),dTes{1}(b),'ko-')
    errorbar(R(b)-R_offset,Tes{2}(b),dTes{2}(b),'ro-')
    ylim([0,6])
    ylabel('$$T_e$$ [eV]','Interpreter','latex','FontSize',12)
    grid on
    box on
        xlim([-5,5])
    
    subplot(3,1,3); hold on
    errorbar(R(b)-R_offset,e_c*Tes{1}(b).*nes{1}(b),e_c*dTes{1}(b).*dnes{1}(b),'ko-')
    errorbar(R(b)-R_offset,e_c*Tes{2}(b).*nes{2}(b),e_c*dTes{2}(b).*dnes{2}(b),'ro-')
    ylim([0,50])
    ylabel('$$P_e$$ [Pa]','Interpreter','latex','FontSize',12)
    xlabel('$$R$$ [cm]','Interpreter','latex','FontSize',12)
    grid on
    box on
        xlim([-5,5])
    
    set(gcf,'color','w')
end

% Plot time traces
figure
for s = 1:length(R)
    subplot(4,4,s)
        h(s) = plot(t_Ni_gf{s},Ni_gf{s},C{s},'lineWidth',2);
        hold on
        plot(t_ech{s}(1:length(ECH{s})),abs(ECH{s})*8e19,C{s},'lineWidth',0.5)
        ylim([0,2e20])
        xlim([4.1,4.8])
        title(['R: ',num2str(R(s)-R_offset)])
        box on
        grid on
end
set(gcf,'color','w')

figure
for s = 1:length(R)
    subplot(4,4,s)
        h(s) = plot(t_Ni_gf{s},Te_gf{s},C{s},'lineWidth',2);
        hold on
        plot(t_ech{s}(1:length(ECH{s})),2*abs(ECH{s}),C{s},'lineWidth',0.5)
        ylim([0,6])
        xlim([4.1,4.8])
        title(['R: ',num2str(R(s)-R_offset)])
        box on
        grid on
end
set(gcf,'color','w')


%% Convert data into Excel
% Save data to Excel spreadsheet:
% Create table to export to EXCEL
if 0
G = [shotlist(b)',PS1(b)',ne{p}(b)',dne{p}(b)',Te{p}(b)',dTe{p}(b)'];
F = {'Shot','R [cm]','ne[m^-3]','dne','Te[eV]','dTe'};
FileName = 'NeTe_Spool_10_5_2018_01_16.xlsx';
xlswrite(FileName,[F;num2cell(G)]);
end