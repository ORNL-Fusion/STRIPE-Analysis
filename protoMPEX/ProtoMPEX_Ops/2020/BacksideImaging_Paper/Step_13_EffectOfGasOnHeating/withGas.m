clear all
close all

ProbeLoc = 'A'; % DLP 10.5
FitShow = 1; 
RawDataShow = 0;

switch ProbeLoc
    case 'A'

    shotlist  = [21000 + [380,381,382,383,384,385,386,387,388,389,390,391,392,393]];
    R         =          [6.5,7.0,7.5,8.0,8.5,9.0,9.5,6.0,5.5,5.0,4.5,4.0,3.5,3.0] ;
% use 0.13 for stdresnorm

% shotlist  = [21000 + [383]];
% R         =          [8.0] ;

% =========================================================================
% % ECH radial scan with no D2 puff
% shotlist  = [21000 + [394,395,396,398,399,400]];
% R         =          [3.0,3.5,4.0,4.5,5.0,5.5] ;
% 
% shotlist  = [21000 + [400,401,402,403,404,405,406]];
% R         =          [5.5,6.0,6.0,6.5,7.0,7.5,8.0] ;
% 
% shotlist  = [21000 + [406,408,409,410]];
% R         =          [8.0,8.5,9.0,9.5] ;
% 
 shotlist  = [21000 + [395,396,398,399,400,401,402,403,404,405,406,408,409]];
 R         =          [3.5,4.0,4.5,5.0,5.5,6.0,6.0,6.5,7.0,7.5,8.0,8.5,9.0] ;
% % the data is very noisy for this data set
% =========================================================================

    SweepType = 'iso'; % "iso" for isolated sweep, "niso" for non-isolated sweep
    AttType = 'Vx1,Ix1'; % Attenuation on the digitized signals
    DLPType = '10MP';
    Config.FitFunction = 2; 
    ChannelType = '5' ;
    Config.tStart = 4.15; % [s]
    Config.tEnd = 5.5;
    Config.Center_V = 0; % Remove offset on V: 1 (yes) 0(no)
    Config.Center_I = 0; % Remove offset on I: 1 (yes) 0(no)
    
case 'B'        
    Config.FitFunction = 2; 

    DLPType = '1MP'; % DLP 1.5 horizontal
    ChannelType = '5' ; % "1" for TARGET_LP,LP_V_RAMP, "2" for other options

    Config.tStart = 4.25; % [s]
    Config.tEnd = 4.7;
    SweepType = 'iso'; % "iso" for isolated sweep, "niso" for non-isolated sweep
    AttType = 'Vx1,Ix1'; % Attenuation on the digitized signals
    Config.Center_V = 0; % Remove offset on V: 1 (yes) 0(no)
    Config.Center_I = 0; % Remove offset on I: 1 (yes) 0(no)
end

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
    case '6'
        DLP = 6.5;
        Config.L_tip = 1.0/1000; % 1.0 as of April 11th 2017
        Config.D_tip = 0.254/1000; % [m]
    case '10'
        DLP = 10.5;
        Config.L_tip = 1.2/1000;
        Config.D_tip = 0.254/1000; % [m]
    case '10MP'
        DLP = 10.5;
        Config.L_tip = 1.8/1000;
        Config.D_tip = 0.254/1000; % [m]
    case '4MP'
        DLP = 4.5;
        Config.L_tip = 1.7/1000;
        Config.D_tip = 0.254/1000; % [m]
    case '1MP'
        DLP = 1.5;
        Config.L_tip = 1.7/1000;
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
end
% #########################################################################
switch AttType
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
Config.AMU = 2; % Ion mass in AMU
Config.AreaType = 1; % 1: Cylindrical + cap
% #########################################################################
% [ni,Te,time,Ifit,Ip,Vp,tm,Vsweep,Isweep,GlitchFlag,SSQres] = DLP_fit_V5_2(Config,shotlist,DataAddress);
[ni,Te,Isat,time,Ifit,Ip,Vp,tm,Vsweep,Isweep,GlitchFlag,SSQres,StdRes,StdResNorm] = DLP_fit_V5_6(Config,shotlist,DataAddress);

%Ni,Te,time,Ifit,Ip,Vp,tm,Vsweep,Isweep,GlitchFlag,SSQres
for s = 1:length(shotlist)
    Ni{s} = 0.5*(ni{s}{1} + ni{s}{2});
    Is{s} = 0.5*(Isat{s}{1} + Isat{s}{2});
end
% #########################################################################
DA{1} = [RootAddress,'PWR_28GHZ']; % Isx
[ECH,t_ech]   = my_mdsvalue_v2(shotlist,DA(1))
% % #########################################################################
DA{1} = [RootAddress,'RF_FWD_PWR'];
[RF,t_rf]   = my_mdsvalue_v2(shotlist,DA(1))
% #########################################################################

%%
% #########################################################################
% PLOT DATA:
% #########################################################################

close all
C = {'k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:'};
t = t_zero(shotlist);

TimePlotStart = 4.1;
TimePlotEnd = 5.5;

figure; 
subplot(2,2,1); hold on
for s = 1:length(shotlist)
    GoodFits{s} = GlitchFlag{s} == 0 & StdResNorm{s}<=0.12 & Ni{s}>0 & Ni{s}<12e19 & Te{s}<=15 & Is{s}>0.01e-3;
            
    t_Ni_gf{s} = time{s}(GoodFits{s});
    Ni_gf{s} = Ni{s}(GoodFits{s});
    
    h(s) = plot(t_Ni_gf{s},Ni_gf{s},C{s},'lineWidth',2);

    hold on
    plot(t_ech{s}(1:length(ECH{s})),ECH{s}*0.1e19,C{s},'lineWidth',0.5)
    
    if s == 2
        continue
    end

end
title('$ n_e $ $ [m^{-3}] $','interpreter','Latex','FontSize',13,'Rotation',0)
legend(h,['DLP ',num2str(DLP)],'location','SouthEast')
ylim([0,1e19])
xlim([TimePlotStart,TimePlotEnd])

subplot(2,2,2); hold on
for s = 1:length(shotlist)
        t_Te_gf{s} = time{s}(GoodFits{s});
        Te_gf{s} = sgolay_t(Te{s}(GoodFits{s}),3,7);
        
        h(s) = plot(t_Te_gf{s},Te_gf{s},C{s},'lineWidth',2);
        
        L{s} = [num2str(shotlist(s)),' ,t=',num2str(t{s}(10:14))];
        plot(t_ech{s}(1:length(ECH{s})),2*ECH{s},C{s},'lineWidth',0.5)

end
legend(h,L,'location','NorthEast')
title('$ T_e $ $ [eV] $','interpreter','Latex','FontSize',13,'Rotation',0)
ylim([0,10])
xlim([TimePlotStart,TimePlotEnd])

subplot(2,2,3); hold on
for s = 1:length(shotlist)
    Pe_gf{s} = e_c.*Ni_gf{s}.*Te_gf{s};
    plot(t_Ni_gf{s},Pe_gf{s},C{s},'lineWidth',2)
end
ylim([0,30])
xlim([TimePlotStart,TimePlotEnd])
title('$ P_e $ $ [Pa] $','interpreter','Latex','FontSize',13,'Rotation',0)

set(findobj('-Property','YTick'),'box','on')
set(gcf,'color','w')

if RawDataShow
    figure; hold on
    for s = 1:length(shotlist)
    plot(tm{s},Vp{s})
    h(s) = plot(tm{s},(Ip{s})*1000);
    ylim([-300,300])
    xlim([4.15,TimePlotEnd])
    grid on
    end
    legend(h,num2str(shotlist'))
    set(gcf,'color','w')
    box on
end

%% Show fits:
Ymax = 150; % mA

if FitShow
    for k = 1:16
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

%%
% #########################################################################
% Steady state values
% #########################################################################
% close all

for s = 1:length(shotlist)
    
    % for D2 puff radial scan
%     rngT{1} = [time{s}>=4.9 & time{s}<=4.95]; 
    rngT{2} = [time{s}>=4.42 & time{s}<=4.43]; 
    rngT{3} = [time{s}>=4.69 & time{s}<=4.711]; 
    
    % for radial scan without D2 puff
    
    rngT{2} = [time{s}>=4.39 & time{s}<=4.54];  
    rngT{3} = [time{s}>=4.71 & time{s}<=4.78]; 

for p = 2:3
        nes{p}(s)    = mean(Ni{s}(rngT{p} & GoodFits{s}));
        dnes{p}(s)   = std(Ni{s}(rngT{p} & GoodFits{s}),1,2);
        Tes{p}(s)    = mean(Te{s}(rngT{p} & GoodFits{s}));
        dTes{p}(s)   = std(Te{s}(rngT{p} & GoodFits{s}),1,2);
    end
end

rOffset = 6.5;
[a1,b1] = sort(R);
[a0,b0] = sort(R)

figure; 
subplot(1,2,1); 
hold on
for p = 2:3
y(p) = errorbar(R(b1)-rOffset,nes{p}(b1),dnes{p}(b1),'LineWidth',1,'color','r','marker','o');
end
set(y(2),'color','k');
% hL{1} = legend([y(3),y(2)],'28 GHz ON','28 GHz OFF');
% set(hL{1},'location','north');
axis('square')
ylim([0,1]*1e19); ylabel('$n_e$ $[m^{-3}]$','Interpreter','Latex','fontsize',13)
xlim([-4,4]); xlabel('$R$ $[cm]$','Interpreter','Latex')
box on
grid on
t = title('$n_e$ $[m^{-3}]$ '); t.Interpreter = 'latex'; t.FontSize = 13;

subplot(1,2,2); hold on
for p = 2:3
k(p) = errorbar(R(b1)-rOffset,Tes{p}(b1),dTes{p}(b1),'LineWidth',1,'color','r','marker','o');
end
set(k(p),'color','k')
hL{2} = legend([k(3),k(2)],'28 GHz ON','28 GHz OFF');
set(hL{2},'location','north');
axis('square')
ylim([0,10]); xlabel('$R$ $[cm]$','Interpreter','Latex')
xlim([-4,4]); ylabel('$T_e$ $[eV]$','Interpreter','Latex','fontsize',13)
box on
grid on
t = title('$T_e$ $[eV]$ '); t.Interpreter = 'latex'; t.FontSize = 13;

set(gcf,'color','w')
set(gcf,'position',[360   174   421   444]);
