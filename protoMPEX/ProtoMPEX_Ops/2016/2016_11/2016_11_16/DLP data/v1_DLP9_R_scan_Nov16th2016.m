% list of shot number and associate timing

close all
clear all

series = 1; 

switch series 
    case 1 % Load data
        load('NiTe_Data_DLP_9_5_v1')
        CMPT = 0;
    case 2 % see data
        CMPT = 1;
    shotlist = 11000 + 400 +[103,104,106];
    tStart = 4 +            [15 ,15 ,15]/100;
    tEnd   = 4 +            [32 ,32 ,32]/100;     
    I_Att     =             [ 2 , 2 , 2];
    case 3 % compute data
        CMPT = 1;
        D = importdata('DLP_9_5_table.xlsx',',',2);
        % tstart: 4, tend:5, shot: 1, y: 2, Iatt: 6
        shotlist = D.data(:,1);
        tStart   = D.data(:,4)/100 + 4;
        tEnd     = D.data(:,5)/100 + 4;
        Rad      = D.data(:,2);
        I_Att    = D.data(:,6);
end
 
if CMPT == 1
% Address:
Stem = '\MPEX::TOP.';
Branch = 'MACHOPS1:';
RootAddress = [Stem,Branch];
DataAddress{1} = [RootAddress,'LP_V_RAMP']; % V
DataAddress{2} = [RootAddress,'TARGET_LP']; % I

% % see raw data:
% [Vs,t_vs]   = my_mdsvalue_v2(shotlist,DataAddress(1));
% [Is,t_Is]   = my_mdsvalue_v2(shotlist,DataAddress(2));

Config.L_tip = 1.2/1000;

Config.FilterDataInput = 1; % Filter input data with savitsky Golay filter order 3 
Config.SGF = 5; % frame size for the SG filter
Config.TimeMode = 2; % Effective time of sweep, (1) start of ramp or (2) mean time of ramp
Config.AMU = 2; % Ion mass in AMU
Config.D_tip = 0.254/1000; % [m]
Config.FitFunction = 2; 
Config.AreaType = 1; % 1: Cylindrical + cap
Config.V_Att = 2;  % Output voltage of DLP box (Voltage) = V_Att*Digitized data 
Config.V_cal = [12.05,0.205];   % Voltage output of DLP = V_cal(1)*Output voltage of DLP box + V_cal(2)
Config.I_cal = [-142.5, 0.015]; % Current output of DLP = (I_cal(2) + Output voltage of DLP box)/Ical(1)

for s = 1:length(shotlist)
    Config.I_Att  = I_Att(s);  % Output voltage of DLP box (Current) = I_att*Digitized data
    Config.tStart = tStart(s); % [s]
    Config.tEnd   = tEnd(s);
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
if 1
    save('NiTe_Data_DLP_9_5_v1')
end
end

%%
close all
C = {'k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m'};
t = t_zero(shotlist);

figure; 
subplot(2,2,1); hold on
for s = 1:length(shotlist)
    plot(time{s},ni{s}{1},C{s},'lineWidth',1);
    plot(time{s},ni{s}{2},C{s},'lineWidth',1);
    h(s) = plot(time{s},Ni{s},C{s},'lineWidth',2);
end
title('$ n_e $ $ [m^{-3}] $','interpreter','Latex','FontSize',13,'Rotation',0)
% legend(h,['DLP ',num2str(DLP)],'location','SouthEast')
%set(gca,'PlotBoxAspectRatio',[1 1 1])
ylim([0,8e19])
xlim([4.15,4.32])

subplot(2,2,2); hold on
for s = 1:length(shotlist)
    h(s) = plot(time{s},Te{s},C{s},'lineWidth',2)
    L{s} = [num2str(shotlist(s)),' ,t=',num2str(t{s}(10:14))];
end
legend(h,L,'location','NorthEast')
title('$ T_e $ $ [eV] $','interpreter','Latex','FontSize',13,'Rotation',0)
ylim([0,12])
xlim([4.15,4.32])

subplot(2,2,3); hold on
for s = 1:length(shotlist)
    plot(time{s},e_c.*Ni{s}.*Te{s},C{s},'lineWidth',2)
end
ylim([0,40])
xlim([4.15,4.32])
title('$ P_e $ $ [Pa] $','interpreter','Latex','FontSize',13,'Rotation',0)

set(findobj('-Property','YTick'),'box','on')
set(gcf,'color','w')

%%
if 1
    figure; hold on
    for s = 1:length(shotlist)
    plot(tm{s},Vp{s})
    h(s) = plot(tm{s},Ip{s}*1000);
    ylim([-100,100])
    xlim([4.15,4.32])
    grid on
    end
    legend(h,num2str(shotlist'))
end

if 1
    
     figure;
    for c = 1:25;
        subplot(5,5,c); hold on
        plot(Vsweep{s}{c},Isweep{s}{c}*1e3,'k')
        plot(Vsweep{s}{c},Ifit{s}{c}*1e3,'r')
        ht = title(['T_e: ',num2str(Te{s}(c))]);
        set(ht,'FontSize',5); grid on
    end
    
    figure;
        for c = 26:37;
        subplot(5,5,c-25); hold on
        plot(Vsweep{s}{c},Isweep{s}{c}*1e3,'k')
        plot(Vsweep{s}{c},Ifit{s}{c}*1e3,'r')
        ht = title(['T_e: ',num2str(Te{s}(c))]);
        set(ht,'FontSize',5); grid on
        end
           
end

%% extract steady state data
if 1
close all

for s = 1:length(shotlist)
    
 rng = find(time{s}>=4.28 & time{s}<=4.305);
 Nea(s) = mean(Ni{s}(rng));
 dNea(s) = std(Ni{s}(rng),1);
 Tea(s) = mean(Te{s}(rng));
 dTea(s) = std(Te{s}(rng),1);
 
end

[n,m] = sort(Rad);
PlotShotNum = 1;

figure;
subplot(2,1,1)
errorbar(Rad(m),Nea(m),dNea(m),'k.-')
grid on
if PlotShotNum
for s = 1:length(shotlist); text(Rad(s),Nea(s) + 0.1e19,num2str(shotlist(s)-11500),'FontSize',5); end
end
text(-6.3,7.7e19,'t = 4.29 to 4.305 sec')
xlabel('R (Probe) [cm]')
ylabel('n_e [m^{-3}]')
title('DLP 9.5')
ylim([0,10e19])

subplot(2,1,2)
errorbar(Rad(m),Tea(m),dTea(m),'k.-')
grid on
if PlotShotNum
for s = 1:length(shotlist); text(Rad(s),Tea(s) + 0.1,num2str(shotlist(s)-11500),'FontSize',5); end
end
ylim([0,5])
xlabel('R (Probe) [cm]')
ylabel('T_e [eV]')
title('DLP 6.5')

set(gcf,'color','w')
end

%% create excel spreadsheet


if 0
D = [(Rad(m)+2.84)',Nea(m)',dNea(m)',Tea(m)',dTea(m)'];
F = {'R [cm]','n_e [m^-3]','dNe [m^-3]','T_e [eV]','dTe [eV]'};
FileName = 'NeTe_Spool_9_5_Nov16th2016.xlsx';
xlswrite(FileName,[F;num2cell(D)]);
end
