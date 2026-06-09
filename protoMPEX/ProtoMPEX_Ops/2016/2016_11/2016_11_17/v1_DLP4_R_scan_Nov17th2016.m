% list of shot number and associate timing

close all
clear all

series = 2; 

% there exist a large difference in Te between shots 115,116 and 118
% According to the logbook, 115 and 116 are not considered in MJ and only
% 118 is considered MJ. that explains why shots 115,116 have much higher
% Te while 118 is much lower in Te

switch series 
    case 1 % Load data
        load('NiTe_Data_DLP_4_5_v1')
        CMPT = 0;
    case 2 % see data
        CMPT = 1;
    shotlist = 11000 + 500 +[92,100];
    tStart = 4 +            [15,27.5]/100;
    tEnd   = 4 +            [32,30]/100;     
    I_Att     =             [ 1,1];
    
    shotlist = 11000 + 500 +[131,132];
    tStart = 4 +            [15,25]/100;
    tEnd   = 4 +            [32,32]/100;     
    I_Att     =             [1,1];
    
        shotlist = 11000 + 500 +[101,102];
    tStart = 4 +            [20,20]/100;
    tEnd   = 4 +            [28.5,28.5]/100;     
    I_Att     =             [1,1];
    
    % look at the differences in Te:
            shotlist = 11000 + 500 +[115,116,118];
    tStart = 4 +            [15,15,15]/100;
    tEnd   = 4 +            [28,28,28]/100;     
    I_Att     =             [1,1];
    
%     % large dTe:
%     shotlist = 11000 + 500 +[97];
%     tStart = 4 +            [15]/100;
%     tEnd   = 4 +            [32]/100;     
%     I_Att     =             [1,1];
    
    rng = find(tStart ~= NaN);
    Rad = 1;
    
    case 3 % compute data
        CMPT = 1;
        D = importdata('DLP_4_5_table.xlsx',',',2);
        % tstart: 4, tend:5, shot: 1, y: 2, Iatt: 6
        shotlist = 11000 + 500 + D.data(:,1);
        tStart   = D.data(:,5)/100 + 4;
        CStart   = D.data(:,8)/100 + 4;
        tEnd     = D.data(:,6)/100 + 4;
        CEnd     = D.data(:,9)/100 + 4;        
        Rad      = D.data(:,3);
        MJ       = D.data(:,7);
        TR2      = D.data(:,2); 
       
        rng = find(tStart ~= NaN & MJ == 1);
end
 
if CMPT == 1
% Address:
Stem = '\MPEX::TOP.';
Branch = 'MACHOPS1:';
RootAddress = [Stem,Branch];
DataAddress{1} = [RootAddress,'INT_4MM_1']; % V
DataAddress{2} = [RootAddress,'INT_4MM_2']; % I

DAddress{1} = [Stem,'FSCOPE','TUBE8']; % I


% % see raw data:
% [Vs,t_vs]   = my_mdsvalue_v2(shotlist,DataAddress(1));
% [Is,t_Is]   = my_mdsvalue_v2(shotlist,DataAddress(2));

Config.L_tip = 1.2/1000;
Config.D_tip = 0.508/1000; % [m], Probe tip for DLP 4.5 is thicker than the other probes in ProtoMPEX

Config.FilterDataInput = 1; % Filter input data with savitsky Golay filter order 3 
Config.SGF = 5; % frame size for the SG filter
Config.TimeMode = 2; % Effective time of sweep, (1) start of ramp or (2) mean time of ramp
Config.AMU = 2; % Ion mass in AMU
Config.FitFunction = 2; 
Config.AreaType = 1; % 1: Cylindrical + cap
Config.V_Att = 1;  % Output voltage of DLP box (Voltage) = V_Att*Digitized data 
Config.V_cal = [(0.46e-3)^-1,0]; % Voltage output of DLP = V_cal(1)*Output voltage of DLP box + V_cal(2) %  2.1739e+03
Config.I_cal = [-1,0]; % Current output of DLP = (I_cal(2) + Output voltage of DLP box)/Ical(1)

for s = 1:length(shotlist(rng))
    Config.I_Att  = 1;  % Output voltage of DLP box (Current) = I_att*Digitized data
    Config.tStart = tStart(rng(s)); % [s]
    Config.tEnd   = tEnd(rng(s));
    [a,b,c,d,f,g,j,p,q] = DLP_fit_V5(Config,shotlist(rng(s)),DataAddress);
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
    save('NiTe_Data_DLP_4_5_v1')
end
end

%%
close all
C = {'k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m'};
t = t_zero(shotlist);

figure; 
subplot(2,2,1); hold on
for s = 1:length(shotlist(rng))
       if isnan(tStart(s))
        continue
        end
    %plot(time{s},ni{s}{1},C{s},'lineWidth',1);
    %plot(time{s},ni{s}{2},C{s},'lineWidth',1);
    h(s) = plot(time{s},Ni{s},C{s},'lineWidth',2);
    shotlist(rng(s))
end
title('$ n_e $ $ [m^{-3}] $','interpreter','Latex','FontSize',13,'Rotation',0)
% legend(h,['DLP ',num2str(DLP)],'location','SouthEast')
%set(gca,'PlotBoxAspectRatio',[1 1 1])
ylim([0,6e19])
xlim([4.15,4.32])

subplot(2,2,2); hold on
for s = 1:length(shotlist(rng))
    h(s) = plot(time{s},Te{s},C{s},'lineWidth',2);
    L{s} = [num2str(shotlist(rng(s))),' ,t=',num2str(t{s}(10:14))];
end
legend(h,L,'location','NorthEast')
title('$ T_e $ $ [eV] $','interpreter','Latex','FontSize',13,'Rotation',0)
ylim([0,12])
xlim([4.15,4.32])

subplot(2,2,3); hold on
for s = 1:length(shotlist(rng))
       if isnan(tStart(s))
        continue
        end
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
    for s = 1:length(shotlist(rng))
    plot(tm{s},Vp{s})
    h(s) = plot(tm{s},Ip{s}*1000);
    ylim([-100,100])
    xlim([4.15,4.32])
    grid on
    end
    legend(h,num2str(shotlist(rng)'))
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
if (series == 1 || series  == 3)
close all

for s = 1:length(shotlist(rng))
    
 rngTime = find(time{s}>=CStart(rng(s)) & time{s}<=CEnd(rng(s)));
 Nea(s) = mean(Ni{s}(rngTime));
 dNea(s) = std(Ni{s}(rngTime),1);
 Tea(s) = mean(Te{s}(rngTime));
 dTea(s) = std(Te{s}(rngTime),1);
 
end

[n,m] = sort(Rad(rng));
PlotShotNum = 0;
PlotTr2 = 1;

figure;
subplot(2,1,1)
errorbar(Rad(rng(m)),Nea(m),dNea(m),'k.-')
grid on
if PlotShotNum
for s = 1:length(shotlist(rng)); text(Rad(rng(s)),Nea(s) + 0.1e19,num2str(shotlist(rng(s))-11500),'FontSize',5); end
end
if PlotTr2
for s = 1:length(shotlist(rng)); text(Rad(rng(s)),Nea(s) + 0.3e19,num2str(TR2(rng(s))),'FontSize',5); end
end
xlabel('R (Probe) [cm]')
ylabel('n_e [m^{-3}]')
title('DLP 4.5')
ylim([0,6e19])
xlim([-8,8])

subplot(2,1,2)
errorbar(Rad(rng(m)),Tea(m),dTea(m),'k.-')
grid on
if PlotShotNum
for s = 1:length(shotlist(rng)); text(Rad(rng(s)),Tea(s) + 0.1,num2str(shotlist(rng(s))-11500),'FontSize',5); end
end
if PlotTr2
for s = 1:length(shotlist(rng)); text(Rad(rng(s)),Tea(s) + 0.6,num2str(TR2(rng(s))),'FontSize',5); end
end
ylim([0,10])
xlabel('R (Probe) [cm]')
ylabel('T_e [eV]')
title('DLP 4.5')
xlim([-8,8])

set(gcf,'color','w')
end

%% create excel spreadsheet


if 1
D = [shotlist(rng(m)),Rad(rng(m)),Nea(m)',dNea(m)',Tea(m)',dTea(m)',TR2(rng(m))];
F = {'Shot','R [cm]','n_e [m^-3]','dNe [m^-3]','T_e [eV]','dTe [eV]','TR2 [A]'};
FileName = 'NeTe_Spool_4_5_Nov17th2016.xlsx';
xlswrite(FileName,[F;num2cell(D)]);
end
