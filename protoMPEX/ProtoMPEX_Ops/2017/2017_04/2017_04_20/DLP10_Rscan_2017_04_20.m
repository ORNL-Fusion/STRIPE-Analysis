close all
clear all

CMPT = 0;
 
if CMPT == 1

shotlist = 14000 + [57,58,59 ,60 ,61,62,63 ,64 ,65,66,67 ,69 ,70,72,73  ,74  ,77,78,79  ,80  ,81,82,83  ,84  ];
r        =         [1 ,1 ,1.5,1.5,2 ,2 ,2.5,2.5,3 ,3 ,3.5,3.5,0 ,0 ,-0.5,-0.5,-1,-1,-1.5,-1.5,-2,-2,-2.5,-2.5];
echPulse =         [1 ,0 ,1  ,0  ,1 ,0 ,1  ,0  ,1 ,0 ,1  ,0  ,1 ,0 ,1   ,0   ,1 ,0 ,1   ,0   ,1 ,0 ,1   ,0   ];
DLPType = '10MP';

Config.tStart = 4.15; % [s]
Config.tEnd = 4.45;

AddressType='s'; % s for standard
CalType='niso'; % niso for not isolated, is for isolated

% Acquiring Ne and Te data
Stem = '\MPEX::TOP.';
Branch = 'MACHOPS1:';
RootAddress = [Stem,Branch];

switch AddressType
    case 's'
DataAddress{1} = [RootAddress,'LP_V_RAMP']; % V
DataAddress{2} = [RootAddress,'TARGET_LP']; % I
Config.V_Att = 2;  % Output voltage of DLP box (Voltage) = V_Att*Digitized data 
Config.I_Att = 5;  % Output voltage of DLP box (Current) = I_att*Digitized data
    case 'n'
DataAddress{1} = [RootAddress,'INT_4MM_1']; % V
DataAddress{2} = [RootAddress,'INT_4MM_2']; % I
Config.V_Att = 1;  % Output voltage of DLP box (Voltage) = V_Att*Digitized data 
Config.I_Att = 1;  % Output voltage of DLP box (Current) = I_att*Digitized data
end

DA{1} = [RootAddress,'PWR_28GHZ']; % Isx
[ECH,t_ech]   = my_mdsvalue_v2(shotlist,DA(1))

DA{1} = [RootAddress,'RF_FWD_PWR'];
[RF,t_rf]   = my_mdsvalue_v2(shotlist,DA(1))

switch CalType 
    case 'iso'
Config.V_cal = [(0.46e-3)^-1,0]; % Voltage output of DLP = V_cal(1)*Output voltage of DLP box + V_cal(2) %  2.1739e+03
Config.I_cal = [-1,0]; % Current output of DLP = (I_cal(2) + Output voltage of DLP box)/Ical(1)
    case 'niso'
Config.V_cal = [12.05,0.205];   % Voltage output of DLP = V_cal(1)*Output voltage of DLP box + V_cal(2)
Config.I_cal = [-142.5, 0.015]; % Current output of DLP = (I_cal(2) + Output voltage of DLP box)/Ical(1)
end

switch DLPType 
    case '4'
DLP = 4.5;
    case '9'
DLP = 9.5;
    case '6'
DLP = 6.5;
Config.L_tip = 1.0/1000; % 1.0 as of April 11th 2017
    case '10'
DLP = 10.5;
Config.L_tip = 1.2/1000;
    case '10MP'
DLP = 10.5;
Config.L_tip = 1.8/1000;
end

Config.FilterDataInput = 1; % Filter input data with savitsky Golay filter order 3 frame size 7
Config.SGF = 7;
Config.TimeMode = 2; % Effective time of sweep, (1) start of ramp or (2) mean time of ramp
Config.AMU = 2; % Ion mass in AMU
Config.D_tip = 0.254/1000; % [m]
Config.FitFunction = 2; 
Config.AreaType = 1; % 1: Cylindrical + cap

[ni,Te,time,Ifit,Ip,Vp,tm,Vsweep,Isweep] = DLP_fit_V5(Config,shotlist,DataAddress);

for s = 1:length(shotlist)
    Ni{s} = 0.5*(ni{s}{1} + ni{s}{2});
end
    save ('DLP10_RscanData.mat')
else
    load('DLP10_RscanData.mat')
end
% note, 


%%
close all
C = {'k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:'};
t = t_zero(shotlist);

TimePlotEnd = 4.48;

figure; 
k = 1;
for s = 1:12
    subplot(4,3,s); hold on
%     plot(time{s},ni{s}{1},C{s},'lineWidth',1);
%     plot(time{s},ni{s}{2},C{s},'lineWidth',1);

    % ECH ON case
    h(k) = plot(time{k},Ni{k},'k','lineWidth',2);
    plot(t_ech{k}(1:length(ECH{k})),ECH{k}*0.5e19,'k','lineWidth',0.5)
    plot(t_rf{k}(1:length(RF{k})),(RF{k}.^2)*6e19,'k','lineWidth',0.5)

    % ECH OFF case
    h(k+1) = plot(time{k+1},Ni{k+1},'r','lineWidth',2);
    
    title(['r = ',num2str(r(k))])
    ylim([0,5e19])
    xlim([4.15,TimePlotEnd])
    
    if r(k) == 0
        rng0 = find(time{s}>=4.35 & time{s}<=4.4);
    else
        rng0 = find(time{s}>=4.35 & time{s}<=4.4);
    end
    % ECH ON case
        ne1(s)  = mean(Ni{k}(rng0));
        dne1(s) = std(Ni{k}(rng0),1,2);
        Te1(s)  = mean(Te{k}(rng0));
        dTe1(s) = std(Te{k}(rng0),1,2);
        r1(s) = r(k);
        
    % ECH OFF case 
        Time_echoff{s} = time{k+1};
        Ni_echoff{s} = Ni{k+1};
        ne0(s)  = mean(Ni_echoff{s}(rng0));
        dne0(s) = std(Ni_echoff{s}(rng0),1,2);
        Te_echoff{s} = Te{k+1};
        Te0(s)  = mean(Te_echoff{s}(rng0));
        dTe0(s) = std(Te_echoff{s}(rng0),1,2);
        r0(s) = r(k+1);
    k = k + 2;
end
set(gcf,'color','w')
box on

figure; 
k = 1;
for s = 1:12
    subplot(4,3,s); hold on
    h(k) = plot(time{k},Te{k},'k','lineWidth',2);
    h(k+1) = plot(time{k+1},Te{k+1},'r','lineWidth',2);
    title(['r = ',num2str(r(k))])
    ylim([0,10])
    xlim([4.15,TimePlotEnd])
    k = k + 2;
end
set(gcf,'color','w')
box on

% subplot(2,2,3); hold on
% for s = 1:length(shotlist)
%     plot(time{s},e_c.*Ni{s}.*Te{s},C{s},'lineWidth',2)
% end
% ylim([0,30])
% xlim([4.15,TimePlotEnd])
% title('$ P_e $ $ [Pa] $','interpreter','Latex','FontSize',13,'Rotation',0)

set(findobj('-Property','YTick'),'box','on')
set(gcf,'color','w')

% GasCompare_2017_04_20

if 1
    figure; hold on
    for s = 1:length(shotlist)
    plot(tm{s},Vp{s})
    h(s) = plot(tm{s},Ip{s}*1000);
    ylim([-100,100])
    xlim([4.15,TimePlotEnd])
    grid on
    end
    legend(h,num2str(shotlist'))
    set(gcf,'color','w')
    box on
end

if 0
    
     figure;
    for c = 1:25;
        subplot(5,5,c); hold on
        plot(Vsweep{s}{c},Isweep{s}{c}*1e3,'k')
        plot(Vsweep{s}{c},Ifit{s}{c}*1e3,'r')
        ht = title(['T_e: ',num2str(Te{s}(c))]);
        set(ht,'FontSize',5); grid on
    end
    
    figure;
        for c = 26:50;
        subplot(5,5,c-25); hold on
        plot(Vsweep{s}{c},Isweep{s}{c}*1e3,'k')
        plot(Vsweep{s}{c},Ifit{s}{c}*1e3,'r')
        ht = title(['T_e: ',num2str(Te{s}(c))]);
        set(ht,'FontSize',5); grid on
        end
        
    figure;
        for c = 51:59;
        subplot(5,5,c-50); hold on
        plot(Vsweep{s}{c},Isweep{s}{c}*1e3,'k')
        plot(Vsweep{s}{c},Ifit{s}{c}*1e3,'r')
        ht = title(['T_e: ',num2str(Te{s}(c))]);
        set(ht,'FontSize',5); grid on
        end  
    
end

%%
[a,b] = sort(r0);

figure; hold on
rng = 1:118;
for s = 1:length(Te_echoff)
plot3(Time_echoff{b(s)}(rng),r0(b(s))*ones(size(Te_echoff{b(s)}(rng))),Te_echoff{b(s)}(rng))
Te2D_echoff(s,:) = Te_echoff{b(s)}(rng);
end
zlim([0,10])

figure; hold on
n = 10; e(1) = errorbar(r0(b),mean(Te2D_echoff(:,n:n+3),2),std(Te2D_echoff(:,n:n+3),1,2));
n = 110; e(2) = errorbar(r0(b),mean(Te2D_echoff(:,n:n+5),2),std(Te2D_echoff(:,n:n+5),1,2));
ylim([0,6])
set(e(1),'color','k')
set(e(2),'color','r')


%%
figure; 
offset = -0.0;
subplot(2,1,1); hold on
h = errorbar(r1(b)+offset,ne1(b),dne1(b),'ro-','LineWidth',1)
h = errorbar(r0(b)+offset,ne0(b),dne0(b),'ko-','LineWidth',1); 
ylim([0,5]*1e19)
xlim([-4,4])
box on
grid on
t = title('$n_e$ $[m^{-3}]$ '); t.Interpreter = 'latex'; t.FontSize = 13;

subplot(2,1,2); hold on
h = errorbar(r1(b)+offset,Te1(b),dTe1(b),'ro-','LineWidth',1)
h = errorbar(r0(b)+offset,Te0(b),dTe0(b),'ko-','LineWidth',1)
ylim([0,7])
xlim([-4,4])
box on
grid on
t = title('$T_e$ $[eV]$ '); t.Interpreter = 'latex'; t.FontSize = 13;

set(gcf,'color','w')
%% Figure for APS poster
close all
figure; 
hold on
h(1) = errorbar(r0(b)-0.5,ne0(b)*1e-19,dne0(b),'ko-','LineWidth',1); 
h(2) = errorbar(r0(b)-0.5,Te0(b),dTe0(b),'ro-','LineWidth',1)
set(gca,'YTick',[0:1:5])
ylim([0,5])
xlim([-4,4])
xlabel('$r$ $[cm]$','Interpreter','latex','FontSize',14)
box on
grid on
t = legend(h,{'$n_e$ $[\times 10^{19}$ $m^{-3}]$ ','$T_e$ $[eV]$ '}); 
t.Interpreter = 'latex'; t.FontSize = 13; t.Location = 'northwest';
% t = title('$T_e$ $[eV]$ '); t.Interpreter = 'latex'; t.FontSize = 13;

set(gcf,'color','w')