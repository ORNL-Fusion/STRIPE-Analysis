% Created 2017_08_13 by JF CANESES

% DLP data for Helicon + ECH
% Based on conditions used on the ECH campaign April 2017
% TR1: 600 A, TR2: 160 A, PS1/PS2: 4000 A
% U fueling, 1.21v:3.85s, 0.43v:4.1s, 0v:4.55s
% 115 kW RF at 13.56 MHz, 300 ms duration

% For this script we are only interested in the DLP 10.5 to show that we
% still have a high density plasma even though the light emission has
% dropped considerable and thus showing that we have a very highly ionized
% plasma column.

close all
clear all

CMPT = 1;
if CMPT == 1;
shotlist        = [13971    ,    17566,    13958];
DLPType         = {'10MP'   ,   '10MP',      '6'};
ChannelType     = {'1'      ,    '1'  ,      '1'};
SweepType       = {'niso'   ,   ' iso',   'niso'};
AttType         = {'Vx2,Ix5','Vx2,Ix5','Vx2,Ix5'};
Center_V        = [1        ,       1 ,        1]; % Remove offset on V: 1 (yes) 0(no)
Center_I        = [0        ,       0 ,        0]; % Remove offset on I: 1 (yes) 0(no)
tStart          = [4.15     ,     4.10,     4.15]; % [s]
tEnd            = [4.45     ,     4.7 ,     4.45];

% #########################################################################
% DLP FITTING ROUTINE
% #########################################################################
for s = 1:length(shotlist)
% #########################################################################
switch DLPType{s} 
    case '4'
        DLP{s} = 4.5;
        Config.L_tip = 1.2/1000; % 1.0 as of April 11th 2017
        Config.D_tip = 2*0.254/1000; % [m]
    case '9'
        DLP{s} = 9.5;
    case '6'
        DLP{s} = 6.5;
        Config.L_tip = 1.0/1000; % 1.0 as of April 11th 2017
        Config.D_tip = 0.254/1000; % [m]
    case '10'
        DLP{s} = 10.5;
        Config.L_tip = 1.2/1000;
        Config.D_tip = 0.254/1000; % [m]
    case '10MP'
        DLP{s} = 10.5;
        Config.L_tip = 1.8/1000;
        Config.D_tip = 0.254/1000; % [m]
    case '4MP'
        DLP{s} = 4.5;
        Config.L_tip = 1.7/1000;
        Config.D_tip = 0.254/1000; % [m]
    case '1MP'
        DLP{s} = 1.5;
        Config.L_tip = 1.7/1000;
        Config.D_tip = 0.254/1000; % [m]
end
% #########################################################################
switch AttType{s}
    case 'Vx2,Ix5'
        Config.V_Att = 2;  % Output voltage of DLP box (Voltage) = V_Att*Digitized data 
        Config.I_Att = 5;  % Output voltage of DLP box (Current) = I_att*Digitized data
    case 'Vx1,Ix1'
        Config.V_Att = 1;  % Output voltage of DLP box (Voltage) = V_Att*Digitized data 
        Config.I_Att = 1;  % Output voltage of DLP box (Current) = I_att*Digitized data
end
% #########################################################################
Stem = '\MPEX::TOP.'; Branch = 'MACHOPS1:'; RootAddress = [Stem,Branch];
switch ChannelType{s}
    case '1'
        DataAddress{1} = [RootAddress,'LP_V_RAMP']; % V
        DataAddress{2} = [RootAddress,'TARGET_LP']; % I
    case '2'
        DataAddress{1} = [RootAddress,'INT_4MM_1']; % V
        DataAddress{2} = [RootAddress,'INT_4MM_2']; % I
    case '3'
        DataAddress{1} = [RootAddress,'GEN_RF_PWR']; % V
        DataAddress{2} = [RootAddress,'INT_4MM_2'];  % I
    case '4'
        DataAddress{1} = [RootAddress,'ICH_LP']; % V
        DataAddress{2} = [RootAddress,'INT_4MM_2'];  % I
end
% #########################################################################
switch SweepType{s}
    case 'iso'
Config.V_cal = [(0.46e-3)^-1,0]; % Voltage output of DLP = V_cal(1)*Output voltage of DLP box + V_cal(2) %  2.1739e+03
Config.I_cal = [-1,0]; % Current output of DLP = (I_cal(2) + Output voltage of DLP box)/Ical(1)
    case 'niso'
Config.V_cal = [12.05,0.205];   % Voltage output of DLP = V_cal(1)*Output voltage of DLP box + V_cal(2)
Config.I_cal = [-142.5, 0.015]; % Current output of DLP = (I_cal(2) + Output voltage of DLP box)/Ical(1)
end
% #########################################################################
Config.Center_V = Center_V(s);
Config.Center_I = Center_I(s);
Config.tStart   = tStart(s);
Config.tEnd     = tEnd(s);
Config.FilterDataInput = 1; % Filter input data with savitsky Golay filter order 3 frame size 7
Config.SGF = 7;
Config.TimeMode = 2; % Effective time of sweep, (1) start of ramp or (2) mean time of ramp
Config.AMU = 2; % Ion mass in AMU
Config.FitFunction = 2; 
Config.AreaType = 1; % 1: Cylindrical + cap
% #########################################################################
[y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11] = DLP_fit_V5_2(Config,shotlist(s),DataAddress);
ni{s}{1}      = y1{1}{1}; ni{s}{2} = y1{1}{2};
Ni{s}         = 0.5*(ni{s}{1} + ni{s}{2});
Te{s}         = y2{1};
time{s}       = y3{1};
Ifit{s}       = y4{1};
Ip{s}         = y5{1};
Vp{s}         = y6{1};
tm{s}         = y7{1};
Vsweep{s}     = y8{1};
Isweep{s}     = y9{1};
GlitchFlag{s} = y10{1};
SSQres{s}     = y11{1};
end
%     save('HeliconECH_71_66_DLPs.mat')
else
%     load('HeliconECH_71_66_DLPs.mat')    
end
%% Test data
close all
C = {'k','b','g','m','c'}; ssqLevel = 0.001;
figure;
subplot(2,1,1); hold on
for s = 1:length(shotlist)
    rng = find(SSQres{s}<ssqLevel);
    rng = find(GlitchFlag{s}==0);
    h(s) = plot(time{s}(rng),Ni{s}(rng),C{s});
end
legend(h,DLPType); ylim([0,7]*1e19); xlim([4.15,4.45])

subplot(2,1,2); hold on
for s = 1:length(shotlist)
    rng = find(time{s}>=4.2 & time{s}<=4.44);
    h(s) = plot(time{s}(rng),Te{s}(rng),C{s});
end
ylim([0,6]); xlim([4.15,4.45]); legend(h,DLPType)
%% Load Color camera and Gas data
% the variables that we need to load from the camera and Gas data are:
load('HeliconECH_66_62_Dataset_1.mat')
%% Plot data:
% Define geometry of plot
close all
TimePlotStart = 4;
TimePlotEnd   = 4.6;
td1 = 4.23;
td2 = 4.43;
%==========================================================================
figure;
f = gcf;
% Define the window white in color and with no menubar. the normalized
% units are relative to the computer screen
set(f,'Menubar','figure','color','w','Units','normalized')
% Set the figure to cover the Right half of the screen
set(f,'Position',[0.4 0.05 0.4 0.65]); % [Left Bottom Width height]
set(f,'Position',[0.6 0.1 0.3 0.65]); % [Left Bottom Width height]
set(f,'Position',[0.6 0.1 0.3 0.5]); % [Left Bottom Width height]
set(f,'Position',[0.6 0.1 0.3 0.5]); % [Left Bottom Width height]

for s = 1:3
    ax(s) = axes;
    set(ax(s),'Units','Normalized','box','on')
    xlim(ax(s),[TimePlotStart,TimePlotEnd])
end
% The default position of any axes in normalized units is given by:
% [0.1300    0.1100    0.7750    0.8150]
% In order allow space for labels we need at least a 0.13 and 0.11 gap
% respectively
SmallOffset = 0.01;
dx = 0.13; 
dy = 0.15;
w = (1-(2*dx))/1;
h = 1.1*((1-(2*dy))/length(ax));
set(ax(1),'Position',[dx (dy + 0*(h-SmallOffset)) w h])
for s = 1:(length(ax)-1)
set(ax(s+1),'Position',[dx (dy + s*(h-SmallOffset)) w h],'XTick',[])
end
%==========================================================================

% RF
n = 3;
ax(n).NextPlot = 'add';
ax(n).YTick = [30:30:120];
ax(n).FontName = 'Times';
ax(n).YLabel.String = '[kW]';
RFsq = RF{1}.^2;
RFnormalized = RFsq./max(RFsq); 
RF_kW = RFnormalized*110; % kW
hRF(1) = plot(trf{1}(1:end-1),RF_kW,'Parent',ax(n),'color','k');
ylim(ax(n),[0,120])
ax(n).YGrid = 'on';
t1 = text(4.01,100,'RF power'); t1.FontName = 'Times';
hRF(1).LineWidth = 1.5;
line(ax(n),[td1,td1],[0,120],'color','r ','linewidth',2,'LineStyle',':')
line(ax(n),[td2,td2],[0,120],'color','r','linewidth',2,'LineStyle',':')

% DLP 6.5, Ne
n = 2;
ax(n).NextPlot = 'add';
ax(n).YTick = [1:1:3];
ax(n).FontName = 'Times';
ylabel(ax(n),'$[m^{-3}]$','FontSize',13); 
% text(ax(n),0.75,0.20,'$n_e$ ${\times}10^{-19}$','Units','normalized','Interpreter','latex')
% text(ax(n),0.55,0.50,'$T_e$','Units','normalized','Interpreter','latex','color','bl')

ax(n).YLabel.Interpreter = 'latex';
ax(n).YLabel.Rotation = 90; ax(n).YLabel.FontSize = 11; 
ax(n).YLabel.HorizontalAlignment = 'center'; ax(n).YLabel.VerticalAlignment = 'bottom';
% s = 2; rng = find(GlitchFlag{s}==0); nskip = 10;
% hne(s) = plot(time{s}(rng),Ni{s}(rng)*1e-19,'Parent',ax(n),'color','k')
s = 3; rng = find(GlitchFlag{s}==0); nskip = 10;
% hne(s) = plot(time{s}(rng),Ni{s}(rng)*1e-19,'Parent',ax(n),'color','k')
% hne(s).LineWidth = 1.5;
s = 3; hne6(s) = plot(time{s}(rng),Ni{s}(rng)*1e-19,'Parent',ax(n),'color','k')
rng = 12:118;
hte6(s) = plot(time{s}(rng),Te{s}(rng),'Parent',ax(n),'color','bl');

hne6(s).LineWidth = 1.5;

Lne = legend([hne6(s),hte6(s)],'$n_e{\times}10^{-19}$','$T_e$ $[eV]$');
Lne.Interpreter = 'latex';
Lne.Box = 'off'; Lne.Location = 'East';
Lne.Position = [0.2    0.4    0.07    0.08];

ylim(ax(n),[0,4])
ax(n).YGrid = 'on';
text(ax(n),0.02,0.8,'Probe C','Units','normalized','Interpreter','latex','FontSize',10)
line(ax(n),[td1,td1],[0,120],'color','r ','linewidth',2,'LineStyle',':')
line(ax(n),[td2,td2],[0,120],'color','r','linewidth',2,'LineStyle',':')


% PG spool 6.5
n = 1;
ax(n).NextPlot = 'add';
ax(n).YTick = [0.1,0.2,0.3];
ax(n).FontName = 'Times';
ylabel(ax(n),'[Pa]','FontSize',11); %
ax(n).YLabel.Interpreter = 'latex';
ax(n).YLabel.Rotation = 90; 
ax(n).YLabel.HorizontalAlignment = 'center'; ax(n).YLabel.VerticalAlignment = 'bottom';
hpgB(1) = plot(tpg2{1}(1:end-1)-0.03,PG3{1}/7.5,'Parent',ax(n),'color','k')
% hpgB(2) = plot(tpg1{1}(1:end-1),PG1{1}/7.5,'Parent',ax(n),'color','k')
hpgB(1).LineWidth = 1.5;
ylim(ax(n),[-0.005,0.35])
ax(n).YGrid = 'on';
text(ax(n),0.02,0.8,'$P3$','Units','normalized','Interpreter','latex','FontSize',10)
line(ax(n),[td1,td1],[0,120],'color','r ','linewidth',2,'LineStyle',':')
line(ax(n),[td2,td2],[0,120],'color','r','linewidth',2,'LineStyle',':')

n = 1;
xlabel(ax(n),'$time$ [s]','Interpreter','latex');
ax(n).XTick = [4:0.1:4.6];


abc = {'(c)','(b)','(a)'};
for s = 1:length(ax)
    ax(s).YLabel.FontSize = 10;
    text(ax(s),1.01,0.85,abc{s},'Units','normalized','Interpreter','latex','FontSize',11)
end
    