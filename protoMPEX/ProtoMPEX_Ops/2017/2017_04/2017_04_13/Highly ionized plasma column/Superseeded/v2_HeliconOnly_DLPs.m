% DLP data for Helicon only
% Based on conditions used on the ECH campaign April 2017
% TR1: 600 A, TR2: 160 A, PS1/PS2: 4000 A
% U fueling, 1.21v:3.85s, 0.43v:4.1s, 0v:4.55s
% 115 kW RF at 13.56 MHz, 300 ms duration

% We first Compute/Load the DLP data from various locations on-axis
% Then, we load the visible camera and Gas data

close all
clear all

CMPT = 0;
if CMPT == 1;
shotlist        = [14265    ,14133    ,13958    , 13971   ];
DLPType         = {'1MP'    ,'4MP'    ,'6'      ,'10MP'   };
ChannelType     = {'4'      ,'2'      ,'1'      ,'1'      };
SweepType       = {'iso'    ,'iso'    ,'niso'   ,'niso'   };
AttType         = {'Vx1,Ix1','Vx1,Ix1','Vx2,Ix5','Vx2,Ix5'};
Center_V        = [1        ,1        ,1        ,1        ]; % Remove offset on V: 1 (yes) 0(no)
Center_I        = [1        ,1        ,0        ,0        ]; % Remove offset on I: 1 (yes) 0(no)
tStart          = [4.15     ,4.15     ,4.15     ,4.15     ]; % [s]
tEnd            = [4.45     ,4.45     ,4.45     ,4.45     ];

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
    %save('HeliconOnly_DLPs.mat')
else
    load('HeliconOnly_DLPs.mat')    
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
load('HeliconOnly_58_62_Dataset_1.mat')
%% Plot data: Ne and Te
% Define geometry of plot
close all
TimePlotStart = 4;
TimePlotEnd   = 4.6;
%==========================================================================
figure;
f = gcf;
% Define the window white in color and with no menubar. the normalized
% units are relative to the computer screen
set(f,'Menubar','figure','color','w','Units','normalized')
% Set the figure to cover the Right half of the screen
set(f,'Position',[0.4 0.05 0.4 0.65]); % [Left Bottom Width height]
set(f,'Position',[0.2 0.1 0.3 0.65]); % [Left Bottom Width height]
set(f,'Position',[0.2 0.1 0.3 0.5]); % [Left Bottom Width height]

for s = 1:5
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
dy = 0.11;
w = (1-(2*dx))/1;
h = 1.1*((1-(2*dy))/length(ax));
set(ax(1),'Position',[dx (dy + 0*(h-SmallOffset)) w h])
for s = 1:(length(ax)-1)
set(ax(s+1),'Position',[dx (dy + s*(h-SmallOffset)) w h],'XTick',[])
end
%==========================================================================

% RF
n = 5;
ax(n).NextPlot = 'add';
ax(n).YTick = [0.5,1];
ax(n).FontName = 'Times';
RFsq = RF{1}.^2;
RFnormalized = RFsq./max(RFsq);
trfnorm = trf{1}(1:end-1);
h(1) = plot(trfnorm,RFnormalized,'Parent',ax(n),'color','k')
ylim(ax(n),[0,1])
% Light Emission from central chamber
n = 5;
ax(n).NextPlot = 'add';
ax(n).YTick = [0.2:0.2:0.8];
ax(n).FontName = 'Times';
ylabel(ax(n),'[a.u]','FontSize',11);
ax(n).YLabel.Interpreter = 'latex';
ax(n).YLabel.Rotation = 90; 
ax(n).YLabel.HorizontalAlignment = 'center'; ax(n).YLabel.VerticalAlignment = 'bottom';
%h(2) = plot(t*1e-3,OnAxisEmission./max(OnAxisEmission),'Parent',ax(n),'color','g')
ylim(ax(n),[0,1])
ax(n).YGrid = 'on';
L5 = legend(h,'$P_{RF}$'); 
L5.Interpreter = 'latex';
L5.Box = 'off'; L5.Location = 'northwest';

% DLP 1.5
n = 4;
ax(n).NextPlot = 'add';
ax(n).YTick = [2:2:6];
ax(n).FontName = 'Times';
ylabel(ax(n),'$n_e$ , $T_e$','FontSize',13);
ax(n).YLabel.Interpreter = 'latex';
ax(n).YLabel.Rotation = 90; 
ax(n).YLabel.HorizontalAlignment = 'center'; ax(n).YLabel.VerticalAlignment = 'bottom';
s = 1; rng = find(GlitchFlag{s}==0); nskip = 10;
hne(1) = plot(time{s}(rng),Ni{s}(rng)*1e-19,'Parent',ax(n),'color','k')
hne(2) = plot(time{s}(rng(nskip:end)),Te{s}(rng(nskip:end)),'Parent',ax(n),'color','r')
ylim(ax(n),[0,8])
Lne = legend(hne,'$n_e$','$T_e$');
Lne.Interpreter = 'latex';
Lne.Box = 'off'; Lne.Location = 'Southwest';
ax(n).YGrid = 'on';
text(ax(n),0.02,0.8,'Probe A','Units','normalized','Interpreter','latex','FontSize',10)

% DLP 4.5
n = 3;
ax(n).NextPlot = 'add';
ax(n).YTick = [2:2:6];
ax(n).FontName = 'Times';
ylabel(ax(n),'$n_e$ , $T_e$','FontSize',13);
ax(n).YLabel.Interpreter = 'latex';
ax(n).YLabel.Rotation = 90; 
ax(n).YLabel.HorizontalAlignment = 'center'; ax(n).YLabel.VerticalAlignment = 'bottom';
s = 2; rng = find(GlitchFlag{s}==0); nskip = 11;
ne4 = Ni{s}(rng);
tne4 = time{s}(rng);
hne(1) = plot(tne4,ne4*1e-19,'Parent',ax(n),'color','k')
hne(2) = plot(time{s}(rng(nskip:end)),Te{s}(rng(nskip:end)),'Parent',ax(n),'color','r')
ylim(ax(n),[0,8])
Lne = legend(hne,'$n_e$','$T_e$'); % \times 10^{-19}
Lne.Interpreter = 'latex';
Lne.Box = 'off'; Lne.Location = 'Southwest';
ax(n).YGrid = 'on';
text(ax(n),0.02,0.8,'Probe B','Units','normalized','Interpreter','latex','FontSize',10)

% DLP 6.5
n = 2;
ax(n).NextPlot = 'add';
ax(n).YTick = [2:2:6];
ax(n).FontName = 'Times';
ylabel(ax(n),'$n_e$ , $T_e$','FontSize',13);
ax(n).YLabel.Interpreter = 'latex';
ax(n).YLabel.Rotation = 90; 
ax(n).YLabel.HorizontalAlignment = 'center'; ax(n).YLabel.VerticalAlignment = 'bottom';
s = 3; rng = find(GlitchFlag{s}==0); nskip = 10;
hne(1) = plot(time{s}(rng),Ni{s}(rng)*1e-19,'Parent',ax(n),'color','k')
hne(2) = plot(time{s}(rng(nskip:end-1)),Te{s}(rng(nskip:end-1)),'Parent',ax(n),'color','r')
ylim(ax(n),[0,8])
Lne = legend(hne,'$n_e$','$T_e$');
Lne.Interpreter = 'latex';
Lne.Box = 'off'; Lne.Location = 'Southwest';
ax(n).YGrid = 'on';
text(ax(n),0.02,0.8,'Probe C','Units','normalized','Interpreter','latex','FontSize',10)

% DLP 10.5
n = 1;
ax(n).NextPlot = 'add';
ax(n).YTick = [2:2:6];
ax(n).FontName = 'Times';
ylabel(ax(n),'$n_e$ , $T_e$','FontSize',13);
xlabel(ax(n),'$time$ [s]','FontSize',11);
ax(n).YLabel.Interpreter = 'latex';ax(n).XLabel.Interpreter = 'latex';
ax(n).YLabel.Rotation = 90; 
ax(n).YLabel.HorizontalAlignment = 'center'; ax(n).YLabel.VerticalAlignment = 'bottom';
s = 4; rng = find(GlitchFlag{s}==0); nskip = 10;
ne10 = Ni{s}(rng);
tne10 = time{s}(rng);
hne(1) = plot(tne10,ne10*1e-19,'Parent',ax(n),'color','k')
hne(2) = plot(time{s}(rng(nskip:end-2)),Te{s}(rng(nskip:end-2)),'Parent',ax(n),'color','r')
ylim(ax(n),[0,8])
Lne = legend(hne,'$n_e$','$T_e$');
Lne.Interpreter = 'latex';
Lne.Box = 'off'; Lne.Location = 'Southwest';
ax(n).YGrid = 'on';
text(ax(n),0.02,0.8,'Probe D','Units','normalized','Interpreter','latex','FontSize',10)

abc = {'(e)','(d)','(c)','(b)','(a)'};
for s = 1:length(ax)
    ax(s).YLabel.FontSize = 11;
    text(ax(s),1.01,0.85,abc{s},'Units','normalized','Interpreter','latex','FontSize',11)
end

if 0
save('DLP_4_10_2017_04_13','RFnormalized','trfnorm','ne4','tne4','ne10','tne10')
end

%% Plot data: Gas data
% Define geometry of plot
%close all
TimePlotStart = 4;
TimePlotEnd   = 4.6;
%==========================================================================
figure
f = gcf;
% Define the window white in color and with no menubar. the normalized
% units are relative to the computer screen
set(f,'Menubar','figure','color','w','Units','normalized')
% Set the figure to cover the Right half of the screen
set(f,'Position',[0.4 0.05 0.4 0.65]); % [Left Bottom Width height]
set(f,'Position',[0.6 0.1 0.3 0.65]); % [Left Bottom Width height]
set(f,'Position',[0.6 0.1 0.3 0.5]); % [Left Bottom Width height]

for s = 1:5
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
dy = 0.11;
w = (1-(2*dx))/1;
h = 1.1*((1-(2*dy))/length(ax));
set(ax(1),'Position',[dx (dy + 0*(h-SmallOffset)) w h])
for s = 1:(length(ax)-1)
set(ax(s+1),'Position',[dx (dy + s*(h-SmallOffset)) w h],'XTick',[])
end
%==========================================================================
% RF
n = 5;
ax(n).NextPlot = 'add';
ax(n).YTick = [0.5,1];
ax(n).FontName = 'Times';
RFsq = RF{1}.^2;
RFnormalized = RFsq./max(RFsq);
trfnorm = trf{1}(1:end-1);
h(1) = plot(trfnorm,RFnormalized,'Parent',ax(n),'color','k')
ylim(ax(n),[0,1])
% Light Emission from central chamber
n = 5;
ax(n).NextPlot = 'add';
ax(n).YTick = [0.5,1];
ax(n).FontName = 'Times';
ylabel(ax(n),'[a.u]','FontSize',11);
ax(n).YLabel.Interpreter = 'latex';
ax(n).YLabel.Rotation = 90; 
ax(n).YLabel.HorizontalAlignment = 'center'; ax(n).YLabel.VerticalAlignment = 'bottom';
%h(2) = plot(t*1e-3,OnAxisEmission./max(OnAxisEmission),'Parent',ax(n),'color','g')
ylim(ax(n),[0,1])
ax(n).YGrid = 'on';
L5 = legend(h(1),'$P_{RF}$'); 
L5.Interpreter = 'latex';
L5.Box = 'off'; L5.Location = 'northwest';

% PG spool 2.5
n = 4;
ax(n).NextPlot = 'add';
ax(n).YTick = [1:1:3];
ax(n).FontName = 'Times';
ylabel(ax(n),'[Pa]','FontSize',11);
ax(n).YLabel.Interpreter = 'latex';
ax(n).YLabel.Rotation = 90; 
ax(n).YLabel.HorizontalAlignment = 'center'; ax(n).YLabel.VerticalAlignment = 'bottom';
hpgA(1) = plot(tpg2{1}(1:end-1),PG2{1}/7.5,'Parent',ax(n),'color','k')
hpgA(2) = plot(tpg2{1}(1:end-1),PG2_GasOnly{1}/7.5,'Parent',ax(n),'color','r')
ylim(ax(n),[-0.05,4])
ax(n).YGrid = 'on';
% LpgA = legend(hpgA,'w/ Plasma','no Plasma');
% LpgA.Interpreter = 'latex';
% LpgA.Box = 'off'; LpgA.Location = 'northwest';
text(ax(n),0.02,0.8,'$P1$','Units','normalized','Interpreter','latex','FontSize',10)

% PG spool 4.5
n = 3;
ax(n).NextPlot = 'add';
ax(n).YTick = [1:1:2];
ax(n).FontName = 'Times';
ylabel(ax(n),'[Pa]','FontSize',11);
ax(n).YLabel.Interpreter = 'latex';
ax(n).YLabel.Rotation = 90; 
ax(n).YLabel.HorizontalAlignment = 'center'; ax(n).YLabel.VerticalAlignment = 'bottom';
hpgA(1) = plot(tpg2{1}(1:end-1),PG4{1}/7.5,'Parent',ax(n),'color','k')
hpgA(2) = plot(tpg2{1}(1:end-1),PG4_GasOnly{1}/7.5,'Parent',ax(n),'color','r')
ylim(ax(n),[-0.05,3])
ax(n).YGrid = 'on';
% LpgA = legend(hpgA,'w/ Plasma','no Plasma');
% LpgA.Interpreter = 'latex';
% LpgA.Box = 'off'; LpgA.Location = 'northwest';
text(ax(n),0.02,0.8,'$P2$','Units','normalized','Interpreter','latex','FontSize',10)

% PG spool 6.5
n = 2;
ax(n).NextPlot = 'add';
ax(n).YTick = [0.2:0.2:0.8];
ax(n).FontName = 'Times';
ylabel(ax(n),'[Pa]','FontSize',11);
ax(n).YLabel.Interpreter = 'latex';
ax(n).YLabel.Rotation = 90; 
ax(n).YLabel.HorizontalAlignment = 'center'; ax(n).YLabel.VerticalAlignment = 'bottom';
hpgB(1) = plot(tpg2{1}(1:end-1),PG3{1}/7.5,'Parent',ax(n),'color','k')
hpgB(2) = plot(tpg2{1}(1:end-1),PG3_GasOnly{1}/7.5,'Parent',ax(n),'color','r')
ylim(ax(n),[-0.01,1])
ax(n).YGrid = 'on';
% LpgB = legend(hpgB,'w/ Plasma','no Plasma');
% LpgB.Interpreter = 'latex';
% LpgB.Box = 'off'; LpgB.Location = 'northwest';
text(ax(n),0.02,0.8,'$P3$','Units','normalized','Interpreter','latex','FontSize',10)
% text(ax(n),0.59,0.62,'${\times}0.5$','Units','normalized','Interpreter'...
%     ,'latex','FontSize',10,'color','r','BackGroundColor','w')


% PG spool 9.5
n = 1;
ax(n).NextPlot = 'add';
ax(n).YTick = [0.1,0.2,0.3];
ax(n).FontName = 'Times';
ylabel(ax(n),'[Pa]','FontSize',11);
xlabel(ax(n),'$time$ [s]','FontSize',12);
ax(n).YLabel.Interpreter = 'latex'; ax(n).XLabel.Interpreter = 'latex';
ax(n).YLabel.Rotation = 90; 
ax(n).YLabel.HorizontalAlignment = 'center'; ax(n).YLabel.VerticalAlignment = 'bottom';
hpgB(1) = plot(tpg1{1}(1:end-1),PG1{1}/7.5,'Parent',ax(n),'color','k')
hpgB(2) = plot(tpg1{1}(1:end-1),PG1_GasOnly{1}/7.5,'Parent',ax(n),'color','r')
ylim(ax(n),[0,0.4])
ax(n).YGrid = 'on';
LpgB = legend(hpgB,'$P_n$ , w/ Plasma','$P_n$ , no Plasma');
LpgB.Interpreter = 'latex';
LpgB.Box = 'off'; LpgB.Location = 'northwest';
LpgB.Position = [0.25    0.1881    0.2394    0.0712];
text(ax(n),0.02,0.8,'$P4$','Units','normalized','Interpreter','latex','FontSize',10)

abc = {'(e)','(d)','(c)','(b)','(a)'};
for s = 1:length(ax)
    ax(s).YLabel.FontSize = 11;
    text(ax(s),1.01,0.85,abc{s},'Units','normalized','Interpreter','latex','FontSize',11)
end
