% Helicon only analsysis

clear all
close all
CMPT = 0;
if CMPT == 1 % 49 seconds
% Based on shot 13962 (probe out) and shot 13958
shotlist = [13962,13958];
% #########################################################################
% 1 - Gathering DLP data: Ne, Te
% #########################################################################
Config.tStart = 4.15; % [s]
Config.tEnd = 4.45;
SweepType = 'niso'; % "iso" for isolated sweep, "niso" for non-isolated sweep
AttType = 'Vx2,Ix5'; % Attenuation on the digitized signals
Config.Center_V = 1; % Remove offset on V: 1 (yes) 0(no)
Config.Center_I = 0; % Remove offset on I: 1 (yes) 0(no)
DLP = 6.5;
Config.L_tip = 1.0/1000; % 1.0 as of April 11th 2017
Config.D_tip = 0.254/1000; % [m]
Config.V_Att = 2;  % Output voltage of DLP box (Voltage) = V_Att*Digitized data 
Config.I_Att = 5;  % Output voltage of DLP box (Current) = I_att*Digitized data
Stem = '\MPEX::TOP.'; Branch = 'MACHOPS1:'; RootAddress = [Stem,Branch];
DataAddress{1} = [RootAddress,'LP_V_RAMP']; % V
DataAddress{2} = [RootAddress,'TARGET_LP']; % I
Config.V_cal = [12.05,0.205];   % Voltage output of DLP = V_cal(1)*Output voltage of DLP box + V_cal(2)
Config.I_cal = [-142.5, 0.015]; % Current output of DLP = (I_cal(2) + Output voltage of DLP box)/Ical(1)
Config.FilterDataInput = 1; % Filter input data with savitsky Golay filter order 3 frame size 7
Config.SGF = 7;
Config.TimeMode = 2; % Effective time of sweep, (1) start of ramp or (2) mean time of ramp
Config.AMU = 2; % Ion mass in AMU
Config.FitFunction = 2; 
Config.AreaType = 1; % 1: Cylindrical + cap
[ne,Te,tne,Ifit,Ip,Vp,tm,Vsweep,Isweep] = DLP_fit_V5_1(Config,shotlist,DataAddress);
for s = 1:length(shotlist)
    Ne{s} = 0.5*(ne{s}{1} + ne{s}{2});
end
% #########################################################################
% 2 - Gather the RF power trace: RF 
% #########################################################################
DA{1} = [RootAddress,'RF_FWD_PWR'];
[RF,trf]   = my_mdsvalue_v2(shotlist,DA(1))
% #########################################################################
% 3 - Gather Pressure gauges: PG1,2,3,4
% #########################################################################
mdsconnect('mpexserver')
address{1} = '\MPEX::TOP.MACHOPS1:PG1'; % PG9.5
address{2} = '\MPEX::TOP.MACHOPS1:PG2'; % PG2.5
address{3} = '\MPEX::TOP.MACHOPS1:PG3'; % PG6.5
address{4} = '\MPEX::TOP.MACHOPS1:PG4'; % PG4.5

% Plasma shots:
[pg1,tpg1] = my_mdsvalue_v2(shotlist(1),address(1));
[pg2,tpg2] = my_mdsvalue_v2(shotlist(1),address(2));
[pg3,tpg3] = my_mdsvalue_v2(shotlist(1),address(3));
[pg4,tpg4] = my_mdsvalue_v2(shotlist(1),address(4));
% Calibration shot:
Calshot = 13968;
[pg1cal,tpg1cal] = my_mdsvalue_v2(Calshot,address(1));
[pg2cal,tpg2cal] = my_mdsvalue_v2(Calshot,address(2));
[pg3cal,tpg3cal] = my_mdsvalue_v2(Calshot,address(3));
[pg4cal,tpg4cal] = my_mdsvalue_v2(Calshot,address(4));
% Gas only shots, no magnets
GasOnlyShot = 13969;
[pg1GasOnly,tpg1GasOnly] = my_mdsvalue_v2(GasOnlyShot,address(1));
[pg2GasOnly,tpg2GasOnly] = my_mdsvalue_v2(GasOnlyShot,address(2));
[pg3GasOnly,tpg3GasOnly] = my_mdsvalue_v2(GasOnlyShot,address(3));
[pg4GasOnly,tpg4GasOnly] = my_mdsvalue_v2(GasOnlyShot,address(4));

for s = 1 :length(shotlist(1))
    PG1{s}         = (pg1{s}-pg1cal{1})*2;
    PG1_GasOnly{s} = (pg1GasOnly{1}-mean(pg1GasOnly{1}(1:100)))*2;
    
    PG2{s} = (pg2{s}-pg2cal{1})*2;
    PG2_GasOnly{s} = (pg2GasOnly{1}-mean(pg2GasOnly{1}(1:100)))*2;
    
    PG3{s} = (pg3{s}-pg3cal{1})*2;
    PG3_GasOnly{s} = (pg3GasOnly{1}-mean(pg3GasOnly{1}(1:100)))*2;
    
    PG4{s} = (pg4{s}-pg4cal{1})*10;
    PG4_GasOnly{s} = (pg4GasOnly{1}-mean(pg4GasOnly{1}(1:100)))*10;
end
T = {'9.5','2.5','6.5','4.5'};
% #########################################################################
% 4 - Loading the visible camera video (Probe out): S
% #########################################################################

FileName{1} = 'slomo_c_13962';
Address = 'C:\Users\nfc\Documents\Proto-MPEX Data Analysis\2017\04 - April\13th - Long RF pulse + ECH\Visible camera selected\';
f = VideoReader([Address,FileName{1},'.mov']);

vidHeight = f.Height;
vidWidth  = f.Width;
S = struct('cdata',zeros(vidHeight,vidWidth,3,'uint8'),...
'colormap',[]);
n = 1;
while hasFrame(f)
    S(n).cdata = readFrame(f);
    n = n+1;
end

if 0
    figure
    set(gcf,'position',[150 150 0.7*vidWidth 0.7*vidHeight]);
    set(gca,'units','pixels');
    set(gca,'position',[0 0 0.5*vidWidth 0.5*vidHeight]);
    movie(S,1,180)
end
    save('HeliconOnly_CCamera_Gas_58_62.mat')
else
    load('HeliconOnly_CCamera_Gas_58_62.mat')
end
%%
%##########################################################################
% 5 - Analyse Video: OnAxisEmission and RadialEmission
% #########################################################################
VidFrameOffset = 163 - 4165;
t = [[1:1:500]-VidFrameOffset]; % 1 ms per frame, since frame rate is 500 fr/s and 500 ms movie
for s = 1:1:500
       R(:,:,s) = S(s).cdata(:,:,1)/3;
       G(:,:,s) = S(s).cdata(:,:,2)/3;
       B(:,:,s) = S(s).cdata(:,:,3)/3;
       fr(:,:,s) = (R(:,:,s)+ G(:,:,s) + B(:,:,s));
   IntensityIntegrated(s) = sum(sum(fr(:,:,s),1),2);
   x1 = 440; x2 = 460;
   y1 = 5  ; y2 = 688-y1;
   rp(:,:,s) = fr(y1:y2,x1:x2,s); 
   RadialEmission(:,s) = mean(rp(:,:,s),2);
   OnAxisEmission(:,s) = mean(mean(rp(280:320,:,s),2),1);
end
%% Load data from Helicon + ECH case, Camera data
D = load('HeliconECH_66_62_Dataset_1.mat');
RadialEmissionECH = D.RadialEmission;
I_end_ECH = D.I_end;
%% Plot data:
close all
TimePlotStart = 4;
TimePlotEnd   = 4.6;
f = gcf;
% Define the window white in color and with no menubar. the normalized
% units are relative to the computer screen
set(f,'Menubar','figure','color','w','Units','normalized')
% Set the figure to cover the Right half of the screen
set(f,'Position',[0.4 0.1 0.5 0.8]); % [Left Bottom Width height]
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

% RF
n = 5;
ax(n).NextPlot = 'add';
ax(n).YTick = [0.5,1];
ax(n).FontName = 'Times';
ax(n).YLabel.String = 'RF [a.u]'
RFsq = RF{1}.^2;
RFnormalized = RFsq./max(RFsq);
plot(trf{1}(1:end-1),RFnormalized,'Parent',ax(n),'color','k')
ylim(ax(n),[0,1])

% Ne and Te
n = 4;
ax(n).NextPlot = 'add';
ax(n).YTick = [1,2,3,4];
ax(n).FontName = 'Times';
ylabel(ax(n),'$n_e$ , $T_e$','FontSize',13);
ax(n).YLabel.Interpreter = 'latex';
ax(n).YLabel.Rotation = 90; 
ax(n).YLabel.HorizontalAlignment = 'center'; ax(n).YLabel.VerticalAlignment = 'bottom';
hne(1) = plot(tne{2},Ne{2}*1e-19,'Parent',ax(n),'color','k')
hne(2) = plot(tne{2}(8:end-2),Te{2}(8:end-2),'Parent',ax(n),'color','r')
ylim(ax(n),[0,5])
Lne = legend(hne,'$n_e \times 10^{-19}$','$T_e$');
Lne.Interpreter = 'latex';
Lne.Box = 'off'; Lne.Location = 'northwest';
ax(n).YGrid = 'on';

% PG spool 4.5
n = 3;
ax(n).NextPlot = 'add';
ax(n).YTick = [5:5:20];
ax(n).FontName = 'Times';
ylabel(ax(n),'$P_n$ [mTorr]','FontSize',11);
ax(n).YLabel.Interpreter = 'latex';
ax(n).YLabel.Rotation = 90; 
ax(n).YLabel.HorizontalAlignment = 'center'; ax(n).YLabel.VerticalAlignment = 'bottom';
hpgA(1) = plot(tpg2{1}(1:end-1),PG4{1},'Parent',ax(n),'color','k')
hpgA(2) = plot(tpg2{1}(1:end-1),PG4_GasOnly{1},'Parent',ax(n),'color','r')
ylim(ax(n),[-0.1,25])
ax(n).YGrid = 'on';
LpgA = legend(hpgA,'w/ Plasma','w/o Plasma');
LpgA.Interpreter = 'latex';
LpgA.Box = 'off'; LpgA.Location = 'northwest';

% PG spool 6.5
n = 2;
ax(n).NextPlot = 'add';
ax(n).YTick = [1:1:7];
ax(n).FontName = 'Times';
ylabel(ax(n),'$P_n$ [mTorr]','FontSize',11);
ax(n).YLabel.Interpreter = 'latex';
ax(n).YLabel.Rotation = 90; 
ax(n).YLabel.HorizontalAlignment = 'center'; ax(n).YLabel.VerticalAlignment = 'bottom';
hpgB(1) = plot(tpg2{1}(1:end-1),PG3{1},'Parent',ax(n),'color','k')
hpgB(2) = plot(tpg2{1}(1:end-1),PG3_GasOnly{1},'Parent',ax(n),'color','r')
ylim(ax(n),[-0.05,8])
ax(n).YGrid = 'on';
LpgB = legend(hpgB,'w/ Plasma','w/o Plasma');
LpgB.Interpreter = 'latex';
LpgB.Box = 'off'; LpgB.Location = 'northwest';

% Light Emission from central chamber
n = 1;
ax(n).NextPlot = 'add';
ax(n).YTick = [0.2:0.2:0.8];
ax(n).FontName = 'Times';
ylabel(ax(n),'$I_{D_2}$ [a.u]','FontSize',11);
ax(n).YLabel.Interpreter = 'latex';
ax(n).YLabel.Rotation = 90; 
ax(n).YLabel.HorizontalAlignment = 'center'; ax(n).YLabel.VerticalAlignment = 'bottom';
plot(t*1e-3,OnAxisEmission./max(OnAxisEmission),'Parent',ax(n),'color','k')
ylim(ax(n),[-0.05,1])
    ax(n).YGrid = 'on';
%%
 % Visible camera images
close all
t1 = 4.177;
t2 = 4.43;
rng = find(t*1e-3 >= t1 & t*1e-3<=t2);

figure
fvid = gcf;
set(fvid,'Menubar','figure','color','w','Units','normalized')
for s = 1:4
    axvid(s) = axes;
    set(axvid(s),'Units','Normalized','box','on')%,'NextPlot','add')
end

I_start_HeliconOnly = (S(rng(1)).cdata); %rgb2gray
imh(1) = image(axvid(1),I_start_HeliconOnly);   %axvid(1) = gca;
% text(axvid(1),0.08,0.9,'a)','Units','normalized','color','w','FontSize',14,'FontName','Times')
% rectangle(axvid(1),'Position',[x1,y1,x2-x1,y2-y1],'EdgeColor','w')
line(axvid(1),[x1,x1],[y1,y2],'Color','w','LineWidth',1,'LineStyle',':')
text(axvid(1),0.1,0.1,'t = 4.17 s','Units','normalized','color','w','FontSize',10,'FontName','Times')

I_end_HeliconOnly = (S(rng(end)).cdata);
imh(2) = image(axvid(2),I_end_HeliconOnly); hold on %axvid(2) = gca;
% text(axvid(2),0.08,0.9,'b)','Units','normalized','color','w','FontSize',14,'FontName','Times')
%rectangle(axvid(2),'Position',[x1,y1,x2-x1,y2-y1],'EdgeColor','w')
line(axvid(2),[x1,x1],[y1,y2],'Color','w','LineWidth',1,'LineStyle',':')
text(axvid(2),0.1,0.1,'t = 4.43 s','Units','normalized','color','w','FontSize',10,'FontName','Times')

imh(3) = image(axvid(3),I_end_ECH); hold on
% text(axvid(3),0.08,0.9,'c)','Units','normalized','color','w','FontSize',14,'FontName','Times')
%rectangle(axvid(2),'Position',[x1,y1,x2-x1,y2-y1],'EdgeColor','w')
line(axvid(3),[x1,x1],[y1,0.85*y2],'Color','w','LineWidth',1,'LineStyle',':')
text(axvid(3),0.1,0.1,'t = 4.43 s (28 GHz)','Units','normalized','color','w','FontSize',10,'FontName','Times')

dn = 10;
RE1 = mean(RadialEmission(:,rng(1:dn)),2);
plot(RE1(end:-1:1),1:length(RE1),'Parent',axvid(4),'color','k');
% text(axvid(4),0.63,0.5,'a)','Units','normalized','color','k','FontSize',12,'FontName','Times')

n2 = length(rng);n1 = n2-dn;
RE2 = mean(RadialEmission(:,rng(n1:n2)),2);
plot(RE2(end:-1:1),1:length(RE2),'Parent',axvid(4),'color','bl');
% text(axvid(4),0.38,0.5,'b)','Units','normalized','color','k','FontSize',10,'FontName','Times')

RE2_ECH = mean(RadialEmissionECH{2}(:,rng(n1:n2)),2);
plot(RE2_ECH(end:-1:1),1:length(RE2_ECH),'Parent',axvid(4),'color','r');
% text(axvid(3),0.28,0.5,'b)','Units','normalized','color','k','FontSize',10,'FontName','Times')

ylim([-10,688-10])

% dx = 0.1;
% dy = 0.1;
% w = (1-2*dx)/4;
% h = (1-2*dy);
% 
% axvid(1).Position = [dx   dy w h];
% axvid(2).Position = [dx+1.1*w dy w h];
% axvid(3).Position = [dx+2.2*w dy w h];
% axvid(4).Position = [dx+3.3*w dy w h];

dx = 0.1;
dy = 0.1;
w = (1-2*dx)/2;
h = (1-2*dy)/2;

axvid(1).Position = [dx     dy+1.1*h  w h]; %(a)
axvid(2).Position = [dx+1*w dy+1.1*h  w h]; %(b)
axvid(3).Position = [dx+1*w dy        w h]; %(c)
axvid(4).Position = [dx     dy        w h];

axis(axvid,'square')
set(axvid,'XTick',[],'YTick',[])

% Save data for access by other scripts
if 0 
save('HeliconOnly_58_62_Dataset_1','RF','trf','PG1','PG2','PG3','PG4','tpg1','tpg2','tpg3','tpg4'...
    ,'PG1_GasOnly','PG2_GasOnly','PG3_GasOnly','PG4_GasOnly','OnAxisEmission','t')
end
%% test abel transform
close all
i = 1:length(RE2_ECH);
figure; hold on
plot(i-340,RE2_ECH)
rng1 = i(340:end);
plot(i(rng1)-340,RE2_ECH(rng1),'r')
rng2 = i(340:-1:1);
plot(i(rng2)-340,RE2_ECH(rng2),'g')

y = @(x,a,b,c) a*( tanh(b*(x-c)) - tanh(b*(x+c)) );

x1 = linspace(-500,500,300); 
y1 = y(x1,18.5,0.025,-149)-1*y(x1,12.8,0.025,-105) + 8*gaussian(x1,0,230);

plot(x1,y1,'r:')

 [ f_rec , X ] = abel_inversion(y1(150:end),4,24,1,0);
%  [ f_rec , X ] = abel_inversion(RE2_ECH(rng1),4,24,1,0);
%%
close all
i = 1:length(RE2);
figure; hold on
plot(i-340,RE2,'m')
rng1 = i(340:end);
plot(i(rng1)-340,RE1(rng1),'r')
rng2 = i(340:-1:1)
plot(i(rng2)-340,RE1(rng2),'g.')

y = @(x,a,b,c) a*( tanh(b*(x-c)) - tanh(b*(x+c)) );

x1 = linspace(-500,500,300); 
y1 = 1.7*(y(x1,17,0.025,-149)-1*y(x1,11.1,0.025,-105) + 8*gaussian(x1,0,370));

plot(x1,y1,'r:')

[ f_rec , X ] = abel_inversion(0.1*y1(150:end),4,24,1,0);
% [ f_rec , X ] = abel_inversion(sgolay_t(RE1(rng1),3,51),4,104,1,0);
% [ f_rec , X ] = abel_inversion(sgolay_t(RE2(340:end),3,51),4,104,1,0);

