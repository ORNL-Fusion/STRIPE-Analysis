% Helicon + ECH analsysis

% For these set of shots, we do no have a DLP 6.5 shot for 300 ms RF+ ECH, we do
% have it for 200 ms. Moreover, we do have DLP 10.5 shots for 300 ms RF +
% ECH
% For this program, we will skip the DLP analysis and only include the
% following:
% RF + ECH trace
% PGA and PGB, these are spool 4.5 and 6.5, w/ and w/o ECH
% Light emission, w/ and w/o ECH

% Shots to use:
% 13985, 300 ms RF, 90 ms ECH at 4.2 sec, DLP 10.5
% 13986, 300 ms RF, no ECH, Probes out (DLP 6.5 at y = 4 cm)

clear all
close all
CMPT = 0;
if CMPT == 1 % 108 sec
    shotlist = 13900 + [86,85]; % ECH OFF and ECH ON
% #########################################################################
% 1 - Gather the RF power trace: RF 
% #########################################################################
Stem = '\MPEX::TOP.'; Branch = 'MACHOPS1:'; RootAddress = [Stem,Branch];
DA{1} = [RootAddress,'RF_FWD_PWR'];
[RF,trf]   = my_mdsvalue_v2(shotlist,DA(1))

DA{1} = [RootAddress,'PWR_28GHZ']; % Isx
[ECH,t_ech]   = my_mdsvalue_v2(shotlist,DA(1))
% #########################################################################
% 2 - Gather Pressure gauges: PG1,2,3,4
% #########################################################################
mdsconnect('mpexserver')
address{1} = '\MPEX::TOP.MACHOPS1:PG1'; % PG9.5
address{2} = '\MPEX::TOP.MACHOPS1:PG2'; % PG2.5
address{3} = '\MPEX::TOP.MACHOPS1:PG3'; % PG6.5
address{4} = '\MPEX::TOP.MACHOPS1:PG4'; % PG4.5

% Plasma shots:
[pg1,tpg1] = my_mdsvalue_v2(shotlist,address(1));
[pg2,tpg2] = my_mdsvalue_v2(shotlist,address(2));
[pg3,tpg3] = my_mdsvalue_v2(shotlist,address(3));
[pg4,tpg4] = my_mdsvalue_v2(shotlist,address(4));
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

for s = 1 :length(shotlist)
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
% 3 - Loading the visible camera video (Probe out): S
% #########################################################################
 Address = 'C:\Users\nfc\Documents\Proto-MPEX Data Analysis\2017\04 - April\13th - Long RF pulse + ECH\Visible camera selected\';
for s = 1:length(shotlist)
    FileName{s} = ['slomo_c_',num2str(shotlist(s))];
    f{s} = VideoReader([Address,FileName{s},'.mov']);
    vidHeight = f{s}.Height;
    vidWidth  = f{s}.Width;
    S{s} = struct('cdata',zeros(vidHeight,vidWidth,3,'uint8'),...
    'colormap',[]);
    n = 1;
    while hasFrame(f{s})
        S{s}(n).cdata = readFrame(f{s});
        n = n+1;
    end
end

if 1
    s = 1;
    figure
    set(gcf,'position',[150 150 0.7*vidWidth 0.7*vidHeight]);
    set(gca,'units','pixels');
    set(gca,'position',[0 0 0.5*vidWidth 0.5*vidHeight]);
    movie(S{s},1,180)
end
    save('HeliconECH_CCamera_Gas_85_86.mat')
else
    load('HeliconECH_CCamera_Gas_85_86.mat')
end
%%
%##########################################################################
% 5 - Analyse Video: OnAxisEmission and RadialEmission
% #########################################################################
VidFrameOffset = 163 - 4165;
t = [[1:1:500]-VidFrameOffset]; % 1 ms per frame, since frame rate is 500 fr/s and 500 ms movie
for k = 1:length(shotlist)
for s = 1:1:500
       R{k}(:,:,s) = S{k}(s).cdata(:,:,1)/3;
       G{k}(:,:,s) = S{k}(s).cdata(:,:,2)/3;
       B{k}(:,:,s) = S{k}(s).cdata(:,:,3)/3;
       fr{k}(:,:,s) = (R{k}(:,:,s)+ G{k}(:,:,s) + B{k}(:,:,s));
   IntensityIntegrated{k}(s) = sum(sum(fr{k}(:,:,s),1),2);
   x1 = 440; x2 = 460;
   y1 = 5  ; y2 = 688-y1;
   rp{k}(:,:,s) = fr{k}(y1:y2,x1:x2,s); 
   RadialEmission{k}(:,s) = mean(rp{k}(:,:,s),2);
   OnAxisEmission{k}(:,s) = mean(mean(rp{k}(280:320,:,s),2),1);
end
end

% The following needs to be updated:

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
N = 4;
for s = 1:N
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

% RF and ECH
n = 4;
ax(n).NextPlot = 'add';
ax(n).YTick = [0.5,1];
ax(n).FontName = 'Times';
ax(n).YLabel.String = 'RF, {\mu}W [a.u]'
RFsq = RF{1}.^2;
RFnormalized = RFsq./max(RFsq); 
ECHnormalized = ECH{2}./max(ECH{2});
plot(trf{1}(1:end-1),RFnormalized,'Parent',ax(n),'color','k')
plot(t_ech{1}(1:end-1),ECHnormalized/7.5,'Parent',ax(n),'color','r')
ylim(ax(n),[0,1])

% % Ne and Te
% n = 3;
% ax(n).NextPlot = 'add';
% ax(n).YTick = [1,2,3,4];
% ax(n).FontName = 'Times';
% ylabel(ax(n),'$n_e$ , $T_e$','FontSize',13);
% ax(n).YLabel.Interpreter = 'latex';
% ax(n).YLabel.Rotation = 90; 
% ax(n).YLabel.HorizontalAlignment = 'center'; ax(n).YLabel.VerticalAlignment = 'bottom';
% hne(1) = plot(tne{2},Ne{2}*1e-19,'Parent',ax(n),'color','k')
% hne(2) = plot(tne{2}(8:end-2),Te{2}(8:end-2),'Parent',ax(n),'color','r')
% ylim(ax(n),[0,5])
% Lne = legend(hne,'$n_e \times 10^{-19}$','$T_e$');
% Lne.Interpreter = 'latex';
% Lne.Box = 'off'; Lne.Location = 'northwest';
% ax(n).YGrid = 'on';

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
hpgA(2) = plot(tpg2{1}(1:end-1),PG4{2},'Parent',ax(n),'color','r')
ylim(ax(n),[-0.1,25])
ax(n).YGrid = 'on';
LpgA = legend(hpgA,'${\mu}W$ ON','${\mu}W$ OFF');
LpgA.Interpreter = 'latex';
LpgA.Box = 'off'; LpgA.Location = 'northwest';

% PG spool 6.5
n = 2;
ax(n).NextPlot = 'add';
ax(n).YTick = [0.5:0.5:2.5];
ax(n).FontName = 'Times';
ylabel(ax(n),'$P_n$ [mTorr]','FontSize',11);
ax(n).YLabel.Interpreter = 'latex';
ax(n).YLabel.Rotation = 90; 
ax(n).YLabel.HorizontalAlignment = 'center'; ax(n).YLabel.VerticalAlignment = 'bottom';
hpgB(1) = plot(tpg2{1}(1:end-1),PG3{1},'Parent',ax(n),'color','k')
hpgB(2) = plot(tpg2{1}(1:end-1),PG3{2},'Parent',ax(n),'color','r')
ylim(ax(n),[-0.05,3])
ax(n).YGrid = 'on';
LpgB = legend(hpgB,'${\mu}W$ ON','${\mu}W$ OFF');
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
hLE(1) = plot(t*1e-3,0.9*OnAxisEmission{1}./max(OnAxisEmission{2}),'Parent',ax(n),'color','k')
hLE(2) = plot(t*1e-3,0.9*OnAxisEmission{2}./max(OnAxisEmission{2}),'Parent',ax(n),'color','r')
ylim(ax(n),[-0.05,1])
    ax(n).YGrid = 'on';

 % Visible camera images
% close all
t1 = 4.177;
t2 = 4.43;
rng = find(t*1e-3 >= t1 & t*1e-3<=t2);

figure
fvid = gcf;
set(fvid,'Menubar','figure','color','w','Units','normalized')
for s = 1:3
    axvid(s) = axes;
    set(axvid(s),'Units','Normalized','box','on')%,'NextPlot','add')
end

I = (S(rng(1)).cdata); %rgb2gray
imh(1) = image(axvid(1),I);   %axvid(1) = gca;
text(axvid(1),0.08,0.9,'a)','Units','normalized','color','w','FontSize',14,'FontName','Times')
% rectangle(axvid(1),'Position',[x1,y1,x2-x1,y2-y1],'EdgeColor','w')
line(axvid(1),[x1,x1],[y1,y2],'Color','w','LineWidth',1,'LineStyle',':')

I = (S(rng(end)).cdata);
imh(2) = image(axvid(2),I); hold on %axvid(2) = gca;
text(axvid(2),0.08,0.9,'b)','Units','normalized','color','w','FontSize',14,'FontName','Times')
%rectangle(axvid(2),'Position',[x1,y1,x2-x1,y2-y1],'EdgeColor','w')
line(axvid(2),[x1,x1],[y1,y2],'Color','w','LineWidth',1,'LineStyle',':')

plot(RadialEmission(end:-1:1,rng(1)),1:length(RadialEmission(:,rng(end))),'Parent',axvid(3),'color','k');
text(axvid(3),0.63,0.5,'a)','Units','normalized','color','k','FontSize',12,'FontName','Times')

plot(RadialEmission(end:-1:1,rng(end)),1:length(RadialEmission(:,rng(end))),'Parent',axvid(3),'color','k');
text(axvid(3),0.24,0.5,'b)','Units','normalized','color','k','FontSize',10,'FontName','Times')
ylim([-10,688-10])


dx = 0.1;
dy = 0.1;
w = (1-2*dx)/3;
h = (1-2*dy);

axvid(1).Position = [dx   dy w h];
axvid(2).Position = [dx+1.1*w dy w h];
axvid(3).Position = [dx+2.2*w dy w h];
axis(axvid,'square')
set(axvid,'XTick',[],'YTick',[])
%%
figure; 
fvid = gcf;
fvid.Color = 'w'
subplot(1,2,1)
surf(S(rng(1)).cdata(:,:,1),'LineStyle','none')
zlim([0,200])
caxis([0,200])
axis('square')
view([0,90])

subplot(1,2,2)
surf(S(rng(end)).cdata(:,:,1),'LineStyle','none')
zlim([0,200])
caxis([0,200])
axis('square')
view([0,90])

figure
surf((S(rng(1)).cdata(:,:,1)-S(rng(end)).cdata(:,:,1)),'LineStyle','none')
% zlim([-100,200])
% caxis([-100,200])
axis('square')
view([0,90])
colormap('jet')





