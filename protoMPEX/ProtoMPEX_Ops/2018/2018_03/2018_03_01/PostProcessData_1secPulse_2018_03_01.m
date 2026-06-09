clear all
close all

DataName = 'RawData_1secPulse_2018_03_01';
load(DataName)
struct2table(AddressTable)
RawData

%% Preview data
PS1 = RawData.MetaData.PS1;
figure;
[~,b] = sort(PS1);

for s = 1:length(PS1)
    subplot(2,2,s)
    plot(RawData.t_I{b(s)}(1:end-1),(RawData.I{b(s)}))
    ylim([-0.2,0.2])
    xlim([4,6])
    title(['PS1: ',num2str(PS1(b(s)))])
end

%% DLP configuration and calculations
Config.tStart          = 4.12;
Config.tEnd            = 5.2 ;
Config.FitFunction     = 2   ;
Config.Center_V        = 0   ;
Config.Center_I        = 0   ;
Config.FilterDataInput = 1   ;
Config.SGF             = 11  ;
Config.TimeMode        = 1   ;
Config.AMU             = 2   ;
Config.AreaType        = 2   ;
Config.L_tip           = 1.8/1000  ; % [m]
Config.D_tip           = 0.254/1000; % [m]
Config.V_Att           = 1   ;
Config.I_Att           = 1   ;
Config.V_cal           = [(0.46e-3)^-1,0]; % Voltage output of DLP
Config.I_cal           = [-1,0]; % Current output of DLP

RawData_DLP = GetRawDataSubset(20000 + [80  ,116],RawData);
RawData_PG  = GetRawDataSubset(20000 + [102 ,105],RawData);

DLPData = DLP_fit_V6(Config,RawData);

%%
C = {'k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:'};
close all
figure
tEnd = 5.9;
tRF = 4.158;

subplot(5,1,1); hold on
for s = [1,2]
    hold on
    % Plasma density
    g = DLPData.Ni_m{s};
    t_g = DLPData.time{s} + 0.006;
    g(t_g>5.149) = 0;
    
    % Electron temperature
    k = DLPData.Te{s};
    t_k = DLPData.time{s} + 0.006;
    k(t_g>5.145) = 0;
    
    GoodFits_ne{s} = DLPData.GlitchFlag{s} == 0 & DLPData.StdResNorm{s}<=0.15...
    & g>-1e17 & g <8e19 & k<=10;

    GoodFits_te{s} = DLPData.GlitchFlag{s} == 0 & DLPData.StdResNorm{s}<=0.15...
    & g>-1e17 & g <8e19 & k<=10;

    ne{s} = sgolay_t(g(GoodFits_ne{s}),3,7);
    t_ne{s} =t_g(GoodFits_ne{s});
   
    te{s} = sgolay_t(k(GoodFits_te{s}),3,7);
    t_te{s} =t_k(GoodFits_te{s});   
    
    
    plot(t_ne{s},ne{s}*1e-19,'LineWidth',2,'Color',C{s})
    plot(t_te{s},te{s}      ,'LineWidth',1,'Color',C{s});

    ylim([0,8])
    xlim([3.7,tEnd])
    line(gca,[tRF,tRF],[0,10],'color','k','linestyle',':','LineWidth',2)

end

subplot(5,1,2); hold on
for s = [1,2,3]
    hold on
    p1{s} = (2/7.5)*(RawData.PG2{s} - mean(RawData.PG2{s}(3500:3800)) );
    t_p1{s} = RawData.t_PG2{s}(1:end-1);
    
    plot(t_p1{s},p1{s},'LineWidth',2,'Color',C{s})
    ylim([0,2.5])
    xlim([3.7,tEnd])
    line(gca,[tRF,tRF],[0,10],'color','k','linestyle',':','LineWidth',2)

end

subplot(5,1,3); hold on
for s = [1,2,3]
    hold on
    p2{s} = (10/7.5)*(RawData.PG4{s} - mean(RawData.PG4{s}(3500:3800)) );
    t_p2{s} = RawData.t_PG4{s}(1:end-1);

    plot(t_p2{s},p2{s},'LineWidth',2,'Color',C{s})
    ylim([0,2.5])
    xlim([3.7,tEnd])
    line(gca,[tRF,tRF],[0,10],'color','k','linestyle',':','LineWidth',2)

end

subplot(5,1,4); hold on
for s = [1,2,3]
    hold on
    p3{s} = (2/7.5)*(RawData.PG3{s} - mean(RawData.PG3{s}(3500:3800)) );
    t_p3{s} = RawData.t_PG3{s}(1:end-1);

    plot(t_p3{s},p3{s},'LineWidth',2,'Color',C{s})
    ylim([0,1])
    xlim([3.7,tEnd])
    line(gca,[tRF,tRF],[0,10],'color','k','linestyle',':','LineWidth',2)

end

subplot(5,1,5); hold on
for s = [1,2,3]
    hold on
    p4{s} = (2/7.5)*(RawData.PG1{s} - mean(RawData.PG1{s}(3500:3800)) );
    t_p4{s} = RawData.t_PG1{s}(1:end-1);

    plot(t_p4{s},p4{s},'LineWidth',2,'Color',C{s})
    ylim([0,0.5])
    xlim([3.7,tEnd])
    line(gca,[tRF,tRF],[0,10],'color','k','linestyle',':','LineWidth',2)

end

set(findobj('-property','YTick'),'box','on')
set(findobj('-property','NextPlot'),'color','w')
set(gcf,'position',[380.3333   73.6667  263.3333  532.0000])

%% Plot for paper:
TimePlotStart = 3.6;
TimePlotEnd   = 6;
%==========================================================================
figure;
f = gcf;
% Define the window white in color and with no menubar. the normalized
% units are relative to the computer screen
set(f,'Menubar','figure','color','w','Units','normalized')
% Set the figure to cover the Right half of the screen
set(f,'Position',[0.6 0.1 0.3 0.75]); % [Left Bottom Width height]

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
dy = 0.15;
w = (1-(2*dx))/1;
h = 1.1*((1-(2*dy))/length(ax));
set(ax(1),'Position',[dx (dy + 0*(h-SmallOffset)) w h])
for s = 1:(length(ax)-1)
set(ax(s+1),'Position',[dx (dy + s*(h-SmallOffset)) w h],'XTick',[])
end

% Ne and Te 
% =========================================================================
n = 5;
ax(n).NextPlot = 'add';
ax(n).YTick = [2:2:8];
ax(n).FontName = 'Times';
ylabel(ax(n),'$[m^{-3}]$','FontSize',11); %
ax(n).YLabel.Interpreter = 'latex';
ax(n).YLabel.Rotation = 90; 
ax(n).YLabel.HorizontalAlignment = 'center'; ax(n).YLabel.VerticalAlignment = 'bottom';

s = 1;
hne = plot(t_ne{s},ne{s}*1e-19,'Parent',ax(n),'color','k');
hne.LineWidth = 2;
hte = plot(t_te{s},te{s},'Parent',ax(n),'color','bl');
hte.LineWidth = 1;

line([tRF,tRF],[0,10],'color','k','linestyle',':','LineWidth',2)

ylim(ax(n),[0,8])
ax(n).YGrid = 'on';
text(ax(n),0.02,0.845,'Probe D','Units','normalized','Interpreter','latex','FontSize',10,'BackgroundColor','none')

Lne = legend([hne,hte],'$n_e{\times}10^{-19}$','$T_e$ $[eV]$');
Lne.Interpreter = 'latex';
Lne.Box = 'off'; Lne.Location = 'East';

% P1 
% =========================================================================
n = 4;
ax(n).NextPlot = 'add';
ax(n).YTick = [1:1:2];
ax(n).FontName = 'Times';
ylabel(ax(n),'[Pa]','FontSize',11); %
ax(n).YLabel.Interpreter = 'latex';
ax(n).YLabel.Rotation = 90; 
ax(n).YLabel.HorizontalAlignment = 'center'; ax(n).YLabel.VerticalAlignment = 'bottom';

s = 2;
hpg1 = plot(t_p1{s},p1{s},'Parent',ax(n),'color','k');
hpg1.LineWidth = 2;

s = 3;
hpg1_0 = plot(t_p1{s},p1{s},'Parent',ax(n),'color','g','LineStyle','-');
hpg1_0.LineWidth = 2;

line(ax(n),[tRF,tRF],[0,5],'color','k','linestyle',':','LineWidth',2)

ylim(ax(n),[0,2.5])
ax(n).YGrid = 'on';
text(ax(n),0.02,0.8,'$P1$','Units','normalized','Interpreter','latex','FontSize',10)

Lpg1 = legend([hpg1,hpg1_0],'$P_n$, w/ Plasma','$P_n$, No Plasma');
Lpg1.Interpreter = 'latex';
Lpg1.Box = 'off'; Lpg1.Location = 'East';

% P2
% =========================================================================
n = 3;
ax(n).NextPlot = 'add';
ax(n).YTick = [1:1:2];
ax(n).FontName = 'Times';
ylabel(ax(n),'[Pa]','FontSize',11); %
ax(n).YLabel.Interpreter = 'latex';
ax(n).YLabel.Rotation = 90; 
ax(n).YLabel.HorizontalAlignment = 'center'; ax(n).YLabel.VerticalAlignment = 'bottom';

s = 2;
hpg1 = plot(t_p2{s},p2{s},'Parent',ax(n),'color','k');
hpg1.LineWidth = 2;

s = 3;
hpg1_0 = plot(t_p2{s},p2{s},'Parent',ax(n),'color','g','LineStyle','-');
hpg1_0.LineWidth = 2;

line(ax(n),[tRF,tRF],[0,5],'color','k','linestyle',':','LineWidth',2)

ylim(ax(n),[0,2.5])
ax(n).YGrid = 'on';
text(ax(n),0.02,0.8,'$P2$','Units','normalized','Interpreter','latex','FontSize',10)

% P3
% =========================================================================
n = 2;
ax(n).NextPlot = 'add';
ax(n).YTick = [0.2:0.2:0.8];
ax(n).FontName = 'Times';
ylabel(ax(n),'[Pa]','FontSize',11); %
ax(n).YLabel.Interpreter = 'latex';
ax(n).YLabel.Rotation = 90; 
ax(n).YLabel.HorizontalAlignment = 'center'; ax(n).YLabel.VerticalAlignment = 'bottom';

s = 2;
hpg1 = plot(t_p3{s},p3{s},'Parent',ax(n),'color','k');
hpg1.LineWidth = 2;

s = 3;
hpg1_0 = plot(t_p3{s},p3{s},'Parent',ax(n),'color','g','LineStyle','-');
hpg1_0.LineWidth = 2;

line(ax(n),[tRF,tRF],[0,5],'color','k','linestyle',':','LineWidth',2)

ylim(ax(n),[0,1])
ax(n).YGrid = 'on';
text(ax(n),0.02,0.8,'$P3$','Units','normalized','Interpreter','latex','FontSize',10)

% P4
% =========================================================================
n = 1;
ax(n).NextPlot = 'add';
ax(n).YTick = [0.1:0.1:0.4];
ax(n).FontName = 'Times';
ylabel(ax(n),'[Pa]','FontSize',11); %
ax(n).YLabel.Interpreter = 'latex';
ax(n).YLabel.Rotation = 90; 
ax(n).YLabel.HorizontalAlignment = 'center'; ax(n).YLabel.VerticalAlignment = 'bottom';

s = 2;
hpg1 = plot(t_p4{s},p4{s},'Parent',ax(n),'color','k');
hpg1.LineWidth = 2;

s = 3;
hpg1_0 = plot(t_p4{s},p4{s},'Parent',ax(n),'color','g','LineStyle','-');
hpg1_0.LineWidth = 2;

line(ax(n),[tRF,tRF],[0,5],'color','k','linestyle',':','LineWidth',2)

ylim(ax(n),[0,0.5])
ax(n).YGrid = 'on';
text(ax(n),0.02,0.8,'$P3$','Units','normalized','Interpreter','latex','FontSize',10)

% Time axis
% =========================================================================
n = 1;
xlabel(ax(n),'$time$ [s]','Interpreter','latex');
ax(n).XTick = [3.5:0.5:TimePlotEnd];

abc = {'(e)','(d)','(c)','(b)','(a)'};
for s = 1:length(ax)
    ax(s).YLabel.FontSize = 10;
    text(ax(s),1.01,0.85,abc{s},'Units','normalized','Interpreter','latex','FontSize',11)
end