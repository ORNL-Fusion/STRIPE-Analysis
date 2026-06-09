clear all
close all

DataName = 'RawData_GasFlowLevelScan_2018_05_01';
load(DataName)
struct2table(AddressTable)
RawData

%% Preview data
FlowLevel = RawData.MetaData.dt;
figure;
[~,b] = sort(FlowLevel);

for s = 1:length(FlowLevel)
    subplot(2,3,s)
    plot(RawData.t_I{b(s)}(1:end-1),(RawData.I{b(s)}))
    ylim([0,0.2])
    xlim([4,5])
    title(['R: ',num2str(FlowLevel(b(s)))])
end

%% DLP configuration and calculations
Config.tStart          = 4.15;
Config.tEnd            = 4.69 ;
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

DLPData = DLP_fit_V6(Config,RawData);

%%
C = {'k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:'};
close all
figure
tStart = 3.6
tEnd   = 4.8

subplot(4,1,3); hold on
for s = [1:4]
    hold on
    g = DLPData.Ni{s}{1};
    t_g = DLPData.time{s};

    g(t_g>4.648) = 0;
    
    GoodFits{s} = DLPData.GlitchFlag{s} == 0 & DLPData.StdResNorm{s}<=0.2...
    & g>-1e17 & g <5.5e19 & DLPData.Te{s}<=20;

    ne{s} = sgolay_t(g(GoodFits{s}),3,7);
    t_ne{s} =t_g(GoodFits{s});

    plot(t_ne{s},ne{s},'LineWidth',2,'Color',C{s})
    ylim([0,6]*1e19)
    xlim([tStart,tEnd])
    
    Isat{s} = sgolay_t(DLPData.Isat{s}{1}(GoodFits{s}),3,13);
    Isat{s}(t_ne{s}>4.648) = 0;
end

clear g f t_g t_f

subplot(4,1,4); hold on
for s = [2,4]
    hold on
    g = DLPData.Te{s};
    t_g = DLPData.time{s};

    g(t_g>4.648) = 0;
    
        GoodFits{s} = DLPData.GlitchFlag{s} == 0 & DLPData.StdResNorm{s}<=0.1...
    & g>-1e17 & g <5.5e19 & DLPData.Te{s}<=10;
    
    te{s} = sgolay_t(g(GoodFits{s}),3,7);

    t_te{s} =t_g(GoodFits{s});
    plot(t_te{s},te{s},'LineWidth',2,'Color',C{s})
    ylim([0,6])
    xlim([tStart,tEnd])
end

subplot(4,1,1); hold on
for s = [1,5]
    hold on
    p1{s} = (2/7.5)*RawData.PG2{s};
    t_p1{s} = RawData.t_PG2{s}(1:end-1);

    plot(t_p1{s},p1{s},'LineWidth',2,'Color',C{s})
    ylim([-0.05,2])
    xlim([tStart,tEnd])
end

subplot(4,1,2); hold on
for s = [1,5]
    hold on
    p3{s} = (2/7.5)*RawData.PG3{s};
    t_p3{s} = RawData.t_PG3{s}(1:end-1);

    if s == 1
    p3{s}(t_p3{s}<4.06) = NaN;
    end
       
    plot(t_p3{s},p3{s}-0.04,'LineWidth',2,'Color',C{s})
    ylim([0,0.3])
    xlim([tStart,tEnd])
end

set(findobj('-property','YTick'),'box','on')
set(findobj('-property','NextPlot'),'color','w')
set(gcf,'position',[380.3333   73.6667  263.3333  532.0000])

%% Plot for paper:
clear ax
C = {'k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:'};
TimePlotStart = 3.60;
TimePlotEnd   = 4.8;
tRF = 4.15;
%==========================================================================
figure;
f = gcf;
% Define the window white in color and with no menubar. the normalized
% units are relative to the computer screen
set(f,'Menubar','figure','color','w','Units','normalized')
% Set the figure to cover the Right half of the screen
set(f,'Position',[0.6 0.1 0.3 0.75]); % [Left Bottom Width height]

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


% PG spool 2.5
n = 3;
ax(n).NextPlot = 'add';
ax(n).YTick = [0.5:0.5:2];
ax(n).FontName = 'Times';
ylabel(ax(n),'[Pa]','FontSize',11); %
ax(n).YLabel.Interpreter = 'latex';
ax(n).YLabel.Rotation = 90; 
ax(n).YLabel.HorizontalAlignment = 'center'; ax(n).YLabel.VerticalAlignment = 'bottom';

s = 1;
hpg1(s) = plot(t_p1{s}-0.012,p1{s}+0.022,'Parent',ax(n),'color','k');
hpg1(s).LineWidth = 2;

s = 5;
hpg1(s) = plot(t_p1{s}-0.012,p1{s}+0.022,'Parent',ax(n),'color','r');
hpg1(s).LineWidth = 2;

line([tRF,tRF],[0,5],'color','k','linestyle',':','LineWidth',2)

ylim(ax(n),[-0.005,2.5])
ax(n).YGrid = 'on';
text(ax(n),0.02,0.8,'$P1$','Units','normalized','Interpreter','latex','FontSize',10)

% PG spool 6.5
n = 2;
ax(n).NextPlot = 'add';
ax(n).YTick = [0.2:0.2:0.8];
ax(n).FontName = 'Times';
ylabel(ax(n),'[Pa]','FontSize',11); %
ax(n).YLabel.Interpreter = 'latex';
ax(n).YLabel.Rotation = 90; 
ax(n).YLabel.HorizontalAlignment = 'center'; ax(n).YLabel.VerticalAlignment = 'bottom';

s = 1;
hpg6(s) = plot(t_p3{s}-0.03,p3{s}-0.045,'Parent',ax(n),'color','k')
hpg6(s).LineWidth = 2;

s = 5;
hpg6(s) = plot(t_p3{s}-0.03,p3{s}-0.045,'Parent',ax(n),'color','r')
hpg6(s).LineWidth = 2;

line(ax(n),[tRF,tRF],[0,5],'color','k','linestyle',':','LineWidth',2)

ylim(ax(n),[-0.005,1])
ax(n).YGrid = 'on';
text(ax(n),0.02,0.8,'$P3$','Units','normalized','Interpreter','latex','FontSize',10)

% Isat DLP
n = 1;
ax(n).NextPlot = 'add';
ax(n).YTick = [20:20:80];
ax(n).FontName = 'Times';
ylabel(ax(n),'$I_+$ $[A]$','FontSize',11); %
ax(n).YLabel.Interpreter = 'latex';
ax(n).YLabel.Rotation = 90; 
ax(n).YLabel.HorizontalAlignment = 'center'; ax(n).YLabel.VerticalAlignment = 'bottom';

s = 2;
hIsat(s) = plot(t_ne{s},Isat{s}*1e3,'Parent',ax(n),'color','k');
hIsat(s).LineWidth = 2;

s = 4;
hIsat(s) = plot(t_ne{s},Isat{s}*1e3,'Parent',ax(n),'color','r');
hIsat(s).LineWidth = 2;
 
line(ax(n),[tRF,tRF],[0,100],'color','k','linestyle',':','LineWidth',2)
ylim(ax(n),[0,100])
ax(n).YGrid = 'on';
text(ax(n),0.02,0.8,'Probe D','Units','normalized','Interpreter','latex','FontSize',10)


n = 1;
xlabel(ax(n),'$time$ [s]','Interpreter','latex');
ax(n).XTick = [3.6:0.2:5];


abc = {'(c)', '(b)','(a)'};
for s = 1:length(ax)
    ax(s).YLabel.FontSize = 10;
    text(ax(s),1.01,0.85,abc{s},'Units','normalized','Interpreter','latex','FontSize',11)
end
set(gcf,'position',[0.6000    0.2787    0.3000    0.5713])