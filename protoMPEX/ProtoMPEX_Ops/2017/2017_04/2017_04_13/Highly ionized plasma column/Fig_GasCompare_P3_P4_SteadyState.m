
mdsconnect('mpexserver')
close all
clear all

shotlist = 20116;

PressureInPascal = 1;
RemoveOffsetStartRF = 1;

address{1} = '\MPEX::TOP.MACHOPS1:PG1'; % PG9.5
address{2} = '\MPEX::TOP.MACHOPS1:PG2'; % PG2.5
address{3} = '\MPEX::TOP.MACHOPS1:PG3'; % PG6.5
address{4} = '\MPEX::TOP.MACHOPS1:PG4'; % PG4.5

if strcmp(address{1},'\MPEX::TOP.MACHOPS1:PG3')
     f = 2;
     T = 'PG6.5';
elseif strcmp(address{1},'\MPEX::TOP.MACHOPS1:PG1')
     f = 2;
     T = 'PG9.5';
elseif strcmp(address{1},'\MPEX::TOP.MACHOPS1:PG4')
     f = 10;
     T = 'PG4.5';
elseif strcmp(address{1},'\MPEX::TOP.MACHOPS1:PG2')
     f = 2;
     T = 'PG2.5';
 end
% f = 2; % use 10 for PG4 and 2 for all other

[PG1,t1] = my_mdsvalue_v2(shotlist,address(1));
[PG2,t2] = my_mdsvalue_v2(shotlist,address(2));
[PG3,t3] = my_mdsvalue_v2(shotlist,address(3));
[PG4,t4] = my_mdsvalue_v2(shotlist,address(4));

% Calibration shot:
Calshot = 20908; % 2018_04_05, 2.0 kA case, 1 sec RF pulse

[PG1cal,t0cal] = my_mdsvalue_v2(Calshot,address(1));
[PG2cal,t1cal] = my_mdsvalue_v2(Calshot,address(4));
[PG3cal,t2cal] = my_mdsvalue_v2(Calshot,address(3));
[PG4cal,t4cal] = my_mdsvalue_v2(Calshot,address(4));


%%
figure;
C = {'k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:'};

if PressureInPascal
    fctr = 0.1333;
else
    fctr = 1;
end

L1 = length(PG1{1});
L2 = length(PG1cal{1});
if L1>L2
    Ldata = L2;
else
    Ldata = L1;
end

for s = 1 :length(shotlist)
    P{1}{s} = ( PG1{s}(1:Ldata)-PG1cal{1}(1:Ldata) )*2*fctr;
    P{2}{s} = ( PG2{s}(1:Ldata)-PG2cal{1}(1:Ldata) )*2*fctr;
    P{3}{s} = ( PG3{s}(1:Ldata)-PG3cal{1}(1:Ldata) )*2*fctr;
    P{4}{s} = ( PG4{s}(1:Ldata)-PG4cal{1}(1:Ldata) )*10*fctr;
    
end
T = {'9.5','2.5','6.5','4.5'};

for p = 1:4
subplot(2,2,p); hold on
for s = 1:length(shotlist)
    
    if s == 31
        s
    end
    TimeP = t1{s}(1:length(P{p}{s}));
    
    if RemoveOffsetStartRF
        rngMean = find(TimeP>3.6 & TimeP <4.0);
        rngMean = find(TimeP>2 & TimeP <3.5);

        Offset{p}{s} = mean(P{p}{s}(rngMean));
    else 
        Offset{p}{s} = 0;
    end
    P{p}{s} = P{p}{s} - Offset{p}{s};
    h(s) = plot(TimeP,P{p}{s},C{s});
    set(h(s),'LineWidth',2)
    title(['PG ',T{p}])
    box on
    grid on
    if PressureInPascal
        ylabel('[Pa]')
    else
        ylabel('mTorr')
    end
    xlim([3.5,6])
    switch T{p}
        case '9.5'
            ylim([-0.5,4]*fctr)
        case '2.5'
            ylim([-0.5,25]*fctr)
        case '6.5'
            ylim([-0.5,8]*fctr)
        case '4.5'
            ylim([-0.5,35]*fctr)  
    end
end
end
legend(h,{num2str(shotlist')},'location','NorthWest')
set(gcf,'color','w')

%% Plot for paper:
TimePlotStart = 3.8;
TimePlotEnd   = 5.4;
tRF = 4.15;
%==========================================================================
figure;
f = gcf;
% Define the window white in color and with no menubar. the normalized
% units are relative to the computer screen
set(f,'Menubar','figure','color','w','Units','normalized')
% Set the figure to cover the Right half of the screen
set(f,'Position',[0.6 0.1 0.3 0.75]); % [Left Bottom Width height]

for s = 1:2
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


% PG spool 10.5
n = 2;
ax(n).NextPlot = 'add';
ax(n).YTick = [0.1:0.1:0.3];
ax(n).FontName = 'Times';
ylabel(ax(n),'[Pa]','FontSize',11); %
ax(n).YLabel.Interpreter = 'latex';
ax(n).YLabel.Rotation = 90; 
ax(n).YLabel.HorizontalAlignment = 'center'; ax(n).YLabel.VerticalAlignment = 'bottom';
for s = 1:length(shotlist)
hpg4(s) = plot(TimeP-0.03,P{1}{s},'Parent',ax(n),'color',C{s});
hpg4(s).LineWidth = 2;
end
line([tRF,tRF],[0,5],'color','k','linestyle',':','LineWidth',2)

ylim(ax(n),[-0.005,0.4])
ax(n).YGrid = 'on';
text(ax(n),0.02,0.8,'$P4$','Units','normalized','Interpreter','latex','FontSize',10)


% PG spool 6.5
n = 1;
ax(n).NextPlot = 'add';
ax(n).YTick = [0.2:0.2:0.8];
ax(n).FontName = 'Times';
ylabel(ax(n),'[Pa]','FontSize',11); %
ax(n).YLabel.Interpreter = 'latex';
ax(n).YLabel.Rotation = 90; 
ax(n).YLabel.HorizontalAlignment = 'center'; ax(n).YLabel.VerticalAlignment = 'bottom';
for s = 1:length(shotlist)
hpg6(s) = plot(TimeP-0.03,P{3}{s},'Parent',ax(n),'color',C{s})
hpg6(s).LineWidth = 2;
end
line(ax(n),[tRF,tRF],[0,5],'color','k','linestyle',':','LineWidth',2)

ylim(ax(n),[-0.005,1])
ax(n).YGrid = 'on';
text(ax(n),0.02,0.8,'$P3$','Units','normalized','Interpreter','latex','FontSize',10)


n = 1;
xlabel(ax(n),'$time$ [s]','Interpreter','latex');
ax(n).XTick = [3.8:0.2:5.4];


abc = {'(b)','(a)'};
for s = 1:length(ax)
    ax(s).YLabel.FontSize = 10;
    text(ax(s),1.01,0.85,abc{s},'Units','normalized','Interpreter','latex','FontSize',11)
end
