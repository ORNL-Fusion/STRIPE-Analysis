close all
clear all

SaveData = 0;
load('SourceMagneticField_2017_12_01.mat');

if SaveData
shotlist = [18000 + [179,180,181,182,184,185,186,187,188,189,190,191,192,193,194,196,198,199,204,205,206,207,208,209,210,212,213,214,216,217,218,219,220,221]]; %
TR2      =          [160,180,200,200,220,220,240,240,240,260,260,260,280,280,300,300,330,330,160,160,160,140,140,120,120,100, 80, 80, 40, 40, 20, 20,  0,  0] ; %

HomeAddress = cd;
cd('\\mpexserver\protompex_data\automated_analysis\GUI')

mdsconnect('mpexserver')
% Inputs
%     shot               single shot number
% Outputs
%     heat_flux          heat flux 3-D array in X,Y,T
%     heat_peak          peak heat flux, 1-D array in T
%     heat_peak_x        x location of peak heat flux, 1-D array in T
%     heat_peak_y        y location of peak heat flux, 1-D array in T
%     heat_power         power to target 1-D array in T
%     heat_radius        radius of heat flux 1-D array in T
%     T                  time 1-D array
%     X                  x 1-D array
%     Y                  Y 1-D array
%     comments           comments
a = 1;
for s = 1:length(shotlist)
   [heat_flux{s},heat_peak{s},heat_peak_x{s}...
       ,heat_peak_y{s},heat_power{s}...
       ,heat_radius{s},T{s},X{s},Y{s}...
       ,comments{s}] = get_IR(shotlist(s));
   
   heat_flux{s} = a*heat_flux{s};
   heat_peak{s} = a*heat_peak{s};
   heat_peak_x{s} = a*heat_peak_x{s};
   heat_peak_y{s} = a*heat_peak_y{s};
   heat_power{s}  = a*heat_power{s};
   
end
cd(HomeAddress)

Stem = '\MPEX::TOP.'; Branch = 'MACHOPS1:'; RootAddress = [Stem,Branch];
DA{1} = [RootAddress,'RF_FWD_PWR'];
[RF,t_rf]   = my_mdsvalue_v2(shotlist,DA(1))

 save('PostProcessData_IR_vs_B0_2017_12_01.mat')
else
 load('PostProcessData_IR_vs_B0_2017_12_01.mat')
end

%% Observe IR video
close all
shotlist = [18000 + [179,180,181,182,184,185,186,187,188,189,190,191,192,193,194,196,198,199,204,205,206,207,208,209,210,212,213,214,216,217,218,219,220,221]]; %
TR2      =          [160,180,200,200,220,220,240,240,240,260,260,260,280,280,300,300,330,330,160,160,160,140,140,120,120,100, 80, 80, 40, 40, 20, 20,  0,  0] ; %

% NOTES
% Dont use the follownit shots:
% 18214: 
% 

p = find(shotlist == [18000+179]); 
[~,b] = sort(TR2);
% p = b(20)

figure; hold on
plot(t_rf{p}(1:end-1),RF{p}/min(RF{p}))
plot(T{p}-0.035,heat_peak{p}/max(heat_peak{p}))

%  The important times are:
%  4.3, 4.4, 4.5, 4.65
tFrame = [4.3,4.4,4.5,4.65];
tFrame = [4.22,4.43,4.55,4.65];

if 1
    figure
    for frame = 10:60;
    imagesc(X{p}*100,Y{p}*100,heat_flux{p}(:,:,frame)*1e-6)
    caxis([0,0.3])
    colormap(flipud(hot))
    axis('square')
    xlim([-3,3])
    ylim([-3,3])
    title(num2str(T{p}(frame)-0.035))
    drawnow
%     pause(0.01)
    end
end
%% Plot steady-state IR heat fluxes vs Bo and RF power
close all
k = 0;
figure
N   = [2,4 ,6 ,7 ,9  ,10 ,12 ,14 ,18 ,19 ,21 ,23 ,27 ,29 ,32 ,34 ];
tr2 = [0,20,40,80,100,120,140,160,180,200,220,240,260,280,300,330];
K = length(tFrame(1:3));
for s = N
    for n = 1:K
        k = k + 1;
        nFrame = find(T{b(s)}-0.035 >= tFrame(n),1);
        HeatFlux{b(s)}{n} = mean(heat_flux{b(s)}(:,:,nFrame-1:nFrame+1),3);
        subplot(length(N),K,k)
        try
        B0_HeatFlux2D{b(s)}{n} = B0Data.B0(b(s));
        HeatFlux2D{b(s)}{n} = HeatFlux{b(s)}{n};
        
        imagesc(X{b(s)}*100,Y{b(s)}*100,HeatFlux2D{b(s)}{n}*1e-6)
        PeakHeatFlux(b(s),n) = max(max(HeatFlux{b(s)}{n}));
        TotalHeat(b(s),n) = mean(heat_power{b(s)}(:,nFrame-1:nFrame+1));

        catch
            continue
        end
        if n==1
            caxis([0,0.9])
        elseif n==2
            caxis([0,0.2])
        else
            caxis([0,0.07])            
        end
         colormap(flipud(hot))
%         colormap(hot)
        axis('square')
        xlim([-3,3])
        ylim([-3,3])
        set(gca,'XTick',[],'YTick',[])
%         title(num2str(T{b(s)}(nFrame)-0.035,4),'FontSize',4)
        title(num2str(TR2(b(s))),'FontSize',6)

    end
end

% Peak Heat flux to target
figure; hold on
for n = 1:K
    plot(B0Data.B0(b(N)),PeakHeatFlux(b(N),n))
end
ylim([0,1.3e6])

% Total Heat to target
figure; hold on
for n = 1:K
    plot(B0Data.B0(b(N)),TotalHeat(b(N),n))
end
ylim([0,700])

%% Plot hand picked IR data
close all
N          = [2 ,4 ,6 ,7 ,9  ,10 ,12 ,14 ,18 ,19 ,21 ,23 ,27 ,29 ,32 ,34 ];
tr2        = [0 ,20,40,80,100,120,140,160,180,200,220,240,260,280,300,330];
ShotEnable = [0 ,1 ,0 ,0 ,1  ,0  ,1  ,0  ,1  ,0  ,1  ,0  ,1  ,0  ,1  ,0  ];
ShotEnable = [0 ,1 ,0 ,0 ,0  ,1  ,0  ,1  ,0  ,1  ,0  ,1  ,0  ,1  ,0  ,1  ];
ShotEnable = [1 ,1 ,0 ,0 ,1  ,0  ,0  ,1  ,0  ,0  ,0  ,1  ,0  ,0  ,0  ,1  ];
% ShotEnable = [1 ,0 ,0 ,1 ,0  ,0  ,0  ,0  ,1  ,0  ,0  ,1  ,0  ,0  ,0  ,1  ];


ShotEnable = [1 ,0 ,0 ,1 ,0  ,0  ,0  ,0  ,1  ,0  ,0  ,0  ,1  ,0  ,0  ,1  ];
Powerlevel = [90,70,50];

SelectedShots = find(ShotEnable > 0);
k = 3*length(SelectedShots)+1;
K = length(tFrame(1:3));
% figure
for s = N(SelectedShots)
    for n = 1:K
        k = k - 1;
        nFrame = find(T{b(s)}-0.035 >= tFrame(n),1);
%        subplot(length(SelectedShots),K,k)
     figure
        try        
            imagesc(X{b(s)}*100,Y{b(s)}*100,HeatFlux2D{b(s)}{n}*1e-6)
        catch
            continue
        end
        if n==1
            caxis([0,0.50])
        elseif n==2
            caxis([0,0.2])
        else
            caxis([0,0.1])            
        end

%          colormap(flipud(hot))
        colormap(hot)
        colorbar
        axis('square')
        xlim([-3,3])
        ylim([-3,3])
        title([num2str(B0Data.B0(b(s))*1000,3),' [mT], ',num2str(Powerlevel(n)),' [kW]'],'FontSize',10)
        set(gcf,'position',[360.3333  373.6667  392.0000  244.0000],'color','w')
        
        set(gca,'XTick',[],'YTick',[])
        set(gca,'YTick',[-2,-1,0,1,2],'XTick',[-2,-1,0,1,2],'FontName','Times')
    end
end

%% Check Quality of IR data
% close all
k = 0;
figure
N = find(TR2(b) == 240);
for s = N
    for n = 1:length(tFrame)
        k = k + 1;
        nFrame = find(T{b(s)}-0.035 >= tFrame(n),1);
%         HeatFlux{b(s)}{n} = mean(heat_flux{b(s)}(:,:,nFrame-3:nFrame),3);
        
        subplot(length(N),3,k)
        try
        imagesc(X{b(s)}*100,Y{b(s)}*100,HeatFlux{b(s)}{n}*1e-6)
        catch
            continue
        end
        caxis([0,0.3])
    %     colormap(flipud(hot))
        colormap(hot)
        axis('square')
        xlim([-3,3])
        ylim([-3,3])
        set(gca,'XTick',[],'YTick',[])
        title(num2str(T{b(s)}(nFrame)-0.035,4),'FontSize',4)
    end
end

%%
% #########################################################################
% Steady state DLP values
% #########################################################################
close all
DLPData = load('PostProcessData_ne_vs_B0_2017_12_01.mat');
for s = 1:length(shotlist)
    GoodFits{s} = DLPData.GlitchFlag{s} == 0 & DLPData.StdResNorm{s}<=0.081...
        & DLPData.Ni{s}>0 & DLPData.Ni{s}<1e21 & DLPData.Te{s}<=4;
end

for s = 1:length(shotlist)
    rngPwr{1} = [DLPData.time{s}>=4.25 & DLPData.time{s}<=4.35]; % 100
    rngPwr{2} = [DLPData.time{s}>=4.40 & DLPData.time{s}<=4.45]; % 80 
    rngPwr{3} = [DLPData.time{s}>=4.50 & DLPData.time{s}<=4.55]; % 60
    rngPwr{4} = [DLPData.time{s}>=4.60 & DLPData.time{s}<=4.65]; % 40
    for p = 1:length(rngPwr)
        ne{p}(s)  = mean(DLPData.Ni{s}(rngPwr{p} & GoodFits{s}));
        dne{p}(s) = std(DLPData.Ni{s}(rngPwr{p} & GoodFits{s}),1,2);
        Te{p}(s)  = mean(DLPData.Te{s}(rngPwr{p} & GoodFits{s}));
        dTe{p}(s) = std(DLPData.Te{s}(rngPwr{p} & GoodFits{s}),1,2);
    end
end

% close all
% [a,b] = sort(TR2);

figure; 
hold on
for p = 1:4
% h(p) = errorbar(TR2(b),ne{p}(b),dne{p}(b),'r.-','LineWidth',1)
h(p) = errorbar(B0Data.B0(b),ne{p}(b),dne{p}(b));
end
h(1).Color = 'k';
h(2).Color = 'r';
h(3).Color = 'g';
h(4).Color = 'm';
h(1).Marker = 'sq';
h(2).Marker = 'o';
h(3).Marker = '^';
h(4).Marker = '*';

for s = 1:4
    h(s).LineStyle = '--';
end

xlim([0,0.1]);  

box on
t = title('$n_e$ $[m^{-3}]$ '); t.Interpreter = 'latex'; t.FontSize = 13;
set(gca,'FontName','Times')
set(gcf,'Position', [360.3333  340.3333  415.3333  277.3333],'color','w')
t = ylabel(gca,'$n_e$ $[m^{-3}]$ '); t.Interpreter = 'latex'; t.FontSize = 14;
t = xlabel(gca,'$B_0$ $[T]$ '); t.Interpreter = 'latex'; t.FontSize = 12;
legend(h,'90 kW','70 kW','50 kW','30 kW')

figure; hold on
for p = 1:4
h(p) = errorbar(B0Data.B0(b),Te{p}(b),dTe{p}(b),'LineWidth',1);
end
h(1).Color = 'k';
h(2).Color = 'r';
h(3).Color = 'g';
h(4).Color = 'm';
h(1).Marker = 'sq';
h(2).Marker = 'o';
h(3).Marker = '^';
h(4).Marker = '*';

ylim([0,5])
xlim([0.02,0.1]); 
box on
t = title('$T_e$ $[eV]$ '); t.Interpreter = 'latex'; t.FontSize = 13;
set(gcf,'color','w')

%% Compare DLP and IR data
close all
figure; hold on
yyaxis left
ax(1) = gca;
hDLP = errorbar(B0Data.B0(b),ne{1}(b),dne{1}(b));
% h(2) = plot(TR2(b(N)),PeakHeatFlux(b(N),1)/(4.2*e_c*2.3*C_s(2.3,2)))
% h(2) = plot(TR2(b(N)),PeakHeatFlux(b(N),1)/(4*e_c*2*C_s(2,2)))
set(ax(1),'XColor','k','YColor','k')
ylim([0,8e19])
xlim([0,0.1])
set(hDLP,'LineStyle','none','Marker','o','MarkerSize',6)
t = ylabel(ax(1),'$n_e$ $[m^{-3}]$ '); t.Interpreter = 'latex'; t.FontSize = 14;

yyaxis right
PlotType = 1;
switch PlotType
    case 1
        hIR = plot(B0Data.B0(b(N)),PeakHeatFlux(b(N),1)*1e-6,'rsq');
        set(hIR,'MarkerFaceColor','r')
        ylim([0,1.5])
        xlim([0,0.1])
    case 2
        hIR = plot(B0Data.B0(b(N)),TotalHeat(b(N),1))
        ylim([0,600])
        xlim([0,0.1])
end
ax(2) = gca;
set(ax(2),'XColor','k','YColor','k','FontName','times')
t = ylabel(ax(2),'$q$ $[MWm^{-2}]$ '); t.Interpreter = 'latex'; t.FontSize = 13;
t = xlabel(ax(2),'$B_0$ $[T]$ '); t.Interpreter = 'latex'; t.FontSize = 13;
box on
set(gcf,'color','w','position',[360.3333  296.3333  480.0000  321.3333])
grid on
l = legend([hDLP,hIR],'DLP - $n_e$','IR - Heat flux');
set(l,'interpreter','Latex','box','on','FontSize',10,'Location','northwest')

%% Comparing IR and DLP derived heat fluxes