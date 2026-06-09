% Figures for "TG suppression/ Edge to core power deposition" paper in
% collaboration with Pawel

clear all
close all
% #########################################################################
% Load data
% #########################################################################
% Details of the data to import:
% =========================================================================
% Figure 1:
% It has its own script
% =========================================================================
% Figure 2:
% It potentially has 4 to 6 parts, these are
% (a) Hollow dIR at 4.2 sec
% (b) Centrally peaked dIR at 4.43 s
% (c) Hollow dIR at 4.43 s
% (d) Radially resolved contour plot in time
% (e) |Br| core and edge from Bdot probe
% =========================================================================
% Figure 3:
% Time evolution of ne at location "A" and "B". Include the normalized RF
% power trace
% =========================================================================
% Figure 4:
% Radial |Bz| from Bdot probe at location "A" for similar conditions as the
% other figures. Gas fueling at z = 2 m, B0 of 0.07 T and 0.75 T at the source and
% Target regions respectively.

% #########################################################################
% Figure 2 data:
% #########################################################################
% Parts (a),(b) and (d)
Dfig2abd = load('IR_data_13962_2017_04_13.mat');
% Parts (c)
Dfig2c   = load('C:\Users\nfc\Documents\ProtoMPEX_Ops\2017\2017_06\2017_06_16\IR_hollow_15083_2017_06_16');
% Parts (e)
Dfig2e = load('C:\Users\nfc\Documents\ProtoMPEX_Ops\2017\2017_06\2017_06_16\BdotData_CoreEdge_2017_06_16');

% #########################################################################
% Figure 3 data:
% #########################################################################
Dfig3 = load('C:\Users\nfc\Documents\ProtoMPEX_Ops\2017\2017_04\2017_04_13\Highly ionized plasma column\DLP_4_10_2017_04_13');

% #########################################################################
% Figure 4 data:
% #########################################################################
Dfig4 = load('C:\Users\nfc\Documents\ProtoMPEX_Ops\2017\2017_01\2017_01_18\BdotData_Bz_TR2_200A_2017_01_18');

%% Figure 2
close all

figure; % =================================================================
subplot(1,3,1)
f = Dfig2abd.dIR_t_4_2s;
contourf(f,20,'lineStyle','none');
caxis([0,1000])
axis('square')
text(20,190,'t = 4.2 s','Color','w','FontName','times','FontSize',11)
colormap(flipud(hot))
ylim([30,200])
xlim([10,190])
set(gca,'YTick',[],'XTick',[])
ylabel('y','FontName','times','FontSize',12)
xlabel('x','FontName','times','FontSize',12)

subplot(1,3,2)
f = Dfig2abd.dIR_t_4_43s;
contourf(f,20,'lineStyle','none');
caxis([0,1000])
axis('square')
text(20,190,'t = 4.43 s','Color','w','FontName','times','FontSize',11)
colormap(flipud(hot))
ylim([30,200])
xlim([10,190])
set(gca,'YTick',[],'XTick',[])
ylabel('y','FontName','times','FontSize',12)
xlabel('x','FontName','times','FontSize',12)

subplot(1,3,3);
f = Dfig2c.IR_15083_t_4_43;
contourf(f,20,'LineStyle','none');
caxis([0,3000])
axis('square')
text(20+10,190-10,'t = 4.43 s','Color','w','FontName','times','FontSize',11)
colormap(flipud(hot))
ylim([30,200]-10)
xlim([10,190]+10)
set(gca,'YTick',[],'XTick',[])
ylabel('y','FontName','times','FontSize',12)
xlabel('x','FontName','times','FontSize',12)
set(gcf,'Color','w')
box on

set(gcf,'color','w')
box on

%%
figure; % =================================================================
fvid = gcf;
set(fvid,'Menubar','figure','color','w','Units','normalized')

subplot(1,3,1); ax(1) = gca;
f = Dfig2abd.dIR_t_4_2s;
contourf(f,20,'lineStyle','none');
caxis([0,1000])
axis('square')
% text(20,190,'t = 4.2 s','Color','w','FontName','times','FontSize',11)
colormap(flipud(hot))
ylim([30,200])
xlim([10,190])
set(gca,'YTick',[],'XTick',[])

subplot(1,3,2); ax(2) = gca;
f = Dfig2abd.dIR_t_4_43s;
contourf(f,20,'lineStyle','none');
caxis([0,1000])
axis('square')
% text(20,190,'t = 4.43 s','Color','w','FontName','times','FontSize',11)
colormap(flipud(hot))
ylim([30,200])
xlim([10,190])
set(gca,'YTick',[],'XTick',[])

subplot(1,3,3); ax(3) = gca;
f = Dfig2c.IR_15083_t_4_43;
contourf(f,20,'LineStyle','none');
caxis([0,3000])
axis('square')
% text(20+10,190-10,'t = 4.43 s','Color','w','FontName','times','FontSize',11)
colormap(flipud(hot))
ylim([30,200]-10)
xlim([10,190]+10)
set(gca,'YTick',[],'XTick',[])
set(gcf,'Color','w')
box on

set(gcf,'color','w')
box on

dx = 0.15;
dy = 0.2;
h = (1-2*dy);
w = (1-2*dx)/3;
 
ax(1).Position = [dx   dy w h];
ax(2).Position = [dx+1.05*w dy w h];
ax(3).Position = [dx+2.10*w dy w h];
axis(ax,'square')

%%
figure; % =================================================================
f = Dfig2abd.RadialHeatFlux;
fx = Dfig2abd.t_rhf;
fy = Dfig2abd.r_rhf;
contourf(fx,fy,f,20,'lineStyle','none');
line([4.2,4.2],[0,400],'color','k','linewidth',2,'linestyle','--')
line([4.434,4.434],[0,400],'color','k','linewidth',2,'linestyle','--')
colormap(flipud(hot))
caxis([0,1000])
xlim([4.12,4.50])
ylim([25,200])
box on
xlabel('t [s]','Interpreter','latex','FontSize',12)
% ylabel('y [a.u]')
set(gca,'YTick',[])
set(gca,'FontName','times','FontSize',11)
set(gca,'PlotBoxAspectRatio',[1 0.3 1],'color','w')
set(gcf,'color','w')

% From Pawel 2017_08_10
% 0.3 mm per pixel

figure; % =================================================================
% Find edge of the plasma in the TG region
rngTime_TG = find(fx>=4.19 & fx<=4.21);
f_TG = f(:,rngTime_TG(1));
[nf_TG_peak,f_TG_peak] = peakseek(f_TG,75,0.5*max(f_TG));
dPixel = 0;

% Find Core of the plasma in the H region
rngTime_H = find(fx>=4.43);
f_H = f(:,rngTime_H(1));
[nf_H_peak,f_H_peak] = peakseek(f_H,75,0.5*max(f_H));

hold on
Intensity2PowerFlux = 0.625/f_H_peak; % MWm^-2/Intensity
H_flux = Intensity2PowerFlux*f(nf_H_peak(1)-dPixel:nf_H_peak(1)+dPixel,:);
TG_flux = Intensity2PowerFlux*f(nf_TG_peak(1)-dPixel:nf_TG_peak(1)+dPixel,:);

ha(1) = plot(fx,mean(H_flux,1),'ksq-'); set(ha(1),'LineWidth',1)
% h(2) = plot(fx,f(nf_TG_peak(1),:),'r.-'); set(h(2),'LineWidth',0.5)
ha(2) = plot(fx,mean(TG_flux,1),'ro-'); set(ha(2),'LineWidth',1)
xlim([4.12,4.50])
ylim(Intensity2PowerFlux*[0,1200])
box on
xlabel('t [s]','Interpreter','latex','FontSize',12)
ylabel('MW/m^{-2}')
% set(gca,'YTick',[])
set(gca,'FontName','times','FontSize',11)
set(gca,'PlotBoxAspectRatio',[1 0.3 1],'color','w')
set(gcf,'color','w')
L = legend(ha,'Core','Edge','location','NorthWest'); L.Box = 'off';
grid on

figure; % =================================================================
hold on
h(1) = plot(Dfig2e.tcore,Dfig2e.Acore,'k'); % x = 6 from the edge
h(2) = plot(Dfig2e.tedge,Dfig2e.Aedge,'r') % x = 11 from the edge

text(4.35,1.0,'Core','Color','k','FontSize',18)
text(4.35,0.35,'Edge','Color','r','FontSize',18)
% xlabel('t [s]')
ylabel('$|B_r|$ [A.U]','Interpreter','latex')

ylim([0,1.5])
xlim([4.12,4.5])
set(gcf,'color','w')
set(gca,'PlotBoxAspectRatio',[1 0.3 1],'color','w')
set(gca,'YTickLabel',[],'FontName','Times','FontSize',12)
set(gca,'XTickLabel',[],'FontName','Times','FontSize',12)
box on

hold on; box on;
%%plot stuff here%%%
set(h(1),'linewidth',3)
set(h(2),'linewidth',3)

set(gca,'fontsize',18,'fontname','times','fontweight','bold','linewidth',1.5)
print(gcf,'C:\Users\nfc\Documents\Proto-MPEX Data Analysis\2017\2017_04\2017_04_13\TG supression stuff\figure3a','-dpng')

%% Figure 3 
close all
figure
hold on
yyaxis left
h(1) = plot(Dfig3.tne4,Dfig3.ne4*1e-19,'k'); set(h(1),'LineWidth',2);
h(2) = plot(Dfig3.tne10,Dfig3.ne10*1e-19,'bl'); set(h(2),'LineWidth',2);
% h(3) = plot(Dfig3.trfnorm,Dfig3.RFnormalized*6.5,'r:'); set(h(3),'LineWidth',1);
xlim([4.12,4.5])
ylim([0,7])
ax1 = gca;

set(gcf,'color','w')
set(ax1,'PlotBoxAspectRatio',[1 0.5 1],'color','w')
set(ax1,'FontName','Times','FontSize',12,'Xcolor','k','Ycolor','k')% ,'YTickLabel',[]
ylabel(ax1,'$n_e$ [$m^{-3}$]$\times$ $10^{19}$','Interpreter','latex','FontSize',17,'color','k','fontweight','bold')

yyaxis right
nskip = 11;
% h(3) = plot(Dfig3.tTe4(nskip:end-2),Dfig3.Te4(nskip:end-2),'k:'); set(h(3),'LineWidth',2);
% h(4) = plot(Dfig3.tTe10(nskip:end-2),Dfig3.Te10(nskip:end-2),'bl:'); set(h(4),'LineWidth',2);
h(3) = plot(Dfig3.tTe4,Dfig3.Te4,'k:'); set(h(3),'LineWidth',2);
h(4) = plot(Dfig3.tTe10,Dfig3.Te10,'bl:'); set(h(4),'LineWidth',2);
xlim([4.12,4.5])
ylim([0,7])
ax2 = gca;
ylabel(ax2,'$T_e$ $[eV]$','Interpreter','latex','FontSize',17,'color','k','fontweight','bold')
set(ax2,'Xcolor','k','Ycolor','k')


xlabel('$t [s]$','Interpreter','Latex','Fontsize',17)
box on
set(gca,'fontsize',13,'fontname','times','fontweight','bold','linewidth',1.5)

legend(ax1,[h(1),h(3)],'n_e','T_e','location','NorthWest','FontSize',18); legend('boxoff')

text(ax1,4.36,6,'Probe A','fontsize',13,'fontname','times','fontweight','bold','color','k')
text(ax1,4.36,1.2,'Probe B','fontsize',13,'fontname','times','fontweight','bold','color','bl')

print(figure(1),'C:\Users\nfc\Documents\Proto-MPEX Data Analysis\2017\2017_04\2017_04_13\TG supression stuff\figure3b','-dpng')
%% Figure 4
close all
figure; 
hold on
Rad = Dfig4.Rad;
BzAmp     = Dfig4.BzAmp;
dBzAmp    = Dfig4.dBzAmp;
ArgBzAmp  = Dfig4.ArgBzAmp;
dArgBzAmp = Dfig4.dArgBzAmp;

h4 = errorbar(Rad,BzAmp,dBzAmp,'ksq-'); h4.MarkerSize = 6;
ylim([0,4])
xlim([-6,6])
set(gca,'PlotBoxAspectRatio',[1.5 1 1])
xlabel('r [cm]','FontName','Times')
ylabel('|B_z| [A.U]','FontName','Times')
box on
set(gca,'FontName','Times','FontSize',11)
set(gcf,'Position',[450 330 470 280],'color','w')


figure; 
subplot(2,1,1)
h4 = errorbar(Rad,BzAmp,dBzAmp,'ksq-'); h4.MarkerSize = 6;
box on
ylim([0,4])
xlim([-6,6])
ylabel('|B_z| [A.U]','FontName','Times')

subplot(2,1,2); hold on
h4 = errorbar(Rad(1:5),ArgBzAmp(1:5)/180,dArgBzAmp(1:5)/180,'ksq-'); h4.MarkerSize = 7;
h4 = errorbar(Rad(6:end),ArgBzAmp(6:end)/180,dArgBzAmp(6:end)/180,'ksq-'); h4.MarkerSize = 7;
xlim([-6,6])
ylim([-3,0])
ylabel('\phi /\pi [Rad]','FontName','Times')
xlabel('r [cm]','FontName','Times')
box on
set(gcf,'color','w','Position',[763.6667  281.0000  370.6667  292.0000])

%% Figure 3 for Highly-ionized plasma.. paper
% close all
ab_x = -15;
ab_y = 195;
figure; % =================================================================
subplot(1,2,1)
f = Dfig2abd.dIR_t_4_2s;
contourf(f,20,'lineStyle','none');
caxis([0,1000])
axis('square')
text(ab_x,ab_y,'a)','Color','k','FontName','times','FontSize',14)
text(20,190,'t = 4.2 s','Color','w','FontName','times','FontSize',11)
text(10,215,'Start of RF pulse','Color','k','FontName','times','FontSize',12)

colormap(flipud(hot))
ylim([30,200])
xlim([10,190])
set(gca,'YTick',[],'XTick',[])
ylabel('y','FontName','times','FontSize',12)
xlabel('x','FontName','times','FontSize',12)

subplot(1,2,2)
f = Dfig2abd.dIR_t_4_43s;
contourf(f,20,'lineStyle','none');
caxis([0,1000])
axis('square')
text(ab_x,ab_y,'b)','Color','k','FontName','times','FontSize',14)
text(20,190,'t = 4.43 s','Color','w','FontName','times','FontSize',11)
text(10,215,'End of RF pulse','Color','k','FontName','times','FontSize',12)
colormap(flipud(hot))
ylim([30,200])
xlim([10,190])
set(gca,'YTick',[],'XTick',[])
ylabel('y','FontName','times','FontSize',12)
xlabel('x','FontName','times','FontSize',12)

set(gcf,'Position',[ 360.3333  197.6667  448.0000  336.0000],'color','w')

