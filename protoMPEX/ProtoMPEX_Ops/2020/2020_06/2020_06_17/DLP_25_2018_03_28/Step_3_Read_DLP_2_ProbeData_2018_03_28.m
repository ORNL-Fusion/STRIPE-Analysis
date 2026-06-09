% Read data:
clear all
close all

load('DLP_2_ProbeData_2018_03_28.mat')


% Plot time-resolved radial data:
% =========================================================================
figure('color','w');
hold on
surf(DLP2.RR,DLP2.TT,DLP2.ne,'LineStyle','none')
hL = line([0,7.5],[1,1]*DLP2.t_rfStart,[1,1]*1e18,'color','r','LineWidth',2);
zlim([-1,10]*1e19)
xlim([0,7.5])
ylim([4.15,4.18])
set(gca,'FontName','Times','FontSize',11)
xlabel('r [cm]','Interpreter','latex','FontSize',14)
ylabel('time [s]','Interpreter','latex','FontSize',14)
zlabel('$n_e$ [m$^{-3}$]','Interpreter','latex','FontSize',14)
title('DLP 2.5 time-resolved radial scan, 2018-03-28','Interpreter','latex','FontSize',14)
view([30,60])
colorbar
hLeg = legend(hL,'Start of RF');
hLeg.Location = 'best'

saveFig = 1;

if saveFig
    saveas(gcf,'Step_3_TimeResolvedRadialScan_DLP_2','tiffn')
end

figure('color','w');
hold on
surf(sqrt(DLP2.Xi),DLP2.TT,DLP2.ne,'LineStyle','none')
hL = line([0,7.5],[1,1]*DLP2.t_rfStart,[1,1]*1e18,'color','r','LineWidth',2);
hXi = line([1,1],[4.15,4.2],[1,1]*1e18,'Color','r','LineWidth',2,'LineStyle','--');
zlim([-1,10]*1e19)
xlim([0,1.3])
ylim([4.15,4.18])
set(gca,'FontName','Times','FontSize',11)
xlabel('$\sqrt{\chi}$','interpreter','latex','Fontsize',12)
ylabel('time [s]','Interpreter','latex','FontSize',14)
zlabel('$n_e$ [m$^{-3}$]','Interpreter','latex','FontSize',14)
title('DLP 2.5 time-resolved radial scan, 2018-03-28','Interpreter','latex','FontSize',14)
view([30,60])
colorbar
hLeg = legend(hL,'Start of RF');
hLeg.Location = 'best';

saveFig = 1;

if saveFig
    saveas(gcf,'Step_3_TimeResolvedRadialScan_DLP_2_Xi','tiffn')
end